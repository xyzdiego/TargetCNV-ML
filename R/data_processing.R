#' Data Processing and Spatial Filtering Module for TargetCNV-ML
#'
#' Pure functions for data extraction, GFF coordinate intersections,
#' spatial rescue, redundancy resolution, and CNV matrix generation.
#'
#' @name data_processing
NULL

#' Limpiar y Estandarizar Identificadores de Captus
#'
#' Extrae el identificador del ortogrupo (OG) y reconstruye una llave universal
#' compatible con el GFF y las tablas de estadisticas de Captus.
#'
#' @param df Data frame con datos crudos de estadisticas de Captus.
#' @param ref_col Nombre de la columna que contiene la referencia (por defecto "ref_name").
#' @param hit_col Nombre de la columna que indica el indice del hit (por defecto "hit").
#'
#' @return Un \code{data.frame} con las columnas adicionales \code{OG_Limpio} e \code{ID_GFF_Compatible}.
#' @export
#' @examples
#' \dontrun{
#' df_limpio <- clean_captus_ids(stats_raw)
#' }
clean_captus_ids <- function(df, ref_col = "ref_name", hit_col = "hit") {
  stopifnot(is.data.frame(df))
  stopifnot(ref_col %in% colnames(df))
  stopifnot(hit_col %in% colnames(df))

  df %>%
    dplyr::mutate(
      OG_Limpio = stringr::str_extract(.data[[ref_col]], "OG[0-9]+"),
      ID_GFF_Compatible = sprintf("%s__%02d", OG_Limpio, as.numeric(.data[[hit_col]]))
    )
}

#' Filtrar Elementos de Anotacion GFF
#'
#' Filtra un objeto de tipo \code{GRanges} importado desde un GFF/GFF3 segun el tipo
#' de caracteristica y limpia el identificador de locus eliminando prefijos y versiones de transcripto.
#'
#' @param gff Objeto \code{GRanges} cargado mediante rtracklayer.
#' @param feature_type Cadena con el tipo de entidad a filtrar (ej. "gene", "protein_match:NUC").
#'
#' @return Un objeto \code{GRanges} filtrado con el metadato \code{Locus_Ref_Limpio}.
#' @export
#' @examples
#' \dontrun{
#' genes_ref <- filter_gff_features(gff_raw, feature_type = "gene")
#' }
filter_gff_features <- function(gff, feature_type = "gene") {
  stopifnot(inherits(gff, "GRanges"))

  filtrado <- gff[gff$type == feature_type]
  if ("ID" %in% colnames(GenomicRanges::mcols(filtrado))) {
    filtrado$Locus_Ref_Limpio <- stringr::str_remove_all(
      as.character(filtrado$ID), "^gene:|\\.[0-9]+$"
    )
  }
  return(filtrado)
}

#' Interseccion Espacial de Coordenadas Genomicas
#'
#' Realiza el cruce espacial estricto entre los hits ensamblados de Captus y las
#' anotaciones genomicas de referencia usando \code{GenomicRanges::findOverlaps}.
#'
#' @param query_gff Objeto \code{GRanges} correspondiente a los contigs de Captus (query).
#' @param subject_gff Objeto \code{GRanges} correspondiente a los genes del genoma de referencia (subject).
#' @param ignore_strand Booleano. Si es FALSE (por defecto), exige concordancia de hebra para maxima precision.
#'
#' @return Un \code{data.frame} con las relaciones espaciales detectadas (\code{Captus_Row}, \code{Ref_Row},
#'   \code{ID_Captus}, \code{OG_Captus}, \code{Locus_Ref_Limpio}).
#' @export
find_spatial_overlaps <- function(query_gff, subject_gff, ignore_strand = FALSE) {
  stopifnot(inherits(query_gff, "GRanges"))
  stopifnot(inherits(subject_gff, "GRanges"))

  hits <- GenomicRanges::findOverlaps(
    query = query_gff,
    subject = subject_gff,
    type = "any",
    ignore.strand = ignore_strand
  )

  df_cruces <- data.frame(
    Captus_Row = S4Vectors::queryHits(hits),
    Ref_Row    = S4Vectors::subjectHits(hits)
  )

  df_cruces$ID_Captus <- as.character(query_gff$Name[df_cruces$Captus_Row])
  df_cruces$OG_Captus <- stringr::str_extract(df_cruces$ID_Captus, "OG[0-9]+")
  df_cruces$Locus_Ref_Limpio <- as.character(subject_gff$Locus_Ref_Limpio[df_cruces$Ref_Row])

  return(df_cruces)
}

#' Clasificacion Espacial Preliminar y Doble Control Biologico
#'
#' Cruza las intersecciones espaciales con el diccionario teorico de homologia
#' para validar emparejamientos de loci y catalogar las secuencias como TP,
#' hit cruzado o posible expansion genomica.
#'
#' @param df_cruces Data frame con los cruces espaciales generado por \code{find_spatial_overlaps}.
#' @param dict_ref Data frame con el diccionario teorico de homologia (ej. Locus_Benoit y OG_Benoit).
#' @param locus_col Nombre de la columna de locus en df_cruces (por defecto "Locus_Ref_Limpio").
#' @param dict_locus_col Nombre de la columna de locus en dict_ref (por defecto "Locus_Benoit").
#' @param dict_og_col Nombre de la columna de ortogrupo en dict_ref (por defecto "OG_Benoit").
#'
#' @return Un \code{data.frame} deduplicado por \code{ID_Captus} con las clasificaciones:
#'   \code{Validacion_Espacial} ("TRUE_POSITIVE" / "FALSE_POSITIVE") y \code{Doble_Control_Biologico}.
#' @export
classify_spatial_hits <- function(df_cruces,
                                  dict_ref,
                                  locus_col = "Locus_Ref_Limpio",
                                  dict_locus_col = "Locus_Benoit",
                                  dict_og_col = "OG_Benoit") {
  stopifnot(is.data.frame(df_cruces))
  stopifnot(is.data.frame(dict_ref))

  df_clasificacion <- df_cruces %>%
    dplyr::mutate(
      ID_Captus = dplyr::if_else(stringr::str_detect(ID_Captus, "__[0-9]+$"),
                                 ID_Captus,
                                 paste0(ID_Captus, "__00"))
    ) %>%
    dplyr::left_join(dict_ref, by = stats::setNames(dict_locus_col, locus_col)) %>%
    dplyr::mutate(
      Validacion_Espacial = dplyr::case_when(
        !is.na(.data[[dict_og_col]]) & OG_Captus == .data[[dict_og_col]] ~ "TRUE_POSITIVE",
        TRUE ~ "FALSE_POSITIVE"
      ),
      Doble_Control_Biologico = dplyr::case_when(
        Validacion_Espacial == "TRUE_POSITIVE" ~ "1_Validado_Benoit_Exacto",
        !is.na(.data[[dict_og_col]]) & OG_Captus != .data[[dict_og_col]] ~ "2_Hit_Cruzado_OG_Vecino",
        is.na(.data[[dict_og_col]]) ~ "3_Expansion_Nueva_Genoma"
      ),
      prioridad = dplyr::if_else(Validacion_Espacial == "TRUE_POSITIVE", 0L, 1L)
    ) %>%
    dplyr::arrange(prioridad) %>%
    dplyr::distinct(ID_Captus, .keep_all = TRUE) %>%
    dplyr::select(ID_Captus, Validacion_Espacial, Doble_Control_Biologico,
                  dplyr::all_of(locus_col))

  return(df_clasificacion)
}

#' Algoritmo de Rescate Espacial por Homologo Tandem
#'
#' Evalua los hits marcados inicialmente como falsos positivos (FP) y rescata aquellos
#' que se ubican en el mismo cromosoma a una distancia sintenica menor o igual al umbral
#' de un Verdadero Positivo ancla para el mismo grupo homologo.
#'
#' @param df_base Data frame con coordenadas y clasificacion estricta inicial.
#' @param og_col Columna del identificador de ortogrupo (por defecto "OG_Captus").
#' @param chrom_col Columna que indica el cromosoma (por defecto "Cromosoma").
#' @param gene_num_col Columna con el indice numerico del gen (por defecto "Num_Gen").
#' @param class_col Columna con la clasificacion estricta (por defecto "Validacion_Estricta").
#' @param max_distance Distancia maxima en numero de genes para considerar un tandem (por defecto 50).
#'
#' @return Un \code{data.frame} actualizado con las columnas \code{Validacion_Espacial} y
#'   \code{Doble_Control_Biologico} tras el rescate espacial.
#' @export
rescue_tandem_duplicates <- function(df_base,
                                     og_col = "OG_Captus",
                                     chrom_col = "Cromosoma",
                                     gene_num_col = "Num_Gen",
                                     class_col = "Validacion_Estricta",
                                     max_distance = 50) {
  stopifnot(is.data.frame(df_base))

  df_rescatado <- df_base %>%
    dplyr::group_by(.data[[og_col]]) %>%
    dplyr::mutate(
      Tiene_Ancla_TP = any(.data[[class_col]] == "TRUE_POSITIVE"),
      Chroms_TP = list(.data[[chrom_col]][.data[[class_col]] == "TRUE_POSITIVE"]),
      Posiciones_TP = list(.data[[gene_num_col]][.data[[class_col]] == "TRUE_POSITIVE"])
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      Es_Tandem = dplyr::if_else(
        .data[[class_col]] == "FALSE_POSITIVE" &
          Tiene_Ancla_TP &
          (.data[[chrom_col]] %in% unlist(Chroms_TP)),
        any(abs(.data[[gene_num_col]] - unlist(Posiciones_TP)[unlist(Chroms_TP) == .data[[chrom_col]]]) <= max_distance),
        FALSE
      ),
      Validacion_Espacial = dplyr::case_when(
        .data[[class_col]] == "TRUE_POSITIVE" ~ "TRUE_POSITIVE",
        Es_Tandem ~ "TRUE_POSITIVE",
        TRUE ~ "FALSE_POSITIVE"
      ),
      Doble_Control_Biologico = dplyr::case_when(
        .data[[class_col]] == "TRUE_POSITIVE" ~ "1_Validado_Benoit_Exacto",
        Es_Tandem ~ "2_Rescatado_Paralogo_Tandem",
        TRUE ~ "3_No_Rescatado_o_Ruido"
      ),
      prioridad = dplyr::case_when(
        Doble_Control_Biologico == "1_Validado_Benoit_Exacto" ~ 0L,
        Doble_Control_Biologico == "2_Rescatado_Paralogo_Tandem" ~ 1L,
        TRUE ~ 2L
      )
    ) %>%
    dplyr::ungroup()

  return(df_rescatado)
}

#' Resolucion de Redundancias Espaciales y Degradacion a Falso Positivo
#'
#' Ordena las secuencias asignadas a un mismo locus genomico segun su calidad
#' de ensamblaje (\code{ref_len_matched}, \code{pct_recovered}) y degrada los hits
#' secundarios redundantes a Falso Positivo para evitar sesgos de conteo.
#'
#' @param df Data frame integrado con estadisticas y validacion espacial.
#' @param locus_col Nombre de la columna del locus de referencia (por defecto "Locus_Ref_Limpio").
#' @param class_col Nombre de la columna de validacion espacial (por defecto "Validacion_Espacial").
#'
#' @return Un \code{data.frame} donde solo el hit principal conserva el estado \code{TRUE_POSITIVE}.
#' @export
resolve_redundancies <- function(df,
                                 locus_col = "Locus_Ref_Limpio",
                                 class_col = "Validacion_Espacial") {
  stopifnot(is.data.frame(df))

  df %>%
    dplyr::group_by(.data[[locus_col]], .data[[class_col]]) %>%
    dplyr::arrange(dplyr::desc(as.numeric(ref_len_matched)),
                   dplyr::desc(as.numeric(pct_recovered)),
                   .by_group = TRUE) %>%
    dplyr::mutate(
      Es_Redundante = dplyr::if_else(
        .data[[class_col]] == "TRUE_POSITIVE" &
          !is.na(.data[[locus_col]]) &
          dplyr::row_number() > 1,
        TRUE,
        FALSE
      ),
      Validacion_Espacial = dplyr::if_else(Es_Redundante, "FALSE_POSITIVE", .data[[class_col]]),
      Doble_Control_Biologico = dplyr::if_else(
        Es_Redundante,
        "5_TP_Redundante_Degradado",
        dplyr::coalesce(.data[["Doble_Control_Biologico"]], "4_Ruido_Intergenico_o_Intron")
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-Es_Redundante)
}

#' Calcular Deltas Locales de Alineamiento por Locus
#'
#' Calcula la diferencia de calidad (\code{delta_wscore} y \code{delta_pct_identity})
#' de cada hit respecto al mejor hit obtenido para el mismo ortogrupo o locus.
#'
#' @param df Data frame con las metricas continuas de alineamiento.
#' @param locus_col Nombre de la columna que agrupa el locus/ortogrupo (por defecto "locus").
#' @param score_col Columna de puntaje (por defecto "score").
#' @param wscore_col Columna de puntaje ponderado (por defecto "wscore").
#' @param identity_col Columna de porcentaje de identidad (por defecto "pct_identity").
#'
#' @return Un \code{data.frame} con las nuevas columnas \code{delta_wscore} y \code{delta_pct_identity}.
#' @export
calculate_locus_deltas <- function(df,
                                   locus_col = "locus",
                                   score_col = "score",
                                   wscore_col = "wscore",
                                   identity_col = "pct_identity") {
  stopifnot(is.data.frame(df))

  df %>%
    dplyr::group_by(.data[[locus_col]]) %>%
    dplyr::arrange(dplyr::desc(.data[[score_col]]), .by_group = TRUE) %>%
    dplyr::mutate(
      delta_wscore = max(.data[[wscore_col]], na.rm = TRUE) - .data[[wscore_col]],
      delta_pct_identity = max(.data[[identity_col]], na.rm = TRUE) - .data[[identity_col]]
    ) %>%
    dplyr::ungroup()
}

#' Control de Fragmentacion Espacial y Extraccion de Matriz de CNV
#'
#' Aplica las reglas biologicas de control de fragmentacion sobre secuencias con veredicto
#' valido de consenso para clasificar copias genicas (Copia Principal, Copia Homologa Confirmada,
#' Fragmentacion de Ensamblaje o Copia Dudosa) y consolida la matriz final de CNV.
#'
#' @param df Data frame que contiene las secuencias candidatas y su veredicto de consenso.
#' @param sample_col Columna que identifica la muestra (por defecto "Sample").
#' @param locus_col Columna del locus u ortogrupo (por defecto "Locus_Ortogrupo").
#' @param cov_threshold Umbral minimo de cobertura recuperada individual para confirmar homologo (por defecto 50.0).
#' @param frag_threshold Umbral maximo de cobertura acumulada para atribuir a fragmentacion (por defecto 120.0).
#' @param verdict_col Columna con el veredicto de consenso (opcional, por defecto "Veredicto_Consenso").
#'
#' @return Una lista con:
#'   \item{cnv_matrix}{Matriz de recuento de copias (Locus x Muestras) en formato ancho.}
#'   \item{cnv_annotated}{Data frame anotado con la clasificacion de cada copia.}
#'   \item{functional_genes}{Data frame filtrado solo con copias funcionales validas.}
#' @export
#' @examples
#' \dontrun{
#' res_cnv <- extract_cnv_matrix(df_resultados_finales, cov_threshold = 50.0, frag_threshold = 120.0)
#' head(res_cnv$cnv_matrix)
#' }
extract_cnv_matrix <- function(df,
                               sample_col = "Sample",
                               locus_col = "Locus_Ortogrupo",
                               cov_threshold = 50.0,
                               frag_threshold = 120.0,
                               verdict_col = "Veredicto_Consenso") {
  stopifnot(is.data.frame(df))
  stopifnot(sample_col %in% colnames(df))
  stopifnot(locus_col %in% colnames(df))
  stopifnot("pct_recovered" %in% colnames(df))

  df_filtrado <- df
  if (verdict_col %in% colnames(df)) {
    df_filtrado <- df_filtrado %>%
      dplyr::filter(.data[[verdict_col]] == "Ortologo_Valido")
  }

  df_cnv_bruto <- df_filtrado %>%
    dplyr::group_by(.data[[sample_col]], .data[[locus_col]]) %>%
    dplyr::arrange(dplyr::desc(as.numeric(pct_recovered)), .by_group = TRUE) %>%
    dplyr::mutate(
      cobertura_acumulada = cumsum(as.numeric(pct_recovered)),
      Clasificacion_Copia = dplyr::case_when(
        dplyr::row_number() == 1 ~ "Copia_Principal",
        dplyr::row_number() > 1 & as.numeric(pct_recovered) >= cov_threshold ~ "Paralogo_Confirmado",
        dplyr::row_number() > 1 & as.numeric(pct_recovered) < cov_threshold & cobertura_acumulada <= frag_threshold ~ "Fragmentacion_de_Ensamblaje",
        TRUE ~ "Paralogo_Dudoso"
      )
    ) %>%
    dplyr::ungroup()

  df_funcionales <- df_cnv_bruto %>%
    dplyr::filter(Clasificacion_Copia %in% c("Copia_Principal", "Paralogo_Confirmado"))

  df_matriz_cnv <- df_funcionales %>%
    dplyr::group_by(.data[[sample_col]], .data[[locus_col]]) %>%
    dplyr::summarise(CNV_Detectado = dplyr::n(), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = dplyr::all_of(sample_col),
      values_from = CNV_Detectado,
      values_fill = 0
    )

  return(list(
    cnv_matrix = df_matriz_cnv,
    cnv_annotated = df_cnv_bruto,
    functional_genes = df_funcionales
  ))
}
