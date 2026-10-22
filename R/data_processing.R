# ==============================================================================
# TargetCNV-ML: Libreria de Procesamiento Espacial y Limpieza de Datos (R/data_processing.R)
# ==============================================================================

#' Limpia y estandariza identificadores de secuencias de CAPTUS
#'
#' Extrae el identificador limpio del ortogrupo (OG) y construye un ID compatible
#' con la nomenclatura de contigs ensamblados en el GFF de CAPTUS (OGxxxx__00).
#' Respeta la logica exacta de Workflow_Central_Sin_Rescate.R.
#'
#' @param df Data frame con las estadisticas crudas de CAPTUS (captus-extract_stats.tsv).
#' @param ref_col Nombre de la columna que contiene la referencia u ortogrupo.
#' @param hit_col Nombre de la columna que contiene el indice numerico del hit.
#' @return Data frame con las columnas adicionales `OG_Limpio` e `ID_GFF_Compatible`.
#' @export
clean_captus_ids <- function(df, ref_col = "ref_name", hit_col = "hit") {
  stopifnot(is.data.frame(df))
  stopifnot(ref_col %in% colnames(df))
  stopifnot(hit_col %in% colnames(df))

  df %>%
    dplyr::mutate(
      OG_Limpio = dplyr::coalesce(
        stringr::str_extract(as.character(.data[[ref_col]]), "OG[0-9]+"),
        as.character(.data[[ref_col]])
      ),
      ID_GFF_Compatible = sprintf("%s__%02d", OG_Limpio, as.numeric(.data[[hit_col]]))
    )
}

#' Filtra caracteristicas especificas de un objeto GFF/GRanges y limpia IDs de genes
#'
#' @param gff Objeto GRanges importado con rtracklayer::import().
#' @param feature_type Tipo de anotacion a conservar (ej. "gene" o "protein_match:NUC").
#' @return Objeto GRanges filtrado con columna `Locus_Ref_Limpio` normalizada.
#' @export
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

#' Limpia prefijos de scaffolds/contigs entre GFF de CAPTUS y GFF de referencia
#'
#' Util cuando CAPTUS agrega prefijos de especie o "contig_" a seqlevels.
#' Solo altera seqlevels si no hay solape directo entre nombres de cromosomas.
#'
#' @param query_gff Objeto GRanges de CAPTUS.
#' @param subject_gff Objeto GRanges de referencia.
#' @param prefix_regex Expresion regular de prefijos a remover.
#' @return Lista con `query_gff` y `subject_gff` armonizados.
#' @export
harmonize_gff_seqlevels <- function(query_gff,
                                    subject_gff,
                                    prefix_regex = "^[Ss]citrullifolium_|^[Ss]torvum_|^contig_") {
  common_levels <- intersect(GenomeInfoDb::seqlevels(query_gff), GenomeInfoDb::seqlevels(subject_gff))
  if (length(common_levels) == 0) {
    GenomeInfoDb::seqlevels(query_gff) <- stringr::str_remove_all(
      GenomeInfoDb::seqlevels(query_gff), prefix_regex
    )
    GenomeInfoDb::seqlevels(subject_gff) <- stringr::str_remove_all(
      GenomeInfoDb::seqlevels(subject_gff), prefix_regex
    )
  }
  list(query_gff = query_gff, subject_gff = subject_gff)
}

#' Ejecuta la interseccion espacial entre contigs de CAPTUS y genes de referencia
#'
#' @param query_gff Objeto GRanges con marcadores funcionales de CAPTUS (protein_match:NUC).
#' @param subject_gff Objeto GRanges con los genes del genoma de referencia.
#' @param ignore_strand Booleano; FALSE exige coincidencia en la misma hebra (criterio estricto).
#' @return Data frame con pares de solape (`Captus_Row`, `Ref_Row`, `ID_Captus`, `OG_Captus`, `Locus_Ref_Limpio`).
#' @export
find_spatial_overlaps <- function(query_gff, subject_gff, ignore_strand = FALSE) {
  stopifnot(inherits(query_gff, "GRanges"))
  stopifnot(inherits(subject_gff, "GRanges"))

  harmonized <- harmonize_gff_seqlevels(query_gff, subject_gff)
  query_gff <- harmonized$query_gff
  subject_gff <- harmonized$subject_gff

  hits <- GenomicRanges::findOverlaps(
    query = query_gff,
    subject = subject_gff,
    type = "any",
    ignore.strand = ignore_strand
  )

  if (length(hits) == 0) {
    return(data.frame(
      Captus_Row = integer(0),
      Ref_Row = integer(0),
      ID_Captus = character(0),
      OG_Captus = character(0),
      Locus_Ref_Limpio = character(0),
      stringsAsFactors = FALSE
    ))
  }

  df_cruces <- data.frame(
    Captus_Row = S4Vectors::queryHits(hits),
    Ref_Row    = S4Vectors::subjectHits(hits),
    stringsAsFactors = FALSE
  )

  df_cruces$ID_Captus <- as.character(query_gff$Name[df_cruces$Captus_Row])
  df_cruces$OG_Captus <- dplyr::coalesce(
    stringr::str_extract(df_cruces$ID_Captus, "OG[0-9]+"),
    stringr::str_remove(df_cruces$ID_Captus, "__[0-9]+$")
  )
  df_cruces$Locus_Ref_Limpio <- as.character(subject_gff$Locus_Ref_Limpio[df_cruces$Ref_Row])

  return(df_cruces)
}

#' Extrae el locus de referencia directamente desde el atributo Target del GFF de CAPTUS
#'
#' Diseñado para genomas altamente fragmentados (ej. S. stellatiglandulosum) donde
#' el atributo Target codifica `<transcript_id>-<OG_ID> <coordenadas>`.
#'
#' @param query_gff Objeto GRanges de CAPTUS (protein_match:NUC).
#' @return Data frame con `ID_Captus`, `OG_Captus`, `Locus_Ref_Limpio` y `Contig_Name`.
#' @export
extract_target_loci_from_gff <- function(query_gff) {
  stopifnot(inherits(query_gff, "GRanges"))
  if (!"Target" %in% colnames(GenomicRanges::mcols(query_gff))) {
    stop("[ERROR] El GFF de CAPTUS no contiene el atributo 'Target' requerido para este modo.", call. = FALSE)
  }

  data.frame(
    ID_Captus = as.character(query_gff$Name),
    OG_Captus = dplyr::coalesce(
      stringr::str_extract(as.character(query_gff$Name), "OG[0-9]+"),
      stringr::str_remove(as.character(query_gff$Name), "__[0-9]+$")
    ),
    Locus_Ref_Limpio = stringr::str_remove(
      stringr::str_extract(as.character(query_gff$Target), "^\\S+"),
      "-OG[0-9]+$"
    ),
    Contig_Name = as.character(GenomicRanges::seqnames(query_gff)),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      Locus_Ref_Limpio = stringr::str_remove(Locus_Ref_Limpio, "\\.[0-9]+$")
    ) %>%
    dplyr::distinct()
}

#' Construye el diccionario de homologia teorica a partir del archivo Excel de referencia
#'
#' Reproduce con exactitud matematica el filtrado de vias metabolicas y ensamblajes
#' de `Workflow_Central_Sin_Rescate.R` cuando se proporciona `Genes_Capture_May_2026.xlsx`,
#' o carga directamente cualquier tabla `.xlsx`, `.tsv` o `.csv` personalizada.
#'
#' @param dict_path Ruta al archivo de diccionario (.xlsx, .tsv, .csv).
#' @param species Especie/muestra objetivo ("Storvum", "Scitrullifolium", "S_stella" o NULL).
#' @return Lista con `dict_df` (columnas `Locus_Benoit`, `OG_Benoit`) y opcionalmente `hog_bridge`.
#' @export
load_reference_dictionary <- function(dict_path, species = NULL) {
  stopifnot(file.exists(dict_path))

  pathways_filtro <- c(
    "Phenylpropanoid", "SGAs", "Anthocyanin", "Flavonoid", "Terpenoids",
    "ASs", "UGT", "Metabolites", "Trichome structure",
    "Glandular Trichome, Leaf development",
    "Trichome structure: Trichome differentiation",
    "Trichome structure: Glandular",
    "Trichome structure: Glandular type VI", "Trichome_structure"
  )

  if (grepl("\\.xlsx$", dict_path, ignore.case = TRUE)) {
    suppressPackageStartupMessages(library(readxl))
    sheets <- readxl::excel_sheets(dict_path)

    # Caso 1: S. stellatiglandulosum (puente HOG -> OG desde CNV_Perez)
    if (!is.null(species) && grepl("stella", species, ignore.case = TRUE) &&
        all(c("CNV_Perez", "Contigs_probe_design", "Sequences_Torva") %in% sheets)) {
      columna_ste <- "S._stellatiglandulosum"
      cnv_ste_raw <- readxl::read_excel(dict_path, sheet = "CNV_Perez") %>%
        dplyr::select(HOG_ID, Probe_ID, dplyr::all_of(columna_ste)) %>%
        dplyr::filter(!is.na(.data[[columna_ste]]) & trimws(.data[[columna_ste]]) != "")
      probes_validos <- cnv_ste_raw$Probe_ID

      contigs_probe <- readxl::read_excel(dict_path, sheet = "Contigs_probe_design") %>%
        dplyr::select(Probe_ID, Contig_name, HOG_ID_Solanum, OG_ID_Pangenomes,
                      Correcte_Contig_Name, Short_name, Pathway, Refseq_RNA, DAVID_Annotation) %>%
        dplyr::filter(Probe_ID %in% probes_validos)

      sequences_ste <- readxl::read_excel(dict_path, sheet = "Sequences_Torva") %>%
        dplyr::filter(stringr::str_detect(Assembly, "stellati"), Probe_ID %in% probes_validos)

      puente_og_hog <- readxl::read_excel(dict_path, sheet = "Contigs_probe_design") %>%
        dplyr::select(OG_ID = OG_ID_Pangenomes, HOG_ID = HOG_ID_Solanum) %>%
        dplyr::filter(!is.na(OG_ID), !is.na(HOG_ID)) %>%
        dplyr::distinct() %>%
        dplyr::add_count(OG_ID, name = "n_hog") %>%
        dplyr::group_by(OG_ID) %>%
        dplyr::slice_max(order_by = n_hog, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::select(OG_ID, HOG_ID)

      dict_ste <- cnv_ste_raw %>%
        dplyr::left_join(contigs_probe, by = "Probe_ID") %>%
        dplyr::left_join(sequences_ste %>% dplyr::select(-Orthogroup_ID), by = "Probe_ID") %>%
        dplyr::filter(trimws(Pathway) %in% pathways_filtro) %>%
        dplyr::left_join(sequences_ste %>% dplyr::select(Seq, Orthogroup_ID), by = "Seq") %>%
        dplyr::filter(!is.na(Contig_ID)) %>%
        dplyr::select(Locus_Benoit = Contig_ID, HOG_Esperado = Orthogroup_ID) %>%
        dplyr::mutate(Locus_Benoit = stringr::str_remove(Locus_Benoit, "\\.[0-9]+$")) %>%
        dplyr::left_join(puente_og_hog, by = c("HOG_Esperado" = "HOG_ID"), relationship = "many-to-many") %>%
        dplyr::transmute(Locus_Benoit = Locus_Benoit, OG_Benoit = dplyr::coalesce(OG_ID, HOG_Esperado)) %>%
        dplyr::distinct()

      return(dict_ste)
    }

    # Caso 2: S. citrullifolium (hoja S._citr + CNV_Benoit)
    if (!is.null(species) && grepl("citr", species, ignore.case = TRUE) &&
        all(c("CNV_Benoit", "Contigs_probe_design", "S._citr") %in% sheets)) {
      cnv_cit <- readxl::read_excel(dict_path, sheet = "CNV_Benoit") %>%
        dplyr::select(OG_ID, Probe_ID, `S._citrullifolium`) %>%
        dplyr::filter(!is.na(`S._citrullifolium`) & trimws(`S._citrullifolium`) != "")
      probes_validos <- cnv_cit$Probe_ID

      contigs_probe <- readxl::read_excel(dict_path, sheet = "Contigs_probe_design") %>%
        dplyr::select(Probe_ID, Contig_name, HOG_ID_Solanum, OG_ID_Pangenomes,
                      Correcte_Contig_Name, Short_name, Pathway, Refseq_RNA, DAVID_Annotation) %>%
        dplyr::filter(Probe_ID %in% probes_validos)

      sequences_cit <- readxl::read_excel(dict_path, sheet = "S._citr") %>%
        dplyr::filter(Probe_ID %in% probes_validos)

      dict_cit <- cnv_cit %>%
        dplyr::select(-OG_ID) %>%
        dplyr::left_join(contigs_probe, by = "Probe_ID") %>%
        dplyr::left_join(sequences_cit %>% dplyr::select(-Orthogroup_ID), by = "Probe_ID") %>%
        dplyr::mutate(Assembly = "S. citrullifolium Benoit et al") %>%
        dplyr::filter(trimws(Pathway) %in% pathways_filtro,
                      Assembly == "S. citrullifolium Benoit et al") %>%
        dplyr::select(-`S._citrullifolium`) %>%
        dplyr::left_join(sequences_cit %>% dplyr::select(Seq, Orthogroup_ID), by = "Seq") %>%
        dplyr::filter(!is.na(Contig_ID)) %>%
        dplyr::select(Locus_Benoit = Contig_ID, OG_Benoit = Orthogroup_ID) %>%
        dplyr::mutate(Locus_Benoit = stringr::str_remove(Locus_Benoit, "\\.[0-9]+$")) %>%
        dplyr::distinct()

      return(dict_cit)
    }

    # Caso 3: S. torvum (experimento central en Workflow_Central_Sin_Rescate.R)
    if (all(c("CNV_Benoit", "Contigs_probe_design", "Sequences_Torva") %in% sheets)) {
      cnv_benoit <- readxl::read_excel(dict_path, sheet = "CNV_Benoit") %>%
        dplyr::select(OG_ID, Probe_ID, dplyr::any_of(c("S._torvum", "S_torvum")))
      col_torvum <- intersect(c("S._torvum", "S_torvum"), colnames(cnv_benoit))[1]
      if (!is.null(col_torvum) && !is.na(col_torvum)) {
        cnv_benoit <- cnv_benoit %>%
          dplyr::filter(!is.na(.data[[col_torvum]]) & trimws(as.character(.data[[col_torvum]])) != "")
      }
      probes_validos <- cnv_benoit$Probe_ID

      contigs_probe <- readxl::read_excel(dict_path, sheet = "Contigs_probe_design") %>%
        dplyr::select(Probe_ID, Contig_name, HOG_ID_Solanum, OG_ID_Pangenomes,
                      Correcte_Contig_Name, Short_name, Pathway, Refseq_RNA, DAVID_Annotation) %>%
        dplyr::filter(Probe_ID %in% probes_validos)

      sequences_torva <- readxl::read_excel(dict_path, sheet = "Sequences_Torva") %>%
        dplyr::filter(Probe_ID %in% probes_validos)

      datos_unificados <- cnv_benoit %>%
        dplyr::select(-OG_ID) %>%
        dplyr::left_join(contigs_probe, by = "Probe_ID") %>%
        dplyr::left_join(sequences_torva %>% dplyr::select(-Orthogroup_ID), by = "Probe_ID")

      dict_benoit <- datos_unificados %>%
        dplyr::filter(trimws(Pathway) %in% pathways_filtro,
                      Assembly == "S. torvum Benoit et al") %>%
        dplyr::select(-dplyr::any_of(c("S._torvum", "S_torvum"))) %>%
        dplyr::left_join(sequences_torva %>% dplyr::select(Seq, Orthogroup_ID), by = "Seq") %>%
        dplyr::filter(!is.na(Contig_ID)) %>%
        dplyr::select(Locus_Benoit = Contig_ID, OG_Benoit = Orthogroup_ID) %>%
        dplyr::mutate(Locus_Benoit = stringr::str_remove(Locus_Benoit, "\\.[0-9]+$")) %>%
        dplyr::distinct()

      return(dict_benoit)
    }

    # Caso 4: Excel con 2 hojas o tabla plana
    if (all(c("CNV_Benoit", "Sequences_Torva") %in% sheets)) {
      cnv_b <- readxl::read_excel(dict_path, sheet = "CNV_Benoit") %>%
        dplyr::select(OG_ID, Probe_ID, dplyr::any_of(c("S._torvum", "S_torvum")))
      seq_t <- readxl::read_excel(dict_path, sheet = "Sequences_Torva")
      dict_ref <- cnv_b %>%
        dplyr::left_join(seq_t, by = "Probe_ID", relationship = "many-to-many") %>%
        dplyr::filter(!is.na(Contig_ID)) %>%
        dplyr::select(Locus_Benoit = Contig_ID, OG_Benoit = Orthogroup_ID) %>%
        dplyr::mutate(Locus_Benoit = stringr::str_remove(Locus_Benoit, "\\.[0-9]+$")) %>%
        dplyr::distinct()
      return(dict_ref)
    }

    return(as.data.frame(readxl::read_excel(dict_path)))
  }

  # Archivo TSV o CSV
  sep_char <- if (grepl("\\.csv$", dict_path, ignore.case = TRUE)) "," else "\t"
  utils::read.delim(dict_path, sep = sep_char, stringsAsFactors = FALSE)
}

#' Clasifica los hits espaciales contrastando contra el diccionario de ortologia
#'
#' @param df_cruces Data frame obtenido con `find_spatial_overlaps()`.
#' @param dict_ref Data frame del diccionario con identificadores de locus y ortogrupo.
#' @param locus_col Columna de locus en `df_cruces`.
#' @param dict_locus_col Columna de locus en `dict_ref`.
#' @param dict_og_col Columna de ortogrupo en `dict_ref`.
#' @return Data frame con `ID_Captus`, `Validacion_Espacial`, `Doble_Control_Biologico` y `Locus_Ref_Limpio`.
#' @export
classify_spatial_hits <- function(df_cruces,
                                  dict_ref,
                                  locus_col = "Locus_Ref_Limpio",
                                  dict_locus_col = "Locus_Benoit",
                                  dict_og_col = "OG_Benoit") {
  stopifnot(is.data.frame(df_cruces))
  stopifnot(is.data.frame(dict_ref))

  if (nrow(df_cruces) == 0) {
    return(data.frame(
      ID_Captus = character(0),
      Validacion_Espacial = character(0),
      Doble_Control_Biologico = character(0),
      Locus_Ref_Limpio = character(0),
      stringsAsFactors = FALSE
    ))
  }

  df_clasificacion <- df_cruces %>%
    dplyr::mutate(
      ID_Captus = dplyr::if_else(
        stringr::str_detect(ID_Captus, "__[0-9]+$"),
        ID_Captus,
        paste0(ID_Captus, "__00")
      )
    ) %>%
    dplyr::left_join(dict_ref, by = stats::setNames(dict_locus_col, locus_col), relationship = "many-to-many") %>%
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
    dplyr::select(
      ID_Captus,
      Validacion_Espacial,
      Doble_Control_Biologico,
      dplyr::all_of(locus_col)
    )

  return(df_clasificacion)
}

#' Resuelve redundancias tecnicas intra-locus degradando hits secundarios a FALSE_POSITIVE
#'
#' Ordena por `ref_len_matched` y `pct_recovered` de mayor a menor dentro de cada
#' locus real validado, conservando el hit principal como TRUE_POSITIVE y
#' reclasificando los hits redundantes como `5_TP_Redundante_Degradado`.
#'
#' @param df Data frame integrado con estadisticas de CAPTUS y validacion espacial.
#' @param locus_col Nombre de la columna con el locus del genoma de referencia.
#' @param class_col Nombre de la columna con la etiqueta `Validacion_Espacial`.
#' @param sample_col Nombre de la columna de muestra (si existe).
#' @return Data frame sin redundancias tecnicas.
#' @export
resolve_redundancies <- function(df,
                                 locus_col = "Locus_Ref_Limpio",
                                 class_col = "Validacion_Espacial",
                                 sample_col = if ("sample_name" %in% colnames(df)) "sample_name" else NULL) {
  stopifnot(is.data.frame(df))

  group_vars <- c(sample_col, locus_col, class_col)

  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::arrange(
      dplyr::desc(as.numeric(ref_len_matched)),
      dplyr::desc(as.numeric(pct_recovered)),
      .by_group = TRUE
    ) %>%
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
    dplyr::select(-dplyr::any_of("Es_Redundante"))
}

#' Calcula diferenciales locales (deltas) intra-locus respecto al hit dominante
#'
#' @param df Data frame con metricas continuas de CAPTUS.
#' @param locus_col Columna que identifica el locus u ortogrupo.
#' @param score_col Columna de puntaje total (`score`).
#' @param wscore_col Columna de puntaje ponderado (`wscore`).
#' @param identity_col Columna de porcentaje de identidad (`pct_identity`).
#' @param sample_col Columna de muestra (si aplica).
#' @return Data frame con las columnas `delta_wscore` y `delta_pct_identity`.
#' @export
calculate_locus_deltas <- function(df,
                                   locus_col = "locus",
                                   score_col = "score",
                                   wscore_col = "wscore",
                                   identity_col = "pct_identity",
                                   sample_col = if ("sample_name" %in% colnames(df)) "sample_name" else NULL) {
  stopifnot(is.data.frame(df))

  group_vars <- c(sample_col, locus_col)

  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::arrange(dplyr::desc(.data[[score_col]]), .by_group = TRUE) %>%
    dplyr::mutate(
      delta_wscore = max(.data[[wscore_col]], na.rm = TRUE) - .data[[wscore_col]],
      delta_pct_identity = max(.data[[identity_col]], na.rm = TRUE) - .data[[identity_col]]
    ) %>%
    dplyr::ungroup()
}

#' Extrae la Matriz de Variacion en el Numero de Copias (CNV) con control de fragmentacion
#'
#' Aplica las reglas biologicas de cobertura individual (`cov_threshold`) y
#' cobertura acumulada (`frag_threshold`) sobre los homologos validados por consenso.
#'
#' @param df Data frame con los veredictos del consenso de 7 modelos.
#' @param sample_col Nombre de la columna de muestra.
#' @param locus_col Nombre de la columna de locus/ortogrupo.
#' @param cov_threshold Cobertura individual minima (%) para confirmar un paralogo genuino (default: 50.0).
#' @param frag_threshold Cobertura acumulada maxima (%) para diagnosticar fragmentacion de ensamblaje (default: 120.0).
#' @param verdict_col Columna con el veredicto de consenso.
#' @return Lista con `cnv_matrix`, `cnv_annotated` y `functional_genes`.
#' @export
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
      dplyr::filter(.data[[verdict_col]] %in% c("Homologo_Valido", "Ortologo_Valido"))
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
    dplyr::filter(Clasificacion_Copia %in% c("Copia_Principal", "Paralogo_Confirmado")) %>%
    dplyr::arrange(.data[[locus_col]])

  if (nrow(df_funcionales) == 0) {
    df_matriz_cnv <- data.frame(Locus_Ortogrupo = character(0), stringsAsFactors = FALSE)
    colnames(df_matriz_cnv)[1] <- locus_col
  } else {
    df_matriz_cnv <- df_funcionales %>%
      dplyr::group_by(.data[[sample_col]], .data[[locus_col]]) %>%
      dplyr::summarise(CNV_Detectado = dplyr::n(), .groups = "drop") %>%
      tidyr::pivot_wider(
        names_from = dplyr::all_of(sample_col),
        values_from = CNV_Detectado,
        values_fill = 0
      )
  }

  return(list(
    cnv_matrix = df_matriz_cnv,
    cnv_annotated = df_cnv_bruto,
    functional_genes = df_funcionales
  ))
}
