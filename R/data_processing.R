# ==============================================================================
# TargetCNV-ML: Libreria de Procesamiento Espacial y Limpieza de Datos (R/data_processing.R)
# ==============================================================================

#' Limpia y estandariza identificadores de secuencias de CAPTUS
#'
#' Extrae el identificador limpio del ortogrupo (OG) y construye un ID compatible
#' con la nomenclatura de contigs ensamblados en el GFF de CAPTUS (OGxxxx__00).
#' Respeta la logica exacta del flujo de referencia.
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
  tryCatch({
    q_levels <- as.character(unique(query_gff@seqnames))
    s_levels <- as.character(unique(subject_gff@seqnames))
    if (length(base::intersect(q_levels, s_levels)) == 0) {
      suppressPackageStartupMessages(library(GenomicRanges))
      seqlevels(query_gff) <- stringr::str_remove_all(seqlevels(query_gff), prefix_regex)
      seqlevels(subject_gff) <- stringr::str_remove_all(seqlevels(subject_gff), prefix_regex)
    }
  }, error = function(e) {
    warning(sprintf("[harmonize_gff_seqlevels] No se pudieron armonizar seqlevels via GenomicRanges (%s); se usara comparacion directa de nombres.", conditionMessage(e)), call. = FALSE)
  })
  list(query_gff = query_gff, subject_gff = subject_gff)
}

#' Ejecuta la interseccion espacial entre contigs de CAPTUS y genes de referencia
#'
#' Incluye respaldo directo sobre los slots S4 (@seqnames, @ranges, @strand) cuando el entorno
#' renv local carece de GenomeInfoDb (evitando el fallo de despacho S4 en intersect/findOverlaps).
#'
#' @param query_gff Objeto GRanges con marcadores funcionales de CAPTUS (protein_match:NUC).
#' @param subject_gff Objeto GRanges con los genes del genoma de referencia.
#' @param ignore_strand Booleano; FALSE exige coincidencia en la misma hebra (criterio estricto).
#' @return Data frame con pares de solape (`Captus_Row`, `Ref_Row`, `ID_Captus`, `OG_Captus`, `Locus_Ref_Limpio`).
#' @export
find_spatial_overlaps <- function(query_gff, subject_gff, ignore_strand = FALSE) {
  stopifnot(inherits(query_gff, "GRanges"))
  stopifnot(inherits(subject_gff, "GRanges"))

  prefix_regex <- "^[Ss]citrullifolium_|^[Ss]torvum_|^contig_"
  harmonized <- harmonize_gff_seqlevels(query_gff, subject_gff, prefix_regex = prefix_regex)
  query_gff <- harmonized$query_gff
  subject_gff <- harmonized$subject_gff

  df_cruces <- tryCatch({
    hits <- GenomicRanges::findOverlaps(
      query = query_gff,
      subject = subject_gff,
      type = "any",
      ignore.strand = ignore_strand
    )
    if (length(hits) == 0) {
      data.frame(Captus_Row = integer(0), Ref_Row = integer(0), stringsAsFactors = FALSE)
    } else {
      data.frame(
        Captus_Row = S4Vectors::queryHits(hits),
        Ref_Row    = S4Vectors::subjectHits(hits),
        stringsAsFactors = FALSE
      )
    }
  }, error = function(e) {
    warning(sprintf("[find_spatial_overlaps] GenomicRanges::findOverlaps fallo (%s); ejecutando cruce espacial directo sobre slots S4 (@seqnames, @ranges, @strand).", conditionMessage(e)), call. = FALSE)
    q_seq <- as.character(query_gff@seqnames)
    s_seq <- as.character(subject_gff@seqnames)
    if (length(base::intersect(unique(q_seq), unique(s_seq))) == 0) {
      q_seq <- stringr::str_remove_all(q_seq, prefix_regex)
      s_seq <- stringr::str_remove_all(s_seq, prefix_regex)
    }
    q_start  <- as.integer(query_gff@ranges@start)
    q_end    <- q_start + as.integer(query_gff@ranges@width) - 1L
    q_strand <- as.character(query_gff@strand)

    s_start  <- as.integer(subject_gff@ranges@start)
    s_end    <- s_start + as.integer(subject_gff@ranges@width) - 1L
    s_strand <- as.character(subject_gff@strand)

    s_by_seq <- split(seq_along(s_seq), s_seq)
    q_hits_list <- vector("list", length(q_seq))

    for (i in seq_along(q_seq)) {
      cand <- s_by_seq[[q_seq[i]]]
      if (is.null(cand) || length(cand) == 0) next
      ov <- cand[q_start[i] <= s_end[cand] & q_end[i] >= s_start[cand]]
      if (length(ov) == 0) next
      if (!ignore_strand && q_strand[i] != "*") {
        ov <- ov[s_strand[ov] == "*" | s_strand[ov] == q_strand[i]]
      }
      if (length(ov) > 0) {
        q_hits_list[[i]] <- data.frame(Captus_Row = rep.int(i, length(ov)), Ref_Row = ov, stringsAsFactors = FALSE)
      }
    }
    res_df <- dplyr::bind_rows(q_hits_list)
    if (nrow(res_df) == 0) {
      data.frame(Captus_Row = integer(0), Ref_Row = integer(0), stringsAsFactors = FALSE)
    } else {
      res_df
    }
  })

  if (nrow(df_cruces) == 0) {
    return(data.frame(
      Captus_Row = integer(0),
      Ref_Row = integer(0),
      ID_Captus = character(0),
      OG_Captus = character(0),
      Locus_Ref_Limpio = character(0),
      stringsAsFactors = FALSE
    ))
  }

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
#' Diseñado para genomas altamente fragmentados a escala de contigs (ej. S. stellatiglandulosum)
#' donde el atributo Target codifica `<transcript_id>-<OG_ID> <coordenadas>` (ej. `g65393.t1-OG0000254 1-148`).
#' Reproduce exactamente las lineas 1343-1353 de Workflow_Central_Sin_Rescate.R.
#'
#' @param query_gff Objeto GRanges de CAPTUS (protein_match:NUC).
#' @return Data frame con `ID_Captus`, `OG_Captus`, `Locus_Ref_Limpio`, `Contig_Name` y `Num_Gen`.
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
    Contig_Name = as.character(query_gff@seqnames),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      Num_Gen = suppressWarnings(as.numeric(stringr::str_extract(Locus_Ref_Limpio, "(?<=^g)[0-9]+")))
    ) %>%
    dplyr::distinct()
}

#' Construye el diccionario de homologia teorica a partir del archivo Excel de referencia
#'
#' Reproduce con exactitud matematica el filtrado de vias metabolicas, hojas de base de datos
#' y ensamblajes de `Workflow_Central_Sin_Rescate.R` para:
#' - `Storvum` (escala cromosomica: `CNV_Benoit` + `Contigs_probe_design` + `Sequences_Torva`)
#' - `Scitrullifolium` (escala cromosomica: `CNV_Benoit` + `Contigs_probe_design` + `S._citr`)
#' - `S_stella` (escala de contigs: `CNV_Perez` + `Contigs_probe_design` + `Sequences_Torva` + Puente `OG_ID <-> HOG_ID`)
#'
#' @param dict_path Ruta al archivo de diccionario (.xlsx, .tsv, .csv).
#' @param species Especie/muestra objetivo ("Storvum", "Scitrullifolium", "S_stella" o NULL).
#' @return Data frame con columnas `Locus_Benoit` y `OG_Benoit` (y atributo `puente_og_hog` en modo contig/HOG).
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

    # Caso 1: S. stellatiglandulosum (Genoma fragmentado a escala de contigs: puente HOG -> OG desde CNV_Perez)
    if (!is.null(species) && grepl("stella", species, ignore.case = TRUE) &&
        all(c("CNV_Perez", "Contigs_probe_design", "Sequences_Torva") %in% sheets)) {
      columna_ste <- "S._stellatiglandulosum"
      cnv_ste_raw <- readxl::read_excel(dict_path, sheet = "CNV_Perez") %>%
        dplyr::select(HOG_ID, Probe_ID, dplyr::all_of(columna_ste)) %>%
        dplyr::filter(!is.na(.data[[columna_ste]]) & trimws(.data[[columna_ste]]) != "")
      probes_validos_ste <- cnv_ste_raw$Probe_ID

      contigs_probe_ste <- readxl::read_excel(dict_path, sheet = "Contigs_probe_design") %>%
        dplyr::select(Probe_ID, Contig_name, HOG_ID_Solanum, OG_ID_Pangenomes,
                      Correcte_Contig_Name, Short_name, Pathway, Refseq_RNA, DAVID_Annotation) %>%
        dplyr::filter(Probe_ID %in% probes_validos_ste)

      cnv_ste <- cnv_ste_raw %>%
        dplyr::left_join(contigs_probe_ste %>% dplyr::select(Probe_ID, OG_ID = OG_ID_Pangenomes), by = "Probe_ID") %>%
        dplyr::select(OG_ID, Probe_ID, dplyr::all_of(columna_ste))

      sequences_ste <- readxl::read_excel(dict_path, sheet = "Sequences_Torva") %>%
        dplyr::filter(stringr::str_detect(Assembly, "stellati")) %>%
        dplyr::filter(Probe_ID %in% probes_validos_ste)

      datos_unificados_ste <- cnv_ste %>%
        dplyr::select(-OG_ID) %>%
        dplyr::left_join(contigs_probe_ste, by = "Probe_ID") %>%
        dplyr::left_join(sequences_ste %>% dplyr::select(-Orthogroup_ID), by = "Probe_ID")

      dict_ste <- datos_unificados_ste %>%
        dplyr::filter(trimws(Pathway) %in% pathways_filtro) %>%
        dplyr::select(-dplyr::all_of(columna_ste)) %>%
        dplyr::left_join(sequences_ste %>% dplyr::select(Seq, Orthogroup_ID), by = "Seq") %>%
        dplyr::filter(!is.na(Contig_ID)) %>%
        dplyr::select(Locus_Benoit = Contig_ID, OG_Benoit = Orthogroup_ID) %>%
        dplyr::mutate(Locus_Benoit = stringr::str_remove(Locus_Benoit, "\\.[0-9]+$")) %>%
        dplyr::distinct()

      # Tabla puente OG <-> HOG exacta de Workflow_Central_Sin_Rescate.R (lineas 1361-1369)
      puente_og_hog <- readxl::read_excel(dict_path, sheet = "Contigs_probe_design") %>%
        dplyr::select(OG_ID = OG_ID_Pangenomes, HOG_ID = HOG_ID_Solanum) %>%
        dplyr::filter(!is.na(OG_ID), !is.na(HOG_ID)) %>%
        dplyr::distinct() %>%
        dplyr::add_count(OG_ID, name = "n_hog") %>%
        dplyr::group_by(OG_ID) %>%
        dplyr::slice_max(order_by = n_hog, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::select(OG_ID, HOG_ID)

      attr(dict_ste, "puente_og_hog") <- puente_og_hog
      attr(dict_ste, "genome_scale")  <- "contig"
      attr(dict_ste, "db_sheets")     <- c("CNV_Perez", "Contigs_probe_design", "Sequences_Torva (Assembly: stellati)")
      return(dict_ste)
    }

    # Caso 2: S. citrullifolium (Genoma a escala cromosomica: hoja S._citr + CNV_Benoit)
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
        dplyr::filter(trimws(Pathway) %in% pathways_filtro) %>%
        dplyr::select(-`S._citrullifolium`) %>%
        dplyr::left_join(sequences_cit %>% dplyr::select(Seq, Orthogroup_ID), by = "Seq") %>%
        dplyr::filter(!is.na(Contig_ID)) %>%
        dplyr::select(Locus_Benoit = Contig_ID, OG_Benoit = Orthogroup_ID) %>%
        dplyr::mutate(Locus_Benoit = stringr::str_remove(Locus_Benoit, "\\.[0-9]+$")) %>%
        dplyr::distinct()

      attr(dict_cit, "genome_scale") <- "chromosome"
      attr(dict_cit, "db_sheets")    <- c("CNV_Benoit", "Contigs_probe_design", "S._citr")
      return(dict_cit)
    }

    # Caso 3: S. torvum (Experimento central a escala cromosomica)
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

      attr(dict_benoit, "genome_scale") <- "chromosome"
      attr(dict_benoit, "db_sheets")    <- c("CNV_Benoit", "Contigs_probe_design", "Sequences_Torva (Assembly: S. torvum Benoit et al)")
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
#' Soporta tanto el cruce directo OG == OG para genomas a escala cromosomica (ej. S. torvum y S. citrullifolium)
#' como el cruce con tabla puente OG <-> HOG para genomas a escala de contigs (ej. S. stellatiglandulosum),
#' reproduciendo exactamente Workflow_Central_Sin_Rescate.R.
#'
#' @param df_cruces Data frame obtenido con `find_spatial_overlaps()` o `extract_target_loci_from_gff()`.
#' @param dict_ref Data frame del diccionario con identificadores de locus y ortogrupo.
#' @param locus_col Columna de locus en `df_cruces`.
#' @param dict_locus_col Columna de locus en `dict_ref`.
#' @param dict_og_col Columna de ortogrupo/HOG en `dict_ref`.
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

  puente_og_hog <- attr(dict_ref, "puente_og_hog")

  # Modo Genoma a Escala de Contigs con Puente OG <-> HOG (Workflow_Central_Sin_Rescate.R lineas 1376-1415)
  if (!is.null(puente_og_hog)) {
    df_clasificacion <- df_cruces %>%
      dplyr::left_join(dict_ref, by = stats::setNames(dict_locus_col, locus_col), relationship = "many-to-many") %>%
      dplyr::left_join(puente_og_hog, by = c("OG_Captus" = "OG_ID"), relationship = "many-to-one") %>%
      dplyr::mutate(
        ID_Captus = dplyr::if_else(
          stringr::str_detect(ID_Captus, "__[0-9]+$"),
          ID_Captus,
          paste0(ID_Captus, "__00")
        ),
        Validacion_Estricta = dplyr::case_when(
          !is.na(.data[[dict_og_col]]) & !is.na(HOG_ID) & HOG_ID == .data[[dict_og_col]] ~ "TRUE_POSITIVE",
          TRUE ~ "FALSE_POSITIVE"
        ),
        Validacion_Espacial = Validacion_Estricta,
        Doble_Control_Biologico = dplyr::case_when(
          Validacion_Estricta == "TRUE_POSITIVE" ~ "1_Validado_Exacto",
          !is.na(.data[[dict_og_col]]) & (is.na(HOG_ID) | HOG_ID != .data[[dict_og_col]]) ~ "3_Hit_Cruzado_OG_Vecino",
          TRUE ~ "4_Expansion_No_Anotada_o_Ruido"
        ),
        prioridad = dplyr::case_when(
          Doble_Control_Biologico == "1_Validado_Exacto" ~ 0L,
          TRUE ~ 1L
        )
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

  # Modo Genoma a Escala Cromosomica con OG_ID directo (Workflow_Central_Sin_Rescate.R lineas 100-117 y 1015-1034)
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

#' Prepara de extremo a extremo el dataset etiquetado (TP/FP) de una especie segun la escala de su genoma
#'
#' Ejecuta el flujo completo de preparacion para una especie a escala cromosomica (`genome_scale = "chromosome"`,
#' modo `overlap` con GFF3 de genes) o a escala de contigs (`genome_scale = "contig"`, modo `target_attr`
#' extrayendo el locus desde el atributo `Target` del GFF de CAPTUS y aplicando el puente `OG_ID <-> HOG_ID`).
#'
#' @param sample_name Nombre exacto de la muestra en `captus_stats` (ej. `"Scitrullifolium"` o `"S_stella"`).
#' @param genome_scale Escala del genoma: `"chromosome"` (escala cromosomica) o `"contig"` (escala de contigs).
#' @param captus_gff Ruta al archivo GFF de contigs de CAPTUS (`protein_match:NUC`).
#' @param captus_stats Ruta al archivo `captus-extract_stats.tsv`.
#' @param dict_ref Ruta al archivo Excel/TSV de diccionario (`Genes_Capture_May_2026.xlsx`).
#' @param gff_ref Ruta al archivo GFF3 de genes del genoma de referencia (obligatorio si `genome_scale == "chromosome"`).
#' @param save_rds Ruta opcional para guardar el data frame preparado (.rds).
#' @return Lista con `data` (data frame preparado), `metadata` (informacion algoritmica y de base de datos).
#' @export
prepare_transfer_species_data <- function(sample_name,
                                          genome_scale = c("chromosome", "contig"),
                                          captus_gff,
                                          captus_stats,
                                          dict_ref,
                                          gff_ref = NULL,
                                          save_rds = NULL) {
  genome_scale <- match.arg(genome_scale)
  stopifnot(file.exists(captus_gff), file.exists(captus_stats), file.exists(dict_ref))
  if (genome_scale == "chromosome" && (is.null(gff_ref) || !file.exists(gff_ref))) {
    stop(sprintf("[ERROR] Para genoma a escala cromosomica (%s) se requiere un archivo gff_ref valido.", sample_name), call. = FALSE)
  }

  suppressPackageStartupMessages(library(rtracklayer))
  gff_captus_raw <- rtracklayer::import(captus_gff)
  gff_captus <- filter_gff_features(gff_captus_raw, feature_type = "protein_match:NUC")

  n_genes_ref <- NA_integer_
  if (genome_scale == "chromosome") {
    gff_ref_raw <- rtracklayer::import(gff_ref)
    genes_ref <- filter_gff_features(gff_ref_raw, feature_type = "gene")
    n_genes_ref <- length(genes_ref)
    df_cruces <- find_spatial_overlaps(query_gff = gff_captus, subject_gff = genes_ref, ignore_strand = FALSE)
  } else {
    df_cruces <- extract_target_loci_from_gff(gff_captus)
  }

  dict_df <- load_reference_dictionary(dict_ref, species = sample_name)
  puente_og_hog <- attr(dict_df, "puente_og_hog")
  db_sheets     <- attr(dict_df, "db_sheets")

  n_ogs_captus <- length(unique(df_cruces$OG_Captus))
  n_ogs_bridge <- if (!is.null(puente_og_hog)) sum(unique(df_cruces$OG_Captus) %in% puente_og_hog$OG_ID) else n_ogs_captus

  col_locus <- intersect(c("Locus_Benoit", "Locus_Ref_Limpio", "Locus", "gene_id"), colnames(dict_df))[1]
  col_og    <- intersect(c("OG_Benoit", "OG_ID", "Orthogroup_ID", "OG"), colnames(dict_df))[1]

  df_clasificacion <- classify_spatial_hits(
    df_cruces = df_cruces,
    dict_ref = dict_df,
    locus_col = "Locus_Ref_Limpio",
    dict_locus_col = col_locus,
    dict_og_col = col_og
  )

  stats_crudo <- utils::read.delim(captus_stats, header = TRUE, stringsAsFactors = FALSE) %>%
    dplyr::filter(sample_name == !!sample_name)
  stats_limpio <- clean_captus_ids(stats_crudo, ref_col = "ref_name", hit_col = "hit")

  stats_integrado <- stats_limpio %>%
    dplyr::left_join(df_clasificacion, by = c("ID_GFF_Compatible" = "ID_Captus")) %>%
    dplyr::filter(marker_type == "NUC") %>%
    dplyr::mutate(
      Validacion_Espacial = tidyr::replace_na(Validacion_Espacial, "FALSE_POSITIVE"),
      Doble_Control_Biologico = tidyr::replace_na(Doble_Control_Biologico, "4_Ruido_Intergenico_o_Intron"),
      ref_len_matched = as.numeric(ref_len_matched),
      pct_recovered = as.numeric(pct_recovered),
      pct_identity = as.numeric(pct_identity),
      score = as.numeric(score),
      wscore = as.numeric(wscore)
    )

  stats_sin_redundancias <- resolve_redundancies(
    df = stats_integrado,
    locus_col = "Locus_Ref_Limpio",
    class_col = "Validacion_Espacial"
  )

  stats_preparado <- calculate_locus_deltas(
    df = stats_sin_redundancias,
    locus_col = "locus",
    score_col = "score",
    wscore_col = "wscore",
    identity_col = "pct_identity"
  )

  attr(stats_preparado, "genome_scale") <- genome_scale
  attr(stats_preparado, "target_mode")  <- if (genome_scale == "chromosome") "overlap" else "target_attr"
  attr(stats_preparado, "db_sheets")    <- db_sheets
  attr(stats_preparado, "n_ogs_captus") <- n_ogs_captus
  attr(stats_preparado, "n_ogs_bridge") <- n_ogs_bridge
  attr(stats_preparado, "n_genes_ref")  <- n_genes_ref
  attr(stats_preparado, "n_hits_gff")   <- nrow(df_cruces)

  if (!is.null(save_rds)) {
    dir.create(dirname(save_rds), recursive = TRUE, showWarnings = FALSE)
    saveRDS(stats_preparado, file = save_rds)
    tsv_out <- sub("\\.rds$", ".tsv", save_rds)
    if (tsv_out != save_rds) {
      utils::write.table(stats_preparado, file = tsv_out, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  }

  return(stats_preparado)
}

#' Extrae la Matriz de Variacion en el Numero de Copias (CNV) con regla de un solo umbral
#'
#' Ordena los homologos validados por cobertura (`pct_recovered`) de mayor a menor
#' dentro de cada muestra y locus, clasificando la primera secuencia como `Copia_Principal`,
#' las adicionales con cobertura >= `cov_threshold` como `Homologo_Confirmado` y las
#' menores a ese umbral como `No_Contado`.
#'
#' @param df Data frame con los veredictos del consenso de 7 modelos.
#' @param sample_col Nombre de la columna de muestra.
#' @param locus_col Nombre de la columna de locus/ortogrupo.
#' @param cov_threshold Cobertura individual minima (%) para confirmar un homologo adicional (default: 50.0).
#' @param verdict_col Columna con el veredicto de consenso.
#' @param ... Argumentos adicionales ignorados por compatibilidad.
#' @return Lista con `cnv_matrix`, `cnv_annotated` y `functional_genes`.
#' @export
extract_cnv_matrix <- function(df,
                               sample_col = "Sample",
                               locus_col = "Locus_Ortogrupo",
                               cov_threshold = 50.0,
                               verdict_col = "Veredicto_Consenso",
                               ...) {
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
      Clasificacion_Copia = dplyr::case_when(
        dplyr::row_number() == 1 ~ "Copia_Principal",
        dplyr::row_number() > 1 & as.numeric(pct_recovered) >= cov_threshold ~ "Homologo_Confirmado",
        TRUE ~ "No_Contado"
      )
    ) %>%
    dplyr::ungroup()

  df_funcionales <- df_cnv_bruto %>%
    dplyr::filter(Clasificacion_Copia %in% c("Copia_Principal", "Homologo_Confirmado")) %>%
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
