#!/usr/bin/env Rscript

# ==============================================================================
# TargetCNV-ML Pipeline
# Script 01: Data Preparation and Spatial Filtering
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

# Ubicacion dinamica del proyecto y carga de modulos
get_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    if (basename(script_dir) == "scripts") dirname(script_dir) else script_dir
  } else {
    getwd()
  }
}

project_root <- get_project_root()
source(file.path(project_root, "R", "data_processing.R"))

# Definicion de argumentos CLI con optparse
option_list <- list(
  make_option(
    c("-g", "--gff_ref"),
    type = "character",
    default = NULL,
    help = "Ruta al archivo GFF/GFF3 del genoma de referencia [Obligatorio]",
    metavar = "FILE"
  ),
  make_option(
    c("-c", "--captus_gff"),
    type = "character",
    default = NULL,
    help = "Ruta al archivo GFF de contigs de CAPTUS (protein_match:NUC) [Obligatorio]",
    metavar = "FILE"
  ),
  make_option(
    c("-s", "--captus_stats"),
    type = "character",
    default = NULL,
    help = "Ruta al archivo TSV de estadisticas de CAPTUS (captus-extract_stats.tsv) [Obligatorio]",
    metavar = "FILE"
  ),
  make_option(
    c("-d", "--dict_ref"),
    type = "character",
    default = NULL,
    help = "Ruta al diccionario de referencia (Excel .xlsx o TSV/CSV) [Opcional]",
    metavar = "FILE"
  ),
  make_option(
    c("-p", "--sample"),
    type = "character",
    default = NULL,
    help = "Nombre exacto de la muestra para filtrar en captus_stats (ej. Storvum) [Opcional]",
    metavar = "STR"
  ),
  make_option(
    c("-o", "--out_dir"),
    type = "character",
    default = file.path(project_root, "data"),
    help = "Directorio donde se guardaran los datos preparados [Por defecto: data/]",
    metavar = "DIR"
  )
)

opt_parser <- OptionParser(
  usage = "Uso: %prog [opciones]",
  option_list = option_list,
  description = "Prepara datos de CAPTUS, ejecuta cruces espaciales con GFF, resuelve redundancias y genera dataset de ML."
)

opt <- parse_args(opt_parser)

# Validacion de argumentos requeridos
if (is.null(opt$gff_ref) || is.null(opt$captus_gff) || is.null(opt$captus_stats)) {
  print_help(opt_parser)
  stop("[ERROR] Debe especificar --gff_ref, --captus_gff y --captus_stats.", call. = FALSE)
}

if (!file.exists(opt$gff_ref)) {
  stop(sprintf("[ERROR] No se encontro el archivo GFF de referencia: %s", opt$gff_ref), call. = FALSE)
}
if (!file.exists(opt$captus_gff)) {
  stop(sprintf("[ERROR] No se encontro el archivo GFF de CAPTUS: %s", opt$captus_gff), call. = FALSE)
}
if (!file.exists(opt$captus_stats)) {
  stop(sprintf("[ERROR] No se encontro el archivo de estadisticas: %s", opt$captus_stats), call. = FALSE)
}
# Un diccionario indicado pero inexistente NO debe degradar en silencio a "todo TRUE_POSITIVE"
if (!is.null(opt$dict_ref) && !file.exists(opt$dict_ref)) {
  stop(sprintf("[ERROR] No se encontro el diccionario de referencia (--dict_ref): %s", opt$dict_ref), call. = FALSE)
}

if (!dir.exists(opt$out_dir)) {
  dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
}

cat("=================================================================\n")
cat(" TargetCNV-ML: Modulo 01 - Preparacion de Datos y Cruce Espacial \n")
cat("=================================================================\n")
cat(sprintf("Genoma de Referencia GFF : %s\n", opt$gff_ref))
cat(sprintf("CAPTUS Contigs GFF       : %s\n", opt$captus_gff))
cat(sprintf("Estadisticas CAPTUS      : %s\n", opt$captus_stats))
cat(sprintf("Directorio de Salida     : %s\n", opt$out_dir))
cat("-----------------------------------------------------------------\n")

# 1. Carga de anotaciones GFF
cat("[1/6] Cargando y filtrando anotaciones GFF...\n")
suppressPackageStartupMessages(library(rtracklayer))

gff_ref_raw <- rtracklayer::import(opt$gff_ref)
genes_ref <- filter_gff_features(gff_ref_raw, feature_type = "gene")
cat(sprintf("  -> Genes de referencia cargados: %d\n", length(genes_ref)))

gff_captus_raw <- rtracklayer::import(opt$captus_gff)
gff_captus <- filter_gff_features(gff_captus_raw, feature_type = "protein_match:NUC")
cat(sprintf("  -> Marcadores funcionales (NUC) CAPTUS: %d\n", length(gff_captus)))

# 2. Interseccion espacial
cat("[2/6] Ejecutando interseccion espacial estricta (GenomicRanges)...\n")
df_cruces <- find_spatial_overlaps(query_gff = gff_captus, subject_gff = genes_ref, ignore_strand = FALSE)
cat(sprintf("  -> Hits espaciales encontrados: %d\n", nrow(df_cruces)))

# 3. Clasificacion con diccionario de referencia (si esta disponible)
cat("[3/6] Evaluando correspondencia con diccionario de referencia...\n")
if (!is.null(opt$dict_ref) && file.exists(opt$dict_ref)) {
  if (grepl("\\.xlsx$", opt$dict_ref, ignore.case = TRUE)) {
    suppressPackageStartupMessages(library(readxl))
    sheets <- readxl::excel_sheets(opt$dict_ref)
    if ("CNV_Benoit" %in% sheets && "Sequences_Torva" %in% sheets) {
      cnv_b <- readxl::read_excel(opt$dict_ref, sheet = "CNV_Benoit") %>%
        dplyr::select(OG_ID, Probe_ID, dplyr::any_of(c("S._torvum", "S_torvum")))
      seq_t <- readxl::read_excel(opt$dict_ref, sheet = "Sequences_Torva")
      dict_ref <- cnv_b %>%
        dplyr::left_join(seq_t, by = "Probe_ID", relationship = "many-to-many") %>%
        dplyr::filter(!is.na(Contig_ID)) %>%
        dplyr::select(Locus_Benoit = Contig_ID, OG_Benoit = Orthogroup_ID) %>%
        dplyr::mutate(Locus_Benoit = stringr::str_remove(Locus_Benoit, "\\.[0-9]+$")) %>%
        dplyr::distinct()
    } else {
      dict_ref <- readxl::read_excel(opt$dict_ref)
    }
  } else {
    dict_ref <- utils::read.delim(opt$dict_ref, stringsAsFactors = FALSE)
  }
  
  col_locus <- intersect(c("Locus_Benoit", "Locus_Ref_Limpio", "Locus", "gene_id"), colnames(dict_ref))[1]
  col_og    <- intersect(c("OG_Benoit", "OG_ID", "Orthogroup_ID", "OG"), colnames(dict_ref))[1]
  
  df_clasificacion <- classify_spatial_hits(
    df_cruces = df_cruces,
    dict_ref = dict_ref,
    locus_col = "Locus_Ref_Limpio",
    dict_locus_col = col_locus,
    dict_og_col = col_og
  )
} else {
  cat("  [AVISO] No se proporciono diccionario (--dict_ref): todos los hits con solape espacial se marcaran\n")
  cat("          como TRUE_POSITIVE. Con esta salida NO es posible entrenar (02_train_models.R requiere TP y FP).\n")
  df_clasificacion <- df_cruces %>%
    dplyr::mutate(
      Validacion_Espacial = "TRUE_POSITIVE",
      Doble_Control_Biologico = "1_Hit_Espacial_Valido"
    ) %>%
    dplyr::distinct(ID_Captus, .keep_all = TRUE)
}

# 4. Integracion con estadisticas de CAPTUS
cat("[4/6] Integrando metricas continuas de alineamiento CAPTUS...\n")
stats_crudo <- utils::read.delim(opt$captus_stats, header = TRUE, stringsAsFactors = FALSE)

if (!is.null(opt$sample)) {
  muestras_disponibles <- unique(stats_crudo$sample_name)
  if (!opt$sample %in% muestras_disponibles) {
    stop(sprintf(
      "[ERROR] La muestra '%s' no existe en %s.\n        Muestras disponibles (%d): %s%s",
      opt$sample, opt$captus_stats, length(muestras_disponibles),
      paste(utils::head(muestras_disponibles, 10), collapse = ", "),
      if (length(muestras_disponibles) > 10) ", ..." else ""
    ), call. = FALSE)
  }
  stats_crudo <- stats_crudo %>% dplyr::filter(sample_name == opt$sample)
  cat(sprintf("  -> Filtrado para la muestra: %s (%d filas)\n", opt$sample, nrow(stats_crudo)))
}

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

# 5. Resolucion de redundancias y calculo de deltas
cat("[5/6] Resolviendo redundancias y computando deltas locales intra-locus...\n")
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

# 6. Exportacion de artefactos
cat("[6/6] Guardando datos preparados...\n")
out_tsv <- file.path(opt$out_dir, "data_prepared.tsv")
out_rds <- file.path(opt$out_dir, "data_prepared.rds")

utils::write.table(stats_preparado, file = out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(stats_preparado, file = out_rds)

cat("-----------------------------------------------------------------\n")
cat(" RESUMEN FINAL DE PREPARACION\n")
cat("-----------------------------------------------------------------\n")
cat(sprintf("Total secuencias analizadas (NUC) : %d\n", nrow(stats_preparado)))
cat(sprintf("Verdaderos Positivos (TP)        : %d\n", sum(stats_preparado$Validacion_Espacial == "TRUE_POSITIVE")))
cat(sprintf("Falsos Positivos & Redundantes   : %d\n", sum(stats_preparado$Validacion_Espacial == "FALSE_POSITIVE")))
cat(sprintf("Archivo TSV generado             : %s\n", out_tsv))
cat(sprintf("Archivo RDS generado             : %s\n", out_rds))
cat("=================================================================\n")
cat("[OK] Proceso completado con exito.\n")
