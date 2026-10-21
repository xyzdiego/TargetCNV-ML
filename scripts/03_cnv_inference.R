#!/usr/bin/env Rscript

# ==============================================================================
# TargetCNV-ML Pipeline
# Script 03: Production CNV Inference & Multi-Model Consensus System
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(caret)
  library(pROC)
})

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
source(file.path(project_root, "R", "model_evaluation.R"))
source(file.path(project_root, "R", "plotting.R"))

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = NULL,
    help = "Ruta al archivo TSV/CSV con datos nuevos de CAPTUS [Obligatorio]",
    metavar = "FILE"
  ),
  make_option(
    c("-m", "--models_dir"),
    type = "character",
    default = file.path(project_root, "models"),
    help = "Directorio con los modelos y umbrales serializados (.rds) [Por defecto: models/]",
    metavar = "DIR"
  ),
  make_option(
    c("-o", "--out_dir"),
    type = "character",
    default = file.path(project_root, "results"),
    help = "Directorio de salida para matrices CNV, tablas y figuras [Por defecto: results/]",
    metavar = "DIR"
  ),
  make_option(
    c("-c", "--cov_threshold"),
    type = "double",
    default = 50.0,
    help = "Umbral de cobertura individual (%%) para confirmar una copia homóloga [Por defecto: 50.0]",
    metavar = "NUM"
  ),
  make_option(
    c("-f", "--frag_threshold"),
    type = "double",
    default = 120.0,
    help = "Umbral de cobertura acumulada (%%) para diagnosticar fragmentacion [Por defecto: 120.0]",
    metavar = "NUM"
  ),
  make_option(
    c("-v", "--min_votes"),
    type = "integer",
    default = 4,
    help = "Minimo de votos afirmativos para consenso (mayoría simple de 7 modelos) [Por defecto: 4]",
    metavar = "INT"
  ),
  make_option(
    c("-t", "--top_loci"),
    type = "integer",
    default = 30,
    help = "Numero de familias/grupos homólogos mas abundantes para graficar en el heatmap [Por defecto: 30]",
    metavar = "INT"
  )
)

opt_parser <- OptionParser(
  usage = "Uso: %prog [opciones]",
  option_list = option_list,
  description = "Aplica el Sistema de Consenso de 7 modelos (ML + Youden + Heuristico) y genera la matriz CNV."
)

opt <- parse_args(opt_parser)

if (is.null(opt$input) || !file.exists(opt$input)) {
  print_help(opt_parser)
  stop(sprintf("[ERROR] Archivo de entrada no especificado o inexistente: %s", opt$input), call. = FALSE)
}

if (!dir.exists(opt$models_dir)) {
  stop(sprintf("[ERROR] El directorio de modelos no existe: %s", opt$models_dir), call. = FALSE)
}

if (!dir.exists(opt$out_dir)) {
  dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
}

cat("=================================================================\n")
cat(" TargetCNV-ML: Modulo 03 - Inferencia de Produccion y Matriz CNV \n")
cat("=================================================================\n")
cat(sprintf("Archivo de Entrada (CAPTUS) : %s\n", opt$input))
cat(sprintf("Directorio de Modelos (.rds): %s\n", opt$models_dir))
cat(sprintf("Directorio de Salida        : %s\n", opt$out_dir))
cat(sprintf("Umbral Cobertura Homólogo   : %.1f%%\n", opt$cov_threshold))
cat(sprintf("Umbral Max. Fragmentacion   : %.1f%%\n", opt$frag_threshold))
cat(sprintf("Votos Minimos Consenso      : %d / 7\n", opt$min_votes))
cat("-----------------------------------------------------------------\n")

# 1. Carga de Modelos y Umbrales
cat("[1/6] Cargando artefactos de Machine Learning y umbrales optimizados...\n")

# Carga de umbrales Youden y modelo heuristico
youden_file <- file.path(opt$models_dir, "youden_thresholds.rds")
heuristic_file <- file.path(opt$models_dir, "heuristic_model.rds")

if (!file.exists(youden_file)) {
  stop(sprintf("[ERROR] No se encontro 'youden_thresholds.rds' en %s", opt$models_dir), call. = FALSE)
}
if (!file.exists(heuristic_file)) {
  stop(sprintf("[ERROR] No se encontro 'heuristic_model.rds' en %s", opt$models_dir), call. = FALSE)
}

youden_thresholds <- readRDS(youden_file)
heuristic_model   <- readRDS(heuristic_file)

# Busqueda y carga de modelos supervisados
posibles_modelos <- list(
  ElasticNet     = file.path(opt$models_dir, "model_elasticnet.rds"),
  NaiveBayes     = file.path(opt$models_dir, "model_naivebayes.rds"),
  SVM_Polinomial = file.path(opt$models_dir, "model_svm_polinomial.rds"),
  Red_Neuronal   = file.path(opt$models_dir, "model_red_neuronal.rds"),
  LDA            = file.path(opt$models_dir, "model_lda.rds"),
  DNN_Profunda   = file.path(opt$models_dir, "model_dnn_profunda.rds")
)

ml_models <- list()
for (m_name in names(posibles_modelos)) {
  m_path <- posibles_modelos[[m_name]]
  if (file.exists(m_path)) {
    ml_models[[m_name]] <- readRDS(m_path)
    cat(sprintf("  -> Modelo cargado: %s\n", m_name))
  }
}

# Si existe un contenedor global 'all_trained_models.rds', usarlo como respaldo
if (length(ml_models) == 0 && file.exists(file.path(opt$models_dir, "all_trained_models.rds"))) {
  ml_models <- readRDS(file.path(opt$models_dir, "all_trained_models.rds"))
  cat(sprintf("  -> Cargados %d modelos desde all_trained_models.rds\n", length(ml_models)))
}

if (length(ml_models) == 0) {
  stop("[ERROR] No se encontro ningun modelo ML serializado en el directorio especificado.", call. = FALSE)
}

# 2. Lectura y Estandarizacion de Datos Nuevos
cat("[2/6] Cargando y preprocesando datos de CAPTUS...\n")
datos_crudos <- utils::read.delim(opt$input, header = TRUE, stringsAsFactors = FALSE)

columnas_requeridas <- c("sample_name", "marker_type", "ref_name", "hit",
                         "pct_recovered", "pct_identity", "score", "wscore")
columnas_faltantes <- setdiff(columnas_requeridas, colnames(datos_crudos))
if (length(columnas_faltantes) > 0) {
  stop(sprintf(
    "[ERROR] El archivo de entrada no tiene el formato de 'captus-extract_stats.tsv'. Faltan columnas: %s",
    paste(columnas_faltantes, collapse = ", ")
  ), call. = FALSE)
}

metricas_numericas <- c("pct_recovered", "pct_identity", "score", "wscore")

df_inferencia <- datos_crudos %>%
  dplyr::filter(marker_type == "NUC") %>%
  clean_captus_ids(ref_col = "ref_name", hit_col = "hit") %>%
  dplyr::mutate(
    Sample = dplyr::coalesce(.data[["sample_name"]], "Muestra_Desconocida"),
    ID_Secuencia = ID_GFF_Compatible,
    Locus_Ortogrupo = dplyr::coalesce(.data[["locus"]], OG_Limpio),
    dplyr::across(dplyr::all_of(metricas_numericas), ~ as.numeric(as.character(.x)))
  ) %>%
  dplyr::filter(stats::complete.cases(dplyr::select(., dplyr::all_of(metricas_numericas)))) %>%
  as.data.frame()

cat(sprintf("  -> Secuencias analizadas: %d procedentes de %d muestra(s)\n",
            nrow(df_inferencia), length(unique(df_inferencia$Sample))))

# 3. Aplicacion del Sistema de Consenso de 7 Modelos
cat("[3/6] Evaluando el Sistema de Consenso de 7 Modelos...\n")

df_eval <- df_inferencia
modelos_votantes <- character()

# A. Modelos ML supervisados (Top 4 disponibles)
top_ml_names <- head(names(ml_models), 4)
if (length(top_ml_names) < 4) {
  warning(sprintf(
    "Solo hay %d modelo(s) ML disponible(s) (se esperan 4): el consenso tendra %d votantes en lugar de 7. Revise --min_votes.",
    length(top_ml_names), length(top_ml_names) + 3
  ), call. = FALSE)
}
for (mod_name in top_ml_names) {
  col_name <- paste0("Pred_ML_", mod_name)
  pred_raw <- stats::predict(ml_models[[mod_name]], newdata = df_inferencia[, metricas_numericas])
  df_eval[[col_name]] <- factor(
    ifelse(as.character(pred_raw) %in% c("TP", "TRUE_POSITIVE"), "TP", "FP"),
    levels = c("TP", "FP")
  )
  modelos_votantes <- c(modelos_votantes, col_name)
  cat(sprintf("  [Voto 1..4] Evaluado %s\n", col_name))
}

# B. Cortes Univariados Youden (Top 2: pct_recovered y pct_identity)
top2_youden_vars <- c("pct_recovered", "pct_identity")
for (y_var in top2_youden_vars) {
  col_name <- paste0("Pred_Youden_", y_var)
  umbral_y <- youden_thresholds[[y_var]]
  df_eval[[col_name]] <- factor(
    ifelse(df_eval[[y_var]] >= umbral_y, "TP", "FP"),
    levels = c("TP", "FP")
  )
  modelos_votantes <- c(modelos_votantes, col_name)
  cat(sprintf("  [Voto 5..6] Evaluado Youden: %s (umbral >= %.2f)\n", y_var, umbral_y))
}

# C. Ensamble Heuristico Ponderado
col_heur <- "Pred_Ens_Heuristico"
df_eval[[col_heur]] <- predict_weighted_heuristic(
  df_predictors = df_eval,
  thresholds = youden_thresholds,
  weights = heuristic_model$weights,
  global_threshold = heuristic_model$global_threshold
)
modelos_votantes <- c(modelos_votantes, col_heur)
cat(sprintf("  [Voto 7] Evaluado Ensamble Heuristico Ponderado (umbral >= %.2f)\n",
            heuristic_model$global_threshold))

# Consolidacion del veredicto multimodelo
df_consenso <- consolidar_consenso(
  df_predictions = df_eval,
  model_cols = modelos_votantes,
  min_votes = opt$min_votes
)

cat("  -> Resumen del Veredicto de Consenso:\n")
print(table(df_consenso$Veredicto_Consenso))

# 4. Control de Fragmentacion Espacial y Extraccion de Matriz CNV
cat("[4/6] Ejecutando control de fragmentacion tecnica y extraccion de CNV...\n")

res_cnv <- extract_cnv_matrix(
  df = df_consenso,
  sample_col = "Sample",
  locus_col = "Locus_Ortogrupo",
  cov_threshold = opt$cov_threshold,
  frag_threshold = opt$frag_threshold,
  verdict_col = "Veredicto_Consenso"
)

df_matriz_cnv      <- res_cnv$cnv_matrix
df_cnv_anotado     <- res_cnv$cnv_annotated
df_funcionales     <- res_cnv$functional_genes

cat(sprintf("  -> Loci/Familias genicas detectadas : %d\n", nrow(df_matriz_cnv)))
cat(sprintf("  -> Copias funcionales confirmadas   : %d\n", nrow(df_funcionales)))
cat("  -> Desglose de clasificacion de copias:\n")
print(table(df_cnv_anotado$Clasificacion_Copia))

# 5. Exportacion de Archivos y Tablas
cat("[5/6] Exportando matrices y catalogos a disco...\n")

file_matriz_tsv   <- file.path(opt$out_dir, "cnv_matrix.tsv")
file_anotado_tsv  <- file.path(opt$out_dir, "cnv_annotated_sequences.tsv")
file_func_tsv     <- file.path(opt$out_dir, "functional_genes_catalog.tsv")
file_consenso_tsv <- file.path(opt$out_dir, "consensus_verdicts.tsv")

utils::write.table(df_matriz_cnv, file = file_matriz_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(df_cnv_anotado, file = file_anotado_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(df_funcionales, file = file_func_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(df_consenso, file = file_consenso_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# 6. Generacion y Guardado de Graficos Diagnosticos
cat("[6/6] Renderizando graficos de publicacion (ggplot2)...\n")

# A. Grafico de Veredicto de Consenso
p_consenso <- plot_consensus_distribution(
  df_consensus = df_consenso,
  verdict_col = "Veredicto_Consenso"
)
tryCatch({
  ggplot2::ggsave(file.path(opt$out_dir, "plot_consensus_distribution.pdf"), plot = p_consenso, width = 7, height = 5)
  ggplot2::ggsave(file.path(opt$out_dir, "plot_consensus_distribution.png"), plot = p_consenso, width = 7, height = 5, device = grDevices::png)
}, error = function(e) {
  warning(sprintf("Aviso al guardar graficos de consenso: %s", e$message))
})

# B. Heatmap de Variacion en el Numero de Copias (CNV)
if (nrow(df_funcionales) > 0) {
  p_heatmap <- plot_cnv_heatmap(
    df_cnv_annotated = df_cnv_anotado,
    sample_col = "Sample",
    locus_col = "Locus_Ortogrupo",
    top_n_loci = opt$top_loci
  )
  tryCatch({
    ggplot2::ggsave(file.path(opt$out_dir, "plot_cnv_heatmap.pdf"), plot = p_heatmap, width = 9, height = 7)
    ggplot2::ggsave(file.path(opt$out_dir, "plot_cnv_heatmap.png"), plot = p_heatmap, width = 9, height = 7, device = grDevices::png)
  }, error = function(e) {
    warning(sprintf("Aviso al guardar graficos de heatmap: %s", e$message))
  })
}

cat("-----------------------------------------------------------------\n")
cat(" VISTA PREVIA DE LA MATRIZ DE CNV:\n")
cat("-----------------------------------------------------------------\n")
print(utils::head(df_matriz_cnv))
cat("=================================================================\n")
cat(sprintf("[OK] Pipeline de inferencia finalizado con exito. Salidas en: %s\n", opt$out_dir))
