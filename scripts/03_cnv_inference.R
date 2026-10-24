#!/usr/bin/env Rscript

# ==============================================================================
# TargetCNV-ML Pipeline
# Script 03: Production CNV Inference, Multi-Model Consensus & Transferability
# ==============================================================================

t_start <- Sys.time()

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
    help = "Ruta al archivo TSV/CSV/RDS con datos nuevos de CAPTUS [Obligatorio]",
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
    help = "Umbral de cobertura individual (%) para confirmar una copia homologa [Por defecto: 50.0]",
    metavar = "NUM"
  ),
  make_option(
    c("-f", "--frag_threshold"),
    type = "double",
    default = NULL,
    help = "[Obsoleto/Ignorado] Mantenido solo por retrocompatibilidad de linea de comandos",
    metavar = "NUM"
  ),
  make_option(
    c("-v", "--min_votes"),
    type = "integer",
    default = 4,
    help = "Minimo de votos afirmativos para consenso (mayoria simple de 7 modelos) [Por defecto: 4]",
    metavar = "INT"
  ),
  make_option(
    c("-t", "--top_loci"),
    type = "integer",
    default = 30,
    help = "Numero de familias/grupos homologos mas abundantes para graficar en el heatmap [Por defecto: 30]",
    metavar = "INT"
  ),
  make_option(
    c("--render_report"),
    action = "store_true",
    default = TRUE,
    help = "Compilar automaticamente el reporte HTML de inferencia (report_template.qmd) [Por defecto: TRUE]"
  ),
  make_option(
    c("--no_render_report"),
    action = "store_false",
    dest = "render_report",
    help = "Omitir la compilacion del reporte HTML de inferencia"
  ),
  make_option(
    c("-r", "--report_file"),
    type = "character",
    default = file.path(project_root, "reporte_inferencia.html"),
    help = "Ruta del reporte HTML de inferencia [Por defecto: reporte_inferencia.html]",
    metavar = "FILE"
  )
)

opt_parser <- OptionParser(
  usage = "Uso: %prog [opciones]",
  option_list = option_list,
  description = "Aplica el Sistema de Consenso de 7 modelos (ML + Youden + Heuristico), evalua transferibilidad si existe Ground Truth y genera la matriz CNV."
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

cat("========================================================\n")
cat(" TargetCNV-ML: Modulo 03 - Inferencia de Produccion y Matriz CNV \n")
cat("========================================================\n")
cat(sprintf("Archivo de Entrada (CAPTUS) : %s\n", opt$input))
cat(sprintf("Directorio de Modelos (.rds): %s\n", opt$models_dir))
cat(sprintf("Directorio de Salida        : %s\n", opt$out_dir))
cat(sprintf("Umbral Cobertura Homologo   : %.1f%%\n", opt$cov_threshold))
cat(sprintf("Votos Minimos Consenso      : %d / 7\n", opt$min_votes))
cat("--------------------------------------------------------\n")

# 1. Carga de Modelos y Umbrales
cat("[1/7] Cargando artefactos de Machine Learning y umbrales optimizados...\n")

youden_file <- file.path(opt$models_dir, "youden_thresholds.rds")
heuristic_file <- file.path(opt$models_dir, "heuristic_model.rds")

if (!file.exists(youden_file)) {
  stop(sprintf("[ERROR] No se encontro 'youden_thresholds.rds' en %s", opt$models_dir), call. = FALSE)
}
if (!file.exists(heuristic_file)) {
  stop(sprintf("[ERROR] No se encontro 'heuristic_model.rds' en %s", opt$models_dir), call. = FALSE)
}

youden_thresholds <- readRDS(youden_file)
heuristic_model <- readRDS(heuristic_file)

# Busqueda y carga de modelos supervisados (priorizando los 4 votantes de consenso)
posibles_modelos <- list(
  ElasticNet     = file.path(opt$models_dir, "model_elasticnet.rds"),
  NaiveBayes     = file.path(opt$models_dir, "model_naivebayes.rds"),
  SVM_Polinomial = file.path(opt$models_dir, "model_svm_polinomial.rds"),
  Red_Neuronal   = file.path(opt$models_dir, "model_red_neuronal.rds"),
  LDA            = file.path(opt$models_dir, "model_lda.rds"),
  MLP            = file.path(opt$models_dir, "model_mlp.rds")
)

ml_models <- list()
for (m_name in names(posibles_modelos)) {
  m_path <- posibles_modelos[[m_name]]
  if (file.exists(m_path)) {
    ml_models[[m_name]] <- readRDS(m_path)
    cat(sprintf("   -> Modelo cargado: %s\n", m_name))
  }
}

if (length(ml_models) == 0 && file.exists(file.path(opt$models_dir, "all_trained_models.rds"))) {
  ml_models <- readRDS(file.path(opt$models_dir, "all_trained_models.rds"))
  cat(sprintf("   -> Cargados %d modelos desde all_trained_models.rds\n", length(ml_models)))
}

if (length(ml_models) == 0) {
  stop("[ERROR] No se encontro ningun modelo ML serializado en el directorio especificado.", call. = FALSE)
}

# 2. Lectura y Estandarizacion de Datos Nuevos
cat("[2/7] Cargando y preprocesando datos de CAPTUS...\n")
if (grepl("\\.rds$", opt$input, ignore.case = TRUE)) {
  datos_crudos <- readRDS(opt$input)
} else {
  sep_char <- if (grepl("\\.csv$", opt$input, ignore.case = TRUE)) "," else "\t"
  datos_crudos <- utils::read.delim(opt$input, sep = sep_char, header = TRUE, stringsAsFactors = FALSE)
}

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
spec_file <- file.path(opt$models_dir, "consensus_spec.rds")
consensus_spec <- if (file.exists(spec_file)) readRDS(spec_file) else NULL

# Detectar si los modelos fueron entrenados con --include_hit (desde consensus_spec.rds o desde el propio modelo caret)
predictores_ml <- if (!is.null(consensus_spec) && !is.null(consensus_spec$predictors)) {
  consensus_spec$predictors
} else {
  mod1 <- ml_models[[1]]
  cols_m1 <- if (!is.null(mod1$trainingData)) setdiff(colnames(mod1$trainingData), ".outcome") else mod1$coefnames
  if (!is.null(cols_m1) && "hit" %in% cols_m1) c(metricas_numericas, "hit") else metricas_numericas
}
columnas_num_todas <- unique(c(metricas_numericas, predictores_ml))

df_inferencia <- datos_crudos %>%
  dplyr::filter(marker_type == "NUC") %>%
  clean_captus_ids(ref_col = "ref_name", hit_col = "hit") %>%
  dplyr::mutate(
    Sample = dplyr::coalesce(as.character(.data[["sample_name"]]), "Muestra_Desconocida"),
    ID_Secuencia = ID_GFF_Compatible,
    Locus_Ortogrupo = if ("locus" %in% colnames(.)) dplyr::coalesce(as.character(.data[["locus"]]), OG_Limpio) else OG_Limpio,
    dplyr::across(dplyr::all_of(columnas_num_todas), ~ as.numeric(as.character(.x)))
  ) %>%
  dplyr::filter(stats::complete.cases(dplyr::select(., dplyr::all_of(columnas_num_todas)))) %>%
  as.data.frame()

cat(sprintf("   -> Secuencias analizadas: %d procedentes de %d muestra(s) (Predictores: %s)\n",
            nrow(df_inferencia), length(unique(df_inferencia$Sample)), paste(predictores_ml, collapse = ", ")))

# 3. Aplicacion del Sistema de Consenso de 7 Modelos
cat("[3/7] Evaluando el Sistema de Consenso de 7 Modelos...\n")
df_eval <- df_inferencia
modelos_votantes <- character()

# A. Modelos ML supervisados (usando consensus_spec si existe, o los 4 canonicos)
voters_objetivo <- if (!is.null(consensus_spec) && !is.null(consensus_spec$ml_voters)) {
  consensus_spec$ml_voters
} else {
  c("ElasticNet", "NaiveBayes", "SVM_Polinomial", "Red_Neuronal")
}
preferidos_ml <- intersect(voters_objetivo, names(ml_models))
top_ml_names <- if (length(preferidos_ml) == 4) preferidos_ml else utils::head(unique(c(preferidos_ml, names(ml_models))), 4)

if (length(top_ml_names) < 4) {
  warning(sprintf(
    "Solo hay %d modelo(s) ML disponible(s) (se esperan 4): el consenso tendra %d votantes en lugar de 7. Revise --models_dir",
    length(top_ml_names), length(top_ml_names) + 3
  ), call. = FALSE)
}

df_x_pred <- df_inferencia[, predictores_ml, drop = FALSE]

for (mod_name in top_ml_names) {
  col_name <- paste0("Pred_ML_", mod_name)
  pred_raw <- stats::predict(ml_models[[mod_name]], newdata = df_x_pred)
  df_eval[[col_name]] <- factor(
    ifelse(as.character(pred_raw) %in% c("TP", "TRUE_POSITIVE"), "TP", "FP"),
    levels = c("TP", "FP")
  )
  modelos_votantes <- c(modelos_votantes, col_name)
  cat(sprintf("   [Voto 1..4] Evaluado %s\n", col_name))
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
  cat(sprintf("   [Voto 5..6] Evaluado Youden: %s (umbral >= %.2f)\n", y_var, umbral_y))
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
cat(sprintf("   [Voto 7] Evaluado Ensamble Heuristico Ponderado (umbral >= %.2f)\n",
            heuristic_model$global_threshold))

# Validar --min_votes frente al total real de modelos votantes
if (is.na(opt$min_votes) || opt$min_votes < 1 || opt$min_votes > length(modelos_votantes)) {
  stop(sprintf(
    "[ERROR] --min_votes (%s) debe estar entre 1 y el total de modelos votantes disponibles (%d).",
    as.character(opt$min_votes), length(modelos_votantes)
  ), call. = FALSE)
}

# Consolidacion del veredicto multimodelo
df_consenso <- consolidar_consenso(
  df_predictions = df_eval,
  model_cols = modelos_votantes,
  min_votes = opt$min_votes
)

cat("   -> Resumen del Veredicto de Consenso:\n")
print(table(df_consenso$Veredicto_Consenso))

# Evaluacion automatica de Transferibilidad Filogenetica (si el input incluye Validacion_Espacial)
if ("Validacion_Espacial" %in% colnames(df_consenso) &&
    length(unique( stats::na.omit(df_consenso$Validacion_Espacial) )) >= 2) {
  cat("\n   [Transferibilidad] Detectada columna Ground Truth ('Validacion_Espacial'). Evaluando desempeno...\n")
  obs_gt <- ifelse(df_consenso$Validacion_Espacial %in% c("TRUE_POSITIVE", "TP"), "TP", "FP")
  transfer_rows <- list()

  for (col_v in modelos_votantes) {
    transfer_rows[[col_v]] <- compute_classification_metrics(
      pred = df_consenso[[col_v]],
      obs = obs_gt,
      positive_class = "TP",
      model_name = col_v
    )
  }
  pred_cons_bin <- ifelse(df_consenso$Veredicto_Consenso == "Homologo_Valido", "TP", "FP")
  transfer_rows[["Consenso_7_Modelos"]] <- compute_classification_metrics(
    pred = pred_cons_bin,
    obs = obs_gt,
    positive_class = "TP",
    model_name = "Consenso_7_Modelos",
    score = df_consenso$Votos_TP
  )
  df_transfer <- dplyr::bind_rows(transfer_rows) %>% dplyr::arrange(dplyr::desc(ROC_AUC), dplyr::desc(F1_Score))
  utils::write.table(
    df_transfer,
    file = file.path(opt$out_dir, "transferability_benchmark.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  cat("   -> Tabla de transferibilidad exportada en: transferability_benchmark.tsv\n")
}

# 4. Control de Fragmentacion Espacial y Extraccion de Matriz CNV
cat("[4/7] Ejecutando control de fragmentacion tecnica y extraccion de CNV...\n")

res_cnv <- extract_cnv_matrix(
  df = df_consenso,
  sample_col = "Sample",
  locus_col = "Locus_Ortogrupo",
  cov_threshold = opt$cov_threshold,
  verdict_col = "Veredicto_Consenso"
)

df_matriz_cnv  <- res_cnv$cnv_matrix
df_cnv_anotado <- res_cnv$cnv_annotated
df_funcionales <- res_cnv$functional_genes

cat(sprintf("   -> Loci/Familias genicas detectadas : %d\n", nrow(df_matriz_cnv)))
cat(sprintf("   -> Copias funcionales confirmadas   : %d\n", nrow(df_funcionales)))
cat("   -> Desglose de clasificacion de copias:\n")
print(table(df_cnv_anotado$Clasificacion_Copia))

# 5. Exportacion de Archivos y Tablas
cat("[5/7] Exportando matrices y catalogos a disco...\n")

file_matriz_tsv   <- file.path(opt$out_dir, "cnv_matrix.tsv")
file_anotado_tsv  <- file.path(opt$out_dir, "cnv_annotated_sequences.tsv")
file_func_tsv     <- file.path(opt$out_dir, "functional_genes_catalog.tsv")
file_consenso_tsv <- file.path(opt$out_dir, "consensus_verdicts.tsv")

utils::write.table(df_matriz_cnv, file = file_matriz_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(df_cnv_anotado, file = file_anotado_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(df_funcionales, file = file_func_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(df_consenso, file = file_consenso_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# 6. Generacion y Guardado de Graficos Diagnosticos
cat("[6/7] Renderizando graficos de publicacion (ggplot2)...\n")

p_consenso <- plot_consensus_distribution(
  df_consensus = df_consenso,
  verdict_col = "Veredicto_Consenso"
)
tryCatch({
  ggplot2::ggsave(file.path(opt$out_dir, "plot_consensus_distribution.pdf"), plot = p_consenso, width = 7, height = 5)
  ggplot2::ggsave(file.path(opt$out_dir, "plot_consensus_distribution.png"), plot = p_consenso, width = 7, height = 5, dpi = 300)
}, error = function(e) {
  warning(sprintf("Aviso al guardar graficos de consenso: %s", e$message))
})

if (nrow(df_funcionales) > 0) {
  p_heatmap <- plot_cnv_heatmap(
    df_cnv_annotated = df_cnv_anotado,
    sample_col = "Sample",
    locus_col = "Locus_Ortogrupo",
    top_n_loci = opt$top_loci
  )
  tryCatch({
    ggplot2::ggsave(file.path(opt$out_dir, "plot_cnv_heatmap.pdf"), plot = p_heatmap, width = 9, height = 7)
    ggplot2::ggsave(file.path(opt$out_dir, "plot_cnv_heatmap.png"), plot = p_heatmap, width = 9, height = 7, dpi = 300)
  }, error = function(e) {
    warning(sprintf("Aviso al guardar graficos de heatmap: %s", e$message))
  })
}

# 7. Renderizado del Reporte HTML de Inferencia (report_template.qmd)
cat("[7/7] Renderizando reporte HTML interactivo de inferencia con Quarto...\n")
qmd_inf <- file.path(project_root, "report_template.qmd")
if (!opt$render_report) {
  cat("   -> Renderizado omitido (--no_render_report).\n")
} else if (!file.exists(qmd_inf)) {
  cat(sprintf("   -> [AVISO] Renderizado omitido: no existe la plantilla %s\n", qmd_inf))
} else if (!nzchar(Sys.which("quarto"))) {
  cat("   -> [AVISO] Renderizado omitido: 'quarto' no esta instalado o no esta en el PATH.\n")
} else {
  out_html <- normalizePath(opt$report_file, mustWork = FALSE)
  dir.create(dirname(out_html), recursive = TRUE, showWarnings = FALSE)
  rendered_in_root <- file.path(dirname(qmd_inf), basename(out_html))

  cmd_quarto <- sprintf(
    "quarto render %s --output %s -P project_root:%s -P out_dir:%s -P consensus_verdicts:%s -P cnv_matrix:%s -P cnv_annotated:%s -P functional_genes:%s -P top_loci:%d",
    shQuote(qmd_inf),
    shQuote(basename(out_html)),
    shQuote(project_root),
    shQuote(normalizePath(opt$out_dir, mustWork = FALSE)),
    shQuote(normalizePath(file_consenso_tsv, mustWork = FALSE)),
    shQuote(normalizePath(file_matriz_tsv, mustWork = FALSE)),
    shQuote(normalizePath(file_anotado_tsv, mustWork = FALSE)),
    shQuote(normalizePath(file_func_tsv, mustWork = FALSE)),
    opt$top_loci
  )
  cat(sprintf("   -> Ejecutando Quarto: %s\n", cmd_quarto))
  res_q <- system(cmd_quarto)
  if (res_q == 0) {
    if (file.exists(rendered_in_root) && normalizePath(rendered_in_root, mustWork = FALSE) != out_html) {
      file.copy(rendered_in_root, out_html, overwrite = TRUE)
      unlink(rendered_in_root)
    }
    cat(sprintf("   -> [OK] Reporte HTML de inferencia generado: %s\n", out_html))
  } else {
    cat(sprintf("   -> [AVISO] Fallo el renderizado con Quarto (codigo de retorno: %d).\n", res_q))
  }
}

elapsed_sec <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

cat("========================================================\n")
cat(sprintf("[OK] Matriz CNV y catalogos generados en: %s (%.2f s)\n", opt$out_dir, elapsed_sec))
cat("========================================================\n")
