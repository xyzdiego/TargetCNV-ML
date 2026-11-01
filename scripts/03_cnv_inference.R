#!/usr/bin/env Rscript

# Módulo 03: Inferencia de variación en el número de copias (CNV), votación de
# consenso multimodelo, calibración empírica del umbral de cobertura y concordancia

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
source(file.path(project_root, "R", "model_specs.R"))
source(file.path(project_root, "R", "data_processing.R"))
source(file.path(project_root, "R", "model_evaluation.R"))
source(file.path(project_root, "R", "plotting.R"))
source(file.path(project_root, "R", "report_helpers.R"))

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = NULL,
    help = "Ruta al archivo TSV/CSV/RDS con estadísticas de extracción de CAPTUS [Obligatorio]",
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
    help = "Directorio de salida para matrices CNV, tablas de concordancia y figuras [Por defecto: results/]",
    metavar = "DIR"
  ),
  make_option(
    c("-c", "--cov_threshold"),
    type = "double",
    default = 50.0,
    help = "Umbral de cobertura individual (%) para confirmar copias homólogas adicionales [Por defecto: 50.0]",
    metavar = "NUM"
  ),
  make_option(
    c("-f", "--frag_threshold"),
    type = "double",
    default = NULL,
    help = "[Obsoleto/Ignorado] Conservado por compatibilidad con llamadas anteriores",
    metavar = "NUM"
  ),
  make_option(
    c("-v", "--min_votes"),
    type = "integer",
    default = 4L,
    help = "Mínimo de votos afirmativos para el consenso (mayoría simple sobre 7 modelos) [Por defecto: 4]",
    metavar = "INT"
  ),
  make_option(
    c("-t", "--top_loci"),
    type = "integer",
    default = 30L,
    help = "Número de ortogrupos más abundantes para representar en el mapa de calor [Por defecto: 30]",
    metavar = "INT"
  ),
  make_option(
    c("--true_cnv"),
    type = "character",
    default = NULL,
    help = "Archivo TSV opcional con conteos reales de copias (columnas: Sample, Locus_Ortogrupo, CNV_Real) para evaluar concordancia",
    metavar = "FILE"
  ),
  make_option(
    c("--render_report"),
    action = "store_true",
    default = TRUE,
    help = "Compilar el reporte HTML de inferencia (report_template.qmd) [Por defecto: TRUE]"
  ),
  make_option(
    c("--no_render_report"),
    action = "store_false",
    dest = "render_report",
    help = "Omitir la compilación del reporte HTML de inferencia"
  ),
  make_option(
    c("-r", "--report_file"),
    type = "character",
    default = file.path(project_root, "reports", "reporte_inferencia.html"),
    help = "Ruta del reporte HTML de inferencia [Por defecto: reports/reporte_inferencia.html]",
    metavar = "FILE"
  ),
  make_option(
    c("--sources"),
    type = "character",
    default = "",
    help = "Ruta opcional al archivo TSV con citas de fuentes genómicas y de homología",
    metavar = "FILE"
  )
)

opt_parser <- OptionParser(
  usage = "Uso: %prog [opciones]",
  option_list = option_list,
  description = "Aplica el consenso de 7 modelos (clasificadores supervisados + cortes de Youden + heurístico ponderado), estima la matriz CNV por ortogrupo, calibra la sensibilidad de --cov_threshold (30%-80%) y evalúa concordancia de copias cuando hay referencia."
)

opt <- parse_args(opt_parser)

if (is.null(opt$input) || !file.exists(opt$input)) {
  print_help(opt_parser)
  stop(sprintf("[ERROR] Archivo de entrada no especificado o inexistente: %s", opt$input), call. = FALSE)
}
if (!dir.exists(opt$models_dir)) {
  stop(sprintf("[ERROR] El directorio de modelos no existe: %s", opt$models_dir), call. = FALSE)
}
if (!is.null(opt$true_cnv) && nzchar(opt$true_cnv) && !file.exists(opt$true_cnv)) {
  stop(sprintf("[ERROR] No se encontró el archivo --true_cnv especificado: %s", opt$true_cnv), call. = FALSE)
}
if (!dir.exists(opt$out_dir)) {
  dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
}

cat("========================================================\n")
cat(" TargetCNV-ML: Módulo 03 - Inferencia de CNV y Consenso Multimodelo\n")
cat("========================================================\n")
cat(sprintf("Archivo de Entrada (CAPTUS) : %s\n", opt$input))
cat(sprintf("Directorio de Modelos (.rds): %s\n", opt$models_dir))
cat(sprintf("Directorio de Salida        : %s\n", opt$out_dir))
cat(sprintf("Umbral Cobertura Homólogo   : %.1f%%\n", opt$cov_threshold))
cat(sprintf("Votos Mínimos Consenso      : %d / 7\n", opt$min_votes))
cat("--------------------------------------------------------\n")

# 1. Carga de modelos y umbrales serializados
cat("[1/8] Cargando clasificadores serializados y umbrales optimizados...\n")
youden_file    <- file.path(opt$models_dir, "youden_thresholds.rds")
heuristic_file <- file.path(opt$models_dir, "heuristic_model.rds")

if (!file.exists(youden_file)) {
  stop(sprintf("[ERROR] No se encontró 'youden_thresholds.rds' en %s. Ejecute primero ./targetcnv.sh train", opt$models_dir), call. = FALSE)
}
if (!file.exists(heuristic_file)) {
  stop(sprintf("[ERROR] No se encontró 'heuristic_model.rds' en %s. Ejecute primero ./targetcnv.sh train", opt$models_dir), call. = FALSE)
}

youden_thresholds <- readRDS(youden_file)
heuristic_model   <- readRDS(heuristic_file)

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
    cat(sprintf("   -> Clasificador cargado: %s\n", m_name))
  }
}

if (length(ml_models) == 0 && file.exists(file.path(opt$models_dir, "all_trained_models.rds"))) {
  ml_models <- normalize_model_names(readRDS(file.path(opt$models_dir, "all_trained_models.rds")))
  cat(sprintf("   -> Cargados %d modelos desde all_trained_models.rds\n", length(ml_models)))
}
if (length(ml_models) == 0) {
  stop("[ERROR] No se encontró ningún modelo supervisado serializado en --models_dir.", call. = FALSE)
}

# 2. Lectura y estandarización de datos de entrada
cat("[2/8] Cargando y estandarizando estadísticas de CAPTUS...\n")
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
    "[ERROR] El archivo de entrada no cumple el formato de 'captus-extract_stats.tsv'. Faltan columnas: %s",
    paste(columnas_faltantes, collapse = ", ")
  ), call. = FALSE)
}

metricas_numericas <- c("pct_recovered", "pct_identity", "score", "wscore")
spec_file <- file.path(opt$models_dir, "consensus_spec.rds")
consensus_spec <- if (file.exists(spec_file)) readRDS(spec_file) else NULL

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

cat(sprintf("   -> Secuencias NUC analizadas: %d procedentes de %d muestra(s) (Predictores: %s)\n",
            nrow(df_inferencia), length(unique(df_inferencia$Sample)), paste(predictores_ml, collapse = ", ")))

# 3. Evaluación del consenso de 7 modelos
cat("[3/8] Aplicando votación de consenso de 7 modelos...\n")
df_eval <- df_inferencia
modelos_votantes <- character()

voters_objetivo <- if (!is.null(consensus_spec) && !is.null(consensus_spec$ml_voters)) {
  consensus_spec$ml_voters
} else {
  c("ElasticNet", "NaiveBayes", "SVM_Polinomial", "Red_Neuronal")
}
preferidos_ml <- intersect(voters_objetivo, names(ml_models))
top_ml_names  <- if (length(preferidos_ml) == 4L) preferidos_ml else utils::head(unique(c(preferidos_ml, names(ml_models))), 4L)

if (length(top_ml_names) < 4L) {
  warning(sprintf(
    "Solo hay %d clasificador(es) ML disponible(s) (se esperan 4): el consenso operará con %d votantes.",
    length(top_ml_names), length(top_ml_names) + 3L
  ), call. = FALSE)
}

df_x_pred <- as.data.frame(df_inferencia[, predictores_ml, drop = FALSE])

for (mod_name in top_ml_names) {
  col_name <- paste0("Pred_ML_", mod_name)
  mod_obj  <- ml_models[[mod_name]]
  if (inherits(mod_obj, "caretStack")) {
    pr_s <- predict_stacking_probs(mod_obj, newdata = df_x_pred)
    pred_raw <- ifelse(pr_s >= 0.5, "TP", "FP")
  } else {
    pred_raw <- stats::predict(mod_obj, newdata = df_x_pred)
  }
  df_eval[[col_name]] <- factor(
    ifelse(as.character(pred_raw) %in% c("TP", "TRUE_POSITIVE"), "TP", "FP"),
    levels = c("TP", "FP")
  )
  modelos_votantes <- c(modelos_votantes, col_name)
  cat(sprintf("   [Voto ML] Evaluado %s\n", col_name))
}

top2_youden_vars <- c("pct_recovered", "pct_identity")
for (y_var in top2_youden_vars) {
  col_name <- paste0("Pred_Youden_", y_var)
  umbral_y <- youden_thresholds[[y_var]]
  df_eval[[col_name]] <- factor(
    ifelse(df_eval[[y_var]] >= umbral_y, "TP", "FP"),
    levels = c("TP", "FP")
  )
  modelos_votantes <- c(modelos_votantes, col_name)
  cat(sprintf("   [Voto Youden] Evaluado %s (umbral >= %.2f)\n", y_var, umbral_y))
}

col_heur <- "Pred_Ens_Heuristico"
df_eval[[col_heur]] <- predict_weighted_heuristic(
  df_predictors = df_eval,
  thresholds = youden_thresholds,
  weights = heuristic_model$weights,
  global_threshold = heuristic_model$global_threshold
)
modelos_votantes <- c(modelos_votantes, col_heur)
cat(sprintf("   [Voto Heurístico] Evaluado ensamble ponderado (umbral >= %.2f)\n",
            heuristic_model$global_threshold))

if (is.na(opt$min_votes) || opt$min_votes < 1L || opt$min_votes > length(modelos_votantes)) {
  stop(sprintf(
    "[ERROR] --min_votes (%s) debe estar entre 1 y el total de modelos votantes disponibles (%d).",
    as.character(opt$min_votes), length(modelos_votantes)
  ), call. = FALSE)
}

df_consenso <- consolidar_consenso(
  df_predictions = df_eval,
  model_cols = modelos_votantes,
  min_votes = opt$min_votes
)

cat("   -> Resumen del veredicto de consenso:\n")
print(table(df_consenso$Veredicto_Consenso))

# Evaluación a nivel de secuencia cuando existe Validacion_Espacial en el input
if ("Validacion_Espacial" %in% colnames(df_consenso) &&
    length(unique(stats::na.omit(df_consenso$Validacion_Espacial))) >= 2L) {
  cat("\n   [Validación] Detectada columna 'Validacion_Espacial'. Calculando métricas de clasificación...\n")
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
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  cat("   -> Métricas de clasificación exportadas en: transferability_benchmark.tsv\n")
}

# 4. Control de fragmentación y extracción de la matriz CNV
cat("[4/8] Ejecutando control de cobertura y extracción de matriz CNV...\n")
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

cat(sprintf("   -> Ortogrupos con copias confirmadas : %d\n", nrow(df_matriz_cnv)))
cat(sprintf("   -> Copias génicas funcionales        : %d\n", nrow(df_funcionales)))
cat("   -> Desglose de clasificación de copias:\n")
print(table(df_cnv_anotado$Clasificacion_Copia))

# 5. Calibración empírica del umbral --cov_threshold (30% a 80%) y concordancia CNV por ortogrupo
cat("[5/8] Evaluando sensibilidad del umbral --cov_threshold (30% a 80%) y concordancia por ortogrupo...\n")
true_cnv_df <- if (!is.null(opt$true_cnv) && nzchar(opt$true_cnv)) {
  utils::read.delim(opt$true_cnv, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
} else {
  NULL
}

df_sens <- evaluate_cnv_threshold_sensitivity(
  df_consenso = df_consenso,
  thresholds = seq(30, 80, by = 5),
  sample_col = "Sample",
  locus_col = "Locus_Ortogrupo",
  verdict_col = "Veredicto_Consenso",
  true_cnv_df = true_cnv_df
)
file_sens_tsv <- file.path(opt$out_dir, "cnv_threshold_sensitivity.tsv")
utils::write.table(df_sens, file = file_sens_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

has_cnv_gt <- !is.null(true_cnv_df) || ("Validacion_Espacial" %in% colnames(df_consenso) &&
  length(unique(stats::na.omit(df_consenso$Validacion_Espacial))) >= 1L)

file_conc_tsv <- ""
if (has_cnv_gt) {
  conc_res <- evaluate_cnv_concordance(
    df_consenso = df_consenso,
    cov_threshold = opt$cov_threshold,
    sample_col = "Sample",
    locus_col = "Locus_Ortogrupo",
    verdict_col = "Veredicto_Consenso",
    true_cnv_df = true_cnv_df
  )
  if (nrow(conc_res$summary) > 0) {
    file_conc_tsv  <- file.path(opt$out_dir, "cnv_concordance_summary.tsv")
    file_pairs_tsv <- file.path(opt$out_dir, "cnv_concordance_pairs.tsv")
    utils::write.table(conc_res$summary, file = file_conc_tsv,  sep = "\t", row.names = FALSE, quote = FALSE)
    utils::write.table(conc_res$pairs,   file = file_pairs_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
    s_row <- conc_res$summary[1, ]
    cat(sprintf(
      "   -> Concordancia CNV por ortogrupo (n = %d): Exacta = %.2f%% | Spearman rho = %.4f | Kappa_w = %.4f | RMSE = %.4f\n",
      s_row$N_Ortogrupos, s_row$Concordancia_Exacta * 100,
      s_row$Spearman_Rho, s_row$Kappa_Cuadratico, s_row$RMSE
    ))
  }
}

# 6. Exportación de matrices y catálogos principales
cat("[6/8] Exportando matrices CNV y catálogos de secuencias...\n")
file_matriz_tsv   <- file.path(opt$out_dir, "cnv_matrix.tsv")
file_anotado_tsv  <- file.path(opt$out_dir, "cnv_annotated_sequences.tsv")
file_func_tsv     <- file.path(opt$out_dir, "functional_genes_catalog.tsv")
file_consenso_tsv <- file.path(opt$out_dir, "consensus_verdicts.tsv")

utils::write.table(df_matriz_cnv,  file = file_matriz_tsv,   sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(df_cnv_anotado, file = file_anotado_tsv,  sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(df_funcionales, file = file_func_tsv,     sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(df_consenso,    file = file_consenso_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# 7. Generación de gráficos diagnósticos
cat("[7/8] Generando figuras diagnósticas (ggplot2)...\n")
p_consenso <- plot_consensus_distribution(
  df_consensus = df_consenso,
  verdict_col = "Veredicto_Consenso"
)
tryCatch({
  ggplot2::ggsave(file.path(opt$out_dir, "plot_consensus_distribution.pdf"), plot = p_consenso, width = 7, height = 5)
  ggplot2::ggsave(file.path(opt$out_dir, "plot_consensus_distribution.png"), plot = p_consenso, width = 7, height = 5, dpi = 300)
}, error = function(e) {
  warning(sprintf("Aviso al guardar gráficos de consenso: %s", e$message))
})

p_sens <- plot_cnv_threshold_sensitivity(df_sens, active_threshold = opt$cov_threshold)
tryCatch({
  ggplot2::ggsave(file.path(opt$out_dir, "plot_cnv_threshold_sensitivity.pdf"), plot = p_sens, width = 7.5, height = 5)
  ggplot2::ggsave(file.path(opt$out_dir, "plot_cnv_threshold_sensitivity.png"), plot = p_sens, width = 7.5, height = 5, dpi = 300)
}, error = function(e) {
  warning(sprintf("Aviso al guardar curva de sensibilidad de cobertura: %s", e$message))
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
    warning(sprintf("Aviso al guardar mapa de calor de CNV: %s", e$message))
  })
}

# 8. Compilación del reporte HTML de inferencia en modo de solo lectura
cat("[8/8] Compilando reporte HTML de inferencia con Quarto...\n")
qmd_inf <- file.path(project_root, "report_template.qmd")
if (!opt$render_report) {
  cat("   -> Compilación omitida (--no_render_report).\n")
} else if (!file.exists(qmd_inf)) {
  cat(sprintf("   -> [AVISO] Compilación omitida: no existe la plantilla %s\n", qmd_inf))
} else if (!nzchar(Sys.which("quarto"))) {
  cat("   -> [AVISO] Compilación omitida: 'quarto' no está disponible en el PATH.\n")
} else {
  out_html <- normalizePath(opt$report_file, mustWork = FALSE)
  dir.create(dirname(out_html), recursive = TRUE, showWarnings = FALSE)
  rendered_in_root <- file.path(dirname(qmd_inf), basename(out_html))

  src_arg <- if (!is.null(opt$sources) && nzchar(trimws(opt$sources))) {
    sprintf(" -P sources_file:%s", shQuote(normalizePath(opt$sources, mustWork = FALSE)))
  } else {
    ""
  }
  conc_arg <- if (nzchar(file_conc_tsv) && file.exists(file_conc_tsv)) {
    sprintf(" -P cnv_concordance:%s", shQuote(normalizePath(file_conc_tsv, mustWork = FALSE)))
  } else {
    ""
  }

  cmd_quarto <- sprintf(
    "quarto render %s --output %s -P project_root:%s -P out_dir:%s -P consensus_verdicts:%s -P cnv_matrix:%s -P cnv_annotated:%s -P functional_genes:%s -P cnv_sensitivity:%s -P cov_threshold:%s -P top_loci:%d%s%s",
    shQuote(qmd_inf),
    shQuote(basename(out_html)),
    shQuote(project_root),
    shQuote(normalizePath(opt$out_dir, mustWork = FALSE)),
    shQuote(normalizePath(file_consenso_tsv, mustWork = FALSE)),
    shQuote(normalizePath(file_matriz_tsv, mustWork = FALSE)),
    shQuote(normalizePath(file_anotado_tsv, mustWork = FALSE)),
    shQuote(normalizePath(file_func_tsv, mustWork = FALSE)),
    shQuote(normalizePath(file_sens_tsv, mustWork = FALSE)),
    as.character(opt$cov_threshold),
    opt$top_loci,
    conc_arg,
    src_arg
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
    cat(sprintf("   -> [AVISO] Falló el renderizado con Quarto (código de retorno: %d).\n", res_q))
  }
}

elapsed_sec <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

cat("========================================================\n")
cat(sprintf("[OK] Matriz CNV y tablas exportadas en: %s (%.2f s)\n", opt$out_dir, elapsed_sec))
cat("========================================================\n")
