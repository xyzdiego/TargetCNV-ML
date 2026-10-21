#!/usr/bin/env Rscript

# ==============================================================================
# TargetCNV-ML Pipeline
# Script 02: Model Training and Cross-Validation with Caret
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
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
source(file.path(project_root, "R", "model_evaluation.R"))
source(file.path(project_root, "R", "data_processing.R"))
source(file.path(project_root, "R", "plotting.R"))

option_list <- list(
  make_option(
    c("-d", "--data"),
    type = "character",
    default = file.path(project_root, "data", "data_prepared.rds"),
    help = "Ruta al archivo de datos preparados (.rds o .tsv) [Obligatorio]",
    metavar = "FILE"
  ),
  make_option(
    c("-m", "--out_models"),
    type = "character",
    default = file.path(project_root, "models"),
    help = "Directorio de destino para guardar los modelos serializados (.rds) [Por defecto: models/]",
    metavar = "DIR"
  ),
  make_option(
    c("-k", "--cv_folds"),
    type = "integer",
    default = 5,
    help = "Numero de particiones (folds) para validacion cruzada k-fold [Por defecto: 5]",
    metavar = "INT"
  ),
  make_option(
    c("-s", "--sampling"),
    type = "character",
    default = "down",
    help = "Estrategia de balanceo de clases: 'down', 'up' o 'none' [Por defecto: down]",
    metavar = "STR"
  ),
  make_option(
    c("--seed"),
    type = "integer",
    default = 1234,
    help = "Semilla de reproducibilidad estocastica [Por defecto: 1234]",
    metavar = "INT"
  ),
  make_option(
    c("-c", "--cores"),
    type = "integer",
    default = max(1, parallel::detectCores() - 1),
    help = "Numero de hilos para entrenamiento paralelo [Por defecto: cores - 1]",
    metavar = "INT"
  ),
  make_option(
    c("--include_hit"),
    action = "store_true",
    default = FALSE,
    help = "Incluir el indice 'hit' de CAPTUS como predictor continuo"
  ),
  make_option(
    c("--render_report"),
    action = "store_true",
    default = TRUE,
    help = "Compilar automaticamente el reporte HTML interactivo tras el entrenamiento [Por defecto: TRUE]"
  ),
  make_option(
    c("--no_render_report"),
    action = "store_false",
    dest = "render_report",
    help = "Omitir la compilacion del reporte HTML interactivo"
  ),
  make_option(
    c("-r", "--report_file"),
    type = "character",
    default = file.path(project_root, "reporte_entrenamiento.html"),
    help = "Ruta del reporte HTML de salida [Por defecto: reporte_entrenamiento.html en la raiz del proyecto]",
    metavar = "FILE"
  )
)

opt_parser <- OptionParser(
  usage = "Uso: %prog [opciones]",
  option_list = option_list,
  description = "Entrena modelos supervisados (ElasticNet, NaiveBayes, SVM, NNet, etc.), optimiza Youden y guarda artefactos .rds."
)

opt <- parse_args(opt_parser)

if (is.null(opt$data) || !file.exists(opt$data)) {
  # Si el default .rds no existe, buscar el .tsv
  tsv_fallback <- sub("\\.rds$", ".tsv", opt$data)
  if (file.exists(tsv_fallback)) {
    opt$data <- tsv_fallback
  } else {
    print_help(opt_parser)
    stop(sprintf("[ERROR] No se encontro el archivo de datos preparados: %s", opt$data), call. = FALSE)
  }
}

if (!dir.exists(opt$out_models)) {
  dir.create(opt$out_models, recursive = TRUE, showWarnings = FALSE)
}

# Validacion de argumentos (mensajes claros antes de iniciar un entrenamiento largo)
if (!opt$sampling %in% c("down", "up", "none")) {
  stop(sprintf("[ERROR] --sampling debe ser 'down', 'up' o 'none' (recibido: '%s').", opt$sampling), call. = FALSE)
}
if (is.na(opt$cv_folds) || opt$cv_folds < 2) {
  stop("[ERROR] --cv_folds debe ser un entero >= 2.", call. = FALSE)
}
if (is.na(opt$cores) || opt$cores < 1) {
  stop("[ERROR] --cores debe ser un entero >= 1.", call. = FALSE)
}
if (opt$sampling == "none") {
  cat("[AVISO] --sampling none: el escenario 'Balanceado' sigue usando down-sampling\n")
  cat("        (el escenario 'Desbalanceado' nunca remuestrea). Use 'up' para sobre-muestreo.\n")
}
if (opt$include_hit) {
  cat("[AVISO] --include_hit agrega 'hit' como predictor: los modelos resultantes ya no son\n")
  cat("        compatibles con 03_cnv_inference.R, que solo evalua las 4 metricas base.\n")
}

cat("=================================================================\n")
cat(" TargetCNV-ML: Modulo 02 - Entrenamiento y Validacion Cruzada    \n")
cat("=================================================================\n")
cat(sprintf("Archivo de datos         : %s\n", opt$data))
cat(sprintf("Directorio modelos (.rds): %s\n", opt$out_models))
cat(sprintf("Folds CV                 : %d\n", opt$cv_folds))
cat(sprintf("Balanceo de Clases       : %s\n", opt$sampling))
cat(sprintf("Cores Paralelos          : %d\n", opt$cores))
cat(sprintf("Semilla                  : %d\n", opt$seed))
cat(sprintf("Compilar Reporte HTML    : %s\n", if (opt$render_report) "SI" else "NO"))
if (opt$render_report) cat(sprintf("Reporte HTML de salida   : %s\n", opt$report_file))
cat("-----------------------------------------------------------------\n")

# 1. Carga de datos
cat("[1/7] Cargando datos preparados...\n")
if (grepl("\\.rds$", opt$data, ignore.case = TRUE)) {
  datos_raw <- readRDS(opt$data)
} else {
  datos_raw <- utils::read.delim(opt$data, stringsAsFactors = FALSE)
}

metricas_base <- c("pct_recovered", "pct_identity", "score", "wscore")
predictores <- if (opt$include_hit && "hit" %in% colnames(datos_raw)) c(metricas_base, "hit") else metricas_base

columnas_faltantes <- setdiff(c("Validacion_Espacial", predictores), colnames(datos_raw))
if (length(columnas_faltantes) > 0) {
  stop(sprintf(
    "[ERROR] El archivo de datos no contiene las columnas requeridas: %s.\n        Genere el archivo con scripts/01_prepare_data.R (con --dict_ref) o revise el formato.",
    paste(columnas_faltantes, collapse = ", ")
  ), call. = FALSE)
}

df_ml <- datos_raw %>%
  dplyr::filter(!is.na(Validacion_Espacial)) %>%
  dplyr::mutate(
    Clase = factor(
      dplyr::if_else(Validacion_Espacial == "TRUE_POSITIVE", "TP", "FP"),
      levels = c("TP", "FP")
    ),
    dplyr::across(dplyr::all_of(predictores), ~ as.numeric(as.character(.x)))
  ) %>%
  dplyr::select(Clase, dplyr::all_of(predictores)) %>%
  stats::na.omit() %>%
  as.data.frame()

cat(sprintf("  -> Muestras validas: %d (TP = %d, FP = %d)\n",
            nrow(df_ml), sum(df_ml$Clase == "TP"), sum(df_ml$Clase == "FP")))

if (min(table(df_ml$Clase)) < opt$cv_folds) {
  stop(sprintf(
    "[ERROR] La clase minoritaria tiene %d muestras (< %d folds). Se requieren TP y FP en el archivo de datos\n        (¿01_prepare_data.R se ejecuto sin --dict_ref?) o reduzca --cv_folds.",
    min(table(df_ml$Clase)), opt$cv_folds
  ), call. = FALSE)
}

# 2. Evaluacion de Cortes Univariados Youden
cat("[2/7] Calculando umbrales univariados optimos de Youden...\n")
youden_eval <- evaluate_youden_predictors(
  df = df_ml,
  response_col = "Clase",
  metric_cols = metricas_base,
  positive_class = "TP"
)
saveRDS(youden_eval$thresholds, file = file.path(opt$out_models, "youden_thresholds.rds"))
cat("  -> Umbrales optimizados:\n")
print(round(youden_eval$thresholds, 2))

# 3. Optimizacion de Ensamble Heuristico Ponderado
cat("[3/7] Optimizando ensamble heuristico (voto ponderado)...\n")
heuristic_res <- optimize_weighted_heuristic(
  df_votes = youden_eval$binary_votes,
  response = df_ml$Clase,
  metric_cols = metricas_base,
  positive_class = "TP"
)
heuristic_obj <- list(
  weights = heuristic_res$weights,
  global_threshold = heuristic_res$global_threshold,
  thresholds = youden_eval$thresholds,
  best_auc = heuristic_res$best_auc
)
saveRDS(heuristic_obj, file = file.path(opt$out_models, "heuristic_model.rds"))
cat("  -> Pesos heurísticos optimizados:\n")
print(heuristic_res$weights)
cat(sprintf("  -> Umbral global: %.2f | AUC: %.4f\n", heuristic_res$global_threshold, heuristic_res$best_auc))

# 4. Configuracion de Paralelizacion
if (opt$cores > 1) {
  suppressPackageStartupMessages(library(doParallel))
  cl <- parallel::makeCluster(opt$cores)
  doParallel::registerDoParallel(cl)
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
  }, add = TRUE)
}

# 5. Definicion de Especificaciones de Modelos
model_specs <- list(
  ElasticNet = list(method = "glmnet", tuneLength = 5),
  NaiveBayes = list(method = "nb", tuneLength = 2),
  SVM_Polinomial = list(method = "svmPoly", tuneLength = 3),
  Red_Neuronal = list(method = "nnet", tuneLength = 3, trace = FALSE, MaxNWts = 5000),
  LDA = list(method = "lda")
)

x_train <- df_ml[, predictores, drop = FALSE]
y_train <- df_ml$Clase

# Funcion auxiliar para entrenar torneo de modelos
train_tournament <- function(control_obj, scenario_label) {
  mods <- list()
  rows <- list()
  for (mod_name in names(model_specs)) {
    spec <- model_specs[[mod_name]]
    cat(sprintf("  [%s] Entrenando: %s ... ", scenario_label, mod_name))
    args_train <- list(
      x = x_train,
      y = y_train,
      method = spec$method,
      trControl = control_obj,
      metric = "ROC",
      preProcess = c("center", "scale")
    )
    if (!is.null(spec$tuneLength)) args_train$tuneLength <- spec$tuneLength
    if (!is.null(spec$tuneGrid))   args_train$tuneGrid   <- spec$tuneGrid
    if (!is.null(spec$trace))      args_train$trace      <- spec$trace
    if (!is.null(spec$MaxNWts))   args_train$MaxNWts   <- spec$MaxNWts

    set.seed(opt$seed)
    mod_fit <- tryCatch({
      do.call(caret::train, args_train)
    }, error = function(e) {
      cat(sprintf("[FALLO: %s]\n", e$message))
      NULL
    })

    if (!is.null(mod_fit)) {
      mods[[mod_name]] <- mod_fit
      pred_cv <- mod_fit$pred
      if (!is.null(pred_cv)) {
        # Metricas sobre predicciones out-of-fold (CV) con el esquema canonico
        rows[[mod_name]] <- compute_classification_metrics(
          pred = pred_cv$pred,
          obs = pred_cv$obs,
          score = pred_cv$TP,
          positive_class = "TP",
          model_name = mod_name
        )
      }
      cat("OK!\n")
    }
  }
  return(list(models = mods, benchmark = dplyr::bind_rows(rows)))
}

# 5A. Entrenar Escenario Desbalanceado (sin remuestreo)
cat("[4/7] Entrenando torneo de modelos (Escenario Desbalanceado)...\n")
fit_control_unbal <- caret::trainControl(
  method = "cv",
  number = opt$cv_folds,
  classProbs = TRUE,
  summaryFunction = caret::twoClassSummary,
  allowParallel = (opt$cores > 1),
  savePredictions = "final",
  sampling = NULL
)
res_unbal <- train_tournament(fit_control_unbal, "Desbalanceado")

# 5B. Entrenar Escenario Balanceado (down-sampling)
cat("[5/7] Entrenando torneo de modelos (Escenario Balanceado)...\n")
fit_control_bal <- caret::trainControl(
  method = "cv",
  number = opt$cv_folds,
  classProbs = TRUE,
  summaryFunction = caret::twoClassSummary,
  allowParallel = (opt$cores > 1),
  savePredictions = "final",
  sampling = if (opt$sampling != "none") opt$sampling else "down"
)
res_bal <- train_tournament(fit_control_bal, "Balanceado")

# Liberar los workers paralelos una vez terminado el entrenamiento
if (opt$cores > 1) {
  try(parallel::stopCluster(cl), silent = TRUE)
  foreach::registerDoSEQ()
}

# El ensamble heuristico no depende del remuestreo: se reporta en ambos escenarios
# para que las dos tablas tengan exactamente las mismas filas y columnas.
df_benchmark_bal   <- dplyr::bind_rows(heuristic_res$benchmark_row, res_bal$benchmark)
df_benchmark_unbal <- dplyr::bind_rows(heuristic_res$benchmark_row, res_unbal$benchmark)
df_benchmark_all   <- dplyr::bind_rows(
  dplyr::mutate(df_benchmark_bal,   Escenario = "Balanceado",    .before = 1),
  dplyr::mutate(df_benchmark_unbal, Escenario = "Desbalanceado", .before = 1)
)

# Serializar artefactos de modelos
saveRDS(res_bal$models, file = file.path(opt$out_models, "models_balanced.rds"))
saveRDS(res_bal$models, file = file.path(opt$out_models, "all_trained_models.rds"))
saveRDS(res_unbal$models, file = file.path(opt$out_models, "models_unbalanced.rds"))

for (mod_name in names(res_bal$models)) {
  saveRDS(res_bal$models[[mod_name]], file = file.path(opt$out_models, sprintf("model_%s.rds", tolower(mod_name))))
}

# Serializar benchmarks
saveRDS(df_benchmark_bal, file = file.path(opt$out_models, "benchmark_models_balanced.rds"))
saveRDS(df_benchmark_bal, file = file.path(opt$out_models, "benchmark_models.rds"))
utils::write.table(df_benchmark_bal, file = file.path(opt$out_models, "benchmark_models.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(df_benchmark_bal, file = file.path(opt$out_models, "benchmark_models_balanced.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

saveRDS(df_benchmark_unbal, file = file.path(opt$out_models, "benchmark_models_unbalanced.rds"))
utils::write.table(df_benchmark_unbal, file = file.path(opt$out_models, "benchmark_models_unbalanced.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(df_benchmark_all, file = file.path(opt$out_models, "benchmark_models_all.rds"))
utils::write.table(df_benchmark_all, file = file.path(opt$out_models, "benchmark_models_all.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Guardar metricas Youden detalladas (mismo esquema de columnas + Umbral e Indice J)
df_youden_bench <- youden_eval$benchmark
df_youden_bench$Umbral <- sapply(metricas_base, function(v) round(youden_eval$thresholds[[v]], 3))
df_youden_bench$Indice_Youden_J <- df_youden_bench$Sensibilidad + df_youden_bench$Especificidad - 1
saveRDS(df_youden_bench, file = file.path(opt$out_models, "youden_benchmark.rds"))
utils::write.table(df_youden_bench, file = file.path(opt$out_models, "youden_benchmark.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Consolidar matrices de confusion
cm_youden <- list()
for (var in metricas_base) {
  pred_y <- factor(ifelse(df_ml[[var]] >= youden_eval$thresholds[[var]], "TP", "FP"), levels = c("TP", "FP"))
  cm_youden[[var]] <- table(Prediccion = pred_y, Referencia = df_ml$Clase)
}

pred_h <- predict_weighted_heuristic(df_ml, youden_eval$thresholds, heuristic_res$weights, heuristic_res$global_threshold)
cm_heur <- table(Prediccion = pred_h, Referencia = df_ml$Clase)

cm_ml_bal <- list()
for (m_name in names(res_bal$models)) {
  p <- res_bal$models[[m_name]]$pred
  if (!is.null(p)) cm_ml_bal[[m_name]] <- table(Prediccion = p$pred, Referencia = p$obs)
}

cm_ml_unbal <- list()
for (m_name in names(res_unbal$models)) {
  p <- res_unbal$models[[m_name]]$pred
  if (!is.null(p)) cm_ml_unbal[[m_name]] <- table(Prediccion = p$pred, Referencia = p$obs)
}

confusion_matrices <- list(
  youden = cm_youden,
  heuristic = cm_heur,
  balanced = cm_ml_bal,
  unbalanced = cm_ml_unbal,
  ml = cm_ml_bal
)
saveRDS(confusion_matrices, file = file.path(opt$out_models, "confusion_matrices.rds"))

# 6. Explicabilidad SHAP (KernelSHAP con ElasticNet Balanceado)
cat("[6/7] Generando explicaciones SHAP (Interpretabilidad con KernelSHAP)...\n")
if ("ElasticNet" %in% names(res_bal$models)) {
  suppressPackageStartupMessages({
    library(kernelshap)
    library(shapviz)
  })
  set.seed(opt$seed)
  bg_idx <- sample(seq_len(nrow(x_train)), min(30, nrow(x_train)))
  shap_calc <- tryCatch({
    kernelshap::kernelshap(
      object = res_bal$models[["ElasticNet"]],
      X = x_train,
      bg_X = x_train[bg_idx, , drop = FALSE],
      pred_fun = function(m, d) predict(m, newdata = d, type = "prob")[, "TP"]
    )
  }, error = function(e) {
    cat(sprintf("  -> [Aviso] Error en calculo SHAP: %s\n", e$message))
    NULL
  })
  if (!is.null(shap_calc)) {
    shap_obj <- shapviz::shapviz(shap_calc)
    saveRDS(shap_obj, file = file.path(opt$out_models, "shap_values.rds"))
    cat("  -> Artefacto SHAP serializado exitosamente.\n")
  }
}

# 7. Renderizado del Reporte HTML
cat("[7/7] Renderizando reporte HTML interactivo con Quarto...\n")
qmd_file <- file.path(project_root, "training_report.qmd")
if (!opt$render_report) {
  cat("  -> Renderizado omitido (--no_render_report).\n")
} else if (!file.exists(qmd_file)) {
  cat(sprintf("  -> [AVISO] Renderizado omitido: no existe la plantilla %s\n", qmd_file))
} else if (!nzchar(Sys.which("quarto"))) {
  cat("  -> [AVISO] Renderizado omitido: 'quarto' no esta instalado o no esta en el PATH.\n")
} else {
  out_html <- normalizePath(opt$report_file, mustWork = FALSE)
  dir.create(dirname(out_html), recursive = TRUE, showWarnings = FALSE)
  # Quarto se ejecuta en el directorio de la plantilla: se entregan rutas absolutas
  cmd_quarto <- sprintf(
    "quarto render \"%s\" -P models_dir:\"%s\" -P data_file:\"%s\" --output-dir \"%s\" --output \"%s\"",
    qmd_file,
    normalizePath(opt$out_models, mustWork = FALSE),
    normalizePath(opt$data, mustWork = FALSE),
    dirname(out_html), basename(out_html)
  )
  cat(sprintf("  -> Ejecutando: %s\n", cmd_quarto))
  status_render <- system(cmd_quarto)
  if (status_render == 0) {
    cat(sprintf("  -> [OK] Reporte HTML compilado exitosamente: %s\n", out_html))
  } else {
    cat(sprintf("  -> [AVISO] Quarto retorno codigo de salida: %d\n", status_render))
  }
}

cat("-----------------------------------------------------------------\n")
options(width = 200)
cat(" BENCHMARK COMPARATIVO FINAL (BALANCEADO):\n")
cat("-----------------------------------------------------------------\n")
print(df_benchmark_bal %>% dplyr::arrange(dplyr::desc(F1_Score)), digits = 3)
cat("-----------------------------------------------------------------\n")
cat(" BENCHMARK COMPARATIVO FINAL (DESBALANCEADO):\n")
cat("-----------------------------------------------------------------\n")
print(df_benchmark_unbal %>% dplyr::arrange(dplyr::desc(F1_Score)), digits = 3)
cat("=================================================================\n")
cat(sprintf("[OK] Todos los modelos han sido serializados en: %s\n", opt$out_models))
