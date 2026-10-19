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

cat("=================================================================\n")
cat(" TargetCNV-ML: Modulo 02 - Entrenamiento y Validacion Cruzada    \n")
cat("=================================================================\n")
cat(sprintf("Archivo de datos         : %s\n", opt$data))
cat(sprintf("Directorio modelos (.rds): %s\n", opt$out_models))
cat(sprintf("Folds CV                 : %d\n", opt$cv_folds))
cat(sprintf("Balanceo de Clases       : %s\n", opt$sampling))
cat(sprintf("Cores Paralelos          : %d\n", opt$cores))
cat(sprintf("Semilla                  : %d\n", opt$seed))
cat("-----------------------------------------------------------------\n")

# 1. Carga de datos
cat("[1/6] Cargando datos preparados...\n")
if (grepl("\\.rds$", opt$data, ignore.case = TRUE)) {
  datos_raw <- readRDS(opt$data)
} else {
  datos_raw <- utils::read.delim(opt$data, stringsAsFactors = FALSE)
}

metricas_base <- c("pct_recovered", "pct_identity", "score", "wscore")
predictores <- if (opt$include_hit && "hit" %in% colnames(datos_raw)) c(metricas_base, "hit") else metricas_base

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

# 2. Evaluacion de Cortes Univariados Youden
cat("[2/6] Calculando umbrales univariados optimos de Youden...\n")
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
cat("[3/6] Optimizando ensamble heuristico (voto ponderado)...\n")
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

# 4. Configuracion de Caret y Paralelizacion
cat("[4/6] Configurando validacion cruzada k-fold con balanceo de clases...\n")
if (opt$cores > 1) {
  suppressPackageStartupMessages(library(doParallel))
  cl <- parallel::makeCluster(opt$cores)
  doParallel::registerDoParallel(cl)
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
  }, add = TRUE)
}

fit_control <- caret::trainControl(
  method = "cv",
  number = opt$cv_folds,
  classProbs = TRUE,
  summaryFunction = caret::twoClassSummary,
  allowParallel = (opt$cores > 1),
  savePredictions = "final",
  sampling = if (opt$sampling != "none") opt$sampling else NULL
)

# 5. Definicion y Entrenamiento de Modelos
cat("[5/6] Entrenando torneo de modelos supervisados...\n")

model_specs <- list(
  ElasticNet = list(method = "glmnet", tuneLength = 5),
  NaiveBayes = list(method = "nb", tuneLength = 2),
  SVM_Polinomial = list(method = "svmPoly", tuneLength = 3),
  Red_Neuronal = list(method = "nnet", tuneLength = 3, trace = FALSE, MaxNWts = 5000),
  LDA = list(method = "lda")
)

x_train <- df_ml[, predictores, drop = FALSE]
y_train <- df_ml$Clase

trained_models <- list()
benchmark_rows <- list()

# Agregar fila heuristica al benchmark
benchmark_rows[["Heuristico"]] <- heuristic_res$benchmark_row

for (mod_name in names(model_specs)) {
  spec <- model_specs[[mod_name]]
  cat(sprintf("  -> Entrenando: %s ... ", mod_name))
  
  args_train <- list(
    x = x_train,
    y = y_train,
    method = spec$method,
    trControl = fit_control,
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
    trained_models[[mod_name]] <- mod_fit
    saveRDS(mod_fit, file = file.path(opt$out_models, sprintf("model_%s.rds", tolower(mod_name))))
    
    # Extraer metricas
    pred_cv <- mod_fit$pred
    if (!is.null(pred_cv)) {
      cm_mod <- caret::confusionMatrix(data = pred_cv$pred, reference = pred_cv$obs, mode = "everything")
      roc_mod <- suppressMessages(pROC::roc(response = pred_cv$obs, predictor = pred_cv$TP, quiet = TRUE))
      prec <- as.numeric(cm_mod$byClass["Precision"])
      benchmark_rows[[mod_name]] <- data.frame(
        Modelo = mod_name,
        Exactitud = as.numeric(cm_mod$overall["Accuracy"]),
        Sensibilidad = as.numeric(cm_mod$byClass["Sensitivity"]),
        Especificidad = as.numeric(cm_mod$byClass["Specificity"]),
        Precision = prec,
        FDR = 1 - prec,
        F1_Score = as.numeric(cm_mod$byClass["F1"]),
        ROC_AUC = as.numeric(roc_mod$auc),
        stringsAsFactors = FALSE
      )
    }
    cat("OK!\n")
  }
}

# Guardar lista completa
saveRDS(trained_models, file = file.path(opt$out_models, "all_trained_models.rds"))

# 6. Consolidacion de Resultados
cat("[6/6] Guardando metricas y resumen de desempeno...\n")
df_benchmark <- dplyr::bind_rows(benchmark_rows)
out_bench_tsv <- file.path(opt$out_models, "benchmark_models.tsv")
utils::write.table(df_benchmark, file = out_bench_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(df_benchmark, file = file.path(opt$out_models, "benchmark_models.rds"))

# Guardar metricas Youden detalladas
df_youden_bench <- youden_eval$benchmark
df_youden_bench$ROC_AUC <- sapply(metricas_base, function(v) as.numeric(youden_eval$roc_objects[[v]]$auc))
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

cm_ml <- list()
for (m_name in names(trained_models)) {
  p <- trained_models[[m_name]]$pred
  if (!is.null(p)) {
    cm_ml[[m_name]] <- table(Prediccion = p$pred, Referencia = p$obs)
  }
}

confusion_matrices <- list(
  youden = cm_youden,
  heuristic = cm_heur,
  ml = cm_ml
)
saveRDS(confusion_matrices, file = file.path(opt$out_models, "confusion_matrices.rds"))

cat("-----------------------------------------------------------------\n")
cat(" BENCHMARK COMPARATIVO FINAL:\n")
cat("-----------------------------------------------------------------\n")
print(df_benchmark %>% dplyr::arrange(dplyr::desc(F1_Score)))
cat("=================================================================\n")
cat(sprintf("[OK] Todos los modelos han sido serializados en: %s\n", opt$out_models))
