#!/usr/bin/env Rscript

# ==============================================================================
# TargetCNV-ML Pipeline
# Script 02: Model Training (Full 18-Model Tournament + Stacking + Youden + SHAP)
# ==============================================================================

t_start <- Sys.time()

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
if (file.exists(file.path(project_root, "R", "model_specs.R"))) source(file.path(project_root, "R", "model_specs.R"))
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
    c("--extended_tournament"),
    action = "store_true",
    default = TRUE,
    help = "Ejecutar el torneo completo de 18 algoritmos del torneo completo [Por defecto: TRUE]"
  ),
  make_option(
    c("--fast_tournament"),
    action = "store_false",
    dest = "extended_tournament",
    help = "Entrenar solo los 5 modelos rapidos principales (ElasticNet, NaiveBayes, SVM_Polinomial, Red_Neuronal, LDA)"
  ),
  make_option(
    c("--train_stacking"),
    action = "store_true",
    default = TRUE,
    help = "Entrenar meta-modelo GLM Stacking (caretEnsemble) combinando los mejores modelos base [Por defecto: TRUE]"
  ),
  make_option(
    c("--no_stacking"),
    action = "store_false",
    dest = "train_stacking",
    help = "Omitir el entrenamiento del meta-modelo GLM Stacking"
  ),
  make_option(
    c("--clean_models"),
    action = "store_true",
    default = FALSE,
    help = "Limpiar los archivos .rds/.tsv previos en --out_models antes de guardar la nueva corrida [Por defecto: FALSE]"
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
  ),
  make_option(
    c("--sources"),
    type = "character",
    default = "",
    help = "Ruta opcional al archivo TSV con las citas de fuentes genomicas y de homologia para el reporte HTML",
    metavar = "FILE"
  )
)

opt_parser <- OptionParser(
  usage = "Uso: %prog [opciones]",
  option_list = option_list,
  description = "Entrena todos los modelos supervisados del torneo completo (RF, GBM, XGBoost, SVM Radial/Lineal/Polinomial, NNet, Bagged NNet, ElasticNet, KNN, CART, NaiveBayes, LDA, DNN/MLP, QDA, Treebag, Stacking), optimiza Youden y genera el reporte HTML."
)

opt <- parse_args(opt_parser)

if (is.null(opt$data) || !file.exists(opt$data)) {
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

if (!opt$sampling %in% c("down", "up", "none")) {
  stop(sprintf("[ERROR] --sampling debe ser 'down', 'up' o 'none' (recibido: '%s').", opt$sampling), call. = FALSE)
}
if (is.na(opt$cv_folds) || opt$cv_folds < 2) {
  stop("[ERROR] --cv_folds debe ser un entero >= 2.", call. = FALSE)
}
if (is.na(opt$cores) || opt$cores < 1) {
  stop("[ERROR] --cores debe ser un entero >= 1.", call. = FALSE)
}

cat("========================================================\n")
cat(" TargetCNV-ML: Modulo 02 - Entrenamiento y Validacion Cruzada  \n")
cat("========================================================\n")
cat(sprintf("Archivo de datos          : %s\n", opt$data))
cat(sprintf("Directorio modelos (.rds) : %s\n", opt$out_models))
cat(sprintf("Folds CV                  : %d\n", opt$cv_folds))
cat(sprintf("Balanceo de Clases        : %s\n", opt$sampling))
cat(sprintf("Torneo de Modelos         : %s\n", if (opt$extended_tournament) "Completo (18 algoritmos ML)" else "Rapido (5 algoritmos ML)"))
cat(sprintf("Meta-Modelo Stacking      : %s\n", if (opt$train_stacking) "SI" else "NO"))
cat(sprintf("Cores Paralelos           : %d\n", opt$cores))
cat(sprintf("Semilla                   : %d\n", opt$seed))
cat(sprintf("Compilar Reporte HTML     : %s\n", if (opt$render_report) "SI" else "NO"))
if (opt$render_report) cat(sprintf("Reporte HTML de salida    : %s\n", opt$report_file))
cat("--------------------------------------------------------\n")

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
    "[ERROR] El archivo de datos no contiene las columnas requeridas: %s.\n Genere el archivo con scripts/01_prepare_data.R",
    paste(columnas_faltantes, collapse = ", ")
  ), call. = FALSE)
}

df_ml <- datos_raw %>%
  dplyr::filter(!is.na(Validacion_Espacial)) %>%
  dplyr::mutate(
    Clase = factor(
      dplyr::if_else(Validacion_Espacial %in% c("TRUE_POSITIVE", "TP"), "TP", "FP"),
      levels = c("TP", "FP")
    ),
    dplyr::across(dplyr::all_of(predictores), ~ as.numeric(as.character(.x)))
  ) %>%
  dplyr::select(Clase, dplyr::all_of(predictores)) %>%
  stats::na.omit() %>%
  as.data.frame()

cat(sprintf("   -> Muestras validas: %d (TP = %d, FP = %d)\n",
            nrow(df_ml), sum(df_ml$Clase == "TP"), sum(df_ml$Clase == "FP")))

if (min(table(df_ml$Clase)) < opt$cv_folds) {
  stop(sprintf(
    "[ERROR] La clase minoritaria tiene %d muestras (< %d folds). Se requieren TP y FP suficientes en el archivo de datos.\n",
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
cat("   -> Umbrales optimizados:\n")
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
cat("   -> Pesos heuristicos optimizados:\n")
print(heuristic_res$weights)
cat(sprintf("   -> Umbral global: %.2f | AUC: %.4f\n", heuristic_res$global_threshold, heuristic_res$best_auc))

# 4. Configuracion de Paralelizacion
cl <- NULL
if (opt$cores > 1) {
  suppressPackageStartupMessages(library(doParallel))
  cl <- parallel::makeCluster(opt$cores)
  doParallel::registerDoParallel(cl)
}

# 5. Especificacion Programatica de Modelos desde R/model_specs.R
model_specs <- get_model_specs(n_obs = nrow(df_ml), extended = opt$extended_tournament)

x_train <- df_ml[, predictores, drop = FALSE]
y_train <- df_ml$Clase

train_tournament <- function(control_obj, scenario_label) {
  mods <- list()
  rows <- list()
  cms  <- list()

  for (mod_name in names(model_specs)) {
    spec <- model_specs[[mod_name]]
    cat(sprintf("   [%s] Entrenando: %-18s ... ", scenario_label, mod_name))
    args_train <- list(
      x = x_train,
      y = y_train,
      method = spec$method,
      trControl = control_obj,
      metric = "ROC"
    )
    if (is.null(spec$no_scale) || !spec$no_scale) {
      args_train$preProcess <- c("center", "scale")
    }
    if (!is.null(spec$tuneLength)) args_train$tuneLength <- spec$tuneLength
    if (!is.null(spec$tuneGrid))   args_train$tuneGrid   <- spec$tuneGrid
    if (!is.null(spec$trace))      args_train$trace      <- spec$trace
    if (!is.null(spec$MaxNWts))    args_train$MaxNWts    <- spec$MaxNWts
    if (!is.null(spec$verbose))    args_train$verbose    <- spec$verbose
    if (!is.null(spec$nthread))    args_train$nthread    <- spec$nthread

    set.seed(opt$seed)
    err_msg <- NULL
    mod_fit <- NULL
    invisible(utils::capture.output({
      mod_fit <- tryCatch({
        suppressWarnings(suppressMessages(do.call(caret::train, args_train)))
      }, error = function(e) {
        err_msg <<- e$message
        NULL
      })
    }, type = "output"))
    if (is.null(mod_fit) && !is.null(err_msg)) {
      cat(sprintf("[OMITIDO: %s]\n", err_msg))
    }

    if (!is.null(mod_fit)) {
      mods[[mod_name]] <- mod_fit
      pred_cv <- mod_fit$pred
      if (!is.null(pred_cv) && nrow(pred_cv) > 0) {
        score_tp <- if ("TP" %in% colnames(pred_cv)) pred_cv$TP else NULL
        rows[[mod_name]] <- compute_classification_metrics(
          pred = pred_cv$pred,
          obs = pred_cv$obs,
          score = score_tp,
          positive_class = "TP",
          model_name = mod_name
        )
        cms[[mod_name]] <- table(
          Prediccion = factor(pred_cv$pred, levels = c("TP", "FP")),
          Referencia = factor(pred_cv$obs,  levels = c("TP", "FP"))
        )
      }
      cat("OK!\n")
    }
  }

  # Meta-Modelo GLM Stacking (caretEnsemble) para el escenario actual
  if (opt$train_stacking && length(mods) >= 2 && requireNamespace("caretEnsemble", quietly = TRUE)) {
    cat(sprintf("   [%s] Entrenando: %-18s ... ", scenario_label, "Ensamble_Stacking"))
    tryCatch({
      bench_tmp <- dplyr::bind_rows(rows) %>% dplyr::arrange(dplyr::desc(ROC_AUC))
      top2_names <- utils::head(bench_tmp$Modelo, 2)
      lista_top <- mods[top2_names]
      class(lista_top) <- c("caretList", "list")
      set.seed(opt$seed)
      ens_stack <- suppressWarnings(suppressMessages(
        caretEnsemble::caretStack(
          lista_top,
          method = "glm",
          metric = "ROC",
          trControl = caret::trainControl(
            method = "cv",
            number = opt$cv_folds,
            index = control_obj$index,
            savePredictions = "final",
            classProbs = TRUE,
            summaryFunction = caret::twoClassSummary
          )
        )
      ))
      mods[["Ensamble_Stacking"]] <- ens_stack

      prob_base_list <- lapply(lista_top, function(m) stats::predict(m, newdata = x_train, type = "prob")[["TP"]])
      df_prob_base <- as.data.frame(prob_base_list)
      colnames(df_prob_base) <- top2_names
      probs_stack <- stats::predict(ens_stack$ens_model, newdata = df_prob_base, type = "prob")[["TP"]]
      preds_stack <- factor(ifelse(probs_stack >= 0.5, "TP", "FP"), levels = c("TP", "FP"))

      rows[["Ensamble_Stacking"]] <- compute_classification_metrics(
        pred = preds_stack,
        obs = y_train,
        score = probs_stack,
        positive_class = "TP",
        model_name = "Ensamble_Stacking"
      )
      cms[["Ensamble_Stacking"]] <- table(Prediccion = preds_stack, Referencia = y_train)
      cat(sprintf("OK! (%s)\n", paste(top2_names, collapse = " + ")))
    }, error = function(e) {
      cat(sprintf("[OMITIDO: %s]\n", e$message))
    })
  }

  return(list(models = mods, benchmark = dplyr::bind_rows(rows), cms = cms))
}

# 5A. Entrenar Escenario Desbalanceado (sin remuestreo, preservando exactamente la semilla 1234)
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

# 5B. Entrenar Escenario Balanceado (respeta --sampling down, up o none y preserva la semilla 1234)
cat("[5/7] Entrenando torneo de modelos (Escenario Balanceado)...\n")
if (opt$sampling == "none") {
  warning("[02_train_models.R] Se especifico --sampling none: el escenario balanceado se entrenara sin remuestreo.", call. = FALSE)
}
fit_control_bal <- caret::trainControl(
  method = "cv",
  number = opt$cv_folds,
  classProbs = TRUE,
  summaryFunction = caret::twoClassSummary,
  allowParallel = (opt$cores > 1),
  savePredictions = "final",
  sampling = if (opt$sampling != "none") opt$sampling else NULL
)
res_bal <- train_tournament(fit_control_bal, "Balanceado")

# Liberar los workers paralelos una vez terminado el entrenamiento
if (opt$cores > 1) {
  try(parallel::stopCluster(cl), silent = TRUE)
  foreach::registerDoSEQ()
}

# Consolidar tablas de benchmark con TODOS los modelos evaluados (Heuristico + Todos los ML + 4 Cortes Youden)
df_benchmark_bal   <- dplyr::bind_rows(heuristic_res$benchmark_row, res_bal$benchmark, youden_eval$benchmark)
df_benchmark_unbal <- dplyr::bind_rows(heuristic_res$benchmark_row, res_unbal$benchmark, youden_eval$benchmark)
df_benchmark_all   <- dplyr::bind_rows(
  dplyr::mutate(df_benchmark_bal,   Escenario = "Balanceado",   .before = 1),
  dplyr::mutate(df_benchmark_unbal, Escenario = "Desbalanceado", .before = 1)
)

# Si se indico --clean_models, limpiar artefactos previos en out_models; de lo contrario, conservar modelos previos no reentrenados en esta corrida
if (isTRUE(opt$clean_models)) {
  cat(sprintf("   [--clean_models] Limpiando archivos .rds y .tsv previos en %s...\n", opt$out_models))
  old_files <- list.files(opt$out_models, pattern = "\\.(rds|tsv)$", full.names = TRUE)
  if (length(old_files) > 0) unlink(old_files)
} else {
  prev_bal_file <- file.path(opt$out_models, "models_balanced.rds")
  if (file.exists(prev_bal_file)) {
    prev_bal <- normalize_model_names(readRDS(prev_bal_file))
    for (pm in setdiff(names(prev_bal), names(res_bal$models))) {
      res_bal$models[[pm]] <- prev_bal[[pm]]
    }
  }
  prev_unbal_file <- file.path(opt$out_models, "models_unbalanced.rds")
  if (file.exists(prev_unbal_file)) {
    prev_unbal <- normalize_model_names(readRDS(prev_unbal_file))
    for (pm in setdiff(names(prev_unbal), names(res_unbal$models))) {
      res_unbal$models[[pm]] <- prev_unbal[[pm]]
    }
  }
}

# Serializar contenedores y modelos individuales (actualizando los modelos entrenados en esta corrida)
saveRDS(res_bal$models,   file = file.path(opt$out_models, "models_balanced.rds"))
saveRDS(res_bal$models,   file = file.path(opt$out_models, "all_trained_models.rds"))
saveRDS(res_unbal$models, file = file.path(opt$out_models, "models_unbalanced.rds"))

for (mod_name in names(res_bal$models)) {
  safe_fn <- gsub("[^a-z0-9_]", "_", tolower(mod_name))
  Safe_fn <- gsub("_+", "_", gsub("^_|_$", "", safe_fn))
  saveRDS(res_bal$models[[mod_name]], file = file.path(opt$out_models, sprintf("model_%s.rds", Safe_fn)))
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

# Guardar metricas Youden detalladas
df_youden_bench <- youden_eval$benchmark
df_youden_bench$Umbral <- sapply(metricas_base, function(v) round(youden_eval$thresholds[[v]], 3))
df_youden_bench$Indice_Youden_J <- df_youden_bench$Sensibilidad + df_youden_bench$Especificidad - 1
saveRDS(df_youden_bench, file = file.path(opt$out_models, "youden_benchmark.rds"))
utils::write.table(df_youden_bench, file = file.path(opt$out_models, "youden_benchmark.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Consolidar matrices de confusion de todos los modelos evaluados
cm_youden <- list()
for (var in metricas_base) {
  pred_y <- factor(ifelse(df_ml[[var]] >= youden_eval$thresholds[[var]], "TP", "FP"), levels = c("TP", "FP"))
  cm_youden[[var]] <- table(Prediccion = pred_y, Referencia = df_ml$Clase)
}

pred_h <- predict_weighted_heuristic(df_ml, youden_eval$thresholds, heuristic_res$weights, heuristic_res$global_threshold)
cm_heur <- table(Prediccion = pred_h, Referencia = df_ml$Clase)

confusion_matrices <- list(
  youden     = cm_youden,
  heuristic  = cm_heur,
  balanced   = res_bal$cms,
  unbalanced = res_unbal$cms,
  ml         = res_bal$cms
)
saveRDS(confusion_matrices, file = file.path(opt$out_models, "confusion_matrices.rds"))

# Guardar manifiesto de consenso y lista de predictores para 03_cnv_inference.R
pref_ml <- intersect(c("ElasticNet", "NaiveBayes", "SVM_Polinomial", "Red_Neuronal"), names(res_bal$models))
if (length(pref_ml) < 4 && nrow(res_bal$benchmark) > 0) {
  ranked_ml <- res_bal$benchmark %>%
    dplyr::filter(Modelo != "Ensamble_Stacking") %>%
    dplyr::arrange(dplyr::desc(ROC_AUC), dplyr::desc(F1_Score)) %>%
    dplyr::pull(Modelo)
  pref_ml <- unique(c(pref_ml, ranked_ml))[seq_len(min(4L, length(unique(c(pref_ml, ranked_ml)))))]
}
consensus_spec <- list(
  predictors = predictores,
  metricas_base = metricas_base,
  ml_voters = pref_ml,
  youden_voters = c("pct_recovered", "pct_identity"),
  n_voters = length(pref_ml) + 2L + 1L,
  seed = opt$seed,
  sampling = opt$sampling
)
saveRDS(consensus_spec, file = file.path(opt$out_models, "consensus_spec.rds"))

# 6. Explicabilidad SHAP (KernelSHAP para Escenario Balanceado y Desbalanceado)
cat("[6/7] Generando explicaciones SHAP (Interpretabilidad con KernelSHAP)...\n")
if ("ElasticNet" %in% names(res_bal$models) || "ElasticNet" %in% names(res_unbal$models)) {
  suppressPackageStartupMessages({
    library(kernelshap)
    library(shapviz)
  })
  set.seed(opt$seed)
  bg_idx <- sample(seq_len(nrow(x_train)), min(30, nrow(x_train)))

  if ("ElasticNet" %in% names(res_bal$models)) {
    shap_calc_bal <- tryCatch({
      kernelshap::kernelshap(
        object = res_bal$models[["ElasticNet"]],
        X = x_train,
        bg_X = x_train[bg_idx, , drop = FALSE],
        pred_fun = function(m, d) predict(m, newdata = d, type = "prob")[, "TP"]
      )
    }, error = function(e) NULL)
    if (!is.null(shap_calc_bal)) {
      shap_obj_bal <- shapviz::shapviz(shap_calc_bal)
      saveRDS(shap_obj_bal, file = file.path(opt$out_models, "shap_values.rds"))
      saveRDS(shap_obj_bal, file = file.path(opt$out_models, "shap_values_balanced.rds"))
      cat("   -> Artefacto SHAP (Balanceado) serializado exitosamente.\n")
    }
  }

  if ("ElasticNet" %in% names(res_unbal$models)) {
    shap_calc_unbal <- tryCatch({
      kernelshap::kernelshap(
        object = res_unbal$models[["ElasticNet"]],
        X = x_train,
        bg_X = x_train[bg_idx, , drop = FALSE],
        pred_fun = function(m, d) predict(m, newdata = d, type = "prob")[, "TP"]
      )
    }, error = function(e) NULL)
    if (!is.null(shap_calc_unbal)) {
      shap_obj_unbal <- shapviz::shapviz(shap_calc_unbal)
      saveRDS(shap_obj_unbal, file = file.path(opt$out_models, "shap_values_unbalanced.rds"))
      cat("   -> Artefacto SHAP (Desbalanceado) serializado exitosamente.\n")
    }
  }
}

# 7. Renderizado del Reporte HTML
cat("[7/7] Renderizando reporte HTML interactivo con Quarto...\n")
qmd_file <- file.path(project_root, "training_report.qmd")
if (!opt$render_report) {
  cat("   -> Renderizado omitido (--no_render_report).\n")
} else if (!file.exists(qmd_file)) {
  cat(sprintf("   -> [AVISO] Renderizado omitido: no existe la plantilla %s\n", qmd_file))
} else if (!nzchar(Sys.which("quarto"))) {
  cat("   -> [AVISO] Renderizado omitido: 'quarto' no esta instalado o no esta en el PATH.\n")
} else {
  out_html <- normalizePath(opt$report_file, mustWork = FALSE)
  dir.create(dirname(out_html), recursive = TRUE, showWarnings = FALSE)
  rendered_in_root <- file.path(dirname(qmd_file), basename(out_html))

  src_arg <- if (!is.null(opt$sources) && nzchar(trimws(opt$sources))) sprintf(" -P sources_file:%s", shQuote(normalizePath(opt$sources, mustWork = FALSE))) else ""
  cmd_quarto <- sprintf(
    "quarto render %s --output %s -P project_root:%s -P data_path:%s -P data_file:%s -P models_dir:%s%s",
    shQuote(qmd_file),
    shQuote(basename(out_html)),
    shQuote(project_root),
    shQuote(normalizePath(opt$data, mustWork = FALSE)),
    shQuote(normalizePath(opt$data, mustWork = FALSE)),
    shQuote(normalizePath(opt$out_models, mustWork = FALSE)),
    src_arg
  )
  cat(sprintf("   -> Ejecutando Quarto: %s\n", cmd_quarto))
  res_q <- system(cmd_quarto)
  if (res_q == 0) {
    if (file.exists(rendered_in_root) && normalizePath(rendered_in_root, mustWork = FALSE) != out_html) {
      file.copy(rendered_in_root, out_html, overwrite = TRUE)
      unlink(rendered_in_root)
    }
    cat(sprintf("   -> [OK] Reporte HTML generado: %s\n", out_html))
  } else {
    cat(sprintf("   -> [AVISO] Fallo el renderizado con Quarto (codigo de retorno: %d).\n", res_q))
  }
}

elapsed_sec <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

cat("--------------------------------------------------------\n")
options(width = 200)
cat(sprintf(" BENCHMARK COMPARATIVO FINAL — BALANCEADO (%d MODELOS):\n", nrow(df_benchmark_bal)))
cat("--------------------------------------------------------\n")
print(df_benchmark_bal %>% dplyr::arrange(dplyr::desc(ROC_AUC), dplyr::desc(F1_Score)), digits = 3)
cat("--------------------------------------------------------\n")
cat(sprintf(" BENCHMARK COMPARATIVO FINAL — DESBALANCEADO (%d MODELOS):\n", nrow(df_benchmark_unbal)))
cat("--------------------------------------------------------\n")
print(df_benchmark_unbal %>% dplyr::arrange(dplyr::desc(ROC_AUC), dplyr::desc(F1_Score)), digits = 3)
cat("========================================================\n")
cat(sprintf("[OK] Todos los modelos serializados en: %s (Tiempo total: %.2f s)\n", opt$out_models, elapsed_sec))
