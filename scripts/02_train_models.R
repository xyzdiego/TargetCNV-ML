#!/usr/bin/env Rscript

# Módulo 02: Diagnóstico de colinealidad, entrenamiento supervisado con validación
# cruzada determinista, optimización de umbrales de Youden y explicabilidad KernelSHAP

t_start <- Sys.time()

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
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
source(file.path(project_root, "R", "model_evaluation.R"))
source(file.path(project_root, "R", "data_processing.R"))
source(file.path(project_root, "R", "plotting.R"))
source(file.path(project_root, "R", "report_helpers.R"))

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
    default = 5L,
    help = "Número de particiones (folds) para validación cruzada k-fold [Por defecto: 5]",
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
    default = 1234L,
    help = "Semilla maestra para particiones y semillas explícitas de caret [Por defecto: 1234]",
    metavar = "INT"
  ),
  make_option(
    c("-c", "--cores"),
    type = "integer",
    default = max(1L, parallel::detectCores() - 1L),
    help = "Número de hilos para entrenamiento paralelo [Por defecto: cores - 1]",
    metavar = "INT"
  ),
  make_option(
    c("--include_hit"),
    action = "store_true",
    default = FALSE,
    help = "Incluir el índice ordinal 'hit' de CAPTUS como predictor continuo adicional"
  ),
  make_option(
    c("--extended_tournament"),
    action = "store_true",
    default = TRUE,
    help = "Ejecutar la comparación completa de 18 clasificadores supervisados [Por defecto: TRUE]"
  ),
  make_option(
    c("--fast_tournament"),
    action = "store_false",
    dest = "extended_tournament",
    help = "Entrenar únicamente los 5 clasificadores base rápidos (ElasticNet, NaiveBayes, SVM_Polinomial, Red_Neuronal, LDA)"
  ),
  make_option(
    c("--train_stacking"),
    action = "store_true",
    default = TRUE,
    help = "Entrenar meta-modelo GLM Stacking (caretEnsemble) sobre los dos mejores clasificadores base [Por defecto: TRUE]"
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
    help = "Eliminar archivos .rds/.tsv previos en --out_models antes de escribir la corrida actual [Por defecto: FALSE]"
  ),
  make_option(
    c("--render_report"),
    action = "store_true",
    default = TRUE,
    help = "Compilar el reporte HTML de entrenamiento (training_report.qmd) en modo de solo lectura [Por defecto: TRUE]"
  ),
  make_option(
    c("--no_render_report"),
    action = "store_false",
    dest = "render_report",
    help = "Omitir la compilación del reporte HTML"
  ),
  make_option(
    c("-r", "--report_file"),
    type = "character",
    default = file.path(project_root, "reports", "reporte_entrenamiento.html"),
    help = "Ruta del reporte HTML de salida [Por defecto: reports/reporte_entrenamiento.html]",
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
  description = "Evalúa colinealidad de Spearman, optimiza cortes univariados de Youden y el ensamble heurístico, entrena los clasificadores supervisados con semillas deterministas en validación cruzada y compila el reporte HTML."
)

opt <- parse_args(opt_parser)

if (is.null(opt$data) || !file.exists(opt$data)) {
  tsv_fallback <- sub("\\.rds$", ".tsv", opt$data)
  if (file.exists(tsv_fallback)) {
    opt$data <- tsv_fallback
  } else {
    print_help(opt_parser)
    stop(sprintf("[ERROR] No se encontró el archivo de datos preparados: %s", opt$data), call. = FALSE)
  }
}

if (!dir.exists(opt$out_models)) {
  dir.create(opt$out_models, recursive = TRUE, showWarnings = FALSE)
}

if (!opt$sampling %in% c("down", "up", "none")) {
  stop(sprintf("[ERROR] --sampling debe ser 'down', 'up' o 'none' (recibido: '%s').", opt$sampling), call. = FALSE)
}
if (is.na(opt$cv_folds) || opt$cv_folds < 2L) {
  stop("[ERROR] --cv_folds debe ser un entero >= 2.", call. = FALSE)
}
if (is.na(opt$cores) || opt$cores < 1L) {
  stop("[ERROR] --cores debe ser un entero >= 1.", call. = FALSE)
}

cat("========================================================\n")
cat(" TargetCNV-ML: Módulo 02 - Entrenamiento y Validación Cruzada\n")
cat("========================================================\n")
cat(sprintf("Archivo de datos          : %s\n", opt$data))
cat(sprintf("Directorio modelos (.rds) : %s\n", opt$out_models))
cat(sprintf("Folds CV                  : %d\n", opt$cv_folds))
cat(sprintf("Balanceo de Clases        : %s\n", opt$sampling))
cat(sprintf("Esquema de Modelos        : %s\n", if (opt$extended_tournament) "Completo (18 clasificadores ML)" else "Rápido (5 clasificadores ML)"))
cat(sprintf("Meta-Modelo Stacking      : %s\n", if (opt$train_stacking) "SI" else "NO"))
cat(sprintf("Núcleos Paralelos         : %d\n", opt$cores))
cat(sprintf("Semilla Determinista      : %d\n", opt$seed))
cat(sprintf("Compilar Reporte HTML     : %s\n", if (opt$render_report) "SI" else "NO"))
if (opt$render_report) cat(sprintf("Reporte HTML de salida    : %s\n", opt$report_file))
cat("--------------------------------------------------------\n")

# 1. Carga y verificación de datos preparados
cat("[1/8] Cargando datos preparados...\n")
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
    "[ERROR] El archivo de datos no contiene las columnas requeridas: %s.\nGenere el archivo con scripts/01_prepare_data.R",
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

cat(sprintf("   -> Secuencias válidas: %d (TP = %d, FP = %d)\n",
            nrow(df_ml), sum(df_ml$Clase == "TP"), sum(df_ml$Clase == "FP")))

if (min(table(df_ml$Clase)) < opt$cv_folds) {
  stop(sprintf(
    "[ERROR] La clase minoritaria tiene %d muestras (< %d folds). Se requieren suficientes registros TP y FP.\n",
    min(table(df_ml$Clase)), opt$cv_folds
  ), call. = FALSE)
}

# 2. Test de colinealidad por rangos de Spearman entre predictores continuos
cat("[2/8] Evaluando colinealidad de Spearman entre predictores...\n")
spearman_res <- compute_spearman_collinearity(df_ml, metric_cols = predictores)
print(
  spearman_res$pairs %>%
    dplyr::mutate(
      Spearman_Rho = round(Spearman_Rho, 4),
      P_Valor = signif(P_Valor, 3)
    ),
  row.names = FALSE
)

# 3. Evaluación de cortes univariados óptimos de Youden
cat("[3/8] Calculando umbrales univariados óptimos de Youden...\n")
youden_eval <- evaluate_youden_predictors(
  df = df_ml,
  response_col = "Clase",
  metric_cols = metricas_base,
  positive_class = "TP"
)
cat("   -> Umbrales optimizados:\n")
print(round(youden_eval$thresholds, 2))

# 4. Optimización del ensamble heurístico ponderado
cat("[4/8] Optimizando ensamble heurístico de votación ponderada...\n")
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
cat("   -> Pesos heurísticos optimizados:\n")
print(heuristic_res$weights)
cat(sprintf("   -> Umbral global: %.2f | AUC: %.4f\n", heuristic_res$global_threshold, heuristic_res$best_auc))

# Configuración de clúster paralelo PSOCK (independiente de la semilla determinista por fold)
cl <- NULL
if (opt$cores > 1L) {
  suppressPackageStartupMessages(library(doParallel))
  cl <- parallel::makeCluster(opt$cores)
  doParallel::registerDoParallel(cl)
}

model_specs <- get_model_specs(n_obs = nrow(df_ml), extended = opt$extended_tournament)
x_train <- as.data.frame(df_ml[, predictores, drop = FALSE])
y_train <- df_ml$Clase

train_tournament <- function(control_obj, scenario_label) {
  mods <- list()
  rows <- list()
  cms  <- list()
  shared_index    <- NULL
  shared_indexOut <- NULL

  for (mod_name in names(model_specs)) {
    spec <- model_specs[[mod_name]]
    cat(sprintf("   [%s] Entrenando: %-25s ... ", scenario_label, mod_name))

    pre_proc <- if (is.null(spec$no_scale) || !spec$no_scale) c("center", "scale") else NULL
    det_ctrl <- build_deterministic_cv_control(
      base_ctrl = control_obj,
      x = x_train,
      y = y_train,
      method = spec$method,
      tuneGrid = spec$tuneGrid,
      tuneLength = spec$tuneLength,
      preProcess = pre_proc,
      seed = opt$seed
    )
    if (is.null(shared_index)) {
      shared_index    <- det_ctrl$trControl$index
      shared_indexOut <- det_ctrl$trControl$indexOut
    }

    args_train <- list(
      x = x_train,
      y = y_train,
      method = spec$method,
      trControl = det_ctrl$trControl,
      tuneGrid = det_ctrl$tuneGrid,
      metric = "ROC"
    )
    if (!is.null(pre_proc))         args_train$preProcess <- pre_proc
    if (!is.null(spec$trace))       args_train$trace      <- spec$trace
    if (!is.null(spec$MaxNWts))     args_train$MaxNWts    <- spec$MaxNWts
    if (!is.null(spec$verbose))     args_train$verbose    <- spec$verbose
    if (!is.null(spec$nthread))     args_train$nthread    <- spec$nthread

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
      cat("OK\n")
    }
  }

  # Meta-modelo GLM Stacking con predicciones fuera de muestra alineadas y semillas explícitas
  if (opt$train_stacking && length(mods) >= 2L && requireNamespace("caretEnsemble", quietly = TRUE)) {
    cat(sprintf("   [%s] Entrenando: %-25s ... ", scenario_label, "Ensamble_Stacking"))
    tryCatch({
      bench_tmp  <- dplyr::bind_rows(rows) %>% dplyr::arrange(dplyr::desc(ROC_AUC), dplyr::desc(F1_Score))
      top2_names <- utils::head(bench_tmp$Modelo, 2L)
      lista_top  <- mods[top2_names]
      class(lista_top) <- c("caretList", "list")

      stack_base_ctrl <- caret::trainControl(
        method = "cv",
        number = opt$cv_folds,
        index = shared_index,
        indexOut = shared_indexOut,
        savePredictions = "final",
        classProbs = TRUE,
        allowParallel = (opt$cores > 1L),
        summaryFunction = caret::twoClassSummary
      )
      dummy_meta_x <- as.data.frame(matrix(0.5, nrow = nrow(x_train), ncol = length(top2_names)))
      colnames(dummy_meta_x) <- top2_names
      det_stack <- build_deterministic_cv_control(
        base_ctrl = stack_base_ctrl,
        x = dummy_meta_x,
        y = y_train,
        method = "glm",
        tuneLength = 1L,
        seed = opt$seed
      )

      set.seed(opt$seed)
      ens_stack <- suppressWarnings(suppressMessages(
        caretEnsemble::caretStack(
          lista_top,
          method = "glm",
          metric = "ROC",
          trControl = det_stack$trControl
        )
      ))
      mods[["Ensamble_Stacking"]] <- ens_stack

      probs_stack <- predict_stacking_probs(ens_stack, newdata = x_train)
      preds_stack <- factor(ifelse(probs_stack >= 0.5, "TP", "FP"), levels = c("TP", "FP"))
      rows[["Ensamble_Stacking"]] <- compute_classification_metrics(
        pred = preds_stack,
        obs = y_train,
        score = probs_stack,
        positive_class = "TP",
        model_name = "Ensamble_Stacking"
      )
      cms[["Ensamble_Stacking"]] <- table(Prediccion = preds_stack, Referencia = y_train)
      cat(sprintf("OK (%s)\n", paste(top2_names, collapse = " + ")))
    }, error = function(e) {
      cat(sprintf("[OMITIDO: %s]\n", e$message))
    })
  }

  list(models = mods, benchmark = dplyr::bind_rows(rows), cms = cms)
}

# 5A. Escenario Desbalanceado (distribución natural de clases)
cat("[5/8] Entrenando clasificadores supervisados (Escenario Desbalanceado)...\n")
fit_control_unbal <- caret::trainControl(
  method = "cv",
  number = opt$cv_folds,
  classProbs = TRUE,
  summaryFunction = caret::twoClassSummary,
  allowParallel = (opt$cores > 1L),
  savePredictions = "final",
  sampling = NULL
)
res_unbal <- train_tournament(fit_control_unbal, "Desbalanceado")

# 5B. Escenario Balanceado (submuestreo/sobremuestreo dentro de cada fold de CV)
cat("[6/8] Entrenando clasificadores supervisados (Escenario Balanceado)...\n")
if (opt$sampling == "none") {
  warning("[02_train_models.R] Se especificó --sampling none: el escenario balanceado se entrenará sin remuestreo.", call. = FALSE)
}
fit_control_bal <- caret::trainControl(
  method = "cv",
  number = opt$cv_folds,
  classProbs = TRUE,
  summaryFunction = caret::twoClassSummary,
  allowParallel = (opt$cores > 1L),
  savePredictions = "final",
  sampling = if (opt$sampling != "none") opt$sampling else NULL
)
res_bal <- train_tournament(fit_control_bal, "Balanceado")

if (opt$cores > 1L && !is.null(cl)) {
  try(parallel::stopCluster(cl), silent = TRUE)
  foreach::registerDoSEQ()
}

# Consolidar tablas de benchmark completas
df_benchmark_bal   <- dplyr::bind_rows(heuristic_res$benchmark_row, res_bal$benchmark, youden_eval$benchmark)
df_benchmark_unbal <- dplyr::bind_rows(heuristic_res$benchmark_row, res_unbal$benchmark, youden_eval$benchmark)
df_benchmark_all   <- dplyr::bind_rows(
  dplyr::mutate(df_benchmark_bal,   Escenario = "Balanceado",    .before = 1),
  dplyr::mutate(df_benchmark_unbal, Escenario = "Desbalanceado", .before = 1)
)

# Gestión de artefactos previos antes de serializar los nuevos resultados
if (isTRUE(opt$clean_models)) {
  cat(sprintf("   [--clean_models] Eliminando archivos .rds y .tsv previos en %s...\n", opt$out_models))
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

# Serializar colinealidad de Spearman, umbrales de Youden y modelo heurístico
saveRDS(spearman_res, file = file.path(opt$out_models, "spearman_collinearity.rds"))
utils::write.table(
  spearman_res$pairs,
  file = file.path(opt$out_models, "spearman_collinearity.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

saveRDS(youden_eval$thresholds, file = file.path(opt$out_models, "youden_thresholds.rds"))
saveRDS(heuristic_obj,          file = file.path(opt$out_models, "heuristic_model.rds"))

# Serializar contenedores de modelos y archivos individuales
saveRDS(res_bal$models,   file = file.path(opt$out_models, "models_balanced.rds"))
saveRDS(res_bal$models,   file = file.path(opt$out_models, "all_trained_models.rds"))
saveRDS(res_unbal$models, file = file.path(opt$out_models, "models_unbalanced.rds"))

for (mod_name in names(res_bal$models)) {
  safe_fn <- gsub("[^a-z0-9_]", "_", tolower(mod_name))
  safe_fn <- gsub("_+", "_", gsub("^_|_$", "", safe_fn))
  saveRDS(res_bal$models[[mod_name]], file = file.path(opt$out_models, sprintf("model_%s.rds", safe_fn)))
}

# Serializar tablas de benchmark
saveRDS(df_benchmark_bal, file = file.path(opt$out_models, "benchmark_models_balanced.rds"))
saveRDS(df_benchmark_bal, file = file.path(opt$out_models, "benchmark_models.rds"))
utils::write.table(df_benchmark_bal, file = file.path(opt$out_models, "benchmark_models.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(df_benchmark_bal, file = file.path(opt$out_models, "benchmark_models_balanced.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

saveRDS(df_benchmark_unbal, file = file.path(opt$out_models, "benchmark_models_unbalanced.rds"))
utils::write.table(df_benchmark_unbal, file = file.path(opt$out_models, "benchmark_models_unbalanced.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(df_benchmark_all, file = file.path(opt$out_models, "benchmark_models_all.rds"))
utils::write.table(df_benchmark_all, file = file.path(opt$out_models, "benchmark_models_all.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Serializar métricas detalladas de Youden
df_youden_bench <- youden_eval$benchmark
df_youden_bench$Umbral <- sapply(metricas_base, function(v) round(youden_eval$thresholds[[v]], 3))
df_youden_bench$Indice_Youden_J <- df_youden_bench$Sensibilidad + df_youden_bench$Especificidad - 1
saveRDS(df_youden_bench, file = file.path(opt$out_models, "youden_benchmark.rds"))
utils::write.table(df_youden_bench, file = file.path(opt$out_models, "youden_benchmark.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Consolidar matrices de confusión
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

# Guardar especificación del sistema de consenso de 7 modelos
pref_ml <- intersect(c("ElasticNet", "NaiveBayes", "SVM_Polinomial", "Red_Neuronal"), names(res_bal$models))
if (length(pref_ml) < 4L && nrow(res_bal$benchmark) > 0) {
  ranked_ml <- res_bal$benchmark %>%
    dplyr::filter(Modelo != "Ensamble_Stacking") %>%
    dplyr::arrange(dplyr::desc(ROC_AUC), dplyr::desc(F1_Score)) %>%
    dplyr::pull(Modelo)
  pref_ml <- unique(c(pref_ml, ranked_ml))[seq_len(min(4L, length(unique(c(pref_ml, ranked_ml)))))]
}
consensus_spec <- list(
  predictors    = predictores,
  metricas_base = metricas_base,
  ml_voters     = pref_ml,
  youden_voters = c("pct_recovered", "pct_identity"),
  n_voters      = length(pref_ml) + 2L + 1L,
  seed          = opt$seed,
  sampling      = opt$sampling
)
saveRDS(consensus_spec, file = file.path(opt$out_models, "consensus_spec.rds"))

# 7. Explicabilidad con KernelSHAP
cat("[7/8] Calculando descomposiciones aditivas SHAP (KernelSHAP)...\n")
if ("ElasticNet" %in% names(res_bal$models) || "ElasticNet" %in% names(res_unbal$models)) {
  suppressPackageStartupMessages({
    library(kernelshap)
    library(shapviz)
  })
  set.seed(opt$seed)
  bg_idx <- sample(seq_len(nrow(x_train)), min(30L, nrow(x_train)))

  if ("ElasticNet" %in% names(res_bal$models)) {
    shap_calc_bal <- tryCatch({
      kernelshap::kernelshap(
        object = res_bal$models[["ElasticNet"]],
        X = x_train,
        bg_X = x_train[bg_idx, , drop = FALSE],
        pred_fun = function(m, d) stats::predict(m, newdata = d, type = "prob")[, "TP"]
      )
    }, error = function(e) NULL)
    if (!is.null(shap_calc_bal)) {
      shap_obj_bal <- shapviz::shapviz(shap_calc_bal)
      saveRDS(shap_obj_bal, file = file.path(opt$out_models, "shap_values.rds"))
      saveRDS(shap_obj_bal, file = file.path(opt$out_models, "shap_values_balanced.rds"))
      cat("   -> Artefacto SHAP (Balanceado) serializado.\n")
    }
  }

  if ("ElasticNet" %in% names(res_unbal$models)) {
    shap_calc_unbal <- tryCatch({
      kernelshap::kernelshap(
        object = res_unbal$models[["ElasticNet"]],
        X = x_train,
        bg_X = x_train[bg_idx, , drop = FALSE],
        pred_fun = function(m, d) stats::predict(m, newdata = d, type = "prob")[, "TP"]
      )
    }, error = function(e) NULL)
    if (!is.null(shap_calc_unbal)) {
      shap_obj_unbal <- shapviz::shapviz(shap_calc_unbal)
      saveRDS(shap_obj_unbal, file = file.path(opt$out_models, "shap_values_unbalanced.rds"))
      cat("   -> Artefacto SHAP (Desbalanceado) serializado.\n")
    }
  }
}

# 8. Compilación de reporte Quarto de solo lectura
cat("[8/8] Compilando reporte HTML de entrenamiento con Quarto...\n")
qmd_file <- file.path(project_root, "training_report.qmd")
if (!opt$render_report) {
  cat("   -> Compilación omitida (--no_render_report).\n")
} else if (!file.exists(qmd_file)) {
  cat(sprintf("   -> [AVISO] Compilación omitida: no existe la plantilla %s\n", qmd_file))
} else if (!nzchar(Sys.which("quarto"))) {
  cat("   -> [AVISO] Compilación omitida: 'quarto' no está disponible en el PATH.\n")
} else {
  out_html <- normalizePath(opt$report_file, mustWork = FALSE)
  dir.create(dirname(out_html), recursive = TRUE, showWarnings = FALSE)
  rendered_in_root <- file.path(dirname(qmd_file), basename(out_html))

  src_arg <- if (!is.null(opt$sources) && nzchar(trimws(opt$sources))) {
    sprintf(" -P sources_file:%s", shQuote(normalizePath(opt$sources, mustWork = FALSE)))
  } else {
    ""
  }
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
    cat(sprintf("   -> [AVISO] Falló el renderizado con Quarto (código de retorno: %d).\n", res_q))
  }
}

elapsed_sec <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

cat("--------------------------------------------------------\n")
options(width = 200)
cat(sprintf(" BENCHMARK COMPARATIVO FINAL - ESCENARIO BALANCEADO (%d METODOS):\n", nrow(df_benchmark_bal)))
cat("--------------------------------------------------------\n")
print(df_benchmark_bal %>% dplyr::arrange(dplyr::desc(ROC_AUC), dplyr::desc(F1_Score)), digits = 3)
cat("--------------------------------------------------------\n")
cat(sprintf(" BENCHMARK COMPARATIVO FINAL - ESCENARIO DESBALANCEADO (%d METODOS):\n", nrow(df_benchmark_unbal)))
cat("--------------------------------------------------------\n")
print(df_benchmark_unbal %>% dplyr::arrange(dplyr::desc(ROC_AUC), dplyr::desc(F1_Score)), digits = 3)
cat("========================================================\n")
cat(sprintf("[OK] Artefactos serializados en: %s (Tiempo total: %.2f s)\n", opt$out_models, elapsed_sec))
