#' Model Evaluation and Consensus Inference Module for TargetCNV-ML
#'
#' Pure functions for confusion matrix metrics extraction, Youden index optimization,
#' heuristic weighted voting, multi-model consensus consolidation, and deployment inference.
#'
#' @name model_evaluation
NULL

#' Evaluar Matriz de Confusion y Metricas de Rendimiento
#'
#' Extrae y estandariza las metricas clave de evaluacion diagnostica (Exactitud,
#' Sensibilidad, Especificidad, Precision, FDR, F1-Score) a partir de una matriz
#' de confusion o de vectores de prediccion y referencia.
#'
#' @param data Vector de predicciones (factor) o tabla de confusion 2x2.
#' @param reference Vector de referencia (factor con clases Ground Truth). Obligatorio si data es vector.
#' @param positive_class Etiqueta de la clase positiva (por defecto "TP").
#' @param model_name Nombre opcional del modelo para registrar en el data frame de salida.
#'
#' @return Un \code{data.frame} de una fila con las columnas: \code{Modelo}, \code{Exactitud},
#'   \code{Sensibilidad}, \code{Especificidad}, \code{Precision}, \code{FDR}, \code{F1_Score}.
#' @export
#' @examples
#' \dontrun{
#' res <- evaluate_confusion_matrix(pred, real, model_name = "ElasticNet")
#' }
evaluate_confusion_matrix <- function(data,
                                      reference = NULL,
                                      positive_class = "TP",
                                      model_name = "Modelo") {
  if (inherits(data, "confusionMatrix")) {
    cm <- data
  } else if (is.table(data) || is.matrix(data)) {
    cm <- list(
      table = data,
      overall = c(Accuracy = sum(diag(data)) / sum(data)),
      byClass = c(
        Sensitivity = data[positive_class, positive_class] / sum(data[, positive_class]),
        Specificity = sum(data[-which(rownames(data) == positive_class), -which(colnames(data) == positive_class)]) /
          sum(data[, -which(colnames(data) == positive_class)]),
        Precision = data[positive_class, positive_class] / sum(data[positive_class, ])
      )
    )
    cm$byClass["F1"] <- 2 * (cm$byClass["Precision"] * cm$byClass["Sensitivity"]) /
      (cm$byClass["Precision"] + cm$byClass["Sensitivity"])
  } else {
    stopifnot(!is.null(reference))
    data_fac <- factor(as.character(data), levels = c("TP", "FP"))
    ref_fac  <- factor(as.character(reference), levels = c("TP", "FP"))
    cm <- caret::confusionMatrix(data = data_fac, reference = ref_fac, mode = "everything")
  }

  acc  <- as.numeric(cm$overall["Accuracy"])
  sens <- as.numeric(cm$byClass["Sensitivity"])
  spec <- as.numeric(cm$byClass["Specificity"])
  prec <- as.numeric(cm$byClass["Precision"])
  f1   <- as.numeric(cm$byClass["F1"])
  fdr  <- 1 - prec

  data.frame(
    Modelo = model_name,
    Exactitud = acc,
    Sensibilidad = sens,
    Especificidad = spec,
    Precision = prec,
    FDR = fdr,
    F1_Score = f1,
    stringsAsFactors = FALSE
  )
}

#' Esquema Canonico de Columnas del Benchmark
#'
#' Orden unico de columnas que deben presentar TODAS las tablas de benchmark
#' (escenario balanceado, desbalanceado, cortes Youden y ensamble heuristico).
#'
#' @format Vectores de caracteres.
#' @export
BENCHMARK_METRIC_COLS <- c(
  "Exactitud", "Exactitud_Balanceada", "Sensibilidad", "Especificidad",
  "Precision", "VPN", "FDR", "F1_Score", "MCC", "Kappa", "ROC_AUC"
)
BENCHMARK_COUNT_COLS <- c("TP", "FP", "FN", "TN")
BENCHMARK_COLS <- c("Modelo", BENCHMARK_METRIC_COLS, BENCHMARK_COUNT_COLS)

#' Calcular el Conjunto Completo de Metricas de Clasificacion Binaria
#'
#' Funcion unica y canonica para construir una fila de benchmark. Calcula desde la
#' matriz de confusion (conteos TP/FP/FN/TN) la Exactitud, Exactitud Balanceada,
#' Sensibilidad, Especificidad, Precision (VPP), VPN, FDR, F1-Score, MCC y Kappa de
#' Cohen; y, si se entrega un puntaje continuo, el ROC-AUC (direccion fija: a mayor
#' puntaje, mayor probabilidad de la clase positiva).
#'
#' Convenciones para casos degenerados: si el modelo no predice ningun positivo,
#' Precision y FDR quedan en \code{NA} (indefinidos), mientras que F1 y MCC valen 0.
#'
#' @param pred Vector (factor o caracter) con la clase predicha.
#' @param obs Vector (factor o caracter) con la clase real (Ground Truth).
#' @param score Vector numerico opcional (probabilidad o puntaje) para calcular el ROC-AUC.
#' @param positive_class Etiqueta de la clase positiva (por defecto \code{"TP"}).
#' @param model_name Nombre del modelo/estrategia a registrar.
#'
#' @return \code{data.frame} de una fila con las columnas \code{BENCHMARK_COLS}.
#' @export
compute_classification_metrics <- function(pred,
                                           obs,
                                           score = NULL,
                                           positive_class = "TP",
                                           model_name = "Modelo") {
  pred_chr <- as.character(pred)
  obs_chr  <- as.character(obs)
  stopifnot(length(pred_chr) == length(obs_chr))
  keep <- !is.na(pred_chr) & !is.na(obs_chr)
  if (!is.null(score)) {
    stopifnot(length(score) == length(obs_chr))
    keep <- keep & !is.na(score)
    score <- as.numeric(score)[keep]
  }
  pred_chr <- pred_chr[keep]
  obs_chr  <- obs_chr[keep]

  is_pos_pred <- pred_chr == positive_class
  is_pos_obs  <- obs_chr == positive_class

  tp <- sum(is_pos_pred & is_pos_obs)
  fp <- sum(is_pos_pred & !is_pos_obs)
  fn <- sum(!is_pos_pred & is_pos_obs)
  tn <- sum(!is_pos_pred & !is_pos_obs)
  n  <- tp + fp + fn + tn

  safe_div <- function(a, b) if (b > 0) a / b else NA_real_

  sens <- safe_div(tp, tp + fn)
  spec <- safe_div(tn, tn + fp)
  prec <- safe_div(tp, tp + fp)
  npv  <- safe_div(tn, tn + fn)
  acc  <- safe_div(tp + tn, n)
  bacc <- mean(c(sens, spec), na.rm = TRUE)
  fdr  <- if (is.na(prec)) NA_real_ else 1 - prec
  f1   <- if ((2 * tp + fp + fn) > 0) 2 * tp / (2 * tp + fp + fn) else 0

  mcc_den <- sqrt(as.numeric(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- if (mcc_den > 0) (as.numeric(tp) * tn - as.numeric(fp) * fn) / mcc_den else 0

  p_obs <- acc
  p_exp <- if (n > 0) (as.numeric(tp + fp) * (tp + fn) + as.numeric(fn + tn) * (fp + tn)) / n^2 else NA_real_
  kappa <- if (!is.na(p_exp) && p_exp < 1) (p_obs - p_exp) / (1 - p_exp) else 0

  auc <- NA_real_
  neg_class <- setdiff(unique(obs_chr), positive_class)[1]
  if (!is.null(score) && any(is_pos_obs) && !is.na(neg_class)) {
    roc_obj <- tryCatch(
      suppressMessages(pROC::roc(
        response = obs_chr, predictor = score,
        levels = c(neg_class, positive_class), direction = "<", quiet = TRUE
      )),
      error = function(e) NULL
    )
    if (!is.null(roc_obj)) auc <- as.numeric(roc_obj$auc)
  }

  data.frame(
    Modelo = model_name,
    Exactitud = acc,
    Exactitud_Balanceada = bacc,
    Sensibilidad = sens,
    Especificidad = spec,
    Precision = prec,
    VPN = npv,
    FDR = fdr,
    F1_Score = f1,
    MCC = mcc,
    Kappa = kappa,
    ROC_AUC = auc,
    TP = as.integer(tp),
    FP = as.integer(fp),
    FN = as.integer(fn),
    TN = as.integer(tn),
    stringsAsFactors = FALSE
  )[, BENCHMARK_COLS]
}

#' Calcular Umbral Optimo segun Indice de Youden
#'
#' Determina el punto de corte continuo optimo para discriminar entre verdaderos
#' y falsos positivos maximizando el indice J de Youden (Sensibilidad + Especificidad - 1).
#'
#' @param response Vector categorico o factor con la clase real.
#' @param predictor Vector numerico continuo con la metrica o probabilidad a evaluar.
#' @param positive_class Nivel correspondiente a la clase positiva (por defecto "TRUE_POSITIVE").
#'
#' @return Una lista con:
#'   \item{threshold}{Punto de corte numerico optimo.}
#'   \item{sensitivity}{Sensibilidad en el umbral optimo.}
#'   \item{specificity}{Especificidad en el umbral optimo.}
#'   \item{auc}{Area bajo la curva ROC (AUC).}
#'   \item{roc_object}{Objeto de clase \code{roc} de pROC.}
#' @export
#' @examples
#' \dontrun{
#' youden_res <- calculate_youden_threshold(df$Validacion, df$pct_recovered)
#' }
calculate_youden_threshold <- function(response, predictor, positive_class = "TRUE_POSITIVE") {
  stopifnot(length(response) == length(predictor))
  val_idx <- !is.na(response) & !is.na(predictor)
  resp <- response[val_idx]
  pred <- predictor[val_idx]

  roc_obj <- suppressMessages(pROC::roc(response = resp, predictor = pred, quiet = TRUE))
  best_coords <- pROC::coords(
    roc_obj, "best",
    ret = c("threshold", "specificity", "sensitivity"),
    best.method = "youden"
  )

  list(
    threshold = as.numeric(best_coords$threshold[1]),
    sensitivity = as.numeric(best_coords$sensitivity[1]),
    specificity = as.numeric(best_coords$specificity[1]),
    auc = as.numeric(roc_obj$auc),
    roc_object = roc_obj
  )
}

#' Evaluar Cortes de Youden sobre Multiples Metricas
#'
#' Aplica de manera iterativa el calculo de cortes optimos de Youden para una
#' lista de metricas continuas, retornando la tabla de votos binarios y el benchmark.
#'
#' @param df Data frame con datos anotados.
#' @param response_col Nombre de la columna de respuesta.
#' @param metric_cols Vector con los nombres de las columnas numericas continuas.
#' @param positive_class Nivel de la clase positiva (por defecto "TRUE_POSITIVE").
#'
#' @return Una lista con:
#'   \item{thresholds}{Vector con los umbrales optimos por metrica.}
#'   \item{binary_votes}{Data frame con valores 1/0 segun cada corte.}
#'   \item{benchmark}{Data frame consolidado con las metricas de rendimiento.}
#'   \item{roc_objects}{Lista con los objetos ROC generados.}
#' @export
evaluate_youden_predictors <- function(df,
                                       response_col,
                                       metric_cols,
                                       positive_class = "TRUE_POSITIVE") {
  stopifnot(is.data.frame(df))
  stopifnot(response_col %in% colnames(df))
  stopifnot(all(metric_cols %in% colnames(df)))

  df_clean <- df %>%
    dplyr::select(dplyr::all_of(c(response_col, metric_cols))) %>%
    stats::na.omit()

  thresholds <- numeric(length(metric_cols))
  names(thresholds) <- metric_cols
  df_votos <- data.frame(matrix(ncol = 0, nrow = nrow(df_clean)))
  benchmark_list <- list()
  roc_list <- list()

  for (var in metric_cols) {
    res_y <- calculate_youden_threshold(
      response = df_clean[[response_col]],
      predictor = df_clean[[var]],
      positive_class = positive_class
    )
    umbral <- res_y$threshold
    thresholds[var] <- umbral
    roc_list[[var]] <- res_y$roc_object

    voto_bin <- ifelse(df_clean[[var]] >= umbral, 1, 0)
    df_votos[[var]] <- voto_bin

    neg_class <- setdiff(unique(as.character(df_clean[[response_col]])), positive_class)[1]
    prediccion <- ifelse(voto_bin == 1, positive_class, neg_class)

    benchmark_list[[var]] <- compute_classification_metrics(
      pred = prediccion,
      obs = df_clean[[response_col]],
      score = df_clean[[var]],
      positive_class = positive_class,
      model_name = paste("Corte:", var)
    )
  }

  list(
    thresholds = thresholds,
    binary_votes = df_votos,
    benchmark = dplyr::bind_rows(benchmark_list),
    roc_objects = roc_list
  )
}

#' Optimizar Ensamble Heuristico por Voto Ponderado
#'
#' Explora la grilla de pesos no negativos para los votos binarios y selecciona
#' la combinacion que maximiza el area bajo la curva ROC (AUC), determinando el
#' umbral global optimo de discriminacion.
#'
#' @param df_votes Data frame con votos binarios (1/0) por metrica.
#' @param response Vector de respuesta con la verdad fundamental.
#' @param metric_cols Columnas a incluir en el modelo ponderado (por defecto todas las de df_votes).
#' @param weight_range Rango de enteros para el grid search (por defecto 0:3).
#' @param positive_class Nivel de la clase positiva (por defecto "TRUE_POSITIVE").
#'
#' @return Una lista con:
#'   \item{weights}{Vector con los mejores pesos normalizados/enteros.}
#'   \item{best_auc}{Maximo AUC alcanzado en la grilla.}
#'   \item{global_threshold}{Umbral optimo global para el puntaje acumulado.}
#'   \item{best_scores}{Vector con los puntajes ponderados resultantes.}
#'   \item{roc_object}{Objeto ROC del modelo ponderado.}
#'   \item{benchmark_row}{Data frame con las metricas del modelo combinado.}
#' @export
optimize_weighted_heuristic <- function(df_votes,
                                        response,
                                        metric_cols = colnames(df_votes),
                                        weight_range = 0:3,
                                        positive_class = "TRUE_POSITIVE") {
  stopifnot(is.data.frame(df_votes))
  stopifnot(length(response) == nrow(df_votes))

  grid_pesos <- expand.grid(rep(list(weight_range), length(metric_cols)))
  colnames(grid_pesos) <- metric_cols
  grid_pesos <- grid_pesos[rowSums(grid_pesos) > 0, , drop = FALSE]

  mat_votos <- as.matrix(df_votes[, metric_cols])
  resp_fac <- factor(response)

  mejor_auc <- 0
  mejores_pesos <- NULL
  mejor_puntaje <- NULL

  for (i in seq_len(nrow(grid_pesos))) {
    pesos_act <- as.numeric(grid_pesos[i, ])
    score_act <- as.vector(mat_votos %*% pesos_act)
    roc_act <- suppressMessages(pROC::roc(response = resp_fac, predictor = score_act, quiet = TRUE))
    if (roc_act$auc > mejor_auc) {
      mejor_auc <- as.numeric(roc_act$auc)
      mejores_pesos <- pesos_act
      mejor_puntaje <- score_act
    }
  }

  names(mejores_pesos) <- metric_cols
  roc_comb <- suppressMessages(pROC::roc(response = resp_fac, predictor = mejor_puntaje, quiet = TRUE))
  umbral_comb <- as.numeric(pROC::coords(roc_comb, "best", ret = "threshold", best.method = "youden")$threshold[1])

  neg_class <- setdiff(levels(resp_fac), positive_class)[1]
  pred_comb <- factor(
    ifelse(mejor_puntaje >= umbral_comb, positive_class, neg_class),
    levels = levels(resp_fac)
  )

  benchmark_row <- compute_classification_metrics(
    pred = pred_comb,
    obs = resp_fac,
    score = mejor_puntaje,
    positive_class = positive_class,
    model_name = "Modelo Combinado (Heuristico)"
  )

  list(
    weights = mejores_pesos,
    best_auc = mejor_auc,
    global_threshold = umbral_comb,
    best_scores = mejor_puntaje,
    roc_object = roc_comb,
    benchmark_row = benchmark_row
  )
}

#' Predecir con Modelo Heuristico Ponderado
#'
#' Aplica los umbrales univariados y el vector de ponderaciones optimizado sobre
#' nuevas observaciones continuas para inferir la clase "TP" o "FP".
#'
#' @param df_predictors Data frame con las columnas continuas requeridas.
#' @param thresholds Vector nombrado con los umbrales univariados.
#' @param weights Vector nombrado con los pesos asignados a cada metrica.
#' @param global_threshold Umbral escalar de decision para el puntaje ponderado.
#'
#' @return Un factor con niveles c("TP", "FP").
#' @export
predict_weighted_heuristic <- function(df_predictors, thresholds, weights, global_threshold) {
  metricas <- names(weights)
  stopifnot(all(metricas %in% colnames(df_predictors)))
  stopifnot(all(metricas %in% names(thresholds)))

  mat_votos <- matrix(0, nrow = nrow(df_predictors), ncol = length(metricas))
  colnames(mat_votos) <- metricas

  for (m in metricas) {
    mat_votos[, m] <- ifelse(df_predictors[[m]] >= thresholds[[m]], 1, 0)
  }

  puntaje <- as.vector(mat_votos %*% weights[metricas])
  factor(ifelse(puntaje >= global_threshold, "TP", "FP"), levels = c("TP", "FP"))
}

#' Consolidar Sistema de Consenso Multimodelo
#'
#' Agrega las decisiones de prediccion provenientes de multiples modelos (Top ML,
#' Cortes Youden y Ensamble Heuristico) y dictamina el veredicto final por regla de mayoria.
#'
#' @param df_predictions Data frame que contiene las predicciones categoricas ("TP" / "FP").
#' @param model_cols Vector con los nombres de las columnas de prediccion a someter a consenso.
#' @param min_votes Umbral minimo de votos TP para validar la secuencia homóloga (por defecto 4).
#'
#' @return El \code{data.frame} de entrada complementado con las columnas \code{Votos_TP}
#'   y \code{Veredicto_Consenso} ("Homologo_Valido" o "Ruido_Descartado").
#' @export
consolidar_consenso <- function(df_predictions, model_cols, min_votes = 4) {
  stopifnot(is.data.frame(df_predictions))
  stopifnot(all(model_cols %in% colnames(df_predictions)))

  # Suma de votos afirmativos por fila
  votos_matrix <- as.matrix(df_predictions[, model_cols, drop = FALSE]) == "TP"
  votos_totales <- rowSums(votos_matrix, na.rm = TRUE)

  df_predictions$Votos_TP <- votos_totales
  df_predictions$Veredicto_Consenso <- dplyr::if_else(
    votos_totales >= min_votes,
    "Homologo_Valido",
    "Ruido_Descartado"
  )

  return(df_predictions)
}

#' Pipeline Completo de Prediccion e Inferencia de CNV
#'
#' Orquesta de extremo a extremo la aplicacion de modelos supervisados (ML), cortes
#' univariados de Youden y ensamble heuristico ponderado sobre un conjunto de datos nuevo,
#' consolidando el veredicto de consenso.
#'
#' @param newdata Data frame con nuevas secuencias candidatas y metricas continuas.
#' @param ml_models Lista nombrada con modelos de \code{caret} entrenados (ej. ElasticNet, SVM, etc.).
#' @param youden_vars Vector con los nombres de las variables para cortes Youden.
#' @param youden_thresholds Vector nombrado con los umbrales optimos correspondientes.
#' @param heuristic_weights Vector nombrado con los pesos del ensamble heuristico.
#' @param heuristic_threshold Umbral escalar del ensamble heuristico.
#' @param min_votes Minimo de votos TP para el veredicto de mayoria (por defecto 4).
#'
#' @return Un \code{data.frame} con las predicciones individuales, los votos totales y el veredicto final.
#' @export
predict_cnv_pipeline <- function(newdata,
                                 ml_models,
                                 youden_vars,
                                 youden_thresholds,
                                 heuristic_weights,
                                 heuristic_threshold,
                                 min_votes = 4) {
  stopifnot(is.data.frame(newdata))
  metricas_req <- names(heuristic_weights)
  stopifnot(all(metricas_req %in% colnames(newdata)))

  df_pred <- newdata
  model_cols <- character()

  # 1. Modelos ML supervisados
  if (!is.null(ml_models) && length(ml_models) > 0) {
    for (m_name in names(ml_models)) {
      col_name <- paste0("Pred_ML_", m_name)
      df_pred[[col_name]] <- stats::predict(ml_models[[m_name]], newdata = newdata[, metricas_req, drop = FALSE])
      # Estandarizar niveles a TP / FP si vienen como TRUE_POSITIVE / FALSE_POSITIVE
      if (all(c("TRUE_POSITIVE", "FALSE_POSITIVE") %in% levels(df_pred[[col_name]]))) {
        df_pred[[col_name]] <- factor(
          ifelse(df_pred[[col_name]] == "TRUE_POSITIVE", "TP", "FP"),
          levels = c("TP", "FP")
        )
      }
      model_cols <- c(model_cols, col_name)
    }
  }

  # 2. Cortes Univariados Youden
  for (y_var in youden_vars) {
    col_name <- paste0("Pred_Youden_", y_var)
    umbral <- youden_thresholds[[y_var]]
    df_pred[[col_name]] <- factor(
      ifelse(df_pred[[y_var]] >= umbral, "TP", "FP"),
      levels = c("TP", "FP")
    )
    model_cols <- c(model_cols, col_name)
  }

  # 3. Ensamble Heuristico Ponderado
  col_ens <- "Pred_Ens_Heuristico"
  df_pred[[col_ens]] <- predict_weighted_heuristic(
    df_predictors = df_pred,
    thresholds = youden_thresholds,
    weights = heuristic_weights,
    global_threshold = heuristic_threshold
  )
  model_cols <- c(model_cols, col_ens)

  # 4. Consolidacion de consenso
  df_final <- consolidar_consenso(
    df_predictions = df_pred,
    model_cols = model_cols,
    min_votes = min_votes
  )

  return(df_final)
}
