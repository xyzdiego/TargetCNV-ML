# ==============================================================================
# TargetCNV-ML: Evaluacion Diagnostica, Cortes de Youden y Consenso (R/model_evaluation.R)
# ==============================================================================

# Cargar especificaciones compartidas si existen en el mismo directorio R/
local({
  specs_candidates <- c(
    file.path("R", "model_specs.R"),
    "model_specs.R"
  )
  for (cand in specs_candidates) {
    if (file.exists(cand)) {
      source(cand, local = FALSE)
      break
    }
  }
})

#' Calcula el esquema canonico de 16 columnas de metricas de clasificacion binaria
#'
#' Estandariza las tablas de benchmark para modelos supervisados (escenarios
#' balanceado y desbalanceado), cortes univariados de Youden y ensamble heuristico.
#' Devuelve NA_real_ cuando el denominador de una metrica es 0 y fija la direccion
#' de la curva ROC como `levels = c("FP", "TP"), direction = "<"`.
#'
#' @param pred Vector de clases predichas ("TP" o "FP").
#' @param obs Vector de clases observadas reales ("TP"/"TRUE_POSITIVE" o "FP"/"FALSE_POSITIVE").
#' @param positive_class Valor que representa la clase positiva en `obs`.
#' @param model_name Nombre del modelo o filtro evaluado.
#' @param score Vector numerico opcional de probabilidades o puntajes continuos para calcular ROC_AUC.
#' @return Data frame con las 16 columnas estandarizadas.
#' @export
compute_classification_metrics <- function(pred,
                                           obs,
                                           positive_class = "TRUE_POSITIVE",
                                           model_name = "Modelo",
                                           score = NULL) {
  pred_clean <- factor(
    dplyr::if_else(as.character(pred) %in% c("TP", "TRUE_POSITIVE", positive_class), "TP", "FP"),
    levels = c("TP", "FP")
  )
  obs_clean <- factor(
    dplyr::if_else(as.character(obs) %in% c("TP", "TRUE_POSITIVE", positive_class), "TP", "FP"),
    levels = c("TP", "FP")
  )

  tp <- sum(pred_clean == "TP" & obs_clean == "TP", na.rm = TRUE)
  fp <- sum(pred_clean == "TP" & obs_clean == "FP", na.rm = TRUE)
  fn <- sum(pred_clean == "FP" & obs_clean == "TP", na.rm = TRUE)
  tn <- sum(pred_clean == "FP" & obs_clean == "FP", na.rm = TRUE)

  total <- tp + fp + fn + tn
  acc     <- if (total > 0) (tp + tn) / total else NA_real_
  sens    <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  spec    <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  prec    <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  vpn     <- if ((tn + fn) > 0) tn / (tn + fn) else NA_real_
  fdr     <- if (!is.na(prec)) 1 - prec else NA_real_
  f1      <- if (!is.na(prec) && !is.na(sens) && (prec + sens) > 0) 2 * (prec * sens) / (prec + sens) else NA_real_
  bal_acc <- if (!is.na(sens) && !is.na(spec)) (sens + spec) / 2 else NA_real_

  num_mcc <- (as.numeric(tp) * as.numeric(tn)) - (as.numeric(fp) * as.numeric(fn))
  den_mcc <- sqrt(as.numeric(tp + fp) * as.numeric(tp + fn) * as.numeric(tn + fp) * as.numeric(tn + fn))
  mcc <- if (!is.na(den_mcc) && den_mcc > 0) num_mcc / den_mcc else NA_real_

  p_o <- acc
  p_e <- if (total > 0) (((as.numeric(tp + fp) * as.numeric(tp + fn)) + (as.numeric(fn + tn) * as.numeric(fp + tn))) / (as.numeric(total)^2)) else NA_real_
  kappa <- if (!is.na(p_e) && (1 - p_e) > 0) (p_o - p_e) / (1 - p_e) else NA_real_

  auc_val <- NA_real_
  if (!is.null(score) && length(score) == length(obs_clean)) {
    tryCatch({
      obs_roc <- factor(as.character(obs_clean), levels = c("FP", "TP"))
      roc_obj <- suppressMessages(
        pROC::roc(
          response = obs_roc,
          predictor = as.numeric(score),
          levels = c("FP", "TP"),
          direction = "<",
          quiet = TRUE
        )
      )
      auc_val <- as.numeric(roc_obj$auc)
    }, error = function(e) {})
  }

  data.frame(
    Modelo = model_name,
    Exactitud = acc,
    Exactitud_Balanceada = bal_acc,
    Sensibilidad = sens,
    Especificidad = spec,
    Precision = prec,
    VPN = vpn,
    FDR = fdr,
    F1_Score = f1,
    MCC = mcc,
    Kappa = kappa,
    ROC_AUC = auc_val,
    TP = as.integer(tp),
    FP = as.integer(fp),
    FN = as.integer(fn),
    TN = as.integer(tn),
    stringsAsFactors = FALSE
  )
}

#' Calcula el umbral optimo de Youden (J = Sensibilidad + Especificidad - 1)
#'
#' Fija explicitamente `levels = c("FP", "TP"), direction = "<"` para garantizar
#' consistencia con la regla de decision `predictor >= threshold -> TP`.
#'
#' @param response Vector con la variable respuesta real.
#' @param predictor Vector numerico continuo de la metrica a evaluar.
#' @param positive_class Etiqueta de la clase positiva.
#' @return Lista con `threshold`, `sensitivity`, `specificity`, `auc` y `roc_object`.
#' @export
calculate_youden_threshold <- function(response, predictor, positive_class = "TRUE_POSITIVE") {
  stopifnot(length(response) == length(predictor))
  val_idx <- !is.na(response) & !is.na(predictor)
  resp_raw <- response[val_idx]
  resp <- factor(
    dplyr::if_else(as.character(resp_raw) %in% c("TP", "TRUE_POSITIVE", positive_class), "TP", "FP"),
    levels = c("FP", "TP")
  )
  pred <- as.numeric(predictor[val_idx])

  roc_obj <- suppressMessages(
    pROC::roc(response = resp, predictor = pred, levels = c("FP", "TP"), direction = "<", quiet = TRUE)
  )
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

#' Evalua multiples predictores continuos mediante el Indice J de Youden
#'
#' @param df Data frame de entrenamiento.
#' @param response_col Nombre de la columna de respuesta binaria.
#' @param metric_cols Vector de nombres de las metricas continuas.
#' @param positive_class Etiqueta de la clase positiva.
#' @return Lista con `thresholds`, `binary_votes`, `benchmark` y `roc_objects`.
#' @export
evaluate_youden_predictors <- function(df, response_col, metric_cols, positive_class = "TRUE_POSITIVE") {
  stopifnot(is.data.frame(df))
  stopifnot(response_col %in% colnames(df))
  stopifnot(all(metric_cols %in% colnames(df)))

  thresh_list <- list()
  votes_df <- data.frame(matrix(ncol = 0, nrow = nrow(df)))
  bench_list <- list()
  roc_list <- list()

  for (m in metric_cols) {
    y_res <- calculate_youden_threshold(df[[response_col]], df[[m]], positive_class = positive_class)
    thresh_list[[m]] <- y_res$threshold
    roc_list[[m]] <- y_res$roc_object

    vote_m <- dplyr::if_else(df[[m]] >= y_res$threshold, "TP", "FP")
    votes_df[[m]] <- vote_m

    m_bench <- compute_classification_metrics(
      pred = vote_m,
      obs = df[[response_col]],
      positive_class = positive_class,
      model_name = paste0("Youden_", m),
      score = df[[m]]
    )
    bench_list[[m]] <- m_bench
  }

  return(list(
    thresholds = unlist(thresh_list),
    binary_votes = votes_df,
    benchmark = dplyr::bind_rows(bench_list),
    roc_objects = roc_list
  ))
}

#' Optimiza por busqueda exhaustiva en grilla los pesos del Ensamble Heuristico Ponderado
#'
#' @param df_votes Data frame con votos binarios ("TP"/"FP") por cada metrica.
#' @param response Vector con las etiquetas reales.
#' @param metric_cols Nombres de las columnas de metricas evaluadas.
#' @param weight_range Rango entero de pesos a explorar (por defecto `0:3`).
#' @param positive_class Etiqueta de la clase positiva.
#' @return Lista con `weights`, `best_auc`, `global_threshold`, `best_scores`, `roc_object` y `benchmark_row`.
#' @export
optimize_weighted_heuristic <- function(df_votes,
                                        response,
                                        metric_cols = colnames(df_votes),
                                        weight_range = 0:3,
                                        positive_class = "TRUE_POSITIVE") {
  grid <- expand.grid(rep(list(weight_range), length(metric_cols)))
  colnames(grid) <- metric_cols
  grid <- grid[rowSums(grid) > 0, ]

  resp_fac <- factor(
    dplyr::if_else(as.character(response) %in% c("TP", "TRUE_POSITIVE", positive_class), "TP", "FP"),
    levels = c("FP", "TP")
  )

  vote_mat <- as.matrix(df_votes[, metric_cols, drop = FALSE] == "TP")
  best_auc <- -1
  best_weights <- NULL
  best_scores <- NULL

  for (i in seq_len(nrow(grid))) {
    w <- as.numeric(grid[i, ])
    scores <- as.vector(vote_mat %*% w)
    roc_i <- suppressMessages(
      pROC::roc(response = resp_fac, predictor = scores, levels = c("FP", "TP"), direction = "<", quiet = TRUE)
    )
    auc_i <- as.numeric(roc_i$auc)

    if (auc_i > best_auc) {
      best_auc <- auc_i
      best_weights <- w
      best_scores <- scores
    }
  }

  names(best_weights) <- metric_cols
  final_y <- calculate_youden_threshold(response, best_scores, positive_class = positive_class)
  final_pred <- dplyr::if_else(best_scores >= final_y$threshold, "TP", "FP")

  bench_row <- compute_classification_metrics(
    pred = final_pred,
    obs = response,
    positive_class = positive_class,
    model_name = "Weighted_Heuristic",
    score = best_scores
  )

  list(
    weights = best_weights,
    best_auc = best_auc,
    global_threshold = final_y$threshold,
    best_scores = best_scores,
    roc_object = final_y$roc_object,
    benchmark_row = bench_row
  )
}

#' Predice clases binarias aplicando el Ensamble Heuristico Ponderado
#'
#' @param df_predictors Data frame con las metricas continuas de nuevas muestras.
#' @param thresholds Vector nombrado de umbrales univariados de Youden.
#' @param weights Vector nombrado de pesos optimos por metrica.
#' @param global_threshold Umbral global de corte sobre la suma ponderada.
#' @return Factor con niveles `c("TP", "FP")`.
#' @export
predict_weighted_heuristic <- function(df_predictors, thresholds, weights, global_threshold) {
  stopifnot(all(names(thresholds) %in% colnames(df_predictors)))
  stopifnot(all(names(weights) %in% names(thresholds)))

  m_names <- names(weights)
  vote_mat <- matrix(0, nrow = nrow(df_predictors), ncol = length(m_names))

  for (j in seq_along(m_names)) {
    m <- m_names[j]
    vote_mat[, j] <- as.numeric(df_predictors[[m]] >= thresholds[m])
  }

  scores <- as.vector(vote_mat %*% weights)
  preds  <- dplyr::if_else(scores >= global_threshold, "TP", "FP")
  return(factor(preds, levels = c("TP", "FP")))
}

#' Consolida la votacion multimodelo (Sistema de Consenso de 7 Modelos)
#'
#' @param df_predictions Data frame con las columnas de prediccion de cada modelo votante.
#' @param model_cols Vector con los nombres de las columnas votantes.
#' @param min_votes Numero minimo de votos afirmativos ("TP") para declarar `Homologo_Valido`.
#' @return Data frame enriquecido con `Votos_TP`, `Proporcion_Acuerdo` y `Veredicto_Consenso`.
#' @export
consolidar_consenso <- function(df_predictions, model_cols, min_votes = 4) {
  stopifnot(is.data.frame(df_predictions))
  stopifnot(all(model_cols %in% colnames(df_predictions)))

  n_models <- length(model_cols)
  votos_matrix <- as.matrix(df_predictions[, model_cols, drop = FALSE]) == "TP"
  votos_totales <- rowSums(votos_matrix, na.rm = TRUE)

  df_predictions$Votos_TP <- as.integer(votos_totales)
  df_predictions$Proporcion_Acuerdo <- round(votos_totales / max(1L, n_models), 3)
  df_predictions$Veredicto_Consenso <- dplyr::if_else(
    votos_totales >= min_votes,
    "Homologo_Valido",
    "Ruido_Descartado"
  )

  return(df_predictions)
}
