# Evaluación diagnóstica de clasificadores, colinealidad de Spearman, consenso y concordancia de CNV

local({
  if (!exists("normalize_model_names", mode = "function")) {
    specs_candidates <- c(file.path("R", "model_specs.R"), "model_specs.R")
    for (cand in specs_candidates) {
      if (file.exists(cand)) {
        source(cand, local = FALSE)
        break
      }
    }
  }
})

#' Calcula el esquema estándar de 16 métricas de clasificación binaria
#'
#' Estandariza las tablas de desempeño para modelos supervisados, cortes univariados
#' de Youden, ensamble heurístico y consenso multimodelo. Devuelve `NA_real_` cuando
#' el denominador de una tasa es cero y fija la orientación ROC en `levels = c("FP", "TP")`.
#'
#' @param pred Vector de clases predichas (`"TP"` o `"FP"`).
#' @param obs Vector de clases observadas reales (`"TP"`/`"TRUE_POSITIVE"` o `"FP"`/`"FALSE_POSITIVE"`).
#' @param positive_class Etiqueta que representa la clase positiva en `obs`.
#' @param model_name Nombre identificador del modelo o regla evaluada.
#' @param score Vector numérico opcional de probabilidades o puntajes continuos para `ROC_AUC`.
#' @return Data frame de una fila con las 16 columnas de evaluación.
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

#' Calcula el punto de corte óptimo de Youden (J = Sensibilidad + Especificidad - 1)
#'
#' @param response Vector con la variable respuesta binaria real.
#' @param predictor Vector numérico continuo de la métrica a evaluar.
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

#' Evalúa múltiples predictores continuos mediante el índice J de Youden
#'
#' @param df Data frame de entrenamiento.
#' @param response_col Nombre de la columna de respuesta binaria.
#' @param metric_cols Vector con los nombres de las métricas continuas.
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

    bench_list[[m]] <- compute_classification_metrics(
      pred = vote_m,
      obs = df[[response_col]],
      positive_class = positive_class,
      model_name = paste0("Youden_", m),
      score = df[[m]]
    )
  }

  list(
    thresholds = unlist(thresh_list),
    binary_votes = votes_df,
    benchmark = dplyr::bind_rows(bench_list),
    roc_objects = roc_list
  )
}

#' Optimiza por búsqueda en grilla los pesos del ensamble heurístico ponderado
#'
#' @param df_votes Data frame con votos binarios (`"TP"`/`"FP"`) por cada métrica.
#' @param response Vector con las etiquetas reales.
#' @param metric_cols Nombres de las columnas de métricas evaluadas.
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

#' Predice clases binarias aplicando el ensamble heurístico ponderado
#'
#' @param df_predictors Data frame con las métricas continuas de nuevas muestras.
#' @param thresholds Vector nombrado de umbrales univariados de Youden.
#' @param weights Vector nombrado de pesos óptimos por métrica.
#' @param global_threshold Umbral global de decisión sobre la suma ponderada.
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
  factor(preds, levels = c("TP", "FP"))
}

#' Extrae probabilidades calibradas de `TP` desde un meta-modelo `caretStack`
#'
#' Alinea la columna de probabilidad de los modelos base con la convención interna
#' de `caretEnsemble::makePredObsMatrix()` (que alimenta al meta-modelo GLM con el
#' segundo nivel del factor `levels = c("TP", "FP")`, es decir `"FP"`), evitando la
#' inversión de probabilidades al predecir sobre nuevos conjuntos de datos.
#'
#' @param ens_stack Objeto de clase `caretStack`.
#' @param newdata Data frame con los predictores continuos de evaluación.
#' @return Vector numérico en `[0, 1]` con la probabilidad estimada de `TP`.
#' @export
predict_stacking_probs <- function(ens_stack, newdata) {
  stopifnot(inherits(ens_stack, "caretStack"))
  base_mods <- ens_stack$models
  pred_names <- names(base_mods)
  x_df <- as.data.frame(newdata)

  td <- ens_stack$ens_model$trainingData
  input_level <- "FP"
  if (!is.null(td) && ".outcome" %in% colnames(td) && all(pred_names %in% colnames(td))) {
    m1 <- pred_names[1]
    mean_tp <- mean(td[[m1]][td$.outcome == "TP"], na.rm = TRUE)
    mean_fp <- mean(td[[m1]][td$.outcome == "FP"], na.rm = TRUE)
    if (!is.na(mean_tp) && !is.na(mean_fp)) {
      input_level <- if (mean_tp < mean_fp) "FP" else "TP"
    }
  } else if (length(base_mods) > 0 && length(base_mods[[1]]$levels) >= 2) {
    input_level <- base_mods[[1]]$levels[2]
  }

  prob_base <- lapply(base_mods, function(bm) {
    pr <- stats::predict(bm, newdata = x_df, type = "prob")
    col_use <- if (input_level %in% colnames(pr)) input_level else colnames(pr)[ncol(pr)]
    as.numeric(pr[[col_use]])
  })
  df_prob_base <- as.data.frame(prob_base, stringsAsFactors = FALSE)
  colnames(df_prob_base) <- pred_names

  meta_probs <- stats::predict(ens_stack$ens_model, newdata = df_prob_base, type = "prob")
  as.numeric(meta_probs[["TP"]])
}

#' Consolida la votación multimodelo del sistema de consenso
#'
#' @param df_predictions Data frame con las columnas de predicción de cada modelo votante.
#' @param model_cols Vector con los nombres de las columnas votantes.
#' @param min_votes Número mínimo de votos afirmativos (`"TP"`) para declarar `Homologo_Valido`.
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

  df_predictions
}

#' Evalúa la colinealidad por rangos de Spearman entre las métricas continuas de CAPTUS
#'
#' Calcula la matriz simétrica de coeficientes de correlación de Spearman (`rho`),
#' los p-valores bilaterales exactos/asintóticos y un resumen tabular por pares para
#' diagnosticar la redundancia estadística entre predictores antes del modelado.
#'
#' @param df Data frame con las variables predictoras numéricas.
#' @param metric_cols Vector de nombres de columnas continuas a evaluar.
#' @return Lista con `cor_matrix`, `p_matrix` y `pairs` (data frame por pares).
#' @export
compute_spearman_collinearity <- function(df,
                                          metric_cols = c("pct_recovered", "pct_identity", "score", "wscore")) {
  stopifnot(is.data.frame(df))
  cols <- intersect(metric_cols, colnames(df))
  stopifnot(length(cols) >= 2L)

  mat_num <- as.data.frame(lapply(df[, cols, drop = FALSE], function(v) as.numeric(as.character(v))))
  k <- length(cols)
  cor_mat <- diag(1.0, nrow = k, ncol = k)
  p_mat   <- matrix(0.0, nrow = k, ncol = k)
  dimnames(cor_mat) <- list(cols, cols)
  dimnames(p_mat)   <- list(cols, cols)

  pair_rows <- list()
  idx <- 0L
  for (i in seq_len(k - 1L)) {
    for (j in seq.int(i + 1L, k)) {
      v1 <- mat_num[[i]]
      v2 <- mat_num[[j]]
      ok <- !is.na(v1) & !is.na(v2)
      n_ok <- sum(ok)
      rho_ij <- NA_real_
      pval_ij <- NA_real_
      if (n_ok >= 3L && stats::sd(v1[ok]) > 0 && stats::sd(v2[ok]) > 0) {
        ct <- suppressWarnings(stats::cor.test(v1[ok], v2[ok], method = "spearman", exact = FALSE))
        rho_ij <- as.numeric(ct$estimate)
        pval_ij <- as.numeric(ct$p.value)
      }
      cor_mat[i, j] <- cor_mat[j, i] <- rho_ij
      p_mat[i, j]   <- p_mat[j, i]   <- pval_ij

      abs_r <- abs(rho_ij)
      nivel <- if (is.na(abs_r)) {
        "Sin varianza"
      } else if (abs_r >= 0.85) {
        "Alta (|rho| >= 0.85)"
      } else if (abs_r >= 0.60) {
        "Moderada (0.60 <= |rho| < 0.85)"
      } else {
        "Baja (|rho| < 0.60)"
      }

      idx <- idx + 1L
      pair_rows[[idx]] <- data.frame(
        Variable_1 = cols[i],
        Variable_2 = cols[j],
        N_Pares = as.integer(n_ok),
        Spearman_Rho = rho_ij,
        P_Valor = pval_ij,
        Colinealidad = nivel,
        stringsAsFactors = FALSE
      )
    }
  }

  pairs_df <- dplyr::bind_rows(pair_rows) %>%
    dplyr::arrange(dplyr::desc(abs(Spearman_Rho)))

  list(
    cor_matrix = cor_mat,
    p_matrix = p_mat,
    pairs = pairs_df
  )
}

#' Calcula el coeficiente Kappa cuadrático ponderado para conteos enteros de copias
#'
#' @param pred_copies Vector entero o numérico con el número de copias estimado por ortogrupo.
#' @param real_copies Vector entero o numérico con el número real de copias por ortogrupo.
#' @return Escalar numérico con el Kappa cuadrático ponderado (`[-1, 1]`).
#' @export
compute_quadratic_weighted_kappa <- function(pred_copies, real_copies) {
  stopifnot(length(pred_copies) == length(real_copies))
  ok <- !is.na(pred_copies) & !is.na(real_copies)
  p_vec <- as.integer(round(as.numeric(pred_copies[ok])))
  r_vec <- as.integer(round(as.numeric(real_copies[ok])))
  n <- length(p_vec)
  if (n == 0L) return(NA_real_)

  min_k <- min(c(0L, p_vec, r_vec))
  max_k <- max(c(1L, p_vec, r_vec))
  states <- seq.int(min_k, max_k)
  m <- length(states)
  if (m <= 1L) return(NA_real_)

  p_fac <- factor(p_vec, levels = states)
  r_fac <- factor(r_vec, levels = states)
  obs_mat <- as.matrix(table(r_fac, p_fac)) / n

  r_marg <- rowSums(obs_mat)
  p_marg <- colSums(obs_mat)
  exp_mat <- outer(r_marg, p_marg)

  span <- as.numeric(max_k - min_k)
  w_mat <- outer(states, states, function(a, b) ((as.numeric(a) - as.numeric(b))^2) / (span^2))

  num <- sum(w_mat * obs_mat)
  den <- sum(w_mat * exp_mat)
  if (den <= 0) {
    return(if (num == 0) 1.0 else NA_real_)
  }
  1.0 - (num / den)
}

#' Extrae el conteo real de copias por muestra y ortogrupo desde un dataset etiquetado
#'
#' Después de resolver redundancias intra-locus, cada registro `TRUE_POSITIVE` corresponde
#' a un locus génico distinto validado espacialmente en el genoma de referencia, incluyendo
#' ortogrupos con cero copias reales (`CNV_Real == 0`).
#'
#' @param df Data frame etiquetado con `Validacion_Espacial`.
#' @param sample_col Columna de muestra (o `"sample_name"` si `"Sample"` no existe).
#' @param locus_col Columna de ortogrupo (o `"locus"` / `"OG_Limpio"`).
#' @param class_col Columna con la verdad de referencia (`"Validacion_Espacial"`).
#' @return Data frame con columnas `Sample`, `Locus_Ortogrupo` y `CNV_Real`.
#' @export
extract_reference_cnv <- function(df,
                                  sample_col = "Sample",
                                  locus_col = "Locus_Ortogrupo",
                                  class_col = "Validacion_Espacial") {
  stopifnot(is.data.frame(df))
  stopifnot(class_col %in% colnames(df))

  s_col <- if (sample_col %in% colnames(df)) {
    sample_col
  } else if ("sample_name" %in% colnames(df)) {
    "sample_name"
  } else {
    NULL
  }
  l_col <- if (locus_col %in% colnames(df)) {
    locus_col
  } else if ("locus" %in% colnames(df)) {
    "locus"
  } else if ("OG_Limpio" %in% colnames(df)) {
    "OG_Limpio"
  } else {
    stop("[ERROR] No se encontró columna de ortogrupo/locus en extract_reference_cnv().", call. = FALSE)
  }

  df_work <- df
  if ("marker_type" %in% colnames(df_work)) {
    df_work <- df_work %>% dplyr::filter(marker_type == "NUC")
  }
  df_work$.Sample_Tmp <- if (!is.null(s_col)) as.character(df_work[[s_col]]) else "Sample_1"
  df_work$.Locus_Tmp  <- as.character(df_work[[l_col]])

  df_work %>%
    dplyr::filter(!is.na(.Locus_Tmp) & nzchar(.Locus_Tmp)) %>%
    dplyr::group_by(Sample = .Sample_Tmp, Locus_Ortogrupo = .Locus_Tmp) %>%
    dplyr::summarise(
      CNV_Real = as.integer(sum(.data[[class_col]] %in% c("TRUE_POSITIVE", "TP"), na.rm = TRUE)),
      .groups = "drop"
    )
}

#' Evalúa la concordancia de número de copias a nivel de ortogrupo (CNV estimado vs. CNV real)
#'
#' Contrasta el número de copias inferido (`CNV_Estimado`) contra el número real de copias
#' (`CNV_Real`) en todos los ortogrupos evaluados (incluyendo ortogrupos con `0` copias),
#' reportando correlación de Spearman (`rho`), Kappa cuadrático ponderado, RMSE, MAE,
#' sesgo medio y tasa de concordancia exacta de copias.
#'
#' @param df_consenso Data frame con veredictos de consenso (`Veredicto_Consenso`) y `pct_recovered`.
#' @param cov_threshold Umbral de cobertura (%) para confirmar copias adicionales (por defecto `50.0`).
#' @param sample_col Columna de muestra en `df_consenso`.
#' @param locus_col Columna de ortogrupo en `df_consenso`.
#' @param verdict_col Columna de veredicto del consenso.
#' @param true_cnv_df Data frame opcional con columnas `Sample`, `Locus_Ortogrupo` y `CNV_Real`.
#' @return Lista con `summary` (métricas globales) y `pairs` (comparación ortogrupo por ortogrupo).
#' @export
evaluate_cnv_concordance <- function(df_consenso,
                                     cov_threshold = 50.0,
                                     sample_col = "Sample",
                                     locus_col = "Locus_Ortogrupo",
                                     verdict_col = "Veredicto_Consenso",
                                     true_cnv_df = NULL) {
  stopifnot(is.data.frame(df_consenso))

  df_work <- df_consenso
  if (!sample_col %in% colnames(df_work)) {
    df_work[[sample_col]] <- if ("sample_name" %in% colnames(df_work)) as.character(df_work$sample_name) else "Sample_1"
  }
  if (!locus_col %in% colnames(df_work)) {
    df_work[[locus_col]] <- if ("locus" %in% colnames(df_work)) as.character(df_work$locus) else as.character(df_work$OG_Limpio)
  }

  ref_df <- if (!is.null(true_cnv_df)) {
    stopifnot(all(c("Sample", "Locus_Ortogrupo", "CNV_Real") %in% colnames(true_cnv_df)))
    true_cnv_df %>%
      dplyr::transmute(
        Sample = as.character(Sample),
        Locus_Ortogrupo = as.character(Locus_Ortogrupo),
        CNV_Real = as.integer(CNV_Real)
      )
  } else if ("Validacion_Espacial" %in% colnames(df_work)) {
    extract_reference_cnv(df_work, sample_col = sample_col, locus_col = locus_col, class_col = "Validacion_Espacial")
  } else {
    stop("[ERROR] evaluate_cnv_concordance() requiere 'Validacion_Espacial' en df_consenso o 'true_cnv_df'.", call. = FALSE)
  }

  cnv_res <- extract_cnv_matrix(
    df = df_work,
    sample_col = sample_col,
    locus_col = locus_col,
    cov_threshold = cov_threshold,
    verdict_col = verdict_col
  )

  pred_counts <- if (nrow(cnv_res$functional_genes) > 0) {
    cnv_res$functional_genes %>%
      dplyr::group_by(
        Sample = as.character(.data[[sample_col]]),
        Locus_Ortogrupo = as.character(.data[[locus_col]])
      ) %>%
      dplyr::summarise(CNV_Estimado = as.integer(dplyr::n()), .groups = "drop")
  } else {
    data.frame(Sample = character(0), Locus_Ortogrupo = character(0), CNV_Estimado = integer(0), stringsAsFactors = FALSE)
  }

  universe <- dplyr::bind_rows(
    ref_df %>% dplyr::select(Sample, Locus_Ortogrupo),
    df_work %>% dplyr::transmute(Sample = as.character(.data[[sample_col]]), Locus_Ortogrupo = as.character(.data[[locus_col]]))
  ) %>%
    dplyr::filter(!is.na(Locus_Ortogrupo) & nzchar(Locus_Ortogrupo)) %>%
    dplyr::distinct()

  pairs_df <- universe %>%
    dplyr::left_join(ref_df, by = c("Sample", "Locus_Ortogrupo")) %>%
    dplyr::left_join(pred_counts, by = c("Sample", "Locus_Ortogrupo")) %>%
    dplyr::filter(!is.na(CNV_Real)) %>%
    dplyr::mutate(
      CNV_Estimado = tidyr::replace_na(as.integer(CNV_Estimado), 0L),
      Diferencia = as.integer(CNV_Estimado - CNV_Real),
      Concordante = (CNV_Estimado == CNV_Real)
    ) %>%
    dplyr::arrange(Sample, Locus_Ortogrupo)

  n_og <- nrow(pairs_df)
  if (n_og == 0L) {
    return(list(summary = data.frame(), pairs = pairs_df))
  }

  diff_vec <- as.numeric(pairs_df$CNV_Estimado - pairs_df$CNV_Real)
  rmse_val <- sqrt(mean(diff_vec^2))
  mae_val  <- mean(abs(diff_vec))
  bias_val <- mean(diff_vec)
  exact_rate <- mean(pairs_df$Concordante)

  rho_val <- NA_real_
  rho_p   <- NA_real_
  if (n_og >= 3L && stats::sd(pairs_df$CNV_Estimado) > 0 && stats::sd(pairs_df$CNV_Real) > 0) {
    ct <- suppressWarnings(stats::cor.test(pairs_df$CNV_Estimado, pairs_df$CNV_Real, method = "spearman", exact = FALSE))
    rho_val <- as.numeric(ct$estimate)
    rho_p   <- as.numeric(ct$p.value)
  }
  kw_val <- compute_quadratic_weighted_kappa(pairs_df$CNV_Estimado, pairs_df$CNV_Real)

  summary_df <- data.frame(
    Umbral_Cobertura = as.numeric(cov_threshold),
    N_Ortogrupos = as.integer(n_og),
    Copias_Reales = as.integer(sum(pairs_df$CNV_Real)),
    Copias_Estimadas = as.integer(sum(pairs_df$CNV_Estimado)),
    Concordancia_Exacta = exact_rate,
    Spearman_Rho = rho_val,
    Spearman_P_Valor = rho_p,
    Kappa_Cuadratico = kw_val,
    RMSE = rmse_val,
    MAE = mae_val,
    Sesgo_Medio = bias_val,
    stringsAsFactors = FALSE
  )

  list(summary = summary_df, pairs = pairs_df)
}

#' Genera la curva de sensibilidad del umbral `--cov_threshold` entre 30% y 80%
#'
#' Evalúa empíricamente cómo varía el conteo de copias (y, cuando hay verdad de referencia,
#' la concordancia exacta, Kappa cuadrático ponderado, Spearman `rho` y RMSE) al desplazar
#' el umbral de cobertura de homólogos adicionales en el intervalo `[30%, 80%]`.
#'
#' @param df_consenso Data frame con veredictos de consenso (`Veredicto_Consenso`) y `pct_recovered`.
#' @param thresholds Vector numérico de umbrales de cobertura a explorar (por defecto `seq(30, 80, by = 5)`).
#' @param sample_col Columna de muestra.
#' @param locus_col Columna de ortogrupo.
#' @param verdict_col Columna de veredicto del consenso.
#' @param true_cnv_df Data frame opcional con `Sample`, `Locus_Ortogrupo` y `CNV_Real`.
#' @return Data frame con una fila por cada umbral evaluado.
#' @export
evaluate_cnv_threshold_sensitivity <- function(df_consenso,
                                               thresholds = seq(30, 80, by = 5),
                                               sample_col = "Sample",
                                               locus_col = "Locus_Ortogrupo",
                                               verdict_col = "Veredicto_Consenso",
                                               true_cnv_df = NULL) {
  stopifnot(is.data.frame(df_consenso))
  df_work <- df_consenso
  if (!sample_col %in% colnames(df_work)) {
    df_work[[sample_col]] <- if ("sample_name" %in% colnames(df_work)) as.character(df_work$sample_name) else "Sample_1"
  }
  if (!locus_col %in% colnames(df_work)) {
    df_work[[locus_col]] <- if ("locus" %in% colnames(df_work)) as.character(df_work$locus) else as.character(df_work$OG_Limpio)
  }

  has_gt <- !is.null(true_cnv_df) || ("Validacion_Espacial" %in% colnames(df_work) &&
    length(unique(stats::na.omit(df_work$Validacion_Espacial))) >= 1L)

  rows <- lapply(thresholds, function(thr) {
    cnv_r <- extract_cnv_matrix(
      df = df_work,
      sample_col = sample_col,
      locus_col = locus_col,
      cov_threshold = thr,
      verdict_col = verdict_col
    )
    ann <- cnv_r$cnv_annotated
    n_prin <- sum(ann$Clasificacion_Copia == "Copia_Principal", na.rm = TRUE)
    n_hom  <- sum(ann$Clasificacion_Copia == "Homologo_Confirmado", na.rm = TRUE)
    n_noc  <- sum(ann$Clasificacion_Copia == "No_Contado", na.rm = TRUE)

    n_multicopy <- if (nrow(cnv_r$functional_genes) > 0) {
      cnv_r$functional_genes %>%
        dplyr::count(.data[[sample_col]], .data[[locus_col]], name = "k") %>%
        dplyr::filter(k >= 2L) %>%
        nrow()
    } else {
      0L
    }

    base_row <- data.frame(
      Umbral_Cobertura = as.numeric(thr),
      Copias_Totales_Estimadas = as.integer(n_prin + n_hom),
      Copias_Principales = as.integer(n_prin),
      Homologos_Confirmados = as.integer(n_hom),
      Fragmentos_No_Contados = as.integer(n_noc),
      Pares_Multicopia = as.integer(n_multicopy),
      stringsAsFactors = FALSE
    )

    if (has_gt) {
      conc <- evaluate_cnv_concordance(
        df_consenso = df_work,
        cov_threshold = thr,
        sample_col = sample_col,
        locus_col = locus_col,
        verdict_col = verdict_col,
        true_cnv_df = true_cnv_df
      )$summary
      if (nrow(conc) > 0L) {
        base_row <- cbind(
          base_row,
          conc[, c("N_Ortogrupos", "Copias_Reales", "Concordancia_Exacta",
                   "Spearman_Rho", "Spearman_P_Valor", "Kappa_Cuadratico", "RMSE", "MAE", "Sesgo_Medio"), drop = FALSE]
        )
      }
    }
    base_row
  })

  dplyr::bind_rows(rows)
}

#' Evalúa la transferibilidad filogenética de todos los modelos sobre una especie receptora
#'
#' Aplica sin reentrenamiento los cortes univariados de Youden, el modelo heurístico
#' ponderado, los clasificadores supervisados (incluyendo `Ensamble_Stacking`) y el
#' consenso de 7 modelos sobre el conjunto etiquetado de una especie receptora.
#'
#' @param df_sp Data frame preparado de la especie receptora con `Validacion_Espacial`.
#' @param models_list Lista nombrada de modelos supervisados (`models_balanced` o `models_unbalanced`).
#' @param youden_thresholds Vector nombrado de umbrales univariados de Youden.
#' @param heuristic_model Lista con los pesos y umbral global del modelo heurístico.
#' @param predictors Vector de nombres de predictores continuos (por defecto las 4 métricas base).
#' @param ml_voters Vector con los 4 clasificadores supervisados que integran el consenso de 7 modelos.
#' @param min_votes Mínimo de votos afirmativos para el consenso (por defecto `4L`).
#' @return Lista con `df_clean`, `df_consensus`, `benchmark`, `cms`, `scores` y `preds`.
#' @export
evaluate_species_transferability <- function(df_sp,
                                             models_list,
                                             youden_thresholds = NULL,
                                             heuristic_model = NULL,
                                             predictors = c("pct_recovered", "pct_identity", "score", "wscore"),
                                             ml_voters = c("ElasticNet", "NaiveBayes", "SVM_Polinomial", "Red_Neuronal"),
                                             min_votes = 4L) {
  stopifnot(is.data.frame(df_sp))
  pred_cols <- intersect(predictors, colnames(df_sp))
  stopifnot(length(pred_cols) > 0L, "Validacion_Espacial" %in% colnames(df_sp))

  models_list <- normalize_model_names(models_list)

  df_clean <- df_sp %>%
    dplyr::filter(!is.na(Validacion_Espacial)) %>%
    dplyr::mutate(
      Sample = if ("sample_name" %in% colnames(.)) as.character(sample_name) else "Sample_1",
      Locus_Ortogrupo = if ("locus" %in% colnames(.)) {
        dplyr::coalesce(as.character(locus), if ("OG_Limpio" %in% colnames(.)) as.character(OG_Limpio) else NA_character_)
      } else if ("OG_Limpio" %in% colnames(.)) {
        as.character(OG_Limpio)
      } else {
        as.character(ref_name)
      },
      Clase_Real = factor(
        dplyr::if_else(Validacion_Espacial %in% c("TRUE_POSITIVE", "TP"), "TP", "FP"),
        levels = c("TP", "FP")
      ),
      dplyr::across(dplyr::all_of(pred_cols), ~ as.numeric(as.character(.x)))
    ) %>%
    dplyr::filter(stats::complete.cases(dplyr::select(., dplyr::all_of(pred_cols)))) %>%
    as.data.frame()

  x_new <- df_clean[, pred_cols, drop = FALSE]
  y_obs <- df_clean$Clase_Real

  rows   <- list()
  cms    <- list()
  scores <- list()
  preds  <- list()

  if (!is.null(youden_thresholds)) {
    for (y_var in names(youden_thresholds)) {
      if (y_var %in% colnames(x_new)) {
        m_label <- paste0("Youden_", y_var)
        thr_val <- as.numeric(youden_thresholds[[y_var]])
        pred_y  <- factor(ifelse(x_new[[y_var]] >= thr_val, "TP", "FP"), levels = c("TP", "FP"))
        rows[[m_label]] <- compute_classification_metrics(
          pred = pred_y, obs = y_obs, score = x_new[[y_var]], positive_class = "TP", model_name = m_label
        )
        cms[[m_label]]    <- table(Prediccion = pred_y, Referencia = y_obs)
        scores[[m_label]] <- as.numeric(x_new[[y_var]])
        preds[[m_label]]  <- pred_y
      }
    }
  }

  if (!is.null(heuristic_model) && !is.null(youden_thresholds)) {
    m_names  <- names(heuristic_model$weights)
    v_mat    <- sapply(m_names, function(m) as.numeric(x_new[[m]] >= youden_thresholds[[m]]))
    scores_h <- as.vector(v_mat %*% heuristic_model$weights)
    pred_h   <- factor(ifelse(scores_h >= heuristic_model$global_threshold, "TP", "FP"), levels = c("TP", "FP"))
    rows[["Weighted_Heuristic"]] <- compute_classification_metrics(
      pred = pred_h, obs = y_obs, score = scores_h, positive_class = "TP", model_name = "Weighted_Heuristic"
    )
    cms[["Weighted_Heuristic"]]    <- table(Prediccion = pred_h, Referencia = y_obs)
    scores[["Weighted_Heuristic"]] <- scores_h
    preds[["Weighted_Heuristic"]]  <- pred_h
  }

  if (!is.null(models_list) && length(models_list) > 0) {
    for (m_name in names(models_list)) {
      mod_obj <- models_list[[m_name]]
      if (inherits(mod_obj, "caretStack")) {
        tryCatch({
          probs_s <- predict_stacking_probs(mod_obj, newdata = x_new)
          preds_s <- factor(ifelse(probs_s >= 0.5, "TP", "FP"), levels = c("TP", "FP"))
          rows[[m_name]] <- compute_classification_metrics(
            pred = preds_s, obs = y_obs, score = probs_s, positive_class = "TP", model_name = m_name
          )
          cms[[m_name]]    <- table(Prediccion = preds_s, Referencia = y_obs)
          scores[[m_name]] <- probs_s
          preds[[m_name]]  <- preds_s
        }, error = function(e) {})
      } else {
        tryCatch({
          pred_c <- stats::predict(mod_obj, newdata = x_new)
          pred_c <- factor(ifelse(as.character(pred_c) %in% c("TP", "TRUE_POSITIVE"), "TP", "FP"), levels = c("TP", "FP"))
          prob_c <- tryCatch(as.numeric(stats::predict(mod_obj, newdata = x_new, type = "prob")[["TP"]]), error = function(e) NULL)
          rows[[m_name]] <- compute_classification_metrics(
            pred = pred_c, obs = y_obs, score = prob_c, positive_class = "TP", model_name = m_name
          )
          cms[[m_name]]    <- table(Prediccion = pred_c, Referencia = y_obs)
          scores[[m_name]] <- prob_c
          preds[[m_name]]  <- pred_c
        }, error = function(e) {})
      }
    }
  }

  voter_candidates <- c(ml_voters, "Youden_pct_recovered", "Youden_pct_identity", "Weighted_Heuristic")
  avail_voters <- intersect(voter_candidates, names(preds))
  df_consensus <- df_clean

  if (length(avail_voters) >= 5L) {
    for (vn in avail_voters) {
      df_consensus[[paste0("Voto_", vn)]] <- as.character(preds[[vn]])
    }
    vote_cols <- paste0("Voto_", avail_voters)
    eff_min_votes <- min(as.integer(min_votes), length(vote_cols))
    df_consensus <- consolidar_consenso(df_consensus, model_cols = vote_cols, min_votes = eff_min_votes)
    pred_cons <- factor(
      ifelse(df_consensus$Veredicto_Consenso == "Homologo_Valido", "TP", "FP"),
      levels = c("TP", "FP")
    )
    rows[["Consenso_7_Modelos"]] <- compute_classification_metrics(
      pred = pred_cons,
      obs = y_obs,
      score = df_consensus$Votos_TP / length(avail_voters),
      positive_class = "TP",
      model_name = "Consenso_7_Modelos"
    )
    cms[["Consenso_7_Modelos"]]    <- table(Prediccion = pred_cons, Referencia = y_obs)
    scores[["Consenso_7_Modelos"]] <- as.numeric(df_consensus$Votos_TP)
    preds[["Consenso_7_Modelos"]]  <- pred_cons
  }

  df_bench <- dplyr::bind_rows(rows) %>%
    dplyr::filter(!is.na(Modelo) & nzchar(as.character(Modelo))) %>%
    dplyr::distinct(Modelo, .keep_all = TRUE) %>%
    dplyr::arrange(dplyr::desc(ROC_AUC), dplyr::desc(F1_Score))

  cms <- cms[intersect(df_bench$Modelo, names(cms))]

  list(
    df_clean = df_clean,
    df_consensus = df_consensus,
    benchmark = df_bench,
    cms = cms,
    scores = scores,
    preds = preds
  )
}
