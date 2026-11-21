test_that("XGBoost, Ensamble_Stacking (ROC_AUC > 0.5) y semillas deterministas funcionan de extremo a extremo", {
  set.seed(1234)
  n_tp <- 45L
  n_fp <- 55L
  df_synth <- data.frame(
    Clase = factor(c(rep("TP", n_tp), rep("FP", n_fp)), levels = c("TP", "FP")),
    pct_recovered = c(rnorm(n_tp, mean = 88, sd = 7), rnorm(n_fp, mean = 36, sd = 12)),
    pct_identity  = c(rnorm(n_tp, mean = 93, sd = 3), rnorm(n_fp, mean = 68, sd = 9)),
    score         = c(rnorm(n_tp, mean = 0.89, sd = 0.05), rnorm(n_fp, mean = 0.39, sd = 0.11)),
    wscore        = c(rnorm(n_tp, mean = 0.84, sd = 0.06), rnorm(n_fp, mean = 0.29, sd = 0.10))
  )
  pred_cols <- c("pct_recovered", "pct_identity", "score", "wscore")
  x_tr <- df_synth[, pred_cols, drop = FALSE]
  y_tr <- df_synth$Clase

  # 1. Test de colinealidad de Spearman
  sp_res <- compute_spearman_collinearity(df_synth, metric_cols = pred_cols)
  expect_equal(dim(sp_res$cor_matrix), c(4L, 4L))
  expect_equal(as.numeric(diag(sp_res$cor_matrix)), rep(1.0, 4L))
  expect_equal(nrow(sp_res$pairs), 6L)

  # 2. Adaptador XGBoost compatible con xgboost >= 3.0 y reproducibilidad secuencial vs paralelo
  xgb_method <- get_safe_xgb_method()
  ctrl_seq <- caret::trainControl(
    method = "cv", number = 4, classProbs = TRUE,
    summaryFunction = caret::twoClassSummary, allowParallel = FALSE,
    savePredictions = "final", sampling = "down"
  )
  det_seq <- build_deterministic_cv_control(ctrl_seq, x_tr, y_tr, xgb_method, tuneLength = 2L, seed = 1234L)
  fit_xgb_seq <- suppressWarnings(caret::train(
    x = x_tr, y = y_tr, method = xgb_method,
    trControl = det_seq$trControl, tuneGrid = det_seq$tuneGrid,
    nthread = 1L, metric = "ROC"
  ))

  probs_xgb_seq <- stats::predict(fit_xgb_seq, newdata = x_tr, type = "prob")[["TP"]]
  expect_true(all(probs_xgb_seq >= 0 & probs_xgb_seq <= 1))
  xgb_metrics <- compute_classification_metrics(
    pred = fit_xgb_seq$pred$pred,
    obs = fit_xgb_seq$pred$obs,
    score = fit_xgb_seq$pred$TP,
    model_name = "XGBoost"
  )
  expect_gt(xgb_metrics$ROC_AUC, 0.75)

  imp_xgb <- caret::varImp(fit_xgb_seq)
  expect_s3_class(imp_xgb, "varImp.train")
  expect_equal(nrow(imp_xgb$importance), 4L)

  if (requireNamespace("doParallel", quietly = TRUE)) {
    cl <- parallel::makeCluster(2L)
    on.exit({
      try(parallel::stopCluster(cl), silent = TRUE)
      foreach::registerDoSEQ()
    }, add = TRUE)
    doParallel::registerDoParallel(cl)

    ctrl_par <- caret::trainControl(
      method = "cv", number = 4, classProbs = TRUE,
      summaryFunction = caret::twoClassSummary, allowParallel = TRUE,
      savePredictions = "final", sampling = "down"
    )
    det_par <- build_deterministic_cv_control(ctrl_par, x_tr, y_tr, xgb_method, tuneLength = 2L, seed = 1234L)
    fit_xgb_par <- suppressWarnings(caret::train(
      x = x_tr, y = y_tr, method = xgb_method,
      trControl = det_par$trControl, tuneGrid = det_par$tuneGrid,
      nthread = 1L, metric = "ROC"
    ))
    probs_xgb_par <- stats::predict(fit_xgb_par, newdata = x_tr, type = "prob")[["TP"]]
    expect_equal(probs_xgb_seq, probs_xgb_par, tolerance = 1e-12)
  }

  # 3. Invariante ROC_AUC > 0.5 en Ensamble_Stacking (tanto CV fuera de muestra como predict_stacking_probs)
  skip_if_not_installed("caretEnsemble")
  det_glmnet <- build_deterministic_cv_control(
    ctrl_seq, x_tr, y_tr, "glmnet", tuneLength = 3L, preProcess = c("center", "scale"), seed = 1234L
  )
  det_lda <- build_deterministic_cv_control(
    ctrl_seq, x_tr, y_tr, "lda", preProcess = c("center", "scale"), seed = 1234L
  )

  m_glmnet <- suppressWarnings(caret::train(
    x = x_tr, y = y_tr, method = "glmnet",
    trControl = det_glmnet$trControl, tuneGrid = det_glmnet$tuneGrid,
    preProcess = c("center", "scale"), metric = "ROC"
  ))
  m_lda <- suppressWarnings(caret::train(
    x = x_tr, y = y_tr, method = "lda",
    trControl = det_lda$trControl, tuneGrid = det_lda$tuneGrid,
    preProcess = c("center", "scale"), metric = "ROC"
  ))

  lista_base <- list(ElasticNet = m_glmnet, LDA = m_lda)
  class(lista_base) <- c("caretList", "list")

  stack_ctrl <- caret::trainControl(
    method = "cv", number = 4,
    index = det_glmnet$trControl$index,
    indexOut = det_glmnet$trControl$indexOut,
    savePredictions = "final", classProbs = TRUE,
    summaryFunction = caret::twoClassSummary
  )
  det_stack <- build_deterministic_cv_control(
    stack_ctrl, data.frame(ElasticNet = rep(0.5, nrow(x_tr)), LDA = rep(0.5, nrow(x_tr))),
    y_tr, method = "glm", tuneLength = 1L, seed = 1234L
  )
  ens_stack <- suppressWarnings(caretEnsemble::caretStack(
    lista_base, method = "glm", metric = "ROC", trControl = det_stack$trControl
  ))

  oof_metrics <- compute_classification_metrics(
    pred = ens_stack$ens_model$pred$pred,
    obs = ens_stack$ens_model$pred$obs,
    score = ens_stack$ens_model$pred$TP,
    model_name = "Ensamble_Stacking_OOF"
  )
  expect_gt(oof_metrics$ROC_AUC, 0.5)
  expect_gt(oof_metrics$ROC_AUC, 0.85)

  probs_new <- predict_stacking_probs(ens_stack, newdata = x_tr)
  preds_new <- factor(ifelse(probs_new >= 0.5, "TP", "FP"), levels = c("TP", "FP"))
  new_metrics <- compute_classification_metrics(
    pred = preds_new, obs = y_tr, score = probs_new, model_name = "Ensamble_Stacking_New"
  )
  expect_gt(new_metrics$ROC_AUC, 0.5)
  expect_gt(new_metrics$ROC_AUC, 0.85)
})

test_that("plot_model_varimp calcula importancias específicas por modelo y evaluate_species_transferability soporta 1 fila", {
  set.seed(1234)
  n_tp <- 40L
  n_fp <- 40L
  df_synth <- data.frame(
    Validacion_Espacial = c(rep("TRUE_POSITIVE", n_tp), rep("FALSE_POSITIVE", n_fp)),
    Clase = factor(c(rep("TP", n_tp), rep("FP", n_fp)), levels = c("TP", "FP")),
    pct_recovered = c(rnorm(n_tp, 90, 5), rnorm(n_fp, 35, 10)),
    pct_identity  = c(rnorm(n_tp, 94, 2), rnorm(n_fp, 70, 8)),
    score         = c(rnorm(n_tp, 0.90, 0.04), rnorm(n_fp, 0.35, 0.10)),
    wscore        = c(rnorm(n_tp, 0.85, 0.05), rnorm(n_fp, 0.25, 0.09)),
    sample_name   = "S1",
    locus         = sprintf("OG%04d", seq_len(n_tp + n_fp)),
    stringsAsFactors = FALSE
  )
  pred_cols <- c("pct_recovered", "pct_identity", "score", "wscore")
  x_tr <- df_synth[, pred_cols, drop = FALSE]
  y_tr <- df_synth$Clase

  ctrl <- caret::trainControl(
    method = "cv", number = 3, classProbs = TRUE,
    summaryFunction = caret::twoClassSummary, allowParallel = FALSE,
    savePredictions = "final"
  )
  det_lda <- build_deterministic_cv_control(ctrl, x_tr, y_tr, "lda", preProcess = c("center", "scale"), seed = 1234L)
  det_nb  <- build_deterministic_cv_control(ctrl, x_tr, y_tr, "nb",  preProcess = c("center", "scale"), seed = 1234L)

  m_lda <- suppressWarnings(caret::train(x = x_tr, y = y_tr, method = "lda", trControl = det_lda$trControl, preProcess = c("center", "scale"), metric = "ROC"))
  m_nb  <- suppressWarnings(caret::train(x = x_tr, y = y_tr, method = "nb",  trControl = det_nb$trControl, tuneGrid = det_nb$tuneGrid, preProcess = c("center", "scale"), metric = "ROC"))

  p_lda <- plot_model_varimp(m_lda, df_data = df_synth, model_name = "LDA")
  p_nb  <- plot_model_varimp(m_nb,  df_data = df_synth, model_name = "NaiveBayes")
  p_num <- plot_model_varimp(c(pct_recovered = 3, pct_identity = 2, score = 1, wscore = 1), model_name = "Heuristic")

  expect_s3_class(p_lda, "ggplot")
  expect_s3_class(p_nb,  "ggplot")
  expect_equal(sum(p_lda$data$Importancia), 100, tolerance = 1e-6)
  expect_equal(sum(p_nb$data$Importancia),  100, tolerance = 1e-6)
  expect_equal(sum(p_num$data$Importancia), 100, tolerance = 1e-6)
  # Ningún peso positivo debe colapsarse a 0% por escalado min-max
  expect_true(all(p_num$data$Importancia > 0))
  # LDA (coeficientes LD1) y NaiveBayes (permutación sobre P(TP)) no deben producir barras idénticas
  v_lda <- p_lda$data$Importancia[match(pred_cols, as.character(p_lda$data$Variable))]
  v_nb  <- p_nb$data$Importancia[match(pred_cols, as.character(p_nb$data$Variable))]
  expect_false(isTRUE(all.equal(v_lda, v_nb, tolerance = 1e-3)))

  # Caso borde: evaluación de transferencia sobre un data frame de 1 sola fila
  y_eval <- evaluate_youden_predictors(df_synth, "Clase", pred_cols, positive_class = "TP")
  h_res  <- optimize_weighted_heuristic(y_eval$binary_votes, df_synth$Clase, pred_cols, positive_class = "TP")
  h_obj  <- list(weights = h_res$weights, global_threshold = h_res$global_threshold)
  res_1row <- evaluate_species_transferability(
    df_sp = df_synth[1, , drop = FALSE],
    models_list = list(LDA = m_lda),
    youden_thresholds = y_eval$thresholds,
    heuristic_model = h_obj,
    predictors = pred_cols
  )
  expect_equal(nrow(res_1row$df_clean), 1L)
  expect_true("Weighted_Heuristic" %in% res_1row$benchmark$Modelo)
  expect_true("LDA" %in% res_1row$benchmark$Modelo)
})


