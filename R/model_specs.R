# Catálogo de modelos supervisados, adaptador xgbTree para xgboost >= 3.0 y semillas CV

# Estandariza alias antiguos de redes neuronales para mantener consistencia en tablas y reportes.
normalize_model_names <- function(mod_list) {
  if (is.null(mod_list) || length(mod_list) == 0) return(mod_list)
  if (is.character(mod_list)) {
    mod_list[mod_list == "DNN_Profunda"]     <- "MLP"
    mod_list[mod_list == "DNN_Regularizada"] <- "MLP_Regularizado"
    mod_list[mod_list %in% c("DNN_Autoencoder", "Autoencoder", "SAE-DNN", "SAE_DNN")] <- "Stacked Autoencoder - DNN"
    return(mod_list)
  }
  nms <- names(mod_list)
  if (!is.null(nms)) {
    nms[nms == "DNN_Profunda"]     <- "MLP"
    nms[nms == "DNN_Regularizada"] <- "MLP_Regularizado"
    nms[nms %in% c("DNN_Autoencoder", "Autoencoder", "SAE-DNN", "SAE_DNN")] <- "Stacked Autoencoder - DNN"
    names(mod_list) <- nms
    mod_list <- mod_list[!duplicated(names(mod_list))]
  }
  mod_list
}

# Envoltorio de xgbTree para xgboost >= 3.0: encapsula el booster ALTREP en una lista
# para que caret pueda adjuntar xNames/obsLevels y usa iterationrange en base 1.
get_safe_xgb_method <- function() {
  xgb_mod <- caret::getModelInfo("xgbTree", regex = FALSE)[[1]]

  unwrap_booster <- function(model_fit) {
    if (is.list(model_fit) && !is.null(model_fit$booster)) model_fit$booster else model_fit
  }

  xgb_mod$fit <- function(x, y, wts, param, lev, last, classProbs, ...) {
    dots <- list(...)
    nthread_val <- if ("nthread" %in% names(dots)) as.integer(dots$nthread) else 1L
    dots$nthread <- NULL

    if (!inherits(x, "xgb.DMatrix")) {
      x <- as.matrix(x)
    }

    if (length(lev) == 2L) {
      y_num <- ifelse(y == lev[1], 1, 0)
      dtrain <- xgboost::xgb.DMatrix(data = x, label = y_num, missing = NA)
      if (!is.null(wts)) xgboost::setinfo(dtrain, "weight", wts)
      xgb_params <- list(
        eta = param$eta,
        max_depth = as.integer(param$max_depth),
        gamma = param$gamma,
        colsample_bytree = param$colsample_bytree,
        min_child_weight = param$min_child_weight,
        subsample = param$subsample,
        nthread = nthread_val,
        objective = "binary:logistic",
        eval_metric = "logloss"
      )
    } else {
      y_num <- as.numeric(y) - 1
      dtrain <- xgboost::xgb.DMatrix(data = x, label = y_num, missing = NA)
      if (!is.null(wts)) xgboost::setinfo(dtrain, "weight", wts)
      xgb_params <- list(
        eta = param$eta,
        max_depth = as.integer(param$max_depth),
        gamma = param$gamma,
        colsample_bytree = param$colsample_bytree,
        min_child_weight = param$min_child_weight,
        subsample = param$subsample,
        nthread = nthread_val,
        num_class = length(lev),
        objective = "multi:softprob"
      )
    }

    bst <- do.call(
      xgboost::xgb.train,
      c(list(params = xgb_params, data = dtrain, nrounds = as.integer(param$nrounds), verbose = 0), dots)
    )
    list(booster = bst)
  }

  xgb_mod$predict <- function(modelFit, newdata, submodels = NULL) {
    bst <- unwrap_booster(modelFit)
    if (!inherits(newdata, "xgb.DMatrix")) {
      newdata <- xgboost::xgb.DMatrix(as.matrix(newdata), missing = NA)
    }
    out <- as.numeric(stats::predict(bst, newdata))
    obs_lev <- modelFit$obsLevels

    if (length(obs_lev) == 2L) {
      res <- ifelse(out >= 0.5, obs_lev[1], obs_lev[2])
    } else {
      out_m <- matrix(out, ncol = length(obs_lev), byrow = TRUE)
      res <- obs_lev[apply(out_m, 1, which.max)]
    }

    if (!is.null(submodels) && nrow(submodels) > 0L) {
      tmp <- vector(mode = "list", length = nrow(submodels) + 1L)
      tmp[[1]] <- res
      for (j in seq_along(submodels$nrounds)) {
        n_it <- max(1L, as.integer(submodels$nrounds[j]))
        p_j <- tryCatch(
          as.numeric(stats::predict(bst, newdata, iterationrange = c(1L, n_it))),
          error = function(e) out
        )
        if (length(obs_lev) == 2L) {
          tmp[[j + 1L]] <- ifelse(p_j >= 0.5, obs_lev[1], obs_lev[2])
        } else {
          p_m <- matrix(p_j, ncol = length(obs_lev), byrow = TRUE)
          tmp[[j + 1L]] <- obs_lev[apply(p_m, 1, which.max)]
        }
      }
      res <- tmp
    }
    res
  }

  xgb_mod$prob <- function(modelFit, newdata, submodels = NULL) {
    bst <- unwrap_booster(modelFit)
    if (!inherits(newdata, "xgb.DMatrix")) {
      newdata <- xgboost::xgb.DMatrix(as.matrix(newdata), missing = NA)
    }
    out <- as.numeric(stats::predict(bst, newdata))
    obs_lev <- modelFit$obsLevels

    if (length(obs_lev) == 2L) {
      out_df <- as.data.frame(cbind(out, 1 - out), stringsAsFactors = FALSE)
      colnames(out_df) <- obs_lev
    } else {
      out_m <- matrix(out, ncol = length(obs_lev), byrow = TRUE)
      colnames(out_m) <- obs_lev
      out_df <- as.data.frame(out_m, stringsAsFactors = FALSE)
    }

    if (!is.null(submodels) && nrow(submodels) > 0L) {
      tmp <- vector(mode = "list", length = nrow(submodels) + 1L)
      tmp[[1]] <- out_df
      for (j in seq_along(submodels$nrounds)) {
        n_it <- max(1L, as.integer(submodels$nrounds[j]))
        p_j <- tryCatch(
          as.numeric(stats::predict(bst, newdata, iterationrange = c(1L, n_it))),
          error = function(e) out
        )
        if (length(obs_lev) == 2L) {
          df_j <- as.data.frame(cbind(p_j, 1 - p_j), stringsAsFactors = FALSE)
          colnames(df_j) <- obs_lev
        } else {
          p_m <- matrix(p_j, ncol = length(obs_lev), byrow = TRUE)
          colnames(p_m) <- obs_lev
          df_j <- as.data.frame(p_m, stringsAsFactors = FALSE)
        }
        tmp[[j + 1L]] <- df_j
      }
      out_df <- tmp
    }
    out_df
  }

  xgb_mod$varImp <- function(object, numTrees = NULL, ...) {
    bst <- unwrap_booster(object)
    x_names <- object$xNames
    imp <- tryCatch(
      xgboost::xgb.importance(feature_names = x_names, model = bst),
      error = function(e) NULL
    )
    if (is.null(imp) || nrow(imp) == 0L) {
      if (is.null(x_names)) return(data.frame(Overall = numeric(0)))
      out_imp <- data.frame(Overall = rep(0, length(x_names)), row.names = x_names)
      return(out_imp)
    }
    imp_df <- as.data.frame(imp)
    all_vars <- if (!is.null(x_names)) x_names else as.character(imp_df$Feature)
    scores <- stats::setNames(rep(0, length(all_vars)), all_vars)
    matched <- intersect(as.character(imp_df$Feature), all_vars)
    scores[matched] <- as.numeric(imp_df$Gain[match(matched, as.character(imp_df$Feature))]) * 100
    data.frame(Overall = as.numeric(scores), row.names = names(scores))
  }

  xgb_mod
}

# Construye particiones (index/indexOut) y semillas explícitas para caret::trainControl
# de modo que los resultados bajo una misma semilla no cambien al variar el número de núcleos.
build_deterministic_cv_control <- function(base_ctrl,
                                           x,
                                           y,
                                           method,
                                           tuneGrid = NULL,
                                           tuneLength = 3L,
                                           preProcess = NULL,
                                           seed = 1234L) {
  x_df <- as.data.frame(x)
  models <- if (is.list(method)) method else caret::getModelInfo(method, regex = FALSE)[[1]]

  set.seed(as.integer(seed))
  rs_seed <- sample.int(.Machine$integer.max, 1L)

  ctrl <- base_ctrl
  if (is.null(ctrl$index)) {
    make_rs_fn <- tryCatch(utils::getFromNamespace("make_resamples", "caret"), error = function(e) NULL)
    if (is.function(make_rs_fn)) {
      ctrl <- withr::with_seed(rs_seed, make_rs_fn(ctrl, outcome = y))
    } else {
      folds <- withr::with_seed(
        rs_seed,
        caret::createFolds(y, k = ctrl$number, list = TRUE, returnTrain = TRUE)
      )
      ctrl$index <- folds
      ctrl$indexOut <- lapply(folds, function(idx) setdiff(seq_along(y), idx))
    }
  }

  grid_df <- tuneGrid
  if (is.null(grid_df)) {
    t_len <- if (!is.null(tuneLength)) as.integer(tuneLength) else if (ctrl$method == "none") 1L else 3L
    if (!is.null(preProcess) &&
        length(models$parameters$parameter) > 1L &&
        all(as.character(models$parameters$parameter) != "parameter")) {
      pp_obj <- caret::preProcess(x_df, method = preProcess)
      grid_df <- models$grid(x = stats::predict(pp_obj, x_df), y = y, len = t_len, search = ctrl$search)
    } else {
      grid_df <- models$grid(x = x_df, y = y, len = t_len, search = ctrl$search)
    }
  }
  grid_df <- as.data.frame(grid_df)
  colnames(grid_df) <- sub("^\\.", "", colnames(grid_df))
  grid_df <- grid_df[!duplicated(grid_df), , drop = FALSE]

  train_info <- if (is.function(models$loop) && nrow(grid_df) > 1L) {
    models$loop(grid_df)
  } else {
    list(loop = grid_df)
  }

  num_rs <- length(ctrl$index)
  if (ctrl$method %in% c("boot632", "optimism_boot", "boot_all")) {
    num_rs <- num_rs + 1L
  }
  n_loop <- max(1L, nrow(train_info$loop))

  raw_seeds <- sample.int(n = 1000000L, size = num_rs * n_loop + 1L)
  seeds_list <- lapply(
    seq(from = 1L, to = num_rs * n_loop, by = n_loop),
    function(idx) as.integer(raw_seeds[idx:(idx + n_loop - 1L)])
  )
  seeds_list[[num_rs + 1L]] <- as.integer(raw_seeds[num_rs * n_loop + 1L])
  ctrl$seeds <- seeds_list

  list(trControl = ctrl, tuneGrid = grid_df)
}

# Especificaciones de los 18 clasificadores (o los 5 modelos base en modo rápido).
get_model_specs <- function(n_obs = 184, extended = TRUE) {
  gbm_grid <- expand.grid(
    interaction.depth = c(1, 3, 5),
    n.trees = c(50, 100, 150),
    shrinkage = c(0.05, 0.1),
    n.minobsinnode = max(2, min(10, floor(n_obs / 15)))
  )
  dnn_grid <- expand.grid(layer1 = c(10, 20), layer2 = c(5, 10), layer3 = c(0, 5))
  dnn_auto_grid <- expand.grid(
    layer1 = c(16, 32),
    layer2 = c(8, 16),
    layer3 = c(4, 8),
    hidden_dropout = c(0, 0.1),
    visible_dropout = c(0)
  )

  if (extended) {
    list(
      RandomForest     = list(method = "rf",               tuneLength = 3, file = "model_randomforest.rds",     no_scale = TRUE),
      GBM              = list(method = "gbm",              tuneGrid = gbm_grid, verbose = FALSE, file = "model_gbm.rds", no_scale = TRUE),
      XGBoost          = list(method = get_safe_xgb_method(), tuneLength = 3, nthread = 1, file = "model_xgboost.rds", no_scale = TRUE),
      SVM_Radial       = list(method = "svmRadial",        tuneLength = 3, file = "model_svm_radial.rds"),
      SVM_Lineal       = list(method = "svmLinear",        tuneLength = 3, file = "model_svm_lineal.rds"),
      SVM_Polinomial   = list(method = "svmPoly",          tuneLength = 3, file = "model_svm_polinomial.rds"),
      Red_Neuronal     = list(method = "nnet",             tuneLength = 3, trace = FALSE, MaxNWts = 5000, file = "model_red_neuronal.rds"),
      Red_Neuronal_Bag = list(method = "avNNet",           tuneLength = 3, trace = FALSE, MaxNWts = 5000, file = "model_red_neuronal_bag.rds"),
      ElasticNet       = list(method = "glmnet",           tuneLength = 5, file = "model_elasticnet.rds"),
      KNN              = list(method = "knn",              tuneLength = 3, file = "model_knn.rds"),
      CART             = list(method = "rpart",            tuneLength = 3, file = "model_cart.rds",             no_scale = TRUE),
      NaiveBayes       = list(method = "nb",               tuneLength = 2, file = "model_naivebayes.rds"),
      LDA              = list(method = "lda",              file = "model_lda.rds"),
      MLP              = list(method = "mlpML",            tuneGrid = dnn_grid, file = "model_mlp.rds"),
      MLP_Regularizado = list(method = "mlpWeightDecayML", tuneLength = 3, file = "model_mlp_regularizado.rds"),
      QDA              = list(method = "qda",              file = "model_qda.rds"),
      Bagging_Arboles  = list(method = "treebag",          file = "model_bagging_arboles.rds",                  no_scale = TRUE),
      `Stacked Autoencoder - DNN` = list(method = "dnn",   tuneGrid = dnn_auto_grid, file = "model_stacked_autoencoder_dnn.rds")
    )
  } else {
    list(
      ElasticNet     = list(method = "glmnet",  tuneLength = 5, file = "model_elasticnet.rds"),
      NaiveBayes     = list(method = "nb",      tuneLength = 2, file = "model_naivebayes.rds"),
      SVM_Polinomial = list(method = "svmPoly", tuneLength = 3, file = "model_svm_polinomial.rds"),
      Red_Neuronal   = list(method = "nnet",    tuneLength = 3, trace = FALSE, MaxNWts = 5000, file = "model_red_neuronal.rds"),
      LDA            = list(method = "lda",     file = "model_lda.rds")
    )
  }
}

# Referencias para que renv detecte los motores que caret carga dinámicamente
.register_caret_engines <- function() {
  pkgs <- c(
    "rmarkdown", "knitr", "yaml", "DT", "plotly", "cowplot", "colorspace",
    "jsonlite", "htmltools", "glmnet", "klaR", "kernlab", "randomForest",
    "gbm", "xgboost", "RSNNS", "deepnet", "ipred", "e1071", "MASS", "nnet",
    "rpart", "rpart.plot", "NeuralNetTools", "caretEnsemble", "kernelshap",
    "shapviz", "withr"
  )
  invisible(vapply(pkgs, requireNamespace, logical(1), quietly = TRUE))
}
