# ==============================================================================
# TargetCNV-ML: Especificaciones de Modelos y Adaptadores (R/model_specs.R)
# ==============================================================================

#' Normaliza los nombres de modelos para mantener nomenclatura uniforme
#'
#' @param mod_list Lista nombrada de modelos o vector de nombres.
#' @return Objeto con nombres estandarizados (`MLP`, `MLP_Regularizado`, `Stacked Autoencoder - DNN`).
#' @export
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

#' Adaptador compatible para XGBoost >= 1.7+ con caret (evita el error por deprecacion de ntreelimit)
#'
#' @return Lista con la definicion del metodo `xgbTree` compatible con `caret::train`.
#' @export
get_safe_xgb_method <- function() {
  xgb_mod <- caret::getModelInfo("xgbTree", regex = FALSE)[[1]]
  xgb_mod$fit <- function(x, y, wts, param, lev, last, classProbs, ...) {
    if (!inherits(x, "xgb.DMatrix")) {
      x <- as.matrix(x)
    }
    if (is.factor(y)) {
      if (length(lev) == 2) {
        y_num <- ifelse(y == lev[1], 1, 0)
        dtrain <- xgboost::xgb.DMatrix(data = x, label = y_num, missing = NA)
        if (!is.null(wts)) xgboost::setinfo(dtrain, "weight", wts)
        out <- xgboost::xgb.train(
          params = list(
            eta = param$eta,
            max_depth = param$max_depth,
            gamma = param$gamma,
            colsample_bytree = param$colsample_bytree,
            min_child_weight = param$min_child_weight,
            subsample = param$subsample,
            objective = "binary:logistic",
            eval_metric = "logloss"
          ),
          data = dtrain,
          nrounds = param$nrounds,
          verbose = 0,
          ...
        )
      } else {
        y_num <- as.numeric(y) - 1
        dtrain <- xgboost::xgb.DMatrix(data = x, label = y_num, missing = NA)
        if (!is.null(wts)) xgboost::setinfo(dtrain, "weight", wts)
        out <- xgboost::xgb.train(
          params = list(
            eta = param$eta,
            max_depth = param$max_depth,
            gamma = param$gamma,
            colsample_bytree = param$colsample_bytree,
            min_child_weight = param$min_child_weight,
            subsample = param$subsample,
            num_class = length(lev),
            objective = "multi:softprob"
          ),
          data = dtrain,
          nrounds = param$nrounds,
          verbose = 0,
          ...
        )
      }
    }
    out
  }
  xgb_mod$predict <- function(modelFit, newdata, submodels = NULL) {
    if (!inherits(newdata, "xgb.DMatrix")) {
      newdata <- xgboost::xgb.DMatrix(as.matrix(newdata), missing = NA)
    }
    out <- stats::predict(modelFit, newdata)
    if (length(modelFit$obsLevels) == 2) {
      res <- ifelse(out >= 0.5, modelFit$obsLevels[1], modelFit$obsLevels[2])
    } else {
      out_m <- matrix(out, ncol = length(modelFit$obsLevels), byrow = TRUE)
      res <- modelFit$obsLevels[apply(out_m, 1, which.max)]
    }
    if (!is.null(submodels)) {
      tmp <- vector(mode = "list", length = nrow(submodels) + 1)
      tmp[[1]] <- res
      for (j in seq(along = submodels$nrounds)) {
        n_it <- submodels$nrounds[j]
        p_j <- tryCatch(
          stats::predict(modelFit, newdata, iterationrange = c(1, n_it)),
          error = function(e) out
        )
        if (length(modelFit$obsLevels) == 2) {
          tmp[[j + 1]] <- ifelse(p_j >= 0.5, modelFit$obsLevels[1], modelFit$obsLevels[2])
        } else {
          p_m <- matrix(p_j, ncol = length(modelFit$obsLevels), byrow = TRUE)
          tmp[[j + 1]] <- modelFit$obsLevels[apply(p_m, 1, which.max)]
        }
      }
      res <- tmp
    }
    res
  }
  xgb_mod$prob <- function(modelFit, newdata, submodels = NULL) {
    if (!inherits(newdata, "xgb.DMatrix")) {
      newdata <- xgboost::xgb.DMatrix(as.matrix(newdata), missing = NA)
    }
    out <- stats::predict(modelFit, newdata)
    if (length(modelFit$obsLevels) == 2) {
      out_df <- as.data.frame(cbind(out, 1 - out), stringsAsFactors = FALSE)
      colnames(out_df) <- modelFit$obsLevels
    } else {
      out_m <- matrix(out, ncol = length(modelFit$obsLevels), byrow = TRUE)
      colnames(out_m) <- modelFit$obsLevels
      out_df <- as.data.frame(out_m, stringsAsFactors = FALSE)
    }
    if (!is.null(submodels)) {
      tmp <- vector(mode = "list", length = nrow(submodels) + 1)
      tmp[[1]] <- out_df
      for (j in seq(along = submodels$nrounds)) {
        n_it <- submodels$nrounds[j]
        p_j <- tryCatch(
          stats::predict(modelFit, newdata, iterationrange = c(1, n_it)),
          error = function(e) out
        )
        if (length(modelFit$obsLevels) == 2) {
          df_j <- as.data.frame(cbind(p_j, 1 - p_j), stringsAsFactors = FALSE)
          colnames(df_j) <- modelFit$obsLevels
        } else {
          p_m <- matrix(p_j, ncol = length(modelFit$obsLevels), byrow = TRUE)
          colnames(p_m) <- modelFit$obsLevels
          df_j <- as.data.frame(p_m, stringsAsFactors = FALSE)
        }
        tmp[[j + 1]] <- df_j
      }
      out_df <- tmp
    }
    out_df
  }
  xgb_mod
}

#' Construye el catalogo de especificaciones de modelos supervisados para caret
#'
#' @param n_obs Numero de observaciones de entrenamiento (para ajustar `n.minobsinnode` en GBM).
#' @param extended Booleano; TRUE devuelve los 18 clasificadores, FALSE devuelve los 5 rapidos.
#' @return Lista nombrada con las configuraciones de entrenamiento.
#' @export
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

# Declaracion explicita de dependencias para escaneo estatico de renv::dependencies() y renv::hydrate()
# Incluye motores dinamicos de caret y paquetes de renderizado Quarto/RMarkdown
if (FALSE) {
  library(rmarkdown)
  library(knitr)
  library(DT)
  library(plotly)
  library(cowplot)
  library(jsonlite)
  library(htmltools)
  library(glmnet)
  library(klaR)
  library(kernlab)
  library(randomForest)
  library(gbm)
  library(xgboost)
  library(RSNNS)
  library(deepnet)
  library(ipred)
  library(e1071)
  library(MASS)
  library(nnet)
  library(rpart)
  library(rpart.plot)
  library(NeuralNetTools)
  library(caretEnsemble)
  library(kernelshap)
  library(shapviz)
}
