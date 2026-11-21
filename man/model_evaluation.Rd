\name{model_evaluation}
\alias{model_evaluation}
\alias{compute_classification_metrics}
\alias{calculate_youden_threshold}
\alias{evaluate_youden_predictors}
\alias{optimize_weighted_heuristic}
\alias{predict_weighted_heuristic}
\alias{predict_stacking_probs}
\alias{consolidar_consenso}
\alias{compute_spearman_collinearity}
\alias{compute_quadratic_weighted_kappa}
\alias{extract_reference_cnv}
\alias{evaluate_cnv_concordance}
\alias{evaluate_cnv_threshold_sensitivity}
\alias{evaluate_species_transferability}
\title{Classification metrics, collinearity, consensus voting, and CNV concordance}
\description{
Statistical evaluation utilities for binary classification benchmarks, Youden
threshold optimization, weighted heuristic voting, Spearman collinearity testing,
multi-model consensus aggregation, orthogroup-level copy-number concordance,
coverage threshold sensitivity sweeps, and cross-species transferability.
}
\usage{
compute_classification_metrics(
  pred,
  obs,
  positive_class = "TRUE_POSITIVE",
  model_name = "Modelo",
  score = NULL
)

calculate_youden_threshold(
  response,
  predictor,
  positive_class = "TRUE_POSITIVE"
)

evaluate_youden_predictors(
  df,
  response_col,
  metric_cols,
  positive_class = "TRUE_POSITIVE"
)

optimize_weighted_heuristic(
  df_votes,
  response,
  metric_cols = colnames(df_votes),
  weight_range = 0:3,
  positive_class = "TRUE_POSITIVE"
)

predict_weighted_heuristic(df_predictors, thresholds, weights, global_threshold)

predict_stacking_probs(ens_stack, newdata)

consolidar_consenso(df_predictions, model_cols, min_votes = 4)

compute_spearman_collinearity(
  df,
  metric_cols = c("pct_recovered", "pct_identity", "score", "wscore")
)

compute_quadratic_weighted_kappa(pred_copies, real_copies)

extract_reference_cnv(
  df,
  sample_col = "Sample",
  locus_col = "Locus_Ortogrupo",
  class_col = "Validacion_Espacial"
)

evaluate_cnv_concordance(
  df_consenso,
  cov_threshold = 50,
  sample_col = "Sample",
  locus_col = "Locus_Ortogrupo",
  verdict_col = "Veredicto_Consenso",
  true_cnv_df = NULL
)

evaluate_cnv_threshold_sensitivity(
  df_consenso,
  thresholds = seq(30, 80, by = 5),
  sample_col = "Sample",
  locus_col = "Locus_Ortogrupo",
  verdict_col = "Veredicto_Consenso",
  true_cnv_df = NULL
)

evaluate_species_transferability(
  df_sp,
  models_list,
  youden_thresholds = NULL,
  heuristic_model = NULL,
  predictors = c("pct_recovered", "pct_identity", "score", "wscore"),
  ml_voters = c("ElasticNet", "NaiveBayes", "SVM_Polinomial", "Red_Neuronal"),
  min_votes = 4L
)
}
\arguments{
\item{pred}{Predicted class vector.}
\item{obs}{Observed ground-truth class vector.}
\item{positive_class}{Positive class label.}
\item{model_name}{Display name for the model.}
\item{score}{Optional numeric score or positive-class probability vector for ROC-AUC.}
\item{response}{Observed binary response vector.}
\item{predictor}{Numeric predictor vector.}
\item{df}{Input data frame.}
\item{response_col}{Name of the binary response column.}
\item{metric_cols}{Character vector of continuous metric column names.}
\item{df_votes}{Data frame of binary \code{"TP"}/\code{"FP"} votes per metric.}
\item{weight_range}{Integer grid of candidate weights.}
\item{df_predictors}{Data frame of continuous predictor columns.}
\item{thresholds}{Named numeric vector of Youden thresholds, or vector of coverage thresholds for sensitivity sweeps.}
\item{weights}{Named numeric vector of heuristic weights.}
\item{global_threshold}{Numeric cutoff on the weighted sum.}
\item{ens_stack}{Fitted \code{caretStack} object.}
\item{newdata}{Data frame of predictors for stacking inference.}
\item{df_predictions}{Data frame containing voter prediction columns.}
\item{model_cols}{Character vector of voter column names.}
\item{min_votes}{Minimum positive votes required to declare \code{"Homologo_Valido"}.}
\item{pred_copies}{Integer or numeric vector of estimated orthogroup copy counts.}
\item{real_copies}{Integer or numeric vector of reference orthogroup copy counts.}
\item{sample_col}{Sample column name.}
\item{locus_col}{Orthogroup/locus column name.}
\item{class_col}{Ground-truth class column name.}
\item{df_consenso}{Data frame containing consensus verdicts and \code{pct_recovered}.}
\item{cov_threshold}{Coverage threshold percentage for counting secondary copies.}
\item{verdict_col}{Consensus verdict column name.}
\item{true_cnv_df}{Optional external data frame with \code{Sample}, \code{Locus_Ortogrupo}, and \code{CNV_Real}.}
\item{df_sp}{Prepared species data frame with \code{Validacion_Espacial}.}
\item{models_list}{Named list of trained supervised models.}
\item{youden_thresholds}{Named numeric vector of Youden thresholds.}
\item{heuristic_model}{List with \code{weights} and \code{global_threshold}.}
\item{predictors}{Character vector of active predictor names.}
\item{ml_voters}{Character vector of supervised classifiers included in the 7-model consensus.}
}
\value{
Summary data frames, numeric scalars, factors, or structured evaluation lists.
}
