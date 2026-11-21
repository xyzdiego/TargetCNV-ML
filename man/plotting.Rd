\name{plotting}
\alias{plotting}
\alias{get_model_color_palette}
\alias{plot_confusion_matrix}
\alias{plot_model_varimp}
\alias{plot_neural_net_topology}
\alias{plot_metric_densities}
\alias{plot_spearman_collinearity}
\alias{plot_cnv_threshold_sensitivity}
\alias{plot_transfer_raincloud}
\alias{plot_model_benchmarks}
\alias{plot_consensus_distribution}
\alias{plot_cnv_heatmap}
\alias{plot_roc_curves}
\title{Diagnostic and summary plots for TargetCNV-ML}
\description{
Generates \code{ggplot2} and \code{NeuralNetTools} figures for confusion matrices,
model-specific variable importance, neural network topologies, metric densities,
Spearman collinearity heatmaps, coverage threshold sensitivity curves, raincloud
score distributions, consensus breakdowns, CNV heatmaps, and ROC curves.
}
\usage{
get_model_color_palette(model_name)

plot_confusion_matrix(
  cm_table,
  title = "Matriz de Confusi\u00f3n",
  subtitle = NULL,
  high_color = NULL,
  low_color = NULL
)

plot_model_varimp(
  model_obj,
  df_data = NULL,
  model_name = "Modelo",
  high_color = NULL
)

plot_neural_net_topology(
  model_obj,
  title = "Arquitectura de la Red Neuronal"
)

plot_metric_densities(
  df,
  response_col = "Clase",
  metric_cols = c("pct_recovered", "pct_identity", "score", "wscore"),
  thresholds = NULL,
  colors = c(
    TRUE_POSITIVE = "#27ae60",
    FALSE_POSITIVE = "#e74c3c",
    TP = "#27ae60",
    FP = "#e74c3c"
  ),
  title = "Separaci\u00f3n entre Verdaderos Positivos (TP) y Falsos Positivos (FP)",
  subtitle = "Distribuci\u00f3n de densidad y puntos de corte de Youden (l\u00ednea punteada)",
  legend_position = "bottom",
  class_col = NULL
)

plot_spearman_collinearity(
  cor_matrix,
  p_matrix = NULL,
  title = "Matriz de Colinealidad por Rangos de Spearman"
)

plot_cnv_threshold_sensitivity(
  df_sens,
  active_threshold = 50,
  title = "Sensibilidad del Conteo de Copias al Umbral de Cobertura (30\% - 80\%)"
)

plot_transfer_raincloud(df_sc, mod_name, sp_name)

plot_model_benchmarks(
  df_benchmarks,
  metric_col = "F1_Score",
  metric_cols = NULL,
  title = "Comparativa de Rendimiento de Modelos",
  high_color = "#27ae60"
)

plot_consensus_distribution(
  df_consensus,
  verdict_col = "Veredicto_Consenso",
  title = "Clasificaci\u00f3n de Secuencias por Consenso Multimodelo",
  subtitle = "Consenso de 7 Modelos (Hom\u00f3logo V\u00e1lido vs. Ruido Descartado)"
)

plot_cnv_heatmap(
  df_cnv_annotated = NULL,
  cnv_matrix = NULL,
  sample_col = "Sample",
  locus_col = "Locus_Ortogrupo",
  top_n_loci = 30,
  title = "Matriz de Variaci\u00f3n en el N\u00famero de Copias (CNV)",
  subtitle = "Conteo de copias funcionales estimadas tras el filtro de cobertura"
)

plot_roc_curves(roc_list, title = "Curvas ROC Comparativas")
}
\arguments{
\item{model_name}{Character model identifier.}
\item{cm_table}{2x2 confusion matrix table.}
\item{title}{Plot title.}
\item{subtitle}{Optional plot subtitle.}
\item{high_color}{Hex color for high values.}
\item{low_color}{Hex color for low values.}
\item{model_obj}{Trained \code{caret::train} object or named numeric vector of weights.}
\item{df_data}{Optional training data frame used for permutation importance fallback.}
\item{df}{Data frame with continuous metrics and binary class labels.}
\item{response_col}{Response column name.}
\item{metric_cols}{Character vector of metric column names.}
\item{thresholds}{Optional named vector of Youden thresholds.}
\item{colors}{Named color vector for \code{TP} and \code{FP} classes.}
\item{legend_position}{Legend placement string.}
\item{class_col}{Optional alias for \code{response_col}.}
\item{cor_matrix}{Symmetric Spearman correlation matrix or list returned by \code{compute_spearman_collinearity()}.}
\item{p_matrix}{Optional matrix of pairwise p-values.}
\item{df_sens}{Data frame returned by \code{evaluate_cnv_threshold_sensitivity()}.}
\item{active_threshold}{Operating coverage threshold highlighted with a vertical dashed line.}
\item{df_sc}{Data frame with \code{Clase_Real} and \code{Puntaje}.}
\item{mod_name}{Model name.}
\item{sp_name}{Short species display name.}
\item{df_benchmarks}{Benchmark summary table.}
\item{metric_col}{Single metric column name for univariate ranking plots.}
\item{df_consensus}{Data frame returned by \code{consolidar_consenso()}.}
\item{verdict_col}{Consensus verdict column name.}
\item{df_cnv_annotated}{Annotated CNV sequence table.}
\item{cnv_matrix}{Wide orthogroup-by-sample CNV count table.}
\item{sample_col}{Sample column name.}
\item{locus_col}{Orthogroup/locus column name.}
\item{top_n_loci}{Maximum number of orthogroups shown in the heatmap.}
\item{roc_list}{Named list of \code{pROC::roc} objects.}
}
\value{
A \code{ggplot} object, a color palette list, or \code{NULL} (invisibly) for base-graphics neural network plots.
}
