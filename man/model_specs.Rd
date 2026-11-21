\name{model_specs}
\alias{model_specs}
\alias{get_safe_xgb_method}
\alias{get_model_specs}
\alias{build_deterministic_cv_control}
\alias{normalize_model_names}
\title{Supervised classifier specifications and deterministic CV controls}
\description{
Defines the 18 supervised classification engines evaluated in TargetCNV-ML,
builds the custom \code{caret} adapter compatible with \code{xgboost >= 3.0},
normalizes legacy model identifiers, and generates explicit per-resample seeds
for deterministic cross-validation independent of worker count.
}
\usage{
get_safe_xgb_method()

get_model_specs(n_obs = 184, extended = TRUE)

build_deterministic_cv_control(
  base_ctrl,
  x,
  y,
  method,
  tuneGrid = NULL,
  tuneLength = 3L,
  preProcess = NULL,
  seed = 1234L
)

normalize_model_names(mod_list)
}
\arguments{
\item{n_obs}{Integer number of training observations used to bound tree and neighborhood grids.}
\item{extended}{Logical; if \code{TRUE} returns all 18 classifiers, otherwise the 5-model fast subset.}
\item{base_ctrl}{Base \code{caret::trainControl} object.}
\item{x}{Data frame of numeric training predictors.}
\item{y}{Two-level factor response (\code{c("TP", "FP")}).}
\item{method}{Character model name or custom \code{caret} method list.}
\item{tuneGrid}{Optional tuning grid data frame.}
\item{tuneLength}{Integer tuning grid resolution when \code{tuneGrid} is \code{NULL}.}
\item{preProcess}{Optional preprocessing vector passed to \code{caret}.}
\item{seed}{Integer master seed used to generate fold indices and per-resample \code{seeds}.}
\item{mod_list}{Named list or character vector of model names to normalize.}
}
\value{
\code{get_safe_xgb_method()} returns a custom \code{caret} model list;
\code{get_model_specs()} returns a named list of classifier specifications;
\code{build_deterministic_cv_control()} returns a list with \code{trControl} and \code{tuneGrid};
\code{normalize_model_names()} returns the input object with canonical model names.
}
