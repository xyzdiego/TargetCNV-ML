#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(testthat)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(caret)
  library(pROC)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
proj_root <- if (length(file_arg) > 0) {
  test_dir_path <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
  if (basename(test_dir_path) == "tests") dirname(test_dir_path) else test_dir_path
} else if (file.exists("DESCRIPTION")) {
  normalizePath(".")
} else {
  normalizePath("..")
}

for (mod_f in c("model_specs.R", "data_processing.R", "model_evaluation.R", "plotting.R", "report_helpers.R")) {
  source(file.path(proj_root, "R", mod_f))
}

Sys.setenv(TARGETCNV_PROJECT_ROOT = proj_root)

testthat::test_dir(
  path = file.path(proj_root, "tests", "testthat"),
  reporter = "summary",
  stop_on_failure = TRUE
)

