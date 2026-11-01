test_that("flujo MWE CLI (01 -> 02 -> 03) se ejecuta de extremo a extremo en directorio aislado", {
  proj_root <- Sys.getenv("TARGETCNV_PROJECT_ROOT", unset = normalizePath("../.."))
  mock_dir  <- file.path(proj_root, "data", "test_data")
  skip_if_not(file.exists(file.path(mock_dir, "mock_captus_extract.tsv")), "Datos mock del MWE no disponibles")

  tmp_base <- file.path(tempdir(), paste0("targetcnv_mwe_", Sys.getpid()))
  prep_dir <- file.path(tmp_base, "prepared")
  mod_dir  <- file.path(tmp_base, "models")
  res_dir  <- file.path(tmp_base, "results")
  on.exit(unlink(tmp_base, recursive = TRUE, force = TRUE), add = TRUE)

  # Paso 01: preparación espacial
  rc1 <- system2(
    "Rscript",
    args = c(
      file.path(proj_root, "scripts", "01_prepare_data.R"),
      "--gff_ref", file.path(mock_dir, "mock_reference.gff3"),
      "--captus_gff", file.path(mock_dir, "mock_captus.gff"),
      "--captus_stats", file.path(mock_dir, "mock_captus_extract.tsv"),
      "--dict_ref", file.path(mock_dir, "mock_dictionary.tsv"),
      "--allow_multi_sample",
      "--out_dir", prep_dir
    ),
    stdout = FALSE, stderr = FALSE
  )
  expect_equal(rc1, 0L)
  expect_true(file.exists(file.path(prep_dir, "data_prepared.rds")))

  df_prep <- readRDS(file.path(prep_dir, "data_prepared.rds"))
  expect_false(any(grepl("^delta_", colnames(df_prep))))
  expect_true(all(c("TRUE_POSITIVE", "FALSE_POSITIVE") %in% df_prep$Validacion_Espacial))

  # Paso 02: entrenamiento rápido determinista con Stacking y test de Spearman
  rc2 <- system2(
    "Rscript",
    args = c(
      file.path(proj_root, "scripts", "02_train_models.R"),
      "--data", file.path(prep_dir, "data_prepared.rds"),
      "--out_models", mod_dir,
      "--cv_folds", "3",
      "--sampling", "down",
      "--seed", "42",
      "--cores", "1",
      "--fast_tournament",
      "--clean_models",
      "--no_render_report"
    ),
    stdout = FALSE, stderr = FALSE
  )
  expect_equal(rc2, 0L)
  expect_true(file.exists(file.path(mod_dir, "consensus_spec.rds")))
  expect_true(file.exists(file.path(mod_dir, "spearman_collinearity.tsv")))
  expect_true(file.exists(file.path(mod_dir, "youden_thresholds.rds")))
  expect_true(file.exists(file.path(mod_dir, "heuristic_model.rds")))

  bench_bal <- readRDS(file.path(mod_dir, "benchmark_models_balanced.rds"))
  if ("Ensamble_Stacking" %in% bench_bal$Modelo) {
    auc_st <- bench_bal$ROC_AUC[bench_bal$Modelo == "Ensamble_Stacking"][1]
    expect_gt(auc_st, 0.5)
  }

  # Paso 03: inferencia de CNV, sensibilidad 30%-80% y concordancia sobre datos etiquetados
  rc3 <- system2(
    "Rscript",
    args = c(
      file.path(proj_root, "scripts", "03_cnv_inference.R"),
      "--input", file.path(prep_dir, "data_prepared.rds"),
      "--models_dir", mod_dir,
      "--out_dir", res_dir,
      "--cov_threshold", "50.0",
      "--min_votes", "4",
      "--top_loci", "10",
      "--no_render_report"
    ),
    stdout = FALSE, stderr = FALSE
  )
  expect_equal(rc3, 0L)
  expect_true(file.exists(file.path(res_dir, "cnv_matrix.tsv")))
  expect_true(file.exists(file.path(res_dir, "cnv_threshold_sensitivity.tsv")))
  expect_true(file.exists(file.path(res_dir, "cnv_concordance_summary.tsv")))

  conc_df <- utils::read.delim(file.path(res_dir, "cnv_concordance_summary.tsv"), stringsAsFactors = FALSE)
  expect_equal(nrow(conc_df), 1L)
  expect_true(is.finite(conc_df$Spearman_Rho[1]) && conc_df$Spearman_Rho[1] > 0)
  expect_true(is.finite(conc_df$Kappa_Cuadratico[1]))
  expect_true(is.finite(conc_df$RMSE[1]) && conc_df$RMSE[1] >= 0)
  expect_true(conc_df$Concordancia_Exacta[1] >= 0 && conc_df$Concordancia_Exacta[1] <= 1)
})
