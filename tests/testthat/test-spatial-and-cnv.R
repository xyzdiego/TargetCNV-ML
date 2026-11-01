test_that("cruce espacial GFF respeta cromosoma armonizado, hebra y deduplicación intra-locus", {
  ref_raw <- GenomicRanges::GRanges(
    seqnames = c("Storvum_Chr01", "Storvum_Chr01", "Storvum_Chr02"),
    ranges   = IRanges::IRanges(start = c(1000, 5000, 2000), end = c(2500, 6500, 3500)),
    strand   = c("+", "-", "+"),
    type     = c("gene", "gene", "gene"),
    ID       = c("gene:Sol01g0010", "gene:Sol01g0020", "gene:Sol02g0010")
  )
  ref_gr <- filter_gff_features(ref_raw, feature_type = "gene")

  captus_gr <- GenomicRanges::GRanges(
    seqnames = c("Chr01", "Chr01", "Chr01", "Chr02"),
    ranges   = IRanges::IRanges(start = c(1100, 1200, 5200, 8000), end = c(2000, 1900, 6100, 9000)),
    strand   = c("+", "+", "+", "+"),
    Name     = c("OG0001__00", "OG0001__01", "OG0002__00", "OG0003__00")
  )

  cruces <- find_spatial_overlaps(query_gff = captus_gr, subject_gff = ref_gr, ignore_strand = FALSE)

  # OG0002__00 está en hebra '+' mientras que Sol01g0020 está en '-', por lo que no debe cruzar
  # OG0003__00 cae en región intergénica distal (8000-9000), tampoco cruza
  expect_equal(nrow(cruces), 2L)
  expect_setequal(cruces$ID_Captus, c("OG0001__00", "OG0001__01"))
  expect_true(all(cruces$Locus_Ref_Limpio == "Sol01g0010"))

  dict_df <- data.frame(
    Locus_Benoit = c("Sol01g0010", "Sol01g0020"),
    OG_Benoit    = c("OG0001", "OG0002"),
    stringsAsFactors = FALSE
  )

  clasif <- classify_spatial_hits(
    df_cruces = cruces,
    dict_ref = dict_df,
    locus_col = "Locus_Ref_Limpio",
    dict_locus_col = "Locus_Benoit",
    dict_og_col = "OG_Benoit"
  )
  expect_true(all(clasif$Validacion_Espacial == "TRUE_POSITIVE"))

  # Dos contigs sobre el mismo locus Sol01g0010: resolve_redundancies conserva solo el de mayor identidad como TP
  df_integrado <- data.frame(
    ID_GFF_Compatible       = c("OG0001__00", "OG0001__01", "OG0003__00"),
    Locus_Ref_Limpio        = c("Sol01g0010", "Sol01g0010", NA_character_),
    Validacion_Espacial     = c("TRUE_POSITIVE", "TRUE_POSITIVE", "FALSE_POSITIVE"),
    Doble_Control_Biologico = c("1_Hit_Espacial_Valido", "1_Hit_Espacial_Valido", "4_Ruido_Intergenico_o_Intron"),
    pct_identity            = c(97.5, 89.2, 61.0),
    ref_len_matched         = c(450L, 210L, 120L),
    pct_recovered           = c(91.0, 42.0, 25.0),
    score                   = c(0.95, 0.72, 0.31),
    wscore                  = c(0.90, 0.58, 0.20),
    stringsAsFactors        = FALSE
  )

  dedup <- resolve_redundancies(df_integrado, locus_col = "Locus_Ref_Limpio", class_col = "Validacion_Espacial")
  expect_equal(sum(dedup$Validacion_Espacial == "TRUE_POSITIVE"), 1L)
  expect_equal(dedup$ID_GFF_Compatible[dedup$Validacion_Espacial == "TRUE_POSITIVE"], "OG0001__00")
  expect_equal(
    dedup$Doble_Control_Biologico[dedup$ID_GFF_Compatible == "OG0001__01"],
    "5_TP_Redundante_Degradado"
  )
  expect_false(any(grepl("^delta_", colnames(dedup))))
})

test_that("extract_cnv_matrix resuelve exactamente los casos borde de conteo de copias", {
  # Caso A: 0 secuencias válidas
  df_empty <- data.frame(
    Sample             = c("S1", "S1"),
    Locus_Ortogrupo    = c("OG1", "OG2"),
    pct_recovered      = c(85.0, 90.0),
    Veredicto_Consenso = c("Ruido_Descartado", "Ruido_Descartado"),
    stringsAsFactors   = FALSE
  )
  res_empty <- extract_cnv_matrix(df_empty, cov_threshold = 50.0)
  expect_equal(nrow(res_empty$cnv_matrix), 0L)
  expect_equal(nrow(res_empty$functional_genes), 0L)

  # Casos B, C y D combinados:
  # - OG_LowPrimary en S1: única secuencia válida con 28% (< 50%) -> debe contarse como 1 Copia_Principal
  # - OG_Boundary en S1: principal 92%, homólogo 50.0% (exacto en el borde -> Homologo_Confirmado),
  #                      fragmento 49.99% (< 50% -> No_Contado) -> conteo exacto = 2 copias
  # - OG_Boundary en S2: ausente -> debe rellenarse con 0 copias en cnv_matrix
  df_edges <- data.frame(
    Sample             = c("S1", "S1", "S1", "S1", "S2"),
    Locus_Ortogrupo    = c("OG_LowPrimary", "OG_Boundary", "OG_Boundary", "OG_Boundary", "OG_LowPrimary"),
    ID_Secuencia       = c("seq1", "seq2", "seq3", "seq4", "seq5"),
    pct_recovered      = c(28.0, 92.0, 50.0, 49.99, 88.0),
    Veredicto_Consenso = rep("Homologo_Valido", 5),
    stringsAsFactors   = FALSE
  )

  res_edges <- extract_cnv_matrix(df_edges, cov_threshold = 50.0)
  mat <- res_edges$cnv_matrix

  row_low <- mat[mat$Locus_Ortogrupo == "OG_LowPrimary", , drop = FALSE]
  row_bnd <- mat[mat$Locus_Ortogrupo == "OG_Boundary", , drop = FALSE]

  expect_equal(as.integer(row_low$S1), 1L)
  expect_equal(as.integer(row_low$S2), 1L)
  expect_equal(as.integer(row_bnd$S1), 2L)
  expect_equal(as.integer(row_bnd$S2), 0L)

  ann_bnd <- res_edges$cnv_annotated %>%
    dplyr::filter(Sample == "S1", Locus_Ortogrupo == "OG_Boundary") %>%
    dplyr::arrange(dplyr::desc(pct_recovered))
  expect_equal(
    ann_bnd$Clasificacion_Copia,
    c("Copia_Principal", "Homologo_Confirmado", "No_Contado")
  )
})

test_that("métricas de concordancia CNV por ortogrupo y curva 30%-80% cumplen invariantes matemáticos", {
  df_sim <- data.frame(
    Sample              = c("S1", "S1", "S1", "S1", "S1"),
    Locus_Ortogrupo     = c("OG1", "OG2", "OG2", "OG3", "OG4"),
    pct_recovered       = c(94.0, 91.0, 62.0, 35.0, 88.0),
    Veredicto_Consenso  = c("Homologo_Valido", "Homologo_Valido", "Homologo_Valido", "Ruido_Descartado", "Homologo_Valido"),
    Validacion_Espacial = c("TRUE_POSITIVE", "TRUE_POSITIVE", "TRUE_POSITIVE", "FALSE_POSITIVE", "TRUE_POSITIVE"),
    stringsAsFactors    = FALSE
  )

  conc <- evaluate_cnv_concordance(df_sim, cov_threshold = 50.0)
  expect_equal(conc$summary$N_Ortogrupos, 4L)
  expect_equal(conc$summary$Concordancia_Exacta, 1.0)
  expect_equal(conc$summary$Spearman_Rho, 1.0, tolerance = 1e-12)
  expect_equal(conc$summary$Kappa_Cuadratico, 1.0, tolerance = 1e-12)
  expect_equal(conc$summary$RMSE, 0.0, tolerance = 1e-12)

  sens <- evaluate_cnv_threshold_sensitivity(df_sim, thresholds = seq(30, 80, by = 10))
  expect_equal(nrow(sens), 6L)
  # Al aumentar el umbral de cobertura, el conteo total estimado nunca puede crecer
  expect_true(all(diff(sens$Copias_Totales_Estimadas) <= 0))
})
