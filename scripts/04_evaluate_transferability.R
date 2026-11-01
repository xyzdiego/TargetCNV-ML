#!/usr/bin/env Rscript

# Módulo 04: Evaluación de transferibilidad filogenética inter-especie,
# concordancia de número de copias a nivel de ortogrupo y calibración (30%-80%)
# del umbral de cobertura sobre genomas a escala cromosómica y de contigs

t_start <- Sys.time()

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(caret)
  library(pROC)
})

get_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    if (basename(script_dir) == "scripts") dirname(script_dir) else script_dir
  } else {
    getwd()
  }
}

project_root <- get_project_root()
source(file.path(project_root, "R", "model_specs.R"))
source(file.path(project_root, "R", "data_processing.R"))
source(file.path(project_root, "R", "model_evaluation.R"))
source(file.path(project_root, "R", "plotting.R"))
source(file.path(project_root, "R", "report_helpers.R"))

option_list <- list(
  make_option(
    c("-m", "--models_dir"),
    type = "character",
    default = file.path(project_root, "models"),
    help = "Directorio con modelos preentrenados de la especie donante [Por defecto: models/]",
    metavar = "DIR"
  ),
  make_option(
    c("--donor_data"),
    type = "character",
    default = file.path(project_root, "data", "data_prepared.rds"),
    help = "Datos preparados de la especie donante (S. torvum) [Por defecto: data/data_prepared.rds]",
    metavar = "FILE"
  ),
  make_option(
    c("--cit_data"),
    type = "character",
    default = file.path(project_root, "data", "citrullifolium", "data_prepared.rds"),
    help = "Datos preparados de S. citrullifolium (escala cromosómica) [Por defecto: data/citrullifolium/data_prepared.rds]",
    metavar = "FILE"
  ),
  make_option(
    c("--ste_data"),
    type = "character",
    default = file.path(project_root, "data", "stellatiglandulosum", "data_prepared.rds"),
    help = "Datos preparados de S. stellatiglandulosum (escala de contigs) [Por defecto: data/stellatiglandulosum/data_prepared.rds]",
    metavar = "FILE"
  ),
  make_option(
    c("--cit_gff_ref"),
    type = "character",
    default = "",
    help = "GFF3 genómico opcional de S. citrullifolium si se requiere preparar desde archivos crudos",
    metavar = "FILE"
  ),
  make_option(
    c("--cit_captus_gff"),
    type = "character",
    default = "",
    help = "GFF de contigs de CAPTUS opcional de S. citrullifolium",
    metavar = "FILE"
  ),
  make_option(
    c("--ste_captus_gff"),
    type = "character",
    default = "",
    help = "GFF de contigs de CAPTUS opcional de S. stellatiglandulosum",
    metavar = "FILE"
  ),
  make_option(
    c("-s", "--captus_stats"),
    type = "character",
    default = "",
    help = "Archivo captus-extract_stats.tsv opcional para preparación bajo demanda desde GFF",
    metavar = "FILE"
  ),
  make_option(
    c("-d", "--dict_ref"),
    type = "character",
    default = "",
    help = "Diccionario de homología opcional (.xlsx o .tsv) para preparación bajo demanda",
    metavar = "FILE"
  ),
  make_option(
    c("--extra_species"),
    type = "character",
    default = "",
    help = "Especies adicionales separadas por ';' con formato 'Muestra|Nombre Cientifico|chromosome_o_contig|ruta_datos.rds'",
    metavar = "STR"
  ),
  make_option(
    c("-c", "--cov_threshold"),
    type = "double",
    default = 50.0,
    help = "Umbral de cobertura individual (%) para evaluación de concordancia CNV [Por defecto: 50.0]",
    metavar = "NUM"
  ),
  make_option(
    c("-v", "--min_votes"),
    type = "integer",
    default = 4L,
    help = "Mínimo de votos afirmativos para el consenso de 7 modelos [Por defecto: 4]",
    metavar = "INT"
  ),
  make_option(
    c("-o", "--out_dir"),
    type = "character",
    default = file.path(project_root, "results"),
    help = "Directorio de salida para artefactos .tsv y .rds de transferencia [Por defecto: results/]",
    metavar = "DIR"
  ),
  make_option(
    c("--render_report"),
    action = "store_true",
    default = TRUE,
    help = "Compilar el reporte HTML de transferencia (transferability_report.qmd) en modo de solo lectura [Por defecto: TRUE]"
  ),
  make_option(
    c("--no_render_report"),
    action = "store_false",
    dest = "render_report",
    help = "Omitir la compilación del reporte HTML de transferencia"
  ),
  make_option(
    c("-r", "--report_file"),
    type = "character",
    default = file.path(project_root, "reports", "reporte_transferencia.html"),
    help = "Ruta del reporte HTML de transferencia [Por defecto: reports/reporte_transferencia.html]",
    metavar = "FILE"
  ),
  make_option(
    c("--sources"),
    type = "character",
    default = "",
    help = "Ruta opcional al archivo TSV con citas de fuentes genómicas y de homología",
    metavar = "FILE"
  )
)

opt_parser <- OptionParser(
  usage = "Uso: %prog [opciones]",
  option_list = option_list,
  description = "Evalúa sin reentrenamiento la transferibilidad filogenética de los modelos de S. torvum sobre especies receptoras, calcula la concordancia de número de copias por ortogrupo (Spearman rho, Kappa cuadrático ponderado, RMSE y concordancia exacta) y calibra --cov_threshold entre 30% y 80%."
)

opt <- parse_args(opt_parser)

if (!dir.exists(opt$models_dir)) {
  stop(sprintf(
    "[ERROR] No existe el directorio de modelos: %s.\nEjecute primero: ./targetcnv.sh train",
    opt$models_dir
  ), call. = FALSE)
}

bal_path   <- file.path(opt$models_dir, "models_balanced.rds")
all_path   <- file.path(opt$models_dir, "all_trained_models.rds")
unbal_path <- file.path(opt$models_dir, "models_unbalanced.rds")
youden_f   <- file.path(opt$models_dir, "youden_thresholds.rds")
heur_f     <- file.path(opt$models_dir, "heuristic_model.rds")

if ((!file.exists(bal_path) && !file.exists(all_path)) || !file.exists(youden_f) || !file.exists(heur_f)) {
  stop(sprintf(
    paste0(
      "[ERROR] Faltan artefactos de entrenamiento en '%s' (models_balanced.rds, youden_thresholds.rds o heuristic_model.rds).\n",
      "        Ejecute primero './targetcnv.sh train' antes de evaluar la transferencia filogenética."
    ),
    opt$models_dir
  ), call. = FALSE)
}

if (!dir.exists(opt$out_dir)) {
  dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
}

cat("========================================================\n")
cat(" TargetCNV-ML: Módulo 04 - Transferencia Filogenética y Concordancia CNV\n")
cat("========================================================\n")
cat(sprintf("Directorio de Modelos       : %s\n", opt$models_dir))
cat(sprintf("Datos Donante (S. torvum)   : %s\n", opt$donor_data))
cat(sprintf("Datos S. citrullifolium     : %s\n", opt$cit_data))
cat(sprintf("Datos S. stellatiglandulosum: %s\n", opt$ste_data))
cat(sprintf("Directorio de Salida        : %s\n", opt$out_dir))
cat(sprintf("Umbral Cobertura Homólogo   : %.1f%%\n", opt$cov_threshold))
cat("--------------------------------------------------------\n")

# 1. Carga en modo de solo lectura de los modelos preentrenados en S. torvum
cat("[1/6] Cargando modelos preentrenados en la especie donante...\n")
models_bal <- normalize_model_names(readRDS(if (file.exists(bal_path)) bal_path else all_path))
models_unbal <- if (file.exists(unbal_path)) normalize_model_names(readRDS(unbal_path)) else models_bal
youden_thresholds <- readRDS(youden_f)
heuristic_model   <- readRDS(heur_f)

spec_file <- file.path(opt$models_dir, "consensus_spec.rds")
consensus_spec <- if (file.exists(spec_file)) readRDS(spec_file) else NULL
predictores <- if (!is.null(consensus_spec$predictors)) {
  consensus_spec$predictors
} else {
  c("pct_recovered", "pct_identity", "score", "wscore")
}
ml_voters <- if (!is.null(consensus_spec$ml_voters)) {
  consensus_spec$ml_voters
} else {
  c("ElasticNet", "NaiveBayes", "SVM_Polinomial", "Red_Neuronal")
}

bench_donor_file <- file.path(opt$models_dir, "benchmark_models_balanced.rds")
bench_donor <- if (file.exists(bench_donor_file)) {
  b <- readRDS(bench_donor_file)
  b$Modelo <- normalize_model_names(as.character(b$Modelo))
  b
} else {
  NULL
}

sources_df <- load_sources_citations(opt$sources, proj_root = project_root)

# Función auxiliar para cargar un dataset preparado o generarlo desde data/reference/
resolve_species_dataset <- function(data_path,
                                    sample_id,
                                    genome_scale,
                                    gff_ref_arg = "",
                                    captus_gff_arg = "") {
  if (!is.null(data_path) && nzchar(data_path) && file.exists(data_path)) {
    df_loaded <- if (grepl("\\.rds$", data_path, ignore.case = TRUE)) {
      readRDS(data_path)
    } else {
      utils::read.delim(data_path, stringsAsFactors = FALSE)
    }
    return(list(df = df_loaded, origin = basename(data_path)))
  }

  tsv_alt <- sub("\\.rds$", ".tsv", data_path)
  if (!is.null(tsv_alt) && nzchar(tsv_alt) && file.exists(tsv_alt)) {
    return(list(df = utils::read.delim(tsv_alt, stringsAsFactors = FALSE), origin = basename(tsv_alt)))
  }

  ref_dir <- file.path(project_root, "data", "reference")
  stats_f <- if (nzchar(opt$captus_stats) && file.exists(opt$captus_stats)) {
    opt$captus_stats
  } else {
    file.path(ref_dir, "captus-extract_stats.tsv")
  }
  dict_f <- if (nzchar(opt$dict_ref) && file.exists(opt$dict_ref)) {
    opt$dict_ref
  } else {
    file.path(ref_dir, "Genes_Capture_May_2026.xlsx")
  }
  cap_f <- if (nzchar(captus_gff_arg) && file.exists(captus_gff_arg)) {
    captus_gff_arg
  } else {
    file.path(ref_dir, paste0(sample_id, "_hit_contigs.gff"))
  }
  ref_f <- if (nzchar(gff_ref_arg) && file.exists(gff_ref_arg)) {
    gff_ref_arg
  } else if (genome_scale == "chromosome") {
    file.path(ref_dir, paste0(sample_id, "_genes.gff3"))
  } else {
    NULL
  }

  if (!file.exists(cap_f) || !file.exists(stats_f) || !file.exists(dict_f)) {
    stop(sprintf(
      "[ERROR] No se encontró el dataset preparado '%s' ni los archivos de referencia para la muestra '%s'.",
      data_path, sample_id
    ), call. = FALSE)
  }

  df_prep <- prepare_transfer_species_data(
    captus_gff = cap_f,
    captus_stats = stats_f,
    dict_ref = dict_f,
    sample_name = sample_id,
    genome_scale = genome_scale,
    gff_ref = ref_f
  )
  list(df = df_prep, origin = basename(cap_f))
}

build_species_info <- function(df_sp, sample_id, sci_name, short_name, genome_scale, origin_file) {
  cit_row <- lookup_source_entry(sources_df, sample_id = sample_id, genome_scale = genome_scale)
  has_cit <- !is.null(cit_row)
  algo_txt <- if (genome_scale == "chromosome") {
    "Intersección espacial estricta por cromosoma y hebra (GenomicRanges::findOverlaps, ignore.strand = FALSE) frente a la anotación GFF3 de referencia, seguida de deduplicación intra-locus."
  } else {
    "Extracción directa de identificadores génicos desde el atributo Target de Scipio/CAPTUS, reconciliación de grupos jerárquicos HOG_ID con ortogrupos pangenómicos OG_ID mediante la tabla puente del diseño de sondas y deduplicación intra-locus."
  }

  list(
    id            = sample_id,
    sci_name      = sci_name,
    short_name    = short_name,
    genome_scale  = genome_scale,
    origin_file   = origin_file,
    has_citation  = has_cit,
    genome_source = if (has_cit) as.character(cit_row$genome_source[1]) else "",
    db_info       = if (has_cit) as.character(cit_row$homology_source[1]) else "",
    algo_info     = algo_txt,
    n_ogs_captus  = attr(df_sp, "n_ogs_captus"),
    n_ogs_bridge  = attr(df_sp, "n_ogs_bridge")
  )
}

# 2. Carga de especies receptoras
cat("[2/6] Cargando conjuntos validados de especies receptoras...\n")
species_specs <- list(
  scitrullifolium = list(
    id = "Scitrullifolium",
    sci_name = "Solanum citrullifolium",
    short_name = "S. citrullifolium",
    genome_scale = "chromosome",
    data_path = opt$cit_data,
    gff_ref = opt$cit_gff_ref,
    captus_gff = opt$cit_captus_gff,
    tsv_out = "transferability_scitrullifolium.tsv"
  ),
  s_stella = list(
    id = "S_stella",
    sci_name = "Solanum stellatiglandulosum",
    short_name = "S. stellatiglandulosum",
    genome_scale = "contig",
    data_path = opt$ste_data,
    gff_ref = "",
    captus_gff = opt$ste_captus_gff,
    tsv_out = "transferability_s_stella.tsv"
  )
)

if (!is.null(opt$extra_species) && nzchar(trimws(opt$extra_species))) {
  entries <- strsplit(opt$extra_species, ";", fixed = TRUE)[[1]]
  for (ent in entries) {
    parts <- trimws(strsplit(ent, "|", fixed = TRUE)[[1]])
    if (length(parts) >= 4L) {
      sp_k <- gsub("[^a-z0-9_]", "_", tolower(parts[1]))
      sc_k <- if (tolower(parts[3]) %in% c("contig", "contigs", "target_attr")) "contig" else "chromosome"
      species_specs[[sp_k]] <- list(
        id = parts[1],
        sci_name = parts[2],
        short_name = parts[2],
        genome_scale = sc_k,
        data_path = parts[4],
        gff_ref = "",
        captus_gff = "",
        tsv_out = sprintf("transferability_%s.tsv", sp_k)
      )
    }
  }
}

# 3. Evaluación de transferibilidad por especie (escenarios balanceado y desbalanceado)
cat("[3/6] Evaluando desempeño de clasificación inter-especie...\n")
transfer_results <- list()
multi_rows <- list()

if (!is.null(bench_donor) && nrow(bench_donor) > 0) {
  multi_rows[["donor"]] <- bench_donor %>%
    dplyr::mutate(
      Especie = "S. torvum (Donante CV)",
      Escala_Genoma = "Escala Cromosómica",
      .before = 1
    )
}

concordance_rows <- list()
concordance_pairs_list <- list()
sensitivity_rows <- list()

# Evaluar también concordancia de copias y curva de sensibilidad en la especie donante si está disponible
if (!is.null(opt$donor_data) && nzchar(opt$donor_data) && file.exists(opt$donor_data)) {
  df_donor_raw <- if (grepl("\\.rds$", opt$donor_data, ignore.case = TRUE)) {
    readRDS(opt$donor_data)
  } else {
    utils::read.delim(opt$donor_data, stringsAsFactors = FALSE)
  }
  eval_donor <- evaluate_species_transferability(
    df_sp = df_donor_raw,
    models_list = models_bal,
    youden_thresholds = youden_thresholds,
    heuristic_model = heuristic_model,
    predictors = predictores,
    ml_voters = ml_voters,
    min_votes = opt$min_votes
  )
  if ("Consenso_7_Modelos" %in% eval_donor$benchmark$Modelo && !is.null(multi_rows[["donor"]])) {
    row_c7_donor <- eval_donor$benchmark %>%
      dplyr::filter(Modelo == "Consenso_7_Modelos") %>%
      dplyr::mutate(Especie = "S. torvum (Donante CV)", Escala_Genoma = "Escala Cromosómica", .before = 1)
    if (!"Consenso_7_Modelos" %in% multi_rows[["donor"]]$Modelo) {
      multi_rows[["donor"]] <- dplyr::bind_rows(multi_rows[["donor"]], row_c7_donor)
    }
  }
  if ("Veredicto_Consenso" %in% colnames(eval_donor$df_consensus)) {
    conc_d <- evaluate_cnv_concordance(
      df_consenso = eval_donor$df_consensus,
      cov_threshold = opt$cov_threshold,
      sample_col = "Sample",
      locus_col = "Locus_Ortogrupo",
      verdict_col = "Veredicto_Consenso"
    )
    if (nrow(conc_d$summary) > 0) {
      concordance_rows[["storvum"]] <- conc_d$summary %>%
        dplyr::mutate(Especie = "S. torvum (Donante)", Escala_Genoma = "Escala Cromosómica", .before = 1)
      concordance_pairs_list[["storvum"]] <- conc_d$pairs %>%
        dplyr::mutate(Especie = "S. torvum (Donante)", Escala_Genoma = "Escala Cromosómica", .before = 1)
    }
    sens_d <- evaluate_cnv_threshold_sensitivity(
      df_consenso = eval_donor$df_consensus,
      thresholds = seq(30, 80, by = 5),
      sample_col = "Sample",
      locus_col = "Locus_Ortogrupo",
      verdict_col = "Veredicto_Consenso"
    )
    sensitivity_rows[["storvum"]] <- sens_d %>%
      dplyr::mutate(Especie = "S. torvum (Donante)", Escala_Genoma = "Escala Cromosómica", .before = 1)
  }
}

for (sp_key in names(species_specs)) {
  sp_cfg <- species_specs[[sp_key]]
  loaded <- resolve_species_dataset(
    data_path = sp_cfg$data_path,
    sample_id = sp_cfg$id,
    genome_scale = sp_cfg$genome_scale,
    gff_ref_arg = sp_cfg$gff_ref,
    captus_gff_arg = sp_cfg$captus_gff
  )
  df_sp <- loaded$df
  sp_info <- build_species_info(
    df_sp = df_sp,
    sample_id = sp_cfg$id,
    sci_name = sp_cfg$sci_name,
    short_name = sp_cfg$short_name,
    genome_scale = sp_cfg$genome_scale,
    origin_file = loaded$origin
  )

  res_bal_sp <- evaluate_species_transferability(
    df_sp = df_sp,
    models_list = models_bal,
    youden_thresholds = youden_thresholds,
    heuristic_model = heuristic_model,
    predictors = predictores,
    ml_voters = ml_voters,
    min_votes = opt$min_votes
  )
  res_unbal_sp <- evaluate_species_transferability(
    df_sp = df_sp,
    models_list = models_unbal,
    youden_thresholds = youden_thresholds,
    heuristic_model = heuristic_model,
    predictors = predictores,
    ml_voters = ml_voters,
    min_votes = opt$min_votes
  )

  sc_lbl <- if (sp_cfg$genome_scale == "chromosome") "Escala Cromosómica" else "Escala de Contigs"
  cat(sprintf(
    "   -> %-24s [%s]: n = %d (TP = %d, FP = %d) | Modelos evaluados = %d\n",
    sp_cfg$short_name, sc_lbl, nrow(res_bal_sp$df_clean),
    sum(res_bal_sp$df_clean$Clase_Real == "TP"),
    sum(res_bal_sp$df_clean$Clase_Real == "FP"),
    nrow(res_bal_sp$benchmark)
  ))

  conc_sp <- NULL
  sens_sp <- NULL
  if ("Veredicto_Consenso" %in% colnames(res_bal_sp$df_consensus)) {
    conc_sp <- evaluate_cnv_concordance(
      df_consenso = res_bal_sp$df_consensus,
      cov_threshold = opt$cov_threshold,
      sample_col = "Sample",
      locus_col = "Locus_Ortogrupo",
      verdict_col = "Veredicto_Consenso"
    )
    if (nrow(conc_sp$summary) > 0) {
      concordance_rows[[sp_key]] <- conc_sp$summary %>%
        dplyr::mutate(Especie = sp_cfg$short_name, Escala_Genoma = sc_lbl, .before = 1)
      concordance_pairs_list[[sp_key]] <- conc_sp$pairs %>%
        dplyr::mutate(Especie = sp_cfg$short_name, Escala_Genoma = sc_lbl, .before = 1)
      s_r <- conc_sp$summary[1, ]
      cat(sprintf(
        "      Concordancia CNV (%s, OGs = %d): Exacta = %.2f%% | Spearman rho = %.4f | Kappa_w = %.4f | RMSE = %.4f\n",
        sp_cfg$short_name, s_r$N_Ortogrupos, s_r$Concordancia_Exacta * 100,
        s_r$Spearman_Rho, s_r$Kappa_Cuadratico, s_r$RMSE
      ))
    }

    sens_sp <- evaluate_cnv_threshold_sensitivity(
      df_consenso = res_bal_sp$df_consensus,
      thresholds = seq(30, 80, by = 5),
      sample_col = "Sample",
      locus_col = "Locus_Ortogrupo",
      verdict_col = "Veredicto_Consenso"
    )
    sensitivity_rows[[sp_key]] <- sens_sp %>%
      dplyr::mutate(Especie = sp_cfg$short_name, Escala_Genoma = sc_lbl, .before = 1)
  }

  utils::write.table(
    res_bal_sp$benchmark,
    file = file.path(opt$out_dir, sp_cfg$tsv_out),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  multi_rows[[sp_key]] <- res_bal_sp$benchmark %>%
    dplyr::mutate(Especie = sp_cfg$short_name, Escala_Genoma = sc_lbl, .before = 1)

  transfer_results[[sp_key]] <- list(
    info = sp_info,
    bal = res_bal_sp,
    unbal = res_unbal_sp,
    cnv_concordance = conc_sp,
    cnv_sensitivity = sens_sp
  )
}

# 4. Exportación de tablas consolidadas y artefacto de solo lectura para Quarto
cat("[4/6] Exportando resultados consolidados de transferibilidad y concordancia CNV...\n")
df_multiespecie_all <- dplyr::bind_rows(multi_rows)
df_cnv_concordance  <- dplyr::bind_rows(concordance_rows)
df_cnv_pairs        <- dplyr::bind_rows(concordance_pairs_list)
df_cnv_sensitivity  <- dplyr::bind_rows(sensitivity_rows)

utils::write.table(
  df_multiespecie_all,
  file = file.path(opt$out_dir, "transferability_all_species.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
if (nrow(df_cnv_concordance) > 0) {
  utils::write.table(
    df_cnv_concordance,
    file = file.path(opt$out_dir, "transferability_cnv_concordance.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}
if (nrow(df_cnv_pairs) > 0) {
  utils::write.table(
    df_cnv_pairs,
    file = file.path(opt$out_dir, "transferability_cnv_pairs.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}
if (nrow(df_cnv_sensitivity) > 0) {
  utils::write.table(
    df_cnv_sensitivity,
    file = file.path(opt$out_dir, "transferability_threshold_sensitivity.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}

# 5. Gráfico comparativo de calibración empírica del umbral --cov_threshold (30% a 80%)
cat("[5/6] Generando curva multiespecie de calibración del umbral --cov_threshold...\n")
if (nrow(df_cnv_sensitivity) > 0) {
  p_tr_sens <- plot_cnv_threshold_sensitivity(
    df_sens = df_cnv_sensitivity,
    active_threshold = opt$cov_threshold,
    title = "Calibración Empírica Multiespecie del Umbral de Cobertura (--cov_threshold)"
  )
  tryCatch({
    ggplot2::ggsave(file.path(opt$out_dir, "plot_transfer_threshold_sensitivity.pdf"), plot = p_tr_sens, width = 8.5, height = 5.2)
    ggplot2::ggsave(file.path(opt$out_dir, "plot_transfer_threshold_sensitivity.png"), plot = p_tr_sens, width = 8.5, height = 5.2, dpi = 300)
  }, error = function(e) {
    warning(sprintf("Aviso al guardar curva de sensibilidad multiespecie: %s", e$message))
  })
}

results_rds_path <- file.path(opt$out_dir, "transferability_results.rds")
transfer_bundle <- list(
  transfer_results = transfer_results,
  df_multiespecie_all = df_multiespecie_all,
  cnv_concordance_summary = df_cnv_concordance,
  cnv_concordance_pairs = df_cnv_pairs,
  cnv_threshold_sensitivity = df_cnv_sensitivity,
  cov_threshold = opt$cov_threshold,
  min_votes = opt$min_votes,
  models_dir = normalizePath(opt$models_dir, mustWork = FALSE)
)
saveRDS(transfer_bundle, file = results_rds_path)

# 6. Compilación del reporte HTML de transferencia en modo de solo lectura
cat("[6/6] Compilando reporte HTML de transferencia con Quarto...\n")
qmd_tr <- file.path(project_root, "transferability_report.qmd")
if (!opt$render_report) {
  cat("   -> Compilación omitida (--no_render_report).\n")
} else if (!file.exists(qmd_tr)) {
  cat(sprintf("   -> [AVISO] Compilación omitida: no existe la plantilla %s\n", qmd_tr))
} else if (!nzchar(Sys.which("quarto"))) {
  cat("   -> [AVISO] Compilación omitida: 'quarto' no está disponible en el PATH.\n")
} else {
  out_html <- normalizePath(opt$report_file, mustWork = FALSE)
  dir.create(dirname(out_html), recursive = TRUE, showWarnings = FALSE)
  rendered_in_root <- file.path(dirname(qmd_tr), basename(out_html))

  src_arg <- if (!is.null(opt$sources) && nzchar(trimws(opt$sources))) {
    sprintf(" -P sources_file:%s", shQuote(normalizePath(opt$sources, mustWork = FALSE)))
  } else {
    ""
  }
  cmd_quarto <- sprintf(
    "quarto render %s --output %s -P project_root:%s -P models_dir:%s -P out_dir:%s -P results_rds:%s%s",
    shQuote(qmd_tr),
    shQuote(basename(out_html)),
    shQuote(project_root),
    shQuote(normalizePath(opt$models_dir, mustWork = FALSE)),
    shQuote(normalizePath(opt$out_dir, mustWork = FALSE)),
    shQuote(normalizePath(results_rds_path, mustWork = FALSE)),
    src_arg
  )
  cat(sprintf("   -> Ejecutando Quarto: %s\n", cmd_quarto))
  res_q <- system(cmd_quarto)
  if (res_q == 0) {
    if (file.exists(rendered_in_root) && normalizePath(rendered_in_root, mustWork = FALSE) != out_html) {
      file.copy(rendered_in_root, out_html, overwrite = TRUE)
      unlink(rendered_in_root)
    }
    cat(sprintf("   -> [OK] Reporte HTML de transferencia generado: %s\n", out_html))
  } else {
    cat(sprintf("   -> [AVISO] Falló el renderizado con Quarto (código de retorno: %d).\n", res_q))
  }
}

elapsed_sec <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

cat("--------------------------------------------------------\n")
if (nrow(df_cnv_concordance) > 0) {
  cat(" CONCORDANCIA DE NUMERO DE COPIAS A NIVEL DE ORTOGRUPO:\n")
  print(
    df_cnv_concordance %>%
      dplyr::mutate(
        Concordancia_Exacta = round(Concordancia_Exacta * 100, 2),
        Spearman_Rho = round(Spearman_Rho, 4),
        Kappa_Cuadratico = round(Kappa_Cuadratico, 4),
        RMSE = round(RMSE, 4)
      ),
    row.names = FALSE
  )
}
cat("========================================================\n")
cat(sprintf("[OK] Artefactos de transferencia guardados en: %s (%.2f s)\n", opt$out_dir, elapsed_sec))

