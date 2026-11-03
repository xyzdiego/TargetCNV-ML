#!/usr/bin/env Rscript
# ==============================================================================
# generate_large_mock.R
# Generador de datos sintéticos de prueba para el Ejemplo Mínimo Funcional (MWE)
# de TargetCNV-ML:
#   - 10 muestras (Sample_01 a Sample_10)
#   - 50 loci genómicos distribuidos en 5 cromosomas
#   - Matriz de variación en el número de copias (0 a 9 copias por locus):
#       * Loci con 0 copias (deleciones o ausencia del gen)
#       * Loci de copia única (1 copia)
#       * Loci con 2 a 5 copias (duplicaciones moderadas)
#       * Loci con 6 a 9 copias (expansiones multigénicas en tándem)
#   - Simulación paramétrica por clase (homólogos válidos vs. artefactos de ruido)
#     sobre las cuatro métricas continuas de CAPTUS (pct_recovered, pct_identity,
#     score y wscore).
# ==============================================================================

set.seed(2026)

# Detectar raíz del proyecto
get_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    if (basename(dirname(script_dir)) == "data") return(dirname(dirname(script_dir)))
    if (basename(script_dir) == "data") return(dirname(script_dir))
    return(script_dir)
  }
  getwd()
}

project_root <- get_project_root()
out_dir <- file.path(project_root, "data", "test_data")

cat(sprintf("Directorio de salida: %s\n", out_dir))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Parámetros de simulación ---
n_samples     <- 10
n_loci        <- 50
n_chromosomes <- 5
samples       <- paste0("Sample_", sprintf("%02d", seq_len(n_samples)))

# Distribución de loci por cromosoma
loci_per_chr <- ceiling(n_loci / n_chromosomes)
chr_names    <- paste0("chr", seq_len(n_chromosomes))

# Generar tabla de loci con coordenadas genómicas
loci_df <- data.frame(
  locus_idx   = seq_len(n_loci),
  og_id       = sprintf("OG%07d", seq_len(n_loci)),
  gene_id     = sprintf("Soltor1_%02dg%06d",
                         rep(seq_len(n_chromosomes), each = loci_per_chr)[seq_len(n_loci)],
                         rep(seq(10, by = 10, length.out = loci_per_chr), n_chromosomes)[seq_len(n_loci)]),
  chr         = rep(chr_names, each = loci_per_chr)[seq_len(n_loci)],
  stringsAsFactors = FALSE
)
loci_df$start  <- (loci_df$locus_idx - 1) * 35000 + 10000
loci_df$end    <- loci_df$start + 4500
loci_df$strand <- sample(c("+", "-"), n_loci, replace = TRUE)

# --- Diseño de la Matriz Biológica de CNV (50 Loci x 10 Muestras) ---
# Matriz estructurada para abarcar desde 0 hasta 9 copias con heterogeneidad entre muestras:
target_cnv <- matrix(0L, nrow = n_loci, ncol = n_samples)
colnames(target_cnv) <- samples
rownames(target_cnv) <- loci_df$og_id

# 1. Loci 1-5: Familias altamente expandidas (4 a 9 copias)
target_cnv[1, ] <- c(6, 7, 6, 8, 7, 6, 9, 7, 8, 6)
target_cnv[2, ] <- c(5, 6, 5, 5, 6, 7, 5, 6, 7, 6)
target_cnv[3, ] <- c(7, 8, 6, 7, 9, 8, 7, 8, 6, 7)
target_cnv[4, ] <- c(4, 5, 6, 4, 5, 6, 5, 6, 4, 5)
target_cnv[5, ] <- c(8, 7, 9, 8, 7, 8, 6, 7, 8, 7)

# 2. Loci 6-15: Expansiones moderadas (2 a 5 copias)
target_cnv[6,  ] <- c(3, 4, 3, 2, 4, 3, 5, 4, 3, 2)
target_cnv[7,  ] <- c(2, 3, 2, 3, 2, 4, 3, 2, 3, 4)
target_cnv[8,  ] <- c(4, 3, 4, 5, 3, 4, 3, 4, 5, 3)
target_cnv[9,  ] <- c(3, 2, 3, 3, 4, 2, 3, 4, 2, 3)
target_cnv[10, ] <- c(2, 2, 3, 2, 2, 3, 4, 2, 3, 2)
target_cnv[11, ] <- c(3, 3, 2, 4, 3, 2, 3, 3, 4, 3)
target_cnv[12, ] <- c(4, 5, 4, 3, 4, 5, 3, 4, 5, 4)
target_cnv[13, ] <- c(2, 3, 4, 2, 3, 2, 4, 3, 2, 3)
target_cnv[14, ] <- c(3, 2, 2, 3, 4, 3, 2, 4, 3, 2)
target_cnv[15, ] <- c(4, 3, 5, 4, 3, 4, 3, 2, 4, 3)

# 3. Loci 16-25: Duplicaciones menores y transiciones (1 a 3 copias)
target_cnv[16, ] <- c(1, 2, 2, 1, 2, 3, 1, 2, 1, 2)
target_cnv[17, ] <- c(2, 1, 1, 2, 2, 1, 3, 2, 1, 1)
target_cnv[18, ] <- c(1, 1, 2, 1, 2, 2, 1, 1, 2, 2)
target_cnv[19, ] <- c(2, 2, 1, 3, 2, 1, 2, 2, 1, 2)
target_cnv[20, ] <- c(3, 2, 2, 1, 2, 2, 1, 3, 2, 1)
target_cnv[21, ] <- c(1, 2, 1, 2, 1, 2, 2, 1, 2, 1)
target_cnv[22, ] <- c(2, 1, 2, 2, 1, 1, 2, 2, 1, 2)
target_cnv[23, ] <- c(1, 1, 2, 1, 2, 2, 1, 2, 2, 1)
target_cnv[24, ] <- c(2, 2, 1, 1, 2, 1, 2, 1, 1, 2)
target_cnv[25, ] <- c(1, 2, 2, 2, 1, 3, 2, 1, 2, 1)

# 4. Loci 26-38: Loci basales de copia única (1 copia mayoritaria, raras 2)
for (li in 26:38) {
  target_cnv[li, ] <- sample(c(1, 1, 1, 1, 1, 1, 1, 2, 1, 2), n_samples)
}

# 5. Loci 39-50: Presencia/ausencia y deleciones marcadas (0, 1 o 2 copias)
target_cnv[39, ] <- c(0, 1, 0, 1, 1, 0, 1, 0, 1, 0)
target_cnv[40, ] <- c(1, 0, 1, 0, 0, 1, 0, 1, 1, 0)
target_cnv[41, ] <- c(0, 0, 1, 0, 1, 0, 0, 1, 0, 1)
target_cnv[42, ] <- c(0, 0, 0, 1, 0, 0, 1, 0, 0, 1)
target_cnv[43, ] <- c(1, 1, 0, 0, 1, 0, 1, 1, 0, 0)
target_cnv[44, ] <- c(0, 1, 0, 1, 0, 1, 0, 0, 1, 0)
target_cnv[45, ] <- c(2, 0, 1, 0, 2, 0, 1, 0, 1, 0)
target_cnv[46, ] <- c(0, 0, 0, 1, 0, 0, 1, 0, 0, 1)
target_cnv[47, ] <- c(1, 0, 1, 0, 0, 0, 1, 0, 1, 0)
target_cnv[48, ] <- c(0, 1, 0, 1, 0, 1, 0, 1, 0, 0)
target_cnv[49, ] <- c(2, 0, 1, 0, 1, 0, 0, 1, 0, 1)
target_cnv[50, ] <- c(0, 0, 1, 0, 0, 1, 0, 0, 1, 0)

cat(sprintf("Matriz de CNV objetivo generada: %d loci x %d muestras\n", nrow(target_cnv), ncol(target_cnv)))
cat("Desglose de valores de copias esperados en la matriz:\n")
print(table(as.vector(target_cnv)))

# --- 1. Generar mock_reference.gff3 ---
cat("Generando mock_reference.gff3 ...\n")
gff3_lines <- c("##gff-version 3")
for (i in seq_len(nrow(loci_df))) {
  gff3_lines <- c(gff3_lines, sprintf(
    "%s\tSimulated\tgene\t%d\t%d\t.\t%s\t.\tID=gene:%s;Name=%s;",
    loci_df$chr[i], loci_df$start[i], loci_df$end[i],
    loci_df$strand[i], loci_df$gene_id[i], loci_df$gene_id[i]
  ))
}
writeLines(gff3_lines, file.path(out_dir, "mock_reference.gff3"))
cat(sprintf("  -> %d genes en %d cromosomas\n", n_loci, n_chromosomes))

# --- 2. Generar mock_captus.gff (contigs de CAPTUS) ---
cat("Generando mock_captus.gff ...\n")

captus_gff_lines <- c("##gff-version 3")
captus_id_counter <- 0

for (i in seq_len(nrow(loci_df))) {
  og   <- loci_df$og_id[i]
  gene <- loci_df$gene_id[i]
  chr  <- loci_df$chr[i]

  max_k <- max(target_cnv[i, ])
  # Incluir en el GFF contigs para todas las copias posibles y 2 contigs adicionales de ruido
  n_gff_hits <- max(max_k + 2, 3)

  for (h in seq_len(n_gff_hits)) {
    captus_id_counter <- captus_id_counter + 1
    hit_idx <- sprintf("%02d", h - 1)

    if (h == 1) {
      # Hit principal (superposición directa con el locus de referencia)
      h_start   <- loci_df$start[i] + 100
      h_end     <- loci_df$end[i] - 150
      score_gff <- 1.00
    } else if (h <= max_k) {
      # Parálogos funcionales confirmados en tándem
      h_start   <- loci_df$end[i] + (h - 1) * 3500 + sample(100:400, 1)
      h_end     <- h_start + 4000
      score_gff <- round(runif(1, 0.88, 0.98), 2)
    } else {
      # Contig de ruido intergénico o hit distal no ortólogo
      h_start   <- loci_df$end[i] + sample(15000:100000, 1)
      h_end     <- h_start + sample(800:1800, 1)
      score_gff <- round(runif(1, 0.20, 0.50), 2)
    }

    name   <- paste0(og, "__", hit_idx)
    target <- sprintf("%s.t1-%s 1-100 101-400", gene, og)
    strand <- loci_df$strand[i]

    captus_gff_lines <- c(captus_gff_lines, sprintf(
      "%s\tScipio:Captus\tprotein_match:NUC\t%d\t%d\t%.2f\t%s\t.\tID=NUC_%04d;Name=%s;Target=%s",
      chr, h_start, h_end, score_gff, strand, captus_id_counter, name, target
    ))
  }
}
writeLines(captus_gff_lines, file.path(out_dir, "mock_captus.gff"))
cat(sprintf("  -> %d contigs protein_match:NUC en GFF\n", captus_id_counter))

# --- 3. Generar mock_dictionary.tsv ---
cat("Generando mock_dictionary.tsv ...\n")
# Anotar los primeros 45 loci en el diccionario de homología
n_known <- 45
dict_df <- data.frame(
  Locus_Benoit = loci_df$gene_id[seq_len(n_known)],
  OG_Benoit    = loci_df$og_id[seq_len(n_known)],
  stringsAsFactors = FALSE
)
write.table(dict_df, file = file.path(out_dir, "mock_dictionary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  -> %d loci anotados en el diccionario\n", n_known))

# --- 4. Generar mock_captus_extract.tsv mediante simulación paramétrica por clase ---
cat("Generando mock_captus_extract.tsv ...\n")

# Simula un vector de métricas continuas de alineamiento según la clase biológica del contig
simulate_contig_metrics <- function(is_functional = TRUE) {
  # Factor latente común por secuencia para reflejar la colinealidad natural entre cobertura, identidad y puntaje
  z_shared <- rnorm(1)
  u_rec <- pnorm(0.75 * z_shared + sqrt(1 - 0.75^2) * rnorm(1))
  u_id  <- pnorm(0.75 * z_shared + sqrt(1 - 0.75^2) * rnorm(1))
  u_sc  <- pnorm(0.80 * z_shared + sqrt(1 - 0.80^2) * rnorm(1))

  if (isTRUE(is_functional)) {
    pct_rec <- round(24 + (99.5 - 24) * qbeta(u_rec, 9.2, 2.0), 2)
    pct_id  <- round(52 + (99.5 - 52) * qbeta(u_id,  9.5, 2.0), 2)
    sc      <- round(0.50 + (0.99 - 0.50) * qbeta(u_sc, 7.5, 1.8), 3)
    wsc     <- round(0.45 + (0.98 - 0.45) * qbeta(u_sc, 7.5, 1.8), 3)
  } else {
    pct_rec <- round(16 + (84 - 16) * qbeta(u_rec, 4.2, 3.2), 2)
    pct_id  <- round(38 + (85 - 38) * qbeta(u_id,  4.5, 3.2), 2)
    sc      <- round(0.14 + (0.66 - 0.14) * qbeta(u_sc, 3.8, 3.8), 3)
    wsc     <- round(0.09 + (0.54 - 0.09) * qbeta(u_sc, 3.8, 3.8), 3)
  }

  list(
    pct_recovered   = pct_rec,
    pct_identity    = pct_id,
    score           = sc,
    wscore          = wsc,
    ref_len_matched = as.integer(round(500 * (pct_rec / 100)))
  )
}

target_total_seqs  <- 1800
total_functional   <- sum(target_cnv)
total_noise_needed <- target_total_seqs - total_functional

noise_matrix <- matrix(1L, nrow = n_loci, ncol = n_samples)
extra_noise <- total_noise_needed - (n_loci * n_samples)
if (extra_noise > 0) {
  set.seed(42)
  extra_slots <- sample(n_loci * n_samples, extra_noise)
  noise_matrix[extra_slots] <- 2L
}

rows <- list()
row_counter <- 0L

for (s_idx in seq_along(samples)) {
  s <- samples[s_idx]

  for (i in seq_len(n_loci)) {
    og       <- loci_df$og_id[i]
    gene     <- loci_df$gene_id[i]
    target_k <- target_cnv[i, s_idx]

    # Contigs correspondientes a copias génicas funcionales
    if (target_k > 0) {
      for (h in seq_len(target_k)) {
        hit_num <- h - 1L
        row_counter <- row_counter + 1L

        m_sim      <- simulate_contig_metrics(is_functional = TRUE)
        hit_len    <- sample(2000:4500, 1)
        cds_len    <- sample(1000:2200, 1)
        intron_len <- sample(100:500, 1)
        flanks_len <- max(0L, hit_len - cds_len - intron_len)

        rows[[row_counter]] <- data.frame(
          sample_name     = s,
          marker_type     = "NUC",
          locus           = og,
          ref_name        = sprintf("%s.1-%s", gene, og),
          ref_coords      = "1-148,149-434",
          ref_type        = "prot",
          ref_len_matched = m_sim$ref_len_matched,
          hit             = hit_num,
          pct_recovered   = m_sim$pct_recovered,
          pct_identity    = m_sim$pct_identity,
          score           = m_sim$score,
          wscore          = m_sim$wscore,
          hit_len         = hit_len,
          cds_len         = cds_len,
          intron_len      = intron_len,
          flanks_len      = flanks_len,
          frameshifts     = "",
          hit_contigs     = 1L,
          hit_l50         = 1L,
          hit_l90         = 1L,
          hit_lg50        = 1L,
          hit_lg90        = 1L,
          ctg_names       = paste0("contig_", sprintf("%02d", sample(1:40, 1))),
          ctg_strands     = sample(c("+", "-"), 1),
          ctg_coords      = sprintf("%d-%d", loci_df$start[i], loci_df$end[i]),
          ctg_avg_depth   = round(rnorm(1, 48, 6), 1),
          stringsAsFactors = FALSE
        )
      }
    }

    # Contigs de ruido intergénico o hibridación cruzada no ortóloga
    n_noise_for_locus <- noise_matrix[i, s_idx]
    for (nh in seq_len(n_noise_for_locus)) {
      noise_hit_num <- target_k + nh - 1L
      row_counter <- row_counter + 1L

      m_sim      <- simulate_contig_metrics(is_functional = FALSE)
      hit_len    <- sample(1200:3000, 1)
      cds_len    <- sample(500:1200, 1)
      intron_len <- sample(100:400, 1)
      flanks_len <- max(0L, hit_len - cds_len - intron_len)

      rows[[row_counter]] <- data.frame(
        sample_name     = s,
        marker_type     = "NUC",
        locus           = og,
        ref_name        = sprintf("%s.1-%s", gene, og),
        ref_coords      = "1-148,149-434",
        ref_type        = "prot",
        ref_len_matched = m_sim$ref_len_matched,
        hit             = noise_hit_num,
        pct_recovered   = m_sim$pct_recovered,
        pct_identity    = m_sim$pct_identity,
        score           = m_sim$score,
        wscore          = m_sim$wscore,
        hit_len         = hit_len,
        cds_len         = cds_len,
        intron_len      = intron_len,
        flanks_len      = flanks_len,
        frameshifts     = "",
        hit_contigs     = 1L,
        hit_l50         = 1L,
        hit_l90         = 1L,
        hit_lg50        = 1L,
        hit_lg90        = 1L,
        ctg_names       = paste0("contig_", sprintf("%02d", sample(1:40, 1))),
        ctg_strands     = sample(c("+", "-"), 1),
        ctg_coords      = sprintf("%d-%d", loci_df$start[i] + 50000, loci_df$end[i] + 52000),
        ctg_avg_depth   = round(rnorm(1, 15, 5), 1),
        stringsAsFactors = FALSE
      )
    }
  }
}

df_tsv <- do.call(rbind, rows)
df_tsv <- df_tsv[order(df_tsv$sample_name, df_tsv$locus, df_tsv$hit), ]
rownames(df_tsv) <- NULL

write.table(df_tsv, file = file.path(out_dir, "mock_captus_extract.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("  -> %d filas totales NUC simuladas | %d muestras | %d loci\n",
            nrow(df_tsv), n_samples, n_loci))

# --- Resumen ---
cat("\n=== RESUMEN DE ARCHIVOS GENERADOS ===\n")
for (fname in c("mock_reference.gff3", "mock_captus.gff",
                "mock_captus_extract.tsv", "mock_dictionary.tsv")) {
  fpath <- file.path(out_dir, fname)
  finfo <- file.info(fpath)
  cat(sprintf("  %-30s  %6.1f KB\n", fname, finfo$size / 1024))
}
cat("\n[OK] Datos sintéticos generados correctamente para el MWE (rango de 0 a 9 copias por locus).\n")
