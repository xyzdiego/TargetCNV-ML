#!/usr/bin/env Rscript
# Genera los datos sintéticos de prueba (MWE) en data/test_data/:
# 10 muestras x 50 ortogrupos en 5 cromosomas con 0 a 9 copias por locus.
# Cada copia funcional en tándem se anota como un gen propio en mock_reference.gff3
# y en mock_dictionary.tsv, mientras que los contigs de ruido caen en regiones intergénicas.

set.seed(2026)

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
project_root <- if (length(file_arg) > 0) {
  sd <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
  if (basename(dirname(sd)) == "data") dirname(dirname(sd)) else if (basename(sd) == "data") dirname(sd) else sd
} else {
  getwd()
}
out_dir <- file.path(project_root, "data", "test_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Generando conjunto sintético en: %s\n", out_dir))

n_samples     <- 10L
n_loci        <- 50L
n_chromosomes <- 5L
samples       <- paste0("Sample_", sprintf("%02d", seq_len(n_samples)))
loci_per_chr  <- ceiling(n_loci / n_chromosomes)
chr_names     <- paste0("chr", seq_len(n_chromosomes))

# Espaciamiento amplio por locus (120 kb) para alojar hasta 9 parálogos en tándem sin solapar vecinos
loci_df <- data.frame(
  locus_idx = seq_len(n_loci),
  og_id     = sprintf("OG%07d", seq_len(n_loci)),
  gene_id   = sprintf(
    "Soltor1_%02dg%06d",
    rep(seq_len(n_chromosomes), each = loci_per_chr)[seq_len(n_loci)],
    rep(seq(10, by = 10, length.out = loci_per_chr), n_chromosomes)[seq_len(n_loci)]
  ),
  chr       = rep(chr_names, each = loci_per_chr)[seq_len(n_loci)],
  stringsAsFactors = FALSE
)
loci_df$slot_in_chr <- ((loci_df$locus_idx - 1L) %% loci_per_chr)
loci_df$start       <- loci_df$slot_in_chr * 120000L + 10000L
loci_df$end         <- loci_df$start + 4500L
loci_df$strand      <- sample(c("+", "-"), n_loci, replace = TRUE)

# Matriz de número de copias objetivo (50 ortogrupos x 10 muestras, rango 0 a 9)
target_cnv <- matrix(0L, nrow = n_loci, ncol = n_samples)
colnames(target_cnv) <- samples
rownames(target_cnv) <- loci_df$og_id

# Familias expandidas (4 a 9 copias)
target_cnv[1, ] <- c(6, 7, 6, 8, 7, 6, 9, 7, 8, 6)
target_cnv[2, ] <- c(5, 6, 5, 5, 6, 7, 5, 6, 7, 6)
target_cnv[3, ] <- c(7, 8, 6, 7, 9, 8, 7, 8, 6, 7)
target_cnv[4, ] <- c(4, 5, 6, 4, 5, 6, 5, 6, 4, 5)
target_cnv[5, ] <- c(8, 7, 9, 8, 7, 8, 6, 7, 8, 7)

# Expansiones moderadas (2 a 5 copias)
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

# Duplicaciones menores (1 a 3 copias)
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

# Copia única mayoritaria (1 o 2 copias)
for (li in 26:38) {
  target_cnv[li, ] <- sample(c(1, 1, 1, 1, 1, 1, 1, 2, 1, 2), n_samples)
}

# Presencia/ausencia (0, 1 o 2 copias)
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

# 1-3. Generar simultáneamente mock_reference.gff3, mock_captus.gff y mock_dictionary.tsv
# Cada copia en tándem (h = 1..max_k) tiene su propio gen anotado y su entrada en el diccionario;
# los hits de ruido (h = max_k+1 .. max_k+2) se ubican en espacio intergénico sin gen de referencia.
gff3_lines       <- c("##gff-version 3")
captus_gff_lines <- c("##gff-version 3")
dict_rows        <- list()
captus_id_counter <- 0L
dict_counter      <- 0L
max_k_vec         <- apply(target_cnv, 1, max)

for (i in seq_len(nrow(loci_df))) {
  og     <- loci_df$og_id[i]
  gene   <- loci_df$gene_id[i]
  chr    <- loci_df$chr[i]
  strand <- loci_df$strand[i]
  max_k  <- max(1L, as.integer(max_k_vec[i]))

  # Anotar el gen principal y sus parálogos en tándem
  for (h in seq_len(max_k)) {
    g_start <- loci_df$start[i] + (h - 1L) * 7000L
    g_end   <- g_start + 4500L
    g_id    <- if (h == 1L) gene else sprintf("%s_p%d", gene, h)

    gff3_lines <- c(gff3_lines, sprintf(
      "%s\tSimulated\tgene\t%d\t%d\t.\t%s\t.\tID=gene:%s;Name=%s;",
      chr, g_start, g_end, strand, g_id, g_id
    ))

    dict_counter <- dict_counter + 1L
    dict_rows[[dict_counter]] <- data.frame(
      Locus_Benoit = g_id,
      OG_Benoit    = og,
      stringsAsFactors = FALSE
    )

    captus_id_counter <- captus_id_counter + 1L
    hit_idx <- sprintf("%02d", h - 1L)
    h_start <- g_start + 120L
    h_end   <- g_end - 120L
    sc_gff  <- if (h == 1L) 1.00 else round(runif(1, 0.88, 0.98), 2)
    name    <- paste0(og, "__", hit_idx)
    target  <- sprintf("%s.t1-%s 1-100 101-400", g_id, og)

    captus_gff_lines <- c(captus_gff_lines, sprintf(
      "%s\tScipio:Captus\tprotein_match:NUC\t%d\t%d\t%.2f\t%s\t.\tID=NUC_%04d;Name=%s;Target=%s",
      chr, h_start, h_end, sc_gff, strand, captus_id_counter, name, target
    ))
  }

  # Dos contigs de ruido intergénico dedicados (índices max_k y max_k + 1) fuera de todos los genes
  for (nh in 1:2) {
    captus_id_counter <- captus_id_counter + 1L
    hit_idx <- sprintf("%02d", max_k + nh - 1L)
    h_start <- loci_df$start[i] + 85000L + (nh - 1L) * 6000L
    h_end   <- h_start + sample(800:1800, 1)
    sc_gff  <- round(runif(1, 0.20, 0.50), 2)
    name    <- paste0(og, "__", hit_idx)
    target  <- sprintf("%s.t1-%s 1-100 101-400", gene, og)

    captus_gff_lines <- c(captus_gff_lines, sprintf(
      "%s\tScipio:Captus\tprotein_match:NUC\t%d\t%d\t%.2f\t%s\t.\tID=NUC_%04d;Name=%s;Target=%s",
      chr, h_start, h_end, sc_gff, strand, captus_id_counter, name, target
    ))
  }
}

writeLines(gff3_lines, file.path(out_dir, "mock_reference.gff3"))
writeLines(captus_gff_lines, file.path(out_dir, "mock_captus.gff"))
dict_df <- do.call(rbind, dict_rows)
write.table(dict_df, file = file.path(out_dir, "mock_dictionary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# 4. Generar mock_captus_extract.tsv (simulación paramétrica por clase)
simulate_contig_metrics <- function(is_functional = TRUE) {
  z_shared <- rnorm(1)
  u_rec <- pnorm(0.75 * z_shared + sqrt(1 - 0.75^2) * rnorm(1))
  u_id  <- pnorm(0.75 * z_shared + sqrt(1 - 0.75^2) * rnorm(1))
  u_sc  <- pnorm(0.80 * z_shared + sqrt(1 - 0.80^2) * rnorm(1))

  if (isTRUE(is_functional)) {
    pct_rec <- round(52 + (99.5 - 52) * qbeta(u_rec, 8.5, 2.0), 2)
    pct_id  <- round(72 + (99.5 - 72) * qbeta(u_id,  9.0, 2.0), 2)
    sc      <- round(0.62 + (0.99 - 0.62) * qbeta(u_sc, 7.5, 1.8), 3)
    wsc     <- round(0.55 + (0.98 - 0.55) * qbeta(u_sc, 7.5, 1.8), 3)
  } else {
    pct_rec <- round(16 + (74 - 16) * qbeta(u_rec, 3.8, 3.5), 2)
    pct_id  <- round(38 + (78 - 38) * qbeta(u_id,  4.0, 3.5), 2)
    sc      <- round(0.14 + (0.58 - 0.14) * qbeta(u_sc, 3.5, 4.0), 3)
    wsc     <- round(0.09 + (0.48 - 0.09) * qbeta(u_sc, 3.5, 4.0), 3)
  }

  list(
    pct_recovered   = pct_rec,
    pct_identity    = pct_id,
    score           = sc,
    wscore          = wsc,
    ref_len_matched = as.integer(round(500 * (pct_rec / 100)))
  )
}

target_total_seqs  <- 1800L
total_functional   <- sum(target_cnv)
total_noise_needed <- target_total_seqs - total_functional

noise_matrix <- matrix(1L, nrow = n_loci, ncol = n_samples)
extra_noise <- total_noise_needed - (n_loci * n_samples)
if (extra_noise > 0L) {
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
    max_k    <- max(1L, as.integer(max_k_vec[i]))

    if (target_k > 0L) {
      for (h in seq_len(target_k)) {
        hit_num <- h - 1L
        row_counter <- row_counter + 1L

        m_sim      <- simulate_contig_metrics(is_functional = TRUE)
        hit_len    <- sample(2000:4500, 1)
        cds_len    <- sample(1000:2200, 1)
        intron_len <- sample(100:500, 1)
        flanks_len <- max(0L, hit_len - cds_len - intron_len)
        g_start    <- loci_df$start[i] + (h - 1L) * 7000L

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
          ctg_strands     = loci_df$strand[i],
          ctg_coords      = sprintf("%d-%d", g_start, g_start + 4500L),
          ctg_avg_depth   = round(rnorm(1, 48, 6), 1),
          stringsAsFactors = FALSE
        )
      }
    }

    # Los contigs de ruido usan siempre los índices reservados max_k y max_k + 1
    n_noise_for_locus <- noise_matrix[i, s_idx]
    for (nh in seq_len(n_noise_for_locus)) {
      noise_hit_num <- max_k + nh - 1L
      row_counter <- row_counter + 1L

      m_sim      <- simulate_contig_metrics(is_functional = FALSE)
      hit_len    <- sample(1200:3000, 1)
      cds_len    <- sample(500:1200, 1)
      intron_len <- sample(100:400, 1)
      flanks_len <- max(0L, hit_len - cds_len - intron_len)
      n_start    <- loci_df$start[i] + 85000L + (nh - 1L) * 6000L

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
        ctg_strands     = loci_df$strand[i],
        ctg_coords      = sprintf("%d-%d", n_start, n_start + 1500L),
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

cat(sprintf("Generados %d genes anotados, %d contigs GFF y %d filas NUC (%d funcionales, %d ruido).\n",
            nrow(dict_df), captus_id_counter, nrow(df_tsv), total_functional, nrow(df_tsv) - total_functional))
