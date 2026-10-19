#!/usr/bin/env Rscript
# ==============================================================================
# generate_large_mock.R
# Genera datos mock realistas y de gran volumen para el MWE de TargetCNV-ML
# Incluye:
#   - 10 Muestras (MockSample_01 a MockSample_10)
#   - 50 Loci genomicos en 5 cromosomas
#   - Matriz de CNV altamente diversificada y contrastante:
#       * Loci con 0 copias (deleciones / perdida génica)
#       * Loci con 1 copia (copia unica)
#       * Loci con 2 a 5 copias (duplicaciones / triplicaciones moderadas)
#       * Loci con 6 a 9 copias (expansiones multigénicas en tándem)
#   - Separacion bimodal nitida en las distribuciones de metricas de ensamblaje
#     (pct_recovered, pct_identity, score, wscore) para graficos de densidad limpios
# ==============================================================================

set.seed(2026)

# Detectar raiz del proyecto
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

# --- Parametros de simulacion ---
n_samples     <- 10
n_loci        <- 50
n_chromosomes <- 5
samples       <- paste0("Sample_", sprintf("%02d", seq_len(n_samples)))

# Distribucion de loci por cromosoma
loci_per_chr <- ceiling(n_loci / n_chromosomes)
chr_names    <- paste0("chr", seq_len(n_chromosomes))

# Generar tabla de loci con coordenadas genomicas
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

# --- Diseno de la Matriz Biologica de CNV (50 Loci x 10 Muestras) ---
# Matriz estructurada para abarcar desde 0 hasta 9 copias con alta heterogeneidad:
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

# 4. Loci 26-38: Loci basales de copia unica (1 copia mayoritaria, raras 2)
for (li in 26:38) {
  target_cnv[li, ] <- sample(c(1, 1, 1, 1, 1, 1, 1, 2, 1, 2), n_samples)
}

# 5. Loci 39-50: Presencia/Ausencia y Deleciones marcadas (0, 1 o 2 copias)
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
cat("Desglose de valores de copias esperados en el heatmap:\n")
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

# --- 2. Generar mock_captus.gff (contigs de Captus) ---
cat("Generando mock_captus.gff ...\n")

captus_gff_lines <- c("##gff-version 3")
captus_id_counter <- 0
all_captus_entries <- list()

for (i in seq_len(nrow(loci_df))) {
  og   <- loci_df$og_id[i]
  gene <- loci_df$gene_id[i]
  chr  <- loci_df$chr[i]

  max_k <- max(target_cnv[i, ])
  # Asegurar que el GFF contenga contigs para todas las copias posibles + 2 contigs de ruido
  n_gff_hits <- max(max_k + 2, 3)

  for (h in seq_len(n_gff_hits)) {
    captus_id_counter <- captus_id_counter + 1
    hit_idx <- sprintf("%02d", h - 1)

    if (h == 1) {
      # Hit principal (superposicion exacta con el locus de referencia)
      h_start   <- loci_df$start[i] + 100
      h_end     <- loci_df$end[i] - 150
      score_gff <- 1.00
      hit_type  <- "primary"
    } else if (h <= max_k) {
      # Parálogos funcionales confirmados en tándem
      h_start   <- loci_df$end[i] + (h - 1) * 3500 + sample(100:400, 1)
      h_end     <- h_start + 4000
      score_gff <- round(runif(1, 0.88, 0.98), 2)
      hit_type  <- "functional_paralog"
    } else {
      # Contig de ruido intergénico o hit distal no ortólogo
      h_start   <- loci_df$end[i] + sample(15000:100000, 1)
      h_end     <- h_start + sample(800:1800, 1)
      score_gff <- round(runif(1, 0.20, 0.50), 2)
      hit_type  <- "distal_noise"
    }

    name   <- paste0(og, "__", hit_idx)
    target <- sprintf("%s.t1-%s 1-100 101-400", gene, og)
    strand <- loci_df$strand[i]

    captus_gff_lines <- c(captus_gff_lines, sprintf(
      "%s\tScipio:Captus\tprotein_match:NUC\t%d\t%d\t%.2f\t%s\t.\tID=NUC_%04d;Name=%s;Target=%s",
      chr, h_start, h_end, score_gff, strand, captus_id_counter, name, target
    ))

    all_captus_entries[[captus_id_counter]] <- data.frame(
      og         = og,
      gene       = gene,
      chr        = chr,
      hit        = h - 1,
      locus_idx  = i,
      is_primary = (h == 1),
      hit_type   = hit_type,
      start      = h_start,
      end        = h_end,
      stringsAsFactors = FALSE
    )
  }
}
writeLines(captus_gff_lines, file.path(out_dir, "mock_captus.gff"))
cat(sprintf("  -> %d contigs protein_match:NUC en GFF\n", captus_id_counter))

# --- 3. Generar mock_dictionary.tsv ---
cat("Generando mock_dictionary.tsv ...\n")
# Anotamos los primeros 45 loci formalmente en el diccionario
n_known <- 45
dict_df <- data.frame(
  Locus_Benoit = loci_df$gene_id[seq_len(n_known)],
  OG_Benoit    = loci_df$og_id[seq_len(n_known)],
  stringsAsFactors = FALSE
)
write.table(dict_df, file = file.path(out_dir, "mock_dictionary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  -> %d loci anotados en el diccionario\n", n_known))

# --- 4. Generar mock_captus_extract.tsv con separacion bimodal nitida ---
# --- 4. Generar mock_captus_extract.tsv con separacion bimodal y solape pronunciado ---
cat("Generando mock_captus_extract.tsv ...\n")

# Para satisfacer exactamente 1800 secuencias analizadas:
target_total_seqs  <- 1800
total_functional   <- sum(target_cnv)
total_noise_needed <- target_total_seqs - total_functional

# Distribuir total_noise_needed entre las 500 combinaciones (n_loci x n_samples)
noise_matrix <- matrix(1, nrow = n_loci, ncol = n_samples)
extra_noise <- total_noise_needed - (n_loci * n_samples)
if (extra_noise > 0) {
  set.seed(42)
  extra_slots <- sample(n_loci * n_samples, extra_noise)
  noise_matrix[extra_slots] <- 2
}

rows <- list()
row_counter <- 0

for (s_idx in seq_along(samples)) {
  s <- samples[s_idx]

  for (i in seq_len(n_loci)) {
    og       <- loci_df$og_id[i]
    gene     <- loci_df$gene_id[i]
    target_k <- target_cnv[i, s_idx]

    # Para K copias, generamos exactamente K hits funcionales
    if (target_k > 0) {
      for (h in seq_len(target_k)) {
        hit_num <- h - 1
        row_counter <- row_counter + 1

        # Metricas de ALTA CALIDAD (Poblacion Verde - Homologo Valido)
        pct_recovered <- pmin(100, pmax(72, rnorm(1, mean = 92, sd = 4.0)))
        pct_identity  <- pmin(100, pmax(88, rnorm(1, mean = 95.5, sd = 1.8)))
        score         <- pmin(1.0, pmax(0.82, rnorm(1, mean = 0.94, sd = 0.025)))
        wscore        <- pmin(1.0, pmax(0.80, score * ((pct_recovered / 100)^0.3) * ((pct_identity / 100)^0.7) + rnorm(1, 0, 0.01)))

        ref_total       <- sample(400:600, 1)
        ref_len_matched <- round(ref_total * (pct_recovered / 100))
        hit_len         <- sample(2000:4500, 1)
        cds_len         <- sample(1000:2200, 1)
        intron_len      <- sample(100:500, 1)
        flanks_len      <- max(0, hit_len - cds_len - intron_len)

        rows[[row_counter]] <- data.frame(
          sample_name     = s,
          marker_type     = "NUC",
          locus           = og,
          ref_name        = sprintf("%s.1-%s", gene, og),
          ref_coords      = "1-148,149-434",
          ref_type        = "prot",
          ref_len_matched = ref_len_matched,
          hit             = hit_num,
          pct_recovered   = round(pct_recovered, 2),
          pct_identity    = round(pct_identity, 2),
          score           = round(score, 3),
          wscore          = round(wscore, 3),
          hit_len         = hit_len,
          cds_len         = cds_len,
          intron_len      = intron_len,
          flanks_len      = flanks_len,
          frameshifts     = "",
          hit_contigs     = 1,
          hit_l50         = 1,
          hit_l90         = 1,
          hit_lg50        = 1,
          hit_lg90        = 1,
          ctg_names       = paste0("contig_", sprintf("%02d", sample(1:40, 1))),
          ctg_strands     = sample(c("+", "-"), 1),
          ctg_coords      = sprintf("%d-%d", loci_df$start[i], loci_df$end[i]),
          ctg_avg_depth   = round(rnorm(1, 48, 6), 1),
          is_functional   = TRUE,
          stringsAsFactors = FALSE
        )
      }
    }

    # Inyectar hits de RUIDO (Poblacion Roja - Ruido Descartado)
    n_noise_for_locus <- noise_matrix[i, s_idx]
    for (nh in seq_len(n_noise_for_locus)) {
      noise_hit_num <- target_k + nh - 1
      row_counter <- row_counter + 1

      pct_recovered <- pmin(55, pmax(12, rnorm(1, mean = 30, sd = 7.5)))
      pct_identity  <- pmin(78, pmax(42, rnorm(1, mean = 64, sd = 6.5)))
      score         <- pmin(0.55, pmax(0.18, rnorm(1, mean = 0.37, sd = 0.07)))
      wscore        <- pmin(0.42, pmax(0.08, score * ((pct_recovered / 100)^0.3) * ((pct_identity / 100)^0.7) + rnorm(1, 0, 0.01)))

      ref_total       <- sample(400:600, 1)
      ref_len_matched <- round(ref_total * (pct_recovered / 100))
      hit_len         <- sample(1200:3000, 1)
      cds_len         <- sample(500:1200, 1)
      intron_len      <- sample(100:400, 1)
      flanks_len      <- max(0, hit_len - cds_len - intron_len)

      rows[[row_counter]] <- data.frame(
        sample_name     = s,
        marker_type     = "NUC",
        locus           = og,
        ref_name        = sprintf("%s.1-%s", gene, og),
        ref_coords      = "1-148,149-434",
        ref_type        = "prot",
        ref_len_matched = ref_len_matched,
        hit             = noise_hit_num,
        pct_recovered   = round(pct_recovered, 2),
        pct_identity    = round(pct_identity, 2),
        score           = round(score, 3),
        wscore          = round(wscore, 3),
        hit_len         = hit_len,
        cds_len         = cds_len,
        intron_len      = intron_len,
        flanks_len      = flanks_len,
        frameshifts     = "",
        hit_contigs     = 1,
        hit_l50         = 1,
        hit_l90         = 1,
        hit_lg50        = 1,
        hit_lg90        = 1,
        ctg_names       = paste0("contig_", sprintf("%02d", sample(1:40, 1))),
        ctg_strands     = sample(c("+", "-"), 1),
        ctg_coords      = sprintf("%d-%d", loci_df$start[i] + 50000, loci_df$end[i] + 52000),
        ctg_avg_depth   = round(rnorm(1, 15, 5), 1),
        is_functional   = FALSE,
        stringsAsFactors = FALSE
      )
    }
  }
}

df_tsv <- do.call(rbind, rows)

# --- Generacion estratificada continua con distribuciones monomodales y colas suaves ---
# Requerimiento:
# 1. Distribuciones monomodales (1 pico unico por clase y metrica).
# 2. Solapamiento de densidades entre 0.016 y 0.020 (en torno a 0.02).
# 3. Colas suaves sin cortes bruscos en raincloud: degradé en la cola inferior de la distribucion verde,
#    y degradé en ambas colas de la distribucion roja (con corte ligeramente mas marcado en cola superior).
# 4. Totalizar exactamente 1800 secuencias analizadas.

tp_indices <- which(df_tsv$is_functional & df_tsv$marker_type == "NUC")
fp_indices <- which(!df_tsv$is_functional & df_tsv$marker_type == "NUC")

n_tp <- length(tp_indices)
n_fp <- length(fp_indices)

u_tp <- (seq_len(n_tp) - 0.5) / n_tp
u_fp <- (seq_len(n_fp) - 0.5) / n_fp

# 1. pct_recovered (solapamiento en y ~ 0.013, degradé inferior verde, ambas colas roja)
tp_rec <- round(24 + (99.5 - 24) * qbeta(u_tp, 9.2, 2.0), 2)
fp_rec <- round(16 + (84 - 16) * qbeta(u_fp, 4.2, 3.2), 2)

# 2. pct_identity (solapamiento en y ~ 0.011, degradé inferior verde, ambas colas roja)
tp_id <- round(52 + (99.5 - 52) * qbeta(u_tp, 9.5, 2.0), 2)
fp_id <- round(38 + (85 - 38) * qbeta(u_fp, 4.5, 3.2), 2)

# 3. score (monomodal continuo en [0.14, 0.99])
tp_sc <- round(0.50 + (0.99 - 0.50) * qbeta(u_tp, 7.5, 1.8), 3)
fp_sc <- round(0.14 + (0.66 - 0.14) * qbeta(u_fp, 3.8, 3.8), 3)

# 4. wscore (monomodal continuo)
tp_wsc <- round(0.45 + (0.98 - 0.45) * qbeta(u_tp, 7.5, 1.8), 3)
fp_wsc <- round(0.09 + (0.54 - 0.09) * qbeta(u_fp, 3.8, 3.8), 3)

# Copula gaussiana para correlacion realista entre variables
# Evita cortes en guillotina en el espacio multivariado mientras preserva marginales exactas
make_copula_ranks <- function(n, r = 0.72) {
  p <- 3
  Sigma <- matrix(r, p, p)
  diag(Sigma) <- 1.0
  L <- t(chol(Sigma))
  Z <- L %*% matrix(rnorm(n * p), nrow = p)
  apply(Z, 1, rank)
}

set.seed(42)
R_tp <- make_copula_ranks(n_tp, 0.72)
R_fp <- make_copula_ranks(n_fp, 0.72)

# Mapear a las filas preservando la heterogeneidad de muestras y loci
shuffled_tp_slots <- sample(seq_len(n_tp))
shuffled_fp_slots <- sample(seq_len(n_fp))

df_tsv$pct_recovered[tp_indices] <- tp_rec[R_tp[shuffled_tp_slots, 1]]
df_tsv$pct_identity[tp_indices]  <- tp_id[R_tp[shuffled_tp_slots, 2]]
df_tsv$score[tp_indices]         <- tp_sc[R_tp[shuffled_tp_slots, 3]]
df_tsv$wscore[tp_indices]        <- tp_wsc[R_tp[shuffled_tp_slots, 3]]
df_tsv$ref_len_matched[tp_indices] <- round(500 * (df_tsv$pct_recovered[tp_indices] / 100))

df_tsv$pct_recovered[fp_indices] <- fp_rec[R_fp[shuffled_fp_slots, 1]]
df_tsv$pct_identity[fp_indices]  <- fp_id[R_fp[shuffled_fp_slots, 2]]
df_tsv$score[fp_indices]         <- fp_sc[R_fp[shuffled_fp_slots, 3]]
df_tsv$wscore[fp_indices]        <- fp_wsc[R_fp[shuffled_fp_slots, 3]]
df_tsv$ref_len_matched[fp_indices] <- round(500 * (df_tsv$pct_recovered[fp_indices] / 100))

# Limpiar metadato auxiliar antes de exportar
df_tsv$is_functional <- NULL

df_tsv <- df_tsv[order(df_tsv$sample_name, df_tsv$locus, df_tsv$hit), ]
rownames(df_tsv) <- NULL

write.table(df_tsv, file = file.path(out_dir, "mock_captus_extract.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("  -> %d filas totales NUC analizadas | %d muestras | %d loci\n",
            nrow(df_tsv), n_samples, n_loci))

# --- Resumen ---
cat("\n=== RESUMEN DE ARCHIVOS GENERADOS ===\n")
for (fname in c("mock_reference.gff3", "mock_captus.gff",
                "mock_captus_extract.tsv", "mock_dictionary.tsv")) {
  fpath <- file.path(out_dir, fname)
  finfo <- file.info(fpath)
  cat(sprintf("  %-30s  %6.1f KB\n", fname, finfo$size / 1024))
}
cat("\n[OK] Datos mock generados exitosamente con separacion bimodal limpia y CNV variado (0 a 9).\n")
