\name{data_processing}
\alias{data_processing}
\alias{clean_captus_ids}
\alias{filter_gff_features}
\alias{harmonize_gff_seqlevels}
\alias{find_spatial_overlaps}
\alias{extract_target_loci_from_gff}
\alias{load_reference_dictionary}
\alias{classify_spatial_hits}
\alias{resolve_redundancies}
\alias{extract_cnv_matrix}
\alias{prepare_transfer_species_data}
\title{Spatial GFF labeling, homology dictionary lookup, and CNV counting}
\description{
Functions for standardizing CAPTUS contig identifiers, intersecting GFF3
annotations against reference genomes (chromosome or contig scale), reconciling
orthogroup dictionaries, removing intra-locus technical redundancies, and
building copy-number matrices.
}
\usage{
clean_captus_ids(df, ref_col = "ref_name", hit_col = "hit")

filter_gff_features(gff, feature_type = "gene")

harmonize_gff_seqlevels(
  query_gff,
  subject_gff,
  prefix_regex = "^[Ss]citrullifolium_|^[Ss]torvum_|^contig_"
)

find_spatial_overlaps(query_gff, subject_gff, ignore_strand = FALSE)

extract_target_loci_from_gff(query_gff)

load_reference_dictionary(dict_path, species = NULL)

classify_spatial_hits(
  df_cruces,
  dict_ref,
  locus_col = "Locus_Ref_Limpio",
  dict_locus_col = "Locus_Benoit",
  dict_og_col = "OG_Benoit"
)

resolve_redundancies(
  df,
  locus_col = "Locus_Ref_Limpio",
  class_col = "Validacion_Espacial",
  sample_col = if ("sample_name" \%in\% colnames(df)) "sample_name" else NULL
)

extract_cnv_matrix(
  df,
  sample_col = "Sample",
  locus_col = "Locus_Ortogrupo",
  cov_threshold = 50,
  verdict_col = "Veredicto_Consenso",
  ...
)

prepare_transfer_species_data(
  sample_name,
  genome_scale = c("chromosome", "contig"),
  captus_gff,
  captus_stats,
  dict_ref,
  gff_ref = NULL,
  save_rds = NULL
)
}
\arguments{
\item{df}{Data frame of CAPTUS extraction statistics or consensus verdicts.}
\item{ref_col}{Column name holding CAPTUS reference identifiers.}
\item{hit_col}{Column name holding CAPTUS hit indices.}
\item{gff}{\code{GRanges} object imported from a GFF/GFF3 file.}
\item{feature_type}{Feature type string to retain (e.g., \code{"gene"} or \code{"protein_match:NUC"}).}
\item{query_gff}{\code{GRanges} object of CAPTUS nuclear contigs.}
\item{subject_gff}{\code{GRanges} object of reference genome genes.}
\item{prefix_regex}{Regular expression stripped from sequence names when seqlevels do not overlap.}
\item{ignore_strand}{Logical; whether strand orientation is ignored in \code{findOverlaps}.}
\item{dict_path}{Path to an Excel (\code{.xlsx}) or delimited (\code{.tsv}/\code{.csv}) homology dictionary.}
\item{species}{Optional species identifier used to pick the matching dictionary sheet.}
\item{df_cruces}{Data frame produced by \code{find_spatial_overlaps()} or \code{extract_target_loci_from_gff()}.}
\item{dict_ref}{Dictionary data frame returned by \code{load_reference_dictionary()}.}
\item{locus_col}{Column name for reference loci or orthogroups.}
\item{dict_locus_col}{Gene locus column name inside \code{dict_ref}.}
\item{dict_og_col}{Orthogroup column name inside \code{dict_ref}.}
\item{class_col}{Column name holding ground-truth labels (\code{"Validacion_Espacial"}).}
\item{sample_col}{Column name identifying samples.}
\item{cov_threshold}{Minimum \code{pct_recovered} percentage to count additional homologous copies.}
\item{verdict_col}{Column name holding consensus verdicts.}
\item{...}{Additional ignored arguments kept for backward compatibility.}
\item{sample_name}{Target sample identifier inside \code{captus_stats}.}
\item{genome_scale}{Reference assembly continuity: \code{"chromosome"} or \code{"contig"}.}
\item{captus_gff}{Path to the CAPTUS contigs GFF file.}
\item{captus_stats}{Path to \code{captus-extract_stats.tsv}.}
\item{gff_ref}{Optional path to reference genome GFF3 (required when \code{genome_scale = "chromosome"}).}
\item{save_rds}{Optional output path to save the prepared data frame as RDS.}
}
\value{
Processed data frames, \code{GRanges} objects, or (for \code{extract_cnv_matrix()}) a list containing
\code{cnv_matrix}, \code{cnv_annotated}, and \code{functional_genes}.
}
