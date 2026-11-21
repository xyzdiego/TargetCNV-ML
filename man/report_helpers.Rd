\name{report_helpers}
\alias{report_helpers}
\alias{read_tsv_safe}
\alias{load_sources_citations}
\alias{lookup_source_entry}
\alias{render_sources_box_html}
\alias{compute_global_fdr_score}
\alias{make_export_dt}
\alias{format_benchmark_dt}
\alias{get_top6_summary_models}
\alias{get_thesis_9_models}
\alias{build_interactive_benchmark_plotly}
\alias{load_readonly_training_artifacts}
\alias{render_top6_summary_cms_grid}
\alias{render_single_model_diagnostics}
\alias{render_single_model_extra}
\alias{build_all_model_tabs_markdown}
\alias{render_shap_waterfall_html}
\alias{render_species_cms_grid}
\alias{render_single_transfer_model}
\alias{build_species_model_tabs_markdown}
\alias{build_multispecies_comparison_plotly}
\alias{render_transfer_multispecies_tabs}
\title{Read-only artifact loaders and HTML helpers for Quarto reports}
\description{
Helper functions invoked by \code{training_report.qmd}, \code{report_template.qmd},
and \code{transferability_report.qmd} to load precomputed artifacts, build
interactive \code{DT} tables and \code{plotly} charts, embed the SHAP Waterfall
explorer, and emit dynamic Quarto tabsets without refitting models.
}
\usage{
read_tsv_safe(path_str)

load_sources_citations(sources_param = "", proj_root = ".")

lookup_source_entry(sources_df, sample_id = NULL, genome_scale = NULL)

render_sources_box_html(
  sources_df,
  sample_ids = NULL,
  default_scale = "chromosome",
  origin_files = "",
  box_class = "formula-box"
)

compute_global_fdr_score(df_row)

make_export_dt(df, caption_txt = NULL, page_len = 20)

format_benchmark_dt(
  df,
  other_df = NULL,
  current_scen = "Balanceado",
  other_scen = "Desbalanceado",
  caption_text = NULL
)

get_top6_summary_models(bench_df)

get_thesis_9_models(bench_df)

build_interactive_benchmark_plotly(
  df_bench,
  context_title = "Escenario Balanceado",
  view_label = "6 Modelos Clave"
)

load_readonly_training_artifacts(models_dir, data_file)

render_top6_summary_cms_grid(
  bench_df,
  all_cms,
  scen_txt = "Escenario Balanceado"
)

render_single_model_diagnostics(
  mod_name,
  models_list,
  bench_df,
  all_cms,
  df_ml_report = NULL,
  heuristic_model = NULL,
  youden_thresholds = NULL,
  scen_label = "Escenario Balanceado"
)

render_single_model_extra(
  mod_name,
  models_list,
  scen_label = "Escenario Balanceado"
)

build_all_model_tabs_markdown(
  scenario = c("balanced", "unbalanced"),
  tr_art
)

render_shap_waterfall_html(shap_obj, df_ml_report, proj_root = ".")

render_species_cms_grid(res_sp, mode = c("top6", "thesis9", "top3_f1"))

render_single_transfer_model(res_sp, mod_name)

build_species_model_tabs_markdown(sp_key, transfer_results)

build_multispecies_comparison_plotly(df_multiespecie_all)

render_transfer_multispecies_tabs(transfer_results, df_multiespecie_all = NULL)
}
\arguments{
\item{path_str}{Path to a TSV file.}
\item{sources_param}{Optional user-specified citation table path.}
\item{proj_root}{Project root directory.}
\item{sources_df}{Citations data frame loaded by \code{load_sources_citations()}.}
\item{sample_id}{Sample identifier or character vector.}
\item{genome_scale}{Fallback assembly scale (\code{"chromosome"} or \code{"contig"}).}
\item{sample_ids}{Sample identifiers present in the rendered dataset.}
\item{default_scale}{Fallback assembly scale string.}
\item{origin_files}{Display string of source file names.}
\item{box_class}{CSS class for the HTML container.}
\item{df_row}{Single-row benchmark data frame.}
\item{df}{Data frame to render.}
\item{caption_txt}{Optional table caption.}
\item{page_len}{Rows displayed per page.}
\item{other_df}{Optional comparison benchmark table from the alternative scenario.}
\item{current_scen}{Label for the active scenario.}
\item{other_scen}{Label for the alternative scenario.}
\item{caption_text}{Optional caption for formatted benchmark tables.}
\item{bench_df}{Benchmark data frame.}
\item{df_bench}{Benchmark data frame for interactive Plotly bars.}
\item{context_title}{Scenario or species title.}
\item{view_label}{Subset description label.}
\item{models_dir}{Directory holding serialized model artifacts.}
\item{data_file}{Path to prepared training dataset.}
\item{all_cms}{Named list of 2x2 confusion matrices.}
\item{scen_txt}{Scenario display string.}
\item{mod_name}{Model name.}
\item{models_list}{Named list of trained models.}
\item{df_ml_report}{Prepared modeling data frame.}
\item{heuristic_model}{Serialized weighted heuristic model list.}
\item{youden_thresholds}{Named numeric vector of Youden thresholds.}
\item{scen_label}{Readable scenario name.}
\item{scenario}{Scenario key (\code{"balanced"} or \code{"unbalanced"}).}
\item{tr_art}{Artifact bundle returned by \code{load_readonly_training_artifacts()}.}
\item{shap_obj}{Serialized \code{shapviz} object.}
\item{res_sp}{Single-species list entry from \code{transferability_results.rds}.}
\item{mode}{Summary grid selection mode (\code{"top6"}, \code{"thesis9"}, or \code{"top3_f1"}).}
\item{sp_key}{Species key inside \code{transfer_results}.}
\item{transfer_results}{Named list of evaluated species results.}
\item{df_multiespecie_all}{Combined multi-species benchmark data frame.}
}
\value{
Data frames, \code{DT::datatable} or \code{plotly} widgets, HTML strings, or invisible \code{NULL} for side-effect rendering helpers.
}
