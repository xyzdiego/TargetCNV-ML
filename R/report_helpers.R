# Funciones auxiliares de solo lectura y presentación interactiva para plantillas Quarto

#' Lee de forma segura un archivo TSV si existe en disco
#'
#' @param path_str Ruta al archivo `.tsv`.
#' @return Data frame o `NULL` si el archivo no existe.
#' @export
read_tsv_safe <- function(path_str) {
  if (!is.null(path_str) && nzchar(path_str) && file.exists(path_str)) {
    utils::read.delim(path_str, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    NULL
  }
}

#' Carga la tabla de fuentes genómicas y citas bibliográficas
#'
#' @param sources_param Ruta opcional proporcionada por parámetro.
#' @param proj_root Directorio raíz del proyecto.
#' @return Data frame con las citas bibliográficas o `NULL`.
#' @export
load_sources_citations <- function(sources_param = "", proj_root = ".") {
  cand <- c(
    if (!is.null(sources_param) && nzchar(trimws(sources_param))) c(sources_param, file.path(proj_root, sources_param)),
    file.path(proj_root, "data", "reference", "sources_citations.tsv"),
    file.path(proj_root, "sources_citations.tsv")
  )
  found <- cand[file.exists(cand)][1]
  if (is.na(found) || !nzchar(found)) return(NULL)
  tryCatch(
    utils::read.delim(found, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
}

#' Busca una entrada de citación por muestra o escala estructural del genoma
#'
#' @param sources_df Data frame cargado con `load_sources_citations()`.
#' @param sample_id Identificador o vector de muestras a buscar.
#' @param genome_scale Escala de respaldo (`"chromosome"` o `"contig"`).
#' @return Data frame de una fila o `NULL`.
#' @export
lookup_source_entry <- function(sources_df, sample_id = NULL, genome_scale = NULL) {
  if (is.null(sources_df) || nrow(sources_df) == 0 || !"sample_id" %in% colnames(sources_df)) return(NULL)
  if (!is.null(sample_id) && length(sample_id) > 0 && any(nzchar(sample_id))) {
    hit <- sources_df[tolower(trimws(sources_df$sample_id)) %in% tolower(trimws(sample_id)), , drop = FALSE]
    if (nrow(hit) > 0) return(hit[1, , drop = FALSE])
  }
  if (!is.null(genome_scale) && nzchar(genome_scale)) {
    def_key <- if (tolower(genome_scale) == "contig") "contig_default" else "chromosome_default"
    hit <- sources_df[tolower(trimws(sources_df$sample_id)) == def_key, , drop = FALSE]
    if (nrow(hit) > 0) return(hit[1, , drop = FALSE])
  }
  NULL
}

#' Genera el bloque HTML de trazabilidad de fuentes y archivos de entrada
#'
#' @param sources_df Data frame de citas bibliográficas.
#' @param sample_ids Identificadores de muestra presentes en el conjunto de datos.
#' @param default_scale Escala de respaldo (`"chromosome"` o `"contig"`).
#' @param origin_files Cadena con los nombres de archivo de origen.
#' @param box_class Clase CSS del contenedor (`"formula-box"` o `"rule-box"`).
#' @return Cadena HTML lista para imprimirse con `cat()`.
#' @export
render_sources_box_html <- function(sources_df,
                                    sample_ids = NULL,
                                    default_scale = "chromosome",
                                    origin_files = "",
                                    box_class = "formula-box") {
  hit <- lookup_source_entry(sources_df, sample_id = sample_ids, genome_scale = default_scale)
  if (!is.null(hit)) {
    paste0(
      sprintf('<div class="%s">', box_class),
      sprintf('<b>Fuente del Genoma de Referencia:</b> %s<br>', as.character(hit$genome_source[1])),
      sprintf('<b>Fuente de Homología y Datos Génicos:</b> %s<br>', as.character(hit$homology_source[1])),
      sprintf('<b>Archivos de Origen:</b> <code>%s</code>', origin_files),
      '</div>\n'
    )
  } else {
    sprintf('<div class="%s"><b>Archivos de Origen:</b> <code>%s</code></div>\n', box_class, origin_files)
  }
}

#' Calcula el índice compuesto con ponderación reforzada sobre el control de FDR
#'
#' @param df_row Fila de una tabla de benchmark con columnas `FDR`, `F1_Score`, `MCC`, `Exactitud_Balanceada` y `ROC_AUC`.
#' @return Puntaje escalar en `[0, 1]`.
#' @export
compute_global_fdr_score <- function(df_row) {
  fdr_v  <- suppressWarnings(as.numeric(df_row$FDR[1]))
  f1_v   <- suppressWarnings(as.numeric(df_row$F1_Score[1]))
  mcc_v  <- suppressWarnings(as.numeric(df_row$MCC[1]))
  bacc_v <- suppressWarnings(as.numeric(df_row$Exactitud_Balanceada[1]))
  auc_v  <- suppressWarnings(as.numeric(df_row$ROC_AUC[1]))

  if (length(fdr_v) == 0  || is.na(fdr_v))  fdr_v  <- 1.0
  if (length(f1_v) == 0   || is.na(f1_v))   f1_v   <- 0.0
  if (length(mcc_v) == 0  || is.na(mcc_v))  mcc_v  <- 0.0
  if (length(bacc_v) == 0 || is.na(bacc_v)) bacc_v <- 0.0
  if (length(auc_v) == 0  || is.na(auc_v))  auc_v  <- 0.0

  0.45 * (1 - fdr_v) + 0.20 * f1_v + 0.15 * mcc_v + 0.10 * bacc_v + 0.10 * auc_v
}

#' Construye una tabla DataTables exportable estándar
#'
#' @param df Data frame a visualizar.
#' @param caption_txt Leyenda opcional de la tabla.
#' @param page_len Número de filas por página.
#' @return Objeto `DT::datatable` o `NULL`.
#' @export
make_export_dt <- function(df, caption_txt = NULL, page_len = 20) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  DT::datatable(
    df,
    rownames = FALSE,
    caption = caption_txt,
    extensions = "Buttons",
    options = list(
      pageLength = page_len,
      lengthMenu = c(10, 20, 50, 100),
      scrollX = TRUE,
      dom = "Bfrtip",
      buttons = list(
        list(extend = "copy",  text = "Copiar Tabla",   exportOptions = list(modifier = list(page = "all"))),
        list(extend = "csv",   text = "Exportar CSV",   exportOptions = list(modifier = list(page = "all"))),
        list(extend = "excel", text = "Exportar Excel", exportOptions = list(modifier = list(page = "all")))
      ),
      language = list(
        search = "Filtrar:",
        zeroRecords = "No se encontraron registros",
        emptyTable = "Sin datos disponibles"
      )
    )
  )
}

#' Formatea una tabla de benchmark con comparación cruzada entre escenarios (Balanceado vs. Desbalanceado)
#'
#' @param df Data frame de benchmark del escenario activo.
#' @param other_df Data frame opcional del escenario alternativo para señalar mejoras con `*`.
#' @param current_scen Etiqueta del escenario actual.
#' @param other_scen Etiqueta del escenario alternativo.
#' @param caption_text Leyenda descriptiva de la tabla.
#' @return Objeto `DT::datatable`.
#' @export
format_benchmark_dt <- function(df,
                                other_df = NULL,
                                current_scen = "Balanceado",
                                other_scen = "Desbalanceado",
                                caption_text = NULL) {
  if (is.null(df) || nrow(df) == 0) return(NULL)

  fmt_fdr_pct <- function(val) {
    if (length(val) == 0 || is.na(val)) "Sin pred. TP (NA)" else sprintf("%.2f%%", as.numeric(val) * 100)
  }

  df_display  <- df
  base_models <- as.character(df_display$Modelo)
  comp_notes  <- rep("Invariante (regla fija)", nrow(df_display))
  flagged_vec <- rep(FALSE, nrow(df_display))

  has_other <- !is.null(other_df) && nrow(other_df) > 0 && "Modelo" %in% colnames(other_df)
  other_models <- if (has_other) as.character(other_df$Modelo) else character(0)

  if (has_other) {
    for (i in seq_len(nrow(df_display))) {
      m_name <- base_models[i]
      if (!is.na(m_name) && !isTRUE(grepl("^Youden_", m_name)) && m_name != "Weighted_Heuristic" && m_name %in% other_models) {
        row_curr  <- df[i, , drop = FALSE]
        row_other <- other_df[which(other_models == m_name)[1], , drop = FALSE]

        score_curr  <- compute_global_fdr_score(row_curr)
        score_other <- compute_global_fdr_score(row_other)
        fdr_curr    <- suppressWarnings(as.numeric(row_curr$FDR[1]))
        fdr_other   <- suppressWarnings(as.numeric(row_other$FDR[1]))

        curr_fdr_na  <- (length(fdr_curr) == 0 || is.na(fdr_curr))
        other_fdr_na <- (length(fdr_other) == 0 || is.na(fdr_other))

        other_is_better <- FALSE
        curr_is_better  <- FALSE

        if (curr_fdr_na && !other_fdr_na) {
          other_is_better <- TRUE
        } else if (!curr_fdr_na && other_fdr_na) {
          curr_is_better <- TRUE
        } else if (!curr_fdr_na && !other_fdr_na) {
          other_is_better <- isTRUE(
            (score_other > score_curr + 1e-4 && fdr_other <= fdr_curr + 0.005) ||
              (fdr_other < fdr_curr - 1e-4 && score_other >= score_curr - 0.002)
          )
          curr_is_better <- !other_is_better && isTRUE(
            score_curr > score_other + 1e-4 || fdr_curr < fdr_other - 1e-4
          )
        } else {
          other_is_better <- isTRUE(score_other > score_curr + 1e-4)
          curr_is_better  <- isTRUE(score_curr > score_other + 1e-4)
        }

        if (isTRUE(other_is_better)) {
          df_display$Modelo[i] <- paste0(m_name, " *")
          flagged_vec[i] <- TRUE
          comp_notes[i] <- sprintf(
            "* Mejor en %s (FDR: %s vs %s | Índice: %.3f vs %.3f)",
            other_scen, fmt_fdr_pct(fdr_other), fmt_fdr_pct(fdr_curr), score_other, score_curr
          )
        } else if (isTRUE(curr_is_better)) {
          comp_notes[i] <- sprintf(
            "Óptimo en %s (FDR: %s vs %s | Índice: %.3f vs %.3f)",
            current_scen, fmt_fdr_pct(fdr_curr), fmt_fdr_pct(fdr_other), score_curr, score_other
          )
        } else {
          comp_notes[i] <- sprintf("Equivalente en ambos (FDR: %s)", fmt_fdr_pct(fdr_curr))
        }
      }
    }
  }

  df_display <- df_display %>%
    dplyr::mutate(Mejor_Version_FDR = comp_notes, .after = "Modelo")

  pct_cols <- intersect(
    c("Exactitud", "Exactitud_Balanceada", "Sensibilidad", "Especificidad", "Precision", "VPN", "FDR"),
    colnames(df_display)
  )
  dec_cols <- intersect(
    c("F1_Score", "MCC", "Kappa", "ROC_AUC", "Umbral", "Indice_Youden_J"),
    colnames(df_display)
  )

  for (col in pct_cols) df_display[[col]] <- round(as.numeric(df_display[[col]]) * 100, 2)
  for (col in dec_cols) df_display[[col]] <- round(as.numeric(df_display[[col]]), 3)

  dt <- make_export_dt(df_display, caption_txt = caption_text, page_len = 50)
  if (length(pct_cols) > 0) dt <- dt %>% DT::formatCurrency(columns = pct_cols, currency = "%", before = FALSE, digits = 2)
  if (length(dec_cols) > 0) dt <- dt %>% DT::formatRound(columns = dec_cols, digits = 3)

  if (any(flagged_vec)) {
    flagged_names <- df_display$Modelo[flagged_vec]
    dt <- dt %>% DT::formatStyle(
      "Modelo",
      color = DT::styleEqual(flagged_names, rep("#b9770e", length(flagged_names))),
      backgroundColor = DT::styleEqual(flagged_names, rep("#fef9e7", length(flagged_names))),
      fontWeight = DT::styleEqual(flagged_names, rep("bold", length(flagged_names)))
    )
  }
  if ("F1_Score" %in% colnames(df_display)) {
    dt <- dt %>% DT::formatStyle(
      "F1_Score",
      background = DT::styleColorBar(c(0, 1), "#a9cce3"),
      backgroundSize = "95% 75%", backgroundRepeat = "no-repeat", backgroundPosition = "center", fontWeight = "bold"
    )
  }
  if ("ROC_AUC" %in% colnames(df_display)) {
    dt <- dt %>% DT::formatStyle(
      "ROC_AUC",
      background = DT::styleColorBar(c(0, 1), "#a3e4d7"),
      backgroundSize = "95% 75%", backgroundRepeat = "no-repeat", backgroundPosition = "center"
    )
  }
  dt
}

#' Selecciona los 6 métodos clave de síntesis (Top 2 Youden + Heurístico + Top 3 ML)
#'
#' @param bench_df Tabla de benchmark con columnas `Modelo`, `ROC_AUC` y `F1_Score`.
#' @return Vector de nombres de modelos ordenados por `ROC_AUC` descendente.
#' @export
get_top6_summary_models <- function(bench_df) {
  if (is.null(bench_df) || nrow(bench_df) == 0) return(character(0))

  top2_youden <- bench_df %>%
    dplyr::filter(grepl("^Youden_", Modelo)) %>%
    dplyr::arrange(dplyr::desc(ROC_AUC), dplyr::desc(F1_Score)) %>%
    dplyr::slice_head(n = 2) %>%
    dplyr::pull(Modelo)

  heur_name <- intersect("Weighted_Heuristic", bench_df$Modelo)

  top3_ml <- bench_df %>%
    dplyr::filter(!grepl("^Youden_", Modelo) & !Modelo %in% c("Weighted_Heuristic", "Ensamble_Stacking", "Consenso_7_Modelos")) %>%
    dplyr::arrange(dplyr::desc(ROC_AUC), dplyr::desc(F1_Score)) %>%
    dplyr::slice_head(n = 3) %>%
    dplyr::pull(Modelo)

  selected_6 <- unique(c(top2_youden, heur_name, top3_ml))
  if (length(selected_6) < 6) {
    extra_models <- bench_df %>%
      dplyr::filter(!Modelo %in% c(selected_6, "Ensamble_Stacking", "Consenso_7_Modelos")) %>%
      dplyr::arrange(dplyr::desc(ROC_AUC), dplyr::desc(F1_Score)) %>%
      dplyr::pull(Modelo)
    selected_6 <- unique(c(selected_6, extra_models))[seq_len(min(6L, length(unique(c(selected_6, extra_models)))))]
  }

  bench_df %>%
    dplyr::filter(Modelo %in% selected_6) %>%
    dplyr::arrange(dplyr::desc(ROC_AUC), dplyr::desc(F1_Score)) %>%
    dplyr::pull(Modelo)
}

#' Selecciona el subconjunto estándar de 9 modelos representativos para contraste inter-especie
#'
#' @param bench_df Tabla de benchmark con columna `Modelo`.
#' @return Vector con hasta 9 modelos ordenados por `ROC_AUC` descendente.
#' @export
get_thesis_9_models <- function(bench_df) {
  if (is.null(bench_df) || nrow(bench_df) == 0) return(character(0))
  canonical_9 <- c(
    "Youden_pct_recovered", "Youden_pct_identity", "Weighted_Heuristic",
    "ElasticNet", "NaiveBayes", "SVM_Polinomial", "MLP", "LDA", "Ensamble_Stacking"
  )
  present_9 <- intersect(canonical_9, bench_df$Modelo)
  bench_df %>%
    dplyr::filter(Modelo %in% present_9) %>%
    dplyr::arrange(dplyr::desc(ROC_AUC), dplyr::desc(F1_Score)) %>%
    dplyr::pull(Modelo)
}

#' Construye el gráfico de barras horizontales interactivo Plotly de benchmark multimétrica
#'
#' @param df_bench Tabla de benchmark.
#' @param context_title Etiqueta del escenario o especie evaluada.
#' @param view_label Subtítulo corto del subconjunto mostrado.
#' @return Objeto `plotly`.
#' @export
build_interactive_benchmark_plotly <- function(df_bench,
                                               context_title = "Escenario Balanceado",
                                               view_label = "6 Modelos Clave") {
  if (is.null(df_bench) || nrow(df_bench) == 0) return(NULL)

  metrics_info <- list(
    list(col = "ROC_AUC",              label = "ROC-AUC",              color = "#27ae60"),
    list(col = "F1_Score",             label = "F1-Score",             color = "#2980b9"),
    list(col = "Exactitud_Balanceada", label = "Exactitud Balanceada", color = "#8e44ad"),
    list(col = "Exactitud",            label = "Exactitud (Accuracy)", color = "#16a085"),
    list(col = "Sensibilidad",         label = "Sensibilidad (TPR)",   color = "#d35400"),
    list(col = "Especificidad",        label = "Especificidad (TNR)",  color = "#2c3e50"),
    list(col = "Precision",            label = "Precisión (PPV)",      color = "#f39c12"),
    list(col = "FDR",                  label = "FDR (Falsos Pos.)",    color = "#c0392b"),
    list(col = "MCC",                  label = "MCC (Matthews)",       color = "#7f8c8d")
  )

  avail_info <- Filter(function(m) m$col %in% colnames(df_bench), metrics_info)
  if (length(avail_info) == 0) return(NULL)

  df_best_to_worst <- df_bench %>%
    dplyr::arrange(dplyr::desc(ROC_AUC), dplyr::desc(F1_Score))
  y_levels_bottom_to_top <- rev(unique(as.character(df_best_to_worst$Modelo)))
  df_sorted <- df_best_to_worst %>%
    dplyr::mutate(Modelo = factor(as.character(Modelo), levels = y_levels_bottom_to_top))

  n_models <- nrow(df_sorted)
  dynamic_height <- if (n_models <= 6) 560 else if (n_models <= 10) 720 else max(950, n_models * 82)

  fig <- plotly::plot_ly(height = dynamic_height)
  n_m <- length(avail_info)

  for (m in avail_info) {
    vals <- round(as.numeric(df_sorted[[m$col]]), 3)
    fig <- fig %>%
      plotly::add_trace(
        y = df_sorted$Modelo, x = vals, name = m$label, type = "bar", orientation = "h",
        marker = list(color = m$color), text = sprintf("%.3f", vals), textposition = "auto",
        hovertemplate = paste0("<b>%{y}</b><br>", m$label, ": %{x:.3f}<extra></extra>")
      )
  }

  buttons_list <- list(
    list(method = "restyle", args = list("visible", as.list(rep(TRUE, n_m))), label = "Vista General (Todas las Métricas)"),
    list(
      method = "restyle",
      args = list("visible", as.list(sapply(avail_info, function(x) x$col %in% c("F1_Score", "ROC_AUC", "Exactitud", "Precision", "Sensibilidad", "Especificidad", "FDR")))),
      label = "7 Métricas Clásicas"
    )
  )

  for (i in seq_len(n_m)) {
    vis_vec <- rep(FALSE, n_m)
    vis_vec[i] <- TRUE
    buttons_list[[length(buttons_list) + 1]] <- list(
      method = "restyle", args = list("visible", as.list(vis_vec)), label = paste("Solo:", avail_info[[i]]$label)
    )
  }

  fig %>%
    plotly::layout(
      title = list(
        text = paste0(
          "<b>Benchmark (", view_label, " - Ordenado por ROC-AUC): ", context_title, "</b><br>",
          "<sup>Seleccione métricas en el menú superior o active/desactive series desde la leyenda inferior</sup>"
        ),
        x = 0.5
      ),
      barmode = "group", bargap = 0.24, bargroupgap = 0.08,
      xaxis = list(title = "Puntaje (0 a 1)", range = c(0, 1.08), tickfont = list(size = 12)),
      yaxis = list(title = "", tickfont = list(size = 12), categoryorder = "array", categoryarray = y_levels_bottom_to_top),
      legend = list(orientation = "h", y = -0.08, x = 0.5, xanchor = "center", font = list(size = 11)),
      margin = list(l = 195, r = 35, t = 95, b = 105),
      updatemenus = list(list(type = "dropdown", active = 0, x = 0.0, y = 1.08, xanchor = "left", yanchor = "top", buttons = buttons_list))
    )
}

#' Carga en modo de solo lectura todos los artefactos precomputados por `02_train_models.R`
#'
#' No ejecuta entrenamiento ni escribe archivos en `models_dir`. Si faltan artefactos
#' de modelos en `models_dir`, devuelve `has_models = FALSE` para mostrar un aviso claro.
#'
#' @param models_dir Directorio con artefactos `.rds` y `.tsv` de entrenamiento.
#' @param data_file Ruta al archivo de datos preparados (`data_prepared.rds` o `.tsv`).
#' @return Lista estructurada con datos, benchmarks, matrices de confusión, SHAP y Spearman.
#' @export
load_readonly_training_artifacts <- function(models_dir, data_file) {
  datos <- if (!is.null(data_file) && nzchar(data_file) && file.exists(data_file)) {
    if (grepl("\\.rds$", data_file, ignore.case = TRUE)) readRDS(data_file) else utils::read.delim(data_file, stringsAsFactors = FALSE)
  } else {
    NULL
  }

  df_ml_report <- if (!is.null(datos) && "Validacion_Espacial" %in% colnames(datos)) {
    pred_cols <- intersect(c("pct_recovered", "pct_identity", "score", "wscore"), colnames(datos))
    datos %>%
      dplyr::filter(!is.na(Validacion_Espacial)) %>%
      dplyr::mutate(
        Clase = factor(ifelse(Validacion_Espacial %in% c("TRUE_POSITIVE", "TP"), "TP", "FP"), levels = c("TP", "FP")),
        dplyr::across(dplyr::all_of(pred_cols), ~ as.numeric(as.character(.x)))
      ) %>%
      dplyr::filter(stats::complete.cases(dplyr::select(., dplyr::all_of(pred_cols)))) %>%
      dplyr::mutate(Fila_ID = dplyr::row_number()) %>%
      as.data.frame()
  } else {
    NULL
  }

  read_rds_opt <- function(fname) {
    fp <- file.path(models_dir, fname)
    if (file.exists(fp)) readRDS(fp) else NULL
  }

  models_bal <- read_rds_opt("models_balanced.rds")
  if (is.null(models_bal)) models_bal <- read_rds_opt("all_trained_models.rds")
  models_bal <- if (!is.null(models_bal)) normalize_model_names(models_bal) else list()

  models_unbal <- read_rds_opt("models_unbalanced.rds")
  models_unbal <- if (!is.null(models_unbal)) normalize_model_names(models_unbal) else list()

  youden_thresholds <- read_rds_opt("youden_thresholds.rds")
  heuristic_model   <- read_rds_opt("heuristic_model.rds")
  cm_list           <- read_rds_opt("confusion_matrices.rds")
  if (is.null(cm_list)) cm_list <- list()

  shap_obj_bal <- read_rds_opt("shap_values_balanced.rds")
  if (is.null(shap_obj_bal)) shap_obj_bal <- read_rds_opt("shap_values.rds")

  spearman_res <- read_rds_opt("spearman_collinearity.rds")
  if (is.null(spearman_res) && !is.null(df_ml_report)) {
    spearman_res <- tryCatch(compute_spearman_collinearity(df_ml_report), error = function(e) NULL)
  }

  load_bench <- function(scen_key) {
    b_rds <- read_rds_opt(sprintf("benchmark_models_%s.rds", scen_key))
    if (is.null(b_rds) && scen_key == "balanced") b_rds <- read_rds_opt("benchmark_models.rds")
    if (is.null(b_rds)) b_rds <- read_tsv_safe(file.path(models_dir, sprintf("benchmark_models_%s.tsv", scen_key)))
    if (!is.null(b_rds) && "Modelo" %in% colnames(b_rds)) {
      b_rds$Modelo <- normalize_model_names(as.character(b_rds$Modelo))
      b_rds <- b_rds %>% dplyr::arrange(dplyr::desc(ROC_AUC), dplyr::desc(F1_Score))
    }
    b_rds
  }

  bench_bal   <- load_bench("balanced")
  bench_unbal <- load_bench("unbalanced")

  assemble_cms <- function(scen_key, bench_df) {
    cms <- list()
    if (!is.null(cm_list$heuristic)) cms[["Weighted_Heuristic"]] <- cm_list$heuristic
    if (!is.null(cm_list[[scen_key]])) {
      sc_cms <- normalize_model_names(cm_list[[scen_key]])
      for (nm in names(sc_cms)) cms[[nm]] <- sc_cms[[nm]]
    }
    if (!is.null(cm_list$youden)) {
      for (yv in names(cm_list$youden)) cms[[paste0("Youden_", yv)]] <- cm_list$youden[[yv]]
    }
    if (!is.null(bench_df) && all(c("Modelo", "TP", "FP", "FN", "TN") %in% colnames(bench_df))) {
      for (i in seq_len(nrow(bench_df))) {
        mn <- as.character(bench_df$Modelo[i])
        if (is.null(cms[[mn]])) {
          mat <- matrix(
            c(bench_df$TP[i], bench_df$FN[i], bench_df$FP[i], bench_df$TN[i]),
            nrow = 2, byrow = FALSE,
            dimnames = list(Prediccion = c("TP", "FP"), Referencia = c("TP", "FP"))
          )
          cms[[mn]] <- as.table(mat)
        }
      }
    }
    if (!is.null(bench_df)) cms <- cms[intersect(bench_df$Modelo, names(cms))]
    cms
  }

  cms_bal   <- assemble_cms("balanced", bench_bal)
  cms_unbal <- assemble_cms("unbalanced", bench_unbal)

  has_models <- (!is.null(bench_bal) && nrow(bench_bal) > 0) || (length(models_bal) > 0)

  list(
    has_models = has_models,
    datos = datos,
    df_ml_report = df_ml_report,
    models_bal = models_bal,
    models_unbal = models_unbal,
    youden_thresholds = youden_thresholds,
    heuristic_model = heuristic_model,
    bench_bal = bench_bal,
    bench_unbal = bench_unbal,
    cms_bal = cms_bal,
    cms_unbal = cms_unbal,
    shap_obj_bal = shap_obj_bal,
    spearman_res = spearman_res
  )
}

#' Dibuja la cuadrícula resumen de 6 matrices de confusión (Top 2 Youden + Heurístico + Top 3 ML)
#'
#' @param bench_df Tabla de benchmark del escenario.
#' @param all_cms Lista nombrada de tablas 2x2 de confusión.
#' @param scen_txt Etiqueta descriptiva del escenario.
#' @export
render_top6_summary_cms_grid <- function(bench_df, all_cms, scen_txt = "Escenario Balanceado") {
  if (is.null(bench_df) || length(all_cms) == 0) {
    cat("> *No hay matrices de confusión disponibles para este escenario.*\n")
    return(invisible(NULL))
  }

  selected_6 <- intersect(get_top6_summary_models(bench_df), names(all_cms))
  plot_list <- lapply(selected_6, function(m_name) {
    tab_m <- all_cms[[m_name]]
    pal   <- get_model_color_palette(m_name)
    acc_m <- sum(diag(tab_m)) / sum(tab_m)
    row_i <- bench_df %>% dplyr::filter(Modelo == m_name)
    role_prefix <- if (grepl("^Youden_", m_name)) "Top Youden: " else if (m_name == "Weighted_Heuristic") "Heurístico: " else "Top ML: "
    sub_txt <- if (nrow(row_i) > 0) {
      sprintf("ROC-AUC: %.3f | F1: %.3f | Acc: %.1f%% | FDR: %.1f%%", row_i$ROC_AUC[1], row_i$F1_Score[1], acc_m * 100, row_i$FDR[1] * 100)
    } else {
      sprintf("Exactitud: %.1f%%", acc_m * 100)
    }
    plot_confusion_matrix(tab_m, title = paste0(role_prefix, m_name), subtitle = sub_txt, high_color = pal$high, low_color = pal$low)
  })

  title_banner <- cowplot::ggdraw() +
    cowplot::draw_label(
      sprintf("Resumen de Matrices Clave (%s): Top 2 Youden + Heurístico + Top 3 ML", scen_txt),
      fontface = "bold", size = 13.5, color = "#2c3e50"
    )
  grid_body <- cowplot::plot_grid(plotlist = plot_list, ncol = 2, align = "hv")
  print(cowplot::plot_grid(title_banner, grid_body, ncol = 1, rel_heights = c(0.06, 1)))
}

#' Renderiza la matriz de confusión e importancia de variables de un modelo individual
#'
#' @param mod_name Nombre del modelo.
#' @param models_list Lista de modelos serializados del escenario.
#' @param bench_df Tabla de benchmark del escenario.
#' @param all_cms Lista de matrices de confusión del escenario.
#' @param df_ml_report Data frame de entrenamiento preparado.
#' @param heuristic_model Objeto del modelo heurístico.
#' @param youden_thresholds Vector de umbrales de Youden.
#' @param scen_label Nombre legible del escenario.
#' @export
render_single_model_diagnostics <- function(mod_name,
                                            models_list,
                                            bench_df,
                                            all_cms,
                                            df_ml_report = NULL,
                                            heuristic_model = NULL,
                                            youden_thresholds = NULL,
                                            scen_label = "Escenario Balanceado") {
  pal <- get_model_color_palette(mod_name)
  row_m <- if (!is.null(bench_df)) bench_df %>% dplyr::filter(Modelo == mod_name) else NULL

  if (!is.null(row_m) && nrow(row_m) > 0) {
    fdr_txt <- if (is.na(row_m$FDR[1])) "NA" else sprintf("%.2f%%", row_m$FDR[1] * 100)
    prec_txt <- if (is.na(row_m$Precision[1])) "NA" else sprintf("%.2f%%", row_m$Precision[1] * 100)
    cat(sprintf(
      "<div class='model-badge' style='background-color:%s;'>%s (%s): F1: %.3f | ROC-AUC: %.3f | Exactitud Bal.: %.2f%% | Sensibilidad: %.2f%% | Especificidad: %.2f%% | Precisión: %s | FDR: %s</div>\n\n",
      pal$high, mod_name, scen_label,
      row_m$F1_Score[1], row_m$ROC_AUC[1],
      row_m$Exactitud_Balanceada[1] * 100,
      row_m$Sensibilidad[1] * 100,
      row_m$Especificidad[1] * 100,
      prec_txt, fdr_txt
    ))
  }

  tab_cm <- all_cms[[mod_name]]
  p_cm <- if (!is.null(tab_cm)) {
    acc_val <- sum(diag(tab_cm)) / sum(tab_cm)
    plot_confusion_matrix(
      tab_cm, title = paste("Matriz de Confusión:", mod_name),
      subtitle = sprintf("%s (Exactitud: %.2f%%)", scen_label, acc_val * 100),
      high_color = pal$high, low_color = pal$low
    )
  } else {
    NULL
  }

  p_imp <- NULL
  if (!is.null(models_list) && mod_name %in% names(models_list)) {
    mod_obj <- models_list[[mod_name]]
    if (inherits(mod_obj, "caretStack")) {
      coefs_s <- tryCatch(stats::coef(mod_obj$ens_model$finalModel), error = function(e) NULL)
      if (!is.null(coefs_s) && length(coefs_s) > 1) {
        p_imp <- plot_model_varimp(model_obj = abs(coefs_s[-1]), model_name = paste0("Pesos Meta-Modelo (", scen_label, ")"), high_color = pal$high)
      }
    } else {
      p_imp <- plot_model_varimp(model_obj = mod_obj, df_data = df_ml_report, model_name = paste0(mod_name, " (", scen_label, ")"), high_color = pal$high)
    }
  } else if (mod_name == "Weighted_Heuristic" && !is.null(heuristic_model)) {
    p_imp <- plot_model_varimp(model_obj = heuristic_model$weights, model_name = "Pesos del Modelo Heurístico", high_color = pal$high)
  } else if (grepl("^Youden_", mod_name) && !is.null(youden_thresholds)) {
    var_y <- sub("^Youden_", "", mod_name)
    vec_y <- stats::setNames(ifelse(names(youden_thresholds) == var_y, 100, 0), names(youden_thresholds))
    p_imp <- plot_model_varimp(
      model_obj = vec_y,
      model_name = paste0("Variable Activa: ", var_y, " (Corte >= ", round(youden_thresholds[[var_y]], 2), ")"),
      high_color = pal$high
    )
  }

  if (!is.null(p_cm) && !is.null(p_imp)) {
    print(cowplot::plot_grid(p_cm, p_imp, ncol = 2, align = "h"))
  } else if (!is.null(p_cm)) {
    print(p_cm)
  } else if (!is.null(p_imp)) {
    print(p_imp)
  }
}

#' Imprime los hiperparámetros óptimos y diagramas estructurales de un modelo entrenado
#'
#' @param mod_name Nombre del modelo.
#' @param models_list Lista de modelos serializados del escenario.
#' @param scen_label Nombre legible del escenario.
#' @export
render_single_model_extra <- function(mod_name, models_list, scen_label = "Escenario Balanceado") {
  if (is.null(models_list) || !mod_name %in% names(models_list)) return(invisible(NULL))
  mod_obj <- models_list[[mod_name]]

  if (mod_name %in% c("Red_Neuronal", "Red_Neuronal_Bag", "MLP", "MLP_Regularizado", "Stacked Autoencoder - DNN", "Autoencoder") && !is.null(mod_obj$finalModel)) {
    plot_neural_net_topology(mod_obj, title = sprintf("Arquitectura de Red: %s (%s)", mod_name, scen_label))
  } else if (mod_name == "CART" && requireNamespace("rpart.plot", quietly = TRUE) && !is.null(mod_obj$finalModel)) {
    tryCatch({
      rpart.plot::rpart.plot(
        mod_obj$finalModel, type = 3, extra = 104, under = TRUE, faclen = 0, cex = 0.85,
        main = sprintf("Árbol de Decisión CART (%s)", scen_label),
        box.palette = c("#27ae60", "#e74c3c")
      )
    }, error = function(e) {})
  }

  if (mod_name == "Red_Neuronal_Bag" && !is.null(mod_obj$finalModel)) {
    n_nets <- if (!is.null(mod_obj$finalModel$model)) length(mod_obj$finalModel$model) else NA
    topo_str <- if (!is.null(mod_obj$finalModel$model[[1]]$n)) paste(mod_obj$finalModel$model[[1]]$n, collapse = " - ") else "4 - size - 1"
    cat(sprintf("Parámetros de Red_Neuronal_Bag (avNNet - %s): sub-redes = %s | topología = %s | size = %s | decay = %.4g\n",
                scen_label, as.character(n_nets), topo_str, as.character(mod_obj$bestTune$size), as.numeric(mod_obj$bestTune$decay)))
  } else if (mod_name == "Red_Neuronal" && !is.null(mod_obj$finalModel)) {
    topo_str <- if (!is.null(mod_obj$finalModel$n)) paste(mod_obj$finalModel$n, collapse = " - ") else "4 - size - 1"
    cat(sprintf("Parámetros de Red_Neuronal (nnet - %s): topología = %s | size = %s | decay = %.4g\n",
                scen_label, topo_str, as.character(mod_obj$bestTune$size), as.numeric(mod_obj$bestTune$decay)))
  } else if (inherits(mod_obj, "caretStack")) {
    cat(sprintf("\nResumen del Meta-Modelo GLM Stacking (%s):\n", scen_label))
    print(mod_obj)
  } else if (mod_name == "ElasticNet" && !is.null(mod_obj$bestTune)) {
    cat(sprintf("\nHiperparámetros óptimos Elastic Net (%s): alpha = %.4g | lambda = %.4g\n", scen_label, mod_obj$bestTune$alpha, mod_obj$bestTune$lambda))
    print(stats::coef(mod_obj$finalModel, s = mod_obj$bestTune$lambda))
  } else if (mod_name == "LDA" && !is.null(mod_obj$finalModel$scaling)) {
    cat(sprintf("\nCoeficientes discriminantes lineales LD1 (%s):\n", scen_label))
    print(mod_obj$finalModel$scaling)
  } else if (!is.null(mod_obj$bestTune)) {
    cat(sprintf("\nHiperparámetros seleccionados en validación cruzada para %s (%s):\n", mod_name, scen_label))
    print(mod_obj$bestTune)
  }
}

#' Genera dinámicamente las pestañas Quarto de diagnóstico para todos los modelos de un escenario
#'
#' @param scenario Clave del escenario (`"balanced"` o `"unbalanced"`).
#' @param tr_art Lista devuelta por `load_readonly_training_artifacts()`.
#' @export
build_all_model_tabs_markdown <- function(scenario = c("balanced", "unbalanced"), tr_art) {
  scenario <- match.arg(scenario)
  bench_df    <- if (scenario == "balanced") tr_art$bench_bal else tr_art$bench_unbal
  all_cms     <- if (scenario == "balanced") tr_art$cms_bal else tr_art$cms_unbal
  models_list <- if (scenario == "balanced") tr_art$models_bal else tr_art$models_unbal
  scen_txt    <- if (scenario == "balanced") "Escenario Balanceado" else "Escenario Desbalanceado"

  ordered_models <- if (!is.null(bench_df) && nrow(bench_df) > 0) intersect(bench_df$Modelo, names(all_cms)) else names(all_cms)
  ordered_models <- c(ordered_models, setdiff(names(all_cms), ordered_models))

  child_Env <- new.env(parent = globalenv())
  child_Env$bench_df          <- bench_df
  child_Env$all_cms           <- all_cms
  child_Env$models_list       <- models_list
  child_Env$scen_txt          <- scen_txt
  child_Env$df_ml_report      <- tr_art$df_ml_report
  child_Env$heuristic_model   <- tr_art$heuristic_model
  child_Env$youden_thresholds <- tr_art$youden_thresholds

  bt3 <- strrep("`", 3)
  child_chunks <- sprintf(
    "#### Resumen: Top 6 (2 Youden + Heurístico + 3 ML)\n\n%s{r top6_%s, fig.width=10, fig.height=10.2, echo=FALSE, warning=FALSE, message=FALSE}\nrender_top6_summary_cms_grid(bench_df, all_cms, scen_txt)\n%s\n",
    bt3, scenario, bt3
  )

  for (idx in seq_along(ordered_models)) {
    m_name <- ordered_models[idx]
    safe_id <- gsub("[^a-zA-Z0-9_]", "_", paste0(scenario, "_", m_name, "_", idx))
    chunk_m <- sprintf(
      "#### %s\n\n%s{r diag_%s, fig.width=9.5, fig.height=4.6, results='asis', echo=FALSE, warning=FALSE, message=FALSE}\nrender_single_model_diagnostics('%s', models_list, bench_df, all_cms, df_ml_report, heuristic_model, youden_thresholds, scen_txt)\n%s\n\n%s{r extra_%s, fig.width=8.5, fig.height=4.8, echo=FALSE, warning=FALSE, message=FALSE}\nrender_single_model_extra('%s', models_list, scen_txt)\n%s\n",
      m_name, bt3, safe_id, m_name, bt3, bt3, safe_id, m_name, bt3
    )
    child_chunks <- c(child_chunks, chunk_m)
  }

  res_md <- knitr::knit_child(text = paste(child_chunks, collapse = "\n\n"), envir = child_Env, quiet = TRUE)
  cat(res_md, sep = "\n")
}

#' Genera el explorador interactivo HTML/Plotly de contribuciones SHAP por secuencia (Waterfall)
#'
#' @param shap_obj Objeto `shapviz` serializado.
#' @param df_ml_report Data frame de entrenamiento con metadatos de locus y clase real.
#' @param proj_root Directorio raíz del proyecto para localizar `assets/waterfall.js`.
#' @export
render_shap_waterfall_html <- function(shap_obj, df_ml_report, proj_root = ".") {
  if (is.null(shap_obj) || is.null(df_ml_report)) {
    cat("> *No hay valores SHAP disponibles para construir el explorador individual.*\n")
    return(invisible(NULL))
  }

  S_mat <- as.matrix(shap_obj$S)
  X_mat <- as.data.frame(shap_obj$X)
  base_val <- as.numeric(shap_obj$baseline)
  n_rows <- min(nrow(S_mat), nrow(df_ml_report))

  wf_records <- lapply(seq_len(n_rows), function(i) {
    row_d <- df_ml_report[i, , drop = FALSE]
    gen_ref <- if ("Locus_Ref_Limpio" %in% colnames(row_d) && !is.na(row_d$Locus_Ref_Limpio) && nzchar(as.character(row_d$Locus_Ref_Limpio))) {
      as.character(row_d$Locus_Ref_Limpio)
    } else if ("locus" %in% colnames(row_d) && !is.na(row_d$locus)) {
      as.character(row_d$locus)
    } else {
      "Sin_Gen_Anotado"
    }
    og_id <- if ("locus" %in% colnames(row_d) && !is.na(row_d$locus)) {
      as.character(row_d$locus)
    } else if ("OG_Limpio" %in% colnames(row_d) && !is.na(row_d$OG_Limpio)) {
      as.character(row_d$OG_Limpio)
    } else {
      "N/A"
    }
    seq_id <- if ("ID_GFF_Compatible" %in% colnames(row_d) && !is.na(row_d$ID_GFF_Compatible)) {
      as.character(row_d$ID_GFF_Compatible)
    } else {
      paste0("Secuencia_", i)
    }
    shap_vec <- as.numeric(S_mat[i, ])
    list(
      row = i, seq_id = seq_id, gene = gen_ref, og = og_id, clase = as.character(row_d$Clase),
      base = round(base_val, 4), fx = round(min(1, max(0, base_val + sum(shap_vec))), 4),
      features = colnames(S_mat), values = round(as.numeric(X_mat[i, ]), 2), shaps = round(shap_vec, 4)
    )
  })

  wf_json <- jsonlite::toJSON(wf_records, auto_unbox = TRUE)
  js_candidates <- c(
    file.path(proj_root, "assets", "waterfall.js"),
    file.path("assets", "waterfall.js")
  )
  js_file <- js_candidates[file.exists(js_candidates)][1]
  js_body <- if (!is.na(js_file) && file.exists(js_file)) {
    paste(readLines(js_file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  } else {
    ""
  }

  fence_open  <- paste0("\n", strrep("`", 3), "{=html}\n")
  fence_close <- paste0("\n", strrep("`", 3), "\n\n")

  cat(paste0(
    fence_open,
    '<div class="wf-controls"><div style="display:flex; flex-wrap:wrap; gap:12px; align-items:center;">',
    '<div><label for="wf-class-filter" style="font-weight:bold; margin-right:6px;">Filtrar por Clase Real:</label>',
    '<select id="wf-class-filter" style="padding:6px 10px; border-radius:5px; border:1px solid #bdc3c7;">',
    '<option value="ALL">Todas las filas (', n_rows, ')</option>',
    '<option value="TP">Solo TP (Verdaderos Positivos)</option>',
    '<option value="FP">Solo FP (Falsos Positivos)</option></select></div>',
    '<div style="flex-grow:1;"><label for="wf-row-select" style="font-weight:bold; margin-right:6px;">Seleccionar Fila / Gen:</label>',
    '<select id="wf-row-select" style="padding:6px 10px; border-radius:5px; border:1px solid #bdc3c7; width:100%; max-width:540px;"></select></div>',
    '<div><button id="wf-prev-btn" type="button" style="padding:6px 12px; border-radius:5px; border:1px solid #2c3e50; background:#fff; cursor:pointer; font-weight:bold;">Anterior</button> ',
    '<button id="wf-next-btn" type="button" style="padding:6px 12px; border-radius:5px; border:1px solid #2c3e50; background:#fff; cursor:pointer; font-weight:bold;">Siguiente</button></div>',
    '</div></div>',
    '<div id="wf-meta-info" class="wf-meta-card"></div>',
    '<div id="wf-plotly-container" style="width:100%; height:460px;"></div>',
    '<script>\nwindow.TARGETCNV_WF_DATA = ', wf_json, ';\n',
    js_body, '\n</script>',
    fence_close
  ))
}

#' Dibuja cuadrículas comparativas de matrices de confusión para transferencia filogenética
#'
#' @param res_sp Elemento de especie dentro de `transfer_results` (`info`, `bal`, `unbal`).
#' @param mode Modo de selección: `"top6"`, `"thesis9"` o `"top3_f1"`.
#' @export
render_species_cms_grid <- function(res_sp, mode = c("top6", "thesis9", "top3_f1")) {
  mode <- match.arg(mode)
  bench_df <- res_sp$bal$benchmark
  all_cms  <- res_sp$bal$cms
  sp_name  <- res_sp$info$short_name
  sc_label <- if (res_sp$info$genome_scale == "chromosome") "Escala Cromosómica" else "Escala de Contigs"

  if (is.null(bench_df) || length(all_cms) == 0) {
    cat("> *No hay matrices de confusión disponibles para esta especie.*\n")
    return(invisible(NULL))
  }

  if (mode == "top6") {
    selected_mods <- intersect(get_top6_summary_models(bench_df), names(all_cms))
    n_cols <- 2
    banner_txt <- sprintf("Resumen de Matrices Clave en %s (%s): Top 2 Youden + Heurístico + Top 3 ML", sp_name, sc_label)
  } else if (mode == "thesis9") {
    selected_mods <- intersect(get_thesis_9_models(bench_df), names(all_cms))
    n_cols <- 3
    banner_txt <- sprintf("Transferencia en %s (%s): 9 Modelos Estandarizados (Cortes Empíricos, ML y Stacking)", sp_name, sc_label)
  } else {
    ml_only <- bench_df %>%
      dplyr::filter(!grepl("^Youden_", Modelo) & !Modelo %in% c("Weighted_Heuristic", "Ensamble_Stacking", "Consenso_7_Modelos")) %>%
      dplyr::arrange(dplyr::desc(F1_Score), dplyr::desc(ROC_AUC)) %>%
      dplyr::slice_head(n = 3) %>%
      dplyr::pull(Modelo)
    selected_mods <- intersect(ml_only, names(all_cms))
    n_cols <- 3
    banner_txt <- sprintf("Top 3 Clasificadores Supervisados en %s (%s - Ordenados por F1-Score)", sp_name, sc_label)
  }

  plot_list <- lapply(selected_mods, function(m_name) {
    tab_m <- all_cms[[m_name]]
    pal   <- get_model_color_palette(m_name)
    acc_m <- sum(diag(tab_m)) / sum(tab_m)
    row_i <- bench_df %>% dplyr::filter(Modelo == m_name)
    fdr_s <- if (nrow(row_i) > 0 && !is.na(row_i$FDR[1])) sprintf("%.1f%%", row_i$FDR[1] * 100) else "NA"
    sub_txt <- if (nrow(row_i) > 0) {
      sprintf("F1: %.3f | Acc: %.1f%% | Sens: %.1f%% | FDR: %s", row_i$F1_Score[1], acc_m * 100, row_i$Sensibilidad[1] * 100, fdr_s)
    } else {
      sprintf("Exactitud: %.1f%%", acc_m * 100)
    }
    plot_confusion_matrix(tab_m, title = m_name, subtitle = sub_txt, high_color = pal$high, low_color = pal$low)
  })

  title_banner <- cowplot::ggdraw() + cowplot::draw_label(banner_txt, fontface = "bold", size = 13, color = "#2c3e50")
  grid_body <- cowplot::plot_grid(plotlist = plot_list, ncol = n_cols, align = "hv")
  print(cowplot::plot_grid(title_banner, grid_body, ncol = 1, rel_heights = c(0.06, 1)))
}

#' Renderiza el diagnóstico individual (matriz + raincloud) de un modelo transferido
#'
#' @param res_sp Elemento de especie dentro de `transfer_results`.
#' @param mod_name Nombre del modelo.
#' @export
render_single_transfer_model <- function(res_sp, mod_name) {
  bench_df <- res_sp$bal$benchmark
  all_cms  <- res_sp$bal$cms
  scores_l <- res_sp$bal$scores
  df_clean <- res_sp$bal$df_clean
  sp_name  <- res_sp$info$short_name
  pal      <- get_model_color_palette(mod_name)

  row_m <- bench_df %>% dplyr::filter(Modelo == mod_name)
  if (nrow(row_m) > 0) {
    fdr_txt <- if (is.na(row_m$FDR[1])) "NA" else sprintf("%.2f%%", row_m$FDR[1] * 100)
    prec_txt <- if (is.na(row_m$Precision[1])) "NA" else sprintf("%.2f%%", row_m$Precision[1] * 100)
    cat(sprintf(
      "<div class='model-badge' style='background-color:%s;'>%s en %s: F1: %.3f | ROC-AUC: %.3f | Exactitud: %.2f%% | Exactitud Bal.: %.2f%% | Sensibilidad: %.2f%% | Especificidad: %.2f%% | Precisión: %s | FDR: %s</div>\n\n",
      pal$high, mod_name, sp_name,
      row_m$F1_Score[1], row_m$ROC_AUC[1],
      row_m$Exactitud[1] * 100, row_m$Exactitud_Balanceada[1] * 100,
      row_m$Sensibilidad[1] * 100, row_m$Especificidad[1] * 100,
      prec_txt, fdr_txt
    ))
  }

  tab_cm <- all_cms[[mod_name]]
  p_cm <- if (!is.null(tab_cm)) {
    acc_val <- sum(diag(tab_cm)) / sum(tab_cm)
    plot_confusion_matrix(
      tab_cm, title = paste0("Matriz de Confusión: ", mod_name),
      subtitle = sprintf("%s (Exactitud: %.2f%%)", sp_name, acc_val * 100),
      high_color = pal$high, low_color = pal$low
    )
  } else {
    NULL
  }

  p_dist <- NULL
  sc_vec <- scores_l[[mod_name]]
  if (!is.null(sc_vec) && length(sc_vec) == nrow(df_clean)) {
    df_sc <- data.frame(Clase_Real = df_clean$Clase_Real, Puntaje = as.numeric(sc_vec)) %>%
      dplyr::filter(!is.na(Puntaje))
    if (nrow(df_sc) > 2) p_dist <- plot_transfer_raincloud(df_sc, mod_name = mod_name, sp_name = sp_name)
  }

  if (!is.null(p_cm) && !is.null(p_dist)) {
    print(cowplot::plot_grid(p_cm, p_dist, ncol = 2, align = "h", rel_widths = c(0.95, 1.05)))
  } else if (!is.null(p_cm)) {
    print(p_cm)
  }
}

#' Genera las subpestañas de matrices y rainclouds por modelo para una especie receptora
#'
#' @param sp_key Clave de la especie en `transfer_results`.
#' @param transfer_results Lista completa leída desde `transferability_results.rds`.
#' @export
build_species_model_tabs_markdown <- function(sp_key, transfer_results) {
  res_sp   <- transfer_results[[sp_key]]
  bench_df <- res_sp$bal$benchmark
  all_cms  <- res_sp$bal$cms

  ordered_models <- unique(c(intersect(bench_df$Modelo, names(all_cms)), names(all_cms)))
  child_Env <- new.env(parent = globalenv())
  child_Env$res_sp <- res_sp

  bt3 <- strrep("`", 3)
  child_chunks <- c(
    sprintf("##### Resumen: Top 6 (2 Youden + Heurístico + 3 ML)\n\n%s{r top6_%s, fig.width=10, fig.height=10.2, echo=FALSE, warning=FALSE, message=FALSE}\nrender_species_cms_grid(res_sp, 'top6')\n%s\n", bt3, sp_key, bt3),
    sprintf("##### Resumen Tesis: 9 Modelos Estándar (3x3)\n\n%s{r thesis9_%s, fig.width=11, fig.height=10.5, echo=FALSE, warning=FALSE, message=FALSE}\nrender_species_cms_grid(res_sp, 'thesis9')\n%s\n", bt3, sp_key, bt3),
    sprintf("##### Top 3 ML Específicos (F1-Score)\n\n%s{r top3f1_%s, fig.width=10.5, fig.height=4.3, echo=FALSE, warning=FALSE, message=FALSE}\nrender_species_cms_grid(res_sp, 'top3_f1')\n%s\n", bt3, sp_key, bt3)
  )

  for (idx in seq_along(ordered_models)) {
    m_name  <- ordered_models[idx]
    safe_id <- gsub("[^a-zA-Z0-9_]", "_", paste0(sp_key, "_", m_name, "_", idx))
    child_chunks <- c(child_chunks, sprintf(
      "##### %s\n\n%s{r mod_%s, fig.width=10.2, fig.height=4.7, results='asis', echo=FALSE, warning=FALSE, message=FALSE}\nrender_single_transfer_model(res_sp, '%s')\n%s\n",
      m_name, bt3, safe_id, m_name, bt3
    ))
  }

  res_md <- knitr::knit_child(text = paste(child_chunks, collapse = "\n\n"), envir = child_Env, quiet = TRUE)
  cat(res_md, sep = "\n")
}

#' Construye el gráfico comparativo multiespecie interactivo en Plotly
#'
#' @param df_multiespecie_all Data frame consolidado con métricas del donante y especies receptoras.
#' @return Objeto `plotly` o `NULL`.
#' @export
build_multispecies_comparison_plotly <- function(df_multiespecie_all) {
  if (is.null(df_multiespecie_all) || nrow(df_multiespecie_all) == 0) return(NULL)
  key_mods <- c(
    "Youden_pct_recovered", "Youden_pct_identity", "Weighted_Heuristic",
    "ElasticNet", "NaiveBayes", "SVM_Polinomial", "XGBoost", "MLP",
    "Red_Neuronal", "LDA", "Ensamble_Stacking", "Consenso_7_Modelos"
  )
  df_comp_sub <- df_multiespecie_all %>%
    dplyr::filter(Modelo %in% key_mods) %>%
    dplyr::mutate(Etiqueta_Especie = paste0(Especie, " [", Escala_Genoma, "]"))

  sp_groups <- unique(df_comp_sub$Etiqueta_Especie)
  palette_sp <- c("#2c3e50", "#16a085", "#d35400", "#8e44ad", "#2980b9")
  metrics_comp <- c("F1_Score", "ROC_AUC", "Exactitud_Balanceada", "Sensibilidad", "Especificidad", "FDR")

  fig_comp <- plotly::plot_ly(height = 680)
  for (m_idx in seq_along(metrics_comp)) {
    m_col <- metrics_comp[m_idx]
    for (s_idx in seq_along(sp_groups)) {
      sp_lbl <- sp_groups[s_idx]
      sub_d <- df_comp_sub %>% dplyr::filter(Etiqueta_Especie == sp_lbl)
      vals_m <- round(as.numeric(sub_d[[m_col]]), 3)
      fig_comp <- fig_comp %>% plotly::add_trace(
        y = sub_d$Modelo, x = vals_m, name = sp_lbl, type = "bar", orientation = "h",
        visible = (m_idx == 1),
        marker = list(color = palette_sp[((s_idx - 1) %% length(palette_sp)) + 1]),
        text = sprintf("%.3f", vals_m), textposition = "auto",
        hovertemplate = paste0("<b>%{y}</b><br>", sp_lbl, "<br>", m_col, ": %{x:.3f}<extra></extra>")
      )
    }
  }

  btns_comp <- lapply(seq_along(metrics_comp), function(m_idx) {
    vis_v <- rep(FALSE, length(metrics_comp) * length(sp_groups))
    idx_start <- (m_idx - 1) * length(sp_groups) + 1
    vis_v[idx_start:(idx_start + length(sp_groups) - 1)] <- TRUE
    list(method = "restyle", args = list("visible", as.list(vis_v)), label = paste("Métrica:", metrics_comp[m_idx]))
  })

  fig_comp %>% plotly::layout(
    title = list(text = "<b>Comparativa Inter-Especie según Escala Estructural del Genoma</b>", x = 0.5),
    barmode = "group",
    xaxis = list(title = "Puntaje (0 a 1)", range = c(0, 1.08)),
    yaxis = list(title = ""),
    legend = list(orientation = "h", y = -0.12, x = 0.5, xanchor = "center"),
    margin = list(l = 190, r = 35, t = 85, b = 110),
    updatemenus = list(list(type = "dropdown", active = 0, x = 0.0, y = 1.10, xanchor = "left", buttons = btns_comp))
  )
}

#' Renderiza las pestañas Quarto multiespecie de solo lectura para `transferability_report.qmd`
#'
#' @param transfer_results Lista precomputada por `scripts/04_evaluate_transferability.R`.
#' @param df_multiespecie_all Tabla comparativa global multiespecie.
#' @export
render_transfer_multispecies_tabs <- function(transfer_results, df_multiespecie_all = NULL) {
  if (is.null(transfer_results) || length(transfer_results) == 0) {
    cat(paste0(
      '<div class="missing-alert">',
      '<b>No se encontraron resultados precomputados de transferencia filogenética.</b><br>',
      'Ejecute primero <code>./targetcnv.sh transfer</code> (o <code>Rscript scripts/04_evaluate_transferability.R</code>) ',
      'para evaluar los modelos entrenados y generar los archivos en <code>results/</code>.',
      '</div>\n'
    ))
    return(invisible(NULL))
  }

  child_Env <- new.env(parent = globalenv())
  child_Env$transfer_results    <- transfer_results
  child_Env$df_multiespecie_all <- df_multiespecie_all

  bt3 <- strrep("`", 3)
  sp_tabs_md <- c("::: {.panel-tabset}\n")

  for (sp_key in names(transfer_results)) {
    res_sp   <- transfer_results[[sp_key]]
    sp_info  <- res_sp$info
    is_chr   <- (sp_info$genome_scale == "chromosome")
    sc_short <- if (is_chr) "Escala Cromosómica" else "Escala de Contigs"
    box_cls  <- if (is_chr) "genome-arch-box" else "genome-arch-box contig"
    badge_cl <- if (is_chr) "chr" else "ctg"

    df_c   <- res_sp$bal$df_clean
    n_tot  <- nrow(df_c)
    n_tp   <- sum(df_c$Clase_Real == "TP", na.rm = TRUE)
    n_fp   <- sum(df_c$Clase_Real == "FP", na.rm = TRUE)
    pct_tp <- round((n_tp / max(1, n_tot)) * 100, 1)
    pct_fp <- round((n_fp / max(1, n_tot)) * 100, 1)

    bridge_txt <- if (!is_chr && !is.null(sp_info$n_ogs_captus) && !is.na(sp_info$n_ogs_captus)) {
      sprintf(" | Cobertura Tabla Puente OG_ID <-> HOG_ID: <b>%d de %d OGs</b>", sp_info$n_ogs_bridge, sp_info$n_ogs_captus)
    } else {
      ""
    }

    src_lines_html <- if (isTRUE(sp_info$has_citation)) {
      paste0(
        sprintf("<b>Fuente del Genoma de Referencia:</b> %s<br>\n", sp_info$genome_source),
        sprintf("<b>Fuente de Homología y Datos Génicos:</b> %s%s<br>\n", sp_info$db_info, bridge_txt),
        sprintf("<b>Archivo de Origen:</b> <code>%s</code><br>\n", sp_info$origin_file)
      )
    } else {
      sprintf("<b>Archivo de Origen:</b> <code>%s</code>%s<br>\n", sp_info$origin_file, bridge_txt)
    }

    header_block <- paste0(
      sprintf("# *%s* (%s)\n\n", sp_info$short_name, sc_short),
      sprintf('<div class="%s">\n', box_cls),
      sprintf('<span class="scale-badge %s">%s</span> <b>Especie Validada:</b> <i>%s</i> (Muestra: <code>%s</code>)<br>\n',
              badge_cl, sc_short, sp_info$sci_name, sp_info$id),
      src_lines_html,
      sprintf("<b>Estrategia de Validación Espacial:</b> %s\n", sp_info$algo_info),
      "</div>\n\n",
      '<div class="kpi-grid">',
      sprintf('<div class="kpi-card"><div class="kpi-title">Secuencias NUC Evaluadas</div><div class="kpi-value">%d</div><div class="kpi-sub">Muestra %s</div></div>', n_tot, sp_info$id),
      sprintf('<div class="kpi-card green"><div class="kpi-title">Verdaderos Positivos (TP)</div><div class="kpi-value">%d (%s%%)</div><div class="kpi-sub">Homólogos validados</div></div>', n_tp, pct_tp),
      sprintf('<div class="kpi-card red"><div class="kpi-title">Falsos Positivos (FP)</div><div class="kpi-value">%d (%s%%)</div><div class="kpi-sub">Ruido o redundancia</div></div>', n_fp, pct_fp),
      sprintf('<div class="kpi-card purple"><div class="kpi-title">Escala del Ensamblaje</div><div class="kpi-value">%s</div><div class="kpi-sub">Modo: %s</div></div>',
              sc_short, if (is_chr) "overlap (GenomicRanges)" else "target_attr + Puente OG/HOG"),
      "</div>\n\n"
    )

    chunk_gt <- paste0(
      sprintf("### Resumen de Secuencias Evaluadas (TP vs. FP) en *%s*\n\n", sp_info$short_name),
      sprintf("%s{r gt_summary_%s, echo=FALSE, warning=FALSE, message=FALSE}\n", bt3, sp_key),
      sprintf("df_sp_curr <- transfer_results[['%s']]$bal$df_clean\n", sp_key),
      "df_res_cl <- df_sp_curr %>%\n",
      "  dplyr::mutate(Clase = ifelse(Clase_Real == 'TP', 'TP (Verdadero Positivo)', 'FP (Falso Positivo)')) %>%\n",
      "  dplyr::group_by(Clase) %>%\n",
      "  dplyr::summarise(\n",
      "    Secuencias = dplyr::n(),\n",
      "    Porcentaje = round((dplyr::n() / nrow(df_sp_curr)) * 100, 2),\n",
      "    `Media % Recuperado` = round(mean(as.numeric(pct_recovered), na.rm = TRUE), 2),\n",
      "    `Media % Identidad`  = round(mean(as.numeric(pct_identity), na.rm = TRUE), 2),\n",
      "    `Media Score`        = round(mean(as.numeric(score), na.rm = TRUE), 2),\n",
      "    `Media WScore`       = round(mean(as.numeric(wscore), na.rm = TRUE), 2),\n",
      "    .groups = 'drop'\n",
      "  ) %>% dplyr::arrange(dplyr::desc(Clase))\n",
      sprintf("DT::datatable(df_res_cl, rownames = FALSE, caption = 'Resumen de secuencias en %s (Total = %d NUC)', options = list(dom = 't', ordering = FALSE)) %%>%% DT::formatCurrency('Porcentaje', currency = '%%', before = FALSE, digits = 2)\n", sp_info$short_name, n_tot),
      sprintf("%s\n\n", bt3)
    )

    chunk_bench <- paste0(
      sprintf("### Benchmark Multimétrica de Transferencia en *%s*\n\n", sp_info$short_name),
      sprintf("%s{r plot_top6_%s, echo=FALSE, warning=FALSE, message=FALSE}\n", bt3, sp_key),
      sprintf("b_full <- transfer_results[['%s']]$bal$benchmark\n", sp_key),
      "df_t6 <- b_full %>% dplyr::filter(Modelo %in% get_top6_summary_models(b_full)) %>% dplyr::arrange(dplyr::desc(ROC_AUC), dplyr::desc(F1_Score))\n",
      sprintf("build_interactive_benchmark_plotly(df_t6, context_title = '%s (%s)', view_label = '6 Modelos Clave')\n", sp_info$short_name, sc_short),
      sprintf("%s\n\n", bt3),
      "<details ontoggle=\"if(this.open){setTimeout(function(){window.dispatchEvent(new Event('resize'));}, 60);}\">\n",
      sprintf("<summary class=\"summary-btn\">Ver Benchmark de los 9 Modelos Estandarizados (%s)</summary>\n\n", sp_info$short_name),
      sprintf("%s{r plot_thesis9_%s, echo=FALSE, warning=FALSE, message=FALSE}\n", bt3, sp_key),
      sprintf("b_full <- transfer_results[['%s']]$bal$benchmark\n", sp_key),
      "df_t9 <- b_full %>% dplyr::filter(Modelo %in% get_thesis_9_models(b_full)) %>% dplyr::arrange(dplyr::desc(ROC_AUC), dplyr::desc(F1_Score))\n",
      sprintf("build_interactive_benchmark_plotly(df_t9, context_title = '%s (%s)', view_label = '9 Modelos Estándar')\n", sp_info$short_name, sc_short),
      sprintf("%s\n\n</details>\n\n", bt3),
      "<details ontoggle=\"if(this.open){setTimeout(function(){window.dispatchEvent(new Event('resize'));}, 60);}\">\n",
      sprintf("<summary class=\"summary-btn\">Ver Benchmark de Todos los Modelos Evaluados (%s)</summary>\n\n", sp_info$short_name),
      sprintf("%s{r plot_all_%s, echo=FALSE, warning=FALSE, message=FALSE}\n", bt3, sp_key),
      sprintf("b_full <- transfer_results[['%s']]$bal$benchmark\n", sp_key),
      sprintf("build_interactive_benchmark_plotly(b_full, context_title = '%s (%s)', view_label = paste0('Todos los ', nrow(b_full), ' Modelos'))\n", sp_info$short_name, sc_short),
      sprintf("%s\n\n</details>\n\n", bt3),
      "<details ontoggle=\"if(this.open){setTimeout(function(){window.dispatchEvent(new Event('resize')); if(window.jQuery && jQuery.fn.dataTable){jQuery.fn.dataTable.tables({visible:true, api:true}).columns.adjust();}}, 60);}\">\n",
      "<summary class=\"summary-btn\">Ver Tabla Completa de Métricas</summary>\n\n",
      "<div class=\"asterisk-note\">Los modelos marcados con <code>*</code> logran menor FDR o mayor índice global al transferirse desde su versión entrenada en el escenario desbalanceado de <i>S. torvum</i>.</div>\n\n",
      sprintf("%s{r table_full_%s, echo=FALSE, warning=FALSE, message=FALSE}\n", bt3, sp_key),
      sprintf("format_benchmark_dt(transfer_results[['%s']]$bal$benchmark, other_df = transfer_results[['%s']]$unbal$benchmark, current_scen = 'Balanceado', other_scen = 'Desbalanceado', caption_text = 'Transferencia en %s (%s)')\n", sp_key, sp_key, sp_info$short_name, sc_short),
      sprintf("%s\n\n</details>\n\n", bt3)
    )

    chunk_cms <- paste0(
      sprintf("### Matrices de Confusión y Distribución de Puntajes (*%s*)\n\n", sp_info$short_name),
      "::: {.panel-tabset}\n\n",
      sprintf("%s{r tabs_cms_%s, results='asis', echo=FALSE, warning=FALSE, message=FALSE}\n", bt3, sp_key),
      sprintf("build_species_model_tabs_markdown('%s', transfer_results)\n", sp_key),
      sprintf("%s\n\n:::\n\n", bt3)
    )

    sp_tabs_md <- c(sp_tabs_md, header_block, chunk_gt, chunk_bench, chunk_cms)
  }

  comp_tab_md <- paste0(
    "# Comparativa Multiespecie\n\n",
    "### Contraste de Desempeño según Escala del Genoma (*S. torvum* vs. Receptoras)\n\n",
    sprintf("%s{r multispecies_comp_plot, echo=FALSE, warning=FALSE, message=FALSE}\n", bt3),
    "build_multispecies_comparison_plotly(df_multiespecie_all)\n",
    sprintf("%s\n\n", bt3),
    "### Tabla Consolidada Multiespecie (`transferability_all_species.tsv`)\n\n",
    sprintf("%s{r multispecies_comp_table, echo=FALSE, warning=FALSE, message=FALSE}\n", bt3),
    "if (!is.null(df_multiespecie_all) && nrow(df_multiespecie_all) > 0) {\n",
    "  df_disp <- df_multiespecie_all\n",
    "  pct_c <- intersect(c('Exactitud', 'Exactitud_Balanceada', 'Sensibilidad', 'Especificidad', 'Precision', 'VPN', 'FDR'), colnames(df_disp))\n",
    "  dec_c <- intersect(c('F1_Score', 'MCC', 'Kappa', 'ROC_AUC'), colnames(df_disp))\n",
    "  for (cc in pct_c) df_disp[[cc]] <- round(as.numeric(df_disp[[cc]]) * 100, 2)\n",
    "  for (cc in dec_c) df_disp[[cc]] <- round(as.numeric(df_disp[[cc]]), 3)\n",
    "  make_export_dt(df_disp, caption_txt = 'Tabla comparativa multiespecie (Donante + Especies Receptoras)', page_len = 30) %>%\n",
    "    DT::formatCurrency(columns = pct_c, currency = '%', before = FALSE, digits = 2)\n",
    "}\n",
    sprintf("%s\n\n", bt3)
  )

  sp_tabs_md <- c(sp_tabs_md, comp_tab_md, ":::\n")
  res_all <- knitr::knit_child(text = paste(sp_tabs_md, collapse = "\n"), envir = child_Env, quiet = TRUE)
  cat(res_all, sep = "\n")
}

