# ==============================================================================
# TargetCNV-ML: Rutinas de Visualizacion Cientifica con ggplot2 (R/plotting.R)
# ==============================================================================

#' Grafica una matriz de confusion 2x2 estilizada con TP en la esquina superior izquierda
#'
#' @param cm_table Tabla o matriz 2x2 de confusion.
#' @param title Titulo principal del grafico.
#' @param subtitle Subtitulo con metricas resumidas (ej. Accuracy y FDR).
#' @param high_color Color hexadecimal para la celda de mayor frecuencia.
#' @return Objeto `ggplot`.
#' @export
plot_confusion_matrix <- function(cm_table,
                                  title = "Matriz de Confusion",
                                  subtitle = NULL,
                                  high_color = "#d35400") {
  stopifnot(is.table(cm_table) || is.matrix(cm_table))

  cm_df <- as.data.frame(as.table(cm_table))
  colnames(cm_df) <- c("Prediccion", "Referencia", "Frecuencia")

  cm_df <- cm_df %>%
    dplyr::mutate(
      Referencia = dplyr::case_when(
        grepl("TRUE|TP", Referencia, ignore.case = TRUE) ~ "TP",
        TRUE ~ "FP"
      ),
      Prediccion = dplyr::case_when(
        grepl("TRUE|TP", Prediccion, ignore.case = TRUE) ~ "TP",
        TRUE ~ "FP"
      ),
      Referencia = factor(Referencia, levels = c("TP", "FP")),
      Prediccion = factor(Prediccion, levels = c("FP", "TP"))
    )

  max_freq <- max(cm_df$Frecuencia, na.rm = TRUE)
  if (is.na(max_freq) || max_freq == 0) max_freq <- 1

  ggplot2::ggplot(cm_df, ggplot2::aes(x = Referencia, y = Prediccion)) +
    ggplot2::geom_tile(ggplot2::aes(fill = Frecuencia), color = "white", linewidth = 1) +
    ggplot2::geom_text(
      ggplot2::aes(label = round(Frecuencia, 0)),
      size = 7,
      fontface = "bold",
      color = ifelse(cm_df$Frecuencia > max_freq / 2, "white", "black")
    ) +
    ggplot2::scale_fill_gradient(low = "#f1f2f6", high = high_color, name = "Frecuencia") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Clase Real (Ground Truth)",
      y = "Prediccion"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 11),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "#7f8c8d", size = 10),
      axis.title = ggplot2::element_text(face = "bold", size = 10),
      legend.position = "none",
      panel.grid = ggplot2::element_blank()
    )
}

#' Grafica las distribuciones de densidad de las metricas continuas (TP vs FP)
#'
#' @param df Data frame con las metricas y la columna de clase.
#' @param response_col Nombre de la columna de respuesta ("Validacion_Espacial" o "Clase").
#' @param metric_cols Vector con los nombres de las metricas continuas.
#' @param thresholds Vector nombrado opcional con los umbrales de Youden para trazar lineas verticales.
#' @param colors Vector nombrado de colores por clase.
#' @param title Titulo del grafico.
#' @param subtitle Subtitulo del grafico.
#' @param legend_position Posicion de la leyenda ("bottom", "right", "none").
#' @return Objeto `ggplot`.
#' @export
plot_metric_densities <- function(df,
                                  response_col = "Clase",
                                  metric_cols = c("pct_recovered", "pct_identity", "score", "wscore"),
                                  thresholds = NULL,
                                  colors = c("TRUE_POSITIVE" = "#27ae60", "FALSE_POSITIVE" = "#e74c3c",
                                             "TP" = "#27ae60", "FP" = "#e74c3c"),
                                  title = "Zonas de Solapamiento y Separacion (TP vs FP)",
                                  subtitle = "Distribucion de densidad para las 4 metricas continuas",
                                  legend_position = "bottom") {
  stopifnot(is.data.frame(df))
  stopifnot(response_col %in% colnames(df))
  stopifnot(all(metric_cols %in% colnames(df)))

  df_long <- df %>%
    dplyr::select(dplyr::all_of(c(response_col, metric_cols))) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(metric_cols),
      names_to = "Metrica",
      values_to = "Valor"
    )

  p <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(x = Valor, fill = .data[[response_col]], color = .data[[response_col]])
  ) +
    ggplot2::geom_density(alpha = 0.45, linewidth = 0.8) +
    ggplot2::facet_wrap(~ Metrica, scales = "free") +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Valor de la Metrica",
      y = "Densidad",
      fill = "Clase",
      color = "Clase"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "#7f8c8d"),
      legend.position = legend_position,
      strip.text = ggplot2::element_text(face = "bold")
    )

  if (!is.null(thresholds) && length(thresholds) > 0) {
    df_thr <- data.frame(
      Metrica = names(thresholds),
      Umbral = as.numeric(thresholds),
      stringsAsFactors = FALSE
    ) %>%
      dplyr::filter(Metrica %in% metric_cols)
    if (nrow(df_thr) > 0) {
      p <- p + ggplot2::geom_vline(
        data = df_thr,
        ggplot2::aes(xintercept = Umbral),
        linetype = "dashed",
        color = "#2c3e50",
        linewidth = 0.8
      )
    }
  }

  return(p)
}

#' Grafica el ranking comparativo de modelos segun una metrica de desempeno
#'
#' @param df_benchmarks Data frame de benchmark con la columna `Modelo`.
#' @param metric_col Metrica a graficar (ej. "F1_Score", "ROC_AUC", "MCC").
#' @param title Titulo del grafico.
#' @param high_color Color hexadecimal para el valor maximo.
#' @return Objeto `ggplot`.
#' @export
plot_model_benchmarks <- function(df_benchmarks,
                                  metric_col = "F1_Score",
                                  title = "Comparativa de Rendimiento de Modelos",
                                  high_color = "#27ae60") {
  stopifnot(is.data.frame(df_benchmarks))
  stopifnot("Modelo" %in% colnames(df_benchmarks))
  stopifnot(metric_col %in% colnames(df_benchmarks))

  df_sorted <- df_benchmarks %>%
    dplyr::filter(!is.na(.data[[metric_col]])) %>%
    dplyr::arrange(.data[[metric_col]]) %>%
    dplyr::mutate(Modelo = factor(Modelo, levels = unique(Modelo)))

  ggplot2::ggplot(df_sorted, ggplot2::aes(x = .data[[metric_col]], y = Modelo, fill = .data[[metric_col]])) +
    ggplot2::geom_col(width = 0.7, color = "#2c3e50", linewidth = 0.2) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.3f", .data[[metric_col]])),
      hjust = -0.1,
      size = 3.5,
      fontface = "bold"
    ) +
    ggplot2::scale_fill_gradient(low = "#d5f5e3", high = high_color) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(title = title, x = metric_col, y = "Modelo / Estrategia") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      legend.position = "none"
    )
}

#' Grafica la distribucion de veredictos del Sistema de Consenso de 7 Modelos
#'
#' @param df_consensus Data frame de salida de `consolidar_consenso()`.
#' @param verdict_col Nombre de la columna con el veredicto ("Veredicto_Consenso").
#' @param title Titulo del grafico.
#' @param subtitle Subtitulo del grafico.
#' @return Objeto `ggplot`.
#' @export
plot_consensus_distribution <- function(df_consensus,
                                        verdict_col = "Veredicto_Consenso",
                                        title = "Clasificacion de Secuencias por Consenso Multimodelo",
                                        subtitle = "Sistema de Consenso de 7 Modelos (Homologo Valido vs Ruido Descartado)") {
  stopifnot(is.data.frame(df_consensus))
  stopifnot(verdict_col %in% colnames(df_consensus))

  df_counts <- df_consensus %>%
    dplyr::count(.data[[verdict_col]], name = "Frecuencia")

  ggplot2::ggplot(df_counts, ggplot2::aes(x = .data[[verdict_col]], y = Frecuencia, fill = .data[[verdict_col]])) +
    ggplot2::geom_col(color = "black", alpha = 0.85, width = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = Frecuencia),
      vjust = -0.4,
      size = 4.5,
      fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "Homologo_Valido" = "#27ae60",
        "Ortologo_Valido" = "#27ae60",
        "Ruido_Descartado" = "#c0392b"
      )
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Veredicto del Ensamble de Consenso",
      y = "Numero de Secuencias",
      fill = "Categoria"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "#7f8c8d"),
      legend.position = "none"
    )
}

#' Grafica el Mapa de Calor (Heatmap) de Variacion en el Numero de Copias (CNV)
#'
#' Soporta tanto el data frame anotado (`df_cnv_annotated` de `extract_cnv_matrix()`)
#' como una matriz ancha (`cnv_matrix`), manteniendo compatibilidad universal.
#'
#' @param df_cnv_annotated Data frame anotado con `Clasificacion_Copia` (opcional).
#' @param cnv_matrix Data frame en formato ancho `Locus x Muestras` (opcional).
#' @param sample_col Columna de muestra en `df_cnv_annotated`.
#' @param locus_col Columna de locus/ortogrupo en `df_cnv_annotated`.
#' @param top_n_loci Numero maximo de loci mas abundantes a representar (default: 30).
#' @param title Titulo del mapa de calor.
#' @param subtitle Subtitulo del mapa de calor.
#' @return Objeto `ggplot`.
#' @export
plot_cnv_heatmap <- function(df_cnv_annotated = NULL,
                             cnv_matrix = NULL,
                             sample_col = "Sample",
                             locus_col = "Locus_Ortogrupo",
                             top_n_loci = 30,
                             title = "Matriz de Variacion en el Numero de Copias (CNV)",
                             subtitle = "Deteccion post-control de fragmentacion tecnica") {
  if (!is.null(df_cnv_annotated) && is.data.frame(df_cnv_annotated) &&
      all(c(sample_col, locus_col) %in% colnames(df_cnv_annotated))) {
    df_valid <- df_cnv_annotated
    if ("Clasificacion_Copia" %in% colnames(df_valid)) {
      df_valid <- df_valid %>%
        dplyr::filter(Clasificacion_Copia %in% c("Copia_Principal", "Paralogo_Confirmado"))
    }
    df_counts <- df_valid %>%
      dplyr::group_by(.data[[sample_col]], .data[[locus_col]]) %>%
      dplyr::summarise(Numero_Copias = dplyr::n(), .groups = "drop")

    top_loci <- df_counts %>%
      dplyr::group_by(.data[[locus_col]]) %>%
      dplyr::summarise(Total = sum(Numero_Copias), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(Total)) %>%
      dplyr::slice_head(n = top_n_loci) %>%
      dplyr::pull(.data[[locus_col]])

    df_long <- df_counts %>%
      dplyr::filter(.data[[locus_col]] %in% top_loci) %>%
      dplyr::rename(Muestra = dplyr::all_of(sample_col), Locus = dplyr::all_of(locus_col))
  } else if (!is.null(cnv_matrix) && is.data.frame(cnv_matrix)) {
    first_col <- colnames(cnv_matrix)[1]
    df_long <- cnv_matrix %>%
      dplyr::mutate(Total_Copias = rowSums(dplyr::across(-dplyr::all_of(first_col)), na.rm = TRUE)) %>%
      dplyr::arrange(dplyr::desc(Total_Copias)) %>%
      dplyr::slice_head(n = top_n_loci) %>%
      dplyr::select(-Total_Copias) %>%
      tidyr::pivot_longer(
        cols = -dplyr::all_of(first_col),
        names_to = "Muestra",
        values_to = "Numero_Copias"
      ) %>%
      dplyr::rename(Locus = dplyr::all_of(first_col))
  } else if (is.data.frame(df_cnv_annotated)) {
    # Caso en que se pase cnv_matrix como primer argumento posicional
    first_col <- colnames(df_cnv_annotated)[1]
    df_long <- df_cnv_annotated %>%
      dplyr::mutate(Total_Copias = rowSums(dplyr::across(-dplyr::all_of(first_col)), na.rm = TRUE)) %>%
      dplyr::arrange(dplyr::desc(Total_Copias)) %>%
      dplyr::slice_head(n = top_n_loci) %>%
      dplyr::select(-Total_Copias) %>%
      tidyr::pivot_longer(
        cols = -dplyr::all_of(first_col),
        names_to = "Muestra",
        values_to = "Numero_Copias"
      ) %>%
      dplyr::rename(Locus = dplyr::all_of(first_col))
  } else {
    stop("[ERROR] Debe proporcionar `df_cnv_annotated` o `cnv_matrix` a plot_cnv_heatmap().", call. = FALSE)
  }

  n_samples <- length(unique(df_long$Muestra))
  hide_x_labels <- n_samples > 25

  ggplot2::ggplot(df_long, ggplot2::aes(x = Muestra, y = Locus, fill = Numero_Copias)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.15) +
    ggplot2::scale_fill_viridis_c(option = "mako", direction = -1, name = "Copias") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Muestras Analizadas",
      y = "Locus / Ortogrupo"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "#7f8c8d"),
      axis.text.x = if (hide_x_labels) ggplot2::element_blank() else ggplot2::element_text(angle = 45, hjust = 1),
      axis.ticks.x = if (hide_x_labels) ggplot2::element_blank() else ggplot2::element_line(),
      panel.grid = ggplot2::element_blank()
    )
}

#' Grafica curvas ROC comparativas a partir de una lista de objetos pROC::roc
#'
#' @param roc_list Lista nombrada de objetos `roc`.
#' @param title Titulo del grafico.
#' @return Objeto `ggplot`.
#' @export
plot_roc_curves <- function(roc_list, title = "Curvas ROC Comparativas") {
  stopifnot(is.list(roc_list))

  df_roc_all <- data.frame()

  for (name in names(roc_list)) {
    roc_obj <- roc_list[[name]]
    if (inherits(roc_obj, "roc")) {
      df_roc <- data.frame(
        Especificidad = roc_obj$specificities,
        Sensibilidad = roc_obj$sensitivities,
        Modelo = paste0(name, " (AUC = ", sprintf("%.3f", as.numeric(roc_obj$auc)), ")"),
        stringsAsFactors = FALSE
      )
      df_roc_all <- dplyr::bind_rows(df_roc_all, df_roc)
    }
  }

  ggplot2::ggplot(df_roc_all, ggplot2::aes(x = 1 - Especificidad, y = Sensibilidad, color = Modelo)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::labs(
      title = title,
      x = "1 - Especificidad (Tasa de Falsos Positivos)",
      y = "Sensibilidad (Tasa de Verdaderos Positivos)"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      legend.position = "right"
    )
}
