#' Visualization and Diagnostic Plotting Module for TargetCNV-ML
#'
#' Pure functions wrapped in ggplot2 and pROC for generating confusion matrix tiles,
#' continuous metric density distributions, model benchmarks, and CNV heatmaps.
#'
#' @name plotting
NULL

#' Graficar Matriz de Confusion en Formato Tile
#'
#' Construye una visualizacion grafica de una matriz de confusion 2x2 utilizando
#' \code{ggplot2::geom_tile}, estandarizando las etiquetas a TP/FP con gradientes de color.
#'
#' @param cm_table Tabla o matriz de confusion (ej. salida de table() o confusionMatrix$table).
#' @param title Titulo del grafico (por defecto "Matriz de Confusion").
#' @param subtitle Subtitulo con metricas de resumen (ej. Accuracy, FDR).
#' @param high_color Color hexadecimal para la frecuencia maxima (por defecto "#d35400").
#'
#' @return Un objeto \code{ggplot2} listo para renderizar o combinar con cowplot.
#' @export
#' @examples
#' \dontrun{
#' plot_cm <- plot_confusion_matrix(cm$table, "ElasticNet", "Acc: 98% | FDR: 2%")
#' }
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

#' Grafico de Densidad para Separacion de Metricas Continuas
#'
#' Genera graficos facetados de densidad con areas transparentes para evaluar
#' zonas de solapamiento y capacidad de discriminacion entre clases (TP vs FP).
#'
#' @param df Data frame que contiene las columnas de predictores y respuesta.
#' @param response_col Nombre de la columna que indica la clase/categoria.
#' @param metric_cols Vector con los nombres de las columnas continuas a graficar.
#' @param colors Vector nombrado de colores para las clases (ej. c("TRUE_POSITIVE" = "#27ae60", ...)).
#' @param title Titulo del grafico.
#'
#' @return Un objeto \code{ggplot2}.
#' @export
plot_metric_densities <- function(df,
                                  response_col,
                                  metric_cols,
                                  colors = c("TRUE_POSITIVE" = "#27ae60", "FALSE_POSITIVE" = "#e74c3c",
                                             "TP" = "#27ae60", "FP" = "#e74c3c"),
                                  title = "Zonas de Solapamiento y Separacion (TP vs FP)") {
  stopifnot(is.data.frame(df))
  stopifnot(response_col %in% colnames(df))
  stopifnot(all(metric_cols %in% colnames(df)))

  df_long <- df %>%
    dplyr::select(dplyr::all_of(c(response_col, metric_cols))) %>%
    stats::na.omit() %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(metric_cols),
      names_to = "Metrica",
      values_to = "Valor"
    )

  resp_levels <- unique(as.character(df_long[[response_col]]))
  sub_colors <- colors[names(colors) %in% resp_levels]

  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Valor, fill = .data[[response_col]])) +
    ggplot2::geom_density(alpha = 0.6) +
    ggplot2::facet_wrap(~Metrica, scales = "free") +
    ggplot2::labs(title = title, x = "Valor de la Metrica", y = "Densidad", fill = "Clase") +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text = ggplot2::element_text(face = "bold")
    )

  if (length(sub_colors) == length(resp_levels)) {
    p <- p + ggplot2::scale_fill_manual(values = sub_colors)
  }

  return(p)
}

#' Grafico Descriptivo de Distribucion (Boxplot o Jitter)
#'
#' Genera comparativas visuales continuas estratificadas por clase para evaluar
#' distribuciones univariadas de ensamblaje.
#'
#' @param df Data frame con datos de secuencias.
#' @param response_col Columna con la clase o validacion espacial.
#' @param metric_cols Vector con las variables numericas.
#' @param type Tipo de representacion: "boxplot" o "jitter" (por defecto "boxplot").
#' @param title Titulo del grafico.
#' @param subtitle Subtitulo opcional.
#'
#' @return Un objeto \code{ggplot2}.
#' @export
plot_descriptive_distributions <- function(df,
                                           response_col,
                                           metric_cols,
                                           type = c("boxplot", "jitter"),
                                           title = "Caracterizacion de Secuencias: TP vs FP",
                                           subtitle = NULL) {
  type <- match.arg(type)
  stopifnot(is.data.frame(df))

  df_long <- df %>%
    dplyr::select(dplyr::all_of(c(response_col, metric_cols))) %>%
    stats::na.omit() %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(metric_cols),
      names_to = "Predictor",
      values_to = "Valor"
    )

  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = .data[[response_col]], y = Valor))

  if (type == "boxplot") {
    p <- p + ggplot2::geom_boxplot(
      ggplot2::aes(fill = .data[[response_col]]),
      alpha = 0.7,
      outlier.colour = "royalblue3",
      outlier.size = 1.5
    )
  } else {
    p <- p + ggplot2::geom_jitter(
      ggplot2::aes(colour = .data[[response_col]]),
      alpha = 0.7
    )
  }

  p <- p +
    ggplot2::facet_wrap(~Predictor, scales = "free_y", ncol = 2) +
    ggplot2::labs(title = title, subtitle = subtitle, x = "", y = "Valor Numerico") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      strip.text = ggplot2::element_text(face = "bold", size = 10),
      axis.text.x = ggplot2::element_text(angle = 15, hjust = 1)
    )

  return(p)
}

#' Benchmark Comparativo de Modelos y Estrategias
#'
#' Genera un grafico de barras horizontales agrupadas comparando multiples
#' modelos o filtros diagnosticos en funcion de sus metricas de rendimiento.
#'
#' @param df_benchmark Data frame que contiene la columna \code{Modelo} y metricas numericas.
#' @param model_col Columna que identifica cada estrategia/modelo (por defecto "Modelo").
#' @param metric_cols Metricas a graficar (por defecto Exactitud, Sensibilidad, Especificidad, Precision, FDR, F1_Score).
#' @param model_order Vector opcional con el orden deseado de los modelos.
#' @param title Titulo del benchmark.
#' @param subtitle Subtitulo explicativo.
#'
#' @return Un objeto \code{ggplot2}.
#' @export
plot_benchmark <- function(df_benchmark,
                           model_col = "Modelo",
                           metric_cols = c("Exactitud", "Sensibilidad", "Especificidad",
                                           "Precision", "FDR", "F1_Score"),
                           model_order = NULL,
                           title = "Benchmark de Estrategias Evaluadas",
                           subtitle = NULL) {
  stopifnot(is.data.frame(df_benchmark))
  stopifnot(model_col %in% colnames(df_benchmark))

  cols_validas <- intersect(metric_cols, colnames(df_benchmark))
  stopifnot(length(cols_validas) > 0)

  df_plot <- df_benchmark %>%
    dplyr::select(dplyr::all_of(c(model_col, cols_validas))) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(cols_validas),
      names_to = "Metrica",
      values_to = "Valor"
    )

  if (!is.null(model_order)) {
    df_plot[[model_col]] <- factor(df_plot[[model_col]], levels = rev(model_order))
  } else {
    df_plot[[model_col]] <- stats::reorder(factor(df_plot[[model_col]]), df_plot$Valor)
  }

  ggplot2::ggplot(df_plot, ggplot2::aes(x = .data[[model_col]], y = Valor, fill = Metrica)) +
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Estrategia Evaluada",
      y = "Puntaje (0 a 1)"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold")
    )
}

#' Curvas ROC Comparativas
#'
#' Grafica multiples curvas ROC sobrepuestas a partir de una lista de objetos \code{pROC::roc},
#' reportando el AUC de cada estrategia en la leyenda.
#'
#' @param roc_list Lista nombrada de objetos \code{roc}.
#' @param title Titulo del grafico (por defecto "Curvas ROC Comparativas").
#'
#' @return Un objeto \code{ggplot2}.
#' @export
plot_roc_comparison <- function(roc_list, title = "Curvas ROC Comparativas") {
  stopifnot(is.list(roc_list))
  stopifnot(length(roc_list) > 0)

  df_roc_total <- data.frame()

  for (mod_name in names(roc_list)) {
    roc_obj <- roc_list[[mod_name]]
    auc_val <- as.numeric(roc_obj$auc)
    etiqueta <- sprintf("%s (AUC = %.3f)", mod_name, auc_val)

    df_temp <- data.frame(
      Especificidad = roc_obj$specificities,
      Sensibilidad = roc_obj$sensitivities,
      Modelo = etiqueta,
      stringsAsFactors = FALSE
    )
    df_roc_total <- rbind(df_roc_total, df_temp)
  }

  ggplot2::ggplot(df_roc_total, ggplot2::aes(x = 1 - Especificidad, y = Sensibilidad, colour = Modelo)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::labs(
      title = title,
      x = "1 - Especificidad (Tasa de Falsos Positivos)",
      y = "Sensibilidad (Tasa de Verdaderos Positivos)",
      colour = "Modelo"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold")
    )
}

#' Heatmap de Variacion en el Numero de Copias (CNV)
#'
#' Genera un mapa de calor para visualizar la abundancia de copias genicas a traves
#' de las muestras para los top ortogrupos con mayor variabilidad o expansion.
#'
#' @param df_cnv_annotated Data frame anotado (salida de \code{extract_cnv_matrix$cnv_annotated}).
#' @param sample_col Columna de la muestra (por defecto "Sample").
#' @param locus_col Columna del locus o grupo homólogo (por defecto "Locus_Ortogrupo").
#' @param copy_classes Clases de copias a considerar (por defecto "Copia_Principal" y copia homóloga confirmada "Paralogo_Confirmado").
#' @param top_n_loci Numero de grupos homólogos más abundantes a mostrar (por defecto 30).
#' @param viridis_option Paleta viridis ("mako", "viridis", "plasma", etc., por defecto "mako").
#' @param title Titulo del grafico.
#' @param subtitle Subtitulo del grafico.
#'
#' @return Un objeto \code{ggplot2}.
#' @export
plot_cnv_heatmap <- function(df_cnv_annotated,
                             sample_col = "Sample",
                             locus_col = "Locus_Ortogrupo",
                             copy_classes = c("Copia_Principal", "Paralogo_Confirmado"),
                             top_n_loci = 30,
                             viridis_option = "mako",
                             title = "Matriz de Variacion Estructural (CNV) - Top Familias",
                             subtitle = "Deteccion post-control de fragmentacion tecnica") {
  stopifnot(is.data.frame(df_cnv_annotated))

  df_funcionales <- df_cnv_annotated %>%
    dplyr::filter(.data[["Clasificacion_Copia"]] %in% copy_classes)

  top_loci <- df_funcionales %>%
    dplyr::group_by(.data[[locus_col]]) %>%
    dplyr::summarise(Total_Copias = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(Total_Copias)) %>%
    dplyr::slice_head(n = top_n_loci) %>%
    dplyr::pull(.data[[locus_col]])

  df_heatmap <- df_funcionales %>%
    dplyr::filter(.data[[locus_col]] %in% top_loci) %>%
    dplyr::group_by(.data[[sample_col]], .data[[locus_col]]) %>%
    dplyr::summarise(Numero_Copias = dplyr::n(), .groups = "drop")

  ggplot2::ggplot(df_heatmap, ggplot2::aes(x = .data[[sample_col]], y = .data[[locus_col]], fill = Numero_Copias)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.1) +
    viridis::scale_fill_viridis(option = viridis_option, direction = -1) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Muestras de Target Capture",
      y = "Locus / Grupo Homólogo",
      fill = "Copias"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

#' Grafico de Barras de Distribucion del Veredicto de Consenso
#'
#' Visualiza el conteo total de secuencias clasificadas como homólogos válidos
#' vs "Ruido_Descartado" por el sistema de consenso multimodelo.
#'
#' @param df_consensus Data frame que contiene la columna del veredicto de consenso.
#' @param verdict_col Columna con la clasificacion (por defecto "Veredicto_Consenso").
#' @param colors Vector nombrado de colores para las categorias.
#' @param title Titulo del grafico.
#' @param subtitle Subtitulo del grafico.
#'
#' @return Un objeto \code{ggplot2}.
#' @export
plot_consensus_distribution <- function(df_consensus,
                                        verdict_col = "Veredicto_Consenso",
                                        colors = c("Homologo_Valido" = "#27ae60",
                                                   "Ortologo_Valido" = "#27ae60",
                                                   "Ruido_Descartado" = "#c0392b"),
                                        title = "Clasificacion de Secuencias por Sistema de Consenso",
                                        subtitle = "Mayoría Simple (>= 4/6 Modelos)") {
  stopifnot(is.data.frame(df_consensus))
  stopifnot(verdict_col %in% colnames(df_consensus))

  ggplot2::ggplot(df_consensus, ggplot2::aes(x = .data[[verdict_col]], fill = .data[[verdict_col]])) +
    ggplot2::geom_bar(color = "black", alpha = 0.8) +
    ggplot2::geom_text(stat = "count", ggplot2::aes(label = ggplot2::after_stat(count)), vjust = -0.5, size = 5) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Veredicto del Ensamble",
      y = "Numero Total de Secuencias"
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold")
    )
}
