# Funciones de visualización diagnóstica con ggplot2 para TargetCNV-ML

#' Paleta de colores distintiva por modelo para matrices de confusión e importancia de variables
#'
#' @param model_name Nombre del modelo o estrategia evaluada.
#' @return Lista con `high` (color principal) y `low` (color suave de fondo).
#' @export
get_model_color_palette <- function(model_name) {
  palettes <- list(
    ElasticNet                  = list(high = "#16a085", low = "#e8f8f5"),
    NaiveBayes                  = list(high = "#2980b9", low = "#ebf5fb"),
    SVM_Polinomial              = list(high = "#d35400", low = "#fdf2e9"),
    Red_Neuronal_Bag            = list(high = "#884ea0", low = "#f5eef8"),
    Red_Neuronal                = list(high = "#c0392b", low = "#fdedec"),
    LDA                         = list(high = "#34495e", low = "#ebedef"),
    Weighted_Heuristic          = list(high = "#8e44ad", low = "#f4ecf7"),
    Modelo_Heuristico           = list(high = "#8e44ad", low = "#f4ecf7"),
    Consenso_7_Modelos          = list(high = "#117864", low = "#e8f8f5"),
    Youden_pct_recovered        = list(high = "#27ae60", low = "#eafaf1"),
    Youden_pct_identity         = list(high = "#229954", low = "#e9f7ef"),
    Youden_score                = list(high = "#138d75", low = "#e8f8f5"),
    Youden_wscore               = list(high = "#117a65", low = "#e8f6f3"),
    Youden                      = list(high = "#27ae60", low = "#eafaf1"),
    RandomForest                = list(high = "#1e8449", low = "#e9f7ef"),
    GBM                         = list(high = "#b9770e", low = "#fef9e7"),
    XGBoost                     = list(high = "#6c3483", low = "#f5eef8"),
    SVM_Radial                  = list(high = "#a04000", low = "#fbeee6"),
    SVM_Lineal                  = list(high = "#ca6f1e", low = "#fef5e7"),
    KNN                         = list(high = "#2471a3", low = "#eaf2f8"),
    CART                        = list(high = "#b7950b", low = "#fef9e7"),
    MLP_Regularizado            = list(high = "#7b241c", low = "#f9ebea"),
    MLP                         = list(high = "#922b21", low = "#f9ebea"),
    DNN_Profunda                = list(high = "#922b21", low = "#f9ebea"),
    DNN_Regularizada            = list(high = "#7b241c", low = "#f9ebea"),
    QDA                         = list(high = "#5d6d7e", low = "#ebedef"),
    Bagging_Arboles             = list(high = "#148f77", low = "#e8f8f5"),
    `Stacked Autoencoder - DNN` = list(high = "#512e5f", low = "#f4ecf7"),
    Autoencoder                 = list(high = "#512e5f", low = "#f4ecf7"),
    DNN_Autoencoder             = list(high = "#512e5f", low = "#f4ecf7"),
    Ensamble_Stacking           = list(high = "#5b2c6f", low = "#f4ecf7")
  )

  clean_name <- trimws(as.character(model_name))
  if (clean_name %in% names(palettes)) {
    return(palettes[[clean_name]])
  }
  keys_by_len <- names(palettes)[order(nchar(names(palettes)), decreasing = TRUE)]
  for (key in keys_by_len) {
    if (grepl(key, clean_name, ignore.case = TRUE)) {
      return(palettes[[key]])
    }
  }
  list(high = "#2980b9", low = "#f1f2f6")
}

#' Grafica una matriz de confusión 2x2 con TP en la esquina superior izquierda
#'
#' @param cm_table Tabla o matriz 2x2 de confusión.
#' @param title Título principal del gráfico.
#' @param subtitle Subtítulo con métricas resumidas.
#' @param high_color Color hexadecimal para la celda de mayor frecuencia.
#' @param low_color Color hexadecimal para las celdas de menor frecuencia.
#' @return Objeto `ggplot`.
#' @export
plot_confusion_matrix <- function(cm_table,
                                  title = "Matriz de Confusión",
                                  subtitle = NULL,
                                  high_color = NULL,
                                  low_color = NULL) {
  stopifnot(is.table(cm_table) || is.matrix(cm_table))

  pal <- get_model_color_palette(title)
  if (is.null(high_color)) high_color <- pal$high
  if (is.null(low_color))  low_color  <- pal$low

  cm_df <- as.data.frame(as.table(cm_table))
  colnames(cm_df) <- c("Prediccion", "Referencia", "Frecuencia")

  cm_df <- cm_df %>%
    dplyr::mutate(
      Referencia = dplyr::if_else(grepl("TRUE|TP", Referencia, ignore.case = TRUE), "TP", "FP"),
      Prediccion = dplyr::if_else(grepl("TRUE|TP", Prediccion, ignore.case = TRUE), "TP", "FP"),
      Referencia = factor(Referencia, levels = c("TP", "FP")),
      Prediccion = factor(Prediccion, levels = c("FP", "TP")),
      Etiqueta_Celda = dplyr::case_when(
        Referencia == "TP" & Prediccion == "TP" ~ "TP (Verdadero Positivo)",
        Referencia == "FP" & Prediccion == "TP" ~ "FP (Falso Positivo)",
        Referencia == "TP" & Prediccion == "FP" ~ "FN (Falso Negativo)",
        TRUE ~ "TN (Verdadero Negativo)"
      )
    )

  max_freq <- max(cm_df$Frecuencia, na.rm = TRUE)
  if (is.na(max_freq) || max_freq == 0) max_freq <- 1

  ggplot2::ggplot(cm_df, ggplot2::aes(x = Referencia, y = Prediccion)) +
    ggplot2::geom_tile(ggplot2::aes(fill = Frecuencia), color = "white", linewidth = 1.2) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%d\n(%s)", round(Frecuencia, 0), Etiqueta_Celda)),
      size = 4.0,
      fontface = "bold",
      color = ifelse(cm_df$Frecuencia > max_freq / 2, "white", "#2c3e50")
    ) +
    ggplot2::scale_fill_gradient(low = low_color, high = high_color, name = "Secuencias") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Clase Real (Referencia)",
      y = "Clase Predicha"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 11.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "#5d6d7e", size = 9.5),
      axis.title = ggplot2::element_text(face = "bold", size = 10),
      axis.text = ggplot2::element_text(face = "bold", size = 10),
      legend.position = "none",
      panel.grid = ggplot2::element_blank()
    )
}

#' Calcula y grafica la importancia relativa de variables de un clasificador
#'
#' @param model_obj Objeto `train` de `caret` o vector numérico nombrado de pesos.
#' @param df_data Data frame opcional con los predictores para respaldo por permutación.
#' @param model_name Nombre del modelo para el título y paleta de color.
#' @param high_color Color principal opcional.
#' @return Objeto `ggplot`.
#' @export
plot_model_varimp <- function(model_obj,
                              df_data = NULL,
                              model_name = "Modelo",
                              high_color = NULL) {
  pal <- get_model_color_palette(model_name)
  if (is.null(high_color)) high_color <- pal$high

  imp_df <- NULL

  if (is.numeric(model_obj) && !is.null(names(model_obj))) {
    imp_df <- data.frame(
      Variable = names(model_obj),
      Importancia = as.numeric(model_obj),
      stringsAsFactors = FALSE
    )
  } else {
    imp_df <- tryCatch({
      vi <- caret::varImp(model_obj, scale = TRUE)$importance
      col_pick <- if ("Overall" %in% colnames(vi)) "Overall" else if ("TP" %in% colnames(vi)) "TP" else colnames(vi)[1]
      data.frame(
        Variable = rownames(vi),
        Importancia = as.numeric(vi[[col_pick]]),
        stringsAsFactors = FALSE
      )
    }, error = function(e) NULL)
  }

  if (is.null(imp_df) || nrow(imp_df) == 0 || all(is.na(imp_df$Importancia))) {
    train_x <- if (is.list(model_obj) && !is.null(model_obj$trainingData)) {
      model_obj$trainingData[, setdiff(colnames(model_obj$trainingData), ".outcome"), drop = FALSE]
    } else if (!is.null(df_data)) {
      df_data[, intersect(c("pct_recovered", "pct_identity", "score", "wscore", "hit"), colnames(df_data)), drop = FALSE]
    } else {
      NULL
    }

    if (!is.null(train_x) && ncol(train_x) > 0) {
      base_probs <- tryCatch(stats::predict(model_obj, newdata = train_x, type = "prob")[, "TP"], error = function(e) NULL)
      if (!is.null(base_probs)) {
        set.seed(1234)
        scores <- sapply(colnames(train_x), function(col_v) {
          diffs <- replicate(5, {
            x_perm <- train_x
            x_perm[[col_v]] <- sample(train_x[[col_v]])
            perm_probs <- stats::predict(model_obj, newdata = x_perm, type = "prob")[, "TP"]
            mean(abs(base_probs - perm_probs), na.rm = TRUE)
          })
          mean(diffs, na.rm = TRUE)
        })
        max_s <- max(scores, na.rm = TRUE)
        if (max_s > 0) scores <- (scores / max_s) * 100
        imp_df <- data.frame(
          Variable = names(scores),
          Importancia = as.numeric(scores),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (is.null(imp_df) || nrow(imp_df) == 0) return(NULL)

  imp_df <- imp_df %>%
    dplyr::mutate(
      Importancia = pmax(0, ifelse(is.na(Importancia), 0, Importancia))
    ) %>%
    dplyr::arrange(Importancia) %>%
    dplyr::mutate(Variable = factor(Variable, levels = unique(Variable)))

  ggplot2::ggplot(imp_df, ggplot2::aes(x = Importancia, y = Variable, fill = Importancia)) +
    ggplot2::geom_col(width = 0.65, color = "#2c3e50", linewidth = 0.25) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.1f", Importancia)),
      hjust = -0.15,
      size = 3.8,
      fontface = "bold",
      color = "#2c3e50"
    ) +
    ggplot2::scale_fill_gradient(low = pal$low, high = high_color) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.2))) +
    ggplot2::labs(
      title = paste("Importancia de Variables:", model_name),
      subtitle = "Peso relativo en la decisión del modelo",
      x = "Importancia Relativa",
      y = "Variable Predictora"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 11.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "#5d6d7e", size = 9.5),
      axis.title = ggplot2::element_text(face = "bold", size = 10),
      axis.text.y = ggplot2::element_text(face = "bold", size = 10),
      legend.position = "none"
    )
}

#' Dibuja la topología de redes neuronales (`nnet`, `avNNet`, `MLP`, `Autoencoder`)
#'
#' @param model_obj Objeto `train` de `caret`.
#' @param title Título del diagrama de arquitectura.
#' @export
plot_neural_net_topology <- function(model_obj,
                                     title = "Arquitectura de la Red Neuronal") {
  if (is.null(model_obj) || is.null(model_obj$finalModel)) {
    return(invisible(NULL))
  }

  fm <- model_obj$finalModel
  bt <- model_obj$bestTune

  if (inherits(fm, "avNNet") && !is.null(fm$model) && length(fm$model) >= 1) {
    if (requireNamespace("NeuralNetTools", quietly = TRUE)) {
      tryCatch({
        NeuralNetTools::plotnet(
          fm$model[[1]],
          alpha_val = 0.75, circle_cex = 4.5, cex_val = 0.85, max_sp = TRUE,
          pos_col = "#16a085", neg_col = "#884ea0", circle_col = "#f5eef8", bord_col = "#884ea0"
        )
        graphics::title(main = title, font.main = 2, cex.main = 1.0)
        return(invisible(NULL))
      }, error = function(e) {})
    }
  }

  if (inherits(fm, "nnet")) {
    if (requireNamespace("NeuralNetTools", quietly = TRUE)) {
      tryCatch({
        NeuralNetTools::plotnet(
          fm,
          alpha_val = 0.75, circle_cex = 4.5, cex_val = 0.85, max_sp = TRUE,
          pos_col = "#16a085", neg_col = "#c0392b", circle_col = "#fdedec", bord_col = "#c0392b"
        )
        graphics::title(main = title, font.main = 2, cex.main = 1.0)
        return(invisible(NULL))
      }, error = function(e) {})
    }
  }

  if (!is.null(bt) && any(c("layer1", "size") %in% colnames(bt))) {
    in_vars <- if (!is.null(model_obj$coefnames)) model_obj$coefnames else c("pct_recovered", "pct_identity", "score", "wscore")
    n_in <- length(in_vars)

    hidden_sizes <- numeric(0)
    if ("size" %in% colnames(bt)) {
      hidden_sizes <- as.numeric(bt$size[1])
    } else {
      for (l_col in c("layer1", "layer2", "layer3")) {
        if (l_col %in% colnames(bt) && !is.na(bt[[l_col]][1]) && bt[[l_col]][1] > 0) {
          hidden_sizes <- c(hidden_sizes, as.numeric(bt[[l_col]][1]))
        }
      }
    }

    layer_names <- c(
      sprintf("Entrada\n(%d vars)", n_in),
      sapply(seq_along(hidden_sizes), function(i) sprintf("Capa Oculta %d\n(%d neuronas)", i, hidden_sizes[i])),
      "Salida\n(TP / FP)"
    )
    layer_counts <- c(n_in, hidden_sizes, 2)
    vis_counts <- pmin(layer_counts, 12)

    nodes_list <- list()
    for (lx in seq_along(vis_counts)) {
      k <- vis_counts[lx]
      y_pos <- if (k == 1) 0 else seq(-1, 1, length.out = k)
      lbls <- if (lx == 1 && k == length(in_vars)) {
        in_vars
      } else if (lx == length(vis_counts)) {
        c("TP", "FP")
      } else {
        paste0("H", lx - 1, "_", seq_len(k))
      }
      nodes_list[[lx]] <- data.frame(
        layer_id = lx,
        layer_label = layer_names[lx],
        x = lx,
        y = y_pos,
        node_label = lbls,
        tipo = if (lx == 1) "Entrada" else if (lx == length(vis_counts)) "Salida" else "Oculta",
        stringsAsFactors = FALSE
      )
    }
    df_nodes <- dplyr::bind_rows(nodes_list)

    edges_list <- list()
    for (lx in seq_len(length(vis_counts) - 1)) {
      n1 <- nodes_list[[lx]]
      n2 <- nodes_list[[lx + 1]]
      grid_e <- expand.grid(i = seq_len(nrow(n1)), j = seq_len(nrow(n2)))
      edges_list[[lx]] <- data.frame(
        x = n1$x[grid_e$i], y = n1$y[grid_e$i],
        xend = n2$x[grid_e$j], yend = n2$y[grid_e$j]
      )
    }
    df_edges <- dplyr::bind_rows(edges_list)

    p_arch <- ggplot2::ggplot() +
      ggplot2::geom_segment(
        data = df_edges,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
        color = "#bdc3c7", alpha = 0.55, linewidth = 0.45
      ) +
      ggplot2::geom_point(
        data = df_nodes,
        ggplot2::aes(x = x, y = y, fill = tipo),
        shape = 21, size = 9, color = "#2c3e50", stroke = 0.9
      ) +
      ggplot2::geom_text(
        data = df_nodes,
        ggplot2::aes(x = x, y = y, label = node_label),
        size = 2.7, fontface = "bold", color = "#2c3e50"
      ) +
      ggplot2::scale_fill_manual(values = c("Entrada" = "#d6eaf8", "Oculta" = "#fdedec", "Salida" = "#d5f5e3")) +
      ggplot2::scale_x_continuous(
        breaks = seq_along(layer_names),
        labels = layer_names,
        expand = ggplot2::expansion(mult = c(0.15, 0.15))
      ) +
      ggplot2::scale_y_continuous(limits = c(-1.25, 1.25)) +
      ggplot2::labs(
        title = title,
        subtitle = sprintf("Estructura de capas: %s nodos", paste(layer_counts, collapse = " -> ")),
        x = "", y = ""
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "#5d6d7e", size = 10),
        axis.text.x = ggplot2::element_text(face = "bold", size = 10, color = "#2c3e50"),
        axis.text.y = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        legend.position = "none"
      )
    print(p_arch)
  }
}

#' Grafica las distribuciones de densidad de las métricas continuas (TP vs. FP)
#'
#' @param df Data frame con las métricas y la columna de clase.
#' @param response_col Nombre de la columna de respuesta (`"Validacion_Espacial"` o `"Clase"`).
#' @param metric_cols Vector con los nombres de las métricas continuas.
#' @param thresholds Vector nombrado opcional con los puntos de corte de Youden.
#' @param colors Vector nombrado de colores por clase.
#' @param title Título del gráfico.
#' @param subtitle Subtítulo del gráfico.
#' @param legend_position Posición de la leyenda (`"bottom"`, `"right"`, `"none"`).
#' @return Objeto `ggplot`.
#' @export
plot_metric_densities <- function(df,
                                  response_col = "Clase",
                                  metric_cols = c("pct_recovered", "pct_identity", "score", "wscore"),
                                  thresholds = NULL,
                                  colors = c("TRUE_POSITIVE" = "#27ae60", "FALSE_POSITIVE" = "#e74c3c",
                                             "TP" = "#27ae60", "FP" = "#e74c3c"),
                                  title = "Separación entre Verdaderos Positivos (TP) y Falsos Positivos (FP)",
                                  subtitle = "Distribución de densidad y puntos de corte de Youden (línea punteada)",
                                  legend_position = "bottom",
                                  class_col = NULL) {
  if (!is.null(class_col) && nzchar(class_col)) response_col <- class_col
  stopifnot(is.data.frame(df))
  stopifnot(response_col %in% colnames(df))
  stopifnot(all(metric_cols %in% colnames(df)))

  df_long <- df %>%
    dplyr::select(dplyr::all_of(c(response_col, metric_cols))) %>%
    dplyr::mutate(
      Clase_Simple = ifelse(.data[[response_col]] %in% c("TRUE_POSITIVE", "TP"), "TP", "FP")
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(metric_cols),
      names_to = "Metrica",
      values_to = "Valor"
    )

  p <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(x = Valor, fill = Clase_Simple, color = Clase_Simple)
  ) +
    ggplot2::geom_density(alpha = 0.45, linewidth = 0.8) +
    ggplot2::facet_wrap(~ Metrica, scales = "free") +
    ggplot2::scale_fill_manual(values = c("TP" = "#27ae60", "FP" = "#e74c3c")) +
    ggplot2::scale_color_manual(values = c("TP" = "#27ae60", "FP" = "#e74c3c")) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Valor de la Métrica",
      y = "Densidad",
      fill = "Clase",
      color = "Clase"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "#5d6d7e"),
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

  p
}

#' Grafica el mapa de calor de colinealidad por rangos de Spearman
#'
#' @param cor_matrix Matriz simétrica de correlaciones de Spearman (`rho`), o lista devuelta por `compute_spearman_collinearity()`.
#' @param p_matrix Matriz opcional de p-valores asociados.
#' @param title Título del mapa de calor.
#' @return Objeto `ggplot`.
#' @export
plot_spearman_collinearity <- function(cor_matrix,
                                       p_matrix = NULL,
                                       title = "Matriz de Colinealidad por Rangos de Spearman") {
  if (is.list(cor_matrix) && !is.data.frame(cor_matrix) && !is.null(cor_matrix$cor_matrix)) {
    if (is.null(p_matrix) && !is.null(cor_matrix$p_matrix)) {
      p_matrix <- cor_matrix$p_matrix
    }
    cor_matrix <- cor_matrix$cor_matrix
  }
  stopifnot(is.matrix(cor_matrix) || is.data.frame(cor_matrix))
  cor_mat <- as.matrix(cor_matrix)
  vars <- colnames(cor_mat)

  grid_df <- expand.grid(Var1 = vars, Var2 = vars, stringsAsFactors = FALSE)
  grid_df$Rho <- as.numeric(cor_mat)
  grid_df$P_Val <- if (!is.null(p_matrix)) as.numeric(as.matrix(p_matrix)) else NA_real_

  grid_df <- grid_df %>%
    dplyr::mutate(
      Var1 = factor(Var1, levels = vars),
      Var2 = factor(Var2, levels = rev(vars)),
      Etiqueta = dplyr::case_when(
        as.character(Var1) == as.character(Var2) ~ "1.000",
        is.na(Rho) ~ "NA",
        !is.na(P_Val) & P_Val < 0.001 ~ sprintf("rho = %.3f\n(p < 0.001)", Rho),
        !is.na(P_Val) ~ sprintf("rho = %.3f\n(p = %.3f)", Rho, P_Val),
        TRUE ~ sprintf("rho = %.3f", Rho)
      )
    )

  ggplot2::ggplot(grid_df, ggplot2::aes(x = Var1, y = Var2, fill = Rho)) +
    ggplot2::geom_tile(color = "white", linewidth = 1.0) +
    ggplot2::geom_text(
      ggplot2::aes(label = Etiqueta),
      size = 3.7,
      fontface = "bold",
      color = ifelse(abs(grid_df$Rho) > 0.65, "white", "#2c3e50")
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#c0392b",
      mid = "#f8f9f9",
      high = "#1f618d",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Spearman rho"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = "Diagnóstico no paramétrico de dependencia monótona entre predictores continuos",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "#5d6d7e", size = 9.8),
      axis.text = ggplot2::element_text(face = "bold", size = 10.5, color = "#2c3e50"),
      panel.grid = ggplot2::element_blank()
    )
}

#' Grafica la curva de sensibilidad del umbral de cobertura (`--cov_threshold`, 30% a 80%)
#'
#' Representa las métricas de concordancia con verdad de referencia (`Concordancia_Exacta`,
#' `Kappa_Cuadratico`, `Spearman_Rho`, `RMSE`) cuando están presentes, o bien la estabilidad
#' del conteo de copias estimadas frente al umbral en muestras no etiquetadas.
#'
#' @param df_sens Data frame generado por `evaluate_cnv_threshold_sensitivity()`.
#' @param active_threshold Umbral de cobertura operativo a resaltar (por defecto `50.0`).
#' @param title Título principal del gráfico.
#' @return Objeto `ggplot`.
#' @export
plot_cnv_threshold_sensitivity <- function(df_sens,
                                           active_threshold = 50.0,
                                           title = "Sensibilidad del Conteo de Copias al Umbral de Cobertura (30% - 80%)") {
  stopifnot(is.data.frame(df_sens), "Umbral_Cobertura" %in% colnames(df_sens))

  has_species <- "Especie" %in% colnames(df_sens)
  has_gt <- all(c("Concordancia_Exacta", "RMSE") %in% colnames(df_sens))

  if (has_gt) {
    metric_map <- c(
      "Concordancia_Exacta" = "Concordancia Exacta",
      "Kappa_Cuadratico"    = "Kappa Cuadrático Ponderado",
      "Spearman_Rho"        = "Spearman rho",
      "RMSE"                = "RMSE (Copias)"
    )
    avail_m <- intersect(names(metric_map), colnames(df_sens))
    keep_cols <- unique(c(if (has_species) "Especie", "Umbral_Cobertura", avail_m))

    df_long <- df_sens[, keep_cols, drop = FALSE] %>%
      tidyr::pivot_longer(cols = dplyr::all_of(avail_m), names_to = "Metrica", values_to = "Valor") %>%
      dplyr::filter(!is.na(Valor)) %>%
      dplyr::mutate(Metrica_Label = factor(metric_map[Metrica], levels = unname(metric_map[avail_m])))

    p <- if (has_species) {
      ggplot2::ggplot(df_long, ggplot2::aes(x = Umbral_Cobertura, y = Valor, color = Especie, group = Especie)) +
        ggplot2::geom_vline(xintercept = active_threshold, linetype = "dashed", color = "#7f8c8d", linewidth = 0.8) +
        ggplot2::geom_line(linewidth = 1.0) +
        ggplot2::geom_point(size = 2.4)
    } else {
      ggplot2::ggplot(df_long, ggplot2::aes(x = Umbral_Cobertura, y = Valor)) +
        ggplot2::geom_vline(xintercept = active_threshold, linetype = "dashed", color = "#7f8c8d", linewidth = 0.8) +
        ggplot2::geom_line(linewidth = 1.0, color = "#1f618d") +
        ggplot2::geom_point(size = 2.4, color = "#1f618d")
    }

    p <- p +
      ggplot2::facet_wrap(~ Metrica_Label, scales = "free_y") +
      ggplot2::scale_x_continuous(breaks = seq(30, 80, by = 10)) +
      ggplot2::labs(
        title = title,
        subtitle = sprintf("Línea vertical punteada: umbral operativo de %.0f%%", active_threshold),
        x = "Umbral de Cobertura para Copias Adicionales (%)",
        y = "Valor de la Métrica",
        color = if (has_species) "Especie" else NULL
      ) +
      ggplot2::theme_bw(base_size = 11.5) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "#5d6d7e"),
        strip.text = ggplot2::element_text(face = "bold"),
        legend.position = if (has_species) "bottom" else "none"
      )
    return(p)
  }

  df_long <- df_sens %>%
    dplyr::select(Umbral_Cobertura, Copias_Totales_Estimadas, Homologos_Confirmados, Fragmentos_No_Contados) %>%
    tidyr::pivot_longer(
      cols = c("Copias_Totales_Estimadas", "Homologos_Confirmados", "Fragmentos_No_Contados"),
      names_to = "Categoria",
      values_to = "Conteo"
    ) %>%
    dplyr::mutate(
      Categoria = dplyr::recode(
        Categoria,
        "Copias_Totales_Estimadas" = "Copias Totales Estimadas",
        "Homologos_Confirmados"    = "Homólogos Adicionales (>= umbral)",
        "Fragmentos_No_Contados"   = "Fragmentos No Contados (< umbral)"
      )
    )

  ggplot2::ggplot(df_long, ggplot2::aes(x = Umbral_Cobertura, y = Conteo, color = Categoria)) +
    ggplot2::geom_vline(xintercept = active_threshold, linetype = "dashed", color = "#7f8c8d", linewidth = 0.8) +
    ggplot2::geom_line(linewidth = 1.0) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::scale_color_manual(values = c(
      "Copias Totales Estimadas"          = "#2c3e50",
      "Homólogos Adicionales (>= umbral)" = "#27ae60",
      "Fragmentos No Contados (< umbral)" = "#d35400"
    )) +
    ggplot2::scale_x_continuous(breaks = seq(30, 80, by = 5)) +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("Estabilidad empírica del conteo frente al umbral (umbral activo: %.0f%%)", active_threshold),
      x = "Umbral de Cobertura Individual (%)",
      y = "Número de Secuencias / Copias",
      color = "Categoría"
    ) +
    ggplot2::theme_bw(base_size = 11.5) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "#5d6d7e"),
      legend.position = "bottom"
    )
}

#' Grafica la distribución híbrida (densidad, boxplot y puntos) de puntajes en transferencia
#'
#' @param df_sc Data frame con columnas `Clase_Real` (`"TP"`/`"FP"`) y `Puntaje`.
#' @param mod_name Nombre del modelo graficado.
#' @param sp_name Nombre corto de la especie receptora.
#' @return Objeto `ggplot`.
#' @export
plot_transfer_raincloud <- function(df_sc, mod_name, sp_name) {
  df_sc <- df_sc %>%
    dplyr::filter(!is.na(Puntaje) & !is.na(Clase_Real)) %>%
    dplyr::mutate(
      Clase_Factor = factor(as.character(Clase_Real), levels = c("FP", "TP")),
      x_num = as.numeric(Clase_Factor)
    )

  if (nrow(df_sc) < 2) return(NULL)

  pal_base <- c("FP" = "#c0392b", "TP" = "#27ae60")
  has_cs   <- requireNamespace("colorspace", quietly = TRUE)

  col_dark <- if (has_cs) {
    stats::setNames(colorspace::darken(pal_base, 0.12, space = "HLS"), names(pal_base))
  } else {
    c("FP" = "#922b21", "TP" = "#1e8449")
  }
  col_fill_dens <- if (has_cs) {
    stats::setNames(colorspace::lighten(pal_base, 0.48), names(pal_base))
  } else {
    stats::setNames(grDevices::adjustcolor(pal_base, alpha.f = 0.55), names(pal_base))
  }
  col_fill_box <- if (has_cs) {
    stats::setNames(colorspace::desaturate(colorspace::lighten(pal_base, 0.80), 0.40), names(pal_base))
  } else {
    stats::setNames(grDevices::adjustcolor(pal_base, alpha.f = 0.22), names(pal_base))
  }

  y_min  <- min(df_sc$Puntaje, na.rm = TRUE)
  y_max  <- max(df_sc$Puntaje, na.rm = TRUE)
  y_span <- if ((y_max - y_min) > 1e-6) (y_max - y_min) else 1.0

  dens_list <- lapply(c("FP", "TP"), function(cl) {
    sub_v <- df_sc$Puntaje[df_sc$Clase_Factor == cl]
    if (length(sub_v) == 0) return(NULL)
    x_b <- if (cl == "FP") 1.0 else 2.0
    bw_val <- if (length(sub_v) > 1 && stats::sd(sub_v) > 1e-6) {
      max(stats::bw.nrd0(sub_v) * 0.65, 0.025 * y_span)
    } else {
      0.04 * y_span
    }
    d_obj <- stats::density(
      sub_v, bw = bw_val, n = 128,
      from = y_min - 0.02 * y_span, to = y_max + 0.02 * y_span
    )
    d_norm <- if (max(d_obj$y) > 0) (d_obj$y / max(d_obj$y)) * 0.36 else rep(0, length(d_obj$y))
    data.frame(
      Clase_Factor = factor(cl, levels = c("FP", "TP")),
      x_poly = c(rep(x_b + 0.15, length(d_obj$x)), rev(x_b + 0.15 + d_norm)),
      y_poly = c(d_obj$x, rev(d_obj$x))
    )
  })
  df_dens <- dplyr::bind_rows(dens_list)

  df_stats <- df_sc %>%
    dplyr::group_by(Clase_Factor, x_num) %>%
    dplyr::summarise(
      mediana = stats::median(Puntaje, na.rm = TRUE),
      n_obs = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      med_label = sprintf("%.2f", mediana),
      n_label = paste("n =", n_obs),
      y_n_pos = y_max + 0.045 * y_span
    )

  set.seed(1)
  df_sc$x_rain <- df_sc$x_num - stats::runif(nrow(df_sc), min = 0.14, max = 0.34)

  ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = df_dens,
      ggplot2::aes(x = x_poly, y = y_poly, group = Clase_Factor, fill = Clase_Factor, color = Clase_Factor),
      linewidth = 0.55, alpha = 0.85
    ) +
    ggplot2::geom_boxplot(
      data = df_sc,
      ggplot2::aes(x = x_num, y = Puntaje, group = Clase_Factor, color = Clase_Factor),
      fill = col_fill_box[as.character(levels(df_sc$Clase_Factor))[sort(unique(df_sc$x_num))]],
      width = 0.20, outlier.shape = NA, linewidth = 0.65
    ) +
    ggplot2::geom_point(
      data = df_sc,
      ggplot2::aes(x = x_rain, y = Puntaje, color = Clase_Factor),
      fill = "white", shape = 21, stroke = 0.42, size = 2.0
    ) +
    ggplot2::geom_point(
      data = df_sc,
      ggplot2::aes(x = x_rain, y = Puntaje, fill = Clase_Factor),
      color = "transparent", shape = 21, stroke = 0.42, size = 2.0, alpha = 0.35
    ) +
    ggplot2::geom_text(
      data = df_stats,
      ggplot2::aes(x = x_num + 0.17, y = mediana, label = med_label, color = Clase_Factor),
      fontface = "bold", size = 3.6, vjust = -0.3
    ) +
    ggplot2::geom_text(
      data = df_stats,
      ggplot2::aes(x = x_num, y = y_n_pos, label = n_label, color = Clase_Factor),
      fontface = "bold", size = 3.5, hjust = 0
    ) +
    ggplot2::coord_flip(
      xlim = c(0.52, 2.62),
      ylim = c(y_min - 0.05 * y_span, y_max + 0.20 * y_span),
      clip = "off"
    ) +
    ggplot2::scale_x_continuous(breaks = c(1, 2), labels = c("FP", "TP")) +
    ggplot2::scale_color_manual(values = col_dark, guide = "none") +
    ggplot2::scale_fill_manual(values = col_fill_dens, guide = "none") +
    ggplot2::labs(
      title = paste0("Separación de Puntaje/Probabilidad: ", mod_name),
      subtitle = paste0("Distribución en ", sp_name, " (TP vs. FP)"),
      x = NULL,
      y = "Probabilidad P(TP) / Puntaje"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(color = col_dark[c("FP", "TP")], face = "bold", size = 13),
      axis.title.x = ggplot2::element_text(face = "bold", size = 10.5, margin = ggplot2::margin(t = 8)),
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle = ggplot2::element_text(color = "#5d6d7e", size = 9.5, margin = ggplot2::margin(0, 0, 10, 0)),
      plot.title.position = "plot",
      plot.margin = ggplot2::margin(12, 18, 10, 12)
    )
}

#' Grafica el ranking comparativo de modelos para una o múltiples métricas
#'
#' @param df_benchmarks Data frame de benchmark con la columna `Modelo`.
#' @param metric_col Métrica única a graficar (ej. `"F1_Score"`).
#' @param metric_cols Vector opcional con varias métricas para vista general.
#' @param title Título del gráfico.
#' @param high_color Color hexadecimal para el valor máximo en modo univariado.
#' @return Objeto `ggplot`.
#' @export
plot_model_benchmarks <- function(df_benchmarks,
                                  metric_col = "F1_Score",
                                  metric_cols = NULL,
                                  title = "Comparativa de Rendimiento de Modelos",
                                  high_color = "#27ae60") {
  stopifnot(is.data.frame(df_benchmarks))
  stopifnot("Modelo" %in% colnames(df_benchmarks))

  if (!is.null(metric_cols) && length(metric_cols) > 1) {
    valid_cols <- intersect(metric_cols, colnames(df_benchmarks))
    stopifnot(length(valid_cols) > 0)

    order_col <- if ("F1_Score" %in% colnames(df_benchmarks)) "F1_Score" else valid_cols[1]
    model_order <- df_benchmarks %>%
      dplyr::arrange(.data[[order_col]]) %>%
      dplyr::pull(Modelo) %>%
      unique()

    df_long <- df_benchmarks %>%
      dplyr::select(Modelo, dplyr::all_of(valid_cols)) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(valid_cols), names_to = "Metrica", values_to = "Valor") %>%
      dplyr::filter(!is.na(Valor)) %>%
      dplyr::mutate(
        Modelo = factor(Modelo, levels = model_order),
        Metrica = factor(Metrica, levels = valid_cols)
      )

    metric_palette <- c(
      "F1_Score" = "#2980b9", "ROC_AUC" = "#27ae60", "Exactitud_Balanceada" = "#8e44ad",
      "Exactitud" = "#16a085", "Sensibilidad" = "#d35400", "Especificidad" = "#2c3e50",
      "Precision" = "#f39c12", "MCC" = "#7f8c8d", "FDR" = "#c0392b"
    )

    return(
      ggplot2::ggplot(df_long, ggplot2::aes(x = Modelo, y = Valor, fill = Metrica)) +
        ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.72) +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(values = metric_palette) +
        ggplot2::scale_y_continuous(limits = c(0, 1.05)) +
        ggplot2::labs(title = title, subtitle = "Todos los modelos evaluados", x = "Modelo / Estrategia", y = "Puntaje (0 a 1)", fill = "Métrica") +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
          plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "#5d6d7e"),
          legend.position = "bottom"
        )
    )
  }

  stopifnot(metric_col %in% colnames(df_benchmarks))

  df_sorted <- df_benchmarks %>%
    dplyr::filter(!is.na(.data[[metric_col]])) %>%
    dplyr::arrange(.data[[metric_col]]) %>%
    dplyr::mutate(Modelo = factor(Modelo, levels = unique(Modelo)))

  ggplot2::ggplot(df_sorted, ggplot2::aes(x = .data[[metric_col]], y = Modelo, fill = .data[[metric_col]])) +
    ggplot2::geom_col(width = 0.7, color = "#2c3e50", linewidth = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", .data[[metric_col]])), hjust = -0.1, size = 3.5, fontface = "bold") +
    ggplot2::scale_fill_gradient(low = "#e8f8f5", high = high_color) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(title = title, x = metric_col, y = "Modelo / Estrategia") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5), legend.position = "none")
}

#' Grafica la distribución de veredictos del sistema de consenso multimodelo
#'
#' @param df_consensus Data frame de salida de `consolidar_consenso()`.
#' @param verdict_col Nombre de la columna con el veredicto (`"Veredicto_Consenso"`).
#' @param title Título del gráfico.
#' @param subtitle Subtítulo del gráfico.
#' @return Objeto `ggplot`.
#' @export
plot_consensus_distribution <- function(df_consensus,
                                        verdict_col = "Veredicto_Consenso",
                                        title = "Clasificación de Secuencias por Consenso Multimodelo",
                                        subtitle = "Consenso de 7 Modelos (Homólogo Válido vs. Ruido Descartado)") {
  stopifnot(is.data.frame(df_consensus))
  stopifnot(verdict_col %in% colnames(df_consensus))

  df_counts <- df_consensus %>%
    dplyr::count(.data[[verdict_col]], name = "Frecuencia")

  ggplot2::ggplot(df_counts, ggplot2::aes(x = .data[[verdict_col]], y = Frecuencia, fill = .data[[verdict_col]])) +
    ggplot2::geom_col(color = "black", alpha = 0.85, width = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = Frecuencia), vjust = -0.4, size = 4.5, fontface = "bold") +
    ggplot2::scale_fill_manual(values = c(
      "Homologo_Valido" = "#27ae60",
      "Ortologo_Valido" = "#27ae60",
      "Ruido_Descartado" = "#c0392b"
    )) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::labs(title = title, subtitle = subtitle, x = "Resultado del Consenso", y = "Número de Secuencias", fill = "Categoría") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "#5d6d7e"),
      legend.position = "none"
    )
}

#' Grafica el mapa de calor de variación en el número de copias (CNV)
#'
#' @param df_cnv_annotated Data frame anotado con `Clasificacion_Copia` (opcional).
#' @param cnv_matrix Data frame en formato ancho `Locus x Muestras` (opcional).
#' @param sample_col Columna de muestra en `df_cnv_annotated`.
#' @param locus_col Columna de locus u ortogrupo en `df_cnv_annotated`.
#' @param top_n_loci Número máximo de loci más abundantes a representar (por defecto `30`).
#' @param title Título del mapa de calor.
#' @param subtitle Subtítulo del mapa de calor.
#' @return Objeto `ggplot`.
#' @export
plot_cnv_heatmap <- function(df_cnv_annotated = NULL,
                             cnv_matrix = NULL,
                             sample_col = "Sample",
                             locus_col = "Locus_Ortogrupo",
                             top_n_loci = 30,
                             title = "Matriz de Variación en el Número de Copias (CNV)",
                             subtitle = "Conteo de copias funcionales estimadas tras el filtro de cobertura") {
  if (!is.null(df_cnv_annotated) && is.data.frame(df_cnv_annotated) &&
      all(c(sample_col, locus_col) %in% colnames(df_cnv_annotated))) {
    all_samples <- sort(unique(as.character(df_cnv_annotated[[sample_col]])))
    df_valid <- df_cnv_annotated
    if ("Clasificacion_Copia" %in% colnames(df_valid)) {
      df_valid <- df_valid %>%
        dplyr::filter(Clasificacion_Copia %in% c("Copia_Principal", "Homologo_Confirmado", "Paralogo_Confirmado"))
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
      dplyr::rename(Muestra = dplyr::all_of(sample_col), Locus = dplyr::all_of(locus_col)) %>%
      tidyr::complete(Muestra = all_samples, Locus = top_loci, fill = list(Numero_Copias = 0))
  } else {
    mat_src <- if (!is.null(cnv_matrix) && is.data.frame(cnv_matrix)) cnv_matrix else df_cnv_annotated
    if (!is.data.frame(mat_src)) {
      stop("[ERROR] Debe proporcionar `df_cnv_annotated` o `cnv_matrix` a plot_cnv_heatmap().", call. = FALSE)
    }
    first_col <- colnames(mat_src)[1]
    num_cols <- setdiff(colnames(mat_src)[sapply(mat_src, is.numeric)], c(first_col, "Total_Copias", "Max_Copias_Muestra"))
    df_long <- mat_src %>%
      dplyr::select(dplyr::all_of(c(first_col, num_cols))) %>%
      dplyr::mutate(Total_Copias = rowSums(dplyr::across(dplyr::all_of(num_cols)), na.rm = TRUE)) %>%
      dplyr::arrange(dplyr::desc(Total_Copias)) %>%
      dplyr::slice_head(n = top_n_loci) %>%
      dplyr::select(-Total_Copias) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(num_cols), names_to = "Muestra", values_to = "Numero_Copias") %>%
      dplyr::mutate(Numero_Copias = tidyr::replace_na(as.numeric(Numero_Copias), 0)) %>%
      dplyr::rename(Locus = dplyr::all_of(first_col))
  }

  hide_x_labels <- length(unique(df_long$Muestra)) > 25

  ggplot2::ggplot(df_long, ggplot2::aes(x = Muestra, y = Locus, fill = Numero_Copias)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.15) +
    ggplot2::scale_fill_viridis_c(option = "mako", direction = -1, name = "Copias") +
    ggplot2::labs(title = title, subtitle = subtitle, x = "Muestras Analizadas", y = "Locus / Ortogrupo") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "#5d6d7e"),
      axis.text.x = if (hide_x_labels) ggplot2::element_blank() else ggplot2::element_text(angle = 45, hjust = 1),
      axis.ticks.x = if (hide_x_labels) ggplot2::element_blank() else ggplot2::element_line(),
      panel.grid = ggplot2::element_blank()
    )
}

#' Grafica curvas ROC comparativas a partir de una lista de objetos `pROC::roc`
#'
#' @param roc_list Lista nombrada de objetos `roc`.
#' @param title Título del gráfico.
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
