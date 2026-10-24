# ==============================================================================
# TargetCNV-ML: Rutinas de Visualizacion Cientifica con ggplot2 (R/plotting.R)
# ==============================================================================

#' Paleta de colores distintiva por modelo para Matrices de Confusion e Importancia de Variables
#'
#' Cada modelo evaluado (tanto ML como Heuristico y Cortes de Youden) recibe un
#' color propio consistente entre su matriz de confusion y su grafico de importancia.
#'
#' @param model_name Nombre del modelo o estrategia.
#' @return Lista con `high` (color principal) y `low` (color suave de gradiente).
#' @export
get_model_color_palette <- function(model_name) {
  palettes <- list(
    ElasticNet            = list(high = "#16a085", low = "#e8f8f5"), # Verde Azulado (Teal)
    NaiveBayes            = list(high = "#2980b9", low = "#ebf5fb"), # Azul Real
    SVM_Polinomial        = list(high = "#d35400", low = "#fdf2e9"), # Naranja Intenso
    Red_Neuronal_Bag      = list(high = "#884ea0", low = "#f5eef8"), # Purpura Bagging NNet
    Red_Neuronal          = list(high = "#c0392b", low = "#fdedec"), # Rojo Carmesi
    LDA                   = list(high = "#34495e", low = "#ebedef"), # Azul Pizarra
    Weighted_Heuristic    = list(high = "#8e44ad", low = "#f4ecf7"), # Morado Amatista
    Modelo_Heuristico     = list(high = "#8e44ad", low = "#f4ecf7"),
    Youden_pct_recovered  = list(high = "#27ae60", low = "#eafaf1"), # Verde Esmeralda
    Youden_pct_identity   = list(high = "#229954", low = "#e9f7ef"), # Verde Bosque
    Youden_score          = list(high = "#138d75", low = "#e8f8f5"), # Verde Mar
    Youden_wscore         = list(high = "#117a65", low = "#e8f6f3"), # Verde Pino
    Youden                = list(high = "#27ae60", low = "#eafaf1"),
    RandomForest          = list(high = "#1e8449", low = "#e9f7ef"),
    GBM                   = list(high = "#b9770e", low = "#fef9e7"),
    XGBoost               = list(high = "#6c3483", low = "#f5eef8"),
    SVM_Radial            = list(high = "#a04000", low = "#fbeee6"),
    SVM_Lineal            = list(high = "#ca6f1e", low = "#fef5e7"),
    KNN                   = list(high = "#2471a3", low = "#eaf2f8"),
    CART                  = list(high = "#b7950b", low = "#fef9e7"),
    MLP_Regularizado      = list(high = "#7b241c", low = "#f9ebea"),
    MLP                   = list(high = "#922b21", low = "#f9ebea"),
    DNN_Profunda          = list(high = "#922b21", low = "#f9ebea"),
    DNN_Regularizada      = list(high = "#7b241c", low = "#f9ebea"),
    QDA                   = list(high = "#5d6d7e", low = "#ebedef"),
    Bagging_Arboles       = list(high = "#148f77", low = "#e8f8f5"),
    `Stacked Autoencoder - DNN` = list(high = "#512e5f", low = "#f4ecf7"),
    Autoencoder           = list(high = "#512e5f", low = "#f4ecf7"),
    DNN_Autoencoder       = list(high = "#512e5f", low = "#f4ecf7"),
    Ensamble_Stacking     = list(high = "#5b2c6f", low = "#f4ecf7")
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
  return(list(high = "#2980b9", low = "#f1f2f6"))
}

#' Grafica una matriz de confusion 2x2 estilizada con TP en la esquina superior izquierda
#'
#' @param cm_table Tabla o matriz 2x2 de confusion.
#' @param title Titulo principal del grafico.
#' @param subtitle Subtitulo con metricas resumidas (ej. Accuracy y FDR).
#' @param high_color Color hexadecimal para la celda de mayor frecuencia.
#' @param low_color Color hexadecimal para las celdas de baja frecuencia.
#' @return Objeto `ggplot`.
#' @export
plot_confusion_matrix <- function(cm_table,
                                  title = "Matriz de Confusion",
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
      Referencia = dplyr::case_when(
        grepl("TRUE|TP", Referencia, ignore.case = TRUE) ~ "TP",
        TRUE ~ "FP"
      ),
      Prediccion = dplyr::case_when(
        grepl("TRUE|TP", Prediccion, ignore.case = TRUE) ~ "TP",
        TRUE ~ "FP"
      ),
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

#' Calcula y grafica la Importancia de Variables de un modelo con su color distintivo
#'
#' @param model_obj Objeto `train` de `caret` o vector nombrado de importancias/pesos.
#' @param df_data Data frame opcional con los predictores y la columna `Clase` (para respaldo).
#' @param model_name Nombre del modelo para el titulo y seleccion de color.
#' @param high_color Color principal opcional (si es NULL, usa la paleta del modelo).
#' @return Objeto `ggplot`.
#' @export
plot_model_varimp <- function(model_obj,
                              df_data = NULL,
                              model_name = "Modelo",
                              high_color = NULL) {
  pal <- get_model_color_palette(model_name)
  if (is.null(high_color)) high_color <- pal$high

  imp_df <- NULL

  # Caso A: vector numerico nombrado (ej. pesos del modelo heuristico o corte univariado)
  if (is.numeric(model_obj) && !is.null(names(model_obj))) {
    imp_df <- data.frame(
      Variable = names(model_obj),
      Importancia = as.numeric(model_obj),
      stringsAsFactors = FALSE
    )
  } else {
    # Caso B: objeto caret::train
    imp_df <- tryCatch({
      vi <- caret::varImp(model_obj, scale = TRUE)$importance
      if ("Overall" %in% colnames(vi)) {
        data.frame(
          Variable = rownames(vi),
          Importancia = as.numeric(vi$Overall),
          stringsAsFactors = FALSE
        )
      } else if ("TP" %in% colnames(vi)) {
        data.frame(
          Variable = rownames(vi),
          Importancia = as.numeric(vi$TP),
          stringsAsFactors = FALSE
        )
      } else {
        data.frame(
          Variable = rownames(vi),
          Importancia = as.numeric(vi[, 1]),
          stringsAsFactors = FALSE
        )
      }
    }, error = function(e) {
      NULL
    })
  }

  # Respaldo basado en sensibilidad de probabilidad predicha (permutacion por variable)
  if (is.null(imp_df) || nrow(imp_df) == 0 || all(is.na(imp_df$Importancia))) {
    train_x <- if (is.list(model_obj) && !is.null(model_obj$trainingData)) {
      model_obj$trainingData[, setdiff(colnames(model_obj$trainingData), ".outcome"), drop = FALSE]
    } else if (!is.null(df_data)) {
      df_data[, setdiff(colnames(df_data), c("Clase", "Validacion_Espacial")), drop = FALSE]
    } else {
      NULL
    }

    if (!is.null(train_x) && ncol(train_x) > 0) {
      base_probs <- tryCatch(
        stats::predict(model_obj, newdata = train_x, type = "prob")[, "TP"],
        error = function(e) NULL
      )
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

  if (is.null(imp_df) || nrow(imp_df) == 0) {
    return(NULL)
  }

  imp_df <- imp_df %>%
    dplyr::mutate(
      Importancia = ifelse(is.na(Importancia), 0, Importancia),
      Importancia = pmax(0, Importancia)
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
      subtitle = "Peso relativo en la decision del modelo",
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

#' Dibuja la arquitectura de Red Neuronal (soporta nnet, avNNet, MLP, MLP_Regularizado y Autoencoder)
#'
#' Para redes de 1 capa oculta (`nnet` y `avNNet`), utiliza `NeuralNetTools::plotnet`
#' sobre el modelo final (o la primera sub-red representativa del ensamble en `avNNet`).
#' Para redes multicapa profundas (`MLP`, `MLP_Regularizado`, `Autoencoder`), dibuja un
#' diagrama de topologia por capas (Entrada -> Capas Ocultas -> Salida) en `ggplot2`.
#'
#' @param model_obj Objeto `train` de `caret`.
#' @param title Titulo del grafico de arquitectura.
#' @export
plot_neural_net_topology <- function(model_obj,
                                     title = "Arquitectura de la Red Neuronal") {
  if (is.null(model_obj) || is.null(model_obj$finalModel)) {
    return(invisible(NULL))
  }

  fm <- model_obj$finalModel
  bt <- model_obj$bestTune

  # Caso 1: Red_Neuronal_Bag (avNNet) -> extraer la primera sub-red nnet del ensamble
  if (inherits(fm, "avNNet") && !is.null(fm$model) && length(fm$model) >= 1) {
    if (requireNamespace("NeuralNetTools", quietly = TRUE)) {
      tryCatch({
        NeuralNetTools::plotnet(
          fm$model[[1]],
          alpha_val = 0.75,
          circle_cex = 4.5,
          cex_val = 0.85,
          max_sp = TRUE,
          pos_col = "#16a085",
          neg_col = "#884ea0",
          circle_col = "#f5eef8",
          bord_col = "#884ea0"
        )
        graphics::title(main = title, font.main = 2, cex.main = 1.0)
        return(invisible(NULL))
      }, error = function(e) {})
    }
  }

  # Caso 2: Red_Neuronal (nnet estandar)
  if (inherits(fm, "nnet")) {
    if (requireNamespace("NeuralNetTools", quietly = TRUE)) {
      tryCatch({
        NeuralNetTools::plotnet(
          fm,
          alpha_val = 0.75,
          circle_cex = 4.5,
          cex_val = 0.85,
          max_sp = TRUE,
          pos_col = "#16a085",
          neg_col = "#c0392b",
          circle_col = "#fdedec",
          bord_col = "#c0392b"
        )
        graphics::title(main = title, font.main = 2, cex.main = 1.0)
        return(invisible(NULL))
      }, error = function(e) {})
    }
  }

  # Caso 3: Redes Multicapa (MLP, MLP_Regularizado, Autoencoder) -> Diagrama de topologia multicapa en ggplot2
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
    # Limitar nodos dibujados visualmente a max 12 por capa para legibilidad, indicando el total en la etiqueta
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
        x    = n1$x[grid_e$i],
        y    = n1$y[grid_e$i],
        xend = n2$x[grid_e$j],
        yend = n2$y[grid_e$j]
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
      ggplot2::scale_fill_manual(
        values = c("Entrada" = "#d6eaf8", "Oculta" = "#fdedec", "Salida" = "#d5f5e3")
      ) +
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
                                  title = "Separacion entre Verdaderos Positivos (TP) y Falsos Positivos (FP)",
                                  subtitle = "Distribucion de densidad y puntos de corte de Youden (linea punteada)",
                                  legend_position = "bottom") {
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
      x = "Valor de la Metrica",
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

  return(p)
}

#' Grafica el ranking comparativo de modelos para una o multiples metricas
#'
#' @param df_benchmarks Data frame de benchmark con la columna `Modelo`.
#' @param metric_col Metrica unica a graficar (ej. "F1_Score").
#' @param metric_cols Vector opcional con varias metricas para vista general multimetrica.
#' @param title Titulo del grafico.
#' @param high_color Color hexadecimal para el valor maximo en modo de una sola metrica.
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
      tidyr::pivot_longer(
        cols = dplyr::all_of(valid_cols),
        names_to = "Metrica",
        values_to = "Valor"
      ) %>%
      dplyr::filter(!is.na(Valor)) %>%
      dplyr::mutate(
        Modelo = factor(Modelo, levels = model_order),
        Metrica = factor(Metrica, levels = valid_cols)
      )

    metric_palette <- c(
      "F1_Score"             = "#2980b9",
      "ROC_AUC"              = "#27ae60",
      "Exactitud_Balanceada" = "#8e44ad",
      "Exactitud"            = "#16a085",
      "Sensibilidad"         = "#d35400",
      "Especificidad"        = "#2c3e50",
      "Precision"            = "#f39c12",
      "MCC"                  = "#7f8c8d",
      "FDR"                  = "#c0392b"
    )

    return(
      ggplot2::ggplot(df_long, ggplot2::aes(x = Modelo, y = Valor, fill = Metrica)) +
        ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.72) +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(values = metric_palette) +
        ggplot2::scale_y_continuous(limits = c(0, 1.05)) +
        ggplot2::labs(
          title = title,
          subtitle = "Todos los modelos evaluados",
          x = "Modelo / Estrategia",
          y = "Puntaje (0 a 1)",
          fill = "Metrica"
        ) +
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
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.3f", .data[[metric_col]])),
      hjust = -0.1,
      size = 3.5,
      fontface = "bold"
    ) +
    ggplot2::scale_fill_gradient(low = "#e8f8f5", high = high_color) +
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
      x = "Resultado del Consenso",
      y = "Numero de Secuencias",
      fill = "Categoria"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "#5d6d7e"),
      legend.position = "none"
    )
}

#' Grafica el Mapa de Calor (Heatmap) de Variacion en el Numero de Copias (CNV)
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
                             subtitle = "Conteo de copias funcionales luego del filtro de fragmentacion") {
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
  } else if (!is.null(cnv_matrix) && is.data.frame(cnv_matrix)) {
    first_col <- colnames(cnv_matrix)[1]
    num_cols <- setdiff(colnames(cnv_matrix)[sapply(cnv_matrix, is.numeric)], c(first_col, "Total_Copias", "Max_Copias_Muestra"))
    df_long <- cnv_matrix %>%
      dplyr::select(dplyr::all_of(c(first_col, num_cols))) %>%
      dplyr::mutate(Total_Copias = rowSums(dplyr::across(dplyr::all_of(num_cols)), na.rm = TRUE)) %>%
      dplyr::arrange(dplyr::desc(Total_Copias)) %>%
      dplyr::slice_head(n = top_n_loci) %>%
      dplyr::select(-Total_Copias) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(num_cols),
        names_to = "Muestra",
        values_to = "Numero_Copias"
      ) %>%
      dplyr::mutate(Numero_Copias = tidyr::replace_na(as.numeric(Numero_Copias), 0)) %>%
      dplyr::rename(Locus = dplyr::all_of(first_col))
  } else if (is.data.frame(df_cnv_annotated)) {
    first_col <- colnames(df_cnv_annotated)[1]
    num_cols <- setdiff(colnames(df_cnv_annotated)[sapply(df_cnv_annotated, is.numeric)], c(first_col, "Total_Copias", "Max_Copias_Muestra"))
    df_long <- df_cnv_annotated %>%
      dplyr::select(dplyr::all_of(c(first_col, num_cols))) %>%
      dplyr::mutate(Total_Copias = rowSums(dplyr::across(dplyr::all_of(num_cols)), na.rm = TRUE)) %>%
      dplyr::arrange(dplyr::desc(Total_Copias)) %>%
      dplyr::slice_head(n = top_n_loci) %>%
      dplyr::select(-Total_Copias) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(num_cols),
        names_to = "Muestra",
        values_to = "Numero_Copias"
      ) %>%
      dplyr::mutate(Numero_Copias = tidyr::replace_na(as.numeric(Numero_Copias), 0)) %>%
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
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "#5d6d7e"),
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
