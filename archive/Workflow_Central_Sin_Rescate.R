## PAQUETES Y RUTAS -----------------------------------------------------
pacman::p_load(
    Biostrings, tidyverse, nnet, GenomicRanges, rtracklayer,
    colorspace, ggsci, pROC, caret, doParallel, gbm, rpart, glmnet,
    rpart.plot, corrplot, conflicted, readxl, cowplot, xgboost,
    NeuralNetTools, klaR, RSNNS, caretEnsemble, DALEX, shapviz,
    deepnet
)

conflicts_prefer(dplyr::first, dplyr::filter, dplyr::select, dplyr::lag,
                 base::as.factor, pROC::coords, colorspace::coords,
                 caret::confusionMatrix, caret::train)

ruta_base <- "~/Documentos/Formacion/Maestria/Tesis/Archivos_Referencia/"
setwd(ruta_base)

## DICCIONARIO TEÓRICO DE REFERENCIA --------------------------------------
archivo <- "Genes_Capture_May_2026.xlsx"

cnv_benoit <- read_excel(archivo, sheet = "CNV_Benoit") %>%
    select(OG_ID, Probe_ID, `S._torvum`) %>%
    filter(!is.na(`S._torvum`) & trimws(`S._torvum`) != "")
probes_validos <- cnv_benoit$Probe_ID

contigs_probe <- read_excel(archivo, sheet = "Contigs_probe_design") %>%
    select(Probe_ID, Contig_name, HOG_ID_Solanum, OG_ID_Pangenomes,
           Correcte_Contig_Name, Short_name, Pathway, Refseq_RNA,
           DAVID_Annotation) %>%
    filter(Probe_ID %in% probes_validos)

sequences_torva <- read_excel(archivo, sheet = "Sequences_Torva") %>%
    filter(Probe_ID %in% probes_validos)

datos_unificados <- cnv_benoit %>% select(-OG_ID) %>%
    left_join(contigs_probe, by = "Probe_ID") %>%
    left_join(sequences_torva %>% select(-Orthogroup_ID), by = "Probe_ID")

pathways_filtro <- c(
    "Phenylpropanoid", "SGAs", "Anthocyanin", "Flavonoid", "Terpenoids",
    "ASs", "UGT", "Metabolites", "Trichome structure",
    "Glandular Trichome, Leaf development",
    "Trichome structure: Trichome differentiation",
    "Trichome structure: Glandular",
    "Trichome structure: Glandular type VI", "Trichome_structure"
)

datos_finales <- datos_unificados %>%
    filter(trimws(Pathway) %in% pathways_filtro) %>%
    filter(Assembly == "S. torvum Benoit et al") %>%
    select(-S._torvum) %>%
    left_join(sequences_torva %>% select(Seq, Orthogroup_ID), by = "Seq") %>%
    distinct() %>%
    left_join(cnv_benoit %>% mutate(Orthogroup_ID = OG_ID) %>%
                  filter(Probe_ID %in% datos_unificados$Probe_ID) %>%
                  select(Orthogroup_ID, S._torvum), by = "Orthogroup_ID") %>% 
    select(c(1,15:16,2:10,11:14))

dict_benoit <- datos_unificados %>%
    filter(trimws(Pathway) %in% pathways_filtro, 
           Assembly == "S. torvum Benoit et al") %>%
    select(-S._torvum) %>%
    left_join(sequences_torva %>% select(Seq, Orthogroup_ID), by = "Seq") %>%
    filter(!is.na(Contig_ID)) %>%
    select(Locus_Benoit = Contig_ID, OG_Benoit = Orthogroup_ID) %>%
    mutate(Locus_Benoit = str_remove(Locus_Benoit, "\\.[0-9]+$")) %>%
    distinct()

## CARGA DE GFFs (REFERENCIA Y CAPTUS) ------------------------------------
path_gff_ref <- paste0(ruta_base, 
                       "Solanum_torvum_genome/Solanum_torvum_genome/",
                       "Storvum_genes.gff3")
gff_ref_completo <- import.gff3(path_gff_ref)
genes_ref <- gff_ref_completo[gff_ref_completo$type == "gene"]
genes_ref$Locus_Ref_Limpio <- str_remove_all(genes_ref$ID, "^gene:|\\.[0-9]+$")

path_gff_captus <- paste0(ruta_base, "Captus_Local/03_extract_metabolismo/",
                          "Storvum__captus-ext/06_assembly_annotated/",
                          "Storvum_hit_contigs.gff")
gff_captus_completo <- import.gff3(path_gff_captus)

# Filtro estricto por marcadores funcionales (Scipio)
gff_captus <- gff_captus_completo[gff_captus_completo$type == "protein_match:NUC"]

## INTERSECCIÓN ESPACIAL Y RELACIÓN BASE ----------------------------------
# type = "any" con hebra obligatoria (ignore.strand = FALSE) para máxima precisión
hits_globales <- findOverlaps(query = gff_captus, subject = genes_ref, 
                              type = "any", ignore.strand = FALSE)

df_cruces <- data.frame(
    Captus_Row = queryHits(hits_globales),
    Ref_Row    = subjectHits(hits_globales)
) %>%
    mutate(
        ID_Captus = gff_captus$Name[Captus_Row],
        OG_Captus = str_extract(ID_Captus, "OG[0-9]+"),
        Locus_Ref_Limpio = genes_ref$Locus_Ref_Limpio[Ref_Row]
    )

## CLASIFICACIÓN ESPACIAL PRELIMINAR ---------------------------------------
df_clasificacion_hits <- df_cruces %>%
    left_join(dict_benoit, by = c("Locus_Ref_Limpio" = "Locus_Benoit")) %>%
    mutate(
        Validacion_Espacial = case_when(
            !is.na(OG_Benoit) & OG_Captus == OG_Benoit ~ "TRUE_POSITIVE",
            TRUE ~ "FALSE_POSITIVE"
        ),
        Doble_Control_Biologico = case_when(
            Validacion_Espacial == "TRUE_POSITIVE" ~ "1_Validado_Benoit_Exacto",
            !is.na(OG_Benoit) & OG_Captus != OG_Benoit ~ "2_Hit_Cruzado_OG_Vecino",
            is.na(OG_Benoit) ~ "3_Expansion_Nueva_Genoma"
        ),
        prioridad = if_else(Validacion_Espacial == "TRUE_POSITIVE", 0L, 1L)
    ) %>%
    arrange(prioridad) %>%
    distinct(ID_Captus, .keep_all = TRUE) %>%
    select(ID_Captus, Validacion_Espacial, Doble_Control_Biologico,
           Locus_Ref_Limpio)

## INTEGRACIÓN DE ESTADÍSTICAS --------------------------------------------
path_stats_captus <- paste0(ruta_base, 
                            "Captus_Local/03_extract_metabolismo/",
                            "captus-extract_stats.tsv")

stats_captus_crudo <- read.delim(path_stats_captus, header = TRUE,
                                 stringsAsFactors = FALSE) %>%
    filter(sample_name == "Storvum") %>%
    group_by(sample_name, marker_type, ref_name) %>%
    mutate(
        OG_Limpio = str_extract(ref_name, "OG[0-9]+"),
        # Reconstrucción fidedigna del ID usando la columna 'hit' nativa
        ID_GFF_Compatible = case_when(
            n() == 1 ~ OG_Limpio,
            n() > 1  ~ sprintf("%s__%02d", OG_Limpio, hit)
        )
    ) %>%
    ungroup()

stats_captus_integrado <- stats_captus_crudo %>%
    left_join(df_clasificacion_hits, 
              by = c("ID_GFF_Compatible" = "ID_Captus")) %>%
    filter(marker_type == "NUC") %>%
    mutate(
        Validacion_Espacial = replace_na(Validacion_Espacial, "FALSE_POSITIVE"),
        Doble_Control_Biologico = replace_na(Doble_Control_Biologico, 
                                             "4_Ruido_Intergenico_o_Intron"),
        ref_len_matched = as.numeric(ref_len_matched),
        pct_recovered = as.numeric(pct_recovered)
    )

## RESOLUCIÓN DE REDUNDANCIAS ---------------------------------------------

stats_captus <- stats_captus_integrado %>%
    group_by(Locus_Ref_Limpio, Validacion_Espacial) %>%
    # Ordenamos de mejor a peor métrica
    arrange(desc(ref_len_matched), desc(pct_recovered), .by_group = TRUE) %>%
    mutate(
        # Solo aplicamos la lógica si es TP y tiene un Locus Real
        Es_Redundante = if_else(
            Validacion_Espacial == "TRUE_POSITIVE" & 
                !is.na(Locus_Ref_Limpio) & row_number() > 1, 
            TRUE, 
            FALSE
        ),
        # Degradación oficial a Falso Positivo si es redundante
        Validacion_Espacial = if_else(Es_Redundante, "FALSE_POSITIVE", 
                                      Validacion_Espacial),
        Doble_Control_Biologico = if_else(Es_Redundante, 
                                          "5_TP_Redundante_Degradado", 
                                          Doble_Control_Biologico)
    ) %>%
    ungroup() %>%
    select(-Es_Redundante)

cat("Total de secuencias procesadas (NUC) :", nrow(stats_captus), "\n")
cat("Verdaderos Positivos Únicos (TP)     :", 
    sum(stats_captus$Validacion_Espacial == "TRUE_POSITIVE"), "\n")
cat("Falsos Positivos & Redundantes (FP)  :", 
    sum(stats_captus$Validacion_Espacial == "FALSE_POSITIVE"), "\n")
cat("\nDesglose de Origen Biológico:\n")
table(stats_captus$Doble_Control_Biologico)

# --- SWITCH: incluir la variable 'hit' como predictor en los modelos de ML --
# TRUE  -> los modelos de ML usan pct_identity, pct_recovered, wscore, score, hit
# FALSE -> los modelos de ML usan solo pct_identity, pct_recovered, wscore, score
incluir_variable_hit <- FALSE

stats_captus <- stats_captus %>%
    group_by(locus) %>%
    # Ordenamos por score para asegurar que el mejor esté arriba
    arrange(desc(score), .by_group = TRUE) %>%
    mutate(
        # Calculamos la diferencia contra el mejor hit del mismo ortogrupo
        delta_wscore = max(wscore, na.rm = TRUE) - wscore,
        delta_pct_identity = max(pct_identity, na.rm = TRUE) - pct_identity
    ) %>%
    ungroup()

# --- Nombres de columnas (explicitos, no por indice -> mas legible y a prueba
#     de que el orden de columnas del tsv cambie en una corrida futura) ------
variable_respuesta <- "Validacion_Espacial"
metricas_numericas <- c("pct_recovered", "pct_identity", "score", "wscore")
col_hit <- "hit"

set.seed(1234)
num_cores <- max(1, parallel::detectCores() - 1)

## UNIFICACION DE METADATA DE GENES (PATHWAY / ANOTACION) -----------------
stats_captus[[variable_respuesta]] <- factor(
    stats_captus[[variable_respuesta]],
    levels = c("TRUE_POSITIVE", "FALSE_POSITIVE")
)

## GRAFICOS DESCRIPTIVOS ----------------------------------------------------
graficar_descriptivo <- function(datos, tipo = c("boxplot", "jitter")) {
    tipo <- match.arg(tipo)
    p <- datos %>%
        select(all_of(c(variable_respuesta, metricas_numericas))) %>%
        pivot_longer(cols = all_of(metricas_numericas),
                     names_to = "Predictor", values_to = "Valor") %>%
        ggplot(aes(x = .data[[variable_respuesta]], y = Valor))
    
    p <- if (tipo == "boxplot") {
        p + geom_boxplot(aes(fill = .data[[variable_respuesta]]),
                         alpha = 0.7, outlier.colour = "royalblue3", 
                         outlier.size = 1.5)
    } else {
        p + geom_jitter(aes(colour = .data[[variable_respuesta]]), alpha = 0.7)
    }
    
    p +
        facet_wrap(~Predictor, scales = "free_y", ncol = 2) +
        scale_fill_npg() + scale_colour_npg() +
        labs(title = "Caracterizacion de Secuencias CAPTUS: TP vs FP",
             subtitle = paste("Contraste de descriptores continuos para",
                              unique(datos$sample_name)),
             x = "", y = "Valor Numerico") +
        theme_bw(base_size = 12) +
        theme(legend.position = "none",
              strip.text = element_text(face = "bold", size = 10),
              axis.text.x = element_text(angle = 15, hjust = 1))
}
print(graficar_descriptivo(stats_captus, "boxplot"))
print(graficar_descriptivo(stats_captus, "jitter"))

stats_captus %>%
    select(all_of(c(variable_respuesta, col_hit, metricas_numericas))) %>%
    pivot_longer(cols = all_of(metricas_numericas),
                 names_to = "Predictor", 
                 values_to = "Valor") %>%
    ggplot(aes(x = .data[[col_hit]], y = Valor)) +
    geom_jitter(aes(colour = .data[[variable_respuesta]]), alpha = 0.4,
                position = position_dodge(width = 0.75)) +
    facet_wrap(~Predictor, scales = "free_y", ncol = 2) +
    scale_colour_npg() +
    labs(title = "Dinamica de la Calidad del Ensamblaje segun Frecuencia de Mapeo",
         subtitle = paste(
             "Contraste TP vs FP estratificado por el numero de hits en",
             unique(stats_captus$sample_name)),
         x = "Numero de Hits Registrados por CAPTUS (Categoria)",
         y = "Valor Numerico", colour = "Validacion Espacial") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom", 
          strip.text = element_text(face = "bold", size = 10))

## GLM (LOGISTICA BINOMIAL / MULTINOMIAL) ----------------------------------
formula_modelo <- as.formula(paste(variable_respuesta, "~", 
                                   paste(metricas_numericas, collapse = " + ")))
num_categorias <- length(levels(stats_captus[[variable_respuesta]]))
if (num_categorias == 2) {
    cat("\nAplicando Regresion Logistica Binomial...\n")
    modelo_glm_simple <- glm(formula_modelo,
                             data = stats_captus, 
                             family = binomial)
    print(summary(modelo_glm_simple))
} else if (num_categorias > 2) {
    cat("\nAplicando Regresion Logistica Multinomial...\n")
    modelo_glm_simple <- multinom(formula_modelo, 
                                  data = stats_captus, trace = FALSE)
    print(summary(modelo_glm_simple))
} else {
    cat("\n¡ADVERTENCIA! Solo hay una categoria de validacion.\n")
}

## CORRELACION Y DISTRIBUCION POR HIT --------------------------------------
corrplot(cor(stats_captus[, metricas_numericas], 
             method = "spearman",
             use = "pairwise.complete.obs"))

stats_captus %>%
    mutate(hit_factor = factor(.data[[col_hit]])) %>%
    ggplot(aes(hit_factor, fill = hit_factor)) +
    geom_bar() + coord_flip() +
    facet_wrap(as.formula(paste("~", variable_respuesta))) +
    theme_bw() + theme(legend.position = "none")

## INDICE DE YOUDEN POR METRICA INDIVIDUAL ---------------------------------
df_corte <- stats_captus %>%
    select(all_of(c(variable_respuesta, metricas_numericas))) %>%
    drop_na()

df_corte %>%
    pivot_longer(cols = all_of(metricas_numericas), 
                 names_to = "Metrica", values_to = "Valor") %>%
    ggplot(aes(x = Valor, fill = .data[[variable_respuesta]])) +
    geom_density(alpha = 0.6) +
    facet_wrap(~Metrica, scales = "free") +
    scale_fill_manual(values = c("#27ae60", "#e74c3c")) +
    labs(title = "Zonas de Solapamiento y Separacion (TP vs FP)",
         x = "Valor de la Metrica", y = "Densidad", fill = "") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom", strip.text = element_text(face = "bold"))

benchmark_youden_list <- list()
matrices_youden_plots <- list()
df_votos_binarios <- data.frame(matrix(ncol = 0, nrow = nrow(df_corte)))

graficar_matriz_confusion <- function(cm_table, titulo, subtitulo,
                                      color_alto = "#d35400") {
    cm_df <- as.data.frame(cm_table)
    colnames(cm_df) <- c("Prediccion", "Referencia", "Frecuencia")
    cm_df <- cm_df %>%
        mutate(
            Referencia = case_when(
                grepl("TRUE|TP", Referencia, ignore.case = TRUE) ~ "TP",
                TRUE ~ "FP"
            ),
            Prediccion = case_when(
                grepl("TRUE|TP", Prediccion, ignore.case = TRUE) ~ "TP",
                TRUE ~ "FP"
            )
        ) %>% 
        mutate(
            # Configuración de niveles para que TP quede en la esquina superior izquierda
            Referencia  = factor(Referencia, levels = c("TP", "FP")),
            Prediccion  = factor(Prediccion, levels = c("FP", "TP")) 
        )
    
    ggplot(cm_df, aes(x = Referencia, y = Prediccion)) +
        geom_tile(aes(fill = Frecuencia), color = "white", linewidth = 1) +
        geom_text(aes(label = round(Frecuencia, 0)), size = 8, fontface = "bold",
                  color = ifelse(cm_df$Frecuencia > max(cm_df$Frecuencia) / 2, 
                                 "white", "black")) +
        scale_fill_gradient(low = "#f1f2f6", high = color_alto, name = "Frecuencia") +
        labs(title = titulo, subtitle = subtitulo, 
             x = "Clase Real (Ground Truth)", 
             y = "Prediccion") +
        theme_minimal(base_size = 12) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
              plot.subtitle = element_text(hjust = 0.5,
                                           color = "#7f8c8d", size = 10),
              axis.title = element_text(face = "bold", size = 10),
              legend.position = "none", panel.grid = element_blank())
}

for (var in metricas_numericas) {
    roc_obj <- roc(response = df_corte[[variable_respuesta]], 
                   predictor = df_corte[[var]], quiet = TRUE)
    best_coords <- pROC::coords(roc_obj, "best",
                                ret = c("threshold", "specificity", 
                                        "sensitivity"), 
                                best.method = "youden")
    umbral <- best_coords$threshold[1]
    
    cat(sprintf("Metrica: %s -> corte sugerido >= %.2f | Sens: %.1f%% | Espec: %.1f%%\n",
                toupper(var), umbral, best_coords$sensitivity[1] * 100, 
                best_coords$specificity[1] * 100))
    
    df_votos_binarios[[var]] <- ifelse(df_corte[[var]] >= umbral, 1, 0)
    predicciones <- factor(ifelse(df_corte[[var]] >= umbral, 
                                  "TRUE_POSITIVE", "FALSE_POSITIVE"),
                           levels = levels(df_corte[[variable_respuesta]]))
    cm <- confusionMatrix(data = predicciones, 
                          reference = df_corte[[variable_respuesta]],
                          mode = "everything")
    precision <- as.numeric(cm$byClass["Precision"])
    
    benchmark_youden_list[[var]] <- data.frame(
        Modelo = paste("Corte:", var),
        Exactitud = as.numeric(cm$overall["Accuracy"]),
        Sensibilidad = as.numeric(cm$byClass["Sensitivity"]),
        Especificidad = as.numeric(cm$byClass["Specificity"]),
        Precision = precision, FDR = 1 - precision,
        F1_Score = as.numeric(cm$byClass["F1"])
    )
    matrices_youden_plots[[var]] <- graficar_matriz_confusion(
        cm$table, paste("Predictor Univariado:", toupper(var)),
        sprintf("Corte: >= %.2f | Acc: %.1f%% | FDR: %.1f%%", umbral, 
                as.numeric(cm$overall["Accuracy"]) * 100, (1 - precision) * 100)
    )
}

## MODELO HEURISTICO COMBINADO (VOTO PONDERADO) ----------------------------
metricas_youden <- c("pct_recovered", "pct_identity", "score", "wscore")
grid_pesos <- expand.grid(w1 = 0:3, w2 = 0:3, w3 = 0:3, w4 = 0:3)
grid_pesos <- grid_pesos[rowSums(grid_pesos) > 0, ]

mejor_auc <- 0
mejores_pesos <- NULL
mejor_puntaje <- NULL
for (i in seq_len(nrow(grid_pesos))) {
    pesos_actuales <- as.numeric(grid_pesos[i, ])
    puntaje_actual <- as.matrix(df_votos_binarios[, metricas_youden]) %*% pesos_actuales
    roc_actual <- suppressMessages(roc(response = df_corte[[variable_respuesta]],
                                       predictor = as.vector(puntaje_actual), 
                                       quiet = TRUE))
    if (roc_actual$auc > mejor_auc) {
        mejor_auc <- roc_actual$auc
        mejores_pesos <- pesos_actuales
        mejor_puntaje <- as.vector(puntaje_actual)
    }
}
cat("\nPesos del modelo heuristico:\n")
print(setNames(mejores_pesos, metricas_youden))

roc_comb <- suppressMessages(roc(response = df_corte[[variable_respuesta]], 
                                 predictor = mejor_puntaje, quiet = TRUE))
umbral_comb <- pROC::coords(roc_comb, "best", ret = "threshold", 
                            best.method = "youden")$threshold[1]
predicciones_comb <- factor(ifelse(mejor_puntaje >= umbral_comb, 
                                   "TRUE_POSITIVE", "FALSE_POSITIVE"),
                            levels = levels(df_corte[[variable_respuesta]]))
cm_comb <- confusionMatrix(data = predicciones_comb, 
                           reference = df_corte[[variable_respuesta]], 
                           mode = "everything")
precision_comb <- as.numeric(cm_comb$byClass["Precision"])

benchmark_youden_list[["Modelo_Combinado"]] <- data.frame(
    Modelo = "Modelo Combinado (Heuristico)",
    Exactitud = as.numeric(cm_comb$overall["Accuracy"]),
    Sensibilidad = as.numeric(cm_comb$byClass["Sensitivity"]),
    Especificidad = as.numeric(cm_comb$byClass["Specificity"]),
    Precision = precision_comb, FDR = 1 - precision_comb,
    F1_Score = as.numeric(cm_comb$byClass["F1"])
)
matrices_youden_plots[["Modelo_Combinado"]] <- graficar_matriz_confusion(
    cm_comb$table, "Modelo Combinado (Heuristico)",
    sprintf("Corte Global: >= %.1f | Acc: %.1f%% | FDR: %.1f%%", umbral_comb, 
            as.numeric(cm_comb$overall["Accuracy"]) * 100, (1 - precision_comb) * 100),
    color_alto = "#8e44ad"
)

matriz_youden_grid <- plot_grid(plotlist = matrices_youden_plots, 
                                ncol = 2, align = "v")
print(plot_grid(ggdraw() + 
                    draw_label(
                        "Matrices de Confusion: Cortes Simples vs. Modelo Combinado", 
                        fontface = "bold", size = 16),
                matriz_youden_grid, ncol = 1, rel_heights = c(0.05, 1)))

benchmark_youden_df <- bind_rows(benchmark_youden_list)
print(
    benchmark_youden_df %>%
        pivot_longer(cols = -Modelo, 
                     names_to = "Metrica", values_to = "Valor") %>%
        ggplot(aes(x = reorder(Modelo, Valor), y = Valor, fill = Metrica)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        coord_flip() + scale_fill_npg() +
        labs(
            title = "Benchmark de Filtros Heuristicos (Indice de Youden)",
            subtitle = "Contraste entre variables aisladas y el Sistema de Voto Ponderado",
            x = "Filtro Utilizado", y = "Puntaje (0 a 1)") +
        theme_bw(base_size = 12) + 
        theme(legend.position = "bottom",
              plot.title = element_text(face = "bold"))
)

# MACHINE LEARNING---------------------------------------------------------

## PREPARACION DE df_ml (respeta el switch incluir_variable_hit) ----------
predictores_ml <- if (incluir_variable_hit) c(metricas_youden, col_hit) else metricas_youden
cat("\nPredictores usados en los modelos de ML:", 
    paste(predictores_ml, collapse = ", "), "\n")

df_ml <- stats_captus %>%
    mutate(
        Clase = factor(if_else(.data[[variable_respuesta]] == "TRUE_POSITIVE", 
                               "TP", "FP"), 
                       levels = c("TP", "FP")),
        across(all_of(predictores_ml), ~ as.numeric(as.character(.x)))
    ) %>%
    select(Clase, all_of(predictores_ml)) %>%
    drop_na() %>%
    as.data.frame()

## ESPECIFICACION PROGRAMATICA DE MODELOS ----------------------------------

gbm_grid <- expand.grid(interaction.depth = c(1, 3, 5), 
                        n.trees = c(50, 100, 150),
                        shrinkage = c(0.05, 0.1), n.minobsinnode = 10)
dnn_grid <- expand.grid(layer1 = c(10, 20), layer2 = c(5, 10), layer3 = c(0, 5))

# Sintoniza el número de neuronas por capa oculta y la tasa de dropout (Stacked Autoencoder)
dnn_auto_grid <- expand.grid(
    layer1 = c(16, 32),
    layer2 = c(8, 16),
    layer3 = c(4, 8),
    hidden_dropout = c(0, 0.1),
    visible_dropout = c(0)
)

especificaciones_modelos <- list(
    RandomForest     = list(method = "rf",                tuneLength = 3),
    GBM              = list(method = "gbm",               tuneGrid = gbm_grid),
    XGBoost          = list(method = "xgbTree",           tuneLength = 3),
    SVM_Radial       = list(method = "svmRadial",         tuneLength = 3),
    SVM_Lineal       = list(method = "svmLinear",         tuneLength = 3),
    SVM_Polinomial   = list(method = "svmPoly",           tuneLength = 3),
    Red_Neuronal     = list(method = "nnet",              tuneLength = 3),
    Red_Neuronal_Bag = list(method = "avNNet",            tuneLength = 3),
    ElasticNet       = list(method = "glmnet",            tuneLength = 3),
    KNN              = list(method = "knn",               tuneLength = 3),
    CART             = list(method = "rpart",             tuneLength = 3),
    NaiveBayes       = list(method = "nb",                tuneLength = 3),
    LDA              = list(method = "lda"),
    DNN_Profunda     = list(method = "mlpML",             tuneGrid = dnn_grid),
    DNN_Regularizada = list(method = "mlpWeightDecayML",  tuneLength = 3),
    QDA              = list(method = "qda"),
    Bagging_Arboles  = list(method = "treebag"),
    DNN_Autoencoder  = list(method = "dnn",               tuneGrid = dnn_auto_grid)
)
metodos_sin_escalado <- c("rf", "gbm", "xgbTree", "rpart", "treebag")

entrenar_modelos <- function(especificaciones, datos, control, semilla = 1234) {
    resultados <- list()
    x_train <- datos[, base::setdiff(names(datos), "Clase"), drop = FALSE]
    y_train <- datos$Clase
    for (nombre in names(especificaciones)) {
        spec <- especificaciones[[nombre]]
        args_train <- list(x = x_train, y = y_train, method = spec$method,
                           trControl = control, metric = "ROC")
        if (!is.null(spec$tuneGrid))   args_train$tuneGrid   <- spec$tuneGrid
        if (!is.null(spec$tuneLength)) args_train$tuneLength <- spec$tuneLength
        if (!spec$method %in% metodos_sin_escalado) args_train$preProcess <- c("center", "scale")
        if (spec$method == "gbm")                   args_train$verbose <- FALSE
        if (spec$method == "xgbTree")               args_train$nthread <- 1
        if (spec$method %in% c("nnet", "avNNet")) { 
            args_train$trace <- FALSE
            args_train$MaxNWts <- 5000 }
        
        set.seed(semilla)
        resultados[[nombre]] <- tryCatch({
            do.call(caret::train, args_train)
        }, error = function(e) {
            cat("\n[!] Alerta: El modelo", nombre,
                "ha fallado y fue omitido. Error:", e$message, "\n")
            return(NULL) # Retorna NULL para este modelo
        })
        if(!is.null(resultados[[nombre]])) {
            cat("Entrenado:", nombre, "\n")
        }
    }
    return(resultados[!sapply(resultados, is.null)])
}

fitControl <- trainControl(method = "cv", number = 5, classProbs = TRUE,
                           summaryFunction = twoClassSummary, 
                           allowParallel = TRUE)
fitControl_balanceado <- trainControl(method = "cv", number = 5,
                                      classProbs = TRUE,
                                      summaryFunction = twoClassSummary, 
                                      allowParallel = TRUE,
                                      savePredictions = "final", 
                                      sampling = "down")

## ENTRENAMIENTO: MODELOS DESBALANCEADOS ----------------------------------
cl <- makePSOCKcluster(num_cores); registerDoParallel(cl)
modelos_desbalanceados <- entrenar_modelos(especificaciones_modelos, df_ml,
                                           fitControl)
stopCluster(cl); registerDoSEQ()

resultados_torneo <- resamples(modelos_desbalanceados)
print(summary(resultados_torneo))
print(bwplot(resultados_torneo, metric = "ROC", 
             main = "Comparacion de Modelos: Area Bajo la Curva (ROC)"))

roc_means_unbal <- summary(resultados_torneo)$statistics$ROC[, "Mean"]
mejor_modelo_unbal <- names(which.max(roc_means_unbal))
print(
    data.frame(Modelo = names(roc_means_unbal), ROC_Medio = roc_means_unbal) %>%
        ggplot(aes(x = reorder(Modelo, ROC_Medio), y = ROC_Medio, fill = Modelo)) +
        geom_col(color = "black", alpha = 0.8, width = 0.7) +
        coord_flip() + scale_fill_npg() +
        labs(title = "Ranking de Modelos (Datos Desbalanceados)",
             subtitle = "Promedio del Area Bajo la Curva (ROC) en Validacion Cruzada",
             x = "Modelo Evaluado", y = "Valor ROC Promedio",
             caption = sprintf("Mejor modelo: %s con ROC de %.4f", 
                               mejor_modelo_unbal, max(roc_means_unbal))) +
        theme_bw(base_size = 12) +
        theme(legend.position = "none", plot.title = element_text(face = "bold"),
              plot.caption = element_text(face = "bold.italic", 
                                          size = 11, color = "#c0392b"))
)

## ARBOL DE REGLAS INTERPRETABLE ------------------------------------------
arbol_reglas <- rpart(Clase ~ ., data = df_ml, method = "class", 
                      control = rpart.control(cp = 0.01))
rpart.plot(arbol_reglas, type = 3, extra = 104, under = TRUE, 
           faclen = 0, cex = 0.85,
           main = "Arbol de Decision: Reglas de Filtrado",
           box.palette = c("#27ae60", "#e74c3c"))

## ENTRENAMIENTO: MODELOS BALANCEADOS (DOWN-SAMPLING) ---------------------
cl <- makePSOCKcluster(num_cores); registerDoParallel(cl)
modelos_balanceados <- entrenar_modelos(especificaciones_modelos, df_ml, 
                                        fitControl_balanceado)
stopCluster(cl); registerDoSEQ()

resultados_torneo_bal <- resamples(modelos_balanceados)

extraer_metricas <- function(modelo, nombre) {
    cm <- confusionMatrix(data = modelo$pred$pred, 
                          reference = modelo$pred$obs, mode = "everything")
    data.frame(
        Modelo = nombre,
        ROC = as.numeric(max(modelo$results$ROC, na.rm = TRUE)),
        Exactitud = as.numeric(cm$overall["Accuracy"]),
        Sensibilidad = as.numeric(cm$byClass["Sensitivity"]),
        Especificidad = as.numeric(cm$byClass["Specificity"]),
        Precision = as.numeric(cm$byClass["Precision"]),
        F1_Score = as.numeric(cm$byClass["F1"])
    )
}
benchmark_df <- purrr::imap_dfr(modelos_balanceados, extraer_metricas)

roc_promedios <- summary(resultados_torneo_bal)$statistics$ROC[, "Mean"]
mejor_modelo_bal <- names(which.max(roc_promedios))
print(
    benchmark_df %>%
        pivot_longer(cols = -Modelo, 
                     names_to = "Metrica", values_to = "Valor") %>%
        ggplot(aes(x = reorder(Modelo, Valor), y = Valor, fill = Metrica)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        coord_flip() + scale_fill_npg() +
        labs(title = "Benchmark de Modelos de Aprendizaje (Balanceados)",
             x = "Modelo Evaluado", y = "Puntaje (0 a 1)",
             caption = sprintf("Mejor modelo general (ROC): %s con %.4f", 
                               mejor_modelo_bal, max(roc_promedios))) +
        theme_bw(base_size = 12) +
        theme(legend.position = "bottom", 
              plot.title = element_text(face = "bold"),
              plot.caption = element_text(face = "bold.italic", 
                                          size = 11, color = "gray35"))
)

# Comparación directa Naive Bayes: desbalanceado vs balanceado
cm_bal <- confusionMatrix(
    data = modelos_balanceados[["NaiveBayes"]]$pred$pred,
    reference = modelos_balanceados[["NaiveBayes"]]$pred$obs,
    mode = "everything"
)

cat(sprintf("ROC  -> Desbalanceado: %.4f | Balanceado: %.4f\n",
            summary(resultados_torneo)$statistics$ROC["NaiveBayes", "Mean"],
            summary(resultados_torneo_bal)$statistics$ROC["NaiveBayes", "Mean"]))

cat(sprintf("Sens -> Desbalanceado: %.4f | Balanceado: %.4f\n",
            summary(resultados_torneo)$statistics$Sens["NaiveBayes", "Mean"],
            as.numeric(cm_bal$byClass["Sensitivity"])))

cat(sprintf("Spec -> Desbalanceado: %.4f | Balanceado: %.4f\n",
            summary(resultados_torneo)$statistics$Spec["NaiveBayes", "Mean"],
            as.numeric(cm_bal$byClass["Specificity"])))

cat(sprintf("Acc  -> Balanceado (Tabla 13): %.4f\n",
            as.numeric(cm_bal$overall["Accuracy"])))

## ATRICES DE CONFUSION: TOP 4 MODELOS -----------------------------------
graficar_matriz <- function(modelo, titulo = "Matriz de Confusion") {
    cm_table <- table(Predict = modelo$pred$pred, Obs = modelo$pred$obs)
    accuracy <- sum(diag(cm_table)) / sum(cm_table)
    graficar_matriz_confusion(cm_table, titulo, 
                              paste("Precision CV (Accuracy):", 
                                    round(accuracy * 100, 2), "%"),
                              color_alto = "#2980b9")
}

top4_nombres <- names(sort(roc_promedios, decreasing = TRUE)[1:4])
plots_top4 <- lapply(top4_nombres, 
                     function(nombre) graficar_matriz(modelos_balanceados[[nombre]], 
                                                      titulo = nombre))
print(plot_grid(
    ggdraw() + draw_label("Top 4 Modelos: Matrices de Confusion", 
                          fontface = "bold", size = 16),
    plot_grid(plotlist = plots_top4, ncol = 2, align = "v"),
    ncol = 1, rel_heights = c(0.1, 1)
))

## INTERPRETABILIDAD ------------------------------------------------------
### Elastic Net: ecuacion e importancia -----------------------------------
modelo_glmnet_bal <- modelos_balanceados[["ElasticNet"]]
mejor_lambda <- modelo_glmnet_bal$bestTune$lambda
print(coef(modelo_glmnet_bal$finalModel, s = mejor_lambda))
plot(varImp(modelo_glmnet_bal, scale = TRUE), 
     main = "Importancia de Variables (Elastic Net)")

### SHAP (interpretabilidad local, Elastic Net como modelo de referencia) ----
x_data <- df_ml %>% select(-Clase)
pred_fun <- function(modelo, datos) predict(modelo, 
                                            newdata = datos, 
                                            type = "prob")[, "TP"]

set.seed(1234)
shap_calc <- kernelshap::kernelshap(object = modelo_glmnet_bal, 
                                    X = x_data, bg_X = x_data, pred_fun = pred_fun)
shap_values <- shapviz(shap_calc)

print(
    sv_importance(shap_values, kind = "bee") +
        labs(title = "Impacto Direccional de las Variables (SHAP)",
             subtitle = "Como cada variable altera la probabilidad de ser Verdadero Positivo") +
        theme_bw(base_size = 12)
)

fila_ejemplo_shap <- 1
print(
    sv_waterfall(shap_values, row_id = fila_ejemplo_shap) +
        labs(title = paste("Anatomia de la Prediccion Individual (Fila", 
                           fila_ejemplo_shap, ")")) +
        theme_minimal()
)

### Topologia de red neuronal (nnet shallow, como referencia visual) ---------
modelo_nn_bal <- modelos_balanceados[["Red_Neuronal"]]
print(summary(modelo_nn_bal$finalModel))
plotnet(modelo_nn_bal$finalModel, alpha_val = 0.7, circle_cex = 5, cex_val = 0.8,
        max_sp = TRUE, pos_col = "darkgreen", neg_col = "darkred",
        main = "Topologia del Perceptron Multicapa (Red Neuronal)")
print(plot(varImp(modelo_nn_bal, scale = TRUE), 
           main = "Importancia del Input en la Red Neuronal"))

## ENSAMBLE DE MODELOS (META-MODELO STACKING) ----------------------------

# Seleccionar los 4 mejores modelos)
top4_nombres <- names(sort(roc_promedios, decreasing = TRUE)[1:2])
lista_top_modelos <- modelos_balanceados[top4_nombres]
class(lista_top_modelos) <- c("caretList", "list")

cat("Modelos base seleccionados para el ensamble:\n", 
    paste(top4_nombres, collapse = ", "), "\n\n")

# Comprobación de Correlación entre los Modelos Base
# Un ensamble es más fuerte cuando los modelos tienen baja correlación en sus errores
resultados_top4 <- resamples(lista_top_modelos)
(correlacion_modelos <- modelCor(resultados_top4))

# Entrenamiento del Meta-Modelo (Se usa un GLM para combinar las predicciones)
set.seed(1234)
ensamble_stack <- caretStack(
    lista_top_modelos,
    method = "glm",
    metric = "ROC",
    trControl = trainControl(
        method = "cv", 
        number = 5, 
        savePredictions = "final", 
        classProbs = TRUE, 
        summaryFunction = twoClassSummary
    )
)

cat("\nResumen del Meta-Modelo (GLM Stacking):\n")
print(ensamble_stack)

# Se generan las predicciones combinadas usando el meta-modelo
predicciones_cruda <- predict(ensamble_stack, 
                              newdata = df_ml %>% select(-Clase))

# Extraer la etiqueta final interpretando la salida
if (is.data.frame(predicciones_cruda) || is.matrix(predicciones_cruda)) {
    # Si devuelve 2 columnas (TP y FP), se extrae el nombre de la columna ganadora
    etiquetas_finales <- colnames(predicciones_cruda)[max.col(predicciones_cruda, 
                                                              ties.method = "first")]
} else if (is.numeric(predicciones_cruda)) {
    # Por si devuelve un vector numérico (probabilidad de TP)
    etiquetas_finales <- ifelse(predicciones_cruda > 0.5, "TP", "FP")
} else {
    # Por si devuelve texto directamente
    etiquetas_finales <- as.character(predicciones_cruda)
}
predicciones_stack <- factor(etiquetas_finales, levels = levels(df_ml$Clase))

# Se contruye la matriz de confusión del ensamble
cm_stack <- confusionMatrix(
    data = predicciones_stack, 
    reference = df_ml$Clase, 
    mode = "everything"
)

# Extraer métricas para comparar contra el mejor modelo individual
exactitud_stack <- as.numeric(cm_stack$overall["Accuracy"])
sensibilidad_stack <- as.numeric(cm_stack$byClass["Sensitivity"])
especificidad_stack <- as.numeric(cm_stack$byClass["Specificity"])

cat("\nDesempeño Final del Meta-Modelo:\n")
cat(sprintf("Exactitud (Accuracy): %.2f%%\n", exactitud_stack * 100))
cat(sprintf("Sensibilidad (TP Rate): %.2f%%\n", sensibilidad_stack * 100))
cat(sprintf("Especificidad (TN Rate): %.2f%%\n", especificidad_stack * 100))

# Graficamos la Matriz de Confusión usando tu función robusta
plot_ensamble <- graficar_matriz_confusion(
    cm_stack$table, 
    titulo = "Meta-Modelo Combinado (GLM Stacking)",
    subtitulo = sprintf("Precisión Global: %.2f%% | Combinación de Top 4", 
                        exactitud_stack * 100),
    color_alto = "#8e44ad" # Un tono morado para diferenciarlo de los modelos base
)

print(plot_ensamble)

## COMPARACION DE ESTRATEGIAS ----------------------------------------------

df_solo_univariados <- benchmark_youden_df %>%
    filter(Modelo != "Modelo Combinado (Heuristico)") %>%
    arrange(desc(F1_Score))

top2_uni_nombres <- df_solo_univariados$Modelo[1:2]
var_cruda_1 <- gsub("Corte: ", "", top2_uni_nombres[1])
var_cruda_2 <- gsub("Corte: ", "", top2_uni_nombres[2])

roc_uni_1 <- suppressMessages(roc(response = df_corte[[variable_respuesta]], 
                                  predictor = df_corte[[var_cruda_1]], 
                                  quiet = TRUE))
roc_uni_2 <- suppressMessages(roc(response = df_corte[[variable_respuesta]], 
                                  predictor = df_corte[[var_cruda_2]], 
                                  quiet = TRUE))

fila_uni_1 <- df_solo_univariados %>% 
    filter(Modelo == top2_uni_nombres[1]) %>% 
    mutate(ROC = as.numeric(roc_uni_1$auc))
fila_uni_2 <- df_solo_univariados %>% 
    filter(Modelo == top2_uni_nombres[2]) %>% 
    mutate(ROC = as.numeric(roc_uni_2$auc))

fila_heuristico <- benchmark_youden_df %>%
    filter(Modelo == "Modelo Combinado (Heuristico)") %>%
    mutate(ROC = as.numeric(roc_comb$auc))

filas_ml <- benchmark_df %>% 
    arrange(desc(ROC)) %>% 
    dplyr::slice(1:2) %>% 
    mutate(FDR = 1 - Precision)
top2_nombres <- filas_ml$Modelo

# Extraemos las probabilidades directas para calcular un ROC real del Ensamble
prob_base_list <- lapply(lista_top_modelos, function(modelo) {
    predict(modelo, newdata = df_ml %>% select(-Clase), 
            type = "prob")[["TP"]]
})
df_prob_base <- as.data.frame(prob_base_list)
colnames(df_prob_base) <- top4_nombres # Nombramos las columnas según el modelo

# Inyectar esas probabilidades directamente al GLM
vector_probs_stack <- predict(ensamble_stack$ens_model, 
                              newdata = df_prob_base, type = "prob")[["TP"]]

# Calcular el ROC real con los números puros
roc_ensamble <- suppressMessages(roc(response = df_ml$Clase, 
                                     predictor = vector_probs_stack, 
                                     quiet = TRUE))

fila_ensamble <- data.frame(
    Modelo = "Ensamble ML (Stacking)",
    Exactitud = exactitud_stack,
    Sensibilidad = sensibilidad_stack,
    Especificidad = especificidad_stack,
    Precision = as.numeric(cm_stack$byClass["Precision"]),
    FDR = 1 - as.numeric(cm_stack$byClass["Precision"]),
    F1_Score = as.numeric(cm_stack$byClass["F1"]),
    ROC = as.numeric(roc_ensamble$auc)
)

df_comparacion_final <- bind_rows(
    fila_uni_1, fila_uni_2, 
    fila_heuristico, 
    filas_ml, 
    fila_ensamble
)

matriz_batalla_grid <- plot_grid(
    matrices_youden_plots[[var_cruda_1]] + 
        labs(title = paste("Univariado 1:", toupper(var_cruda_1))),
    matrices_youden_plots[[var_cruda_2]] + 
        labs(title = paste("Univariado 2:", toupper(var_cruda_2))),
    matrices_youden_plots[["Modelo_Combinado"]],
    plot_ensamble + labs(title = "Ensamble ML (Stacking)"),
    graficar_matriz(modelos_balanceados[[top2_nombres[1]]], 
                    titulo = paste("ML 1:", top2_nombres[1])),
    graficar_matriz(modelos_balanceados[[top2_nombres[2]]], 
                    titulo = paste("ML 2:", top2_nombres[2])),
    ncol = 2, align = "v"
)

print(plot_grid(
    ggdraw() + draw_label("Cortes Empíricos vs. Machine Learning vs. Ensambles", 
                          fontface = "bold", size = 16),
    matriz_batalla_grid, ncol = 1, rel_heights = c(0.04, 1)
))

orden_modelos <- c(
    "Ensamble ML (Stacking)", 
    top2_nombres[1], top2_nombres[2], 
    "Modelo Combinado (Heuristico)", 
    top2_uni_nombres[1], top2_uni_nombres[2]
)

print(
    df_comparacion_final %>%
        pivot_longer(cols = -Modelo, names_to = "Metrica", values_to = "Valor") %>%
        mutate(Modelo = factor(Modelo, levels = rev(orden_modelos))) %>%
        ggplot(aes(x = Modelo, y = Valor, fill = Metrica)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        coord_flip() + scale_fill_npg() +
        labs(title = "Benchmark Definitivo de Estrategias de Clasificación",
             subtitle = "Comparativa: Univariados Top 2 | Heurístico | ML Top 2 | Meta-Modelo Ensamble",
             x = "Estrategia Evaluada", y = "Puntaje (0 a 1)") +
        theme_bw(base_size = 12) + 
        theme(legend.position = "bottom", 
              plot.title = element_text(face = "bold"))
)




#  EXPERIMENTO DE TRANSFERIBILIDAD: S. scitrullifolium

## 1. DICCIONARIO TEÓRICO DE REFERENCIA (CITRULLIFOLIUM) ----------------------
archivo <- "Genes_Capture_May_2026.xlsx"

cnv_citrulli <- read_excel(archivo, sheet = "CNV_Benoit") %>%
    select(OG_ID, Probe_ID, `S._citrullifolium`) %>%
    filter(!is.na(`S._citrullifolium`) & trimws(`S._citrullifolium`) != "")
probes_validos_cit <- cnv_citrulli$Probe_ID

contigs_probe_cit <- read_excel(archivo, sheet = "Contigs_probe_design") %>%
    select(Probe_ID, Contig_name, HOG_ID_Solanum, OG_ID_Pangenomes,
           Correcte_Contig_Name, Short_name, Pathway, Refseq_RNA,
           DAVID_Annotation) %>%
    filter(Probe_ID %in% probes_validos_cit)

sequences_cit <- read_excel(archivo, sheet = "S._citr") %>%
    filter(Probe_ID %in% probes_validos_cit)

datos_unificados_cit <- cnv_citrulli %>% select(-OG_ID) %>%
    left_join(contigs_probe_cit, by = "Probe_ID") %>%
    left_join(sequences_cit %>% select(-Orthogroup_ID), by = "Probe_ID") %>%
    mutate(Assembly = "S. citrullifolium Benoit et al")
    
dict_citrulli <- datos_unificados_cit %>%
    filter(trimws(Pathway) %in% pathways_filtro, 
           Assembly == "S. citrullifolium Benoit et al") %>%
    select(-`S._citrullifolium`) %>%
    left_join(sequences_cit %>% select(Seq, Orthogroup_ID), by = "Seq") %>%
    filter(!is.na(Contig_ID)) %>%
    select(Locus_Cit = Contig_ID, OG_Cit = Orthogroup_ID) %>%
    mutate(Locus_Cit = str_remove(Locus_Cit, "\\.[0-9]+$")) %>%
    distinct()

datos_finales_cit <- datos_unificados_cit %>%
    filter(trimws(Pathway) %in% pathways_filtro) %>%
    filter(Assembly == "S. citrullifolium Benoit et al") %>%
    select(-`S._citrullifolium`) %>%
    left_join(sequences_cit %>% select(Seq, Orthogroup_ID), by = "Seq") %>%
    distinct() %>%
    left_join(
        cnv_citrulli %>% 
            mutate(Orthogroup_ID = OG_ID) %>%
            filter(Probe_ID %in% datos_unificados_cit$Probe_ID) %>%
            select(Orthogroup_ID, `S._citrullifolium`), 
        by = "Orthogroup_ID"
    ) %>%
    select(c(1,15,2:10,11:14))

## 2. CARGA DE GFFs (REFERENCIA Y CAPTUS) ------------------------------------
path_gff_ref_cit <- paste0(ruta_base, 
                           "Solanum_citrullifolium_genome/Solanum_citrullifolium_genome/",
                           "Scitrullifolium_genes.gff3")
gff_ref_completo_cit <- import.gff3(path_gff_ref_cit)
genes_ref_cit <- gff_ref_completo_cit[gff_ref_completo_cit$type == "gene"]
genes_ref_cit$Locus_Ref_Limpio <- str_remove_all(genes_ref_cit$ID, "^gene:|\\.[0-9]+$")

path_gff_captus_cit <- paste0(ruta_base,
                              "Captus_Local/03_extract_citr_metabolism_complete/",
                              "Scitrullifolium__captus-ext/06_assembly_annotated/",
                              "Scitrullifolium_hit_contigs.gff")
gff_captus_completo_cit <- import.gff3(path_gff_captus_cit)

# Filtro estricto por marcadores funcionales (Scipio)
gff_captus_cit <- gff_captus_completo_cit[gff_captus_completo_cit$type == "protein_match:NUC"]

# Escudo de nombres: Limpiamos prefijos añadidos por CAPTUS para asegurar el cruce
seqlevels(gff_captus_cit) <- str_remove_all(seqlevels(gff_captus_cit), "^[Ss]citrullifolium_|^contig_")
seqlevels(genes_ref_cit)  <- str_remove_all(seqlevels(genes_ref_cit), "^[Ss]citrullifolium_|^contig_")

## 3. INTERSECCIÓN ESPACIAL Y RELACIÓN BASE ----------------------------------
hits_globales_cit <- findOverlaps(query = gff_captus_cit, subject = genes_ref_cit, 
                                  type = "any", ignore.strand = FALSE)

if (length(hits_globales_cit) == 0) {
    stop("Cero overlaps. Verifica los nombres de los scaffolds (seqnames) entre ambos GFFs.")
} else {
    cat(sprintf("\nÉxito: Se encontraron %d solapamientos físicos.\n", length(hits_globales_cit)))
}

df_cruces_cit <- data.frame(
    Captus_Row = queryHits(hits_globales_cit),
    Ref_Row    = subjectHits(hits_globales_cit)
) %>%
    mutate(
        ID_Captus = gff_captus_cit$Name[Captus_Row],
        OG_Captus = str_extract(ID_Captus, "OG[0-9]+"),
        Locus_Ref_Limpio = genes_ref_cit$Locus_Ref_Limpio[Ref_Row]
    )

## 4. CLASIFICACIÓN ESPACIAL PRELIMINAR ---------------------------------------
df_clasificacion_hits_cit <- df_cruces_cit %>%
    left_join(dict_citrulli, by = c("Locus_Ref_Limpio" = "Locus_Cit")) %>%
    mutate(
        Validacion_Espacial = case_when(
            !is.na(OG_Cit) & OG_Captus == OG_Cit ~ "TRUE_POSITIVE",
            TRUE ~ "FALSE_POSITIVE"
        ),
        Doble_Control_Biologico = case_when(
            Validacion_Espacial == "TRUE_POSITIVE" ~ "1_Validado_Benoit_Exacto",
            !is.na(OG_Cit) & OG_Captus != OG_Cit ~ "2_Hit_Cruzado_OG_Vecino",
            is.na(OG_Cit) ~ "3_Expansion_Nueva_Genoma"
        ),
        prioridad = if_else(Validacion_Espacial == "TRUE_POSITIVE", 0L, 1L),
        ID_Captus = if_else(str_detect(ID_Captus, "__[0-9]+$"),
                              ID_Captus,
                              paste0(ID_Captus, "__00"))
    ) %>%
    arrange(prioridad) %>%
    distinct(ID_Captus, .keep_all = TRUE) %>%
    select(ID_Captus, Validacion_Espacial, Doble_Control_Biologico, Locus_Ref_Limpio)

## 5. INTEGRACIÓN DE ESTADÍSTICAS Y REDUNDANCIAS ------------------------------
# Leemos el TSV general y filtramos por esta especie
stats_captus_crudo_cit <- read.delim(path_stats_captus, header = TRUE, 
                                     stringsAsFactors = FALSE) %>%
    filter(sample_name == "Scitrullifolium") %>%
    group_by(sample_name, marker_type, ref_name) %>%
    mutate(
        OG_Limpio = str_extract(ref_name, "OG[0-9]+"),
        ID_GFF_Compatible = sprintf("%s__%02d", OG_Limpio, hit)
    ) %>%
    ungroup()

stats_captus_cit <- stats_captus_crudo_cit %>%
    left_join(df_clasificacion_hits_cit, 
              by = c("ID_GFF_Compatible" = "ID_Captus")) %>% 
    filter(marker_type == "NUC") %>%
    mutate(
        Validacion_Espacial = replace_na(Validacion_Espacial, "FALSE_POSITIVE"),
        Doble_Control_Biologico = replace_na(Doble_Control_Biologico, 
                                             "4_Ruido_Intergenico_o_Intron"),
        ref_len_matched = as.numeric(ref_len_matched),
        pct_recovered = as.numeric(pct_recovered)
    ) %>%
    group_by(Locus_Ref_Limpio, Validacion_Espacial) %>%
    arrange(desc(ref_len_matched), desc(pct_recovered), .by_group = TRUE) %>%
    mutate(
        Es_Redundante = if_else(
            Validacion_Espacial == "TRUE_POSITIVE" & !is.na(Locus_Ref_Limpio) & 
                row_number() > 1, 
            TRUE, FALSE
        ),
        Validacion_Espacial = if_else(Es_Redundante,
                                      "FALSE_POSITIVE", Validacion_Espacial),
        Doble_Control_Biologico = if_else(Es_Redundante, 
                                          "5_TP_Redundante_Degradado", 
                                          Doble_Control_Biologico),
    ) %>%
    ungroup() %>%
    select(-Es_Redundante) %>%
    # Cálculo de los Deltas
    group_by(locus) %>%
    arrange(desc(score), .by_group = TRUE) %>%
    mutate(
        delta_wscore = max(wscore, na.rm = TRUE) - wscore,
        delta_pct_identity = max(pct_identity, na.rm = TRUE) - pct_identity
    ) %>%
    ungroup()

cat("\nResumen S. citrullifolium:\n")
cat("Total de secuencias procesadas (NUC) :", nrow(stats_captus_cit), "\n")
cat("Verdaderos Positivos Únicos (TP)     :", sum(stats_captus_cit$Validacion_Espacial == "TRUE_POSITIVE"), "\n")
cat("Falsos Positivos & Redundantes (FP)  :", sum(stats_captus_cit$Validacion_Espacial == "FALSE_POSITIVE"), "\n")
cat("\nDesglose de Origen Biológico (S. citrullifolium):\n")
print(table(stats_captus_cit$Doble_Control_Biologico))

cat("\nDistribución de hits (S. citrullifolium):\n")
print(
    stats_captus_cit %>% 
        group_by(hit, Validacion_Espacial) %>% 
        count() %>%
        pivot_wider(names_from = Validacion_Espacial, values_from = n, values_fill = 0)
)


## 6. PRUEBA DE TRANSFERENCIA: PREDICCIÓN EN S. CITRULLIFOLIUM ----------------
cat("\n--- INICIANDO PRUEBA DE TRANSFERENCIA FILOGENÉTICA ---\n")

# 1. Preparar la matriz de datos de la nueva especie (S. citrullifolium)
df_prueba_cit <- stats_captus_cit %>%
    mutate(
        Clase_Real = factor(if_else(Validacion_Espacial == "TRUE_POSITIVE", "TP", "FP"),
                            levels = c("TP", "FP")),
        across(all_of(metricas_numericas), ~ as.numeric(as.character(.x)))
    ) %>%
    drop_na(all_of(metricas_numericas)) %>%
    as.data.frame()

# CRÍTICO: Aislamos los predictores exactos para evitar que kernlab colapse
df_prueba_x <- df_prueba_cit[, metricas_numericas]

cat(sprintf("Evaluando %d secuencias en la nueva especie...\n", nrow(df_prueba_cit)))

# 2. Predecir utilizando los modelos base entrenados en S. torvum (usando solo df_prueba_x)
pred_nb_cit <- predict(modelos_balanceados[["NaiveBayes"]], newdata = df_prueba_x)
pred_svm_cit <- predict(modelos_balanceados[["SVM_Polinomial"]], newdata = df_prueba_x)
pred_en_cit <- predict(modelos_balanceados[["ElasticNet"]], newdata = df_prueba_x)

# 3. Predecir utilizando el Meta-Modelo (Stacking)
prob_base_cit <- data.frame(
    NaiveBayes = predict(modelos_balanceados[["NaiveBayes"]], newdata = df_prueba_x, type = "prob")[["TP"]],
    SVM_Polinomial = predict(modelos_balanceados[["SVM_Polinomial"]], newdata = df_prueba_x, type = "prob")[["TP"]]
)

# El GLM supervisor calcula la probabilidad final
pred_stack_prob_cit <- predict(ensamble_stack$ens_model, newdata = prob_base_cit, type = "prob")[["TP"]]
pred_stack_cit <- factor(ifelse(pred_stack_prob_cit > 0.5, "TP", "FP"), levels = c("TP", "FP"))

# 4. Generar las Matrices de Confusión de la Transferencia
cm_nb_cit <- confusionMatrix(data = pred_nb_cit, reference = df_prueba_cit$Clase_Real, mode = "everything")
cm_svm_cit <- confusionMatrix(data = pred_svm_cit, reference = df_prueba_cit$Clase_Real, mode = "everything")
cm_en_cit <- confusionMatrix(data = pred_en_cit, reference = df_prueba_cit$Clase_Real, mode = "everything")
cm_stack_cit <- confusionMatrix(data = pred_stack_cit, reference = df_prueba_cit$Clase_Real, mode = "everything")

# 5. Graficar los resultados
plot_nb_cit <- graficar_matriz_confusion(cm_nb_cit$table, "Naive Bayes",
                                         sprintf("Acc: %.1f%% | FDR: %.1f%%", 
                                                 cm_nb_cit$overall["Accuracy"] * 100, 
                                                 (1 - cm_nb_cit$byClass["Precision"]) * 100))

plot_svm_cit <- graficar_matriz_confusion(cm_svm_cit$table, "SVM Polinomial",
                                          sprintf("Acc: %.1f%% | FDR: %.1f%%", 
                                                  cm_svm_cit$overall["Accuracy"] * 100, 
                                                  (1 - cm_svm_cit$byClass["Precision"]) * 100))

plot_en_cit <- graficar_matriz_confusion(cm_en_cit$table, "ElasticNet",
                                         sprintf("Acc: %.1f%% | FDR: %.1f%%", 
                                                 cm_en_cit$overall["Accuracy"] * 100, 
                                                 (1 - cm_en_cit$byClass["Precision"]) * 100))

plot_stack_cit <- graficar_matriz_confusion(cm_stack_cit$table, "Meta-Modelo Stacking",
                                            sprintf("Acc: %.1f%% | FDR: %.1f%%", 
                                                    cm_stack_cit$overall["Accuracy"] * 100, 
                                                    (1 - cm_stack_cit$byClass["Precision"]) * 100), 
                                            color_alto = "#8e44ad")

print(plot_grid(
    ggdraw() + draw_label("Prueba de Transferencia: S. citrullifolium", fontface = "bold", size = 16),
    plot_grid(plot_nb_cit, plot_en_cit, plot_svm_cit, plot_stack_cit, ncol = 2, align = "v"),
    ncol = 1, rel_heights = c(0.1, 1)
))

# 6. Imprimir resumen en consola
cat("\nResultados Comparativos en la Nueva Especie:\n")
cat(sprintf("ElasticNet     : F1=%.3f | Acc=%.1f%% | Sens=%.1f%% | FDR=%.1f%%\n", 
            cm_en_cit$byClass["F1"], cm_en_cit$overall["Accuracy"]*100, cm_en_cit$byClass["Sensitivity"]*100, (1 - cm_en_cit$byClass["Precision"])*100))
cat(sprintf("Naive Bayes    : F1=%.3f | Acc=%.1f%% | Sens=%.1f%% | FDR=%.1f%%\n", 
            cm_nb_cit$byClass["F1"], cm_nb_cit$overall["Accuracy"]*100, cm_nb_cit$byClass["Sensitivity"]*100, (1 - cm_nb_cit$byClass["Precision"])*100))
cat(sprintf("SVM Polinomial : F1=%.3f | Acc=%.1f%% | Sens=%.1f%% | FDR=%.1f%%\n", 
            cm_svm_cit$byClass["F1"], cm_svm_cit$overall["Accuracy"]*100, cm_svm_cit$byClass["Sensitivity"]*100, (1 - cm_svm_cit$byClass["Precision"])*100))
cat(sprintf("Stacking (GLM) : F1=%.3f | Acc=%.1f%% | Sens=%.1f%% | FDR=%.1f%%\n", 
            cm_stack_cit$byClass["F1"], cm_stack_cit$overall["Accuracy"]*100, cm_stack_cit$byClass["Sensitivity"]*100, (1 - cm_stack_cit$byClass["Precision"])*100))

## 7. INCLUSIÓN DE CORTES YOUDEN Y HEURÍSTICO EN LA TRANSFERENCIA -----------
cat("\n--- EVALUANDO CORTES EMPÍRICOS EN LA NUEVA ESPECIE ---\n")

# 1. Rescatamos los umbrales exactos aprendidos en S. torvum (df_corte)
umbrales_torvum <- list()
for (var in metricas_numericas) {
    roc_obj_temp <- suppressMessages(roc(response = df_corte[[variable_respuesta]], 
                                         predictor = df_corte[[var]], quiet = TRUE))
    umbrales_torvum[[var]] <- pROC::coords(roc_obj_temp, "best", 
                                           ret = "threshold", best.method = "youden")$threshold[1]
}

# 2. Aplicamos estos cortes univariados a los datos nuevos de S. citrullifolium
pred_youden_cit <- list()
for (var in metricas_numericas) {
    pred_binaria <- ifelse(df_prueba_cit[[var]] >= umbrales_torvum[[var]], "TP", "FP")
    pred_youden_cit[[var]] <- factor(pred_binaria, levels = c("TP", "FP"))
}

# 3. Aplicamos el Ensamble Heurístico (Voto Ponderado) aprendido en S. torvum
votos_cit <- data.frame(matrix(ncol = 0, nrow = nrow(df_prueba_cit)))
for (var in metricas_numericas) {
    votos_cit[[var]] <- ifelse(df_prueba_cit[[var]] >= umbrales_torvum[[var]], 1, 0)
}
puntaje_cit <- as.matrix(votos_cit[, metricas_numericas]) %*% mejores_pesos
pred_heuristico_cit <- factor(ifelse(puntaje_cit >= umbral_comb, "TP", "FP"), levels = c("TP", "FP"))

# 4. Generamos las Matrices de Confusión (Top 2 Univariados y el Heurístico)
cm_y_rec_cit <- confusionMatrix(data = pred_youden_cit[["pct_recovered"]], 
                                reference = df_prueba_cit$Clase_Real, mode = "everything")
cm_y_iden_cit <- confusionMatrix(data = pred_youden_cit[["pct_identity"]], 
                                 reference = df_prueba_cit$Clase_Real, mode = "everything")
cm_heuristico_cit <- confusionMatrix(data = pred_heuristico_cit, 
                                     reference = df_prueba_cit$Clase_Real, mode = "everything")

# 5. Imprimimos la Gran Tabla Comparativa de Transferencia en Consola
cat("\n========================================================================\n")
cat("RESULTADOS DE TRANSFERENCIA (S. torvum -> S. citrullifolium)\n")
cat("========================================================================\n")
cat(sprintf("Youden (pct_recovered): F1=%.3f | Acc=%.1f%% | Sens=%.1f%% | FDR=%.1f%%\n", 
            cm_y_rec_cit$byClass["F1"], cm_y_rec_cit$overall["Accuracy"]*100, 
            cm_y_rec_cit$byClass["Sensitivity"]*100, (1 - cm_y_rec_cit$byClass["Precision"])*100))

cat(sprintf("Youden (pct_identity) : F1=%.3f | Acc=%.1f%% | Sens=%.1f%% | FDR=%.1f%%\n", 
            cm_y_iden_cit$byClass["F1"], cm_y_iden_cit$overall["Accuracy"]*100, 
            cm_y_iden_cit$byClass["Sensitivity"]*100, (1 - cm_y_iden_cit$byClass["Precision"])*100))

cat(sprintf("Ensamble Heurístico   : F1=%.3f | Acc=%.1f%% | Sens=%.1f%% | FDR=%.1f%%\n", 
            cm_heuristico_cit$byClass["F1"], cm_heuristico_cit$overall["Accuracy"]*100, 
            cm_heuristico_cit$byClass["Sensitivity"]*100, (1 - cm_heuristico_cit$byClass["Precision"])*100))
cat("------------------------------------------------------------------------\n")
cat(sprintf("ML: ElasticNet        : F1=%.3f | Acc=%.1f%% | Sens=%.1f%% | FDR=%.1f%%\n", 
            cm_en_cit$byClass["F1"], cm_en_cit$overall["Accuracy"]*100, 
            cm_en_cit$byClass["Sensitivity"]*100, (1 - cm_en_cit$byClass["Precision"])*100))

cat(sprintf("ML: Naive Bayes       : F1=%.3f | Acc=%.1f%% | Sens=%.1f%% | FDR=%.1f%%\n", 
            cm_nb_cit$byClass["F1"], cm_nb_cit$overall["Accuracy"]*100, 
            cm_nb_cit$byClass["Sensitivity"]*100, (1 - cm_nb_cit$byClass["Precision"])*100))

cat(sprintf("ML: SVM Polinomial    : F1=%.3f | Acc=%.1f%% | Sens=%.1f%% | FDR=%.1f%%\n", 
            cm_svm_cit$byClass["F1"], cm_svm_cit$overall["Accuracy"]*100, 
            cm_svm_cit$byClass["Sensitivity"]*100, (1 - cm_svm_cit$byClass["Precision"])*100))

cat(sprintf("Meta-Modelo (Stacking): F1=%.3f | Acc=%.1f%% | Sens=%.1f%% | FDR=%.1f%%\n", 
            cm_stack_cit$byClass["F1"], cm_stack_cit$overall["Accuracy"]*100, 
            cm_stack_cit$byClass["Sensitivity"]*100, (1 - cm_stack_cit$byClass["Precision"])*100))
cat("========================================================================\n")

# --- 8. GENERACIÓN DEL GRÁFICO COMPLETO DE TRANSFERENCIA ---

# Generar gráficos individuales para los univariados y el heurístico
plot_y_rec_cit <- graficar_matriz_confusion(cm_y_rec_cit$table, "Youden: PCT_RECOVERED",
                                            sprintf("Acc: %.1f%% | FDR: %.1f%%", 
                                                    cm_y_rec_cit$overall["Accuracy"] * 100, 
                                                    (1 - cm_y_rec_cit$byClass["Precision"]) * 100))

plot_y_iden_cit <- graficar_matriz_confusion(cm_y_iden_cit$table, "Youden: PCT_IDENTITY",
                                             sprintf("Acc: %.1f%% | FDR: %.1f%%", 
                                                     cm_y_iden_cit$overall["Accuracy"] * 100, 
                                                     (1 - cm_y_iden_cit$byClass["Precision"]) * 100))

plot_heuristico_cit <- graficar_matriz_confusion(cm_heuristico_cit$table, "Modelo Heurístico Combinado",
                                                 sprintf("Acc: %.1f%% | FDR: %.1f%%", 
                                                         cm_heuristico_cit$overall["Accuracy"] * 100, 
                                                         (1 - cm_heuristico_cit$byClass["Precision"]) * 100))

# Ensamblar la grilla completa con todos los modelos
matriz_transferencia_grid <- plot_grid(
    plot_y_rec_cit, plot_y_iden_cit, 
    plot_heuristico_cit, plot_en_cit, 
    plot_nb_cit, plot_svm_cit, 
    plot_stack_cit, 
    ncol = 2, align = "v"
)

# Imprimir el gráfico final con título
print(plot_grid(
    ggdraw() + draw_label("Transferencia en S. citrullifolium: Cortes Empíricos vs. ML vs. Stacking", 
                          fontface = "bold", size = 16),
    matriz_transferencia_grid, ncol = 1, rel_heights = c(0.04, 1)
))





# ============================================================================
# EXPERIMENTO DE TRANSFERIBILIDAD 2: S. stellatiglandulosum (Fragmentado)
# ============================================================================

## 1. DICCIONARIO TEÓRICO DE REFERENCIA (PUENTE HOG -> OG) -------------------
archivo <- "Genes_Capture_May_2026.xlsx"
columna_ste <- "S._stellatiglandulosum" 

# 1. Leemos CNV_Perez usando HOG_ID
cnv_ste_raw <- read_excel(archivo, sheet = "CNV_Perez") %>%
    select(HOG_ID, Probe_ID, all_of(columna_ste)) %>%
    filter(!is.na(.data[[columna_ste]]) & trimws(.data[[columna_ste]]) != "")
probes_validos_ste <- cnv_ste_raw$Probe_ID

# 2. Usamos Contigs_probe_design como tabla puente
contigs_probe_ste <- read_excel(archivo, sheet = "Contigs_probe_design") %>%
    select(Probe_ID, Contig_name, HOG_ID_Solanum, OG_ID_Pangenomes,
           Correcte_Contig_Name, Short_name, Pathway, Refseq_RNA,
           DAVID_Annotation) %>%
    filter(Probe_ID %in% probes_validos_ste)

# 3. Reconstruimos cnv_ste inyectando el OG_ID estándar para no romper el pipeline
cnv_ste <- cnv_ste_raw %>%
    left_join(contigs_probe_ste %>% select(Probe_ID, OG_ID = OG_ID_Pangenomes), 
              by = "Probe_ID") %>%
    select(OG_ID, Probe_ID, all_of(columna_ste))

# 4. Leemos secuencias desde Sequences_Torva filtrando por el Assembly de Perez
sequences_ste <- read_excel(archivo, sheet = "Sequences_Torva") %>%
    filter(str_detect(Assembly, "stellati")) %>%
    filter(Probe_ID %in% probes_validos_ste)

# 5. Consolidamos el diccionario unificado
datos_unificados_ste <- cnv_ste %>% select(-OG_ID) %>%
    left_join(contigs_probe_ste, by = "Probe_ID") %>%
    left_join(sequences_ste %>% select(-Orthogroup_ID), by = "Probe_ID")

dict_ste <- datos_unificados_ste %>%
    filter(trimws(Pathway) %in% pathways_filtro) %>%
    select(-all_of(columna_ste)) %>%
    left_join(sequences_ste %>% select(Seq, Orthogroup_ID), by = "Seq") %>%
    filter(!is.na(Contig_ID)) %>%
    select(Locus_Ste = Contig_ID, OG_Ste = Orthogroup_ID) %>%
    mutate(Locus_Ste = str_remove(Locus_Ste, "\\.[0-9]+$")) %>%
    distinct()

## 2. CARGA DEL GFF DE CAPTUS -------------------------------------------------
path_gff_captus_ste <- paste0(ruta_base,
                              "Captus_Local/03_extract_stella_metabolism_complete/",
                              "S_stella__captus-ext/06_assembly_annotated/",
                              "S_stella_hit_contigs.gff")
gff_captus_completo_ste <- rtracklayer::import(path_gff_captus_ste, format = "gff3")
gff_captus_ste <- gff_captus_completo_ste[gff_captus_completo_ste$type == "protein_match:NUC"]

## 3. EXTRACCION DIRECTA DEL LOCUS DE REFERENCIA DESDE EL ATRIBUTO TARGET ----
# Target tiene forma "g65393.t1-OG0000254 1-148 149-448" (transcript_id +
# "-" + OG_ID + coordenadas locales del alineamiento, separadas por espacio).
# Cada hit puede tener varias filas (una por exon del alineamiento) que
# comparten el mismo Name/Target -- se colapsan con distinct() al final.
df_cruces_ste <- data.frame(
    ID_Captus = gff_captus_ste$Name,
    OG_Captus = str_extract(gff_captus_ste$Name, "OG[0-9]+"),
    Locus_Ref_Limpio = str_remove(str_extract(gff_captus_ste$Target, "^\\S+"),
                                  "-OG[0-9]+$"),
    Contig_Name = as.character(seqnames(gff_captus_ste))
) %>%
    mutate(
        Num_Gen = as.numeric(str_extract(Locus_Ref_Limpio, "(?<=^g)[0-9]+"))
    ) %>%
    distinct()

cat(sprintf("Filas unicas tras colapsar exones del mismo hit: %d\n",
            nrow(df_cruces_ste)))

## 4. CLASIFICACION (reemplazo de la mutate de validacion) ------------------

# Tabla puente OG <-> HOG, construida desde el mismo Excel de referencia
puente_og_hog <- read_excel(archivo, sheet = "Contigs_probe_design") %>%
    select(OG_ID = OG_ID_Pangenomes, HOG_ID = HOG_ID_Solanum) %>%
    filter(!is.na(OG_ID), !is.na(HOG_ID)) %>%
    distinct() %>%
    add_count(OG_ID, name = "n_hog") %>%
    group_by(OG_ID) %>%
    slice_max(order_by = n_hog, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(OG_ID, HOG_ID)

# Diagnóstico rápido: confirmar que el puente cubre los OGs de Captus
cat(sprintf("OGs de Captus cubiertos por el puente: %d de %d\n",
            sum(unique(df_cruces_ste$OG_Captus) %in% puente_og_hog$OG_ID),
            length(unique(df_cruces_ste$OG_Captus))))

df_clasificacion_base_ste <- df_cruces_ste %>%
    # 1. Join contra el diccionario para traer el HOG esperado del locus
    left_join(dict_ste, by = c("Locus_Ref_Limpio" = "Locus_Ste")) %>%
    # 2. Join contra el puente para saber a qué HOG equivale el OG de Captus
    left_join(puente_og_hog, by = c("OG_Captus" = "OG_ID"),
              relationship = "many-to-one") %>%
    mutate(
        ID_Captus = if_else(str_detect(ID_Captus, "__[0-9]+$"),
                            ID_Captus,
                            paste0(ID_Captus, "__00")),
        # CORRECCION: el formato es "g65393.t1" -- el numero va entre "g" y "."
        Num_Gen = as.numeric(str_extract(Locus_Ref_Limpio, "(?<=^g)[0-9]+")),
        # Ahora comparamos HOG vs HOG, en el mismo espacio de nombres
        Validacion_Estricta = case_when(
            !is.na(OG_Ste) & !is.na(HOG_ID) & HOG_ID == OG_Ste ~ "TRUE_POSITIVE",
            TRUE ~ "FALSE_POSITIVE"
        )
    )

cat("\nValidacion Estricta (puente OG->HOG):\n")
print(table(df_clasificacion_base_ste$Validacion_Estricta))

# Clasificación Directa Sin Rescate por Vecindad (Criterio Conservador)
df_clasificacion_ste <- df_clasificacion_base_ste %>%
    mutate(
        Validacion_Espacial = Validacion_Estricta,
        Doble_Control_Biologico = case_when(
            Validacion_Estricta == "TRUE_POSITIVE" ~ "1_Validado_Exacto",
            !is.na(OG_Ste) & HOG_ID != OG_Ste     ~ "3_Hit_Cruzado_OG_Vecino",
            TRUE                                   ~ "4_Expansion_No_Anotada_o_Ruido"
        ),
        prioridad = case_when(
            Doble_Control_Biologico == "1_Validado_Exacto" ~ 0L,
            TRUE                                           ~ 1L
        )
    ) %>%
    arrange(prioridad) %>%
    distinct(ID_Captus, .keep_all = TRUE) %>%
    select(ID_Captus, Validacion_Espacial, Doble_Control_Biologico, Locus_Ref_Limpio)

# Auditoría en consola
cat("\nClasificación S. stellatiglandulosum (Sin Rescate por Vecindad):\n")
print(table(df_clasificacion_ste$Validacion_Espacial))
print(table(df_clasificacion_ste$Doble_Control_Biologico))

## 5. INTEGRACIÓN DE ESTADÍSTICAS Y RESOLUCIÓN DE REDUNDANCIAS ---------------
stats_captus_crudo_ste <- read.delim(path_stats_captus, header = TRUE, 
                                     stringsAsFactors = FALSE) %>%
    filter(sample_name == "S_stella") %>%  # CORRECCIÓN: Usamos el nombre exacto de la muestra
    group_by(sample_name, marker_type, ref_name) %>%
    mutate(
        OG_Limpio = str_extract(ref_name, "OG[0-9]+"),
        ID_GFF_Compatible = sprintf("%s__%02d", OG_Limpio, hit)
    ) %>%
    ungroup()

stats_captus_ste <- stats_captus_crudo_ste %>%
    left_join(df_clasificacion_ste, 
              by = c("ID_GFF_Compatible" = "ID_Captus")) %>% 
    filter(marker_type == "NUC") %>%
    mutate(
        Validacion_Espacial = replace_na(Validacion_Espacial, "FALSE_POSITIVE"),
        Doble_Control_Biologico = replace_na(Doble_Control_Biologico, 
                                             "4_Ruido_Intergenico_o_Intron"),
        ref_len_matched = as.numeric(ref_len_matched),
        pct_recovered = as.numeric(pct_recovered)
    ) %>%
    group_by(Locus_Ref_Limpio, Validacion_Espacial) %>%
    arrange(desc(ref_len_matched), desc(pct_recovered), .by_group = TRUE) %>%
    mutate(
        Es_Redundante = if_else(
            Validacion_Espacial == "TRUE_POSITIVE" & !is.na(Locus_Ref_Limpio) & 
                row_number() > 1, 
            TRUE, FALSE
        ),
        Validacion_Espacial = if_else(Es_Redundante, "FALSE_POSITIVE", Validacion_Espacial),
        Doble_Control_Biologico = if_else(Es_Redundante, 
                                          "5_TP_Redundante_Degradado", 
                                          Doble_Control_Biologico)
    ) %>%
    ungroup() %>%
    select(-Es_Redundante) %>%
    group_by(locus) %>%
    arrange(desc(score), .by_group = TRUE) %>%
    mutate(
        delta_wscore = max(wscore, na.rm = TRUE) - wscore,
        delta_pct_identity = max(pct_identity, na.rm = TRUE) - pct_identity
    ) %>%
    ungroup()

cat("\nResumen Final S. stellatiglandulosum (Sin Rescate por Vecindad):\n")
cat("Total procesadas (NUC) :", nrow(stats_captus_ste), "\n")
cat("TP Únicos              :", sum(stats_captus_ste$Validacion_Espacial == "TRUE_POSITIVE"), "\n")
cat("FP & Redundantes       :", sum(stats_captus_ste$Validacion_Espacial == "FALSE_POSITIVE"), "\n")

## 6. PRUEBA DE TRANSFERENCIA: PREDICCIÓN EN S. STELLATIGLANDULOSUM ----------
cat("\n--- INICIANDO PREDICCIÓN EN GENOMA FRAGMENTADO ---\n")

df_prueba_ste <- stats_captus_ste %>%
    mutate(
        Clase_Real = factor(if_else(Validacion_Espacial == "TRUE_POSITIVE", "TP", "FP"),
                            levels = c("TP", "FP")),
        across(all_of(metricas_numericas), ~ as.numeric(as.character(.x)))
    ) %>%
    drop_na(all_of(metricas_numericas)) %>%
    as.data.frame()

# CRÍTICO: Aislar predictores exactos para SVM
df_prueba_x_ste <- df_prueba_ste[, metricas_numericas]

pred_nb_ste <- predict(modelos_balanceados[["NaiveBayes"]], newdata = df_prueba_x_ste)
pred_svm_ste <- predict(modelos_balanceados[["SVM_Polinomial"]], newdata = df_prueba_x_ste)
pred_en_ste <- predict(modelos_balanceados[["ElasticNet"]], newdata = df_prueba_x_ste)

prob_base_ste <- data.frame(
    NaiveBayes = predict(modelos_balanceados[["NaiveBayes"]], newdata = df_prueba_x_ste, type = "prob")[["TP"]],
    SVM_Polinomial = predict(modelos_balanceados[["SVM_Polinomial"]], newdata = df_prueba_x_ste, type = "prob")[["TP"]]
)
pred_stack_prob_ste <- predict(ensamble_stack$ens_model, newdata = prob_base_ste, type = "prob")[["TP"]]
pred_stack_ste <- factor(ifelse(pred_stack_prob_ste > 0.5, "TP", "FP"), levels = c("TP", "FP"))

cm_nb_ste <- confusionMatrix(data = pred_nb_ste, reference = df_prueba_ste$Clase_Real, mode = "everything")
cm_svm_ste <- confusionMatrix(data = pred_svm_ste, reference = df_prueba_ste$Clase_Real, mode = "everything")
cm_en_ste <- confusionMatrix(data = pred_en_ste, reference = df_prueba_ste$Clase_Real, mode = "everything")
cm_stack_ste <- confusionMatrix(data = pred_stack_ste, reference = df_prueba_ste$Clase_Real, mode = "everything")

## 7. CORTES YOUDEN Y HEURÍSTICO ---------------------------------------------
pred_youden_ste <- list()
for (var in metricas_numericas) {
    pred_binaria <- ifelse(df_prueba_ste[[var]] >= umbrales_torvum[[var]], "TP", "FP")
    pred_youden_ste[[var]] <- factor(pred_binaria, levels = c("TP", "FP"))
}

votos_ste <- data.frame(matrix(ncol = 0, nrow = nrow(df_prueba_ste)))
for (var in metricas_numericas) {
    votos_ste[[var]] <- ifelse(df_prueba_ste[[var]] >= umbrales_torvum[[var]], 1, 0)
}
puntaje_ste <- as.matrix(votos_ste[, metricas_numericas]) %*% mejores_pesos
pred_heuristico_ste <- factor(ifelse(puntaje_ste >= umbral_comb, "TP", "FP"), levels = c("TP", "FP"))

cm_y_rec_ste <- confusionMatrix(data = pred_youden_ste[["pct_recovered"]], 
                                reference = df_prueba_ste$Clase_Real, mode = "everything")
cm_y_iden_ste <- confusionMatrix(data = pred_youden_ste[["pct_identity"]], 
                                 reference = df_prueba_ste$Clase_Real, mode = "everything")
cm_heuristico_ste <- confusionMatrix(data = pred_heuristico_ste, 
                                     reference = df_prueba_ste$Clase_Real, mode = "everything")

## 8. IMPRESIÓN Y GRÁFICO FINAL (S. stellatiglandulosum) ---------------------
cat("\n========================================================================\n")
cat("RESULTADOS DE TRANSFERENCIA (S. torvum -> S. stellatiglandulosum)\n")
cat("========================================================================\n")
cat(sprintf("Youden (pct_recovered): F1=%.3f | Acc=%.1f%% | Sens=%.1f%% | FDR=%.1f%%\n", 
            cm_y_rec_ste$byClass["F1"], cm_y_rec_ste$overall["Accuracy"]*100, 
            cm_y_rec_ste$byClass["Sensitivity"]*100, (1 - cm_y_rec_ste$byClass["Precision"])*100))
cat(sprintf("Youden (pct_identity) : F1=%.3f | Acc=%.1f%% | Sens=%.1f%% | FDR=%.1f%%\n", 
            cm_y_iden_ste$byClass["F1"], cm_y_iden_ste$overall["Accuracy"]*100, 
            cm_y_iden_ste$byClass["Sensitivity"]*100, (1 - cm_y_iden_ste$byClass["Precision"])*100))
cat(sprintf("Ensamble Heurístico   : F1=%.3f | Acc=%.1f%% | Sens=%.1f%% | FDR=%.1f%%\n", 
            cm_heuristico_ste$byClass["F1"], cm_heuristico_ste$overall["Accuracy"]*100, 
            cm_heuristico_ste$byClass["Sensitivity"]*100, (1 - cm_heuristico_ste$byClass["Precision"])*100))
cat("------------------------------------------------------------------------\n")
cat(sprintf("ML: ElasticNet        : F1=%.3f | Acc=%.1f%% | Sens=%.1f%% | FDR=%.1f%%\n", 
            cm_en_ste$byClass["F1"], cm_en_ste$overall["Accuracy"]*100, 
            cm_en_ste$byClass["Sensitivity"]*100, (1 - cm_en_ste$byClass["Precision"])*100))
cat(sprintf("ML: Naive Bayes       : F1=%.3f | Acc=%.1f%% | Sens=%.1f%% | FDR=%.1f%%\n", 
            cm_nb_ste$byClass["F1"], cm_nb_ste$overall["Accuracy"]*100, 
            cm_nb_ste$byClass["Sensitivity"]*100, (1 - cm_nb_ste$byClass["Precision"])*100))
cat(sprintf("ML: SVM Polinomial    : F1=%.3f | Acc=%.1f%% | Sens=%.1f%% | FDR=%.1f%%\n", 
            cm_svm_ste$byClass["F1"], cm_svm_ste$overall["Accuracy"]*100, 
            cm_svm_ste$byClass["Sensitivity"]*100, (1 - cm_svm_ste$byClass["Precision"])*100))
cat(sprintf("Meta-Modelo (Stacking): F1=%.3f | Acc=%.1f%% | Sens=%.1f%% | FDR=%.1f%%\n", 
            cm_stack_ste$byClass["F1"], cm_stack_ste$overall["Accuracy"]*100, 
            cm_stack_ste$byClass["Sensitivity"]*100, (1 - cm_stack_ste$byClass["Precision"])*100))
cat("========================================================================\n")

plot_y_rec_ste <- graficar_matriz_confusion(cm_y_rec_ste$table, "Youden: PCT_RECOVERED",
                                            sprintf("Acc: %.1f%% | FDR: %.1f%%", cm_y_rec_ste$overall["Accuracy"] * 100, (1 - cm_y_rec_ste$byClass["Precision"]) * 100))
plot_y_iden_ste <- graficar_matriz_confusion(cm_y_iden_ste$table, "Youden: PCT_IDENTITY",
                                             sprintf("Acc: %.1f%% | FDR: %.1f%%", cm_y_iden_ste$overall["Accuracy"] * 100, (1 - cm_y_iden_ste$byClass["Precision"]) * 100))
plot_heuristico_ste <- graficar_matriz_confusion(cm_heuristico_ste$table, "Modelo Heurístico",
                                                 sprintf("Acc: %.1f%% | FDR: %.1f%%", cm_heuristico_ste$overall["Accuracy"] * 100, (1 - cm_heuristico_ste$byClass["Precision"]) * 100))
plot_nb_ste <- graficar_matriz_confusion(cm_nb_ste$table, "Naive Bayes",
                                         sprintf("Acc: %.1f%% | FDR: %.1f%%", cm_nb_ste$overall["Accuracy"] * 100, (1 - cm_nb_ste$byClass["Precision"]) * 100))
plot_svm_ste <- graficar_matriz_confusion(cm_svm_ste$table, "SVM Polinomial",
                                          sprintf("Acc: %.1f%% | FDR: %.1f%%", cm_svm_ste$overall["Accuracy"] * 100, (1 - cm_svm_ste$byClass["Precision"]) * 100))
plot_en_ste <- graficar_matriz_confusion(cm_en_ste$table, "ElasticNet",
                                         sprintf("Acc: %.1f%% | FDR: %.1f%%", cm_en_ste$overall["Accuracy"] * 100, (1 - cm_en_ste$byClass["Precision"]) * 100))
plot_stack_ste <- graficar_matriz_confusion(cm_stack_ste$table, "Meta-Modelo Stacking",
                                            sprintf("Acc: %.1f%% | FDR: %.1f%%", cm_stack_ste$overall["Accuracy"] * 100, (1 - cm_stack_ste$byClass["Precision"]) * 100), color_alto = "#8e44ad")

matriz_transferencia_grid_ste <- plot_grid(
    plot_y_rec_ste, plot_y_iden_ste, 
    plot_heuristico_ste, plot_en_ste, 
    plot_nb_ste, plot_svm_ste, 
    plot_stack_ste, 
    ncol = 2, align = "v"
)

print(plot_grid(
    ggdraw() + draw_label("Transferencia en S. stellatiglandulosum: Genoma Fragmentado", 
                          fontface = "bold", size = 16),
    matriz_transferencia_grid_ste, ncol = 1, rel_heights = c(0.04, 1)
))

# --- 9. GENERATIO TABULAE COMPARATIVAE (BENCHMARK) PRO S. STELLATIGLANDULOSUM ---

fila_y_rec_ste <- data.frame(
    Modelo = "Corte: pct_recovered",
    Exactitud = as.numeric(cm_y_rec_ste$overall["Accuracy"]),
    Sensibilidad = as.numeric(cm_y_rec_ste$byClass["Sensitivity"]),
    Especificidad = as.numeric(cm_y_rec_ste$byClass["Specificity"]),
    Precision = as.numeric(cm_y_rec_ste$byClass["Precision"]),
    FDR = 1 - as.numeric(cm_y_rec_ste$byClass["Precision"]),
    F1_Score = as.numeric(cm_y_rec_ste$byClass["F1"])
)

fila_y_iden_ste <- data.frame(
    Modelo = "Corte: pct_identity",
    Exactitud = as.numeric(cm_y_iden_ste$overall["Accuracy"]),
    Sensibilidad = as.numeric(cm_y_iden_ste$byClass["Sensitivity"]),
    Especificidad = as.numeric(cm_y_iden_ste$byClass["Specificity"]),
    Precision = as.numeric(cm_y_iden_ste$byClass["Precision"]),
    FDR = 1 - as.numeric(cm_y_iden_ste$byClass["Precision"]),
    F1_Score = as.numeric(cm_y_iden_ste$byClass["F1"])
)

fila_heuristico_ste <- data.frame(
    Modelo = "Modelo Heurístico Combinado",
    Exactitud = as.numeric(cm_heuristico_ste$overall["Accuracy"]),
    Sensibilidad = as.numeric(cm_heuristico_ste$byClass["Sensitivity"]),
    Especificidad = as.numeric(cm_heuristico_ste$byClass["Specificity"]),
    Precision = as.numeric(cm_heuristico_ste$byClass["Precision"]),
    FDR = 1 - as.numeric(cm_heuristico_ste$byClass["Precision"]),
    F1_Score = as.numeric(cm_heuristico_ste$byClass["F1"])
)

fila_en_ste <- data.frame(
    Modelo = "ElasticNet",
    Exactitud = as.numeric(cm_en_ste$overall["Accuracy"]),
    Sensibilidad = as.numeric(cm_en_ste$byClass["Sensitivity"]),
    Especificidad = as.numeric(cm_en_ste$byClass["Specificity"]),
    Precision = as.numeric(cm_en_ste$byClass["Precision"]),
    FDR = 1 - as.numeric(cm_en_ste$byClass["Precision"]),
    F1_Score = as.numeric(cm_en_ste$byClass["F1"])
)

fila_nb_ste <- data.frame(
    Modelo = "NaiveBayes",
    Exactitud = as.numeric(cm_nb_ste$overall["Accuracy"]),
    Sensibilidad = as.numeric(cm_nb_ste$byClass["Sensitivity"]),
    Especificidad = as.numeric(cm_nb_ste$byClass["Specificity"]),
    Precision = as.numeric(cm_nb_ste$byClass["Precision"]),
    FDR = 1 - as.numeric(cm_nb_ste$byClass["Precision"]),
    F1_Score = as.numeric(cm_nb_ste$byClass["F1"])
)

fila_svm_ste <- data.frame(
    Modelo = "SVM_Polinomial",
    Exactitud = as.numeric(cm_svm_ste$overall["Accuracy"]),
    Sensibilidad = as.numeric(cm_svm_ste$byClass["Sensitivity"]),
    Especificidad = as.numeric(cm_svm_ste$byClass["Specificity"]),
    Precision = as.numeric(cm_svm_ste$byClass["Precision"]),
    FDR = 1 - as.numeric(cm_svm_ste$byClass["Precision"]),
    F1_Score = as.numeric(cm_svm_ste$byClass["F1"])
)

fila_stack_ste <- data.frame(
    Modelo = "Ensamble ML (Stacking)",
    Exactitud = as.numeric(cm_stack_ste$overall["Accuracy"]),
    Sensibilidad = as.numeric(cm_stack_ste$byClass["Sensitivity"]),
    Especificidad = as.numeric(cm_stack_ste$byClass["Specificity"]),
    Precision = as.numeric(cm_stack_ste$byClass["Precision"]),
    FDR = 1 - as.numeric(cm_stack_ste$byClass["Precision"]),
    F1_Score = as.numeric(cm_stack_ste$byClass["F1"])
)

df_comparacion_final_ste <- bind_rows(
    fila_y_rec_ste, fila_y_iden_ste,
    fila_heuristico_ste,
    fila_en_ste, fila_nb_ste, fila_svm_ste,
    fila_stack_ste
)

orden_modelos_ste <- c(
    "Ensamble ML (Stacking)",
    "NaiveBayes", "SVM_Polinomial", "ElasticNet",
    "Modelo Heurístico Combinado",
    "Corte: pct_identity", "Corte: pct_recovered"
)

print(
    df_comparacion_final_ste %>%
        pivot_longer(cols = -Modelo, names_to = "Metrica", values_to = "Valor") %>%
        mutate(Modelo = factor(Modelo, levels = rev(orden_modelos_ste))) %>%
        ggplot(aes(x = Modelo, y = Valor, fill = Metrica)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        coord_flip() + scale_fill_npg() +
        labs(title = "Benchmark Definitivo de Estrategias de Clasificación",
             subtitle = "Transferencia en S. stellatiglandulosum (Genoma Fragmentado)",
             x = "Estrategia Evaluada", y = "Puntaje (0 a 1)") +
        theme_bw(base_size = 12) +
        theme(legend.position = "bottom",
              plot.title = element_text(face = "bold"))
)

# --- EXPLORACIÓN: LOS 3 MEJORES MODELOS ESPECÍFICOS PARA S. STELLATIGLANDULOSUM ---
cat("\n--- BÚSQUEDA DINÁMICA DE LOS MEJORES MODELOS EN EL GENOMA FRAGMENTADO ---\n")

# 1. Predecir con TODOS los modelos balanceados y extraer sus métricas
resultados_exploracion <- list()

for (nombre_modelo in names(modelos_balanceados)) {
    modelo <- modelos_balanceados[[nombre_modelo]]
    
    # Generar predicción
    pred_temp <- predict(modelo, newdata = df_prueba_x_ste)
    
    # Calcular matriz de confusión
    cm_temp <- confusionMatrix(data = pred_temp, 
                               reference = df_prueba_ste$Clase_Real, 
                               mode = "everything")
    
    # Extraer métricas clave (manejando posibles NA en F1 si no hay predicciones TP)
    f1_val <- as.numeric(cm_temp$byClass["F1"])
    f1_val <- ifelse(is.na(f1_val), 0, f1_val) # Reemplazar NA por 0 para poder ordenar
    
    resultados_exploracion[[nombre_modelo]] <- data.frame(
        Modelo = nombre_modelo,
        F1_Score = f1_val,
        Exactitud = as.numeric(cm_temp$overall["Accuracy"]),
        Sensibilidad = as.numeric(cm_temp$byClass["Sensitivity"]),
        Precision = as.numeric(cm_temp$byClass["Precision"]),
        FDR = 1 - as.numeric(cm_temp$byClass["Precision"])
    )
}

# 2. Consolidar y ordenar de mejor a peor según F1-Score
df_exploracion <- bind_rows(resultados_exploracion) %>%
    arrange(desc(F1_Score))

cat("\nRanking completo de modelos sobre S. stellatiglandulosum:\n")
print(df_exploracion)

# 3. Extraer los nombres de los 3 mejores
top3_ste_nombres <- df_exploracion$Modelo[1:3]
cat(sprintf("\nLos 3 mejores algoritmos para lidiar con fragmentación son:\n1. %s\n2. %s\n3. %s\n", 
            top3_ste_nombres[1], top3_ste_nombres[2], top3_ste_nombres[3]))

# 4. Graficar los 3 mejores para comparativa visual
plots_top3_ste <- list()
for (nombre in top3_ste_nombres) {
    pred_temp <- predict(modelos_balanceados[[nombre]], newdata = df_prueba_x_ste)
    cm_temp <- confusionMatrix(data = pred_temp, reference = df_prueba_ste$Clase_Real, mode = "everything")
    
    plots_top3_ste[[nombre]] <- graficar_matriz_confusion(
        cm_temp$table, 
        titulo = nombre,
        subtitulo = sprintf("Acc: %.1f%% | Sens: %.1f%% | FDR: %.1f%%", 
                            as.numeric(cm_temp$overall["Accuracy"]) * 100, 
                            as.numeric(cm_temp$byClass["Sensitivity"]) * 100,
                            (1 - as.numeric(cm_temp$byClass["Precision"])) * 100)
    )
}

# Imprimir la grilla gráfica
print(plot_grid(
    ggdraw() + draw_label("Top 3 Algoritmos Específicos para S. stellatiglandulosum", 
                          fontface = "bold", size = 15),
    plot_grid(plotlist = plots_top3_ste, ncol = 3, align = "v"),
    ncol = 1, rel_heights = c(0.1, 1)
))


# ============================================================================
# ACTUALIZACIÓN: 9 MODELOS ESTANDARIZADOS PARA AMBOS ESCENARIOS
# ============================================================================

## A. ESCENARIO 1: S. citrullifolium (Genoma Cromosómico, Distante) -----------

# 1. Nuevas predicciones
pred_dnn_cit <- predict(modelos_balanceados[["DNN_Profunda"]], newdata = df_prueba_x)
pred_lda_cit <- predict(modelos_balanceados[["LDA"]], newdata = df_prueba_x)

cm_dnn_cit <- confusionMatrix(data = pred_dnn_cit, reference = df_prueba_cit$Clase_Real, mode = "everything")
cm_lda_cit <- confusionMatrix(data = pred_lda_cit, reference = df_prueba_cit$Clase_Real, mode = "everything")

# 2. Nuevos gráficos de matrices
plot_dnn_cit <- graficar_matriz_confusion(cm_dnn_cit$table, "DNN Profunda",
                                          sprintf("Acc: %.1f%% | FDR: %.1f%%", cm_dnn_cit$overall["Accuracy"] * 100, (1 - cm_dnn_cit$byClass["Precision"]) * 100))
plot_lda_cit <- graficar_matriz_confusion(cm_lda_cit$table, "LDA",
                                          sprintf("Acc: %.1f%% | FDR: %.1f%%", cm_lda_cit$overall["Accuracy"] * 100, (1 - cm_lda_cit$byClass["Precision"]) * 100))

# 3. Cuadrícula de 9 matrices (3x3)
matriz_9_cit <- plot_grid(
    plot_y_rec_cit, plot_y_iden_cit, plot_heuristico_cit, 
    plot_en_cit, plot_nb_cit, plot_svm_cit, 
    plot_dnn_cit, plot_lda_cit, plot_stack_cit, 
    ncol = 3, align = "v"
)

print(plot_grid(
    ggdraw() + draw_label("Transferencia S. citrullifolium (9 Modelos Estándar)", fontface = "bold", size = 16),
    matriz_9_cit, ncol = 1, rel_heights = c(0.05, 1)
))

# 4. Actualización del Benchmark (S. citrullifolium)
fila_dnn_cit <- data.frame(Modelo = "DNN_Profunda", Exactitud = as.numeric(cm_dnn_cit$overall["Accuracy"]),
                           Sensibilidad = as.numeric(cm_dnn_cit$byClass["Sensitivity"]), Especificidad = as.numeric(cm_dnn_cit$byClass["Specificity"]),
                           Precision = as.numeric(cm_dnn_cit$byClass["Precision"]), FDR = 1 - as.numeric(cm_dnn_cit$byClass["Precision"]), F1_Score = as.numeric(cm_dnn_cit$byClass["F1"]))

fila_lda_cit <- data.frame(Modelo = "LDA", Exactitud = as.numeric(cm_lda_cit$overall["Accuracy"]),
                           Sensibilidad = as.numeric(cm_lda_cit$byClass["Sensitivity"]), Especificidad = as.numeric(cm_lda_cit$byClass["Specificity"]),
                           Precision = as.numeric(cm_lda_cit$byClass["Precision"]), FDR = 1 - as.numeric(cm_lda_cit$byClass["Precision"]), F1_Score = as.numeric(cm_lda_cit$byClass["F1"]))

# Extraer filas existentes de df_comparacion_final (asumiendo que df_comparacion_final original tiene los otros 7)
# Se reconstruye el df manualmente para asegurar limpieza:
df_comp_9_cit <- bind_rows(
    data.frame(Modelo = "Corte: pct_recovered", Exactitud = as.numeric(cm_y_rec_cit$overall["Accuracy"]), Sensibilidad = as.numeric(cm_y_rec_cit$byClass["Sensitivity"]), Especificidad = as.numeric(cm_y_rec_cit$byClass["Specificity"]), Precision = as.numeric(cm_y_rec_cit$byClass["Precision"]), FDR = 1 - as.numeric(cm_y_rec_cit$byClass["Precision"]), F1_Score = as.numeric(cm_y_rec_cit$byClass["F1"])),
    data.frame(Modelo = "Corte: pct_identity", Exactitud = as.numeric(cm_y_iden_cit$overall["Accuracy"]), Sensibilidad = as.numeric(cm_y_iden_cit$byClass["Sensitivity"]), Especificidad = as.numeric(cm_y_iden_cit$byClass["Specificity"]), Precision = as.numeric(cm_y_iden_cit$byClass["Precision"]), FDR = 1 - as.numeric(cm_y_iden_cit$byClass["Precision"]), F1_Score = as.numeric(cm_y_iden_cit$byClass["F1"])),
    data.frame(Modelo = "Modelo Heurístico Combinado", Exactitud = as.numeric(cm_heuristico_cit$overall["Accuracy"]), Sensibilidad = as.numeric(cm_heuristico_cit$byClass["Sensitivity"]), Especificidad = as.numeric(cm_heuristico_cit$byClass["Specificity"]), Precision = as.numeric(cm_heuristico_cit$byClass["Precision"]), FDR = 1 - as.numeric(cm_heuristico_cit$byClass["Precision"]), F1_Score = as.numeric(cm_heuristico_cit$byClass["F1"])),
    data.frame(Modelo = "ElasticNet", Exactitud = as.numeric(cm_en_cit$overall["Accuracy"]), Sensibilidad = as.numeric(cm_en_cit$byClass["Sensitivity"]), Especificidad = as.numeric(cm_en_cit$byClass["Specificity"]), Precision = as.numeric(cm_en_cit$byClass["Precision"]), FDR = 1 - as.numeric(cm_en_cit$byClass["Precision"]), F1_Score = as.numeric(cm_en_cit$byClass["F1"])),
    data.frame(Modelo = "NaiveBayes", Exactitud = as.numeric(cm_nb_cit$overall["Accuracy"]), Sensibilidad = as.numeric(cm_nb_cit$byClass["Sensitivity"]), Especificidad = as.numeric(cm_nb_cit$byClass["Specificity"]), Precision = as.numeric(cm_nb_cit$byClass["Precision"]), FDR = 1 - as.numeric(cm_nb_cit$byClass["Precision"]), F1_Score = as.numeric(cm_nb_cit$byClass["F1"])),
    data.frame(Modelo = "SVM_Polinomial", Exactitud = as.numeric(cm_svm_cit$overall["Accuracy"]), Sensibilidad = as.numeric(cm_svm_cit$byClass["Sensitivity"]), Especificidad = as.numeric(cm_svm_cit$byClass["Specificity"]), Precision = as.numeric(cm_svm_cit$byClass["Precision"]), FDR = 1 - as.numeric(cm_svm_cit$byClass["Precision"]), F1_Score = as.numeric(cm_svm_cit$byClass["F1"])),
    fila_dnn_cit, fila_lda_cit,
    data.frame(Modelo = "Ensamble ML (Stacking)", Exactitud = as.numeric(cm_stack_cit$overall["Accuracy"]), Sensibilidad = as.numeric(cm_stack_cit$byClass["Sensitivity"]), Especificidad = as.numeric(cm_stack_cit$byClass["Specificity"]), Precision = as.numeric(cm_stack_cit$byClass["Precision"]), FDR = 1 - as.numeric(cm_stack_cit$byClass["Precision"]), F1_Score = as.numeric(cm_stack_cit$byClass["F1"]))
)

orden_9_modelos <- c("Ensamble ML (Stacking)", "DNN_Profunda", "LDA", "SVM_Polinomial", "NaiveBayes", "ElasticNet", "Modelo Heurístico Combinado", "Corte: pct_identity", "Corte: pct_recovered")

print(
    df_comp_9_cit %>%
        pivot_longer(cols = -Modelo, names_to = "Metrica", values_to = "Valor") %>%
        mutate(Modelo = factor(Modelo, levels = rev(orden_9_modelos))) %>%
        ggplot(aes(x = Modelo, y = Valor, fill = Metrica)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        coord_flip() + scale_fill_npg() +
        labs(title = "Benchmark de Estrategias: S. citrullifolium (9 Modelos)", x = "Estrategia Evaluada", y = "Puntaje (0 a 1)") +
        theme_bw(base_size = 12) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
)

## B. ESCENARIO 2: S. stellatiglandulosum (Genoma Fragmentado) ----------------

# 1. Nuevas predicciones
pred_dnn_ste <- predict(modelos_balanceados[["DNN_Profunda"]], newdata = df_prueba_x_ste)
pred_lda_ste <- predict(modelos_balanceados[["LDA"]], newdata = df_prueba_x_ste)

cm_dnn_ste <- confusionMatrix(data = pred_dnn_ste, reference = df_prueba_ste$Clase_Real, mode = "everything")
cm_lda_ste <- confusionMatrix(data = pred_lda_ste, reference = df_prueba_ste$Clase_Real, mode = "everything")

# 2. Nuevos gráficos de matrices
plot_dnn_ste <- graficar_matriz_confusion(cm_dnn_ste$table, "DNN Profunda",
                                          sprintf("Acc: %.1f%% | FDR: %.1f%%", cm_dnn_ste$overall["Accuracy"] * 100, (1 - cm_dnn_ste$byClass["Precision"]) * 100))
plot_lda_ste <- graficar_matriz_confusion(cm_lda_ste$table, "LDA",
                                          sprintf("Acc: %.1f%% | FDR: %.1f%%", cm_lda_ste$overall["Accuracy"] * 100, (1 - cm_lda_ste$byClass["Precision"]) * 100))

# 3. Cuadrícula de 9 matrices (3x3)
matriz_9_ste <- plot_grid(
    plot_y_rec_ste, plot_y_iden_ste, plot_heuristico_ste, 
    plot_en_ste, plot_nb_ste, plot_svm_ste, 
    plot_dnn_ste, plot_lda_ste, plot_stack_ste, 
    ncol = 3, align = "v"
)

print(plot_grid(
    ggdraw() + draw_label("Transferencia S. stellatiglandulosum (9 Modelos Estándar)", fontface = "bold", size = 16),
    matriz_9_ste, ncol = 1, rel_heights = c(0.05, 1)
))

# 4. Actualización del Benchmark (S. stellatiglandulosum)
fila_dnn_ste <- data.frame(Modelo = "DNN_Profunda", Exactitud = as.numeric(cm_dnn_ste$overall["Accuracy"]),
                           Sensibilidad = as.numeric(cm_dnn_ste$byClass["Sensitivity"]), Especificidad = as.numeric(cm_dnn_ste$byClass["Specificity"]),
                           Precision = as.numeric(cm_dnn_ste$byClass["Precision"]), FDR = 1 - as.numeric(cm_dnn_ste$byClass["Precision"]), F1_Score = as.numeric(cm_dnn_ste$byClass["F1"]))

fila_lda_ste <- data.frame(Modelo = "LDA", Exactitud = as.numeric(cm_lda_ste$overall["Accuracy"]),
                           Sensibilidad = as.numeric(cm_lda_ste$byClass["Sensitivity"]), Especificidad = as.numeric(cm_lda_ste$byClass["Specificity"]),
                           Precision = as.numeric(cm_lda_ste$byClass["Precision"]), FDR = 1 - as.numeric(cm_lda_ste$byClass["Precision"]), F1_Score = as.numeric(cm_lda_ste$byClass["F1"]))

df_comp_9_ste <- bind_rows(
    data.frame(Modelo = "Corte: pct_recovered", Exactitud = as.numeric(cm_y_rec_ste$overall["Accuracy"]), Sensibilidad = as.numeric(cm_y_rec_ste$byClass["Sensitivity"]), Especificidad = as.numeric(cm_y_rec_ste$byClass["Specificity"]), Precision = as.numeric(cm_y_rec_ste$byClass["Precision"]), FDR = 1 - as.numeric(cm_y_rec_ste$byClass["Precision"]), F1_Score = as.numeric(cm_y_rec_ste$byClass["F1"])),
    data.frame(Modelo = "Corte: pct_identity", Exactitud = as.numeric(cm_y_iden_ste$overall["Accuracy"]), Sensibilidad = as.numeric(cm_y_iden_ste$byClass["Sensitivity"]), Especificidad = as.numeric(cm_y_iden_ste$byClass["Specificity"]), Precision = as.numeric(cm_y_iden_ste$byClass["Precision"]), FDR = 1 - as.numeric(cm_y_iden_ste$byClass["Precision"]), F1_Score = as.numeric(cm_y_iden_ste$byClass["F1"])),
    data.frame(Modelo = "Modelo Heurístico Combinado", Exactitud = as.numeric(cm_heuristico_ste$overall["Accuracy"]), Sensibilidad = as.numeric(cm_heuristico_ste$byClass["Sensitivity"]), Especificidad = as.numeric(cm_heuristico_ste$byClass["Specificity"]), Precision = as.numeric(cm_heuristico_ste$byClass["Precision"]), FDR = 1 - as.numeric(cm_heuristico_ste$byClass["Precision"]), F1_Score = as.numeric(cm_heuristico_ste$byClass["F1"])),
    data.frame(Modelo = "ElasticNet", Exactitud = as.numeric(cm_en_ste$overall["Accuracy"]), Sensibilidad = as.numeric(cm_en_ste$byClass["Sensitivity"]), Especificidad = as.numeric(cm_en_ste$byClass["Specificity"]), Precision = as.numeric(cm_en_ste$byClass["Precision"]), FDR = 1 - as.numeric(cm_en_ste$byClass["Precision"]), F1_Score = as.numeric(cm_en_ste$byClass["F1"])),
    data.frame(Modelo = "NaiveBayes", Exactitud = as.numeric(cm_nb_ste$overall["Accuracy"]), Sensibilidad = as.numeric(cm_nb_ste$byClass["Sensitivity"]), Especificidad = as.numeric(cm_nb_ste$byClass["Specificity"]), Precision = as.numeric(cm_nb_ste$byClass["Precision"]), FDR = 1 - as.numeric(cm_nb_ste$byClass["Precision"]), F1_Score = as.numeric(cm_nb_ste$byClass["F1"])),
    data.frame(Modelo = "SVM_Polinomial", Exactitud = as.numeric(cm_svm_ste$overall["Accuracy"]), Sensibilidad = as.numeric(cm_svm_ste$byClass["Sensitivity"]), Especificidad = as.numeric(cm_svm_ste$byClass["Specificity"]), Precision = as.numeric(cm_svm_ste$byClass["Precision"]), FDR = 1 - as.numeric(cm_svm_ste$byClass["Precision"]), F1_Score = as.numeric(cm_svm_ste$byClass["F1"])),
    fila_dnn_ste, fila_lda_ste,
    data.frame(Modelo = "Ensamble ML (Stacking)", Exactitud = as.numeric(cm_stack_ste$overall["Accuracy"]), Sensibilidad = as.numeric(cm_stack_ste$byClass["Sensitivity"]), Especificidad = as.numeric(cm_stack_ste$byClass["Specificity"]), Precision = as.numeric(cm_stack_ste$byClass["Precision"]), FDR = 1 - as.numeric(cm_stack_ste$byClass["Precision"]), F1_Score = as.numeric(cm_stack_ste$byClass["F1"]))
)

print(
    df_comp_9_ste %>%
        pivot_longer(cols = -Modelo, names_to = "Metrica", values_to = "Valor") %>%
        mutate(Modelo = factor(Modelo, levels = rev(orden_9_modelos))) %>%
        ggplot(aes(x = Modelo, y = Valor, fill = Metrica)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        coord_flip() + scale_fill_npg() +
        labs(title = "Benchmark de Estrategias: S. stellatiglandulosum (9 Modelos)", x = "Estrategia Evaluada", y = "Puntaje (0 a 1)") +
        theme_bw(base_size = 12) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
)

# Imprimir DataFrames en consola para tener los datos de la tesis actualizados
cat("\nTabla Final: S. citrullifolium\n")
print(df_comp_9_cit %>% select(Modelo, F1_Score, Exactitud, Sensibilidad, FDR))

cat("\nTabla Final: S. stellatiglandulosum\n")
print(df_comp_9_ste %>% select(Modelo, F1_Score, Exactitud, Sensibilidad, FDR))








# 1. CARGA Y PREPROCESAMIENTO DE NUEVOS DATOS --------------------------------
path_nuevos_datos <- paste0(ruta_base, "CAPTUS_Muestras/captus-extract_stats.tsv") 

if(file.exists(path_nuevos_datos)){
    datos_nuevos_crudos <- read.delim(path_nuevos_datos, header = TRUE, 
                                      stringsAsFactors = FALSE)
} else {
    cat("\n[!] Archivo no encontrado. Usando S. citrullifolium como demostración.\n")
    datos_nuevos_crudos <- stats_captus_crudo_cit 
}

# Formateo y limpieza estructural
df_inferencia <- datos_nuevos_crudos %>%
    filter(marker_type == "NUC") %>%
    mutate(
        OG_Limpio = str_extract(ref_name, "OG[0-9]+"),
        ID_Secuencia = case_when(
            n() == 1 ~ OG_Limpio,
            n() > 1  ~ sprintf("%s__%02d", OG_Limpio, hit)
        ),
        across(all_of(metricas_numericas), ~ as.numeric(as.character(.x)))
    ) %>%
    drop_na(all_of(metricas_numericas)) %>% 
    as.data.frame()

df_predictores_nuevos <- df_inferencia %>% select(all_of(metricas_numericas))

cat(sprintf("Procesando %d secuencias candidatas de %s...\n", 
            nrow(df_inferencia), unique(df_inferencia$sample_name)[1]))

df_resultados_finales <- data.frame(
    Sample = df_inferencia$sample_name,
    ID_Secuencia = df_inferencia$ID_Secuencia,
    Locus_Ortogrupo = df_inferencia$locus,
    pct_recovered = df_inferencia$pct_recovered
)

# 2. APLICACIÓN DE LOS 4 MEJORES MODELOS DE MACHINE LEARNING -----------------
# 2 de S. citrullifolium (ElasticNet, NaiveBayes) + 2 de S. stellatiglandulosum (DNN_Profunda, SVM_Polinomial)
modelos_consenso <- c("ElasticNet", "NaiveBayes", "DNN_Profunda", "SVM_Polinomial")
cat("Aplicando los 4 modelos estadísticos base...\n")

for(mod_name in modelos_consenso) {
    df_resultados_finales[[paste0("Pred_ML_", mod_name)]] <- predict(
        modelos_balanceados[[mod_name]], newdata = df_predictores_nuevos
    )
}

# 3. APLICACIÓN DE LOS 2 CORTES HEURÍSTICOS (YOUDEN) -------------------------
cat("Aplicando los 2 cortes heurísticos univariados (Youden)...\n")
vars_youden <- c("pct_recovered", "pct_identity")

for (var in vars_youden) {
    roc_obj_temp <- suppressMessages(roc(df_corte[[variable_respuesta]], df_corte[[var]], quiet = TRUE))
    umbral_temp <- pROC::coords(roc_obj_temp, "best", ret = "threshold", best.method = "youden")$threshold[1]
    
    pred_binaria <- ifelse(df_inferencia[[var]] >= umbral_temp, 1, 0)
    df_resultados_finales[[paste0("Pred_Youden_", var)]] <- factor(ifelse(pred_binaria == 1, "TP", "FP"), levels = c("TP", "FP"))
}

# 4. APLICACIÓN DEL ENSAMBLE HEURÍSTICO COMBINADO ----------------------------
cat("Aplicando Ensamble Heurístico Combinado...\n")
votos_nuevos <- data.frame(matrix(ncol = 0, nrow = nrow(df_inferencia)))

for (var in metricas_numericas) {
    roc_obj_temp <- suppressMessages(roc(df_corte[[variable_respuesta]], df_corte[[var]], quiet = TRUE))
    umbral_temp <- pROC::coords(roc_obj_temp, "best", ret = "threshold", best.method = "youden")$threshold[1]
    votos_nuevos[[var]] <- ifelse(df_inferencia[[var]] >= umbral_temp, 1, 0)
}

puntaje_nuevos <- as.matrix(votos_nuevos[, metricas_numericas]) %*% mejores_pesos
df_resultados_finales$Pred_Ens_Heuristico <- factor(ifelse(puntaje_nuevos >= umbral_comb, "TP", "FP"), levels = c("TP", "FP"))

# 5. GENERACIÓN DEL CATÁLOGO (CONSENSO DE 7 ESTRATEGIAS) ---------------------
cat("Calculando veredicto de consenso (Umbral: >= 4/7)...\n")
df_resultados_finales <- df_resultados_finales %>%
    mutate(
        Votos_TP = (Pred_Ens_Heuristico == "TP") + 
            (Pred_ML_ElasticNet == "TP") +
            (Pred_ML_NaiveBayes == "TP") +
            (Pred_ML_DNN_Profunda == "TP") +
            (Pred_ML_SVM_Polinomial == "TP") +
            (Pred_Youden_pct_recovered == "TP") +
            (Pred_Youden_pct_identity == "TP"),
        Veredicto_Consenso = if_else(Votos_TP >= 4, "Ortologo_Valido", "Ruido_Descartado")
    )

# 6. FILTRADO BIOLÓGICO Y EXTRACCIÓN DE CNV ----------------------------------
cat("Ejecutando control de fragmentación técnica y extracción de matriz CNV...\n")
df_cnv_bruto <- df_resultados_finales %>%
    filter(Veredicto_Consenso == "Ortologo_Valido") %>%
    group_by(Sample, Locus_Ortogrupo) %>%
    arrange(desc(pct_recovered), .by_group = TRUE) %>%
    mutate(
        cobertura_acumulada = cumsum(pct_recovered),
        Clasificacion_Copia = case_when(
            row_number() == 1 ~ "Copia_Principal",
            row_number() > 1 & pct_recovered >= 50.0 ~ "Paralogo_Confirmado",
            row_number() > 1 & pct_recovered < 50.0 & cobertura_acumulada <= 120.0 ~ "Fragmentacion_de_Ensamblaje",
            TRUE ~ "Paralogo_Dudoso"
        )
    ) %>%
    ungroup()

df_matriz_cnv <- df_cnv_bruto %>%
    filter(Clasificacion_Copia %in% c("Copia_Principal", "Paralogo_Confirmado")) %>%
    group_by(Sample, Locus_Ortogrupo) %>%
    summarise(CNV_Detectado = n(), .groups = "drop") %>%
    pivot_wider(names_from = Sample, values_from = CNV_Detectado, values_fill = 0)

catalogo_genes_funcionales <- df_cnv_bruto %>% 
    filter(Clasificacion_Copia %in% c("Copia_Principal", "Paralogo_Confirmado")) %>%
    arrange(Locus_Ortogrupo)

cat("\n=======================================================\n")
cat("RESULTADOS DEL DESPLIEGUE FINAL (7 MODELOS)\n")
cat("=======================================================\n")
print(table(df_resultados_finales$Veredicto_Consenso))
cat(sprintf("\nTotal de copias funcionales rescatadas: %d\n", nrow(catalogo_genes_funcionales)))

# 7. GENERACIÓN DE GRÁFICOS (FIGURAS 27 Y 28) --------------------------------
library(ggplot2)
library(viridis)

# Gráfico A: Distribución del Consenso
plot_consenso <- ggplot(df_resultados_finales, aes(x = Veredicto_Consenso, fill = Veredicto_Consenso)) +
    geom_bar(color = "black", alpha = 0.8) +
    scale_fill_manual(values = c("Ortologo_Valido" = "#27ae60", "Ruido_Descartado" = "#c0392b")) +
    geom_text(stat='count', aes(label=..count..), vjust=-0.5, size=5) +
    theme_bw(base_size = 14) +
    labs(title = "Clasificación de Secuencias en Muestras Silvestres",
         subtitle = "Sistema de Consenso (Mayoría: >= 4/7 Modelos)",
         x = "Veredicto del Ensamble",
         y = "Número de Secuencias", fill = "Categoría") +
    theme(legend.position = "none")

# Gráfico B: Heatmap de Variación Estructural (Top 30)
top_loci_cnv <- df_matriz_cnv %>%
    rowwise() %>%
    mutate(Total_Copias = sum(c_across(-Locus_Ortogrupo))) %>%
    ungroup() %>%
    arrange(desc(Total_Copias)) %>%
    slice_head(n = 30) %>%
    pull(Locus_Ortogrupo)

df_heatmap_cnv <- df_cnv_bruto %>%
    filter(Clasificacion_Copia %in% c("Copia_Principal", "Paralogo_Confirmado")) %>%
    filter(Locus_Ortogrupo %in% top_loci_cnv) %>%
    group_by(Sample, Locus_Ortogrupo) %>%
    summarise(Numero_Copias = n(), .groups = "drop")

plot_heatmap <- ggplot(df_heatmap_cnv, aes(x = Sample, y = Locus_Ortogrupo, fill = Numero_Copias)) +
    geom_tile(color = "white", size = 0.1) +
    scale_fill_viridis(option = "mako", direction = -1, breaks = c(1, 2, 3, 4, 5)) +
    theme_minimal(base_size = 12) +
    labs(title = "Matriz de CNV - Top 30 Familias Expandidas",
         subtitle = "Filtrado biológico y control de fragmentación",
         x = "Muestras Silvestres", y = "Locus / Ortogrupo", fill = "Copias") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank())

print(plot_consenso)
print(plot_heatmap)






























## 8. DESPLIEGUE: INFERENCIA DE CNV EN ESPECIES DE TARGET CAPTURE ------------

# 1. CARGA Y PREPROCESAMIENTO DE NUEVOS DATOS
# Apuntamos directamente a la lectura del workbook Excel
path_nuevos_datos <- paste0(ruta_base, 
                            "CAPTUS_Muestras/captus-extract_stats.tsv") 

if(file.exists(path_nuevos_datos)){
    datos_nuevos_crudos <- read.delim(path_nuevos_datos, header = TRUE, 
                                      stringsAsFactors = FALSE)
} else {
    cat("\n[!] Archivo Excel no encontrado. Usando S. citrullifolium puro como demostración.\n")
    datos_nuevos_crudos <- stats_captus_crudo_cit 
}

# Formateo y limpieza estructural
df_inferencia <- datos_nuevos_crudos %>%
    filter(marker_type == "NUC") %>%
    mutate(
        OG_Limpio = str_extract(ref_name, "OG[0-9]+"),
        ID_Secuencia = case_when(
            n() == 1 ~ OG_Limpio,
            n() > 1  ~ sprintf("%s__%02d", OG_Limpio, hit)
        ),
        across(all_of(metricas_numericas), ~ as.numeric(as.character(.x)))
    ) %>%
    drop_na(all_of(metricas_numericas)) %>% 
    as.data.frame()

# Matriz pura para alimentar los modelos
df_predictores_nuevos <- df_inferencia %>% select(all_of(metricas_numericas))

cat(sprintf("Procesando %d secuencias candidatas de %s...\n", 
            nrow(df_inferencia), unique(df_inferencia$sample_name)[1]))

# Tabla base de resultados (incluimos pct_recovered directo para evitar joins posteriores)
df_resultados_finales <- data.frame(
    Sample = df_inferencia$sample_name,
    ID_Secuencia = df_inferencia$ID_Secuencia,
    Locus_Ortogrupo = df_inferencia$locus,
    pct_recovered = df_inferencia$pct_recovered
)


# 2. APLICACIÓN DE LOS 3 MEJORES MODELOS DE MACHINE LEARNING
cat("Aplicando los 3 mejores modelos ML individuales...\n")
mejores_3_ml <- names(sort(roc_promedios, decreasing = TRUE)[1:3])

for(mod_name in mejores_3_ml) {
    df_resultados_finales[[paste0("Pred_ML_", mod_name)]] <- predict(
        modelos_balanceados[[mod_name]], newdata = df_predictores_nuevos
    )
}

# 3. APLICACIÓN DE LOS 2 MEJORES MODELOS HEURÍSTICOS (UNIVARIADOS - YOUDEN)
cat("Aplicando los 2 mejores cortes heurísticos (Youden)...\n")
# Extraemos el nombre limpio de la variable desde top2_youden
mejores_2_youden_vars <- str_replace(top2_youden, "Youden_", "") 

for (var in mejores_2_youden_vars) {
    roc_obj_temp <- suppressMessages(roc(df_corte[[variable_respuesta]], df_corte[[var]], quiet = TRUE))
    umbral_temp <- pROC::coords(roc_obj_temp, "best", ret = "threshold", best.method = "youden")$threshold[1]
    
    pred_binaria <- ifelse(df_inferencia[[var]] >= umbral_temp, 1, 0)
    df_resultados_finales[[paste0("Pred_Youden_", var)]] <- factor(ifelse(pred_binaria == 1, "TP", "FP"), levels = c("TP", "FP"))
}


# 4. APLICACIÓN DEL ENSAMBLE HEURÍSTICO (VOTO PONDERADO)
cat("Aplicando Ensamble Heurístico Combinado...\n")
votos_nuevos <- data.frame(matrix(ncol = 0, nrow = nrow(df_inferencia)))

for (var in metricas_numericas) {
    roc_obj_temp <- suppressMessages(roc(df_corte[[variable_respuesta]], df_corte[[var]], quiet = TRUE))
    umbral_temp <- pROC::coords(roc_obj_temp, "best", ret = "threshold", best.method = "youden")$threshold[1]
    votos_nuevos[[var]] <- ifelse(df_inferencia[[var]] >= umbral_temp, 1, 0)
}

puntaje_nuevos <- as.matrix(votos_nuevos[, metricas_numericas]) %*% mejores_pesos
df_resultados_finales$Pred_Ens_Heuristico <- factor(ifelse(puntaje_nuevos >= umbral_comb, "TP", "FP"), levels = c("TP", "FP"))


# 5. GENERACIÓN DEL CATÁLOGO DE ORTÓLOGOS (CONSENSO DE 6 MODELOS)
cat("Calculando veredicto de consenso...\n")
# Consenso: Si 4 de los 6 modelos dictaminan TP (mayoría), se considera Ortólogo Válido
df_resultados_finales <- df_resultados_finales %>%
    mutate(
        Votos_TP = (Pred_Ens_Heuristico == "TP") + 
            (.[[paste0("Pred_ML_", mejores_3_ml[1])]] == "TP") +
            (.[[paste0("Pred_ML_", mejores_3_ml[2])]] == "TP") +
            (.[[paste0("Pred_ML_", mejores_3_ml[3])]] == "TP") +
            (.[[paste0("Pred_Youden_", mejores_2_youden_vars[1])]] == "TP") +
            (.[[paste0("Pred_Youden_", mejores_2_youden_vars[2])]] == "TP"),
        Veredicto_Consenso = if_else(Votos_TP >= 4, "Ortologo_Valido", "Ruido_Descartado")
    )


# 6. FILTRADO BIOLÓGICO, CONTROL DE FRAGMENTACIÓN Y EXTRACCIÓN DE CNV
cat("Ejecutando rescate de parálogos sin GFF y extracción de matriz CNV...\n")

df_cnv_bruto <- df_resultados_finales %>%
    filter(Veredicto_Consenso == "Ortologo_Valido") %>%
    group_by(Sample, Locus_Ortogrupo) %>%
    # Ordenamos por cobertura (de mayor a menor) para evaluar las acumulaciones
    arrange(desc(pct_recovered), .by_group = TRUE) %>%
    mutate(
        # LÓGICA DE RESCATE SIN GFF:
        cobertura_acumulada = cumsum(pct_recovered),
        
        Clasificacion_Copia = case_when(
            row_number() == 1 ~ "Copia_Principal",
            row_number() > 1 & pct_recovered >= 50.0 ~ "Paralogo_Confirmado",
            row_number() > 1 & pct_recovered < 50.0 & cobertura_acumulada <= 120.0 ~ "Fragmentacion_de_Ensamblaje",
            TRUE ~ "Paralogo_Dudoso"
        )
    ) %>%
    ungroup()

# Extraemos el número de copias final (CNV) por Ortogrupo para visualización y exportación
df_matriz_cnv <- df_cnv_bruto %>%
    filter(Clasificacion_Copia %in% c("Copia_Principal", "Paralogo_Confirmado")) %>%
    group_by(Sample, Locus_Ortogrupo) %>%
    summarise(
        CNV_Detectado = n(),
        .groups = "drop"
    ) %>%
    pivot_wider(names_from = Sample, values_from = CNV_Detectado, values_fill = 0)

catalogo_genes_funcionales <- df_cnv_bruto %>% 
    filter(Clasificacion_Copia %in% c("Copia_Principal", "Paralogo_Confirmado")) %>%
    arrange(Locus_Ortogrupo)

# Resumen de salida
cat("\n=======================================================\n")
cat("RESULTADOS DEL DESPLIEGUE FINAL (INFERENCIA)\n")
cat("=======================================================\n")
print(table(df_resultados_finales$Veredicto_Consenso))

cat(sprintf("\n¡Listo! Se han rescatado %d copias génicas funcionales en la nueva especie.\n", 
            nrow(catalogo_genes_funcionales)))

cat("\n--- MATRIZ FINAL DE CNV (Vista Previa) ---\n")
print(head(df_matriz_cnv))

# ============================================================================
# 7. GENERACIÓN DE GRÁFICOS PARA LA SECCIÓN 4.5
# ============================================================================
library(ggplot2)
library(viridis)
library(tidyr)

# Gráfico A: Distribución del Veredicto de Consenso
plot_consenso <- ggplot(df_resultados_finales, aes(x = Veredicto_Consenso, fill = Veredicto_Consenso)) +
    geom_bar(color = "black", alpha = 0.8) +
    scale_fill_manual(values = c("Ortologo_Valido" = "#27ae60", "Ruido_Descartado" = "#c0392b")) +
    geom_text(stat='count', aes(label=..count..), vjust=-0.5, size=5) +
    theme_bw(base_size = 14) +
    labs(title = "Clasificación de Secuencias en 105 Muestras Silvestres",
         subtitle = "Sistema de Consenso (Mayoría Simple: >= 4/6 Modelos)",
         x = "Veredicto del Ensamble",
         y = "Número Total de Secuencias",
         fill = "Categoría") +
    theme(legend.position = "none")

print(plot_consenso)

# Gráfico B: Heatmap de Variación en el Número de Copias (Top 30 Loci más variables)
# Identificamos los 30 ortogrupos con mayor acumulación de copias para facilitar la visualización
top_loci_cnv <- df_matriz_cnv %>%
    rowwise() %>%
    mutate(Total_Copias = sum(c_across(-Locus_Ortogrupo))) %>%
    ungroup() %>%
    arrange(desc(Total_Copias)) %>%
    slice_head(n = 30) %>%
    pull(Locus_Ortogrupo)

df_heatmap_cnv <- df_cnv_bruto %>%
    filter(Clasificacion_Copia %in% c("Copia_Principal", "Paralogo_Confirmado")) %>%
    filter(Locus_Ortogrupo %in% top_loci_cnv) %>%
    group_by(Sample, Locus_Ortogrupo) %>%
    summarise(Numero_Copias = n(), .groups = "drop")

plot_heatmap <- ggplot(df_heatmap_cnv, aes(x = Sample, y = Locus_Ortogrupo, fill = Numero_Copias)) +
    geom_tile(color = "white", size = 0.1) +
    scale_fill_viridis(option = "mako", direction = -1, breaks = c(1, 2, 3, 4, 5)) +
    theme_minimal(base_size = 12) +
    labs(title = "Matriz de Variación Estructural (CNV) - Top 30 Familias",
         subtitle = "Detección post-control de fragmentación técnica",
         x = "Muestras de Target Capture",
         y = "Locus / Ortogrupo",
         fill = "Copias") +
    theme(axis.text.x = element_blank(), # Oculta los nombres de las 105 muestras por espacio
          axis.ticks.x = element_blank(),
          panel.grid = element_blank())

print(plot_heatmap)
