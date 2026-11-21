# Declaración de variables de evaluación no estándar (dplyr / ggplot2) para R CMD check
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".", ".Locus_Tmp", ".Sample_Tmp",
    "Assembly", "CNV_Detectado", "CNV_Estimado", "CNV_Real", "Categoria",
    "Clase_Factor", "Clase_Real", "Clase_Simple", "Clasificacion_Copia",
    "Conteo", "Contig_ID", "Contig_name", "Copias_Totales_Estimadas",
    "Correcte_Contig_Name", "DAVID_Annotation", "Doble_Control_Biologico",
    "Es_Redundante", "Escala_Genoma", "Especie", "Especificidad", "Etiqueta",
    "Etiqueta_Celda", "Etiqueta_Especie", "F1_Score", "Fragmentos_No_Contados",
    "Frecuencia", "Freq", "HOG_ID", "HOG_ID_Solanum", "HOG_Limpio",
    "Homologos_Confirmados", "ID_Captus", "ID_GFF_Compatible", "Importancia",
    "Label", "Locus", "Locus_Benoit", "Locus_Dict_Limpio", "Locus_Ortogrupo",
    "Locus_Ref_Limpio", "Metrica", "Metrica_Label", "Modelo", "Muestra",
    "Numero_Copias", "OG_Captus", "OG_Dict_Limpio", "OG_ID", "OG_ID_Pangenomes",
    "OG_Limpio", "Orthogroup_ID", "P_Val", "Pathway", "Prediccion", "Probe_ID",
    "Puntaje", "ROC_AUC", "Referencia", "Refseq_RNA", "Rho", "S._citrullifolium",
    "Sample", "Sensibilidad", "Seq", "Short_name", "Spearman_Rho", "Total",
    "Total_Copias", "Umbral", "Umbral_Cobertura", "Validacion_Espacial",
    "Validacion_Estricta", "Valor", "Var1", "Var2", "Variable",
    "Veredicto_Consenso", "k", "locus", "marker_type", "med_label", "mediana",
    "n_hog", "n_label", "n_obs", "node_label", "pct_identity", "pct_recovered",
    "prioridad", "ref_len_matched", "ref_name", "sample_name", "score", "tipo",
    "wscore", "x", "x_num", "x_poly", "x_rain", "xend", "y", "y_n_pos",
    "y_poly", "yend"
  ))
}
