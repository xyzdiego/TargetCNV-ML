# ============================================================================
# DIAGNOSTICO: comparar el formato de ID en las 3 fuentes que deben cruzar
#              para S. stellatiglandulosum
#
# TP = 0 / FP = 100% con "4_Expansion_No_Anotada_o_Ruido" en TODAS las filas
# significa que el join contra dict_ste fallo por completo (namespace
# distinto), no que la clasificacion sea parcialmente incorrecta.
#
# Este script no modifica nada, solo imprime muestras para comparar formatos.
# ============================================================================

cat("==================================================================\n")
cat("1) Locus_Ste (diccionario armado desde el Excel)\n")
cat("==================================================================\n")
print(head(dict_ste$Locus_Ste, 15))
cat(sprintf("Total de valores unicos en dict_ste$Locus_Ste: %d\n",
            length(unique(dict_ste$Locus_Ste))))

cat("\n==================================================================\n")
cat("2) gene_id crudo del GTF (BRAKER3), sin contig ni fuente agregados\n")
cat("==================================================================\n")
print(head(unique(exones_ref_ste$gene_id), 15))

cat("\n==================================================================\n")
cat("3) transcript_id crudo del GTF (por si dict_ste usa nivel transcrito,\n")
cat("   no nivel gen)\n")
cat("==================================================================\n")
print(head(unique(exones_ref_ste$transcript_id), 15))

cat("\n==================================================================\n")
cat("4) Atributo Target del GFF de Captus (revela el ID de referencia que\n")
cat("   Scipio/Captus uso internamente para el alineamiento)\n")
cat("==================================================================\n")
if ("Target" %in% names(mcols(gff_captus_ste))) {
    print(head(gff_captus_ste$Target, 15))
} else {
    cat("La columna 'Target' no se parseo como mcol independiente.\n")
    cat("Revisa el atributo crudo:\n")
    print(head(gff_captus_completo_ste$Attributes, 15))  # fallback si existe
}

cat("\n==================================================================\n")
cat("5) Interseccion directa: cuantos valores de Locus_Ste aparecen TAL\n")
cat("   CUAL dentro del gene_id crudo del GTF (sin contig/fuente)\n")
cat("==================================================================\n")
cat(sprintf("Coinciden exactas (Locus_Ste == gene_id crudo): %d de %d\n",
            sum(dict_ste$Locus_Ste %in% unique(exones_ref_ste$gene_id)),
            length(unique(dict_ste$Locus_Ste))))
cat(sprintf("Coinciden exactas (Locus_Ste == transcript_id crudo): %d de %d\n",
            sum(dict_ste$Locus_Ste %in% unique(exones_ref_ste$transcript_id)),
            length(unique(dict_ste$Locus_Ste))))
