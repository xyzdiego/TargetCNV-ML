# ============================================================================
# FIX DEFINITIVO: usar el atributo Target de Captus como llave directa
#
# DIAGNOSTICO CONFIRMADO:
#   dict_ste$Locus_Ste == transcript_id crudo del GTF, 51/51 coincidencias
#   exactas (NO gene_id a nivel gen). Y el atributo Target del GFF de
#   Captus ("g65393.t1-OG0000254 1-148 149-448") ya contiene ese mismo
#   transcript_id, seguido del OG con el que se diseño la sonda.
#
#   Conclusion: Captus/Scipio ya resuelve, hit por hit, a que transcrito
#   de referencia alineo. No hace falta reconstruir genes desde el GTF ni
#   hacer findOverlaps espacial -- eso solo introducia la ambiguedad
#   AUGUSTUS/gmst sin necesidad. Se reemplaza el cruce espacial por
#   extraccion directa del Target.
#
#   (El GTF y genes_ref_ste ya no son necesarios para este paso de
#   clasificacion TP/FP; consérvalos solo si los usas en otro punto del
#   pipeline, p. ej. si necesitas los modelos genicos completos para algo
#   distinto al etiquetado TP/FP.)
# ============================================================================

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

## 4. CLASIFICACION Y RESCATE POR VECINDAD EN CONTIG -------------------------
df_clasificacion_base_ste <- df_cruces_ste %>%
    left_join(dict_ste, by = c("Locus_Ref_Limpio" = "Locus_Ste")) %>%
    mutate(
        ID_Captus = if_else(str_detect(ID_Captus, "__[0-9]+$"),
                            ID_Captus,
                            paste0(ID_Captus, "__00")),
        Num_OG_Captus = str_extract(OG_Captus, "[0-9]+"),
        Num_OG_Ref = str_extract(OG_Ste, "[0-9]+"),
        Validacion_Estricta = case_when(
            !is.na(Num_OG_Ref) & Num_OG_Captus == Num_OG_Ref ~ "TRUE_POSITIVE",
            TRUE ~ "FALSE_POSITIVE"
        )
    )

cat("\nValidacion Estricta (antes del rescate por tandem):\n")
print(table(df_clasificacion_base_ste$Validacion_Estricta))

# Rescate por vecindad en tandem -- misma logica que ya tenias, ahora sobre
# Contig_Name/Num_Gen derivados directamente del propio hit de Captus
df_clasificacion_rescate_ste <- df_clasificacion_base_ste %>%
    group_by(OG_Captus) %>%
    mutate(
        Tiene_Ancla_TP = any(Validacion_Estricta == "TRUE_POSITIVE"),
        Contigs_TP = list(Contig_Name[Validacion_Estricta == "TRUE_POSITIVE"]),
        Posiciones_TP = list(Num_Gen[Validacion_Estricta == "TRUE_POSITIVE"])
    ) %>%
    rowwise() %>%
    mutate(
        Es_Tandem = if_else(
            Validacion_Estricta == "FALSE_POSITIVE" &
                Tiene_Ancla_TP &
                (Contig_Name %in% unlist(Contigs_TP)),
            any(abs(Num_Gen - unlist(Posiciones_TP)[unlist(Contigs_TP) == Contig_Name]) <= 50),
            FALSE
        ),
        Validacion_Espacial = case_when(
            Validacion_Estricta == "TRUE_POSITIVE" ~ "TRUE_POSITIVE",
            Es_Tandem ~ "TRUE_POSITIVE",
            TRUE ~ "FALSE_POSITIVE"
        ),
        Doble_Control_Biologico = case_when(
            Validacion_Estricta == "TRUE_POSITIVE" ~ "1_Validado_Exacto",
            Es_Tandem ~ "2_Rescatado_Paralogo_Tandem",
            !is.na(Num_OG_Ref) & Num_OG_Captus != Num_OG_Ref ~ "3_Hit_Cruzado_OG_Vecino",
            TRUE ~ "4_Expansion_No_Anotada_o_Ruido"
        ),
        prioridad = case_when(
            Doble_Control_Biologico == "1_Validado_Exacto" ~ 0L,
            Doble_Control_Biologico == "2_Rescatado_Paralogo_Tandem" ~ 1L,
            TRUE ~ 2L
        )
    ) %>%
    ungroup() %>%
    arrange(prioridad) %>%
    distinct(ID_Captus, .keep_all = TRUE) %>%
    select(ID_Captus, Validacion_Espacial, Doble_Control_Biologico, Locus_Ref_Limpio)

cat("\nAuditoria de Rescate Espacial (S. stellatiglandulosum) -- via Target:\n")
print(table(df_clasificacion_rescate_ste$Doble_Control_Biologico))

# A partir de aqui, tu Seccion 5 (integracion de estadisticas / TSV) y
# en adelante NO necesitan cambios: siguen consumiendo df_clasificacion_rescate_ste
# con las mismas columnas de siempre (ID_Captus, Validacion_Espacial,
# Doble_Control_Biologico, Locus_Ref_Limpio).
