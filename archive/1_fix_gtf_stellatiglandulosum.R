# ============================================================================
# FIX: Cruce de coordenadas GFF (Captus) vs GTF (BRAKER3) para
#      S. stellatiglandulosum
#
# PROBLEMA RAIZ:
#   El GTF de referencia contiene modelos genicos de DOS fuentes distintas
#   (columna "source": AUGUSTUS y gmst). Cada predictor numera sus genes de
#   forma independiente y reiniciando el contador (g1, g2, g3...), muchas
#   veces por contig. Esto es normal en salidas de BRAKER3/TSEBRA: se
#   conserva la fuente original de cada transcrito seleccionado.
#
#   Tu Locus_Ref_Limpio original era paste0(contig, "_", gene_id), lo cual
#   NO es unico cuando AUGUSTUS y gmst comparten "g1" en el mismo contig:
#   dos genes fisicamente distintos terminan fusionados bajo el mismo ID
#   en split()+reduce(), lo que corrompe el join contra dict_ste y en
#   cascada la clasificacion TRUE_POSITIVE/FALSE_POSITIVE y el rescate por
#   tandem.
#
# FIX: incluir la fuente (source) en el ID unico, y propagar Contig_Name
#      directamente desde seqnames() en vez de reconstruirlo con regex
#      sobre el ID compuesto (que ahora incluye la fuente).
# ============================================================================

## 2. CARGA DE GTF (REFERENCIA BRAKER3) Y GFF CAPTUS ------------------------
path_gtf_ref_ste <- paste0(ruta_base,
                           "Solanum_stellatiglandulosum_genome/braker3_Stella.gtf")

gtf_ref_completo_ste <- rtracklayer::import(path_gtf_ref_ste, format = "gtf")
exones_ref_ste <- gtf_ref_completo_ste[gtf_ref_completo_ste$type == "exon"]

# --- DIAGNOSTICO: cuantificar colisiones de gene_id entre fuentes ----------
# (corre esto ANTES del fix para confirmar el diagnostico; deberia dar > 0)
colisiones_ste <- exones_ref_ste %>%
    as.data.frame() %>%
    distinct(seqnames, source, gene_id) %>%
    count(seqnames, gene_id) %>%
    filter(n > 1)
cat(sprintf("Colisiones de gene_id entre fuentes (AUGUSTUS/gmst): %d\n",
            nrow(colisiones_ste)))

# --- CORRECCION CLAVE: la fuente entra al ID unico -------------------------
exones_ref_ste$Fuente <- as.character(exones_ref_ste$source)
exones_ref_ste$Locus_Ref_Limpio <- paste(
    as.character(seqnames(exones_ref_ste)),
    exones_ref_ste$Fuente,
    str_remove_all(exones_ref_ste$gene_id, "^gene:|\\.[a-zA-Z0-9]+$"),
    sep = "_"
)

# Agrupamos y fusionamos coordenadas espaciales por este ID ya desambiguado
genes_list <- split(exones_ref_ste, exones_ref_ste$Locus_Ref_Limpio)
genes_ref_ste <- unlist(GenomicRanges::reduce(genes_list))
genes_ref_ste$Locus_Ref_Limpio <- names(genes_ref_ste)
names(genes_ref_ste) <- NULL

cat(sprintf("Genes de referencia unicos tras el fix: %d\n",
            length(unique(genes_ref_ste$Locus_Ref_Limpio))))

# --- Carga y preparacion del GFF de CAPTUS (sin cambios) -------------------
path_gff_captus_ste <- paste0(ruta_base,
                              "Captus_Local/03_extract_stella_metabolism_complete/",
                              "S_stella__captus-ext/06_assembly_annotated/",
                              "S_stella_hit_contigs.gff")
gff_captus_completo_ste <- rtracklayer::import(path_gff_captus_ste, format = "gff3")
gff_captus_ste <- gff_captus_completo_ste[gff_captus_completo_ste$type == "protein_match:NUC"]

seqlevels(gff_captus_ste) <- str_remove_all(seqlevels(gff_captus_ste), "^[Ss]stellati[a-z]*_|^contig_")
seqlevels(genes_ref_ste)  <- str_remove_all(seqlevels(genes_ref_ste), "^[Ss]stellati[a-z]*_|^contig_")


## 3. INTERSECCION ESPACIAL Y RELACION BASE ----------------------------------
hits_globales_ste <- findOverlaps(query = gff_captus_ste, subject = genes_ref_ste,
                                  type = "any", ignore.strand = FALSE)

if (length(hits_globales_ste) == 0) {
    stop("Cero overlaps en S. stellatiglandulosum. Verifica los seqnames.")
} else {
    cat(sprintf("\nExito: Se encontraron %d solapamientos fisicos en S. stellatiglandulosum.\n",
                length(hits_globales_ste)))
}

df_cruces_ste <- data.frame(
    Captus_Row = queryHits(hits_globales_ste),
    Ref_Row    = subjectHits(hits_globales_ste)
) %>%
    mutate(
        ID_Captus = gff_captus_ste$Name[Captus_Row],
        OG_Captus = str_extract(ID_Captus, "OG[0-9]+"),
        Locus_Ref_Limpio = genes_ref_ste$Locus_Ref_Limpio[Ref_Row],
        # NUEVO: Contig_Name real, tomado directamente de seqnames del
        # GRanges de referencia -- ya no se reconstruye con regex sobre
        # el string compuesto (que ahora incluye tambien la fuente).
        Contig_Name = as.character(seqnames(genes_ref_ste))[Ref_Row]
    )


## 4. CLASIFICACION ESPACIAL PRELIMINAR Y RESCATE EN CONTIGS -----------------
df_clasificacion_base_ste <- df_cruces_ste %>%
    left_join(dict_ste, by = c("Locus_Ref_Limpio" = "Locus_Ste")) %>%
    mutate(
        ID_Captus = if_else(str_detect(ID_Captus, "__[0-9]+$"),
                            ID_Captus,
                            paste0(ID_Captus, "__00")),

        # Contig_Name ya viene propagado desde la Seccion 3 (df_cruces_ste).
        # Solo extraemos el numero de gen, que sigue al final del ID y no
        # se ve afectado por el prefijo de fuente insertado en el medio.
        Num_Gen = as.numeric(str_extract(Locus_Ref_Limpio, "(?<=_g)[0-9]+$")),

        Num_OG_Captus = str_extract(OG_Captus, "[0-9]+"),
        Num_OG_Ref = str_extract(OG_Ste, "[0-9]+"),

        Validacion_Estricta = case_when(
            !is.na(Num_OG_Ref) & Num_OG_Captus == Num_OG_Ref ~ "TRUE_POSITIVE",
            TRUE ~ "FALSE_POSITIVE"
        )
    )

# El resto del bloque de rescate por vecindad en tandem (group_by(OG_Captus),
# Tiene_Ancla_TP, Es_Tandem, etc.) NO cambia -- sigue usando Contig_Name y
# Num_Gen, que ahora llegan correctamente desambiguados por fuente.
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

cat("\nAuditoria de Rescate Espacial (S. stellatiglandulosum) -- post-fix:\n")
print(table(df_clasificacion_rescate_ste$Doble_Control_Biologico))
