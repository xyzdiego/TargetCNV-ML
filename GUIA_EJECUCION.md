# Guía de Ejecución de TargetCNV-ML

[English Execution Guide](EXECUTION_GUIDE.md) | [README en Español](README_ES.md) | [English README](README.md)

Este manual detalla la configuración del entorno, los parámetros de línea de comandos y los flujos experimentales de TargetCNV-ML desarrollados en la tesis de Maestría en Bioinformática para estimar variaciones en el número de copias génicas (CNV) a partir de datos de captura dirigida (*Target Capture* / *Hyb-Seq*).

Los mensajes de consola, encabezados `.tsv` y reportes HTML compilados con Quarto dentro de `reports/` (`reporte_entrenamiento.html`, `reporte_inferencia.html` y `reporte_transferencia.html`; empaquetados en el repositorio como `reports/reportes_html.zip`) se generan en español.

## 1. Requisitos e Instalación

### Software base
- Linux (Ubuntu 22.04/24.04, Debian 12+, Arch/CachyOS), macOS o Windows mediante WSL2.
- **R 4.6.2** y **Bioconductor 3.23** para restaurar exactamente `renv.lock` (compatible con **R $\ge 4.5$ / Bioconductor $\ge 3.22$** por el paquete `Seqinfo`).
- **Quarto CLI** (`1.10.20` o $\ge 1.4$) para compilar las plantillas `.qmd`.
- Utilidades de compilación y terminal: `bash`, `curl` o `wget`, `unzip`, `md5sum`, `gcc`, `g++`, `gfortran`, `cmake` y `make`.

### 1.1 Restauración con `renv`

```bash
git clone https://github.com/xyzdiego/TargetCNV-ML.git
cd TargetCNV-ML
chmod +x targetcnv.sh run_test_pipeline.sh scripts/*.R

./targetcnv.sh restore
./targetcnv.sh doctor
```

**Integridad de `renv.lock`:** use siempre `./targetcnv.sh restore` (o `./targetcnv.sh setup`), que llama a `renv::restore(prompt = FALSE)` sin alterar el lockfile. Evite `renv::hydrate()` o `renv::snapshot()` durante la instalación, pues enlazan librerías globales del sistema y sobrescriben `renv.lock`.

### 1.2 Entornos Pixi, Apptainer y Paquete R

```bash
# Pixi (conda-forge + bioconda + extras CRAN)
pixi install
pixi run setup-cran
pixi run doctor

# Contenedor Apptainer / Singularity (congela R 4.6.2 + Bioconductor 3.23 + Quarto 1.10.20)
apptainer build targetcnv-ml.sif Apptainer.def
apptainer run targetcnv-ml.sif doctor
```

Al incluir `DESCRIPTION`, `NAMESPACE` y páginas de ayuda en `man/`, el paquete `TargetCNVML` (60 funciones exportadas en `R/`) puede cargarse con `pkgload::load_all(".")` o instalarse con `R CMD INSTALL .`, además de ejecutarse desde los scripts CLI en `scripts/`.

## 2. Separación entre Cómputo y Reportes Quarto

Todo ajuste de modelos, validación cruzada, votación de consenso y evaluación de transferibilidad ocurre únicamente en `scripts/01_prepare_data.R` a `scripts/04_evaluate_transferability.R`. Las tres plantillas Quarto (`training_report.qmd`, `report_template.qmd` y `transferability_report.qmd`) trabajan en modo de **solo lectura**, apoyadas en `assets/report.css` y `assets/waterfall.js`. Si faltan modelos en `models/`, el reporte muestra un aviso pidiendo ejecutar `./targetcnv.sh train` sin lanzar entrenamiento interno.

## 3. Referencia de Comandos CLI

### 3.0 Envoltorio `./targetcnv.sh`

| Subcomando | Ejemplo | Descripción |
| :--- | :--- | :--- |
| `restore` / `setup` | `./targetcnv.sh restore` | Restaura la librería local desde `renv.lock`. |
| `doctor` | `./targetcnv.sh doctor` | Comprueba `Rscript`, `quarto`, `curl`, `unzip`, `md5sum` y paquetes de R. |
| `test` | `./targetcnv.sh test [--no-report]` | Corre el MWE aislado sobre `data/test_data/`. |
| `zenodo` | `./targetcnv.sh zenodo [--out-dir .]` | Descarga `TargetCNV_Models.zip`, valida MD5 y comprueba los modelos extraídos. |
| `prepare` | `./targetcnv.sh prepare [opciones de 01]` | Pasa opciones a `scripts/01_prepare_data.R`. |
| `train` | `./targetcnv.sh train [opciones de 02]` | Pasa opciones a `scripts/02_train_models.R`. |
| `infer` | `./targetcnv.sh infer -i stats.tsv --auto-zenodo` | Ejecuta `scripts/03_cnv_inference.R` (descarga modelos si `models/` está vacío). |
| `transfer` | `./targetcnv.sh transfer [--sources ...]` | Ejecuta `scripts/04_evaluate_transferability.R` y compila el reporte de transferencia. |
| `report` | `./targetcnv.sh report [--models-dir models]` | Recompila `reports/reporte_entrenamiento.html` desde artefactos existentes. |
| `reproduce-storvum` | `./targetcnv.sh reproduce-storvum --cores 4` | Corre `01` -> `02` -> `03` sobre *Solanum torvum*. |
| `run-all` | `./targetcnv.sh run-all -g ref.gff3 -c cap.gff -s stats.tsv -d dict.xlsx -i nuevos.tsv` | Corre `01` -> `02` -> `03` en una especie nueva generando reportes separados. |

**Banderas que eliminan o sobrescriben archivos previos:**
- `--clean_models` en `02_train_models.R`: borra los `.rds` y `.tsv` existentes en `--out_models` antes de escribir el nuevo entrenamiento. No está activa por defecto ni es invocada por `./targetcnv.sh reproduce-storvum`. Úsela antes de un reentrenamiento completo desde cero; evítela con `--fast_tournament` si desea conservar los 18 modelos previos.
- `--sync-root` en `./targetcnv.sh test` o `./run_test_pipeline.sh`: copia las salidas sintéticas del MWE sobre `models/` y `results/` en la raíz. Nunca la use cuando `models/` contenga modelos reales o descargados de Zenodo.

### 3.1 Módulo 01: Preparación de Datos y Etiquetado Espacial (`01_prepare_data.R`)

Asigna la etiqueta `TRUE_POSITIVE` o `FALSE_POSITIVE` a cada contig nuclear de CAPTUS (`protein_match:NUC`):
- **Escala cromosómica (`--genome_scale chromosome` / `--target_mode overlap`):** armoniza prefijos de cromosoma y cruza coordenadas y hebra (`GenomicRanges::findOverlaps`, `ignore.strand = FALSE`) frente a anotaciones `gene` en `--gff_ref`.
- **Escala de contigs (`--genome_scale contig` / `--target_mode target_attr`):** extrae el gen objetivo desde el atributo `Target` del GFF de CAPTUS y traduce grupos `HOG_ID` a ortogrupos `OG_ID` usando la tabla puente de diseño de sondas.
- **Deduplicación intra-locus:** dentro de cada muestra y gen de referencia, `resolve_redundancies()` conserva como `TRUE_POSITIVE` solo el contig de mayor longitud emparejada (`ref_len_matched`) y cobertura (`pct_recovered`), degradando fragmentos extra del mismo locus a `FALSE_POSITIVE` (`5_TP_Redundante_Degradado`).

```bash
Rscript scripts/01_prepare_data.R \
  --gff_ref data/reference/Storvum_genes.gff3 \
  --captus_gff data/reference/Storvum_hit_contigs.gff \
  --captus_stats data/reference/captus-extract_stats.tsv \
  --dict_ref data/reference/Genes_Capture_May_2026.xlsx \
  --sample Storvum \
  --genome_scale chromosome \
  --out_dir data/
```

| Corto | Largo | Tipo | Requerido | Por Defecto | Descripción |
| :--- | :--- | :--- | :---: | :--- | :--- |
| `-g` | `--gff_ref` | Archivo | Sí\* | `NULL` | GFF/GFF3 de genes de referencia (\*obligatorio en modo `chromosome`). |
| `-c` | `--captus_gff` | Archivo | **Sí** | `NULL` | GFF de contigs ensamblados por CAPTUS (`protein_match:NUC`). |
| `-s` | `--captus_stats` | Archivo | **Sí** | `NULL` | Tabla `captus-extract_stats.tsv` de CAPTUS. |
| `-d` | `--dict_ref` | Archivo | No | `NULL` | Diccionario de homología (`.xlsx`, `.tsv` o `.csv`). |
| `-p` | `--sample` | Texto | No | `NULL` | Muestra específica a filtrar (`Storvum`, `Scitrullifolium`, `S_stella`). |
| | `--genome_scale` | Texto | No | `chromosome` | `'chromosome'` (cruce GFF3) o `'contig'` (atributo `Target` + puente `OG_ID`/`HOG_ID`). |
| | `--target_mode` | Texto | No | `overlap` | Alias técnico: `'overlap'` o `'target_attr'`. |
| | `--allow_multi_sample` | Bandera | No | `FALSE` | Permite procesar tablas multimuestra sin `--sample` (solo MWE sintético). |
| `-o` | `--out_dir` | Carpeta | No | `data/` | Carpeta de salida para `data_prepared.rds` y `data_prepared.tsv`. |

### 3.2 Módulo 02: Colinealidad, Entrenamiento Determinista y SHAP (`02_train_models.R`)

Evalúa la redundancia de predictores y entrena los 24 métodos sobre `pct_recovered`, `pct_identity`, `score` y `wscore`:
1. **Colinealidad de Spearman:** exporta `models/spearman_collinearity.tsv` y `.rds` (en *S. torvum*, `score` vs. `wscore` $\rho = 0.9977$; `pct_recovered` vs. `pct_identity` $\rho = 0.7776$).
2. **Cortes de Youden (4) y `Weighted_Heuristic`:** maximiza $J = \text{Sensibilidad} + \text{Especificidad} - 1$ por variable y busca pesos enteros (`0:3`) de máxima ROC-AUC.
3. **Validación cruzada determinista (18 clasificadores ML + `Ensamble_Stacking`):** pasa la lista explícita `seeds` a `caret::trainControl()` a partir de `--seed`, manteniendo resultados idénticos con distinto número de núcleos en escenarios **Balanceado** (`down`) y **Desbalanceado**, más explicabilidad KernelSHAP y compilación de `reports/reporte_entrenamiento.html`.

```bash
Rscript scripts/02_train_models.R \
  --data data/data_prepared.rds \
  --out_models models/ \
  --cv_folds 5 \
  --sampling down \
  --cores 4 \
  --seed 1234 \
  --sources data/reference/sources_citations.tsv \
  --report_file reports/reporte_entrenamiento.html
```

| Corto | Largo | Tipo | Requerido | Por Defecto | Descripción |
| :--- | :--- | :--- | :---: | :--- | :--- |
| `-d` | `--data` | Archivo | **Sí** | `data/data_prepared.rds` | Datos etiquetados del Módulo 01. |
| `-m` | `--out_models` | Carpeta | No | `models/` | Directorio de salida para modelos `.rds`, tablas `.tsv` y SHAP. |
| `-k` | `--cv_folds` | Entero | No | `5` | Pliegues de validación cruzada ($\ge 2$). |
| `-s` | `--sampling` | Texto | No | `down` | Balanceo en el escenario balanceado (`'down'`, `'up'`, `'none'`). |
| `-c` | `--cores` | Entero | No | `detectCores() - 1` | Núcleos PSOCK paralelos. |
| | `--seed` | Entero | No | `1234` | Semilla maestra para particiones, rejillas y `seeds` por pliegue. |
| | `--extended_tournament` / `--fast_tournament` | Bandera | No | `TRUE` | Entrena los 18 modelos ML (por defecto) o el subconjunto rápido de 5 modelos. |
| | `--train_stacking` / `--no_stacking` | Bandera | No | `TRUE` | Ajusta (por defecto) u omite `Ensamble_Stacking`. |
| | `--clean_models` | Bandera | No | `FALSE` | Destructiva: elimina `.rds` y `.tsv` previos en `--out_models`. |
| | `--include_hit` | Bandera | No | `FALSE` | Añade `hit` de CAPTUS como quinto predictor (registrado en `consensus_spec.rds`; se recomienda `FALSE` entre lotes). |
| | `--render_report` / `--no_render_report` | Bandera | No | `TRUE` | Compila (por defecto) u omite `reports/reporte_entrenamiento.html`. |
| `-r` | `--report_file` | Archivo | No | `reports/reporte_entrenamiento.html` | Ruta del HTML de salida. |
| | `--sources` | Archivo | No | `""` | Tabla TSV con citas de genomas (`data/reference/sources_citations.tsv`). |

Las tablas de benchmark comparten 16 columnas (`Modelo`, `Exactitud`, `Exactitud_Balanceada`, `Sensibilidad`, `Especificidad`, `Precision`, `VPN`, `FDR`, `F1_Score`, `MCC`, `Kappa`, `ROC_AUC`, `TP`, `FP`, `FN`, `TN`). En el HTML, `*` señala los modelos con mayor índice ponderado por FDR en el escenario alternativo:

$$\text{Índice Global}_{\text{FDR}} = 0.45(1 - \text{FDR}) + 0.20\,\text{F1} + 0.15\,\text{MCC} + 0.10\,\text{Exactitud Bal.} + 0.10\,\text{ROC-AUC}$$

### 3.3 Módulo 03: Inferencia de CNV y Sensibilidad (`03_cnv_inference.R`)

Aplica el comité de consenso ($\ge 4/7$ votos, cargando dinámicamente los votantes registrados en `consensus_spec.rds`, o por defecto `ElasticNet`, `NaiveBayes`, `SVM_Polinomial`, `Red_Neuronal`, `Youden_pct_recovered`, `Youden_pct_identity` y `Weighted_Heuristic`). Las secuencias `Homologo_Valido` se ordenan por `pct_recovered` descendente en cada muestra y ortogrupo:
- `Copia_Principal`: secuencia de mayor cobertura del ortogrupo en esa muestra.
- `Homologo_Confirmado`: copia secundaria con `pct_recovered >= --cov_threshold` (`50.0%` por defecto).
- `No_Contado`: fragmento por debajo de `--cov_threshold`, retenido en el catálogo anotado pero excluido del conteo final.

Además exporta la curva de sensibilidad entre `30%` y `80%` (`cnv_threshold_sensitivity.tsv`). Dado que `captus-extract_stats.tsv` reporta la cobertura escalar total sin coordenadas por fragmento sobre la sonda, esta regla puede subestimar copias reales con recuperación menor al `50%` o sumar por separado dos fragmentos no solapantes que superen individualmente ese umbral.

```bash
Rscript scripts/03_cnv_inference.R \
  --input data/CAPTUS_Muestras/captus-extract_stats.tsv \
  --models_dir models/ \
  --out_dir results/ \
  --cov_threshold 50.0 \
  --min_votes 4 \
  --top_loci 30 \
  --sources data/reference/sources_citations.tsv
```

| Corto | Largo | Tipo | Requerido | Por Defecto | Descripción |
| :--- | :--- | :--- | :---: | :--- | :--- |
| `-i` | `--input` | Archivo | **Sí** | `NULL` | Tabla `.tsv`, `.csv` o `.rds` de extracción CAPTUS (104 muestras en el conjunto poblacional de referencia). |
| `-m` | `--models_dir` | Carpeta | No | `models/` | Directorio con modelos `.rds`. |
| `-o` | `--out_dir` | Carpeta | No | `results/` | Carpeta de salida para matrices, sensibilidad y gráficos. |
| `-c` | `--cov_threshold` | Decimal | No | `50.0` | Umbral mínimo de `pct_recovered` (%) para confirmar copias adicionales. |
| `-v` | `--min_votes` | Entero | No | `4` | Votos mínimos (sobre 7) para declarar `Homologo_Valido`. |
| `-t` | `--top_loci` | Entero | No | `30` | Número de ortogrupos más abundantes en el mapa de calor. |
| | `--true_cnv` | Archivo | No | `""` | Tabla opcional de conteos reales ($\text{CNV}^{\text{real}}_{\text{OG}}$) para evaluar concordancia. |
| | `--render_report` / `--no_render_report` | Bandera | No | `TRUE` | Compila (por defecto) u omite `reports/reporte_inferencia.html`. |
| `-r` | `--report_file` | Archivo | No | `reports/reporte_inferencia.html` | Ruta del reporte HTML. |
| | `--sources` | Archivo | No | `""` | Citas de fuentes genómicas (`data/reference/sources_citations.tsv`). |

Salidas en `<out_dir>/`: `cnv_matrix.tsv`, `cnv_annotated_sequences.tsv`, `functional_genes_catalog.tsv`, `consensus_verdicts.tsv`, `cnv_threshold_sensitivity.tsv`, `cnv_concordance_summary.tsv` y `cnv_concordance_pairs.tsv` (cuando existen etiquetas espaciales o `--true_cnv`), gráficos PDF/PNG y `reports/reporte_inferencia.html`.

### 3.4 Módulo 04: Transferibilidad Filogenética (`04_evaluate_transferability.R`)

Evalúa sin reentrenamiento los 24 métodos de *S. torvum* y el consenso de 7 modelos sobre *S. citrullifolium* (escala cromosómica) y *S. stellatiglandulosum* (escala de contigs), calcula la concordancia CNV por ortogrupo (Spearman $\rho$, Kappa cuadrático $\kappa_w$, RMSE, MAE y tasa exacta) y barre `--cov_threshold` entre `30%` y `80%`, escribiendo `results/transferability_*.tsv|.rds` y compilando `reports/reporte_transferencia.html`.

```bash
./targetcnv.sh transfer --sources data/reference/sources_citations.tsv
```

### 3.5 Recompilación Manual de Reportes (Solo Lectura)

```bash
./targetcnv.sh report

quarto render report_template.qmd \
  -P consensus_verdicts:"results/consensus_verdicts.tsv" \
  -P cnv_matrix:"results/cnv_matrix.tsv" \
  -P cnv_annotated:"results/cnv_annotated_sequences.tsv" \
  -P cnv_sensitivity:"results/cnv_threshold_sensitivity.tsv" \
  -P sources_file:"data/reference/sources_citations.tsv" \
  -P top_loci:30 \
  --output "reporte_inferencia.html" && mv -f reporte_inferencia.html reports/

./targetcnv.sh transfer --sources data/reference/sources_citations.tsv
```

## 4. Flujo 1: Batería de Pruebas y Ejemplo Mínimo Funcional (MWE)

```bash
Rscript tests/testthat.R
./targetcnv.sh test
```

`tests/testthat.R` comprueba intersección GFF con armonización de prefijos y hebra, determinismo de `build_deterministic_cv_control()` entre 1 y 2 núcleos, ajuste de `XGBoost`, invariante `ROC_AUC > 0.5` en `Ensamble_Stacking`, cálculo específico de importancia en `plot_model_varimp()`, casos borde de `extract_cnv_matrix()` y corrida completa `01` -> `02` -> `03` en directorio temporal.

## 5. Flujo 2: Modelos Preentrenados en Zenodo

Los modelos de *S. torvum* están depositados en Zenodo ([doi:10.5281/zenodo.23463987](https://doi.org/10.5281/zenodo.23463987), `TargetCNV_Models.zip`).

```bash
./targetcnv.sh zenodo
./targetcnv.sh infer -i data/CAPTUS_Muestras/captus-extract_stats.tsv --auto-zenodo
```

Para publicar una versión actualizada tras reentrenar:
1. Empaquetar y calcular el hash: `zip -r TargetCNV_Models.zip models/ && md5sum TargetCNV_Models.zip`.
2. Crear una nueva versión en Zenodo bajo el DOI `10.5281/zenodo.23463987` y subir `TargetCNV_Models.zip`.
3. Actualizar `ZENODO_URL` y `ZENODO_MD5` en la cabecera de `targetcnv.sh`.

## 6. Flujo 3: Especies Nuevas y Transferencia Inter-Especie

### 6.1 Corrida completa sobre una especie nueva (`run-all`)

```bash
./targetcnv.sh run-all \
  -g ruta/a/nuevo_genoma.gff3 \
  -c ruta/a/nuevos_contigs_captus.gff \
  -s ruta/a/captus-extract_stats.tsv \
  -d ruta/a/diccionario_homologia.xlsx \
  -p MuestraReferencia \
  -i ruta/a/captus_poblacion.tsv \
  --out_data data/nueva_especie \
  -m models/nueva_especie \
  -o results/nueva_especie \
  -r reports/reporte_entrenamiento_nueva_especie.html \
  --infer_report_file reports/reporte_inferencia_nueva_especie.html
```

### 6.2 Transferibilidad en *S. citrullifolium* y *S. stellatiglandulosum*

```bash
Rscript scripts/01_prepare_data.R \
  --gff_ref data/reference/Scitrullifolium_genes.gff3 \
  --captus_gff data/reference/Scitrullifolium_hit_contigs.gff \
  --captus_stats data/reference/captus-extract_stats.tsv \
  --dict_ref data/reference/Genes_Capture_May_2026.xlsx \
  --sample Scitrullifolium \
  --genome_scale chromosome \
  --out_dir data/citrullifolium/

Rscript scripts/01_prepare_data.R \
  --captus_gff data/reference/S_stella_hit_contigs.gff \
  --captus_stats data/reference/captus-extract_stats.tsv \
  --dict_ref data/reference/Genes_Capture_May_2026.xlsx \
  --sample S_stella \
  --genome_scale contig \
  --out_dir data/stellatiglandulosum/

./targetcnv.sh transfer --sources data/reference/sources_citations.tsv
```

#### Concordancia de copias por ortogrupo (`--cov_threshold 50.0%`, Consenso de 7 Modelos)

| Especie | Escala del Genoma | Ortogrupos ($N$) | Copias Reales | Copias Estimadas | Concordancia Exacta (%) | $\rho$ de Spearman | $\kappa_w$ Cuadrático | RMSE |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| ***S. torvum*** (Donante, en muestra) | Escala Cromosómica | 58 | 87 | 79 | **79.31** | **0.5983** ($p = 7.0 \times 10^{-7}$) | **0.5269** | 0.8710 |
| ***S. citrullifolium*** (Transferencia) | Escala Cromosómica | 55 | 56 | 66 | **70.91** | **0.4653** ($p = 3.4 \times 10^{-4}$) | **0.5506** | 0.5394 |
| ***S. stellatiglandulosum*** (Transferencia) | Escala de Contigs | 58 | 19 | 78 | 24.14 | 0.0831 ($p = 0.535$) | 0.0742 | 1.3834 |

## 7. Flujo 4: Reproducción del Experimento en *Solanum torvum*

```bash
./targetcnv.sh reproduce-storvum --cores 4
```

### Resultados representativos en validación cruzada ($K = 5$, `--seed 1234`, $N = 184$: 87 TP vs. 97 FP)

#### Escenario Balanceado (*Down-Sampling*)

| Modelo | Exactitud Bal. (%) | Sensibilidad (%) | Especificidad (%) | Precisión (%) | FDR (%) | F1 Score | MCC | Kappa | ROC-AUC |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Weighted_Heuristic** | 88.86 | 83.91 | 93.81 | 92.41 | 7.59 | **0.880** | 0.784 | 0.781 | 0.921 |
| **Ensamble_Stacking** | 88.29 | 82.76 | 93.81 | 92.31 | 7.69 | **0.873** | 0.774 | 0.770 | **0.932** |
| **Red_Neuronal (`nnet`)** | 87.31 | 83.91 | 90.72 | 89.02 | 10.98 | **0.864** | 0.750 | 0.749 | 0.927 |
| **SVM_Polinomial (`svmPoly`)** | 87.14 | 80.46 | 93.81 | 92.11 | 7.89 | **0.859** | 0.753 | 0.748 | **0.934** |
| **LDA \*** | 86.74 | 82.76 | 90.72 | 88.89 | 11.11 | **0.857** | 0.739 | 0.737 | 0.922 |
| **NaiveBayes (`nb`)** | 86.56 | 79.31 | 93.81 | 92.00 | 8.00 | **0.852** | 0.743 | 0.736 | **0.934** |
| **ElasticNet (`glmnet`) \*** | 85.25 | 83.91 | 86.60 | 84.88 | 15.12 | **0.844** | 0.706 | 0.705 | 0.930 |
| **XGBoost \*** | 84.50 | 79.31 | 89.69 | 87.34 | 12.66 | **0.831** | 0.696 | 0.693 | 0.920 |

#### Escenario Desbalanceado (Proporción original: 87 TP vs. 97 FP)

| Modelo | Exactitud Bal. (%) | Sensibilidad (%) | Especificidad (%) | Precisión (%) | FDR (%) | F1 Score | MCC | Kappa | ROC-AUC |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Weighted_Heuristic** | 88.86 | 83.91 | 93.81 | 92.41 | 7.59 | **0.880** | 0.784 | 0.781 | 0.921 |
| **Ensamble_Stacking** | 88.35 | 83.91 | 92.78 | 91.25 | 8.75 | **0.874** | 0.772 | 0.770 | **0.935** |
| **LDA** | 87.26 | 82.76 | 91.75 | 90.00 | 10.00 | **0.862** | 0.750 | 0.748 | 0.922 |
| **NaiveBayes (`nb`)** | 86.56 | 79.31 | 93.81 | 92.00 | 8.00 | **0.852** | 0.743 | 0.736 | **0.934** |
| **ElasticNet (`glmnet`)** | 85.77 | 83.91 | 87.63 | 85.88 | 14.12 | **0.849** | 0.716 | 0.716 | 0.930 |
| **Red_Neuronal (`nnet`) \*** | 84.56 | 80.46 | 88.66 | 86.42 | 13.58 | **0.833** | 0.695 | 0.694 | 0.921 |
| **XGBoost** | 84.56 | 80.46 | 88.66 | 86.42 | 13.58 | **0.833** | 0.695 | 0.694 | 0.918 |
| **SVM_Polinomial (`svmPoly`) \*** | 83.03 | 90.80 | 75.26 | 76.70 | 23.30 | **0.832** | 0.664 | 0.654 | 0.930 |

*(\* indica mejor puntaje en el índice global ponderado por FDR en el escenario alternativo; las tablas completas de los 24 métodos se exportan en `models/benchmark_models_balanced.tsv` y `models/benchmark_models_unbalanced.tsv`).*

## 8. Citación y Licencia

- **Genomas cromosómicos (*S. torvum* y *S. citrullifolium*):** Solanum Pan-Genomics (https://solpangenomics.com/dist/pages/downloads/index.php), SGN (https://solgenomics.net/ftp/genomes/Solanum_pangenomics/) y NCBI BioProject `PRJNA1073673` (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1073673); Benoit, M. *et al.* (2025), *Nature*, 640(8057): 135-145, [doi:10.1038/s41586-025-08619-6](https://doi.org/10.1038/s41586-025-08619-6).
- **Genoma a escala de contigs (*S. stellatiglandulosum*):** NCBI BioProject `PRJNA1070281` (https://www.ncbi.nlm.nih.gov/bioproject/?term=txid3132602%5BOrganism:noexp%5D); Pérez-Mesa, P. A. & Roda, F. (2025), *Current Opinion in Plant Biology*, 85: 102727, [doi:10.1016/j.pbi.2025.102727](https://doi.org/10.1016/j.pbi.2025.102727); Pérez-Mesa, P. A., Pardo, R. A., Alseekh, S. *et al.* (2026), *Plant Molecular Biology*, 116: 67, [doi:10.1007/s11103-026-01718-6](https://doi.org/10.1007/s11103-026-01718-6).

```bibtex
@software{targetcnv_ml_2026,
  author       = {Benítez Duarte, Diego Andrés},
  title        = {{TargetCNV-ML: Machine Learning and Spatial Inference for Copy Number Variations in Target Capture}},
  month        = jul,
  year         = 2026,
  publisher    = {Zenodo},
  version      = {v1.1.0},
  doi          = {10.5281/zenodo.23463987},
  url          = {https://doi.org/10.5281/zenodo.23463987}
}

@mastersthesis{benitez_duarte_thesis_2026,
  author       = {Benítez Duarte, Diego Andrés},
  title        = {{Desarrollo de un flujo de trabajo bioinformático para la Detección de Variación en el Número de Copias (CNV) de genes a partir de secuencias genómicas de lectura corta}},
  school       = {Universidad Nacional de Colombia},
  year         = 2026,
  address      = {Colombia},
  note         = {Tesis de Maestría en Bioinformática}
}
```

Distribuido bajo la Licencia MIT ([LICENSE](LICENSE)). Los conjuntos derivados y modelos preentrenados en Zenodo se distribuyen bajo CC BY 4.0.
