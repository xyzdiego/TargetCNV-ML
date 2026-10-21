# Guía Exhaustiva de Ejecución y Operación del Pipeline TargetCNV-ML

TargetCNV-ML es un software bioinformático modular y ejecutable por línea de comandos (CLI) diseñado para la detección de **Variaciones en el Número de Copias (CNV)** y la resolución de **expansiones génicas homólogas** a partir de secuencias obtenidas mediante Captura de Genes (*Target Capture*).

Esta guía detalla de punta a punta:
1. Todos los comandos CLI, sus banderas (*flags*), argumentos y valores posibles.
2. El procedimiento para reproducir el pipeline con los datos de prueba simulados (*Minimal Working Example*).
3. El procedimiento para descargar y utilizar los modelos preentrenados oficiales archivados en **Zenodo** sin necesidad de reentrenar.
4. El procedimiento para entrenar el sistema con **nuevos datos genómicos** (nuevas especies o paneles de captura).
5. La reproducción exacta del experimento de referencia con ***Solanum torvum***.
6. El registro de correcciones y mejoras de experiencia de usuario (UX/UI) aplicadas bajo el criterio de mantener una reproducibilidad del 100% con respecto al código de investigación original (`Workflow_Central_Sin_Rescate.R`).

---

## Índice

- [1. Requisitos del Sistema e Instalación](#1-requisitos-del-sistema-e-instalación)
- [2. Estructura de Directorios y Archivos](#2-estructura-de-directorios-y-archivos)
- [3. Referencia Completa de Módulos CLI](#3-referencia-completa-de-módulos-cli)
  - [3.1 Módulo 01: Preparación de Datos (`01_prepare_data.R`)](#31-módulo-01-preparación-de-datos-01_prepare_datar)
  - [3.2 Módulo 02: Entrenamiento y Benchmark (`02_train_models.R`)](#32-módulo-02-entrenamiento-y-benchmark-02_train_modelsr)
  - [3.3 Módulo 03: Inferencia de CNV y Consenso (`03_cnv_inference.R`)](#33-módulo-03-inferencia-de-cnv-y-consenso-03_cnv_inferencer)
  - [3.4 Compilación Manual del Reporte HTML Quarto](#34-compilación-manual-del-reporte-html-quarto)
- [4. Flujo 1: Ejecución del Ejemplo Mínimo Funcional (MWE / Test Data)](#4-flujo-1-ejecución-del-ejemplo-mínimo-funcional-mwe--test-data)
- [5. Flujo 2: Descarga y Uso de Modelos Preentrenados desde Zenodo](#5-flujo-2-descarga-y-uso-de-modelos-preentrenados-desde-zenodo)
- [6. Flujo 3: Entrenamiento con Nuevos Datos Genómicos](#6-flujo-3-entrenamiento-con-nuevos-datos-genómicos)
- [7. Flujo 4: Reproducción del Experimento de Referencia (*Solanum torvum*)](#7-flujo-4-reproducción-del-experimento-de-referencia-solanum-torvum)
- [8. Mejoras de UX/UI y Correcciones Aplicadas (Preservando Resultados)](#8-mejoras-de-uxui-y-correcciones-aplicadas-preservando-resultados)

---

## 1. Requisitos del Sistema e Instalación

### Requisitos de Software
- **Sistema Operativo:** Linux (Ubuntu 20.04+, Debian 11+, Arch/CachyOS, Fedora) o macOS.
- **R:** $\ge 4.0.0$ (recomendado R 4.3+ o 4.6).
- **Quarto CLI:** $\ge 1.3+$ (necesario para compilar el reporte interactivo `.html`).
- **Herramientas de sistema:** `curl` o `wget`, `unzip`, `bash`, `tar`.

### Instalación y Restauración del Entorno Aislado (`renv`)

El proyecto utiliza `renv` para congelar las versiones exactas de las librerías estadísticas, bioinformáticas y de visualización.

```bash
# 1. Clonar el repositorio
git clone https://github.com/xyzdiego/TargetCNV-ML.git
cd TargetCNV-ML

# 2. Restaurar el entorno de paquetes R (autónomo y reproducible)
Rscript -e "renv::restore()"
```

> [!NOTE]
> Las dependencias principales gestionadas por `renv` incluyen: `caret`, `glmnet`, `klaR`, `pROC`, `GenomicRanges`, `rtracklayer`, `optparse`, `kernelshap`, `shapviz`, `NeuralNetTools`, `plotly`, `DT`, `crosstalk`, `dplyr` y `tidyr`.

---

## 2. Estructura de Directorios y Archivos

```text
TargetCNV-ML/
├── .Rprofile                    # Activación automática de renv al iniciar R
├── .gitignore                   # Reglas de exclusión para datos crudos y artefactos pesados
├── renv.lock                    # Archivo de bloqueo con versiones exactas de R y paquetes
├── training_report.qmd          # Plantilla interactiva Quarto para el reporte de entrenamiento
├── run_test_pipeline.sh         # Script ejecutor del Minimal Working Example (MWE)
├── GUIA_EJECUCION.md            # Este documento técnico
├── README.md                    # Documentación general en inglés
├── README_es.md                 # Documentación general en español
│
├── scripts/                     # Scripts ejecutables desde la terminal (CLI)
│   ├── 01_prepare_data.R        # Paso 1: Intersección espacial y anotación Ground Truth
│   ├── 02_train_models.R        # Paso 2: Torneo ML, optimización Youden y reporte
│   └── 03_cnv_inference.R       # Paso 3: Consenso de 7 modelos y matriz de CNV
│
├── R/                           # Librerías internas modulares y funciones puras
│   ├── data_processing.R        # Procesamiento GFF, diccionarios y conteo de CNV
│   ├── model_evaluation.R       # Métricas diagnósticas, cortes Youden, regla de consenso
│   └── plotting.R               # Gráficos ggplot2/plotly (densidades, matrices, benchmarks)
│
├── data/
│   ├── reference/               # Archivos genómicos de referencia de S. torvum
│   ├── CAPTUS_Muestras/         # Estadísticas de CAPTUS para muestras biológicas
│   ├── test_data/               # Datos simulados livianos para pruebas (< 25 KB)
│   ├── data_prepared.rds        # Salida procesada lista para entrenamiento (S. torvum)
│   └── data_prepared.tsv        # Versión tabular legible de los datos procesados
│
├── models/                      # Modelos ML serializados (.rds) y tablas de benchmark (.tsv)
└── results/                     # Matrices de CNV, tablas de consenso y gráficos (.pdf/.png)
```

---

## 3. Referencia Completa de Módulos CLI

### 3.1 Módulo 01: Preparación de Datos (`01_prepare_data.R`)

Realiza la intersección genómica espacial de coordenadas mediante `GenomicRanges`, contrasta las secuencias contra un diccionario de control biológico u ortología, resuelve redundancias técnicas de ensamblaje degradando hits secundarios a falso positivo y calcula los diferenciales locales de calidad (*deltas* intra-locus).

#### Comando de Ejemplo
```bash
./scripts/01_prepare_data.R \
  --gff_ref data/reference/Storvum_genes.gff3 \
  --captus_gff data/reference/Storvum_hit_contigs.gff \
  --captus_stats data/reference/captus-extract_stats.tsv \
  --dict_ref data/reference/Genes_Capture_May_2026.xlsx \
  --sample Storvum \
  --out_dir data/
```

#### Catálogo Completo de Banderas (Flags)

| Bandera Corta | Bandera Larga | Tipo | Requerido | Valor por Defecto | Descripción y Opciones Válidas |
| :--- | :--- | :--- | :---: | :--- | :--- |
| `-g` | `--gff_ref` | Archivo | **SÍ** | *Ninguno* | Ruta al archivo GFF/GFF3 con la anotación de genes del genoma de referencia (debe contener registros de tipo `gene`). |
| `-c` | `--captus_gff` | Archivo | **SÍ** | *Ninguno* | Ruta al archivo GFF de cóntigos ensamblados por CAPTUS (características `protein_match:NUC`). |
| `-s` | `--captus_stats` | Archivo | **SÍ** | *Ninguno* | Archivo tabular `captus-extract_stats.tsv` con las métricas de ensamblaje generadas por CAPTUS. |
| `-d` | `--dict_ref` | Archivo | No | `NULL` | Diccionario de correspondencia ortológica (`.xlsx`, `.tsv` o `.csv`). Si se omite, los hits se marcan heurísticamente como `TRUE_POSITIVE`, lo cual es útil para exploración pero **no** para entrenamiento supervisado. |
| `-p` | `--sample` | Cadena | No | `NULL` | Nombre exacto de la muestra a extraer dentro de `captus_stats` (ej. `Storvum`). Si se especifica, el script valida que la muestra exista; si se omite, procesa todas las filas. |
| `-o` | `--out_dir` | Carpeta | No | `data/` | Directorio de salida donde se escribirán los datos preparados. |
| `-h` | `--help` | Bandera | No | - | Muestra el menú de ayuda con la descripción de los argumentos. |

#### Salidas Generadas
- `<out_dir>/data_prepared.rds`: Data frame serializado en binario nativo de R con todas las variables calculadas (`Validacion_Espacial`, `pct_recovered`, `pct_identity`, `score`, `wscore`, `delta_wscore`, etc.).
- `<out_dir>/data_prepared.tsv`: Versión en texto plano separado por tabulaciones.

---

### 3.2 Módulo 02: Entrenamiento y Benchmark (`02_train_models.R`)

Calcula los umbrales univariados óptimos mediante el **Índice J de Youden**, optimiza por búsqueda exhaustiva un **ensamble heurístico ponderado**, y entrena mediante `caret` un torneo de 5 algoritmos supervisados (`ElasticNet`, `NaiveBayes`, `SVM Polinomial`, `Red Neuronal Artificial`, `LDA`) bajo dos escenarios en paralelo:
1. **Escenario Balanceado:** Aplica submuestreo de la clase mayoritaria (*down-sampling*) en cada partición de validación cruzada.
2. **Escenario Desbalanceado:** Entrena sobre la distribución natural de clases sin remuestreo.

Adicionalmente, computa la interpretabilidad direccional mediante **KernelSHAP** sobre el modelo Elastic Net y compila automáticamente el reporte interactivo Quarto `reporte_entrenamiento.html`.

#### Comando de Ejemplo
```bash
./scripts/02_train_models.R \
  --data data/data_prepared.rds \
  --out_models models/ \
  --cv_folds 5 \
  --sampling down \
  --cores 4 \
  --seed 1234 \
  --report_file reporte_entrenamiento.html
```

#### Catálogo Completo de Banderas (Flags)

| Bandera Corta | Bandera Larga | Tipo | Requerido | Valor por Defecto | Descripción y Opciones Válidas |
| :--- | :--- | :--- | :---: | :--- | :--- |
| `-d` | `--data` | Archivo | **SÍ** | `data/data_prepared.rds` | Ruta al archivo de datos preparados (`.rds` o `.tsv`). Valida que existan las columnas `Validacion_Espacial` y las métricas predictoras. |
| `-m` | `--out_models` | Carpeta | No | `models/` | Directorio donde se almacenarán los modelos `.rds` y las tablas de métricas `.tsv`. |
| `-k` | `--cv_folds` | Entero | No | `5` | Número de particiones (*folds*) para validación cruzada k-fold ($\ge 2$). El script verifica que la clase minoritaria cuente con suficientes observaciones. |
| `-s` | `--sampling` | Cadena | No | `down` | Estrategia de balanceo para el escenario balanceado: `'down'` (submuestreo), `'up'` (sobremuestreo) o `'none'`. |
| `-c` | `--cores` | Entero | No | `detectCores() - 1` | Número de núcleos de CPU para paralelizar el entrenamiento mediante `doParallel`. Al finalizar, el cluster se libera automáticamente. |
| | `--seed` | Entero | No | `1234` | Semilla del generador pseudoaleatorio para garantizar reproducibilidad matemática exacta en la partición y el ajuste estocástico. |
| | `--include_hit` | Bandera | No | `FALSE` | Si se activa, incorpora la variable ordinal entera `hit` como predictor continuo. *(Nota: desactivado por defecto para garantizar compatibilidad con el módulo de inferencia).* |
| | `--render_report` | Bandera | No | `TRUE` | Compila automáticamente el reporte HTML interactivo tras culminar el entrenamiento. |
| | `--no_render_report` | Bandera | No | - | Desactiva la compilación del reporte Quarto (útil para ejecuciones en entornos sin Quarto o servidores headless). |
| `-r` | `--report_file` | Archivo | No | `reporte_entrenamiento.html` | Ruta completa y nombre del archivo HTML que compilará Quarto. |
| `-h` | `--help` | Bandera | No | - | Muestra la ayuda en consola. |

#### Esquema Unificado de Métricas Generadas

Tanto el escenario balanceado como el desbalanceado, los cortes univariados de Youden y el ensamble heurístico comparten el mismo conjunto de columnas estandarizadas:

1. `Modelo`: Nombre del algoritmo o estrategia evaluada.
2. `Exactitud` (*Accuracy*): $\frac{TP + TN}{TP + FP + FN + TN}$.
3. `Exactitud_Balanceada` (*Balanced Accuracy*): $\frac{\text{Sensibilidad} + \text{Especificidad}}{2}$.
4. `Sensibilidad` (*Recall / TPR*): $\frac{TP}{TP + FN}$.
5. `Especificidad` (*TNR*): $\frac{TN}{TN + FP}$.
6. `Precision` (*VPP*): $\frac{TP}{TP + FP}$.
7. `VPN` (*Valor Predictivo Negativo*): $\frac{TN}{TN + FN}$.
8. `FDR` (*False Discovery Rate*): $1 - \text{Precisión} = \frac{FP}{TP + FP}$.
9. `F1_Score`: $2 \times \frac{\text{Precisión} \times \text{Sensibilidad}}{\text{Precisión} + \text{Sensibilidad}}$.
10. `MCC` (*Coeficiente de Correlación de Matthews*): Medida de calidad robusta frente a desbalances.
11. `Kappa`: Coeficiente Kappa de Cohen sobre concordancia observada vs esperada.
12. `ROC_AUC`: Área bajo la curva ROC calculada a partir de las probabilidades o puntuaciones continuas.
13. `TP`, `FP`, `FN`, `TN`: Conteos enteros de la matriz de confusión.

#### Salidas Generadas en `<out_models>/`
- `all_trained_models.rds` / `models_balanced.rds`: Lista con los 5 modelos entrenados bajo el escenario balanceado.
- `models_unbalanced.rds`: Lista con los modelos bajo el escenario desbalanceado.
- `model_elasticnet.rds`, `model_naivebayes.rds`, `model_svm_polinomial.rds`, `model_red_neuronal.rds`, `model_lda.rds`: Modelos individuales balanceados.
- `youden_thresholds.rds`: Vector con los puntos de corte univariados óptimos para cada métrica.
- `heuristic_model.rds`: Pesos óptimos y umbral global del ensamble ponderado.
- `confusion_matrices.rds`: Tablas 2x2 de confusión para todos los modelos y cortes.
- `shap_values.rds`: Valores de explicabilidad calculados con KernelSHAP.
- `benchmark_models_balanced.tsv` y `.rds`: Tabla comparativa con todas las métricas en el escenario balanceado.
- `benchmark_models_unbalanced.tsv` y `.rds`: Tabla comparativa con todas las métricas en el escenario desbalanceado.
- `benchmark_models_all.tsv` y `.rds`: Tabla consolidada con la columna `Escenario`.
- `youden_benchmark.tsv` y `.rds`: Rendimiento de cada corte univariado junto con su `Umbral` e `Indice_Youden_J`.

---

### 3.3 Módulo 03: Inferencia de CNV y Consenso (`03_cnv_inference.R`)

Ejecuta la inferencia diagnóstica en producción sobre muestras nuevas secuenciadas por captura de genes. Integra el **Sistema de Consenso de 7 Modelos**:
- **4 Votos de Machine Learning:** Los 4 modelos supervisados con mejor rendimiento (`Red_Neuronal`, `SVM_Polinomial`, `ElasticNet`, `NaiveBayes`).
- **2 Votos Univariados de Youden:** Cortes óptimos en `pct_recovered` y `pct_identity`.
- **1 Voto del Ensamble Heurístico Ponderado:** Combinación lineal ponderada de las 4 métricas base.

**Regla de Consenso:** Mayoría simple ($\ge 4$ votos de 7) clasifica la secuencia como `Homologo_Valido`.

Posteriormente, aplica las **Reglas Biológicas de Control de Fragmentación**:
- Si la secuencia es la primera del locus: `Copia_Principal`.
- Si es una copia adicional y su `pct_recovered` $\ge \text{cov\_threshold}$ (50%): `Paralogo_Confirmado`.
- Si su cobertura es $< 50\%$ pero la suma acumulada de cobertura en el locus $\le \text{frag\_threshold}$ (120%): `Fragmentacion_de_Ensamblaje` (artefacto de ensamblador, no se suma como duplicación biológica).
- En otro caso: `Paralogo_Dudoso`.

Finalmente, extrae la **Matriz de CNV** (`Locus` $\times$ `Muestras`) y cataloga las copias funcionales.

#### Comando de Ejemplo
```bash
./scripts/03_cnv_inference.R \
  --input data/CAPTUS_Muestras/captus-extract_stats.tsv \
  --models_dir models/ \
  --out_dir results/ \
  --cov_threshold 50.0 \
  --frag_threshold 120.0 \
  --min_votes 4 \
  --top_loci 30
```

#### Catálogo Completo de Banderas (Flags)

| Bandera Corta | Bandera Larga | Tipo | Requerido | Valor por Defecto | Descripción y Opciones Válidas |
| :--- | :--- | :--- | :---: | :--- | :--- |
| `-i` | `--input` | Archivo | **SÍ** | *Ninguno* | Archivo TSV/CSV con métricas de CAPTUS para las muestras que se desean diagnosticar. Valida la presencia de columnas clave. |
| `-m` | `--models_dir` | Carpeta | No | `models/` | Carpeta que contiene los modelos serializados (`youden_thresholds.rds`, `heuristic_model.rds` y los modelos supervisados). |
| `-o` | `--out_dir` | Carpeta | No | `results/` | Directorio de destino para tablas diagnósticas, matrices y figuras. |
| `-c` | `--cov_threshold` | Decimal | No | `50.0` | Porcentaje mínimo individual de cobertura recuperada (`pct_recovered`) requerido para validar una copia homóloga biológica genuina. |
| `-f` | `--frag_threshold` | Decimal | No | `120.0` | Porcentaje acumulado de cobertura en el locus que delimita un artefacto técnico de fragmentación. |
| `-v` | `--min_votes` | Entero | No | `4` | Número mínimo de votos positivos para validar un homólogo (por defecto mayoría simple 4 de 7). Rango sugerido: 4 a 6. |
| `-t` | `--top_loci` | Entero | No | `30` | Número de familias génicas u ortogrupos más abundantes a graficar en el mapa de calor (*heatmap*). |
| `-h` | `--help` | Bandera | No | - | Muestra la ayuda en consola. |

#### Salidas Generadas en `<out_dir>/`
- `cnv_matrix.tsv`: Matriz de variación en el número de copias (filas: ortogrupos/loci, columnas: muestras, valores: número de copias funcionales detectadas).
- `cnv_annotated_sequences.tsv`: Tabla con cada secuencia candidata, sus predicciones individuales en los 7 modelos, votos obtenidos, veredicto de consenso y clasificación funcional de la copia.
- `functional_genes_catalog.tsv`: Subconjunto filtrado exclusivamente con las copias funcionales confirmadas (`Copia_Principal` y `Paralogo_Confirmado`).
- `consensus_verdicts.tsv`: Resumen de veredictos por secuencia.
- `plot_consensus_distribution.pdf` y `.png`: Gráfico de barras con la distribución de veredictos de consenso (`Homologo_Valido` vs `Ruido_Descartado`).
- `plot_cnv_heatmap.pdf` y `.png`: Mapa de calor de abundancia de CNV para los principales loci.

---

### 3.4 Compilación Manual del Reporte HTML Quarto

El script `02_train_models.R` compila este reporte de forma automática al finalizar el entrenamiento. Si desea actualizar o recompilar el reporte `reporte_entrenamiento.html` tras modificar parámetros visuales o texto sin volver a entrenar los modelos:

```bash
quarto render training_report.qmd \
  -P models_dir:"models" \
  -P data_file:"data/data_prepared.rds" \
  --output "reporte_entrenamiento.html"
```

El reporte resultante es un archivo HTML autocontenido (con CSS, JavaScript de DataTables, Plotly y mathjax embebidos) que puede abrirse directamente en cualquier navegador web sin necesidad de conexión a internet ni de un servidor Shiny.

---

## 4. Flujo 1: Ejecución del Ejemplo Mínimo Funcional (MWE / Test Data)

El repositorio incluye un conjunto de datos sintético liviano en `data/test_data/` para validar la instalación en segundos.

### Ejecución Básica (Aislada, sin tocar modelos de producción)
```bash
chmod +x run_test_pipeline.sh
./run_test_pipeline.sh
```

El script ejecutará de forma secuencial:
1. `01_prepare_data.R` sobre `data/test_data/mock_*.gff3`, guardando en `data/test_data/prepared/`.
2. `02_train_models.R` sobre los datos simulados con validación cruzada de 3 pliegues, guardando en `data/test_data/models/` y generando `data/test_data/reporte_entrenamiento_test.html`.
3. `03_cnv_inference.R` sobre `data/test_data/mock_captus_extract.tsv`, guardando en `data/test_data/results/`.

### Opción `--sync-root` (Sincronización hacia la raíz)
Si desea que las salidas de prueba se copien a las carpetas raíz `models/` y `results/`:
```bash
./run_test_pipeline.sh --sync-root
```
> [!CAUTION]
> El flag `--sync-root` sobrescribirá los modelos y tablas de producción presentes en `models/` con los modelos entrenados a partir de los datos simulados. Úselo únicamente para pruebas de integración continua.

---

## 5. Flujo 2: Descarga y Uso de Modelos Preentrenados desde Zenodo

Siguiendo las buenas prácticas de ciencia abierta, los modelos serializados entrenados con el genoma de referencia de *Solanum torvum* se encuentran archivados y respaldados con un Identificador de Objeto Digital (DOI) en Zenodo.

- **DOI Oficial:** [`10.5281/zenodo.23463987`](https://doi.org/10.5281/zenodo.23463987)
- **Registro Zenodo:** [https://zenodo.org/records/23463987](https://zenodo.org/records/23463987)
- **Archivo:** `TargetCNV_Models.zip` (Tamaño: ~382 KB)
- **Hash MD5 Oficial:** `13ec7b35711158079f962a96314eb671`

Si no dispone de datos de referencia con anotación GFF para entrenar nuevos modelos y desea realizar inferencia de CNV directamente sobre muestras biológicas nuevas secuenciadas con el panel de captura de Solanáceas:

### Paso 1: Descargar el paquete de modelos
```bash
# Descarga directa vía curl
curl -L -o TargetCNV_Models.zip https://zenodo.org/api/records/23463987/files/TargetCNV_Models.zip/content

# Alternativamente vía wget
# wget -O TargetCNV_Models.zip https://zenodo.org/api/records/23463987/files/TargetCNV_Models.zip/content
```

### Paso 2: Verificar la integridad del archivo (Checksum MD5)
```bash
md5sum TargetCNV_Models.zip
# La salida debe coincidir exactamente con:
# 13ec7b35711158079f962a96314eb671  TargetCNV_Models.zip
```

### Paso 3: Descomprimir en la raíz del proyecto
El archivo zip ya contiene el directorio `models/` con todos los clasificadores preentrenados y umbrales:
```bash
unzip -o TargetCNV_Models.zip
```

Verifique los modelos extraídos:
```bash
ls -la models/
# Debe listar:
# all_trained_models.rds
# benchmark_models.rds
# confusion_matrices.rds
# heuristic_model.rds
# model_elasticnet.rds
# model_lda.rds
# model_naivebayes.rds
# model_red_neuronal.rds
# model_svm_polinomial.rds
# youden_benchmark.rds
# youden_thresholds.rds
```

### Paso 4: Ejecutar inferencia diagnóstica directa (Paso 3)
```bash
./scripts/03_cnv_inference.R \
  --input data/CAPTUS_Muestras/captus-extract_stats.tsv \
  --models_dir models \
  --out_dir results \
  --cov_threshold 50.0 \
  --frag_threshold 120.0 \
  --min_votes 4 \
  --top_loci 30
```

---

## 6. Flujo 3: Entrenamiento con Nuevos Datos Genómicos

Si desea calibrar y entrenar TargetCNV-ML para una **nueva especie**, un **nuevo genoma de referencia** o un **nuevo panel de captura génica**:

### 6.1 Requisitos de los Archivos de Entrada

1. **Anotación GFF3 del Genoma de Referencia (`--gff_ref`):**
   - Archivo GFF/GFF3 estándar con anotación de coordenadas génicas (registros con la 3ª columna igual a `gene`).
   - El atributo `ID` en la 9ª columna debe contener el identificador del locus (ej. `ID=gene:Solyc01g005000` o `ID=Soltor1_01g000010.1`).
2. **Cóntigos Ensamblados por CAPTUS en GFF (`--captus_gff`):**
   - Archivo GFF generado por el módulo de alineamiento o extracción de CAPTUS (con registros `protein_match:NUC`).
   - El atributo `Name` debe registrar la relación ortológica del cóntigo (ej. `OG0000254__00`).
3. **Tabla de Estadísticas de CAPTUS (`--captus_stats`):**
   - El archivo `captus-extract_stats.tsv` generado por CAPTUS.
   - Columnas mínimas requeridas: `sample_name`, `marker_type`, `locus`, `ref_name`, `hit`, `pct_recovered`, `pct_identity`, `score`, `wscore`.
4. **Diccionario de Homología Teórica / Control Biológico (`--dict_ref`):**
   - Para que el entrenamiento sea posible, el conjunto de datos debe contener tanto **Verdaderos Positivos (TP)** como **Falsos Positivos o Ruido (FP)**.
   - El diccionario puede ser un archivo Excel (`.xlsx`), TSV o CSV que contenga:
     - Una columna que identifique el gen de referencia (ej. `Locus_Benoit`, `Locus`, `gene_id` o `Locus_Ref_Limpio`).
     - Una columna que identifique el ortogrupo teórico al que pertenece (ej. `OG_Benoit`, `OG_ID`, `Orthogroup_ID` u `OG`).
   - Un hit de CAPTUS cuya intersección espacial en el genoma coincida con el ortogrupo teórico se clasifica como `TRUE_POSITIVE`. Si cae en una coordenada de otro ortogrupo o fuera de objetivo, se cataloga como `FALSE_POSITIVE`.

### 6.2 Paso 1: Preparar los Nuevos Datos
```bash
./scripts/01_prepare_data.R \
  --gff_ref ruta/a/nuevo_genoma.gff3 \
  --captus_gff ruta/a/nuevos_contigs_captus.gff \
  --captus_stats ruta/a/captus-extract_stats.tsv \
  --dict_ref ruta/a/diccionario_homologia.xlsx \
  --sample NombreMuestraReferencia \
  --out_dir data/nueva_especie/
```

### 6.3 Paso 2: Entrenar el Torneo de Modelos y Generar el Reporte HTML
```bash
./scripts/02_train_models.R \
  --data data/nueva_especie/data_prepared.rds \
  --out_models models/nueva_especie/ \
  --cv_folds 5 \
  --sampling down \
  --cores 4 \
  --seed 1234 \
  --report_file reporte_entrenamiento_nueva_especie.html
```

El script:
- Verificará que el número de TP y FP sea suficiente para los folds solicitados.
- Optimizará los puntos de corte univariados de Youden para la nueva especie.
- Calibrará los pesos del ensamble heurístico.
- Ajustará los hiperparámetros de los 5 algoritmos de ML mediante validación cruzada.
- Generará las explicaciones KernelSHAP.
- Compilará el reporte HTML `reporte_entrenamiento_nueva_especie.html`.

### 6.4 Paso 3: Aplicar Inferencia de CNV con los Nuevos Modelos
```bash
./scripts/03_cnv_inference.R \
  --input ruta/a/captus_muestras_estudio.tsv \
  --models_dir models/nueva_especie/ \
  --out_dir results/nueva_especie/ \
  --cov_threshold 50.0 \
  --frag_threshold 120.0 \
  --min_votes 4 \
  --top_loci 30
```

---

## 7. Flujo 4: Reproducción del Experimento de Referencia (*Solanum torvum*)

A continuación se muestran los comandos exactos para reproducir los resultados de tesis sobre los datos biológicos reales de *Solanum torvum*:

### Paso 1: Extracción de Características e Intersección Espacial
```bash
./scripts/01_prepare_data.R \
  --gff_ref data/reference/Storvum_genes.gff3 \
  --captus_gff data/reference/Storvum_hit_contigs.gff \
  --captus_stats data/reference/captus-extract_stats.tsv \
  --dict_ref data/reference/Genes_Capture_May_2026.xlsx \
  --sample Storvum \
  --out_dir data/
```
*Salida:* `data/data_prepared.rds` conteniendo 184 secuencias válidas NUC (87 Verdaderos Positivos y 97 Falsos Positivos/Ruido).

### Paso 2: Calibración, Torneo de Machine Learning y Reporte HTML
```bash
./scripts/02_train_models.R \
  --data data/data_prepared.rds \
  --out_models models \
  --cv_folds 5 \
  --sampling down \
  --cores 4 \
  --seed 1234 \
  --report_file reporte_entrenamiento.html
```

#### Rendimiento Obtenido en Validación Cruzada (K=5, Semilla 1234)

##### Escenario Balanceado (*Down-Sampling*)
| Modelo | Exactitud Bal. (%) | Sensibilidad (%) | Especificidad (%) | Precisión (%) | FDR (%) | F1 Score | MCC | Kappa | ROC-AUC |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Modelo Heurístico** | 88.86 | 83.91 | 93.81 | 92.41 | 7.59 | **0.880** | 0.784 | 0.781 | 0.921 |
| **Red Neuronal (`nnet`)** | 87.32 | 83.91 | 90.72 | 88.98 | 11.02 | **0.864** | 0.750 | 0.749 | 0.927 |
| **SVM Polinomial (`svmPoly`)** | 87.14 | 80.46 | 93.81 | 92.11 | 7.89 | **0.859** | 0.753 | 0.748 | **0.934** |
| **LDA** | 86.74 | 82.76 | 90.72 | 88.89 | 11.11 | **0.857** | 0.739 | 0.737 | 0.922 |
| **Naive Bayes (`nb`)** | 86.56 | 79.31 | 93.81 | 92.00 | 8.00 | **0.852** | 0.743 | 0.736 | **0.934** |
| **Elastic Net (`glmnet`)** | 85.28 | 83.91 | 86.60 | 84.88 | 15.12 | **0.844** | 0.706 | 0.705 | 0.930 |

##### Escenario Desbalanceado (Distribución Natural 87 TP vs 97 FP)
| Modelo | Exactitud Bal. (%) | Sensibilidad (%) | Especificidad (%) | Precisión (%) | FDR (%) | F1 Score | MCC | Kappa | ROC-AUC |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Modelo Heurístico** | 88.86 | 83.91 | 93.81 | 92.41 | 7.59 | **0.880** | 0.784 | 0.781 | 0.921 |
| **LDA** | 87.26 | 82.76 | 91.75 | 90.00 | 10.00 | **0.862** | 0.750 | 0.748 | 0.922 |
| **Naive Bayes (`nb`)** | 86.56 | 79.31 | 93.81 | 92.00 | 8.00 | **0.852** | 0.743 | 0.736 | **0.934** |
| **Elastic Net (`glmnet`)** | 85.76 | 83.91 | 87.63 | 85.88 | 14.12 | **0.849** | 0.716 | 0.716 | 0.930 |
| **Red Neuronal (`nnet`)** | 84.60 | 80.46 | 88.66 | 86.42 | 13.58 | **0.833** | 0.695 | 0.694 | 0.921 |
| **SVM Polinomial (`svmPoly`)** | 83.04 | 90.80 | 75.26 | 76.70 | 23.30 | **0.832** | 0.664 | 0.654 | 0.930 |

### Paso 3: Inferencia Multimuestra y Extracción de CNV
```bash
./scripts/03_cnv_inference.R \
  --input data/CAPTUS_Muestras/captus-extract_stats.tsv \
  --models_dir models \
  --out_dir results \
  --cov_threshold 50.0 \
  --frag_threshold 120.0 \
  --min_votes 4 \
  --top_loci 30
```

---

## 8. Mejoras de UX/UI y Correcciones Aplicadas (Preservando Resultados)

En consonancia con el requerimiento de **no alterar los resultados numéricos ni la lógica biológica de `Workflow_Central_Sin_Rescate.R`**, se realizaron intervenciones estrictamente en la arquitectura de software, validaciones preventivas, experiencia de usuario y estandarización visual:

### 1. Estandarización Universal del Formato de Salida de Benchmarks
- **Problema previo:** El archivo `benchmark_models_balanced.tsv` y `benchmark_models_unbalanced.tsv` presentaban discrepancias estructurales: uno omitía el modelo heurístico, faltaban métricas diagnósticas (VPN, MCC, Kappa, Exactitud Balanceada) y el ROC-AUC del ensamble figuraba como `NA`.
- **Solución implementada:** Se diseñó la función canónica `compute_classification_metrics()` en `R/model_evaluation.R`. Ahora, todas las tablas de benchmark (`balanceado`, `desbalanceado`, `youden_benchmark` y el consolidado `benchmark_models_all.tsv`) contienen las mismas 16 columnas en idéntico orden, permitiendo comparaciones tabulares directas y análisis reproducibles.

### 2. Aislamiento del Ejemplo Mínimo Funcional (`run_test_pipeline.sh`)
- **Problema previo:** La ejecución del script de pruebas `run_test_pipeline.sh` copiaba indiscriminadamente los modelos y resultados simulados de juguete a `models/` y `results/`, sobrescribiendo los modelos reales entrenados con *S. torvum*.
- **Solución implementada:** Se encapsuló la suite de pruebas para que opere exclusivamente dentro de `data/test_data/prepared/`, `data/test_data/models/` y `data/test_data/results/`. La copia hacia las carpetas de la raíz ahora es una opción explícita y controlada mediante la bandera `--sync-root`.

### 3. Validación Temprana de Argumentos y Mensajes de Diagnóstico
- **En `01_prepare_data.R`:** Si el usuario especifica una muestra con `--sample` que no existe en el archivo de estadísticas de CAPTUS, el script se detiene de inmediato informando el error y listando las muestras disponibles en el archivo, evitando generar silenciosamente una tabla vacía.
- **En `01_prepare_data.R`:** Si se especifica un archivo en `--dict_ref` pero la ruta es incorrecta, el script falla explícitamente en lugar de continuar silenciosamente marcando todos los datos como verdaderos positivos.
- **En `02_train_models.R`:** Se añadió la opción `--report_file` para definir el destino del reporte interactivo sin pisar reportes previos, se agregaron validaciones de parámetros (`--sampling`, `--cv_folds`, `--cores`) y se valida que la clase minoritaria cuente con más muestras que el número de folds antes de iniciar el cómputo intensivo.
- **En `03_cnv_inference.R`:** Se verifica que el archivo de entrada contenga todas las columnas obligatorias de CAPTUS antes de cargar los modelos en memoria.

### 4. Experiencia Visual Interactiva en Quarto (`training_report.qmd`)
- Se implementó un formateador único (`format_benchmark_dt`) que estandariza las tablas dinámicas de DataTables en Quarto:
  - Formato porcentual a 2 decimales para Exactitud, Sensibilidad, Especificidad, Precisión, VPN y FDR.
  - Formato decimal a 3 decimales para F1-Score, MCC, Kappa y ROC-AUC.
  - Barras horizontales de gradiente de color para una lectura visual intuitiva del `F1_Score` (azul claro) y `ROC_AUC` (verde claro).
  - Manejo elegante de valores no definidos mediante guion (`–`).
- Se dinamizaron los textos del resumen metodológico para que lean los conteos reales del dataset en memoria en lugar de valores estáticos, manteniendo plena validez si se reentrena con otras especies.

### 5. Corrección de Referencias a Modelos en Zenodo
- Se corrigió en la documentación el nombre y formato del archivo empaquetado en Zenodo: el artefacto oficial es `TargetCNV_Models.zip` (que se descomprime directamente con `unzip` conservando la carpeta `models/`), en lugar de la referencia errónea a un archivo `.tar.gz`.
