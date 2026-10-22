# Guía de Ejecución de TargetCNV-ML

**Idiomas / Languages:** [Español (GUIA\_EJECUCION.md)](https://www.google.com/search?q=GUIA_EJECUCION.md) | [English (EXECUTION\_GUIDE.md)](https://www.google.com/search?q=EXECUTION_GUIDE.md)  
**Resumen General:** [README en Español (README\_es.md)](https://www.google.com/search?q=README_es.md) | [English README (README.md)](README.md)

**TargetCNV-ML** es una herramienta bioinformática de línea de comandos (CLI) que detecta **Variaciones en el Número de Copias (CNV)** e identifica **expansiones de genes homólogos** a partir de secuencias obtenidas mediante captura de genes (*Target Capture*).

En esta guía encontrarás:
1. Requisitos del sistema, instalación y verificación del entorno (`./targetcnv.sh doctor`).
2. Estructura de carpetas y archivos del repositorio.
3. Uso del script principal en una línea (`./targetcnv.sh`) y descripción de todos los parámetros de los tres módulos (`01_prepare_data.R`, `02_train_models.R`, `03_cnv_inference.R`) y del reporte HTML (`training_report.qmd`).
4. **Flujo 1:** Prueba rápida con los datos simulados incluidos en el repositorio (*Minimal Working Example* / MWE).
5. **Flujo 2:** Descarga y uso de los modelos preentrenados publicados en **Zenodo** ([`DOI: 10.5281/zenodo.23463987`](https://doi.org/10.5281/zenodo.23463987)) sin necesidad de volver a entrenar.
6. **Flujo 3:** Entrenamiento con nuevos datos genómicos y pruebas de **transferibilidad entre especies** (tanto en genomas ensamblados a nivel de cromosoma como *Solanum citrullifolium*, como en genomas fragmentados como *Solanum stellatiglandulosum*).
7. **Flujo 4:** Reproducción paso a paso del experimento de referencia con ***Solanum torvum***.
8. Resumen de mejoras técnicas y de interfaz aplicadas conservando exactamente los mismos resultados del código original (`Workflow_Central_Sin_Rescate.R`).

## 1. Requisitos del Sistema e Instalación

### Requisitos de Software
- **Sistema Operativo:** Linux (Ubuntu 20.04+, Debian 11+, Arch/CachyOS, Fedora), macOS o Windows (mediante WSL2).
- **R:** $\ge 4.0.0$ (recomendado R 4.3+ o 4.6).
- **Quarto CLI:** $\ge 1.3+$ (necesario para generar el reporte interactivo `.html`).
- **Herramientas de sistema:** `bash`, `curl` o `wget`, `unzip`, `md5sum`, `tar`.

### Instalación y Restauración del Entorno de R (`renv`)

El proyecto usa `renv` para instalar exactamente las mismas versiones de los paquetes estadísticos y bioinformáticos con los que se desarrolló el pipeline.

```bash
# 1. Clonar el repositorio
git clone https://github.com/xyzdiego/TargetCNV-ML.git
cd TargetCNV-ML

# 2. Dar permisos de ejecución a los scripts y restaurar las librerías de R
chmod +x targetcnv.sh run_test_pipeline.sh scripts/*.R
Rscript -e "renv::restore()"

# 3. Comprobar que todas las herramientas y paquetes estén listos
./targetcnv.sh doctor
```

> [!NOTE]
> Los paquetes principales gestionados por `renv` son: `caret`, `glmnet`, `klaR`, `pROC`, `GenomicRanges`, `rtracklayer`, `optparse`, `readxl`, `kernelshap`, `shapviz`, `NeuralNetTools`, `plotly`, `DT`, `cowplot`, `dplyr`, `tidyr` y `stringr`.

## 2. Estructura de Directorios y Archivos

```text
TargetCNV-ML/
├── .Rprofile                    # Activa renv automáticamente al abrir R
├── .gitignore                   # Reglas para ignorar archivos temporales o pesados
├── renv.lock                    # Registro de versiones exactas de R y paquetes
├── targetcnv.sh                 # Script principal en Bash para ejecutar flujos en 1 línea
├── run_test_pipeline.sh         # Script que ejecuta el ejemplo de prueba (MWE)
├── training_report.qmd          # Plantilla Quarto para el reporte HTML de entrenamiento
├── GUIA_EJECUCION.md            # Guía de ejecución en español (este archivo)
├── EXECUTION_GUIDE.md           # Guía de ejecución en inglés
├── README.md                    # Resumen general en inglés
├── README_es.md                 # Resumen general en español
│
├── scripts/                     # Scripts ejecutables desde la terminal (CLI)
│   ├── 01_prepare_data.R        # Paso 1: Cruce espacial GFF y etiquetado TP/FP
│   ├── 02_train_models.R        # Paso 2: Entrenamiento de modelos, cortes de Youden y reporte
│   └── 03_cnv_inference.R       # Paso 3: Consenso de 7 modelos y matriz de CNV
│
├── R/                           # Funciones modulares de R
│   ├── data_processing.R        # Lectura de GFF, diccionarios y conteo de copias (CNV)
│   ├── model_evaluation.R       # Cálculo de métricas, índice de Youden y regla de consenso
│   └── plotting.R               # Funciones para generar gráficos en ggplot2 y plotly
│
├── data/
│   ├── reference/               # Archivos genómicos de referencia de S. torvum
│   ├── CAPTUS_Muestras/         # Estadísticas de CAPTUS para muestras biológicas
│   ├── test_data/               # Datos simulados ligeros para pruebas (< 25 KB)
│   ├── data_prepared.rds        # Datos procesados listos para entrenar (S. torvum)
│   └── data_prepared.tsv        # Versión en tabla de texto de los datos procesados
│
├── models/                      # Modelos guardados (.rds) y tablas de métricas (.tsv)
└── results/                     # Matrices de CNV, tablas de consenso y gráficos (.pdf/.png)
```

## 3. Referencia de Comandos CLI

### 3.0 Ejecución Rápida con `targetcnv.sh`

El script `./targetcnv.sh` reúne las tareas más comunes para no tener que escribir rutas largas en la terminal:

| Subcomando | Ejemplo | Qué hace |
| :--- | :--- | :--- |
| `doctor` | `./targetcnv.sh doctor` | Revisa si están instalados `Rscript`, `quarto`, `curl`, `unzip`, `md5sum` y los paquetes de R. |
| `test` | `./targetcnv.sh test [--sync-root] [--no-report]` | Corre el ejemplo de prueba con datos simulados dentro de `data/test_data/`. |
| `zenodo` | `./targetcnv.sh zenodo [--out-dir .]` | Descarga `TargetCNV_Models.zip` desde Zenodo, comprueba su código MD5 y extrae la carpeta `models/`. |
| `infer` | `./targetcnv.sh infer -i stats.tsv --auto-zenodo` | Calcula la matriz de CNV en muestras nuevas; con `--auto-zenodo` descarga los modelos si `models/` está vacío. |
| `reproduce-storvum` | `./targetcnv.sh reproduce-storvum --cores 4` | Corre de principio a fin los pasos 01, 02 y 03 con los datos reales de *Solanum torvum*. |
| `run-all` | `./targetcnv.sh run-all -g ref.gff3 -c cap.gff -s stats.tsv -d dict.xlsx -i nuevos.tsv` | Ejecuta en orden los pasos `01_prepare_data.R` $\rightarrow$ `02_train_models.R` $\rightarrow$ `03_cnv_inference.R`. |
| `prepare` | `./targetcnv.sh prepare [opciones de 01]` | Envía los argumentos directamente a `scripts/01_prepare_data.R`. |
| `train` | `./targetcnv.sh train [opciones de 02]` | Envía los argumentos directamente a `scripts/02_train_models.R`. |
| `report` | `./targetcnv.sh report --models-dir models --data data/data_prepared.rds` | Vuelve a generar el reporte HTML con Quarto sin reentrenar los modelos. |

### 3.1 Módulo 01: Preparación de Datos (`01_prepare_data.R`)

Este script identifica qué cóntigos ensamblados por CAPTUS corresponden a genes reales del genoma de referencia y cuáles son ruido o copias redundantes:
1. Cruza las coordenadas genómicas de los cóntigos de CAPTUS (`protein_match:NUC`) con los genes del archivo `.gff3` de referencia usando `GenomicRanges` (o lee el gen directamente del atributo `Target` si se usa `--target_mode target_attr` en genomas muy fragmentados).
2. Compara cada cruce con el diccionario de ortología para etiquetar los aciertos como `TRUE_POSITIVE` y los alineamientos en otros genes o fuera de objetivo como `FALSE_POSITIVE`.
3. Si varios cóntigos caen sobre el mismo gen de referencia, conserva el cóntigo de mayor longitud y cobertura como `TRUE_POSITIVE` y reclasifica los restantes como `FALSE_POSITIVE` (`5_TP_Redundante_Degradado`).
4. Calcula la diferencia de calidad de cada hit con respecto al mejor hit de su mismo locus (`delta_wscore` y `delta_pct_identity`).

#### Ejemplo de Uso
```bash
./scripts/01_prepare_data.R \
  --gff_ref data/reference/Storvum_genes.gff3 \
  --captus_gff data/reference/Storvum_hit_contigs.gff \
  --captus_stats data/reference/captus-extract_stats.tsv \
  --dict_ref data/reference/Genes_Capture_May_2026.xlsx \
  --sample Storvum \
  --target_mode overlap \
  --out_dir data/
```

#### Parámetros Disponibles

| Corto | Largo | Tipo | Requerido | Valor por Defecto | Descripción |
| :--- | :--- | :--- | :---: | :--- | :--- |
| `-g` | `--gff_ref` | Archivo | **SÍ*** | *Ninguno* | Archivo GFF/GFF3 con la anotación de genes del genoma de referencia (registros `gene`). *Requerido en modo `overlap`. |
| `-c` | `--captus_gff` | Archivo | **SÍ** | *Ninguno* | Archivo GFF de cóntigos ensamblados por CAPTUS (registros `protein_match:NUC`). |
| `-s` | `--captus_stats` | Archivo | **SÍ** | *Ninguno* | Archivo `captus-extract_stats.tsv` con las métricas numéricas de alineamiento de CAPTUS. |
| `-d` | `--dict_ref` | Archivo | No | `NULL` | Tabla de correspondencia entre genes y ortogrupos (`.xlsx`, `.tsv` o `.csv`). Compatible con el archivo multi-hoja `Genes_Capture_May_2026.xlsx` o tablas simples. |
| `-p` | `--sample` | Texto | No | `NULL` | Nombre de la muestra a filtrar dentro de `captus_stats` (ej. `Storvum`, `Scitrullifolium`, `S_stella`). |
| | `--target_mode` | Texto | No | `overlap` | `'overlap'`: cruza coordenadas con `--gff_ref`. `'target_attr'`: extrae el gen de referencia desde el atributo `Target` del GFF de CAPTUS (útil en genomas fragmentados). |
| `-o` | `--out_dir` | Carpeta | No | `data/` | Carpeta donde se guardarán los archivos preparados. |
| `-h` | `--help` | Bandera | No | - | Muestra la ayuda en pantalla. |

#### Archivos de Salida
- `<out_dir>/data_prepared.rds`: Tabla en formato nativo de R con las etiquetas `Validacion_Espacial`, `Doble_Control_Biologico` y las métricas numéricas.
- `<out_dir>/data_prepared.tsv`: La misma tabla en texto separado por tabulaciones.

### 3.2 Módulo 02: Entrenamiento y Evaluación de Modelos (`02_train_models.R`)

Con los datos etiquetados del Paso 1, este script entrena y compara diferentes métodos para separar homólogos reales (`TP`) de ruido (`FP`):
1. **Cortes individuales (Youden):** Calcula para cada una de las 4 métricas (`pct_recovered`, `pct_identity`, `score`, `wscore`) el umbral que maximiza el índice $J$ de Youden ($\text{Sensibilidad} + \text{Especificidad} - 1$).
2. **Modelo heurístico ponderado:** Prueba combinaciones de pesos enteros (de 0 a 3) sobre los votos de las 4 métricas para encontrar la suma ponderada con mayor área bajo la curva ROC (AUC).
3. **Modelos de Machine Learning (`caret`):** Entrena 5 algoritmos (`ElasticNet`, `NaiveBayes`, `SVM_Polinomial`, `Red_Neuronal`, `LDA`) mediante validación cruzada en dos escenarios:
   - **Balanceado:** Iguala el número de `TP` y `FP` en cada partición mediante submuestreo (*down-sampling*).
   - **Desbalanceado:** Entrena con la proporción original de los datos.
4. **Explicabilidad (KernelSHAP) y Reporte:** Calcula la contribución de cada variable sobre las predicciones del modelo Elastic Net y genera el reporte interactivo `reporte_entrenamiento.html`.

#### Ejemplo de Uso
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

#### Parámetros Disponibles

| Corto | Largo | Tipo | Requerido | Valor por Defecto | Descripción |
| :--- | :--- | :--- | :---: | :--- | :--- |
| `-d` | `--data` | Archivo | **SÍ** | `data/data_prepared.rds` | Archivo de datos preparados (`.rds` o `.tsv`) generado en el Paso 1. |
| `-m` | `--out_models` | Carpeta | No | `models/` | Carpeta donde se guardarán los modelos `.rds` y las tablas de resultados `.tsv`. |
| `-k` | `--cv_folds` | Entero | No | `5` | Número de particiones para la validación cruzada ($\ge 2$). |
| `-s` | `--sampling` | Texto | No | `down` | Método de balanceo para el escenario balanceado: `'down'` (submuestreo), `'up'` (sobremuestreo) o `'none'`. |
| `-c` | `--cores` | Entero | No | `detectCores() - 1` | Número de núcleos de procesador para entrenar en paralelo. |
| | `--seed` | Entero | No | `1234` | Semilla aleatoria para obtener resultados exactamente reproducibles. |
| | `--extended_tournament` | Bandera | No | `FALSE` | Entrena los 18 algoritmos evaluados en el script de investigación original (añade `RandomForest`, `GBM`, `XGBoost`, `DNN`, etc.). |
| | `--train_stacking` | Bandera | No | `FALSE` | Entrena un meta-modelo GLM (`caretStack`) que combina las probabilidades de los dos mejores modelos (`model_stacking.rds`). |
| | `--include_hit` | Bandera | No | `FALSE` | Añade la columna `hit` como variable predictora. *(Desactivado por defecto para mantener compatibilidad con `03_cnv_inference.R`).* |
| | `--render_report` | Bandera | No | `TRUE` | Genera automáticamente el reporte HTML al terminar el entrenamiento. |
| | `--no_render_report` | Bandera | No | - | No genera el reporte HTML (útil en servidores sin Quarto). |
| `-r` | `--report_file` | Archivo | No | `reporte_entrenamiento.html` | Nombre y ruta del archivo HTML del reporte. |
| `-h` | `--help` | Bandera | No | - | Muestra la ayuda en pantalla. |

#### Métricas Incluidas en las Tablas de Resultados (16 Columnas)

Todas las tablas de comparación (`benchmark_models_balanced.tsv`, `benchmark_models_unbalanced.tsv`, `youden_benchmark.tsv` y `benchmark_models_all.tsv`) incluyen las mismas 16 columnas:

1. `Modelo`: Nombre del clasificador o regla de corte.
2. `Exactitud` (*Accuracy*): Proporción total de aciertos $\frac{TP + TN}{TP + FP + FN + TN}$.
3. `Exactitud_Balanceada`: Promedio entre sensibilidad y especificidad $\frac{\text{Sensibilidad} + \text{Especificidad}}{2}$.
4. `Sensibilidad` (*Recall*): Proporción de verdaderos positivos detectados $\frac{TP}{TP + FN}$.
5. `Especificidad`: Proporción de falsos positivos correctamente descartados $\frac{TN}{TN + FP}$.
6. `Precision` (*Valor Predictivo Positivo*): Proporción de predicciones positivas que son realmente verdaderas $\frac{TP}{TP + FP}$.
7. `VPN` (*Valor Predictivo Negativo*): $\frac{TN}{TN + FN}$.
8. `FDR` (*Tasa de Falsos Descubrimientos*): $1 - \text{Precisión} = \frac{FP}{TP + FP}$.
9. `F1_Score`: Media armónica entre precisión y sensibilidad $2 \times \frac{\text{Precisión} \times \text{Sensibilidad}}{\text{Precisión} + \text{Sensibilidad}}$.
10. `MCC` (*Coeficiente de Correlación de Matthews*): Medida de desempeño entre $-1$ y $+1$ que toma en cuenta las cuatro celdas de la matriz de confusión.
11. `Kappa`: Índice Kappa de Cohen de concordancia más allá del azar.
12. `ROC_AUC`: Área bajo la curva ROC.
13. `TP`, `FP`, `FN`, `TN`: Conteos de la matriz de confusión.

#### Archivos Generados en `<out_models>/`
- `all_trained_models.rds` y `models_balanced.rds`: Lista con los modelos entrenados en el escenario balanceado.
- `models_unbalanced.rds`: Lista con los modelos entrenados sin balanceo.
- `model_elasticnet.rds`, `model_naivebayes.rds`, `model_svm_polinomial.rds`, `model_red_neuronal.rds`, `model_lda.rds`: Modelos individuales listos para usar en el Paso 3.
- `youden_thresholds.rds`: Puntos de corte óptimos de cada métrica.
- `heuristic_model.rds`: Pesos y umbral del modelo heurístico ponderado.
- `confusion_matrices.rds`: Matrices de confusión 2x2 de todos los métodos.
- `shap_values.rds`: Valores KernelSHAP para los gráficos de importancia de variables.
- `benchmark_models_balanced.tsv` / `.rds`, `benchmark_models_unbalanced.tsv` / `.rds`, `benchmark_models_all.tsv` / `.rds` y `youden_benchmark.tsv` / `.rds`.

### 3.3 Módulo 03: Inferencia de CNV y Consenso (`03_cnv_inference.R`)

Este script toma la tabla de estadísticas de CAPTUS de muestras nuevas, filtra el ruido mediante votación de modelos y cuenta las copias génicas por locus y muestra:

1. **Sistema de consenso de 7 modelos:** Cada secuencia recibe 7 evaluaciones independientes:
   - 4 votos de modelos de Machine Learning (`ElasticNet`, `NaiveBayes`, `SVM_Polinomial`, `Red_Neuronal`).
   - 2 votos de cortes univariados de Youden (`pct_recovered` y `pct_identity`).
   - 1 voto del modelo heurístico ponderado.
   Si la secuencia suma **4 o más votos positivos** (`--min_votes 4`), se clasifica como `Homologo_Valido`; de lo contrario, se descarta como `Ruido_Descartado`.

2. **Reglas de control de fragmentación:** Dentro de cada muestra y locus, las secuencias aceptadas como `Homologo_Valido` se ordenan de mayor a menor cobertura (`pct_recovered`):
   - La primera secuencia (mayor cobertura) se registra como `Copia_Principal`.
   - Cualquier secuencia adicional con cobertura individual $\ge 50\%$ (`--cov_threshold 50.0`) se confirma como duplicación real (`Paralogo_Confirmado`).
   - Si una secuencia adicional tiene cobertura $< 50\%$ y la cobertura acumulada en el locus es $\le 120\%$ (`--frag_threshold 120.0`), se marca como `Fragmentacion_de_Ensamblaje` y **no** se suma como copia extra.
   - Si supera el 120% acumulado pero con fragmentos cortos ($< 50\%$), se marca como `Paralogo_Dudoso`.

3. **Conteo final de CNV:** La matriz `cnv_matrix.tsv` cuenta únicamente las copias funcionales confirmadas (`Copia_Principal` + `Paralogo_Confirmado`). Si el archivo de entrada incluye la columna `Validacion_Espacial`, el script también calcula el desempeño de transferencia entre especies y guarda `transferability_benchmark.tsv`.

#### Ejemplo de Uso
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

#### Parámetros Disponibles

| Corto | Largo | Tipo | Requerido | Valor por Defecto | Descripción |
| :--- | :--- | :--- | :---: | :--- | :--- |
| `-i` | `--input` | Archivo | **SÍ** | *Ninguno* | Archivo `.tsv`, `.csv` o `.rds` con las métricas de CAPTUS de las muestras a analizar. |
| `-m` | `--models_dir` | Carpeta | No | `models/` | Carpeta donde se encuentran los archivos `.rds` de los modelos y umbrales. |
| `-o` | `--out_dir` | Carpeta | No | `results/` | Carpeta donde se guardarán las tablas de CNV y los gráficos. |
| `-c` | `--cov_threshold` | Decimal | No | `50.0` | Cobertura mínima individual (`pct_recovered` en %) para confirmar un parálogo real. |
| `-f` | `--frag_threshold` | Decimal | No | `120.0` | Cobertura acumulada máxima (%) para considerar fragmentos cortos como parte de una misma copia fragmentada. |
| `-v` | `--min_votes` | Entero | No | `4` | Mínimo de votos positivos (de 7) para aceptar una secuencia como `Homologo_Valido`. |
| `-t` | `--top_loci` | Entero | No | `30` | Cantidad de loci con más copias que se mostrarán en el mapa de calor. |
| `-h` | `--help` | Bandera | No | - | Muestra la ayuda en pantalla. |

#### Archivos Generados en `<out_dir>/`
- `cnv_matrix.tsv`: Matriz de copias por gen (filas: loci/ortogrupos, columnas: muestras).
- `cnv_annotated_sequences.tsv`: Detalle de todas las secuencias válidas con su cobertura acumulada y su categoría (`Copia_Principal`, `Paralogo_Confirmado`, `Fragmentacion_de_Ensamblaje` o `Paralogo_Dudoso`).
- `functional_genes_catalog.tsv`: Lista filtrada únicamente con las copias funcionales (`Copia_Principal` y `Paralogo_Confirmado`).
- `consensus_verdicts.tsv`: Tabla con los votos de cada uno de los 7 modelos para todas las secuencias evaluadas.
- `transferability_benchmark.tsv`: *(Se genera automáticamente si el archivo de entrada tiene la columna `Validacion_Espacial`)* Métricas de rendimiento de los modelos al aplicarlos sobre otra especie.
- `plot_consensus_distribution.pdf` y `.png`: Gráfico de barras con el conteo de secuencias aceptadas y descartadas.
- `plot_cnv_heatmap.pdf` y `.png`: Mapa de calor con el número de copias en los loci más variables.

### 3.4 Compilación Manual del Reporte HTML en Quarto

El script `02_train_models.R` genera el reporte automáticamente. Si deseas volver a compilar `reporte_entrenamiento.html` sin repetir el entrenamiento:

```bash
# En una línea con targetcnv.sh:
./targetcnv.sh report

# O directamente con Quarto:
quarto render training_report.qmd \
  -P project_root:"." \
  -P models_dir:"models" \
  -P data_file:"data/data_prepared.rds" \
  --output "reporte_entrenamiento.html"
```

## 4. Flujo 1: Prueba Rápida con Datos Simulados (MWE)

El repositorio incluye un conjunto pequeño de datos sintéticos en `data/test_data/` para comprobar en pocos segundos que todo funciona correctamente.

### Ejecución Aislada (sin modificar las carpetas `models/` ni `results/` de la raíz)
```bash
./targetcnv.sh test
# Equivalente a: ./run_test_pipeline.sh
```

Este comando ejecuta los tres pasos dentro de `data/test_data/`:
1. `01_prepare_data.R` con los archivos `mock_*.gff3`, guardando en `data/test_data/prepared/`.
2. `02_train_models.R` con validación cruzada de 3 pliegues, guardando en `data/test_data/models/` y generando `data/test_data/reporte_entrenamiento_test.html`.
3. `03_cnv_inference.R` sobre `data/test_data/mock_captus_extract.tsv`, guardando la matriz de CNV y los gráficos en `data/test_data/results/`.

### Opciones Adicionales de Prueba
- `./targetcnv.sh test --no-report`: Ejecuta toda la prueba omitiendo la compilación de Quarto para terminar más rápido.
- `./targetcnv.sh test --sync-root`: Copia los modelos y resultados de prueba a las carpetas principales `models/` y `results/`. *(Atención: esto reemplaza los modelos que tengas guardados en `models/`).*

## 5. Flujo 2: Descarga y Uso de Modelos Preentrenados desde Zenodo

Los modelos entrenados con el genoma de referencia de *Solanum torvum* están depositados en Zenodo:
- **DOI:** [`10.5281/zenodo.23463987`](https://doi.org/10.5281/zenodo.23463987)
- **Enlace directo al registro:** [https://zenodo.org/records/23463987](https://www.google.com/search?q=https://zenodo.org/records/23463987)
- **Archivo:** `TargetCNV_Models.zip` (~382 KB)
- **Código MD5:** `13ec7b35711158079f962a96314eb671`

Si quieres analizar muestras nuevas directamente sin entrenar modelos desde cero, tienes dos opciones:

### Opción A: Descarga y Análisis Automático en 1 Línea
```bash
# Si la carpeta models/ está vacía, descarga el zip de Zenodo, verifica el MD5 y corre la inferencia:
./targetcnv.sh infer -i data/CAPTUS_Muestras/captus-extract_stats.tsv --auto-zenodo
```

### Opción B: Paso a Paso
```bash
# Paso 1: Descargar, comprobar MD5 y descomprimir los modelos con targetcnv.sh
./targetcnv.sh zenodo

# (O manualmente con curl, md5sum y unzip)
curl -L -o TargetCNV_Models.zip https://zenodo.org/api/records/23463987/files/TargetCNV_Models.zip/content
md5sum TargetCNV_Models.zip
# Debe mostrar: 13ec7b35711158079f962a96314eb671  TargetCNV_Models.zip
unzip -o TargetCNV_Models.zip

# Paso 2: Ejecutar la inferencia de CNV sobre tus muestras
./scripts/03_cnv_inference.R \
  --input data/CAPTUS_Muestras/captus-extract_stats.tsv \
  --models_dir models \
  --out_dir results \
  --cov_threshold 50.0 \
  --frag_threshold 120.0 \
  --min_votes 4 \
  --top_loci 30
```

## 6. Flujo 3: Entrenamiento con Nuevos Datos y Transferibilidad entre Especies

### 6.1 Entrenar y Analizar una Nueva Especie en 1 Línea
Si cuentas con el genoma de referencia anotado (`.gff3`), los cóntigos de CAPTUS (`.gff`), la tabla `captus-extract_stats.tsv` y una tabla de correspondencia entre genes y ortogrupos:

```bash
./targetcnv.sh run-all \
  -g ruta/a/nuevo_genoma.gff3 \
  -c ruta/a/nuevos_contigs_captus.gff \
  -s ruta/a/captus-extract_stats.tsv \
  -d ruta/a/diccionario_homologia.xlsx \
  -p NombreMuestraReferencia \
  -i ruta/a/captus_muestras_poblacion.tsv \
  --out_data data/nueva_especie \
  -m models/nueva_especie \
  -o results/nueva_especie \
  -r reporte_entrenamiento_nueva_especie.html
```

### 6.2 Pruebas de Transferibilidad en Otras Especies (*S. citrullifolium* y *S. stellatiglandulosum*)

También puedes evaluar qué tan bien funcionan los modelos entrenados en *S. torvum* cuando se aplican sobre otras especies con distinta distancia evolutiva o calidad de ensamblaje:

#### Caso A: Genoma ensamblado a nivel de cromosomas (*Solanum citrullifolium*)
```bash
# 1. Preparar los datos de referencia de S. citrullifolium (armoniza automáticamente prefijos de cromosomas)
./scripts/01_prepare_data.R \
  --gff_ref ruta/a/Scitrullifolium_genes.gff3 \
  --captus_gff ruta/a/Scitrullifolium_hit_contigs.gff \
  --captus_stats ruta/a/captus-extract_stats.tsv \
  --dict_ref data/reference/Genes_Capture_May_2026.xlsx \
  --sample Scitrullifolium \
  --out_dir data/citrullifolium/

# 2. Evaluar los modelos de S. torvum sobre S. citrullifolium
./scripts/03_cnv_inference.R \
  --input data/citrullifolium/data_prepared.rds \
  --models_dir models/ \
  --out_dir results/citrullifolium/
```

#### Caso B: Genoma altamente fragmentado (*Solanum stellatiglandulosum*, puente HOG $\leftrightarrow$ OG y atributo `Target`)
```bash
# 1. Preparar los datos leyendo el locus directamente del atributo Target del GFF (--target_mode target_attr)
./scripts/01_prepare_data.R \
  --captus_gff ruta/a/S_stella_hit_contigs.gff \
  --captus_stats ruta/a/captus-extract_stats.tsv \
  --dict_ref data/reference/Genes_Capture_May_2026.xlsx \
  --sample S_stella \
  --target_mode target_attr \
  --out_dir data/stellatiglandulosum/

# 2. Evaluar la transferencia y guardar la tabla comparativa en results/stellatiglandulosum/transferability_benchmark.tsv
./scripts/03_cnv_inference.R \
  --input data/stellatiglandulosum/data_prepared.rds \
  --models_dir models/ \
  --out_dir results/stellatiglandulosum/
```

## 7. Flujo 4: Reproducción del Experimento de Referencia (*Solanum torvum*)

### En 1 Sola Línea
```bash
./targetcnv.sh reproduce-storvum --cores 4
```

### Paso a Paso con los Scripts Individuales
```bash
# Paso 1: Cruce espacial y preparación de datos (genera 184 secuencias NUC: 87 TP y 97 FP)
./scripts/01_prepare_data.R \
  --gff_ref data/reference/Storvum_genes.gff3 \
  --captus_gff data/reference/Storvum_hit_contigs.gff \
  --captus_stats data/reference/captus-extract_stats.tsv \
  --dict_ref data/reference/Genes_Capture_May_2026.xlsx \
  --sample Storvum \
  --out_dir data/

# Paso 2: Entrenamiento de modelos, cortes de Youden, KernelSHAP y reporte HTML
./scripts/02_train_models.R \
  --data data/data_prepared.rds \
  --out_models models \
  --cv_folds 5 \
  --sampling down \
  --cores 4 \
  --seed 1234 \
  --report_file reporte_entrenamiento.html

# Paso 3: Inferencia en muestras biológicas y obtención de la matriz de CNV
./scripts/03_cnv_inference.R \
  --input data/CAPTUS_Muestras/captus-extract_stats.tsv \
  --models_dir models \
  --out_dir results \
  --cov_threshold 50.0 \
  --frag_threshold 120.0 \
  --min_votes 4 \
  --top_loci 30
```

### Resultados en Validación Cruzada (K=5, Semilla 1234)

#### Escenario Balanceado (*Down-Sampling*)
| Modelo | Exactitud Bal. (%) | Sensibilidad (%) | Especificidad (%) | Precisión (%) | FDR (%) | F1 Score | MCC | Kappa | ROC-AUC |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Modelo Heurístico** | 88.86 | 83.91 | 93.81 | 92.41 | 7.59 | **0.880** | 0.784 | 0.781 | 0.921 |
| **Red Neuronal (`nnet`)** | 87.32 | 83.91 | 90.72 | 88.98 | 11.02 | **0.864** | 0.750 | 0.749 | 0.927 |
| **SVM Polinomial (`svmPoly`)** | 87.14 | 80.46 | 93.81 | 92.11 | 7.89 | **0.859** | 0.753 | 0.748 | **0.934** |
| **LDA** | 86.74 | 82.76 | 90.72 | 88.89 | 11.11 | **0.857** | 0.739 | 0.737 | 0.922 |
| **Naive Bayes (`nb`)** | 86.56 | 79.31 | 93.81 | 92.00 | 8.00 | **0.852** | 0.743 | 0.736 | **0.934** |
| **Elastic Net (`glmnet`)** | 85.28 | 83.91 | 86.60 | 84.88 | 15.12 | **0.844** | 0.706 | 0.705 | 0.930 |

#### Escenario Desbalanceado (Proporción Original: 87 TP vs. 97 FP)
| Modelo | Exactitud Bal. (%) | Sensibilidad (%) | Especificidad (%) | Precisión (%) | FDR (%) | F1 Score | MCC | Kappa | ROC-AUC |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Modelo Heurístico** | 88.86 | 83.91 | 93.81 | 92.41 | 7.59 | **0.880** | 0.784 | 0.781 | 0.921 |
| **LDA** | 87.26 | 82.76 | 91.75 | 90.00 | 10.00 | **0.862** | 0.750 | 0.748 | 0.922 |
| **Naive Bayes (`nb`)** | 86.56 | 79.31 | 93.81 | 92.00 | 8.00 | **0.852** | 0.743 | 0.736 | **0.934** |
| **Elastic Net (`glmnet`)** | 85.76 | 83.91 | 87.63 | 85.88 | 14.12 | **0.849** | 0.716 | 0.716 | 0.930 |
| **Red Neuronal (`nnet`)** | 84.60 | 80.46 | 88.66 | 86.42 | 13.58 | **0.833** | 0.695 | 0.694 | 0.921 |
| **SVM Polinomial (`svmPoly`)** | 83.04 | 90.80 | 75.26 | 76.70 | 23.30 | **0.832** | 0.664 | 0.654 | 0.930 |

## 8. Resumen de Mejoras y Correcciones Aplicadas (Conservando los Resultados)

Para mejorar la facilidad de uso sin cambiar los resultados numéricos ni las reglas biológicas de `Workflow_Central_Sin_Rescate.R`, se realizaron los siguientes ajustes:

1. **Script unificado `targetcnv.sh`:** Permite correr pruebas, revisar dependencias, descargar modelos de Zenodo verificando su MD5, o ejecutar los tres módulos en una sola línea de Bash.
2. **Corrección en la lectura de diccionarios Excel (`load_reference_dictionary`):** Se reemplazó la llamada `readxl::excel_sheet()` (que producía error porque la función en `readxl` se llama `readxl::read_excel()`) y se integró el filtrado exacto de vías metabólicas y ensamblajes de `Workflow_Central_Sin_Rescate.R` para *S. torvum*, *S. citrullifolium* y *S. stellatiglandulosum* (incluyendo la tabla puente `HOG_ID <-> OG_ID` y el modo `--target_mode target_attr`).
3. **Formato único en tablas de métricas y funciones de graficado:** Todas las tablas de resultados comparten las mismas 16 columnas calculadas con `compute_classification_metrics()`, se estandarizó el formato `OGxxxx__00` con doble guion bajo en `clean_captus_ids()`, y se ajustó `plot_cnv_heatmap()` en `R/plotting.R` para aceptar tanto tablas de secuencias anotadas como matrices de conteo.
4. **Plantilla `training_report.qmd` incluida en la raíz:** Se incorporó el archivo `training_report.qmd` (que no estaba en la rama principal de GitHub) con tablas interactivas `DT`, gráficos `plotly`, matrices de confusión y gráficos `shapviz`, asegurando además que el archivo HTML se guarde en la ruta indicada por `--report_file`.
5. **Opciones para torneo de 18 modelos, GLM Stacking y transferibilidad:** Se añadieron las opciones `--extended_tournament` y `--train_stacking` en `02_train_models.R`, y la generación automática de `transferability_benchmark.tsv` en `03_cnv_inference.R` cuando los datos de entrada incluyen etiquetas de referencia.
