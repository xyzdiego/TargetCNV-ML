# TargetCNV-ML: Detección de Variaciones en el Número de Copias en Datos de Captura de Genes

**Idiomas / Languages:** [Español (README\_es.md)](https://www.google.com/search?q=README_es.md) | [English (README.md)](README.md)  
**Guías de Ejecución:** [Guía de Ejecución (GUIA\_EJECUCION.md)](https://www.google.com/search?q=GUIA_EJECUCION.md) | [Execution Guide (EXECUTION\_GUIDE.md)](https://www.google.com/search?q=EXECUTION_GUIDE.md)

[](https://doi.org/10.5281/zenodo.23463987)





**TargetCNV-ML** es un pipeline bioinformático de línea de comandos (CLI) en R y Bash diseñado para detectar **Variaciones en el Número de Copias (CNV)** e identificar **expansiones de familias génicas** a partir de ensamblajes de secuenciación por captura de genes (*Target Capture* / *Hyb-Seq*).

## ¿Qué problema resuelve TargetCNV-ML?

La captura de genes permite secuenciar cientos de loci de interés a bajo costo en estudios de filogenómica, evolución y mejoramiento genético. Sin embargo, estimar el número de copias génicas directamente desde los cóntigos ensamblados *de novo* presenta dos dificultades frecuentes:

1. **Alineamientos cruzados y ruido fuera de objetivo:** Las sondas de captura pueden hibridar con parálogos cercanos, pseudogenes o regiones intrónicas/intergénicas. Esto genera cóntigos que se alinean de forma ambigua y aumentan los falsos positivos.
2. **Fragmentación técnica del ensamblaje:** Los ensambladores *de novo* suelen dividir un mismo gen largo en dos o más fragmentos cortos. Si se cuentan todos los cóntigos por igual, estos fragmentos se confunden con duplicaciones génicas reales.

**TargetCNV-ML** resuelve ambos puntos en tres etapas:
- **Cruce espacial de coordenadas (`GenomicRanges`):** Compara las coordenadas de los cóntigos ensamblados (`protein_match:NUC`) con la anotación del genoma de referencia (`.gff3`) y una tabla de ortología de referencia para distinguir genes reales (`TRUE_POSITIVE`) de alineamientos cruzados, ruido o copias técnicas redundantes (`FALSE_POSITIVE`).
- **Clasificación por consenso de 7 modelos:** Combina **4 modelos supervisados de Machine Learning** (`ElasticNet`, `NaiveBayes`, `SVM_Polinomial`, `Red_Neuronal`), **2 puntos de corte univariados basados en el índice de Youden** (`pct_recovered`, `pct_identity`) y **1 modelo heurístico ponderado**. Una secuencia se acepta como homólogo válido si recibe al menos 4 votos positivos de los 7 métodos.
- **Control de fragmentación por cobertura:** Distingue copias adicionales completas (`Paralogo_Confirmado`, cobertura individual $\ge 50\%$) de fragmentos técnicos de un mismo gen (`Fragmentacion_de_Ensamblaje`, cobertura individual $< 50\%$ y cobertura acumulada en el locus $\le 120\%$).

## Estructura del Repositorio

```text
TargetCNV-ML/
├── targetcnv.sh               # Script principal en Bash para correr los flujos en una línea
├── run_test_pipeline.sh       # Script de prueba con datos simulados (MWE)
├── training_report.qmd        # Plantilla Quarto para el reporte HTML interactivo
├── README.md                  # Resumen general en inglés
├── README_es.md               # Resumen general en español
├── EXECUTION_GUIDE.md         # Guía de ejecución en inglés
├── GUIA_EJECUCION.md          # Guía de ejecución en español
├── renv.lock                  # Versiones exactas de los paquetes de R
├── .Rprofile                  # Activación automática del entorno renv al iniciar R
├── LICENSE                    # Licencia MIT
│
├── scripts/                   # Scripts ejecutables desde la terminal (CLI)
│   ├── 01_prepare_data.R      # Paso 1: Cruce espacial GFF, tabla de referencia y deltas
│   ├── 02_train_models.R      # Paso 2: Cortes de Youden, entrenamiento ML, SHAP y reporte HTML
│   └── 03_cnv_inference.R     # Paso 3: Votación de 7 modelos y generación de matriz CNV
│
├── R/                         # Funciones internas modulares en R
│   ├── data_processing.R      # Filtrado GFF, solape espacial, lectura de diccionarios y CNV
│   ├── model_evaluation.R     # Métricas de clasificación, cortes de Youden y consenso
│   └── plotting.R             # Gráficos de densidades, matrices de confusión y mapas de calor
│
├── data/                      # Datos de entrada y conjunto de prueba
│   ├── reference/             # GFF3 de S. torvum, cóntigos CAPTUS y diccionario de genes
│   ├── CAPTUS_Muestras/       # Estadísticas de alineamiento CAPTUS para muestras biológicas
│   ├── test_data/             # Datos sintéticos ligeros para pruebas rápidas (< 25 KB)
│   └── data_prepared.rds      # Datos procesados de S. torvum (184 secuencias NUC)
│
├── models/                    # Modelos entrenados (.rds) y tablas de métricas (.tsv)
└── results/                   # Matrices de CNV (.tsv), catálogos de genes y figuras (.pdf/.png)
```

## Requisitos e Instalación

### Requisitos del Sistema
- **Sistema Operativo:** Linux, macOS o Windows (mediante WSL2).
- **R:** $\ge 4.0.0$ (se recomienda R 4.3+ o 4.6).
- **Quarto CLI:** $\ge 1.3+$ (opcional, solo necesario para compilar el reporte interactivo `reporte_entrenamiento.html`).
- **Utilidades de terminal:** `bash`, `curl` o `wget`, `unzip`, `md5sum`.

### Pasos de Instalación

```bash
# 1. Clonar el repositorio
git clone https://github.com/xyzdiego/TargetCNV-ML.git
cd TargetCNV-ML

# 2. Restaurar los paquetes de R definidos en renv.lock
Rscript -e "renv::restore()"

# 3. Dar permisos de ejecución, revisar el entorno y correr la prueba rápida (MWE)
chmod +x targetcnv.sh run_test_pipeline.sh scripts/*.R
./targetcnv.sh doctor
./targetcnv.sh test
```

## Ejecución Rápida en Una Línea (`./targetcnv.sh`)

El script `./targetcnv.sh` permite ejecutar cualquiera de las tareas del proyecto con un solo comando:

| Tarea | Comando en Terminal |
| :--- | :--- |
| **Verificar herramientas y paquetes de R** | `./targetcnv.sh doctor` |
| **Ejecutar la prueba rápida con datos simulados (MWE)** | `./targetcnv.sh test` |
| **Descargar modelos preentrenados desde Zenodo** | `./targetcnv.sh zenodo` |
| **Predecir CNV en muestras nuevas (descargando modelos si faltan)** | `./targetcnv.sh infer -i data/CAPTUS_Muestras/captus-extract_stats.tsv --auto-zenodo` |
| **Reproducir el experimento de *Solanum torvum*** | `./targetcnv.sh reproduce-storvum --cores 4` |
| **Correr los 3 pasos completos con una nueva especie** | `./targetcnv.sh run-all -g ref.gff3 -c captus.gff -s stats.tsv -d dict.xlsx -p Muestra -i pob_stats.tsv` |
| **Volver a generar el reporte HTML de entrenamiento** | `./targetcnv.sh report` |

## Flujo de Trabajo Paso a Paso (CLI)

También puedes ejecutar cada uno de los tres scripts de R por separado:

### Paso 1: Preparación de Datos (`01_prepare_data.R`)
Cruza las coordenadas de los cóntigos de CAPTUS (`protein_match:NUC`) con los genes del genoma de referencia (`gene`), compara los resultados con el diccionario de ortología, marca como `FALSE_POSITIVE` los hits secundarios redundantes de un mismo gen y calcula las diferencias de calidad dentro de cada locus (`delta_wscore`, `delta_pct_identity`).

```bash
./scripts/01_prepare_data.R \
  --gff_ref data/reference/Storvum_genes.gff3 \
  --captus_gff data/reference/Storvum_hit_contigs.gff \
  --captus_stats data/reference/captus-extract_stats.tsv \
  --dict_ref data/reference/Genes_Capture_May_2026.xlsx \
  --sample Storvum \
  --out_dir data/
```

### Paso 2: Entrenamiento y Evaluación de Modelos (`02_train_models.R`)
Calcula los puntos de corte de Youden para cada métrica, busca la mejor combinación de pesos para el modelo heurístico, entrena los clasificadores de Machine Learning con validación cruzada de 5 particiones (con y sin balanceo de clases), calcula los valores KernelSHAP para explicar el peso de cada variable y genera el reporte HTML.

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
*(Opciones adicionales: `--extended_tournament` entrena los 18 clasificadores explorados en el estudio original; `--train_stacking` ajusta un modelo GLM que combina los mejores clasificadores).*

### Paso 3: Inferencia por Consenso y Matriz de CNV (`03_cnv_inference.R`)
Evalúa nuevas muestras con los 7 modelos, aplica la regla de votación por mayoría ($\ge 4/7$ votos), filtra los fragmentos de ensamblaje y guarda la matriz de CNV (`Locus x Muestras`), el catálogo de copias funcionales y los gráficos de resumen.

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

## Resultados de Referencia (*Solanum torvum*)

Desempeño sobre 184 secuencias `NUC` de *Solanum torvum* (87 Verdaderos Positivos y 97 Falsos Positivos) en validación cruzada de 5 particiones (`seed = 1234`, escenario balanceado con submuestreo):

| Modelo / Método | Exactitud Bal. (%) | Sensibilidad (%) | Especificidad (%) | Precisión (%) | FDR (%) | F1-Score | MCC | ROC-AUC |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Modelo Heurístico Ponderado** | 88.86 | 83.91 | 93.81 | 92.41 | 7.59 | **0.880** | 0.784 | 0.921 |
| **Red Neuronal (`nnet`)** | 87.32 | 83.91 | 90.72 | 88.98 | 11.02 | **0.864** | 0.750 | 0.927 |
| **SVM Polinomial (`svmPoly`)** | 87.14 | 80.46 | 93.81 | 92.11 | 7.89 | **0.859** | 0.753 | **0.934** |
| **LDA** | 86.74 | 82.76 | 90.72 | 88.89 | 11.11 | **0.857** | 0.739 | 0.922 |
| **Naive Bayes (`nb`)** | 86.56 | 79.31 | 93.81 | 92.00 | 8.00 | **0.852** | 0.743 | **0.934** |
| **Elastic Net (`glmnet`)** | 85.28 | 83.91 | 86.60 | 84.88 | 15.12 | **0.844** | 0.706 | 0.930 |

## Modelos Preentrenados en Zenodo

Los modelos calibrados con *Solanum torvum* están publicados en Zenodo:
- **DOI:** [`10.5281/zenodo.23463987`](https://doi.org/10.5281/zenodo.23463987)
- **Archivo:** `TargetCNV_Models.zip` (`MD5: 13ec7b35711158079f962a96314eb671`)
- **Descarga automática:** Ejecuta `./targetcnv.sh zenodo` para descargar el archivo, comprobar el MD5 y descomprimir los modelos en `models/`.

Para consultar todas las opciones y flujos de trabajo paso a paso, revisa la [Guía de Ejecución (GUIA\_EJECUCION.md)](https://www.google.com/search?q=GUIA_EJECUCION.md).

## Citación y Licencia

Si utilizas **TargetCNV-ML** en tu investigación, por favor cita:

```bibtex
@software{TargetCNV_ML_2026,
  author       = {Diego Torres and Contributors},
  title        = {TargetCNV-ML: Machine Learning & Spatial Inference for Copy Number Variations in Target Capture},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.23463987},
  url          = {https://github.com/xyzdiego/TargetCNV-ML}
}
```

Distribuido bajo la **Licencia MIT** — consulta el archivo [LICENSE](https://www.google.com/search?q=LICENSE) para más detalles.
