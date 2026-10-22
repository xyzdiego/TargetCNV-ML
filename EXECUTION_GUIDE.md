# TargetCNV-ML Execution Guide

**Languages / Idiomas:** [English (EXECUTION\_GUIDE.md)](https://www.google.com/search?q=EXECUTION_GUIDE.md) | [Español (GUIA\_EJECUCION.md)](https://www.google.com/search?q=GUIA_EJECUCION.md)  
**General Overview:** [English README (README.md)](README.md) | [README en Español (README\_es.md)](https://www.google.com/search?q=README_es.md)

**TargetCNV-ML** is a command-line (CLI) bioinformatics pipeline that detects **Copy Number Variations (CNV)** and identifies **homologous gene expansions** from target capture (*Hyb-Seq* / *Target Enrichment*) sequencing assemblies.

In this guide you will find:
1. System requirements, installation steps, and environment check (`./targetcnv.sh doctor`).
2. Directory and file layout of the repository.
3. One-line execution with `./targetcnv.sh` and a complete reference of all options for the three CLI modules (`01_prepare_data.R`, `02_train_models.R`, `03_cnv_inference.R`) and the HTML report (`training_report.qmd`).
4. **Workflow 1:** Running the quick test with the included synthetic dataset (*Minimal Working Example* / MWE).
5. **Workflow 2:** Downloading and running the pretrained models hosted on **Zenodo** ([`DOI: 10.5281/zenodo.23463987`](https://doi.org/10.5281/zenodo.23463987)) without retraining.
6. **Workflow 3:** Training models with new genomic data and testing **cross-species transferability** (on chromosome-level assemblies like *Solanum citrullifolium* and fragmented assemblies like *Solanum stellatiglandulosum*).
7. **Workflow 4:** Reproducing the ***Solanum torvum*** reference experiment step by step.
8. Summary of technical and user-interface improvements that keep 100% of the original numerical results from `Workflow_Central_Sin_Rescate.R`.

## 1. System Requirements & Installation

### Software Requirements
- **Operating System:** Linux (Ubuntu 20.04+, Debian 11+, Arch/CachyOS, Fedora), macOS, or Windows (via WSL2).
- **R:** $\ge 4.0.0$ (R 4.3+ or 4.6 recommended).
- **Quarto CLI:** $\ge 1.3+$ (required to build the interactive `.html` report).
- **System tools:** `bash`, `curl` or `wget`, `unzip`, `md5sum`, `tar`.

### Installing & Restoring the R Environment (`renv`)

The project uses `renv` so that R installs the exact same package versions used to develop and validate the pipeline.

```bash
# 1. Clone the repository
git clone https://github.com/xyzdiego/TargetCNV-ML.git
cd TargetCNV-ML

# 2. Make scripts executable and restore R packages
chmod +x targetcnv.sh run_test_pipeline.sh scripts/*.R
Rscript -e "renv::restore()"

# 3. Check that all tools and R libraries are ready
./targetcnv.sh doctor
```

> [!NOTE]
> The main R libraries managed by `renv` are: `caret`, `glmnet`, `klaR`, `pROC`, `GenomicRanges`, `rtracklayer`, `optparse`, `readxl`, `kernelshap`, `shapviz`, `NeuralNetTools`, `plotly`, `DT`, `cowplot`, `dplyr`, `tidyr`, and `stringr`.

## 2. Directory & File Structure

```text
TargetCNV-ML/
├── .Rprofile                    # Activates renv automatically when R starts
├── .gitignore                   # Rules to ignore temporary or large files
├── renv.lock                    # Exact versions of R and installed packages
├── targetcnv.sh                 # Main Bash script to run workflows in 1 line
├── run_test_pipeline.sh         # Script that runs the test dataset (MWE)
├── training_report.qmd          # Quarto template for the HTML training report
├── EXECUTION_GUIDE.md           # Execution Guide in English (this file)
├── GUIA_EJECUCION.md            # Execution Guide in Spanish
├── README.md                    # Main overview in English
├── README_es.md                 # Main overview in Spanish
│
├── scripts/                     # Command-line R scripts (CLI)
│   ├── 01_prepare_data.R        # Step 1: GFF coordinate overlap & TP/FP labeling
│   ├── 02_train_models.R        # Step 2: Model training, Youden cutoffs & HTML report
│   └── 03_cnv_inference.R       # Step 3: 7-model consensus vote & CNV matrix
│
├── R/                           # Internal modular R functions
│   ├── data_processing.R        # GFF reading, dictionaries & CNV copy counting
│   ├── model_evaluation.R       # Performance metrics, Youden index & consensus rule
│   └── plotting.R               # Functions to generate ggplot2 and plotly figures
│
├── data/
│   ├── reference/               # Reference genomic files for S. torvum
│   ├── CAPTUS_Muestras/         # CAPTUS alignment stats for biological samples
│   ├── test_data/               # Small synthetic dataset for testing (< 25 KB)
│   ├── data_prepared.rds        # Processed dataset ready for training (S. torvum)
│   └── data_prepared.tsv        # Tab-separated text version of the processed dataset
│
├── models/                      # Saved models (.rds) and benchmark tables (.tsv)
└── results/                     # CNV matrices, consensus tables & plots (.pdf/.png)
```

## 3. CLI Command Reference

### 3.0 Quick Execution with `targetcnv.sh`

The `./targetcnv.sh` helper script groups the most common operations into simple subcommands:

| Subcommand | Example | What it does |
| :--- | :--- | :--- |
| `doctor` | `./targetcnv.sh doctor` | Checks whether `Rscript`, `quarto`, `curl`, `unzip`, `md5sum`, and all R packages are installed. |
| `test` | `./targetcnv.sh test [--sync-root] [--no-report]` | Runs the synthetic test dataset inside `data/test_data/`. |
| `zenodo` | `./targetcnv.sh zenodo [--out-dir .]` | Downloads `TargetCNV_Models.zip` from Zenodo, checks its MD5 hash, and extracts `models/`. |
| `infer` | `./targetcnv.sh infer -i stats.tsv --auto-zenodo` | Estimates the CNV matrix on new samples; with `--auto-zenodo`, downloads the models first if `models/` is empty. |
| `reproduce-storvum` | `./targetcnv.sh reproduce-storvum --cores 4` | Runs Steps 01, 02, and 03 from start to finish on the *Solanum torvum* dataset. |
| `run-all` | `./targetcnv.sh run-all -g ref.gff3 -c cap.gff -s stats.tsv -d dict.xlsx -i new.tsv` | Runs `01_prepare_data.R` $\rightarrow$ `02_train_models.R` $\rightarrow$ `03_cnv_inference.R` in sequence. |
| `prepare` | `./targetcnv.sh prepare [01 options]` | Passes all arguments directly to `scripts/01_prepare_data.R`. |
| `train` | `./targetcnv.sh train [02 options]` | Passes all arguments directly to `scripts/02_train_models.R`. |
| `report` | `./targetcnv.sh report --models-dir models --data data/data_prepared.rds` | Rebuilds the HTML report with Quarto without retraining the models. |

### 3.1 Module 01: Data Preparation (`01_prepare_data.R`)

This script identifies which CAPTUS contigs correspond to real genes in the reference genome and which ones are noise or redundant fragments:
1. Matches the genomic coordinates of CAPTUS contigs (`protein_match:NUC`) against `gene` features in the reference `.gff3` file using `GenomicRanges` (or reads the reference gene ID directly from the GFF `Target` attribute when `--target_mode target_attr` is used for fragmented assemblies).
2. Compares each match against the reference orthology dictionary to label true matches as `TRUE_POSITIVE` and cross-homolog or off-target alignments as `FALSE_POSITIVE`.
3. When multiple contigs map to the same reference gene as `TRUE_POSITIVE`, keeps the longest, highest-coverage contig as `TRUE_POSITIVE` and reclassifies the secondary fragments as `FALSE_POSITIVE` (`5_TP_Redundante_Degradado`).
4. Computes the quality difference of each hit relative to the best hit in the same locus (`delta_wscore` and `delta_pct_identity`).

#### Example Usage
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

#### Available Options

| Short | Long | Type | Required | Default | Description |
| :--- | :--- | :--- | :---: | :--- | :--- |
| `-g` | `--gff_ref` | File | **YES*** | *None* | Reference genome GFF/GFF3 annotation file (must contain `gene` entries). *Required in `overlap` mode. |
| `-c` | `--captus_gff` | File | **YES** | *None* | CAPTUS assembled contigs GFF file (`protein_match:NUC` entries). |
| `-s` | `--captus_stats` | File | **YES** | *None* | CAPTUS `captus-extract_stats.tsv` alignment metrics table. |
| `-d` | `--dict_ref` | File | No | `NULL` | Gene-to-orthogroup correspondence table (`.xlsx`, `.tsv`, or `.csv`). Supports the multi-sheet `Genes_Capture_May_2026.xlsx` file or simple tables. |
| `-p` | `--sample` | Text | No | `NULL` | Sample name to filter inside `captus_stats` (e.g., `Storvum`, `Scitrullifolium`, `S_stella`). |
| | `--target_mode` | Text | No | `overlap` | `'overlap'`: matches coordinates against `--gff_ref`. `'target_attr'`: extracts the reference gene ID from the `Target` attribute of the CAPTUS GFF (for fragmented assemblies). |
| `-o` | `--out_dir` | Folder | No | `data/` | Output directory where prepared files will be saved. |
| `-h` | `--help` | Flag | No | - | Shows the help menu. |

#### Output Files
- `<out_dir>/data_prepared.rds`: R data frame containing `Validacion_Espacial`, `Doble_Control_Biologico`, and all alignment metrics.
- `<out_dir>/data_prepared.tsv`: Tab-separated text version of the same table.

### 3.2 Module 02: Model Training & Evaluation (`02_train_models.R`)

Using the labeled dataset from Step 1, this script trains and compares classification methods to separate real homologs (`TP`) from noise (`FP`):
1. **Univariate cutoffs (Youden):** For each of the 4 continuous metrics (`pct_recovered`, `pct_identity`, `score`, `wscore`), finds the cutoff that maximizes Youden's $J$ statistic ($\text{Sensitivity} + \text{Specificity} - 1$).
2. **Weighted heuristic model:** Tests integer weight combinations (0 to 3) across the 4 binary Youden votes to find the weighted sum that maximizes the area under the ROC curve (AUC).
3. **Machine Learning models (`caret`):** Trains 5 classifiers (`ElasticNet`, `NaiveBayes`, `SVM_Polinomial`, `Red_Neuronal`, `LDA`) using cross-validation under two scenarios:
   - **Balanced:** Balances `TP` and `FP` counts inside each fold using down-sampling.
   - **Unbalanced:** Trains on the original class proportions without resampling.
4. **Explainability (KernelSHAP) & Report:** Computes KernelSHAP values on the Elastic Net model to show how each metric shifts the probability of `TP`, and compiles `reporte_entrenamiento.html`.

#### Example Usage
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

#### Available Options

| Short | Long | Type | Required | Default | Description |
| :--- | :--- | :--- | :---: | :--- | :--- |
| `-d` | `--data` | File | **YES** | `data/data_prepared.rds` | Prepared dataset file (`.rds` or `.tsv`) generated in Step 1. |
| `-m` | `--out_models` | Folder | No | `models/` | Folder where `.rds` models and `.tsv` benchmark tables will be saved. |
| `-k` | `--cv_folds` | Integer | No | `5` | Number of cross-validation folds ($\ge 2$). |
| `-s` | `--sampling` | Text | No | `down` | Class-balancing method for the balanced scenario: `'down'` (down-sampling), `'up'` (up-sampling), or `'none'`. |
| `-c` | `--cores` | Integer | No | `detectCores() - 1` | Number of CPU cores used for parallel training. |
| | `--seed` | Integer | No | `1234` | Random seed for reproducible results. |
| | `--extended_tournament` | Flag | No | `FALSE` | Trains all 18 classifiers tested in the original study (adds `RandomForest`, `GBM`, `XGBoost`, `DNN`, etc.). |
| | `--train_stacking` | Flag | No | `FALSE` | Fits a GLM stacking meta-model (`caretStack`) combining the top two classifiers (`model_stacking.rds`). |
| | `--include_hit` | Flag | No | `FALSE` | Adds the `hit` column as a predictor. *(Off by default to keep compatibility with `03_cnv_inference.R`).* |
| | `--render_report` | Flag | No | `TRUE` | Automatically compiles the HTML report when training finishes. |
| | `--no_render_report` | Flag | No | - | Skips HTML report generation (useful on servers without Quarto). |
| `-r` | `--report_file` | File | No | `reporte_entrenamiento.html` | Output path and filename for the HTML report. |
| `-h` | `--help` | Flag | No | - | Shows the help menu. |

#### Metrics Included in the Benchmark Tables (16 Columns)

All benchmark tables (`benchmark_models_balanced.tsv`, `benchmark_models_unbalanced.tsv`, `youden_benchmark.tsv`, and `benchmark_models_all.tsv`) share the same 16 columns:

1. `Modelo`: Name of the classifier or cutoff rule.
2. `Exactitud` (*Accuracy*): Overall proportion of correct predictions $\frac{TP + TN}{TP + FP + FN + TN}$.
3. `Exactitud_Balanceada` (*Balanced Accuracy*): Average of sensitivity and specificity $\frac{\text{Sensitivity} + \text{Specificity}}{2}$.
4. `Sensibilidad` (*Recall*): Proportion of true positives detected $\frac{TP}{TP + FN}$.
5. `Especificidad`: Proportion of false positives correctly rejected $\frac{TN}{TN + FP}$.
6. `Precision` (*Positive Predictive Value*): Proportion of positive calls that are truly positive $\frac{TP}{TP + FP}$.
7. `VPN` (*Negative Predictive Value*): $\frac{TN}{TN + FN}$.
8. `FDR` (*False Discovery Rate*): $1 - \text{Precision} = \frac{FP}{TP + FP}$.
9. `F1_Score`: Harmonic mean of precision and sensitivity $2 \times \frac{\text{Precision} \times \text{Sensitivity}}{\text{Precision} + \text{Sensitivity}}$.
10. `MCC` (*Matthews Correlation Coefficient*): Correlation metric between $-1$ and $+1$ that accounts for all four confusion matrix cells.
11. `Kappa`: Cohen's Kappa agreement score.
12. `ROC_AUC`: Area under the ROC curve.
13. `TP`, `FP`, `FN`, `TN`: Confusion matrix integer counts.

#### Output Files in `<out_models>/`
- `all_trained_models.rds` & `models_balanced.rds`: List of models trained in the balanced scenario.
- `models_unbalanced.rds`: List of models trained without class balancing.
- `model_elasticnet.rds`, `model_naivebayes.rds`, `model_svm_polinomial.rds`, `model_red_neuronal.rds`, `model_lda.rds`: Individual models used in Step 3.
- `youden_thresholds.rds`: Optimal cutoff values for each metric.
- `heuristic_model.rds`: Weights and cutoff for the weighted heuristic model.
- `confusion_matrices.rds`: $2 \times 2$ confusion matrices for all evaluated methods.
- `shap_values.rds`: KernelSHAP values used for variable importance plots.
- `benchmark_models_balanced.tsv` / `.rds`, `benchmark_models_unbalanced.tsv` / `.rds`, `benchmark_models_all.tsv` / `.rds`, and `youden_benchmark.tsv` / `.rds`.

### 3.3 Module 03: CNV Inference & Consensus (`03_cnv_inference.R`)

This script takes a CAPTUS statistics table from new samples, removes assembly noise by combining 7 models, and counts functional gene copies per locus and sample:

1. **7-Model consensus vote:** Each sequence is evaluated by 7 methods:
   - 4 Machine Learning classifiers (`ElasticNet`, `NaiveBayes`, `SVM_Polinomial`, `Red_Neuronal`).
   - 2 univariate Youden cutoffs (`pct_recovered` and `pct_identity`).
   - 1 weighted heuristic model.
   If a sequence receives **4 or more positive votes** (`--min_votes 4`), it is accepted as `Homologo_Valido`; otherwise it is discarded as `Ruido_Descartado`.

2. **Fragmentation control rules:** Within each sample and locus, sequences accepted as `Homologo_Valido` are sorted from highest to lowest coverage (`pct_recovered`):
   - The first sequence (highest coverage) is labeled `Copia_Principal` (primary copy).
   - Any additional sequence with individual coverage $\ge 50\%$ (`--cov_threshold 50.0`) is labeled `Paralogo_Confirmado` (confirmed paralog).
   - If an additional sequence has coverage $< 50\%$ and the cumulative coverage in the locus is $\le 120\%$ (`--frag_threshold 120.0`), it is labeled `Fragmentacion_de_Ensamblaje` (assembly fragment) and is **not** counted as an extra gene copy.
   - If cumulative coverage exceeds 120% with short fragments ($< 50\%$), it is labeled `Paralogo_Dudoso` (dubious paralog).

3. **Final CNV count:** The matrix `cnv_matrix.tsv` counts only confirmed functional copies (`Copia_Principal` + `Paralogo_Confirmado`). If the input file includes a `Validacion_Espacial` column, the script also evaluates cross-species transferability and saves `transferability_benchmark.tsv`.

#### Example Usage
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

#### Available Options

| Short | Long | Type | Required | Default | Description |
| :--- | :--- | :--- | :---: | :--- | :--- |
| `-i` | `--input` | File | **YES** | *None* | Input `.tsv`, `.csv`, or `.rds` file with CAPTUS metrics for the samples to analyze. |
| `-m` | `--models_dir` | Folder | No | `models/` | Folder containing the `.rds` model and threshold files. |
| `-o` | `--out_dir` | Folder | No | `results/` | Output folder for CNV tables and plots. |
| `-c` | `--cov_threshold` | Double | No | `50.0` | Minimum individual coverage (`pct_recovered` in %) required to confirm a real paralog. |
| `-f` | `--frag_threshold` | Double | No | `120.0` | Maximum cumulative coverage (%) to treat short contigs as fragments of a single gene copy. |
| `-v` | `--min_votes` | Integer | No | `4` | Minimum positive votes (out of 7) required to accept a sequence as `Homologo_Valido`. |
| `-t` | `--top_loci` | Integer | No | `30` | Number of top loci with the highest copy counts to show in the heatmap. |
| `-h` | `--help` | Flag | No | - | Shows the help menu. |

#### Output Files in `<out_dir>/`
- `cnv_matrix.tsv`: Gene copy number matrix (rows: loci/orthogroups, columns: samples).
- `cnv_annotated_sequences.tsv`: Table of all valid sequences showing cumulative coverage and copy category (`Copia_Principal`, `Paralogo_Confirmado`, `Fragmentacion_de_Ensamblaje`, or `Paralogo_Dudoso`).
- `functional_genes_catalog.tsv`: Filtered list containing only confirmed functional copies (`Copia_Principal` and `Paralogo_Confirmado`).
- `consensus_verdicts.tsv`: Full table showing the individual votes of all 7 models for every evaluated sequence.
- `transferability_benchmark.tsv`: *(Generated automatically when `--input` contains `Validacion_Espacial`)* Classification metrics when applying the models to another species.
- `plot_consensus_distribution.pdf` & `.png`: Bar plot showing accepted versus discarded sequences.
- `plot_cnv_heatmap.pdf` & `.png`: Heatmap of copy numbers across the most variable loci.

### 3.4 Manual Compilation of the Quarto HTML Report

`02_train_models.R` builds the HTML report automatically. To rebuild `reporte_entrenamiento.html` without retraining:

```bash
# In one line with targetcnv.sh:
./targetcnv.sh report

# Or directly with Quarto:
quarto render training_report.qmd \
  -P project_root:"." \
  -P models_dir:"models" \
  -P data_file:"data/data_prepared.rds" \
  --output "reporte_entrenamiento.html"
```

## 4. Workflow 1: Quick Test with Synthetic Data (MWE)

The repository includes a small synthetic dataset in `data/test_data/` so you can verify your installation in a few seconds.

### Isolated Run (Does Not Modify Root `models/` or `results/`)
```bash
./targetcnv.sh test
# Equivalent to: ./run_test_pipeline.sh
```

This command runs all three steps inside `data/test_data/`:
1. `01_prepare_data.R` on `data/test_data/mock_*.gff3`, saving to `data/test_data/prepared/`.
2. `02_train_models.R` using 3-fold cross-validation, saving to `data/test_data/models/` and generating `data/test_data/reporte_entrenamiento_test.html`.
3. `03_cnv_inference.R` on `data/test_data/mock_captus_extract.tsv`, saving the CNV matrix and plots in `data/test_data/results/`.

### Additional Test Options
- `./targetcnv.sh test --no-report`: Runs the test skipping Quarto report compilation for faster checks.
- `./targetcnv.sh test --sync-root`: Copies the test models and results to the root `models/` and `results/` folders. *(Note: this overwrites any existing models in `models/`).*

## 5. Workflow 2: Downloading & Using Pretrained Models from Zenodo

The models trained on the *Solanum torvum* reference genome are archived on Zenodo:
- **DOI:** [`10.5281/zenodo.23463987`](https://doi.org/10.5281/zenodo.23463987)
- **Record URL:** [https://zenodo.org/records/23463987](https://www.google.com/search?q=https://zenodo.org/records/23463987)
- **File:** `TargetCNV_Models.zip` (~382 KB)
- **MD5 Checksum:** `13ec7b35711158079f962a96314eb671`

If you want to run CNV inference on new samples right away without training new models, you have two options:

### Option A: Automatic Download & Inference in 1 Line
```bash
# If models/ is empty, downloads the zip from Zenodo, checks the MD5 hash, and runs inference:
./targetcnv.sh infer -i data/CAPTUS_Muestras/captus-extract_stats.tsv --auto-zenodo
```

### Option B: Step-by-Step
```bash
# Step 1: Download, verify MD5, and unpack the models with targetcnv.sh
./targetcnv.sh zenodo

# (Or manually with curl, md5sum, and unzip)
curl -L -o TargetCNV_Models.zip https://zenodo.org/api/records/23463987/files/TargetCNV_Models.zip/content
md5sum TargetCNV_Models.zip
# Should print: 13ec7b35711158079f962a96314eb671  TargetCNV_Models.zip
unzip -o TargetCNV_Models.zip

# Step 2: Run CNV inference on your samples
./scripts/03_cnv_inference.R \
  --input data/CAPTUS_Muestras/captus-extract_stats.tsv \
  --models_dir models \
  --out_dir results \
  --cov_threshold 50.0 \
  --frag_threshold 120.0 \
  --min_votes 4 \
  --top_loci 30
```

## 6. Workflow 3: Training on New Data & Cross-Species Transferability

### 6.1 Train and Predict on a New Species in 1 Line
If you have an annotated reference genome (`.gff3`), CAPTUS contigs (`.gff`), the `captus-extract_stats.tsv` table, and a gene-to-orthogroup reference table:

```bash
./targetcnv.sh run-all \
  -g path/to/new_genome.gff3 \
  -c path/to/new_captus_contigs.gff \
  -s path/to/captus-extract_stats.tsv \
  -d path/to/orthology_dictionary.xlsx \
  -p ReferenceSampleName \
  -i path/to/population_captus_stats.tsv \
  --out_data data/new_species \
  -m models/new_species \
  -o results/new_species \
  -r training_report_new_species.html
```

### 6.2 Cross-Species Transferability Tests (*S. citrullifolium* & *S. stellatiglandulosum*)

You can also test how well the models trained on *S. torvum* perform when applied to other species with different evolutionary distances or assembly qualities:

#### Case A: Chromosome-Level Assembly (*Solanum citrullifolium*)
```bash
# 1. Prepare S. citrullifolium reference data (automatically harmonizes chromosome prefixes)
./scripts/01_prepare_data.R \
  --gff_ref path/to/Scitrullifolium_genes.gff3 \
  --captus_gff path/to/Scitrullifolium_hit_contigs.gff \
  --captus_stats path/to/captus-extract_stats.tsv \
  --dict_ref data/reference/Genes_Capture_May_2026.xlsx \
  --sample Scitrullifolium \
  --out_dir data/citrullifolium/

# 2. Evaluate the S. torvum models on S. citrullifolium
./scripts/03_cnv_inference.R \
  --input data/citrullifolium/data_prepared.rds \
  --models_dir models/ \
  --out_dir results/citrullifolium/
```

#### Case B: Highly Fragmented Assembly (*Solanum stellatiglandulosum*, HOG $\leftrightarrow$ OG Bridge & `Target` Attribute)
```bash
# 1. Prepare reference data by reading the locus directly from the GFF Target attribute (--target_mode target_attr)
./scripts/01_prepare_data.R \
  --captus_gff path/to/S_stella_hit_contigs.gff \
  --captus_stats path/to/captus-extract_stats.tsv \
  --dict_ref data/reference/Genes_Capture_May_2026.xlsx \
  --sample S_stella \
  --target_mode target_attr \
  --out_dir data/stellatiglandulosum/

# 2. Evaluate transferability and save metrics to results/stellatiglandulosum/transferability_benchmark.tsv
./scripts/03_cnv_inference.R \
  --input data/stellatiglandulosum/data_prepared.rds \
  --models_dir models/ \
  --out_dir results/stellatiglandulosum/
```

## 7. Workflow 4: Reproducing the Reference Experiment (*Solanum torvum*)

### In 1 Single Line
```bash
./targetcnv.sh reproduce-storvum --cores 4
```

### Step-by-Step with Individual Scripts
```bash
# Step 1: Coordinate matching and data preparation (produces 184 NUC sequences: 87 TP and 97 FP)
./scripts/01_prepare_data.R \
  --gff_ref data/reference/Storvum_genes.gff3 \
  --captus_gff data/reference/Storvum_hit_contigs.gff \
  --captus_stats data/reference/captus-extract_stats.tsv \
  --dict_ref data/reference/Genes_Capture_May_2026.xlsx \
  --sample Storvum \
  --out_dir data/

# Step 2: Model training, Youden cutoffs, KernelSHAP, and HTML report
./scripts/02_train_models.R \
  --data data/data_prepared.rds \
  --out_models models \
  --cv_folds 5 \
  --sampling down \
  --cores 4 \
  --seed 1234 \
  --report_file reporte_entrenamiento.html

# Step 3: Multi-sample inference and CNV matrix generation
./scripts/03_cnv_inference.R \
  --input data/CAPTUS_Muestras/captus-extract_stats.tsv \
  --models_dir models \
  --out_dir results \
  --cov_threshold 50.0 \
  --frag_threshold 120.0 \
  --min_votes 4 \
  --top_loci 30
```

### Cross-Validation Results (K=5, Seed 1234)

#### Balanced Scenario (*Down-Sampling*)
| Model | Balanced Acc. (%) | Sensitivity (%) | Specificity (%) | Precision (%) | FDR (%) | F1 Score | MCC | Kappa | ROC-AUC |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Weighted Heuristic** | 88.86 | 83.91 | 93.81 | 92.41 | 7.59 | **0.880** | 0.784 | 0.781 | 0.921 |
| **Neural Network (`nnet`)** | 87.32 | 83.91 | 90.72 | 88.98 | 11.02 | **0.864** | 0.750 | 0.749 | 0.927 |
| **Polynomial SVM (`svmPoly`)** | 87.14 | 80.46 | 93.81 | 92.11 | 7.89 | **0.859** | 0.753 | 0.748 | **0.934** |
| **LDA** | 86.74 | 82.76 | 90.72 | 88.89 | 11.11 | **0.857** | 0.739 | 0.737 | 0.922 |
| **Naive Bayes (`nb`)** | 86.56 | 79.31 | 93.81 | 92.00 | 8.00 | **0.852** | 0.743 | 0.736 | **0.934** |
| **Elastic Net (`glmnet`)** | 85.28 | 83.91 | 86.60 | 84.88 | 15.12 | **0.844** | 0.706 | 0.705 | 0.930 |

#### Unbalanced Scenario (Original Proportions: 87 TP vs. 97 FP)
| Model | Balanced Acc. (%) | Sensitivity (%) | Specificity (%) | Precision (%) | FDR (%) | F1 Score | MCC | Kappa | ROC-AUC |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Weighted Heuristic** | 88.86 | 83.91 | 93.81 | 92.41 | 7.59 | **0.880** | 0.784 | 0.781 | 0.921 |
| **LDA** | 87.26 | 82.76 | 91.75 | 90.00 | 10.00 | **0.862** | 0.750 | 0.748 | 0.922 |
| **Naive Bayes (`nb`)** | 86.56 | 79.31 | 93.81 | 92.00 | 8.00 | **0.852** | 0.743 | 0.736 | **0.934** |
| **Elastic Net (`glmnet`)** | 85.76 | 83.91 | 87.63 | 85.88 | 14.12 | **0.849** | 0.716 | 0.716 | 0.930 |
| **Neural Network (`nnet`)** | 84.60 | 80.46 | 88.66 | 86.42 | 13.58 | **0.833** | 0.695 | 0.694 | 0.921 |
| **Polynomial SVM (`svmPoly`)** | 83.04 | 90.80 | 75.26 | 76.70 | 23.30 | **0.832** | 0.664 | 0.654 | 0.930 |

## 8. Summary of Improvements & Fixes (Preserving Results)

To make the software easier to run without changing any numerical results or biological rules from `Workflow_Central_Sin_Rescate.R`, the following updates were applied:

1. **Unified `targetcnv.sh` script:** Lets you run tests, check dependencies, download and verify Zenodo models via MD5, or execute all three pipeline steps in a single Bash command.
2. **Fixed Excel dictionary reader (`load_reference_dictionary`):** Replaced the invalid `readxl::excel_sheet()` call with `readxl::read_excel()` and included the exact metabolic pathway and assembly filters from `Workflow_Central_Sin_Rescate.R` for *S. torvum*, *S. citrullifolium*, and *S. stellatiglandulosum* (including the `HOG_ID <-> OG_ID` bridge table and `--target_mode target_attr`).
3. **Consistent benchmark tables and plotting functions:** All benchmark tables share the same 16 columns calculated by `compute_classification_metrics()`, sequence IDs use the standard double-underscore format (`OGxxxx__00`) in `clean_captus_ids()`, and `plot_cnv_heatmap()` in `R/plotting.R` accepts both annotated sequence tables and wide count matrices.
4. **Added `training_report.qmd` to the repository root:** Included the interactive Quarto report template (missing from the `main` branch on GitHub) with `DT` tables, `plotly` charts, confusion matrices, and `shapviz` plots, and fixed output file placement when `--report_file` points to a subfolder.
5. **Support for the 18-model tournament, GLM Stacking, and transferability tables:** Added `--extended_tournament` and `--train_stacking` to `02_train_models.R`, and automatic generation of `transferability_benchmark.tsv` in `03_cnv_inference.R` when the input dataset contains reference labels.
