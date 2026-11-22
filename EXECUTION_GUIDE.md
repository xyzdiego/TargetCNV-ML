# TargetCNV-ML Execution Guide

[Guía de Ejecución en Español](GUIA_EJECUCION.md) | [English README](README.md) | [README en Español](README_ES.md)

This guide covers environment setup, CLI arguments, Zenodo model management, and experimental workflows for TargetCNV-ML. Console logs, TSV column names, and the rendered Quarto reports in `reports/` (`reporte_entrenamiento.html`, `reporte_inferencia.html`, `reporte_transferencia.html`; bundled in `reports/reportes_html.zip`) are produced in Spanish to match the Master's thesis.

## 1. Environment Setup

### Requirements
- Linux (tested on Ubuntu 22.04/24.04, Debian 12+, Arch/CachyOS), macOS, or Windows with WSL2.
- **R 4.6.2** with **Bioconductor 3.23** for exact `renv.lock` parity (compatible with **R $\ge 4.5$ / Bioconductor $\ge 3.22$** due to `Seqinfo`).
- **Quarto CLI** (`1.10.20` or $\ge 1.4$) to render HTML reports.
- System utilities and compilers: `bash`, `curl` or `wget`, `unzip`, `md5sum`, `gcc`, `g++`, `gfortran`, `cmake`, `make`.

### 1.1 Restoring with `renv`

```bash
git clone https://github.com/xyzdiego/TargetCNV-ML.git
cd TargetCNV-ML
chmod +x targetcnv.sh run_test_pipeline.sh scripts/*.R

./targetcnv.sh restore
./targetcnv.sh doctor
```

**Lockfile integrity:** Always install with `./targetcnv.sh restore` (or `./targetcnv.sh setup`), which runs `renv::restore(prompt = FALSE)` against `renv.lock`. Do not call `renv::hydrate()` or `renv::snapshot()` during user setup, because hydrating pulls unpinned host packages and snapshots overwrite `renv.lock`.

### 1.2 Pixi, Apptainer, and R Package Installation

```bash
# Pixi (conda-forge + bioconda + CRAN extras)
pixi install
pixi run setup-cran
pixi run doctor

# Apptainer / Singularity container (bundles R 4.6.2 + Quarto 1.10.20)
apptainer build targetcnv-ml.sif Apptainer.def
apptainer run targetcnv-ml.sif doctor
```

The repository is structured as a valid R package (`DESCRIPTION`, `NAMESPACE`, `man/`), so the 60 exported functions under `R/` can be loaded either by `scripts/*.R` via `source()` or directly with `pkgload::load_all(".")` / `R CMD INSTALL .`.

## 2. Read-Only Reporting Design

Computation lives exclusively in `scripts/01_prepare_data.R` through `scripts/04_evaluate_transferability.R`. The three Quarto documents (`training_report.qmd`, `report_template.qmd`, `transferability_report.qmd`) only read precomputed `.rds` and `.tsv` files from `models/` and `results/`, styled with `assets/report.css` and `assets/waterfall.js`. If `models/` is empty, the training and transferability reports print a notice pointing to `./targetcnv.sh train` instead of fitting models on the fly.

## 3. Command-Line Reference

### 3.0 Wrapper Script (`./targetcnv.sh`)

| Subcommand | Example | Action |
| :--- | :--- | :--- |
| `restore` / `setup` | `./targetcnv.sh restore` | Restores R packages from `renv.lock`. |
| `doctor` | `./targetcnv.sh doctor` | Checks `Rscript`, `quarto`, `curl`, `unzip`, `md5sum`, and R packages. |
| `test` | `./targetcnv.sh test [--no-report]` | Runs the MWE on `data/test_data/`. |
| `zenodo` | `./targetcnv.sh zenodo [--out-dir .]` | Fetches `TargetCNV_Models.zip`, verifies MD5, and checks extracted files. |
| `prepare` | `./targetcnv.sh prepare [01 args]` | Calls `scripts/01_prepare_data.R`. |
| `train` | `./targetcnv.sh train [02 args]` | Calls `scripts/02_train_models.R`. |
| `infer` | `./targetcnv.sh infer -i stats.tsv --auto-zenodo` | Calls `scripts/03_cnv_inference.R` (downloads models if missing). |
| `transfer` | `./targetcnv.sh transfer [--sources ...]` | Calls `scripts/04_evaluate_transferability.R` and compiles the transfer report. |
| `report` | `./targetcnv.sh report [--models-dir models]` | Re-renders `reports/reporte_entrenamiento.html` from saved artifacts. |
| `reproduce-storvum` | `./targetcnv.sh reproduce-storvum --cores 4` | Runs `01` -> `02` -> `03` on the *Solanum torvum* dataset. |
| `run-all` | `./targetcnv.sh run-all -g ref.gff3 -c cap.gff -s stats.tsv -d dict.xlsx -i new.tsv` | Runs `01` -> `02` -> `03` on a new dataset with separate HTML reports. |

**Flags that delete or overwrite outputs:**
- `--clean_models` (`02_train_models.R`) removes all existing `.rds` and `.tsv` files in `--out_models` before writing new ones. It is disabled by default and is not called by `reproduce-storvum`. Use it before a full re-train; skip it with `--fast_tournament` if you want to keep previously trained 18-model artifacts.
- `--sync-root` (`./targetcnv.sh test` / `./run_test_pipeline.sh`) copies synthetic MWE outputs over root `models/` and `results/`. Avoid it whenever real or Zenodo models are stored in `models/`.

### 3.1 Module 01: Spatial Labeling (`01_prepare_data.R`)

Labels each nuclear CAPTUS hit (`protein_match:NUC`) as `TRUE_POSITIVE` or `FALSE_POSITIVE`:
- **Chromosome scale (`--genome_scale chromosome` / `--target_mode overlap`):** harmonizes sequence prefixes and intersects coordinates and strand (`GenomicRanges::findOverlaps`, `ignore.strand = FALSE`) against `gene` features in `--gff_ref`.
- **Contig scale (`--genome_scale contig` / `--target_mode target_attr`):** extracts reference gene IDs from the CAPTUS GFF `Target` attribute and maps hierarchical `HOG_ID` groups to `OG_ID` orthogroups via the probe-design bridge table.
- **Intra-locus deduplication:** within each sample and reference locus, `resolve_redundancies()` keeps only the positive hit with greatest `ref_len_matched` and `pct_recovered` as `TRUE_POSITIVE`; extra hits on that same gene become `FALSE_POSITIVE` (`5_TP_Redundante_Degradado`).

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

| Short | Long | Type | Required | Default | Description |
| :--- | :--- | :--- | :---: | :--- | :--- |
| `-g` | `--gff_ref` | File | Yes\* | `NULL` | Reference genome GFF/GFF3 `gene` file (\*required in `chromosome` mode). |
| `-c` | `--captus_gff` | File | **Yes** | `NULL` | CAPTUS assembled contigs GFF (`protein_match:NUC`). |
| `-s` | `--captus_stats` | File | **Yes** | `NULL` | CAPTUS `captus-extract_stats.tsv` alignment metrics table. |
| `-d` | `--dict_ref` | File | No | `NULL` | Homology dictionary (`.xlsx`, `.tsv`, or `.csv`). |
| `-p` | `--sample` | String | No | `NULL` | Sample name filter inside `captus-extract_stats.tsv` (`Storvum`, `Scitrullifolium`, `S_stella`). |
| | `--genome_scale` | String | No | `chromosome` | `'chromosome'` (GFF3 overlap) or `'contig'` (`Target` + `OG_ID`/`HOG_ID` bridge). |
| | `--target_mode` | String | No | `overlap` | Explicit mode alias: `'overlap'` or `'target_attr'`. |
| | `--allow_multi_sample` | Flag | No | `FALSE` | Allows multi-sample input without `--sample` (intended for MWE synthetic data). |
| `-o` | `--out_dir` | Dir | No | `data/` | Output directory for `data_prepared.rds` and `data_prepared.tsv`. |

### 3.2 Module 02: Collinearity, Deterministic CV, and SHAP (`02_train_models.R`)

Benchmarks 24 methods over `pct_recovered`, `pct_identity`, `score`, and `wscore`:
1. **Spearman collinearity:** exports `models/spearman_collinearity.tsv` and `.rds` (in *S. torvum*, `score` vs. `wscore` $\rho = 0.9977$; `pct_recovered` vs. `pct_identity` $\rho = 0.7776$).
2. **Youden cutoffs (4) and `Weighted_Heuristic`:** maximizes $J = \text{Sensitivity} + \text{Specificity} - 1$ per metric and searches integer weights (`0:3`) for the highest-AUC combination.
3. **Deterministic CV (18 ML classifiers + `Ensamble_Stacking`):** constructs explicit per-resample `seeds` in `caret::trainControl()` from `--seed` (identical across `--cores`) under **Balanced** (`down`-sampling) and **Unbalanced** regimes, followed by KernelSHAP calculation and rendering of `reports/reporte_entrenamiento.html`.

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

| Short | Long | Type | Required | Default | Description |
| :--- | :--- | :--- | :---: | :--- | :--- |
| `-d` | `--data` | File | **Yes** | `data/data_prepared.rds` | Labeled dataset from Module 01. |
| `-m` | `--out_models` | Dir | No | `models/` | Output folder for serialized `.rds` models, benchmarks, and SHAP objects. |
| `-k` | `--cv_folds` | Integer | No | `5` | Number of stratified CV folds ($\ge 2$). |
| `-s` | `--sampling` | String | No | `down` | Balancing method for the balanced scenario (`'down'`, `'up'`, `'none'`). |
| `-c` | `--cores` | Integer | No | `detectCores() - 1` | PSOCK worker count. |
| | `--seed` | Integer | No | `1234` | Master seed controlling folds, tuning grids, and explicit `seeds`. |
| | `--extended_tournament` / `--fast_tournament` | Flag | No | `TRUE` | Train all 18 ML classifiers (default) or the 5-model quick subset. |
| | `--train_stacking` / `--no_stacking` | Flag | No | `TRUE` | Fit (default) or skip `Ensamble_Stacking`. |
| | `--clean_models` | Flag | No | `FALSE` | Destructive: deletes existing `.rds`/`.tsv` files in `--out_models` first. |
| | `--include_hit` | Flag | No | `FALSE` | Adds CAPTUS `hit` rank as a fifth predictor (recorded in `consensus_spec.rds`; default `FALSE` is recommended across datasets). |
| | `--render_report` / `--no_render_report` | Flag | No | `TRUE` | Compile (default) or skip `reports/reporte_entrenamiento.html`. |
| `-r` | `--report_file` | File | No | `reports/reporte_entrenamiento.html` | Output HTML path. |
| | `--sources` | File | No | `""` | Optional TSV with genome/literature citations (`data/reference/sources_citations.tsv`). |

Benchmark tables contain 16 columns (`Modelo`, `Exactitud`, `Exactitud_Balanceada`, `Sensibilidad`, `Especificidad`, `Precision`, `VPN`, `FDR`, `F1_Score`, `MCC`, `Kappa`, `ROC_AUC`, `TP`, `FP`, `FN`, `TN`). Models tagged with `*` in the HTML tables reach a higher FDR-weighted composite index in the alternative regime:

$$\text{Composite Score}_{\text{FDR}} = 0.45(1 - \text{FDR}) + 0.20\,\text{F1} + 0.15\,\text{MCC} + 0.10\,\text{Balanced Acc.} + 0.10\,\text{ROC-AUC}$$

### 3.3 Module 03: CNV Inference and Coverage Sensitivity (`03_cnv_inference.R`)

Applies 7-model consensus voting ($\ge 4/7$ votes, dynamically reading active voters from `consensus_spec.rds` and falling back to `ElasticNet`, `NaiveBayes`, `SVM_Polinomial`, `Red_Neuronal`, `Youden_pct_recovered`, `Youden_pct_identity`, and `Weighted_Heuristic`). Valid sequences (`Homologo_Valido`) are sorted by descending `pct_recovered` within each sample and orthogroup:
- `Copia_Principal`: top-coverage valid sequence in that orthogroup and sample.
- `Homologo_Confirmado`: secondary valid hit with `pct_recovered >= --cov_threshold` (default `50.0%`).
- `No_Contado`: secondary valid hit below `--cov_threshold`, kept in the annotated table but excluded from copy counts.

The script also sweeps `--cov_threshold` from `30%` to `80%` (`cnv_threshold_sensitivity.tsv`). Because `captus-extract_stats.tsv` provides scalar recovery percentages rather than per-fragment alignments along the reference protein, this cutoff can under-count true copies recovered below `50%` or count two non-overlapping fragments if both individually exceed `50%`.

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

| Short | Long | Type | Required | Default | Description |
| :--- | :--- | :--- | :---: | :--- | :--- |
| `-i` | `--input` | File | **Yes** | `NULL` | `.tsv`, `.csv`, or `.rds` CAPTUS extraction table (104 samples in the reference dataset). |
| `-m` | `--models_dir` | Dir | No | `models/` | Folder holding trained `.rds` artifacts. |
| `-o` | `--out_dir` | Dir | No | `results/` | Destination folder for matrices, sensitivity tables, and plots. |
| `-c` | `--cov_threshold` | Float | No | `50.0` | Minimum `pct_recovered` (%) to count secondary hits as `Homologo_Confirmado`. |
| `-v` | `--min_votes` | Integer | No | `4` | Minimum votes (out of 7) for `Homologo_Valido`. |
| `-t` | `--top_loci` | Integer | No | `30` | Number of top orthogroups shown in the CNV heatmap. |
| | `--true_cnv` | File | No | `""` | Optional reference copy-number table ($\text{CNV}^{\text{real}}_{\text{OG}}$) for concordance metrics. |
| | `--render_report` / `--no_render_report` | Flag | No | `TRUE` | Compile (default) or skip `reports/reporte_inferencia.html`. |
| `-r` | `--report_file` | File | No | `reports/reporte_inferencia.html` | Output HTML path. |
| | `--sources` | File | No | `""` | Optional genome citations TSV (`data/reference/sources_citations.tsv`). |

Main outputs written to `<out_dir>/`: `cnv_matrix.tsv`, `cnv_annotated_sequences.tsv`, `functional_genes_catalog.tsv`, `consensus_verdicts.tsv`, `cnv_threshold_sensitivity.tsv`, `cnv_concordance_summary.tsv` (when labels or `--true_cnv` exist), PDF/PNG figures, and `reports/reporte_inferencia.html`.

### 3.4 Module 04: Cross-Species Transferability (`04_evaluate_transferability.R`)

Evaluates all 24 pre-trained *S. torvum* methods and the 7-model consensus without re-training on recipient species (*S. citrullifolium* at chromosome scale and *S. stellatiglandulosum* at contig scale), calculates orthogroup copy-number concordance ($\widehat{\text{CNV}}_{\text{OG}}$ vs. $\text{CNV}^{\text{real}}_{\text{OG}}$), sweeps `--cov_threshold` across `30%` to `80%`, exports `results/transferability_*.tsv|.rds`, and renders `reports/reporte_transferencia.html`.

```bash
./targetcnv.sh transfer --sources data/reference/sources_citations.tsv
```

## 4. Test Suite and Minimal Working Example

```bash
Rscript tests/testthat.R
./targetcnv.sh test
```

## 5. Pre-Trained Models on Zenodo

Pre-trained *S. torvum* artifacts are stored at [https://doi.org/10.5281/zenodo.23463987](https://doi.org/10.5281/zenodo.23463987) (`TargetCNV_Models.zip`).

```bash
./targetcnv.sh zenodo
./targetcnv.sh infer -i data/CAPTUS_Muestras/captus-extract_stats.tsv --auto-zenodo
```

To update the Zenodo archive after a full re-train:
1. Run `zip -r TargetCNV_Models.zip models/ && md5sum TargetCNV_Models.zip`.
2. Create a new version under DOI `10.5281/zenodo.23463987` and upload the updated zip.
3. Update `ZENODO_URL` and `ZENODO_MD5` at the top of `targetcnv.sh`.

## 6. New Species and Cross-Species Workflows

### 6.1 End-to-end `run-all` on a New Species

```bash
./targetcnv.sh run-all \
  -g path/to/new_genome.gff3 \
  -c path/to/new_captus_contigs.gff \
  -s path/to/captus-extract_stats.tsv \
  -d path/to/homology_dictionary.xlsx \
  -p ReferenceSample \
  -i path/to/population_captus_stats.tsv \
  --out_data data/new_species \
  -m models/new_species \
  -o results/new_species \
  -r reports/training_report_new_species.html \
  --infer_report_file reports/inference_report_new_species.html
```

### 6.2 *S. citrullifolium* and *S. stellatiglandulosum* Transferability

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

#### Orthogroup Copy-Number Concordance (`--cov_threshold 50.0%`, 7-Model Consensus)

| Species | Genome Scale | Orthogroups ($N$) | Reference Copies | Estimated Copies | Exact Concordance (%) | Spearman's $\rho$ | Quadratic $\kappa_w$ | RMSE |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| ***S. torvum*** (Donor, in-sample) | Chromosome Scale | 58 | 87 | 79 | **79.31** | **0.5983** ($p = 7.0 \times 10^{-7}$) | **0.5269** | 0.8710 |
| ***S. citrullifolium*** (Transfer) | Chromosome Scale | 55 | 56 | 66 | **70.91** | **0.4653** ($p = 3.4 \times 10^{-4}$) | **0.5506** | 0.5394 |
| ***S. stellatiglandulosum*** (Transfer) | Contig Scale | 58 | 19 | 78 | 24.14 | 0.0831 ($p = 0.535$) | 0.0742 | 1.3834 |

## 7. Reproducing the *Solanum torvum* Experiment

```bash
./targetcnv.sh reproduce-storvum --cores 4
```

### Representative Cross-Validation Benchmarks ($K = 5$, `--seed 1234`, $N = 184$: 87 TP vs. 97 FP)

#### Balanced Scenario (*Down-Sampling*)

| Model | Balanced Acc. (%) | Sensitivity (%) | Specificity (%) | Precision (%) | FDR (%) | F1 Score | MCC | Kappa | ROC-AUC |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Weighted_Heuristic** | 88.86 | 83.91 | 93.81 | 92.41 | 7.59 | **0.880** | 0.784 | 0.781 | 0.921 |
| **Ensamble_Stacking** | 88.29 | 82.76 | 93.81 | 92.31 | 7.69 | **0.873** | 0.774 | 0.770 | **0.932** |
| **Red_Neuronal (`nnet`)** | 87.31 | 83.91 | 90.72 | 89.02 | 10.98 | **0.864** | 0.750 | 0.749 | 0.927 |
| **SVM_Polinomial (`svmPoly`)** | 87.14 | 80.46 | 93.81 | 92.11 | 7.89 | **0.859** | 0.753 | 0.748 | **0.934** |
| **LDA \*** | 86.74 | 82.76 | 90.72 | 88.89 | 11.11 | **0.857** | 0.739 | 0.737 | 0.922 |
| **NaiveBayes (`nb`)** | 86.56 | 79.31 | 93.81 | 92.00 | 8.00 | **0.852** | 0.743 | 0.736 | **0.934** |
| **ElasticNet (`glmnet`) \*** | 85.25 | 83.91 | 86.60 | 84.88 | 15.12 | **0.844** | 0.706 | 0.705 | 0.930 |
| **XGBoost \*** | 84.50 | 79.31 | 89.69 | 87.34 | 12.66 | **0.831** | 0.696 | 0.693 | 0.920 |

#### Unbalanced Scenario (Original Class Counts: 87 TP vs. 97 FP)

| Model | Balanced Acc. (%) | Sensitivity (%) | Specificity (%) | Precision (%) | FDR (%) | F1 Score | MCC | Kappa | ROC-AUC |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Weighted_Heuristic** | 88.86 | 83.91 | 93.81 | 92.41 | 7.59 | **0.880** | 0.784 | 0.781 | 0.921 |
| **Ensamble_Stacking** | 88.35 | 83.91 | 92.78 | 91.25 | 8.75 | **0.874** | 0.772 | 0.770 | **0.935** |
| **LDA** | 87.26 | 82.76 | 91.75 | 90.00 | 10.00 | **0.862** | 0.750 | 0.748 | 0.922 |
| **NaiveBayes (`nb`)** | 86.56 | 79.31 | 93.81 | 92.00 | 8.00 | **0.852** | 0.743 | 0.736 | **0.934** |
| **ElasticNet (`glmnet`)** | 85.77 | 83.91 | 87.63 | 85.88 | 14.12 | **0.849** | 0.716 | 0.716 | 0.930 |
| **Red_Neuronal (`nnet`) \*** | 84.56 | 80.46 | 88.66 | 86.42 | 13.58 | **0.833** | 0.695 | 0.694 | 0.921 |
| **XGBoost** | 84.56 | 80.46 | 88.66 | 86.42 | 13.58 | **0.833** | 0.695 | 0.694 | 0.918 |
| **SVM_Polinomial (`svmPoly`) \*** | 83.03 | 90.80 | 75.26 | 76.70 | 23.30 | **0.832** | 0.664 | 0.654 | 0.930 |

## 8. Data Sources and License

- **Chromosome-scale genomes (*S. torvum*, *S. citrullifolium*):** Solanum Pan-Genomics (https://solpangenomics.com/dist/pages/downloads/index.php), SGN (https://solgenomics.net/ftp/genomes/Solanum_pangenomics/), NCBI BioProject `PRJNA1073673` (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1073673); Benoit, M. *et al.* (2025), *Nature*, 640(8057): 135-145, [doi:10.1038/s41586-025-08619-6](https://doi.org/10.1038/s41586-025-08619-6).
- **Contig-scale genome (*S. stellatiglandulosum*):** NCBI BioProject `PRJNA1070281` (https://www.ncbi.nlm.nih.gov/bioproject/?term=txid3132602%5BOrganism:noexp%5D); Pérez-Mesa, P. A. & Roda, F. (2025), *Current Opinion in Plant Biology*, 85: 102727, [doi:10.1016/j.pbi.2025.102727](https://doi.org/10.1016/j.pbi.2025.102727); Pérez-Mesa, P. A., Pardo, R. A., Alseekh, S. *et al.* (2026), *Plant Molecular Biology*, 116: 67, [doi:10.1007/s11103-026-01718-6](https://doi.org/10.1007/s11103-026-01718-6).

Source code is licensed under the MIT License ([LICENSE](LICENSE)). Derived datasets and pre-trained models on Zenodo are licensed under CC BY 4.0.
