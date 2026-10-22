# TargetCNV-ML: Copy Number Variation Detection from Target Capture Data

**Languages / Idiomas:** [English (README.md)](README.md) | [Español (README\_es.md)](https://www.google.com/search?q=README_es.md)  
**Execution Guides:** [Execution Guide (EXECUTION\_GUIDE.md)](https://www.google.com/search?q=EXECUTION_GUIDE.md) | [Guía de Ejecución (GUIA\_EJECUCION.md)](https://www.google.com/search?q=GUIA_EJECUCION.md)

[](https://doi.org/10.5281/zenodo.23463987)





**TargetCNV-ML** is a command-line (CLI) bioinformatics pipeline written in R and Bash to estimate **Copy Number Variations (CNV)** and identify **gene family expansions** from target capture (*Hyb-Seq* / *Target Enrichment*) sequencing assemblies.

## Why TargetCNV-ML?

Target capture is widely used in phylogenomics, evolutionary biology, and plant breeding because it sequences specific regions of interest at a fraction of whole-genome sequencing costs. However, counting gene copies from *de novo* target capture assemblies presents two practical problems:

1. **Off-target and cross-homolog alignments:** Capture probes can hybridize with closely related paralogs, pseudogenes, or flanking non-coding regions. As a result, assembled contigs often align to multiple similar loci, increasing false positives.
2. **Assembly fragmentation:** *De novo* assemblers frequently break a single long gene into two or more short contigs. Simple contig counting mistakes these fragments for extra gene copies.

**TargetCNV-ML** addresses both issues in three steps:
- **Genomic coordinate matching (`GenomicRanges`):** Compares assembled contigs (`protein_match:NUC`) against a reference genome annotation (`.gff3`) and a reference orthology table to label true gene hits (`TRUE_POSITIVE`) versus off-target or redundant hits (`FALSE_POSITIVE`).
- **7-Model Consensus Classification:** Combines **4 supervised machine learning models** (`ElasticNet`, `NaiveBayes`, `SVM_Polinomial`, `Red_Neuronal`), **2 univariate Youden index cutoffs** (`pct_recovered`, `pct_identity`), and **1 weighted heuristic score**. A sequence is accepted as a valid homolog if at least 4 of the 7 methods vote positive.
- **Coverage-based fragmentation filter:** Separates complete additional gene copies (`Paralogo_Confirmado`, individual coverage $\ge 50\%$) from fragmented contigs of the same gene (`Fragmentacion_de_Ensamblaje`, individual coverage $< 50\%$ and cumulative locus coverage $\le 120\%$).

## Repository Structure

```text
TargetCNV-ML/
├── targetcnv.sh               # Main Bash helper script to run workflows in one line
├── run_test_pipeline.sh       # Test script using the small synthetic dataset (MWE)
├── training_report.qmd        # Quarto template for the interactive HTML training report
├── README.md                  # Main overview in English
├── README_es.md               # Main overview in Spanish
├── EXECUTION_GUIDE.md         # Execution Guide in English
├── GUIA_EJECUCION.md          # Execution Guide in Spanish
├── renv.lock                  # Exact R package versions used by the pipeline
├── .Rprofile                  # Activates renv automatically when starting R
├── LICENSE                    # MIT License
│
├── scripts/                   # Command-line R scripts
│   ├── 01_prepare_data.R      # Step 1: Coordinate matching, reference table join & local deltas
│   ├── 02_train_models.R      # Step 2: Threshold calculation, ML training, SHAP & HTML report
│   └── 03_cnv_inference.R     # Step 3: 7-model consensus voting & CNV matrix generation
│
├── R/                         # Internal R functions
│   ├── data_processing.R      # GFF filtering, spatial overlaps, dictionary loading & CNV rules
│   ├── model_evaluation.R     # Classification metrics, Youden cutoffs & consensus voting
│   └── plotting.R             # Plots for densities, confusion matrices, benchmarks & heatmaps
│
├── data/                      # Input files and test data
│   ├── reference/             # S. torvum reference GFF3, CAPTUS contigs & reference gene list
│   ├── CAPTUS_Muestras/       # CAPTUS alignment statistics for biological samples
│   ├── test_data/             # Small synthetic dataset for quick testing (< 25 KB)
│   └── data_prepared.rds      # Prepared S. torvum training dataset (184 NUC sequences)
│
├── models/                    # Saved R models (.rds) and performance tables (.tsv)
└── results/                   # Output CNV matrices (.tsv), sequence catalogs & plots (.pdf/.png)
```

## Requirements & Installation

### System Requirements
- **Operating System:** Linux, macOS, or Windows (via WSL2).
- **R:** $\ge 4.0.0$ (R 4.3+ or 4.6 recommended).
- **Quarto CLI:** $\ge 1.3+$ (optional, only required to generate the HTML report `reporte_entrenamiento.html`).
- **Standard command-line tools:** `bash`, `curl` or `wget`, `unzip`, `md5sum`.

### Installation Steps

```bash
# 1. Clone the repository
git clone https://github.com/xyzdiego/TargetCNV-ML.git
cd TargetCNV-ML

# 2. Restore the R package environment from renv.lock
Rscript -e "renv::restore()"

# 3. Make scripts executable, check dependencies, and run the test dataset
chmod +x targetcnv.sh run_test_pipeline.sh scripts/*.R
./targetcnv.sh doctor
./targetcnv.sh test
```

## Quick Execution in One Line (`./targetcnv.sh`)

You can run any workflow with a single command using `./targetcnv.sh`:

| Task | Command |
| :--- | :--- |
| **Check installed tools and R packages** | `./targetcnv.sh doctor` |
| **Run the test dataset (MWE)** | `./targetcnv.sh test` |
| **Download pretrained models from Zenodo** | `./targetcnv.sh zenodo` |
| **Predict CNVs on new samples (auto-downloading models if needed)** | `./targetcnv.sh infer -i data/CAPTUS_Muestras/captus-extract_stats.tsv --auto-zenodo` |
| **Reproduce the *Solanum torvum* experiment** | `./targetcnv.sh reproduce-storvum --cores 4` |
| **Run all 3 steps on a new species** | `./targetcnv.sh run-all -g ref.gff3 -c captus.gff -s stats.tsv -d dict.xlsx -p SampleName -i pop_stats.tsv` |
| **Rebuild the HTML training report** | `./targetcnv.sh report` |

## Step-by-Step CLI Workflow

If you prefer to run each module separately to customize parameters:

### Step 1: Data Preparation (`01_prepare_data.R`)
Finds coordinate overlaps between CAPTUS contigs (`protein_match:NUC`) and reference genes (`gene`), checks matches against the reference orthology dictionary, marks secondary hits in the same gene as `FALSE_POSITIVE`, and computes within-locus differences (`delta_wscore`, `delta_pct_identity`).

```bash
./scripts/01_prepare_data.R \
  --gff_ref data/reference/Storvum_genes.gff3 \
  --captus_gff data/reference/Storvum_hit_contigs.gff \
  --captus_stats data/reference/captus-extract_stats.tsv \
  --dict_ref data/reference/Genes_Capture_May_2026.xlsx \
  --sample Storvum \
  --out_dir data/
```

### Step 2: Model Training & Evaluation (`02_train_models.R`)
Calculates optimal Youden index cutoffs for each alignment metric, finds the best weights for the heuristic combination, trains the machine learning classifiers using 5-fold cross-validation (both balanced via down-sampling and unbalanced), computes KernelSHAP variable importance, and generates the HTML report.

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
*(Optional flags: `--extended_tournament` trains all 18 classifiers from the original study; `--train_stacking` fits a GLM stacking model on top of the best classifiers).*

### Step 3: Consensus Prediction & CNV Matrix (`03_cnv_inference.R`)
Applies the 7-model consensus vote to new samples, filters out assembly fragments, and exports the copy number matrix (`Locus x Samples`), the catalog of functional copies, and summary plots.

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

## Reference Results (*Solanum torvum*)

Results on 184 `NUC` sequences from *Solanum torvum* (87 True Positives and 97 False Positives) using 5-fold cross-validation (`seed = 1234`, balanced with down-sampling):

| Model / Method | Balanced Acc. (%) | Sensitivity (%) | Specificity (%) | Precision (%) | FDR (%) | F1-Score | MCC | ROC-AUC |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Weighted Heuristic** | 88.86 | 83.91 | 93.81 | 92.41 | 7.59 | **0.880** | 0.784 | 0.921 |
| **Neural Network (`nnet`)** | 87.32 | 83.91 | 90.72 | 88.98 | 11.02 | **0.864** | 0.750 | 0.927 |
| **Polynomial SVM (`svmPoly`)** | 87.14 | 80.46 | 93.81 | 92.11 | 7.89 | **0.859** | 0.753 | **0.934** |
| **LDA** | 86.74 | 82.76 | 90.72 | 88.89 | 11.11 | **0.857** | 0.739 | 0.922 |
| **Naive Bayes (`nb`)** | 86.56 | 79.31 | 93.81 | 92.00 | 8.00 | **0.852** | 0.743 | **0.934** |
| **Elastic Net (`glmnet`)** | 85.28 | 83.91 | 86.60 | 84.88 | 15.12 | **0.844** | 0.706 | 0.930 |

## Pretrained Models on Zenodo

The models trained on *Solanum torvum* are available on Zenodo:
- **DOI:** [`10.5281/zenodo.23463987`](https://doi.org/10.5281/zenodo.23463987)
- **File:** `TargetCNV_Models.zip` (`MD5: 13ec7b35711158079f962a96314eb671`)
- **Download command:** Run `./targetcnv.sh zenodo` to download the file, check its MD5 hash, and unpack the models into `models/`.

For full details on all flags and workflows, see the [Execution Guide (EXECUTION\_GUIDE.md)](https://www.google.com/search?q=EXECUTION_GUIDE.md).

## Citation & License

If you use **TargetCNV-ML** in your work, please cite:

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

Released under the **MIT License** — see [LICENSE](https://www.google.com/search?q=LICENSE) for details.
