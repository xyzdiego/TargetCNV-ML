#!/usr/bin/env bash

# Wrapper CLI para TargetCNV-ML (preparación, entrenamiento, inferencia, transferencia y entorno).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

ZENODO_URL="https://zenodo.org/api/records/23463987/files/TargetCNV_Models.zip/content"
ZENODO_ZIP="TargetCNV_Models.zip"
ZENODO_MD5="13ec7b35711158079f962a96314eb671"

resolve_default_sources_file() {
  if [ -f "data/reference/sources_citations.tsv" ]; then
    echo "data/reference/sources_citations.tsv"
  elif [ -f "$SCRIPT_DIR/data/reference/sources_citations.tsv" ]; then
    echo "$SCRIPT_DIR/data/reference/sources_citations.tsv"
  elif [ -f "sources_citations.tsv" ]; then
    echo "sources_citations.tsv"
  else
    echo "$SCRIPT_DIR/data/reference/sources_citations.tsv"
  fi
}

resolve_abs_path() {
  local target="$1"
  local dir_part base_part
  if [ -d "$target" ]; then
    (cd "$target" && pwd -P)
  else
    dir_part="$(dirname "$target")"
    base_part="$(basename "$target")"
    if [ -d "$dir_part" ]; then
      printf '%s/%s\n' "$(cd "$dir_part" && pwd -P)" "$base_part"
    else
      case "$target" in
        /*) printf '%s\n' "$target" ;;
        *)  printf '%s/%s\n' "$(pwd -P)" "$target" ;;
      esac
    fi
  fi
}

show_help() {
  cat <<'EOF'
TargetCNV-ML - Interfaz de línea de comandos

Uso:
  ./targetcnv.sh <subcomando> [opciones]

Subcomandos:
  restore (alias: setup)
      Restaura la librería de R exactamente desde renv.lock (renv::restore).
  doctor
      Comprueba R, paquetes requeridos, Quarto y utilidades del sistema.
  test [--sync-root] [--no-report]
      Ejecuta el Ejemplo Mínimo Funcional (MWE) sobre los datos en data/test_data/.
      --sync-root sobrescribe models/ y results/ en la raíz del proyecto.
  zenodo [--out-dir DIR]
      Descarga los modelos preentrenados (DOI: 10.5281/zenodo.23463987) y valida MD5.
  prepare [--genome_scale chromosome|contig] [opciones de 01_prepare_data.R]
      Ejecuta el módulo 01 (intersección espacial GFF o atributo Target).
  train [opciones de 02_train_models.R]
      Ejecuta el módulo 02 (colinealidad Spearman, CV determinista, SHAP y reporte).
  infer -i <captus_stats.tsv> [--auto-zenodo] [opciones de 03_cnv_inference.R]
      Ejecuta el módulo 03 (consenso de 7 modelos, matriz CNV y curva 30%-80%).
  transfer [--models-dir DIR] [--cit-data FILE] [--ste-data FILE] [--output FILE]
      Ejecuta el módulo 04 (transferibilidad inter-especie y concordancia CNV).
  reproduce-storvum [--cores N]
      Ejecuta secuencialmente los módulos 01 -> 02 -> 03 sobre Solanum torvum.
  run-all -g <ref.gff3> -c <captus.gff> -s <stats.tsv> -d <dict> [-i <infer.tsv>]
      Ejecuta los módulos 01, 02 y 03 sobre un conjunto de datos nuevo.
  report [--models-dir DIR] [--data FILE] [--output FILE]
      Compila training_report.qmd en modo de solo lectura.

Ejemplos:
  ./targetcnv.sh restore
  ./targetcnv.sh doctor
  ./targetcnv.sh test
  ./targetcnv.sh infer -i data/CAPTUS_Muestras/captus-extract_stats.tsv --auto-zenodo
  ./targetcnv.sh transfer
  ./targetcnv.sh reproduce-storvum --cores 4
EOF
}

cmd_restore() {
  echo "[restore] Restaurando entorno desde $SCRIPT_DIR/renv.lock..."
  if [ ! -f "$SCRIPT_DIR/renv.lock" ]; then
    echo "[ERROR] No se encontró renv.lock en $SCRIPT_DIR." >&2
    exit 1
  fi
  TARGETCNV_PROJ="$SCRIPT_DIR" Rscript -e '
    if (!requireNamespace("renv", quietly = TRUE)) {
      stop("El paquete renv no está disponible. Instale renv antes de ejecutar restore.", call. = FALSE)
    }
    renv::restore(project = Sys.getenv("TARGETCNV_PROJ"), prompt = FALSE)
  '
  echo "[restore] Librería restaurada conforme a renv.lock."
}

cmd_dev_hydrate() {
  local update_lock=0
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --update-lock) update_lock=1; shift ;;
      -h|--help)
        echo "Uso interno: ./targetcnv.sh dev-hydrate [--update-lock]"
        echo "Comando reservado para mantenimiento del entorno; para instalar use './targetcnv.sh restore'."
        exit 0 ;;
      *) echo "[ERROR] Opción desconocida en dev-hydrate: $1" >&2; exit 1 ;;
    esac
  done

  echo "[dev-hydrate] Hidratación local desde librerías del sistema (solo desarrollo)."
  local user_r_libs
  user_r_libs="$(Rscript --vanilla -e 'cat(paste(.libPaths(), collapse = ":"))')"
  TARGETCNV_PROJ="$SCRIPT_DIR" USER_R_LIBS="$user_r_libs" UPDATE_LOCK_FLAG="$update_lock" Rscript -e '
    proj <- Sys.getenv("TARGETCNV_PROJ")
    ext_libs <- strsplit(Sys.getenv("USER_R_LIBS"), ":", fixed = TRUE)[[1]]
    ext_libs <- ext_libs[nzchar(ext_libs) & dir.exists(ext_libs)]
    renv::hydrate(project = proj, sources = ext_libs, prompt = FALSE)
    if (identical(Sys.getenv("UPDATE_LOCK_FLAG"), "1")) {
      renv::snapshot(project = proj, prompt = FALSE)
      cat("[dev-hydrate] renv.lock actualizado (--update-lock).\n")
    } else {
      cat("[dev-hydrate] Librería local hidratada sin modificar renv.lock.\n")
    }
  '
}

cmd_doctor() {
  echo "[doctor] Verificando utilidades del sistema y paquetes de R..."
  local status=0

  for tool in Rscript curl unzip md5sum; do
    if command -v "$tool" >/dev/null 2>&1; then
      printf "  [OK]    %-10s %s\n" "$tool" "$(command -v "$tool")"
    else
      printf "  [FALTA] %-10s\n" "$tool"
      status=1
    fi
  done

  if command -v quarto >/dev/null 2>&1; then
    printf "  [OK]    %-10s %s\n" "quarto" "$(quarto --version 2>/dev/null || echo 'instalado')"
  else
    echo "  [AVISO] quarto no detectado (requerido solo para compilar reportes HTML)"
  fi

  Rscript -e '
    pkgs <- c("optparse", "dplyr", "tidyr", "stringr", "ggplot2", "caret", "glmnet",
              "klaR", "pROC", "GenomicRanges", "Seqinfo", "rtracklayer", "readxl", "kernelshap",
              "shapviz", "DT", "plotly", "cowplot", "NeuralNetTools", "caretEnsemble",
              "kernlab", "xgboost", "gbm", "randomForest", "RSNNS", "deepnet", "ipred",
              "doParallel", "knitr", "rmarkdown", "testthat", "withr")
    missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing_pkgs) == 0) {
      cat("  [OK]    Todos los paquetes requeridos de R están instalados.\n")
    } else {
      cat("  [FALTA] Paquetes ausentes:", paste(missing_pkgs, collapse = ", "), "\n")
      cat("          Ejecute: ./targetcnv.sh restore\n")
      quit(status = 1)
    }
  ' || status=1

  if [ "$status" -eq 0 ]; then
    echo "[doctor] Entorno verificado."
  else
    echo "[doctor] Corrija los elementos pendientes antes de correr el flujo." >&2
  fi
  return "$status"
}

cmd_zenodo() {
  local out_dir="."
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --out-dir) out_dir="$2"; shift 2 ;;
      -h|--help)
        echo "Uso: ./targetcnv.sh zenodo [--out-dir DIR]"
        exit 0 ;;
      *) echo "[ERROR] Opción desconocida en zenodo: $1" >&2; exit 1 ;;
    esac
  done

  mkdir -p "$out_dir"
  local zip_path="${out_dir%/}/${ZENODO_ZIP}"
  rm -f "$zip_path"

  echo "[zenodo] Descargando modelos preentrenados (DOI: 10.5281/zenodo.23463987)..."
  if command -v curl >/dev/null 2>&1; then
    if ! curl -L --fail --retry 3 --connect-timeout 30 --progress-bar -o "$zip_path" "$ZENODO_URL"; then
      echo "[ERROR] Falló la descarga con curl desde Zenodo." >&2
      rm -f "$zip_path"
      exit 1
    fi
  elif command -v wget >/dev/null 2>&1; then
    if ! wget --tries=3 --timeout=30 -O "$zip_path" "$ZENODO_URL"; then
      echo "[ERROR] Falló la descarga con wget desde Zenodo." >&2
      rm -f "$zip_path"
      exit 1
    fi
  else
    echo "[ERROR] Se requiere curl o wget para descargar desde Zenodo." >&2
    exit 1
  fi

  if [ ! -s "$zip_path" ]; then
    echo "[ERROR] El archivo descargado está vacío." >&2
    rm -f "$zip_path"
    exit 1
  fi

  local actual_md5=""
  if command -v md5sum >/dev/null 2>&1; then
    actual_md5="$(md5sum "$zip_path" | awk '{print $1}')"
  elif command -v md5 >/dev/null 2>&1; then
    actual_md5="$(md5 -q "$zip_path")"
  fi

  if [ -n "$actual_md5" ]; then
    if [ "$actual_md5" = "$ZENODO_MD5" ]; then
      echo "  [OK] MD5 verificado (${actual_md5})"
    else
      echo "[ERROR] MD5 inválido: esperado ${ZENODO_MD5}, recibido ${actual_md5}" >&2
      rm -f "$zip_path"
      exit 1
    fi
  else
    echo "  [AVISO] No se encontró md5sum ni md5; se omite la validación MD5."
  fi

  if ! unzip -o "$zip_path" -d "$out_dir"; then
    echo "[ERROR] Falló la descompresión de $zip_path." >&2
    rm -f "$zip_path"
    exit 1
  fi
  rm -f "$zip_path"

  local target_models_dir="${out_dir%/}/models"
  if [ ! -f "$target_models_dir/youden_thresholds.rds" ] || [ ! -f "$target_models_dir/heuristic_model.rds" ]; then
    echo "[ERROR] Faltan youden_thresholds.rds o heuristic_model.rds en ${target_models_dir}/" >&2
    exit 1
  fi

  local n_ind_models
  n_ind_models="$(find "$target_models_dir" -maxdepth 1 -name 'model_*.rds' | wc -l | tr -d ' ')"
  if [ -f "$target_models_dir/consensus_spec.rds" ] && [ "$n_ind_models" -ge 18 ]; then
    echo "  [OK] Inventario completo: consensus_spec.rds y ${n_ind_models} modelos en ${target_models_dir}/"
  else
    echo "  [AVISO] Se extrajeron ${n_ind_models} archivos model_*.rds (falta consensus_spec.rds o el paquete de 18 modelos)."
    echo "          Para regenerar los 18 clasificadores localmente ejecute: ./targetcnv.sh train --clean_models"
  fi
}

cmd_infer() {
  local auto_zenodo=0
  local models_dir="models"
  local pass_args=()

  while [[ $# -gt 0 ]]; do
    case "$1" in
      --auto-zenodo)
        auto_zenodo=1
        shift ;;
      -m|--models_dir)
        models_dir="${2:-models}"
        pass_args+=("$1" "$models_dir")
        shift 2 ;;
      --sources)
        if [[ $# -ge 2 && "${2:-}" != -* ]]; then
          pass_args+=("--sources" "$2")
          shift 2
        else
          pass_args+=("--sources" "$(resolve_default_sources_file)")
          shift 1
        fi ;;
      *)
        pass_args+=("$1")
        shift ;;
    esac
  done

  if [ "$auto_zenodo" -eq 1 ]; then
    if [ ! -f "$models_dir/youden_thresholds.rds" ] || [ ! -f "$models_dir/heuristic_model.rds" ]; then
      echo "[infer] No se detectaron modelos en '$models_dir'; descargando desde Zenodo..."
      cmd_zenodo --out-dir "."
    fi
  fi

  Rscript "$SCRIPT_DIR/scripts/03_cnv_inference.R" "${pass_args[@]}"
}

cmd_reproduce_storvum() {
  local cores=4
  while [[ $# -gt 0 ]]; do
    case "$1" in
      -c|--cores) cores="$2"; shift 2 ;;
      -h|--help)
        echo "Uso: ./targetcnv.sh reproduce-storvum [--cores N]"
        exit 0 ;;
      *) echo "[ERROR] Opción desconocida: $1" >&2; exit 1 ;;
    esac
  done

  echo "[reproduce-storvum 1/3] Preparación espacial (S. torvum)..."
  Rscript "$SCRIPT_DIR/scripts/01_prepare_data.R" \
    --gff_ref data/reference/Storvum_genes.gff3 \
    --captus_gff data/reference/Storvum_hit_contigs.gff \
    --captus_stats data/reference/captus-extract_stats.tsv \
    --dict_ref data/reference/Genes_Capture_May_2026.xlsx \
    --sample Storvum \
    --out_dir data/

  echo "[reproduce-storvum 2/3] Entrenamiento determinista y SHAP..."
  Rscript "$SCRIPT_DIR/scripts/02_train_models.R" \
    --data data/data_prepared.rds \
    --out_models models/ \
    --cv_folds 5 \
    --sampling down \
    --cores "$cores" \
    --seed 1234 \
    --report_file reports/reporte_entrenamiento.html

  echo "[reproduce-storvum 3/3] Inferencia multimuestra y matriz CNV..."
  Rscript "$SCRIPT_DIR/scripts/03_cnv_inference.R" \
    --input data/CAPTUS_Muestras/captus-extract_stats.tsv \
    --models_dir models/ \
    --out_dir results/ \
    --cov_threshold 50.0 \
    --min_votes 4 \
    --top_loci 30 \
    --report_file reports/reporte_inferencia.html

  echo "[reproduce-storvum] Experimento de Solanum torvum completado."
}

cmd_run_all() {
  local gff_ref=""
  local captus_gff=""
  local captus_stats=""
  local dict_ref=""
  local sample_name=""
  local target_mode="overlap"
  local infer_input=""
  local out_data="data"
  local out_models="models"
  local out_results="results"
  local cv_folds=5
  local sampling="down"
  local cores=4
  local seed=1234
  local min_votes=4
  local cov_threshold="50.0"
  local top_loci=30
  local report_file="reports/reporte_entrenamiento.html"
  local infer_report_file=""

  while [[ $# -gt 0 ]]; do
    case "$1" in
      -g|--gff_ref) gff_ref="$2"; shift 2 ;;
      -c|--captus_gff) captus_gff="$2"; shift 2 ;;
      -s|--captus_stats) captus_stats="$2"; shift 2 ;;
      -d|--dict_ref) dict_ref="$2"; shift 2 ;;
      -p|--sample) sample_name="$2"; shift 2 ;;
      --target_mode) target_mode="$2"; shift 2 ;;
      --genome_scale)
        if [[ "$2" == "contig" || "$2" == "contigs" ]]; then target_mode="target_attr"; else target_mode="overlap"; fi
        shift 2 ;;
      -i|--infer_input) infer_input="$2"; shift 2 ;;
      --out_data) out_data="$2"; shift 2 ;;
      -m|--out_models) out_models="$2"; shift 2 ;;
      -o|--out_results) out_results="$2"; shift 2 ;;
      -k|--cv_folds) cv_folds="$2"; shift 2 ;;
      --sampling) sampling="$2"; shift 2 ;;
      --cores) cores="$2"; shift 2 ;;
      --seed) seed="$2"; shift 2 ;;
      -v|--min_votes) min_votes="$2"; shift 2 ;;
      --cov_threshold) cov_threshold="$2"; shift 2 ;;
      -t|--top_loci) top_loci="$2"; shift 2 ;;
      -r|--report_file) report_file="$2"; shift 2 ;;
      --infer_report_file) infer_report_file="$2"; shift 2 ;;
      -h|--help)
        echo "Uso: ./targetcnv.sh run-all -g <ref.gff3> -c <captus.gff> -s <stats.tsv> -d <dict> [-p <sample>] [-i <infer.tsv>]"
        exit 0 ;;
      *) echo "[ERROR] Argumento desconocido en run-all: $1" >&2; exit 1 ;;
    esac
  done

  if [ -z "$captus_gff" ] || [ -z "$captus_stats" ] || [ -z "$dict_ref" ]; then
    echo "[ERROR] 'run-all' requiere -c (--captus_gff), -s (--captus_stats) y -d (--dict_ref)." >&2
    exit 1
  fi
  if [ -z "$infer_input" ]; then
    infer_input="$captus_stats"
  fi
  if [ -z "$infer_report_file" ]; then
    if [[ "$report_file" != "reports/reporte_entrenamiento.html" ]]; then
      if [[ "$report_file" == *entrenamiento* ]]; then
        infer_report_file="${report_file/entrenamiento/inferencia}"
      elif [[ "$report_file" == *training* ]]; then
        infer_report_file="${report_file/training/inference}"
      else
        infer_report_file="${report_file%.html}_inferencia.html"
      fi
    elif [ -n "$sample_name" ]; then
      infer_report_file="reports/reporte_inferencia_${sample_name}.html"
    elif [[ "$out_results" != "results" ]]; then
      infer_report_file="reports/reporte_inferencia_$(basename "$out_results").html"
    else
      infer_report_file="reports/reporte_inferencia.html"
    fi
  fi

  local prep_args=(--captus_gff "$captus_gff" --captus_stats "$captus_stats" --target_mode "$target_mode" --out_dir "$out_data")
  [ -n "$gff_ref" ] && prep_args+=(--gff_ref "$gff_ref")
  [ -n "$dict_ref" ] && prep_args+=(--dict_ref "$dict_ref")
  [ -n "$sample_name" ] && prep_args+=(--sample "$sample_name")

  echo "[run-all 1/3] 01_prepare_data.R..."
  Rscript "$SCRIPT_DIR/scripts/01_prepare_data.R" "${prep_args[@]}"

  echo "[run-all 2/3] 02_train_models.R..."
  Rscript "$SCRIPT_DIR/scripts/02_train_models.R" \
    --data "$out_data/data_prepared.rds" \
    --out_models "$out_models" \
    --cv_folds "$cv_folds" \
    --sampling "$sampling" \
    --cores "$cores" \
    --seed "$seed" \
    --report_file "$report_file"

  echo "[run-all 3/3] 03_cnv_inference.R..."
  Rscript "$SCRIPT_DIR/scripts/03_cnv_inference.R" \
    --input "$infer_input" \
    --models_dir "$out_models" \
    --out_dir "$out_results" \
    --cov_threshold "$cov_threshold" \
    --min_votes "$min_votes" \
    --top_loci "$top_loci" \
    --report_file "$infer_report_file"

  echo "[run-all] Salidas: $out_data/data_prepared.rds | $out_models/ | $out_results/cnv_matrix.tsv"
  echo "[run-all] Reportes: $report_file | $infer_report_file"
}

cmd_report() {
  local models_dir="models"
  local data_file="data/data_prepared.rds"
  local sources_file=""
  local output_file="reports/reporte_entrenamiento.html"

  while [[ $# -gt 0 ]]; do
    case "$1" in
      -m|--models-dir) models_dir="$2"; shift 2 ;;
      -d|--data) data_file="$2"; shift 2 ;;
      --sources)
        if [[ $# -ge 2 && "${2:-}" != -* ]]; then
          sources_file="$2"
          shift 2
        else
          sources_file="$(resolve_default_sources_file)"
          shift 1
        fi ;;
      -o|--output) output_file="$2"; shift 2 ;;
      -h|--help)
        echo "Uso: ./targetcnv.sh report [--models-dir models] [--data data/data_prepared.rds] [--sources FILE] [--output reports/reporte_entrenamiento.html]"
        exit 0 ;;
      *) echo "[ERROR] Opción desconocida en report: $1" >&2; exit 1 ;;
    esac
  done

  local out_base abs_models abs_data abs_sources abs_out
  out_base="$(basename "$output_file")"
  mkdir -p "$(dirname "$output_file")"
  abs_models="$(resolve_abs_path "$models_dir")"
  abs_data="$(resolve_abs_path "$data_file")"
  abs_out="$(resolve_abs_path "$output_file")"
  abs_sources=""
  if [ -n "$sources_file" ]; then
    abs_sources="$(resolve_abs_path "$sources_file")"
  fi

  # Si SCRIPT_DIR es de solo lectura (ej. contenedor Apptainer), compilar en un directorio temporal
  if [ -w "$SCRIPT_DIR" ]; then
    quarto render "$SCRIPT_DIR/training_report.qmd" \
      -P project_root:"$SCRIPT_DIR" \
      -P models_dir:"$abs_models" \
      -P data_path:"$abs_data" \
      -P data_file:"$abs_data" \
      -P sources_file:"$abs_sources" \
      --output "$out_base"

    if [ -f "$SCRIPT_DIR/$out_base" ] && [ "$abs_out" != "$SCRIPT_DIR/$out_base" ]; then
      mv -f "$SCRIPT_DIR/$out_base" "$abs_out"
    fi
  else
    local tmp_render
    tmp_render="$(mktemp -d)"
    trap 'rm -rf "$tmp_render"' RETURN
    cp "$SCRIPT_DIR/training_report.qmd" "$tmp_render/"
    cp -r "$SCRIPT_DIR/assets" "$tmp_render/"
    quarto render "$tmp_render/training_report.qmd" \
      -P project_root:"$SCRIPT_DIR" \
      -P models_dir:"$abs_models" \
      -P data_path:"$abs_data" \
      -P data_file:"$abs_data" \
      -P sources_file:"$abs_sources" \
      --output "$out_base"
    mv -f "$tmp_render/$out_base" "$abs_out"
  fi
  echo "[report] Reporte compilado: ${output_file}"
}

cmd_transfer() {
  local models_dir="models"
  local donor_data="data/data_prepared.rds"
  local cit_data="data/citrullifolium/data_prepared.rds"
  local ste_data="data/stellatiglandulosum/data_prepared.rds"
  local cit_gff_ref=""
  local cit_captus_gff=""
  local ste_captus_gff=""
  local captus_stats=""
  local dict_ref=""
  local sources_file=""
  local extra_species=""
  local cov_threshold="50.0"
  local min_votes="4"
  local out_dir="results"
  local output_file="reports/reporte_transferencia.html"
  local no_report=0

  while [[ $# -gt 0 ]]; do
    case "$1" in
      -m|--models-dir|--models_dir) models_dir="$2"; shift 2 ;;
      --donor-data|--donor_data) donor_data="$2"; shift 2 ;;
      --cit-data|--cit_data) cit_data="$2"; shift 2 ;;
      --ste-data|--ste_data) ste_data="$2"; shift 2 ;;
      --cit-gff-ref|--cit_gff_ref) cit_gff_ref="$2"; shift 2 ;;
      --cit-captus-gff|--cit_captus_gff) cit_captus_gff="$2"; shift 2 ;;
      --ste-captus-gff|--ste_captus_gff) ste_captus_gff="$2"; shift 2 ;;
      -s|--stats|--captus_stats) captus_stats="$2"; shift 2 ;;
      -d|--dict|--dict_ref) dict_ref="$2"; shift 2 ;;
      -c|--cov-threshold|--cov_threshold) cov_threshold="$2"; shift 2 ;;
      -v|--min-votes|--min_votes) min_votes="$2"; shift 2 ;;
      --no-report|--no_render_report) no_report=1; shift ;;
      --sources)
        if [[ $# -ge 2 && "${2:-}" != -* ]]; then
          sources_file="$2"
          shift 2
        else
          sources_file="$(resolve_default_sources_file)"
          shift 1
        fi ;;
      --extra-species|--extra_species) extra_species="$2"; shift 2 ;;
      --out-dir|--out_dir) out_dir="$2"; shift 2 ;;
      -o|--output|--report_file) output_file="$2"; shift 2 ;;
      -h|--help)
        echo "Uso: ./targetcnv.sh transfer [--models-dir models] [--cit-data FILE] [--ste-data FILE] [--cov-threshold 50.0] [--out-dir results] [--output reports/reporte_transferencia.html]"
        exit 0 ;;
      *) echo "[ERROR] Opción desconocida en transfer: $1" >&2; exit 1 ;;
    esac
  done

  local tr_args=(
    --models_dir "$models_dir"
    --donor_data "$donor_data"
    --cit_data "$cit_data"
    --ste_data "$ste_data"
    --cov_threshold "$cov_threshold"
    --min_votes "$min_votes"
    --out_dir "$out_dir"
    --report_file "$output_file"
  )
  [ -n "$cit_gff_ref" ]    && tr_args+=(--cit_gff_ref "$cit_gff_ref")
  [ -n "$cit_captus_gff" ] && tr_args+=(--cit_captus_gff "$cit_captus_gff")
  [ -n "$ste_captus_gff" ] && tr_args+=(--ste_captus_gff "$ste_captus_gff")
  [ -n "$captus_stats" ]   && tr_args+=(--captus_stats "$captus_stats")
  [ -n "$dict_ref" ]       && tr_args+=(--dict_ref "$dict_ref")
  [ -n "$sources_file" ]   && tr_args+=(--sources "$sources_file")
  [ -n "$extra_species" ]  && tr_args+=(--extra_species "$extra_species")
  [ "$no_report" -eq 1 ]   && tr_args+=(--no_render_report)

  Rscript "$SCRIPT_DIR/scripts/04_evaluate_transferability.R" "${tr_args[@]}"
}

SUBCOMMAND="${1:-help}"
shift || true

case "$SUBCOMMAND" in
  restore|setup)
    cmd_restore "$@" ;;
  dev-hydrate)
    cmd_dev_hydrate "$@" ;;
  snapshot)
    echo "[AVISO] El antiguo subcomando 'snapshot' fue desactivado para proteger renv.lock." >&2
    echo "        Para instalar dependencias use: ./targetcnv.sh restore" >&2
    echo "        Para desarrollo interno use:    ./targetcnv.sh dev-hydrate [--update-lock]" >&2
    exit 1 ;;
  test|mwe)
    bash "$SCRIPT_DIR/run_test_pipeline.sh" "$@" ;;
  zenodo|download-models)
    cmd_zenodo "$@" ;;
  infer|inference)
    cmd_infer "$@" ;;
  transfer|transferability|report-transfer)
    cmd_transfer "$@" ;;
  prepare|prep)
    Rscript "$SCRIPT_DIR/scripts/01_prepare_data.R" "$@" ;;
  train)
    Rscript "$SCRIPT_DIR/scripts/02_train_models.R" "$@" ;;
  run-all|all)
    cmd_run_all "$@" ;;
  reproduce-storvum|storvum)
    cmd_reproduce_storvum "$@" ;;
  report)
    cmd_report "$@" ;;
  doctor|check)
    cmd_doctor "$@" ;;
  -h|--help|help)
    show_help ;;
  *)
    echo "[ERROR] Subcomando desconocido: '$SUBCOMMAND'" >&2
    show_help
    exit 1 ;;
esac
