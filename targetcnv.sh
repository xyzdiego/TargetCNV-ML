#!/usr/bin/env bash

# TargetCNV-ML: Pipeline wrapper script (targetcnv.sh)
# Interfaz de línea de comandos para preparación espacial, entrenamiento,
# inferencia de CNV, transferencia filogenética y gestión reproducible del entorno.

set -euo pipefail

GREEN='\033[0;32m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BOLD='\033[1m'
NC='\033[0m'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

ZENODO_URL="https://zenodo.org/api/records/23463987/files/TargetCNV_Models.zip/content"
ZENODO_ZIP="TargetCNV_Models.zip"
ZENODO_MD5="13ec7b35711158079f962a96314eb671"

resolve_default_sources_file() {
  if [ -f "$SCRIPT_DIR/data/reference/sources_citations.tsv" ]; then
    echo "data/reference/sources_citations.tsv"
  elif [ -f "$SCRIPT_DIR/sources_citations.tsv" ]; then
    echo "sources_citations.tsv"
  else
    echo "data/reference/sources_citations.tsv"
  fi
}

print_banner() {
  echo -e "${BOLD}${CYAN}======================================================================${NC}"
  echo -e "${BOLD}${CYAN}  TargetCNV-ML | Pipeline CLI Wrapper for Spatial ML and CNV Calling  ${NC}"
  echo -e "${BOLD}${CYAN}======================================================================${NC}"
}

show_help() {
  print_banner
  echo -e "
${BOLD}USO / USAGE:${NC}
  ./targetcnv.sh <subcomando> [opciones]

${BOLD}SUBCOMANDOS DISPONIBLES:${NC}
  ${GREEN}restore${NC} (alias: ${GREEN}setup${NC})
      Restaura de forma segura la librería de R exactamente desde renv.lock
      mediante renv::restore(prompt = FALSE), sin alterar el lockfile.

  ${GREEN}doctor${NC}
      Verifica la disponibilidad de R, paquetes requeridos, Quarto y utilidades
      del sistema antes de ejecutar el flujo de trabajo.

  ${GREEN}test${NC} [--sync-root] [--no-report]
      Ejecuta el Ejemplo Mínimo Funcional (MWE) con datos simulados en data/test_data/.
      Nota: --sync-root sobrescribe models/ y results/ en la raíz del repositorio.

  ${GREEN}zenodo${NC} [--out-dir DIR]
      Descarga los modelos preentrenados desde Zenodo (DOI: 10.5281/zenodo.23463987),
      verifica la suma MD5, limpia archivos corruptos si falla y valida el contenido.

  ${GREEN}prepare${NC} [--genome_scale chromosome|contig] [opciones de 01_prepare_data.R]
      Ejecuta el Módulo 01 (intersección espacial GenomicRanges o atributo Target).

  ${GREEN}train${NC} [opciones de 02_train_models.R]
      Ejecuta el Módulo 02 (colinealidad de Spearman, Youden, clasificadores ML con
      semillas deterministas, KernelSHAP y reporte HTML de solo lectura).

  ${GREEN}infer${NC} -i <captus_stats.tsv> [--auto-zenodo] [opciones de 03_cnv_inference.R]
      Ejecuta el Módulo 03 (consenso de 7 modelos, matriz CNV, calibración 30%-80%
      de --cov_threshold y concordancia por ortogrupo si existe referencia).

  ${GREEN}transfer${NC} [--models-dir DIR] [--cit-data FILE] [--ste-data FILE] [--output FILE]
      Ejecuta el Módulo 04 (04_evaluate_transferability.R) sobre especies receptoras
      y compila el reporte de transferencia en modo de solo lectura.

  ${GREEN}reproduce-storvum${NC} [--cores N]
      Reproduce secuencialmente el experimento completo con Solanum torvum (01 -> 02 -> 03).

  ${GREEN}run-all${NC} -g <ref.gff3> -c <captus.gff> -s <stats.tsv> -d <dict> [-i <infer.tsv>]
      Ejecuta de principio a fin los módulos 01, 02 y 03 para un nuevo conjunto de datos.

  ${GREEN}report${NC} [--models-dir DIR] [--data FILE] [--output FILE]
      Compila training_report.qmd en modo de solo lectura sin reentrenar modelos.

${BOLD}EJEMPLOS:${NC}
  ./targetcnv.sh restore
  ./targetcnv.sh doctor
  ./targetcnv.sh test
  ./targetcnv.sh infer -i data/CAPTUS_Muestras/captus-extract_stats.tsv --auto-zenodo
  ./targetcnv.sh transfer
  ./targetcnv.sh reproduce-storvum --cores 4
"
}

cmd_restore() {
  print_banner
  echo -e "\n${BOLD}${BLUE}[RESTORE] Restaurando entorno reproducible desde renv.lock...${NC}"
  if [ ! -f "$SCRIPT_DIR/renv.lock" ]; then
    echo -e "${RED}[ERROR] No se encontró renv.lock en $SCRIPT_DIR.${NC}"
    exit 1
  fi
  Rscript -e '
    if (!requireNamespace("renv", quietly = TRUE)) {
      stop("El paquete renv no está disponible. Instale renv antes de ejecutar restore.", call. = FALSE)
    }
    renv::restore(prompt = FALSE)
  '
  echo -e "${BOLD}${GREEN}[OK] Librería restaurada conforme a renv.lock.${NC}\n"
}

cmd_dev_hydrate() {
  local update_lock=0
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --update-lock) update_lock=1; shift ;;
      -h|--help)
        echo "Uso interno de desarrollo: ./targetcnv.sh dev-hydrate [--update-lock]"
        echo "AVISO: Comando reservado para mantenimiento del entorno de desarrollo."
        echo "       Para instalar dependencias como usuario final, use: ./targetcnv.sh restore"
        exit 0 ;;
      *) echo -e "${RED}[ERROR] Opción desconocida en dev-hydrate: $1${NC}"; exit 1 ;;
    esac
  done

  print_banner
  echo -e "\n${YELLOW}[AVISO DE DESARROLLO] Ejecutando hidratación local desde librerías del sistema.${NC}"
  echo -e "${YELLOW}No utilice este comando para reproducir el entorno; use './targetcnv.sh restore'.${NC}"
  USER_R_LIBS="$(Rscript --vanilla -e 'cat(paste(.libPaths(), collapse = ":"))')"
  UPDATE_LOCK_FLAG="$update_lock" Rscript -e "
    ext_libs <- strsplit('$USER_R_LIBS', ':')[[1]]
    ext_libs <- ext_libs[nzchar(ext_libs) & dir.exists(ext_libs)]
    renv::hydrate(sources = ext_libs, prompt = FALSE)
    if (identical(Sys.getenv('UPDATE_LOCK_FLAG'), '1')) {
      renv::snapshot(prompt = FALSE)
      cat('[OK] renv.lock actualizado explícitamente (--update-lock).\n')
    } else {
      cat('[OK] Librería local hidratada sin modificar renv.lock.\n')
    }
  "
}

cmd_doctor() {
  print_banner
  echo -e "\n${BOLD}${BLUE}[DOCTOR] Verificando utilidades del sistema y dependencias de R...${NC}\n"
  local status=0

  for tool in Rscript curl unzip md5sum; do
    if command -v "$tool" >/dev/null 2>&1; then
      echo -e "  ${GREEN}[OK]${NC}    Herramienta detectada : ${BOLD}$tool${NC} ($(command -v "$tool"))"
    else
      echo -e "  ${RED}[FALTA]${NC} Herramienta ausente   : ${BOLD}$tool${NC}"
      status=1
    fi
  done

  if command -v quarto >/dev/null 2>&1; then
    echo -e "  ${GREEN}[OK]${NC}    Quarto CLI            : ${BOLD}$(quarto --version 2>/dev/null || echo 'Instalado')${NC}"
  else
    echo -e "  ${YELLOW}[AVISO]${NC} Quarto CLI no detectado (necesario únicamente para compilar reportes HTML)"
  fi

  echo -e "\n${BOLD}${BLUE}[DOCTOR] Verificando paquetes requeridos en R...${NC}"
  Rscript -e '
    pkgs <- c("optparse", "dplyr", "tidyr", "stringr", "ggplot2", "caret", "glmnet",
              "klaR", "pROC", "GenomicRanges", "rtracklayer", "readxl", "kernelshap",
              "shapviz", "DT", "plotly", "cowplot", "NeuralNetTools", "caretEnsemble",
              "kernlab", "xgboost", "gbm", "randomForest", "RSNNS", "deepnet", "ipred",
              "doParallel", "knitr", "rmarkdown", "testthat", "withr")
    missing_pkgs <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) == 0) {
      cat("  [OK] Todos los paquetes requeridos de R están instalados.\n")
    } else {
      cat("  [FALTA] Paquetes ausentes:", paste(missing_pkgs, collapse = ", "), "\n")
      cat("          Ejecute: ./targetcnv.sh restore\n")
      quit(status = 1)
    }
  ' || status=1

  echo ""
  if [ "$status" -eq 0 ]; then
    echo -e "${BOLD}${GREEN}[OK] Entorno verificado correctamente.${NC}\n"
  else
    echo -e "${BOLD}${YELLOW}[AVISO] Corrija los elementos pendientes antes de ejecutar el pipeline completo.${NC}\n"
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
      *) echo -e "${RED}[ERROR] Opción desconocida en zenodo: $1${NC}"; exit 1 ;;
    esac
  done

  print_banner
  echo -e "\n${BOLD}${BLUE}[ZENODO] Descargando modelos preentrenados (DOI: 10.5281/zenodo.23463987)...${NC}"
  rm -f "$ZENODO_ZIP"

  if command -v curl >/dev/null 2>&1; then
    if ! curl -L --fail --retry 3 --connect-timeout 30 --progress-bar -o "$ZENODO_ZIP" "$ZENODO_URL"; then
      echo -e "${RED}[ERROR] Falló la descarga con curl desde Zenodo. Limpiando archivo parcial.${NC}"
      rm -f "$ZENODO_ZIP"
      exit 1
    fi
  elif command -v wget >/dev/null 2>&1; then
    if ! wget --tries=3 --timeout=30 -O "$ZENODO_ZIP" "$ZENODO_URL"; then
      echo -e "${RED}[ERROR] Falló la descarga con wget desde Zenodo. Limpiando archivo parcial.${NC}"
      rm -f "$ZENODO_ZIP"
      exit 1
    fi
  else
    echo -e "${RED}[ERROR] Se requiere curl o wget para descargar artefactos desde Zenodo.${NC}"
    exit 1
  fi

  if [ ! -s "$ZENODO_ZIP" ]; then
    echo -e "${RED}[ERROR] El archivo descargado está vacío o corrupto.${NC}"
    rm -f "$ZENODO_ZIP"
    exit 1
  fi

  if command -v md5sum >/dev/null 2>&1; then
    echo -e "${BOLD}${BLUE}[ZENODO] Verificando suma de comprobación MD5...${NC}"
    local actual_md5
    actual_md5="$(md5sum "$ZENODO_ZIP" | awk '{print $1}')"
    if [ "$actual_md5" = "$ZENODO_MD5" ]; then
      echo -e "  ${GREEN}[OK] Suma MD5 verificada (${actual_md5})${NC}"
    else
      echo -e "  ${RED}[ERROR] Suma MD5 inválida: esperado ${ZENODO_MD5}, recibido ${actual_md5}${NC}"
      rm -f "$ZENODO_ZIP"
      exit 1
    fi
  elif command -v md5 >/dev/null 2>&1; then
    echo -e "${BOLD}${BLUE}[ZENODO] Verificando suma de comprobación MD5...${NC}"
    local actual_md5
    actual_md5="$(md5 -q "$ZENODO_ZIP")"
    if [ "$actual_md5" = "$ZENODO_MD5" ]; then
      echo -e "  ${GREEN}[OK] Suma MD5 verificada (${actual_md5})${NC}"
    else
      echo -e "  ${RED}[ERROR] Suma MD5 inválida: esperado ${ZENODO_MD5}, recibido ${actual_md5}${NC}"
      rm -f "$ZENODO_ZIP"
      exit 1
    fi
  else
    echo -e "  ${YELLOW}[AVISO] No se encontró md5sum ni md5: se omitió la validación criptográfica.${NC}"
  fi

  mkdir -p "$out_dir"
  if ! unzip -o "$ZENODO_ZIP" -d "$out_dir"; then
    echo -e "${RED}[ERROR] Falló la descompresión de $ZENODO_ZIP. Eliminando archivo corrupto.${NC}"
    rm -f "$ZENODO_ZIP"
    exit 1
  fi
  rm -f "$ZENODO_ZIP"

  local target_models_dir="${out_dir}/models"
  if [ ! -f "$target_models_dir/youden_thresholds.rds" ] || [ ! -f "$target_models_dir/heuristic_model.rds" ]; then
    echo -e "${RED}[ERROR] El paquete extraído no contiene youden_thresholds.rds o heuristic_model.rds en ${target_models_dir}/${NC}"
    exit 1
  fi

  local n_ind_models
  n_ind_models="$(find "$target_models_dir" -maxdepth 1 -name 'model_*.rds' | wc -l | tr -d ' ')"
  if [ -f "$target_models_dir/consensus_spec.rds" ] && [ "$n_ind_models" -ge 18 ]; then
    echo -e "  ${GREEN}[OK] Inventario verificado: consensus_spec.rds presente y ${n_ind_models} modelos serializados.${NC}"
  else
    echo -e "  ${YELLOW}[AVISO] El archivo de Zenodo contiene ${n_ind_models} archivos model_*.rds (falta consensus_spec.rds o el conjunto completo de 18+1 modelos actualizados).${NC}"
    echo -e "  ${YELLOW}        Para regenerar localmente los 18 clasificadores (incluyendo XGBoost), Stacking y consensus_spec.rds, ejecute: ./targetcnv.sh train --clean_models${NC}"
  fi

  echo -e "${GREEN}[OK] Modelos disponibles en: ${BOLD}${target_models_dir}/${NC}\n"
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
      echo -e "${YELLOW}[AUTO-ZENODO] No se detectaron modelos en '$models_dir'. Descargando desde Zenodo...${NC}"
      cmd_zenodo --out-dir "."
    fi
  fi

  Rscript ./scripts/03_cnv_inference.R "${pass_args[@]}"
}

cmd_reproduce_storvum() {
  local cores=4
  while [[ $# -gt 0 ]]; do
    case "$1" in
      -c|--cores) cores="$2"; shift 2 ;;
      -h|--help)
        echo "Uso: ./targetcnv.sh reproduce-storvum [--cores N]"
        exit 0 ;;
      *) echo -e "${RED}[ERROR] Opción desconocida: $1${NC}"; exit 1 ;;
    esac
  done

  print_banner
  echo -e "\n${BOLD}${BLUE}[STORVUM 1/3] Preparación de datos e intersección espacial...${NC}"
  Rscript ./scripts/01_prepare_data.R \
    --gff_ref data/reference/Storvum_genes.gff3 \
    --captus_gff data/reference/Storvum_hit_contigs.gff \
    --captus_stats data/reference/captus-extract_stats.tsv \
    --dict_ref data/reference/Genes_Capture_May_2026.xlsx \
    --sample Storvum \
    --out_dir data/

  echo -e "\n${BOLD}${BLUE}[STORVUM 2/3] Entrenamiento supervisado determinista, colinealidad, SHAP y reporte HTML...${NC}"
  Rscript ./scripts/02_train_models.R \
    --data data/data_prepared.rds \
    --out_models models/ \
    --cv_folds 5 \
    --sampling down \
    --cores "$cores" \
    --seed 1234 \
    --report_file reports/reporte_entrenamiento.html

  echo -e "\n${BOLD}${BLUE}[STORVUM 3/3] Inferencia multimuestra y extracción de matriz CNV...${NC}"
  Rscript ./scripts/03_cnv_inference.R \
    --input data/CAPTUS_Muestras/captus-extract_stats.tsv \
    --models_dir models/ \
    --out_dir results/ \
    --cov_threshold 50.0 \
    --min_votes 4 \
    --top_loci 30 \
    --report_file reports/reporte_inferencia.html

  echo -e "\n${BOLD}${GREEN}[OK] Experimento de referencia con Solanum torvum reproducido correctamente.${NC}\n"
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
      *) echo -e "${RED}[ERROR] Argumento desconocido en run-all: $1${NC}"; exit 1 ;;
    esac
  done

  if [ -z "$captus_gff" ] || [ -z "$captus_stats" ] || [ -z "$dict_ref" ]; then
    echo -e "${RED}[ERROR] 'run-all' requiere al menos -c (--captus_gff), -s (--captus_stats) y -d (--dict_ref).${NC}"
    echo -e "${YELLOW}        El diccionario (-d) es indispensable para clasificar TP y FP antes del entrenamiento.${NC}"
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

  print_banner
  local prep_args=(--captus_gff "$captus_gff" --captus_stats "$captus_stats" --target_mode "$target_mode" --out_dir "$out_data")
  [ -n "$gff_ref" ] && prep_args+=(--gff_ref "$gff_ref")
  [ -n "$dict_ref" ] && prep_args+=(--dict_ref "$dict_ref")
  [ -n "$sample_name" ] && prep_args+=(--sample "$sample_name")

  echo -e "\n${BOLD}${BLUE}[PASO 1/3] Ejecutando 01_prepare_data.R...${NC}"
  Rscript ./scripts/01_prepare_data.R "${prep_args[@]}"

  echo -e "\n${BOLD}${BLUE}[PASO 2/3] Ejecutando 02_train_models.R...${NC}"
  Rscript ./scripts/02_train_models.R \
    --data "$out_data/data_prepared.rds" \
    --out_models "$out_models" \
    --cv_folds "$cv_folds" \
    --sampling "$sampling" \
    --cores "$cores" \
    --seed "$seed" \
    --report_file "$report_file"

  echo -e "\n${BOLD}${BLUE}[PASO 3/3] Ejecutando 03_cnv_inference.R...${NC}"
  Rscript ./scripts/03_cnv_inference.R \
    --input "$infer_input" \
    --models_dir "$out_models" \
    --out_dir "$out_results" \
    --cov_threshold "$cov_threshold" \
    --min_votes "$min_votes" \
    --top_loci "$top_loci" \
    --report_file "$infer_report_file"

  echo -e "\n${BOLD}${GREEN}[OK] Flujo completo finalizado.${NC}"
  echo -e "  - Datos preparados      : ${CYAN}${out_data}/data_prepared.rds${NC}"
  echo -e "  - Modelos               : ${CYAN}${out_models}/${NC}"
  echo -e "  - Matriz CNV            : ${CYAN}${out_results}/cnv_matrix.tsv${NC}"
  echo -e "  - Reporte Entrenamiento : ${CYAN}${report_file}${NC}"
  echo -e "  - Reporte Inferencia    : ${CYAN}${infer_report_file}${NC}\n"
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
        echo "Uso: ./targetcnv.sh report [--models-dir models] [--data data/data_prepared.rds] [--sources data/reference/sources_citations.tsv] [--output reports/reporte_entrenamiento.html]"
        exit 0 ;;
      *) echo -e "${RED}[ERROR] Opción desconocida en report: $1${NC}"; exit 1 ;;
    esac
  done

  local out_base
  out_base="$(basename "$output_file")"
  mkdir -p "$(dirname "$output_file")"

  quarto render training_report.qmd \
    -P project_root:"$SCRIPT_DIR" \
    -P models_dir:"$models_dir" \
    -P data_path:"$data_file" \
    -P data_file:"$data_file" \
    -P sources_file:"$sources_file" \
    --output "$out_base"

  if [ -f "$SCRIPT_DIR/$out_base" ] && [ "$(cd "$(dirname "$output_file")" && pwd)/$out_base" != "$SCRIPT_DIR/$out_base" ]; then
    mv -f "$SCRIPT_DIR/$out_base" "$output_file"
  fi
  echo -e "${BOLD}${GREEN}[OK] Reporte compilado en modo de solo lectura: ${CYAN}${output_file}${NC}"
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
        echo "Uso: ./targetcnv.sh transfer [--models-dir models] [--cit-data data/citrullifolium/data_prepared.rds] [--ste-data data/stellatiglandulosum/data_prepared.rds] [--cov-threshold 50.0] [--out-dir results] [--output reports/reporte_transferencia.html]"
        exit 0 ;;
      *) echo -e "${RED}[ERROR] Opción desconocida en transfer: $1${NC}"; exit 1 ;;
    esac
  done

  print_banner
  echo -e "\n${BOLD}${BLUE}[TRANSFER] Ejecutando evaluación de transferibilidad y concordancia CNV (04_evaluate_transferability.R)...${NC}"

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

  Rscript ./scripts/04_evaluate_transferability.R "${tr_args[@]}"
}

SUBCOMMAND="${1:-help}"
shift || true

case "$SUBCOMMAND" in
  restore|setup)
    cmd_restore "$@" ;;
  dev-hydrate)
    cmd_dev_hydrate "$@" ;;
  snapshot)
    echo -e "${YELLOW}[AVISO] El antiguo subcomando 'snapshot' sobrescribía renv.lock con paquetes globales del sistema y ha sido desactivado para proteger la reproducibilidad.${NC}"
    echo -e "        - Para instalar el entorno reproducible desde renv.lock, ejecute: ${BOLD}./targetcnv.sh restore${NC}"
    echo -e "        - Si desarrolla el paquete internamente, utilice: ${BOLD}./targetcnv.sh dev-hydrate [--update-lock]${NC}"
    exit 1 ;;
  test|mwe)
    bash ./run_test_pipeline.sh "$@" ;;
  zenodo|download-models)
    cmd_zenodo "$@" ;;
  infer|inference)
    cmd_infer "$@" ;;
  transfer|transferability|report-transfer)
    cmd_transfer "$@" ;;
  prepare|prep)
    Rscript ./scripts/01_prepare_data.R "$@" ;;
  train)
    Rscript ./scripts/02_train_models.R "$@" ;;
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
    echo -e "${RED}[ERROR] Subcomando desconocido: '$SUBCOMMAND'${NC}"
    show_help
    exit 1 ;;
esac
