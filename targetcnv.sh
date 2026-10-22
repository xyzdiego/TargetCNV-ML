#!/usr/bin/env bash

# ==============================================================================
# TargetCNV-ML: Unified Command-Line Orchestrator (targetcnv.sh)
# Ejecuta flujos completos o individuales en una sola linea de Bash
# ==============================================================================

set -euo pipefail

# Colores ANSI
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

print_banner() {
  echo -e "${BOLD}${CYAN}======================================================================${NC}"
  echo -e "${BOLD}${CYAN}  TargetCNV-ML | Automated CLI & Spatial ML Pipeline for CNV Calling  ${NC}"
  echo -e "${BOLD}${CYAN}======================================================================${NC}"
}

show_help() {
  print_banner
  echo -e "

${BOLD}USO / USAGE:${NC}
  ./targetcnv.sh <subcomando> [opciones]

${BOLD}SUBCOMANDOS RAPIDOS (1 LINEA EN BASH):${NC}
  ${GREEN}test${NC} [--sync-root]
      Ejecuta el Ejemplo Minimo Funcional (MWE) con datos simulados en segundos.

  ${GREEN}zenodo${NC} [--out-dir DIR]
      Descarga automaticamente los modelos preentrenados oficiales desde Zenodo
      (DOI: 10.5281/zenodo.23463987), verifica el hash MD5 y descomprime en models/.

  ${GREEN}infer${NC} -i <captus_stats.tsv> [--auto-zenodo] [opciones de 03_cnv_inference.R]
      Ejecuta inferencia de CNV con el consenso de 7 modelos sobre muestras nuevas.
      Si se pasa --auto-zenodo y models/ esta vacio, descarga los modelos primero.

  ${GREEN}reproduce-storvum${NC} [--cores N]
      Reproduce en 1 comando el experimento de referencia completo con Solanum torvum
      (Paso 01 -> Paso 02 -> Paso 03 + Reporte HTML Quarto).

  ${GREEN}run-all${NC} -g <ref.gff3> -c <captus.gff> -s <stats.tsv> [-d <dict.xlsx>] [-i <infer.tsv>] [opciones]
      Ejecuta el pipeline completo de extremo a extremo (01 -> 02 -> 03) en una sola orden.

  ${GREEN}prepare${NC} [opciones de 01_prepare_data.R]
      Ejecuta directamente el Modulo 01 (Interseccion espacial y preparacion).

  ${GREEN}train${NC} [opciones de 02_train_models.R]
      Ejecuta directamente el Modulo 02 (Torneo ML, Youden, SHAP y reporte Quarto).

  ${GREEN}report${NC} [--models-dir DIR] [--data FILE] [--output FILE]
      Recompila manualmente el reporte interactivo HTML con Quarto sin reentrenar.

  ${GREEN}doctor${NC}
      Verifica la instalacion de R, renv, Quarto, herramientas del sistema y paquetes.

${BOLD}EJEMPLOS / EXAMPLES:${NC}
  # 1. Validar instalacion con el MWE:
  ./targetcnv.sh test

  # 2. Descargar modelos de Zenodo e inferir CNV en muestras nuevas en 1 linea:
  ./targetcnv.sh infer -i data/CAPTUS_Muestras/captus-extract_stats.tsv --auto-zenodo

  # 3. Reproducir experimento completo de Solanum torvum:
  ./targetcnv.sh reproduce-storvum --cores 4

  # 4. Entrenar e inferir con una nueva especie en 1 sola linea:
  ./targetcnv.sh run-all \\
    -g ruta/genoma.gff3 -c ruta/contigs.gff -s ruta/captus_stats.tsv \\
    -d ruta/diccionario.xlsx -p MiMuestra -i ruta/muestras_poblacion.tsv
"
}

cmd_doctor() {
  print_banner
  echo -e "\n${BOLD}${BLUE}[DOCTOR] Verificando entorno del sistema y dependencias R...${NC}\n"
  local status=0

  for tool in Rscript curl unzip md5sum; do
    if command -v "$tool" >/dev/null 2>&1; then
      echo -e "  ${GREEN}✓${NC} Herramienta encontrada : ${BOLD}$tool${NC} ($(command -v "$tool"))"
    else
      echo -e "  ${RED}✗${NC} Herramienta faltante   : ${BOLD}$tool${NC}"
      status=1
    fi
  done

  if command -v quarto >/dev/null 2>&1; then
    echo -e "  ${GREEN}✓${NC} Quarto CLI             : ${BOLD}$(quarto --version 2>/dev/null || echo 'Instalado')${NC}"
  else
    echo -e "  ${YELLOW}!${NC} Quarto CLI no detectado (opcional: necesario solo para compilar el reporte HTML)"
  fi

  echo -e "\n${BOLD}${BLUE}[DOCTOR] Verificando paquetes clave en R...${NC}"
  Rscript -e '
    pkgs <- c("optparse", "dplyr", "tidyr", "stringr", "ggplot2", "caret", "glmnet",
              "klaR", "pROC", "GenomicRanges", "rtracklayer", "readxl", "kernelshap",
              "shapviz", "DT", "plotly")
    missing_pkgs <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) == 0) {
      cat("  \033[0;32m✓ Todos los paquetes principales de R estan instalados.\033[0m\n")
    } else {
      cat("  \033[1;33m! Paquetes faltantes:\033[0m", paste(missing_pkgs, collapse = ", "), "\n")
      cat("    Ejecute: Rscript -e \"renv::restore()\"\n")
      quit(status = 1)
    }
  ' || status=1

  echo ""
  if [ "$status" -eq 0 ]; then
    echo -e "${BOLD}${GREEN}[OK] El entorno esta listo para ejecutar TargetCNV-ML.${NC}\n"
  else
    echo -e "${BOLD}${YELLOW}[AVISO] Revise los elementos marcados arriba antes de ejecutar en produccion.${NC}\n"
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
      *) echo -e "${RED}[ERROR] Opcion desconocida en zenodo: $1${NC}"; exit 1 ;;
    esac
  done

  print_banner
  echo -e "\n${BOLD}${BLUE}[ZENODO] Descargando modelos preentrenados oficiales (DOI: 10.5281/zenodo.23463987)...${NC}"
  if command -v curl >/dev/null 2>&1; then
    curl -L --fail --progress-bar -o "$ZENODO_ZIP" "$ZENODO_URL"
  elif command -v wget >/dev/null 2>&1; then
    wget -O "$ZENODO_ZIP" "$ZENODO_URL"
  else
    echo -e "${RED}[ERROR] Se requiere curl o wget para descargar desde Zenodo.${NC}"
    exit 1
  fi

  if command -v md5sum >/dev/null 2>&1; then
    echo -e "${BOLD}${BLUE}[ZENODO] Verificando suma de comprobacion MD5...${NC}"
    local actual_md5
    actual_md5="$(md5sum "$ZENODO_ZIP" | awk '{print $1}')"
    if [ "$actual_md5" = "$ZENODO_MD5" ]; then
      echo -e "  ${GREEN}✓ Hash MD5 verificado (${actual_md5})${NC}"
    else
      echo -e "  ${RED}✗ Error de integridad MD5: esperado ${ZENODO_MD5}, recibido ${actual_md5}${NC}"
      exit 1
    fi
  fi

  mkdir -p "$out_dir"
  unzip -o "$ZENODO_ZIP" -d "$out_dir"
  echo -e "${GREEN}✓ Modelos extraidos exitosamente en: ${BOLD}${out_dir}/models/${NC}\n"
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
        models_dir="$2"
        pass_args+=("$1" "$2")
        shift 2 ;;
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
      *) echo -e "${RED}[ERROR] Opcion desconocida: $1${NC}"; exit 1 ;;
    esac
  done

  print_banner
  echo -e "\n${BOLD}${BLUE}[STORVUM 1/3] Preparacion de datos e interseccion espacial...${NC}"
  Rscript ./scripts/01_prepare_data.R \
    --gff_ref data/reference/Storvum_genes.gff3 \
    --captus_gff data/reference/Storvum_hit_contigs.gff \
    --captus_stats data/reference/captus-extract_stats.tsv \
    --dict_ref data/reference/Genes_Capture_May_2026.xlsx \
    --sample Storvum \
    --out_dir data/

  echo -e "\n${BOLD}${BLUE}[STORVUM 2/3] Entrenamiento de modelos, Youden, SHAP y Reporte HTML...${NC}"
  Rscript ./scripts/02_train_models.R \
    --data data/data_prepared.rds \
    --out_models models/ \
    --cv_folds 5 \
    --sampling down \
    --cores "$cores" \
    --seed 1234 \
    --report_file reporte_entrenamiento.html

  echo -e "\n${BOLD}${BLUE}[STORVUM 3/3] Inferencia multimuestra y extraccion de matriz CNV...${NC}"
  Rscript ./scripts/03_cnv_inference.R \
    --input data/CAPTUS_Muestras/captus-extract_stats.tsv \
    --models_dir models/ \
    --out_dir results/ \
    --cov_threshold 50.0 \
    --frag_threshold 120.0 \
    --min_votes 4 \
    --top_loci 30

  echo -e "\n${BOLD}${GREEN}✓ Experimento de referencia Solanum torvum reproducido con exito.${NC}\n"
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
  local frag_threshold="120.0"
  local top_loci=30
  local report_file="reporte_entrenamiento.html"

  while [[ $# -gt 0 ]]; do
    case "$1" in
      -g|--gff_ref) gff_ref="$2"; shift 2 ;;
      -c|--captus_gff) captus_gff="$2"; shift 2 ;;
      -s|--captus_stats) captus_stats="$2"; shift 2 ;;
      -d|--dict_ref) dict_ref="$2"; shift 2 ;;
      -p|--sample) sample_name="$2"; shift 2 ;;
      --target_mode) target_mode="$2"; shift 2 ;;
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
      --frag_threshold) frag_threshold="$2"; shift 2 ;;
      -t|--top_loci) top_loci="$2"; shift 2 ;;
      -r|--report_file) report_file="$2"; shift 2 ;;
      -h|--help)
        echo "Uso: ./targetcnv.sh run-all -g <ref.gff3> -c <captus.gff> -s <stats.tsv> [-d <dict>] [-p <sample>] [-i <infer.tsv>]"
        exit 0 ;;
      *) echo -e "${RED}[ERROR] Argumento desconocido en run-all: $1${NC}"; exit 1 ;;
    esac
  done

  if [ -z "$captus_gff" ] || [ -z "$captus_stats" ]; then
    echo -e "${RED}[ERROR] 'run-all' requiere al menos -c (--captus_gff) y -s (--captus_stats).${NC}"
    exit 1
  fi
  if [ -z "$infer_input" ]; then
    infer_input="$captus_stats"
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
    --frag_threshold "$frag_threshold" \
    --min_votes "$min_votes" \
    --top_loci "$top_loci"

  echo -e "\n${BOLD}${GREEN}✓ Pipeline completo finalizado exitosamente.${NC}"
  echo -e "  • Datos preparados : ${CYAN}${out_data}/data_prepared.rds${NC}"
  echo -e "  • Modelos          : ${CYAN}${out_models}/${NC}"
  echo -e "  • Resultados CNV   : ${CYAN}${out_results}/cnv_matrix.tsv${NC}"
  echo -e "  • Reporte HTML     : ${CYAN}${report_file}${NC}\n"
}

cmd_report() {
  local models_dir="models"
  local data_file="data/data_prepared.rds"
  local output_file="reporte_entrenamiento.html"

  while [[ $# -gt 0 ]]; do
    case "$1" in
      -m|--models-dir) models_dir="$2"; shift 2 ;;
      -d|--data) data_file="$2"; shift 2 ;;
      -o|--output) output_file="$2"; shift 2 ;;
      -h|--help)
        echo "Uso: ./targetcnv.sh report [--models-dir models] [--data data/data_prepared.rds] [--output reporte_entrenamiento.html]"
        exit 0 ;;
      *) echo -e "${RED}[ERROR] Opcion desconocida en report: $1${NC}"; exit 1 ;;
    esac
  done

  quarto render training_report.qmd \
    -P project_root:"$SCRIPT_DIR" \
    -P models_dir:"$models_dir" \
    -P data_path:"$data_file" \
    -P data_file:"$data_file" \
    --output "$(basename "$output_file")"
}

# Enrutador principal de subcomandos
SUBCOMMAND="${1:-help}"
shift || true

case "$SUBCOMMAND" in
  test|mwe)
    bash ./run_test_pipeline.sh "$@" ;;
  zenodo|download-models)
    cmd_zenodo "$@" ;;
  infer|inference)
    cmd_infer "$@" ;;
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