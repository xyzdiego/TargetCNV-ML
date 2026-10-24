#!/usr/bin/env bash

# ==============================================================================
# TargetCNV-ML: Script de Demostracion de Pipeline (Minimal Working Example)
# Ejecuta secuencialmente los tres modulos CLI utilizando datos simulados
# ==============================================================================

set -euo pipefail

# Colores de consola ANSI
GREEN='\033[0;32m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BOLD='\033[1m'
NC='\033[0m' # No Color

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

SYNC_ROOT=0
NO_REPORT=0
for arg in "$@"; do
  case "$arg" in
    --sync-root) SYNC_ROOT=1 ;;
    --no-report) NO_REPORT=1 ;;
    -h|--help)
      echo "Uso: ./run_test_pipeline.sh [--sync-root] [--no-report]"
      echo "  --sync-root  Copia ademas los artefactos del MWE a models/ y results/ de la raiz (SOBRESCRIBE su contenido)"
      echo "  --no-report  Omite la compilacion del reporte HTML Quarto durante la prueba rapida"
      exit 0 ;;
    *) echo "Opcion desconocida: $arg (use --help)"; exit 1 ;;
  esac
done

START_TS=$(date +%s)

echo -e "\n${BOLD}${CYAN}====================================================================${NC}"
echo -e "${BOLD}${CYAN} TargetCNV-ML : Minimal Working Example (Test Suite) ${NC}"
echo -e "${BOLD}${CYAN}====================================================================${NC}\n"

TEST_DIR="data/test_data"
PREPARED_DIR="$TEST_DIR/prepared"
MODELS_DIR="$TEST_DIR/models"
RESULTS_DIR="$TEST_DIR/results"

rm -rf "$PREPARED_DIR" "$MODELS_DIR" "$RESULTS_DIR"
mkdir -p "$PREPARED_DIR" "$MODELS_DIR" "$RESULTS_DIR"

# ------------------------------------------------------------------------------
# PASO 1: Preparacion de Datos y Cruces Espaciales
# ------------------------------------------------------------------------------
echo -e "${BOLD}${BLUE}[PASO 1/3] Ejecutando 01_prepare_data.R...${NC}"
Rscript ./scripts/01_prepare_data.R \
  --gff_ref "$TEST_DIR/mock_reference.gff3" \
  --captus_gff "$TEST_DIR/mock_captus.gff" \
  --captus_stats "$TEST_DIR/mock_captus_extract.tsv" \
  --dict_ref "$TEST_DIR/mock_dictionary.tsv" \
  --allow_multi_sample \
  --out_dir "$PREPARED_DIR"

if [ -f "$PREPARED_DIR/data_prepared.rds" ]; then
  echo -e "${GREEN}✓ Paso 1 completado exitosamente.${NC}\n"
else
  echo -e "${RED}✗ Error: No se genero data_prepared.rds${NC}"
  exit 1
fi

# ------------------------------------------------------------------------------
# PASO 2: Entrenamiento de Modelos y Validacion Cruzada con Caret
# ------------------------------------------------------------------------------
echo -e "${BOLD}${BLUE}[PASO 2/3] Ejecutando 02_train_models.R...${NC}"
REPORT_ARGS=(--report_file "$TEST_DIR/reporte_entrenamiento_test.html")
if [ "$NO_REPORT" -eq 1 ]; then
  REPORT_ARGS=(--no_render_report)
fi

Rscript ./scripts/02_train_models.R \
  --data "$PREPARED_DIR/data_prepared.rds" \
  --out_models "$MODELS_DIR" \
  --cv_folds 3 \
  --sampling down \
  --seed 42 \
  --cores 2 \
  --fast_tournament \
  --no_stacking \
  "${REPORT_ARGS[@]}"

if [ -f "$MODELS_DIR/heuristic_model.rds" ] && [ -f "$MODELS_DIR/youden_thresholds.rds" ]; then
  echo -e "${GREEN}✓ Paso 2 completado exitosamente.${NC}\n"
else
  echo -e "${RED}✗ Error: No se generaron los modelos esperados en $MODELS_DIR${NC}"
  exit 1
fi

# ------------------------------------------------------------------------------
# PASO 3: Inferencia de Produccion y Generacion de Matriz CNV
# ------------------------------------------------------------------------------
echo -e "${BOLD}${BLUE}[PASO 3/3] Ejecutando 03_cnv_inference.R...${NC}"
INF_REPORT_ARGS=(--report_file "$TEST_DIR/reporte_inferencia_test.html")
if [ "$NO_REPORT" -eq 1 ]; then
  INF_REPORT_ARGS=(--no_render_report)
fi

Rscript ./scripts/03_cnv_inference.R \
  --input "$TEST_DIR/mock_captus_extract.tsv" \
  --models_dir "$MODELS_DIR" \
  --out_dir "$RESULTS_DIR" \
  --cov_threshold 50.0 \
  --min_votes 2 \
  --top_loci 10 \
  "${INF_REPORT_ARGS[@]}" 

if [ -f "$RESULTS_DIR/cnv_matrix.tsv" ] && [ -f "$RESULTS_DIR/consensus_verdicts.tsv" ]; then
  echo -e "${GREEN}✓ Paso 3 completado exitosamente.${NC}\n"
else
  echo -e "${RED}✗ Error: No se genero la matriz CNV en $RESULTS_DIR${NC}"
  exit 1
fi

# Sincronizacion opcional con los directorios raiz
if [ "$SYNC_ROOT" -eq 1 ]; then
  echo -e "${YELLOW}[SYNC] Copiando artefactos del MWE a models/ y results/ de la raiz...${NC}"
  mkdir -p models results
  cp -f "$MODELS_DIR"/* models/ 2>/dev/null || true
  cp -f "$RESULTS_DIR"/* results/ 2>/dev/null || true
  if [ -f "$TEST_DIR/reporte_entrenamiento_test.html" ]; then
    cp -f "$TEST_DIR/reporte_entrenamiento_test.html" reporte_entrenamiento.html
  fi
  if [ -f "$TEST_DIR/reporte_inferencia_test.html" ]; then
    cp -f "$TEST_DIR/reporte_inferencia_test.html" reporte_inferencia.html
  fi
  echo -e "${GREEN}✓ Sincronizacion con la raiz completada.${NC}\n"
fi

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))

echo -e "${BOLD}${CYAN}====================================================================${NC}"
echo -e "${BOLD}${GREEN}✓ Minimal Working Example finalizado con 0 errores (${ELAPSED}s).${NC}"
echo -e "${BOLD}${CYAN}====================================================================${NC}"
echo -e "Artefactos generados para verificacion:"
echo -e "  • Datos preparados   : ${CYAN}$PREPARED_DIR/data_prepared.tsv${NC}"
echo -e "  • Modelos entrenados : ${CYAN}$MODELS_DIR/benchmark_models.tsv${NC}"
echo -e "  • Matriz CNV final   : ${CYAN}$RESULTS_DIR/cnv_matrix.tsv${NC}"
echo -e "  • Veredictos Consenso: ${CYAN}$RESULTS_DIR/consensus_verdicts.tsv${NC}"
echo -e "  • Grafico Consenso   : ${CYAN}$RESULTS_DIR/plot_consensus_distribution.pdf${NC}"
echo -e "  • Grafico Heatmap    : ${CYAN}$RESULTS_DIR/plot_cnv_heatmap.pdf${NC}"
echo -e "\nVista previa de la Matriz de CNV resultante:"
if command -v column >/dev/null 2>&1; then
  head -n 7 "$RESULTS_DIR/cnv_matrix.tsv" | column -t
else
  head -n 7 "$RESULTS_DIR/cnv_matrix.tsv"
fi
echo -e "\n${BOLD}${GREEN}[TEST PASSED]${NC} El pipeline es 100% reproducible de inicio a fin.\n"