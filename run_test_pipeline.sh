#!/usr/bin/env bash

# TargetCNV-ML: Ejemplo Mínimo Funcional (Minimal Working Example)
# Ejecuta secuencialmente los módulos 01, 02 y 03 sobre el conjunto sintético de validación.

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

SYNC_ROOT=0
NO_REPORT=0
for arg in "$@"; do
  case "$arg" in
    --sync-root) SYNC_ROOT=1 ;;
    --no-report) NO_REPORT=1 ;;
    -h|--help)
      echo "Uso: ./run_test_pipeline.sh [--sync-root] [--no-report]"
      echo "  --sync-root  Copia los artefactos del MWE a models/ y results/ en la raíz (operación destructiva: sobrescribe archivos existentes)"
      echo "  --no-report  Omite la compilación de reportes HTML Quarto durante la prueba rápida"
      exit 0 ;;
    *) echo "Opción desconocida: $arg (use --help)"; exit 1 ;;
  esac
done

START_TS=$(date +%s)

echo -e "\n${BOLD}${CYAN}====================================================================${NC}"
echo -e "${BOLD}${CYAN} TargetCNV-ML : Minimal Working Example (MWE) ${NC}"
echo -e "${BOLD}${CYAN}====================================================================${NC}\n"

TEST_DIR="data/test_data"
PREPARED_DIR="$TEST_DIR/prepared"
MODELS_DIR="$TEST_DIR/models"
RESULTS_DIR="$TEST_DIR/results"

rm -rf "$PREPARED_DIR" "$MODELS_DIR" "$RESULTS_DIR"
mkdir -p "$PREPARED_DIR" "$MODELS_DIR" "$RESULTS_DIR"

# Paso 1: Preparación de datos y cruce espacial
echo -e "${BOLD}${BLUE}[PASO 1/3] Ejecutando 01_prepare_data.R...${NC}"
Rscript ./scripts/01_prepare_data.R \
  --gff_ref "$TEST_DIR/mock_reference.gff3" \
  --captus_gff "$TEST_DIR/mock_captus.gff" \
  --captus_stats "$TEST_DIR/mock_captus_extract.tsv" \
  --dict_ref "$TEST_DIR/mock_dictionary.tsv" \
  --allow_multi_sample \
  --out_dir "$PREPARED_DIR"

if [ -f "$PREPARED_DIR/data_prepared.rds" ]; then
  echo -e "${GREEN}[OK] Paso 1 completado.${NC}\n"
else
  echo -e "${RED}[ERROR] No se generó data_prepared.rds${NC}"
  exit 1
fi

REPORTS_DIR="reports"
mkdir -p "$REPORTS_DIR"

# Paso 2: Entrenamiento supervisado y validación cruzada determinista
echo -e "${BOLD}${BLUE}[PASO 2/3] Ejecutando 02_train_models.R...${NC}"
REPORT_ARGS=(--report_file "$REPORTS_DIR/reporte_entrenamiento_test.html")
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

if [ -f "$MODELS_DIR/heuristic_model.rds" ] && [ -f "$MODELS_DIR/youden_thresholds.rds" ] && [ -f "$MODELS_DIR/spearman_collinearity.tsv" ]; then
  echo -e "${GREEN}[OK] Paso 2 completado.${NC}\n"
else
  echo -e "${RED}[ERROR] Faltan artefactos esperados en $MODELS_DIR${NC}"
  exit 1
fi

# Paso 3: Inferencia de CNV y calibración de cobertura
echo -e "${BOLD}${BLUE}[PASO 3/3] Ejecutando 03_cnv_inference.R...${NC}"
INF_REPORT_ARGS=(--report_file "$REPORTS_DIR/reporte_inferencia_test.html")
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

if [ -f "$RESULTS_DIR/cnv_matrix.tsv" ] && [ -f "$RESULTS_DIR/consensus_verdicts.tsv" ] && [ -f "$RESULTS_DIR/cnv_threshold_sensitivity.tsv" ]; then
  echo -e "${GREEN}[OK] Paso 3 completado.${NC}\n"
else
  echo -e "${RED}[ERROR] No se generó la matriz CNV o la tabla de sensibilidad en $RESULTS_DIR${NC}"
  exit 1
fi

if [ "$SYNC_ROOT" -eq 1 ]; then
  echo -e "${YELLOW}[AVISO] Copiando artefactos del MWE a models/, results/ y reports/ (--sync-root)...${NC}"
  mkdir -p models results "$REPORTS_DIR"
  cp -f "$MODELS_DIR"/* models/ 2>/dev/null || true
  cp -f "$RESULTS_DIR"/* results/ 2>/dev/null || true
  if [ -f "$REPORTS_DIR/reporte_entrenamiento_test.html" ]; then
    cp -f "$REPORTS_DIR/reporte_entrenamiento_test.html" "$REPORTS_DIR/reporte_entrenamiento.html"
  fi
  if [ -f "$REPORTS_DIR/reporte_inferencia_test.html" ]; then
    cp -f "$REPORTS_DIR/reporte_inferencia_test.html" "$REPORTS_DIR/reporte_inferencia.html"
  fi
  echo -e "${GREEN}[OK] Sincronización con la raíz completada.${NC}\n"
fi

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))

echo -e "${BOLD}${CYAN}====================================================================${NC}"
echo -e "${BOLD}${GREEN}[OK] Ejemplo Mínimo Funcional finalizado sin errores (${ELAPSED}s).${NC}"
echo -e "${BOLD}${CYAN}====================================================================${NC}"
echo -e "Artefactos generados para verificación:"
echo -e "  - Datos preparados     : ${CYAN}$PREPARED_DIR/data_prepared.tsv${NC}"
echo -e "  - Colinealidad Spearman: ${CYAN}$MODELS_DIR/spearman_collinearity.tsv${NC}"
echo -e "  - Modelos entrenados   : ${CYAN}$MODELS_DIR/benchmark_models.tsv${NC}"
echo -e "  - Matriz CNV final     : ${CYAN}$RESULTS_DIR/cnv_matrix.tsv${NC}"
echo -e "  - Sensibilidad umbral  : ${CYAN}$RESULTS_DIR/cnv_threshold_sensitivity.tsv${NC}"
echo -e "  - Veredictos consenso  : ${CYAN}$RESULTS_DIR/consensus_verdicts.tsv${NC}"
echo -e "\nVista previa de la matriz de CNV resultante:"
if command -v column >/dev/null 2>&1; then
  head -n 7 "$RESULTS_DIR/cnv_matrix.tsv" | column -t
else
  head -n 7 "$RESULTS_DIR/cnv_matrix.tsv"
fi
echo -e "\n${BOLD}${GREEN}[TEST PASSED]${NC} Ejecución reproducible verificada de principio a fin.\n"
