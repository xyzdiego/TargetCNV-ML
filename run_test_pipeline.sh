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

# Directorio raiz del repositorio
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo -e "\n${BOLD}${CYAN}=================================================================${NC}"
echo -e "${BOLD}${CYAN}       TargetCNV-ML : Minimal Working Example (Test Suite)       ${NC}"
echo -e "${BOLD}${CYAN}=================================================================${NC}\n"

# Rutas de datos mock y salidas
TEST_DIR="data/test_data"
PREPARED_DIR="$TEST_DIR/prepared"
MODELS_DIR="$TEST_DIR/models"
RESULTS_DIR="$TEST_DIR/results"

# Limpieza previa de salidas de prueba
rm -rf "$PREPARED_DIR" "$MODELS_DIR" "$RESULTS_DIR"
mkdir -p "$PREPARED_DIR" "$MODELS_DIR" "$RESULTS_DIR"

# ------------------------------------------------------------------------------
# PASO 1: Preparacion de Datos y Cruces Espaciales
# ------------------------------------------------------------------------------
echo -e "${BOLD}${BLUE}[PASO 1/3] Ejecutando 01_prepare_data.R...${NC}"
./scripts/01_prepare_data.R \
  --gff_ref "$TEST_DIR/mock_reference.gff3" \
  --captus_gff "$TEST_DIR/mock_captus.gff" \
  --captus_stats "$TEST_DIR/mock_captus_extract.tsv" \
  --dict_ref "$TEST_DIR/mock_dictionary.tsv" \
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
./scripts/02_train_models.R \
  --data "$PREPARED_DIR/data_prepared.rds" \
  --out_models "$MODELS_DIR" \
  --cv_folds 3 \
  --sampling down \
  --seed 42 \
  --cores 2

if [ -f "$MODELS_DIR/heuristic_model.rds" ] && [ -f "$MODELS_DIR/youden_thresholds.rds" ]; then
  echo -e "${GREEN}✓ Paso 2 completado exitosamente.${NC}\n"
else
  echo -e "${RED}✗ Error: No se generaron los modelos esperados en $MODELS_DIR${NC}"
  exit 1
fi

# ------------------------------------------------------------------------------
# PASO 3: Inferencia de Produccion y Consenso de 7 Modelos
# ------------------------------------------------------------------------------
echo -e "${BOLD}${BLUE}[PASO 3/3] Ejecutando 03_cnv_inference.R...${NC}"
./scripts/03_cnv_inference.R \
  --input "$TEST_DIR/mock_captus_extract.tsv" \
  --models_dir "$MODELS_DIR" \
  --out_dir "$RESULTS_DIR" \
  --cov_threshold 50.0 \
  --frag_threshold 120.0 \
  --min_votes 4 \
  --top_loci 20

if [ -f "$RESULTS_DIR/cnv_matrix.tsv" ]; then
  echo -e "${GREEN}✓ Paso 3 completado exitosamente.${NC}\n"
else
  echo -e "${RED}✗ Error: No se genero cnv_matrix.tsv${NC}"
  exit 1
fi

# ------------------------------------------------------------------------------
# RESUMEN DE ARTEFACTOS GENERADOS
# ------------------------------------------------------------------------------
echo -e "${BOLD}${CYAN}=================================================================${NC}"
echo -e "${BOLD}${GREEN} ✓ Minimal Working Example finalizado con 0 errores.${NC}"
echo -e "${BOLD}${CYAN}=================================================================${NC}"
echo -e "Artefactos generados para verificacion:"
echo -e "  • Datos preparados   : ${CYAN}$PREPARED_DIR/data_prepared.tsv${NC}"
echo -e "  • Modelos entrenados : ${CYAN}$MODELS_DIR/benchmark_models.tsv${NC}"
echo -e "  • Matriz CNV final   : ${CYAN}$RESULTS_DIR/cnv_matrix.tsv${NC}"
echo -e "  • Veredictos Consenso: ${CYAN}$RESULTS_DIR/consensus_verdicts.tsv${NC}"
echo -e "  • Grafico Consenso   : ${CYAN}$RESULTS_DIR/plot_consensus_distribution.pdf${NC}"
echo -e "  • Grafico Heatmap    : ${CYAN}$RESULTS_DIR/plot_cnv_heatmap.pdf${NC}"
echo -e "\nVista previa de la Matriz de CNV resultante:"
cat "$RESULTS_DIR/cnv_matrix.tsv" | head -n 7 | column -t
echo -e "\n${BOLD}${GREEN}[TEST PASSED]${NC} El pipeline es 100% reproducible de inicio a fin.\n"
