#!/usr/bin/env bash

# Ejemplo Mínimo Funcional (MWE) de TargetCNV-ML: ejecuta 01 -> 02 -> 03 sobre data/test_data/.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

SYNC_ROOT=0
NO_REPORT=0
for arg in "$@"; do
  case "$arg" in
    --sync-root) SYNC_ROOT=1 ;;
    --no-report) NO_REPORT=1 ;;
    -h|--help)
      echo "Uso: ./run_test_pipeline.sh [--sync-root] [--no-report]"
      echo "  --sync-root  Copia las salidas del MWE a models/ y results/ en la raíz (sobrescribe archivos previos)"
      echo "  --no-report  Omite la compilación de los reportes HTML Quarto"
      exit 0 ;;
    *) echo "[ERROR] Opción desconocida: $arg (use --help)" >&2; exit 1 ;;
  esac
done

START_TS=$(date +%s)

TEST_DIR="$SCRIPT_DIR/data/test_data"
if [ -w "$TEST_DIR" ]; then
  WORK_BASE="$TEST_DIR"
  REPORTS_DIR="$SCRIPT_DIR/reports"
else
  # Compatibilidad con imágenes de contenedor en sistema de archivos de solo lectura
  WORK_BASE="$PWD/targetcnv_mwe_out"
  REPORTS_DIR="$WORK_BASE/reports"
fi

PREPARED_DIR="$WORK_BASE/prepared"
MODELS_DIR="$WORK_BASE/models"
RESULTS_DIR="$WORK_BASE/results"

rm -rf "$PREPARED_DIR" "$MODELS_DIR" "$RESULTS_DIR"
mkdir -p "$PREPARED_DIR" "$MODELS_DIR" "$RESULTS_DIR" "$REPORTS_DIR"

echo "[MWE 1/3] Ejecutando 01_prepare_data.R..."
Rscript "$SCRIPT_DIR/scripts/01_prepare_data.R" \
  --gff_ref "$TEST_DIR/mock_reference.gff3" \
  --captus_gff "$TEST_DIR/mock_captus.gff" \
  --captus_stats "$TEST_DIR/mock_captus_extract.tsv" \
  --dict_ref "$TEST_DIR/mock_dictionary.tsv" \
  --allow_multi_sample \
  --out_dir "$PREPARED_DIR"

[ -f "$PREPARED_DIR/data_prepared.rds" ] || { echo "[ERROR] No se generó data_prepared.rds" >&2; exit 1; }

echo "[MWE 2/3] Ejecutando 02_train_models.R..."
REPORT_ARGS=(--report_file "$REPORTS_DIR/reporte_entrenamiento_test.html")
if [ "$NO_REPORT" -eq 1 ]; then
  REPORT_ARGS=(--no_render_report)
fi

Rscript "$SCRIPT_DIR/scripts/02_train_models.R" \
  --data "$PREPARED_DIR/data_prepared.rds" \
  --out_models "$MODELS_DIR" \
  --cv_folds 3 \
  --sampling down \
  --seed 42 \
  --cores 2 \
  --fast_tournament \
  --no_stacking \
  "${REPORT_ARGS[@]}"

if [ ! -f "$MODELS_DIR/heuristic_model.rds" ] || [ ! -f "$MODELS_DIR/youden_thresholds.rds" ] || [ ! -f "$MODELS_DIR/spearman_collinearity.tsv" ]; then
  echo "[ERROR] Faltan artefactos esperados en $MODELS_DIR" >&2
  exit 1
fi

echo "[MWE 3/3] Ejecutando 03_cnv_inference.R..."
INF_REPORT_ARGS=(--report_file "$REPORTS_DIR/reporte_inferencia_test.html")
if [ "$NO_REPORT" -eq 1 ]; then
  INF_REPORT_ARGS=(--no_render_report)
fi

Rscript "$SCRIPT_DIR/scripts/03_cnv_inference.R" \
  --input "$TEST_DIR/mock_captus_extract.tsv" \
  --models_dir "$MODELS_DIR" \
  --out_dir "$RESULTS_DIR" \
  --cov_threshold 50.0 \
  --min_votes 2 \
  --top_loci 10 \
  "${INF_REPORT_ARGS[@]}"

if [ ! -f "$RESULTS_DIR/cnv_matrix.tsv" ] || [ ! -f "$RESULTS_DIR/consensus_verdicts.tsv" ] || [ ! -f "$RESULTS_DIR/cnv_threshold_sensitivity.tsv" ]; then
  echo "[ERROR] Faltan salidas esperadas en $RESULTS_DIR" >&2
  exit 1
fi

if [ "$SYNC_ROOT" -eq 1 ]; then
  echo "[MWE] Sincronizando salidas con models/, results/ y reports/ (--sync-root)..."
  mkdir -p models results "$REPORTS_DIR"
  cp -f "$MODELS_DIR"/* models/ 2>/dev/null || true
  cp -f "$RESULTS_DIR"/* results/ 2>/dev/null || true
  if [ -f "$REPORTS_DIR/reporte_entrenamiento_test.html" ]; then
    cp -f "$REPORTS_DIR/reporte_entrenamiento_test.html" "$REPORTS_DIR/reporte_entrenamiento.html"
  fi
  if [ -f "$REPORTS_DIR/reporte_inferencia_test.html" ]; then
    cp -f "$REPORTS_DIR/reporte_inferencia_test.html" "$REPORTS_DIR/reporte_inferencia.html"
  fi
fi

ELAPSED=$(( $(date +%s) - START_TS ))
echo ""
echo "[MWE OK] Flujo de prueba completado en ${ELAPSED}s."
echo "  datos:        $PREPARED_DIR/data_prepared.tsv"
echo "  spearman:     $MODELS_DIR/spearman_collinearity.tsv"
echo "  benchmark:    $MODELS_DIR/benchmark_models.tsv"
echo "  cnv_matrix:   $RESULTS_DIR/cnv_matrix.tsv"
echo "  sensibilidad: $RESULTS_DIR/cnv_threshold_sensitivity.tsv"
echo ""
echo "Vista previa de cnv_matrix.tsv:"
if command -v column >/dev/null 2>&1; then
  head -n 7 "$RESULTS_DIR/cnv_matrix.tsv" | column -t
else
  head -n 7 "$RESULTS_DIR/cnv_matrix.tsv"
fi
