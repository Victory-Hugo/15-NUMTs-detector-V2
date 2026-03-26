#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash pipe/01-run_numts_frequency.sh --config conf/Config.yaml
EOF
}

if [[ $# -ne 2 || "$1" != "--config" ]]; then
  usage >&2
  exit 1
fi

CONFIG_PATH="$2"
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

source "$PROJECT_DIR/script/load_config.sh"
load_config "$CONFIG_PATH"

CONDA_ROOT="${HOME}/miniconda3"
CONDA_SH="${CONDA_ROOT}/etc/profile.d/conda.sh"
if [[ ! -f "$CONDA_SH" ]]; then
  echo "conda.sh not found: $CONDA_SH" >&2
  exit 1
fi

set +u
source "$CONDA_SH"
conda activate "${ENV_CONDA_ENV}"
set -u

PYTHON_BIN="$(command -v python)"
if [[ -z "$PYTHON_BIN" ]]; then
  echo "Python not found in env: ${ENV_CONDA_ENV}" >&2
  exit 1
fi

export PYTHONPATH="${PROJECT_DIR}/python:${PYTHONPATH:-}"

mkdir -p \
  "${PROJECT_DIR}/output/00-qc" \
  "${PROJECT_DIR}/output/cache" \
  "${PROJECT_DIR}/output/01-sample-events" \
  "${PROJECT_DIR}/output/02-catalog" \
  "${PROJECT_DIR}/output/03-annotation" \
  "${PROJECT_DIR}/output/04-main-results/tables" \
  "${PROJECT_DIR}/output/04-main-results/plots" \
  "${PROJECT_DIR}/output/05-mt-breakpoint/tables" \
  "${PROJECT_DIR}/output/05-mt-breakpoint/plots" \
  "${PROJECT_DIR}/output/06-stratified"

LENGTH_SUMMARY_CACHE_TSV="${PROJECT_DIR}/output/cache/length_bed.full_region_summary.tsv"
STEP1_REQUIRED_FILES=(
  "${PROJECT_DIR}/output/00-qc/01-input_summary.tsv"
  "${PROJECT_DIR}/output/00-qc/02-sample_inclusion.tsv"
  "${PROJECT_DIR}/output/00-qc/03-sample_exclusion.tsv"
  "${PROJECT_DIR}/output/00-qc/04-event_traceability.tsv"
  "${PROJECT_DIR}/output/00-qc/05-input_validation_report.md"
  "${PROJECT_DIR}/output/00-qc/_length_bed_region_summary.tsv"
)

step1_cache_ready=true
for required_file in "${STEP1_REQUIRED_FILES[@]}"; do
  if [[ ! -s "$required_file" ]]; then
    step1_cache_ready=false
    break
  fi
done

if [[ "$step1_cache_ready" == true ]]; then
  echo "[1/9] validate inputs (cached 00-qc)"
else
  echo "[1/9] validate inputs"
  "$PYTHON_BIN" "${PROJECT_DIR}/python/01_validate_inputs.py" \
    --breakpoint-tsv "${INPUT_BREAKPOINT_TSV}" \
    --cluster-tsv "${INPUT_CLUSTER_TSV}" \
    --length-bed-gz "${INPUT_LENGTH_BED_GZ}" \
    --meta-tsv "${INPUT_META_TSV}" \
    --meta-id-col "${QC_META_ID_COL}" \
    --meta-qc-col "${QC_META_QC_COL}" \
    --meta-qc-pass-value "${QC_META_QC_PASS_VALUE}" \
    --output-dir "${PROJECT_DIR}/output/00-qc" \
    --length-summary-cache-tsv "${LENGTH_SUMMARY_CACHE_TSV}"
fi

echo "[2/9] build sample event table"
"$PYTHON_BIN" "${PROJECT_DIR}/python/02_build_sample_event_table.py" \
  --breakpoint-tsv "${INPUT_BREAKPOINT_TSV}" \
  --cluster-tsv "${INPUT_CLUSTER_TSV}" \
  --length-summary-tsv "${PROJECT_DIR}/output/00-qc/_length_bed_region_summary.tsv" \
  --sample-inclusion-tsv "${PROJECT_DIR}/output/00-qc/02-sample_inclusion.tsv" \
  --output-tsv "${PROJECT_DIR}/output/01-sample-events/sample_numt_events.tsv"

echo "[3/9] build distinct catalog"
"$PYTHON_BIN" "${PROJECT_DIR}/python/03_build_distinct_catalog.py" \
  --sample-events-tsv "${PROJECT_DIR}/output/01-sample-events/sample_numt_events.tsv" \
  --sample-inclusion-tsv "${PROJECT_DIR}/output/00-qc/02-sample_inclusion.tsv" \
  --nuclear-merge-gap-bp "${CATALOG_NUCLEAR_MERGE_GAP_BP}" \
  --mt-merge-gap-bp "${CATALOG_MT_MERGE_GAP_BP}" \
  --distinct-id-prefix "${CATALOG_DISTINCT_ID_PREFIX}" \
  --prevalent-min "${FREQUENCY_PREVALENT_MIN}" \
  --common-min "${FREQUENCY_COMMON_MIN}" \
  --rare-min "${FREQUENCY_RARE_MIN}" \
  --catalog-tsv "${PROJECT_DIR}/output/02-catalog/distinct_numt_catalog.tsv" \
  --sample-map-tsv "${PROJECT_DIR}/output/02-catalog/distinct_numt_sample_map.tsv"

echo "[4/9] annotate catalog"
"$PYTHON_BIN" "${PROJECT_DIR}/python/04_annotate_catalog.py" \
  --sample-events-tsv "${PROJECT_DIR}/output/01-sample-events/sample_numt_events.tsv" \
  --catalog-tsv "${PROJECT_DIR}/output/02-catalog/distinct_numt_catalog.tsv" \
  --sample-map-tsv "${PROJECT_DIR}/output/02-catalog/distinct_numt_sample_map.tsv" \
  --length-summary-tsv "${PROJECT_DIR}/output/00-qc/_length_bed_region_summary.tsv" \
  --mt-gene-table "${INPUT_MT_GENE_TABLE}" \
  --sample-events-output-tsv "${PROJECT_DIR}/output/03-annotation/sample_numt_events.annotated.tsv" \
  --catalog-output-tsv "${PROJECT_DIR}/output/03-annotation/distinct_numt_catalog.annotated.tsv"

echo "[5/9] compute main results"
"$PYTHON_BIN" "${PROJECT_DIR}/python/05_compute_main_results.py" \
  --catalog-tsv "${PROJECT_DIR}/output/03-annotation/distinct_numt_catalog.annotated.tsv" \
  --sample-map-tsv "${PROJECT_DIR}/output/02-catalog/distinct_numt_sample_map.tsv" \
  --sample-events-tsv "${PROJECT_DIR}/output/03-annotation/sample_numt_events.annotated.tsv" \
  --sample-inclusion-tsv "${PROJECT_DIR}/output/00-qc/02-sample_inclusion.tsv" \
  --output-dir "${PROJECT_DIR}/output/04-main-results/tables"

echo "[6/9] compute mt breakpoint results"
"$PYTHON_BIN" "${PROJECT_DIR}/python/06_compute_mt_breakpoint_results.py" \
  --catalog-tsv "${PROJECT_DIR}/output/03-annotation/distinct_numt_catalog.annotated.tsv" \
  --output-dir "${PROJECT_DIR}/output/05-mt-breakpoint/tables"

echo "[7/9] compute stratified results"
"$PYTHON_BIN" "${PROJECT_DIR}/python/07_compute_stratified_results.py" \
  --catalog-tsv "${PROJECT_DIR}/output/03-annotation/distinct_numt_catalog.annotated.tsv" \
  --sample-map-tsv "${PROJECT_DIR}/output/02-catalog/distinct_numt_sample_map.tsv" \
  --sample-events-tsv "${PROJECT_DIR}/output/03-annotation/sample_numt_events.annotated.tsv" \
  --sample-inclusion-tsv "${PROJECT_DIR}/output/00-qc/02-sample_inclusion.tsv" \
  --meta-tsv "${INPUT_META_TSV}" \
  --meta-id-col "${QC_META_ID_COL}" \
  --geography-col "${STRATIFY_GEOGRAPHY_COL}" \
  --ethnicity-col "${STRATIFY_ETHNICITY_COL}" \
  --sex-col "${STRATIFY_SEX_COL}" \
  --prevalent-min "${FREQUENCY_PREVALENT_MIN}" \
  --common-min "${FREQUENCY_COMMON_MIN}" \
  --rare-min "${FREQUENCY_RARE_MIN}" \
  --min-group-n-warn "${STRATIFY_MIN_GROUP_N_WARN}" \
  --output-dir "${PROJECT_DIR}/output/06-stratified"

echo "[8/9] render plots"
"$PYTHON_BIN" "${PROJECT_DIR}/python/08_render_plots.py" \
  --main-results-dir "${PROJECT_DIR}/output/04-main-results" \
  --mt-results-dir "${PROJECT_DIR}/output/05-mt-breakpoint" \
  --stratified-dir "${PROJECT_DIR}/output/06-stratified" \
  --format-main "${PLOT_FORMAT_MAIN}" \
  --format-aux "${PLOT_FORMAT_AUX}" \
  --dpi "${PLOT_DPI}" \
  --font-family "${PLOT_FONT_FAMILY}"

echo "[9/9] validate outputs"
"$PYTHON_BIN" "${PROJECT_DIR}/python/09_validate_outputs.py" \
  --main-results-dir "${PROJECT_DIR}/output/04-main-results/tables" \
  --mt-results-dir "${PROJECT_DIR}/output/05-mt-breakpoint/tables" \
  --stratified-dir "${PROJECT_DIR}/output/06-stratified" \
  --sample-events-tsv "${PROJECT_DIR}/output/03-annotation/sample_numt_events.annotated.tsv" \
  --report-path "${PROJECT_DIR}/output/00-qc/06-output_validation_report.md"

echo "Pipeline completed"
