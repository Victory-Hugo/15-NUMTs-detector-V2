#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
CONFIG_PATH="$PROJECT_ROOT/conf/Config.yaml"
PIPELINE_LOG="$PROJECT_ROOT/tmp/logs/pipeline.log"

mkdir -p "$PROJECT_ROOT/tmp/logs"
: > "$PIPELINE_LOG"

log() {
    local message="$1"
    printf '[%s] %s\n' "$(date '+%F %T')" "$message" | tee -a "$PIPELINE_LOG"
}

require_var() {
    local var_name="$1"
    if [[ -z "${!var_name:-}" ]]; then
        log "Missing required config value: $var_name"
        exit 1
    fi
}

resolve_path() {
    local raw_path="$1"
    if [[ "$raw_path" = /* ]]; then
        printf '%s\n' "$raw_path"
    else
        printf '%s\n' "$PROJECT_ROOT/$raw_path"
    fi
}

run_with_log() {
    local log_name="$1"
    shift
    local log_path="$PROJECT_ROOT/tmp/logs/${log_name}.log"
    : > "$log_path"
    log "Start $log_name"
    "$@" > >(tee -a "$log_path") 2> >(tee -a "$log_path" >&2)
    log "Done $log_name"
}

convert_svg_to_png() {
    local input_svg="$1"
    local output_png="$2"
    local log_path="$PROJECT_ROOT/tmp/logs/01_plot_ideogram.log"

    if [[ ! -f "$input_svg" ]]; then
        log "SVG not found, skip PNG conversion: $input_svg"
        return 1
    fi

    log "Convert SVG to PNG: $(basename "$input_svg")"
    rsvg-convert -a --dpi-x 300 --dpi-y 300 -o "$output_png" "$input_svg" >>"$log_path" 2>&1
}

if [[ ! -f "$CONFIG_PATH" ]]; then
    log "Config file not found: $CONFIG_PATH"
    exit 1
fi

source "$PROJECT_ROOT/script/load_config.sh" "$CONFIG_PATH"

require_var CFG_INPUT_BREAKPOINT_TSV_GZ
require_var CFG_INPUT_DISTINCT_NUMT_CLUSTER_INPUT_TSV_GZ
require_var CFG_INPUT_MERGE_BED_TSV_GZ
require_var CFG_INPUT_META_TSV
require_var CFG_META_ID_COL
require_var CFG_META_QC_COL
require_var CFG_META_QC_PASS_VALUE
require_var CFG_REFERENCE_KARYOTYPE_TXT
require_var CFG_ANALYSIS_FREQUENCY_CLUSTER_GAP_BP
require_var CFG_ANALYSIS_DISTINCT_NUMT_MIN_SAMPLE_SUPPORTS
require_var CFG_ANALYSIS_DISTINCT_NUMT_PRIMARY_MIN_SAMPLE_SUPPORT
require_var CFG_ANALYSIS_IDEOGRAM_BIN_SIZE_BP
require_var CFG_RUNTIME_CONDA_SH
require_var CFG_RUNTIME_CONDA_ENV
require_var CFG_RUNTIME_THREADS
require_var CFG_RUNTIME_CHUNK_ROWS
require_var CFG_OUTPUT_TMP_DIR
require_var CFG_OUTPUT_OUT_DIR
require_var CFG_OUTPUT_QC_SUBDIR
require_var CFG_OUTPUT_TABLE_SUBDIR
require_var CFG_OUTPUT_FIGURE_PYTHON_SUBDIR
require_var CFG_OUTPUT_FIGURE_R_SUBDIR

if [[ ! -f "$CFG_RUNTIME_CONDA_SH" ]]; then
    log "conda.sh not found: $CFG_RUNTIME_CONDA_SH"
    exit 1
fi

source "$CFG_RUNTIME_CONDA_SH"
if [[ "${CONDA_DEFAULT_ENV:-}" != "$CFG_RUNTIME_CONDA_ENV" ]]; then
    conda activate "$CFG_RUNTIME_CONDA_ENV"
fi

PYTHON_BIN=$(command -v python)
RSCRIPT_BIN=$(command -v Rscript)

if [[ -z "$PYTHON_BIN" || -z "$RSCRIPT_BIN" ]]; then
    log "python or Rscript not found after activating $CFG_RUNTIME_CONDA_ENV"
    exit 1
fi

log "Checking Python dependencies"
"$PYTHON_BIN" -c "import pandas, numpy, matplotlib, seaborn"
log "Checking R dependencies"
"$RSCRIPT_BIN" -e "pkgs <- c('RIdeogram','ggplot2','readr','dplyr','tidyr','forcats','patchwork','scales','svglite'); miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; if (length(miss) > 0) { stop(paste('Missing R packages:', paste(miss, collapse=', '))) }"

TMP_DIR=$(resolve_path "$CFG_OUTPUT_TMP_DIR")
OUT_DIR=$(resolve_path "$CFG_OUTPUT_OUT_DIR")
QC_OUT_DIR="$OUT_DIR/$CFG_OUTPUT_QC_SUBDIR"
TABLE_OUT_DIR="$OUT_DIR/$CFG_OUTPUT_TABLE_SUBDIR"
PYTHON_FIG_OUT_DIR="$OUT_DIR/$CFG_OUTPUT_FIGURE_PYTHON_SUBDIR"
R_FIG_OUT_DIR="$OUT_DIR/$CFG_OUTPUT_FIGURE_R_SUBDIR"
THREADS="$CFG_RUNTIME_THREADS"
CHUNK_ROWS="$CFG_RUNTIME_CHUNK_ROWS"
KARYOTYPE_TXT="$CFG_REFERENCE_KARYOTYPE_TXT"

mkdir -p "$TMP_DIR" "$OUT_DIR" "$TMP_DIR/logs" "$QC_OUT_DIR" "$TABLE_OUT_DIR" "$PYTHON_FIG_OUT_DIR" "$R_FIG_OUT_DIR"

PRIMARY_MIN_SUPPORT="$CFG_ANALYSIS_DISTINCT_NUMT_PRIMARY_MIN_SAMPLE_SUPPORT"
PRIMARY_FREQ_CLASS_TSV="$TABLE_OUT_DIR/4-numt-frequency-class-summary.min-support-${PRIMARY_MIN_SUPPORT}.tsv"

run_with_log "01_filter_qc_samples" \
    "$PYTHON_BIN" "$PROJECT_ROOT/python/01_filter_qc_samples.py" \
    --meta-tsv "$CFG_INPUT_META_TSV" \
    --breakpoint-tsv-gz "$CFG_INPUT_BREAKPOINT_TSV_GZ" \
    --merge-bed-tsv-gz "$CFG_INPUT_MERGE_BED_TSV_GZ" \
    --meta-id-col "$CFG_META_ID_COL" \
    --meta-qc-col "$CFG_META_QC_COL" \
    --meta-qc-pass-value "$CFG_META_QC_PASS_VALUE" \
    --pass-sample-tsv "$TMP_DIR/1-pass-samples.tsv" \
    --breakpoint-out-gz "$TMP_DIR/2-confident_breakpoints.pass.tsv.gz" \
    --merge-bed-out-gz "$TMP_DIR/3-merge_bed.pass.tsv.gz" \
    --summary-out "$QC_OUT_DIR/1-qc-filter-summary.tsv" \
    --threads "$THREADS"

# 从 QC summary 读取通过质控的样本数作为 frequency denominator
QC_SUMMARY_TSV="$QC_OUT_DIR/1-qc-filter-summary.tsv"
FREQUENCY_DENOMINATOR=$(awk -F'\t' '$1=="meta_pass_samples"{print $2}' "$QC_SUMMARY_TSV")
if [[ -z "$FREQUENCY_DENOMINATOR" ]]; then
    log "Failed to read meta_pass_samples from $QC_SUMMARY_TSV"
    exit 1
fi
log "Frequency denominator (meta_pass_samples): $FREQUENCY_DENOMINATOR"

run_with_log "02_build_length_table" \
    "$PYTHON_BIN" "$PROJECT_ROOT/python/02_build_length_table.py" \
    --input-gz "$TMP_DIR/3-merge_bed.pass.tsv.gz" \
    --output-tsv "$TABLE_OUT_DIR/2-numt-length-by-region.tsv" \
    --summary-tsv "$TABLE_OUT_DIR/2-numt-length-summary.tsv" \
    --chunk-rows "$CHUNK_ROWS" \
    --temp-dir "$TMP_DIR/length_buckets" &
pid_length=$!

run_with_log "03_build_event_table" \
    "$PYTHON_BIN" "$PROJECT_ROOT/python/03_build_event_table.py" \
    --input-gz "$TMP_DIR/2-confident_breakpoints.pass.tsv.gz" \
    --events-out "$TABLE_OUT_DIR/3-numt-events.tsv" \
    --ideogram-out "$TABLE_OUT_DIR/3-numt-events-for-ideogram.tsv" \
    --chr-count-out "$TABLE_OUT_DIR/3-numt-chromosome-counts.tsv" \
    --chunk-rows "$CHUNK_ROWS" &
pid_events=$!

wait "$pid_events"

run_with_log "04_cluster_frequency" \
    "$PYTHON_BIN" "$PROJECT_ROOT/python/04_cluster_frequency.py" \
    --cluster-input-file "$CFG_INPUT_DISTINCT_NUMT_CLUSTER_INPUT_TSV_GZ" \
    --cluster-prefix "$TMP_DIR/4-mt_disc_breakpoint_input" \
    --cluster-out "$TABLE_OUT_DIR/4-numt-frequency-by-cluster.tsv" \
    --class-summary-out "$TABLE_OUT_DIR/4-numt-frequency-class-summary.tsv" \
    --support-summary-out "$TABLE_OUT_DIR/4-numt-support-summary.tsv" \
    --top-out "$TABLE_OUT_DIR/4-numt-top-recurrent-clusters.tsv" \
    --cluster-gap-bp "$CFG_ANALYSIS_FREQUENCY_CLUSTER_GAP_BP" \
    --denominator "$FREQUENCY_DENOMINATOR" \
    --min-supports "$CFG_ANALYSIS_DISTINCT_NUMT_MIN_SAMPLE_SUPPORTS" \
    --primary-min-support "$CFG_ANALYSIS_DISTINCT_NUMT_PRIMARY_MIN_SAMPLE_SUPPORT" \
    --threads "$THREADS"

wait "$pid_length"

run_with_log "05_plot_python_figures" \
    "$PYTHON_BIN" "$PROJECT_ROOT/python/05_plot_python_figures.py" \
    --length-tsv "$TABLE_OUT_DIR/2-numt-length-by-region.tsv" \
    --chr-count-tsv "$TABLE_OUT_DIR/3-numt-chromosome-counts.tsv" \
    --freq-cluster-tsv "$TABLE_OUT_DIR/4-numt-frequency-by-cluster.tsv" \
    --freq-class-tsv "$TABLE_OUT_DIR/4-numt-frequency-class-summary.tsv" \
    --support-summary-tsv "$TABLE_OUT_DIR/4-numt-support-summary.tsv" \
    --relative-frequency-tsv "$TABLE_OUT_DIR/5-numt-relative-frequency-percentage.tsv" \
    --out-prefix "$PYTHON_FIG_OUT_DIR/5-numt-statistical-summary" &
pid_python_plots=$!

run_with_log "01_plot_ideogram" \
    "$RSCRIPT_BIN" "$PROJECT_ROOT/src/01_plot_ideogram.R" \
    --input "$TABLE_OUT_DIR/3-numt-events.tsv" \
    --marker-input "$TABLE_OUT_DIR/4-numt-frequency-by-cluster.tsv" \
    --output-marker "$R_FIG_OUT_DIR/3-numt-ideogram-marker.svg" \
    --output-heatmap "$R_FIG_OUT_DIR/3-numt-ideogram-heatmap.svg" \
    --output-heatmap-1kb "$R_FIG_OUT_DIR/3-numt-ideogram-heatmap-1kb.svg" \
    --karyotype "$KARYOTYPE_TXT" \
    --bin-size "$CFG_ANALYSIS_IDEOGRAM_BIN_SIZE_BP" \
    --frequency-denominator "$FREQUENCY_DENOMINATOR" &
pid_ideogram=$!

run_with_log "02_plot_cluster_distribution" \
    "$RSCRIPT_BIN" "$PROJECT_ROOT/src/02_plot_cluster_distribution.R" \
    --input "$TABLE_OUT_DIR/4-numt-frequency-by-cluster.tsv" \
    --output-pdf "$R_FIG_OUT_DIR/4-numt-cluster-distribution-combined.pdf" \
    --output-png "$R_FIG_OUT_DIR/4-numt-cluster-distribution-combined.png" \
    --output-svg "$R_FIG_OUT_DIR/4-numt-cluster-distribution-combined.svg" \
    --frequency-denominator "$FREQUENCY_DENOMINATOR" &
pid_cluster_r=$!

run_with_log "06_mtdna_length_frequency" \
    "$PYTHON_BIN" "$PROJECT_ROOT/python/06_mtdna_length_frequency.py" \
    --merge-bed-gz "$TMP_DIR/3-merge_bed.pass.tsv.gz" \
    --cluster-detail-tsv "$TMP_DIR/4-mt_disc_breakpoint_input.min-support-1.allCluster.tsv" \
    --freq-cluster-tsv "$TABLE_OUT_DIR/4-numt-frequency-by-cluster.tsv" \
    --output-tsv "$TABLE_OUT_DIR/6-numt-mtdna-length-by-cluster.tsv" \
    --output-pdf "$PYTHON_FIG_OUT_DIR/6-numt-mtdna-length-frequency.pdf" &
pid_mtdna_scatter=$!

wait "$pid_python_plots"
wait "$pid_ideogram"
wait "$pid_cluster_r"
wait "$pid_mtdna_scatter"

convert_svg_to_png "$R_FIG_OUT_DIR/3-numt-ideogram-marker.svg" "$R_FIG_OUT_DIR/3-numt-ideogram-marker.png"
convert_svg_to_png "$R_FIG_OUT_DIR/3-numt-ideogram-heatmap.svg" "$R_FIG_OUT_DIR/3-numt-ideogram-heatmap.png"
convert_svg_to_png "$R_FIG_OUT_DIR/3-numt-ideogram-heatmap-1kb.svg" "$R_FIG_OUT_DIR/3-numt-ideogram-heatmap-1kb.png"

log "Pipeline finished successfully"
