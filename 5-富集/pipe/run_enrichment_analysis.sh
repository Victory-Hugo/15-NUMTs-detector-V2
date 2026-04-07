#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
CONFIG_PATH="${1:-$PROJECT_ROOT/conf/Config.yaml}"

source "$PROJECT_ROOT/script/load_config.sh" "$CONFIG_PATH"

require_var() {
    local var_name="$1"
    if [[ -z "${!var_name:-}" ]]; then
        echo "错误: 配置缺少必要字段: $var_name" >&2
        exit 1
    fi
}

resolve_path() {
    local base_dir="$1"
    local raw_path="$2"
    if [[ "$raw_path" = /* ]]; then
        printf '%s\n' "$raw_path"
    else
        printf '%s\n' "$base_dir/$raw_path"
    fi
}

resolve_tool() {
    local base_dir="$1"
    local raw_tool="$2"
    if [[ "$raw_tool" = /* ]]; then
        printf '%s\n' "$raw_tool"
    elif [[ "$raw_tool" == */* ]]; then
        printf '%s\n' "$base_dir/$raw_tool"
    else
        printf '%s\n' "$raw_tool"
    fi
}

require_var CFG_PROJECT_BASE_DIR
require_var CFG_PROJECT_BREAKPOINT_TSV_GZ
require_var CFG_PROJECT_CLUSTER_INPUT_TSV_GZ
require_var CFG_PROJECT_CLUSTER_DETAIL_TSV
require_var CFG_PROJECT_FREQUENCY_CLUSTER_TSV
require_var CFG_PATHS_PYTHON_DIR
require_var CFG_PATHS_REGIONS_DIR
require_var CFG_PATHS_OUTPUT_DIR
require_var CFG_PATHS_REPORT_DIR
require_var CFG_TOOLS_PYTHON_BIN
require_var CFG_TOOLS_BEDTOOLS_BIN
require_var CFG_RUNTIME_TARGET_BED_GLOB
require_var CFG_RUNTIME_FREQUENCY_CLASSES
require_var CFG_RUNTIME_FLANK_BP
require_var CFG_RUNTIME_SIMULATION_RUNS
require_var CFG_RUNTIME_RANDOM_SEED
require_var CFG_RUNTIME_JOBS
require_var CFG_RUNTIME_GENOME_LENGTH

BASE_DIR=$(resolve_path "$PROJECT_ROOT" "$CFG_PROJECT_BASE_DIR")
PYTHON_DIR=$(resolve_path "$BASE_DIR" "$CFG_PATHS_PYTHON_DIR")
REGIONS_DIR=$(resolve_path "$BASE_DIR" "$CFG_PATHS_REGIONS_DIR")
OUTPUT_DIR=$(resolve_path "$BASE_DIR" "$CFG_PATHS_OUTPUT_DIR")
REPORT_DIR=$(resolve_path "$BASE_DIR" "$CFG_PATHS_REPORT_DIR")
FREQUENCY_CLUSTER_TSV=$(resolve_path "$BASE_DIR" "$CFG_PROJECT_FREQUENCY_CLUSTER_TSV")
BREAKPOINT_TSV_GZ=$(resolve_path "$BASE_DIR" "$CFG_PROJECT_BREAKPOINT_TSV_GZ")
CLUSTER_INPUT_TSV_GZ=$(resolve_path "$BASE_DIR" "$CFG_PROJECT_CLUSTER_INPUT_TSV_GZ")
CLUSTER_DETAIL_TSV=$(resolve_path "$BASE_DIR" "$CFG_PROJECT_CLUSTER_DETAIL_TSV")
PYTHON_BIN=$(resolve_tool "$BASE_DIR" "$CFG_TOOLS_PYTHON_BIN")
BEDTOOLS_BIN=$(resolve_tool "$BASE_DIR" "$CFG_TOOLS_BEDTOOLS_BIN")
TARGET_BED_GLOB="$CFG_RUNTIME_TARGET_BED_GLOB"
FREQUENCY_CLASSES="$CFG_RUNTIME_FREQUENCY_CLASSES"
FLANK_BP="$CFG_RUNTIME_FLANK_BP"
SIMULATION_RUNS="$CFG_RUNTIME_SIMULATION_RUNS"
RANDOM_SEED="$CFG_RUNTIME_RANDOM_SEED"
JOBS="$CFG_RUNTIME_JOBS"
GENOME_LENGTH="$CFG_RUNTIME_GENOME_LENGTH"

for input_file in "$BREAKPOINT_TSV_GZ" "$CLUSTER_INPUT_TSV_GZ" "$CLUSTER_DETAIL_TSV" "$FREQUENCY_CLUSTER_TSV"; do
    if [[ ! -f "$input_file" ]]; then
        echo "错误: 输入文件不存在: $input_file" >&2
        exit 1
    fi
done

if [[ ! -d "$REGIONS_DIR" ]]; then
    echo "错误: 目标区域BED目录不存在: $REGIONS_DIR" >&2
    exit 1
fi

if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
    echo "错误: 找不到 Python 解释器: $PYTHON_BIN" >&2
    exit 1
fi

if ! command -v "$BEDTOOLS_BIN" >/dev/null 2>&1; then
    echo "错误: 找不到 bedtools: $BEDTOOLS_BIN" >&2
    exit 1
fi

PREPARED_DIR="$OUTPUT_DIR/1-prepared-breakpoint-flanks"
ENRICHMENT_DIR="$OUTPUT_DIR/2-enrichment"
SUMMARY_DIR="$OUTPUT_DIR/3-summary"
FIGURE_DIR="$OUTPUT_DIR/4-figures"
LOG_DIR="$OUTPUT_DIR/logs"
TASK_MANIFEST="$PREPARED_DIR/enrichment_task_manifest.tsv"
COMBINED_SUMMARY="$SUMMARY_DIR/enrichment_summary.tsv"
REPORT_MD="$REPORT_DIR/enrichment_report.md"

rm -rf "$PREPARED_DIR" "$ENRICHMENT_DIR" "$SUMMARY_DIR" "$FIGURE_DIR" "$LOG_DIR" "$REPORT_DIR"
mkdir -p "$PREPARED_DIR" "$ENRICHMENT_DIR" "$SUMMARY_DIR" "$FIGURE_DIR" "$LOG_DIR" "$REPORT_DIR"

echo "=== Step 1: prepare nuclear breakpoint flank BED files and task manifest ==="
"$PYTHON_BIN" "$PYTHON_DIR/prepare_breakpoint_flanks.py" \
    --breakpoint-tsv-gz "$BREAKPOINT_TSV_GZ" \
    --cluster-input-tsv-gz "$CLUSTER_INPUT_TSV_GZ" \
    --cluster-detail-tsv "$CLUSTER_DETAIL_TSV" \
    --frequency-cluster-tsv "$FREQUENCY_CLUSTER_TSV" \
    --regions-dir "$REGIONS_DIR" \
    --target-bed-glob "$TARGET_BED_GLOB" \
    --frequency-classes "$FREQUENCY_CLASSES" \
    --flank-bp "$FLANK_BP" \
    --prepared-dir "$PREPARED_DIR" \
    --enrichment-dir "$ENRICHMENT_DIR" \
    --task-manifest "$TASK_MANIFEST" \
    2>&1 | tee "$LOG_DIR/1-prepare_breakpoint_flanks.log"

TASK_COUNT=$(tail -n +2 "$TASK_MANIFEST" | wc -l)
if [[ "$TASK_COUNT" -le 0 ]]; then
    echo "错误: 没有生成可运行的富集任务。" >&2
    exit 1
fi

echo "=== Step 2: run full enrichment tasks with shell job control ==="
echo "任务数: $TASK_COUNT"
echo "并行任务数: $JOBS"
echo "每任务模拟次数: $SIMULATION_RUNS"

run_one_task() {
    local task_id="$1"
    local region_id="$2"
    local region_name="$3"
    local target_bed="$4"
    local frequency_class="$5"
    local breakpoint_bed="$6"
    local task_output_dir="$7"
    local task_log="$task_output_dir/stdout.log"

    mkdir -p "$task_output_dir"
    {
        echo "[$(date '+%F %T')] Start $task_id"
        "$PYTHON_BIN" "$PYTHON_DIR/run_one_breakpoint_enrichment.py" \
            --task-id "$task_id" \
            --region-id "$region_id" \
            --region-name "$region_name" \
            --target-bed "$target_bed" \
            --frequency-class "$frequency_class" \
            --breakpoint-bed "$breakpoint_bed" \
            --output-dir "$task_output_dir" \
            --bedtools-bin "$BEDTOOLS_BIN" \
            --simulation-runs "$SIMULATION_RUNS" \
            --random-seed "$RANDOM_SEED" \
            --genome-length "$GENOME_LENGTH" \
            --flank-bp "$FLANK_BP"
        echo "[$(date '+%F %T')] Done $task_id"
    } >"$task_log" 2>&1
}

while IFS=$'\t' read -r task_id region_id region_name target_bed frequency_class breakpoint_bed task_output_dir; do
    run_one_task "$task_id" "$region_id" "$region_name" "$target_bed" "$frequency_class" "$breakpoint_bed" "$task_output_dir" &
    while [[ "$(jobs -rp | wc -l)" -ge "$JOBS" ]]; do
        wait -n
    done
done < <(tail -n +2 "$TASK_MANIFEST")

wait

echo "=== Step 3: collect enrichment summaries ==="
"$PYTHON_BIN" "$PYTHON_DIR/collect_enrichment_results.py" \
    --task-manifest "$TASK_MANIFEST" \
    --output-tsv "$COMBINED_SUMMARY" \
    2>&1 | tee "$LOG_DIR/3-collect_enrichment_results.log"

echo "=== Step 4: plot enrichment summaries ==="
"$PYTHON_BIN" "$PYTHON_DIR/plot_enrichment_results.py" \
    --summary-tsv "$COMBINED_SUMMARY" \
    --output-dir "$FIGURE_DIR" \
    2>&1 | tee "$LOG_DIR/4-plot_enrichment_results.log"

echo "=== Step 5: write Markdown report ==="
"$PYTHON_BIN" "$PYTHON_DIR/write_enrichment_report.py" \
    --summary-tsv "$COMBINED_SUMMARY" \
    --prepared-summary "$PREPARED_DIR/prepared_breakpoint_flanks_summary.tsv" \
    --output-md "$REPORT_MD" \
    2>&1 | tee "$LOG_DIR/5-write_enrichment_report.log"

echo "=== Pipeline finished successfully ==="
echo "Combined summary: $COMBINED_SUMMARY"
echo "Figures: $FIGURE_DIR"
echo "Report: $REPORT_MD"
