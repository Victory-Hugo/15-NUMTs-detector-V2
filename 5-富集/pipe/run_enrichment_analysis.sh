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
require_var CFG_PATHS_SRC_DIR
require_var CFG_PATHS_REGIONS_DIR
require_var CFG_PATHS_OUTPUT_DIR
require_var CFG_PATHS_REPORT_DIR
require_var CFG_TOOLS_PYTHON_BIN
require_var CFG_TOOLS_RSCRIPT_BIN
require_var CFG_TOOLS_BEDTOOLS_BIN
require_var CFG_RUNTIME_TARGET_BED_GLOB
require_var CFG_RUNTIME_FREQUENCY_CLASSES
require_var CFG_RUNTIME_FLANK_BP
require_var CFG_RUNTIME_SIMULATION_RUNS
require_var CFG_RUNTIME_RANDOM_SEED
require_var CFG_RUNTIME_JOBS
require_var CFG_RUNTIME_GENOME_LENGTH
require_var CFG_MTDNA_MT_GENOME_LENGTH
require_var CFG_MTDNA_MT_REGIONS_DIR

BASE_DIR=$(resolve_path "$PROJECT_ROOT" "$CFG_PROJECT_BASE_DIR")
PYTHON_DIR=$(resolve_path "$BASE_DIR" "$CFG_PATHS_PYTHON_DIR")
SRC_DIR=$(resolve_path "$BASE_DIR" "$CFG_PATHS_SRC_DIR")
REGIONS_DIR=$(resolve_path "$BASE_DIR" "$CFG_PATHS_REGIONS_DIR")
OUTPUT_DIR=$(resolve_path "$BASE_DIR" "$CFG_PATHS_OUTPUT_DIR")
REPORT_DIR=$(resolve_path "$BASE_DIR" "$CFG_PATHS_REPORT_DIR")
FREQUENCY_CLUSTER_TSV=$(resolve_path "$BASE_DIR" "$CFG_PROJECT_FREQUENCY_CLUSTER_TSV")
BREAKPOINT_TSV_GZ=$(resolve_path "$BASE_DIR" "$CFG_PROJECT_BREAKPOINT_TSV_GZ")
CLUSTER_INPUT_TSV_GZ=$(resolve_path "$BASE_DIR" "$CFG_PROJECT_CLUSTER_INPUT_TSV_GZ")
CLUSTER_DETAIL_TSV=$(resolve_path "$BASE_DIR" "$CFG_PROJECT_CLUSTER_DETAIL_TSV")
PYTHON_BIN=$(resolve_tool "$BASE_DIR" "$CFG_TOOLS_PYTHON_BIN")
RSCRIPT_BIN=$(resolve_tool "$BASE_DIR" "$CFG_TOOLS_RSCRIPT_BIN")
BEDTOOLS_BIN=$(resolve_tool "$BASE_DIR" "$CFG_TOOLS_BEDTOOLS_BIN")
TARGET_BED_GLOB="$CFG_RUNTIME_TARGET_BED_GLOB"
FREQUENCY_CLASSES="$CFG_RUNTIME_FREQUENCY_CLASSES"
FLANK_BP="$CFG_RUNTIME_FLANK_BP"
SIMULATION_RUNS="$CFG_RUNTIME_SIMULATION_RUNS"
RANDOM_SEED="$CFG_RUNTIME_RANDOM_SEED"
JOBS="$CFG_RUNTIME_JOBS"
GENOME_LENGTH="$CFG_RUNTIME_GENOME_LENGTH"
MT_GENOME_LENGTH="$CFG_MTDNA_MT_GENOME_LENGTH"
MT_REGIONS_DIR=$(resolve_path "$BASE_DIR" "$CFG_MTDNA_MT_REGIONS_DIR")
MT_GENES_TXT="$BASE_DIR/data/mtDNA-region/MitoGenes.txt"
MT_FINE_REGIONS_DIR="$MT_REGIONS_DIR/fine_grained"

MT_PREPARED_DIR="$OUTPUT_DIR/1-mt-breakpoint-flanks"
MT_ENRICHMENT_DIR="$OUTPUT_DIR/2-mt-enrichment"
MT_SUMMARY_DIR="$OUTPUT_DIR/3-mt-summary"
MT_FIGURE_DIR="$OUTPUT_DIR/4-mt-figures"
MT_TASK_MANIFEST="$MT_PREPARED_DIR/mt_enrichment_task_manifest.tsv"
MT_COMBINED_SUMMARY="$MT_SUMMARY_DIR/mt_enrichment_summary.tsv"
MT_FINE_PREPARED_DIR="$OUTPUT_DIR/1-mt-breakpoint-flanks-fine"
MT_FINE_ENRICHMENT_DIR="$OUTPUT_DIR/2-mt-enrichment-fine"
MT_FINE_SUMMARY_DIR="$OUTPUT_DIR/3-mt-summary-fine"
MT_FINE_FIGURE_DIR="$OUTPUT_DIR/4-mt-figures-fine"
MT_FINE_TASK_MANIFEST="$MT_FINE_PREPARED_DIR/mt_enrichment_task_manifest.tsv"
MT_FINE_COMBINED_SUMMARY="$MT_FINE_SUMMARY_DIR/mt_enrichment_summary.fine.tsv"
MT_FINE_REPORT_MD="$REPORT_DIR/mt_enrichment_report.fine.md"
CENTROMERE_BED="$REGIONS_DIR/15-Centromeres.bed"
TELOMERE_BED="$BASE_DIR/data/landmark-beds/Telomeres.bed"

for input_file in "$BREAKPOINT_TSV_GZ" "$CLUSTER_INPUT_TSV_GZ" "$CLUSTER_DETAIL_TSV" "$FREQUENCY_CLUSTER_TSV"; do
    if [[ ! -f "$input_file" ]]; then
        echo "错误: 输入文件不存在: $input_file" >&2
        exit 1
    fi
done

if [[ ! -f "$CENTROMERE_BED" ]]; then
    echo "错误: 着丝粒BED文件不存在: $CENTROMERE_BED" >&2
    exit 1
fi

if [[ ! -f "$TELOMERE_BED" ]]; then
    echo "错误: 端粒BED文件不存在: $TELOMERE_BED" >&2
    exit 1
fi

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

if ! command -v "$RSCRIPT_BIN" >/dev/null 2>&1; then
    echo "错误: 找不到 Rscript: $RSCRIPT_BIN" >&2
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

rm -rf "$PREPARED_DIR" "$ENRICHMENT_DIR" "$SUMMARY_DIR" "$FIGURE_DIR" "$LOG_DIR" "$REPORT_DIR" \
       "$MT_PREPARED_DIR" "$MT_ENRICHMENT_DIR" "$MT_SUMMARY_DIR" "$MT_FIGURE_DIR" \
       "$MT_FINE_PREPARED_DIR" "$MT_FINE_ENRICHMENT_DIR" "$MT_FINE_SUMMARY_DIR" "$MT_FINE_FIGURE_DIR" \
       "$MT_FINE_REGIONS_DIR"
mkdir -p "$PREPARED_DIR" "$ENRICHMENT_DIR" "$SUMMARY_DIR" "$FIGURE_DIR" "$LOG_DIR" "$REPORT_DIR" \
         "$MT_PREPARED_DIR" "$MT_ENRICHMENT_DIR" "$MT_SUMMARY_DIR" "$MT_FIGURE_DIR" \
         "$MT_FINE_PREPARED_DIR" "$MT_FINE_ENRICHMENT_DIR" "$MT_FINE_SUMMARY_DIR" "$MT_FINE_FIGURE_DIR"

echo "=== Step 0: generate mtDNA target category BED files ==="
"$PYTHON_BIN" "$PYTHON_DIR/prepare_mt_target_beds.py" \
    --mito-genes-txt "$MT_GENES_TXT" \
    --output-dir "$MT_REGIONS_DIR" \
    --fine-output-dir "$MT_FINE_REGIONS_DIR" \
    2>&1 | tee "$LOG_DIR/0-prepare_mt_target_beds.log"

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

echo "=== Step 3b: collect nuclear enrichment summaries (unstratified, all class only) ==="
"$PYTHON_BIN" "$PYTHON_DIR/collect_enrichment_results.py" \
    --task-manifest "$TASK_MANIFEST" \
    --output-tsv "$SUMMARY_DIR/enrichment_summary.all_only.tsv" \
    --filter-class "all" \
    2>&1 | tee "$LOG_DIR/3b-collect_nuclear_all_only.log"

echo "=== Step 4: plot enrichment summaries ==="
"$PYTHON_BIN" "$PYTHON_DIR/plot_enrichment_results.py" \
    --summary-tsv "$COMBINED_SUMMARY" \
    --output-dir "$FIGURE_DIR" \
    2>&1 | tee "$LOG_DIR/4-plot_enrichment_results.log"

"$RSCRIPT_BIN" "$SRC_DIR/03_plot_enrichment_pheatmap.R" \
    --summary-tsv "$COMBINED_SUMMARY" \
    --output-dir "$FIGURE_DIR" \
    --frequency-classes "all,ultra-rare,rare,low-frequency,common" \
    2>&1 | tee "$LOG_DIR/4-plot_enrichment_pheatmap.log"

echo "=== Step 4b: plot nuclear unstratified (all class only) figures ==="
"$PYTHON_BIN" "$PYTHON_DIR/plot_enrichment_results.py" \
    --summary-tsv "$SUMMARY_DIR/enrichment_summary.all_only.tsv" \
    --output-dir "$FIGURE_DIR/all_only" \
    2>&1 | tee "$LOG_DIR/4b-plot_nuclear_all_only.log"

"$RSCRIPT_BIN" "$SRC_DIR/03_plot_enrichment_pheatmap.R" \
    --summary-tsv "$SUMMARY_DIR/enrichment_summary.all_only.tsv" \
    --output-dir "$FIGURE_DIR/all_only" \
    --frequency-classes "all" \
    2>&1 | tee "$LOG_DIR/4b-plot_nuclear_all_only_pheatmap.log"

echo "=== Step 5: write Markdown report ==="
"$PYTHON_BIN" "$PYTHON_DIR/write_enrichment_report.py" \
    --summary-tsv "$COMBINED_SUMMARY" \
    --prepared-summary "$PREPARED_DIR/prepared_breakpoint_flanks_summary.tsv" \
    --output-md "$REPORT_MD" \
    2>&1 | tee "$LOG_DIR/5-write_enrichment_report.log"

echo "=== Step 5b: plot nuclear breakpoint distances to centromeres and telomeres ==="
"$RSCRIPT_BIN" "$SRC_DIR/04_plot_centromere_telomere_distance.R" \
    --breakpoint-bed "$PREPARED_DIR/nuclear_breakpoint_flanks.all_records.bed" \
    --centromere-bed "$CENTROMERE_BED" \
    --telomere-bed "$TELOMERE_BED" \
    --output-dir "$FIGURE_DIR" \
    --frequency-classes "$FREQUENCY_CLASSES" \
    2>&1 | tee "$LOG_DIR/5b-plot_centromere_telomere_distance.log"

echo "=== Step 6: prepare mtDNA breakpoint flank BED files ==="
"$PYTHON_BIN" "$PYTHON_DIR/prepare_mt_breakpoint_flanks.py" \
    --breakpoint-tsv-gz "$BREAKPOINT_TSV_GZ" \
    --cluster-input-tsv-gz "$CLUSTER_INPUT_TSV_GZ" \
    --cluster-detail-tsv "$CLUSTER_DETAIL_TSV" \
    --frequency-cluster-tsv "$FREQUENCY_CLUSTER_TSV" \
    --mt-regions-dir "$MT_REGIONS_DIR" \
    --target-bed-glob "$TARGET_BED_GLOB" \
    --frequency-classes "$FREQUENCY_CLASSES" \
    --flank-bp "$FLANK_BP" \
    --prepared-dir "$MT_PREPARED_DIR" \
    --enrichment-dir "$MT_ENRICHMENT_DIR" \
    --task-manifest "$MT_TASK_MANIFEST" \
    2>&1 | tee "$LOG_DIR/6-prepare_mt_breakpoint_flanks.log"

MT_TASK_COUNT=$(tail -n +2 "$MT_TASK_MANIFEST" | wc -l)
if [[ "$MT_TASK_COUNT" -le 0 ]]; then
    echo "错误: 没有生成可运行的 mtDNA 富集任务。" >&2
    exit 1
fi

echo "=== Step 7: run mtDNA enrichment tasks with shell job control ==="
echo "mtDNA任务数: $MT_TASK_COUNT"

run_one_mt_task() {
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
        "$PYTHON_BIN" "$PYTHON_DIR/run_one_mt_enrichment.py" \
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
            --genome-length "$MT_GENOME_LENGTH" \
            --flank-bp "$FLANK_BP"
        echo "[$(date '+%F %T')] Done $task_id"
    } >"$task_log" 2>&1
}

while IFS=$'\t' read -r task_id region_id region_name target_bed frequency_class breakpoint_bed task_output_dir; do
    run_one_mt_task "$task_id" "$region_id" "$region_name" "$target_bed" "$frequency_class" "$breakpoint_bed" "$task_output_dir" &
    while [[ "$(jobs -rp | wc -l)" -ge "$JOBS" ]]; do
        wait -n
    done
done < <(tail -n +2 "$MT_TASK_MANIFEST")

wait

echo "=== Step 8: collect mtDNA enrichment summaries ==="
"$PYTHON_BIN" "$PYTHON_DIR/collect_enrichment_results.py" \
    --task-manifest "$MT_TASK_MANIFEST" \
    --output-tsv "$MT_COMBINED_SUMMARY" \
    2>&1 | tee "$LOG_DIR/8-collect_mt_enrichment_results.log"

"$PYTHON_BIN" "$PYTHON_DIR/collect_enrichment_results.py" \
    --task-manifest "$MT_TASK_MANIFEST" \
    --output-tsv "$MT_SUMMARY_DIR/mt_enrichment_summary.all_only.tsv" \
    --filter-class "all" \
    2>&1 | tee "$LOG_DIR/8b-collect_mt_all_only.log"

echo "=== Step 9: plot mtDNA enrichment summaries ==="
"$PYTHON_BIN" "$PYTHON_DIR/plot_enrichment_results.py" \
    --summary-tsv "$MT_COMBINED_SUMMARY" \
    --output-dir "$MT_FIGURE_DIR" \
    2>&1 | tee "$LOG_DIR/9-plot_mt_enrichment_results.log"

"$PYTHON_BIN" "$PYTHON_DIR/plot_enrichment_results.py" \
    --summary-tsv "$MT_SUMMARY_DIR/mt_enrichment_summary.all_only.tsv" \
    --output-dir "$MT_FIGURE_DIR/all_only" \
    2>&1 | tee "$LOG_DIR/9b-plot_mt_all_only.log"

"$RSCRIPT_BIN" "$SRC_DIR/03_plot_enrichment_pheatmap.R" \
    --summary-tsv "$MT_SUMMARY_DIR/mt_enrichment_summary.all_only.tsv" \
    --output-dir "$MT_FIGURE_DIR/all_only" \
    --frequency-classes "all" \
    2>&1 | tee "$LOG_DIR/9b-plot_mt_all_only_pheatmap.log"

"$RSCRIPT_BIN" "$SRC_DIR/03_plot_enrichment_pheatmap.R" \
    --summary-tsv "$MT_COMBINED_SUMMARY" \
    --output-dir "$MT_FIGURE_DIR" \
    --frequency-classes "all,ultra-rare,rare,low-frequency,common" \
    2>&1 | tee "$LOG_DIR/9c-plot_mt_pheatmap.log"

echo "=== Step 10: write mtDNA Markdown report ==="
"$PYTHON_BIN" "$PYTHON_DIR/write_enrichment_report.py" \
    --summary-tsv "$MT_COMBINED_SUMMARY" \
    --prepared-summary "$MT_PREPARED_DIR/prepared_breakpoint_flanks_summary.tsv" \
    --output-md "$REPORT_DIR/mt_enrichment_report.md" \
    2>&1 | tee "$LOG_DIR/10-write_mt_enrichment_report.log"

echo "=== Step 11: prepare fine-grained mtDNA breakpoint flank BED files ==="
"$PYTHON_BIN" "$PYTHON_DIR/prepare_mt_breakpoint_flanks.py" \
    --breakpoint-tsv-gz "$BREAKPOINT_TSV_GZ" \
    --cluster-input-tsv-gz "$CLUSTER_INPUT_TSV_GZ" \
    --cluster-detail-tsv "$CLUSTER_DETAIL_TSV" \
    --frequency-cluster-tsv "$FREQUENCY_CLUSTER_TSV" \
    --mt-regions-dir "$MT_FINE_REGIONS_DIR" \
    --target-bed-glob "$TARGET_BED_GLOB" \
    --frequency-classes "$FREQUENCY_CLASSES" \
    --flank-bp "$FLANK_BP" \
    --prepared-dir "$MT_FINE_PREPARED_DIR" \
    --enrichment-dir "$MT_FINE_ENRICHMENT_DIR" \
    --task-manifest "$MT_FINE_TASK_MANIFEST" \
    2>&1 | tee "$LOG_DIR/11-prepare_mt_breakpoint_flanks_fine.log"

MT_FINE_TASK_COUNT=$(tail -n +2 "$MT_FINE_TASK_MANIFEST" | wc -l)
if [[ "$MT_FINE_TASK_COUNT" -le 0 ]]; then
    echo "错误: 没有生成可运行的细粒度 mtDNA 富集任务。" >&2
    exit 1
fi

echo "=== Step 12: run fine-grained mtDNA enrichment tasks with shell job control ==="
echo "细粒度mtDNA任务数: $MT_FINE_TASK_COUNT"

while IFS=$'\t' read -r task_id region_id region_name target_bed frequency_class breakpoint_bed task_output_dir; do
    run_one_mt_task "$task_id" "$region_id" "$region_name" "$target_bed" "$frequency_class" "$breakpoint_bed" "$task_output_dir" &
    while [[ "$(jobs -rp | wc -l)" -ge "$JOBS" ]]; do
        wait -n
    done
done < <(tail -n +2 "$MT_FINE_TASK_MANIFEST")

wait

echo "=== Step 13: collect fine-grained mtDNA enrichment summaries ==="
"$PYTHON_BIN" "$PYTHON_DIR/collect_enrichment_results.py" \
    --task-manifest "$MT_FINE_TASK_MANIFEST" \
    --output-tsv "$MT_FINE_COMBINED_SUMMARY" \
    2>&1 | tee "$LOG_DIR/13-collect_mt_enrichment_results_fine.log"

"$PYTHON_BIN" "$PYTHON_DIR/collect_enrichment_results.py" \
    --task-manifest "$MT_FINE_TASK_MANIFEST" \
    --output-tsv "$MT_FINE_SUMMARY_DIR/mt_enrichment_summary.fine.all_only.tsv" \
    --filter-class "all" \
    2>&1 | tee "$LOG_DIR/13b-collect_mt_all_only_fine.log"

echo "=== Step 14: plot fine-grained mtDNA enrichment summaries ==="
"$PYTHON_BIN" "$PYTHON_DIR/plot_enrichment_results.py" \
    --summary-tsv "$MT_FINE_COMBINED_SUMMARY" \
    --output-dir "$MT_FINE_FIGURE_DIR" \
    2>&1 | tee "$LOG_DIR/14-plot_mt_enrichment_results_fine.log"

"$PYTHON_BIN" "$PYTHON_DIR/plot_enrichment_results.py" \
    --summary-tsv "$MT_FINE_SUMMARY_DIR/mt_enrichment_summary.fine.all_only.tsv" \
    --output-dir "$MT_FINE_FIGURE_DIR/all_only" \
    2>&1 | tee "$LOG_DIR/14b-plot_mt_all_only_fine.log"

"$RSCRIPT_BIN" "$SRC_DIR/03_plot_enrichment_pheatmap.R" \
    --summary-tsv "$MT_FINE_SUMMARY_DIR/mt_enrichment_summary.fine.all_only.tsv" \
    --output-dir "$MT_FINE_FIGURE_DIR/all_only" \
    --frequency-classes "all" \
    2>&1 | tee "$LOG_DIR/14b-plot_mt_all_only_pheatmap_fine.log"

"$RSCRIPT_BIN" "$SRC_DIR/03_plot_enrichment_pheatmap.R" \
    --summary-tsv "$MT_FINE_COMBINED_SUMMARY" \
    --output-dir "$MT_FINE_FIGURE_DIR" \
    --frequency-classes "all,ultra-rare,rare,low-frequency,common" \
    2>&1 | tee "$LOG_DIR/14c-plot_mt_pheatmap_fine.log"

echo "=== Step 15: write fine-grained mtDNA Markdown report ==="
"$PYTHON_BIN" "$PYTHON_DIR/write_enrichment_report.py" \
    --summary-tsv "$MT_FINE_COMBINED_SUMMARY" \
    --prepared-summary "$MT_FINE_PREPARED_DIR/prepared_breakpoint_flanks_summary.tsv" \
    --output-md "$MT_FINE_REPORT_MD" \
    2>&1 | tee "$LOG_DIR/15-write_mt_enrichment_report_fine.log"

echo "=== Pipeline finished successfully ==="
echo "Combined summary: $COMBINED_SUMMARY"
echo "Figures: $FIGURE_DIR"
echo "Report: $REPORT_MD"
echo "Fine-grained mtDNA summary: $MT_FINE_COMBINED_SUMMARY"
echo "Fine-grained mtDNA figures: $MT_FINE_FIGURE_DIR"
echo "Fine-grained mtDNA report: $MT_FINE_REPORT_MD"
