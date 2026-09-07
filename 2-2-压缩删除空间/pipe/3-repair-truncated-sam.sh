#!/usr/bin/env bash
# pipe/3-repair-truncated-sam.sh
# 修复步骤 1 中因 SAM 截断而压缩失败的样本：保留可解析部分，转为规范 BAM。
# 样本数量少，生成单个 LSF 作业串行处理。
#
# 用法：
#   bash pipe/3-repair-truncated-sam.sh                    # dry-run，只报告可解析记录数
#   bash pipe/3-repair-truncated-sam.sh --dry-run false --submit
#   bash pipe/3-repair-truncated-sam.sh --collect
set -euo pipefail

SCRIPT_NAME="$(basename "$0")"
STEP_NAME="3-repair-truncated-sam"
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

CONFIG_PATH="${PROJECT_DIR}/conf/${STEP_NAME}.yaml"
OVERRIDE_DRY_RUN=""
DO_SUBMIT=""
DO_COLLECT="false"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)  CONFIG_PATH="$2";      shift 2 ;;
        --dry-run) OVERRIDE_DRY_RUN="$2"; shift 2 ;;
        --submit)  DO_SUBMIT="true";      shift 1 ;;
        --collect) DO_COLLECT="true";     shift 1 ;;
        *) echo "[ERROR] Unknown argument: $1" >&2; exit 2 ;;
    esac
done

# shellcheck source=../script/load_config.sh
source "${PROJECT_DIR}/script/load_config.sh" "$CONFIG_PATH"

INPUT_ROOT="$CFG_PATHS_INPUT_ROOT"
OUTPUT_RESULT="$CFG_PATHS_OUTPUT_RESULT"
TEMP_DIR="$CFG_PATHS_TEMP"
LOG_DIR="$CFG_PATHS_LOG"
SAMTOOLS_BIN="$CFG_TOOLS_SAMTOOLS_BIN"
SAMPLE_IDS="$CFG_SAMPLES_IDS"
SOURCE_SUFFIX="$CFG_COMPRESS_SOURCE_SUFFIX"
TARGET_SUFFIX="$CFG_COMPRESS_TARGET_SUFFIX"
THREADS="$CFG_COMPRESS_THREADS"
DRY_RUN="${OVERRIDE_DRY_RUN:-$CFG_RUNTIME_DRY_RUN}"
KEEP_SAM="$CFG_RUNTIME_KEEP_SAM"
QUEUE="$CFG_LSF_QUEUE"
CORES="$CFG_LSF_CORES"
WALLTIME="$CFG_LSF_WALLTIME"
JOB_PREFIX="$CFG_LSF_JOB_PREFIX"
SUBMIT="${DO_SUBMIT:-$CFG_LSF_SUBMIT}"

WORKER="${PROJECT_DIR}/script/repair_sam_sample.sh"
JOB_DIR="${TEMP_DIR}/jobs"
REPORT_DIR="${TEMP_DIR}/reports"
MERGED_REPORT="${OUTPUT_RESULT}/3-1-repair_sam_report.tsv"

mkdir -p "$OUTPUT_RESULT" "$TEMP_DIR" "$LOG_DIR" "$JOB_DIR" "$REPORT_DIR"

log() { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }

if [[ "$DO_COLLECT" == "true" ]]; then
    {
        printf 'sample_id\taction\tcategory\tfiles\tbytes\n'
        find "$REPORT_DIR" -maxdepth 1 -type f -name '*.tsv' -exec cat {} +
    } > "$MERGED_REPORT"
    log "报告：$MERGED_REPORT"
    column -t -s $'\t' "$MERGED_REPORT" 2>/dev/null || cat "$MERGED_REPORT"
    exit 0
fi

[[ -d "$INPUT_ROOT" ]] || { echo "[ERROR] input_root 不存在: $INPUT_ROOT" >&2; exit 1; }
[[ -x "$WORKER"     ]] || { echo "[ERROR] worker 不可执行: $WORKER" >&2; exit 1; }
[[ -x "$SAMTOOLS_BIN" ]] || { echo "[ERROR] samtools 不可执行: $SAMTOOLS_BIN" >&2; exit 1; }

log "配置    : $CONFIG_PATH"
log "待修复  : $SAMPLE_IDS"
log "dry_run : $DRY_RUN"
log "keep_sam: $KEEP_SAM"

job_name="${JOB_PREFIX}_0000"
job_file="${JOB_DIR}/${job_name}.lsf"

{
cat <<LSF
#!/bin/bash
#BSUB -J ${job_name}
#BSUB -q ${QUEUE}
#BSUB -n ${CORES}
#BSUB -R "span[hosts=1]"
#BSUB -W ${WALLTIME}
#BSUB -o ${LOG_DIR}/${job_name}.%J.out
#BSUB -e ${LOG_DIR}/${job_name}.%J.err

set -uo pipefail
export OMP_NUM_THREADS="\${LSB_DJOB_NUMPROC:-${CORES}}"

echo "[\$(date '+%F %T')] start ${job_name}"
LSF
for sid in $SAMPLE_IDS; do
cat <<LSF
"${WORKER}" \\
    --sample-dir "${INPUT_ROOT}/${sid}" \\
    --report "${REPORT_DIR}/${sid}.tsv" \\
    --dry-run '${DRY_RUN}' \\
    --keep-sam '${KEEP_SAM}' \\
    --samtools '${SAMTOOLS_BIN}' \\
    --source-suffix '${SOURCE_SUFFIX}' \\
    --target-suffix '${TARGET_SUFFIX}' \\
    --threads '${THREADS}'
echo "  ${sid} rc=\$?"
LSF
done
cat <<LSF
echo "[\$(date '+%F %T')] done ${job_name}"
exit 0
LSF
} > "$job_file"
chmod +x "$job_file"

log "已生成作业脚本：$job_file"

if [[ "$SUBMIT" != "true" ]]; then
    log "未提交（--submit 未指定）"
    exit 0
fi

command -v bsub >/dev/null 2>&1 || { echo "[ERROR] 未找到 bsub" >&2; exit 1; }
bsub < "$job_file" > "${job_file}.submit.log" 2>&1
log "提交: $(head -n1 "${job_file}.submit.log")"
log "作业结束后执行：bash pipe/${SCRIPT_NAME} --collect"
