#!/usr/bin/env bash
# pipe/2-clean-vardetect-output.sh
# 清理 1-5-Vardetect 输出目录：删除可由 *aln.fa + *_numts.bed 再生的逐碱基 TSV
# 及可秒级重建的拼接 FASTA。将样本切块后生成 LSF 作业脚本，按需提交（不在前台执行清理）。
#
# 用法：
#   bash pipe/2-clean-vardetect-output.sh                    # 生成作业脚本（默认 dry-run，不提交）
#   bash pipe/2-clean-vardetect-output.sh --submit
#   bash pipe/2-clean-vardetect-output.sh --dry-run false --limit 50 --submit
#   bash pipe/2-clean-vardetect-output.sh --collect          # 合并样本报告
set -euo pipefail

SCRIPT_NAME="$(basename "$0")"
STEP_NAME="2-clean-vardetect-output"
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

CONFIG_PATH="${PROJECT_DIR}/conf/${STEP_NAME}.yaml"
OVERRIDE_DRY_RUN=""
OVERRIDE_LIMIT=""
DO_SUBMIT=""
DO_COLLECT="false"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)  CONFIG_PATH="$2";      shift 2 ;;
        --dry-run) OVERRIDE_DRY_RUN="$2"; shift 2 ;;
        --limit)   OVERRIDE_LIMIT="$2";   shift 2 ;;
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
PARALLEL_BIN="$CFG_TOOLS_PARALLEL_BIN"
DELETE_PATTERNS="$CFG_DELETE_PATTERNS"
KEEP_PATTERNS="$CFG_KEEP_PATTERNS"
DRY_RUN="${OVERRIDE_DRY_RUN:-$CFG_RUNTIME_DRY_RUN}"
LIMIT="${OVERRIDE_LIMIT:-$CFG_RUNTIME_LIMIT}"
JOBS="$CFG_RUNTIME_JOBS"
SAMPLES_PER_JOB="$CFG_RUNTIME_SAMPLES_PER_JOB"
QUEUE="$CFG_LSF_QUEUE"
CORES="$CFG_LSF_CORES"
WALLTIME="$CFG_LSF_WALLTIME"
JOB_PREFIX="$CFG_LSF_JOB_PREFIX"
SUBMIT="${DO_SUBMIT:-$CFG_LSF_SUBMIT}"

WORKER="${PROJECT_DIR}/script/clean_vardetect_sample.sh"
CHUNK_DIR="${TEMP_DIR}/chunks"
JOB_DIR="${TEMP_DIR}/jobs"
REPORT_DIR="${TEMP_DIR}/reports"
SAMPLE_LIST="${TEMP_DIR}/samples.txt"
MERGED_REPORT="${OUTPUT_RESULT}/2-1-clean_vardetect_report.tsv"
SUMMARY_REPORT="${OUTPUT_RESULT}/2-2-clean_vardetect_summary.tsv"

mkdir -p "$OUTPUT_RESULT" "$TEMP_DIR" "$LOG_DIR" "$CHUNK_DIR" "$JOB_DIR" "$REPORT_DIR"

# 模式清单落盘（每行一个）。作业脚本以文件路径传给 worker，
# 避免 GNU parallel 重新解析命令行时把带空格的模式串拆散。
DELETE_PATTERNS_FILE="${TEMP_DIR}/delete_patterns.txt"
KEEP_PATTERNS_FILE="${TEMP_DIR}/keep_patterns.txt"
printf '%s\n' $DELETE_PATTERNS > "$DELETE_PATTERNS_FILE"
printf '%s\n' $KEEP_PATTERNS   > "$KEEP_PATTERNS_FILE"

log() { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }

if [[ "$DO_COLLECT" == "true" ]]; then
    log "合并样本报告：$REPORT_DIR"
    {
        printf 'sample_id\taction\tcategory\tfiles\tbytes\n'
        find "$REPORT_DIR" -maxdepth 1 -type f -name '*.tsv' -exec cat {} +
    } > "$MERGED_REPORT"

    # 先落盘再拆表头排序：head/tail 直接读文件，避免管道缓冲吞掉 tail 的输入。
    awk -F'\t' 'NR>1 && $2!="status" {f[$2"\t"$3]+=$4; b[$2"\t"$3]+=$5}
        END {
            printf "action\tcategory\tfiles\tbytes\tgib\n"
            for (k in f) printf "%s\t%d\t%d\t%.3f\n", k, f[k], b[k], b[k]/1073741824
        }' "$MERGED_REPORT" > "${SUMMARY_REPORT}.raw"
    {
        head -n 1 "${SUMMARY_REPORT}.raw"
        tail -n +2 "${SUMMARY_REPORT}.raw" | LC_ALL=C sort -k1,1 -k4,4gr
    } > "$SUMMARY_REPORT"
    rm -f "${SUMMARY_REPORT}.raw"

    log "长表：$MERGED_REPORT"
    log "汇总：$SUMMARY_REPORT"
    column -t -s $'\t' "$SUMMARY_REPORT" 2>/dev/null || cat "$SUMMARY_REPORT"
    exit 0
fi

[[ -d "$INPUT_ROOT" ]] || { echo "[ERROR] input_root 不存在: $INPUT_ROOT" >&2; exit 1; }
[[ -x "$WORKER"     ]] || { echo "[ERROR] worker 不可执行: $WORKER" >&2; exit 1; }
[[ -x "$PARALLEL_BIN" ]] || { echo "[ERROR] parallel 不可执行: $PARALLEL_BIN" >&2; exit 1; }

log "配置        : $CONFIG_PATH"
log "输入根目录  : $INPUT_ROOT"
log "dry_run     : $DRY_RUN"
log "limit       : $LIMIT"

find "$INPUT_ROOT" -mindepth 1 -maxdepth 1 -type d ! -name 'logs' ! -name 'log' \
     -printf '%p\n' 2>/dev/null | LC_ALL=C sort > "$SAMPLE_LIST"

TOTAL_SAMPLES=$(wc -l < "$SAMPLE_LIST")
if [[ "$LIMIT" != "0" ]]; then
    head -n "$LIMIT" "$SAMPLE_LIST" > "${SAMPLE_LIST}.use"
else
    cp "$SAMPLE_LIST" "${SAMPLE_LIST}.use"
fi
USE_SAMPLES=$(wc -l < "${SAMPLE_LIST}.use")
log "样本总数    : $TOTAL_SAMPLES（本次处理 $USE_SAMPLES）"
[[ "$USE_SAMPLES" -gt 0 ]] || { echo "[ERROR] 没有可处理的样本" >&2; exit 1; }

rm -f "${CHUNK_DIR}"/chunk_*.txt "${JOB_DIR}"/*.lsf
split -l "$SAMPLES_PER_JOB" -d -a 4 "${SAMPLE_LIST}.use" "${CHUNK_DIR}/chunk_"
for f in "${CHUNK_DIR}"/chunk_*; do [[ "$f" == *.txt ]] || mv "$f" "${f}.txt"; done

N_JOBS=0
for chunk in "${CHUNK_DIR}"/chunk_*.txt; do
    idx="$(basename "$chunk" .txt)"; idx="${idx#chunk_}"
    job_name="${JOB_PREFIX}_${idx}"
    job_file="${JOB_DIR}/${job_name}.lsf"

    cat > "$job_file" <<LSF
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

CHUNK="${chunk}"
WORKER="${WORKER}"
REPORT_DIR="${REPORT_DIR}"
PARALLEL_BIN="${PARALLEL_BIN}"

echo "[\$(date '+%F %T')] start ${job_name}: \$(wc -l < "\$CHUNK") samples"

"\$PARALLEL_BIN" --no-run-if-empty -j ${JOBS} \\
    "\$WORKER" \\
        --sample-dir {} \\
        --delete-patterns-file "${DELETE_PATTERNS_FILE}" \\
        --keep-patterns-file "${KEEP_PATTERNS_FILE}" \\
        --report "\$REPORT_DIR/{/}.tsv" \\
        --dry-run '${DRY_RUN}' \\
    :::: "\$CHUNK"
rc=\$?

echo "[\$(date '+%F %T')] done ${job_name}: parallel exit=\$rc"
exit 0
LSF
    chmod +x "$job_file"
    N_JOBS=$((N_JOBS + 1))
done

log "已生成 $N_JOBS 个 LSF 作业脚本：$JOB_DIR"

if [[ "$SUBMIT" != "true" ]]; then
    log "未提交（--submit 未指定）。检查无误后执行："
    log "  bash pipe/${SCRIPT_NAME} --dry-run ${DRY_RUN} --limit ${LIMIT} --submit"
    exit 0
fi

command -v bsub >/dev/null 2>&1 || { echo "[ERROR] 未找到 bsub" >&2; exit 1; }
SUBMITTED=0
for job_file in "${JOB_DIR}"/*.lsf; do
    if bsub < "$job_file" > "${job_file}.submit.log" 2>&1; then
        SUBMITTED=$((SUBMITTED + 1))
        log "提交 $(basename "$job_file"): $(head -n1 "${job_file}.submit.log")"
    else
        echo "[ERROR] 提交失败: $job_file" >&2
        cat "${job_file}.submit.log" >&2
    fi
done
log "已提交 $SUBMITTED / $N_JOBS 个作业"
log "作业结束后执行：bash pipe/${SCRIPT_NAME} --collect"
