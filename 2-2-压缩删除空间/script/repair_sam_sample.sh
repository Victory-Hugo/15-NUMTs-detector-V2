#!/usr/bin/env bash
# script/repair_sam_sample.sh
# 把一个截断的 SAM 中可解析的部分转为规范 BAM，丢弃损坏尾部。
#
# 做法：samtools view -h 读到出错行即停止，但此前已输出的记录是完整有效的；
#       把这段管道输出转成 BAM。转换后校验 BAM 记录数 > 0 才替换原 SAM。
#
# 报告：sample_id \t action \t category \t files \t bytes
set -uo pipefail

SAMPLE_DIR=""
REPORT=""
DRY_RUN="true"
KEEP_SAM="false"
SAMTOOLS_BIN="samtools"
SOURCE_SUFFIX="_all_regions.bwa.sam"
TARGET_SUFFIX="_all_regions.bwa.bam"
THREADS="2"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample-dir)    SAMPLE_DIR="$2";    shift 2 ;;
        --report)        REPORT="$2";        shift 2 ;;
        --dry-run)       DRY_RUN="$2";       shift 2 ;;
        --keep-sam)      KEEP_SAM="$2";      shift 2 ;;
        --samtools)      SAMTOOLS_BIN="$2";  shift 2 ;;
        --source-suffix) SOURCE_SUFFIX="$2"; shift 2 ;;
        --target-suffix) TARGET_SUFFIX="$2"; shift 2 ;;
        --threads)       THREADS="$2";       shift 2 ;;
        *) echo "[ERROR] Unknown argument: $1" >&2; exit 2 ;;
    esac
done

if [[ -z "$SAMPLE_DIR" || ! -d "$SAMPLE_DIR" ]]; then
    echo "[ERROR] Invalid --sample-dir: $SAMPLE_DIR" >&2
    exit 2
fi
[[ -n "$REPORT" ]] || { echo "[ERROR] Missing --report" >&2; exit 2; }

SAMPLE_ID="$(basename "$SAMPLE_DIR")"
SAM_PATH="${SAMPLE_DIR}/${SAMPLE_ID}${SOURCE_SUFFIX}"
BAM_PATH="${SAMPLE_DIR}/${SAMPLE_ID}${TARGET_SUFFIX}"

emit() {
    printf '%s\t%s\t%s\t%s\t%s\n' "$SAMPLE_ID" "$1" "$2" "$3" "$4" >> "$REPORT"
}

: > "$REPORT"

if [[ ! -s "$SAM_PATH" ]]; then
    emit "skipped_empty_or_missing" "all_regions.bwa.sam" 1 0
    emit "status" "skipped" 1 0
    exit 0
fi

if [[ -s "$BAM_PATH" ]]; then
    emit "skipped_bam_exists" "all_regions.bwa.bam" 1 0
    emit "status" "skipped" 1 0
    exit 0
fi

SAM_BYTES=$(stat -c '%s' "$SAM_PATH" 2>/dev/null || echo 0)

# 可解析记录数（samtools 在出错行停止，之前的行已计入）
READABLE=$("$SAMTOOLS_BIN" view "$SAM_PATH" 2>/dev/null | wc -l)

if [[ "$DRY_RUN" == "true" ]]; then
    emit "would_repair" "all_regions.bwa.sam" 1 "$SAM_BYTES"
    emit "readable_records" "all_regions.bwa.sam" "$READABLE" 0
    emit "status" "dry_run" 1 0
    exit 0
fi

if [[ "$READABLE" -le 0 ]]; then
    emit "repair_failed_no_readable_records" "all_regions.bwa.sam" 1 "$SAM_BYTES"
    emit "status" "repair_failed" 1 0
    exit 1
fi

TMP_BAM="${BAM_PATH}.tmp"
# view -h 的退出码必然非零（尾部损坏），因此只看产物是否有效，不看退出码
"$SAMTOOLS_BIN" view -h "$SAM_PATH" 2>/dev/null \
    | "$SAMTOOLS_BIN" view -@ "$THREADS" -b -o "$TMP_BAM" - 2>/dev/null

if [[ ! -s "$TMP_BAM" ]]; then
    rm -f -- "$TMP_BAM"
    emit "repair_failed_empty_bam" "all_regions.bwa.sam" 1 "$SAM_BYTES"
    emit "status" "repair_failed" 1 0
    exit 1
fi

if ! "$SAMTOOLS_BIN" quickcheck "$TMP_BAM" 2>/dev/null; then
    rm -f -- "$TMP_BAM"
    emit "repair_failed_quickcheck" "all_regions.bwa.sam" 1 "$SAM_BYTES"
    emit "status" "repair_failed" 1 0
    exit 1
fi

N_BAM=$("$SAMTOOLS_BIN" view -c "$TMP_BAM" 2>/dev/null || echo 0)
if [[ "$N_BAM" -ne "$READABLE" ]]; then
    rm -f -- "$TMP_BAM"
    emit "repair_failed_count_mismatch" "all_regions.bwa.sam" 1 "$SAM_BYTES"
    emit "status" "repair_failed" 1 0
    exit 1
fi

mv -f -- "$TMP_BAM" "$BAM_PATH"
BAM_BYTES=$(stat -c '%s' "$BAM_PATH" 2>/dev/null || echo 0)
emit "repaired_from_sam" "all_regions.bwa.sam" 1 "$SAM_BYTES"
emit "repaired_to_bam"   "all_regions.bwa.bam" 1 "$BAM_BYTES"
emit "recovered_records" "all_regions.bwa.bam" "$N_BAM" 0

if [[ "$KEEP_SAM" == "true" ]]; then
    emit "kept_original_sam" "all_regions.bwa.sam" 1 "$SAM_BYTES"
else
    rm -f -- "$SAM_PATH"
    emit "deleted_original_sam" "all_regions.bwa.sam" 1 "$SAM_BYTES"
fi

emit "status" "ok" 1 0
exit 0
