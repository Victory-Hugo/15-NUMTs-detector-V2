#!/usr/bin/env bash
# script/clean_numts_sample.sh
# 清理单个 1-NUMTs 样本目录：删除可再生的下游文件，并把 bwa.sam 压缩为 BAM。
#
# 安全约定：
#   1. keep 优先于 delete —— 任何匹配保留清单的文件绝不会被删除。
#   2. 既不匹配 keep 也不匹配 delete 的文件记为 unknown，原样保留。
#   3. SAM 先转 BAM 并校验记录数，一致后才删除原 SAM。
#
# 报告为 tidy long 格式，每行一条：sample_id \t action \t category \t files \t bytes
set -uo pipefail

SAMPLE_DIR=""
DELETE_PATTERNS_FILE=""
KEEP_PATTERNS_FILE=""
REPORT=""
DRY_RUN="true"
SAMTOOLS_BIN="samtools"
COMPRESS_ENABLED="true"
SOURCE_SUFFIX="_all_regions.bwa.sam"
TARGET_SUFFIX="_all_regions.bwa.bam"
COMPRESS_THREADS="2"
VERIFY="true"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample-dir)         SAMPLE_DIR="$2";        shift 2 ;;
        --delete-patterns-file) DELETE_PATTERNS_FILE="$2"; shift 2 ;;
        --keep-patterns-file)   KEEP_PATTERNS_FILE="$2";   shift 2 ;;
        --report)             REPORT="$2";            shift 2 ;;
        --dry-run)            DRY_RUN="$2";           shift 2 ;;
        --samtools)           SAMTOOLS_BIN="$2";      shift 2 ;;
        --compress)           COMPRESS_ENABLED="$2";  shift 2 ;;
        --source-suffix)      SOURCE_SUFFIX="$2";     shift 2 ;;
        --target-suffix)      TARGET_SUFFIX="$2";     shift 2 ;;
        --compress-threads)   COMPRESS_THREADS="$2";  shift 2 ;;
        --verify)             VERIFY="$2";            shift 2 ;;
        *) echo "[ERROR] Unknown argument: $1" >&2; exit 2 ;;
    esac
done

if [[ -z "$SAMPLE_DIR" || ! -d "$SAMPLE_DIR" ]]; then
    echo "[ERROR] Invalid --sample-dir: $SAMPLE_DIR" >&2
    exit 2
fi
if [[ -z "$REPORT" ]]; then
    echo "[ERROR] Missing --report" >&2
    exit 2
fi

SAMPLE_ID="$(basename "$SAMPLE_DIR")"

# 模式清单从文件读取（每行一个），避免 GNU parallel 重新解析命令行时拆散带空格的参数。
load_patterns() {
    local file="$1" line
    [[ -n "$file" && -f "$file" ]] || return 0
    while IFS= read -r line; do
        line="${line%$'\r'}"
        [[ -z "$line" || "$line" == \#* ]] && continue
        printf '%s\n' "$line"
    done < "$file"
}

mapfile -t DELETE_ARR < <(load_patterns "$DELETE_PATTERNS_FILE")
mapfile -t KEEP_ARR   < <(load_patterns "$KEEP_PATTERNS_FILE")

if [[ ${#DELETE_ARR[@]} -eq 0 ]]; then
    echo "[ERROR] 删除模式清单为空: $DELETE_PATTERNS_FILE" >&2
    exit 2
fi

# 归类：返回 keep / delete / unknown。keep 优先。
classify() {
    local name="$1" pat
    for pat in "${KEEP_ARR[@]}"; do
        [[ -n "$pat" && "$name" == $pat ]] && { echo "keep"; return; }
    done
    for pat in "${DELETE_ARR[@]}"; do
        [[ -n "$pat" && "$name" == $pat ]] && { echo "delete"; return; }
    done
    echo "unknown"
}

# 把文件名归入一个可读的统计类别，便于报告聚合。
categorize() {
    local name="$1"
    case "$name" in
        *_all_regions.fasta.cap.*)  echo "all_regions.fasta.cap" ;;
        *_all_regions.fasta)        echo "all_regions.fasta" ;;
        *_all_regions.psl)          echo "all_regions.psl" ;;
        *_all_regions.bwa.sam)      echo "all_regions.bwa.sam" ;;
        *_all_regions.bwa.bam)      echo "all_regions.bwa.bam" ;;
        *.mt.disc.sam)              echo "mt.disc.sam" ;;
        *.mt.split.sam)             echo "mt.split.sam" ;;
        *.Breakpoints.old.tsv)      echo "Breakpoints.old.tsv" ;;
        *.cluster.old.tsv)          echo "cluster.old.tsv" ;;
        *.AllBreakpoints.tsv)       echo "AllBreakpoints.tsv" ;;
        *.ConfidentBreakpoints.tsv) echo "ConfidentBreakpoints.tsv" ;;
        *.breakpointINPUT.tsv)      echo "breakpointINPUT.tsv" ;;
        *.cluster.tsv)              echo "cluster.tsv" ;;
        *.cluster.summary.tsv)      echo "cluster.summary.tsv" ;;
        *.psl)                      echo "region.psl" ;;
        *)                          echo "other" ;;
    esac
}

declare -A ACT_FILES ACT_BYTES

record() {
    local action="$1" category="$2" bytes="$3"
    local key="${action}|${category}"
    ACT_FILES["$key"]=$(( ${ACT_FILES["$key"]:-0} + 1 ))
    ACT_BYTES["$key"]=$(( ${ACT_BYTES["$key"]:-0} + bytes ))
}

STATUS="ok"

# ------------------------------------------------------------------ #
# 阶段 1：删除可再生文件。先删再压缩，以便在磁盘吃紧时腾出写 BAM 的空间。
# ------------------------------------------------------------------ #
while IFS=$'\t' read -r fsize fname; do
    [[ -z "$fname" ]] && continue
    verdict="$(classify "$fname")"
    category="$(categorize "$fname")"

    case "$verdict" in
        delete)
            if [[ "$DRY_RUN" == "true" ]]; then
                record "would_delete" "$category" "$fsize"
            elif rm -f -- "${SAMPLE_DIR}/${fname}"; then
                record "deleted" "$category" "$fsize"
            else
                record "delete_failed" "$category" "$fsize"
                STATUS="delete_failed"
            fi
            ;;
        keep)    record "kept"    "$category" "$fsize" ;;
        unknown) record "unknown" "$category" "$fsize" ;;
    esac
done < <(find "$SAMPLE_DIR" -maxdepth 1 -type f -printf '%s\t%f\n' 2>/dev/null)

# ------------------------------------------------------------------ #
# 阶段 2：SAM -> BAM
# ------------------------------------------------------------------ #
if [[ "$COMPRESS_ENABLED" == "true" ]]; then
    while IFS= read -r sam_path; do
        [[ -z "$sam_path" ]] && continue
        bam_path="${sam_path%${SOURCE_SUFFIX}}${TARGET_SUFFIX}"
        sam_bytes=$(stat -c '%s' "$sam_path" 2>/dev/null || echo 0)

        if [[ "$DRY_RUN" == "true" ]]; then
            record "would_compress" "all_regions.bwa.sam" "$sam_bytes"
            continue
        fi

        if [[ -s "$bam_path" ]]; then
            record "compress_skipped_exists" "all_regions.bwa.bam" 0
            continue
        fi

        tmp_bam="${bam_path}.tmp"
        if ! "$SAMTOOLS_BIN" view -@ "$COMPRESS_THREADS" -b -o "$tmp_bam" "$sam_path" 2>/dev/null; then
            rm -f -- "$tmp_bam"
            record "compress_failed" "all_regions.bwa.sam" "$sam_bytes"
            STATUS="compress_failed"
            continue
        fi

        if [[ "$VERIFY" == "true" ]]; then
            n_sam=$("$SAMTOOLS_BIN" view -c "$sam_path" 2>/dev/null || echo -1)
            n_bam=$("$SAMTOOLS_BIN" view -c "$tmp_bam" 2>/dev/null || echo -2)
            if [[ "$n_sam" != "$n_bam" ]]; then
                rm -f -- "$tmp_bam"
                record "verify_failed" "all_regions.bwa.sam" "$sam_bytes"
                STATUS="verify_failed"
                continue
            fi
        fi

        mv -f -- "$tmp_bam" "$bam_path"
        bam_bytes=$(stat -c '%s' "$bam_path" 2>/dev/null || echo 0)
        rm -f -- "$sam_path"
        record "compressed_from_sam" "all_regions.bwa.sam" "$sam_bytes"
        record "compressed_to_bam"   "all_regions.bwa.bam" "$bam_bytes"
    done < <(find "$SAMPLE_DIR" -maxdepth 1 -type f -name "*${SOURCE_SUFFIX}" 2>/dev/null)
fi

# ------------------------------------------------------------------ #
# 写报告（tidy long）
# ------------------------------------------------------------------ #
{
    for key in "${!ACT_FILES[@]}"; do
        printf '%s\t%s\t%s\t%s\t%s\n' \
            "$SAMPLE_ID" "${key%%|*}" "${key##*|}" "${ACT_FILES[$key]}" "${ACT_BYTES[$key]}"
    done
    printf '%s\t%s\t%s\t%s\t%s\n' "$SAMPLE_ID" "status" "$STATUS" 1 0
} > "$REPORT"

[[ "$STATUS" == "ok" ]] || exit 1
exit 0
