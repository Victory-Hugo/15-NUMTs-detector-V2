#!/usr/bin/env bash
# script/clean_vardetect_sample.sh
# 清理单个 1-5-Vardetect 样本目录：删除可由 *aln.fa + *_numts.bed 再生的逐碱基 TSV
# 及可秒级重建的拼接 FASTA。递归处理 alnHuman/ alnHumanChimp/ assembly/ psl*/ 子目录。
#
# 安全约定：
#   1. keep 优先于 delete —— 任何匹配保留清单的文件绝不会被删除。
#   2. 既不匹配 keep 也不匹配 delete 的文件记为 unknown，原样保留。
#
# 报告为 tidy long 格式：sample_id \t action \t category \t files \t bytes
set -uo pipefail

SAMPLE_DIR=""
DELETE_PATTERNS_FILE=""
KEEP_PATTERNS_FILE=""
REPORT=""
DRY_RUN="true"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample-dir)      SAMPLE_DIR="$2";      shift 2 ;;
        --delete-patterns-file) DELETE_PATTERNS_FILE="$2"; shift 2 ;;
        --keep-patterns-file)   KEEP_PATTERNS_FILE="$2";   shift 2 ;;
        --report)          REPORT="$2";          shift 2 ;;
        --dry-run)         DRY_RUN="$2";         shift 2 ;;
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

categorize() {
    local name="$1"
    case "$name" in
        *.humanMTaln.fa)                echo "humanMTaln.fa" ;;
        *.humanchimpMTaln.fa)           echo "humanchimpMTaln.fa" ;;
        *.humanMTaln.fa.numt.tsv)       echo "humanMTaln.numt.tsv" ;;
        *.humanchimpMTaln.fa.numt.tsv)  echo "humanchimpMTaln.numt.tsv" ;;
        *.numtVar.tsv)                  echo "numtVar.tsv" ;;
        *.numtVarFilterPos.tsv)         echo "numtVarFilterPos.tsv" ;;
        *.full.tsv)                     echo "full.tsv" ;;
        *.numtDhumanChimp.sum.tsv)      echo "numtDhumanChimp.sum.tsv" ;;
        *.numtDhumanChimp.tsv)          echo "numtDhumanChimp.tsv" ;;
        *_all_regions.filtered.fasta)   echo "filtered.fasta" ;;
        *_all_regions.humanchimpMT.fasta) echo "humanchimpMT.fasta" ;;
        *_all_regions.humanMT.fasta)    echo "humanMT.fasta" ;;
        *_numts.bed)                    echo "numts.bed" ;;
        *.human.psl)                    echo "human.psl" ;;
        *.chimp.psl)                    echo "chimp.psl" ;;
        *.fasta.cap|*.fasta.cap.*)      echo "fasta.cap" ;;
        fastaFiles.list)                echo "fastaFiles.list" ;;
        core.[0-9]*)                    echo "cap3.coredump" ;;
        *)                              echo "other" ;;
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

while IFS=$'\t' read -r fsize frel; do
    [[ -z "$frel" ]] && continue
    fname="$(basename "$frel")"
    verdict="$(classify "$fname")"
    category="$(categorize "$fname")"

    case "$verdict" in
        delete)
            if [[ "$DRY_RUN" == "true" ]]; then
                record "would_delete" "$category" "$fsize"
            elif rm -f -- "${SAMPLE_DIR}/${frel}"; then
                record "deleted" "$category" "$fsize"
            else
                record "delete_failed" "$category" "$fsize"
                STATUS="delete_failed"
            fi
            ;;
        keep)    record "kept"    "$category" "$fsize" ;;
        unknown) record "unknown" "$category" "$fsize" ;;
    esac
done < <(find "$SAMPLE_DIR" -type f -printf '%s\t%P\n' 2>/dev/null)

{
    for key in "${!ACT_FILES[@]}"; do
        printf '%s\t%s\t%s\t%s\t%s\n' \
            "$SAMPLE_ID" "${key%%|*}" "${key##*|}" "${ACT_FILES[$key]}" "${ACT_BYTES[$key]}"
    done
    printf '%s\t%s\t%s\t%s\t%s\n' "$SAMPLE_ID" "status" "$STATUS" 1 0
} > "$REPORT"

[[ "$STATUS" == "ok" ]] || exit 1
exit 0
