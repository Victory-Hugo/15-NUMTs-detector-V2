#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: $0 --input <circos_10K.txt> --links-output <links.txt> --mt-output <mtFragments.txt>" >&2
}

INPUT_FILE=""
LINKS_OUTPUT=""
MT_OUTPUT=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)
            INPUT_FILE="$2"
            shift 2
            ;;
        --links-output)
            LINKS_OUTPUT="$2"
            shift 2
            ;;
        --mt-output)
            MT_OUTPUT="$2"
            shift 2
            ;;
        *)
            usage
            exit 1
            ;;
    esac
done

if [[ -z "$INPUT_FILE" || -z "$LINKS_OUTPUT" || -z "$MT_OUTPUT" ]]; then
    usage
    exit 1
fi

if [[ ! -f "$INPUT_FILE" ]]; then
    echo "[ERROR] Input file not found: $INPUT_FILE" >&2
    exit 1
fi

mkdir -p "$(dirname "$LINKS_OUTPUT")" "$(dirname "$MT_OUTPUT")"

awk '
NF == 0 {
    next
}

NF != 6 {
    printf "[ERROR] Expected 6 columns at line %d, got %d\n", NR, NF > "/dev/stderr"
    exit 1
}

{
    nuclear_chr = $1
    nuclear_start = $2 + 0
    nuclear_end = $3 + 0
    mt_chr = $4
    mt_start = $5 + 0
    mt_end = $6 + 0

    if (nuclear_start <= nuclear_end) {
        norm_start = nuclear_start
        norm_end = nuclear_end
    } else {
        norm_start = nuclear_end
        norm_end = nuclear_start
    }

    printf "%s\t%d\t%d\t%s\t%d\t%d\n", nuclear_chr, norm_start, norm_end, mt_chr, mt_start, mt_end > links_file
    printf "%s\t%d\t%d\n", mt_chr, mt_start, mt_end > mt_file
}
' links_file="$LINKS_OUTPUT" mt_file="$MT_OUTPUT" "$INPUT_FILE"

echo "[OK] Wrote links file: $LINKS_OUTPUT"
echo "[OK] Wrote mt fragments file: $MT_OUTPUT"

