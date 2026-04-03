#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
PYTHON_SCRIPT="$PROJECT_DIR/python/filter_merge_bed_confident.py"

CONFIDENT_GZ=""
MERGE_BED_GZ=""
OUTPUT_GZ=""
OVERWRITE=0

usage() {
  cat <<'EOF'
使用方法:
  2-filter-merge-bed-confident.sh \
    --confident-gz /path/to/3-confident_breakpoints.tsv.gz \
    --merge-bed-gz /path/to/8-merge_bed.tsv.gz \
    --output-gz /path/to/9-merge_bed_confident.tsv.gz \
    [--overwrite]
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --confident-gz)
      CONFIDENT_GZ="$2"
      shift 2
      ;;
    --merge-bed-gz)
      MERGE_BED_GZ="$2"
      shift 2
      ;;
    --output-gz)
      OUTPUT_GZ="$2"
      shift 2
      ;;
    --overwrite)
      OVERWRITE=1
      shift 1
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if [[ -z "$CONFIDENT_GZ" || -z "$MERGE_BED_GZ" || -z "$OUTPUT_GZ" ]]; then
  usage >&2
  exit 1
fi

if [[ ! -f "$PYTHON_SCRIPT" ]]; then
  echo "Python script not found: $PYTHON_SCRIPT" >&2
  exit 1
fi

if [[ ! -f "$CONFIDENT_GZ" ]]; then
  echo "Confident breakpoint file not found: $CONFIDENT_GZ" >&2
  exit 1
fi

if [[ ! -f "$MERGE_BED_GZ" ]]; then
  echo "Merge bed file not found: $MERGE_BED_GZ" >&2
  exit 1
fi

CMD=(
  python3 "$PYTHON_SCRIPT"
  --confident-gz "$CONFIDENT_GZ"
  --merge-bed-gz "$MERGE_BED_GZ"
  --output-gz "$OUTPUT_GZ"
)

if [[ "$OVERWRITE" -eq 1 ]]; then
  CMD+=(--overwrite)
fi

"${CMD[@]}"
