#!/usr/bin/env bash
set -euo pipefail

SCRIPT_NAME=$(basename "$0")
PROJECT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
LOG_DIR="$PROJECT_DIR/log"
LOG_FILE="$LOG_DIR/${SCRIPT_NAME}.log"

mkdir -p "$LOG_DIR"

INPUT_FILE="/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/2-整理tsv文件/output/7-mt_disc_cluster.tsv"
OUTPUT_DIR="/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/2-整理tsv文件/output/"
OUTPUT_FILE="$OUTPUT_DIR/mt.disc.sam.cluster.concat.circos.txt"
FORCE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input)
      INPUT_FILE="$2"
      shift 2
      ;;
    --output-dir)
      OUTPUT_DIR="$2"
      OUTPUT_FILE="$OUTPUT_DIR/mt.disc.sam.cluster.concat.circos.txt"
      shift 2
      ;;
    --output)
      OUTPUT_FILE="$2"
      shift 2
      ;;
    --force)
      FORCE=1
      shift 1
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

INPUT_FILE="${INPUT_FILE//$'\r'/}"
OUTPUT_DIR="${OUTPUT_DIR//$'\r'/}"
OUTPUT_FILE="${OUTPUT_FILE//$'\r'/}"

if [[ -f "$LOG_FILE" ]] && grep -q "SUCCESS" "$LOG_FILE" && [[ "$FORCE" -eq 0 ]]; then
  exit 0
fi

exec > >(tee -a "$LOG_FILE") 2>&1

echo "[$(date '+%F %T')] START"

if [[ ! -f "$INPUT_FILE" ]]; then
  echo "Input file not found: $INPUT_FILE" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

python "$PROJECT_DIR/python/mt.cluster.concat→circos.py" \
  "$INPUT_FILE" \
  "$OUTPUT_FILE"

echo "[$(date '+%F %T')] SUCCESS"
