#!/usr/bin/env bash
set -euo pipefail

SCRIPT_NAME=$(basename "$0")
PROJECT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
LOG_DIR="$PROJECT_DIR/log"
LOG_FILE="$LOG_DIR/${SCRIPT_NAME}.log"

mkdir -p "$LOG_DIR"

CONFIG_PATH="$PROJECT_DIR/conf/merge_config.json"
JOBS=4
FORCE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG_PATH="$2"
      shift 2
      ;;
    --jobs)
      JOBS="$2"
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

CONFIG_PATH="${CONFIG_PATH//$'\r'/}"

if [[ -f "$LOG_FILE" ]] && grep -q "SUCCESS" "$LOG_FILE" && [[ "$FORCE" -eq 0 ]]; then
  exit 0
fi

exec > >(tee -a "$LOG_FILE") 2>&1

echo "[$(date '+%F %T')] START"

auto_parallel_cmd="parallel"
if ! command -v parallel >/dev/null 2>&1; then
  auto_parallel_cmd="xargs"
fi

export PYTHONPATH="$PROJECT_DIR/python:${PYTHONPATH:-}"

TMP_DIR="$PROJECT_DIR/tmp"
mkdir -p "$TMP_DIR"
KEYS_FILE="$TMP_DIR/merge_keys.txt"
FORCE_FLAG=""
if [[ "$FORCE" -eq 1 ]]; then
  FORCE_FLAG="--force"
fi

if [[ "$auto_parallel_cmd" == "parallel" ]]; then
  python -m merge_tables --config "$CONFIG_PATH" --list-keys > "$KEYS_FILE"
  parallel -j "$JOBS" --arg-file "$KEYS_FILE" \
    python -m merge_tables --config "$CONFIG_PATH" --suffix-key {} $FORCE_FLAG
else
  python -m merge_tables --config "$CONFIG_PATH" --list-keys > "$KEYS_FILE"
  xargs -I {} -P "$JOBS" python -m merge_tables --config "$CONFIG_PATH" --suffix-key {} $FORCE_FLAG \
    < "$KEYS_FILE"
fi

echo "[$(date '+%F %T')] SUCCESS"
