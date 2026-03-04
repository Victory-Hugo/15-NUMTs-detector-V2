#!/usr/bin/env bash
set -euo pipefail

SCRIPT_NAME=$(basename "$0")
PROJECT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
LOG_DIR="$PROJECT_DIR/log"
LOG_FILE="$LOG_DIR/${SCRIPT_NAME}.log"

mkdir -p "$LOG_DIR"

CONFIG_PATH="$PROJECT_DIR/conf/merge_config.json"
JOBS=""
FORCE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG_PATH="$2"
      shift 2
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

CONFIG_PATH="${CONFIG_PATH//$'\r'/}"

if [[ ! -f "$CONFIG_PATH" ]]; then
  echo "Config file not found: $CONFIG_PATH" >&2
  exit 1
fi

mapfile -t CONFIG_RUNTIME < <(
  python "$PROJECT_DIR/python/merge_runtime_config.py" --config "$CONFIG_PATH"
)

JOBS="${CONFIG_RUNTIME[0]}"
FORCE="${CONFIG_RUNTIME[1]}"

if ! [[ "$JOBS" =~ ^[0-9]+$ ]] || [[ "$JOBS" -lt 1 ]]; then
  echo "Invalid jobs value in config: $JOBS" >&2
  exit 1
fi

if ! [[ "$FORCE" =~ ^[01]$ ]]; then
  echo "Invalid force value in config: $FORCE" >&2
  exit 1
fi

if [[ -f "$LOG_FILE" ]] && grep -q "SUCCESS" "$LOG_FILE" && [[ "$FORCE" -eq 0 ]]; then
  echo "[$(date '+%F %T')] SKIP (force=0 and existing SUCCESS in $LOG_FILE)"
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
LISTS_DIR="$TMP_DIR/merge_filelists"
FORCE_FLAG=""
if [[ "$FORCE" -eq 1 ]]; then
  FORCE_FLAG="--force"
fi

echo "[$(date '+%F %T')] STEP1 build merge file lists"
python "$PROJECT_DIR/python/build_merge_filelists.py" \
  --config "$CONFIG_PATH" \
  --keys-file "$KEYS_FILE" \
  --lists-dir "$LISTS_DIR"

echo "[$(date '+%F %T')] STEP2 merge using prepared txt lists"
if [[ "$auto_parallel_cmd" == "parallel" ]]; then
  parallel -j "$JOBS" --arg-file "$KEYS_FILE" \
    python -m merge_tables --config "$CONFIG_PATH" --suffix-key {} --files-list "$LISTS_DIR/{}.txt" $FORCE_FLAG
else
  xargs -I {} -P "$JOBS" python -m merge_tables --config "$CONFIG_PATH" --suffix-key {} --files-list "$LISTS_DIR/{}.txt" $FORCE_FLAG \
    < "$KEYS_FILE"
fi

echo "[$(date '+%F %T')] SUCCESS"
