#!/usr/bin/env bash
set -euo pipefail

SCRIPT_NAME=$(basename "$0")
PROJECT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
CONFIG_PATH="$PROJECT_DIR/conf/2-1-filter-merge-bed-confident.yaml"

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

  或:
  2-filter-merge-bed-confident.sh \
    [--config /path/to/2-filter-merge-bed-confident.yaml] \
    [--overwrite]
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG_PATH="$2"
      shift 2
      ;;
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

CONFIG_PATH="${CONFIG_PATH//$'\r'/}"

if [[ ! -f "$CONFIG_PATH" ]]; then
  echo "Config file not found: $CONFIG_PATH" >&2
  exit 1
fi

read_config_value() {
  local section="$1"
  local key="$2"
  awk -v target_section="$section" -v target_key="$key" '
    function trim(s) {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", s)
      return s
    }
    function unquote(s, first, last) {
      s = trim(s)
      sub(/[[:space:]]+#.*$/, "", s)
      s = trim(s)
      first = substr(s, 1, 1)
      last = substr(s, length(s), 1)
      if ((first == "'"'"'" && last == "'"'"'") || (first == "\"" && last == "\"")) {
        s = substr(s, 2, length(s) - 2)
      }
      return s
    }
    /^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
    /^[^[:space:]][^:]*:[[:space:]]*$/ {
      current_section = $0
      sub(/:.*/, "", current_section)
      next
    }
    current_section == target_section && /^[[:space:]]{2}[^[:space:]][^:]*:[[:space:]]*/ {
      line = $0
      key_name = line
      sub(/^[[:space:]]+/, "", key_name)
      sub(/:.*/, "", key_name)
      if (key_name == target_key) {
        sub(/^[^:]*:[[:space:]]*/, "", line)
        print unquote(line)
        found = 1
        exit
      }
    }
    END { if (!found) exit 2 }
  ' "$CONFIG_PATH"
}

normalize_bool() {
  case "${1,,}" in
    1|true|yes|y) echo 1 ;;
    0|false|no|n|"") echo 0 ;;
    *) echo "Invalid boolean value: $1" >&2; exit 1 ;;
  esac
}

resolve_path() {
  local path_value="$1"
  local base_dir="$2"
  if [[ "$path_value" = /* ]]; then
    printf '%s\n' "$path_value"
  else
    printf '%s/%s\n' "$base_dir" "$path_value"
  fi
}

require_value() {
  local name="$1"
  local value="$2"
  if [[ -z "$value" ]]; then
    echo "Missing required config value: $name" >&2
    exit 1
  fi
}

BASE_DIR=$(read_config_value project base_dir)
PYTHON_DIR_RAW=$(read_config_value paths python_dir)
LOG_DIR_RAW=$(read_config_value paths log_dir)
PYTHON_BIN=$(read_config_value tools python_bin)
CONFIG_OVERWRITE=$(normalize_bool "$(read_config_value runtime overwrite)")

require_value "project.base_dir" "$BASE_DIR"
require_value "paths.python_dir" "$PYTHON_DIR_RAW"
require_value "paths.log_dir" "$LOG_DIR_RAW"
require_value "tools.python_bin" "$PYTHON_BIN"

PYTHON_DIR=$(resolve_path "$PYTHON_DIR_RAW" "$BASE_DIR")
LOG_DIR=$(resolve_path "$LOG_DIR_RAW" "$BASE_DIR")

if [[ -z "$CONFIDENT_GZ" ]]; then
  CONFIDENT_GZ=$(resolve_path "$(read_config_value paths confident_gz)" "$BASE_DIR")
fi
if [[ -z "$MERGE_BED_GZ" ]]; then
  MERGE_BED_GZ=$(resolve_path "$(read_config_value paths merge_bed_gz)" "$BASE_DIR")
fi
if [[ -z "$OUTPUT_GZ" ]]; then
  OUTPUT_GZ=$(resolve_path "$(read_config_value paths output_gz)" "$BASE_DIR")
fi
if [[ "$OVERWRITE" -eq 0 ]]; then
  OVERWRITE="$CONFIG_OVERWRITE"
fi

require_value "paths.confident_gz" "$CONFIDENT_GZ"
require_value "paths.merge_bed_gz" "$MERGE_BED_GZ"
require_value "paths.output_gz" "$OUTPUT_GZ"

PYTHON_SCRIPT="$PYTHON_DIR/filter_merge_bed_confident.py"

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

mkdir -p "$LOG_DIR" "$(dirname "$OUTPUT_GZ")"
LOG_FILE="$LOG_DIR/${SCRIPT_NAME}.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "[$(date '+%F %T')] START"
echo "Config: $CONFIG_PATH"
echo "Confident: $CONFIDENT_GZ"
echo "Merge bed: $MERGE_BED_GZ"
echo "Output:    $OUTPUT_GZ"
echo "Overwrite: $OVERWRITE"

CMD=(
  "$PYTHON_BIN" "$PYTHON_SCRIPT"
  --confident-gz "$CONFIDENT_GZ"
  --merge-bed-gz "$MERGE_BED_GZ"
  --output-gz "$OUTPUT_GZ"
)

if [[ "$OVERWRITE" -eq 1 ]]; then
  CMD+=(--overwrite)
fi

"${CMD[@]}"

echo "[$(date '+%F %T')] SUCCESS"
