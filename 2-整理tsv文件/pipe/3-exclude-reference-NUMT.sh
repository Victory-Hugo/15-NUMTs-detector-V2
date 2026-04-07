#!/usr/bin/env bash
set -euo pipefail

SCRIPT_NAME=$(basename "$0")
PROJECT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
CONFIG_PATH="$PROJECT_DIR/conf/3-exclude-reference.yaml"
FORCE_OVERWRITE=0

usage() {
  cat <<'EOF'
使用方法:
  3-exclude-reference-NUMT.sh [--config /path/to/3-exclude-reference.yaml] [--overwrite]

功能:
  调用 python/filter_exclude_reference_numts.py，从严格阈值结果中排除参考 NUMTs 重叠记录。
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG_PATH="$2"
      shift 2
      ;;
    --overwrite)
      FORCE_OVERWRITE=1
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
INPUT_DIR_RAW=$(read_config_value paths input_dir)
OUTPUT_DIR_RAW=$(read_config_value paths output_dir)
REF_FILE_RAW=$(read_config_value paths ref_file)
PYTHON_BIN=$(read_config_value tools python_bin)
EXTENSION=$(read_config_value runtime extension)
OVERWRITE=$(read_config_value runtime overwrite)

require_value "project.base_dir" "$BASE_DIR"
require_value "paths.python_dir" "$PYTHON_DIR_RAW"
require_value "paths.log_dir" "$LOG_DIR_RAW"
require_value "paths.input_dir" "$INPUT_DIR_RAW"
require_value "paths.output_dir" "$OUTPUT_DIR_RAW"
require_value "paths.ref_file" "$REF_FILE_RAW"
require_value "tools.python_bin" "$PYTHON_BIN"
require_value "runtime.extension" "$EXTENSION"

PYTHON_DIR=$(resolve_path "$PYTHON_DIR_RAW" "$BASE_DIR")
LOG_DIR=$(resolve_path "$LOG_DIR_RAW" "$BASE_DIR")
INPUT_DIR=$(resolve_path "$INPUT_DIR_RAW" "$BASE_DIR")
OUTPUT_DIR=$(resolve_path "$OUTPUT_DIR_RAW" "$BASE_DIR")
REF_FILE=$(resolve_path "$REF_FILE_RAW" "$BASE_DIR")
OVERWRITE=$(normalize_bool "$OVERWRITE")
if [[ "$FORCE_OVERWRITE" -eq 1 ]]; then
  OVERWRITE=1
fi

if ! [[ "$EXTENSION" =~ ^[0-9]+$ ]]; then
  echo "Invalid runtime.extension: $EXTENSION" >&2
  exit 1
fi

PYTHON_SCRIPT="$PYTHON_DIR/filter_exclude_reference_numts.py"
if [[ ! -f "$PYTHON_SCRIPT" ]]; then
  echo "Python script not found: $PYTHON_SCRIPT" >&2
  exit 1
fi
if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Input directory not found: $INPUT_DIR" >&2
  exit 1
fi
if [[ ! -f "$REF_FILE" ]]; then
  echo "Reference NUMTs file not found: $REF_FILE" >&2
  exit 1
fi

mkdir -p "$LOG_DIR" "$OUTPUT_DIR"
LOG_FILE="$LOG_DIR/${SCRIPT_NAME}.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "[$(date '+%F %T')] START"
echo "Config: $CONFIG_PATH"
echo "Input:  $INPUT_DIR"
echo "Output: $OUTPUT_DIR"
echo "Ref:    $REF_FILE"
echo "Extension: $EXTENSION"
echo "Overwrite: $OVERWRITE"

CMD=(
  "$PYTHON_BIN" "$PYTHON_SCRIPT"
  --input-dir "$INPUT_DIR"
  --output-dir "$OUTPUT_DIR"
  --ref-file "$REF_FILE"
  --extension "$EXTENSION"
)
if [[ "$OVERWRITE" -eq 1 ]]; then
  CMD+=(--overwrite)
fi

"${CMD[@]}"

echo "[$(date '+%F %T')] SUCCESS"
