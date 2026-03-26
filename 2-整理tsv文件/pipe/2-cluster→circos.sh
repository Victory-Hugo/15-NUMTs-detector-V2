#!/usr/bin/env bash
set -euo pipefail

SCRIPT_NAME=$(basename "$0")
PROJECT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
LOG_DIR="$PROJECT_DIR/log"
LOG_FILE="$LOG_DIR/${SCRIPT_NAME}.log"

mkdir -p "$LOG_DIR"

CONFIG_PATH="$PROJECT_DIR/conf/merge_config.yaml"
INPUT_FILE=""
OUTPUT_DIR=""
OUTPUT_FILE=""
FORCE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG_PATH="$2"
      shift 2
      ;;
    --input)
      INPUT_FILE="$2"
      shift 2
      ;;
    --output-dir)
      OUTPUT_DIR="$2"
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

CONFIG_PATH="${CONFIG_PATH//$'\r'/}"
INPUT_FILE="${INPUT_FILE//$'\r'/}"
OUTPUT_DIR="${OUTPUT_DIR//$'\r'/}"
OUTPUT_FILE="${OUTPUT_FILE//$'\r'/}"

if [[ ! -f "$CONFIG_PATH" ]]; then
  echo "Config file not found: $CONFIG_PATH" >&2
  exit 1
fi

parse_force() {
  case "${1,,}" in
    1|true|yes|y) echo 1 ;;
    0|false|no|n|"") echo 0 ;;
    *) echo "Invalid force value in config: $1" >&2; exit 1 ;;
  esac
}

read_config() {
  awk '
    function trim(s) {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", s)
      return s
    }
    function unquote(s, first, last) {
      s = trim(s)
      sub(/[[:space:]]+#.*$/, "", s)
      first = substr(s, 1, 1)
      last = substr(s, length(s), 1)
      if ((first == "'"'"'" && last == "'"'"'") || (first == "\"" && last == "\"")) {
        s = substr(s, 2, length(s) - 2)
      }
      return s
    }
    /^[[:space:]]*#/ || /^[[:space:]]*$/ {
      next
    }
    /^[^[:space:]][^:]*:[[:space:]]*.*$/ {
      top_key = $0
      sub(/:.*/, "", top_key)
      line = $0
      sub(/^[^:]*:[[:space:]]*/, "", line)
      if (top_key == "output_dir") {
        print "output_dir\t" unquote(line)
      } else if (top_key == "force") {
        print "force\t" unquote(line)
      } else if (top_key == "circos") {
        in_circos = 1
      } else {
        in_circos = 0
      }
      next
    }
    in_circos && /^[[:space:]]{2}enabled:[[:space:]]*/ {
      line = $0
      sub(/^[^:]*:[[:space:]]*/, "", line)
      print "circos.enabled\t" unquote(line)
      next
    }
    in_circos && /^[[:space:]]{2}input_name:[[:space:]]*/ {
      line = $0
      sub(/^[^:]*:[[:space:]]*/, "", line)
      print "circos.input_name\t" unquote(line)
      next
    }
    in_circos && /^[[:space:]]{2}output_name:[[:space:]]*/ {
      line = $0
      sub(/^[^:]*:[[:space:]]*/, "", line)
      print "circos.output_name\t" unquote(line)
      next
    }
  ' "$CONFIG_PATH"
}

CONFIG_OUTPUT_DIR=""
CONFIG_FORCE=""
CIRCOS_ENABLED="true"
CIRCOS_INPUT_NAME=""
CIRCOS_OUTPUT_NAME=""

while IFS=$'\t' read -r key value; do
  case "$key" in
    output_dir) CONFIG_OUTPUT_DIR="$value" ;;
    force) CONFIG_FORCE="$value" ;;
    circos.enabled) CIRCOS_ENABLED="$value" ;;
    circos.input_name) CIRCOS_INPUT_NAME="$value" ;;
    circos.output_name) CIRCOS_OUTPUT_NAME="$value" ;;
  esac
done < <(read_config)

if [[ -z "$CONFIG_OUTPUT_DIR" ]]; then
  echo "Missing output_dir in $CONFIG_PATH" >&2
  exit 1
fi

if [[ -z "$CONFIG_FORCE" ]]; then
  CONFIG_FORCE=0
fi

if [[ -z "$FORCE" ]]; then
  FORCE=$(parse_force "$CONFIG_FORCE")
fi

case "${CIRCOS_ENABLED,,}" in
  1|true|yes|y) ;;
  0|false|no|n|"")
    echo "[$(date '+%F %T')] SKIP (circos.enabled=false)"
    exit 0
    ;;
  *)
    echo "Invalid circos.enabled value in config: $CIRCOS_ENABLED" >&2
    exit 1
    ;;
esac

if [[ -z "$OUTPUT_DIR" ]]; then
  OUTPUT_DIR="$CONFIG_OUTPUT_DIR"
fi

if [[ -z "$INPUT_FILE" ]]; then
  if [[ -z "$CIRCOS_INPUT_NAME" ]]; then
    echo "Missing circos.input_name in $CONFIG_PATH" >&2
    exit 1
  fi
  INPUT_FILE="$OUTPUT_DIR/$CIRCOS_INPUT_NAME"
fi

if [[ -z "$OUTPUT_FILE" ]]; then
  if [[ -z "$CIRCOS_OUTPUT_NAME" ]]; then
    echo "Missing circos.output_name in $CONFIG_PATH" >&2
    exit 1
  fi
  OUTPUT_FILE="$OUTPUT_DIR/$CIRCOS_OUTPUT_NAME"
fi

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

python3 "$PROJECT_DIR/python/mt.cluster.concat→circos.py" \
  "$INPUT_FILE" \
  "$OUTPUT_FILE"

echo "[$(date '+%F %T')] SUCCESS"
