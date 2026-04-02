#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
CONFIG_PATH="${1:-$PROJECT_ROOT/conf/Config.yaml}"

if [[ ! -f "$CONFIG_PATH" ]]; then
    echo "[ERROR] Config file not found: $CONFIG_PATH" >&2
    exit 1
fi

source "$PROJECT_ROOT/script/load_config.sh" "$CONFIG_PATH"
bash "$PROJECT_ROOT/script/check_env.sh" "$PROJECT_ROOT" "$CONFIG_PATH"

require_var() {
    local var_name="$1"
    if [[ -z "${!var_name:-}" ]]; then
        echo "[ERROR] Missing required config value: $var_name" >&2
        exit 1
    fi
}

require_var CFG_PATHS_INPUT_FILE
require_var CFG_PATHS_REFERENCE_DIR
require_var CFG_PATHS_OUTPUT_DIR
require_var CFG_PATHS_RUNTIME_DIR
require_var CFG_PATHS_CIRCOS_TEMPLATE
require_var CFG_PATHS_RENDERED_CIRCOS_CONFIG
require_var CFG_PATHS_LINKS_FILE
require_var CFG_PATHS_MT_FRAGMENTS_FILE
require_var CFG_PATHS_PNG_OUTPUT
require_var CFG_PATHS_SVG_OUTPUT
require_var CFG_TOOLS_CONDA_BIN
require_var CFG_TOOLS_CIRCOS_BIN
require_var CFG_TOOLS_CIRCOS_ENV_NAME
require_var CFG_TOOLS_CIRCOS_ETC_DIR

INPUT_FILE="$PROJECT_ROOT/$CFG_PATHS_INPUT_FILE"
REFERENCE_DIR="$PROJECT_ROOT/$CFG_PATHS_REFERENCE_DIR"
OUTPUT_DIR="$PROJECT_ROOT/$CFG_PATHS_OUTPUT_DIR"
RUNTIME_DIR="$PROJECT_ROOT/$CFG_PATHS_RUNTIME_DIR"
TEMPLATE_PATH="$PROJECT_ROOT/$CFG_PATHS_CIRCOS_TEMPLATE"
RENDERED_CIRCOS_CONFIG="$PROJECT_ROOT/$CFG_PATHS_RENDERED_CIRCOS_CONFIG"
LINKS_FILE="$PROJECT_ROOT/$CFG_PATHS_LINKS_FILE"
MT_FRAGMENTS_FILE="$PROJECT_ROOT/$CFG_PATHS_MT_FRAGMENTS_FILE"
PNG_OUTPUT="$PROJECT_ROOT/$CFG_PATHS_PNG_OUTPUT"
SVG_OUTPUT="$PROJECT_ROOT/$CFG_PATHS_SVG_OUTPUT"
CONDA_BIN="$CFG_TOOLS_CONDA_BIN"
CIRCOS_BIN="$CFG_TOOLS_CIRCOS_BIN"
CIRCOS_ENV_NAME="$CFG_TOOLS_CIRCOS_ENV_NAME"
CIRCOS_ETC_DIR="$CFG_TOOLS_CIRCOS_ETC_DIR"

mkdir -p "$OUTPUT_DIR" "$RUNTIME_DIR" "$RUNTIME_DIR/data"

if [[ ! -d "$REFERENCE_DIR" ]]; then
    echo "[ERROR] Reference directory not found: $REFERENCE_DIR" >&2
    exit 1
fi

bash "$PROJECT_ROOT/script/build_runtime_inputs.sh" \
    --input "$INPUT_FILE" \
    --links-output "$LINKS_FILE" \
    --mt-output "$MT_FRAGMENTS_FILE"

sed \
    -e "s|__CIRCOS_ETC_DIR__|$CIRCOS_ETC_DIR|g" \
    -e "s|__PROJECT_ROOT__|$PROJECT_ROOT|g" \
    "$TEMPLATE_PATH" > "$RENDERED_CIRCOS_CONFIG"

"$CONDA_BIN" run -n "$CIRCOS_ENV_NAME" "$CIRCOS_BIN" -conf "$RENDERED_CIRCOS_CONFIG" -outputdir "$OUTPUT_DIR"

if [[ ! -f "$PNG_OUTPUT" ]]; then
    echo "[ERROR] Missing expected PNG output: $PNG_OUTPUT" >&2
    exit 1
fi

if [[ ! -f "$SVG_OUTPUT" ]]; then
    echo "[ERROR] Missing expected SVG output: $SVG_OUTPUT" >&2
    exit 1
fi

echo "[OK] Pipeline finished."
echo "[OK] PNG: $PNG_OUTPUT"
echo "[OK] SVG: $SVG_OUTPUT"
