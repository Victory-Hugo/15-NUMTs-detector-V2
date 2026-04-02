#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${1:?missing project root}"
CONFIG_PATH="${2:?missing config path}"

source "$PROJECT_ROOT/script/load_config.sh" "$CONFIG_PATH"

required_files=(
    "$PROJECT_ROOT/$CFG_PATHS_INPUT_FILE"
    "$PROJECT_ROOT/$CFG_PATHS_CIRCOS_TEMPLATE"
    "$PROJECT_ROOT/$CFG_PATHS_REFERENCE_DIR/karyotype.human.hg38.txt"
    "$PROJECT_ROOT/$CFG_PATHS_REFERENCE_DIR/MitoGenes.txt"
    "$PROJECT_ROOT/$CFG_PATHS_REFERENCE_DIR/MitoGenes_Dloop.txt"
    "$PROJECT_ROOT/$CFG_PATHS_REFERENCE_DIR/MitoGenes_codingGenes.txt"
    "$PROJECT_ROOT/$CFG_PATHS_REFERENCE_DIR/MitoGenes_rRNA.txt"
    "$PROJECT_ROOT/$CFG_PATHS_REFERENCE_DIR/MitoGenes_tRNA.txt"
    "$PROJECT_ROOT/$CFG_PATHS_REFERENCE_DIR/MitoGenes_color.txt"
    "$PROJECT_ROOT/conf/confs/custom_colors.conf"
    "$PROJECT_ROOT/conf/confs/ideogram_NUMTs.conf"
    "$PROJECT_ROOT/conf/confs/karyotype_NUMTs.conf"
    "$PROJECT_ROOT/conf/confs/links_allBreakpoints.soma.conf"
    "$PROJECT_ROOT/conf/confs/ticks_NUMTs.conf"
)

for path in "${required_files[@]}"; do
    if [[ ! -f "$path" ]]; then
        echo "[ERROR] Required file not found: $path" >&2
        exit 1
    fi
done

if [[ ! -x "$CFG_TOOLS_CONDA_BIN" ]]; then
    echo "[ERROR] Conda executable not found or not executable: $CFG_TOOLS_CONDA_BIN" >&2
    exit 1
fi

if [[ ! -d "$CFG_TOOLS_CIRCOS_ETC_DIR" ]]; then
    echo "[ERROR] Circos etc directory not found: $CFG_TOOLS_CIRCOS_ETC_DIR" >&2
    exit 1
fi

echo "[OK] Environment check passed."
