#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONFIG_PATH="${PROJECT_ROOT}/conf/Config.yaml"

source "${PROJECT_ROOT}/script/load_config.sh" "${CONFIG_PATH}"

CONDA_BIN="${CFG_RUNTIME_CONDA_BASE}/condabin/conda"

if [[ ! -x "${CONDA_BIN}" ]]; then
    echo "[ERROR] 未找到 conda: ${CONDA_BIN}" >&2
    exit 1
fi

if ! "${CONDA_BIN}" env list | awk '{print $1}' | grep -Fxq "${CFG_RUNTIME_CONDA_ENV}"; then
    echo "[ERROR] 未找到 conda 环境: ${CFG_RUNTIME_CONDA_ENV}" >&2
    exit 1
fi

if [[ ! -f "${CFG_PATHS_REFERENCE}" ]]; then
    echo "[ERROR] 未找到参考基因组: ${CFG_PATHS_REFERENCE}" >&2
    exit 1
fi

if [[ ! -f "${CFG_PATHS_SAMPLE_LIST}" ]]; then
    echo "[ERROR] 未找到样本列表: ${CFG_PATHS_SAMPLE_LIST}" >&2
    exit 1
fi

for tool_path in "${CFG_TOOLS_PYTHON}" "${CFG_TOOLS_SAMTOOLS}" "${CFG_TOOLS_SAMBLASTER}" "${CFG_TOOLS_BLAT}"; do
    if ! "${CONDA_BIN}" run -n "${CFG_RUNTIME_CONDA_ENV}" bash -lc "command -v '${tool_path}' >/dev/null"; then
        echo "[ERROR] 环境 ${CFG_RUNTIME_CONDA_ENV} 中缺少工具: ${tool_path}" >&2
        exit 1
    fi
done

while IFS= read -r bam_path || [[ -n "${bam_path}" ]]; do
    [[ -z "${bam_path}" ]] && continue
    if [[ ! -f "${bam_path}" ]]; then
        echo "[ERROR] 未找到 BAM 文件: ${bam_path}" >&2
        exit 1
    fi
    if [[ ! -f "${bam_path}.bai" ]]; then
        echo "[ERROR] 未找到 BAM 索引: ${bam_path}.bai" >&2
        exit 1
    fi
done < "${CFG_PATHS_SAMPLE_LIST}"

echo "[INFO] 环境检查通过"
