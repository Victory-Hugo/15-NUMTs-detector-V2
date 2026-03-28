#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONFIG_PATH="${PROJECT_ROOT}/conf/Config.yaml"

source "${PROJECT_ROOT}/script/load_config.sh" "${CONFIG_PATH}"

CONDA_SH="${CFG_RUNTIME_CONDA_BASE}/etc/profile.d/conda.sh"
if [[ ! -f "${CONDA_SH}" ]]; then
    echo "[ERROR] 未找到 conda 初始化脚本: ${CONDA_SH}" >&2
    exit 1
fi

# shellcheck disable=SC1090
source "${CONDA_SH}"
conda activate "${CFG_RUNTIME_CONDA_ENV}"

mkdir -p "${CFG_PATHS_OUTPUT_ROOT}" "${CFG_PATHS_OUTPUT_ROOT}/logs" "${CFG_PATHS_OUTPUT_ROOT}/multi_sample"

"${PROJECT_ROOT}/script/check_env.sh"

PIPELINE_LOG="${CFG_PATHS_OUTPUT_ROOT}/logs/pipeline.log"
SUCCESS_LOG="${CFG_PATHS_OUTPUT_ROOT}/logs/success.log"
FAILURE_LOG="${CFG_PATHS_OUTPUT_ROOT}/logs/failure.log"
MERGED_BREAKPOINT_INPUT="${CFG_PATHS_OUTPUT_ROOT}/multi_sample/merged.breakpointINPUT.tsv"
GROUP_SCRIPT="${PROJECT_ROOT}/python/groupNumtCluster_fromMultipleSamples.py"
SAMPLE_WRAPPER="${PROJECT_ROOT}/pipe/NUMTs_detection.sh"

: > "${PIPELINE_LOG}"
: > "${SUCCESS_LOG}"
: > "${FAILURE_LOG}"
: > "${MERGED_BREAKPOINT_INPUT}"

log_message() {
    local level="$1"
    local message="$2"
    local timestamp
    timestamp="$(date '+%Y-%m-%d %H:%M:%S')"
    echo "[${timestamp}] [${level}] ${message}" | tee -a "${PIPELINE_LOG}"
}

sample_failures=0

while IFS= read -r bam_path || [[ -n "${bam_path}" ]]; do
    [[ -z "${bam_path}" ]] && continue
    sample_name="$(basename "${bam_path}")"
    sample_id="${sample_name%.bam}"

    log_message "INFO" "开始处理样本 ${sample_id}"
    if bash "${SAMPLE_WRAPPER}" "${bam_path}" "${CFG_PATHS_OUTPUT_ROOT}" "${CFG_PATHS_REFERENCE}"; then
        echo -e "$(date '+%Y-%m-%d %H:%M:%S')\t${sample_id}\tSUCCESS" >> "${SUCCESS_LOG}"
        log_message "INFO" "样本 ${sample_id} 处理成功"
    else
        echo -e "$(date '+%Y-%m-%d %H:%M:%S')\t${sample_id}\tFAILED" >> "${FAILURE_LOG}"
        log_message "ERROR" "样本 ${sample_id} 处理失败"
        sample_failures=$((sample_failures + 1))
    fi
done < "${CFG_PATHS_SAMPLE_LIST}"

if [[ "${sample_failures}" -ne 0 ]]; then
    log_message "ERROR" "共有 ${sample_failures} 个样本失败，停止多样本汇总"
    exit 1
fi

while IFS= read -r bam_path || [[ -n "${bam_path}" ]]; do
    [[ -z "${bam_path}" ]] && continue
    sample_name="$(basename "${bam_path}")"
    sample_id="${sample_name%.bam}"
    breakpoint_file="${CFG_PATHS_OUTPUT_ROOT}/${sample_id}/${sample_id}.mt.disc.sam.breakpointINPUT.tsv"
    if [[ -s "${breakpoint_file}" ]]; then
        cat "${breakpoint_file}" >> "${MERGED_BREAKPOINT_INPUT}"
    fi
done < "${CFG_PATHS_SAMPLE_LIST}"

if [[ ! -s "${MERGED_BREAKPOINT_INPUT}" ]]; then
    log_message "WARNING" "所有样本的 breakpointINPUT.tsv 均为空，跳过多样本聚类"
    exit 0
fi

log_message "INFO" "开始执行多样本聚类"
"${CFG_TOOLS_PYTHON}" "${GROUP_SCRIPT}" --input-file "${MERGED_BREAKPOINT_INPUT}"
log_message "INFO" "主流程完成"
