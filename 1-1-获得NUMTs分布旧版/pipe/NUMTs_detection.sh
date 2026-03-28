#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONFIG_PATH="${PROJECT_ROOT}/conf/Config.yaml"

if [[ $# -ne 3 ]]; then
    echo "[ERROR] 用法: $0 <bam> <output_dir> <reference>" >&2
    exit 1
fi

INPUT_BAM="$1"
OUTPUT_ROOT="$2"
REFERENCE="$3"

if [[ ! -f "${INPUT_BAM}" ]]; then
    echo "[ERROR] BAM 文件不存在: ${INPUT_BAM}" >&2
    exit 1
fi

if [[ ! -f "${REFERENCE}" ]]; then
    echo "[ERROR] 参考基因组不存在: ${REFERENCE}" >&2
    exit 1
fi

source "${PROJECT_ROOT}/script/load_config.sh" "${CONFIG_PATH}"

CONDA_SH="${CFG_RUNTIME_CONDA_BASE}/etc/profile.d/conda.sh"
if [[ -f "${CONDA_SH}" && "${CONDA_DEFAULT_ENV:-}" != "${CFG_RUNTIME_CONDA_ENV}" ]]; then
    # shellcheck disable=SC1090
    source "${CONDA_SH}"
    conda activate "${CFG_RUNTIME_CONDA_ENV}"
fi

PYTHON_BIN="${CFG_TOOLS_PYTHON}"
SAMTOOLS_BIN="${CFG_TOOLS_SAMTOOLS}"
SAMBLASTER_BIN="${CFG_TOOLS_SAMBLASTER}"
BLAT_BIN="${CFG_TOOLS_BLAT}"
THREADS="${CFG_RUNTIME_THREADS}"
SORT_MEM="${CFG_RUNTIME_SORT_MEM}"

for required_cmd in "${PYTHON_BIN}" "${SAMTOOLS_BIN}" "${SAMBLASTER_BIN}" "${BLAT_BIN}"; do
    if ! command -v "${required_cmd}" >/dev/null 2>&1; then
        echo "[ERROR] 未找到命令: ${required_cmd}" >&2
        exit 1
    fi
done

SAMPLE_NAME="$(basename "${INPUT_BAM}")"
SAMPLE_ID="${SAMPLE_NAME%.bam}"
SAMPLE_OUTPUT_DIR="${OUTPUT_ROOT}/${SAMPLE_ID}"
LOG_DIR="${OUTPUT_ROOT}/logs"
mkdir -p "${SAMPLE_OUTPUT_DIR}" "${LOG_DIR}"

SAMPLE_LOG="${LOG_DIR}/${SAMPLE_ID}.log"
exec > >(tee -a "${SAMPLE_LOG}") 2>&1

CLUSTER_SCRIPT="${PROJECT_ROOT}/python/searchNumtCluster_fromDiscordantReads.py"
BREAKPOINT_SCRIPT="${PROJECT_ROOT}/python/searchBreakpoint_fromblatoutputs.py"

INPUT_DISC="${SAMPLE_OUTPUT_DIR}/${SAMPLE_ID}.mt.disc.sam"
INPUT_SPLIT="${SAMPLE_OUTPUT_DIR}/${SAMPLE_ID}.mt.split.sam"
BREAKPOINT_INPUT="${INPUT_DISC}.breakpointINPUT.tsv"

echo "[INFO] 开始处理样本: ${SAMPLE_ID}"
echo "[INFO] 输出目录: ${SAMPLE_OUTPUT_DIR}"

if [[ -s "${INPUT_DISC}" && -e "${INPUT_SPLIT}" ]]; then
    echo "[INFO] 已存在 disc/split 结果，跳过提取步骤"
else
    "${SAMTOOLS_BIN}" view -@"${THREADS}" -h -F 2 "${INPUT_BAM}" \
        | grep -e '^@' -e 'MT' -e 'chrM' \
        | "${SAMTOOLS_BIN}" sort -@"${THREADS}" -m "${SORT_MEM}" -n - \
        | "${SAMTOOLS_BIN}" view -h - \
        | "${SAMBLASTER_BIN}" --ignoreUnmated -e -d "${INPUT_DISC}" -s "${INPUT_SPLIT}" -o /dev/null
fi

"${PYTHON_BIN}" "${CLUSTER_SCRIPT}" \
    --sample-id "${SAMPLE_ID}" \
    --wgs-bam "${INPUT_BAM}" \
    --input-disc "${INPUT_DISC}"

echo "[INFO] 聚类步骤完成，开始搜索断点"

if [[ ! -s "${BREAKPOINT_INPUT}" ]]; then
    echo "[INFO] breakpointINPUT.tsv 为空，样本 ${SAMPLE_ID} 无候选区域"
    exit 0
fi

while IFS=$'\t' read -r sample_id cluster_no input_dis input_split input_wgs chr start end || [[ -n "${sample_id:-}" ]]; do
    [[ -z "${sample_id:-}" ]] && continue

    region="${chr}:${start}-${end}"
    output_prefix="${SAMPLE_OUTPUT_DIR}/${sample_id}_${chr}.${start}.${end}"

    "${SAMTOOLS_BIN}" view "${input_wgs}" "${region}" > "${output_prefix}.sam"

    awk '$6 !~ /150M|149M|148M|149S|148S/' "${output_prefix}.sam" | cut -f1,10 > "${output_prefix}.fasta"
    perl -pi -e 's/^/>/g' "${output_prefix}.fasta"
    perl -pi -e 's/\t/\n/g' "${output_prefix}.fasta"

    "${BLAT_BIN}" "${REFERENCE}" "${output_prefix}.fasta" "${output_prefix}.psl"
    "${PYTHON_BIN}" "${BREAKPOINT_SCRIPT}" \
        --input-psl "${output_prefix}.psl" \
        --sample-id "${sample_id}" \
        --chromosome "${chr}" \
        --start "${start}" \
        --end "${end}" \
        --output-prefix "${output_prefix}"

    rm -f "${output_prefix}.fasta" "${output_prefix}.sam"
done < "${BREAKPOINT_INPUT}"

echo "[INFO] 样本 ${SAMPLE_ID} 处理完成"
