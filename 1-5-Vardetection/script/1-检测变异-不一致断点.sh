#! /bin/bash

################################################################################
## 该脚本用于检测来自不一致和分裂读取的NUMT变异
## 需要安装CAP3、blat和clustalo
################################################################################

# set -euo pipefail # 注释掉 set -e 来允许单个样本失败而不影响整个批次的处理

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONF_FILE="${1:-${SCRIPT_DIR}/../conf/VarDetection_fromDiscSplitReads.conf}"

if [ ! -f "${CONF_FILE}" ]; then
    echo "error: config not found: ${CONF_FILE}"
    exit 1
fi

# shellcheck disable=SC1090
. "${CONF_FILE}"

: "${PYTHON_BIN:=python3}"
: "${CAP3_BIN:=cap3}"
: "${BLAT_BIN:=blat}"
: "${CLUSTALO_BIN:=clustalo}"

tool_exists() {
    local tool="$1"
    [ -x "${tool}" ] || command -v "${tool}" > /dev/null 2>&1
}

if [ ! -f "${LIST_PATH}" ]; then
    echo "error: list_path not found: ${LIST_PATH}"
    exit 1
fi

process_sample() {
    local sample_id="$1"
    local input_dir="$2"

    if [ -z "${sample_id}" ] || [ -z "${input_dir}" ]; then
        echo "warning: skip empty sample entry"
        return 0
    fi

    if [ ! -d "${input_dir}" ]; then
        echo "error: input dir not found: ${input_dir}"
        return 1
    fi

    local output_dir="${OUTPUT_BASE}/${sample_id}"
    local logs_dir="${LOGS_DIR}"
    local success_log="${logs_dir}/success.log"
    local success_lock="${logs_dir}/success.lock"
    local output_dir_cap3="${output_dir}/assembly"
    local outputdir_psl_human="${output_dir}/pslHuman"
    local outputdir_psl_chimp="${output_dir}/pslChimp"
    local outputdir_aln_human="${output_dir}/alnHuman"
    local outputdir_aln_humanchimp="${output_dir}/alnHumanChimp"

    local fasta_list="${output_dir}/fastaFiles.list"
    local numts_bed="${output_dir}/${sample_id}_numts.bed"

    mkdir -p "${output_dir}" "${output_dir_cap3}" "${outputdir_psl_human}" "${outputdir_psl_chimp}" "${outputdir_aln_human}" "${outputdir_aln_humanchimp}" "${logs_dir}"

    if [ -f "${success_log}" ] && grep -Fxq "${sample_id}" "${success_log}"; then
        echo "skip ${sample_id}: already completed"
        return 0
    fi

    if [ ! -f "${fasta_list}" ]; then
        find "${input_dir}" -name "*.fasta" > "${fasta_list}"
    fi

    while IFS= read -r contigINPUT; do
        if [ -z "${contigINPUT}" ]; then
            continue
        fi
        echo "processing ${contigINPUT} ..."

        local sampleIndex
        local contigBasename
        sampleIndex="$(basename "${contigINPUT}" .fasta)"
        contigBasename="$(basename "${contigINPUT}")"

        ### assembly ###
        local cap3_out
        local assembled_fasta
        cap3_out="${output_dir_cap3}/${contigBasename}.cap"

        # Check sequence count to avoid CAP3 memory issues
        local seq_count
        seq_count=$(grep -c "^>" "${contigINPUT}" || echo "0")
        local max_seq_for_cap3="${MAX_SEQ_FOR_CAP3:-500}"

        if [ "${seq_count}" -gt "${max_seq_for_cap3}" ]; then
            echo "INFO: Skipping CAP3 assembly for ${contigINPUT}"
            echo "  Reason: Input contains ${seq_count} sequences (threshold: ${max_seq_for_cap3})"
            echo "  CAP3 may cause memory overflow or segmentation fault with large inputs"
            echo "  Impact: Using original reads instead of assembled contigs"
            echo "  Note: This does NOT affect variant detection or age inference results"
            echo "  Note: Output may contain redundant variant calls from overlapping reads"
            assembled_fasta="${contigINPUT}"
        elif tool_exists "${CAP3_BIN}"; then
            (cd "${output_dir_cap3}" && "${CAP3_BIN}" "${contigINPUT}" > "${cap3_out}") || {
                echo "warning: CAP3 failed, fallback to original fasta"
                assembled_fasta="${contigINPUT}"
            }
            assembled_fasta="${cap3_out}.contigs"
            if [ ! -f "${assembled_fasta}" ]; then
                echo "warning: CAP3 output missing, fallback to original fasta"
                assembled_fasta="${contigINPUT}"
            fi
        else
            assembled_fasta="${contigINPUT}"
        fi

        ### filter fasta to avoid empty/invalid sequences ###
        local filtered_fasta
        filtered_fasta="${output_dir}/${sampleIndex}.filtered.fasta"
        "${PYTHON_BIN}" "${FILTER_FASTA_PY}" --input "${assembled_fasta}" --output "${filtered_fasta}" --min-acgt "${MIN_ACGT}"
        if [ ! -s "${filtered_fasta}" ]; then
            echo "warning: filtered fasta empty, skip ${contigINPUT}"
            continue
        fi

        ### reference alignment - psl file ###
        local psl_human
        local psl_chimp
        psl_human="${outputdir_psl_human}/${sampleIndex}.human.psl"
        psl_chimp="${outputdir_psl_chimp}/${sampleIndex}.chimp.psl"
        if tool_exists "${BLAT_BIN}"; then
            "${BLAT_BIN}" "${REF_HUMAN}" "${filtered_fasta}" "${psl_human}"
            "${BLAT_BIN}" "${REF_CHIMP}" "${filtered_fasta}" "${psl_chimp}"
        else
            echo "warning: blat not found, skip PSL generation"
        fi

        ### build numts.bed from psl ###
        if [ -f "${psl_human}" ]; then
            "${PYTHON_BIN}" "${PSL_TO_BED_PY}" "${psl_human}" "${numts_bed}"
        else
            echo "error: PSL file not found: ${psl_human}"
            return 1
        fi

        ### multiple alignment ###
        local humanMT_fasta
        local humanchimpMT_fasta
        humanMT_fasta="${output_dir}/${sampleIndex}.humanMT.fasta"
        humanchimpMT_fasta="${output_dir}/${sampleIndex}.humanchimpMT.fasta"

        cat "${REF_HUMAN}" "${contigINPUT}" > "${humanMT_fasta}"
        cat "${REF_HUMANCHIMP}" "${contigINPUT}" > "${humanchimpMT_fasta}"

        local alnHuman
        local alnHumanChimp
        alnHuman="${outputdir_aln_human}/${sampleIndex}"
        alnHumanChimp="${outputdir_aln_humanchimp}/${sampleIndex}"
        if tool_exists "${CLUSTALO_BIN}"; then
            "${CLUSTALO_BIN}" -i "${humanMT_fasta}" -o "${alnHuman}.humanMTaln.fa" --outfmt=fa --force
            "${CLUSTALO_BIN}" -i "${humanchimpMT_fasta}" -o "${alnHumanChimp}.humanchimpMTaln.fa" --outfmt=fa --force
        else
            echo "warning: clustalo not found, skip alignment"
            continue
        fi

        ### variant detection ###
        "${PYTHON_BIN}" "${PY_HUMAN}" "${alnHuman}.humanMTaln.fa" "${numts_bed}"
        "${PYTHON_BIN}" "${PY_HUMANCHIMP}" "${alnHumanChimp}.humanchimpMTaln.fa" "${numts_bed}"
    done < "${fasta_list}"

    (
        flock -x 200
        if [ -f "${success_log}" ] && grep -Fxq "${sample_id}" "${success_log}"; then
            exit 0
        fi
        echo "${sample_id}" >> "${success_log}"
    ) 200>"${success_lock}"
}

if [ "${JOBS}" -le 1 ]; then
    while IFS=$'\t' read -r sample_id input_dir; do
        if [ -z "${sample_id}" ] || [[ "${sample_id}" =~ ^# ]]; then
            continue
        fi
        process_sample "${sample_id}" "${input_dir}"
    done < "${LIST_PATH}"
else
    pids=()
    while IFS=$'\t' read -r sample_id input_dir; do
        if [ -z "${sample_id}" ] || [[ "${sample_id}" =~ ^# ]]; then
            continue
        fi
        process_sample "${sample_id}" "${input_dir}" &
        pids+=("$!")
        if [ "${#pids[@]}" -ge "${JOBS}" ]; then
            # wait -n requires Bash 4.3+; use a portable polling loop instead
            while true; do
                new_pids=()
                for pid in "${pids[@]}"; do
                    if kill -0 "${pid}" 2>/dev/null; then
                        new_pids+=("${pid}")
                    fi
                done
                pids=("${new_pids[@]}")
                if [ "${#pids[@]}" -lt "${JOBS}" ]; then
                    break
                fi
                sleep 1
            done
        fi
    done < "${LIST_PATH}"
    wait
fi
