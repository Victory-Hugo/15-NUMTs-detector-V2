#! /bin/bash

################################################################################
## This script detects NUMT variants from discordance and split reads
## CAP3, blat and clustalo need to be installed
################################################################################

set -euo pipefail


CONF_FILE="/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/1-5-Vardetection/conf/VarDetection_fromDiscSplitReads.conf"

if [ ! -f "${CONF_FILE}" ]; then
    echo "error: config not found: ${CONF_FILE}"
    exit 1
fi

# shellcheck disable=SC1090
. "${CONF_FILE}"

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
        if command -v cap3 &> /dev/null; then
            (cd "${output_dir_cap3}" && cap3 "${contigINPUT}" > "${cap3_out}")
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
        python3 "${FILTER_FASTA_PY}" --input "${assembled_fasta}" --output "${filtered_fasta}" --min-acgt "${MIN_ACGT}"
        if [ ! -s "${filtered_fasta}" ]; then
            echo "warning: filtered fasta empty, skip ${contigINPUT}"
            continue
        fi

        ### reference alignment - psl file ###
        local psl_human
        local psl_chimp
        psl_human="${outputdir_psl_human}/${sampleIndex}.human.psl"
        psl_chimp="${outputdir_psl_chimp}/${sampleIndex}.chimp.psl"
        if command -v blat &> /dev/null; then
            blat "${REF_HUMAN}" "${filtered_fasta}" "${psl_human}"
            blat "${REF_CHIMP}" "${filtered_fasta}" "${psl_chimp}"
        else
            echo "warning: blat not found, skip PSL generation"
        fi

        ### build numts.bed from psl ###
        if [ -f "${psl_human}" ]; then
            python3 "${PSL_TO_BED_PY}" "${psl_human}" "${numts_bed}"
        else
            echo "error: PSL file not found: ${psl_human}"
            return 1
        fi

        ### multiple alignment ###
        local humanMT_fasta
        local humanchimpMT_fasta
        humanMT_fasta="${output_dir}/${sampleIndex}.humanMT.fasta"
        humanchimpMT_fasta="${output_dir}/${sampleIndex}.humanchimpMT.fasta"

        cat "${REF_HUMAN}" "${filtered_fasta}" > "${humanMT_fasta}"
        cat "${REF_HUMANCHIMP}" "${filtered_fasta}" > "${humanchimpMT_fasta}"

        local alnHuman
        local alnHumanChimp
        alnHuman="${outputdir_aln_human}/${sampleIndex}"
        alnHumanChimp="${outputdir_aln_humanchimp}/${sampleIndex}"
        if command -v clustalo &> /dev/null; then
            clustalo -i "${humanMT_fasta}" -o "${alnHuman}.humanMTaln.fa" --outfmt=fa --force
            clustalo -i "${humanchimpMT_fasta}" -o "${alnHumanChimp}.humanchimpMTaln.fa" --outfmt=fa --force
        else
            echo "warning: clustalo not found, skip alignment"
            continue
        fi

        ### variant detection ###
        python3 "${PY_HUMAN}" "${alnHuman}.humanMTaln.fa" "${numts_bed}"
        python3 "${PY_HUMANCHIMP}" "${alnHumanChimp}.humanchimpMTaln.fa" "${numts_bed}"
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
            wait -n
            new_pids=()
            for pid in "${pids[@]}"; do
                if kill -0 "${pid}" 2>/dev/null; then
                    new_pids+=("${pid}")
                fi
            done
            pids=("${new_pids[@]}")
        fi
    done < "${LIST_PATH}"
    wait
fi
