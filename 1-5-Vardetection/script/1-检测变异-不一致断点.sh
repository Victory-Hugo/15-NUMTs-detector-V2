#! /bin/bash

################################################################################
## 该脚本用于检测来自不一致和分裂读取的NUMT变异
## 需要安装CAP3、blat和clustalo
################################################################################

set -u -o pipefail

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
: "${CLUSTALO_THREADS:=1}"
: "${STEP_TIMEOUT_SEC:=0}"
: "${FAIL_LOG:=${LOGS_DIR}/fail.log}"

tool_exists() {
    local tool="$1"
    [ -x "${tool}" ] || command -v "${tool}" > /dev/null 2>&1
}

timestamp() {
    date '+%Y-%m-%d %H:%M:%S'
}

log_msg() {
    local level="$1"
    local sample_id="$2"
    local contig_id="$3"
    shift 3
    local prefix
    prefix="$(timestamp) [${level}] [sample:${sample_id}]"
    if [ -n "${contig_id}" ]; then
        prefix="${prefix} [contig:${contig_id}]"
    fi
    echo "${prefix} $*"
}

run_with_optional_timeout() {
    if [ "${STEP_TIMEOUT_SEC}" -gt 0 ] && command -v timeout > /dev/null 2>&1; then
        timeout --preserve-status "${STEP_TIMEOUT_SEC}" "$@"
    else
        "$@"
    fi
}

record_sample_status() {
    local lock_file="$1"
    local status_file="$2"
    local sample_id="$3"
    local message="$4"

    (
        flock -x 200
        if [ -f "${status_file}" ] && grep -Fxq "${sample_id}" "${status_file}"; then
            exit 0
        fi
        if [ -n "${message}" ]; then
            printf '%s\t%s\n' "${sample_id}" "${message}" >> "${status_file}"
        else
            printf '%s\n' "${sample_id}" >> "${status_file}"
        fi
    ) 200>"${lock_file}"
}

require_nonempty_file() {
    local path="$1"
    local description="$2"
    if [ ! -s "${path}" ]; then
        echo "error: ${description} missing or empty: ${path}" >&2
        return 1
    fi
    return 0
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
    local fail_log="${FAIL_LOG}"
    local fail_lock="${logs_dir}/fail.lock"
    local output_dir_cap3="${output_dir}/assembly"
    local outputdir_psl_human="${output_dir}/pslHuman"
    local outputdir_psl_chimp="${output_dir}/pslChimp"
    local outputdir_aln_human="${output_dir}/alnHuman"
    local outputdir_aln_humanchimp="${output_dir}/alnHumanChimp"

    local fasta_list="${output_dir}/fastaFiles.list"
    local numts_bed="${output_dir}/${sample_id}_numts.bed"

    mkdir -p "${output_dir}" "${output_dir_cap3}" "${outputdir_psl_human}" "${outputdir_psl_chimp}" "${outputdir_aln_human}" "${outputdir_aln_humanchimp}" "${logs_dir}"

    if [ -f "${success_log}" ] && grep -Fxq "${sample_id}" "${success_log}"; then
        log_msg "INFO" "${sample_id}" "" "skip already completed"
        return 0
    fi

    if [ ! -f "${fasta_list}" ]; then
        find "${input_dir}" -name "*.fasta" > "${fasta_list}"
    fi

    if [ ! -s "${fasta_list}" ]; then
        log_msg "ERROR" "${sample_id}" "" "fasta list missing or empty: ${fasta_list}"
        record_sample_status "${fail_lock}" "${fail_log}" "${sample_id}" "empty_fasta_list"
        return 1
    fi

    exec 3< "${fasta_list}"
    while IFS= read -r contigINPUT <&3; do
        if [ -z "${contigINPUT}" ]; then
            continue
        fi
        if [ ! -f "${contigINPUT}" ]; then
            log_msg "ERROR" "${sample_id}" "" "contig input not found: ${contigINPUT}"
            exec 3<&-
            record_sample_status "${fail_lock}" "${fail_log}" "${sample_id}" "missing_contig_input"
            return 1
        fi

        local sampleIndex
        local contigBasename
        sampleIndex="$(basename "${contigINPUT}" .fasta)"
        contigBasename="$(basename "${contigINPUT}")"
        log_msg "INFO" "${sample_id}" "${sampleIndex}" "processing ${contigINPUT}"

        ### assembly ###
        local cap3_out
        local assembled_fasta
        cap3_out="${output_dir_cap3}/${contigBasename}.cap"

        # Check sequence count to avoid CAP3 memory issues
        local seq_count
        seq_count=$(grep -c "^>" "${contigINPUT}" || echo "0")
        local max_seq_for_cap3="${MAX_SEQ_FOR_CAP3:-500}"

        if [ "${seq_count}" -gt "${max_seq_for_cap3}" ]; then
            log_msg "INFO" "${sample_id}" "${sampleIndex}" "skip CAP3: ${seq_count} sequences exceed threshold ${max_seq_for_cap3}"
            assembled_fasta="${contigINPUT}"
        elif tool_exists "${CAP3_BIN}"; then
            if (cd "${output_dir_cap3}" && run_with_optional_timeout "${CAP3_BIN}" "${contigINPUT}" > "${cap3_out}" < /dev/null); then
                if [ -s "${cap3_out}.contigs" ]; then
                    assembled_fasta="${cap3_out}.contigs"
                else
                    log_msg "WARNING" "${sample_id}" "${sampleIndex}" "CAP3 completed but contigs output missing; fallback to original fasta"
                    assembled_fasta="${contigINPUT}"
                fi
            else
                log_msg "WARNING" "${sample_id}" "${sampleIndex}" "CAP3 failed; fallback to original fasta"
                assembled_fasta="${contigINPUT}"
            fi
        else
            log_msg "WARNING" "${sample_id}" "${sampleIndex}" "CAP3 not found; fallback to original fasta"
            assembled_fasta="${contigINPUT}"
        fi

        ### filter fasta to avoid empty/invalid sequences ###
        local filtered_fasta
        filtered_fasta="${output_dir}/${sampleIndex}.filtered.fasta"
        if ! run_with_optional_timeout "${PYTHON_BIN}" "${FILTER_FASTA_PY}" --input "${assembled_fasta}" --output "${filtered_fasta}" --min-acgt "${MIN_ACGT}" < /dev/null; then
            log_msg "ERROR" "${sample_id}" "${sampleIndex}" "filter_fasta_nonempty.py failed"
            exec 3<&-
            record_sample_status "${fail_lock}" "${fail_log}" "${sample_id}" "filter_fasta_failed:${sampleIndex}"
            return 1
        fi
        if [ ! -s "${filtered_fasta}" ]; then
            log_msg "WARNING" "${sample_id}" "${sampleIndex}" "filtered fasta empty; skip contig"
            continue
        fi

        ### reference alignment - psl file ###
        local psl_human
        local psl_chimp
        psl_human="${outputdir_psl_human}/${sampleIndex}.human.psl"
        psl_chimp="${outputdir_psl_chimp}/${sampleIndex}.chimp.psl"
        if ! tool_exists "${BLAT_BIN}"; then
            log_msg "ERROR" "${sample_id}" "${sampleIndex}" "blat not found"
            exec 3<&-
            record_sample_status "${fail_lock}" "${fail_log}" "${sample_id}" "missing_blat"
            return 1
        fi
        if ! run_with_optional_timeout "${BLAT_BIN}" "${REF_HUMAN}" "${filtered_fasta}" "${psl_human}" < /dev/null; then
            log_msg "ERROR" "${sample_id}" "${sampleIndex}" "blat failed for human reference"
            exec 3<&-
            record_sample_status "${fail_lock}" "${fail_log}" "${sample_id}" "blat_human_failed:${sampleIndex}"
            return 1
        fi
        if ! run_with_optional_timeout "${BLAT_BIN}" "${REF_CHIMP}" "${filtered_fasta}" "${psl_chimp}" < /dev/null; then
            log_msg "ERROR" "${sample_id}" "${sampleIndex}" "blat failed for chimp reference"
            exec 3<&-
            record_sample_status "${fail_lock}" "${fail_log}" "${sample_id}" "blat_chimp_failed:${sampleIndex}"
            return 1
        fi

        ### build numts.bed from psl ###
        if ! require_nonempty_file "${psl_human}" "human PSL file"; then
            log_msg "ERROR" "${sample_id}" "${sampleIndex}" "PSL file not found: ${psl_human}"
            exec 3<&-
            record_sample_status "${fail_lock}" "${fail_log}" "${sample_id}" "missing_psl:${sampleIndex}"
            return 1
        fi
        if ! run_with_optional_timeout "${PYTHON_BIN}" "${PSL_TO_BED_PY}" "${psl_human}" "${numts_bed}" < /dev/null; then
            log_msg "ERROR" "${sample_id}" "${sampleIndex}" "psl_to_numts_bed.py failed"
            exec 3<&-
            record_sample_status "${fail_lock}" "${fail_log}" "${sample_id}" "psl_to_bed_failed:${sampleIndex}"
            return 1
        fi
        if ! require_nonempty_file "${numts_bed}" "NUMTs BED"; then
            log_msg "ERROR" "${sample_id}" "${sampleIndex}" "NUMTs BED missing or empty"
            exec 3<&-
            record_sample_status "${fail_lock}" "${fail_log}" "${sample_id}" "missing_numts_bed:${sampleIndex}"
            return 1
        fi

        ### multiple alignment ###
        local humanMT_fasta
        local humanchimpMT_fasta
        humanMT_fasta="${output_dir}/${sampleIndex}.humanMT.fasta"
        humanchimpMT_fasta="${output_dir}/${sampleIndex}.humanchimpMT.fasta"

        if ! cat "${REF_HUMAN}" "${contigINPUT}" > "${humanMT_fasta}"; then
            log_msg "ERROR" "${sample_id}" "${sampleIndex}" "failed to build humanMT fasta"
            exec 3<&-
            record_sample_status "${fail_lock}" "${fail_log}" "${sample_id}" "build_human_fasta_failed:${sampleIndex}"
            return 1
        fi
        if ! cat "${REF_HUMANCHIMP}" "${contigINPUT}" > "${humanchimpMT_fasta}"; then
            log_msg "ERROR" "${sample_id}" "${sampleIndex}" "failed to build humanchimpMT fasta"
            exec 3<&-
            record_sample_status "${fail_lock}" "${fail_log}" "${sample_id}" "build_humanchimp_fasta_failed:${sampleIndex}"
            return 1
        fi

        local alnHuman
        local alnHumanChimp
        alnHuman="${outputdir_aln_human}/${sampleIndex}"
        alnHumanChimp="${outputdir_aln_humanchimp}/${sampleIndex}"
        if ! tool_exists "${CLUSTALO_BIN}"; then
            log_msg "ERROR" "${sample_id}" "${sampleIndex}" "clustalo not found"
            exec 3<&-
            record_sample_status "${fail_lock}" "${fail_log}" "${sample_id}" "missing_clustalo"
            return 1
        fi
        if ! run_with_optional_timeout "${CLUSTALO_BIN}" -i "${humanMT_fasta}" -o "${alnHuman}.humanMTaln.fa" --outfmt=fa --force --threads="${CLUSTALO_THREADS}" < /dev/null; then
            log_msg "ERROR" "${sample_id}" "${sampleIndex}" "clustalo failed for human alignment"
            exec 3<&-
            record_sample_status "${fail_lock}" "${fail_log}" "${sample_id}" "clustalo_human_failed:${sampleIndex}"
            return 1
        fi
        if ! require_nonempty_file "${alnHuman}.humanMTaln.fa" "human alignment"; then
            log_msg "ERROR" "${sample_id}" "${sampleIndex}" "human alignment missing or empty"
            exec 3<&-
            record_sample_status "${fail_lock}" "${fail_log}" "${sample_id}" "missing_human_alignment:${sampleIndex}"
            return 1
        fi
        if ! run_with_optional_timeout "${CLUSTALO_BIN}" -i "${humanchimpMT_fasta}" -o "${alnHumanChimp}.humanchimpMTaln.fa" --outfmt=fa --force --threads="${CLUSTALO_THREADS}" < /dev/null; then
            log_msg "ERROR" "${sample_id}" "${sampleIndex}" "clustalo failed for human+chimp alignment"
            exec 3<&-
            record_sample_status "${fail_lock}" "${fail_log}" "${sample_id}" "clustalo_humanchimp_failed:${sampleIndex}"
            return 1
        fi
        if ! require_nonempty_file "${alnHumanChimp}.humanchimpMTaln.fa" "human+chimp alignment"; then
            log_msg "ERROR" "${sample_id}" "${sampleIndex}" "human+chimp alignment missing or empty"
            exec 3<&-
            record_sample_status "${fail_lock}" "${fail_log}" "${sample_id}" "missing_humanchimp_alignment:${sampleIndex}"
            return 1
        fi

        ### variant detection ###
        if ! run_with_optional_timeout "${PYTHON_BIN}" "${PY_HUMAN}" "${alnHuman}.humanMTaln.fa" "${numts_bed}" < /dev/null; then
            log_msg "ERROR" "${sample_id}" "${sampleIndex}" "generateVariantTable.Human.py failed"
            exec 3<&-
            record_sample_status "${fail_lock}" "${fail_log}" "${sample_id}" "human_variant_failed:${sampleIndex}"
            return 1
        fi
        if ! require_nonempty_file "${alnHuman}.humanMTaln.fa.numt.tsv" "human numt output"; then
            log_msg "ERROR" "${sample_id}" "${sampleIndex}" "human numt output missing or empty"
            exec 3<&-
            record_sample_status "${fail_lock}" "${fail_log}" "${sample_id}" "missing_human_numt:${sampleIndex}"
            return 1
        fi
        if ! run_with_optional_timeout "${PYTHON_BIN}" "${PY_HUMANCHIMP}" "${alnHumanChimp}.humanchimpMTaln.fa" "${numts_bed}" < /dev/null; then
            log_msg "ERROR" "${sample_id}" "${sampleIndex}" "generateVariantTable.HumanChimp.py failed"
            exec 3<&-
            record_sample_status "${fail_lock}" "${fail_log}" "${sample_id}" "humanchimp_variant_failed:${sampleIndex}"
            return 1
        fi
        if ! require_nonempty_file "${alnHumanChimp}.humanchimpMTaln.fa.numtDhumanChimp.sum.tsv" "human+chimp summary output"; then
            log_msg "ERROR" "${sample_id}" "${sampleIndex}" "human+chimp summary output missing or empty"
            exec 3<&-
            record_sample_status "${fail_lock}" "${fail_log}" "${sample_id}" "missing_humanchimp_summary:${sampleIndex}"
            return 1
        fi
        log_msg "INFO" "${sample_id}" "${sampleIndex}" "contig completed"
    done
    exec 3<&-

    record_sample_status "${success_lock}" "${success_log}" "${sample_id}" ""
    log_msg "INFO" "${sample_id}" "" "sample completed"
    return 0
}

wait_for_active_pids() {
    local pid
    local status=0
    for pid in "$@"; do
        if ! wait "${pid}"; then
            status=1
        fi
    done
    return "${status}"
}

overall_status=0

if [ "${JOBS}" -le 1 ]; then
    exec 4< "${LIST_PATH}"
    while IFS=$'\t' read -r sample_id input_dir <&4; do
        if [ -z "${sample_id}" ] || [[ "${sample_id}" =~ ^# ]]; then
            continue
        fi
        if ! process_sample "${sample_id}" "${input_dir}"; then
            overall_status=1
        fi
    done
    exec 4<&-
else
    pids=()
    all_pids=()
    exec 4< "${LIST_PATH}"
    while IFS=$'\t' read -r sample_id input_dir <&4; do
        if [ -z "${sample_id}" ] || [[ "${sample_id}" =~ ^# ]]; then
            continue
        fi
        process_sample "${sample_id}" "${input_dir}" &
        pids+=("$!")
        all_pids+=("$!")
        if [ "${#pids[@]}" -ge "${JOBS}" ]; then
            if ! wait_for_active_pids "${pids[@]}"; then
                overall_status=1
            fi
            pids=()
        fi
    done
    exec 4<&-
    if [ "${#pids[@]}" -gt 0 ]; then
        if ! wait_for_active_pids "${pids[@]}"; then
            overall_status=1
        fi
    fi
fi

exit "${overall_status}"
