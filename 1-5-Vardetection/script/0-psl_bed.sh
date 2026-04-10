#! /bin/bash

set -u -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONF_FILE="${1:-${SCRIPT_DIR}/../conf/0-psl_bed.conf}"

if [ ! -f "${CONF_FILE}" ]; then
    echo "error: config not found: ${CONF_FILE}" >&2
    exit 1
fi

# shellcheck disable=SC1090
. "${CONF_FILE}"

: "${LIST_PATH:?error: LIST_PATH is not set}"
: "${PSL_TO_BED_PY:?error: PSL_TO_BED_PY is not set}"
: "${PYTHON_BIN:=python3}"
: "${PSL_SUFFIX:=_all_regions.psl}"
: "${PSL_HUMAN_SUBDIR:=pslHuman}"
: "${PSL_HUMAN_SUFFIX:=_all_regions.human.psl}"
: "${BED_SUFFIX:=_numts.bed}"

if [ ! -f "${LIST_PATH}" ]; then
    echo "error: list_path not found: ${LIST_PATH}" >&2
    exit 1
fi

if [ ! -f "${PSL_TO_BED_PY}" ]; then
    echo "error: psl_to_bed script not found: ${PSL_TO_BED_PY}" >&2
    exit 1
fi

process_sample() {
    local sample_id="$1"
    local sample_dir="$2"
    local psl_file
    local bed_file
    local alt_psl_file

    if [ -z "${sample_id}" ] || [ -z "${sample_dir}" ]; then
        echo "warning: skip empty sample entry" >&2
        return 0
    fi

    if [ ! -d "${sample_dir}" ]; then
        echo "error: sample dir not found for ${sample_id}: ${sample_dir}" >&2
        return 1
    fi

    psl_file="${sample_dir}/${sample_id}${PSL_SUFFIX}"
    alt_psl_file="${sample_dir}/${PSL_HUMAN_SUBDIR}/${sample_id}${PSL_HUMAN_SUFFIX}"
    bed_file="${sample_dir}/${sample_id}${BED_SUFFIX}"

    if [ -s "${psl_file}" ]; then
        :
    elif [ -s "${alt_psl_file}" ]; then
        psl_file="${alt_psl_file}"
    else
        echo "error: psl missing or empty for ${sample_id}: ${psl_file} or ${alt_psl_file}" >&2
        return 1
    fi

    if ! "${PYTHON_BIN}" "${PSL_TO_BED_PY}" "${psl_file}" "${bed_file}"; then
        echo "error: failed to build bed for ${sample_id}" >&2
        return 1
    fi

    if [ ! -s "${bed_file}" ]; then
        echo "error: generated bed missing or empty for ${sample_id}: ${bed_file}" >&2
        return 1
    fi

    echo "ok: ${bed_file}"
    return 0
}

main() {
    local has_error=0
    local sample_id
    local sample_dir

    while IFS=$'\t' read -r sample_id sample_dir _; do
        if [ -z "${sample_id}${sample_dir}" ]; then
            continue
        fi

        if ! process_sample "${sample_id}" "${sample_dir}"; then
            has_error=1
        fi
    done < "${LIST_PATH}"

    return "${has_error}"
}

main "$@"
