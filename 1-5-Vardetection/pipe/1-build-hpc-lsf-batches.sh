#!/bin/bash

set -euo pipefail

# ==============================================================================
# 下面这一段是本脚本使用的固定路径配置区。
# 如果你把整个项目复制到别的目录，优先修改这里，不要在脚本中间到处找路径。
# ==============================================================================

# 项目根目录
PIPELINE_ROOT="/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/1-5-Vardetection"

# pipe 脚本目录
SCRIPT_DIR="${PIPELINE_ROOT}/pipe"

# 默认使用的主配置文件
DEFAULT_BASE_CONF="${PIPELINE_ROOT}/conf/VarDetection_fromDiscSplitReads.conf"

# 实际执行分析流程的主脚本
RUN_SCRIPT="${PIPELINE_ROOT}/script/1-检测变异-不一致断点.sh"

# 默认的批次输出目录
DEFAULT_BATCH_ROOT="${PIPELINE_ROOT}/conf/hpc_batches"

# ==============================================================================
# 下面这一段是运行参数默认值，不是输入文件路径。
# ==============================================================================

DEFAULT_LINES_PER_CHUNK=1
DEFAULT_JOB_PREFIX="Vardetect"
DEFAULT_QUEUE="normal"
DEFAULT_WALL_TIME="2000:00"

BASE_CONF=""
LINES_PER_CHUNK="${DEFAULT_LINES_PER_CHUNK}"
BATCH_ID=""
JOB_PREFIX="${DEFAULT_JOB_PREFIX}"
QUEUE="${DEFAULT_QUEUE}"
WALL_TIME="${DEFAULT_WALL_TIME}"
BATCH_ROOT_BASE="${DEFAULT_BATCH_ROOT}"
GENERATED_DIRS=()

usage() {
    cat <<EOF
Usage:
  bash pipe/1-build-hpc-lsf-batches.sh --base-conf conf/VarDetection_fromDiscSplitReads.conf [options]

Options:
  --base-conf PATH         Base config file. Default: ${DEFAULT_BASE_CONF}
  --lines-per-chunk N      Number of samples per split file. Default: ${DEFAULT_LINES_PER_CHUNK}
  --batch-id ID            Batch directory name. Default: batch_YYYYmmdd_HHMMSS
  --job-prefix PREFIX      LSF job name prefix. Default: ${DEFAULT_JOB_PREFIX}
  --queue QUEUE            LSF queue name. Default: ${DEFAULT_QUEUE}
  --wall-time HH:MM        LSF wall time. Default: ${DEFAULT_WALL_TIME}
  --batch-root PATH        Directory to store generated batches. Default: ${DEFAULT_BATCH_ROOT}
  --help                   Show this help message
EOF
}

cleanup_on_error() {
    local exit_code=$?
    if [ "${exit_code}" -ne 0 ] && [ "${#GENERATED_DIRS[@]}" -gt 0 ]; then
        for dir in "${GENERATED_DIRS[@]}"; do
            if [ -d "${dir}" ]; then
                rm -rf "${dir}"
            fi
        done
    fi
    exit "${exit_code}"
}

trap cleanup_on_error EXIT

abs_path() {
    local path="$1"
    if [ -d "${path}" ]; then
        (cd "${path}" && pwd)
    else
        local dir
        local base
        dir="$(dirname "${path}")"
        base="$(basename "${path}")"
        (cd "${dir}" && printf '%s/%s\n' "$(pwd)" "${base}")
    fi
}

while [ "$#" -gt 0 ]; do
    case "$1" in
        --base-conf)
            BASE_CONF="$2"
            shift 2
            ;;
        --lines-per-chunk)
            LINES_PER_CHUNK="$2"
            shift 2
            ;;
        --batch-id)
            BATCH_ID="$2"
            shift 2
            ;;
        --job-prefix)
            JOB_PREFIX="$2"
            shift 2
            ;;
        --queue)
            QUEUE="$2"
            shift 2
            ;;
        --wall-time)
            WALL_TIME="$2"
            shift 2
            ;;
        --batch-root)
            BATCH_ROOT_BASE="$2"
            shift 2
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            echo "error: unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [ -z "${BASE_CONF}" ]; then
    BASE_CONF="${DEFAULT_BASE_CONF}"
fi

if [ ! -f "${BASE_CONF}" ]; then
    echo "error: base config not found: ${BASE_CONF}" >&2
    exit 1
fi

BASE_CONF="$(abs_path "${BASE_CONF}")"

if ! [[ "${LINES_PER_CHUNK}" =~ ^[0-9]+$ ]] || [ "${LINES_PER_CHUNK}" -lt 1 ]; then
    echo "error: --lines-per-chunk must be an integer >= 1" >&2
    exit 1
fi

if [ ! -f "${RUN_SCRIPT}" ]; then
    echo "error: run script not found: ${RUN_SCRIPT}" >&2
    exit 1
fi

# shellcheck disable=SC1090
. "${BASE_CONF}"

if [ -z "${LIST_PATH:-}" ]; then
    echo "error: LIST_PATH is not defined in ${BASE_CONF}" >&2
    exit 1
fi

LIST_PATH="$(abs_path "${LIST_PATH}")"

if [ ! -f "${LIST_PATH}" ]; then
    echo "error: list path not found: ${LIST_PATH}" >&2
    exit 1
fi

if [ -z "${BATCH_ID}" ]; then
    BATCH_ID="batch_$(date '+%Y%m%d_%H%M%S')"
fi

mkdir -p "${BATCH_ROOT_BASE}"
BATCH_ROOT_BASE="$(abs_path "${BATCH_ROOT_BASE}")"
BATCH_ROOT="${BATCH_ROOT_BASE}/${BATCH_ID}"
SPLIT_DIR="${BATCH_ROOT}/split"
CONF_DIR="${BATCH_ROOT}/conf"
LSF_DIR="${BATCH_ROOT}/lsf"
LOG_DIR="${BATCH_ROOT}/logs"
README_FILE="${BATCH_ROOT}/README.submit.txt"
CLEAN_MASTER_LIST="${BATCH_ROOT}/cleaned_master_list.txt"

if [ -e "${BATCH_ROOT}" ]; then
    echo "error: batch directory already exists: ${BATCH_ROOT}" >&2
    exit 1
fi

mkdir -p "${SPLIT_DIR}" "${CONF_DIR}" "${LSF_DIR}" "${LOG_DIR}"
GENERATED_DIRS+=("${BATCH_ROOT}")

awk '
BEGIN {
    valid = 0
}
{
    line = $0
    gsub(/\r$/, "", line)
    if (line ~ /^[[:space:]]*$/) {
        next
    }
    if (line ~ /^[[:space:]]*#/) {
        next
    }
    n = split(line, fields, "\t")
    if (n < 2 || fields[1] == "" || fields[2] == "") {
        printf("error: malformed sample list at line %d: %s\n", NR, line) > "/dev/stderr"
        exit 1
    }
    print line
    valid++
}
END {
    if (valid == 0) {
        print "error: no valid sample records found after filtering LIST_PATH" > "/dev/stderr"
        exit 2
    }
}
' "${LIST_PATH}" > "${CLEAN_MASTER_LIST}"

split -l "${LINES_PER_CHUNK}" -d -a 4 --additional-suffix=.txt \
    "${CLEAN_MASTER_LIST}" "${SPLIT_DIR}/list_"

mapfile -t SPLIT_FILES < <(find "${SPLIT_DIR}" -maxdepth 1 -type f -name 'list_*.txt' | sort)

if [ "${#SPLIT_FILES[@]}" -eq 0 ]; then
    echo "error: split did not generate any files" >&2
    exit 1
fi

task_count=0

for split_file in "${SPLIT_FILES[@]}"; do
    split_file="$(abs_path "${split_file}")"
    split_name="$(basename "${split_file}" .txt)"
    suffix="${split_name#list_}"
    conf_name="VarDetection_fromDiscSplitReads_${suffix}.conf"
    lsf_name="${JOB_PREFIX}_${suffix}.lsf"
    job_name="${JOB_PREFIX}_${suffix}"
    conf_path="${CONF_DIR}/${conf_name}"
    lsf_path="${LSF_DIR}/${lsf_name}"

    awk -v list_path="${split_file}" '
    BEGIN {
        list_updated = 0
        jobs_updated = 0
    }
    /^LIST_PATH=/ {
        print "LIST_PATH=\"" list_path "\""
        list_updated = 1
        next
    }
    /^JOBS=/ {
        print "JOBS=1"
        jobs_updated = 1
        next
    }
    {
        print
    }
    END {
        if (!list_updated) {
            print "LIST_PATH=\"" list_path "\""
        }
        if (!jobs_updated) {
            print "JOBS=1"
        }
    }
    ' "${BASE_CONF}" > "${conf_path}"

    cat > "${lsf_path}" <<EOF
#!/bin/bash
#BSUB -J ${job_name}
#BSUB -q ${QUEUE}
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -W ${WALL_TIME}
#BSUB -o ${LOG_DIR}/${job_name}.out
#BSUB -e ${LOG_DIR}/${job_name}.err

set -euo pipefail

PIPELINE_ROOT="${PIPELINE_ROOT}"
CONF_FILE="${conf_path}"
RUN_SCRIPT="${RUN_SCRIPT}"

bash "\${RUN_SCRIPT}" "\${CONF_FILE}"
EOF

    task_count=$((task_count + 1))
done

cat > "${README_FILE}" <<EOF
Batch ID: ${BATCH_ID}
Generated at: $(date '+%Y-%m-%d %H:%M:%S %z')
Base conf: ${BASE_CONF}
Original list: ${LIST_PATH}
Batches root: ${BATCH_ROOT_BASE}
Lines per chunk: ${LINES_PER_CHUNK}
Task count: ${task_count}

Files:
- split/: per-task sample lists
- conf/: per-task config files with LIST_PATH overridden and JOBS=1
- lsf/: manual submission scripts
- logs/: LSF stdout/stderr targets

Manual submission:
cd "${LSF_DIR}"
for f in *.lsf; do bsub < "\$f"; done

Single job example:
bsub < ${JOB_PREFIX}_0000.lsf
EOF

echo "batch generated: ${BATCH_ROOT}"
echo "task count: ${task_count}"
echo "lsf dir: ${LSF_DIR}"

trap - EXIT
