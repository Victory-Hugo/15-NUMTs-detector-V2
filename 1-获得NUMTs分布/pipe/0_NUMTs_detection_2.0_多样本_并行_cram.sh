#! /bin/bash
#*##############################################################################
#* 本脚本用于从全基因组测序 CRAM 文件中检测 NUMTs（多样本并行）
#* 优先使用 GNU parallel；若不可用则退回 xargs 并行
#*##############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
GENERAL_CONF="/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/1-获得NUMTs分布/confs/0_NUMTs_detection_2.0_多样本_并行.conf"
RUN_MODE=${1:-""}
WORKER_INPUT_CRAM=${2:-""}
WORKER_SAMPLE_ID=${3:-""}
if [[ ! -f "$GENERAL_CONF" ]]; then
  echo "[ERROR] 未找到配置文件: $GENERAL_CONF" >&2
  exit 1
fi
source "$GENERAL_CONF"

# Ensure a valid working directory even if the caller's cwd was deleted.
if ! pwd >/dev/null 2>&1; then
  if [[ -n "${BASE_DIR:-}" && -d "${BASE_DIR}" ]]; then
    cd "$BASE_DIR"
  else
    cd /
  fi
fi

## ===============================  配置  ==================================== ##
: "${SAMPLE_LIST:?missing SAMPLE_LIST in general.conf}"
: "${OUTPUT_ROOT:?missing OUTPUT_ROOT in general.conf}"
: "${BASE_DIR:?missing BASE_DIR in general.conf}"
: "${CONFIG_DIR:?missing CONFIG_DIR in general.conf}"
: "${DEFAULT_CONFIG:?missing DEFAULT_CONFIG in general.conf}"
: "${REF_GRCh:?missing REF_GRCh in general.conf}"
: "${CLUSTER_SCRIPT:?missing CLUSTER_SCRIPT in general.conf}"
: "${SAM2PSL_SCRIPT:?missing SAM2PSL_SCRIPT in general.conf}"
: "${BREAKPOINT_SCRIPT:?missing BREAKPOINT_SCRIPT in general.conf}"
: "${THREADS:?missing THREADS in general.conf}"
: "${BWA_MEM_THREADS:?missing BWA_MEM_THREADS in general.conf}"
: "${SAMTOOLS_VIEW_MEM:?missing SAMTOOLS_VIEW_MEM in general.conf}"
: "${SAMTOOLS_SORT_MEM:?missing SAMTOOLS_SORT_MEM in general.conf}"
: "${PARALLEL_JOBS:?missing PARALLEL_JOBS in general.conf}"
: "${LOG_SUBDIR:?missing LOG_SUBDIR in general.conf}"
: "${OLD_LOG_NAME:?missing OLD_LOG_NAME in general.conf}"
###############################################################################

# CRAM 特有配置：需要指定参考基因组（必须）
: "${REF_CRAM:?missing REF_CRAM in general.conf (required for CRAM files)}"

export NUMTS_CONFIG="${NUMTS_CONFIG:-$DEFAULT_CONFIG}"
echo ">>> 使用配置: $NUMTS_CONFIG"
echo ">>> 样本列表: $SAMPLE_LIST"
echo ">>> 输出根目录: $OUTPUT_ROOT"
echo ">>> 并行数量: $PARALLEL_JOBS"
echo ">>> CRAM 参考基因组: $REF_CRAM"

if [[ "$LOG_SUBDIR" = /* ]]; then
  LOG_DIR="$LOG_SUBDIR"
else
  LOG_DIR="${OUTPUT_ROOT}/${LOG_SUBDIR}"
fi
mkdir -p "$LOG_DIR"

if [[ "$OLD_LOG_NAME" = /* ]]; then
  OLD_LOG="$OLD_LOG_NAME"
else
  OLD_LOG="${LOG_DIR}/${OLD_LOG_NAME}"
fi

SUCCESS_LOG="${LOG_DIR}/success.log"
FAILURE_LOG="${LOG_DIR}/failure.log"
LOCK_FILE="${LOG_DIR}/.log.lock"

log_success() {
  local SAMPLE_ID=$1
  local TIMESTAMP
  TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')
  (
    flock -x 200
    echo -e "${TIMESTAMP}\t${SAMPLE_ID}\tSUCCESS" >> "$SUCCESS_LOG"
  ) 200>"$LOCK_FILE"
}

log_failure() {
  local SAMPLE_ID=$1
  local STEP=$2
  local ERROR_MSG=$3
  local TIMESTAMP
  TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')
  (
    flock -x 200
    echo -e "${TIMESTAMP}\t${SAMPLE_ID}\tFAILED\t${STEP}\t${ERROR_MSG}" >> "$FAILURE_LOG"
  ) 200>"$LOCK_FILE"
}

log_old_status() {
  local ITEM=$1
  local STATUS=$2
  local SAMPLE_ID=$3
  local TIMESTAMP
  TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')
  (
    flock -x 200
    echo "[${TIMESTAMP}] ${SAMPLE_ID} - ${ITEM} - ${STATUS}" >> "$OLD_LOG"
  ) 200>"$LOCK_FILE"
}

# 从配置读取线粒体染色体名称
MITOCHONDRIAL_CHR=$(
  python3 - "$NUMTS_CONFIG" <<'PY'
import json
import sys

with open(sys.argv[1], "r") as f:
    cfg = json.load(f)
print(cfg.get("cluster", {}).get("mt_chrom", ""))
PY
)
if [[ -z "$MITOCHONDRIAL_CHR" ]]; then
  echo "[ERROR] 无法从配置读取 mt_chrom: $NUMTS_CONFIG" >&2
  exit 1
fi

process_sample() {
  local INPUT_CRAM=$1
  local SAMPLE_ID=$2
  if [[ -z "$INPUT_CRAM" || -z "$SAMPLE_ID" ]]; then
    echo "[WARN] 空样本行，跳过"
    return 0
  fi
  local OUTPUT_DIR="${OUTPUT_ROOT}/${SAMPLE_ID}"
  local CURRENT_STEP=""

  handle_error() {
    local exit_code=$?
    log_failure "$SAMPLE_ID" "$CURRENT_STEP" "Exit code: ${exit_code}"
    echo "[ERROR] ${SAMPLE_ID} 在步骤 ${CURRENT_STEP} 失败，退出码: ${exit_code}" >&2
    return 1
  }

  trap 'handle_error' ERR

  mkdir -p "$OUTPUT_DIR"

  local DISC_SAM=${OUTPUT_DIR}/${SAMPLE_ID}.mt.disc.sam
  local SPLIT_SAM=${OUTPUT_DIR}/${SAMPLE_ID}.mt.split.sam
  local BREAK_IN=${DISC_SAM}.breakpointINPUT.tsv
  local CLUSTER_OLD=${DISC_SAM}.cluster.old.tsv
  local COMBINED_FASTA=${OUTPUT_DIR}/${SAMPLE_ID}_all_regions.fasta
  local BWA_SAM=${OUTPUT_DIR}/${SAMPLE_ID}_all_regions.bwa.sam
  local ALL_PSL=${OUTPUT_DIR}/${SAMPLE_ID}_all_regions.psl

  echo "=============================================="
  echo ">>> 开始处理样本: $SAMPLE_ID"
  echo ">>> CRAM文件: $INPUT_CRAM"
  echo ">>> 输出目录: $OUTPUT_DIR"
  echo "=============================================="

  ###########################################################################
  ## 1. 提取 线粒体 相关 read 并生成 disc/split SAM
  ###########################################################################
  CURRENT_STEP="step1_disc_split_SAM"
  if [[ ! -s "$DISC_SAM" ]]; then
    echo ">>> step1  生成 disc/split SAM"
    samtools view -@${THREADS} -m"${SAMTOOLS_VIEW_MEM}" -h -F2 -T "$REF_CRAM" "$INPUT_CRAM" \
      | grep -e '^@' -e "${MITOCHONDRIAL_CHR}" \
      | samtools sort -@${THREADS} -m"${SAMTOOLS_SORT_MEM}" -n - \
      | samtools view -h - \
      | samblaster --ignoreUnmated -e \
          -d "$DISC_SAM" \
          -s "$SPLIT_SAM" \
          -o /dev/null
  fi

  ###########################################################################
  ## 2. 运行聚类脚本，生成 breakpointINPUT.tsv
  ###########################################################################
  CURRENT_STEP="step2_clustering"
  if [[ ! -s "$BREAK_IN" || ! -s "$CLUSTER_OLD" ]]; then
    echo ">>> step2  运行聚类脚本"
    python3 "$CLUSTER_SCRIPT" "$DISC_SAM" "$SAMPLE_ID" "$INPUT_CRAM"
  fi
  if [[ -s "$CLUSTER_OLD" ]]; then
    log_old_status "cluster.old.tsv" "exists" "$SAMPLE_ID"
  else
    log_old_status "cluster.old.tsv" "missing" "$SAMPLE_ID"
  fi

  ###########################################################################
  ## 3. 合并所有 region 的序列到一个 FASTA
  ###########################################################################
  CURRENT_STEP="step3_generate_FASTA"
  if [[ ! -s "$COMBINED_FASTA" ]]; then
    echo ">>> step3  生成合并 FASTA"
    > "$COMBINED_FASTA"
    grep -v '^$' "$BREAK_IN" | while IFS=$'\t' read -r sampleID cluster_no \
            disFile splitFile wgsFile chr start end; do
        chr_clean=$(echo "$chr" | sed "s/[(),' ]//g")
        PREFIX="${sampleID}_${chr_clean}.${start}.${end}"

        samtools view -T "$REF_CRAM" "$wgsFile" "${chr_clean}:${start}-${end}" \
          | awk '$6 !~ /150M|149M|148M|149S|148S/' \
          | cut -f1,10 \
          | while IFS=$'\t' read -r readID seq; do
                printf ">%s|%s\n%s\n" "$PREFIX" "$readID" "$seq" >> "$COMBINED_FASTA"
            done
    done
  fi

  if [[ ! -s "$COMBINED_FASTA" ]]; then
    echo "[WARN] ${SAMPLE_ID}: 没有候选区域，跳过后续步骤"
    log_success "$SAMPLE_ID"
    trap - ERR
    return 0
  fi

  ###########################################################################
  ## 4. bwa-mem 全量比对
  ###########################################################################
  CURRENT_STEP="step4_bwa_mem"
  if [[ ! -s "$BWA_SAM" ]]; then
    echo ">>> step4  bwa-mem 全量比对"
    if [[ ! -f "${REF_GRCh}.bwt" ]]; then
      echo "    building bwa index ..."
      bwa index "$REF_GRCh"
    fi
    bwa mem -M -t${BWA_MEM_THREADS} "$REF_GRCh" "$COMBINED_FASTA" > "$BWA_SAM"
  fi

  ###########################################################################
  ## 5. SAM → PSL
  ###########################################################################
  CURRENT_STEP="step5_SAM_to_PSL"
  if [[ ! -s "$ALL_PSL" ]]; then
    echo ">>> step5  SAM → PSL"
    python3 "$SAM2PSL_SCRIPT" "$BWA_SAM" "${REF_GRCh}.fai" "$ALL_PSL"
  fi

  ###########################################################################
  ## 6. 按 PREFIX 拆 PSL 并调用断点脚本
  ###########################################################################
  CURRENT_STEP="step6_breakpoint_detection"
  echo ">>> step6  拆分 PSL 并逐条调用断点脚本"
  local PSL_HEADER
  PSL_HEADER=$(head -n5 "$ALL_PSL")

  grep -v '^$' "$BREAK_IN" | while IFS=$'\t' read -r sampleID cluster_no \
          disFile splitFile wgsFile chr start end; do
      chr_clean=$(echo "$chr" | sed "s/[(),' ]//g")
      PREFIX="${sampleID}_${chr_clean}.${start}.${end}"
      PSL_PART="${OUTPUT_DIR}/${PREFIX}.psl"
      LEGACY_BP="${OUTPUT_DIR}/${PREFIX}.Breakpoints.old.tsv"

      if [[ -s "$LEGACY_BP" ]]; then
          log_old_status "${PREFIX}.Breakpoints.old.tsv" "exists" "$SAMPLE_ID"
          continue
      fi

      {   printf "%s\n" "$PSL_HEADER"
          awk -F'\t' -v q="$PREFIX" '
              {
                split($10, a, "|")
                if (a[1]==q) print
              }' "$ALL_PSL"
      } > "$PSL_PART"

      if [[ $(wc -l < "$PSL_PART") -le 5 ]]; then
          echo "  [WARN] $PREFIX 没在 PSL 找到比对，跳过"
          log_old_status "${PREFIX}.Breakpoints.old.tsv" "no_psl_match" "$SAMPLE_ID"
          continue
      fi
      python3 "$BREAKPOINT_SCRIPT" \
          "$PSL_PART" \
          "$sampleID" \
          "$chr_clean" \
          "$start" \
          "$end" \
          "${OUTPUT_DIR}/${PREFIX}"
      if [[ -s "$LEGACY_BP" ]]; then
          log_old_status "${PREFIX}.Breakpoints.old.tsv" "generated" "$SAMPLE_ID"
      else
          log_old_status "${PREFIX}.Breakpoints.old.tsv" "missing_after_run" "$SAMPLE_ID"
      fi
  done

  log_success "$SAMPLE_ID"
  echo ">>> 样本 ${SAMPLE_ID} 处理完成"
  trap - ERR
}

if [[ "$RUN_MODE" == "--worker" ]]; then
  if [[ -z "$WORKER_INPUT_CRAM" ]]; then
    echo "[ERROR] worker 缺少 CRAM 路径" >&2
    exit 1
  fi
  if [[ -z "$WORKER_SAMPLE_ID" ]]; then
    WORKER_SAMPLE_ID=${WORKER_INPUT_CRAM##*/}
    WORKER_SAMPLE_ID=${WORKER_SAMPLE_ID%.cram}
  fi
  process_sample "$WORKER_INPUT_CRAM" "$WORKER_SAMPLE_ID"
  exit $?
fi

PENDING_LIST=$(mktemp)
trap 'rm -f "$PENDING_LIST"' EXIT

if [[ ! -s "$SAMPLE_LIST" ]]; then
  echo "[ERROR] 样本列表为空: $SAMPLE_LIST" >&2
  exit 1
fi

if [[ -s "$SUCCESS_LOG" ]]; then
  awk -F'\t' '
    NR==FNR {
      if ($3=="SUCCESS" && $2!="") done[$2]=1
      next
    }
    $0!~/^#/ && $0!~/^$/ {
      line=$0
      sub(/\r$/, "", line)
      n=split(line, a, "\t")
      cram=a[1]; sid=a[2]
      sub(/\r/, "", cram)
      sub(/\r/, "", sid)
      if (cram=="") next
      if (sid=="") {
        sid=cram
        sub(/^.*\//, "", sid)
        sub(/\.cram$/, "", sid)
      }
      if (!done[sid]) print cram "\t" sid
    }
  ' "$SUCCESS_LOG" "$SAMPLE_LIST" > "$PENDING_LIST"
else
  awk -F'\t' '
    $0!~/^#/ && $0!~/^$/ {
      line=$0
      sub(/\r$/, "", line)
      n=split(line, a, "\t")
      cram=a[1]; sid=a[2]
      sub(/\r/, "", cram)
      sub(/\r/, "", sid)
      if (cram=="") next
      if (sid=="") {
        sid=cram
        sub(/^.*\//, "", sid)
        sub(/\.cram$/, "", sid)
      }
      print cram "\t" sid
    }
  ' "$SAMPLE_LIST" > "$PENDING_LIST"
fi

PENDING_COUNT=$(wc -l < "$PENDING_LIST" | tr -d ' ')
echo ">>> 待处理样本数: $PENDING_COUNT"
if [[ "$PENDING_COUNT" -eq 0 ]]; then
  echo ">>> 没有待处理样本"
  exit 0
fi

if command -v parallel >/dev/null 2>&1; then
  echo ">>> 使用 GNU parallel 并行"
  parallel --no-run-if-empty --colsep '\t' -j "${PARALLEL_JOBS}" \
    "${SCRIPT_DIR}/0_NUMTs_detection_2.0_多样本_并行_cram.sh" --worker {1} {2} \
    :::: "$PENDING_LIST"
else
  echo ">>> 未找到 GNU parallel，使用 xargs 并行"
  xargs -P "${PARALLEL_JOBS}" -n 2 -a "$PENDING_LIST" \
    "${SCRIPT_DIR}/0_NUMTs_detection_2.0_多样本_并行_cram.sh" --worker
fi

echo "=== ALL DONE ==="
