#! /bin/bash
#*##############################################################################
#* 本脚本用于从全基因组测序 BAM 文件中检测 NUMTs
#* 运行该流程需要安装 samtools、samblaster 和 bwa
#* samtools 可在 http://www.htslib.org/download/ 下载
#* samblaster 可在 https://github.com/GregoryFaust/samblaster 下载
#*##############################################################################

#*##############################################################################
#* NUMTs 检测流程 —— 一步式版本（BWA-MEM → SAM → PSL → 断点）
#* ❶ 预处理（disc/split）  ❷ 聚类脚本  ❸ 合并 FASTA（header=PREFIX|readID）
#* ❹ bwa-mem 全量比对      ❺ SAM→PSL   ❻ 按 PREFIX 拆 PSL + 断点脚本
#*##############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
GENERAL_CONF="/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/1-获得NUMTs分布/confs/0_NUMTs_detection_2.0.conf"
if [[ ! -f "$GENERAL_CONF" ]]; then
  echo "[ERROR] 未找到配置文件: $GENERAL_CONF" >&2
  exit 1
fi
source "$GENERAL_CONF"


## ===============================  配置  ==================================== ##
: "${INPUT_BAM:?missing INPUT_BAM in general.conf}"
: "${OUTPUT_DIR:?missing OUTPUT_DIR in general.conf}"
: "${BASE_DIR:?missing BASE_DIR in general.conf}"
: "${CONFIG_DIR:?missing CONFIG_DIR in general.conf}"
: "${DEFAULT_CONFIG:?missing DEFAULT_CONFIG in general.conf}"
: "${REF_GRCh:?missing REF_GRCh in general.conf}"
: "${CLUSTER_SCRIPT:?missing CLUSTER_SCRIPT in general.conf}"
: "${SAM2PSL_SCRIPT:?missing SAM2PSL_SCRIPT in general.conf}"
: "${BREAKPOINT_SCRIPT:?missing BREAKPOINT_SCRIPT in general.conf}"
: "${THREADS:?missing THREADS in general.conf}"
: "${BWA_MEM_THREADS:?missing BWA_MEM_THREADS in general.conf}"
: "${LOG_SUBDIR:?missing LOG_SUBDIR in general.conf}"
: "${OLD_LOG_NAME:?missing OLD_LOG_NAME in general.conf}"
###############################################################################

SAMPLE_ID=${INPUT_BAM##*/}; SAMPLE_ID=${SAMPLE_ID%.bam}

mkdir -p "$OUTPUT_DIR"

DISC_SAM=${OUTPUT_DIR}/${SAMPLE_ID}.mt.disc.sam
SPLIT_SAM=${OUTPUT_DIR}/${SAMPLE_ID}.mt.split.sam
BREAK_IN=${DISC_SAM}.breakpointINPUT.tsv          # 聚类脚本会生成
CLUSTER_OLD=${DISC_SAM}.cluster.old.tsv
COMBINED_FASTA=${OUTPUT_DIR}/${SAMPLE_ID}_all_regions.fasta
BWA_SAM=${OUTPUT_DIR}/${SAMPLE_ID}_all_regions.bwa.sam
ALL_PSL=${OUTPUT_DIR}/${SAMPLE_ID}_all_regions.psl
LOG_DIR="${OUTPUT_DIR}/${LOG_SUBDIR}"
OLD_LOG="${LOG_DIR}/${OLD_LOG_NAME}"
mkdir -p "$LOG_DIR"
###############################################################################
## 0. 决定使用的配置文件
###############################################################################
export NUMTS_CONFIG="${NUMTS_CONFIG:-$DEFAULT_CONFIG}"
echo ">>> step0  请注意你使用了: $NUMTS_CONFIG"
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

log_old_status() {
  local ITEM=$1
  local STATUS=$2
  local TIMESTAMP
  TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')
  echo "[${TIMESTAMP}] ${SAMPLE_ID} - ${ITEM} - ${STATUS}" >> "$OLD_LOG"
}
###############################################################################
## 1. 提取 线粒体 相关 read 并生成 disc/split SAM
###############################################################################
if [[ ! -s "$DISC_SAM" ]]; then
  echo ">>> step1  生成 disc/split SAM"
  samtools view -@${THREADS} -m"${SAMTOOLS_VIEW_MEM}" -h -F2 "$INPUT_BAM" \
    | grep -e '^@' -e "${MITOCHONDRIAL_CHR}" \
    | samtools sort -@${THREADS} -m"${SAMTOOLS_SORT_MEM}" -n - \
    | samtools view -h - \
    | samblaster --ignoreUnmated -e \
        -d "$DISC_SAM" \
        -s "$SPLIT_SAM" \
        -o /dev/null
fi

###############################################################################
## 2. 运行聚类脚本，生成 breakpointINPUT.tsv
###############################################################################
if [[ ! -s "$BREAK_IN" || ! -s "$CLUSTER_OLD" ]]; then
  echo ">>> step2  运行聚类脚本"
  python3 "$CLUSTER_SCRIPT" "$DISC_SAM" "$SAMPLE_ID" "$INPUT_BAM"
fi
if [[ -s "$CLUSTER_OLD" ]]; then
  log_old_status "cluster.old.tsv" "exists"
else
  log_old_status "cluster.old.tsv" "missing"
fi

###############################################################################
## 3. 合并所有 region 的序列到一个 FASTA（header=PREFIX|readID）
###############################################################################
if [[ ! -s "$COMBINED_FASTA" ]]; then
  echo ">>> step3  生成合并 FASTA"
  > "$COMBINED_FASTA"
  grep -v '^$' "$BREAK_IN" | while IFS=$'	' read -r sampleID cluster_no \
          disFile splitFile wgsBAM chr start end; do
      chr_clean=$(echo "$chr" | sed "s/[(),' ]//g")
      PREFIX="${sampleID}_${chr_clean}.${start}.${end}"

      samtools view "$wgsBAM" "${chr_clean}:${start}-${end}" \
        | awk '$6 !~ /150M|149M|148M|149S|148S/' \
        | cut -f1,10 \
        | while IFS=$'	' read -r readID seq; do
              printf ">%s|%s\n%s\n" "$PREFIX" "$readID" "$seq" >> "$COMBINED_FASTA"
          done
  done
fi

###############################################################################
## 4. bwa-mem 全量比对
###############################################################################
if [[ ! -s "$BWA_SAM" ]]; then
  echo ">>> step4  bwa-mem 全量比对"
  if [[ ! -f "${REF_GRCh}.bwt" ]]; then
    echo "    building bwa index ..."
    bwa index "$REF_GRCh"
  fi
  bwa mem -M -t${BWA_MEM_THREADS} "$REF_GRCh" "$COMBINED_FASTA" > "$BWA_SAM"
fi

###############################################################################
## 5. SAM → PSL
###############################################################################
if [[ ! -s "$ALL_PSL" ]]; then
  echo ">>> step5  SAM → PSL"
  python3 "$SAM2PSL_SCRIPT" "$BWA_SAM" "${REF_GRCh}.fai" "$ALL_PSL"
fi

###############################################################################
## 6. 按 PREFIX 拆 PSL 并调用断点脚本
###############################################################################
echo ">>> step6  拆分 PSL 并逐条调用断点脚本"
PSL_HEADER=$(head -n5 "$ALL_PSL")

grep -v '^$' "$BREAK_IN" | while IFS=$'\t' read -r sampleID cluster_no \
        disFile splitFile wgsBAM chr start end; do
    chr_clean=$(echo "$chr" | sed "s/[(),' ]//g")
    PREFIX="${sampleID}_${chr_clean}.${start}.${end}"
    PSL_PART="${OUTPUT_DIR}/${PREFIX}.psl"
    LEGACY_BP="${OUTPUT_DIR}/${PREFIX}.Breakpoints.old.tsv"

    if [[ -s "$LEGACY_BP" ]]; then
        log_old_status "${PREFIX}.Breakpoints.old.tsv" "exists"
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
        log_old_status "${PREFIX}.Breakpoints.old.tsv" "no_psl_match"
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
        log_old_status "${PREFIX}.Breakpoints.old.tsv" "generated"
    else
        log_old_status "${PREFIX}.Breakpoints.old.tsv" "missing_after_run"
    fi
done
echo "=== PIPELINE FINISHED ==="
