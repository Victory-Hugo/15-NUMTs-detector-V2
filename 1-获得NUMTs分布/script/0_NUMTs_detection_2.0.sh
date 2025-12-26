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


## ===============================  配置  ==================================== ##
INPUT_BAM=/mnt/d/5-NCBI-Reference/3-Human/example/HG001.bam  #? 输入文件
OUTPUT_DIR=/mnt/d/5-NCBI-Reference/3-Human/example/output/HG001/                           #? 输出目录

BASE_DIR='/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/1-获得NUMTs分布/' # 脚本基础目录
CONFIG_DIR="${BASE_DIR}/confs"
DEFAULT_CONFIG="${CONFIG_DIR}/GRCH38.json"

REF_GRCh=/mnt/d/5-NCBI-Reference/3-Human/Reference/GRCh38_latest_genomic.fna
CLUSTER_SCRIPT=$BASE_DIR/script/0_1_找聚类.py
SAM2PSL_SCRIPT=$BASE_DIR/script/0_2_sam2psl.py  # 前面提供的脚本
BREAKPOINT_SCRIPT=$BASE_DIR/script/0_3_找断点.py

THREADS=16
###############################################################################

SAMPLE_ID=${INPUT_BAM##*/}; SAMPLE_ID=${SAMPLE_ID%.bam}

mkdir -p "$OUTPUT_DIR"

DISC_SAM=${OUTPUT_DIR}/${SAMPLE_ID}.mt.disc.sam
SPLIT_SAM=${OUTPUT_DIR}/${SAMPLE_ID}.mt.split.sam
BREAK_IN=${DISC_SAM}.breakpointINPUT.tsv          # 聚类脚本会生成
COMBINED_FASTA=${OUTPUT_DIR}/${SAMPLE_ID}_all_regions.fasta
BWA_SAM=${OUTPUT_DIR}/${SAMPLE_ID}_all_regions.bwa.sam
ALL_PSL=${OUTPUT_DIR}/${SAMPLE_ID}_all_regions.psl
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
###############################################################################
## 1. 提取 线粒体 相关 read 并生成 disc/split SAM
###############################################################################
if [[ ! -s "$DISC_SAM" ]]; then
  echo ">>> step1  生成 disc/split SAM"
  samtools view -@${THREADS} -m32G -h -F2 "$INPUT_BAM" \
    | grep -e '^@' -e "${MITOCHONDRIAL_CHR}" \
    | samtools sort -@${THREADS} -m32G -n - \
    | samtools view -h - \
    | samblaster --ignoreUnmated -e \
        -d "$DISC_SAM" \
        -s "$SPLIT_SAM" \
        -o /dev/null
fi

###############################################################################
## 2. 运行聚类脚本，生成 breakpointINPUT.tsv
###############################################################################
if [[ ! -s "$BREAK_IN" ]]; then
  echo ">>> step2  运行聚类脚本"
  python3 "$CLUSTER_SCRIPT" "$DISC_SAM" "$SAMPLE_ID" "$INPUT_BAM"
fi

###############################################################################
## 3. 合并所有 region 的序列到一个 FASTA（header=PREFIX|readID）
###############################################################################
echo ">>> step3  生成合并 FASTA"
> "$COMBINED_FASTA"
grep -v '^$' "$BREAK_IN" | while IFS=$'\t' read -r sampleID cluster_no \
        disFile splitFile wgsBAM chr start end; do
    chr_clean=$(echo "$chr" | sed "s/[(),' ]//g")
    PREFIX="${sampleID}_${chr_clean}.${start}.${end}"

    samtools view "$wgsBAM" "${chr_clean}:${start}-${end}" \
      | awk '$6 !~ /150M|149M|148M|149S|148S/' \
      | cut -f1,10 \
      | while IFS=$'\t' read -r readID seq; do
            printf ">%s|%s\n%s\n" "$PREFIX" "$readID" "$seq" >> "$COMBINED_FASTA"
        done
done

###############################################################################
## 4. bwa-mem 全量比对
###############################################################################
echo ">>> step4  bwa-mem 全量比对"
if [[ ! -f "${REF_GRCh}.bwt" ]]; then
  echo "    building bwa index ..."
  bwa index "$REF_GRCh"
fi
bwa mem -M -t${THREADS} "$REF_GRCh" "$COMBINED_FASTA" > "$BWA_SAM"

###############################################################################
## 5. SAM → PSL
###############################################################################
echo ">>> step5  SAM → PSL"
python3 "$SAM2PSL_SCRIPT" "$BWA_SAM" "${REF_GRCh}.fai" "$ALL_PSL"

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

    {   printf "%s\n" "$PSL_HEADER"
        awk -F'\t' -v q="$PREFIX" '
            {
              split($10, a, "|")
              if (a[1]==q) print
            }' "$ALL_PSL"
    } > "$PSL_PART"

    if [[ $(wc -l < "$PSL_PART") -le 5 ]]; then
        echo "  [WARN] $PREFIX 没在 PSL 找到比对，跳过"
        continue
    fi
    python3 "$BREAKPOINT_SCRIPT" \
        "$PSL_PART" \
        "$sampleID" \
        "$chr_clean" \
        "$start" \
        "$end" \
        "${OUTPUT_DIR}/${PREFIX}"
done
echo "=== PIPELINE FINISHED ==="
