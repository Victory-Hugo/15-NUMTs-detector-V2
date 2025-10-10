#! /bin/bash
################################################################################
## This script detects NUMTs from whole genome sequencing BAM files
## samtools, samblaster and blat need to be installed to run the pipeline
## samtools can be downloaded at http://www.htslib.org/download/
## samblaster can be downloaded at https://github.com/GregoryFaust/samblaster
## blat can be downloaded at http://hgdownload.soe.ucsc.edu/admin/exe/
################################################################################
## load modules on HPC
# module load samblaster/0.1.24
# module load samtools 
# module load Sambamba/0.6.6
###############################################################################
# NUMTs detection pipeline  ——  One-shot version (BWA-MEM → SAM → PSL → breakpoints)
# ❶ 预处理（disc/split）  ❷ 聚类脚本  ❸ 合并 FASTA（header=PREFIX|readID）
# ❹ bwa-mem 全量比对      ❺ SAM→PSL   ❻ 按 PREFIX 拆 PSL + 断点脚本
###############################################################################
set -euo pipefail

## ===============================  配置  ==================================== ##
MITOCHONDRIAL_CHR=chrM  # 线粒体染色体名称（MT|chrM）
INPUT_BAM=/mnt/d/9_Reference/HG001.bam
OUTPUT_DIR=/mnt/f/OneDrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/1-获得NUMTs分布/output/

BASE_DIR='/mnt/f/OneDrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/1-获得NUMTs分布/' # 脚本基础目录
REF_GRCh38=$BASE_DIR/download/GRCh38_latest_genomic.fna
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
## 1. 提取 MT 相关 read 并生成 disc/split SAM
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
if [[ ! -f "${REF_GRCh38}.bwt" ]]; then
  echo "    building bwa index ..."
  bwa index "$REF_GRCh38"
fi
bwa mem -M -t${THREADS} "$REF_GRCh38" "$COMBINED_FASTA" > "$BWA_SAM"

###############################################################################
## 5. SAM → PSL
###############################################################################
echo ">>> step5  SAM → PSL"
python3 "$SAM2PSL_SCRIPT" "$BWA_SAM" "${REF_GRCh38}.fai" "$ALL_PSL"

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
        "chr$chr_clean" \
        "$start" \
        "$end" \
        "${OUTPUT_DIR}/${PREFIX}"
done
echo "=== PIPELINE FINISHED ==="