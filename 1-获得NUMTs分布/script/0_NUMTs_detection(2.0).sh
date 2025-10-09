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

INPUT_BAM=/home/luolintao/S20_NUMTs-detection-1.0/data/Y21100000492606.deduped.bam
OUTPUT_DIR=/home/luolintao/S20_NUMTs-detection-1.0/output

BASE_DIR='/home/luolintao/S20_NUMTs-detection-1.0/' # 脚本基础目录
REF_GRCh38=$BASE_DIR/download/GRCh38_latest_genomic.fna
CLUSTER_SCRIPT=$BASE_DIR/script/0_1_找聚类.py
SAM2PSL_SCRIPT=$BASE_DIR/script/0_2_sam2psl.py  # 前面提供的脚本
BREAKPOINT_SCRIPT=$BASE_DIR/script/0_2_找断点.py

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
  samtools view -@${THREADS} -m4G -h -F2 "$INPUT_BAM" \
    | grep -e '^@' -e 'MT' \
    | samtools sort -@${THREADS} -m4G -n - \
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











































# # 使用制表符作为分隔符读取文件的每一行，并将其字段分配给相应的变量
# while IFS=$'\t' read -r sampleID cluster_no disFile splitFile wgsBAM chr start end; do
#     # 去除 chr 字段中的多余字符（如括号和引号）
#     chr=$(echo $chr | sed "s/[(),' ]//g")
#     # 输出当前行的各个字段，方便调试
#     echo "$sampleID $cluster_no $disFile $splitFile $wgsBAM $chr $start $end"    
#     # 构建区域字符串，用于 samtools 提取指定区域
#     REGION="${chr}:${start}-${end}"    
#     # 构建输出文件的前缀路径
#     OUTPUT="${OUTPUT_DIR}/${sampleID}_${chr}.${start}.${end}"    
#     #使用 samtools 从 BAM 文件中提取指定区域的数据，并保存为 SAM 格式
#     samtools view "${wgsBAM}" "${REGION}" > "${OUTPUT}.sam"    
#     # 使用 awk 过滤掉特定 CIGAR 值的行，并提取第1和第10列，将结果保存为 FASTA 格式
#     awk '$6 !~ /150M|149M|148M|149S|148S/' "${OUTPUT}.sam" | cut -f1,10 > "${OUTPUT}.fasta"
#     # 使用 perl 脚本在每行前添加 '>' 符号，符合 FASTA 格式要求
#     perl -pi -e 's/^/>/g' "${OUTPUT}.fasta"    
#     # 使用 perl 脚本将制表符替换为换行符，符合 FASTA 格式要求
#     perl -pi -e 's/\t/\n/g' "${OUTPUT}.fasta"    
#     # 使用 BLAT 工具将 FASTA 文件比对到参考基因组，结果保存为 PSL 格式
#     /home/luolintao/S21_blat-35.1/bin/blat "${REF_GRCh38}" "${OUTPUT}.fasta" "${OUTPUT}.psl"    
#     # 运行 Python 脚本处理 BLAT 生成的 PSL 文件，提取断点信息
#     python3 "${BREAKPOINT_SCRIPT}" "${OUTPUT}.psl" "${sampleID}" "${chr}" "${start}" "${end}" "${OUTPUT}"
#     # 删除临时生成的 FASTA 文件
#     rm "${OUTPUT}.fasta"    
#     # 删除临时生成的 SAM 文件
#     rm "${OUTPUT}.sam"
# done <<< "$filelines"