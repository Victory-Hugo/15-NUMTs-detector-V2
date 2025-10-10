#!/bin/bash

# 从 missing_ids.list 获取所有 BAM 路径
mapfile -t BAM_FILES < /public/home/heguanglin/09.luolintao/NUMTs-detection-1.0/script_Gansu/missing_ids.list
TOTAL_FILES=${#BAM_FILES[@]}

# 分成 200 个子脚本
NUM_SCRIPTS=200
# 每个脚本至少分到的文件数（向下取整）
FILES_PER_SCRIPT=$(( TOTAL_FILES / NUM_SCRIPTS ))
# 余数文件，前 REMAINDER 个脚本各多分 1 个
REMAINDER=$(( TOTAL_FILES % NUM_SCRIPTS ))

for (( i=0; i<NUM_SCRIPTS; i++ )); do
    # 计算本脚本分到的文件区间
    # start = i * FILES_PER_SCRIPT + min(i, REMAINDER)
    START=$(( i * FILES_PER_SCRIPT + ( i < REMAINDER ? i : REMAINDER ) ))
    # count = FILES_PER_SCRIPT + (i < REMAINDER ? 1 : 0)
    COUNT=$(( FILES_PER_SCRIPT + ( i < REMAINDER ? 1 : 0 ) ))
    # 取出子集
    ASSIGNED=( "${BAM_FILES[@]:$START:$COUNT}" )

    # 生成子脚本
    SCRIPT="script_${i}.sh"
    cat > "$SCRIPT" << 'EOF'
#!/bin/bash
#SBATCH -p normal
#SBATCH -J NUMT_$i
#SBATCH --time=2000:00:00
#SBATCH -n 1
#SBATCH -o NUMT_out_$i.log
#SBATCH -e NUMT_err_$i.log

set -euo pipefail

# 本脚本处理的 BAM 文件列表
BAM_FILES=(
EOF
    # 把 ASSIGNED 数组写进文件
    for bam in "${ASSIGNED[@]}"; do
        echo "    $bam" >> "$SCRIPT"
    done

    # 剩下的流水线模板
    cat >> "$SCRIPT" << 'EOF'
)

# 每个 BAM 依次跑
for INPUT_BAM in "${BAM_FILES[@]}"; do
    SAMPLE_ID=$(basename "$INPUT_BAM" .bam)
    OUTPUT_DIR="/home/luolintao/S20_NUMTs-detection-1.0/output/${SAMPLE_ID}"
    mkdir -p "$OUTPUT_DIR"

    BASE_DIR='/home/luolintao/S20_NUMTs-detection-1.0/'
    REF_GRCh38=$BASE_DIR/download/GRCh38_latest_genomic.fna
    CLUSTER_SCRIPT=$BASE_DIR/script/0_1_找聚类.py
    SAM2PSL_SCRIPT=$BASE_DIR/script/0_2_sam2psl.py
    BREAKPOINT_SCRIPT=$BASE_DIR/script/0_3_找断点.py
    THREADS=16

    DISC_SAM=${OUTPUT_DIR}/${SAMPLE_ID}.mt.disc.sam
    SPLIT_SAM=${OUTPUT_DIR}/${SAMPLE_ID}.mt.split.sam
    BREAK_IN=${DISC_SAM}.breakpointINPUT.tsv
    COMBINED_FASTA=${OUTPUT_DIR}/${SAMPLE_ID}_all_regions.fasta
    BWA_SAM=${OUTPUT_DIR}/${SAMPLE_ID}_all_regions.bwa.sam
    ALL_PSL=${OUTPUT_DIR}/${SAMPLE_ID}_all_regions.psl

    echo ">>> [$SAMPLE_ID] Step 1: extract MT reads"
    if [[ ! -s "$DISC_SAM" ]]; then
      samtools view -@${THREADS} -m4G -h -F2 "$INPUT_BAM" \
        | grep -e '^@' -e 'MT' \
        | samtools sort -@${THREADS} -m4G -n - \
        | samtools view -h - \
        | samblaster --ignoreUnmated -e \
            -d "$DISC_SAM" \
            -s "$SPLIT_SAM" \
            -o /dev/null
    fi

    echo ">>> [$SAMPLE_ID] Step 2: cluster breakpoint"
    if [[ ! -s "$BREAK_IN" ]]; then
      python3 "$CLUSTER_SCRIPT" "$DISC_SAM" "$SAMPLE_ID" "$INPUT_BAM"
    fi

    echo ">>> [$SAMPLE_ID] Step 3: merge FASTA"
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

    echo ">>> [$SAMPLE_ID] Step 4: bwa-mem"
    if [[ ! -f "${REF_GRCh38}.bwt" ]]; then
      bwa index "$REF_GRCh38"
    fi
    bwa mem -M -t${THREADS} "$REF_GRCh38" "$COMBINED_FASTA" > "$BWA_SAM"

    echo ">>> [$SAMPLE_ID] Step 5: SAM→PSL"
    python3 "$SAM2PSL_SCRIPT" "$BWA_SAM" "${REF_GRCh38}.fai" "$ALL_PSL"

    echo ">>> [$SAMPLE_ID] Step 6: split PSL & call breakpoint"
    PSL_HEADER=$(head -n5 "$ALL_PSL")
    grep -v '^$' "$BREAK_IN" | while IFS=$'\t' read -r sampleID cluster_no \
            disFile splitFile wgsBAM chr start end; do
        chr_clean=$(echo "$chr" | sed "s/[(),' ]//g")
        PREFIX="${sampleID}_${chr_clean}.${start}.${end}"
        PSL_PART="${OUTPUT_DIR}/${PREFIX}.psl"

        { printf "%s\n" "$PSL_HEADER"
          awk -F'\t' -v q="$PREFIX" '${
              split($10,a,"|");
              if(a[1]==q) print
          }' "$ALL_PSL"
        } > "$PSL_PART"

        if [[ $(wc -l < "$PSL_PART") -le 5 ]]; then
            echo "  [WARN] $PREFIX no alignments, skip"
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

    echo "=== [$SAMPLE_ID] PIPELINE FINISHED ==="
done

EOF

    # 确保子脚本可执行
    chmod +x "$SCRIPT"
done
