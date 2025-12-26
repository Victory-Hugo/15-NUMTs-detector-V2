#!/bin/bash
#*##############################################################################
#* 服务器版本 - 批量 NUMTs 检测流程
#* 读取样本列表，将任务分配到 10 个子脚本中并行处理
#* 用法: bash 服务器版本.sh [样本列表文件] [输出根目录] [子脚本数量]
#*##############################################################################

set -euo pipefail

## =============================== 参数解析 ================================== ##
SAMPLE_LIST=${1:-"/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/1-获得NUMTs分布/confs/input_bam.path.txt"}
OUTPUT_ROOT=${2:-"/mnt/d/5-NCBI-Reference/3-Human/example/output"}
NUM_SCRIPTS=${3:-10}

## =============================== 基础配置 ================================== ##
BASE_DIR='/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/1-获得NUMTs分布/'
CONFIG_DIR="${BASE_DIR}/confs"
DEFAULT_CONFIG="${CONFIG_DIR}/GRCH38.json"  # 默认使用 GRCh38

REF_GRCh=/mnt/d/5-NCBI-Reference/3-Human/Reference/GRCh38_latest_genomic.fna
CLUSTER_SCRIPT=$BASE_DIR/script/0_1_找聚类.py
SAM2PSL_SCRIPT=$BASE_DIR/script/0_2_sam2psl.py
BREAKPOINT_SCRIPT=$BASE_DIR/script/0_3_找断点.py

THREADS=16

## =============================== 脚本生成目录 ============================== ##
SCRIPT_OUTPUT_DIR="${BASE_DIR}/pipe/generated_scripts"
mkdir -p "$SCRIPT_OUTPUT_DIR"

## =============================== 读取样本列表 ============================== ##
mapfile -t SAMPLES < <(grep -v '^#' "$SAMPLE_LIST" | grep -v '^$')
TOTAL_SAMPLES=${#SAMPLES[@]}

if [[ $TOTAL_SAMPLES -eq 0 ]]; then
    echo "[ERROR] 样本列表为空: $SAMPLE_LIST"
    exit 1
fi

echo "=========================================="
echo "总样本数: $TOTAL_SAMPLES"
echo "子脚本数: $NUM_SCRIPTS"
echo "输出根目录: $OUTPUT_ROOT"
echo "使用配置: $DEFAULT_CONFIG"
echo "=========================================="

## =============================== 计算分配 ================================== ##
# 如果样本数少于脚本数，调整脚本数
if [[ $TOTAL_SAMPLES -lt $NUM_SCRIPTS ]]; then
    NUM_SCRIPTS=$TOTAL_SAMPLES
    echo "[INFO] 样本数少于脚本数，调整为 $NUM_SCRIPTS 个脚本"
fi

SAMPLES_PER_SCRIPT=$(( (TOTAL_SAMPLES + NUM_SCRIPTS - 1) / NUM_SCRIPTS ))

## =============================== 生成子脚本 ================================ ##
generate_subscript() {
    local SCRIPT_ID=$1
    local START_IDX=$2
    local END_IDX=$3
    local SCRIPT_FILE="${SCRIPT_OUTPUT_DIR}/batch_${SCRIPT_ID}.sh"
    
    cat > "$SCRIPT_FILE" << 'SCRIPT_HEADER'
#!/bin/bash
#*##############################################################################
#* 自动生成的批量处理脚本
#*##############################################################################

set -euo pipefail

SCRIPT_HEADER

    cat >> "$SCRIPT_FILE" << SCRIPT_CONFIG
## =============================== 配置 ====================================== ##
BASE_DIR='${BASE_DIR}'
CONFIG_DIR="\${BASE_DIR}/confs"
export NUMTS_CONFIG="${DEFAULT_CONFIG}"

REF_GRCh=${REF_GRCh}
CLUSTER_SCRIPT=\$BASE_DIR/script/0_1_找聚类.py
SAM2PSL_SCRIPT=\$BASE_DIR/script/0_2_sam2psl.py
BREAKPOINT_SCRIPT=\$BASE_DIR/script/0_3_找断点.py

THREADS=${THREADS}
OUTPUT_ROOT="${OUTPUT_ROOT}"

## =============================== 日志配置 ================================== ##
LOG_DIR="${OUTPUT_ROOT}/logs"
mkdir -p "\$LOG_DIR"
SUCCESS_LOG="\${LOG_DIR}/success.log"
FAILURE_LOG="\${LOG_DIR}/failure.log"
LOCK_FILE="\${LOG_DIR}/.log.lock"

# 使用 flock 防止并行写入竞争
log_success() {
    local SAMPLE_ID=\$1
    local TIMESTAMP=\$(date '+%Y-%m-%d %H:%M:%S')
    (
        flock -x 200
        echo "[\${TIMESTAMP}] \${SAMPLE_ID} - SUCCESS" >> "\$SUCCESS_LOG"
    ) 200>"\$LOCK_FILE"
}

log_failure() {
    local SAMPLE_ID=\$1
    local STEP=\$2
    local ERROR_MSG=\$3
    local TIMESTAMP=\$(date '+%Y-%m-%d %H:%M:%S')
    (
        flock -x 200
        echo "[\${TIMESTAMP}] \${SAMPLE_ID} - FAILED at \${STEP}: \${ERROR_MSG}" >> "\$FAILURE_LOG"
    ) 200>"\$LOCK_FILE"
}

echo ">>> 使用配置: \$NUMTS_CONFIG"
echo ">>> 日志目录: \$LOG_DIR"

# 从配置读取线粒体染色体名称
MITOCHONDRIAL_CHR=\$(
  python3 - "\$NUMTS_CONFIG" <<'PY'
import json
import sys

with open(sys.argv[1], "r") as f:
    cfg = json.load(f)
print(cfg.get("cluster", {}).get("mt_chrom", ""))
PY
)
if [[ -z "\$MITOCHONDRIAL_CHR" ]]; then
  echo "[ERROR] 无法从配置读取 mt_chrom: \$NUMTS_CONFIG" >&2
  exit 1
fi

###############################################################################
## 处理函数
###############################################################################
process_sample() {
    local INPUT_BAM=\$1
    local SAMPLE_ID=\$2
    local OUTPUT_DIR="\${OUTPUT_ROOT}/\${SAMPLE_ID}"
    local CURRENT_STEP=""
    
    echo "=============================================="
    echo ">>> 开始处理样本: \$SAMPLE_ID"
    echo ">>> BAM文件: \$INPUT_BAM"
    echo ">>> 输出目录: \$OUTPUT_DIR"
    echo "=============================================="
    
    # 错误处理函数
    handle_error() {
        local exit_code=\$?
        log_failure "\$SAMPLE_ID" "\$CURRENT_STEP" "Exit code: \$exit_code"
        echo "[ERROR] \$SAMPLE_ID 在步骤 \$CURRENT_STEP 失败，退出码: \$exit_code"
        return 1
    }
    
    # 设置错误陷阱
    trap 'handle_error' ERR
    
    mkdir -p "\$OUTPUT_DIR"
    
    DISC_SAM=\${OUTPUT_DIR}/\${SAMPLE_ID}.mt.disc.sam
    SPLIT_SAM=\${OUTPUT_DIR}/\${SAMPLE_ID}.mt.split.sam
    BREAK_IN=\${DISC_SAM}.breakpointINPUT.tsv
    COMBINED_FASTA=\${OUTPUT_DIR}/\${SAMPLE_ID}_all_regions.fasta
    BWA_SAM=\${OUTPUT_DIR}/\${SAMPLE_ID}_all_regions.bwa.sam
    ALL_PSL=\${OUTPUT_DIR}/\${SAMPLE_ID}_all_regions.psl
    
    ###########################################################################
    ## 1. 提取线粒体相关 read 并生成 disc/split SAM
    ###########################################################################
    CURRENT_STEP="step1_disc_split_SAM"
    if [[ ! -s "\$DISC_SAM" ]]; then
      echo ">>> step1  生成 disc/split SAM"
      samtools view -@\${THREADS} -m32G -h -F2 "\$INPUT_BAM" \\
        | grep -e '^@' -e "\${MITOCHONDRIAL_CHR}" \\
        | samtools sort -@\${THREADS} -m32G -n - \\
        | samtools view -h - \\
        | samblaster --ignoreUnmated -e \\
            -d "\$DISC_SAM" \\
            -s "\$SPLIT_SAM" \\
            -o /dev/null
    fi
    
    ###########################################################################
    ## 2. 运行聚类脚本，生成 breakpointINPUT.tsv
    ###########################################################################
    CURRENT_STEP="step2_clustering"
    if [[ ! -s "\$BREAK_IN" ]]; then
      echo ">>> step2  运行聚类脚本"
      python3 "\$CLUSTER_SCRIPT" "\$DISC_SAM" "\$SAMPLE_ID" "\$INPUT_BAM"
    fi
    
    ###########################################################################
    ## 3. 合并所有 region 的序列到一个 FASTA
    ###########################################################################
    CURRENT_STEP="step3_generate_FASTA"
    echo ">>> step3  生成合并 FASTA"
    > "\$COMBINED_FASTA"
    if [[ -s "\$BREAK_IN" ]]; then
        grep -v '^\$' "\$BREAK_IN" | while IFS=\$'\\t' read -r sampleID cluster_no \\
                disFile splitFile wgsBAM chr start end; do
            chr_clean=\$(echo "\$chr" | sed "s/[(),' ]//g")
            PREFIX="\${sampleID}_\${chr_clean}.\${start}.\${end}"
        
            samtools view "\$wgsBAM" "\${chr_clean}:\${start}-\${end}" \\
              | awk '\$6 !~ /150M|149M|148M|149S|148S/' \\
              | cut -f1,10 \\
              | while IFS=\$'\\t' read -r readID seq; do
                    printf ">%s|%s\\n%s\\n" "\$PREFIX" "\$readID" "\$seq" >> "\$COMBINED_FASTA"
                done
        done
    fi
    
    # 检查 FASTA 是否为空
    if [[ ! -s "\$COMBINED_FASTA" ]]; then
        echo "[WARN] \$SAMPLE_ID: 没有找到候选区域，跳过后续步骤"
        log_success "\$SAMPLE_ID"
        trap - ERR
        return 0
    fi
    
    ###########################################################################
    ## 4. bwa-mem 全量比对
    ###########################################################################
    CURRENT_STEP="step4_bwa_mem"
    echo ">>> step4  bwa-mem 全量比对"
    if [[ ! -f "\${REF_GRCh}.bwt" ]]; then
      echo "    building bwa index ..."
      bwa index "\$REF_GRCh"
    fi
    bwa mem -M -t\${THREADS} "\$REF_GRCh" "\$COMBINED_FASTA" > "\$BWA_SAM"
    
    ###########################################################################
    ## 5. SAM → PSL
    ###########################################################################
    CURRENT_STEP="step5_SAM_to_PSL"
    echo ">>> step5  SAM → PSL"
    python3 "\$SAM2PSL_SCRIPT" "\$BWA_SAM" "\${REF_GRCh}.fai" "\$ALL_PSL"
    
    ###########################################################################
    ## 6. 按 PREFIX 拆 PSL 并调用断点脚本
    ###########################################################################
    CURRENT_STEP="step6_breakpoint_detection"
    echo ">>> step6  拆分 PSL 并逐条调用断点脚本"
    PSL_HEADER=\$(head -n5 "\$ALL_PSL")
    
    grep -v '^\$' "\$BREAK_IN" | while IFS=\$'\\t' read -r sampleID cluster_no \\
            disFile splitFile wgsBAM chr start end; do
        chr_clean=\$(echo "\$chr" | sed "s/[(),' ]//g")
        PREFIX="\${sampleID}_\${chr_clean}.\${start}.\${end}"
        PSL_PART="\${OUTPUT_DIR}/\${PREFIX}.psl"
    
        {   printf "%s\\n" "\$PSL_HEADER"
            awk -F'\\t' -v q="\$PREFIX" '
                {
                  split(\$10, a, "|")
                  if (a[1]==q) print
                }' "\$ALL_PSL"
        } > "\$PSL_PART"
    
        if [[ \$(wc -l < "\$PSL_PART") -le 5 ]]; then
            echo "  [WARN] \$PREFIX 没在 PSL 找到比对，跳过"
            continue
        fi
        python3 "\$BREAKPOINT_SCRIPT" \\
            "\$PSL_PART" \\
            "\$sampleID" \\
            "\$chr_clean" \\
            "\$start" \\
            "\$end" \\
            "\${OUTPUT_DIR}/\${PREFIX}"
    done
    
    # 所有步骤完成，记录成功
    log_success "\$SAMPLE_ID"
    echo ">>> 样本 \$SAMPLE_ID 处理完成"
    
    # 清除错误陷阱
    trap - ERR
}

###############################################################################
## 样本列表
###############################################################################
SCRIPT_CONFIG

    # 添加样本处理调用
    for (( i=START_IDX; i<END_IDX && i<TOTAL_SAMPLES; i++ )); do
        IFS=$'\t' read -r BAM_PATH SAMPLE_NAME <<< "${SAMPLES[$i]}"
        echo "process_sample \"$BAM_PATH\" \"$SAMPLE_NAME\"" >> "$SCRIPT_FILE"
    done
    
    cat >> "$SCRIPT_FILE" << 'SCRIPT_FOOTER'

echo "=== BATCH SCRIPT FINISHED ==="
SCRIPT_FOOTER

    chmod +x "$SCRIPT_FILE"
    echo "生成脚本: $SCRIPT_FILE (样本 $((START_IDX+1)) - $END_IDX)"
}

## =============================== 生成所有子脚本 ============================ ##
echo ""
echo ">>> 开始生成子脚本..."

for (( i=0; i<NUM_SCRIPTS; i++ )); do
    START_IDX=$((i * SAMPLES_PER_SCRIPT))
    END_IDX=$((START_IDX + SAMPLES_PER_SCRIPT))
    
    # 确保不超过总样本数
    if [[ $END_IDX -gt $TOTAL_SAMPLES ]]; then
        END_IDX=$TOTAL_SAMPLES
    fi
    
    # 如果起始索引已超过样本数，跳过
    if [[ $START_IDX -ge $TOTAL_SAMPLES ]]; then
        break
    fi
    
    generate_subscript $((i+1)) $START_IDX $END_IDX
done

## =============================== 生成主运行脚本 ============================ ##
MASTER_SCRIPT="${SCRIPT_OUTPUT_DIR}/run_all.sh"
cat > "$MASTER_SCRIPT" << MASTER_HEADER
#!/bin/bash
#*##############################################################################
#* 主运行脚本 - 可选择串行或并行执行所有子脚本
#* 用法: 
#*   串行: bash run_all.sh
#*   并行: bash run_all.sh parallel
#*##############################################################################

cd "${SCRIPT_OUTPUT_DIR}"

MODE=\${1:-"serial"}

if [[ "\$MODE" == "parallel" ]]; then
    echo ">>> 并行模式: 同时运行所有子脚本"
MASTER_HEADER

for (( i=0; i<NUM_SCRIPTS; i++ )); do
    START_IDX=$((i * SAMPLES_PER_SCRIPT))
    if [[ $START_IDX -ge $TOTAL_SAMPLES ]]; then
        break
    fi
    echo "    nohup bash batch_$((i+1)).sh > batch_$((i+1)).log 2>&1 &" >> "$MASTER_SCRIPT"
done

cat >> "$MASTER_SCRIPT" << 'MASTER_MID'
    echo ">>> 所有子脚本已在后台启动，可通过 tail -f batch_*.log 查看进度"
else
    echo ">>> 串行模式: 依次运行所有子脚本"
MASTER_MID

for (( i=0; i<NUM_SCRIPTS; i++ )); do
    START_IDX=$((i * SAMPLES_PER_SCRIPT))
    if [[ $START_IDX -ge $TOTAL_SAMPLES ]]; then
        break
    fi
    echo "    echo \">>> 运行 batch_$((i+1)).sh\"" >> "$MASTER_SCRIPT"
    echo "    bash batch_$((i+1)).sh" >> "$MASTER_SCRIPT"
done

cat >> "$MASTER_SCRIPT" << 'MASTER_FOOTER'
fi

echo "=== ALL DONE ==="
MASTER_FOOTER

chmod +x "$MASTER_SCRIPT"

echo ""
echo "=========================================="
echo "脚本生成完成！"
echo "生成目录: $SCRIPT_OUTPUT_DIR"
echo ""
echo "使用方法:"
echo "  串行运行: bash ${MASTER_SCRIPT}"
echo "  并行运行: bash ${MASTER_SCRIPT} parallel"
echo ""
echo "或者单独运行某个子脚本:"
echo "  bash ${SCRIPT_OUTPUT_DIR}/batch_1.sh"
echo "=========================================="
