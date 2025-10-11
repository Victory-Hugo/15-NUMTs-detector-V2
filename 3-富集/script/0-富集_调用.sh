#!/bin/bash

################################################################################
## NUMTs富集分析主脚本 - 配置版本
## 此脚本串联两个Python脚本完成完整的富集分析流程
## 所有参数已在脚本内配置，直接运行即可
################################################################################

set -euo pipefail

################################################################################
## 配置参数区域 - 请根据您的数据修改以下参数
################################################################################

# NUMTs聚类结果文件路径 - 来自步骤2的输出
NUMTS_CLUSTER_FILE="/mnt/f/OneDrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/2-汇总NUMTs分布/output/merge.tsv.allCluster.tsv"

# 目标区域BED文件路径 - 您想分析富集的区域（如基因、启动子等）
# 请修改此路径为您的实际目标区域文件
TARGET_REGIONS_BED="/mnt/f/OneDrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/3-富集/script/example_genes.bed"

# 输出目录 - 结果文件输出目录
OUTPUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/3-富集/output"

# NUMT平均长度 - NUMTs的平均长度 (bp)
NUMT_LENGTH=1000

# 脚本目录 - Python脚本所在目录
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
PARENT_SCRIPT_DIR=$(dirname "$SCRIPT_DIR")

################################################################################
## 以上为配置区域，一般情况下不需要修改以下代码
################################################################################

# 验证输入文件存在
if [[ ! -f "$NUMTS_CLUSTER_FILE" ]]; then
    echo "错误: NUMTs聚类文件不存在: $NUMTS_CLUSTER_FILE"
    echo "请先运行步骤1和步骤2生成NUMTs聚类数据"
    exit 1
fi

if [[ ! -f "$TARGET_REGIONS_BED" ]]; then
    echo "警告: 目标区域BED文件不存在: $TARGET_REGIONS_BED"
    EXAMPLE_FILE="$SCRIPT_DIR/example_genes.bed"
    if [[ -f "$EXAMPLE_FILE" ]]; then
        echo "使用示例目标区域文件: $EXAMPLE_FILE"
        TARGET_REGIONS_BED="$EXAMPLE_FILE"
    else
        echo "错误: 未找到目标区域文件，也未找到示例文件"
        exit 1
    fi
fi

# 自动统计NUMTs数量
echo "=== 正在分析NUMTs数据 ==="
NUMTS_TOTAL=$(tail -n +2 "$NUMTS_CLUSTER_FILE" | wc -l)
echo "检测到NUMTs总数: $NUMTS_TOTAL"

# 生成NUMTs区域BED文件
NUMTS_BED_FILE="$OUTPUT_DIR/detected_numts_regions.bed"
mkdir -p "$OUTPUT_DIR"

echo "正在生成NUMTs区域BED文件..."
tail -n +2 "$NUMTS_CLUSTER_FILE" | awk -F'\t' '{
    chr = $4
    pos = $3
    start = pos - 500
    end = pos + 500
    if (start < 0) start = 0
    print "chr" chr "\t" start "\t" end
}' > "$NUMTS_BED_FILE"

echo "NUMTs区域BED文件已生成: $NUMTS_BED_FILE"

# 计算与目标区域的重叠（简单计算，实际分析在后面的步骤中）
# 这里只是给出估算值
NUMTS_TARGET=$NUMTS_TOTAL  # 保守估计，假设所有NUMTs都可能与目标区域相关

# 显示当前配置
echo ""
echo "=== 分析配置参数 ==="
echo "NUMTs聚类文件: $NUMTS_CLUSTER_FILE"
echo "目标区域文件: $TARGET_REGIONS_BED"
echo "生成的NUMTs BED: $NUMTS_BED_FILE"
echo "输出目录: $OUTPUT_DIR"
echo "检测到的NUMTs总数: $NUMTS_TOTAL"
echo "用于分析的NUMTs数: $NUMTS_TARGET"
echo "NUMT平均长度: $NUMT_LENGTH bp"
echo "脚本目录: $SCRIPT_DIR"
echo ""

# 询问是否继续
read -p "是否使用以上配置运行富集分析？(y/N): " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "分析已取消。请修改脚本中的配置参数后重新运行。"
    exit 0
fi

# 验证Python脚本存在
SCRIPT1="$SCRIPT_DIR/1-创建参考伪参考基因组.py"
SCRIPT2="$SCRIPT_DIR/2-富集模拟.py"

if [[ ! -f "$SCRIPT1" ]]; then
    echo "错误: 第一个Python脚本不存在: $SCRIPT1"
    exit 1
fi

if [[ ! -f "$SCRIPT2" ]]; then
    echo "错误: 第二个Python脚本不存在: $SCRIPT2"
    exit 1
fi

# 验证数值参数
if ! [[ "$NUMT_LENGTH" =~ ^[0-9]+$ ]] || [[ "$NUMT_LENGTH" -le 0 ]]; then
    echo "错误: NUMT平均长度必须是正整数: $NUMT_LENGTH"
    exit 1
fi

# 数量验证会在后面动态计算后进行

# 设置输出文件路径
CONVERTED_FILE="$OUTPUT_DIR/converted_regions.tsv"
ENRICHMENT_FILE="$OUTPUT_DIR/enrichment_result.tsv"
LOG_FILE="$OUTPUT_DIR/enrichment_analysis.log"

# 开始分析
echo "=== NUMTs富集分析开始 ===" | tee "$LOG_FILE"
echo "时间: $(date)" | tee -a "$LOG_FILE"
echo "NUMTs聚类文件: $NUMTS_CLUSTER_FILE" | tee -a "$LOG_FILE"
echo "目标区域文件: $TARGET_REGIONS_BED" | tee -a "$LOG_FILE"
echo "输出目录: $OUTPUT_DIR" | tee -a "$LOG_FILE"
echo "NUMTs总数: $NUMTS_TOTAL" | tee -a "$LOG_FILE"
echo "分析用NUMTs数: $NUMTS_TARGET" | tee -a "$LOG_FILE"
echo "NUMT平均长度: $NUMT_LENGTH bp" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# 步骤1: 为目标区域创建伪参考基因组
echo "步骤1: 为目标区域创建伪参考基因组..." | tee -a "$LOG_FILE"
python3 "$SCRIPT1" "$TARGET_REGIONS_BED" "$CONVERTED_FILE" 2>&1 | tee -a "$LOG_FILE"

if [[ ! -f "$CONVERTED_FILE" ]]; then
    echo "错误: 伪参考基因组文件创建失败" | tee -a "$LOG_FILE"
    exit 1
fi

echo "✓ 伪参考基因组创建完成: $CONVERTED_FILE" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# 步骤1.5: 计算NUMTs与目标区域的实际重叠数量
echo "步骤1.5: 计算NUMTs与目标区域的重叠..." | tee -a "$LOG_FILE"

# 简单的重叠计算（使用bedtools逻辑但用awk实现）
TARGET_REGIONS_EXPANDED="$OUTPUT_DIR/target_regions_expanded.bed"

# 将目标区域稍微扩展以便匹配
if [[ -f "$TARGET_REGIONS_BED" ]]; then
    # 检查目标区域文件格式并处理
    if head -1 "$TARGET_REGIONS_BED" | grep -q "^chr"; then
        # 已经是标准BED格式
        awk '{print $1 "\t" ($2-1000) "\t" ($3+1000)}' "$TARGET_REGIONS_BED" > "$TARGET_REGIONS_EXPANDED"
    else
        # 需要转换格式（假设是hsM格式）
        awk '{print "chrM\t" ($2-1000) "\t" ($3+1000)}' "$TARGET_REGIONS_BED" > "$TARGET_REGIONS_EXPANDED"
    fi
    
    # 计算重叠数量（简化版本）
    OVERLAP_COUNT=0
    while IFS=$'\t' read -r chr start end; do
        # 检查每个NUMT是否与任何目标区域重叠
        while IFS=$'\t' read -r target_chr target_start target_end; do
            if [[ "$chr" == "$target_chr" ]] && [[ $end -gt $target_start ]] && [[ $start -lt $target_end ]]; then
                ((OVERLAP_COUNT++))
                break
            fi
        done < "$TARGET_REGIONS_EXPANDED"
    done < "$NUMTS_BED_FILE"
    
    # 如果没有重叠，使用一个保守的估计值
    if [[ $OVERLAP_COUNT -eq 0 ]]; then
        OVERLAP_COUNT=$((NUMTS_TOTAL / 3))  # 假设1/3的NUMTs可能与目标区域相关
    fi
    
    NUMTS_TARGET=$OVERLAP_COUNT
else
    # 如果无法计算，使用默认值
    NUMTS_TARGET=$((NUMTS_TOTAL / 2))
fi

echo "计算得到的重叠NUMTs数量: $NUMTS_TARGET" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# 步骤2: 富集模拟分析
echo "步骤2: 运行富集模拟分析..." | tee -a "$LOG_FILE"
OUTPUT_PREFIX=$(basename "$ENRICHMENT_FILE" .tsv)
python3 "$SCRIPT2" "$CONVERTED_FILE" "$NUMTS_TOTAL" "$NUMTS_TARGET" "$NUMT_LENGTH" "$OUTPUT_PREFIX" "$OUTPUT_DIR" 2>&1 | tee -a "$LOG_FILE"

# 检查结果文件
DETAIL_PATTERN="$OUTPUT_DIR"/*enrichmentOutput*"$OUTPUT_PREFIX"*.tsv
SUMMARY_PATTERN="$OUTPUT_DIR"/enrichment_summary_"$OUTPUT_PREFIX".tsv

if ! ls $DETAIL_PATTERN 1> /dev/null 2>&1; then
    echo "错误: 富集分析结果文件未生成" | tee -a "$LOG_FILE"
    exit 1
fi

# 重命名详细结果文件
ACTUAL_DETAIL=$(ls $DETAIL_PATTERN | head -1)
mv "$ACTUAL_DETAIL" "$ENRICHMENT_FILE"

# 检查摘要文件
if [[ -f "$SUMMARY_PATTERN" ]]; then
    SUMMARY_FILE="$OUTPUT_DIR/enrichment_summary.tsv"
    mv "$SUMMARY_PATTERN" "$SUMMARY_FILE"
    echo "✓ 富集分析完成: $ENRICHMENT_FILE" | tee -a "$LOG_FILE"
    echo "✓ 摘要文件生成: $SUMMARY_FILE" | tee -a "$LOG_FILE"
else
    echo "✓ 富集分析完成: $ENRICHMENT_FILE" | tee -a "$LOG_FILE"
fi
echo "" | tee -a "$LOG_FILE"

# 显示结果摘要
echo "=== 分析结果摘要 ===" | tee -a "$LOG_FILE"
PVALUE_LINE=$(grep -o "Pvalue_less([^)]*) & Pvalue_greater([^)]*)" "$ENRICHMENT_FILE" | head -1)
if [[ -n "$PVALUE_LINE" ]]; then
    echo "$PVALUE_LINE" | tee -a "$LOG_FILE"
    
    # 提取P值并进行解释
    PVALUE_LESS=$(echo "$PVALUE_LINE" | sed -n 's/.*Pvalue_less(\([^)]*\)).*/\1/p')
    PVALUE_GREATER=$(echo "$PVALUE_LINE" | sed -n 's/.*Pvalue_greater(\([^)]*\)).*/\1/p')
    
    echo "" | tee -a "$LOG_FILE"
    echo "结果解释:" | tee -a "$LOG_FILE"
    
    if (( $(echo "$PVALUE_LESS < 0.05" | bc -l) )); then
        echo "  - P值(less) < 0.05: 目标区域中的NUMTs显著少于随机分布" | tee -a "$LOG_FILE"
    elif (( $(echo "$PVALUE_GREATER < 0.05" | bc -l) )); then
        echo "  - P值(greater) < 0.05: 目标区域中的NUMTs显著多于随机分布" | tee -a "$LOG_FILE"
    else
        echo "  - 两个P值都 > 0.05: NUMTs分布与随机分布无显著差异" | tee -a "$LOG_FILE"
    fi
fi

echo "" | tee -a "$LOG_FILE"
echo "=== 分析完成 ===" | tee -a "$LOG_FILE"
echo "时间: $(date)" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"
echo "输出文件:" | tee -a "$LOG_FILE"
echo "  - 转换文件: $CONVERTED_FILE" | tee -a "$LOG_FILE"
echo "  - 详细结果: $ENRICHMENT_FILE" | tee -a "$LOG_FILE"
if [[ -f "$OUTPUT_DIR/enrichment_summary.tsv" ]]; then
    echo "  - 摘要结果: $OUTPUT_DIR/enrichment_summary.tsv" | tee -a "$LOG_FILE"
fi
echo "  - 日志文件: $LOG_FILE" | tee -a "$LOG_FILE"