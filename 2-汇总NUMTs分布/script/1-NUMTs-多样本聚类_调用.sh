#!/bin/bash

################################################################################
## 多样本 NUMTs 聚类分析调用脚本
## 功能：合并多个样本的 breakpointINPUT.tsv 文件，进行跨样本聚类分析
################################################################################

# 设置脚本在遇到错误时停止执行
set -e

#? 示例路径
# BASE_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/"
# INPUT_DIR="${BASE_DIR}/2-汇总NUMTs分布/data"
# OUTPUT_DIR="${BASE_DIR}/2-汇总NUMTs分布/output"
# SCRIPT_PATH="${BASE_DIR}/2-汇总NUMTs分布/script/1-NUMTs-多样本聚类.py"

BASE_DIR="/mnt/f/12_甘肃藏族_NUMT/1-所有的NUMTs/"
INPUT_DIR="${BASE_DIR}/data"
OUTPUT_DIR="${BASE_DIR}/output"
SCRIPT_PATH="/mnt/f/OneDrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/2-汇总NUMTs分布/script/1-NUMTs-多样本聚类.py"

# 创建输出目录（如果不存在）
mkdir -p "$OUTPUT_DIR"

# 定义合并后的文件名
MERGED_FILE="$OUTPUT_DIR/merge.tsv"

echo "========================================="
echo "多样本 NUMTs 聚类分析开始"
echo "========================================="
echo "输入目录: $INPUT_DIR"
echo "输出目录: $OUTPUT_DIR"
echo "处理脚本: $SCRIPT_PATH"
echo

# 检查输入目录是否存在
if [ ! -d "$INPUT_DIR" ]; then
    echo "错误: 输入目录不存在: $INPUT_DIR"
    exit 1
fi

# 检查Python脚本是否存在
if [ ! -f "$SCRIPT_PATH" ]; then
    echo "错误: Python脚本不存在: $SCRIPT_PATH"
    exit 1
fi

# 查找所有的 *mt.disc.sam.breakpointINPUT.tsv 文件
echo "正在查找 *mt.disc.sam.breakpointINPUT.tsv 文件..."
INPUT_FILES=($(find "$INPUT_DIR" -name "*mt.disc.sam.breakpointINPUT.tsv" -type f))

# 检查是否找到文件
if [ ${#INPUT_FILES[@]} -eq 0 ]; then
    echo "警告: 在 $INPUT_DIR 中未找到任何 *mt.disc.sam.breakpointINPUT.tsv 文件"
    exit 1
fi

echo "找到 ${#INPUT_FILES[@]} 个输入文件:"
for file in "${INPUT_FILES[@]}"; do
    echo "  - $(basename "$file")"
done
echo

# 合并所有输入文件
echo "正在合并输入文件到: $MERGED_FILE"
# 清空或创建合并文件
> "$MERGED_FILE"

# 逐个添加文件内容
for file in "${INPUT_FILES[@]}"; do
    if [ -f "$file" ] && [ -s "$file" ]; then
        echo "  正在处理: $(basename "$file")"
        cat "$file" >> "$MERGED_FILE"
    else
        echo "  警告: 文件为空或不存在: $(basename "$file")"
    fi
done

# 检查合并文件是否成功创建且非空
if [ ! -f "$MERGED_FILE" ] || [ ! -s "$MERGED_FILE" ]; then
    echo "错误: 合并文件创建失败或为空: $MERGED_FILE"
    exit 1
fi

echo "合并完成! 合并文件包含 $(wc -l < "$MERGED_FILE") 行数据"
echo

# 运行Python聚类分析脚本
echo "正在运行多样本聚类分析..."
echo "命令: python3 \"$SCRIPT_PATH\" \"$MERGED_FILE\""

# 切换到输出目录运行脚本（确保输出文件在正确位置）
cd "$OUTPUT_DIR"
python3 "$SCRIPT_PATH" "$MERGED_FILE"

# 检查输出文件是否生成
OUTPUT_FILES=(
    "${MERGED_FILE}.allCluster.tsv"
    "${MERGED_FILE}.allCluster.sum.tsv"
)

echo
echo "检查输出文件..."
for output_file in "${OUTPUT_FILES[@]}"; do
    if [ -f "$output_file" ]; then
        echo "✓ 成功生成: $(basename "$output_file") ($(wc -l < "$output_file") 行)"
    else
        echo "✗ 未生成: $(basename "$output_file")"
    fi
done

echo
echo "========================================="
echo "多样本 NUMTs 聚类分析完成!"
echo "========================================="
echo "输出文件位置:"
echo "  - 合并输入文件: $MERGED_FILE"
echo "  - 详细聚类结果: ${MERGED_FILE}.allCluster.tsv"
echo "  - 聚类汇总统计: ${MERGED_FILE}.allCluster.sum.tsv"
echo

# 显示简要统计信息
if [ -f "${MERGED_FILE}.allCluster.sum.tsv" ]; then
    echo "聚类统计概览:"
    echo "  总聚类数: $(tail -n +2 "${MERGED_FILE}.allCluster.sum.tsv" | wc -l)"
    echo "  涉及染色体: $(tail -n +2 "${MERGED_FILE}.allCluster.tsv" | cut -f4 | sort -u | tr '\n' ' ')"
fi

echo "分析完成!"
