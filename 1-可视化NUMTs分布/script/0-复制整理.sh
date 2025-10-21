#!/bin/bash
#BSUB -J Example                    # 作业名称（可修改）
#BSUB -q q6240                      # 使用的队列（q6240: 48核/节点；q8358: 64核/节点）
#BSUB -n 8                         # 请求 8 个 CPU 核心
#BSUB -R "span[hosts=1]"            # 所有核心集中在同一节点
#BSUB -W 04:00                      # 最大运行时间 4 小时，可调整
#BSUB -o %J.Example.out             # 标准输出日志文件 (%J 为 Job ID)
#BSUB -e %J.Example.err             # 错误输出日志文件

# -------------------------------
# 环境设置
# -------------------------------
# 定义彩色输出函数
echo_black()  { echo -e "\033[1;30m$*\033[0m"; }
echo_red()    { echo -e "\033[1;31m$*\033[0m"; }
echo_green()  { echo -e "\033[1;32m$*\033[0m"; }
echo_yellow() { echo -e "\033[1;33m$*\033[0m"; }
echo_blue()   { echo -e "\033[1;34m$*\033[0m"; }
echo_magenta(){ echo -e "\033[1;35m$*\033[0m"; }  
echo_cyan()   { echo -e "\033[1;36m$*\033[0m"; } 
echo_white()  { echo -e "\033[37m$*\033[0m"; }

# 使用示例
echo_cyan "=== 任务开始执行 ==="
echo_green "黑色" "红色" "绿色" "黄色" "蓝色" "洋红" "青色" "白色"

# -------------------------------
# 输入输出文件（请根据需要修改）
# -------------------------------
INPUT="/share/home/grp-wangyf/luolintao/1-luolintao/2-NUMTs/1-整理/script"
OUTPUT="/share/home/grp-wangyf/luolintao/1-luolintao/2-NUMTs/1-整理/output/"

mkdir -p ${OUTPUT}
# -------------------------------
# 执行命令
# -------------------------------
echo_yellow "切换到输入目录: ${INPUT}"
cd "${INPUT}" || {
    echo_red "错误：无法切换到目录 ${INPUT}"
    exit 1
}

echo_cyan "开始查找并复制 .tsv 文件..."

# 方法1：使用 find 结合 xargs 进行并行处理
echo_magenta "使用方法1: find + xargs 并行复制"
find . -type f -name "*.tsv" -print0 | xargs -0 -P 8 -I {} cp {} "${OUTPUT}"

# 方法2：使用 find 的 -exec 结合 + 号进行批量处理（效率更高）
# echo_magenta "使用方法2: find -exec + 批量复制"
# find . -type f -name "*.tsv" -exec cp -t "${OUTPUT}" {} +

# 方法3：使用 while 循环逐文件处理（单线程，作为备选）
# echo_magenta "使用方法3: while 循环逐文件复制"
# find . -type f -name "*.tsv" | while read -r file; do
#     cp "$file" "${OUTPUT}"
# done

echo_green "文件复制完成。"

# 验证复制结果
echo_cyan "验证复制结果..."
count_original=$(find . -type f -name "*.tsv" | wc -l)
count_copied=$(find "${OUTPUT}" -type f -name "*.tsv" | wc -l)

echo_blue "原始目录中找到 .tsv 文件数量: ${count_original}"
echo_blue "目标目录中复制的 .tsv 文件数量: ${count_copied}"

if [ "${count_original}" -eq "${count_copied}" ]; then
    echo_green "✓ 文件复制验证成功，数量匹配"
else
    echo_red "✗ 文件复制验证失败，数量不匹配"
fi

echo_blue "任务已经全部完成。"
