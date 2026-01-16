#!/bin/bash
#*===================================================================
#* 脚本名称: 1-liftover.sh
#* 描述:
#* 本脚本用于执行NUMTs的liftover操作，依赖于配置文件和Python脚本。主要流程包括读取配置、检查输入输出路径、执行liftover操作，并记录日志。支持通过--config参数指定配置文件路径。
#* 使用方法:
#*     bash 1-liftover.sh [--config 配置文件路径]
#* 参数说明:
#*     --config    指定配置文件路径（可选，默认路径为../conf/config_cluster_summary.yaml）
#* 主要流程:
#* 1. 设置基础目录、日志目录、临时目录和Python脚本目录。
#* 2. 解析命令行参数，支持自定义配置文件路径。
#* 3. 通过Python脚本读取配置文件中的input.path和output.path。
#* 4. 检查输入文件是否存在，输出路径是否有效。
#* 5. 检查日志文件，若已成功处理则跳过。
#* 6. 调用liftover_numts.py脚本进行liftover操作，并根据结果记录日志。
#* 依赖:
#* - conda环境中的ucsc-liftover工具
#* - Python脚本: config_reader.py, liftover_numts.py
#* 日志说明:
#* - 日志文件记录每次操作的状态（SUCCESS/FAIL/DONE），便于追踪和重复执行时跳过已完成任务。
#*===================================================================


#* conda install -c bioconda ucsc-liftover 安装
set -euo pipefail

script_name=$(basename "$0")
base_dir=$(cd "$(dirname "$0")/.." && pwd)
log_dir="$base_dir/log"
tmp_dir="$base_dir/tmp"
py_dir="$base_dir/python"

mkdir -p "$log_dir" "$tmp_dir"

config_path="$base_dir/conf/config_cluster_summary.yaml"

while [ $#* -gt 0 ]; do
    case "$1" in
        --config)
            config_path="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 1
            ;;
    esac

done

log_file="$log_dir/${script_name}.log"

input_path=$(python "$py_dir/config_reader.py" --config "$config_path" --key input.path)
output_path=$(python "$py_dir/config_reader.py" --config "$config_path" --key output.path)

if [ -z "$input_path" ] || [ -z "$output_path" ]; then
    echo "input.path or output.path missing in config" >&2
    exit 1
fi

if [ ! -f "$input_path" ]; then
    echo "input file not found: $input_path" >&2
    exit 1
fi

if [ -f "$log_file" ]; then
    if awk -F'\t' -v in_path="$input_path" -v out_path="$output_path" '$1=="SUCCESS" && $2==in_path && $3==out_path {found=1} END{exit !found}' "$log_file"; then
        printf "SKIP\t%s\t%s\n" "$script_name" "already succeeded" >> "$log_file"
        exit 0
    fi
fi

set +e
python "$py_dir/liftover_numts.py" --input "$input_path" --output "$output_path" --config "$config_path"
status=$?
set -e
ts=$(date -Iseconds)
if [ $status -eq 0 ]; then
    printf "SUCCESS\t%s\t%s\t%s\n" "$input_path" "$output_path" "$ts" >> "$log_file"
else
    printf "FAIL\t%s\t%s\t%s\n" "$input_path" "$output_path" "$ts" >> "$log_file"
    exit $status
fi

printf "DONE\t%s\t%s\n" "$script_name" "completed" >> "$log_file"
