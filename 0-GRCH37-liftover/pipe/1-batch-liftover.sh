#!/bin/bash
#*===================================================================
#* 脚本名称: 1-batch-liftover.sh
#* 描述:
#* 批量将 NUMTs TSV 文件从 GRCh37 liftover 到 GRCh38。
#* 从 conf/Config.yaml 读取所有配置，通过 python/batch_liftover.py 调度。
#* 使用方法:
#*     bash pipe/1-batch-liftover.sh [--config 配置文件路径]
#* 参数:
#*     --config    指定配置文件路径（可选，默认为 ../conf/Config.yaml）
#*===================================================================

set -euo pipefail

script_name=$(basename "$0")
base_dir=$(cd "$(dirname "$0")/.." && pwd)
py_dir="$base_dir/python"
config_path="$base_dir/conf/Config.yaml"

while [ $# -gt 0 ]; do
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

python_bin=$(python "$py_dir/config_reader.py" --config "$config_path" --key tools.python_bin 2>/dev/null || echo "python")
log_dir_rel=$(python "$py_dir/config_reader.py" --config "$config_path" --key runtime.log_dir 2>/dev/null || echo "log")
log_dir="$base_dir/$log_dir_rel"
mkdir -p "$log_dir"

log_file="$log_dir/${script_name}.log"
ts=$(date -Iseconds)

echo "[$ts] START $script_name config=$config_path" | tee -a "$log_file"

set +e
"$python_bin" "$py_dir/batch_liftover.py" --config "$config_path"
status=$?
set -e

ts=$(date -Iseconds)
if [ $status -eq 0 ]; then
    printf "SUCCESS\t%s\t%s\n" "$script_name" "$ts" | tee -a "$log_file"
else
    printf "FAIL\t%s\t%s\tstatus=%d\n" "$script_name" "$ts" "$status" | tee -a "$log_file"
    exit $status
fi
