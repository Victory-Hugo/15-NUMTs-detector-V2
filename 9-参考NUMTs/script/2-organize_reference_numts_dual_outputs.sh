#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "$0")" && pwd)
base_dir=$(cd "$script_dir/.." && pwd)
config_path="$base_dir/conf/config.yaml"

while [[ $# -gt 0 ]]; do
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

python_bin=$(python "$base_dir/python/config_reader.py" --config "$config_path" --key runtime.python_bin)
strict_output=$("$python_bin" "$base_dir/python/config_reader.py" --config "$config_path" --key dual_output.nuc_only.strict_path)
strict_summary=$("$python_bin" "$base_dir/python/config_reader.py" --config "$config_path" --key dual_output.nuc_only.strict_summary_path)
all_output=$("$python_bin" "$base_dir/python/config_reader.py" --config "$config_path" --key dual_output.nuc_only.all_sources_path)
all_summary=$("$python_bin" "$base_dir/python/config_reader.py" --config "$config_path" --key dual_output.nuc_only.all_sources_summary_path)

"$base_dir/script/1-organize_reference_numts.sh" \
  --config "$config_path" \
  --dedup-mode nuc_only \
  --output "$strict_output" \
  --summary "$strict_summary"

"$base_dir/script/1-organize_reference_numts.sh" \
  --config "$config_path" \
  --dedup-mode nuc_only \
  --output "$all_output" \
  --summary "$all_summary" \
  --ignore-source-filter
