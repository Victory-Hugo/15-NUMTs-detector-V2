#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "$0")" && pwd)
base_dir=$(cd "$script_dir/.." && pwd)
config_path="$base_dir/conf/config.yaml"
dedup_mode="event"
output_override=""
summary_override=""
ignore_source_filter="false"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      config_path="$2"
      shift 2
      ;;
    --dedup-mode)
      dedup_mode="$2"
      shift 2
      ;;
    --output)
      output_override="$2"
      shift 2
      ;;
    --summary)
      summary_override="$2"
      shift 2
      ;;
    --ignore-source-filter)
      ignore_source_filter="true"
      shift
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

python_bin=$(python "$base_dir/python/config_reader.py" --config "$config_path" --key runtime.python_bin)
input_path=$("$python_bin" "$base_dir/python/config_reader.py" --config "$config_path" --key input.path)
output_path=$("$python_bin" "$base_dir/python/config_reader.py" --config "$config_path" --key output.path)
summary_path=$("$python_bin" "$base_dir/python/config_reader.py" --config "$config_path" --key runtime.summary_path)
tmp_dir=$("$python_bin" "$base_dir/python/config_reader.py" --config "$config_path" --key runtime.tmp_dir)
liftover_bin=$("$python_bin" "$base_dir/python/config_reader.py" --config "$config_path" --key liftover.bin)
min_match=$("$python_bin" "$base_dir/python/config_reader.py" --config "$config_path" --key liftover.min_match)
hg18_chain=$("$python_bin" "$base_dir/python/config_reader.py" --config "$config_path" --key liftover.hg18_to_hg38_chain)
hg19_chain=$("$python_bin" "$base_dir/python/config_reader.py" --config "$config_path" --key liftover.hg19_to_hg38_chain)
hs1_chain=$("$python_bin" "$base_dir/python/config_reader.py" --config "$config_path" --key liftover.hs1_to_hg38_chain)
max_distance=$("$python_bin" "$base_dir/python/config_reader.py" --config "$config_path" --key dedup.max_breakpoint_distance)
mapfile -t include_sources < <("$python_bin" "$base_dir/python/config_reader.py" --config "$config_path" --key filters.include_sources 2>/dev/null || true)

if [[ "$dedup_mode" != "event" && "$dedup_mode" != "nuc_only" ]]; then
  echo "Unsupported dedup mode: $dedup_mode" >&2
  exit 1
fi

if [[ -n "$output_override" ]]; then
  output_path="$output_override"
fi

if [[ -n "$summary_override" ]]; then
  summary_path="$summary_override"
fi

if [[ ! -f "$input_path" ]]; then
  echo "Missing input file: $input_path" >&2
  exit 1
fi

mkdir -p "$(dirname "$output_path")" "$tmp_dir"

"$base_dir/script/0-download_chain_files.sh" --config "$config_path"

python_args=(
  --input "$input_path"
  --output "$output_path"
  --summary "$summary_path"
  --tmp-dir "$tmp_dir"
  --liftover-bin "$liftover_bin"
  --hg18-chain "$hg18_chain"
  --hg19-chain "$hg19_chain"
  --hs1-chain "$hs1_chain"
  --min-match "$min_match"
  --max-distance "$max_distance"
  --dedup-mode "$dedup_mode"
)

for source in "${include_sources[@]}"; do
  [[ -n "$source" ]] || continue
  python_args+=(--include-source "$source")
done

if [[ "$ignore_source_filter" == "true" ]]; then
  python_args+=(--ignore-source-filter)
fi

"$python_bin" "$base_dir/python/organize_reference_numts.py" "${python_args[@]}"
