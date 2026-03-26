#!/usr/bin/env bash
set -euo pipefail

suffix_key="$1"
project_tmp_dir="$2"

# shellcheck disable=SC1090
source "$project_tmp_dir/merge_group_runner.env"

output_name=$(awk -F '\t' -v key="$suffix_key" '$1 == key { print $3; exit }' "$MERGE_GROUPS_FILE")
if [[ -z "$output_name" ]]; then
  echo "suffix_key not found in config: $suffix_key" >&2
  exit 1
fi

merge_one_group "$suffix_key" "$LISTS_DIR/${suffix_key}.txt" "$output_name"
