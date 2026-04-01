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

download_chain() {
  local path_key="$1"
  local url_key="$2"
  local dst
  local url

  dst=$("$python_bin" "$base_dir/python/config_reader.py" --config "$config_path" --key "$path_key")
  url=$("$python_bin" "$base_dir/python/config_reader.py" --config "$config_path" --key "$url_key")

  if [[ -f "$dst" ]]; then
    echo "[download] exists: $dst"
    return 0
  fi

  mkdir -p "$(dirname "$dst")"
  echo "[download] $url -> $dst"
  curl -L --fail --output "$dst" "$url"
}

download_chain "liftover.hg18_to_hg38_chain" "liftover.hg18_to_hg38_chain_url"
download_chain "liftover.hg19_to_hg38_chain" "liftover.hg19_to_hg38_chain_url"
download_chain "liftover.hs1_to_hg38_chain" "liftover.hs1_to_hg38_chain_url"
