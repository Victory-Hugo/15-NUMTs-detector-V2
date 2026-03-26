#!/usr/bin/env bash

set -euo pipefail

load_config() {
  local config_path="$1"
  if [[ ! -f "$config_path" ]]; then
    echo "Config file not found: $config_path" >&2
    return 1
  fi

  local current_section=""
  while IFS= read -r raw_line || [[ -n "$raw_line" ]]; do
    local line="${raw_line%%#*}"
    line="${line%$'\r'}"
    [[ -z "${line//[[:space:]]/}" ]] && continue

    if [[ "$line" =~ ^([A-Za-z0-9_]+):[[:space:]]*$ ]]; then
      current_section="${BASH_REMATCH[1]}"
      continue
    fi

    if [[ "$line" =~ ^[[:space:]]{2}([A-Za-z0-9_\.]+):[[:space:]]*(.*)$ ]]; then
      local key="${BASH_REMATCH[1]}"
      local value="${BASH_REMATCH[2]}"
      value="${value#\"}"
      value="${value%\"}"
      value="${value#\'}"
      value="${value%\'}"
      local var_name
      var_name="$(printf '%s_%s' "$current_section" "$key" | tr '[:lower:].' '[:upper:]_')"
      printf -v "$var_name" '%s' "$value"
      export "$var_name"
    fi
  done < "$config_path"
}
