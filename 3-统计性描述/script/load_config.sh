#!/usr/bin/env bash

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    echo "[ERROR] Source this script instead: source script/load_config.sh <config_path>" >&2
    exit 1
fi

_config_path="${1:-}"
if [[ -z "$_config_path" ]]; then
    echo "[ERROR] Missing config path. Usage: source script/load_config.sh <config_path>" >&2
    return 1
fi

if [[ ! -f "$_config_path" ]]; then
    echo "[ERROR] Config file not found: $_config_path" >&2
    return 1
fi

mapfile -t _config_exports < <(
    awk '
    function trim(value) {
        sub(/^[[:space:]]+/, "", value)
        sub(/[[:space:]]+$/, "", value)
        return value
    }

    function unquote(value, first_char, last_char) {
        first_char = substr(value, 1, 1)
        last_char = substr(value, length(value), 1)
        if ((first_char == "\"" && last_char == "\"") || (first_char == "'"'"'" && last_char == "'"'"'")) {
            return substr(value, 2, length(value) - 2)
        }
        return value
    }

    function shq(value) {
        gsub(/\047/, "'\''\"'\''\"'\''", value)
        return "'"'"'" value "'"'"'"
    }

    /^[[:space:]]*($|#)/ {
        next
    }

    match($0, /^([A-Za-z0-9_-]+):[[:space:]]*$/, groups) {
        section = groups[1]
        next
    }

    match($0, /^  ([A-Za-z0-9_-]+):[[:space:]]*(.+)[[:space:]]*$/, groups) {
        if (section == "") {
            printf "[ERROR] Invalid YAML layout near line %d: %s\n", NR, $0 > "/dev/stderr"
            exit 1
        }

        key = groups[1]
        value = trim(groups[2])

        if (value ~ /^[^"'\'']+[[:space:]]+#/) {
            sub(/[[:space:]]+#.*/, "", value)
            value = trim(value)
        }

        value = unquote(value)
        var_name = toupper(section "_" key)
        gsub(/-/, "_", var_name)
        printf "export CFG_%s=%s\n", var_name, shq(value)
        next
    }

    {
        printf "[ERROR] Unsupported YAML syntax in %s at line %d: %s\n", FILENAME, NR, $0 > "/dev/stderr"
        exit 1
    }
    ' "$_config_path"
)

for _config_export in "${_config_exports[@]}"; do
    eval "$_config_export"
done

resolve_cfg_placeholders() {
    local max_passes=10
    local pass var_name var_value placeholder_key replacement changed found normalized_key
    local -a cfg_vars

    mapfile -t cfg_vars < <(compgen -A variable | grep '^CFG_')

    for ((pass = 1; pass <= max_passes; pass++)); do
        changed=0

        for var_name in "${cfg_vars[@]}"; do
            var_value="${!var_name}"

            while [[ "$var_value" =~ \{([A-Za-z0-9_-]+)\} ]]; do
                placeholder_key="${BASH_REMATCH[1]}"
                normalized_key="$(printf '%s' "$placeholder_key" | tr '[:lower:]-' '[:upper:]_')"
                replacement=""
                found=0

                for _candidate_var in "${cfg_vars[@]}"; do
                    if [[ "${_candidate_var#CFG_}" == "$normalized_key" || "${_candidate_var#CFG_}" == *"_$normalized_key" ]]; then
                        replacement="${!_candidate_var}"
                        found=1
                        break
                    fi
                done

                if [[ "$found" -eq 0 ]]; then
                    break
                fi

                var_value="${var_value//\{$placeholder_key\}/$replacement}"
                changed=1
            done

            printf -v "$var_name" '%s' "$var_value"
            export "$var_name"
        done

        if [[ "$changed" -eq 0 ]]; then
            return 0
        fi
    done

    return 1
}

if ! resolve_cfg_placeholders; then
    echo "[ERROR] Failed to fully resolve config placeholders within 10 passes: $_config_path" >&2
    return 1
fi

unset _config_path
unset _candidate_var
unset _config_export
unset _config_exports
unset cfg_vars
unset found
unset max_passes
unset normalized_key
unset pass
unset placeholder_key
unset replacement
unset changed
unset var_name
unset var_value
