#!/usr/bin/env bash
# pipe/1-merge-vardetection-tsv.sh
# 将多样本 Vardetection 输出的三类 TSV 文件分别合并为单个结果文件。
# 用法：bash pipe/1-merge-vardetection-tsv.sh [--config <yaml>] [--force]
set -euo pipefail

SCRIPT_NAME=$(basename "$0")
PROJECT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
SCRIPT_DIR="$PROJECT_DIR/script"
LOG_DIR="$PROJECT_DIR/log"
LOG_FILE="$LOG_DIR/${SCRIPT_NAME}.log"
MERGE_HELPER="$SCRIPT_DIR/merge_tsv_group.sh"

mkdir -p "$LOG_DIR"

CONFIG_PATH="$PROJECT_DIR/conf/Config.yaml"
FORCE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG_PATH="$2"
      shift 2
      ;;
    --force)
      FORCE=1
      shift 1
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

CONFIG_PATH="${CONFIG_PATH//$'\r'/}"

if [[ ! -f "$CONFIG_PATH" ]]; then
  echo "Config file not found: $CONFIG_PATH" >&2
  exit 1
fi

if [[ ! -x "$MERGE_HELPER" ]]; then
  echo "Merge helper not found or not executable: $MERGE_HELPER" >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# 解析 YAML（仅支持本项目所需的平铺结构与 suffix_groups 分组）
# ---------------------------------------------------------------------------
parse_config() {
  awk '
    function emit_group() {
      if (current_key != "") {
        if (pattern == "" || output_name == "") {
          print "Invalid suffix_groups entry for key: " current_key > "/dev/stderr"
          exit 1
        }
        print "GROUP\t" current_key "\t" pattern "\t" output_name
      }
    }
    function trim(s) {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", s)
      return s
    }
    function unquote(s, first, last) {
      s = trim(s)
      sub(/[[:space:]]+#.*$/, "", s)
      first = substr(s, 1, 1)
      last = substr(s, length(s), 1)
      if ((first == "'"'"'" && last == "'"'"'") || (first == "\"" && last == "\"")) {
        s = substr(s, 2, length(s) - 2)
      }
      return s
    }
    /^[[:space:]]*#/ || /^[[:space:]]*$/ {
      next
    }
    /^[^[:space:]][^:]*:[[:space:]]*.*$/ {
      if (in_suffix) {
        emit_group()
        current_key = ""
        pattern = ""
        output_name = ""
      }
      top_key = $0
      sub(/:.*/, "", top_key)
      value = $0
      sub(/^[^:]*:[[:space:]]*/, "", value)

      if (top_key == "input_dir" || top_key == "output_dir" || top_key == "tmp_dir" || top_key == "jobs" || top_key == "force") {
        print "META\t" top_key "\t" unquote(value)
      } else if (top_key == "suffix_groups") {
        in_suffix = 1
      } else {
        in_suffix = 0
      }
      next
    }
    in_suffix && /^[[:space:]]{2}[^[:space:]][^:]*:[[:space:]]*$/ {
      emit_group()
      current_key = $0
      gsub(/^[[:space:]]+/, "", current_key)
      sub(/:.*/, "", current_key)
      pattern = ""
      output_name = ""
      next
    }
    in_suffix && current_key != "" && /^[[:space:]]{4}pattern:[[:space:]]*/ {
      line = $0
      sub(/^[^:]*:[[:space:]]*/, "", line)
      pattern = unquote(line)
      next
    }
    in_suffix && current_key != "" && /^[[:space:]]{4}output_name:[[:space:]]*/ {
      line = $0
      sub(/^[^:]*:[[:space:]]*/, "", line)
      output_name = unquote(line)
      next
    }
    END {
      if (in_suffix) {
        emit_group()
      }
    }
  ' "$CONFIG_PATH"
}

normalize_force() {
  case "${1,,}" in
    1|true|yes|y) echo 1 ;;
    0|false|no|n|"") echo 0 ;;
    *) echo "Invalid force value in config: $1" >&2; exit 1 ;;
  esac
}

MERGE_META_FILE=$(mktemp)
MERGE_GROUPS_FILE=$(mktemp)
trap 'rm -f "$MERGE_META_FILE" "$MERGE_GROUPS_FILE"' EXIT

parse_config | while IFS=$'\t' read -r record_type key value extra; do
  if [[ "$record_type" == "META" ]]; then
    printf '%s\t%s\n' "$key" "$value"
  elif [[ "$record_type" == "GROUP" ]]; then
    printf '%s\t%s\t%s\n' "$key" "$value" "$extra"
  fi
done | awk -F '\t' '
  $1 == "input_dir" || $1 == "output_dir" || $1 == "tmp_dir" || $1 == "jobs" || $1 == "force" {
    print > meta_file
    next
  }
  {
    print > groups_file
  }
' meta_file="$MERGE_META_FILE" groups_file="$MERGE_GROUPS_FILE"

INPUT_DIR=""
OUTPUT_DIR=""
TMP_DIR=""
JOBS=""

while IFS=$'\t' read -r key value; do
  case "$key" in
    input_dir)  INPUT_DIR="$value"  ;;
    output_dir) OUTPUT_DIR="$value" ;;
    tmp_dir)    TMP_DIR="$value"    ;;
    jobs)       JOBS="$value"       ;;
    force)
      if [[ "$FORCE" -ne 1 ]]; then
        FORCE=$(normalize_force "$value")
      fi
      ;;
  esac
done < "$MERGE_META_FILE"

if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" || -z "$TMP_DIR" || -z "$JOBS" ]]; then
  echo "Missing required config fields (input_dir / output_dir / tmp_dir / jobs) in $CONFIG_PATH" >&2
  exit 1
fi

if ! [[ "$JOBS" =~ ^[0-9]+$ ]] || [[ "$JOBS" -lt 1 ]]; then
  echo "Invalid jobs value in config: $JOBS" >&2
  exit 1
fi

if ! [[ "$FORCE" =~ ^[01]$ ]]; then
  echo "Invalid force value: $FORCE" >&2
  exit 1
fi

if [[ -f "$LOG_FILE" ]] && grep -q "SUCCESS" "$LOG_FILE" && [[ "$FORCE" -eq 0 ]]; then
  echo "[$(date '+%F %T')] SKIP (force=0 and existing SUCCESS in $LOG_FILE)"
  exit 0
fi

exec > >(tee -a "$LOG_FILE") 2>&1

echo "[$(date '+%F %T')] START  config=$CONFIG_PATH"

# 优先使用 GNU parallel；否则降级为 xargs
auto_parallel_cmd="parallel"
if ! command -v parallel >/dev/null 2>&1; then
  auto_parallel_cmd="xargs"
fi

mkdir -p "$OUTPUT_DIR" "$TMP_DIR"
KEYS_FILE="$TMP_DIR/merge_keys.txt"
LISTS_DIR="$TMP_DIR/merge_filelists"
MERGE_RUNNER="$TMP_DIR/merge_group_runner.sh"
mkdir -p "$LISTS_DIR"

# ---------------------------------------------------------------------------
# STEP1：为每个分组构建待合并文件列表
# ---------------------------------------------------------------------------
build_merge_lists() {
  local mergeable_count=0

  : > "$KEYS_FILE"
  rm -f "$LISTS_DIR"/*.txt

  while IFS=$'\t' read -r suffix_key pattern output_name; do
    local file_list_path="$LISTS_DIR/${suffix_key}.txt"
    LC_ALL=C find "$INPUT_DIR" -type f -name "$pattern" | LC_ALL=C sort > "$file_list_path"

    local file_count
    file_count=$(awk 'END { print NR + 0 }' "$file_list_path")

    if [[ "$file_count" -gt 0 ]]; then
      printf '%s\n' "$suffix_key" >> "$KEYS_FILE"
      mergeable_count=$((mergeable_count + 1))
    fi

    printf '%s\t%s\t%s\n' "$suffix_key" "$file_count" "$file_list_path"
  done < <(LC_ALL=C sort -t $'\t' -k1,1 "$MERGE_GROUPS_FILE")

  if [[ "$mergeable_count" -eq 0 ]]; then
    echo "No input files were discovered for any suffix group under: $INPUT_DIR" >&2
    exit 1
  fi
}

# ---------------------------------------------------------------------------
# 调用 merge helper 合并单个分组
# ---------------------------------------------------------------------------
merge_one_group() {
  local suffix_key="$1"
  local files_list="$2"
  local output_name="$3"
  local final_output="$OUTPUT_DIR/$output_name"

  bash "$MERGE_HELPER" \
    --suffix-key "$suffix_key" \
    --files-list "$files_list" \
    --output-path "$final_output" \
    --force "$FORCE" \
    --source-mode "none"
}

# 生成可被 parallel/xargs 调用的 runner 脚本（需要将函数导出到子进程）
cat > "$MERGE_RUNNER" <<'RUNNER_EOF'
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
RUNNER_EOF
chmod +x "$MERGE_RUNNER"

cat > "$TMP_DIR/merge_group_runner.env" <<ENV_EOF
OUTPUT_DIR=$(printf '%q' "$OUTPUT_DIR")
FORCE=$(printf '%q' "$FORCE")
LISTS_DIR=$(printf '%q' "$LISTS_DIR")
MERGE_GROUPS_FILE=$(printf '%q' "$MERGE_GROUPS_FILE")
MERGE_HELPER=$(printf '%q' "$MERGE_HELPER")
$(declare -f merge_one_group)
ENV_EOF

echo "[$(date '+%F %T')] STEP1 build merge file lists"
build_merge_lists

echo "[$(date '+%F %T')] STEP2 merge using prepared file lists"
if [[ "$auto_parallel_cmd" == "parallel" ]]; then
  parallel -j "$JOBS" --arg-file "$KEYS_FILE" "$MERGE_RUNNER" {} "$TMP_DIR"
else
  xargs -I {} -P "$JOBS" "$MERGE_RUNNER" {} "$TMP_DIR" < "$KEYS_FILE"
fi

echo "[$(date '+%F %T')] SUCCESS"
echo "[$(date '+%F %T')] Output: $OUTPUT_DIR"
