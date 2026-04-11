#!/usr/bin/env bash
# script/merge_tsv_group.sh
# 通用 TSV 合并辅助脚本：自动检测共同 header 并去重，支持原子写入。
# 由 pipe/1-merge-vardetection-tsv.sh 调用，不应直接运行。
set -euo pipefail

SUFFIX_KEY=""
FILES_LIST=""
OUTPUT_PATH=""
FORCE="0"
SOURCE_MODE="none"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --suffix-key)   SUFFIX_KEY="$2";   shift 2 ;;
    --files-list)   FILES_LIST="$2";   shift 2 ;;
    --output-path)  OUTPUT_PATH="$2";  shift 2 ;;
    --force)        FORCE="$2";        shift 2 ;;
    --source-mode)  SOURCE_MODE="$2";  shift 2 ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

if [[ -z "$SUFFIX_KEY" || -z "$FILES_LIST" || -z "$OUTPUT_PATH" ]]; then
  echo "Usage: bash script/merge_tsv_group.sh --suffix-key <key> --files-list <txt> --output-path <tsv> [--force 0|1] [--source-mode none]" >&2
  exit 1
fi

if [[ ! -f "$FILES_LIST" ]]; then
  echo "File list not found: $FILES_LIST" >&2
  exit 1
fi

if ! [[ "$FORCE" =~ ^[01]$ ]]; then
  echo "Invalid force value: $FORCE" >&2
  exit 1
fi

if [[ "$SOURCE_MODE" != "none" ]]; then
  echo "Unsupported source mode: $SOURCE_MODE" >&2
  exit 1
fi

if [[ -f "$OUTPUT_PATH" && "$FORCE" -eq 0 ]]; then
  echo "[$(date '+%F %T')] SKIP existing output: $OUTPUT_PATH"
  exit 0
fi

mkdir -p "$(dirname "$OUTPUT_PATH")"
TEMP_OUTPUT="${OUTPUT_PATH}.tmp"

awk -v list_file="$FILES_LIST" -v out_file="$TEMP_OUTPUT" '
  function fail(message) {
    print message > "/dev/stderr"
    exit 1
  }

  function trim_cr(text) {
    gsub(/\r/, "", text)
    return text
  }

  {
    current_line = trim_cr($0)
    if (current_line == "" || current_line ~ /^#/) {
      next
    }

    file = current_line
    files[++file_total] = file

    status = (getline current_first < file)
    if (status < 0) {
      fail("Listed file does not exist: " file)
    }
    close(file)

    if (status == 0) {
      current_first = ""
    }

    current_first = trim_cr(current_first)
    first_line[file] = current_first
    counts[current_first]++
  }

  END {
    if (file_total == 0) {
      fail("No files found in list: " list_file)
    }

    # 取出现次数最多的首行作为 header（需至少出现 2 次才视为共同 header）
    header = ""
    header_count = 0
    for (line in counts) {
      if (counts[line] > header_count) {
        header = line
        header_count = counts[line]
      }
    }
    if (header_count < 2) {
      header = ""
    }

    if (header != "") {
      print header > out_file
    }

    for (i = 1; i <= file_total; i++) {
      file = files[i]
      first = first_line[file]
      skip_first = (header != "" && first == header) ? 1 : 0

      if (header != "" && first != header && first ~ /^#/) {
        fail("Header mismatch: expected [" header "] but got [" first "] in " file)
      }

      close(file)
      line_no = 0
      while ((getline line < file) > 0) {
        line_no++
        if (skip_first && line_no == 1) {
          continue
        }
        print trim_cr(line) > out_file
      }
      close(file)
    }
  }
' "$FILES_LIST"

mv -f "$TEMP_OUTPUT" "$OUTPUT_PATH"
echo "[$(date '+%F %T')] DONE  $SUFFIX_KEY -> $OUTPUT_PATH"
