#!/usr/bin/env bash
set -euo pipefail

SUFFIX_KEY=""
FILES_LIST=""
OUTPUT_PATH=""
FORCE="0"
SOURCE_MODE="none"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --suffix-key)
      SUFFIX_KEY="$2"
      shift 2
      ;;
    --files-list)
      FILES_LIST="$2"
      shift 2
      ;;
    --output-path)
      OUTPUT_PATH="$2"
      shift 2
      ;;
    --force)
      FORCE="$2"
      shift 2
      ;;
    --source-mode)
      SOURCE_MODE="$2"
      shift 2
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

if [[ -z "$SUFFIX_KEY" || -z "$FILES_LIST" || -z "$OUTPUT_PATH" ]]; then
  echo "Usage: bash script/merge_tsv_group.sh --suffix-key <key> --files-list <txt> --output-path <tsv> [--force 0|1] [--source-mode none|breakpoint_region]" >&2
  exit 1
fi

if [[ ! -f "$FILES_LIST" ]]; then
  echo "No files found in list: $FILES_LIST" >&2
  exit 1
fi

if ! [[ "$FORCE" =~ ^[01]$ ]]; then
  echo "Invalid force value: $FORCE" >&2
  exit 1
fi

if [[ "$SOURCE_MODE" != "none" && "$SOURCE_MODE" != "breakpoint_region" ]]; then
  echo "Unsupported source mode: $SOURCE_MODE" >&2
  exit 1
fi

if [[ -f "$OUTPUT_PATH" && "$FORCE" -eq 0 ]]; then
  echo "[$(date '+%F %T')] SKIP existing output: $OUTPUT_PATH"
  exit 0
fi

mkdir -p "$(dirname "$OUTPUT_PATH")"
TEMP_OUTPUT="${OUTPUT_PATH}.tmp"

awk -v list_file="$FILES_LIST" -v out_file="$TEMP_OUTPUT" -v source_mode="$SOURCE_MODE" '
  function fail(message) {
    print message > "/dev/stderr"
    exit 1
  }

  function trim_cr(text) {
    gsub(/\r/, "", text)
    return text
  }

  function file_basename(path, parts, count) {
    count = split(path, parts, "/")
    return parts[count]
  }

  function source_header(mode) {
    if (mode == "breakpoint_region") {
      return "source_file_basename\tsource_sample_id\tsource_region_chr\tsource_region_start\tsource_region_end\tsource_region_key"
    }
    return ""
  }

  function build_source_metadata(file_path, mode, base_name, sample_id, region_chr, region_start, region_end, region_key, matches) {
    if (mode == "none") {
      return ""
    }

    base_name = file_basename(file_path)
    if (mode == "breakpoint_region") {
      if (match(base_name, /^(.*)_(chr[^.]+)\.([0-9]+)\.([0-9]+)\.(chr[^.]+)\.(AllBreakpoints|ConfidentBreakpoints)\.tsv$/, matches)) {
        sample_id = matches[1]
        region_chr = matches[2]
        region_start = matches[3]
        region_end = matches[4]
        region_key = sample_id "_" region_chr "." region_start "." region_end
        return base_name "\t" sample_id "\t" region_chr "\t" region_start "\t" region_end "\t" region_key
      }
      fail("Unable to parse breakpoint region from filename: " base_name)
    }

    fail("Unsupported source mode in awk: " mode)
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

    extra_header = source_header(source_mode)
    if (extra_header != "" && header == "") {
      fail("Source metadata requires a shared header, but none was detected for " list_file)
    }

    if (header != "") {
      if (extra_header != "") {
        print header "\t" extra_header > out_file
      } else {
        print header > out_file
      }
    }

    for (i = 1; i <= file_total; i++) {
      file = files[i]
      first = first_line[file]
      skip_first = (header != "" && first == header) ? 1 : 0

      if (header != "" && first != header && first ~ /^#/) {
        fail("Header mismatch: expected " header " but got " first " in " file)
      }

      metadata = build_source_metadata(file, source_mode)

      close(file)
      line_no = 0
      while ((getline line < file) > 0) {
        line_no++
        if (skip_first && line_no == 1) {
          continue
        }
        line = trim_cr(line)
        if (metadata != "") {
          print line "\t" metadata > out_file
        } else {
          print line > out_file
        }
      }
      close(file)
    }
  }
' "$FILES_LIST"

mv -f "$TEMP_OUTPUT" "$OUTPUT_PATH"
