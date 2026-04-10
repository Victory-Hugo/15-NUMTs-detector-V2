#!/usr/bin/env python3
import argparse
import csv
import os
import re
import subprocess
import sys
import uuid
from typing import Dict, List, Optional, Tuple, Union

import yaml


def load_config(path: str) -> Dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def normalize_chr(raw_chr: str, chr_prefix: str) -> str:
    chr_val = raw_chr.strip()
    if chr_val.lower().startswith(chr_prefix.lower()):
        return chr_val
    return f"{chr_prefix}{chr_val}"


def is_mt_chr(raw_chr: str, mt_chroms: List[str]) -> bool:
    return raw_chr.strip() in mt_chroms


def parse_pos(raw_pos: str) -> int:
    # Handles integer-like strings or floats stored as strings
    return int(float(raw_pos))


def parse_pos_list(raw_list: str) -> List[int]:
    return [int(x) for x in re.findall(r"\d+", raw_list)]


def format_pos_list(values: List[int]) -> str:
    return "[" + ", ".join(str(v) for v in values) + "]"


def parse_cluster_id(raw_value: str) -> Optional[Tuple[str, int, int, str]]:
    parts = raw_value.split("_")
    if len(parts) < 4:
        return None
    chr_part = parts[0]
    try:
        start = parse_pos(parts[1])
        end = parse_pos(parts[2])
    except ValueError:
        return None
    rest = "_".join(parts[3:])
    return chr_part, start, end, rest


def resolve_index(index: Optional[int], index_from_end: Optional[int], row_len: int) -> int:
    if index is not None:
        return index
    if index_from_end is not None:
        return row_len - index_from_end
    raise ValueError("Column index not configured.")


def _run_liftover_cmd(
    liftover_bin: str,
    chain: str,
    min_match: Optional[float],
    bed_in: str,
    bed_out: str,
    bed_unmapped: str,
) -> None:
    cmd = [liftover_bin]
    if min_match is not None:
        cmd.append(f"-minMatch={min_match}")
    cmd.extend([bed_in, chain, bed_out, bed_unmapped])
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            "liftOver failed:\n"
            f"STDOUT: {result.stdout}\n"
            f"STDERR: {result.stderr}"
        )


def _read_lifted_pos(bed_out: str) -> Dict[str, Tuple[int, int]]:
    lifted_pos: Dict[str, Tuple[int, int]] = {}
    with open(bed_out, "r", encoding="utf-8") as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            lifted_pos[fields[3]] = (int(fields[1]) + 1, int(fields[2]))
    return lifted_pos


def _cleanup(*paths: str) -> None:
    for path in paths:
        if os.path.exists(path):
            os.remove(path)


def run_dual_pos(input_path: str, output_path: str, config: Dict) -> Dict:
    """无表头双位置列格式：chr 在固定列，有效位置在 pos_col_index_a 或 pos_col_index_b（另一列为 -1）。"""
    delimiter = config["input"]["delimiter"]
    chr_col_idx: int = config["input"]["chr_col_index"]
    pos_col_a: int = config["input"]["pos_col_index_a"]
    pos_col_b: int = config["input"]["pos_col_index_b"]
    chr_prefix: str = config["input"]["chr_prefix"]
    mt_chroms: List[str] = config["input"]["mt_chroms"]
    liftover_bin: str = config["liftover"]["bin"]
    chain: str = config["liftover"]["chain"]
    min_match: Optional[float] = config["liftover"].get("min_match")
    tmp_dir: str = config["runtime"]["tmp_dir"]
    os.makedirs(tmp_dir, exist_ok=True)

    run_id = uuid.uuid4().hex
    bed_in = os.path.join(tmp_dir, f"liftover_in_{run_id}.bed")
    bed_out = os.path.join(tmp_dir, f"liftover_out_{run_id}.bed")
    bed_unmapped = os.path.join(tmp_dir, f"liftover_unmapped_{run_id}.bed")

    rows: List[List[str]] = []
    bed_lines: List[str] = []
    mt_row_ids: set = set()
    pos_a_submitted = 0
    pos_b_submitted = 0

    with open(input_path, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter=delimiter)
        for idx, row in enumerate(reader):
            if not row:
                continue
            rows.append(row)
            raw_chr = row[chr_col_idx]
            if is_mt_chr(raw_chr, mt_chroms):
                mt_row_ids.add(idx)
                continue
            chrom = normalize_chr(raw_chr, chr_prefix)
            val_a = row[pos_col_a].strip()
            if val_a not in ("-1", ""):
                bed_lines.append(f"{chrom}\t{parse_pos(val_a) - 1}\t{parse_pos(val_a)}\t{idx}:A")
                pos_a_submitted += 1
            val_b = row[pos_col_b].strip()
            if val_b not in ("-1", ""):
                bed_lines.append(f"{chrom}\t{parse_pos(val_b) - 1}\t{parse_pos(val_b)}\t{idx}:B")
                pos_b_submitted += 1

    lifted_pos: Dict[str, Tuple[int, int]] = {}
    if bed_lines:
        with open(bed_in, "w", encoding="utf-8") as f:
            f.write("\n".join(bed_lines) + "\n")
        _run_liftover_cmd(liftover_bin, chain, min_match, bed_in, bed_out, bed_unmapped)
        lifted_pos = _read_lifted_pos(bed_out)

    pos_a_lifted = sum(1 for k in lifted_pos if k.endswith(":A"))
    pos_b_lifted = sum(1 for k in lifted_pos if k.endswith(":B"))
    nuclear_rows = len(rows) - len(mt_row_ids)

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter=delimiter)
        for idx, row in enumerate(rows):
            if idx in mt_row_ids:
                writer.writerow(row)
                continue
            if f"{idx}:A" in lifted_pos:
                row[pos_col_a] = str(lifted_pos[f"{idx}:A"][0])
            if f"{idx}:B" in lifted_pos:
                row[pos_col_b] = str(lifted_pos[f"{idx}:B"][0])
            writer.writerow(row)

    _cleanup(bed_in, bed_out, bed_unmapped)
    return {
        "format": "dual_pos",
        "total_rows": len(rows),
        "mt_rows": len(mt_row_ids),
        "nuclear_rows": nuclear_rows,
        "pos_a_submitted": pos_a_submitted,
        "pos_a_lifted": pos_a_lifted,
        "pos_a_failed": pos_a_submitted - pos_a_lifted,
        "pos_b_submitted": pos_b_submitted,
        "pos_b_lifted": pos_b_lifted,
        "pos_b_failed": pos_b_submitted - pos_b_lifted,
        "on_fail_action": "keep",
        "rows_output": len(rows),
    }


def run_encoded_key(input_path: str, output_path: str, config: Dict) -> Dict:
    """无表头编码键格式：核基因组坐标编码于指定列，格式 {SampleID}_{chr}.{start}.{end}|{readname}。"""
    import re as _re
    KEY_PATTERN = _re.compile(r'^(.+)_([^_.]+)\.(\d+)\.(\d+)\|(.+)$')

    delimiter = config["input"]["delimiter"]
    key_col: int = config["input"]["key_col_index"]
    chr_prefix: str = config["input"]["chr_prefix"]
    mt_chroms: List[str] = config["input"]["mt_chroms"]
    liftover_bin: str = config["liftover"]["bin"]
    chain: str = config["liftover"]["chain"]
    min_match: Optional[float] = config["liftover"].get("min_match")
    tmp_dir: str = config["runtime"]["tmp_dir"]
    os.makedirs(tmp_dir, exist_ok=True)

    run_id = uuid.uuid4().hex
    bed_in = os.path.join(tmp_dir, f"liftover_in_{run_id}.bed")
    bed_out = os.path.join(tmp_dir, f"liftover_out_{run_id}.bed")
    bed_unmapped = os.path.join(tmp_dir, f"liftover_unmapped_{run_id}.bed")

    rows: List[List[str]] = []
    bed_lines: List[str] = []
    mt_row_ids: set = set()
    key_meta: Dict[int, Tuple[str, str, int, int, str]] = {}  # idx -> (sample_id, raw_chr, start, end, rest)

    with open(input_path, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter=delimiter)
        for idx, row in enumerate(reader):
            if not row:
                continue
            rows.append(row)
            m = KEY_PATTERN.match(row[key_col])
            if not m:
                continue
            sample_id, raw_chr, start_str, end_str, rest = m.groups()
            start, end = int(start_str), int(end_str)
            if is_mt_chr(raw_chr, mt_chroms):
                mt_row_ids.add(idx)
                continue
            chrom = normalize_chr(raw_chr, chr_prefix)
            key_meta[idx] = (sample_id, raw_chr, start, end, rest)
            bed_lines.append(f"{chrom}\t{start - 1}\t{end}\t{idx}")

    lifted_pos: Dict[str, Tuple[int, int]] = {}
    if bed_lines:
        with open(bed_in, "w", encoding="utf-8") as f:
            f.write("\n".join(bed_lines) + "\n")
        _run_liftover_cmd(liftover_bin, chain, min_match, bed_in, bed_out, bed_unmapped)
        lifted_pos = _read_lifted_pos(bed_out)

    keys_parsed = len(key_meta)
    keys_lifted = sum(1 for idx in key_meta if str(idx) in lifted_pos)
    keys_failed = keys_parsed - keys_lifted

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter=delimiter)
        for idx, row in enumerate(rows):
            row_key = str(idx)
            if idx in mt_row_ids or idx not in key_meta:
                writer.writerow(row)
                continue
            sample_id, raw_chr, orig_start, orig_end, rest = key_meta[idx]
            if row_key in lifted_pos:
                new_start, new_end = lifted_pos[row_key]
                row[key_col] = f"{sample_id}_{raw_chr}.{new_start}.{new_end}|{rest}"
            writer.writerow(row)

    _cleanup(bed_in, bed_out, bed_unmapped)
    return {
        "format": "encoded_key",
        "total_rows": len(rows),
        "mt_rows": len(mt_row_ids),
        "nuclear_rows": len(rows) - len(mt_row_ids),
        "keys_parsed": keys_parsed,
        "keys_lifted": keys_lifted,
        "keys_failed": keys_failed,
        "on_fail_action": "keep",
        "rows_output": len(rows),
    }


def run(
    input_path: str,
    output_path: str,
    config_path: Optional[str] = None,
    config: Optional[Dict] = None,
) -> None:
    if config is None:
        config = load_config(config_path)

    fmt = config["input"].get("format")
    if fmt == "dual_pos":
        return run_dual_pos(input_path, output_path, config)
    if fmt == "encoded_key":
        return run_encoded_key(input_path, output_path, config)

    delimiter = config["input"]["delimiter"]
    has_header = config["input"]["has_header"]
    chr_col = config["input"].get("chr_col")
    pos_col = config["input"].get("pos_col")
    start_col = config["input"].get("start_col")
    end_col = config["input"].get("end_col")
    reads_list_col = config["input"].get("reads_list_col")
    cluster_id_col = config["input"].get("cluster_id_col")
    chr_col_index = config["input"].get("chr_col_index")
    pos_col_index = config["input"].get("pos_col_index")
    start_col_index = config["input"].get("start_col_index")
    end_col_index = config["input"].get("end_col_index")
    chr_col_index_from_end = config["input"].get("chr_col_index_from_end")
    pos_col_index_from_end = config["input"].get("pos_col_index_from_end")
    start_col_index_from_end = config["input"].get("start_col_index_from_end")
    end_col_index_from_end = config["input"].get("end_col_index_from_end")
    cluster_id_col_index = config["input"].get("cluster_id_col_index")
    cluster_id_col_index_from_end = config["input"].get("cluster_id_col_index_from_end")
    chr_prefix = config["input"]["chr_prefix"]
    mt_chroms = config["input"]["mt_chroms"]

    liftover_bin = config["liftover"]["bin"]
    chain = config["liftover"]["chain"]
    min_match = config["liftover"].get("min_match")

    on_fail = config["output"]["on_fail"]

    tmp_dir = config["runtime"]["tmp_dir"]
    os.makedirs(tmp_dir, exist_ok=True)

    run_id = uuid.uuid4().hex
    bed_in = os.path.join(tmp_dir, f"liftover_in_{run_id}.bed")
    bed_out = os.path.join(tmp_dir, f"liftover_out_{run_id}.bed")
    bed_unmapped = os.path.join(tmp_dir, f"liftover_unmapped_{run_id}.bed")

    rows: List[Union[Dict[str, str], List[str]]] = []
    bed_lines: List[str] = []
    mt_row_ids = set()
    fieldnames: List[str] = []
    index_meta: Dict[str, int] = {}
    list_positions: Dict[str, List[int]] = {}
    cluster_meta: Dict[str, Tuple[str, str]] = {}

    with open(input_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=delimiter) if has_header else csv.reader(f, delimiter=delimiter)
        if has_header:
            fieldnames = reader.fieldnames or []
            if not fieldnames:
                raise ValueError("Header row is empty; cannot locate columns.")
            for idx, row in enumerate(reader):
                rows.append(row)
                if not chr_col:
                    raise ValueError("input.chr_col is required when has_header=true.")
                raw_chr = row[chr_col]
                raw_pos = row[pos_col] if pos_col else ""
                if is_mt_chr(raw_chr, mt_chroms):
                    mt_row_ids.add(idx)
                    continue
                if pos_col:
                    pos = parse_pos(raw_pos)
                    chrom = normalize_chr(raw_chr, chr_prefix)
                    bed_lines.append(f"{chrom}\t{pos - 1}\t{pos}\t{idx}")
                else:
                    if not start_col or not end_col:
                        raise ValueError("input.start_col and input.end_col required when pos_col is missing.")
                    start = parse_pos(row[start_col])
                    end = parse_pos(row[end_col])
                    chrom = normalize_chr(raw_chr, chr_prefix)
                    bed_lines.append(f"{chrom}\t{start - 1}\t{end}\t{idx}")
                if reads_list_col:
                    raw_list = row[reads_list_col]
                    pos_list = parse_pos_list(raw_list)
                    list_positions[str(idx)] = pos_list
                    chrom = normalize_chr(raw_chr, chr_prefix)
                    for j, p in enumerate(pos_list):
                        bed_lines.append(f"{chrom}\t{p - 1}\t{p}\t{idx}:L{j}")
                if cluster_id_col:
                    cluster_value = row[cluster_id_col]
                    parsed = parse_cluster_id(cluster_value)
                    if parsed:
                        chr_token, c_start, c_end, rest = parsed
                        cluster_meta[str(idx)] = (chr_token, rest)
                        chrom = normalize_chr(chr_token, chr_prefix)
                        bed_lines.append(f"{chrom}\t{c_start - 1}\t{c_end}\t{idx}:C")
        else:
            for idx, row in enumerate(reader):
                if not row:
                    continue
                if not index_meta:
                    row_len = len(row)
                    if chr_col_index is not None or chr_col_index_from_end is not None:
                        index_meta["chr"] = resolve_index(chr_col_index, chr_col_index_from_end, row_len)
                    if pos_col_index is not None or pos_col_index_from_end is not None:
                        index_meta["pos"] = resolve_index(pos_col_index, pos_col_index_from_end, row_len)
                    elif start_col_index is not None or start_col_index_from_end is not None:
                        index_meta["start"] = resolve_index(start_col_index, start_col_index_from_end, row_len)
                        index_meta["end"] = resolve_index(end_col_index, end_col_index_from_end, row_len)
                    if cluster_id_col_index is not None or cluster_id_col_index_from_end is not None:
                        index_meta["cluster_id"] = resolve_index(
                            cluster_id_col_index, cluster_id_col_index_from_end, row_len
                        )
                rows.append(row)
                raw_chr = row[index_meta["chr"]] if "chr" in index_meta else ""
                cluster_value = row[index_meta["cluster_id"]] if "cluster_id" in index_meta else ""
                parsed_cluster = parse_cluster_id(cluster_value) if cluster_value else None
                if not raw_chr and parsed_cluster:
                    raw_chr = parsed_cluster[0]
                if is_mt_chr(raw_chr, mt_chroms):
                    mt_row_ids.add(idx)
                    continue
                if raw_chr:
                    chrom = normalize_chr(raw_chr, chr_prefix)
                    if "pos" in index_meta:
                        pos = parse_pos(row[index_meta["pos"]])
                        bed_lines.append(f"{chrom}\t{pos - 1}\t{pos}\t{idx}")
                    elif "start" in index_meta and "end" in index_meta:
                        start = parse_pos(row[index_meta["start"]])
                        end = parse_pos(row[index_meta["end"]])
                        bed_lines.append(f"{chrom}\t{start - 1}\t{end}\t{idx}")
                if parsed_cluster:
                    chr_token, c_start, c_end, rest = parsed_cluster
                    cluster_meta[str(idx)] = (chr_token, rest)
                    chrom = normalize_chr(chr_token, chr_prefix)
                    bed_lines.append(f"{chrom}\t{c_start - 1}\t{c_end}\t{idx}:C")

    # 统计 bed 提交量：核心坐标 + reads列表 + cluster_id
    nuclear_bed_count = sum(1 for l in bed_lines if ":" not in l.split("\t")[3])
    reads_list_bed_count = sum(1 for l in bed_lines if ":L" in l.split("\t")[3])
    cluster_bed_count = sum(1 for l in bed_lines if l.split("\t")[3].endswith(":C"))

    if bed_lines:
        with open(bed_in, "w", encoding="utf-8") as f:
            f.write("\n".join(bed_lines))
            f.write("\n")

        cmd = [liftover_bin]
        if min_match is not None:
            cmd.append(f"-minMatch={min_match}")
        cmd.extend([bed_in, chain, bed_out, bed_unmapped])

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(
                "liftOver failed:\n"
                f"STDOUT: {result.stdout}\n"
                f"STDERR: {result.stderr}"
            )

        lifted_pos: Dict[str, Tuple[int, int]] = {}
        with open(bed_out, "r", encoding="utf-8") as f:
            for line in f:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 4:
                    continue
                row_id = fields[3]
                start = int(fields[1])
                end = int(fields[2])
                lifted_pos[row_id] = (start + 1, end)

        failed_ids = set()
        with open(bed_unmapped, "r", encoding="utf-8") as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 4:
                    continue
                failed_ids.add(fields[3])
    else:
        lifted_pos = {}
        failed_ids = set()

    nuclear_lifted = sum(1 for k in lifted_pos if ":" not in k)
    nuclear_failed = sum(1 for k in failed_ids if ":" not in k)
    reads_lifted = sum(1 for k in lifted_pos if ":L" in k)
    cluster_lifted = sum(1 for k in lifted_pos if k.endswith(":C"))

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    if has_header:
        with open(output_path, "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter=delimiter)
            writer.writeheader()
            for idx, row in enumerate(rows):
                row_key = str(idx)
                if idx in mt_row_ids:
                    writer.writerow(row)
                    continue
                if row_key in lifted_pos:
                    new_start, new_end = lifted_pos[row_key]
                    if pos_col:
                        row[pos_col] = str(new_start)
                    else:
                        row[start_col] = str(new_start)
                        row[end_col] = str(new_end)
                if reads_list_col and row_key in list_positions:
                    updated_list: List[int] = []
                    for j, orig in enumerate(list_positions[row_key]):
                        key = f"{row_key}:L{j}"
                        if key in lifted_pos:
                            updated_list.append(lifted_pos[key][0])
                        else:
                            updated_list.append(orig)
                    row[reads_list_col] = format_pos_list(updated_list)
                if cluster_id_col and row_key in cluster_meta:
                    key = f"{row_key}:C"
                    if key in lifted_pos:
                        new_start, new_end = lifted_pos[key]
                        chr_token, rest = cluster_meta[row_key]
                        row[cluster_id_col] = f"{chr_token}_{new_start}_{new_end}_{rest}"
                writer.writerow(row)
                continue
                if row_key in failed_ids:
                    if on_fail == "drop":
                        continue
                    if on_fail == "na":
                        if pos_col:
                            row[pos_col] = ""
                        else:
                            row[start_col] = ""
                            row[end_col] = ""
                    writer.writerow(row)
                    continue
                writer.writerow(row)
    else:
        with open(output_path, "w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f, delimiter=delimiter)
            for idx, row in enumerate(rows):
                row_key = str(idx)
                if idx in mt_row_ids:
                    writer.writerow(row)
                    continue
                if row_key in lifted_pos:
                    new_start, new_end = lifted_pos[row_key]
                    if "pos" in index_meta:
                        row[index_meta["pos"]] = str(new_start)
                    else:
                        row[index_meta["start"]] = str(new_start)
                        row[index_meta["end"]] = str(new_end)
                if "cluster_id" in index_meta and row_key in cluster_meta:
                    key = f"{row_key}:C"
                    if key in lifted_pos:
                        new_start, new_end = lifted_pos[key]
                        chr_token, rest = cluster_meta[row_key]
                        row[index_meta["cluster_id"]] = f"{chr_token}_{new_start}_{new_end}_{rest}"
                writer.writerow(row)
                continue
                if row_key in failed_ids:
                    if on_fail == "drop":
                        continue
                    if on_fail == "na":
                        if "pos" in index_meta:
                            row[index_meta["pos"]] = ""
                        else:
                            row[index_meta["start"]] = ""
                            row[index_meta["end"]] = ""
                    writer.writerow(row)
                    continue
                writer.writerow(row)

    for path in (bed_in, bed_out, bed_unmapped):
        if os.path.exists(path):
            os.remove(path)

    stats: Dict = {
        "format": "standard",
        "total_rows": len(rows),
        "mt_rows": len(mt_row_ids),
        "nuclear_rows": len(rows) - len(mt_row_ids),
        "nuclear_bed_submitted": nuclear_bed_count,
        "nuclear_lifted": nuclear_lifted,
        "nuclear_failed": nuclear_failed,
        "reads_list_positions_submitted": reads_list_bed_count,
        "reads_list_positions_lifted": reads_lifted,
        "cluster_ids_submitted": cluster_bed_count,
        "cluster_ids_lifted": cluster_lifted,
        "on_fail_action": on_fail,
        "rows_output": len(rows),
        "has_header": has_header,
        "has_reads_list": reads_list_col is not None,
        "has_cluster_id": (cluster_id_col is not None or "cluster_id" in index_meta),
    }
    return stats


def main() -> None:
    parser = argparse.ArgumentParser(description="LiftOver NUMTs breakpoint positions to GRCh38.")
    parser.add_argument("--input", required=True, help="Input TSV file")
    parser.add_argument("--output", required=True, help="Output TSV file")
    parser.add_argument("--config", required=True, help="YAML config path")
    args = parser.parse_args()
    run(args.input, args.output, args.config)


if __name__ == "__main__":
    main()
