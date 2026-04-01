#!/usr/bin/env python3
import argparse
import csv
import math
import os
import statistics
import subprocess
import uuid
from collections import Counter, defaultdict
from typing import Dict, Iterable, List, Optional, Tuple


OUTPUT_COLUMNS = [
    "Chr",
    "Mt Start",
    "Mt End",
    "Mt fragment length",
    "Nuc Start",
    "Nuc End",
    "Chr fragment length",
    "Difference",
    "Source",
]

VALID_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y"}


class UnionFind:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))
        self.rank = [0] * size

    def find(self, node: int) -> int:
        while self.parent[node] != node:
            self.parent[node] = self.parent[self.parent[node]]
            node = self.parent[node]
        return node

    def union(self, left: int, right: int) -> None:
        root_left = self.find(left)
        root_right = self.find(right)
        if root_left == root_right:
            return
        if self.rank[root_left] < self.rank[root_right]:
            root_left, root_right = root_right, root_left
        self.parent[root_right] = root_left
        if self.rank[root_left] == self.rank[root_right]:
            self.rank[root_left] += 1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Clean, liftover, and deduplicate reference NUMTs.")
    parser.add_argument("--input", required=True, help="Input TSV path")
    parser.add_argument("--output", required=True, help="Output TSV path")
    parser.add_argument("--summary", required=True, help="Summary report path")
    parser.add_argument("--tmp-dir", required=True, help="Temporary directory")
    parser.add_argument("--liftover-bin", required=True, help="liftOver executable")
    parser.add_argument("--hg18-chain", required=True, help="hg18 -> hg38 chain")
    parser.add_argument("--hg19-chain", required=True, help="hg19 -> hg38 chain")
    parser.add_argument("--hs1-chain", required=True, help="hs1 -> hg38 chain")
    parser.add_argument("--min-match", type=float, default=0.95, help="liftOver -minMatch value")
    parser.add_argument("--max-distance", type=int, default=1000, help="Dedup distance threshold")
    parser.add_argument(
        "--include-source",
        action="append",
        default=[],
        help="Restrict processing to listed Source values. May be passed multiple times.",
    )
    parser.add_argument(
        "--ignore-source-filter",
        action="store_true",
        help="Ignore --include-source values and process all sources in the input file.",
    )
    parser.add_argument(
        "--dedup-mode",
        choices=["event", "nuc_only"],
        default="event",
        help="Dedup rule: event requires both nuc and mt breakpoints; nuc_only uses only nuclear breakpoints.",
    )
    return parser.parse_args()


def parse_int(value: str) -> Optional[int]:
    text = (value or "").strip()
    if not text:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def parse_positive_int(value: str) -> Optional[int]:
    parsed = parse_int(value)
    if parsed is None or parsed <= 0:
        return None
    return parsed


def format_int(value: Optional[int]) -> str:
    return "" if value is None else str(value)


def normalize_chr(raw_chr: str) -> Optional[str]:
    text = (raw_chr or "").strip()
    if not text:
        return None
    if text.lower().startswith("chr"):
        text = text[3:]
    text = text.upper()
    if text == "M" or text == "MT":
        return None
    return text if text in VALID_CHROMS else None


def normalize_reference(raw_value: str) -> Optional[str]:
    text = (raw_value or "").strip().lower()
    if text == "hg18":
        return "hg18"
    if text == "hg19" or text == "grch37":
        return "hg19"
    if text == "hg38" or text == "grch38":
        return "hg38"
    if text == "t2t-chr13":
        return "hs1"
    return None


def resolve_interval(
    start: Optional[int],
    end: Optional[int],
    length: Optional[int],
) -> Tuple[Optional[int], Optional[int], Optional[int]]:
    if start is not None and end is not None:
        if end < start:
            start, end = end, start
    elif start is not None and length is not None:
        end = start + length - 1
    elif end is not None and length is not None:
        start = end - length + 1
        if start <= 0:
            start = None
            end = None
    elif start is not None:
        end = start
    elif end is not None:
        start = end

    if start is not None and end is not None:
        length = end - start + 1
    elif length is not None and length <= 0:
        length = None

    return start, end, length


def rounded_median(values: Iterable[int]) -> int:
    median_value = statistics.median(values)
    return int(math.floor(float(median_value) + 0.5))


def read_and_clean_rows(
    input_path: str,
    include_sources: Optional[set[str]],
    stats: Counter,
) -> List[Dict[str, Optional[int]]]:
    rows: List[Dict[str, Optional[int]]] = []

    with open(input_path, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row_number, raw_row in enumerate(reader, start=2):
            stats["input_rows"] += 1

            source = (raw_row.get("Source", "") or "").strip()
            if include_sources and source not in include_sources:
                stats["drop_excluded_source"] += 1
                continue

            ref_genome = normalize_reference(raw_row.get("Ref_Genome", ""))
            if ref_genome is None:
                stats["drop_unknown_ref"] += 1
                continue

            chrom = normalize_chr(raw_row.get("Chr", ""))
            if chrom is None:
                stats["drop_invalid_chr"] += 1
                continue

            nuc_start = parse_positive_int(raw_row.get("Nuc Start", ""))
            nuc_end = parse_positive_int(raw_row.get("Nuc End", ""))
            nuc_len = parse_positive_int(raw_row.get("Chr fragment length", ""))
            nuc_start, nuc_end, nuc_len = resolve_interval(nuc_start, nuc_end, nuc_len)
            if nuc_start is None or nuc_end is None or nuc_len is None:
                stats["drop_invalid_nuclear_interval"] += 1
                continue

            mt_start = parse_positive_int(raw_row.get("Mt Start", ""))
            mt_end = parse_positive_int(raw_row.get("Mt End", ""))
            mt_len = parse_positive_int(raw_row.get("Mt fragment length", ""))
            mt_start, mt_end, mt_len = resolve_interval(mt_start, mt_end, mt_len)

            if mt_start is None or mt_end is None:
                mt_start = None
                mt_end = None
                if mt_len is None or mt_len <= 0:
                    mt_len = None
                stats["rows_missing_mt_breakpoint"] += 1

            difference = abs(nuc_len - mt_len) if nuc_len is not None and mt_len is not None else None

            cleaned = {
                "row_number": row_number,
                "ref_genome": ref_genome,
                "Chr": chrom,
                "Mt Start": mt_start,
                "Mt End": mt_end,
                "Mt fragment length": mt_len,
                "Nuc Start": nuc_start,
                "Nuc End": nuc_end,
                "Chr fragment length": nuc_len,
                "Difference": difference,
                "Source": source,
                "liftover_status": "not_needed" if ref_genome == "hg38" else "pending",
            }
            rows.append(cleaned)
            stats[f"rows_ref_{ref_genome}"] += 1

    stats["clean_rows"] = len(rows)
    return rows


def write_bed(path: str, entries: List[Tuple[str, str, int, int]]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        for row_id, chrom, start, end in entries:
            handle.write(f"chr{chrom}\t{start - 1}\t{end}\t{row_id}\n")


def run_liftover_group(
    rows: List[Dict[str, Optional[int]]],
    source_ref: str,
    chain_path: str,
    liftover_bin: str,
    min_match: float,
    tmp_dir: str,
    stats: Counter,
) -> None:
    targets = [row for row in rows if row["ref_genome"] == source_ref]
    if not targets:
        return
    if not os.path.exists(chain_path):
        raise FileNotFoundError(f"Missing chain file for {source_ref}: {chain_path}")

    run_id = uuid.uuid4().hex
    bed_in = os.path.join(tmp_dir, f"{source_ref}_{run_id}.in.bed")
    bed_out = os.path.join(tmp_dir, f"{source_ref}_{run_id}.out.bed")
    bed_unmapped = os.path.join(tmp_dir, f"{source_ref}_{run_id}.unmapped.bed")

    entries = []
    for index, row in enumerate(targets):
        entries.append((f"{source_ref}:{index}", row["Chr"], row["Nuc Start"], row["Nuc End"]))
    write_bed(bed_in, entries)

    cmd = [liftover_bin, f"-minMatch={min_match}", bed_in, chain_path, bed_out, bed_unmapped]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"liftOver failed for {source_ref}\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )

    lifted: Dict[str, List[Tuple[str, int, int]]] = defaultdict(list)
    if os.path.exists(bed_out):
        with open(bed_out, "r", encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 4:
                    continue
                target_chr = normalize_chr(fields[0])
                if target_chr is None:
                    continue
                lifted[fields[3]].append((target_chr, int(fields[1]) + 1, int(fields[2])))

    unmapped = set()
    if os.path.exists(bed_unmapped):
        with open(bed_unmapped, "r", encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 4:
                    continue
                unmapped.add(fields[3])

    for index, row in enumerate(targets):
        row_id = f"{source_ref}:{index}"
        if row_id in unmapped:
            row["liftover_status"] = "unmapped"
            stats[f"liftover_{source_ref}_unmapped"] += 1
            continue
        mapped_segments = lifted.get(row_id, [])
        if len(mapped_segments) != 1:
            row["liftover_status"] = "ambiguous"
            stats[f"liftover_{source_ref}_ambiguous"] += 1
            continue
        mapped_chr, mapped_start, mapped_end = mapped_segments[0]
        row["Chr"] = mapped_chr
        row["Nuc Start"] = mapped_start
        row["Nuc End"] = mapped_end
        row["Chr fragment length"] = mapped_end - mapped_start + 1
        row["Difference"] = (
            abs(row["Chr fragment length"] - row["Mt fragment length"])
            if row["Mt fragment length"] is not None
            else None
        )
        row["ref_genome"] = "hg38"
        row["liftover_status"] = "success"
        stats[f"liftover_{source_ref}_success"] += 1

    for path in (bed_in, bed_out, bed_unmapped):
        if os.path.exists(path):
            os.remove(path)


def deduplicate_rows(
    rows: List[Dict[str, Optional[int]]],
    max_distance: int,
    dedup_mode: str,
    stats: Counter,
) -> List[Dict[str, Optional[int]]]:
    eligible = []
    passthrough = []
    for row in rows:
        if row["ref_genome"] != "hg38":
            stats["drop_not_converted_to_hg38"] += 1
            continue
        if dedup_mode == "event" and (row["Mt Start"] is None or row["Mt End"] is None):
            passthrough.append(row)
            stats["passthrough_missing_mt"] += 1
            continue
        eligible.append(row)

    stats["eligible_for_dedup"] = len(eligible)

    merged_rows: List[Dict[str, Optional[int]]] = []
    by_chr: Dict[str, List[Dict[str, Optional[int]]]] = defaultdict(list)
    for row in eligible:
        by_chr[row["Chr"]].append(row)

    for chrom, chrom_rows in by_chr.items():
        chrom_rows.sort(
            key=lambda item: (
                item["Nuc Start"],
                item["Nuc End"],
                item["Mt Start"] if item["Mt Start"] is not None else 10 ** 12,
                item["Mt End"] if item["Mt End"] is not None else 10 ** 12,
                item["Source"],
            )
        )
        uf = UnionFind(len(chrom_rows))
        right = 0
        for left in range(len(chrom_rows)):
            while right < len(chrom_rows) and chrom_rows[right]["Nuc Start"] - chrom_rows[left]["Nuc Start"] <= max_distance:
                right += 1
            for other in range(left + 1, right):
                first = chrom_rows[left]
                second = chrom_rows[other]
                if abs(first["Nuc End"] - second["Nuc End"]) > max_distance:
                    continue
                if dedup_mode == "event":
                    if abs(first["Mt Start"] - second["Mt Start"]) > max_distance:
                        continue
                    if abs(first["Mt End"] - second["Mt End"]) > max_distance:
                        continue
                uf.union(left, other)

        groups: Dict[int, List[Dict[str, Optional[int]]]] = defaultdict(list)
        for index, row in enumerate(chrom_rows):
            groups[uf.find(index)].append(row)

        for group_rows in groups.values():
            if len(group_rows) == 1:
                merged_rows.append(group_rows[0])
                continue

            stats["dedup_merged_events"] += 1
            stats["dedup_removed_rows"] += len(group_rows) - 1
            merged = build_merged_row(group_rows, dedup_mode)
            merged_rows.append(merged)

    merged_rows.extend(passthrough)
    return merged_rows


def build_merged_row(group_rows: List[Dict[str, Optional[int]]], dedup_mode: str) -> Dict[str, Optional[int]]:
    chrom = group_rows[0]["Chr"]
    nuc_start = rounded_median(row["Nuc Start"] for row in group_rows)
    nuc_end = rounded_median(row["Nuc End"] for row in group_rows)

    mt_start_values = [row["Mt Start"] for row in group_rows if row["Mt Start"] is not None]
    mt_end_values = [row["Mt End"] for row in group_rows if row["Mt End"] is not None]

    mt_start: Optional[int]
    mt_end: Optional[int]
    mt_length: Optional[int]
    if mt_start_values and mt_end_values:
        mt_start = rounded_median(mt_start_values)
        mt_end = rounded_median(mt_end_values)
        if mt_end < mt_start:
            mt_start, mt_end = mt_end, mt_start
        mt_length = mt_end - mt_start + 1
    else:
        mt_start = None
        mt_end = None
        mt_length = None

    if nuc_end < nuc_start:
        nuc_start, nuc_end = nuc_end, nuc_start

    nuc_length = nuc_end - nuc_start + 1
    sources = sorted({row["Source"] for row in group_rows if row["Source"]})

    return {
        "Chr": chrom,
        "Mt Start": mt_start,
        "Mt End": mt_end,
        "Mt fragment length": mt_length,
        "Nuc Start": nuc_start,
        "Nuc End": nuc_end,
        "Chr fragment length": nuc_length,
        "Difference": abs(nuc_length - mt_length) if mt_length is not None else None,
        "Source": ";".join(sources),
        "ref_genome": "hg38",
        "liftover_status": f"merged_{dedup_mode}",
    }


def sort_key(row: Dict[str, Optional[int]]) -> Tuple[int, int, int, int, int, str]:
    chrom = row["Chr"]
    chrom_rank = 23 if chrom == "X" else 24 if chrom == "Y" else int(chrom)
    mt_start = row["Mt Start"] if row["Mt Start"] is not None else 10 ** 12
    mt_end = row["Mt End"] if row["Mt End"] is not None else 10 ** 12
    return (
        chrom_rank,
        row["Nuc Start"],
        row["Nuc End"],
        mt_start,
        mt_end,
        row["Source"],
    )


def write_output(rows: List[Dict[str, Optional[int]]], output_path: str) -> None:
    final_rows = sorted(rows, key=sort_key)
    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in final_rows:
            output_row = {
                "Chr": row["Chr"],
                "Mt Start": format_int(row["Mt Start"]),
                "Mt End": format_int(row["Mt End"]),
                "Mt fragment length": format_int(row["Mt fragment length"]),
                "Nuc Start": format_int(row["Nuc Start"]),
                "Nuc End": format_int(row["Nuc End"]),
                "Chr fragment length": format_int(row["Chr fragment length"]),
                "Difference": format_int(row["Difference"]),
                "Source": row["Source"],
            }
            writer.writerow(output_row)


def write_summary(summary_path: str, stats: Counter, output_rows: int) -> None:
    lines = ["metric\tvalue"]
    for key in sorted(stats):
        lines.append(f"{key}\t{stats[key]}")
    lines.append(f"output_rows\t{output_rows}")
    with open(summary_path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))
        handle.write("\n")


def ensure_dependencies(args: argparse.Namespace) -> None:
    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input file not found: {args.input}")
    if not os.path.exists(args.liftover_bin):
        raise FileNotFoundError(f"liftOver binary not found: {args.liftover_bin}")
    os.makedirs(args.tmp_dir, exist_ok=True)
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    os.makedirs(os.path.dirname(args.summary), exist_ok=True)


def main() -> None:
    args = parse_args()
    ensure_dependencies(args)

    stats: Counter = Counter()
    include_sources = None if args.ignore_source_filter else (set(args.include_source) if args.include_source else None)
    rows = read_and_clean_rows(args.input, include_sources, stats)
    run_liftover_group(rows, "hg18", args.hg18_chain, args.liftover_bin, args.min_match, args.tmp_dir, stats)
    run_liftover_group(rows, "hg19", args.hg19_chain, args.liftover_bin, args.min_match, args.tmp_dir, stats)
    run_liftover_group(rows, "hs1", args.hs1_chain, args.liftover_bin, args.min_match, args.tmp_dir, stats)

    final_rows = deduplicate_rows(rows, args.max_distance, args.dedup_mode, stats)
    write_output(final_rows, args.output)
    write_summary(args.summary, stats, len(final_rows))

    print(f"Output written to: {args.output}")
    print(f"Summary written to: {args.summary}")
    print(f"Final rows: {len(final_rows)}")


if __name__ == "__main__":
    main()
