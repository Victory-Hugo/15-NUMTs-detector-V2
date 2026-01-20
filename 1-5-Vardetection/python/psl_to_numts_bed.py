#!/usr/bin/env python3

import argparse
import sys
from typing import Dict, Tuple

MT_TARGETS = {"chrM", "MT", "NC_012920.1"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Parse PSL and write a 4-column numts.bed TSV for mtDNA best hits."
        )
    )
    parser.add_argument("psl", help="Input PSL file")
    parser.add_argument("out", help="Output numts.bed (TSV, 4 columns)")
    return parser.parse_args()


def is_header_line(line: str) -> bool:
    stripped = line.strip()
    if not stripped:
        return True
    return stripped.startswith("psLayout") or stripped.startswith("match")


def best_hit_key(matches: int, mismatches: int) -> Tuple[int, int]:
    # Higher matches is better; lower mismatches is better.
    return (matches, -mismatches)


def main() -> int:
    args = parse_args()

    best_hits: Dict[str, Tuple[int, int, str, int, int]] = {}
    # qName -> (matches, mismatches, tName, tStart, tEnd)

    with open(args.psl, "r", encoding="utf-8") as f:
        for raw in f:
            if is_header_line(raw):
                continue
            parts = raw.rstrip("\n").split()
            if len(parts) < 17:
                continue

            try:
                matches = int(parts[0])
                mismatches = int(parts[1])
                q_name = parts[9]
                t_name = parts[13]
                t_start = int(parts[15])
                t_end = int(parts[16])
            except (ValueError, IndexError):
                continue

            if t_name not in MT_TARGETS:
                continue

            current = best_hits.get(q_name)
            if current is None:
                best_hits[q_name] = (matches, mismatches, t_name, t_start, t_end)
                continue

            cur_matches, cur_mismatches, _, _, _ = current
            if best_hit_key(matches, mismatches) > best_hit_key(cur_matches, cur_mismatches):
                best_hits[q_name] = (matches, mismatches, t_name, t_start, t_end)

    with open(args.out, "w", encoding="utf-8") as out_f:
        for q_name, (_, _, _, t_start, t_end) in sorted(best_hits.items()):
            numt_start = t_start + 1
            numt_end = t_end
            out_f.write(f"{q_name}\t{numt_start}\t{numt_end}\tbest_hit_to_mtDNA\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
