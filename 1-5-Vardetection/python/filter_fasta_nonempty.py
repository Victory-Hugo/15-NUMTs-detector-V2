#!/usr/bin/env python3

import argparse
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter FASTA records; keep sequences with >= min A/C/G/T bases."
    )
    parser.add_argument("--input", required=True, help="Input FASTA")
    parser.add_argument("--output", required=True, help="Output FASTA")
    parser.add_argument("--min-acgt", type=int, default=1, help="Minimum A/C/G/T count")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    kept = 0

    with open(args.input, "r", encoding="utf-8") as f_in, open(
        args.output, "w", encoding="utf-8"
    ) as f_out:
        header = None
        seq_parts = []

        def flush():
            nonlocal kept
            if header is None:
                return
            seq = "".join(seq_parts).strip()
            if not seq:
                return
            seq_upper = seq.upper()
            acgt_count = sum(seq_upper.count(b) for b in "ACGT")
            if acgt_count < args.min_acgt:
                return
            f_out.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                f_out.write(seq[i : i + 60] + "\n")
            kept += 1

        for line in f_in:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        flush()

    return 0 if kept > 0 else 2


if __name__ == "__main__":
    sys.exit(main())
