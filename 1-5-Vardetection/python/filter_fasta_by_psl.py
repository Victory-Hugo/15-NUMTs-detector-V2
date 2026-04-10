#!/usr/bin/env python3

import argparse
import sys

MT_TARGETS = {"chrM", "MT", "NC_012920.1"}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter FASTA records to keep only sequences with mtDNA hits in a PSL file."
    )
    parser.add_argument("--psl", required=True, help="Input PSL file (from BLAT vs mtDNA)")
    parser.add_argument("--fasta", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Output filtered FASTA file")
    return parser.parse_args()


def is_header_line(line: str) -> bool:
    stripped = line.strip()
    if not stripped:
        return True
    return stripped.startswith("psLayout") or stripped.startswith("match")


def load_hit_names(psl_path: str) -> set:
    hit_names = set()
    with open(psl_path, "r", encoding="utf-8") as f:
        for raw in f:
            if is_header_line(raw):
                continue
            parts = raw.rstrip("\n").split()
            if len(parts) < 17:
                continue
            try:
                t_name = parts[13]
                q_name = parts[9]
            except IndexError:
                continue
            if t_name in MT_TARGETS:
                hit_names.add(q_name)
    return hit_names


def filter_fasta(fasta_path: str, hit_names: set, output_path: str) -> int:
    kept = 0
    with open(fasta_path, "r", encoding="utf-8") as f_in, open(
        output_path, "w", encoding="utf-8"
    ) as f_out:
        header = None
        seq_parts = []
        current_key = None
        write_current = False

        def flush():
            nonlocal kept
            if header is None:
                return
            if write_current and seq_parts:
                seq = "".join(seq_parts)
                f_out.write(f">{header}\n")
                for i in range(0, len(seq), 60):
                    f_out.write(seq[i : i + 60] + "\n")
                kept += 1

        for line in f_in:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                header = line[1:]
                current_key = header.split()[0] if header else ""
                write_current = current_key in hit_names
                seq_parts = []
            else:
                if write_current:
                    seq_parts.append(line.strip())
        flush()

    return kept


def main() -> int:
    args = parse_args()

    try:
        hit_names = load_hit_names(args.psl)
    except OSError as e:
        print(f"error: cannot read PSL file: {e}", file=sys.stderr)
        return 1

    if not hit_names:
        print(
            f"warning: no MT hits found in PSL (targets: {sorted(MT_TARGETS)})",
            file=sys.stderr,
        )
        # Write empty output file
        try:
            open(args.output, "w").close()
        except OSError:
            pass
        return 2

    try:
        kept = filter_fasta(args.fasta, hit_names, args.output)
    except OSError as e:
        print(f"error: {e}", file=sys.stderr)
        return 1

    if kept == 0:
        print("warning: no sequences matched MT hits", file=sys.stderr)
        return 2

    return 0


if __name__ == "__main__":
    sys.exit(main())
