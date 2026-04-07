#!/usr/bin/env python3
"""Create pseudo-reference coordinates for target BED regions.

Example:
    python create_pseudo_reference.py --input-bed genes.bed --output-tsv converted.tsv
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)

CHR_OFFSETS = {
    "chr1": 0,
    "chr2": 248956422,
    "chr3": 491149951,
    "chr4": 689445510,
    "chr5": 879660065,
    "chr6": 1061198324,
    "chr7": 1232004303,
    "chr8": 1391350276,
    "chr9": 1536488912,
    "chr10": 1674883629,
    "chr11": 1808681051,
    "chr12": 1943767673,
    "chr13": 2077042982,
    "chr14": 2191407310,
    "chr15": 2298451028,
    "chr16": 2400442217,
    "chr17": 2490780562,
    "chr18": 2574038003,
    "chr19": 2654411288,
    "chr20": 2713028904,
    "chr21": 2777473071,
    "chr22": 2824183054,
    "chrX": 2875001522,
    "chrY": 3031042417,
}


def run(input_bed: str | Path, output_tsv: str | Path | None = None) -> dict[str, str | int]:
    """Convert a BED file to pseudo-reference coordinates."""
    input_path = Path(input_bed)
    if output_tsv is None:
        output_path = Path(f"{input_path}.convert.noGap.tsv")
    else:
        output_path = Path(output_tsv)

    if not input_path.is_file():
        raise FileNotFoundError(f"Input BED file not found: {input_path}")

    df_raw = pd.read_csv(input_path, names=["chr", "start", "end"], sep="\t", usecols=[0, 1, 2])
    df = df_raw[df_raw["chr"].isin(CHR_OFFSETS)].copy()
    offsets = df["chr"].map(CHR_OFFSETS)
    df["refStart_noChr"] = df["start"] + offsets
    df["refEnd_noChr"] = df["end"] + offsets

    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False, sep="\t")
    log.info("转换完成，输出文件: %s", output_path)
    log.info("处理了 %d 个区域", len(df))
    return {"output_tsv": str(output_path), "regions": int(len(df))}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Create pseudo-reference coordinates from a BED file.")
    parser.add_argument("legacy_input_bed", nargs="?", help="Legacy positional input BED path.")
    parser.add_argument("legacy_output_tsv", nargs="?", help="Legacy positional output TSV path.")
    parser.add_argument("--input-bed", help="Input BED file with chr/start/end columns.")
    parser.add_argument("--output-tsv", help="Output TSV path. Defaults to <input>.convert.noGap.tsv.")
    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = build_parser().parse_args(argv)
    input_bed = args.input_bed or args.legacy_input_bed
    output_tsv = args.output_tsv or args.legacy_output_tsv
    if input_bed is None:
        raise SystemExit("Missing input BED. Use --input-bed <file> or the legacy positional argument.")
    run(input_bed=input_bed, output_tsv=output_tsv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
