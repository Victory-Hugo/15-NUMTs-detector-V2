#!/usr/bin/env python3
"""Generate mtDNA target BED files from MitoGenes.txt."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)

MT_CATEGORY_MAP: dict[str, set[str]] = {
    "Protein_coding": {"ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
                       "COX1", "COX2", "COX3", "ATP6", "ATP8", "CYTB"},
    "rRNA": {"RNR1", "RNR2"},
    "Dloop": {"D-loop-A", "D-loop-B"},
    # tRNA: gene_name.startswith("TRN") — handled programmatically
}
COARSE_CATEGORY_ORDER = ["Protein_coding", "tRNA", "rRNA", "Dloop"]


def load_mito_genes(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chr", "start", "end", "gene_name"],
        dtype={"chr": str, "gene_name": str},
    )
    df["chr"] = df["chr"].str.replace("hsM", "chrM", regex=False)
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["start", "end"]).copy()
    df["start"] = df["start"].astype(int) - 1  # 1-based → 0-based
    df["end"] = df["end"].astype(int)
    df = df.loc[df["end"] > df["start"]].copy()
    return df


def classify_gene(gene_name: str) -> str | None:
    for category, gene_set in MT_CATEGORY_MAP.items():
        if gene_name in gene_set:
            return category
    if gene_name.startswith("TRN"):
        return "tRNA"
    return None


def write_beds_by_group(
    df: pd.DataFrame,
    output_dir: str | Path,
    group_column: str,
    output_order: list[str] | None = None,
) -> dict[str, Path]:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    result: dict[str, Path] = {}
    groups = output_order if output_order is not None else sorted(df[group_column].dropna().astype(str).unique().tolist())
    for group_name in groups:
        cat_df = df.loc[df[group_column] == group_name, ["chr", "start", "end", "gene_name"]].copy()
        cat_df = cat_df.sort_values(["chr", "start"]).reset_index(drop=True)
        bed_path = output_path / f"mt_{group_name}.bed"
        cat_df.to_csv(bed_path, sep="\t", header=False, index=False)
        result[group_name] = bed_path
        log.info("Wrote %d rows to %s", len(cat_df), bed_path)
    return result


def write_category_beds(
    mito_genes_txt: str | Path,
    output_dir: str | Path,
    fine_output_dir: str | Path | None = None,
) -> dict[str, dict[str, Path]]:
    df = load_mito_genes(mito_genes_txt)
    df["category"] = df["gene_name"].apply(classify_gene)
    unclassified = df.loc[df["category"].isna(), "gene_name"].tolist()
    if unclassified:
        log.warning("Unclassified genes (skipped): %s", unclassified)
    df = df.dropna(subset=["category"]).copy()

    outputs = {
        "coarse": write_beds_by_group(
            df=df,
            output_dir=output_dir,
            group_column="category",
            output_order=COARSE_CATEGORY_ORDER,
        )
    }

    if fine_output_dir is not None:
        fine_df = df.copy()
        fine_df["fine_category"] = fine_df["gene_name"].astype(str)
        outputs["fine"] = write_beds_by_group(
            df=fine_df,
            output_dir=fine_output_dir,
            group_column="fine_category",
        )

    return outputs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generate mtDNA target category BED files.")
    parser.add_argument("--mito-genes-txt", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--fine-output-dir", default=None)
    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = build_parser().parse_args(argv)
    write_category_beds(args.mito_genes_txt, args.output_dir, args.fine_output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
