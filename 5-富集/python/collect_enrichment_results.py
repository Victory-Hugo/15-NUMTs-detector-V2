#!/usr/bin/env python3
"""Collect per-task enrichment summaries into one table."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)


def run(task_manifest: str | Path, output_tsv: str | Path) -> dict[str, str | int]:
    manifest_path = Path(task_manifest)
    output_path = Path(output_tsv)
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Task manifest not found: {manifest_path}")

    manifest = pd.read_csv(manifest_path, sep="\t")
    frames = []
    missing = []
    for row in manifest.itertuples(index=False):
        summary_path = Path(row.output_dir) / "summary.tsv"
        if summary_path.is_file():
            frames.append(pd.read_csv(summary_path, sep="\t"))
        else:
            missing.append(str(summary_path))
    if missing:
        raise FileNotFoundError("Missing task summaries:\n" + "\n".join(missing))
    if not frames:
        raise ValueError(f"No task summaries listed in {manifest_path}")

    combined = pd.concat(frames, ignore_index=True)
    combined = combined.sort_values(["region_name", "frequency_class"]).reset_index(drop=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_path, sep="\t", index=False)
    log.info("Collected %d task summaries into %s", len(combined), output_path)
    return {"output_tsv": str(output_path), "rows": int(len(combined))}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Collect per-task enrichment summaries.")
    parser.add_argument("--task-manifest", required=True)
    parser.add_argument("--output-tsv", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = build_parser().parse_args(argv)
    run(**vars(args))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
