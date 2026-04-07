#!/usr/bin/env python3
"""Run NUMTs enrichment permutation simulation.

Example:
    python enrichment_simulation.py --input-tsv converted.tsv --numts-total 100 \
        --numts-target 10 --numt-length 200 --output-prefix output --output-dir output
"""

from __future__ import annotations

import argparse
import logging
import os
import random
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


def run(
    input_tsv: str | Path,
    numts_total: int,
    numts_target: int,
    numt_length: int,
    output_prefix: str,
    output_dir: str | Path | None = None,
    simulation_runs: int = 10000,
) -> dict[str, str | float | int]:
    """Run permutation enrichment analysis and write detail and summary files."""
    try:
        import pyranges as pr
    except ImportError as exc:
        raise RuntimeError("未安装pyranges，请先安装后再运行。") from exc

    input_path = Path(input_tsv)
    if not input_path.is_file():
        raise FileNotFoundError(f"Input TSV file not found: {input_path}")
    if numts_total <= 0:
        raise ValueError("--numts-total must be a positive integer.")
    if numt_length <= 0:
        raise ValueError("--numt-length must be a positive integer.")
    if simulation_runs <= 1:
        raise ValueError("--simulation-runs must be greater than 1.")

    output_path = Path(output_dir) if output_dir is not None else input_path.parent
    output_path.mkdir(parents=True, exist_ok=True)

    log.info("使用PyRanges进行高效区间重叠计算")
    numts_target_perc = numts_target / numts_total
    df_ref = pd.read_csv(input_path, sep="\t")
    df_ref["chr"] = "chr"
    df_ref = df_ref[["chr", "refStart_noChr", "refEnd_noChr"]]
    df_ref_pr = pr.PyRanges(
        df_ref.rename(
            columns={"chr": "Chromosome", "refStart_noChr": "Start", "refEnd_noChr": "End"}
        )
    )

    freq_list = []
    for i in range(1, simulation_runs):
        random_start = random.sample(range(1, 3088269832), numts_total)
        random_numts = pd.DataFrame(random_start, columns=["randomStart"])
        random_numts["randomStart"] = random_numts["randomStart"].astype(int)
        random_numts["randomEnd"] = random_numts["randomStart"] + numt_length
        random_numts["chr"] = "chr"
        random_numts = random_numts[["chr", "randomStart", "randomEnd"]]

        random_numts_pr = pr.PyRanges(
            random_numts.rename(
                columns={"chr": "Chromosome", "randomStart": "Start", "randomEnd": "End"}
            )
        )
        overlap_count = len(random_numts_pr.join(df_ref_pr).df.drop_duplicates(subset=["Start", "End"]))
        freq_list.append(overlap_count / numts_total)

        if i % 100 == 0:
            log.info("  模拟进度: %d/%d (%.1f%%)", i, simulation_runs - 1, i * 100 / (simulation_runs - 1))

    df_pro = pd.DataFrame(freq_list, columns=["Percentage"])
    pvalue_less = 1 - (len(df_pro[df_pro["Percentage"] > numts_target_perc]) / simulation_runs)
    pvalue_greater = 1 - (len(df_pro[df_pro["Percentage"] < numts_target_perc]) / simulation_runs)
    df_pro.columns = [f"Pvalue_less({pvalue_less}) & Pvalue_greater({pvalue_greater})"]

    log.info("Pvalue(less): %s", pvalue_less)
    log.info("Pvalue(greater): %s", pvalue_greater)

    detail_file = output_path / f"{os.path.basename(input_path)}.enrichmentOutput.{output_prefix}.tsv"
    summary_file = output_path / f"enrichment_summary_{output_prefix}.tsv"
    df_pro.to_csv(detail_file, index=True, sep="\t")

    summary_df = pd.DataFrame(
        {
            "Parameter": [
                "Total_NUMTs",
                "Target_Region_NUMTs",
                "NUMT_Length_bp",
                "Simulation_Runs",
                "Observed_Percentage",
                "P_value_less",
                "P_value_greater",
                "Significance_less_than_random",
                "Significance_greater_than_random",
            ],
            "Value": [
                numts_total,
                numts_target,
                numt_length,
                simulation_runs - 1,
                f"{numts_target_perc:.4f}",
                f"{pvalue_less:.6f}",
                f"{pvalue_greater:.6f}",
                "Yes" if pvalue_less < 0.05 else "No",
                "Yes" if pvalue_greater < 0.05 else "No",
            ],
            "Description": [
                "Total number of NUMTs in the genome",
                "Number of NUMTs in target regions",
                "Average length of NUMTs",
                "Number of permutation simulations",
                "Observed percentage of NUMTs in target regions",
                "Probability that random distribution has fewer overlaps",
                "Probability that random distribution has more overlaps",
                "Is target region depleted of NUMTs?",
                "Is target region enriched for NUMTs?",
            ],
        }
    )
    summary_df.to_csv(summary_file, index=False, sep="\t")

    simulation_percentages = df_pro.iloc[:, 0].values
    log.info("富集分析完成！")
    log.info("  - 详细结果: %s", detail_file)
    log.info("  - 摘要结果: %s", summary_file)
    log.info("")
    log.info("=== 模拟结果统计 ===")
    log.info("随机模拟中重叠百分比统计:")
    log.info("  平均值: %.6f", np.mean(simulation_percentages))
    log.info("  标准差: %.6f", np.std(simulation_percentages))
    log.info("  最小值: %.6f", np.min(simulation_percentages))
    log.info("  最大值: %.6f", np.max(simulation_percentages))
    log.info("  观察值: %.6f", numts_target_perc)
    log.info("=====================")

    return {
        "detail_file": str(detail_file),
        "summary_file": str(summary_file),
        "pvalue_less": float(pvalue_less),
        "pvalue_greater": float(pvalue_greater),
        "simulation_runs": int(simulation_runs),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run NUMTs enrichment permutation simulation.")
    parser.add_argument("legacy_input_tsv", nargs="?", help="Legacy positional converted TSV input.")
    parser.add_argument("legacy_numts_total", nargs="?", type=int, help="Legacy positional total NUMTs.")
    parser.add_argument("legacy_numts_target", nargs="?", type=int, help="Legacy positional target NUMTs.")
    parser.add_argument("legacy_numt_length", nargs="?", type=int, help="Legacy positional NUMT length.")
    parser.add_argument("legacy_output_prefix", nargs="?", help="Legacy positional output prefix.")
    parser.add_argument("legacy_output_dir", nargs="?", help="Legacy positional output directory.")
    parser.add_argument("--input-tsv", help="Converted pseudo-reference TSV from create_pseudo_reference.py.")
    parser.add_argument("--numts-total", type=int, help="Total number of NUMTs.")
    parser.add_argument("--numts-target", type=int, help="Number of NUMTs overlapping target regions.")
    parser.add_argument("--numt-length", type=int, help="Average NUMT length in bp.")
    parser.add_argument("--output-prefix", help="Output prefix.")
    parser.add_argument("--output-dir", help="Output directory.")
    parser.add_argument("--simulation-runs", type=int, default=10000, help="Permutation loop upper bound.")
    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = build_parser().parse_args(argv)
    input_tsv = args.input_tsv or args.legacy_input_tsv
    numts_total = args.numts_total if args.numts_total is not None else args.legacy_numts_total
    numts_target = args.numts_target if args.numts_target is not None else args.legacy_numts_target
    numt_length = args.numt_length if args.numt_length is not None else args.legacy_numt_length
    output_prefix = args.output_prefix or args.legacy_output_prefix
    output_dir = args.output_dir or args.legacy_output_dir

    missing = [
        name
        for name, value in {
            "--input-tsv": input_tsv,
            "--numts-total": numts_total,
            "--numts-target": numts_target,
            "--numt-length": numt_length,
            "--output-prefix": output_prefix,
        }.items()
        if value is None
    ]
    if missing:
        raise SystemExit(f"Missing required argument(s): {', '.join(missing)}")

    run(input_tsv, numts_total, numts_target, numt_length, output_prefix, output_dir, args.simulation_runs)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
