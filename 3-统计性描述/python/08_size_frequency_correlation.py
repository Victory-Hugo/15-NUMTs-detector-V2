#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""计算 NUMT 长度与群体频率的统计相关性（Pearson + Spearman），输出摘要 TSV。"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import List, Optional

import pandas as pd
from scipy import stats

LOG = logging.getLogger(__name__)

_FILTER_DESC = "is_dloop_artifact==False, dropna(numt_length_bp, frequency)"


def configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="计算 NUMT 长度与群体频率的 Pearson 和 Spearman 相关性，输出摘要 TSV。"
    )
    parser.add_argument(
        "--input-tsv",
        required=True,
        help="6-numt-mtdna-length-by-cluster.tsv（step 06 输出）",
    )
    parser.add_argument(
        "--output-tsv",
        required=True,
        help="相关性统计摘要输出路径（TSV）",
    )
    return parser


def load_and_filter(input_tsv: str) -> pd.DataFrame:
    LOG.info("加载输入文件: %s", input_tsv)
    df = pd.read_csv(input_tsv, sep="\t")
    df["numt_length_bp"] = pd.to_numeric(df["numt_length_bp"], errors="coerce")
    df["frequency"] = pd.to_numeric(df["frequency"], errors="coerce")
    df = df[df["is_dloop_artifact"] == False].dropna(
        subset=["numt_length_bp", "frequency"]
    )
    LOG.info("过滤后剩余 %d 行", len(df))
    return df


def compute_correlations(df: pd.DataFrame) -> pd.DataFrame:
    n = len(df)
    x = df["numt_length_bp"].values
    y = df["frequency"].values

    pearson_r, pearson_p = stats.pearsonr(x, y)
    spearman_r, spearman_p = stats.spearmanr(x, y)

    LOG.info(
        "Pearson  r=%.4f  p=%.3e | Spearman rho=%.4f  p=%.3e  (n=%d)",
        pearson_r,
        pearson_p,
        spearman_r,
        spearman_p,
        n,
    )

    rows = [
        {
            "test": "pearson",
            "variable_x": "numt_length_bp",
            "variable_y": "frequency",
            "n": n,
            "statistic": round(float(pearson_r), 6),
            "pvalue": float(pearson_p),
            "significant_alpha_0.05": bool(pearson_p < 0.05),
            "filter": _FILTER_DESC,
        },
        {
            "test": "spearman",
            "variable_x": "numt_length_bp",
            "variable_y": "frequency",
            "n": n,
            "statistic": round(float(spearman_r), 6),
            "pvalue": float(spearman_p),
            "significant_alpha_0.05": bool(spearman_p < 0.05),
            "filter": _FILTER_DESC,
        },
    ]
    return pd.DataFrame(rows)


def run(input_tsv: str, output_tsv: str) -> int:
    df = load_and_filter(input_tsv)
    result = compute_correlations(df)
    Path(output_tsv).parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_tsv, sep="\t", index=False)
    LOG.info("相关性摘要已保存: %s", output_tsv)
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    configure_logging()
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(input_tsv=args.input_tsv, output_tsv=args.output_tsv)


if __name__ == "__main__":
    raise SystemExit(main())
