#!/usr/bin/env python3
"""旧版多样本 NUMTs 聚类脚本的 Python 3 兼容实现。"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd


LOG = logging.getLogger(__name__)

INPUT_COLUMNS = ["sampleID", "Cluster_No", "disFile", "splitFile", "wgsBAM", "chr", "start", "end"]
DETAIL_COLUMNS = ["GroupID", "Index", "POS", "CHR", "mergedClusterID"]
SUMMARY_COLUMNS = ["mergedClusterID", "GroupID", "Index", "POS", "CHR", "min", "max"]


def build_parser() -> argparse.ArgumentParser:
    """构建命令行参数解析器。"""
    parser = argparse.ArgumentParser(description="对多个样本的旧版 breakpointINPUT.tsv 进行聚类汇总。")
    parser.add_argument("legacy_args", nargs="*", help="旧版位置参数: input_file")
    parser.add_argument("--input-file", dest="input_file", help="合并后的 breakpointINPUT.tsv 文件")
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="日志级别",
    )
    return parser


def configure_logging(log_level: str) -> None:
    """初始化日志。"""
    logging.basicConfig(level=getattr(logging, log_level), format="[%(levelname)s] %(message)s")


def resolve_arguments(args: argparse.Namespace) -> str:
    """统一解析命名参数和旧版位置参数。"""
    input_file = args.input_file
    if args.legacy_args:
        if len(args.legacy_args) != 1:
            raise ValueError("旧版调用方式要求 1 个位置参数: input_file")
        input_file = input_file or args.legacy_args[0]
    if not input_file:
        raise ValueError("缺少必要参数: input_file")
    return input_file


def cluster(data: Sequence[int], maxgap: int) -> list[list[int]]:
    """按旧版规则对位置列表聚类。"""
    if not data:
        return []
    sorted_data = sorted(int(value) for value in data)
    groups = [[sorted_data[0]]]
    for value in sorted_data[1:]:
        if abs(value - groups[-1][-1]) <= maxgap:
            groups[-1].append(value)
        else:
            groups.append([value])
    return groups


def write_empty_outputs(input_file: str) -> None:
    """写出空聚类结果。"""
    pd.DataFrame(columns=DETAIL_COLUMNS).to_csv(
        input_file + ".allCluster.tsv",
        sep="\t",
        header=True,
        index=False,
    )
    pd.DataFrame(columns=SUMMARY_COLUMNS).to_csv(
        input_file + ".allCluster.sum.tsv",
        sep="\t",
        header=True,
        index=True,
    )


def run(input_file: str) -> int:
    """执行旧版多样本聚类逻辑。"""
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"输入文件不存在: {input_file}")

    if input_path.stat().st_size == 0:
        LOG.warning("输入文件为空，输出空结果")
        write_empty_outputs(input_file)
        return 0

    df0 = pd.read_csv(
        input_file,
        sep="\t",
        header=None,
        names=INPUT_COLUMNS,
        dtype=str,
        engine="python",
    )
    if df0.empty:
        LOG.warning("输入文件为空表，输出空结果")
        write_empty_outputs(input_file)
        return 0

    df = df0.drop_duplicates(["sampleID", "chr", "start"]).copy()
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["start", "end"]).copy()

    if df.empty:
        LOG.warning("输入文件缺少可用坐标，输出空结果")
        write_empty_outputs(input_file)
        return 0

    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["position"] = (df["start"] + df["end"]) / 2
    df = df.sort_values(["chr", "position"]).copy()
    df = df.drop_duplicates(["chr", "position"]).copy()

    detail_frames: list[pd.DataFrame] = []
    for cluster_id, myclusters in df.groupby("chr"):
        positions = myclusters["position"].astype(int).tolist()
        sub_clusters = cluster(positions, maxgap=1000)
        for group_id, sub_cluster in enumerate(sub_clusters):
            group_frame = pd.DataFrame(
                {
                    "GroupID": [group_id] * len(sub_cluster),
                    "Index": list(range(len(sub_cluster))),
                    "POS": sub_cluster,
                    "CHR": [cluster_id] * len(sub_cluster),
                }
            )
            detail_frames.append(group_frame)

    if detail_frames:
        output1 = pd.concat(detail_frames, ignore_index=True)
        output1["mergedClusterID"] = output1["GroupID"].astype(str) + "_" + output1["CHR"].astype(str)
        output2_1 = output1.groupby(["mergedClusterID"]).count().reset_index()
        output2_2 = output1.groupby(["mergedClusterID"]).POS.agg(["min", "max"]).reset_index()
        output2 = pd.merge(output2_1, output2_2, how="outer", on=["mergedClusterID"])
    else:
        output1 = pd.DataFrame(columns=DETAIL_COLUMNS)
        output2 = pd.DataFrame(columns=SUMMARY_COLUMNS)

    output1.to_csv(input_file + ".allCluster.tsv", sep="\t", header=True, index=False)
    output2.to_csv(input_file + ".allCluster.sum.tsv", sep="\t", header=True, index=True)
    LOG.info("多样本聚类结果已写出: %s.allCluster.tsv", input_file)
    return 0


def main(argv: Iterable[str] | None = None) -> int:
    """命令行入口。"""
    parser = build_parser()
    args = parser.parse_args(argv)
    configure_logging(args.log_level)
    try:
        input_file = resolve_arguments(args)
        return run(input_file=input_file)
    except Exception as exc:  # pragma: no cover - CLI 兜底
        LOG.error("%s", exc)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
