#!/usr/bin/env python3
"""旧版 NUMTs 聚类脚本的 Python 3 兼容实现。"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd


LOG = logging.getLogger(__name__)

SAM_COLUMNS = [
    "QNAME",
    "FLAG",
    "RNAME",
    "POS",
    "MAPQ",
    "CIGAR",
    "RNEXT",
    "PNEXT",
    "TLEN",
    "SEQ",
    "QUAL",
    "SM",
    "RG",
    "NM",
    "BC",
    "OC",
    "ZX",
    "ZY",
    "SA",
]

CLUSTER_SUMMARY_COLUMNS = [
    "IndividualID",
    "Cluster_ID",
    "Cluster_No",
    "subCluster_No",
    "size",
]

BREAKPOINT_COLUMNS = [
    "IndividualID",
    "Cluster_No",
    "disFile",
    "splitFile",
    "wgsBAM",
    "chr",
    "start",
    "end",
]


def build_parser() -> argparse.ArgumentParser:
    """构建命令行参数解析器。"""
    parser = argparse.ArgumentParser(
        description="从 mtDNA discordant SAM 文件中提取旧版 NUMTs 聚类结果。"
    )
    parser.add_argument("legacy_args", nargs="*", help="旧版位置参数: sampleID wgsBAM input_disc")
    parser.add_argument("--sample-id", dest="sample_id", help="样本 ID")
    parser.add_argument("--wgs-bam", dest="wgs_bam", help="原始 WGS BAM 文件路径")
    parser.add_argument("--input-disc", dest="input_disc", help="discordant SAM 文件路径")
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="日志级别",
    )
    return parser


def configure_logging(log_level: str) -> None:
    """初始化日志。"""
    logging.basicConfig(
        level=getattr(logging, log_level),
        format="[%(levelname)s] %(message)s",
    )


def resolve_arguments(args: argparse.Namespace) -> tuple[str, str, str]:
    """统一解析命名参数和旧版位置参数。"""
    sample_id = args.sample_id
    wgs_bam = args.wgs_bam
    input_disc = args.input_disc

    if args.legacy_args:
        if len(args.legacy_args) != 3:
            raise ValueError("旧版调用方式要求 3 个位置参数: sampleID wgsBAM input_disc")
        legacy_sample_id, legacy_wgs_bam, legacy_input_disc = args.legacy_args
        sample_id = sample_id or legacy_sample_id
        wgs_bam = wgs_bam or legacy_wgs_bam
        input_disc = input_disc or legacy_input_disc

    if not sample_id or not wgs_bam or not input_disc:
        raise ValueError("缺少必要参数: sample_id, wgs_bam, input_disc")

    return sample_id, wgs_bam, input_disc


def cluster(data: Sequence[int], maxgap: int) -> list[list[int]]:
    """按旧版规则对位置列表进行聚类。"""
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


def read_discordant_table(input_disc: str) -> pd.DataFrame:
    """读取 discordant SAM 表。"""
    rows: list[list[str | None]] = []
    with open(input_disc, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            if not raw_line or raw_line.startswith("@"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < len(SAM_COLUMNS):
                fields.extend([None] * (len(SAM_COLUMNS) - len(fields)))
            rows.append(fields[: len(SAM_COLUMNS)])

    if not rows:
        return pd.DataFrame(columns=SAM_COLUMNS)

    table = pd.DataFrame(rows, columns=SAM_COLUMNS)

    table["RNAME"] = table["RNAME"].fillna("")
    table["RNEXT"] = table["RNEXT"].fillna("")
    table["POS"] = pd.to_numeric(table["POS"], errors="coerce")
    table["MAPQ"] = pd.to_numeric(table["MAPQ"], errors="coerce")
    table["PNEXT"] = pd.to_numeric(table["PNEXT"], errors="coerce")
    table = table.dropna(subset=["POS", "MAPQ", "PNEXT"])
    table["POS"] = table["POS"].astype(int)
    table["MAPQ"] = table["MAPQ"].astype(int)
    table["PNEXT"] = table["PNEXT"].astype(int)
    return table


def build_empty_outputs(input_disc: str) -> None:
    """生成空结果文件。"""
    pd.DataFrame(columns=SAM_COLUMNS + ["subCluster_No", "Cluster_No", "Cluster_ID", "IndividualID", "clusterID"]).to_csv(
        input_disc + ".cluster.tsv",
        sep="\t",
        header=True,
        index=False,
    )
    pd.DataFrame(columns=CLUSTER_SUMMARY_COLUMNS).to_csv(
        input_disc + ".cluster.summary.tsv",
        sep="\t",
        header=False,
        index=True,
    )
    Path(input_disc + ".breakpointINPUT.tsv").write_text("", encoding="utf-8")


def ensure_summary_frame() -> pd.DataFrame:
    """创建空的 summary 表。"""
    return pd.DataFrame(columns=CLUSTER_SUMMARY_COLUMNS)


def build_breakpoint_input(
    summary_frame: pd.DataFrame,
    input_disc: str,
    wgs_bam: str,
) -> pd.DataFrame:
    """根据旧版逻辑生成 breakpointINPUT.tsv。"""
    if summary_frame.empty:
        return pd.DataFrame(columns=BREAKPOINT_COLUMNS)

    output2 = summary_frame.copy()
    output2["disFile"] = input_disc
    output2["splitFile"] = input_disc.replace("disc", "split")
    output2["wgsBAM"] = wgs_bam

    cluster_split = output2["Cluster_ID"].str.split("_", expand=True)
    cluster_split.columns = ["chr", "start", "end", "chrM", "mtstart", "mtend"]

    output3 = pd.concat(
        [
            output2.drop(columns=["Cluster_ID", "subCluster_No", "size"]),
            cluster_split[["chr", "start", "end"]],
        ],
        axis=1,
    )
    output3 = output3.drop_duplicates(["chr", "start", "end"]).copy()
    output3["start"] = pd.to_numeric(output3["start"], errors="coerce").fillna(0).astype(int) - 500
    output3["end"] = pd.to_numeric(output3["end"], errors="coerce").fillna(0).astype(int) + 500 + 150
    output3["Cluster_No"] = pd.to_numeric(output3["Cluster_No"], errors="coerce").fillna(0).astype(int)
    return output3[BREAKPOINT_COLUMNS]


def run(sample_id: str, wgs_bam: str, input_disc: str) -> int:
    """执行旧版聚类逻辑。"""
    input_path = Path(input_disc)
    if not input_path.exists():
        raise FileNotFoundError(f"输入文件不存在: {input_disc}")
    if not Path(wgs_bam).exists():
        raise FileNotFoundError(f"WGS BAM 文件不存在: {wgs_bam}")

    normalized_sample_id = sample_id.replace(".mt.disc", "")
    df0 = read_discordant_table(input_disc)
    LOG.info("读取到 %s 条 discordant 记录", len(df0))

    if df0.empty:
        LOG.warning("输入 discordant 文件为空，输出空结果")
        build_empty_outputs(input_disc)
        return 0

    df = df0[~df0["RNAME"].str.contains(r"Un_|random|\.", regex=True)].copy()
    df1 = df[
        (df["MAPQ"] > 0)
        & (
            df["RNAME"].isin(["MT", "chrM"])
            | df["RNEXT"].isin(["MT", "chrM"])
        )
    ].copy()
    df3 = df1[~df1["RNAME"].isin(["chrM", "MT"])].copy()
    df3 = df3.sort_values(["RNAME", "POS"]).copy()

    if df3.empty:
        LOG.warning("未检测到满足旧版规则的核基因组 discordant 记录")
        build_empty_outputs(input_disc)
        return 0

    output_frames: list[pd.DataFrame] = []

    for cluster_id, myclusters in df3.groupby("RNAME"):
        pos_clusters = cluster(myclusters["POS"].astype(int).tolist(), maxgap=500)
        for pos_values in pos_clusters:
            if len(pos_values) < 2:
                continue
            df_cluster = df3[df3["POS"].isin(pos_values)].copy()
            df_cluster_pair_mt = df[df["QNAME"].isin(df_cluster["QNAME"])].copy()
            mt_clusters = cluster(df_cluster["PNEXT"].astype(int).tolist(), maxgap=500)
            subcluster_frames: list[pd.DataFrame] = []
            for mt_values in mt_clusters:
                df_cluster_pair_mt_out = df_cluster_pair_mt[
                    df_cluster_pair_mt["PNEXT"].isin(mt_values)
                ].copy()
                if df_cluster_pair_mt_out.empty:
                    continue
                df_cluster_pair_mt_out["subCluster_No"] = len(mt_values)
                df_cluster_pair_mt_out["Cluster_No"] = len(pos_values)
                df_cluster_pair_mt_out["Cluster_ID"] = (
                    f"{cluster_id}_{int(df_cluster['POS'].min())}_{int(df_cluster['POS'].max())}"
                    f"_MTboth_{int(df_cluster_pair_mt_out['PNEXT'].min())}_{int(df_cluster_pair_mt_out['PNEXT'].max())}"
                )
                subcluster_frames.append(df_cluster_pair_mt_out)
            if subcluster_frames:
                output_frames.append(pd.concat(subcluster_frames, ignore_index=True))

    if output_frames:
        output1 = pd.concat(output_frames, ignore_index=True)
        output1["IndividualID"] = normalized_sample_id
        output1["clusterID"] = output1["RNAME"].astype(str) + "_" + output1["Cluster_No"].astype(str)
        output2 = (
            output1.groupby(["IndividualID", "Cluster_ID", "Cluster_No", "subCluster_No"])
            .size()
            .to_frame("size")
            .reset_index()
        )
    else:
        output1 = pd.DataFrame(columns=SAM_COLUMNS + ["subCluster_No", "Cluster_No", "Cluster_ID", "IndividualID", "clusterID"])
        output2 = ensure_summary_frame()

    output1.to_csv(input_disc + ".cluster.tsv", sep="\t", header=True, index=False)
    output2.to_csv(input_disc + ".cluster.summary.tsv", sep="\t", header=False, index=True)
    build_breakpoint_input(output2, input_disc, wgs_bam).to_csv(
        input_disc + ".breakpointINPUT.tsv",
        sep="\t",
        header=False,
        index=False,
    )
    LOG.info("聚类结果已写出: %s.cluster.tsv", input_disc)
    return 0


def main(argv: Iterable[str] | None = None) -> int:
    """命令行入口。"""
    parser = build_parser()
    args = parser.parse_args(argv)
    configure_logging(args.log_level)
    try:
        sample_id, wgs_bam, input_disc = resolve_arguments(args)
        return run(sample_id=sample_id, wgs_bam=wgs_bam, input_disc=input_disc)
    except Exception as exc:  # pragma: no cover - CLI 兜底
        LOG.error("%s", exc)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
