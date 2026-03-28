#!/usr/bin/env python3
"""旧版断点搜索脚本的 Python 3 兼容实现。"""

from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path
from typing import Iterable

import pandas as pd


LOG = logging.getLogger(__name__)

PSL_COLUMNS = [
    "match",
    "misMatch",
    "repMatch",
    "Ns",
    "QgapCount",
    "QgapBases",
    "TgapCount",
    "TgapBases",
    "strand",
    "Qname",
    "Qsize",
    "Qstart",
    "Qend",
    "Tname",
    "Tsize",
    "Tstart",
    "Tend",
    "blockCount",
    "blockSizes",
    "qStarts",
    "tStarts",
]

OUTPUT_COLUMNS = ["pointGroup", "Group", "chr", "Tend", "strand", "readsCount", "Tstart", "sampleID"]


def build_parser() -> argparse.ArgumentParser:
    """构建命令行参数解析器。"""
    parser = argparse.ArgumentParser(description="根据旧版 PSL 结果搜索 NUMTs 断点。")
    parser.add_argument("legacy_args", nargs="*", help="旧版位置参数: input_psl sample_id chr start end output_prefix")
    parser.add_argument("--input-psl", dest="input_psl", help="输入 PSL 文件")
    parser.add_argument("--sample-id", dest="sample_id", help="样本 ID")
    parser.add_argument("--chromosome", dest="chromosome", help="候选区域染色体")
    parser.add_argument("--start", dest="start", type=int, help="候选区域起点")
    parser.add_argument("--end", dest="end", type=int, help="候选区域终点")
    parser.add_argument("--output-prefix", dest="output_prefix", help="输出文件前缀")
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


def resolve_arguments(args: argparse.Namespace) -> tuple[str, str, str, int, int, str]:
    """统一解析命名参数和旧版位置参数。"""
    input_psl = args.input_psl
    sample_id = args.sample_id
    chromosome = args.chromosome
    start = args.start
    end = args.end
    output_prefix = args.output_prefix

    if args.legacy_args:
        if len(args.legacy_args) != 6:
            raise ValueError("旧版调用方式要求 6 个位置参数: input_psl sample_id chr start end output_prefix")
        legacy_input_psl, legacy_sample_id, legacy_chr, legacy_start, legacy_end, legacy_output_prefix = args.legacy_args
        input_psl = input_psl or legacy_input_psl
        sample_id = sample_id or legacy_sample_id
        chromosome = chromosome or legacy_chr
        start = start if start is not None else int(legacy_start)
        end = end if end is not None else int(legacy_end)
        output_prefix = output_prefix or legacy_output_prefix

    if not input_psl or not sample_id or not chromosome or start is None or end is None or not output_prefix:
        raise ValueError("缺少必要参数: input_psl, sample_id, chromosome, start, end, output_prefix")

    return input_psl, sample_id, chromosome, int(start), int(end), output_prefix


def f_nu(row: pd.Series) -> str:
    """按旧版规则标记核基因组断点方向。"""
    mismatch_len = 3
    read_len = 150
    if row["strand"] == "+" and row["Qend"] >= (read_len - mismatch_len):
        return "nu_Tstart_Bright"
    if row["strand"] == "+" and row["Qstart"] <= mismatch_len:
        return "nu_Tend_Bleft"
    if row["strand"] == "-":
        return "nu_NegStrand"
    return "nu_useLess"


def f_mt(row: pd.Series) -> str:
    """按旧版规则标记线粒体断点方向。"""
    mismatch_len = 3
    read_len = 150
    if row["strand"] == "+" and row["Qend"] >= (read_len - mismatch_len):
        return "mt_Tstart"
    if row["strand"] == "+" and row["Qstart"] <= mismatch_len:
        return "mt_Tend"
    if row["strand"] == "-" and row["Qend"] >= (read_len - mismatch_len):
        return "mt_Tend"
    if row["strand"] == "-" and row["Qstart"] <= mismatch_len:
        return "mt_Tstart"
    return "mt_useLess"


def read_psl_table(input_psl: str) -> pd.DataFrame:
    """读取 PSL 文件。"""
    if Path(input_psl).stat().st_size == 0:
        return pd.DataFrame(columns=PSL_COLUMNS)

    try:
        table = pd.read_csv(
            input_psl,
            skiprows=5,
            sep="\t",
            names=PSL_COLUMNS,
            dtype=str,
        )
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=PSL_COLUMNS)

    if table.empty:
        return table

    numeric_columns = [
        "match",
        "misMatch",
        "Qstart",
        "Qend",
        "Tstart",
        "Tend",
    ]
    for column in numeric_columns:
        table[column] = pd.to_numeric(table[column], errors="coerce")
    table = table.dropna(subset=numeric_columns)
    return table


def normalize_target_name(target_name: str) -> str:
    """把 RefSeq accession 名称映射为旧版脚本可识别的染色体名称。"""
    if pd.isna(target_name):
        return ""

    text = str(target_name).strip()
    if text in {"chrM", "MT", "chrX", "X", "chrY", "Y"}:
        return text

    match = re.match(r"^NC_0*([0-9]+)\.\d+$", text)
    if not match:
        return text

    chrom_number = int(match.group(1))
    if 1 <= chrom_number <= 22:
        return f"chr{chrom_number}"
    if chrom_number == 23:
        return "chrX"
    if chrom_number == 24:
        return "chrY"
    if chrom_number == 12920:
        return "chrM"
    return text


def group_counts(table: pd.DataFrame, group_columns: list[str]) -> pd.DataFrame:
    """执行旧版 groupby 统计。"""
    if table.empty:
        output_columns = group_columns + ["readsCount"]
        return pd.DataFrame(columns=output_columns)
    return table.groupby(group_columns).size().reset_index(name="readsCount")


def empty_output_frame() -> pd.DataFrame:
    """生成空输出表。"""
    return pd.DataFrame(columns=OUTPUT_COLUMNS)


def run(input_psl: str, sample_id: str, chromosome: str, start: int, end: int, output_prefix: str) -> int:
    """执行旧版断点搜索逻辑。"""
    if not Path(input_psl).exists():
        raise FileNotFoundError(f"输入 PSL 文件不存在: {input_psl}")

    mismatch_len = 3
    read_len = 150
    hg37_chr = chromosome.replace("chr", "")

    df_input = read_psl_table(input_psl)
    if df_input.empty:
        LOG.warning("PSL 文件没有可用记录，输出空结果")
        empty_output_frame().to_csv(output_prefix + ".Breakpoints.tsv", sep="\t", header=False, index=False)
        return 0

    df_input["matchLEN"] = df_input["Tend"] - df_input["Tstart"]
    df_input["Tname_normalized"] = df_input["Tname"].map(normalize_target_name)
    df2 = df_input[
        (df_input["matchLEN"] < (read_len - 10))
        & (df_input["misMatch"] <= mismatch_len)
    ].copy()

    df_mt = df2[df2["Tname_normalized"].isin(["chrM", "MT"])].copy()
    df_nu = df2[
        (df2["Tname_normalized"].isin([chromosome, hg37_chr]))
        & (df2["Tstart"] > start)
        & (df2["Tend"] < end)
    ].copy()

    df_nu["pointGroup"] = df_nu.apply(f_nu, axis=1) if not df_nu.empty else pd.Series(dtype=str)
    df_mt["pointGroup"] = df_mt.apply(f_mt, axis=1) if not df_mt.empty else pd.Series(dtype=str)
    df_nu["chr"] = chromosome
    df_nu["Group"] = df_nu["pointGroup"].str.replace(r"_T.*B", "", regex=True)
    df_mt["chr"] = "chrM"

    df_nu_left = df_nu[df_nu["pointGroup"] == "nu_Tend_Bleft"].copy()
    df_nu_left_g = group_counts(df_nu_left, ["pointGroup", "Group", "chr", "Tend", "strand"])
    df_nu_right = df_nu[df_nu["pointGroup"] == "nu_Tstart_Bright"].copy()
    df_nu_right_g = group_counts(df_nu_right, ["pointGroup", "Group", "chr", "Tstart", "strand"])

    df_nu_mt = df_nu[df_nu["Qname"].isin(df_mt["Qname"].tolist())].copy()
    df_nu_left_mt = df_nu_mt[df_nu_mt["pointGroup"] == "nu_Tend_Bleft"].copy()
    df_nu_left_g_mt = group_counts(df_nu_left_mt, ["pointGroup", "Group", "chr", "Tend", "strand"])
    df_nu_right_mt = df_nu_mt[df_nu_mt["pointGroup"] == "nu_Tstart_Bright"].copy()
    df_nu_right_g_mt = group_counts(df_nu_right_mt, ["pointGroup", "Group", "chr", "Tstart", "strand"])

    df_mt_tend = df_mt[df_mt["pointGroup"] == "mt_Tend"].copy()
    df_mt_tend_g = group_counts(df_mt_tend, ["pointGroup", "chr", "Tend", "strand"])
    df_mt_tstart = df_mt[df_mt["pointGroup"] == "mt_Tstart"].copy()
    df_mt_tstart_g = group_counts(df_mt_tstart, ["pointGroup", "chr", "Tstart", "strand"])

    df_mt_tend_left = df_mt_tend[df_mt_tend["Qname"].isin(df_nu_left["Qname"].tolist())].copy()
    df_mt_tend_right = df_mt_tend[df_mt_tend["Qname"].isin(df_nu_right["Qname"].tolist())].copy()
    df_mt_tstart_left = df_mt_tstart[df_mt_tstart["Qname"].isin(df_nu_left["Qname"].tolist())].copy()
    df_mt_tstart_right = df_mt_tstart[df_mt_tstart["Qname"].isin(df_nu_right["Qname"].tolist())].copy()

    df_mt_left_tstart_g = group_counts(df_mt_tstart_left, ["pointGroup", "chr", "Tstart", "strand"])
    df_mt_left_tstart_g["Group"] = "mtLeft"
    df_mt_left_tend_g = group_counts(df_mt_tend_left, ["pointGroup", "chr", "Tend", "strand"])
    df_mt_left_tend_g["Group"] = "mtLeft"
    df_mt_right_tstart_g = group_counts(df_mt_tstart_right, ["pointGroup", "chr", "Tstart", "strand"])
    df_mt_right_tstart_g["Group"] = "mtRight"
    df_mt_right_tend_g = group_counts(df_mt_tend_right, ["pointGroup", "chr", "Tend", "strand"])
    df_mt_right_tend_g["Group"] = "mtRight"

    df_numt_conf = pd.concat(
        [
            df_nu_left_g_mt,
            df_nu_right_g_mt,
            df_mt_left_tstart_g,
            df_mt_left_tend_g,
            df_mt_right_tstart_g,
            df_mt_right_tend_g,
        ],
        ignore_index=True,
        sort=False,
    )

    if df_numt_conf.empty:
        output = empty_output_frame()
    else:
        output = df_numt_conf.copy()
        if "Tstart" not in output.columns:
            output["Tstart"] = pd.NA
        if "Tend" not in output.columns:
            output["Tend"] = pd.NA
        output["sampleID"] = sample_id
        output["Tstart"] = pd.to_numeric(output["Tstart"], errors="coerce").fillna(-1).astype(int)
        output["Tend"] = pd.to_numeric(output["Tend"], errors="coerce").fillna(-1).astype(int)
        output = output.reindex(columns=OUTPUT_COLUMNS)

    output.to_csv(output_prefix + ".Breakpoints.tsv", sep="\t", header=False, index=False)
    LOG.info("断点结果已写出: %s.Breakpoints.tsv", output_prefix)
    return 0


def main(argv: Iterable[str] | None = None) -> int:
    """命令行入口。"""
    parser = build_parser()
    args = parser.parse_args(argv)
    configure_logging(args.log_level)
    try:
        resolved = resolve_arguments(args)
        return run(*resolved)
    except Exception as exc:  # pragma: no cover - CLI 兜底
        LOG.error("%s", exc)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
