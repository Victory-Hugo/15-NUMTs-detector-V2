import logging
from pathlib import Path

import pandas as pd


log = logging.getLogger(__name__)

INPUT_NAME = "6-numt-mtdna-length-by-cluster.tsv"
OUTPUT_NAME = "Circos_tile.txt"
REQUIRED_COLUMNS = [
    "is_dloop_artifact",
    "cluster_min_mt_start",
    "cluster_max_mt_end",
]
DEFAULT_INPUT = (
    "/mnt/l/20-NUMTs/6-NUMTs频率分布描述/2-严格阈值-v1/output/02-tables/"
    + INPUT_NAME
)
DEFAULT_OUTPUT = (
    "/mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/4-可视化NUMTs分布/"
    "0-准备输入数据/output/"
    + OUTPUT_NAME
)


def _normalize_false_mask(series: pd.Series) -> pd.Series:
    text = series.astype(str).str.strip().str.lower()
    return text.eq("false")


def _build_hsm_region(start_value: object, end_value: object) -> str:
    start = int(float(start_value))
    end = int(float(end_value))
    return f"hsM {start} {end}"


def run(input_path: str, output_path: str) -> dict:
    input_file = Path(input_path)
    output_file = Path(output_path)

    if not input_file.exists():
        raise FileNotFoundError(f"输入文件不存在: {input_file}")

    df = pd.read_csv(input_file, sep="\t", dtype=object)

    missing_columns = [column for column in REQUIRED_COLUMNS if column not in df.columns]
    if missing_columns:
        raise ValueError(f"缺少必要列: {', '.join(missing_columns)}")

    filtered = df.loc[_normalize_false_mask(df["is_dloop_artifact"])].copy()
    filtered = filtered.dropna(subset=["cluster_min_mt_start", "cluster_max_mt_end"]).copy()

    filtered["cluster_min_mt_start"] = pd.to_numeric(
        filtered["cluster_min_mt_start"], errors="coerce"
    )
    filtered["cluster_max_mt_end"] = pd.to_numeric(
        filtered["cluster_max_mt_end"], errors="coerce"
    )
    filtered = filtered.dropna(subset=["cluster_min_mt_start", "cluster_max_mt_end"]).copy()

    filtered["mt_region_hsm"] = filtered.apply(
        lambda row: _build_hsm_region(
            row["cluster_min_mt_start"], row["cluster_max_mt_end"]
        ),
        axis=1,
    )

    output_file.parent.mkdir(parents=True, exist_ok=True)
    filtered["mt_region_hsm"].to_csv(output_file, index=False, header=False)

    log.info("Wrote %s rows to %s", len(filtered), output_file)
    return {
        "input_path": str(input_file),
        "output_path": str(output_file),
        "row_count": len(filtered),
    }


def build_parser():
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "保留 is_dloop_artifact == FALSE 的行，去除 mt 起止为空的行，"
            "输出无表头的 hsM 起始 终止 单列文本。"
        )
    )
    parser.add_argument("--input-path", default=DEFAULT_INPUT, help="输入 TSV 文件路径")
    parser.add_argument("--output-path", default=DEFAULT_OUTPUT, help="输出文本文件路径")
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="日志级别",
    )
    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    result = run(input_path=args.input_path, output_path=args.output_path)
    logging.info(
        "处理完成: input=%s output=%s rows=%s",
        Path(result["input_path"]),
        Path(result["output_path"]),
        result["row_count"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
