#!/usr/bin/env python3
"""Write a bilingual Markdown report for NUMT breakpoint enrichment results."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)

REGION_LABELS = {
    "1-启动子区域.bed": "promoter regions",
    "2-蛋白编码区域.bed": "protein-coding regions",
    "3-大规模片段重复区域.bed": "segmental duplication regions",
    "4-去除片段重复区域的人类STR区域.bed": "STR regions after removing segmental duplications",
    "05-CpG_islands.bed": "CpG islands",
    "06-Microsat.bed": "microsatellites",
    "07-Start_codon.bed": "start codons",
    "08-snRNA.bed": "snRNA regions",
    "09-FuncElems.bed": "functional elements",
    "10-Retroposon.bed": "retroposons",
    "11-CDS.bed": "CDS regions",
    "12-Stop_codon.bed": "stop codons",
    "13-SINE.bed": "SINEs",
    "14-rmsk-DNA.bed": "RepeatMasker DNA elements",
    "15-Centromeres.bed": "centromeres",
    "16-Satellite.bed": "satellites",
    "17-LTR.bed": "LTR elements",
    "18-Exon.bed": "exons",
    "19-UTR.bed": "UTRs",
    "20-LINE.bed": "LINEs",
    "21-Simple_Repeats.bed": "simple repeats",
    "22-Genomic_superdups.bed": "genomic superduplications",
    "23-snoRNA_miRNA.bed": "snoRNA/miRNA regions",
    "24-Gene.bed": "gene regions",
    "25-Intron.bed": "introns",
    "26-Regulatory_elements.bed": "regulatory elements",
}


def english_region_name(region_name: str) -> str:
    return REGION_LABELS.get(str(region_name), str(region_name))


def format_empirical_p(p_value: float, simulation_runs: int) -> str:
    if p_value == 0:
        return f"<{1 / simulation_runs:.4g}"
    return f"{p_value:.4g}"


def format_result_sentence(row: pd.Series, english: bool = False) -> str:
    region = row["region_name"]
    cls = row["frequency_class"]
    observed = float(row["observed_percentage"]) * 100
    p_less = float(row["p_value_less"])
    p_greater = float(row["p_value_greater"])
    simulation_runs = int(row["simulation_runs"])
    if p_greater < 0.05:
        direction_cn = "显著富集"
        direction_en = "significantly enriched"
        p_value = p_greater
    elif p_less < 0.05:
        direction_cn = "显著减少"
        direction_en = "significantly depleted"
        p_value = p_less
    else:
        direction_cn = "未见显著偏离随机分布"
        direction_en = "not significantly different from random expectation"
        p_value = min(p_less, p_greater)
    formatted_p = format_empirical_p(p_value, simulation_runs)
    if english:
        return f"In {english_region_name(region)}, {cls} NUMT breakpoint flanks were {direction_en} (observed overlap {observed:.2f}%, empirical P={formatted_p})."
    return f"`{region}` 中 `{cls}` NUMT breakpoint flanks {direction_cn}（观测重叠比例 {observed:.2f}%，经验 P={formatted_p}）。"


def run(summary_tsv: str | Path, prepared_summary: str | Path, output_md: str | Path) -> dict[str, str | int]:
    summary_path = Path(summary_tsv)
    prepared_path = Path(prepared_summary)
    output_path = Path(output_md)
    if not summary_path.is_file():
        raise FileNotFoundError(f"Summary TSV not found: {summary_path}")
    if not prepared_path.is_file():
        raise FileNotFoundError(f"Prepared NUMT summary not found: {prepared_path}")

    df = pd.read_csv(summary_path, sep="\t")
    prepared = pd.read_csv(prepared_path, sep="\t")
    strata_cols = ["frequency_class", "breakpoint_count", "cluster_count"]
    significant = df.loc[(df["p_value_less"] < 0.05) | (df["p_value_greater"] < 0.05)].copy()
    significant["best_p"] = significant[["p_value_less", "p_value_greater"]].min(axis=1)
    significant = significant.sort_values(["best_p", "region_name", "frequency_class"])

    lines: list[str] = []
    lines.append("# NUMT nuclear breakpoint flank enrichment report")
    lines.append("")
    lines.append("## 中文摘要")
    lines.append("")
    lines.append(
        "本分析按照已有研究的富集思路，使用 confident nuclear breakpoints "
        "上下游各 100 bp 的 flanks 作为分析区间，并按 `all`、`common`、`low-frequency`、"
        "`rare` 和 `ultra-rare` 分层，对目标基因组区域进行全量 permutation 富集分析。"
    )
    lines.append("")
    lines.append("### Breakpoint flank 分层数量")
    lines.append("")
    lines.append(prepared[strata_cols].to_markdown(index=False))
    lines.append("")
    lines.append("### 主要结果")
    lines.append("")
    if significant.empty:
        lines.append("在经验 P < 0.05 阈值下，未观察到显著富集或显著减少的目标区域/NUMT 分层组合。")
    else:
        for _, row in significant.iterrows():
            lines.append(f"- {format_result_sentence(row, english=False)}")
    lines.append("")
    lines.append("简短结论：NUMT 核基因组断点 flanks 在不同功能区域和频率层级中呈现可检测的区域差异；显著结果应优先结合重复序列背景、可比对性和公共 NUMT 注释进一步解释。")
    lines.append("")
    lines.append("## English Summary")
    lines.append("")
    lines.append(
        "Following the nuclear-genome enrichment strategy in Figure 7D of the paper, "
        "this analysis used 100 bp flanks on both sides of confident nuclear NUMT breakpoints "
        "and performed full permutation-based enrichment testing across target genomic regions, "
        "stratified by `all`, `common`, `low-frequency`, `rare`, and `ultra-rare` NUMTs."
    )
    lines.append("")
    lines.append("### Breakpoint flank strata counts")
    lines.append("")
    lines.append(prepared[strata_cols].to_markdown(index=False))
    lines.append("")
    lines.append("### Main results")
    lines.append("")
    if significant.empty:
        lines.append("No target-region/frequency-class combination showed significant enrichment or depletion at empirical P < 0.05.")
    else:
        for _, row in significant.iterrows():
            lines.append(f"- {format_result_sentence(row, english=True)}")
    lines.append("")
    lines.append(
        "Concise conclusion: NUMT nuclear breakpoint flanks show region- and frequency-class-specific distributional patterns. "
        "Significant enrichment or depletion signals should be interpreted together with repeat context, mappability, "
        "and known NUMT annotations before assigning biological mechanism."
    )
    lines.append("")
    lines.append("## Full result table")
    lines.append("")
    table_cols = [
        "region_name",
        "frequency_class",
        "breakpoints_total",
        "breakpoints_target",
        "observed_percentage",
        "p_value_less",
        "p_value_greater",
    ]
    lines.append(df[table_cols].to_markdown(index=False))
    lines.append("")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines), encoding="utf-8")
    log.info("Wrote report: %s", output_path)
    return {"output_md": str(output_path), "significant_results": int(len(significant))}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Write bilingual Markdown enrichment report.")
    parser.add_argument("--summary-tsv", required=True)
    parser.add_argument("--prepared-summary", required=True)
    parser.add_argument("--output-md", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = build_parser().parse_args(argv)
    run(**vars(args))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
