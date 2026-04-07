#!/usr/bin/env python3
"""Prepare GENCODE-derived hg38 BED files for NUMT enrichment."""

from __future__ import annotations

import argparse
import logging
import re
import subprocess
from pathlib import Path

log = logging.getLogger(__name__)

PRIMARY_CHROMS = {f"chr{idx}" for idx in range(1, 23)} | {"chrX", "chrY"}
OUTPUTS = {
    "start_codon": "07-Start_codon.bed",
    "snRNA": "08-snRNA.bed",
    "CDS": "11-CDS.bed",
    "stop_codon": "12-Stop_codon.bed",
    "exon": "18-Exon.bed",
    "UTR": "19-UTR.bed",
    "snoRNA_miRNA": "23-snoRNA_miRNA.bed",
    "gene": "24-Gene.bed",
    "intron": "25-Intron.bed",
}


def parse_attributes(raw_attributes: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for match in re.finditer(r'([A-Za-z0-9_]+) "([^"]*)"', raw_attributes):
        attrs[match.group(1)] = match.group(2)
    return attrs


def open_text(path: Path):
    if path.suffix == ".gz":
        import gzip

        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def write_raw_beds(gencode_gtf_gz: Path, tmp_dir: Path) -> dict[str, Path]:
    handles = {name: (tmp_dir / f"{name}.gencode.raw.bed").open("w", encoding="utf-8") for name in OUTPUTS}
    try:
        with open_text(gencode_gtf_gz) as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue
                chrom, _, feature, start_raw, end_raw, _, _, _, attributes_raw = fields
                if chrom not in PRIMARY_CHROMS:
                    continue
                start = int(start_raw) - 1
                end = int(end_raw)
                if end <= start:
                    continue
                attrs = parse_attributes(attributes_raw)
                gene_type = attrs.get("gene_type") or attrs.get("gene_biotype") or attrs.get("transcript_type", "")
                bed_line = f"{chrom}\t{start}\t{end}\n"

                if feature in {"start_codon", "CDS", "stop_codon", "exon", "UTR", "gene"}:
                    handles[feature].write(bed_line)
                if feature == "gene" and gene_type == "snRNA":
                    handles["snRNA"].write(bed_line)
                if feature == "gene" and gene_type in {"snoRNA", "miRNA"}:
                    handles["snoRNA_miRNA"].write(bed_line)
    finally:
        for handle in handles.values():
            handle.close()
    return {name: tmp_dir / f"{name}.gencode.raw.bed" for name in OUTPUTS}


def sort_merge(input_bed: Path, output_bed: Path, bedtools_bin: str) -> None:
    sort_proc = subprocess.Popen(["sort", "-k1,1", "-k2,2n", str(input_bed)], stdout=subprocess.PIPE, text=True)
    with output_bed.open("w", encoding="utf-8") as out_handle:
        merge_proc = subprocess.run([bedtools_bin, "merge", "-i", "-"], stdin=sort_proc.stdout, stdout=out_handle, text=True)
    if sort_proc.stdout is not None:
        sort_proc.stdout.close()
    sort_return = sort_proc.wait()
    if sort_return != 0:
        raise RuntimeError(f"sort failed for {input_bed}")
    if merge_proc.returncode != 0:
        raise RuntimeError(f"bedtools merge failed for {input_bed}")


def build_introns(gene_bed: Path, exon_bed: Path, intron_bed: Path, bedtools_bin: str) -> None:
    subtract_proc = subprocess.Popen(
        [bedtools_bin, "subtract", "-a", str(gene_bed), "-b", str(exon_bed)],
        stdout=subprocess.PIPE,
        text=True,
    )
    sort_proc = subprocess.Popen(["sort", "-k1,1", "-k2,2n"], stdin=subtract_proc.stdout, stdout=subprocess.PIPE, text=True)
    if subtract_proc.stdout is not None:
        subtract_proc.stdout.close()
    with intron_bed.open("w", encoding="utf-8") as out_handle:
        merge_proc = subprocess.run([bedtools_bin, "merge", "-i", "-"], stdin=sort_proc.stdout, stdout=out_handle, text=True)
    if sort_proc.stdout is not None:
        sort_proc.stdout.close()
    if subtract_proc.wait() != 0:
        raise RuntimeError("bedtools subtract failed while building introns")
    if sort_proc.wait() != 0:
        raise RuntimeError("sort failed while building introns")
    if merge_proc.returncode != 0:
        raise RuntimeError("bedtools merge failed while building introns")


def run(gencode_gtf_gz: str | Path, tmp_dir: str | Path, output_dir: str | Path, bedtools_bin: str) -> dict[str, str]:
    gtf_path = Path(gencode_gtf_gz)
    tmp_path = Path(tmp_dir)
    output_path = Path(output_dir)
    if not gtf_path.is_file():
        raise FileNotFoundError(f"GENCODE GTF not found: {gtf_path}")
    tmp_path.mkdir(parents=True, exist_ok=True)
    output_path.mkdir(parents=True, exist_ok=True)

    raw_paths = write_raw_beds(gtf_path, tmp_path)
    for name, raw_path in raw_paths.items():
        if name == "intron":
            continue
        sort_merge(raw_path, output_path / OUTPUTS[name], bedtools_bin)
    build_introns(
        gene_bed=output_path / OUTPUTS["gene"],
        exon_bed=output_path / OUTPUTS["exon"],
        intron_bed=output_path / OUTPUTS["intron"],
        bedtools_bin=bedtools_bin,
    )
    log.info("Wrote GENCODE-derived BED files to %s", output_path)
    return {name: str(output_path / filename) for name, filename in OUTPUTS.items()}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prepare GENCODE-derived article region BED files.")
    parser.add_argument("--gencode-gtf-gz", required=True)
    parser.add_argument("--tmp-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--bedtools-bin", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = build_parser().parse_args(argv)
    run(**vars(args))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
