#!/usr/bin/env python3
"""Filter NUMTs TSV.GZ tables to a selected sample ID set.

The module is importable via run(...) and executable as a CLI. It streams all
large TSV.GZ files line by line; only the ID list is held in memory.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import gzip
import logging
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Set

LOG = logging.getLogger(__name__)


EXPECTED_FILES = (
    "1-all_breakpoints.tsv.gz",
    "2-breakpoints_old.tsv.gz",
    "3-confident_breakpoints.tsv.gz",
    "4-mt_disc_breakpoint_input.tsv.gz",
    "5-mt_disc_cluster_old.tsv.gz",
    "6-mt_disc_cluster_summary.tsv.gz",
    "7-mt_disc_cluster.tsv.gz",
    "8-merge_bed.tsv.gz",
)

OPTIONAL_FILES = (
    "9-merge_bed_confident.tsv.gz",
)

REGION_QUERY_RE = re.compile(r"^(?P<sample>.+)_chr[^.]+\.\d+\.\d+(?:\||$)")


@dataclass
class FileStats:
    input_path: Path
    output_path: Path
    total_rows: int = 0
    kept_rows: int = 0
    dropped_rows: int = 0
    has_header: bool = False


def open_gzip_text(path: Path, mode: str):
    kwargs = {"encoding": "utf-8", "newline": ""}
    if "w" in mode:
        kwargs["compresslevel"] = 1
    return gzip.open(path, mode, **kwargs)


def split_tsv(line: str) -> List[str]:
    return line.rstrip("\n").split("\t")


def build_header_index(header_line: str) -> Dict[str, int]:
    return {name: idx for idx, name in enumerate(split_tsv(header_line))}


def load_sample_ids(id_file: Path) -> Set[str]:
    sample_ids: Set[str] = set()
    with id_file.open("r", encoding="utf-8") as handle:
        for line_no, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split()
            if not fields:
                continue
            sample_id = fields[0]
            if not sample_id:
                raise ValueError(f"Empty sample ID at {id_file}:{line_no}")
            sample_ids.add(sample_id)
    if not sample_ids:
        raise ValueError(f"No sample IDs found in {id_file}")
    return sample_ids


def ensure_outputs_can_be_written(output_dir: Path, overwrite: bool) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    known_outputs = EXPECTED_FILES + OPTIONAL_FILES
    existing = [output_dir / name for name in known_outputs if (output_dir / name).exists()]
    if existing and not overwrite:
        existing_text = "\n".join(str(path) for path in existing)
        raise FileExistsError("Output exists; use --overwrite to replace:\n" + existing_text)


def validate_input_dir(input_dir: Path) -> Dict[str, Path]:
    if not input_dir.is_dir():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    paths: Dict[str, Path] = {}
    missing: List[str] = []
    for name in EXPECTED_FILES:
        path = input_dir / name
        if path.exists():
            paths[name] = path
        else:
            missing.append(str(path))
    if missing:
        raise FileNotFoundError("Missing required inputs:\n" + "\n".join(missing))

    for name in OPTIONAL_FILES:
        path = input_dir / name
        if path.exists():
            paths[name] = path
        else:
            LOG.warning("Optional input missing; skipping: %s", path)

    known = set(EXPECTED_FILES + OPTIONAL_FILES)
    unexpected = sorted(path.name for path in input_dir.glob("*.tsv.gz") if path.name not in known)
    if unexpected:
        raise ValueError("Unexpected TSV.GZ files in input directory:\n" + "\n".join(unexpected))
    return paths


def get_index(header_idx: Dict[str, int], path: Path, candidates: Iterable[str]) -> int:
    for name in candidates:
        if name in header_idx:
            return header_idx[name]
    raise ValueError(f"Missing any of columns {tuple(candidates)} in {path}")


def sample_from_fixed_column(parts: List[str], column_index: int, path: Path, line_no: int) -> str:
    try:
        return parts[column_index]
    except IndexError as exc:
        raise ValueError(f"Expected at least {column_index + 1} columns at {path}:{line_no}") from exc


def sample_from_region_query(parts: List[str], path: Path, line_no: int) -> str:
    query_name = sample_from_fixed_column(parts, 0, path, line_no)
    match = REGION_QUERY_RE.match(query_name)
    if not match:
        raise ValueError(f"Cannot parse sample ID from query name at {path}:{line_no}: {query_name!r}")
    return match.group("sample")


def filter_no_header(
    input_path: Path,
    output_path: Path,
    sample_ids: Set[str],
    sample_getter: Callable[[List[str], Path, int], str],
) -> FileStats:
    stats = FileStats(input_path=input_path, output_path=output_path, has_header=False)
    tmp_path = output_path.with_suffix(output_path.suffix + ".tmp")
    try:
        with open_gzip_text(input_path, "rt") as src, open_gzip_text(tmp_path, "wt") as dst:
            for line_no, line in enumerate(src, start=1):
                if not line.strip():
                    continue
                stats.total_rows += 1
                parts = split_tsv(line)
                sample_id = sample_getter(parts, input_path, line_no)
                if sample_id in sample_ids:
                    dst.write(line)
                    stats.kept_rows += 1
                else:
                    stats.dropped_rows += 1
        os.replace(tmp_path, output_path)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise
    return stats


def filter_with_header(
    input_path: Path,
    output_path: Path,
    sample_ids: Set[str],
    sample_columns: Iterable[str],
) -> FileStats:
    stats = FileStats(input_path=input_path, output_path=output_path, has_header=True)
    tmp_path = output_path.with_suffix(output_path.suffix + ".tmp")
    try:
        with open_gzip_text(input_path, "rt") as src, open_gzip_text(tmp_path, "wt") as dst:
            header = src.readline()
            if not header:
                raise ValueError(f"Empty input: {input_path}")
            dst.write(header)
            sample_idx = get_index(build_header_index(header), input_path, sample_columns)
            for line_no, line in enumerate(src, start=2):
                if not line.strip():
                    continue
                stats.total_rows += 1
                parts = split_tsv(line)
                sample_id = sample_from_fixed_column(parts, sample_idx, input_path, line_no)
                if sample_id in sample_ids:
                    dst.write(line)
                    stats.kept_rows += 1
                else:
                    stats.dropped_rows += 1
        os.replace(tmp_path, output_path)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise
    return stats


def filter_one_file(input_path: Path, output_path: Path, sample_ids: Set[str]) -> FileStats:
    name = input_path.name
    if name in {"1-all_breakpoints.tsv.gz", "3-confident_breakpoints.tsv.gz"}:
        return filter_with_header(input_path, output_path, sample_ids, ("sampleID", "source_sample_id"))
    if name == "5-mt_disc_cluster_old.tsv.gz":
        return filter_with_header(input_path, output_path, sample_ids, ("IndividualID",))
    if name == "7-mt_disc_cluster.tsv.gz":
        return filter_with_header(input_path, output_path, sample_ids, ("IndividualID",))
    if name == "2-breakpoints_old.tsv.gz":
        return filter_no_header(
            input_path, output_path, sample_ids, lambda parts, path, line_no: sample_from_fixed_column(parts, 7, path, line_no)
        )
    if name == "4-mt_disc_breakpoint_input.tsv.gz":
        return filter_no_header(
            input_path, output_path, sample_ids, lambda parts, path, line_no: sample_from_fixed_column(parts, 0, path, line_no)
        )
    if name == "6-mt_disc_cluster_summary.tsv.gz":
        return filter_no_header(
            input_path, output_path, sample_ids, lambda parts, path, line_no: sample_from_fixed_column(parts, 1, path, line_no)
        )
    if name in {"8-merge_bed.tsv.gz", "9-merge_bed_confident.tsv.gz"}:
        return filter_no_header(input_path, output_path, sample_ids, sample_from_region_query)
    raise ValueError(f"Unsupported input file: {input_path}")


def run(
    id_file: str,
    input_dir: str,
    output_dir: str,
    overwrite: bool = False,
    jobs: int = 1,
) -> Dict[str, FileStats]:
    id_path = Path(id_file).expanduser().resolve()
    input_path = Path(input_dir).expanduser().resolve()
    output_path = Path(output_dir).expanduser().resolve()

    if not id_path.exists():
        raise FileNotFoundError(f"Sample ID file not found: {id_path}")
    if jobs < 1:
        raise ValueError("--jobs must be >= 1")

    sample_ids = load_sample_ids(id_path)
    input_paths = validate_input_dir(input_path)
    ensure_outputs_can_be_written(output_path, overwrite)

    LOG.info("Loaded %d sample IDs from %s", len(sample_ids), id_path)
    LOG.info("Input directory: %s", input_path)
    LOG.info("Output directory: %s", output_path)
    LOG.info("Parallel jobs: %d", jobs)

    results: Dict[str, FileStats] = {}
    max_workers = min(jobs, len(input_paths))
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_name = {
            executor.submit(filter_one_file, input_paths[name], output_path / name, sample_ids): name
            for name in input_paths
        }
        for future in as_completed(future_to_name):
            name = future_to_name[future]
            stats = future.result()
            results[name] = stats
            LOG.info(
                "%s input_rows=%d kept_rows=%d dropped_rows=%d output=%s",
                name,
                stats.total_rows,
                stats.kept_rows,
                stats.dropped_rows,
                stats.output_path,
            )
    return results


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Filter all-sample NUMTs TSV.GZ outputs to selected sample IDs."
    )
    parser.add_argument("--id-file", required=True, help="One-column sample ID list.")
    parser.add_argument("--input-dir", required=True, help="Directory with all-sample TSV.GZ inputs.")
    parser.add_argument("--output-dir", required=True, help="Directory for filtered TSV.GZ outputs.")
    parser.add_argument("--jobs", type=int, default=1, help="Number of TSV.GZ files to filter in parallel.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs.")
    return parser


def main(argv=None) -> int:
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    parser = build_parser()
    args = parser.parse_args(argv)
    run(
        id_file=args.id_file,
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        overwrite=args.overwrite,
        jobs=args.jobs,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
