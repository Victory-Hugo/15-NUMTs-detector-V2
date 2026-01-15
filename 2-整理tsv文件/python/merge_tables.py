import argparse
import json
import logging
import os
import shutil
from collections import Counter
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

LOGGER = logging.getLogger(__name__)


def _read_first_line(file_path: Path) -> str:
    with file_path.open("r", encoding="utf-8") as handle:
        first_line = handle.readline()
    return first_line.rstrip("\n")


def _iter_data_lines(file_path: Path, skip_first: bool) -> Iterable[str]:
    with file_path.open("r", encoding="utf-8") as handle:
        if skip_first:
            _ = handle.readline()
        for line in handle:
            yield line


def _detect_header(files: List[Path]) -> Tuple[str, int]:
    first_lines = [_read_first_line(file_path) for file_path in files]
    counter = Counter(first_lines)
    header, count = counter.most_common(1)[0]
    if count < 2:
        return "", 0
    return header, count


def _collect_files(input_dir: Path, pattern: str) -> List[Path]:
    return sorted(input_dir.rglob(pattern))


def run(
    config_path: str,
    suffix_key: str,
    output_path: Optional[str] = None,
    force: bool = False,
) -> str:
    with open(config_path, "r", encoding="utf-8") as handle:
        config = json.load(handle)

    suffix_groups = config.get("suffix_groups", {})
    if suffix_key not in suffix_groups:
        raise KeyError(f"suffix_key not found in config: {suffix_key}")

    group_conf = suffix_groups[suffix_key]
    pattern = group_conf["pattern"]

    input_dir = Path(config["input_dir"]).resolve()
    output_dir = Path(config["output_dir"]).resolve()
    tmp_dir = Path(config["tmp_dir"]).resolve()
    output_name = group_conf["output_name"]

    output_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    final_output = Path(output_path).resolve() if output_path else output_dir / output_name

    if final_output.exists() and not force:
        raise FileExistsError(f"Output exists; use --force to overwrite: {final_output}")

    files = _collect_files(input_dir, pattern)
    if not files:
        raise FileNotFoundError(f"No files match pattern {pattern} under {input_dir}")

    header, header_count = _detect_header(files)
    if header_count == 0:
        LOGGER.info("No shared header detected; concatenate without header")

    temp_output = final_output.with_suffix(final_output.suffix + ".tmp")

    LOGGER.info("Merging %d files for key %s", len(files), suffix_key)
    LOGGER.info("Output: %s", final_output)

    with temp_output.open("w", encoding="utf-8") as out_handle:
        if header:
            out_handle.write(header + "\n")
        for file_path in files:
            first_line = _read_first_line(file_path)
            skip_first = bool(header) and first_line == header
            if header and first_line != header and first_line.startswith("#"):
                raise ValueError(
                    "Header mismatch: "
                    f"expected {header} but got {first_line} in {file_path}"
                )
            for line in _iter_data_lines(file_path, skip_first=skip_first):
                out_handle.write(line)

    try:
        os.replace(temp_output, final_output)
    except OSError as exc:
        if exc.errno != 18:
            raise
        shutil.move(str(temp_output), str(final_output))
    return str(final_output)


def _list_keys(config_path: str) -> List[str]:
    with open(config_path, "r", encoding="utf-8") as handle:
        config = json.load(handle)
    return sorted(config.get("suffix_groups", {}).keys())


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Merge TSV files by suffix group.")
    parser.add_argument("--config", required=True, help="Path to JSON config file")
    parser.add_argument("--suffix-key", help="Key in config.suffix_groups to merge")
    parser.add_argument(
        "--output-path",
        help="Override output path; default is config.output_dir/output_name",
    )
    parser.add_argument("--force", action="store_true", help="Overwrite output")
    parser.add_argument(
        "--list-keys",
        action="store_true",
        help="List available suffix keys and exit",
    )
    return parser


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    parser = build_arg_parser()
    args = parser.parse_args()

    if args.list_keys:
        keys = _list_keys(args.config)
        output = "\n".join(keys)
        os.write(1, (output + "\n").encode("utf-8"))
        return

    if not args.suffix_key:
        raise ValueError("--suffix-key is required unless --list-keys is used")

    run(
        config_path=args.config,
        suffix_key=args.suffix_key,
        output_path=args.output_path,
        force=args.force,
    )


if __name__ == "__main__":
    main()
