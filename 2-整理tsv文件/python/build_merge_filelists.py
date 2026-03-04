#!/usr/bin/env python3
import argparse
import json
from pathlib import Path
from typing import Dict, List, Tuple


def collect_files(input_dir: Path, pattern: str) -> List[Path]:
    return sorted(input_dir.rglob(pattern))


def build_lists(
    config_path: Path,
    keys_file: Path,
    lists_dir: Path,
) -> Tuple[Path, Path]:
    with config_path.open("r", encoding="utf-8") as handle:
        config = json.load(handle)

    input_dir = Path(config["input_dir"]).expanduser().resolve()
    suffix_groups: Dict[str, Dict[str, str]] = config.get("suffix_groups", {})
    if not suffix_groups:
        raise ValueError("No suffix_groups found in config")

    lists_dir.mkdir(parents=True, exist_ok=True)
    keys_file.parent.mkdir(parents=True, exist_ok=True)

    mergeable_keys: List[str] = []
    summary_lines: List[str] = []

    for suffix_key in sorted(suffix_groups):
        group_conf = suffix_groups[suffix_key]
        pattern = group_conf["pattern"]
        files = collect_files(input_dir, pattern)

        file_list_path = lists_dir / f"{suffix_key}.txt"
        with file_list_path.open("w", encoding="utf-8") as out_handle:
            for file_path in files:
                out_handle.write(str(file_path) + "\n")

        if files:
            mergeable_keys.append(suffix_key)

        summary_lines.append(f"{suffix_key}\t{len(files)}\t{file_list_path}")

    with keys_file.open("w", encoding="utf-8") as keys_handle:
        for key in mergeable_keys:
            keys_handle.write(key + "\n")

    if not mergeable_keys:
        raise FileNotFoundError("No input files were discovered for any suffix group")

    for line in summary_lines:
        print(line)
    print(f"MERGE_KEYS\t{len(mergeable_keys)}\t{keys_file}")
    print(f"MERGE_LIST_DIR\t{lists_dir}")

    return keys_file, lists_dir


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build per-group txt file lists for two-step TSV merge workflow."
    )
    parser.add_argument("--config", required=True, help="Path to JSON config file")
    parser.add_argument("--keys-file", required=True, help="Output path for merge_keys.txt")
    parser.add_argument(
        "--lists-dir",
        required=True,
        help="Directory to write per-group txt file lists",
    )
    args = parser.parse_args()

    build_lists(
        config_path=Path(args.config).expanduser().resolve(),
        keys_file=Path(args.keys_file).expanduser().resolve(),
        lists_dir=Path(args.lists_dir).expanduser().resolve(),
    )


if __name__ == "__main__":
    main()
