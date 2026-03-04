#!/usr/bin/env python3
import argparse
import json
import sys


def parse_force(raw_value):
    if isinstance(raw_value, bool):
        return 1 if raw_value else 0
    if isinstance(raw_value, int):
        if raw_value in (0, 1):
            return raw_value
        raise ValueError(f"Invalid force value: {raw_value!r}")
    if isinstance(raw_value, str):
        normalized = raw_value.strip().lower()
        if normalized in {"1", "true", "yes", "y"}:
            return 1
        if normalized in {"0", "false", "no", "n", ""}:
            return 0
    raise ValueError(f"Invalid force value: {raw_value!r}")


def main():
    parser = argparse.ArgumentParser(
        description="Read jobs/force runtime settings from merge config JSON."
    )
    parser.add_argument("--config", required=True, help="Path to merge config JSON")
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as handle:
        config = json.load(handle)

    jobs = config.get("jobs", 4)
    if not isinstance(jobs, int) or jobs < 1:
        raise ValueError(f"Invalid jobs value: {jobs!r}")

    force = parse_force(config.get("force", False))

    # Keep output stable for shell parsing via mapfile.
    sys.stdout.write(f"{jobs}\n{force}\n")


if __name__ == "__main__":
    main()
