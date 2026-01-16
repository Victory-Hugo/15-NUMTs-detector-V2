#!/usr/bin/env python3
import argparse
import sys
from typing import Any

import yaml


def get_value(cfg: dict, key: str) -> Any:
    cur: Any = cfg
    for part in key.split("."):
        if not isinstance(cur, dict) or part not in cur:
            raise KeyError(key)
        cur = cur[part]
    return cur


def main() -> None:
    parser = argparse.ArgumentParser(description="Read a value from YAML config.")
    parser.add_argument("--config", required=True, help="YAML config path")
    parser.add_argument("--key", required=True, help="Dotted key path")
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    try:
        value = get_value(cfg, args.key)
    except KeyError:
        print("", end="")
        sys.exit(1)

    if isinstance(value, list):
        print("\n".join(str(v) for v in value))
    else:
        print(value)


if __name__ == "__main__":
    main()
