#!/usr/bin/env python3
"""解析与基础计算工具。"""

from __future__ import annotations

import ast
import math
import re
from typing import Iterable, List, Sequence

import numpy as np

LEFT_GROUP = "nu_Tend_Bleft"
RIGHT_GROUP = "nu_Tstart_Bright"
MT_START_GROUP = "mt_Tstart"
MT_END_GROUP = "mt_Tend"


def parse_int_list(text: object) -> List[int]:
    """把 cluster 表中的字符串列表解析成整数列表。"""
    if text is None or (isinstance(text, float) and math.isnan(text)):
        return []
    if isinstance(text, list):
        return [int(x) for x in text]
    value = str(text).strip()
    if value == "":
        return []
    try:
        parsed = ast.literal_eval(value)
    except (SyntaxError, ValueError) as exc:
        raise ValueError(f"无法解析整数列表: {value}") from exc
    if not isinstance(parsed, (list, tuple)):
        raise ValueError(f"列表字段不是 list/tuple: {value}")
    result = []
    for item in parsed:
        result.append(int(item))
    return result


def normalize_chr_name(value: object) -> str:
    """尽量保留 chr 风格，不主动猜测复杂命名。"""
    text = str(value).strip()
    if text == "":
        return text
    if text.startswith("chr"):
        return text
    if text in {"M", "MT"}:
        return "chrM"
    return f"chr{text}"


def build_region_key(sample_id: object, chrom: object, start: object, end: object) -> str:
    return f"{sample_id}_{chrom}.{int(start)}.{int(end)}"


def parse_length_region_key(raw_key: object) -> dict:
    """解析长度 bed 第一列中的区域键。"""
    token = str(raw_key).split("|", 1)[0]
    match = re.match(r"^(?P<sample_id>.+)_(?P<chr>chr[^.]+)\.(?P<start>\d+)\.(?P<end>\d+)$", token)
    if not match:
        raise ValueError(f"无法解析长度 bed 区域键: {raw_key}")
    return {
        "sample_id": match.group("sample_id"),
        "chr": match.group("chr"),
        "start": int(match.group("start")),
        "end": int(match.group("end")),
        "region_key": token,
    }


def weighted_median(values: Sequence[float], weights: Sequence[float] | None = None) -> float | None:
    """计算加权中位数。"""
    if values is None or len(values) == 0:
        return None
    arr = np.asarray(values, dtype=float)
    if weights is None:
        weights_arr = np.ones(arr.shape[0], dtype=float)
    else:
        weights_arr = np.asarray(weights, dtype=float)
    order = np.argsort(arr)
    arr = arr[order]
    weights_arr = weights_arr[order]
    total = weights_arr.sum()
    if total <= 0:
        weights_arr = np.ones(arr.shape[0], dtype=float)
        total = weights_arr.sum()
    cumulative = np.cumsum(weights_arr)
    idx = int(np.searchsorted(cumulative, total / 2.0, side="left"))
    idx = min(max(idx, 0), arr.shape[0] - 1)
    return float(arr[idx])


def single_linkage_labels(values: Iterable[float], gap: int) -> List[int]:
    """对排序后的一维数值做单链聚类标签。"""
    labels: List[int] = []
    current_label = -1
    previous = None
    for value in values:
        if previous is None or float(value) - float(previous) > gap:
            current_label += 1
        labels.append(current_label)
        previous = float(value)
    return labels


def safe_midpoint(start: float | int | None, end: float | int | None) -> float | None:
    if start is None or end is None or np.isnan(start) or np.isnan(end):
        return None
    return (float(start) + float(end)) / 2.0


def classify_frequency(
    carrier_count: int,
    carrier_frequency: float,
    prevalent_min: float,
    common_min: float,
    rare_min: float,
) -> str:
    """按约定顺序计算频率类别。"""
    if int(carrier_count) == 1:
        return "private"
    if float(carrier_frequency) > prevalent_min:
        return "prevalent"
    if common_min <= float(carrier_frequency) <= prevalent_min:
        return "common"
    if rare_min <= float(carrier_frequency) < common_min:
        return "rare"
    return "ultra_rare"

