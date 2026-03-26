#!/usr/bin/env python3
"""长度与 mtDNA 区间注释。"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from .catalog_utils import finalize_sample_event_table

LOG = logging.getLogger(__name__)
MTDNA_LENGTH = 16569


def annotate_sample_events(sample_events: pd.DataFrame, length_summary: pd.DataFrame) -> pd.DataFrame:
    return finalize_sample_event_table(sample_events, length_summary)


def _interval_overlap(start1: float, end1: float, start2: int, end2: int) -> int:
    left = max(int(round(start1)), start2)
    right = min(int(round(end1)), end2)
    return max(0, right - left + 1)


def annotate_mt_interval(start: float, end: float, mt_genes: pd.DataFrame) -> dict:
    overlaps = []
    for row in mt_genes.itertuples(index=False):
        overlap = _interval_overlap(start, end, row.POSSTART, row.POSEND)
        if overlap > 0:
            overlaps.append((row.GENENAME, row.GENE_CLASS, overlap, row.POSSTART, row.POSEND))
    if not overlaps:
        return {
            "mt_primary_gene": "Unannotated",
            "mt_primary_gene_class": "Unannotated",
            "mt_overlapping_genes": "",
            "mt_breakpoint_bin_start": int(round(start)) if not pd.isna(start) else np.nan,
            "mt_breakpoint_bin_end": int(round(end)) if not pd.isna(end) else np.nan,
        }
    overlaps.sort(key=lambda item: (-item[2], item[3], item[4], item[0]))
    primary = overlaps[0]
    return {
        "mt_primary_gene": primary[0],
        "mt_primary_gene_class": primary[1],
        "mt_overlapping_genes": ",".join(item[0] for item in overlaps),
        "mt_breakpoint_bin_start": int(round(start)),
        "mt_breakpoint_bin_end": int(round(end)),
    }


def annotate_catalog(
    catalog: pd.DataFrame,
    sample_map: pd.DataFrame,
    sample_events_annotated: pd.DataFrame,
    mt_genes: pd.DataFrame,
) -> pd.DataFrame:
    event_lengths = sample_map.merge(
        sample_events_annotated[
            [
                "sample_event_id",
                "mt_length_for_main_bp",
                "mt_breakpoint_span_bp",
                "mt_source_span_bp",
                "nuclear_breakpoint_span_bp",
            ]
        ],
        on="sample_event_id",
        how="left",
    )
    medians = (
        event_lengths.groupby("distinct_numt_id", as_index=False)
        .agg(
            median_mt_length_for_main_bp=("mt_length_for_main_bp", "median"),
            median_mt_breakpoint_span_bp=("mt_breakpoint_span_bp", "median"),
            median_mt_source_union_span_bp=("mt_source_span_bp", "median"),
            median_nuclear_breakpoint_span_bp=("nuclear_breakpoint_span_bp", "median"),
        )
        .copy()
    )
    annotated = catalog.merge(medians, on="distinct_numt_id", how="left", suffixes=("", "_updated")).copy()
    if "median_mt_length_for_main_bp_updated" in annotated.columns:
        if "median_mt_length_for_main_bp" in annotated.columns:
            annotated["median_mt_length_for_main_bp"] = annotated["median_mt_length_for_main_bp_updated"].fillna(
                annotated["median_mt_length_for_main_bp"]
            )
        else:
            annotated["median_mt_length_for_main_bp"] = annotated["median_mt_length_for_main_bp_updated"]
        annotated = annotated.drop(columns=["median_mt_length_for_main_bp_updated"])
    if "median_mt_breakpoint_span_bp_updated" in annotated.columns:
        if "median_mt_breakpoint_span_bp" in annotated.columns:
            annotated["median_mt_breakpoint_span_bp"] = annotated["median_mt_breakpoint_span_bp_updated"].fillna(
                annotated["median_mt_breakpoint_span_bp"]
            )
        else:
            annotated["median_mt_breakpoint_span_bp"] = annotated["median_mt_breakpoint_span_bp_updated"]
        annotated = annotated.drop(columns=["median_mt_breakpoint_span_bp_updated"])
    if "median_mt_source_union_span_bp_updated" in annotated.columns:
        if "median_mt_source_union_span_bp" in annotated.columns:
            annotated["median_mt_source_union_span_bp"] = annotated["median_mt_source_union_span_bp_updated"].fillna(
                annotated["median_mt_source_union_span_bp"]
            )
        else:
            annotated["median_mt_source_union_span_bp"] = annotated["median_mt_source_union_span_bp_updated"]
        annotated = annotated.drop(columns=["median_mt_source_union_span_bp_updated"])
    annotated["median_mt_source_span_bp"] = annotated["median_mt_source_union_span_bp"]
    if "median_nuclear_breakpoint_span_bp_updated" in annotated.columns:
        if "median_nuclear_breakpoint_span_bp" in annotated.columns:
            annotated["median_nuclear_breakpoint_span_bp"] = annotated[
                "median_nuclear_breakpoint_span_bp_updated"
            ].fillna(annotated["median_nuclear_breakpoint_span_bp"])
        else:
            annotated["median_nuclear_breakpoint_span_bp"] = annotated["median_nuclear_breakpoint_span_bp_updated"]
        annotated = annotated.drop(columns=["median_nuclear_breakpoint_span_bp_updated"])

    mt_annotations = [
        annotate_mt_interval(row.distinct_mt_start, row.distinct_mt_end, mt_genes)
        for row in annotated.itertuples(index=False)
    ]
    mt_annotation_df = pd.DataFrame.from_records(mt_annotations)
    annotated = pd.concat([annotated.reset_index(drop=True), mt_annotation_df], axis=1)
    return annotated
