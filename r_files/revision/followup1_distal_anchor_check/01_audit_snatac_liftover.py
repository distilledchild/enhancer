#!/usr/bin/env python3
"""Audit retention of source rn6 snATAC peaks in the local rn7 peak catalog."""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import h5py
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
REVISION_DIR = SCRIPT_DIR.parent
REVISION_MAIN_DIR = REVISION_DIR / "revision_main"
SOURCE_H5 = SCRIPT_DIR / "inputs/public/GSM5820551_filtered_peak_bc_matrix.h5"
LOCAL_RN7 = (
    REVISION_MAIN_DIR
    / "inputs/data/Duttke2022_snATAC_peaks_rn7.narrowPeak"
)
RESULTS_DIR = SCRIPT_DIR / "results"


def decode(values):
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values]


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    if not SOURCE_H5.exists() or not LOCAL_RN7.exists():
        missing = [str(p) for p in (SOURCE_H5, LOCAL_RN7) if not p.exists()]
        raise FileNotFoundError("Missing required input(s): " + ", ".join(missing))

    with h5py.File(SOURCE_H5, "r") as handle:
        features = handle["matrix/features"]
        source_names = decode(features["name"][:])

    local = pd.read_csv(LOCAL_RN7, sep="\t", header=None, comment="#")
    if local.shape[1] < 4:
        raise ValueError("The local rn7 narrowPeak must contain at least four columns.")
    local = local.iloc[:, :10].copy()
    local.columns = [
        "chr",
        "start0",
        "end0",
        "source_rn6_peak",
        "score",
        "strand",
        "signal",
        "pvalue",
        "qvalue",
        "peak",
    ][: local.shape[1]]

    source_set = set(source_names)
    local_names = local["source_rn6_peak"].astype(str)
    local_set = set(local_names)
    canonical_pattern = re.compile(r"^chr(?:[1-9]|1[0-9]|20|X|Y|M)$")
    duplicate_interval_count = int(
        local.duplicated(subset=["chr", "start0", "end0"]).sum()
    )
    noncanonical_count = int(
        (~local["chr"].astype(str).map(lambda x: bool(canonical_pattern.match(x)))).sum()
    )
    matched_source_count = len(source_set & local_set)
    source_not_retained = sorted(source_set - local_set)
    local_not_in_source = sorted(local_set - source_set)

    summary = pd.DataFrame(
        [
            ("source_rn6_h5_features", len(source_names), "GSM5820551 H5 feature rows"),
            ("source_rn6_unique_features", len(source_set), "Unique H5 feature names"),
            ("local_rn7_peak_rows", len(local), "Local lifted narrowPeak rows"),
            ("local_rn7_unique_source_names", len(local_set), "Unique source names retained in rn7 file"),
            ("source_features_retained", matched_source_count, "Source names represented in local rn7 file"),
            ("source_features_not_retained", len(source_not_retained), "Source names absent from local rn7 file"),
            ("local_names_not_in_source_h5", len(local_not_in_source), "Unexpected local source names"),
            ("duplicate_rn7_intervals", duplicate_interval_count, "Duplicate chr/start/end rows"),
            ("noncanonical_rn7_intervals", noncanonical_count, "Intervals outside chr1-20/X/Y/M"),
            (
                "source_feature_retention_percent",
                round(100.0 * matched_source_count / len(source_set), 6),
                "Name-based retention; not a substitute for a liftOver mapping log",
            ),
        ],
        columns=["metric", "value", "definition"],
    )
    summary.to_csv(RESULTS_DIR / "snatac_rn6_to_rn7_audit.tsv", sep="\t", index=False)

    with gzip.open(
        RESULTS_DIR / "snatac_source_rn6_features_not_retained.tsv.gz", "wt"
    ) as handle:
        handle.write("source_rn6_peak\n")
        handle.writelines(f"{name}\n" for name in source_not_retained)

    if local_not_in_source:
        raise RuntimeError(
            f"Found {len(local_not_in_source)} local peak names absent from the source H5."
        )
    if duplicate_interval_count or noncanonical_count:
        raise RuntimeError("The local rn7 peak catalog failed interval QC.")

    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
