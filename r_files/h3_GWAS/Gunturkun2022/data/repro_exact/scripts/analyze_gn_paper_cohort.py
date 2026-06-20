#!/usr/bin/env python3
"""Derive paper trait coverage and candidate cohort lists from GeneNetwork data."""

from __future__ import annotations

import collections
import csv
import sys
from pathlib import Path


MISSING = {"", "NA", "NaN", "nan", "x", "X", None}
NOIT = list(range(10418, 10424))
OFT = list(range(10424, 10430))
SIT = list(range(10430, 10441))
BEHAVIOR = NOIT + OFT + SIT
COVARIATES = [10443, 10444, 10445, 10446]


def numeric_value(row: dict[str, str], trait_id: int) -> str | None:
    value = row.get(f"HSR_{trait_id}", "")
    if value in MISSING:
        return None
    try:
        float(value)
    except (TypeError, ValueError):
        return None
    return value


def load_metadata(path: Path) -> dict[int, dict[str, str]]:
    meta: dict[int, dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("trait_id"):
                meta[int(row["trait_id"])] = row
    return meta


def write_counts(rows: list[dict[str, str]], meta: dict[int, dict[str, str]], out_dir: Path) -> None:
    fields = ["trait_id", "name", "nonmissing", "sex_0_female", "sex_1_male", "center_values", "batch_values"]
    with (out_dir / "behavior_trait_nonmissing_counts.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for trait_id in BEHAVIOR + COVARIATES:
            sample_rows = [row for row in rows if numeric_value(row, trait_id) is not None]
            sex_counts = collections.Counter(numeric_value(row, 10443) for row in sample_rows)
            center_counts = collections.Counter(numeric_value(row, 10446) for row in sample_rows)
            batch_counts = collections.Counter(numeric_value(row, 10444) for row in sample_rows)
            writer.writerow(
                {
                    "trait_id": trait_id,
                    "name": meta.get(trait_id, {}).get("name", ""),
                    "nonmissing": len(sample_rows),
                    "sex_0_female": sex_counts.get("0.0", 0),
                    "sex_1_male": sex_counts.get("1.0", 0),
                    "center_values": ";".join(f"{k}:{v}" for k, v in sorted(center_counts.items())),
                    "batch_values": ";".join(
                        f"{int(float(k))}:{v}" for k, v in sorted(batch_counts.items(), key=lambda kv: float(kv[0]))
                    ),
                }
            )


def sample_set(rows: list[dict[str, str]], traits: list[int], mode: str) -> set[str]:
    if mode == "any":
        return {row["id"] for row in rows if any(numeric_value(row, trait_id) is not None for trait_id in traits)}
    if mode == "all":
        return {row["id"] for row in rows if all(numeric_value(row, trait_id) is not None for trait_id in traits)}
    raise ValueError(mode)


def write_set_summary(rows: list[dict[str, str]], out_dir: Path) -> dict[str, set[str]]:
    sets = {
        "NOIT_any": sample_set(rows, NOIT, "any"),
        "NOIT_all": sample_set(rows, NOIT, "all"),
        "OFT_any": sample_set(rows, OFT, "any"),
        "OFT_all": sample_set(rows, OFT, "all"),
        "OFT_totaldistance_10424": sample_set(rows, [10424], "all"),
        "SIT_any": sample_set(rows, SIT, "any"),
        "SIT_all": sample_set(rows, SIT, "all"),
        "ALL23_any": sample_set(rows, BEHAVIOR, "any"),
        "ALL23_all": sample_set(rows, BEHAVIOR, "all"),
    }
    sets["OFT_any_and_NOIT_any_and_SIT_any"] = sets["OFT_any"] & sets["NOIT_any"] & sets["SIT_any"]
    sets["OFT10424_and_NOIT_any_and_SIT_any"] = sets["OFT_totaldistance_10424"] & sets["NOIT_any"] & sets["SIT_any"]

    row_by_id = {row["id"]: row for row in rows}
    fields = ["set_name", "n", "female_0", "male_1", "center_values", "batch_values"]
    with (out_dir / "paper_cohort_candidate_sets.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for name, ids in sets.items():
            selected = [row_by_id[sample_id] for sample_id in ids]
            sex_counts = collections.Counter(numeric_value(row, 10443) for row in selected)
            center_counts = collections.Counter(numeric_value(row, 10446) for row in selected)
            batch_counts = collections.Counter(numeric_value(row, 10444) for row in selected)
            writer.writerow(
                {
                    "set_name": name,
                    "n": len(ids),
                    "female_0": sex_counts.get("0.0", 0),
                    "male_1": sex_counts.get("1.0", 0),
                    "center_values": ";".join(f"{k}:{v}" for k, v in sorted(center_counts.items())),
                    "batch_values": ";".join(
                        f"{int(float(k))}:{v}" for k, v in sorted(batch_counts.items(), key=lambda kv: float(kv[0]))
                    ),
                }
            )
    return sets


def write_candidate_cohort(rows: list[dict[str, str]], meta: dict[int, dict[str, str]], ids: set[str], out_dir: Path) -> None:
    fields = ["sample_id"] + [f"{trait_id}:{meta.get(trait_id, {}).get('name', '')}" for trait_id in BEHAVIOR + COVARIATES]
    with (out_dir / "cohort_oft_totaldistance_10424_n1246.samples.txt").open("w", encoding="utf-8") as handle:
        for sample_id in sorted(ids):
            handle.write(sample_id + "\n")
    with (out_dir / "cohort_oft_totaldistance_10424_n1246.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(fields)
        for row in sorted((row for row in rows if row["id"] in ids), key=lambda r: r["id"]):
            writer.writerow([row["id"]] + [numeric_value(row, trait_id) or "" for trait_id in BEHAVIOR + COVARIATES])


def main() -> int:
    if len(sys.argv) != 4:
        print("Usage: analyze_gn_paper_cohort.py HSNIH-PalmerPublish.csv METADATA.tsv OUTPUT_DIR", file=sys.stderr)
        return 2
    csv_path = Path(sys.argv[1])
    metadata_path = Path(sys.argv[2])
    out_dir = Path(sys.argv[3])
    out_dir.mkdir(parents=True, exist_ok=True)

    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    meta = load_metadata(metadata_path)

    write_counts(rows, meta, out_dir)
    sets = write_set_summary(rows, out_dir)
    candidate = sets["OFT_totaldistance_10424"]
    write_candidate_cohort(rows, meta, candidate, out_dir)
    selected = [row for row in rows if row["id"] in candidate]
    sex_counts = collections.Counter(numeric_value(row, 10443) for row in selected)
    matches_paper = len(candidate) == 1246 and sex_counts.get("0.0", 0) == 620 and sex_counts.get("1.0", 0) == 626

    summary = [
        f"source_csv\t{csv_path}",
        f"source_metadata\t{metadata_path}",
        f"rows\t{len(rows)}",
        f"candidate_set\tOFT_totaldistance_10424",
        f"candidate_n\t{len(candidate)}",
        f"candidate_female_0\t{sex_counts.get('0.0', 0)}",
        f"candidate_male_1\t{sex_counts.get('1.0', 0)}",
        f"candidate_matches_paper_n_and_sex\t{str(matches_paper).lower()}",
        "paper_reported_n\t1246",
        "paper_reported_female_male\t620 female / 626 male",
    ]
    (out_dir / "paper_cohort_analysis_summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("\n".join(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
