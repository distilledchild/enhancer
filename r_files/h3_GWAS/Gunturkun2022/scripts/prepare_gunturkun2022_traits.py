#!/usr/bin/env python3
"""Prepare Gunturkun 2022 HS rat phenotypes for GCTA.

Inputs:
  - GeneNetwork HSNIH-PalmerPublish phenotype matrix CSV
  - PLINK .fam from HS_genotypes_v4

Outputs:
  - trait-specific GCTA phenotype files: FID IID PHENO
  - keep files for each trait
  - a union keep file for all target traits
  - metadata TSV describing sample counts and applied covariate adjustment

This intentionally makes a reproducible local approximation of the paper's
phenotype preprocessing: sex-stratified rank inverse-normal transform, optional
categorical covariate residualization when covariates explain >2% variance, and
a final rank inverse-normal transform.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np


TRAITS = list(range(10418, 10441))
COVARIATES = {
    "sex": "HSR_10443",
    "batch_within_center": "HSR_10444",
    "coat_color": "HSR_10445",
    "center": "HSR_10446",
}


def inv_norm_cdf(p: float) -> float:
    """Acklam's inverse normal CDF approximation."""
    if not 0.0 < p < 1.0:
        raise ValueError("p must be in (0, 1)")

    a = [
        -3.969683028665376e01,
        2.209460984245205e02,
        -2.759285104469687e02,
        1.383577518672690e02,
        -3.066479806614716e01,
        2.506628277459239e00,
    ]
    b = [
        -5.447609879822406e01,
        1.615858368580409e02,
        -1.556989798598866e02,
        6.680131188771972e01,
        -1.328068155288572e01,
    ]
    c = [
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e00,
        -2.549732539343734e00,
        4.374664141464968e00,
        2.938163982698783e00,
    ]
    d = [
        7.784695709041462e-03,
        3.224671290700398e-01,
        2.445134137142996e00,
        3.754408661907416e00,
    ]

    plow = 0.02425
    phigh = 1.0 - plow
    if p < plow:
        q = math.sqrt(-2.0 * math.log(p))
        return (
            (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5])
            / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
        )
    if p <= phigh:
        q = p - 0.5
        r = q * q
        return (
            (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5])
            * q
            / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0)
        )
    q = math.sqrt(-2.0 * math.log(1.0 - p))
    return -(
        (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5])
        / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
    )


def parse_float(value: str) -> float | None:
    if value in {"", "x", "NA", "NaN", "nan", "-9"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def rank_inverse_normal(values: np.ndarray) -> np.ndarray:
    n = len(values)
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(n, dtype=float)
    i = 0
    while i < n:
        j = i + 1
        while j < n and values[order[j]] == values[order[i]]:
            j += 1
        avg_rank = (i + 1 + j) / 2.0
        ranks[order[i:j]] = avg_rank
        i = j
    probs = (ranks - 0.5) / n
    return np.array([inv_norm_cdf(float(p)) for p in probs], dtype=float)


def sex_stratified_rint(values: np.ndarray, sex: np.ndarray) -> np.ndarray:
    out = np.empty_like(values, dtype=float)
    for s in sorted(set(sex.tolist())):
        idx = np.where(sex == s)[0]
        if len(idx) < 3:
            out[idx] = rank_inverse_normal(values[idx])
        else:
            out[idx] = rank_inverse_normal(values[idx])
    return out


def categorical_design(rows: list[dict[str, str]], covariates: list[str]) -> np.ndarray:
    cols = [np.ones(len(rows), dtype=float)]
    for cov in covariates:
        raw = [rows[i].get(COVARIATES[cov], "x") for i in range(len(rows))]
        levels = sorted({v for v in raw if v not in {"", "x", "NA", "nan"}})
        for level in levels[1:]:
            cols.append(np.array([1.0 if v == level else 0.0 for v in raw], dtype=float))
    return np.column_stack(cols)


def residualize_if_needed(y: np.ndarray, rows: list[dict[str, str]], covariates: list[str]) -> tuple[np.ndarray, str, float]:
    if not covariates:
        return y, "none", 0.0
    x0 = np.ones((len(y), 1), dtype=float)
    rss0 = float(np.sum((y - x0 @ np.linalg.lstsq(x0, y, rcond=None)[0]) ** 2))
    x = categorical_design(rows, covariates)
    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    fitted = x @ beta
    rss1 = float(np.sum((y - fitted) ** 2))
    r2 = 0.0 if rss0 == 0 else max(0.0, 1.0 - rss1 / rss0)
    if r2 <= 0.02 or x.shape[1] <= 1:
        return y, "none_r2_le_0.02", r2
    residual = y - fitted + float(np.mean(y))
    return residual, "+".join(covariates), r2


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phenotype-csv", required=True)
    parser.add_argument("--fam", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    pheno_csv = Path(args.phenotype_csv)
    fam_path = Path(args.fam)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fam_iids = set()
    with fam_path.open() as handle:
        for line in handle:
            parts = line.split()
            if len(parts) >= 2:
                fam_iids.add(parts[1])

    with pheno_csv.open(newline="") as handle:
        reader = csv.DictReader(handle)
        matrix_rows = [row for row in reader if row["id"] in fam_iids]

    union_ids: set[str] = set()
    metadata_rows = []
    for trait in TRAITS:
        col = f"HSR_{trait}"
        usable_rows = []
        raw_values = []
        sex_values = []
        for row in matrix_rows:
            y = parse_float(row.get(col, "x"))
            sex = parse_float(row.get(COVARIATES["sex"], "x"))
            if y is None or sex is None:
                continue
            usable_rows.append(row)
            raw_values.append(y)
            sex_values.append(sex)

        if not usable_rows:
            continue

        y0 = np.array(raw_values, dtype=float)
        sex = np.array(sex_values, dtype=float)
        y1 = sex_stratified_rint(y0, sex)
        y2, covar_used, covar_r2 = residualize_if_needed(
            y1, usable_rows, ["batch_within_center", "coat_color", "center"]
        )
        y_final = rank_inverse_normal(y2)

        pheno_path = out_dir / f"trait_{trait}.pheno"
        keep_path = out_dir / f"trait_{trait}.keep"
        with pheno_path.open("w") as pheno, keep_path.open("w") as keep:
            for row, value in zip(usable_rows, y_final):
                iid = row["id"]
                pheno.write(f"0\t{iid}\t{value:.10g}\n")
                keep.write(f"0\t{iid}\n")
                union_ids.add(iid)

        metadata_rows.append(
            {
                "trait": str(trait),
                "n": str(len(usable_rows)),
                "covariates_used": covar_used,
                "covariate_r2": f"{covar_r2:.6g}",
                "pheno": str(pheno_path),
                "keep": str(keep_path),
            }
        )

    union_path = out_dir / "target_traits_union.keep"
    with union_path.open("w") as handle:
        for iid in sorted(union_ids):
            handle.write(f"0\t{iid}\n")

    meta_path = out_dir / "trait_metadata.tsv"
    with meta_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["trait", "n", "covariates_used", "covariate_r2", "pheno", "keep"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(metadata_rows)

    print(f"wrote {len(metadata_rows)} trait phenotype files")
    print(f"union keep: {union_path} ({len(union_ids)} IDs)")
    print(f"metadata: {meta_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
