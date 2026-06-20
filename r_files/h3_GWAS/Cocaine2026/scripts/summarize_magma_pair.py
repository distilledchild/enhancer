#!/usr/bin/env python3
"""Summarize cMAGMA vs H-MAGMA gene results for one phenotype."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def read_genes(path: Path, method: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", engine="python")
    df["METHOD"] = method
    df["KEY"] = (
        df["CHR"].astype(str)
        + ":"
        + df["START"].astype(str)
        + ":"
        + df["STOP"].astype(str)
    )
    df["BONF_0_05"] = df["P"] < (0.05 / len(df))
    df["BH_FDR"] = df["P"].rank(method="first") / len(df)
    ordered = df.sort_values("P").copy()
    m = len(ordered)
    q = ordered["P"] * m / pd.Series(range(1, m + 1), index=ordered.index)
    ordered["BH_Q"] = q[::-1].cummin()[::-1].clip(upper=1)
    df = df.drop(columns=["BH_FDR"]).merge(
        ordered[["GENE", "KEY", "BH_Q"]], on=["GENE", "KEY"], how="left"
    )
    return df


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trait", required=True)
    parser.add_argument("--cmagma", required=True, type=Path)
    parser.add_argument("--hmagma", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    c = read_genes(args.cmagma, "cMAGMA")
    h = read_genes(args.hmagma, "H-MAGMA_strict_Duttke_Telese_nonTSS")

    rows = []
    for label, threshold in [
        ("BH_FDR_0.05", lambda d: d["BH_Q"] < 0.05),
        ("BH_FDR_0.10", lambda d: d["BH_Q"] < 0.10),
        ("Bonferroni_0.05", lambda d: d["BONF_0_05"]),
    ]:
        c_sig = c[threshold(c)]
        h_sig = h[threshold(h)]
        c_keys = set(c_sig["KEY"])
        h_only = h_sig[~h_sig["KEY"].isin(c_keys)].copy()
        rows.append(
            {
                "trait": args.trait,
                "criterion": label,
                "cmagma_significant": len(c_sig),
                "hmagma_significant": len(h_sig),
                "hmagma_only_vs_cmagma_sig": len(h_only),
            }
        )
        h_only.sort_values(["P", "BH_Q"]).to_csv(
            args.out_dir / f"{args.trait}.{label}.hmagma_only.tsv",
            sep="\t",
            index=False,
        )

    summary = pd.DataFrame(rows)
    summary.to_csv(args.out_dir / f"{args.trait}.summary.tsv", sep="\t", index=False)

    merged = h.merge(
        c[["KEY", "GENE", "NSNPS", "ZSTAT", "P", "BH_Q", "BONF_0_05"]].rename(
            columns={
                "GENE": "CMAGMA_GENE",
                "NSNPS": "CMAGMA_NSNPS",
                "ZSTAT": "CMAGMA_ZSTAT",
                "P": "CMAGMA_P",
                "BH_Q": "CMAGMA_BH_Q",
                "BONF_0_05": "CMAGMA_BONF_0_05",
            }
        ),
        on="KEY",
        how="left",
    ).rename(
        columns={
            "GENE": "HMAGMA_GENE",
            "NSNPS": "HMAGMA_NSNPS",
            "ZSTAT": "HMAGMA_ZSTAT",
            "P": "HMAGMA_P",
            "BH_Q": "HMAGMA_BH_Q",
            "BONF_0_05": "HMAGMA_BONF_0_05",
        }
    )
    merged["NSNP_DELTA_H_MINUS_C"] = merged["HMAGMA_NSNPS"] - merged["CMAGMA_NSNPS"]
    merged.sort_values(["HMAGMA_P", "HMAGMA_BH_Q"]).to_csv(
        args.out_dir / f"{args.trait}.hmagma_with_cmagma_match.tsv",
        sep="\t",
        index=False,
    )

    print(summary.to_string(index=False))
    print("\nTop H-MAGMA rows:")
    cols = [
        "HMAGMA_GENE",
        "KEY",
        "HMAGMA_NSNPS",
        "CMAGMA_NSNPS",
        "NSNP_DELTA_H_MINUS_C",
        "HMAGMA_ZSTAT",
        "CMAGMA_ZSTAT",
        "HMAGMA_P",
        "CMAGMA_P",
        "HMAGMA_BH_Q",
        "CMAGMA_BH_Q",
        "HMAGMA_BONF_0_05",
        "CMAGMA_BONF_0_05",
    ]
    print(merged.sort_values("HMAGMA_P")[cols].head(20).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
