#!/usr/bin/env python3
"""Convert GCTA MLMA output to MAGMA p-value input."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mlma", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    mlma = Path(args.mlma)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    with mlma.open() as src, out.open("w", newline="") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"No header in {mlma}")
        p_col = "p" if "p" in reader.fieldnames else "P"
        snp_col = "SNP"
        writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
        writer.writerow(["SNP", "p"])
        kept = 0
        for row in reader:
            snp = row.get(snp_col, "")
            p = row.get(p_col, "")
            if not snp or not p or p == "NA":
                continue
            try:
                p_value = float(p)
            except ValueError:
                continue
            if not math.isfinite(p_value) or p_value <= 0 or p_value > 1:
                continue
            writer.writerow([snp, p])
            kept += 1
    print(f"wrote {kept} SNP p-values to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
