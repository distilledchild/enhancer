#!/usr/bin/env python3
"""Summarize exact-public 10424 GCTA output against paper Table 2 peaks."""

from __future__ import annotations

import csv
import heapq
import math
import sys
from pathlib import Path


EXPECTED = {
    "chr10:94549701": 7.286,
    "chr11:33359859": 8.268,
}


def neglog10(p: float) -> float:
    return -math.log10(p) if p > 0 else float("inf")


def read_rows(mlma: Path):
    with mlma.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            try:
                row["_p"] = float(row["p"])
                row["_bp"] = int(row["bp"])
                row["_chr"] = str(row["Chr"])
            except (KeyError, ValueError):
                continue
            if math.isfinite(row["_p"]) and row["_p"] > 0:
                row["_mlog10p"] = neglog10(row["_p"])
                yield row


def write_rows(rows, path: Path, fields):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> int:
    if len(sys.argv) != 3:
        print("Usage: summarize_exact_10424.py MLMA OUT_DIR", file=sys.stderr)
        return 2

    mlma = Path(sys.argv[1])
    out_dir = Path(sys.argv[2])
    out_dir.mkdir(parents=True, exist_ok=True)

    top_heap = []
    expected_hits = {}
    windows = {snp: [] for snp in EXPECTED}
    n = 0
    for row in read_rows(mlma):
        n += 1
        record = {
            "Chr": row["Chr"],
            "SNP": row["SNP"],
            "bp": row["bp"],
            "A1": row["A1"],
            "A2": row["A2"],
            "Freq": row["Freq"],
            "b": row["b"],
            "se": row["se"],
            "p": row["p"],
            "minus_log10p": f"{row['_mlog10p']:.9g}",
        }
        heapq.heappush(top_heap, (row["_mlog10p"], n, record))
        if len(top_heap) > 200:
            heapq.heappop(top_heap)
        if row["SNP"] in EXPECTED:
            expected_hits[row["SNP"]] = record | {"paper_minus_log10p": str(EXPECTED[row["SNP"]])}
        for snp in EXPECTED:
            chrom, bp = snp.replace("chr", "").split(":")
            if row["_chr"] == chrom and abs(row["_bp"] - int(bp)) <= 2_000_000:
                windows[snp].append(record | {"paper_peak": snp, "distance_bp": str(row["_bp"] - int(bp))})

    top_rows = [record for _, _, record in sorted(top_heap, reverse=True)]
    top_rows.sort(key=lambda row: float(row["minus_log10p"]), reverse=True)
    fields = ["Chr", "SNP", "bp", "A1", "A2", "Freq", "b", "se", "p", "minus_log10p"]
    write_rows(top_rows[:50], out_dir / "trait_10424.exact_public.top50.tsv", fields)

    expected_rows = []
    for snp, paper in EXPECTED.items():
        if snp in expected_hits:
            row = expected_hits[snp]
            row["paper_minus_log10p"] = str(paper)
            row["delta_minus_log10p"] = f"{float(row['minus_log10p']) - paper:.9g}"
        else:
            chrom, bp = snp.replace("chr", "").split(":")
            row = {
                "Chr": chrom,
                "SNP": snp,
                "bp": bp,
                "p": "",
                "minus_log10p": "",
                "paper_minus_log10p": str(paper),
                "delta_minus_log10p": "",
            }
        expected_rows.append(row)
    write_rows(expected_rows, out_dir / "trait_10424.exact_public.expected_peaks.tsv",
               fields + ["paper_minus_log10p", "delta_minus_log10p"])

    window_rows = []
    for snp, rows in windows.items():
        rows.sort(key=lambda row: float(row["minus_log10p"]), reverse=True)
        window_rows.extend(rows[:20])
    write_rows(window_rows, out_dir / "trait_10424.exact_public.expected_peak_windows_top20.tsv",
               ["paper_peak", "distance_bp"] + fields)

    summary_lines = [
        f"mlma\t{mlma}",
        f"valid_rows\t{n}",
        "expected_peaks\t" + ",".join(EXPECTED),
    ]
    for row in expected_rows:
        summary_lines.append(
            f"{row['SNP']}\tobserved_minus_log10p={row.get('minus_log10p', '')}\t"
            f"paper_minus_log10p={row.get('paper_minus_log10p', '')}\t"
            f"delta={row.get('delta_minus_log10p', '')}"
        )
    if top_rows:
        top = top_rows[0]
        summary_lines.append(f"top_hit\t{top['SNP']}\tminus_log10p={top['minus_log10p']}\tp={top['p']}")
    (out_dir / "trait_10424.exact_public.summary.txt").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    print("\n".join(summary_lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
