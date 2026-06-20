#!/usr/bin/env python3
"""Audit local round8 genotype files against paper requirements and Table 2 SNPs."""

from __future__ import annotations

import csv
import hashlib
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path


def sha256_file(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def file_manifest(prefix: Path) -> list[dict[str, str]]:
    rows = []
    for ext in (".bed", ".bim", ".fam", ".sample"):
        path = prefix.with_suffix(ext)
        stat = path.stat()
        rows.append(
            {
                "path": str(path),
                "filename": path.name,
                "bytes": str(stat.st_size),
                "mtime_epoch": str(int(stat.st_mtime)),
                "sha256": sha256_file(path),
            }
        )
    return rows


def write_manifest(rows: list[dict[str, str]], out_dir: Path) -> None:
    fields = ["filename", "bytes", "mtime_epoch", "sha256", "path"]
    with (out_dir / "round8_file_manifest.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def audit_bim(bim: Path, table2_tsv: Path, out_dir: Path) -> dict[str, int]:
    chr_counts = Counter()
    chr_min = defaultdict(lambda: None)
    chr_max = defaultdict(lambda: None)
    index: dict[str, tuple[str, str, str, str, str, str]] = {}
    with bim.open(encoding="utf-8") as handle:
        for line in handle:
            chrom, snp, cm, bp, a1, a2 = line.rstrip("\n").split("\t")
            chr_counts[chrom] += 1
            bp_int = int(bp)
            chr_min[chrom] = bp_int if chr_min[chrom] is None else min(chr_min[chrom], bp_int)
            chr_max[chrom] = bp_int if chr_max[chrom] is None else max(chr_max[chrom], bp_int)
            index[snp] = (chrom, snp, cm, bp, a1, a2)

    with (out_dir / "round8_chr_distribution.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["chromosome", "snp_count", "min_bp", "max_bp"])
        for chrom in sorted(chr_counts, key=lambda value: int(value) if value.isdigit() else value):
            writer.writerow([chrom, chr_counts[chrom], chr_min[chrom], chr_max[chrom]])

    table2_rows = []
    with table2_tsv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        last_test = ""
        for row in reader:
            if row.get("Test"):
                last_test = row["Test"]
            row["Test"] = last_test
            snp = row.get("Top SNP", "")
            match = index.get(snp)
            table2_rows.append(
                {
                    "test": row.get("Test", ""),
                    "trait": row.get("Trait", ""),
                    "top_snp": snp,
                    "minus_log10p": row.get("−log10P", ""),
                    "present_in_round8": "yes" if match else "no",
                    "round8_chrom": match[0] if match else "",
                    "round8_snp_id": match[1] if match else "",
                    "round8_bp": match[3] if match else "",
                    "round8_a1": match[4] if match else "",
                    "round8_a2": match[5] if match else "",
                }
            )
    fields = [
        "test",
        "trait",
        "top_snp",
        "minus_log10p",
        "present_in_round8",
        "round8_chrom",
        "round8_snp_id",
        "round8_bp",
        "round8_a1",
        "round8_a2",
    ]
    with (out_dir / "table2_snps_in_round8.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(table2_rows)

    return {"bim_snps": sum(chr_counts.values()), "table2_snps": len(table2_rows), "table2_present": sum(row["present_in_round8"] == "yes" for row in table2_rows)}


def audit_samples(prefix: Path, cohort_samples: Path, out_dir: Path) -> dict[str, int]:
    fam_ids = []
    with prefix.with_suffix(".fam").open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 2:
                fam_ids.append(fields[1])
    cohort_ids = [line.strip() for line in cohort_samples.read_text(encoding="utf-8").splitlines() if line.strip()]
    fam_set = set(fam_ids)
    cohort_set = set(cohort_ids)
    missing = sorted(cohort_set - fam_set)
    extra = sorted(fam_set - cohort_set)
    (out_dir / "cohort_oft10424_missing_from_round8_fam.txt").write_text("\n".join(missing) + ("\n" if missing else ""), encoding="utf-8")
    (out_dir / "round8_fam_not_in_cohort_oft10424.txt").write_text("\n".join(extra) + ("\n" if extra else ""), encoding="utf-8")
    with (out_dir / "round8_cohort_overlap_summary.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["round8_fam_samples", len(fam_ids)])
        writer.writerow(["cohort_oft10424_samples", len(cohort_ids)])
        writer.writerow(["cohort_present_in_round8_fam", len(cohort_set & fam_set)])
        writer.writerow(["cohort_missing_from_round8_fam", len(missing)])
        writer.writerow(["round8_fam_not_in_cohort", len(extra)])
    return {
        "fam_samples": len(fam_ids),
        "cohort_samples": len(cohort_ids),
        "cohort_present": len(cohort_set & fam_set),
        "cohort_missing": len(missing),
    }


def main() -> int:
    if len(sys.argv) != 5:
        print("Usage: audit_round8_genotype.py ROUND8_PREFIX TABLE2.tsv COHORT_SAMPLES.txt OUTPUT_DIR", file=sys.stderr)
        return 2

    prefix = Path(sys.argv[1])
    table2_tsv = Path(sys.argv[2])
    cohort_samples = Path(sys.argv[3])
    out_dir = Path(sys.argv[4])
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest = file_manifest(prefix)
    write_manifest(manifest, out_dir)
    bim_stats = audit_bim(prefix.with_suffix(".bim"), table2_tsv, out_dir)
    sample_stats = audit_samples(prefix, cohort_samples, out_dir)

    summary = [
        f"round8_prefix\t{prefix}",
        f"bim_snps\t{bim_stats['bim_snps']}",
        f"fam_samples\t{sample_stats['fam_samples']}",
        f"cohort_oft10424_samples\t{sample_stats['cohort_samples']}",
        f"cohort_present_in_round8_fam\t{sample_stats['cohort_present']}",
        f"cohort_missing_from_round8_fam\t{sample_stats['cohort_missing']}",
        f"table2_snps\t{bim_stats['table2_snps']}",
        f"table2_present_in_round8\t{bim_stats['table2_present']}",
    ]
    (out_dir / "round8_audit_summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("\n".join(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
