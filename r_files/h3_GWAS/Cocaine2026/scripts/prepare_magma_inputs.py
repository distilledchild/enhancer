#!/usr/bin/env python3
"""Extract selected cocaine GWAS MLMA files and build MAGMA p-value inputs."""

from __future__ import annotations

import argparse
import csv
import math
import re
import zipfile
from pathlib import Path


DEFAULT_TRAITS = [
    "regressedlr_sha_mean_to_01_03",
    "regressedlr_pc1_lga",
    "regressedlr_lga_total_intake",
    "regressedlr_shock_03_calculated",
]


def trait_member_pattern(trait: str) -> re.Pattern[str]:
    return re.compile(rf"^gwas/{re.escape(trait)}_chrgwas(.+)\.mlma$")


def chromosome_sort_key(chromosome: str) -> tuple[int, int | str]:
    if chromosome.isdigit():
        return (0, int(chromosome))
    if chromosome == "x":
        return (1, chromosome)
    if chromosome == "mt":
        return (2, chromosome)
    return (3, chromosome)


def iter_trait_rows(zip_path: Path, trait: str):
    pattern = trait_member_pattern(trait)
    with zipfile.ZipFile(zip_path) as archive:
        members: list[tuple[str, str]] = []
        for name in archive.namelist():
            match = pattern.match(name)
            if match:
                members.append((match.group(1), name))
        if not members:
            raise SystemExit(f"No MLMA members found for {trait}")
        members.sort(key=lambda item: chromosome_sort_key(item[0]))
        for chromosome, member in members:
            with archive.open(member) as handle:
                text = (line.decode("utf-8") for line in handle)
                reader = csv.DictReader(text, delimiter="\t")
                if reader.fieldnames is None:
                    raise SystemExit(f"No header found in {member}")
                for row in reader:
                    yield chromosome, member, row


def write_trait_input(zip_path: Path, trait: str, out_dir: Path) -> tuple[int, int]:
    out_path = out_dir / f"{trait}.magma_pval.tsv"
    seen: set[str] = set()
    kept = 0
    skipped = 0
    with out_path.open("w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["SNP", "p"])
        for _chromosome, _member, row in iter_trait_rows(zip_path, trait):
            snp = row.get("SNP", "")
            p_raw = row.get("p", row.get("P", ""))
            if not snp or snp in seen or not p_raw or p_raw == "NA":
                skipped += 1
                continue
            try:
                p_value = float(p_raw)
            except ValueError:
                skipped += 1
                continue
            if not math.isfinite(p_value) or p_value <= 0 or p_value > 1:
                skipped += 1
                continue
            writer.writerow([snp, p_raw])
            seen.add(snp)
            kept += 1
    return kept, skipped


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--zip", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--traits", nargs="*", default=DEFAULT_TRAITS)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for trait in args.traits:
        kept, skipped = write_trait_input(args.zip, trait, args.out_dir)
        print(f"{trait}\tkept={kept}\tskipped={skipped}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
