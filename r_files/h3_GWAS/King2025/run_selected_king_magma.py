#!/usr/bin/env python3
from __future__ import annotations

import math
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


TRAITS = [
    "crf_ny_incentive_value_index",
    "crf_ny_lever_presses",
    "pavca_ny_d5_response_bias",
]

DATA_ROOT = Path("/Volumes/external_1000GB_all/playground/enhancer/data/King_2025")
OUT_ROOT = Path("/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/King2025")
LOCAL_ROOT = Path("/Users/pete/Desktop/playground/enhancer")

RAW_MLMA_ROOT = DATA_ROOT / "gwas_summary_stats/raw_chrgwas_mlma"
TRAIT_N_PATH = DATA_ROOT / "gwas_summary_stats/trait_n_from_heritability.tsv"
QTL_PATH = DATA_ROOT / "ucsd_bb5030313v/s3_downloads/results__qtls__finalqtl.csv"
INTEGRITY_PATH = DATA_ROOT / "gwas_summary_stats/king_download_integrity_summary.tsv"

MAGMA = Path("/Volumes/external_1000GB_all/playground/enhancer/tools/magma")
BFILE = Path("/Volumes/external_1000GB_all/playground/enhancer/data/hs_data/v4/HS_genotypes_v4")
ANNOT_CMAGMA = LOCAL_ROOT / "r_files/h3_GWAS/nicotine/hs.exonpro.ONLY.annot"
ANNOT_HMAGMA = LOCAL_ROOT / "r_files/h3_GWAS/nicotine/hs.hic.duttke_telese_atac.non_tss_promoter_filtered.annot"
COUNTS_PATH = LOCAL_ROOT / "r_files/h3_GWAS/nicotine/cmagma_vs_hmagma_duttke_telese_non_tss_promoter_snp_counts.csv"

MAP_SOURCE_RN7_BED = OUT_ROOT / "pavca_ny_d5_index.rn7.bed"
MAP_PATH = DATA_ROOT / "gwas_summary_stats/king_round8_to_rn7_snp_map.tsv"

LOG_DIR = OUT_ROOT / "logs"
REPORT_DIR = OUT_ROOT / "reports"


def log(message: str) -> None:
    print(time.strftime("[%Y-%m-%d %H:%M:%S]"), message, flush=True)


def bh_adjust(pvals: pd.Series | np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(adj, 1.0)
    return out


def run_command(cmd: list[str], stdout_path: Path | None = None) -> None:
    log("RUN " + " ".join(cmd))
    t0 = time.time()
    if stdout_path:
        stdout_path.parent.mkdir(parents=True, exist_ok=True)
        with stdout_path.open("w") as out:
            proc = subprocess.run(cmd, stdout=out, stderr=subprocess.STDOUT, text=True)
    else:
        proc = subprocess.run(cmd, text=True)
    elapsed = time.time() - t0
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed with status {proc.returncode}: {' '.join(cmd)}")
    log(f"DONE elapsed={elapsed:.1f}s")


def count_lines(path: Path) -> int:
    with path.open() as handle:
        return sum(1 for _ in handle)


def ensure_round8_to_rn7_map() -> None:
    if MAP_PATH.exists():
        log(f"Using existing SNP map: {MAP_PATH}")
        return
    if not MAP_SOURCE_RN7_BED.exists():
        raise FileNotFoundError(
            f"Missing {MAP_SOURCE_RN7_BED}. Run the pavca_ny_d5_index liftOver first."
        )
    log(f"Building SNP map from {MAP_SOURCE_RN7_BED}")
    tmp = MAP_PATH.with_suffix(".tmp")
    rows = 0
    with MAP_SOURCE_RN7_BED.open() as inp, tmp.open("w") as out:
        out.write("original_snp\trn7_snp\n")
        for line in inp:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, _start, end, original_snp, *_ = line.rstrip("\n").split("\t")
            rn7_snp = f"{chrom.removeprefix('chr')}:{end}"
            out.write(f"{original_snp}\t{rn7_snp}\n")
            rows += 1
    tmp.replace(MAP_PATH)
    log(f"Wrote {MAP_PATH} rows={rows:,}")


def load_snp_map() -> dict[str, str]:
    ensure_round8_to_rn7_map()
    log("Loading SNP map into memory")
    mapping: dict[str, str] = {}
    with MAP_PATH.open() as handle:
        next(handle)
        for line in handle:
            old, new = line.rstrip("\n").split("\t")
            mapping[old] = new
    log(f"Loaded SNP map entries={len(mapping):,}")
    return mapping


def verify_trait_source(trait: str, trait_n: pd.DataFrame, qtls: pd.DataFrame) -> dict[str, object]:
    rows = 0
    files = 0
    min_p = math.inf
    min_snp = None
    for chrom in range(1, 21):
        path = RAW_MLMA_ROOT / trait / f"regressedlr_{trait}_chrgwas{chrom}.mlma"
        if not path.exists():
            raise FileNotFoundError(path)
        files += 1
        for chunk in pd.read_csv(path, sep="\t", usecols=["SNP", "p"], chunksize=500_000):
            rows += len(chunk)
            idx = chunk["p"].idxmin()
            p = float(chunk.loc[idx, "p"])
            if p < min_p:
                min_p = p
                min_snp = str(chunk.loc[idx, "SNP"])
    tn = trait_n.loc[trait_n["trait"].eq(trait)].iloc[0]
    qtrait = qtls.loc[qtls["trait"].eq(trait)].copy()
    qtop = qtrait.sort_values("p", ascending=False).iloc[0].to_dict() if len(qtrait) else {}
    return {
        "trait": trait,
        "N": int(float(tn["n"])),
        "heritability": float(tn["heritability"]),
        "h2_se": float(tn["h2_se"]),
        "raw_mlma_files": files,
        "raw_mlma_rows": rows,
        "top_snp_from_mlma": min_snp,
        "top_p_from_mlma": min_p,
        "top_minus_log10p_from_mlma": -math.log10(min_p),
        "public_report_top_snp": qtop.get("SNP"),
        "public_report_top_minus_log10p": qtop.get("p"),
        "public_report_qtls": len(qtrait),
        "public_report_5pct_qtls": int((qtrait["significance_level"].astype(str) == "5%").sum()) if len(qtrait) else 0,
        "public_report_10pct_qtls": int((qtrait["significance_level"].astype(str) == "10%").sum()) if len(qtrait) else 0,
        "top_snp_match_report": bool(qtop and qtop.get("SNP") == min_snp),
    }


def create_rn7_pval(trait: str, snp_map: dict[str, str]) -> dict[str, object]:
    out = OUT_ROOT / f"{trait}_rn7.magma.tsv"
    tmp = out.with_suffix(".tmp")
    values: dict[str, float] = {}
    rows = 0
    mapped_rows = 0
    duplicate_rn7 = 0
    for chrom in range(1, 21):
        path = RAW_MLMA_ROOT / trait / f"regressedlr_{trait}_chrgwas{chrom}.mlma"
        for chunk in pd.read_csv(path, sep="\t", usecols=["SNP", "p"], chunksize=500_000):
            rows += len(chunk)
            for original_snp, pval in zip(chunk["SNP"].astype(str), chunk["p"]):
                rn7_snp = snp_map.get(original_snp)
                if rn7_snp is None or pd.isna(pval):
                    continue
                mapped_rows += 1
                p = float(pval)
                prev = values.get(rn7_snp)
                if prev is None:
                    values[rn7_snp] = p
                else:
                    duplicate_rn7 += 1
                    if p < prev:
                        values[rn7_snp] = p
    with tmp.open("w") as handle:
        handle.write("SNP\tp\n")
        for snp, p in values.items():
            handle.write(f"{snp}\t{p:.12g}\n")
    tmp.replace(out)
    return {
        "trait": trait,
        "raw_rows": rows,
        "mapped_rows": mapped_rows,
        "unique_rn7_snp_ids": len(values),
        "duplicate_rn7_ids_collapsed": duplicate_rn7,
        "pval_path": str(out),
    }


def parse_magma_log(path: Path) -> dict[str, object]:
    text = path.read_text(errors="ignore")
    valid_match = re.search(
        r"valid SNP p-values for\s+(\d+)\s+SNPs in data\s+\(([^)]+)\)",
        text,
    )
    gene_match = re.search(r"found\s+(\d+)\s+genes containing valid SNPs", text)
    return {
        "valid_snp_pvalues_in_data": int(valid_match.group(1)) if valid_match else None,
        "valid_snp_pvalue_percent": valid_match.group(2) if valid_match else None,
        "genes_containing_valid_snps": int(gene_match.group(1)) if gene_match else None,
    }


def run_magma_for_trait(trait: str, n: int) -> dict[str, object]:
    pval = OUT_ROOT / f"{trait}_rn7.magma.tsv"
    outputs = {}
    for method, annot, suffix in [
        ("cMAGMA", ANNOT_CMAGMA, "cmagma"),
        ("strict H-MAGMA", ANNOT_HMAGMA, "hmagma_strict_duttke_telese"),
    ]:
        out_prefix = OUT_ROOT / f"{trait}_{suffix}"
        log_path = LOG_DIR / f"{trait}_{suffix}.stdout.log"
        for extra in [".genes.out", ".genes.raw", ".log"]:
            target = Path(str(out_prefix) + extra)
            if target.exists():
                target.unlink()
        if log_path.exists():
            log_path.unlink()
        cmd = [
            str(MAGMA),
            "--bfile",
            str(BFILE),
            "--pval",
            str(pval),
            "use=SNP,p",
            f"N={n}",
            "--gene-annot",
            str(annot),
            "--out",
            str(out_prefix),
        ]
        run_command(cmd, log_path)
        outputs[method] = {
            "genes_out": str(Path(str(out_prefix) + ".genes.out")),
            "log": str(log_path),
            **parse_magma_log(log_path),
        }
    return outputs


def summarize_trait(
    trait: str,
    source_info: dict[str, object],
    pval_info: dict[str, object],
    magma_info: dict[str, object],
    counts: pd.DataFrame,
    integrity: dict[str, str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    cm = pd.read_csv(OUT_ROOT / f"{trait}_cmagma.genes.out", sep=r"\s+")
    hm = pd.read_csv(OUT_ROOT / f"{trait}_hmagma_strict_duttke_telese.genes.out", sep=r"\s+")
    for df in (cm, hm):
        df["BH_FDR"] = bh_adjust(df["P"])
        df["coord_key"] = df["CHR"].astype(str) + ":" + df["START"].astype(str) + ":" + df["STOP"].astype(str)
        n = len(df)
        df["BONF_0_05"] = df["P"] < 0.05 / n
        df["BONF_0_10"] = df["P"] < 0.10 / n
        df["BH_0_05"] = df["BH_FDR"] < 0.05
        df["BH_0_10"] = df["BH_FDR"] < 0.10

    hm = hm.merge(
        counts[
            [
                "ensg",
                "gene_name",
                "n_snp_cmagma",
                "n_snp_hmagma_filtered",
                "added_snps_filtered",
                "noise_snps_removed",
            ]
        ],
        how="left",
        left_on="GENE",
        right_on="ensg",
    )
    c_for_join = cm[
        [
            "coord_key",
            "GENE",
            "NSNPS",
            "NPARAM",
            "ZSTAT",
            "P",
            "BH_FDR",
            "BONF_0_05",
            "BONF_0_10",
            "BH_0_05",
            "BH_0_10",
        ]
    ].rename(
        columns={
            "GENE": "cmagma_gene_id",
            "NSNPS": "cmagma_NSNP",
            "NPARAM": "cmagma_NPARAM",
            "ZSTAT": "cmagma_ZSTAT",
            "P": "cmagma_P",
            "BH_FDR": "cmagma_BH_FDR",
            "BONF_0_05": "cmagma_BONF_0_05",
            "BONF_0_10": "cmagma_BONF_0_10",
            "BH_0_05": "cmagma_BH_0_05",
            "BH_0_10": "cmagma_BH_0_10",
        }
    )
    h_for_join = hm[
        [
            "coord_key",
            "GENE",
            "gene_name",
            "CHR",
            "START",
            "STOP",
            "NSNPS",
            "NPARAM",
            "N",
            "ZSTAT",
            "P",
            "BH_FDR",
            "BONF_0_05",
            "BONF_0_10",
            "BH_0_05",
            "BH_0_10",
            "n_snp_cmagma",
            "n_snp_hmagma_filtered",
            "added_snps_filtered",
            "noise_snps_removed",
        ]
    ].rename(
        columns={
            "GENE": "hmagma_gene_id",
            "NSNPS": "hmagma_NSNP",
            "NPARAM": "hmagma_NPARAM",
            "ZSTAT": "hmagma_ZSTAT",
            "P": "hmagma_P",
            "BH_FDR": "hmagma_BH_FDR",
            "BONF_0_05": "hmagma_BONF_0_05",
            "BONF_0_10": "hmagma_BONF_0_10",
            "BH_0_05": "hmagma_BH_0_05",
            "BH_0_10": "hmagma_BH_0_10",
        }
    )
    comp = h_for_join.merge(c_for_join, how="left", on="coord_key")
    comp["delta_NSNP_magma_out"] = comp["hmagma_NSNP"] - comp["cmagma_NSNP"]
    comp["delta_ZSTAT"] = comp["hmagma_ZSTAT"] - comp["cmagma_ZSTAT"]
    comp["h_minus_c_log10P"] = (-np.log10(comp["hmagma_P"])) - (-np.log10(comp["cmagma_P"]))
    for label in ["BONF_0_05", "BONF_0_10", "BH_0_05", "BH_0_10"]:
        comp[f"hmagma_only_{label}"] = comp[f"hmagma_{label}"].eq(True) & ~comp[f"cmagma_{label}"].eq(True)

    summary_rows = []
    for method, df in [("cMAGMA", cm), ("strict H-MAGMA", hm)]:
        n = len(df)
        for criterion, threshold, mask in [
            ("tested_genes", np.nan, np.ones(n, dtype=bool)),
            ("Bonferroni P < 0.05/n", 0.05 / n, df["P"] < 0.05 / n),
            ("Bonferroni P < 0.10/n", 0.10 / n, df["P"] < 0.10 / n),
            ("BH-FDR < 0.05", 0.05, df["BH_FDR"] < 0.05),
            ("BH-FDR < 0.10", 0.10, df["BH_FDR"] < 0.10),
        ]:
            summary_rows.append(
                {
                    "trait": trait,
                    "criterion": criterion,
                    "method": method,
                    "n_genes": n,
                    "threshold": threshold,
                    "n_significant": int(mask.sum()),
                }
            )
    for criterion, col in [
        ("Bonferroni P < 0.05/n", "hmagma_only_BONF_0_05"),
        ("Bonferroni P < 0.10/n", "hmagma_only_BONF_0_10"),
        ("BH-FDR < 0.05", "hmagma_only_BH_0_05"),
        ("BH-FDR < 0.10", "hmagma_only_BH_0_10"),
    ]:
        summary_rows.append(
            {
                "trait": trait,
                "criterion": criterion,
                "method": "H-MAGMA-only vs cMAGMA coordinate match",
                "n_genes": len(comp),
                "threshold": np.nan,
                "n_significant": int(comp[col].sum()),
            }
        )
    summary = pd.DataFrame(summary_rows)

    sig_cols = [
        "hmagma_gene_id",
        "gene_name",
        "CHR",
        "START",
        "STOP",
        "hmagma_NSNP",
        "cmagma_NSNP",
        "delta_NSNP_magma_out",
        "hmagma_ZSTAT",
        "cmagma_ZSTAT",
        "delta_ZSTAT",
        "hmagma_P",
        "hmagma_BH_FDR",
        "cmagma_P",
        "cmagma_BH_FDR",
        "n_snp_cmagma",
        "n_snp_hmagma_filtered",
        "added_snps_filtered",
        "noise_snps_removed",
        "hmagma_BONF_0_05",
        "cmagma_BONF_0_05",
        "hmagma_BH_0_05",
        "cmagma_BH_0_05",
        "hmagma_BH_0_10",
        "cmagma_BH_0_10",
    ]
    comp_sorted = comp.sort_values(["hmagma_P", "hmagma_gene_id"])
    top_h = comp_sorted.head(50).copy()
    comp.to_csv(REPORT_DIR / f"{trait}_all_gene_comparison.tsv", sep="\t", index=False)
    cm.to_csv(REPORT_DIR / f"{trait}_cmagma_with_bh.tsv", sep="\t", index=False)
    hm.to_csv(REPORT_DIR / f"{trait}_hmagma_with_bh.tsv", sep="\t", index=False)
    top_h[sig_cols].to_csv(REPORT_DIR / f"{trait}_top50_hmagma_with_cmagma_nsnp.tsv", sep="\t", index=False)

    honly_key_frames = []
    for name, col in [
        ("bonferroni05", "hmagma_only_BONF_0_05"),
        ("bonferroni10", "hmagma_only_BONF_0_10"),
        ("bh_fdr05", "hmagma_only_BH_0_05"),
        ("bh_fdr10", "hmagma_only_BH_0_10"),
    ]:
        df = comp_sorted.loc[comp_sorted[col]].copy()
        df[sig_cols].to_csv(REPORT_DIR / f"{trait}_hmagma_only_{name}.tsv", sep="\t", index=False)
        tmp = df[sig_cols].copy()
        tmp.insert(0, "criterion", name)
        tmp.insert(0, "trait", trait)
        honly_key_frames.append(tmp)
    honly_key = pd.concat(honly_key_frames, ignore_index=True) if honly_key_frames else pd.DataFrame()

    md = []
    md.append(f"# King 2025 {trait} cMAGMA/H-MAGMA report")
    md.append("")
    md.append("## Input verification")
    md.append(f"- Public King MLMA files: 20 autosomes, {source_info['raw_mlma_rows']:,} SNP rows.")
    md.append(
        f"- Top SNP check: MLMA top `{source_info['top_snp_from_mlma']}`, "
        f"-log10(p) {source_info['top_minus_log10p_from_mlma']:.6g}; "
        f"public report top `{source_info['public_report_top_snp']}`, "
        f"-log10(p) {source_info['public_report_top_minus_log10p']:.6g}; "
        f"match = {source_info['top_snp_match_report']}."
    )
    md.append(
        f"- Trait N used for MAGMA: {source_info['N']}; "
        f"h2={source_info['heritability']:.6g}, SE={source_info['h2_se']:.6g}."
    )
    md.append(
        f"- Public QTLs: total={source_info['public_report_qtls']}, "
        f"5%={source_info['public_report_5pct_qtls']}, "
        f"10%={source_info['public_report_10pct_qtls']}."
    )
    md.append(
        f"- rn7 p-value mapping: raw rows={pval_info['raw_rows']:,}, "
        f"mapped rows={pval_info['mapped_rows']:,}, "
        f"unique rn7 SNP IDs={pval_info['unique_rn7_snp_ids']:,}, "
        f"duplicates collapsed={pval_info['duplicate_rn7_ids_collapsed']}."
    )
    md.append(
        f"- MAGMA matched SNP p-values: cMAGMA={magma_info['cMAGMA']['valid_snp_pvalues_in_data']:,} "
        f"({magma_info['cMAGMA']['valid_snp_pvalue_percent']}), "
        f"H-MAGMA={magma_info['strict H-MAGMA']['valid_snp_pvalues_in_data']:,} "
        f"({magma_info['strict H-MAGMA']['valid_snp_pvalue_percent']})."
    )
    md.append("")
    md.append("## Correction summary")
    md.append(summary.to_markdown(index=False))
    md.append("")
    md.append("## H-MAGMA-only genes")
    for label, name in [
        ("Bonferroni 0.05", "bonferroni05"),
        ("Bonferroni 0.10", "bonferroni10"),
        ("BH-FDR < 0.05", "bh_fdr05"),
        ("BH-FDR < 0.10", "bh_fdr10"),
    ]:
        df = pd.read_csv(REPORT_DIR / f"{trait}_hmagma_only_{name}.tsv", sep="\t")
        md.append(f"### {label}")
        md.append(df.to_markdown(index=False) if len(df) else "None.")
        md.append("")
    md.append("## Top H-MAGMA genes with cMAGMA NSNP/P")
    md.append(top_h[sig_cols].head(25).to_markdown(index=False))
    (REPORT_DIR / f"{trait}_magma_report.md").write_text("\n".join(md) + "\n")
    return summary, honly_key


def main() -> int:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    trait_n = pd.read_csv(TRAIT_N_PATH, sep="\t")
    qtls = pd.read_csv(QTL_PATH)
    counts = pd.read_csv(COUNTS_PATH)
    integrity_df = pd.read_csv(INTEGRITY_PATH, sep="\t")
    integrity = dict(zip(integrity_df.iloc[:, 0], integrity_df.iloc[:, 1]))
    snp_map = load_snp_map()

    all_summary = []
    all_honly = []
    run_meta = []
    for trait in TRAITS:
        log(f"=== {trait} ===")
        source_info = verify_trait_source(trait, trait_n, qtls)
        log(
            f"Verified {trait}: rows={source_info['raw_mlma_rows']:,}, "
            f"top={source_info['top_snp_from_mlma']} match={source_info['top_snp_match_report']}"
        )
        pval_info = create_rn7_pval(trait, snp_map)
        log(
            f"Created rn7 pval for {trait}: unique={pval_info['unique_rn7_snp_ids']:,}, "
            f"mapped={pval_info['mapped_rows']:,}"
        )
        magma_info = run_magma_for_trait(trait, int(source_info["N"]))
        summary, honly = summarize_trait(trait, source_info, pval_info, magma_info, counts, integrity)
        all_summary.append(summary)
        all_honly.append(honly)
        run_meta.append({**source_info, **pval_info})
        log(f"Finished {trait}")

    summary_all = pd.concat(all_summary, ignore_index=True)
    honly_all = pd.concat(all_honly, ignore_index=True) if all_honly else pd.DataFrame()
    meta_all = pd.DataFrame(run_meta)
    summary_all.to_csv(REPORT_DIR / "selected_1_3_traits_magma_correction_summary.tsv", sep="\t", index=False)
    honly_all.to_csv(REPORT_DIR / "selected_1_3_traits_hmagma_only_key_genes.tsv", sep="\t", index=False)
    meta_all.to_csv(REPORT_DIR / "selected_1_3_traits_input_verification.tsv", sep="\t", index=False)

    md = ["# King 2025 Selected Phenotypes cMAGMA/H-MAGMA Summary", ""]
    md.append("## Traits")
    md.append(meta_all[["trait", "N", "heritability", "h2_se", "public_report_qtls", "public_report_5pct_qtls", "public_report_10pct_qtls", "top_snp_from_mlma", "top_minus_log10p_from_mlma"]].to_markdown(index=False))
    md.append("")
    md.append("## Correction Summary")
    md.append(summary_all.to_markdown(index=False))
    md.append("")
    md.append("## H-MAGMA-only Key Genes")
    if len(honly_all):
        md.append(honly_all.to_markdown(index=False))
    else:
        md.append("None.")
    (REPORT_DIR / "selected_1_3_traits_summary.md").write_text("\n".join(md) + "\n")
    log(f"Wrote combined report: {REPORT_DIR / 'selected_1_3_traits_summary.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
