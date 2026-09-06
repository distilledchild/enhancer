# Cocaine2026 GWAS Summary-Stat Analysis

This folder stores lightweight code and documentation for the cocaine
self-administration GWAS follow-up analysis. Large UCSD downloads, MAGMA
inputs, logs, and result tables are intentionally kept on the external SSD and
ignored by git.

## Source

- Dataset: UCSD Digital Collections object `bb2334903t`
- Data DOI: `10.6075/J0QN675X`
- Study: genome-wide association study of cocaine self-administration behavior
  in Heterogeneous Stock rats.

## Local External Work Root

```text
/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Cocaine2026
```

## What Was Run

The four-trait MAGMA analysis below used the previous annotation; it is not a
MAGMA rerun with the revised 2001-bp loop set.

The public GWAS result archive contained chromosome-split GCTA `.mlma` files
for cocaine-related phenotypes. Four key addiction-relevant phenotypes were
selected and converted to MAGMA p-value input:

- `regressedlr_sha_mean_to_01_03` - short-access early timeout
- `regressedlr_pc1_lga` - long-access PC1 addiction-like behavior
- `regressedlr_lga_total_intake` - long-access total intake
- `regressedlr_shock_03_calculated` - shock compulsivity

For each phenotype, both cMAGMA and strict H-MAGMA were run using the project
annotations:

- cMAGMA: `hs.exonpro.ONLY.annot`
- strict H-MAGMA: Duttke/Telese ATAC-filtered Hi-C annotation with
  TSS/promoter-overlap peaks excluded

## Result Summary

Across the four selected phenotypes, no cMAGMA or H-MAGMA gene survived
BH-FDR or Bonferroni correction. The dataset is therefore useful as a
negative/limited-result screen and as a public summary-stat workflow example,
but it is not currently a main discovery dataset for H-MAGMA novelty.

## Lightweight Scripts Kept In Git

```text
scripts/prepare_magma_inputs.py
scripts/summarize_magma_pair.py
scripts/scan_mlma_peaks.awk
scripts/screen_extended_traits.R
scripts/validate_extended_trait_screen.R
```

Generated data and output folders are ignored:

```text
data/
magma_inputs/
results/
logs/
reports/
```

## Expanded Trait Screen (2026-09-06)

All 29 archived traits (638 MLMA files) were screened at the user-requested
lead-SNP threshold of `-log10(P) >= 7`. **None passed**; the maximum was 6.249806
for `regressedlr_lga_mean_11_14`. No new H-MAGMA run was launched and the
threshold was not relaxed. All-row P-value checks, independent R re-reading
of the strongest trait and both PR traits, and seven source-report lead-SNP
checks passed.

The archive does not contain the proposed sensitization/extinction/reinstatement
traits. Its two PR variables measure active lever presses, not breakpoint.
Threshold 7 is stricter than the paper's 5% threshold of 5.58; this screen does
not establish that gene-level MAGMA results would all be nonsignificant.

- [Review report and reproduction steps](/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/Cocaine2026/EXTENDED_TRAIT_SCREEN_REPORT.md)
- [All 29 trait results](/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Cocaine2026/extended_trait_screen_2001bp/all_29_trait_peak_screen.tsv)

Generated screen outputs reside in
`/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Cocaine2026/extended_trait_screen_2001bp/`.
The directory name refers to the intended downstream loop context; this
genome-wide screen itself did not filter SNPs using the loop annotation.
