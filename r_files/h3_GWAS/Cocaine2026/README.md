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
```

Generated data and output folders are ignored:

```text
data/
magma_inputs/
results/
logs/
reports/
```
