# Kuhn2025 Heroin/Nociception GWAS Reproduction And MAGMA Analysis

This folder stores lightweight documentation for the Kuhn et al. 2025 heroin
vulnerability/nociception GWAS reproduction and downstream cMAGMA/H-MAGMA
analysis. Large summary-stat, MAGMA, and result files are kept on the external
SSD and ignored by git.

## Source

- Study: Kuhn et al. 2025, Molecular Psychiatry
- Title: Genome-wide association study reveals multiple loci for nociception
  and opioid consumption behaviors associated with heroin vulnerability in
  outbred rats
- DOI: `10.1038/s41380-025-02922-4`
- PMID: `40000848`

## Local External Work Root

```text
/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Kuhn2025
```

## What Was Reproduced / Validated

Local GWAS reproduction checks were generated for key phenotypes:

- total heroin consumption
- break point (`bp`)
- escalation of heroin intake over 12h
- baseline tail flick (`tail_flick_bl1`)

For these checks, the local outputs were compared to the paper/report values:

- top SNP
- chromosome and base-pair position
- allele frequency
- beta
- standard error
- p-value / `-log10(p)`
- sample size
- MLMA row count

The checked phenotypes matched the reported values to rounding precision, with
only tiny p-value or `-log10(p)` differences where expected from formatting.

## Downstream Analysis

For each selected phenotype, MAGMA p-value inputs were generated and analyzed
with:

- cMAGMA annotation: `hs.exonpro.ONLY.annot`
- strict H-MAGMA annotation: Duttke/Telese ATAC-filtered Hi-C annotation with
  TSS/promoter-overlap peaks excluded

The current downstream gene-level results do not provide a strong H-MAGMA-only
novelty story, but the dataset is a strong GWAS reproduction/validation example.

## Key External Result Files

```text
kuhn_total_heroin_reproduction_compact.tsv
kuhn_*_gwas_reproduction_check.csv
kuhn_*_magma_hmagma_summary.csv
regressedlr_*_cmagma.genes.out
regressedlr_*_hmagma_strict_duttke_telese.genes.out
```

These result files are intentionally ignored by git.
