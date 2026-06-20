# Lara 2024 Selected cMAGMA/H-MAGMA Results

Date: 2026-06-18

Selected phenotypes were limited to the three GWAS traits with significant QTLs in the UCSD GWAS report and paper: `auc`, `exponential_k`, and `dd_indiff_2`.

## Source Links

- GitHub phenotype/pipeline repository: https://github.com/Palmer-Lab-UCSD/HSrat_delaydiscounting
- UCSD Digital Collections GWAS data/report object: https://library.ucsd.edu/dc/object/bb7761581p
- DOI listed for the deposited GWAS data/results/report: https://doi.org/10.6075/J0BG2PDS
- Publication DOI: https://doi.org/10.1111/gbb.12909
- PMID/PMCID: 39119916 / PMC11310854

## Provenance Confirmation

- UCSD GWAS report QTL table: `auc` top SNP `20:32221020`, `dd_indiff_2` top SNP `14:26702994`, `exponential_k` top SNP `20:32221020`.
- Local `.mlma` files reproduce the report rounded values:
  - `auc`: Freq 0.828526, beta 0.360774, se 0.0759683, -log10p 5.689487.
  - `dd_indiff_2`: Freq 0.791057, beta 0.343372, se 0.073712, -log10p 5.496434.
  - `exponential_k`: Freq 0.828526, beta -0.352334, se 0.0760059, -log10p 5.448707.
- PMC article text reports the same loci and rounded -log10(p) values for chromosome 20 and 14.

Verdict: CONFIRMED for downstream use. This is not a raw-genotype GWAS re-run, but the deposited `.mlma` files match the official UCSD GWAS report and publication-level reported loci/statistics to rounding precision.

## Correction Summary

| trait | method | genes | top gene | symbol | top P | Bonf P | BH q | Bonf<0.05 | Bonf<0.10 | FDR<0.05 | FDR<0.10 |
|---|---:|---:|---|---|---:|---:|---:|---:|---:|---:|---:|
| auc | cMAGMA | 18643 | 20:32238545:32239613 |  | 2.4055e-06 | 0.044845737 | 0.044845737 | 1 | 1 | 1 | 1 |
| auc | H-MAGMA_strict | 18963 | ENSRNOG00000066051 |  | 2.4055e-06 | 0.045615497 | 0.045615497 | 1 | 1 | 1 | 1 |
| exponential_k | cMAGMA | 18643 | 20:32238545:32239613 |  | 6.476e-05 | 1 | 0.74022964 | 0 | 0 | 0 | 0 |
| exponential_k | H-MAGMA_strict | 18963 | ENSRNOG00000066051 |  | 6.476e-05 | 1 | 0.7529354 | 0 | 0 | 0 | 0 |
| dd_indiff_2 | cMAGMA | 18643 | 14:26368277:27105860 |  | 0.00016329 | 1 | 0.73626101 | 0 | 0 | 0 | 0 |
| dd_indiff_2 | H-MAGMA_strict | 18963 | ENSRNOG00000030149 | Adgrl3 | 0.00016329 | 1 | 0.72635562 | 0 | 0 | 0 | 0 |

Interpretation:

- `auc` survives both Bonferroni 0.05 and BH FDR 0.05 in cMAGMA and H-MAGMA.
- `exponential_k` and `dd_indiff_2` point to the expected paper loci but do not survive gene-level Bonferroni or BH FDR 0.10 in these MAGMA annotations.
- `dd_indiff_2` H-MAGMA top gene maps to `Adgrl3`, matching the chromosome 14 paper locus.
- The chromosome 20 top result lies in the paper's Slc35f1-associated GWAS locus. The top ENSRNOG identifier from the local annotation did not resolve cleanly to a current Ensembl REST display name, so the coordinate/locus evidence should be used for that locus label.

## Output Files

- `/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Lara2024_delaydiscounting/reports/selected_correction_summary.tsv`
- `/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Lara2024_delaydiscounting/reports/selected_top100_gene_level_with_corrections.tsv`
- `/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Lara2024_delaydiscounting/reports/selected_gene_level_with_corrections.tsv`
- `/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Lara2024_delaydiscounting/reports/selected_candidate_locus_hits_with_corrections.tsv`
