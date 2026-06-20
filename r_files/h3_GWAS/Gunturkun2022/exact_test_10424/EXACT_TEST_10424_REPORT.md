# Exact-Public 10424 Sanity Test

Date: 2026-06-18

## Inputs

Genotype:

`/Volumes/external_1000GB_all/playground/enhancer/data/King_2025/bb15123938_2_1/round8_unpruned`

Phenotype:

`10424 open_field_totaldistance`

Reconstructed cohort:

`/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022/data/repro_exact/genenetwork/cohort_oft_totaldistance_10424_n1246.samples.txt`

Prepared phenotype:

`data/trait_10424.exact_public.pheno`

## Public-Data Preprocessing

The available public data allowed:

1. sex-stratified rank inverse-normal transform,
2. per-covariate screen for `batch` and `color`,
3. residualization of selected covariates,
4. final rank inverse-normal transform.

The paper also mentions age, but the per-sample age covariate was not found in the public GeneNetwork metadata.

Selected covariates:

- `batch`: p = `1.33027205100112e-33`, R2 = `0.142547370677957`
- `color`: p = `1.89537483190892e-24`, R2 = `0.0902133645869209`

N:

- total: `1246`
- female: `620`
- male: `626`

## GCTA Run

GRM runtime:

- `22.39 sec`

MLMA-LOCO runtime:

- `4 minutes 18 sec`

Output:

`results/trait_10424.exact_public.loco.mlma`

Valid SNP rows:

`3,513,489`

## Paper Peak Comparison

| SNP | Paper -log10P | Observed -log10P | Delta |
| --- | ---: | ---: | ---: |
| `chr10:94549701` | `7.286` | `5.4242281` | `-1.8617719` |
| `chr11:33359859` | `8.268` | `6.93788829` | `-1.33011171` |

## Chromosome Tops

Only chromosomes 10 and 11 crossed the paper threshold `-log10P > 5.609`.

| Chr | Top SNP | -log10P |
| --- | --- | ---: |
| 10 | `chr10:95040064` | `5.69526826` |
| 11 | `chr11:33359859` | `6.93788829` |

## Interpretation

This run partially reproduces Table 2:

- The chr11 peak is at the exact paper SNP and is the genome-wide top locus.
- The chr10 locus is present and crosses the paper threshold, but the local top SNP is about `490 kb` from the paper peak.
- Both loci are weaker than the paper-reported p-values.

Therefore this is not yet a perfect reproduction. The remaining likely causes are missing per-sample age covariate, exact author-side covariate selection details, or differences between public GeneNetwork-processed values and the C-GORD analysis input.
