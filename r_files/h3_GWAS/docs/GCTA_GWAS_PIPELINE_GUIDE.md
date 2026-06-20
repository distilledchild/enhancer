# GCTA GWAS Pipeline Guide

This guide documents the practical HS rat GWAS workflow used in this project.
The goal is to make the pipeline learnable and resume/CV-defensible: not just
"I ran a command," but "I prepared phenotypes, ran/validated GCTA GWAS, checked
published concordance, and converted SNP-level results for downstream MAGMA and
H-MAGMA analysis."

## Current Evidence Level By Dataset

| Dataset | What was done locally | Best wording |
|---|---|---|
| Kuhn 2025 heroin/nociception | Local GCTA-style reproduction checks for key phenotypes; top SNP, chr/bp, freq, beta, SE, -log10P, N, and MLMA row counts match the paper/report to rounding precision. | Near-exact GWAS reproduction/validation. |
| King 2025 cue-reactivity/sign-tracking | Public GWAS report validation; N, heritability, QTL count, top SNP, and top -log10P matched for selected key phenotypes; summary stats converted and analyzed with cMAGMA/H-MAGMA. | Near-exact report-level validation plus downstream GWAS summary-stat analysis. |
| Lara 2024 delay discounting | Deposited `.mlma` files match the UCSD report and article-level top loci/statistics to rounding precision; MAGMA-ready inputs generated and analyzed. | Deposited GWAS summary-stat provenance confirmed; not a raw genotype GWAS rerun. |
| Gunturkun 2022 open field/social behavior | Local GCTA pipeline assembled for phenotype preparation, GRM, MLMA-LOCO, and MAGMA conversion; exact reproduction is limited by missing age/preprocessing details. | Practical GWAS generation pipeline; not exact paper reproduction yet. |

Use this distinction in CV/resume wording. Do not overstate King/Lara as raw
genotype-level reruns unless the raw phenotype/genotype preprocessing was
actually reproduced end to end.

## Conceptual Pipeline

1. Collect inputs.
   - PLINK genotype reference: `.bed`, `.bim`, `.fam`.
   - Phenotype table: one row per animal/sample.
   - Covariates: sex, batch, cohort, center, coat color, age, etc. if available.
   - Publication/report values for reproducibility checks.

2. Match sample IDs.
   - Match phenotype IDs to `.fam` sample IDs.
   - Create `keep` files for the analyzed sample subset.
   - Confirm sample size `N` against the paper/report.

3. Process phenotype values.
   - Check missingness and outliers.
   - Apply rank inverse-normal transformation if the source pipeline does so.
   - Regress covariates when appropriate.
   - Export a GCTA phenotype file.

4. Build the genomic relationship matrix.
   - Use autosomes only for rat: `--autosome-num 20 --autosome`.
   - Apply MAF filter, usually `--maf 0.01`.
   - Use the same sample subset as the target phenotype or a valid union set.

5. Run GCTA MLMA or MLMA-LOCO.
   - MLMA-LOCO is preferred for genome-wide scans because each chromosome is
     tested while leaving that chromosome out of the GRM.
   - Output is `.mlma` or `.loco.mlma`.

6. Validate GWAS output.
   - Check top SNP, chr, bp, allele frequency, beta, SE, p-value, and -log10P.
   - Check N and row counts.
   - Check accepted SNP count if the report gives it.
   - Tiny rounding differences are expected; direction and locus should match.

7. Convert `.mlma` to downstream input.
   - MAGMA input requires at minimum:

```text
SNP    p
```

   - Keep the exact N used for the phenotype.

8. Run cMAGMA and H-MAGMA.
   - cMAGMA: conventional local/exonic/proximal annotation.
   - H-MAGMA: Hi-C-informed annotation.
   - Strict project version: Duttke/Telese ATAC-filtered, TSS/promoter-overlap
     peaks excluded for noncoding regulatory assignment.

## Local Toolchain

GCTA wrapper:

```text
/Volumes/external_1000GB_all/playground/enhancer/tools/gcta-1.95.2-macOS-arm64/gcta-1.95.2-macOS-arm64/run.sh
```

GCTA binary:

```text
/Volumes/external_1000GB_all/playground/enhancer/tools/gcta-1.95.2-macOS-arm64/gcta-1.95.2-macOS-arm64/bin/gcta64
```

MAGMA binary:

```text
/Volumes/external_1000GB_all/playground/enhancer/tools/magma
```

HS rat genotype reference:

```text
/Volumes/external_1000GB_all/playground/enhancer/data/hs_data/v4/HS_genotypes_v4
```

cMAGMA annotation:

```text
/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/nicotine/hs.exonpro.ONLY.annot
```

Strict H-MAGMA annotation:

```text
/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/nicotine/hs.hic.duttke_telese_atac.non_tss_promoter_filtered.annot
```

MLMA-to-MAGMA converter:

```text
/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022/scripts/mlma_to_magma.py
```

## Example GCTA Flow

The clearest reusable implementation is the Gunturkun pipeline:

```text
/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022/scripts/run_gcta_gwas.sh
```

Core commands:

```bash
THREADS=10 /Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022/scripts/run_gcta_gwas.sh prepare

THREADS=10 /Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022/scripts/run_gcta_gwas.sh full-grm

THREADS=10 /Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022/scripts/run_gcta_gwas.sh full-trait 10424
```

Internally, the full-trait mode does:

```bash
gcta64 \
  --bfile HS_genotypes_v4 \
  --autosome-num 20 \
  --autosome \
  --pheno trait_<id>.pheno \
  --grm target_traits_union_grm \
  --mlma-loco \
  --out trait_<id> \
  --thread-num 10
```

Then it converts:

```bash
python3 mlma_to_magma.py \
  --mlma trait_<id>.loco.mlma \
  --out trait_<id>.magma.tsv
```

## Reproducibility Audit Checklist

For each phenotype, create or inspect a reproduction-check table with:

- `SNP`
- `Chr`
- `bp`
- `Freq`
- `beta`
- `SE`
- `p`
- `-log10(p)`
- `N`
- number of SNP rows in `.mlma`
- accepted SNP count, if reported
- genome-wide threshold, if reported

Interpretation of match quality:

| Result | Meaning |
|---|---|
| All values match to rounding | Strong reproduction/validation. |
| Top SNP/locus matches, rounded stats match, tiny p difference | Acceptable; likely rounding or formatting. |
| Same locus but different top SNP | Partial reproduction; LD or filtering difference possible. |
| Different locus/top SNP | Not reproduced; inspect phenotype preprocessing, covariates, sample subset, build, and filters. |

## Existing Reproducibility Files

Kuhn 2025:

```text
/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Kuhn2025/*_gwas_reproduction_check.csv
/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Kuhn2025/kuhn_total_heroin_reproduction_compact.tsv
```

King 2025:

```text
/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/King2025/reports/selected_1_3_traits_summary.md
/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/King2025/reports/selected_1_3_traits_input_verification.tsv
```

Lara 2024:

```text
/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Lara2024_delaydiscounting/reports/LARA2024_MAGMA_READINESS_REPORT.md
/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Lara2024_delaydiscounting/reports/SELECTED_MAGMA_CORRECTION_REPORT.md
```

Gunturkun 2022:

```text
/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022/README_GWAS_PIPELINE.md
/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022/EXACT_REPRO_AUDIT.md
```

## Learning Plan

### Level 1: Run The Pipeline

Goal: be able to explain every input/output file.

- Identify `.bed/.bim/.fam`, phenotype, covariates, GRM, `.mlma`, and MAGMA
  input files.
- Run one smoke test trait.
- Confirm `.mlma` columns: `Chr SNP bp A1 A2 Freq b se p`.
- Convert `.mlma` to `SNP p`.

### Level 2: Reproduce Or Validate A Published GWAS

Goal: be able to say whether the output agrees with a paper.

- Compare top SNP and top -log10P.
- Compare N and SNP row counts.
- Compare beta/SE/frequency if reported.
- Write a small reproduction table.

### Level 3: Interpret The GWAS

Goal: move from "the file exists" to "I understand the signal."

- How many genome-wide significant loci?
- Are top hits coding/proximal, intergenic, or regulatory-looking?
- Are loci broad or single sharp peaks?
- Is the sample size/power enough for gene-level correction?
- Does the signal make biological sense for the phenotype?

### Level 4: Downstream Gene-Level Analysis

Goal: connect GWAS to cMAGMA/H-MAGMA.

- Run cMAGMA and strict H-MAGMA.
- Apply BH-FDR and Bonferroni.
- Compare H-MAGMA-only genes to cMAGMA.
- Check whether H-MAGMA signal is caused only by more SNPs.
- Inspect candidate loci with GWAS, Hi-C loop, ATAC peak, and gene body.

### Level 5: Resume/CV Narrative

Goal: describe the work without overstating it.

Strong wording:

```text
Reproduced and validated HS rat GWAS analyses using GCTA MLMA/MLMA-LOCO and
SNP-level summary statistics, benchmarking local results against published
top loci, effect sizes, sample sizes, and report-level QTL statistics.
```

Downstream wording:

```text
Integrated SNP-level GWAS summary statistics with cMAGMA and Hi-C-informed
H-MAGMA annotations to prioritize regulatory-contact genes from addiction- and
behavior-relevant rat GWAS datasets.
```

Avoid:

```text
Independently reproduced every GWAS from raw data
```

unless raw phenotype, covariate, genotype, and model settings were fully
reconstructed and checked.
