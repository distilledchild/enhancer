# GWAS Result Interpretation Guide

This guide explains how to analyze GWAS and cMAGMA/H-MAGMA outputs after the
pipeline has run. The purpose is to build real analytical fluency for CV,
resume, fellowship, and paper discussions.

## First Question: What Kind Of Result Do We Have?

After a GWAS or summary-stat analysis, classify the result into one of four
groups.

| Type | Meaning | Example interpretation |
|---|---|---|
| Strong GWAS + strong gene-level result | The GWAS locus is strong and MAGMA/H-MAGMA identifies corrected genes. | Good candidate for a paper figure. |
| Strong GWAS + weak gene-level result | SNP-level signal exists, but gene aggregation does not survive correction. | Useful GWAS, but weak for H-MAGMA novelty. |
| Weak GWAS + weak gene-level result | Top SNPs and top genes are only nominal. | Not a main discovery dataset. |
| H-MAGMA-only corrected genes | H-MAGMA identifies corrected genes missed by cMAGMA. | Core evidence for the project. |

## GWAS-Level Interpretation

Inspect the `.mlma` file first.

Typical GCTA MLMA columns:

```text
Chr SNP bp A1 A2 Freq b se p
```

### What To Check

1. Sample size.
   - Does N match the paper/report?
   - If N is small, expect weak power.

2. Top SNP.
   - What is the smallest p-value?
   - What is `-log10(p)`?
   - Does it pass the paper's genome-wide threshold?

3. Locus structure.
   - Is the top SNP in a coding gene?
   - Is it intergenic/noncoding?
   - Is the signal a broad LD block or a sharp peak?

4. Effect direction and uncertainty.
   - Is beta positive or negative?
   - Is SE reasonable?
   - Is frequency extreme?

5. Biological fit.
   - Is the phenotype addiction-, reward-, impulsivity-, anxiety-, or cortical-control relevant?
   - Does the locus contain plausible genes?

## Reproducibility Interpretation

Use the paper/report as a benchmark.

High-confidence reproduction:

- same top SNP
- same chr/bp
- same N
- same frequency after rounding
- same beta/SE after rounding
- same -log10P after rounding
- same or very similar p-value
- same MLMA row count or accepted SNP count

Kuhn is currently in this category for the checked phenotypes.

Report-level validation:

- same public report N
- same heritability
- same QTL count
- same top SNP
- same top -log10P
- downstream summary stats usable for cMAGMA/H-MAGMA

King is currently in this category.

Deposited summary-stat provenance confirmation:

- official `.mlma` files match UCSD report and article values
- not necessarily a raw genotype/phenotype rerun

Lara is currently in this category.

## MAGMA / H-MAGMA Interpretation

MAGMA output:

```text
*.genes.out
```

Important columns:

| Column | Meaning |
|---|---|
| `GENE` | Gene ID or coordinate ID. |
| `CHR`, `START`, `STOP` | Gene genomic interval. |
| `NSNPS` | Number of SNPs assigned to that gene. |
| `N` | Sample size used. |
| `ZSTAT` | Gene-level Z statistic. |
| `P` | Gene-level p-value. |

### Correction

Use both:

- Bonferroni: stringent.
- BH-FDR: preferred for H-MAGMA-style discovery, especially when the GWAS has
  fewer than 20 genome-wide significant loci.

Project threshold convention:

```text
Primary: BH-FDR < 0.05
Exploratory/low-hit GWAS: BH-FDR < 0.10
Strict check: Bonferroni < 0.05
```

## cMAGMA vs H-MAGMA

The key comparison is not simply whether H-MAGMA has more significant genes.
The real questions are:

1. Does H-MAGMA identify genes that cMAGMA misses?
2. Are those genes significant after correction?
3. Did H-MAGMA add biologically meaningful distal SNPs?
4. Are the added SNPs inside Duttke/Telese ATAC peaks?
5. Were TSS/promoter-overlap peaks excluded?
6. Is the signal stronger because of regulatory assignment, not just SNP count?

### H-MAGMA-Only Definition

For this project, compare by genomic coordinate when cMAGMA and H-MAGMA use
different gene IDs:

```text
KEY = CHR:START:STOP
```

An H-MAGMA-only gene is:

```text
significant in H-MAGMA under the chosen threshold
and not significant in cMAGMA at the same coordinate threshold
```

If cMAGMA has no row for the coordinate because no valid local SNPs were
available, record that explicitly.

## NSNP Interpretation

`NSNPS` is critical.

| Pattern | Interpretation |
|---|---|
| H NSNP = c NSNP and P is identical | H-MAGMA added nothing for this gene. |
| H NSNP > c NSNP and P improves | H-MAGMA-added SNPs may strengthen the gene signal. |
| H NSNP > c NSNP but P worsens | Added SNPs may add noise. |
| cMAGMA missing but H-MAGMA significant | Strong possible regulatory-contact rescue, but verify SNP assignments carefully. |

For a paper-level candidate, do not stop at `H NSNP > c NSNP`.
Check the actual added SNPs, their p-values, ATAC overlap, and loop assignment.

## Dataset-Specific Lessons So Far

### King 2025

Best current main dataset.

- Selected phenotypes match public report values.
- H-MAGMA-only genes exist under BH-FDR and Bonferroni criteria.
- Repeated candidates across selected phenotypes were extracted.
- `Lamtor1` is currently a strong biological lead because of the
  Ragulator/mTORC1 axis and reward/cue-learning relevance.

Key files:

```text
/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/King2025/reports/selected_1_3_traits_summary.md
/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/King2025/reports/selected_1_3_traits_hmagma_only_key_genes.tsv
/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/King2025/reports/king_1_3_repeated_genes_biology_annotation.tsv
```

### Kuhn 2025

Excellent GWAS reproduction exercise, weak H-MAGMA novelty so far.

- Key GWAS phenotypes reproduce nearly exactly.
- Gene-level cMAGMA/H-MAGMA does not currently provide a strong H-MAGMA-only
  novelty story for the tested heroin/nociception traits.
- This is still valuable for resume wording because it demonstrates GWAS
  reproduction and validation.

Key files:

```text
/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Kuhn2025/*_gwas_reproduction_check.csv
/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Kuhn2025/*_magma_hmagma_summary.csv
```

### Lara 2024

Good public summary-stat analysis, limited H-MAGMA novelty.

- Official `.mlma` files match report/article-level statistics.
- `auc` survives gene-level correction in both cMAGMA and H-MAGMA, but it is
  not H-MAGMA-specific.
- `exponential_k` and `dd_indiff_2` do not survive gene-level correction.

Key files:

```text
/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Lara2024_delaydiscounting/reports/SELECTED_MAGMA_CORRECTION_REPORT.md
/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Lara2024_delaydiscounting/reports/selected_correction_summary.tsv
```

### Cocaine 2026

Good dataset-hunt lesson, weak gene-level discovery.

- Four key addiction-relevant phenotypes were run through cMAGMA/H-MAGMA.
- No cMAGMA or H-MAGMA gene survived BH-FDR or Bonferroni correction.
- This likely reflects modest sample size and weak gene-level signal rather
  than a broken pipeline.

Key file:

```text
/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Cocaine2026/reports/cocaine2026_key_phenotypes_cmagma_hmagma_summary.tsv
```

## How To Write A Result Interpretation

Use this template for each phenotype:

```text
Phenotype:
Sample size:
GWAS source:
GWAS validation:
Top SNP/locus:
Top SNP p-value / -log10P:
Published match:

cMAGMA:
  tested genes:
  top gene:
  top P:
  BH-FDR significant genes:
  Bonferroni significant genes:

H-MAGMA:
  tested genes:
  top gene:
  top P:
  BH-FDR significant genes:
  Bonferroni significant genes:

H-MAGMA-only:
  number:
  strongest gene:
  NSNP cMAGMA -> H-MAGMA:
  evidence from added SNPs:
  ATAC/loop support:

Interpretation:
  Is this a main discovery, supporting result, negative result, or QC result?
```

## What Counts As A Resume-Worthy GWAS Analysis?

Minimum:

- Ran or validated a GWAS pipeline.
- Understood genotype/phenotype inputs.
- Produced or inspected SNP-level summary statistics.
- Matched output to a published report.

Better:

- Performed phenotype preprocessing/QC.
- Ran GCTA MLMA/MLMA-LOCO.
- Checked sample size, top SNP, effect size, p-value, and row counts.
- Converted summary stats for MAGMA.

Strong:

- Compared cMAGMA vs H-MAGMA.
- Applied multiple-testing correction.
- Identified H-MAGMA-only genes.
- Tested whether signals are explained by added SNP count.
- Built biological/locus-level interpretation.

Best:

- Integrated GWAS, Hi-C, ATAC-seq, gene annotation, and biology.
- Produced locus figures.
- Connected candidate genes/pathways to addiction-relevant phenotypes.

## Safe CV / Resume Bullet Options

Option 1:

```text
Reproduced and validated HS rat GWAS analyses using GCTA MLMA/MLMA-LOCO,
benchmarking local outputs against published SNP-level loci, sample sizes,
effect estimates, and report-level QTL statistics.
```

Option 2:

```text
Processed and interpreted SNP-level GWAS summary statistics from multiple HS
rat behavioral genetics datasets, including heroin vulnerability, cue-reactivity,
and delay discounting phenotypes.
```

Option 3:

```text
Integrated GWAS summary statistics with cMAGMA and frontal-cortex Hi-C-informed
H-MAGMA annotations to prioritize regulatory-contact candidate genes in
addiction- and behavior-relevant rat phenotypes.
```

Option 4:

```text
Performed reproducibility audits of published GWAS results by comparing local
or deposited MLMA outputs to reported top SNPs, effect sizes, sample sizes,
heritability estimates, and QTL counts.
```

Avoid:

```text
Discovered novel addiction genes from all GWAS datasets
```

because that overstates the current evidence. King supports a stronger
H-MAGMA-specific story; Kuhn, Lara, and cocaine are currently better described
as validation, specificity, or limited-result analyses.
