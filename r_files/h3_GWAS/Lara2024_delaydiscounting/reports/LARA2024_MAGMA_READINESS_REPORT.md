# Lara 2024 Delay Discounting MAGMA Readiness Report

Date: 2026-06-18

## Source

- GitHub phenotype and processing repo: https://github.com/Palmer-Lab-UCSD/HSrat_delaydiscounting
- UCSD Digital Collections object: https://library.ucsd.edu/dc/object/bb7761581p
- DOI listed by the GitHub repo for GWAS data/results/report: https://doi.org/10.6075/J0BG2PDS
- Publication: Lara et al. 2024, Genes Brain Behav. doi: 10.1111/gbb.12909

## Local Data

Work root:

```text
/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Lara2024_delaydiscounting
```

Downloaded UCSD files:

```text
data/ucsd_1_1.csv      data dictionary
data/ucsd_2_1.html     GWAS report
data/ucsd_3_1.csv      delay discounting traits
data/ucsd_4_1.zip      GWAS results archive
```

Extracted GWAS summary stats:

```text
data/gwas_results/gwas_dd_results/*.loco.mlma
```

There are 12 autosomal/LOCO GCTA MLMA summary-stat files. Each has 5,326,108 SNP rows and columns:

```text
Chr SNP bp A1 A2 Freq b se p
```

## MAGMA Input Conversion

Converted all 12 `.loco.mlma` files to MAGMA p-value inputs:

```text
magma_inputs/*.magma.tsv
```

Each output has:

```text
SNP p
```

and 5,326,108 SNP p-values.

Conversion script used:

```text
/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022/scripts/mlma_to_magma.py
```

Phenotype sample size from `data/ucsd_3_1.csv`:

```text
N = 629 for all regressedlr_* traits
sex counts: F = 319, M = 310
```

## Compatibility Checks

MAGMA bfile used:

```text
/Volumes/external_1000GB_all/playground/enhancer/data/hs_data/v4/HS_genotypes_v4
```

MAGMA binary:

```text
/Volumes/external_1000GB_all/playground/enhancer/tools/magma
MAGMA version v1.08 mac
```

For `regressedlr_auc.magma.tsv`:

```text
p-value SNPs: 5,326,108
SNPs found in v4 BIM: 5,288,132
match rate: 99.287% of Lara p-value SNPs
```

cMAGMA annotation:

```text
/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/nicotine/hs.exonpro.ONLY.annot
annotation genes: 21,235
annotation unique SNPs: 289,346
annotation SNPs overlapping Lara AUC p-values: 199,733
annotation overlap rate: 69.029%
```

H-MAGMA strict annotation:

```text
/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/nicotine/hs.hic.duttke_telese_atac.non_tss_promoter_filtered.annot
annotation genes: 21,383
annotation unique SNPs: 320,281
annotation SNPs overlapping Lara AUC p-values: 222,804
annotation overlap rate: 69.565%
```

## Smoke Test Results

Trait tested:

```text
regressedlr_auc
N = 629
```

cMAGMA smoke run:

```text
out prefix: results/regressedlr_auc_cmagma
valid SNP p-values in bfile: 5,288,132
gene definitions read: 21,235
genes with valid SNPs in genotype data: 18,643
elapsed: 00:09:20
result files:
  results/regressedlr_auc_cmagma.genes.out
  results/regressedlr_auc_cmagma.genes.raw
```

H-MAGMA strict smoke run:

```text
out prefix: results/regressedlr_auc_hmagma_strict_duttke_telese
valid SNP p-values in bfile: 5,288,132
gene definitions read: 21,383
genes with valid SNPs in genotype data: 18,963
elapsed: 00:11:11
result files:
  results/regressedlr_auc_hmagma_strict_duttke_telese.genes.out
  results/regressedlr_auc_hmagma_strict_duttke_telese.genes.raw
```

Top AUC cMAGMA gene-level hit:

```text
GENE                  CHR  START     STOP      NSNPS  N  ZSTAT   P
20:32238545:32239613  20   32238545  32239613  4      629 4.5729  2.4055e-06
```

Top AUC H-MAGMA strict gene-level hit:

```text
GENE              CHR  START     STOP      NSNPS  N  ZSTAT   P
ENSRNOG00000066051 20  32238545  32239613  4      629 4.5729  2.4055e-06
```

## Conclusion

The Lara 2024 public GWAS summary stats are directly usable for local MAGMA-style cMAGMA and H-MAGMA runs after a simple `.mlma` to `SNP p` conversion. The AUC smoke test completed successfully for both cMAGMA and H-MAGMA strict annotations.

The only caveat is annotation coverage: the Lara p-value SNPs match the v4 BIM very well, but the existing cMAGMA/H-MAGMA annotation SNP sets overlap about 69% with the Lara summary-stat SNPs. This is sufficient to run, but not a perfect annotation coverage match.

## Batch Template

Sequential full batch template:

```bash
BASE=/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Lara2024_delaydiscounting
MAGMA=/Volumes/external_1000GB_all/playground/enhancer/tools/magma
BFILE=/Volumes/external_1000GB_all/playground/enhancer/data/hs_data/v4/HS_genotypes_v4
CM_ANN=/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/nicotine/hs.exonpro.ONLY.annot
HM_ANN=/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/nicotine/hs.hic.duttke_telese_atac.non_tss_promoter_filtered.annot

mkdir -p "$BASE/results" "$BASE/logs"

for pval in "$BASE"/magma_inputs/*.magma.tsv; do
  trait=$(basename "$pval" .magma.tsv)

  "$MAGMA" \
    --bfile "$BFILE" \
    --pval "$pval" use=SNP,p N=629 \
    --gene-annot "$CM_ANN" \
    --out "$BASE/results/${trait}_cmagma" \
    > "$BASE/logs/${trait}_cmagma.stdout.log" 2>&1

  "$MAGMA" \
    --bfile "$BFILE" \
    --pval "$pval" use=SNP,p N=629 \
    --gene-annot "$HM_ANN" \
    --out "$BASE/results/${trait}_hmagma_strict_duttke_telese" \
    > "$BASE/logs/${trait}_hmagma_strict_duttke_telese.stdout.log" 2>&1
done
```
