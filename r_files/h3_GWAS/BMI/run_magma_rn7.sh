#!/bin/bash
MAGMA='/Volumes/external_1000GB_all/playground/enhancer/tools/magma'
BFILE='/Volumes/external_1000GB_all/playground/enhancer/data/hs_data/v4/HS_genotypes_v4'
ANNOT_HMAGMA='/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/nicotine/hs.hic.duttke_telese_atac.non_tss_promoter_filtered.annot'
ANNOT_CMAGMA='/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/nicotine/hs.exonpro.ONLY.annot'
OUTDIR='/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/BMI'
PVAL_FILE="$OUTDIR/results.body_weight_g_rn7.magma.tsv"
BASE="results.body_weight_g_rn7"

echo "Running H-MAGMA..."
"$MAGMA" --bfile "$BFILE" \
  --pval "$PVAL_FILE" use=SNP,p N=3173 \
  --gene-annot "$ANNOT_HMAGMA" \
  --out "$OUTDIR/${BASE}_hmagma"

echo "Running cMAGMA..."
"$MAGMA" --bfile "$BFILE" \
  --pval "$PVAL_FILE" use=SNP,p N=3173 \
  --gene-annot "$ANNOT_CMAGMA" \
  --out "$OUTDIR/${BASE}_cmagma"
