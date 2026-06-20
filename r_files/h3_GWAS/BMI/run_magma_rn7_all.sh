#!/bin/bash
set -euo pipefail

MAGMA="/Volumes/external_1000GB_all/playground/enhancer/tools/magma"
BFILE="/Volumes/external_1000GB_all/playground/enhancer/data/hs_data/v4/HS_genotypes_v4"
ANNOT_HMAGMA="/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/nicotine/hs.hic.duttke_telese_atac.non_tss_promoter_filtered.annot"
ANNOT_CMAGMA="/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/nicotine/hs.exonpro.ONLY.annot"
OUTDIR="/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/BMI"
N=3173

phenotypes=(
  "results.bmi_w_tail"
  "results.bmi_wo_tail"
  "results.body_weight_g"
  "results.length_w_tail_cm"
  "results.length_wo_tail_cm"
)

for base in "${phenotypes[@]}"; do
  pval_file="$OUTDIR/${base}_rn7.magma.tsv"
  out_base="$OUTDIR/${base}_rn7"

  echo "== $base H-MAGMA =="
  "$MAGMA" \
    --bfile "$BFILE" \
    --pval "$pval_file" use=SNP,p N="$N" \
    --gene-annot "$ANNOT_HMAGMA" \
    --out "${out_base}_hmagma"

  echo "== $base cMAGMA =="
  "$MAGMA" \
    --bfile "$BFILE" \
    --pval "$pval_file" use=SNP,p N="$N" \
    --gene-annot "$ANNOT_CMAGMA" \
    --out "${out_base}_cmagma"
done
