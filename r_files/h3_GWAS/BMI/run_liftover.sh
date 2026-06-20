#!/bin/bash
INFILE="/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/BMI/results.body_weight_g.csv"
BED_RN6="/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/BMI/results.body_weight_g.rn6.bed"
BED_RN7="/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/BMI/results.body_weight_g.rn7.bed"
UNMAPPED="/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/BMI/results.body_weight_g.unmapped.bed"
MAGMA_TSV="/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/BMI/results.body_weight_g_rn7.magma.tsv"

LIFTOVER="/Volumes/external_1000GB_all/playground/enhancer/tools/liftOver"
CHAIN="/Users/pete/Desktop/playground/enhancer/data/rn6ToRn7.over.chain"

echo "Creating BED with P-values..."
awk -F',' 'NR>1 {
  chr=$1; bp=$3; snp=$2; pval=$9;
  if (chr !~ /^chr/) { chr = "chr" chr }
  print chr "\t" (bp-1) "\t" bp "\t" snp "\t" pval
}' "$INFILE" > "$BED_RN6"

echo "Running liftOver..."
"$LIFTOVER" "$BED_RN6" "$CHAIN" "$BED_RN7" "$UNMAPPED"

echo "Creating new MAGMA TSV..."
echo -e "SNP\tp" > "$MAGMA_TSV"
awk '{
  chr=$1; sub(/^chr/, "", chr);
  pos=$3; pval=$5;
  print chr ":" pos "\t" pval
}' "$BED_RN7" >> "$MAGMA_TSV"

echo "Done! Final MAGMA TSV line count:"
wc -l "$MAGMA_TSV"
