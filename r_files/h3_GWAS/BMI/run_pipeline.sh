#!/bin/bash
if [ -z "$1" ]; then
  echo "Usage: $0 <phenotype_base>"
  echo "Example: $0 results.bmi_wo_tail"
  exit 1
fi

BASE=$1
DIR="/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/BMI"
INFILE="$DIR/${BASE}.csv"
OUTBED="$DIR/${BASE}.rn6.bed"
BED_RN7="$DIR/${BASE}.rn7.bed"
UNMAPPED="$DIR/${BASE}.unmapped.bed"
MAGMA_TSV="$DIR/${BASE}_rn7.magma.tsv"

LIFTOVER="/Volumes/external_1000GB_all/playground/enhancer/tools/liftOver"
CHAIN="/Users/pete/Desktop/playground/enhancer/data/rn6ToRn7.over.chain"

echo "=== Processing ${BASE} ==="

echo "[1/3] Extracting BED format from CSV..."
awk -F',' 'NR>1 {
  chr=$1; bp=$3; snp=$2; pval=$9;
  if (chr !~ /^chr/) { chr = "chr" chr }
  print chr "\t" (bp-1) "\t" bp "\t" snp "\t" pval
}' "$INFILE" > "$OUTBED"

echo "[2/3] Running liftOver..."
"$LIFTOVER" "$OUTBED" "$CHAIN" "$BED_RN7" "$UNMAPPED"

echo "[3/3] Creating new MAGMA TSV..."
echo -e "SNP\tp" > "$MAGMA_TSV"
awk '{
  chr=$1; sub(/^chr/, "", chr);
  pos=$3; pval=$5;
  print chr ":" pos "\t" pval
}' "$BED_RN7" >> "$MAGMA_TSV"

echo "Done! Final MAGMA TSV line count:"
wc -l "$MAGMA_TSV"

MAGMA='/Volumes/external_1000GB_all/playground/enhancer/tools/magma'
BFILE='/Volumes/external_1000GB_all/playground/enhancer/data/hs_data/v4/HS_genotypes_v4'
ANNOT_HMAGMA="/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/nicotine/hs.hic.duttke_telese_atac.non_tss_promoter_filtered.annot"
ANNOT_CMAGMA="/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/nicotine/hs.exonpro.ONLY.annot"

echo "=== Running H-MAGMA ==="
"$MAGMA" --bfile "$BFILE" \
  --pval "$MAGMA_TSV" use=SNP,p N=3173 \
  --gene-annot "$ANNOT_HMAGMA" \
  --out "$DIR/${BASE}_rn7_hmagma"

echo "=== Running cMAGMA ==="
"$MAGMA" --bfile "$BFILE" \
  --pval "$MAGMA_TSV" use=SNP,p N=3173 \
  --gene-annot "$ANNOT_CMAGMA" \
  --out "$DIR/${BASE}_rn7_cmagma"
