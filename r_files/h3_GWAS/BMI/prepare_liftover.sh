#!/bin/bash
INFILE="/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/BMI/results.body_weight_g.csv"
OUTBED="/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/BMI/results.body_weight_g.rn6.bed"

echo "Extracting BED format from CSV..."
# Chr,SNP,bp,A1,A2,Freq,b,se,p
# BED: chr start end snp_id
awk -F',' 'NR>1 {
  chr=$1;
  bp=$3;
  snp=$2;
  # Make sure chr has "chr" prefix
  if (chr !~ /^chr/) { chr = "chr" chr }
  start = bp - 1;
  end = bp;
  print chr "\t" start "\t" end "\t" snp
}' "$INFILE" > "$OUTBED"

echo "BED file created. Line count:"
wc -l "$OUTBED"
