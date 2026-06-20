#!/bin/bash
set -euo pipefail

OUTDIR="/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/BMI"
LIFTOVER="/Volumes/external_1000GB_all/playground/enhancer/tools/liftOver"
CHAIN="/Users/pete/Desktop/playground/enhancer/data/rn6ToRn7.over.chain"
MAPPING_SOURCE="results.body_weight_g"
MAP_BED="$OUTDIR/rn6_to_rn7_liftover_map.bed"
MAP_UNMAPPED="$OUTDIR/rn6_to_rn7_liftover_unmapped.bed"
MAP_DUPLICATES="$OUTDIR/rn6_to_rn7_duplicate_ids.tsv"

phenotypes=(
  "results.bmi_w_tail"
  "results.bmi_wo_tail"
  "results.body_weight_g"
  "results.length_w_tail_cm"
  "results.length_wo_tail_cm"
)

echo "Creating one rn6 BED map from $MAPPING_SOURCE..."
awk -F',' 'BEGIN { OFS="\t" }
  NR > 1 && $1 != "" && $3 != "" && $2 != "" {
    chr = $1
    bp = $3
    snp = $2
    if (chr !~ /^chr/) {
      chr = "chr" chr
    }
    if (bp > 0) {
      print chr, bp - 1, bp, snp
    }
  }' "$OUTDIR/${MAPPING_SOURCE}.csv" > "$OUTDIR/${MAPPING_SOURCE}.rn6.bed"

echo "Running rn6 -> rn7 liftOver once..."
"$LIFTOVER" "$OUTDIR/${MAPPING_SOURCE}.rn6.bed" "$CHAIN" "$MAP_BED" "$MAP_UNMAPPED"

echo "Writing duplicate rn7 coordinate audit..."
awk -v dup="$MAP_DUPLICATES" 'BEGIN { OFS="\t" }
  {
    chr = $1
    sub(/^chr/, "", chr)
    rn7 = chr ":" $3
    count[rn7]++
    line[NR] = $0
    rn7_id[NR] = rn7
  }
  END {
    print "rn7_snp", "source_chr", "rn7_start", "rn7_end", "source_snp" > dup
    for (i = 1; i <= NR; i++) {
      if (count[rn7_id[i]] > 1) {
        split(line[i], a, "\t")
        chr = a[1]
        sub(/^chr/, "", chr)
        print rn7_id[i], chr, a[2], a[3], a[4] > dup
      }
    }
  }' "$MAP_BED"

for base in "${phenotypes[@]}"; do
  infile="$OUTDIR/${base}.csv"
  magma_tsv="$OUTDIR/${base}_rn7.magma.tsv"

  echo "Creating deduplicated rn7 MAGMA TSV for $base..."
  awk -v map="$MAP_BED" 'BEGIN {
      FS = OFS = "\t"
      while ((getline < map) > 0) {
        chr = $1
        sub(/^chr/, "", chr)
        rn7 = chr ":" $3
        source = $4
        rn7_count[rn7]++
        source_to_rn7[source] = rn7
      }
      close(map)
      FS = ","
      OFS = "\t"
      print "SNP", "p"
    }
    NR > 1 {
      source = $2
      pval = $9
      rn7 = source_to_rn7[source]
      if (rn7 != "" && rn7_count[rn7] == 1 && pval != "") {
        print rn7, pval
      }
    }' "$infile" > "$magma_tsv"

  wc -l "$magma_tsv"
done

wc -l "$OUTDIR/${MAPPING_SOURCE}.rn6.bed" "$MAP_BED" "$MAP_UNMAPPED" "$MAP_DUPLICATES"
