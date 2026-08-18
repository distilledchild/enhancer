# ATAC matched-null validation

Run `atac_validation.R` after completing `../revision_main/`.

The analysis reads the revised pooled resource from
`../revision_main/results/`, coordinate-normalized annotations from
`../revision_main/cache_data/`, and chromosome sizes from
`../revision_main/inputs/data/tracks/rn7.chrom.sizes`. It writes permutation
tables and figures to `results/`.

The matched-null analysis tests whether candidate distal anchors overlap
TSS-excluded ATAC peaks more often than anchors matched by chromosome,
HiCCUPS resolution, anchor side and width, and loop-distance bin. ATAC overlap
is interpreted as open-chromatin support, not functional enhancer validation.
