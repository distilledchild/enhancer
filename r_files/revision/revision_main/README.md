# Revised pooled loop analysis

Run these files in order:

1. `promoter_enhancer_interaction_resubmit_coord_prep.R`
2. `promoter_enhancer_interaction_resubmit.R`

The first script audits BED/BEDPE and GTF coordinate conventions, creates
strand-aware TSS and promoter annotations, and writes reusable R objects to
`cache_data/`. The second script builds the revised pooled HiCCUPS call
resource and evidence-layered annotations in `results/`.

## Inputs

`inputs/data/` contains the rn7 GTF, EPD files, ATAC peaks, predicted CTCF
motifs, chromosome sizes, and library-complexity table. Full-depth HiCCUPS
calls are stored under:

`inputs/hic/2023A/hic30_w_sb_options/<sample>/hiccups_5k10k25k/merged_loops.bedpe`

Both scripts locate the shared `r_files/funcs.R` by walking up from their own
directory, so the same Google Drive tree can be run on another computer.
