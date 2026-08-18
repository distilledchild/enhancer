# Resubmission analysis layout

The reviewer-driven analyses are organized into four reproducible modules:

1. `revision_main/`: coordinate normalization and the revised pooled loop resource.
2. `ATAC_validation/`: matched-null validation of TSS-excluded ATAC overlap.
3. `downsampling/`: full-depth, 140M-contact, and 250M-contact sensitivity analyses.
4. `hrdp_genotype_diversity/`: SNP-based genetic diversity analysis of the ten strains.

Run `revision_main` first because its coordinate cache and output tables are
inputs to the ATAC and downsampling modules. The modules share
`r_files/funcs.R`, which is synchronized to the Google Drive `r_files` root.

Each module writes rebuildable files to its own `results/` directory. Large or
source-controlled inputs are stored under that module's `inputs/` directory in
the Google Drive copy and are not deleted by the post-commit synchronization.

## Software requirements

Use a recent R installation with `tidyverse`, `data.table`, `fs`, `digest`,
`scales`, `ggrepel`, `patchwork`, `GenomicRanges`, `GenomeInfoDb`, `IRanges`,
`S4Vectors`, `rtracklayer`, `AnnotationDbi`, `clusterProfiler`, and
`org.Rn.eg.db`. The genotype module additionally uses PLINK 2.0 and downloads
the official Apple Silicon binary when needed; on other platforms, set
`PLINK2_BIN` to a local PLINK 2.0 executable.

The Google Drive copy preserves the same `r_files/revision/` structure and
includes the shared `r_files/funcs.R`. Collaborators should download the four
module directories and `funcs.R` together so the scripts can resolve all
relative paths without using the original author's home directory.
