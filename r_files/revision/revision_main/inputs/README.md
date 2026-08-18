# Revision-main inputs

This directory mirrors the source-data layout expected by the coordinate
preparation script:

- `data/`: GTF, EPD, ATAC, CTCF-motif, chromosome-size, and QC inputs.
- `hic/2023A/hic30_w_sb_options/`: full-depth merged HiCCUPS loop calls.

The post-commit hook preserves this directory and does not recopy its large
files after every commit.
