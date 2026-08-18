# Downsampling validation bundle

This directory documents the files required to rerun the R-based 140M and
matched full-depth/140M/250M depth-sensitivity analyses. Raw Juicer
`merged_nodups.txt` contact files are not included because they are retained on
ISAAC and are required only to regenerate the downsampled contact sets.

## Run from this directory

- `downsampling_140M.R`: standalone all-ten-library full-depth/140M
  validation.
- `downsampling_140M_250M_comparison.R`: matched seven-library
  full-depth/140M/250M comparison.
- `downsampling_depth_functions.R`: shared validation functions.
- `../../../funcs.R`: shared HiCCUPS parsing, pooling, and approximate-locus
  functions.
- `build_downsampling_input_manifests.R`: input
  inventory and provenance checks.

The `r_files/revision/downsampling/140M/` and `250M/` directories retain the
HPC downsampling, `.hic` creation, and Colab HiCCUPS execution scripts.

## Required R-validation inputs

### Full-depth calls

For each of the ten libraries:

`inputs/hic/2023A/hic30_w_sb_options/<sample>/hiccups_5k10k25k/merged_loops.bedpe`

### Downsampled calls and provenance

The following tree is copied without changing sample or result-directory
names:

`inputs/juicer_downsample_q30_140M_250M/`

- `140M/`: all ten libraries.
- `250M/`: the seven libraries with at least 250M source contacts.
- Each eligible sample contains `downsampling_qc.tsv`,
  `hic_creation_qc.tsv`, the downsampled `.hic`, and a
  `*.hiccups.5k10k25k/` directory containing `merged_loops.bedpe` and the
  5/10/25-kb `postprocessed_pixels_*.bedpe` files.

### Coordinate-normalized annotations

`../revision_main/cache_data/`

- `df.ensembl.transcript.coordinate.normalized.rds`
- `df.promoter.annotation.coordinate.normalized.rds`
- `gr.atac.rds`

## Default outputs

- `results/full_140M/`
- `results/full_140M_250M/`

These output directories are rebuildable and are synchronized by the Git
post-commit hook. Large source files under `inputs/` are staged separately and
are protected from deletion by the hook.

## Post-commit synchronization

The repository post-commit hook synchronizes this code and its rebuildable
results to the matching Google Drive directory. The large files under
`inputs/` are staged once and protected from deletion by the hook.

## Scope

The copied bundle reproduces the R validation from completed HiCCUPS calls. To
recreate the contact downsampling itself, use the retained HPC scripts together
with the original ISAAC `merged_nodups.txt`, `inter_30.txt`, restriction-site,
chromosome-size, and Juicer Tools inputs.
