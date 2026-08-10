# Sequencing-Depth Sensitivity Analysis Summary

## Reviewer concern

The reviewer identified sequencing depth as a major technical driver of the
between-library differences in Hi-C loop yield and requested normalization to
a common number of contacts. The reviewer also questioned whether pairwise
loop sharing could support strain-level biological conclusions without
replication.

## Purpose

The downsampling analyses test whether loop-call yield and downstream
annotations remain stable after equalizing usable contact depth. They do not
create biological replication, establish a universal minimum sequencing
depth, or validate strain-specific loop differences.

Usable contacts were sampled from duplicate-removed Juicer
`merged_nodups.txt` files. Both reads were required to have MAPQ >=30, and
intra-fragment pairs were excluded. Downsampled `.hic` files were rebuilt and
HiCCUPS was rerun at 5, 10, and 25 kb with the same settings used across the
depth conditions.

## Analysis designs

### All-ten-library 140M analysis

All ten libraries were retained by downsampling each to exactly 140 million
valid contacts. This analysis directly evaluates normalization across the
complete sample set, but the low target causes substantial loss of individual
calls, particularly at 5 kb.

### Matched seven-library depth comparison

Only seven libraries had at least 250 million valid contacts. Full-depth,
140M, and 250M results were therefore compared within the same seven
libraries: `592BB`, `607`, `D765A`, `DA21A`, `DA68A`, `DBA9A`, and `DE8BA`.
The other three libraries were excluded from this matched comparison only
because their original contact counts were below 250M.

The all-ten and matched-seven estimates answer different questions and must
not be combined as if they came from one comparison.

## Results

### Standalone all-ten 140M analysis

| Metric | Full depth | 140M |
|---|---:|---:|
| Mean loop calls per library | 5,753.7 | 3,062.2 |
| Loop-count CV | 0.335 | 0.163 |
| Original source depth vs loop count, Spearman rho | 0.794 | -0.115 |
| Pooled exact call records | 31,021 | 14,414 |
| Full-depth exact calls recovered | - | 35.25% |

These values are valid for the ten-library analysis. The post-downsampling
correlation uses original source depth as the predictor and asks whether a
residual association remains after equalization; it does not prove that every
technical depth effect was removed.

### Matched seven-library comparison

| Condition | Mean loop calls | Loop-count CV | Original source depth vs loop count, rho | P value |
|---|---:|---:|---:|---:|
| Full depth | 6,740.6 | 0.1813 | 0.3929 | 0.383 |
| 140M | 3,146.0 | 0.1859 | -0.4643 | 0.294 |
| 250M | 5,358.1 | 0.1454 | -0.4286 | 0.337 |

Within the matched seven libraries, 140M did not reduce loop-count CV relative
to full depth. The 250M CV was descriptively lower, but the sample size is
seven and none of the depth correlations was statistically significant.

| Query condition | Exact full-depth recovery | Approximate-locus full-depth recovery | Exact pooled Jaccard |
|---|---:|---:|---:|
| 140M | 28.26% | 48.28% | 0.244 |
| 250M | 59.27% | 74.26% | 0.501 |

The higher 250M target retained substantially more full-depth information.
Approximate loci are a positional sensitivity analysis constructed jointly
across conditions; exact resolution-specific calls remain the primary unit.

Regulatory-category composition was also closer to full depth at 250M than at
140M. Total category-composition variation was 2.18 percentage points for
full-depth versus 250M and 7.00 percentage points for full-depth versus 140M.
Gene presence was more stable at 250M, but gene-loop count stability remained
moderate (full-depth versus 250M Spearman rho = 0.609 for exact-call counts).

## Defensible interpretation

1. Sequencing depth materially affects loop yield and individual call recovery.
2. The 140M all-ten analysis is useful for complete-cohort normalization but is
   an aggressive downsampling level with strong information loss.
3. The matched seven-library comparison shows that 250M preserves substantially
   more full-depth calls and annotation composition than 140M.
4. Neither analysis creates biological replication or supports claims of
   conserved topology, consensus loops, strain-specific divergence, or genetic
   effects on looping.
5. Recurrent multi-library calls represent technical recurrence, not proof of
   biological conservation.
6. Gene ranking and GO enrichment remain exploratory because their stability is
   moderate and each strain is represented by one library.
7. Results must be reported by resolution because 5 kb recovery is especially
   depth sensitive.

## Files to audit

- `r_files/downsampling_140M.R`: standalone ten-library 140M analysis.
- `r_files/downsampling_140M_250M_comparison.R`: matched seven-library
  full-depth/140M/250M comparison.
- `r_files/downsampling_depth_functions.R`: shared comparison, recovery,
  provenance, category, and gene-stability functions.
- `r_files/revision/downsampling/build_downsampling_input_manifests.R`: input
  inventory and provenance manifest generation.
- `r_files/revision/downsampling/140M/`: 140M downsampling, `.hic` creation, and
  HiCCUPS execution scripts.
- `r_files/revision/downsampling/250M/`: 250M downsampling, `.hic` creation, and
  HiCCUPS execution scripts.
- `r_files/downsampling_140M_250M_outputs/`: rebuildable matched-comparison
  output tables.
