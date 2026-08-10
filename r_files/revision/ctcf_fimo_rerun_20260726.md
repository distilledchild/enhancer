# CTCF FIMO rerun (2026-07-26)

## Purpose

Rebuild the CTCF sequence-motif annotation independently of promoter-enhancer
loop selection. CTCF is treated only as predicted structural sequence evidence,
not as CTCF binding or direct evidence of enhancer function.

## Reference and software

- Genome FASTA: `/Volumes/external_1000GB_all/playground/enhancer/data/rn7chr.fa`
- Assembly: mRatBN7.2/rn7, 22 sequences, 2,633,473,415 bases
- MEME Suite environment: `/Volumes/external_1000GB_all/playground/enhancer/tools/conda-envs/meme-suite`
- FIMO/Tomtom version: 5.5.9
- Primary motif: JASPAR 2022 CTCF `MA0139.1`, width 19 bp
- Background: order-0 frequencies estimated from the same rn7 FASTA

## Motif redundancy audit

The historical 100-motif collection was compared against itself with Tomtom
(Pearson distance, minimum overlap 8). Of 9,900 possible non-self directed
comparisons, 9,253 were significant at q <= 0.05. `MA0139.1` was significant
against 96 of the other 99 motifs. Therefore, counts across the historical
collection cannot be interpreted as counts of independent CTCF sites. The
historical collection is retained only for legacy/sensitivity comparisons.

## Primary FIMO scan

FIMO scanned the full rn7 genome with `MA0139.1`, the FASTA-derived background,
a reporting threshold of p <= 1e-4, and `--max-stored-scores 2000000`.

- Reported windows at p <= 1e-4: 1,191,942
- Windows at q <= 0.05: 20,653
- Windows at q <= 0.01: 1,071
- Overlap-collapsed loci at q <= 0.05: 20,641
- Overlap-collapsed loci at q <= 0.01: 1,071

The old combined FIMO file contained 75,870 `MA0139.1` rows because the default
stored-score limit dynamically truncated its nominal p <= 1e-4 output. Its
q <= 0.05 and q <= 0.01 counts were exactly 20,653 and 1,071, respectively,
matching the new scan and validating the rerun.

## Loop-level structural annotation

The current coordinate-normalized pooled resource contains 31,021 unique loop
records (not the 31,019 records from the older 58,992-loop copied input). CTCF
motif loci were annotated at each anchor but were not used to retain or discard
loops.

At q <= 0.05:

- Any predicted motif at both anchors: 5,245/31,021 (16.9%)
- Predicted convergent orientation (`+` at left, `-` at right): 3,894/31,021 (12.6%)
- Revised putative-regulatory loops with motifs at both anchors: 1,870/10,469 (17.9%)
- Revised putative-regulatory loops with predicted convergent orientation: 1,382/10,469 (13.2%)

At q <= 0.01:

- Any predicted motif at both anchors: 61/31,021 (0.2%)
- Predicted convergent orientation: 40/31,021 (0.1%)
- Revised putative-regulatory loops with motifs at both anchors: 25/10,469 (0.2%)

The q <= 0.01 tier is too sparse to serve as the primary structural annotation.
The q <= 0.05 tier is the practical primary annotation, with q <= 0.01 reported
as a strict sensitivity analysis. Neither tier is a P-E interaction filter.

## Reproducible analysis files

- `ctcf_fimo_postprocess.R`: filters FIMO results and collapses overlapping shifted windows
- `ctcf_fimo_loop_annotation.R`: annotates loop anchors and predicted motif orientation
- Full FIMO directory: `/Volumes/external_1000GB_all/playground/enhancer/data/ctcf/fimo_rerun_20260726`

