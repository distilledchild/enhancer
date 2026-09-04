# TSS +/-2-kb promoter-window sensitivity analysis (4,001-bp window; 50-bp ATAC overlap)

## Purpose

The production analysis defines promoter-associated anchors by direct overlap
with Ensembl or EPD TSS +/-1-kb windows. This sensitivity analysis uses TSS
+/-2-kb windows, matching the window described in Cao et al. (2026).

Both linked definitions were changed together:

- promoter-anchor assignment: TSS +/-2 kb (4,001 bp, 1-based closed);
- exclusion from distal ATAC evidence: TSS +/-2 kb;
- minimum non-TSS ATAC overlap remained 50 bp;
- the denominator remained 31,021 pooled HiCCUPS call records shorter than 2 Mb.

The production +/-1-kb outputs were not overwritten.

## Main results

| Metric | +/-1 kb | +/-2 kb | Difference |
|---|---:|---:|---:|
| No direct promoter/TSS | 14,678 (47.3%) | 13,579 (43.8%) | -1,099 (-3.54 pp) |
| Exactly one promoter anchor | 12,295 (39.6%) | 12,859 (41.5%) | +564 (+1.82 pp) |
| Promoter at both anchors | 4,048 (13.0%) | 4,583 (14.8%) | +535 (+1.72 pp) |
| Any promoter-associated loop | 16,343 (52.7%) | 17,442 (56.2%) | +1,099 (+3.54 pp) |
| Putative regulatory | 13,376 (43.1%) | 14,210 (45.8%) | +834 (+2.69 pp) |
| Promoter-associated without distal ATAC | 2,967 (9.6%) | 3,232 (10.4%) | +265 (+0.85 pp) |

Putative-regulatory membership transitions:

- retained putative: 13,147;
- gained putative: 1,063;
- lost putative: 229;
- neither definition: 16,582.

Thus, 98.3% of the primary putative set was retained
(13,147 / 13,376). The Jaccard similarity between the two putative sets was
91.1% (13,147 / 14,439). Major-category assignments were unchanged for 95.2%
of all loops (29,528 / 31,021).

Putative-regulatory percentages increased at every HiCCUPS resolution:

| Resolution | +/-1 kb | +/-2 kb | Difference |
|---|---:|---:|---:|
| 5 kb | 29.1% | 32.3% | +3.2 pp |
| 10 kb | 39.4% | 42.2% | +2.8 pp |
| 25 kb | 54.0% | 56.3% | +2.3 pp |

## Gene-level loop-count sensitivity

Gene-level counts include only direction-supported assignments: genes at the
promoter anchor of single-promoter putative loops, and genes on the supported
promoter side or sides of dual-promoter putative loops. Each distinct pooled
`loop_id` is counted once per Ensembl gene. Calls at different HiCCUPS
resolutions remain distinct pooled loop records.

The 5-, 10-, and 25-kb columns refer to HiCCUPS loop-calling resolution, not
the ATAC-overlap threshold. The minimum ATAC overlap is 50 bp in this analysis.

| Metric | +/-1 kb | +/-2 kb | Difference |
|---|---:|---:|---:|
| Direction-supported gene-loop pairs | 19,292 | 20,878 | +1,586 |
| Unique genes with at least one putative loop | 8,775 | 9,114 | +339 |
| Genes with at least 2 putative loops | 4,631 | 4,947 | +316 |
| Genes with at least 5 putative loops | 830 | 988 | +158 |
| Genes with at least 10 putative loops | 41 | 56 | +15 |

Among the union of genes found under either definition, 7,397 retained the
same loop count, 830 increased, 321 decreased, 566 appeared only under the
+/-2-kb definition, and 227 appeared only under the +/-1-kb definition. The
8,548 shared genes had a Spearman loop-count correlation of 0.922 and a
Pearson correlation of 0.929. Top-gene overlap was 80% for the top 25, 76% for
the top 50, and 75% for the top 100.

Consolidated gene-list output:

- `outputs/promoter-window-anchor-sensitivity/promoter_window_anchor_sensitivity_atac50bp_pooled_5k_10k_25k_gene_lists.xlsx`
  contains the 2,001-, 3,001-, and 4,001-bp promoter-window results for both
  full-anchor interval and one-base anchor-midpoint matching. Each of its six
  sheets contains `ensembl_gene_id`, `gene_symbol`, and `n`, where `n` is the
  number of distinct pooled 5-, 10-, and 25-kb loop records assigned to the
  gene. Superseded standalone promoter-window gene-list files are removed only
  after this workbook is created successfully.

Other gene-level sensitivity files retained here:

- `promoter_window_2001bp_atac_overlap_50bp_vs_promoter_window_4001bp_atac_overlap_50bp_gene_loop_count_comparison.tsv`: side-by-side
  +/-1-kb versus +/-2-kb counts, rank changes, and membership categories;
- `promoter_window_2001bp_atac_overlap_50bp_vs_promoter_window_4001bp_atac_overlap_50bp_gained_putative_gene_loop_counts.tsv`: gene
  counts contributed by the 1,063 gained putative loops;
- `promoter_window_2001bp_atac_overlap_50bp_vs_promoter_window_4001bp_atac_overlap_50bp_lost_putative_gene_loop_counts.tsv`: gene counts
  contributed by the 229 lost putative loops;
- `promoter_window_2001bp_atac_overlap_50bp_vs_promoter_window_4001bp_atac_overlap_50bp_gene_loop_count_summary.tsv` and
  `promoter_window_2001bp_atac_overlap_50bp_vs_promoter_window_4001bp_atac_overlap_50bp_gene_rank_stability.tsv`: compact sensitivity
  summaries.

All TSV exports also record the promoter-window width and minimum ATAC-overlap
threshold in leading columns. The one-ID-per-line text file carries both values
in its filename so it remains directly usable as a gene-list input.

## Interpretation

Doubling the promoter/TSS flank produces the expected moderate increase in
promoter-associated and putative-regulatory calls, but it does not change the
qualitative conclusion or the resolution trend. The primary +/-1-kb result is
therefore robust to the Cao et al. +/-2-kb definition. The +/-2-kb result can be
reported as a sensitivity analysis; changing the primary definition is not
required by these results alone.

The high shared-gene correlation and substantial top-rank overlap also show
that the broad gene-level loop-count pattern is preserved, although membership
and exact counts near the promoter-window boundary are definition-sensitive.

The 229 lost putative calls are expected consequences of applying the wider
definition consistently. Some previously distal ATAC signal becomes part of
the expanded TSS exclusion, and some single-promoter loops become dual-promoter
loops without a completed promoter-to-opposite-residual-ATAC direction.

## Reproduction

Run from the directory containing the `enhancer` repository:

```bash
env \
  PROMOTER_WINDOW_FLANK_BP=2000 \
  ATAC_MINIMUM_OVERLAP_BP=50 \
  PROMOTER_WINDOW_SENSITIVITY_ONLY=1 \
  RESUBMIT_OUTPUT_DIR=/Users/pete/Desktop/playground/enhancer/r_files/revision/revision_main/sensitivity_promoter_window_2kb \
  Rscript /Users/pete/Desktop/playground/enhancer/r_files/revision/revision_main/01_promoter_enhancer_interaction_resubmit_loop_annotation.R
```

The default remains 1,000 bp when `PROMOTER_WINDOW_FLANK_BP` is unset. A
default-window regression run reproduced every locked production count exactly.
