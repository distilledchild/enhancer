# 140M-Contact Downsampling Analysis Summary

## Completion Status

- Completed libraries: 10/10 (`592BB`, `607`, `74AA`, `A2DB`, `D765A`, `DA08A`, `DA21A`, `DA68A`, `DBA9A`, and `DE8BA`)
- Common depth: exactly 140,000,000 MAPQ >=30 valid Hi-C contacts per library
- Rationale for target: immediately below the lowest original usable-contact count (141,993,330), allowing all ten libraries to be retained
- Rebuilt maps: 5, 10, and 25 kb `.hic` resolutions
- Loop calling: HiCCUPS at 5, 10, and 25 kb, followed by the same `<2 Mb` restriction
- Analysis run: 2026-07-30 19:26:28 CDT

## Main Results

| Metric | Full depth | Downsampled to 140M | Interpretation |
|---|---:|---:|---|
| Mean total loop calls per library | 5,753.7 | 3,062.2 | Fewer calls at the common lower depth |
| Across-library loop-count CV | 0.335 | 0.163 | 51.5% reduction in between-library dispersion |
| Across-library loop-count range | 6,017 | 1,360 | 77.4% reduction |
| Source depth vs total loop calls, Spearman rho | 0.794 | -0.115 | The original positive depth association was removed |
| Pooled exact call records | 31,021 | 14,414 | Full-depth union contains many depth-sensitive calls |
| Full-depth pooled calls exactly recovered | - | 10,936 (35.25%) | Individual exact-call recovery is limited |
| Pooled exact Jaccard | - | 0.317 | Partial overlap between full and downsampled unions |

## Resolution-Specific Pooled Recovery

| Resolution | Full-depth pooled calls | Downsampled pooled calls | Exact shared calls | Full-depth calls recovered |
|---|---:|---:|---:|---:|
| 5 kb | 6,522 | 879 | 777 | 11.91% |
| 10 kb | 11,977 | 5,022 | 3,867 | 32.29% |
| 25 kb | 12,522 | 8,513 | 6,292 | 50.25% |
| All | 31,021 | 14,414 | 10,936 | 35.25% |

The median sample-level exact recovery was 22.01% overall. A secondary HiCCUPS-radius matching analysis increased median recovery to 26.25%; it is a sensitivity analysis and does not redefine the exact resource unit.

## Regulatory-Annotation Stability

Broad category proportions were more stable than exact call identity.

| Category | Full-depth proportion | Downsampled proportion | Change |
|---|---:|---:|---:|
| No direct promoter/TSS | 47.32% | 42.66% | -4.66 percentage points |
| Promoter-promoter-compatible | 13.05% | 15.05% | +2.01 percentage points |
| Putative regulatory with opposite-anchor non-TSS ATAC | 33.75% | 37.26% | +3.51 percentage points |
| Single promoter without opposite-anchor non-TSS ATAC | 5.89% | 5.03% | -0.86 percentage points |

All 10,936 exactly shared calls retained the same annotation category because category assignment depends on fixed genomic annotations rather than sequencing depth.

## Gene-Level Stability

- Full-depth genes: 6,420
- Downsampled genes: 4,203
- Shared genes: 3,757
- Gene-presence Jaccard: 0.547
- Zero-filled loop-count/rank Spearman rho: 0.485

Gene-level stability is moderate rather than strong. Gene ranking and GO analysis must remain exploratory and should not be presented as validation or as a stable strain-specific biological result.

## Source-Depth Influence

- Source contact depth versus pooled calls supported per library: rho 0.794 at full depth and -0.115 after downsampling.
- Full-depth calls supported only by above-median-depth libraries had 13.65% exact recovery in the downsampled pool.
- Other full-depth calls had 57.96% exact recovery.

This confirms that the original full-depth pooled resource was materially influenced by high-depth libraries.

## Defensible Interpretation

The analysis directly addresses the request to normalize libraries to comparable contact depth. It shows that downsampling substantially reduces technical loop-count variability and removes the positive association between original contact depth and loop yield. It also shows that many individual exact calls and gene rankings are depth-sensitive.

Therefore:

1. Do not claim conserved chromatin topology, strain-specific biological divergence, a consensus strain map, or genetic-distance effects.
2. Retain the 31,021 full-depth exact call records as a pooled union resource only if source-library support and HiCCUPS provenance are preserved.
3. Report the 14,414-call downsampled union and recovery metrics as a prominent sensitivity analysis.
4. Treat broad annotation categories as comparatively stable, but individual loop and gene results as depth-sensitive.
5. State explicitly that one library per strain prevents biological strain inference; downsampling controls depth but does not create replication.

## Files

- Analysis script: `r_files/downsampling_140M.R`
- Output directory: `r_files/downsampling_140M_outputs/`
- Key tables: `downsampling_loop_count_dispersion.tsv`, `downsampling_pooled_exact_recovery.tsv`, `downsampling_category_composition_change.tsv`, `downsampling_gene_count_stability_summary.tsv`, and `downsampling_pool_influence_correlations.tsv`

