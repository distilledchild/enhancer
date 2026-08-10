# NAR Resubmission Analysis Review

> **Authoritative status as of 2026-07-30**
> This document supersedes earlier interim reviews. It integrates the reviewer
> comments, the current reanalysis code, the official outputs, and the analytical
> principles finalized during subsequent discussions.

## 1. Final Study Positioning

The study should be reframed as:

> **A pooled rat frontal cortex chromatin-loop annotation resource with layered structural and putative regulatory evidence.**

The central principle is not to present one final high-confidence promoter-enhancer
set as a biological ground truth. Instead, the pooled Hi-C call resource is annotated
with promoter/TSS evidence, ATAC-seq accessibility, predicted CTCF motif intervals,
resolution, loop distance, HiCCUPS QC metrics, and source-library support.

The following claims should be avoided:

- Validated promoter-enhancer interactions
- High-confidence P-E interactions
- 31,021 unique or independent biological loops
- Consensus loops or conserved strain topology
- CTCF binding sites or unique CTCF sites
- Strain-specific biological differences or genetic-distance associations

## 2. Current Analysis Units and Fixed Counts

| Analysis unit or category | Count | Final interpretation |
|---|---:|---|
| Exact pooled HiCCUPS call records under 2 Mb | 31,021 | Primary resource and provenance unit; not the number of independent physical loops |
| Approximate cross-resolution loci | 20,335 | Gene-count sensitivity unit only |
| Direct promoter/TSS annotation at exactly one anchor | 12,295 | Directionally assignable promoter-side contacts |
| Direct promoter/TSS annotation at both anchors | 4,048 | Promoter-promoter-compatible contacts retained separately |
| No direct promoter/TSS annotation at either anchor | 14,678 | Unassigned or lower-regulatory-evidence resource calls |
| Single-promoter calls with at least 50 bp of opposite-anchor TSS-excluded ATAC overlap | 10,469 | Open-chromatin-supported putative promoter-to-distal-element contacts |
| Single-promoter calls without that ATAC support | 1,826 | Retained in the resource but excluded from the ATAC-supported subset |
| Genes assigned in the main putative set | 6,420 | Exact-call gene-count and exploratory GO input |
| Genes in the canonical-TSS sensitivity analysis | 5,712 | Annotation sensitivity input |

### Effect of changing the HiCCUPS input

| Processing step | Previous non-`sb` input | Current `sb` input | Difference |
|---|---:|---:|---:|
| Sample-level calls across 10 libraries | 58,992 | **59,000** | +8 |
| After exact coordinate/resolution pooling | 31,773 | **31,778** | +5 |
| After applying the `<2 Mb` filter | 31,019 | **31,021** | +2 |

If `12,295 + 1,094 = 13,389` is reported, it must be described only as the
**directionally assignable set**. It is not uniformly an ATAC-supported P-E set:
only 10,469 of the 12,295 single-promoter calls have opposite-anchor,
TSS-excluded ATAC support.

### Secondary annotation of 4,048 dual-promoter calls

| Category | Count | Definition and interpretation |
|---|---:|---|
| One-direction stringent asymmetric mixed-regulatory evidence | 1,094 | Only one promoter/TSS side has promoter-ATAC support, while the opposite anchor contains separate residual non-TSS ATAC; only one directional combination is complete |
| One-direction broader mixed promoter-regulatory evidence | 629 | Both promoter/TSS windows are ATAC-positive, but residual non-TSS ATAC occurs at only one anchor; only one directional combination is complete |
| Both directions supported | 1,184 | Both promoter-to-residual-ATAC combinations are possible; bidirectional or directionally ambiguous |
| Neither direction supported | 1,141 | No complete promoter-to-residual-ATAC direction; retained as promoter-promoter-compatible or lower-support contacts |

These four categories preserve dual-promoter calls as secondary annotations. They do
not reclassify the opposite anchor as a validated enhancer or establish regulatory
directionality.

## 3. Status of Reviewer-Facing Concerns

| Reviewer-facing issue | Status | Current conclusion |
|---|---|---|
| BED/BEDPE/GTF coordinate systems and true TSS definition | Resolved | All inputs use rn7, chromosome prefixes, and 1-based inclusive intervals; TSSs are strand-aware transcript 5-prime starts |
| Midpoint-based promoter/TSS assignment | Resolved | Direct evidence uses overlap with the complete anchor interval |
| Gene-body or last-exon containment filter | Resolved | Removed from selection and retained only as a descriptive flag |
| Forced P-E direction for dual-promoter calls | Resolved | All calls are retained with secondary ATAC-pattern categories |
| CTCF >=6 dominating P-E selection | Resolved | Removed from selection, categories, confidence labels, and GO inputs; CTCF is a predicted motif annotation only |
| Treating 31,021 calls as independent biological loops | Resolved | Estimand is explicitly exact, resolution-specific call records; 20,335 approximate loci provide sensitivity analysis |
| Loss of HiCCUPS provenance and QC fields | Resolved | Observed, expected, FDR, numCollapsed, centroid, radius, and source support are retained |
| Lack of an ATAC matched-null analysis | Resolved | A 1,000-permutation matched Hi-C-anchor null and resolution-stratified sensitivity analyses are complete |
| Overinterpretation of ATAC as enhancer validation | Managed through claim restriction | ATAC is used only as orthogonal open-chromatin annotation |
| Cross-resolution instability in gene ranking | Resolved analytically | Exact-call main analysis plus approximate-locus and canonical-TSS sensitivities |
| Previous Top-54 and GO overinterpretation | Resolved | Top-54 analysis removed; GO is exploratory and uses complete gene-count inputs |
| Sequencing-depth dependence | Resolved analytically; remains a study limitation | All ten libraries were downsampled to 140M MAPQ >=30 contacts, `.hic` files were rebuilt, and HiCCUPS was rerun at 5/10/25 kb; results confirm substantial depth sensitivity and preclude strain-specific inference |
| Strain-level genetic-distance analysis | Excluded from scope | No strain-specific claim is made, and replication plus a validated VCF-based distance matrix are unavailable |
| Lack of sample-matched functional validation | Residual data limitation | External ATAC cannot validate enhancer activity or target-gene regulation |
| Gene-biotype scope | Decision required | The resource preserves all biotypes; whether to add a protein-coding-only GO sensitivity remains to be decided |
| Sample, sex, and anatomical metadata | Incomplete | Must be completed in the revised Methods and response |
| Production workflow and official output freeze | Resolved | Legacy production dependencies were removed; official outputs, checksums, and source provenance are synchronized |

## 4. Resolved Analytical Components

### 4.1 Coordinate System, TSS Definition, and Promoter Assignment

- All analytical coordinates are normalized to rn7, chromosome-prefixed,
  1-based inclusive intervals.
- Ensembl transcript TSSs are strand-aware 5-prime transcript starts:
  `transcript_start` on the plus strand and `transcript_end` on the minus strand.
- Ensembl true TSS and EPD promoter annotations are both retained.
- Direct evidence is defined by overlap with complete `x1-x2` or `y1-y2` anchor
  intervals rather than anchor midpoints.
- All direct loop-anchor-gene assignments are preserved; one gene is not forced per loop.
- Inward-facing proximity within 10 kb is a secondary tier, while 200 kb proximity is exploratory.
- Transcript, last-exon, or whole-gene containment is an annotation flag, not a selection filter.

This component is scientifically defensible and should be described explicitly to reviewers.

### 4.2 Role of CTCF

The analysis retains results from the original 100-PWM `fimo.4.tsv` scan. The reader
preserves motif ID, alternative ID, score, p-value, q-value, strand, matched sequence,
and multiplicity metadata.

The permitted description of the 3,191,859 records is:

> **3,191,859 exact-coordinate-deduplicated predicted CTCF motif intervals generated from 100 PWM models.**

Required interpretation rules are:

- Do not call them independent CTCF loci, unique CTCF sites, or in-vivo binding events.
- Do not interpret overlapping PWM predictions as separate biological sites.
- Do not use `CTCF >=6` for loop retention, P-E direction, major categories,
  confidence labels, or GO gene-set selection.
- Retain CTCF only as predicted structural sequence annotation.
- For figures, use wording such as `predicted CTCF motif intervals showed positional concentration near loop anchors`.
- Without a formal matched background, use `positional concentration` or `density peaks`, not `enrichment`.
- A separate MA0139.1, q-value <=0.05 analysis may be reported as optional sensitivity,
  but it is not required to replace the original 100-PWM annotation.

The CTCF issue is resolved in the conceptual framing, terminology, code, and official
outputs. No official category or gene set now depends on a CTCF threshold.

### 4.3 Exact Calls and Approximate Loci

The primary unit is defined as:

> **31,021 exact, resolution-specific pooled HiCCUPS call records under 2 Mb.**

Source-library calls are pooled only when they have the same resolution, chromosome,
and anchor coordinates. Supporting libraries and strains are retained. This unit is
therefore a reproducible HiCCUPS call record, not an asserted independent physical
chromatin-loop locus.

The 20,335 approximate cross-resolution loci group nearby or overlapping calls across
resolutions. They are used only to assess gene-count robustness:

- 31,021 exact calls: primary resource and provenance
- 20,335 approximate loci: gene-count sensitivity
- Main gene metric: number of exact loop calls per gene
- Same-resolution shifted-call clustering: optional, not required for the current design

When several exact calls form one approximate locus, the locus is not represented by a
single chosen 5, 10, or 25 kb call. It functions as a cluster identifier linking all
member calls; gene assignments are collapsed to unique locus-gene pairs for sensitivity.

HiCCUPS observed counts, four expected/FDR metrics, numCollapsed, centroid, and radius
are retained as source-call provenance and descriptive QC. They are not used to create
an arbitrary new confidence filter.

### 4.4 ATAC Matched-Null Analysis

The final analysis uses the Duttke et al. 2022 rat prefrontal-cortex snATAC peak catalog,
aligned to rn7 as `Duttke2022_snATAC_peaks_rn7.narrowPeak`. The Yuan et al. SD-rat ATAC
dataset is not used in the final analysis.

For the 12,295 single-promoter calls, the opposite anchor was selected without using
the ATAC outcome. Controls were matched from pooled Hi-C anchors by chromosome,
resolution/anchor width, anchor side, and loop-distance bin. Exact-stratum controls
were available for 12,294 of 12,295 calls; one call required distance-bin relaxation.

Results from 1,000 permutations are:

| Resolution | Observed >=50 bp | Matched null | Absolute difference | Enrichment ratio | Empirical P |
|---|---:|---:|---:|---:|---:|
| 5 kb | 77.04% | 73.28% | 3.75 percentage points | 1.051 | 0.001 |
| 10 kb | 83.57% | 81.35% | 2.22 percentage points | 1.027 | 0.001 |
| 25 kb | 89.49% | 88.66% | 0.83 percentage points | 1.009 | 0.019 |
| All | 85.15% | 83.30% | 1.85 percentage points | 1.022 | 0.001 |

The null-model concern is resolved, but the effect is modest and smallest at 25 kb.
The appropriate conclusion is:

> ATAC-seq provides modest, resolution-dependent open-chromatin support beyond matched Hi-C-anchor expectations.

ATAC cannot be used as functional validation of enhancer activity, target-gene
regulation, or P-E interactions. Sensitivities for at least 1 bp, at least 50 bp,
anchor-overlap fraction, and resolution strata should be reported together.

### 4.5 Gene Ranking and GO Analysis

- The previous Top-54 analysis is removed and should not appear in the revised manuscript.
- Loop-gene pairs are deduplicated before exact-call counts are calculated.
- Approximate-locus counts and canonical-TSS-only assignments are sensitivity analyses.
- Exact-call versus approximate-locus gene-rank Spearman rho is 0.844.
- All-annotation versus canonical-TSS rank Spearman rho is 0.735 for exact calls and 0.701 for approximate loci.
- The main GO input contains 6,420 of 10,925 observable genes; the canonical sensitivity contains 5,712 of 9,878 genes.
- GO provides exploratory biological context only. It does not validate an individual loop, enhancer, target gene, or regulatory direction.

## 5. Partially Resolved Issues and Residual Risks

### 5.1 Sequencing Depth and Library Replication

Technical downsampling is complete for all ten libraries. Each MAPQ >=30 contact set
was uniformly downsampled to exactly 140 million contacts, a target just below the
lowest original usable-contact count (141,993,330). `.hic` files were rebuilt and
HiCCUPS was rerun at 5, 10, and 25 kb before applying the same `<2 Mb` restriction.

The results confirm that sequencing depth was a major technical driver:

- The across-library coefficient of variation in total loop calls fell from 0.335 at
  full depth to 0.163 after downsampling, a 51.5% reduction.
- The across-library range fell from 6,017 to 1,360 calls, a 77.4% reduction.
- The Spearman correlation between original contact depth and total loop calls changed
  from 0.794 at full depth to -0.115 after downsampling.
- The corresponding correlation with the number of pooled calls supported by each
  library also changed from 0.794 to -0.115.

Downsampling also demonstrates substantial individual-call sensitivity:

- The downsampled pooled union contained 14,414 exact, resolution-specific calls,
  compared with 31,021 at full depth.
- Exactly 10,936 full-depth pooled calls were recovered (35.25%; pooled Jaccard 0.317).
- Pooled exact recovery was 11.91%, 32.29%, and 50.25% at 5, 10, and 25 kb,
  respectively.
- Median sample-level exact recovery was 22.01% overall; a HiCCUPS-radius matching
  sensitivity increased this to 26.25%.
- Full-depth calls supported only by above-median-depth libraries had 13.65% exact
  recovery, versus 57.96% for the remaining calls.

Broad annotation composition was more stable than exact call identity. The four main
category proportions changed by at most 4.66 percentage points, and every exactly
shared call retained the same category. However, gene-level stability was only
moderate: 3,757 genes were shared between 6,420 full-depth and 4,203 downsampled genes
(presence Jaccard 0.547), and the zero-filled loop-count/rank Spearman correlation was
0.485. Gene ranking and GO therefore remain exploratory.

Each strain is still represented by one Hi-C library. Downsampling controls contact
depth but cannot create biological replication or distinguish strain effects from
other library-specific effects. The revised manuscript must therefore remove claims
of conserved topology, biological divergence, consensus strain architecture, or
genetic-distance associations. The 31,021 full-depth pooled records may remain the
primary union resource because they preserve available calls and source provenance,
but the 14,414-call downsampled union, resolution-stratified recovery, and source-depth
effects must be reported as a prominent sensitivity analysis and limitation.

### 5.2 External ATAC and Functional-Evidence Limitations

The Duttke ATAC dataset is useful orthogonal rat-PFC annotation, but it is not matched
to the Hi-C libraries by strain, sex, age, cell composition, or experimental preparation.
Combining aggregate ATAC with bulk Hi-C can therefore annotate one call with signals
arising from different cell populations.

Sample-matched expression, H3K27ac/H3K4me1, CTCF ChIP-seq, and perturbation assays are
not available. These are genuine limitations and should be stated directly in the
Discussion and response letter.

### 5.3 Resolution Dependence

Anchors at 5, 10, and 25 kb have different widths. Wider anchors have more opportunity
to overlap a promoter/TSS or ATAC peak. The matched-null analysis controls much of this
for ATAC, but direct promoter/TSS assignment remains resolution-dependent.

Report the following separately by resolution:

- Direct promoter/TSS overlap rate
- Opposite-anchor TSS-excluded ATAC overlap rate
- Putative category count and proportion
- Matched-null effect size

### 5.4 Gene-Biotype Scope

The earlier Ensembl branch was restricted to `Ensembl_canonical + protein_coding`,
whereas the revised resource preserves all annotated biotypes. Retaining all biotypes
is appropriate for the resource. For downstream GO analysis, the observable gene
universe must be stated, and the team should decide whether to add a protein-coding-only
GO sensitivity analysis.

### 5.5 Methods Metadata

The revised Methods and response must provide:

- Biological sample and library counts, including the strain represented by each library
- Sex of each sample
- Anatomical boundaries of the frontal-cortex dissection or an accompanying collection figure
- A primary reference supporting the description of HRDP as genetically diverse
- Strain, tissue, age/sex, and coordinate-conversion provenance for the external ATAC dataset

### 5.6 Production Pipeline and Official Outputs

Previously identified inconsistencies between scientific decisions and code or outputs
have been addressed as follows:

- CTCF-based categories were removed. The mutually exclusive release contains 10,469
  putative regulatory calls, 4,048 promoter-promoter-compatible calls, 1,826
  single-promoter calls without opposite-anchor ATAC support, and 14,678 calls without
  direct promoter/TSS support, totaling 31,021.
- `CTCF >=6` was removed from loop selection, categories, highlighted subsets,
  directional assignment, and GO gene-set construction.
- The old-versus-new comparison script is no longer required by the production main script.
- Legacy candidate RDS, legacy final-loop CSV, legacy rn7 EPD BED, and start-codon RDS
  are no longer required for coordinate preparation.
- The coordinate cache is rebuilt from ten original HiCCUPS BEDPE files, full FIMO output,
  ATAC narrowPeak, Ensembl GTF, EPD rn6 source/coordinate/mapping/chain files, and depth-QC input.
- The 1,000-permutation ATAC outputs are frozen with seed `20260727`, run metadata,
  matched-control QC, session information, and figures.
- The ten-library 140M-contact analysis is complete and stored in
  `r_files/downsampling_140M_outputs/`, including sample status, resolution-stratified
  recovery, category composition, gene stability, source-depth influence, figures,
  run metadata, and session information.
- Progress documents use the current 31,021-call release and revised categories.
- `resubmit_output_manifest.tsv` records file sizes and SHA-256 hashes for 51 official files.
- `source_data_versions.tsv` records path, modification time, SHA-256, assembly,
  provenance, and analysis release for 22 required inputs and scripts; no source or hash is missing.

The manuscript and response letter must use only the counts and tables from this
release. Before submission, the downsampling script and compact outputs must be added
to the final checksum manifest and code release.

## 6. Wording Rules for the Manuscript and Response Letter

### Permitted wording

- Pooled rat frontal cortex chromatin-loop annotation resource
- Exact, resolution-specific pooled HiCCUPS call records
- Open-chromatin-supported putative promoter-to-distal-element contacts
- Promoter-promoter-compatible contacts
- Predicted CTCF motif intervals
- Layered structural and putative regulatory evidence
- Exploratory functional interpretation
- Library support or number of supporting libraries

### Wording to avoid

- Unique loops, independent loops, or consensus loops
- High-confidence P-E interactions
- Validated enhancers or validated target genes
- CTCF binding sites or CTCF-mediated P-E loops
- ATAC validation of enhancer function
- GO validation of loop biology
- Strain-specific regulation, conservation, or genetic effects

### Safe core result statements

> We assembled 31,021 exact, resolution-specific pooled HiCCUPS call records under 2 Mb and annotated them with strand-aware promoter/TSS assignments, open-chromatin support, predicted CTCF motif intervals, resolution, distance, HiCCUPS QC metrics, and source-library support.

> Among 12,295 calls with direct promoter/TSS annotation at exactly one anchor, 10,469 showed at least 50 bp of TSS-excluded ATAC overlap at the opposite anchor and were designated open-chromatin-supported putative promoter-to-distal-element contacts.

> Matched Hi-C-anchor permutations indicated statistically significant but modest ATAC support, so accessibility was interpreted as orthogonal annotation rather than functional validation.

> Equalizing all ten libraries to 140 million MAPQ-filtered contacts reduced the across-library coefficient of variation in loop-call counts from 0.335 to 0.163 and removed the positive association between original sequencing depth and loop yield (Spearman rho 0.794 at full depth versus -0.115 after downsampling). Exact pooled-call recovery was 35.25%, demonstrating that individual loop calls remain depth-sensitive even though broad annotation-category proportions were comparatively stable.

## 7. Remaining Priorities Before Resubmission

1. Integrate the completed downsampling methods, resolution-stratified results, and figures into the manuscript and supplement.
2. Add the downsampling script and compact outputs to the final Git release, checksum manifest, and source-provenance table.
3. Decide whether to add a protein-coding-only GO sensitivity analysis while retaining all biotypes in the resource.
4. Complete Methods metadata for samples, sex, strains, dissection, HRDP references, and ATAC provenance.
5. Rewrite the manuscript, response letter, figures, and supplements using the current release counts and framing.
6. Perform a final consistency audit for prohibited terminology and superseded counts.

## 8. Objective Final Assessment

The reanalysis is substantially stronger than the original submission. Coordinate and
TSS normalization, interval-based direct assignment, retention of dual-promoter calls,
removal of gene-body containment filtering, explicit exact-call estimands, preservation
of HiCCUPS provenance, the ATAC matched-null analysis, and approximate-locus plus
canonical-TSS sensitivity analyses are scientifically defensible.

The completed downsampling directly addresses the reviewer's normalization request and
shows why the revised claim restrictions are necessary. It strengthens the technical
defensibility of the pooled resource, but it does not rescue strain-specific biological
inference or establish that every full-depth call is robust. The remaining work is to
integrate these results, complete metadata and optional biotype sensitivity, freeze the
release, and rewrite the manuscript under the constrained resource framing.

Acceptance remains uncertain because individual loop and gene results are depth-sensitive,
biological replication and sample-matched functional validation are absent, and editorial
assessment of novelty remains a material risk.

## 9. Reference Files

- Main analysis: `r_files/revision/promoter_enhancer_interaction_resubmit.R`
- Coordinate preparation: `r_files/revision/promoter_enhancer_interaction_resubmit_01_coord_prep.R`
- Shared functions: `r_files/funcs.R`
- ATAC matched-null analysis: `r_files/atac_validation/atac_validation.R`
- Depth-normalization analysis: `r_files/downsampling_140M.R`
- Depth-normalization outputs: `r_files/downsampling_140M_outputs/`
- Official outputs: `r_files/revision/resubmit_outputs/`
- Reviewer comments: `Obsidian Vault/research/enhancer/revision/00_reviewer-comment.md`
