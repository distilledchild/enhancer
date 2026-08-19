################################################################################
# Legacy baseline analyses for explicit old-vs-new comparison
#
# This comparison-only file is sourced by promoter_enhancer_interaction_resubmit.R after the
# revised assignment, CTCF-sensitivity, and downstream objects are available.
# It intentionally shares the caller environment and is not a standalone run.
################################################################################

################################################################################
# 8-1. Legacy promoter/TSS gene-annotation baseline
################################################################################

# Reproduce the previous one-assignment-per-loop promoter/TSS evidence from
# the coordinate-normalized legacy object prepared in Section 1-6. This is
# retained only as a baseline for the later old-vs-new comparison and will not
# define the revised direct or proximal assignment tiers.

df.promoter.tss.evidence <- df.promoter.tss.candidate %>%
  transmute(
    loop_id,
    promoter_or_tss_support,
    promoter_support,
    tss_support,
    promoter_tss_component = component,
    promoter_tss_distance = as.integer(distance),
    promoter_tss_gene_id = gene_id,
    promoter_tss_gene_name = gene_name,
    promoter_anchor_side,
    candidate_regulatory_anchor_side,
    promoter_tss_where = WHERE,
    promoter_tss_classification = classification,
    passes_promoter_tss_200kb = (
      promoter_or_tss_support &
        promoter_tss_distance < 200000L
    )
  )

################################################################################
# 8-2. Legacy ATAC open-chromatin baseline
#
# The candidate-regulatory anchor is the loop anchor opposite the selected
# promoter/TSS anchor. ATAC overlap is supportive evidence of accessibility.
# It is not direct evidence that the anchor functions as an enhancer.
################################################################################

# Use the normalized ATAC object prepared in Section 1-5.
gr.atac.union <- GenomicRanges::reduce(gr.atac)

gr.promoter.anchor <- create_anchor_granges(
  df.promoter.tss.candidate,
  "promoter"
)
gr.candidate.regulatory.anchor <- create_anchor_granges(
  df.promoter.tss.candidate,
  "candidate_regulatory"
)

# Verify one promoter and one candidate-regulatory anchor per legacy loop.
assert_analysis_condition(
  length(gr.promoter.anchor) == nrow(df.promoter.tss.candidate) &&
    length(gr.candidate.regulatory.anchor) ==
      nrow(df.promoter.tss.candidate),
  paste0(
    "Anchor construction did not preserve one promoter and one ",
    "candidate-regulatory anchor per candidate loop."
  )
)

promoter.atac.hit <- countOverlaps(
  gr.promoter.anchor,
  gr.atac.union,
  minoverlap = 50
) > 0

candidate.regulatory.atac.hit <- countOverlaps(
  gr.candidate.regulatory.anchor,
  gr.atac.union,
  minoverlap = 50
) > 0

# Reproduce the previous non-TSS ATAC definition with the normalized legacy
# start-codon-derived object from Section 1-6. The revised ATAC calculation will
# instead use true TSS coordinates generated in Step 2.
gr.legacy.start.codon <- GRanges(
  seqnames = df.legacy.start.codon.annotation$chr,
  ranges = IRanges(
    start = df.legacy.start.codon.annotation$start,
    end = df.legacy.start.codon.annotation$end
  ),
  strand = df.legacy.start.codon.annotation$strand
)

gr.legacy.start.codon.region <- promoters(
  gr.legacy.start.codon,
  upstream = 1000,
  downstream = 1000
) %>%
  GenomicRanges::reduce(ignore.strand = TRUE)

common.seqlevels.atac.tss <- intersect(
  seqlevels(gr.atac.union),
  seqlevels(gr.legacy.start.codon.region)
)

gr.atac.for.tss <- keepSeqlevels(
  gr.atac.union,
  common.seqlevels.atac.tss,
  pruning.mode = "coarse"
)
gr.legacy.start.codon.region.for.atac <- keepSeqlevels(
  gr.legacy.start.codon.region,
  common.seqlevels.atac.tss,
  pruning.mode = "coarse"
)

gr.atac.non.tss <- GenomicRanges::setdiff(
  gr.atac.for.tss,
  gr.legacy.start.codon.region.for.atac,
  ignore.strand = TRUE
)

tss.anchor.hit <- findOverlaps(
  gr.candidate.regulatory.anchor,
  gr.legacy.start.codon.region,
  ignore.strand = TRUE
)

tss.bp.by.candidate.anchor <- numeric(
  length(gr.candidate.regulatory.anchor)
)

# Accumulate legacy TSS-overlap bases only when overlaps are present.
if (length(tss.anchor.hit) > 0) {
  gr.tss.anchor.overlap <- pintersect(
    gr.candidate.regulatory.anchor[queryHits(tss.anchor.hit)],
    gr.legacy.start.codon.region[subjectHits(tss.anchor.hit)],
    ignore.strand = TRUE
  )

  tss.bp.sum <- rowsum(
    width(gr.tss.anchor.overlap),
    group = queryHits(tss.anchor.hit),
    reorder = FALSE
  )

  tss.bp.by.candidate.anchor[
    as.integer(rownames(tss.bp.sum))
  ] <- tss.bp.sum[, 1]
}

candidate.anchor.has.non.tss.fragment <- (
  tss.bp.by.candidate.anchor <
    width(gr.candidate.regulatory.anchor)
)

candidate.regulatory.non.tss.atac.hit <- (
  candidate.anchor.has.non.tss.fragment &
    (
      countOverlaps(
        gr.candidate.regulatory.anchor,
        gr.atac.non.tss,
        minoverlap = 50
      ) > 0
    )
)

df.atac.evidence <- tibble(
  loop_id = gr.promoter.anchor$loop_id,
  atac_promoter_anchor = promoter.atac.hit,
  atac_candidate_regulatory_anchor = candidate.regulatory.atac.hit,
  candidate_regulatory_anchor_has_non_tss_fragment =
    candidate.anchor.has.non.tss.fragment,
  candidate_regulatory_anchor_non_tss_bp =
    width(gr.candidate.regulatory.anchor) -
      tss.bp.by.candidate.anchor,
  atac_candidate_regulatory_anchor_non_tss =
    candidate.regulatory.non.tss.atac.hit
)

################################################################################
# 8-3. Previous final set and legacy loop-level evidence table
################################################################################

df.current.final.loop <- readr::read_csv(
  legacy.current.final.loop.file,
  show_col_types = FALSE
)
check_required_columns(
  df.current.final.loop,
  c("loop.id", "category"),
  "df.current.final.loop"
)

df.current.final.loop <- df.current.final.loop %>%
  dplyr::rename(loop_id = loop.id) %>%
  transmute(
    loop_id,
    current_final_set = TRUE,
    current_final_category = category
  )

df.loop.evidence <- df.loop.universe %>%
  left_join(
    df.ctcf.evidence %>% dplyr::select(-resolution),
    by = "loop_id"
  ) %>%
  left_join(df.promoter.tss.evidence, by = "loop_id") %>%
  left_join(df.atac.evidence, by = "loop_id") %>%
  left_join(df.current.final.loop, by = "loop_id") %>%
  mutate(
    promoter_or_tss_support = coalesce(
      promoter_or_tss_support,
      FALSE
    ),
    promoter_support = coalesce(promoter_support, FALSE),
    tss_support = coalesce(tss_support, FALSE),
    passes_promoter_tss_200kb = coalesce(
      passes_promoter_tss_200kb,
      FALSE
    ),
    atac_promoter_anchor = coalesce(
      atac_promoter_anchor,
      FALSE
    ),
    atac_candidate_regulatory_anchor = coalesce(
      atac_candidate_regulatory_anchor,
      FALSE
    ),
    candidate_regulatory_anchor_has_non_tss_fragment = coalesce(
      candidate_regulatory_anchor_has_non_tss_fragment,
      FALSE
    ),
    candidate_regulatory_anchor_non_tss_bp = coalesce(
      candidate_regulatory_anchor_non_tss_bp,
      0
    ),
    atac_candidate_regulatory_anchor_non_tss = coalesce(
      atac_candidate_regulatory_anchor_non_tss,
      FALSE
    ),
    current_final_set = coalesce(current_final_set, FALSE),
    current_final_category = replace_na(
      current_final_category,
      "not_current_final"
    ),
    ensembl_gene_id = str_extract(
      promoter_tss_gene_id,
      "ENSRNOG[0-9]+"
    ),
    gene_symbol = promoter_tss_gene_name,
    putative_regulatory_support = (
      passes_promoter_tss_200kb &
        atac_candidate_regulatory_anchor_non_tss
    ),
    putative_regulatory_support_any_atac = (
      passes_promoter_tss_200kb &
        atac_candidate_regulatory_anchor
    ),
    strongest_candidate_support = (
      putative_regulatory_support &
        passes_ctcf_ge6_both_anchors
    ),
    proposed_major_category = case_when(
      putative_regulatory_support ~
        "putative_regulatory_promoter_TSS_nonTSS_ATAC_supported",
      passes_ctcf_ge6_both_anchors ~
        "structural_CTCF_supported_not_ATAC_regulatory",
      TRUE ~ "lower_support_or_uncertain"
    ),
    proposed_detailed_category = case_when(
      strongest_candidate_support ~
        "strongest_candidate_PE_promoter_TSS_nonTSS_ATAC_CTCF",
      putative_regulatory_support &
        !passes_ctcf_ge6_both_anchors ~
        "putative_regulatory_promoter_TSS_nonTSS_ATAC_no_strong_CTCF",
      promoter_or_tss_support &
        atac_candidate_regulatory_anchor &
        !atac_candidate_regulatory_anchor_non_tss &
        passes_ctcf_ge6_both_anchors ~
        "CTCF_promoter_TSS_ATAC_signal_TSS_sensitive",
      promoter_or_tss_support &
        atac_candidate_regulatory_anchor &
        !atac_candidate_regulatory_anchor_non_tss ~
        "promoter_TSS_ATAC_signal_TSS_sensitive",
      promoter_or_tss_support &
        passes_ctcf_ge6_both_anchors ~
        "structural_CTCF_promoter_TSS_no_nonTSS_ATAC",
      passes_ctcf_ge6_both_anchors ~
        "structural_CTCF_only_no_promoter_TSS_nonTSS_ATAC",
      promoter_or_tss_support ~
        "promoter_TSS_only_no_nonTSS_ATAC_no_strong_CTCF",
      ctcf_both_anchors_any ~
        "weak_structural_CTCF_both_anchors_below_ge6",
      TRUE ~ "low_support_or_uncertain"
    )
  ) %>%
  arrange(chr1, start1, end1, chr2, start2, end2)

# Confirm that loop-level joins preserve the pooled resource row count.
assert_analysis_row_count(
  df.loop.evidence,
  nrow(df.loop.universe),
  "Loop-level joins changed the pooled loop-resource row count."
)

# Require the legacy evidence table to contain one unique row per loop.
assert_analysis_unique_key(
  df.loop.evidence,
  "loop_id",
  "The loop-level evidence table is not one row per loop."
)

# Confirm that the conservative subset remains inside its parent category.
assert_analysis_condition(
  !any(
    df.loop.evidence$strongest_candidate_support &
      df.loop.evidence$proposed_major_category !=
        "putative_regulatory_promoter_TSS_nonTSS_ATAC_supported"
  ),
  paste0(
    "The conservative highlighted subset is not contained in ",
    "the putative-regulatory category."
  )
)

# Cross-tabulate the previous and revised mutually exclusive labels on the same
# rebuilt pooled universe. This table is descriptive and makes category migration
# auditable; it does not imply that either label is experimental ground truth.
df.legacy.revised.category.crosswalk <-
  df.loop.evidence %>%
  dplyr::select(loop_id, resolution, proposed_major_category) %>%
  left_join(
    df.revised.loop.evidence %>%
      dplyr::select(loop_id, revised_major_category),
    by = "loop_id"
  ) %>%
  count(
    resolution,
    proposed_major_category,
    revised_major_category,
    name = "n_loops"
  ) %>%
  group_by(resolution, proposed_major_category) %>%
  mutate(
    pct_within_legacy_category = round(
      100 * n_loops / sum(n_loops),
      1
    )
  ) %>%
  ungroup() %>%
  arrange(
    resolution,
    proposed_major_category,
    desc(n_loops),
    revised_major_category
  )

################################################################################
# 9. Legacy filtering audit and category summaries
################################################################################

n.all.distinct <- n_distinct(
  df.loop.distinct$loop_id
)
n.lt2mb <- nrow(df.loop.universe)
n.ctcf.any.both <- count_true(df.loop.evidence$ctcf_both_anchors_any)
n.ctcf.ge6.both <- count_true(
  df.loop.evidence$passes_ctcf_ge6_both_anchors
)
n.promoter.tss <- count_true(
  df.loop.evidence$promoter_or_tss_support
)
n.promoter.tss.200kb <- count_true(
  df.loop.evidence$passes_promoter_tss_200kb
)
n.current.final <- count_true(df.loop.evidence$current_final_set)
n.candidate.any.atac <- count_true(
  df.loop.evidence$putative_regulatory_support_any_atac
)
n.candidate.non.tss.atac <- count_true(
  df.loop.evidence$putative_regulatory_support
)
n.strongest.candidate <- count_true(
  df.loop.evidence$strongest_candidate_support
)

df.filtering.logic.audit <- tribble(
  ~filter_step,
  ~code_source,
  ~evidence_type,
  ~input_n,
  ~output_n,
  ~filter_condition,
  ~interpretation,
  ~revision_action,
  "raw_distinct_hiccups_loops",
  "enhancer_promoter_interaction_figures.R: loop preprocessing",
  "technical",
  n.all.distinct,
  n.all.distinct,
  "Load distinct HiCCUPS loops across all 10 libraries.",
  "Pooled physical loop universe before revision filters.",
  "Retain as the starting resource universe.",
  "loop_distance_lt_2mb",
  "enhancer_promoter_interaction.R: Section 1",
  "technical",
  n.all.distinct,
  n.lt2mb,
  "loop distance < 2,000,000 bp",
  "Technical scope filter; not regulatory evidence.",
  "Retain and describe as a loop-length scope threshold.",
  "ctcf_any_motif_both_anchors",
  "enhancer_promoter_interaction.R: Section 2",
  "CTCF-driven",
  n.lt2mb,
  n.ctcf.any.both,
  "At least one predicted CTCF motif at both loop anchors.",
  "Structural CTCF support; not direct enhancer evidence.",
  "Use as structural annotation.",
  "ctcf_ge6_both_anchors_previous_filter",
  "enhancer_promoter_interaction.R: CTCF threshold filtering",
  "CTCF-driven",
  n.lt2mb,
  n.ctcf.ge6.both,
  "At least 6 predicted CTCF motifs at both loop anchors.",
  "Strong structural motif support under the previous threshold.",
  "Show threshold sensitivity; do not use alone as regulatory proof.",
  "promoter_or_tss_directional_annotation_200kb",
  "enhancer_promoter_interaction.R: directional gene assignment",
  "promoter/TSS-driven",
  n.lt2mb,
  n.promoter.tss.200kb,
  "Selected promoter/TSS is within 200 kb of its assigned anchor.",
  "Target-gene annotation support; not enhancer activity evidence.",
  "Retain as promoter/TSS evidence.",
  "current_final_loop_set",
  "enhancer_promoter_interaction.R: previous final intersection",
  "combined_CTCF_promoter_TSS",
  n.lt2mb,
  n.current.final,
  "CTCF >=6 at both anchors plus promoter/TSS annotation.",
  "Structurally strong but overclaimed if called validated P-E interactions.",
  "Reframe as the previous CTCF-supported promoter/TSS set.",
  "candidate_regulatory_anchor_ATAC_any",
  "atac_validation.R: paired anchor overlap",
  "ATAC-driven",
  n.promoter.tss,
  n.candidate.any.atac,
  "Opposite anchor overlaps an ATAC peak by at least 50 bp.",
  "Open-chromatin support at the candidate-regulatory anchor.",
  "Use as supportive, not functional, evidence.",
  "candidate_regulatory_anchor_nonTSS_ATAC",
  "atac_validation.R: TSS-exclusion sensitivity analysis",
  "ATAC-driven",
  n.promoter.tss,
  n.candidate.non.tss.atac,
  "Opposite anchor overlaps ATAC signal remaining after TSS +/-1 kb removal.",
  "Conservative accessibility support not explained only by known TSS signal.",
  "Use to define the putative-regulatory evidence category."
)

df.loop.evidence.summary <- bind_rows(
  df.loop.evidence %>%
    count(proposed_major_category, name = "n_loops") %>%
    mutate(
      summary_type = "proposed_major_category",
      category = proposed_major_category
    ) %>%
    dplyr::select(summary_type, category, n_loops),
  df.loop.evidence %>%
    count(proposed_detailed_category, name = "n_loops") %>%
    mutate(
      summary_type = "proposed_detailed_category",
      category = proposed_detailed_category
    ) %>%
    dplyr::select(summary_type, category, n_loops),
  df.loop.evidence %>%
    count(current_final_set, name = "n_loops") %>%
    mutate(
      summary_type = "current_final_set",
      category = as.character(current_final_set)
    ) %>%
    dplyr::select(summary_type, category, n_loops),
  df.loop.evidence %>%
    count(
      current_final_category,
      proposed_major_category,
      name = "n_loops"
    ) %>%
    mutate(
      summary_type = "current_category_by_new_major_category",
      category = paste(
        current_final_category,
        proposed_major_category,
        sep = " | "
      )
    ) %>%
    dplyr::select(summary_type, category, n_loops)
) %>%
  arrange(summary_type, desc(n_loops), category)

df.loop.category.by.resolution <- bind_rows(
  df.loop.evidence %>%
    count(
      resolution,
      proposed_major_category,
      name = "n_loops"
    ) %>%
    group_by(resolution) %>%
    mutate(
      pct_loops = round(100 * n_loops / sum(n_loops), 1)
    ) %>%
    ungroup(),
  df.loop.evidence %>%
    count(proposed_major_category, name = "n_loops") %>%
    mutate(
      resolution = "ALL",
      pct_loops = round(100 * n_loops / sum(n_loops), 1)
    )
) %>%
  arrange(resolution, desc(n_loops))

category.count.check <- df.loop.evidence %>%
  count(proposed_major_category, name = "n_loops")

# Validate that the three legacy major categories partition all pooled loops.
assert_analysis_condition(
  sum(category.count.check$n_loops) == n.lt2mb &&
    nrow(category.count.check) == 3L,
  paste0(
    "The three major categories do not form a complete, mutually ",
    "exclusive partition of the pooled resource."
  )
)

df.category.definition <- tribble(
  ~resource_group,
  ~required_evidence,
  ~excluded_or_not_required,
  ~allowed_interpretation,
  "pooled_loop_annotation_resource",
  "Distinct HiCCUPS loop with distance <2 Mb.",
  "No CTCF, promoter/TSS, or ATAC evidence is required.",
  "Pooled rat frontal cortex chromatin-loop annotation resource.",
  "putative_regulatory_promoter_TSS_nonTSS_ATAC_supported",
  "Promoter/TSS assignment within 200 kb and non-TSS ATAC support at the opposite anchor.",
  "Strong CTCF support is not required.",
  "Putative regulatory loop with layered annotation and accessibility support.",
  "structural_CTCF_supported_not_ATAC_regulatory",
  "At least 6 predicted CTCF motifs at both anchors.",
  "Does not satisfy promoter/TSS plus non-TSS ATAC criteria.",
  "Structurally supported loop; not direct evidence of enhancer function.",
  "lower_support_or_uncertain",
  "Member of pooled loop resource.",
  "Neither of the two higher-evidence category definitions is satisfied.",
  "Lower-support or uncertain loop retained for resource completeness.",
  "highlighted_conservative_subset",
  "Promoter/TSS assignment, opposite-anchor non-TSS ATAC, and CTCF >=6 at both anchors.",
  "This is nested within the putative-regulatory category.",
  "More conservative highlighted subset; not a separate fourth category."
)


################################################################################
# 10-2. Legacy CTCF threshold sensitivity baseline
################################################################################

create_ctcf_sensitivity_row <- function(df, threshold, resolution.label) {
  passes.threshold <- (
    df$ctcf_count_anchor1 >= threshold &
      df$ctcf_count_anchor2 >= threshold
  )
  n.passes <- sum(passes.threshold)

  tibble(
    resolution = resolution.label,
    ctcf_min_threshold = threshold,
    n_loops_universe = nrow(df),
    n_passes_ctcf_thr = n.passes,
    pct_passes_ctcf_thr = round(100 * mean(passes.threshold), 1),
    n_ctcf_and_promoter_tss = sum(
      passes.threshold &
        df$promoter_or_tss_support
    ),
    pct_ctcf_and_promoter_tss_of_ctcf = round(
      100 * sum(
        passes.threshold &
          df$promoter_or_tss_support
      ) / max(n.passes, 1),
      1
    ),
    n_ctcf_promoter_tss_nonTSS_ATAC = sum(
      passes.threshold &
        df$putative_regulatory_support
    ),
    pct_ctcf_promoter_tss_nonTSS_ATAC_of_ctcf = round(
      100 * sum(
        passes.threshold &
          df$putative_regulatory_support
      ) / max(n.passes, 1),
      1
    ),
    n_putative_regulatory = sum(df$putative_regulatory_support)
  )
}

df.ctcf.sensitivity.by.threshold <- map_dfr(
  ctcf.thresholds,
  function(threshold.i) {
    bind_rows(
      map_dfr(
        sort(unique(df.loop.evidence$resolution)),
        function(resolution.i) {
          create_ctcf_sensitivity_row(
            df.loop.evidence %>%
              filter(resolution == resolution.i),
            threshold.i,
            resolution.i
          )
        }
      ),
      create_ctcf_sensitivity_row(
        df.loop.evidence,
        threshold.i,
        "ALL"
      )
    )
  }
) %>%
  arrange(ctcf_min_threshold, resolution)

df.ctcf.threshold.summary <- bind_rows(
  df.loop.evidence %>%
    group_by(resolution) %>%
    summarise(
      n_loops = n(),
      n_ctcf_any_both = sum(ctcf_both_anchors_any),
      pct_ctcf_any_both = percent_true(ctcf_both_anchors_any),
      n_ctcf_ge6_both = sum(passes_ctcf_ge6_both_anchors),
      pct_ctcf_ge6_both = percent_true(
        passes_ctcf_ge6_both_anchors
      ),
      median_min_ctcf_count = stats::median(
        ctcf_min_anchor_count
      ),
      q1_min_ctcf_count = stats::quantile(
        ctcf_min_anchor_count,
        0.25
      ),
      q3_min_ctcf_count = stats::quantile(
        ctcf_min_anchor_count,
        0.75
      ),
      .groups = "drop"
    ),
  df.loop.evidence %>%
    summarise(
      resolution = "ALL",
      n_loops = n(),
      n_ctcf_any_both = sum(ctcf_both_anchors_any),
      pct_ctcf_any_both = percent_true(ctcf_both_anchors_any),
      n_ctcf_ge6_both = sum(passes_ctcf_ge6_both_anchors),
      pct_ctcf_ge6_both = percent_true(
        passes_ctcf_ge6_both_anchors
      ),
      median_min_ctcf_count = stats::median(
        ctcf_min_anchor_count
      ),
      q1_min_ctcf_count = stats::quantile(
        ctcf_min_anchor_count,
        0.25
      ),
      q3_min_ctcf_count = stats::quantile(
        ctcf_min_anchor_count,
        0.75
      )
    )
) %>%
  arrange(resolution)

# The resolution-adjusted rule holds the motif density constant relative to
# the legacy 5-kb >=6 rule: 5 kb >=6, 10 kb >=12, and 25 kb >=30 motifs at
# both anchors. This is a technical sensitivity analysis of motif burden, not
# a biological threshold for CTCF binding or enhancer activity.
summarise_ctcf_resolution_adjustment <- function(df, resolution.label) {
  fixed.pass <- df$passes_ctcf_ge6_both_anchors
  adjusted.pass <- df$passes_ctcf_resolution_adjusted_both_anchors
  putative.pass <- df$putative_regulatory_support
  adjusted.rule <- if (resolution.label == "ALL") {
    "5K>=6; 10K>=12; 25K>=30"
  } else {
    paste0(
      resolution.label,
      ">=",
      unique(df$ctcf_resolution_adjusted_min_threshold)
    )
  }

  tibble(
    resolution = resolution.label,
    resolution_adjusted_rule = adjusted.rule,
    n_loops_universe = nrow(df),
    n_fixed_ge6_both = sum(fixed.pass),
    pct_fixed_ge6_both = round(100 * mean(fixed.pass), 1),
    n_resolution_adjusted_both = sum(adjusted.pass),
    pct_resolution_adjusted_both = round(
      100 * mean(adjusted.pass),
      1
    ),
    n_pass_both_rules = sum(fixed.pass & adjusted.pass),
    n_fixed_ge6_only = sum(fixed.pass & !adjusted.pass),
    n_resolution_adjusted_only = sum(!fixed.pass & adjusted.pass),
    n_fail_both_rules = sum(!fixed.pass & !adjusted.pass),
    n_putative_regulatory = sum(putative.pass),
    n_putative_regulatory_fixed_ge6 = sum(
      putative.pass & fixed.pass
    ),
    pct_putative_regulatory_fixed_ge6 = round(
      100 * sum(putative.pass & fixed.pass) /
        max(sum(putative.pass), 1),
      1
    ),
    n_putative_regulatory_resolution_adjusted = sum(
      putative.pass & adjusted.pass
    ),
    pct_putative_regulatory_resolution_adjusted = round(
      100 * sum(putative.pass & adjusted.pass) /
        max(sum(putative.pass), 1),
      1
    )
  )
}

df.ctcf.resolution.adjusted.summary <- bind_rows(
  map_dfr(
    sort(unique(df.loop.evidence$resolution)),
    function(resolution.i) {
      summarise_ctcf_resolution_adjustment(
        df.loop.evidence %>%
          filter(resolution == resolution.i),
        resolution.i
      )
    }
  ),
  summarise_ctcf_resolution_adjustment(
    df.loop.evidence,
    "ALL"
  )
) %>%
  arrange(resolution)

df.ctcf.resolution.adjusted.loop.comparison <- df.loop.evidence %>%
  transmute(
    loop_id,
    resolution,
    ctcf_count_anchor1,
    ctcf_count_anchor2,
    ctcf_min_anchor_count,
    fixed_min_threshold = 6L,
    passes_ctcf_ge6_both_anchors,
    resolution_adjusted_min_threshold =
      ctcf_resolution_adjusted_min_threshold,
    passes_ctcf_resolution_adjusted_both_anchors,
    putative_regulatory_support,
    fixed_vs_adjusted_status = case_when(
      passes_ctcf_ge6_both_anchors &
        passes_ctcf_resolution_adjusted_both_anchors ~
        "passes_both_rules",
      passes_ctcf_ge6_both_anchors &
        !passes_ctcf_resolution_adjusted_both_anchors ~
        "passes_fixed_ge6_only",
      !passes_ctcf_ge6_both_anchors &
        passes_ctcf_resolution_adjusted_both_anchors ~
        "passes_resolution_adjusted_only",
      TRUE ~ "fails_both_rules"
    )
  ) %>%
  arrange(resolution, loop_id)

################################################################################
# 11. Legacy ATAC support by resolution and paired-anchor comparison
################################################################################

df.atac.candidate <- df.loop.evidence %>%
  filter(passes_promoter_tss_200kb)

df.atac.support.by.resolution <- bind_rows(
  map_dfr(
    sort(unique(df.atac.candidate$resolution)),
    function(resolution.i) {
      summarise_atac_support(
        df.atac.candidate %>%
          filter(resolution == resolution.i),
        resolution.i
      )
    }
  ),
  summarise_atac_support(df.atac.candidate, "ALL")
) %>%
  arrange(resolution)

# Preserve the original one-row output name for downstream compatibility.
df.atac.support.summary <- df.atac.support.by.resolution %>%
  filter(resolution == "ALL") %>%
  dplyr::select(-resolution)

df.atac.paired.anchor.mcnemar <- bind_rows(
  map_dfr(
    sort(unique(df.atac.candidate$resolution)),
    function(resolution.i) {
      summarise_mcnemar(
        df.atac.candidate %>%
          filter(resolution == resolution.i),
        resolution.i
      )
    }
  ),
  summarise_mcnemar(df.atac.candidate, "ALL")
) %>%
  arrange(resolution)


################################################################################
# 12-2. Legacy resource tables and downstream gene summaries
################################################################################

resource.columns <- c(
  "loop_id",
  "chr1", "start1", "end1",
  "chr2", "start2", "end2",
  "resolution", "loop_distance",
  "proposed_major_category",
  "proposed_detailed_category",
  "putative_regulatory_support",
  "strongest_candidate_support",
  "current_final_set",
  "current_final_category",
  "ctcf_count_anchor1",
  "ctcf_count_anchor2",
  "ctcf_min_anchor_count",
  "ctcf_both_anchors_any",
  "passes_ctcf_ge6_both_anchors",
  "promoter_or_tss_support",
  "promoter_support",
  "tss_support",
  "promoter_tss_component",
  "promoter_tss_distance",
  "promoter_anchor_side",
  "candidate_regulatory_anchor_side",
  "gene_symbol",
  "ensembl_gene_id",
  "atac_promoter_anchor",
  "atac_candidate_regulatory_anchor",
  "candidate_regulatory_anchor_has_non_tss_fragment",
  "candidate_regulatory_anchor_non_tss_bp",
  "atac_candidate_regulatory_anchor_non_tss"
)

resource.tables <- list(
  pooled_loop_annotation_resource = df.loop.evidence,
  highlighted_promoter_TSS_nonTSS_ATAC_CTCF_supported_subset =
    df.loop.evidence %>%
      filter(strongest_candidate_support),
  putative_regulatory_promoter_TSS_nonTSS_ATAC_loops =
    df.loop.evidence %>%
      filter(putative_regulatory_support),
  structural_CTCF_supported_loops =
    df.loop.evidence %>%
      filter(
        proposed_major_category ==
          "structural_CTCF_supported_not_ATAC_regulatory"
      ),
  lower_support_or_uncertain_loops =
    df.loop.evidence %>%
      filter(
        proposed_major_category ==
          "lower_support_or_uncertain"
      ),
  all_loop_evidence = df.loop.evidence,
  # Backward-compatible alias for the previous output name.
  main_strongest_candidate_PE_loops =
    df.loop.evidence %>%
      filter(strongest_candidate_support)
)

df.gene.count.by.set <- bind_rows(
  df.loop.evidence %>%
    filter(strongest_candidate_support) %>%
    count(
      gene_set = "main_strongest_candidate_PE_loops",
      gene_symbol,
      ensembl_gene_id,
      sort = TRUE
    ),
  df.loop.evidence %>%
    filter(putative_regulatory_support) %>%
    count(
      gene_set =
        "putative_regulatory_promoter_TSS_nonTSS_ATAC_loops",
      gene_symbol,
      ensembl_gene_id,
      sort = TRUE
    ),
  df.loop.evidence %>%
    filter(current_final_set) %>%
    count(
      gene_set = "current_final_CTCF_promoter_TSS_loops",
      gene_symbol,
      ensembl_gene_id,
      sort = TRUE
    ),
  df.loop.evidence %>%
    filter(
      proposed_major_category ==
        "structural_CTCF_supported_not_ATAC_regulatory"
    ) %>%
    count(
      gene_set = "structural_CTCF_supported_loops",
      gene_symbol,
      ensembl_gene_id,
      sort = TRUE
    )
) %>%
  filter(!is.na(gene_symbol) | !is.na(ensembl_gene_id)) %>%
  arrange(gene_set, desc(n), gene_symbol)

df.gene.count.threshold.summary <- df.gene.count.by.set %>%
  group_by(gene_set) %>%
  summarise(
    n_genes = n(),
    n_ensembl_genes = n_distinct(
      ensembl_gene_id,
      na.rm = TRUE
    ),
    max_interactions_per_gene = max(n),
    n_genes_ge_5 = sum(n >= 5),
    n_genes_ge_10 = sum(n >= 10),
    n_genes_ge_11 = sum(n >= 11),
    n_genes_ge_12 = sum(n >= 12),
    .groups = "drop"
  ) %>%
  arrange(gene_set)

go.gene.sets <- list(
  strongest_ge10 = prepare_go_gene_set(
    df.gene.count.by.set,
    "main_strongest_candidate_PE_loops",
    10
  ),
  strongest_ge11 = prepare_go_gene_set(
    df.gene.count.by.set,
    "main_strongest_candidate_PE_loops",
    11
  ),
  putative_ge10 = prepare_go_gene_set(
    df.gene.count.by.set,
    "putative_regulatory_promoter_TSS_nonTSS_ATAC_loops",
    10
  ),
  putative_ge11 = prepare_go_gene_set(
    df.gene.count.by.set,
    "putative_regulatory_promoter_TSS_nonTSS_ATAC_loops",
    11
  ),
  current_final_ge10 = prepare_go_gene_set(
    df.gene.count.by.set,
    "current_final_CTCF_promoter_TSS_loops",
    10
  ),
  current_final_ge11 = prepare_go_gene_set(
    df.gene.count.by.set,
    "current_final_CTCF_promoter_TSS_loops",
    11
  )
)

go.universe <- df.loop.evidence %>%
  filter(
    passes_promoter_tss_200kb,
    !is.na(ensembl_gene_id)
  ) %>%
  distinct(ensembl_gene_id) %>%
  pull(ensembl_gene_id)

df.go.result <- map_dfr(
  names(go.gene.sets),
  function(set.name) {
    run_go_enrichment(
      go.gene.sets[[set.name]],
      go.universe,
      set.name,
      output.dir
    )
  }
)
