# lintr: disable
if (basename(getwd()) != "enhancer" && dir.exists("enhancer")) {
  setwd("./enhancer")
}
getwd()
funcs.file <- "./funcs_enhancer.R"
source(funcs.file)
list2env(resolve_enhancer_analysis_paths(funcs.file), envir = environment())

library("tidyverse")
library("GenomicRanges")
library("GenomeInfoDb")

options(tibble.width = Inf, tibble.print_max = Inf, tibble.max_extra_cols = Inf, scipen = 999)

################################################################################
# 12. Revised resource tables and downstream gene summaries
################################################################################

# Keep the one-row-per-loop resource separate from the many-row loop-gene
# annotation. A pooled loop is never duplicated merely because multiple genes or
# transcript isoforms overlap an anchor.
df.revised.loop.gene.resource <- df.direct.gene.assignment.position.flags %>%
  left_join(
    df.loop.evidence %>%
      dplyr::select(
        loop_id, resolution, chr1, start1, end1, chr2, start2, end2, loop_distance, n_direct_anchor_sides,
        revised_putative_regulatory_support, revised_single_promoter_distal_atac_support,
        revised_dual_promoter_directional_distal_atac_support, revised_promoter_promoter_compatible,
        revised_single_promoter_without_opposite_atac, revised_dual_promoter_without_directional_distal_atac,
        revised_promoter_associated_without_distal_atac, revised_no_direct_promoter_tss,
        supported_regulatory_direction, revised_major_category, revised_detailed_category, predicted_ctcf_motif_interval_count_anchor1,
        predicted_ctcf_motif_interval_count_anchor2, predicted_ctcf_motif_intervals_both_anchors, predicted_ctcf_motif_annotation_class, candidate_anchor_non_tss_atac_ge50
      ),
    by = c("loop_id", "resolution")
  ) %>%
  mutate(
    ensembl_gene_id = gene_id,
    gene_symbol = gene_name,
    is_direction_supported_putative_gene_assignment =
      revised_putative_regulatory_support &
      (
        n_direct_anchor_sides == 1L |
          supported_regulatory_direction == "both_directions_supported" |
          (supported_regulatory_direction == "anchor1_promoter_to_anchor2_regulatory" & anchor_side == "anchor1") |
          (supported_regulatory_direction == "anchor2_promoter_to_anchor1_regulatory" & anchor_side == "anchor2")
      ),
    revised_gene_assignment_role = case_when(
      is_direction_supported_putative_gene_assignment & n_direct_anchor_sides == 1L ~
        "single_promoter_gene_in_putative_regulatory_loop",
      is_direction_supported_putative_gene_assignment & n_direct_anchor_sides == 2L ~
        "direction_supported_promoter_gene_in_both_anchor_promoter_loop",
      revised_promoter_promoter_compatible & anchor_side == "anchor1" ~ "promoter_TSS_gene_at_anchor1_in_promoter_promoter_loop",
      revised_promoter_promoter_compatible & anchor_side == "anchor2" ~ "promoter_TSS_gene_at_anchor2_in_promoter_promoter_loop",
      revised_single_promoter_without_opposite_atac ~ "direct_promoter_TSS_gene_without_opposite_ATAC_support",
      TRUE ~ "direct_promoter_TSS_gene_in_unassigned_loop"
    ),
    candidate_regulatory_anchor_side = if_else(
      is_direction_supported_putative_gene_assignment,
      opposite_anchor_side,
      NA_character_
    )
  ) %>%
  arrange(loop_id, anchor_side, ensembl_gene_id)

# Verify that direct loop-gene assignments remain complete and unique.
assert_analysis_condition(
  nrow(df.revised.loop.gene.resource) == nrow(df.direct.promoter.tss.gene.assignment) &&
    !anyDuplicated(df.revised.loop.gene.resource$assignment_id),
  "Revised loop-gene resource did not preserve unique direct assignments."
)

# Require every revised putative-regulatory loop to retain a direct gene.
assert_analysis_condition(
  n_distinct(df.revised.loop.gene.resource$loop_id[
    df.revised.loop.gene.resource$is_direction_supported_putative_gene_assignment
  ]) == sum(df.loop.evidence$revised_putative_regulatory_support),
  "Not every revised putative-regulatory loop has a direct gene assignment."
)

# Build a separate secondary inward <=10-kb loop-gene resource if available from 01_2_.
df.revised.secondary.loop.gene.resource <- tibble(
  assignment_id = character(),
  loop_id = character(),
  resolution = character(),
  anchor_side = character(),
  opposite_anchor_side = character(),
  ensembl_gene_id = character(),
  gene_symbol = character(),
  revised_gene_assignment_role = character()
)

if (exists("df.proximal.gene.assignment.position.flags") && exists("df.secondary.inward.proximal.gene.assignment")) {
  df.revised.secondary.loop.gene.resource <- df.proximal.gene.assignment.position.flags %>%
    filter(is_secondary_inward_proximal_candidate_10kb) %>%
    left_join(
      df.loop.evidence %>%
        dplyr::select(
          loop_id,
          resolution,
          chr1,
          start1,
          end1,
          chr2,
          start2,
          end2,
          loop_distance,
          revised_major_category,
          revised_detailed_category,
          predicted_ctcf_motif_interval_count_anchor1,
          predicted_ctcf_motif_interval_count_anchor2,
          predicted_ctcf_motif_annotation_class,
          candidate_anchor_non_tss_atac_ge50
        ),
      by = c("loop_id", "resolution")
    ) %>%
    mutate(
      ensembl_gene_id = gene_id,
      gene_symbol = gene_name,
      revised_gene_assignment_role = "secondary_inward_proximal_10kb_promoter_TSS_candidate"
    ) %>%
    arrange(loop_id, anchor_side, ensembl_gene_id)

  assert_analysis_condition(
    nrow(df.revised.secondary.loop.gene.resource) ==
      nrow(df.secondary.inward.proximal.gene.assignment) &&
      !anyDuplicated(df.revised.secondary.loop.gene.resource$assignment_id),
    "Revised secondary inward loop-gene resource lost unique assignments."
  )
}

df.dual.primary.rep.input <- if (exists("df.dual.primary.representative.gene.assignment")) {
  df.dual.primary.representative.gene.assignment %>%
    transmute(
      gene_set = "revised_dual_promoter_representative_anchor_sensitivity",
      loop_id,
      resolution,
      ensembl_gene_id = gene_id,
      gene_symbol = gene_name
    )
} else {
  tibble(
    gene_set = character(),
    loop_id = character(),
    resolution = character(),
    ensembl_gene_id = character(),
    gene_symbol = character()
  )
}

# Generate explicit loop-gene memberships. Counting is performed after exact
# loop/gene deduplication, so multiple annotation sources and promoter-promoter
# assignments at both anchors cannot inflate a gene's interaction count.
df.revised.gene.loop.membership <- bind_rows(
  df.revised.loop.gene.resource %>%
    transmute(
      gene_set = "revised_direct_promoter_TSS_all",
      loop_id,
      resolution,
      ensembl_gene_id,
      gene_symbol
    ),
  df.revised.loop.gene.resource %>%
    filter(is_direction_supported_putative_gene_assignment) %>%
    transmute(
      gene_set = "revised_putative_regulatory",
      loop_id,
      resolution,
      ensembl_gene_id,
      gene_symbol
    ),
  df.revised.loop.gene.resource %>%
    filter(revised_promoter_associated_without_distal_atac) %>%
    transmute(
      gene_set = "revised_promoter_associated_without_distal_ATAC_support",
      loop_id,
      resolution,
      ensembl_gene_id,
      gene_symbol
    ),
  df.revised.loop.gene.resource %>%
    filter(revised_promoter_promoter_compatible) %>%
    transmute(
      gene_set = "revised_promoter_promoter_compatible",
      loop_id,
      resolution,
      ensembl_gene_id,
      gene_symbol
    ),
  df.revised.secondary.loop.gene.resource %>%
    transmute(
      gene_set = "revised_secondary_inward_proximal_10kb",
      loop_id,
      resolution,
      ensembl_gene_id,
      gene_symbol
    ),
  df.dual.primary.rep.input
) %>%
  filter(!is.na(ensembl_gene_id)) %>%
  distinct(gene_set, loop_id, ensembl_gene_id, .keep_all = TRUE) %>%
  arrange(gene_set, loop_id, ensembl_gene_id)

df.revised.gene.count.by.set <- summarise_revised_gene_loop_counts(df.revised.gene.loop.membership)

# Reuse the legacy approach-2 interface on deduplicated revised loop/gene
# objects. All pooled loops are already <2 Mb; genes with at least two distinct
# loops are exported as ready-to-paste Ensembl-ID inputs for g:Profiler.
df.approach2.loop.gene.input <- df.revised.gene.loop.membership %>%
  left_join(df.loop.distinct.2mb %>% dplyr::select(loop_id, loop_distance), by = "loop_id") %>%
  transmute(gene_set, gene_id = ensembl_gene_id, gene_symbol, loop.id = loop_id, distance = loop_distance) %>%
  distinct(gene_set, gene_id, loop.id, .keep_all = TRUE)

approach2.gene.set.names <- sort(unique(df.approach2.loop.gene.input$gene_set))
approach2.results <- set_names(
  map(approach2.gene.set.names, function(gene.set.i) {
    approach_2nd_analyze_loops_by_threshold(
      df = df.approach2.loop.gene.input %>% filter(gene_set == gene.set.i) %>% dplyr::select(gene_id, loop.id, distance),
      threshold_distance = 2000000L,
      top_n_genes = 80L,
      print_top_n = 50L,
      minimum_interactions = 2L,
      create_plot = FALSE,
      print_results = FALSE
    )
  }),
  approach2.gene.set.names
)

df.approach2.gene.symbol.lookup <- df.approach2.loop.gene.input %>%
  group_by(gene_set, gene_id) %>%
  summarise(
    gene_symbol = {
      symbol.values <- sort(unique(na.omit(gene_symbol)))
      if (length(symbol.values) == 0L) NA_character_ else symbol.values[[1]]
    },
    .groups = "drop"
  )

df.approach2.gene.loop.count <- imap_dfr(approach2.results, function(result.i, gene.set.i) result.i$gene_loop_count %>% mutate(gene_set = gene.set.i, .before = 1)) %>%
  left_join(df.approach2.gene.symbol.lookup, by = c("gene_set", "gene_id")) %>%
  arrange(gene_set, desc(n), gene_symbol, gene_id)

df.approach2.multiple.interaction.genes <- imap_dfr(approach2.results, function(result.i, gene.set.i) result.i$multiple_interaction_genes %>% mutate(gene_set = gene.set.i, .before = 1)) %>%
  left_join(df.approach2.gene.symbol.lookup, by = c("gene_set", "gene_id")) %>%
  arrange(gene_set, desc(n), gene_symbol, gene_id)

df.approach2.multiple.interaction.summary <- imap_dfr(approach2.results, function(result.i, gene.set.i) result.i$summary %>% mutate(gene_set = gene.set.i, .before = 1)) %>%
  arrange(gene_set)

df.revised.gene.count.threshold.summary <- df.revised.gene.count.by.set %>%
  group_by(gene_set) %>%
  summarise(
    n_genes = n(),
    max_interactions_per_gene = max(n),
    n_genes_ge_5 = sum(n >= 5L),
    n_genes_ge_10 = sum(n >= 10L),
    n_genes_ge_11 = sum(n >= 11L),
    n_genes_ge_12 = sum(n >= 12L),
    .groups = "drop"
  ) %>%
  arrange(gene_set)

df.revised.gene.count.threshold.detail <- tidyr::crossing(
  gene_set = sort(unique(df.revised.gene.count.by.set$gene_set)),
  minimum_interactions = c(1L, 5L, 10L, 11L, 12L)
) %>%
  left_join(
    df.revised.gene.count.by.set %>%
      tidyr::crossing(minimum_interactions = c(1L, 5L, 10L, 11L, 12L)) %>%
      filter(n >= minimum_interactions) %>%
      count(gene_set, minimum_interactions, name = "n_genes"),
    by = c("gene_set", "minimum_interactions")
  ) %>%
  mutate(n_genes = coalesce(n_genes, 0L)) %>%
  arrange(gene_set, minimum_interactions)

revised.go.gene.sets <- list(
  revised_putative_all = prepare_go_gene_set(
    df.revised.gene.count.by.set,
    "revised_putative_regulatory",
    1L
  )
)

assert_analysis_condition(
  nrow(revised.go.gene.sets$revised_putative_all) > 0L &&
    !anyDuplicated(revised.go.gene.sets$revised_putative_all$ensembl_gene_id),
  "The primary exact-call GO input is empty or contains duplicate Ensembl gene IDs."
)

# The primary GO background contains every Ensembl gene observable by the
# all-transcript plus EPD direct promoter/TSS assignment process.
revised.go.universe <- df.revised.loop.gene.resource %>%
  filter(!is.na(ensembl_gene_id)) %>%
  distinct(ensembl_gene_id) %>%
  pull(ensembl_gene_id)

df.revised.downstream.analysis.definition <- tribble(
  ~analysis_item, ~definition,
  "pooled_loop_resource", paste0(
    "One row per each of 31,021 pooled loops; gene multiplicity does not ",
    "duplicate loop records."
  ),
  "loop_gene_resource", paste0(
    "One row per direct loop-anchor-Ensembl-gene assignment, preserving all ",
    "overlapping genes without nearest-gene selection."
  ),
  "secondary_inward_loop_gene_resource", paste0(
    "Separate loop-anchor-gene candidates within 10 kb toward the loop ",
    "interior, excluding the same loop-anchor-gene from the primary direct ",
    "tier; this resource is not merged into the primary GO input."
  ),
  "interaction_count", paste0(
    "Number of distinct exact HiCCUPS loop calls per stable Ensembl gene ",
    "within each revised evidence set; duplicate annotation sources and ",
    "anchors are collapsed."
  ),
  "approximate_loop_locus_sensitivity", paste0(
    "A sensitivity-only count in which qualifying overlapping ",
    "cross-resolution loop calls are grouped into approximate loci; the ",
    "31,021-call pooled resource remains unchanged."
  ),
  "canonical_TSS_sensitivity", paste0(
    "A separate exploratory GO sensitivity uses only direct Ensembl-canonical ",
    "transcript TSS assignments; it does not replace the all-transcript plus ",
    "EPD exact-call primary analysis."
  ),
  "GO_universe", paste0(
    "The primary analysis uses all Ensembl genes observable by the revised ",
    "direct assignment process; the canonical-TSS sensitivity uses the ",
    "corresponding canonical-TSS-observable background."
  ),
  "GO_interpretation", paste0(
    "GO enrichment is an exploratory functional interpretation of gene-set ",
    "overrepresentation; it does not validate individual loops, target genes, ",
    "enhancers, or regulatory direction."
  ),
  "top_gene_thresholds", paste0(
    ">=10 and >=11 interaction analyses are sensitivity summaries retained for ",
    "comparability; no single threshold is treated as biological validation."
  )
)

revised.resource.tables <- list(
  revised_pooled_loop_annotation_resource = df.loop.evidence,
  revised_direct_loop_gene_assignments = df.revised.loop.gene.resource,
  revised_secondary_inward_10kb_loop_gene_assignments = df.revised.secondary.loop.gene.resource,
  revised_secondary_inward_10kb_loops = if (exists("df.secondary.inward.proximal.loop.summary")) {
    df.secondary.inward.proximal.loop.summary %>% filter(has_any_secondary_inward_proximal_10kb)
  } else {
    tibble()
  },
  revised_putative_regulatory_loops = df.loop.evidence %>% filter(revised_putative_regulatory_support),
  revised_promoter_associated_without_distal_atac_loops =
    df.loop.evidence %>% filter(revised_promoter_associated_without_distal_atac),
  revised_promoter_promoter_compatible_loops = df.loop.evidence %>% filter(revised_promoter_promoter_compatible),
  revised_single_promoter_without_opposite_atac_loops = df.loop.evidence %>% filter(revised_single_promoter_without_opposite_atac),
  revised_no_direct_promoter_tss_loops = df.loop.evidence %>% filter(revised_no_direct_promoter_tss)
)

################################################################################
# 12-1. HiCCUPS provenance and gene-ranking sensitivity analyses
#
# Exact HiCCUPS calls remain the main pooled resource. Nearby calls at different
# resolutions are grouped only for a conservative ranking-sensitivity analysis.
# Gene ranks are also recalculated using Ensembl-canonical TSS assignments only.
################################################################################

# Preserve every source-library HiCCUPS quality field for provenance and QC.
df.hiccups.sample.loop.quality <- df.sample.loop.1based %>%
  dplyr::select(sample, strain, source_file, sample_loop_id, loop_id, resolution, resolution_bp, chr1, start1, end1, chr2, start2, end2, loop_distance, passes_lt2mb, starts_with("hiccups_"))

df.hiccups.quality.field.definition <- tribble(
  ~field_family, ~resubmission_use,
  "observed", paste0(
    "Retained source contact count for provenance and descriptive QC; not a ",
    "promoter-enhancer selection filter."
  ),
  "expectedBL_expectedDonut_expectedH_expectedV", paste0(
    "Retained HiCCUPS local-background estimates for provenance and ",
    "descriptive QC; not recombined across libraries as inferential evidence."
  ),
  "fdrBL_fdrDonut_fdrH_fdrV", paste0(
    "Retained caller-reported source FDR fields; summarized descriptively ",
    "across supporting libraries and not used as a new P-E filter."
  ),
  "numCollapsed", paste0(
    "Retained number of collapsed HiCCUPS pixels/calls as source provenance; ",
    "not used as a P-E filter."
  ),
  "centroid1_centroid2_radius", paste0(
    "Retained source geometry; pooled centroid medians support the ",
    "all-resolution canonical-loop-locus sensitivity analysis."
  )
)

# Construct an all-resolution canonical-locus sensitivity without changing the
# exact pooled HiCCUPS call IDs or their source-library provenance.
approximate.loop.locus.analysis <- build_approximate_loop_loci(df.loop.distinct.2mb)
df.approximate.loop.locus.edge <- approximate.loop.locus.analysis$edge
df.approximate.loop.locus.map <- approximate.loop.locus.analysis$map
df.approximate.loop.locus.summary <- approximate.loop.locus.analysis$summary
df.approximate.loop.locus.method <- approximate.loop.locus.analysis$method

# Repeat canonicalization with half-sized distance tolerances to quantify how
# strongly the locus count depends on the HiCCUPS-derived midpoint thresholds.
strict.approximate.loop.locus.analysis <- build_approximate_loop_loci(
  df.loop.distinct.2mb,
  merge.distance.bp = c("5K" = 10000L, "10K" = 10000L, "25K" = 25000L)
)
df.approximate.loop.locus.threshold.sensitivity <- bind_rows(
  df.approximate.loop.locus.method %>% mutate(tolerance_set = "HiCCUPS_default_20kb_20kb_50kb", .before = 1),
  strict.approximate.loop.locus.analysis$method %>% mutate(tolerance_set = "half_distance_10kb_10kb_25kb", .before = 1)
)
rm(strict.approximate.loop.locus.analysis)

# Select the primary main-set loop-gene membership based on all retained TSS and
# EPD evidence, then derive a stricter Ensembl-canonical-TSS-only sensitivity.
main.putative.loop.ids <- df.loop.evidence %>%
  filter(revised_putative_regulatory_support) %>%
  pull(loop_id)

df.main.all.annotation.gene.loop.membership <- df.revised.gene.loop.membership %>%
  filter(gene_set == "revised_putative_regulatory") %>%
  dplyr::select(loop_id, resolution, ensembl_gene_id, gene_symbol) %>%
  distinct(loop_id, ensembl_gene_id, .keep_all = TRUE)

df.main.putative.promoter.anchor <- df.revised.loop.gene.resource %>%
  filter(is_direction_supported_putative_gene_assignment) %>%
  distinct(loop_id, anchor_side)

df.canonical.tss.main.loop.gene.membership <- df.direct.true.tss.anchor.overlap %>%
  filter(coalesce(is_ensembl_canonical, FALSE)) %>%
  inner_join(df.main.putative.promoter.anchor, by = c("loop_id", "anchor_side")) %>%
  transmute(loop_id, resolution, ensembl_gene_id = gene_id, gene_symbol = gene_name, transcript_id, transcript_id_versioned, true_tss_id = str_remove(annotation_id, ":TSS_pm1kb$")) %>%
  distinct(loop_id, ensembl_gene_id, .keep_all = TRUE) %>%
  arrange(loop_id, ensembl_gene_id)

# Canonical-TSS sensitivity records must be a subset of the all-annotation main
# loop-gene membership, never a new source of loop or gene assignments.
assert_analysis_condition(
  nrow(anti_join(df.canonical.tss.main.loop.gene.membership, df.main.all.annotation.gene.loop.membership, by = c("loop_id", "ensembl_gene_id"))) == 0L,
  "Canonical-TSS sensitivity produced a loop-gene assignment outside the main set."
)

# Build the canonical-TSS-only gene set and its matching observable background.
revised.go.gene.sets$revised_putative_canonical_TSS_only <- df.canonical.tss.main.loop.gene.membership %>%
  group_by(ensembl_gene_id) %>%
  summarise(
    gene_symbol = {
      symbol.values <- sort(unique(na.omit(gene_symbol)))
      if (length(symbol.values) == 0L) NA_character_ else symbol.values[1]
    },
    n = n_distinct(loop_id),
    .groups = "drop"
  ) %>%
  arrange(desc(n), gene_symbol, ensembl_gene_id)

revised.go.canonical.universe <- df.direct.true.tss.anchor.overlap %>%
  filter(coalesce(is_ensembl_canonical, FALSE), !is.na(gene_id)) %>%
  distinct(gene_id) %>%
  pull(gene_id)

revised.go.universes <- list(
  revised_putative_all = revised.go.universe,
  revised_putative_canonical_TSS_only = revised.go.canonical.universe
)

# Run primary and canonical-TSS sensitivity GO analyses with matched universes.
df.revised.go.input.summary <- imap_dfr(
  revised.go.gene.sets,
  function(df.gene, set.name) {
    universe.i <- revised.go.universes[[set.name]]
    tibble(
      gene_set = set.name,
      analysis_role = if_else(
        set.name == "revised_putative_all",
        "primary_exact_call_exploratory_GO",
        "canonical_TSS_only_exploratory_GO_sensitivity"
      ),
      n_input_ensembl_genes = n_distinct(df.gene$ensembl_gene_id),
      n_universe_ensembl_genes = length(unique(universe.i))
    )
  }
)

df.revised.go.result <- map_dfr(
  names(revised.go.gene.sets),
  function(set.name) {
    run_go_enrichment(
      revised.go.gene.sets[[set.name]],
      revised.go.universes[[set.name]],
      set.name,
      output.dir
    )
  }
)

df.revised.go.significance.summary <- df.revised.go.input.summary %>%
  left_join(
    df.revised.go.result %>%
      group_by(gene_set) %>%
      summarise(n_GO_BP_terms_tested = n(), n_GO_BP_terms_FDR_lt_0_05 = sum(p.adjust < 0.05, na.rm = TRUE), minimum_adjusted_p = min(p.adjust, na.rm = TRUE), .groups = "drop"),
    by = "gene_set"
  ) %>%
  mutate(
    n_GO_BP_terms_tested = coalesce(n_GO_BP_terms_tested, 0L),
    n_GO_BP_terms_FDR_lt_0_05 = coalesce(n_GO_BP_terms_FDR_lt_0_05, 0L)
  ) %>%
  arrange(gene_set)

df.main.all.annotation.gene.loop.locus.count <- summarise_gene_loop_locus_counts(df.main.all.annotation.gene.loop.membership, df.approximate.loop.locus.map, evidence.definition = "all_transcript_TSS_plus_EPD_primary_assignment")
df.main.canonical.tss.gene.loop.locus.count <- summarise_gene_loop_locus_counts(df.canonical.tss.main.loop.gene.membership, df.approximate.loop.locus.map, evidence.definition = "Ensembl_canonical_transcript_TSS_only")

# Place all four count metrics on one row per gene for rank-correlation
# sensitivity without selecting a fixed number of top genes.
df.gene.rank.sensitivity.detail <- df.main.all.annotation.gene.loop.locus.count %>%
  transmute(ensembl_gene_id, all_annotation_gene_symbol = gene_symbol, all_annotation_n_exact_loop_calls = n_exact_loop_calls, all_annotation_n_approximate_loop_loci = n_approximate_loop_loci) %>%
  full_join(
    df.main.canonical.tss.gene.loop.locus.count %>%
      transmute(ensembl_gene_id, canonical_tss_gene_symbol = gene_symbol, canonical_tss_n_exact_loop_calls = n_exact_loop_calls, canonical_tss_n_approximate_loop_loci = n_approximate_loop_loci),
    by = "ensembl_gene_id"
  ) %>%
  mutate(
    gene_symbol = coalesce(all_annotation_gene_symbol, canonical_tss_gene_symbol),
    across(ends_with("loop_calls"), ~ replace_na(.x, 0L)),
    across(ends_with("loop_loci"), ~ replace_na(.x, 0L)),
    rank_all_annotation_exact_calls = min_rank(dplyr::desc(all_annotation_n_exact_loop_calls)),
    rank_all_annotation_approximate_loci = min_rank(dplyr::desc(all_annotation_n_approximate_loop_loci)),
    rank_canonical_tss_exact_calls = min_rank(dplyr::desc(canonical_tss_n_exact_loop_calls)),
    rank_canonical_tss_approximate_loci = min_rank(dplyr::desc(canonical_tss_n_approximate_loop_loci))
  ) %>%
  dplyr::select(ensembl_gene_id, gene_symbol, all_annotation_n_exact_loop_calls, all_annotation_n_approximate_loop_loci, canonical_tss_n_exact_loop_calls, canonical_tss_n_approximate_loop_loci, starts_with("rank_")) %>%
  arrange(rank_all_annotation_exact_calls, ensembl_gene_id)

df.gene.rank.sensitivity.correlation <- bind_rows(
  summarise_gene_count_spearman(df.gene.rank.sensitivity.detail, "all_annotation_n_exact_loop_calls", "all_annotation_n_approximate_loop_loci", "all_annotation_exact_calls_vs_approximate_loci"),
  summarise_gene_count_spearman(df.gene.rank.sensitivity.detail, "all_annotation_n_exact_loop_calls", "canonical_tss_n_exact_loop_calls", "all_annotation_vs_canonical_TSS_exact_calls"),
  summarise_gene_count_spearman(df.gene.rank.sensitivity.detail, "all_annotation_n_approximate_loop_loci", "canonical_tss_n_approximate_loop_loci", "all_annotation_vs_canonical_TSS_approximate_loci")
)

#
# End of 02_promoter_enhancer_interaction_resubmit_gene_resources.R
