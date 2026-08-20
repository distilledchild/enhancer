# lintr: disable
setwd("./enhancer")
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
    df.revised.loop.evidence %>%
      dplyr::select(
        loop_id, resolution, chr1, start1, end1, chr2, start2, end2, loop_distance, n_direct_anchor_sides,
        revised_putative_regulatory_support, revised_promoter_promoter_compatible, revised_single_promoter_without_opposite_atac,
        revised_no_direct_promoter_tss, revised_major_category, revised_detailed_category, predicted_ctcf_motif_interval_count_anchor1,
        predicted_ctcf_motif_interval_count_anchor2, predicted_ctcf_motif_intervals_both_anchors, predicted_ctcf_motif_annotation_class, candidate_anchor_non_tss_atac_ge50
      ),
    by = c("loop_id", "resolution")
  ) %>%
  mutate(
    ensembl_gene_id = gene_id,
    gene_symbol = gene_name,
    revised_gene_assignment_role = case_when(
      revised_putative_regulatory_support ~ "direct_promoter_TSS_gene_in_putative_regulatory_loop",
      revised_promoter_promoter_compatible & anchor_side == "anchor1" ~ "promoter_TSS_gene_at_anchor1_in_promoter_promoter_loop",
      revised_promoter_promoter_compatible & anchor_side == "anchor2" ~ "promoter_TSS_gene_at_anchor2_in_promoter_promoter_loop",
      revised_single_promoter_without_opposite_atac ~ "direct_promoter_TSS_gene_without_opposite_ATAC_support",
      TRUE ~ "direct_promoter_TSS_gene_in_unassigned_loop"
    ),
    candidate_regulatory_anchor_side = if_else(revised_putative_regulatory_support, opposite_anchor_side, NA_character_)
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
  n_distinct(df.revised.loop.gene.resource$loop_id[df.revised.loop.gene.resource$revised_putative_regulatory_support]) == sum(df.revised.loop.evidence$revised_putative_regulatory_support),
  "Not every revised putative-regulatory loop has a direct gene assignment."
)

# Build a separate secondary inward <=10-kb loop-gene resource. It is exported
# and countable but is not merged into the primary direct GO input.
df.revised.secondary.loop.gene.resource <- df.proximal.gene.assignment.position.flags %>%
  filter(is_secondary_inward_proximal_candidate_10kb) %>%
  left_join(
    df.revised.loop.evidence %>%
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

# Verify that secondary inward loop-gene assignments remain complete and unique.
assert_analysis_condition(
  nrow(df.revised.secondary.loop.gene.resource) ==
    nrow(df.secondary.inward.proximal.gene.assignment) &&
    !anyDuplicated(df.revised.secondary.loop.gene.resource$assignment_id),
  "Revised secondary inward loop-gene resource lost unique assignments."
)

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
    filter(revised_putative_regulatory_support) %>%
    transmute(
      gene_set = "revised_putative_regulatory",
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
  df.dual.primary.representative.gene.assignment %>%
    transmute(
      gene_set = paste0(
        "revised_dual_promoter_representative_anchor_",
        "sensitivity"
      ),
      loop_id,
      resolution,
      ensembl_gene_id = gene_id,
      gene_symbol = gene_name
    )
) %>%
  filter(!is.na(ensembl_gene_id)) %>%
  distinct(gene_set, loop_id, ensembl_gene_id, .keep_all = TRUE) %>%
  arrange(gene_set, loop_id, ensembl_gene_id)

df.revised.gene.count.by.set <- summarise_revised_gene_loop_counts(df.revised.gene.loop.membership)

# Reuse the legacy approach-2 interface on deduplicated revised loop/gene
# objects. All pooled loops are already <2 Mb; genes with at least two distinct
# loops are exported as ready-to-paste Ensembl-ID inputs for g:Profiler.
df.approach2.loop.gene.input <- df.revised.gene.loop.membership %>%
  left_join(df.loop.universe %>% dplyr::select(loop_id, loop_distance), by = "loop_id") %>%
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

df.revised.gene.count.threshold.detail <- crossing(
  gene_set = sort(unique(df.revised.gene.count.by.set$gene_set)),
  minimum_interactions = c(1L, 5L, 10L, 11L, 12L)
) %>%
  left_join(
    df.revised.gene.count.by.set %>%
      crossing(minimum_interactions = c(1L, 5L, 10L, 11L, 12L)) %>%
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

assert_analysis_row_count(
  revised.go.gene.sets$revised_putative_all,
  6420L,
  "The primary exact-call GO input no longer contains 6,420 genes."
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
  revised_pooled_loop_annotation_resource = df.revised.loop.evidence,
  revised_direct_loop_gene_assignments = df.revised.loop.gene.resource,
  revised_secondary_inward_10kb_loop_gene_assignments = df.revised.secondary.loop.gene.resource,
  revised_secondary_inward_10kb_loops = df.secondary.inward.proximal.loop.summary %>% filter(has_any_secondary_inward_proximal_10kb),
  revised_putative_regulatory_loops = df.revised.loop.evidence %>% filter(revised_putative_regulatory_support),
  revised_promoter_promoter_compatible_loops = df.revised.loop.evidence %>% filter(revised_promoter_promoter_compatible),
  revised_single_promoter_without_opposite_atac_loops = df.revised.loop.evidence %>% filter(revised_single_promoter_without_opposite_atac),
  revised_no_direct_promoter_tss_loops = df.revised.loop.evidence %>% filter(revised_no_direct_promoter_tss)
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
approximate.loop.locus.analysis <- build_approximate_loop_loci(df.loop.universe)
df.approximate.loop.locus.edge <- approximate.loop.locus.analysis$edge
df.approximate.loop.locus.map <- approximate.loop.locus.analysis$map
df.approximate.loop.locus.summary <- approximate.loop.locus.analysis$summary
df.approximate.loop.locus.method <- approximate.loop.locus.analysis$method

# Repeat canonicalization with half-sized distance tolerances to quantify how
# strongly the locus count depends on the HiCCUPS-derived midpoint thresholds.
strict.approximate.loop.locus.analysis <- build_approximate_loop_loci(
  df.loop.universe,
  merge.distance.bp = c("5K" = 10000L, "10K" = 10000L, "25K" = 25000L)
)
df.approximate.loop.locus.threshold.sensitivity <- bind_rows(
  df.approximate.loop.locus.method %>% mutate(tolerance_set = "HiCCUPS_default_20kb_20kb_50kb", .before = 1),
  strict.approximate.loop.locus.analysis$method %>% mutate(tolerance_set = "half_distance_10kb_10kb_25kb", .before = 1)
)
rm(strict.approximate.loop.locus.analysis)

# Select the primary main-set loop-gene membership based on all retained TSS and
# EPD evidence, then derive a stricter Ensembl-canonical-TSS-only sensitivity.
main.putative.loop.ids <- df.revised.loop.evidence %>%
  filter(revised_putative_regulatory_support) %>%
  pull(loop_id)

df.main.all.annotation.gene.loop.membership <- df.revised.gene.loop.membership %>%
  filter(gene_set == "revised_putative_regulatory") %>%
  dplyr::select(loop_id, resolution, ensembl_gene_id, gene_symbol) %>%
  distinct(loop_id, ensembl_gene_id, .keep_all = TRUE)

df.main.putative.promoter.anchor <- df.revised.atac.unambiguous.orientation %>%
  filter(loop_id %in% main.putative.loop.ids) %>%
  transmute(loop_id, anchor_side = promoter_anchor_side)

df.canonical.tss.main.loop.gene.membership <- df.primary.direct.true.tss.anchor.overlap %>%
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

revised.go.canonical.universe <- df.primary.direct.true.tss.anchor.overlap %>%
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

################################################################################
# 12-2. Revised Figure 5a-c positional-density analysis
#
# This rebuild follows the original Figure 5 transformation in
# enhancer_promoter_interaction.R and
# enhancer_promoter_interaction_for_submission.R: one loop length is added on
# both sides, and the two anchor midpoints map to relative positions 0 and 1.
# Panel a reports predicted CTCF motif intervals as structural sequence
# annotation. Panel b uses strand-aware Ensembl transcript TSS coordinates
# rather than start-codon intervals. Panel c uses the coordinate-normalized EPD
# rn7 promoter annotation reconstructed in Step 1. Identical genomic sites are
# counted once so transcript or annotation duplication does not change density
# weighting.
################################################################################

figure5.resolution.colors <- c(
  "5K" = "#a6cee3",
  "10K" = "#1f78b4",
  "25K" = "#1f3a93"
)

# Define the expanded loop window and retain the unrounded anchor midpoints for
# an exact relative-position transformation: anchor1 = 0 and anchor2 = 1.
df.figure5.loop.window <- df.loop.universe %>%
  transmute(
    loop_id,
    chr = chr1,
    resolution = factor(resolution, levels = names(figure5.resolution.colors)),
    anchor1_midpoint = (start1 + end1) / 2,
    anchor2_midpoint = (start2 + end2) / 2,
    anchor_midpoint_distance = anchor2_midpoint - anchor1_midpoint,
    expanded_start_unclipped = anchor1_midpoint - anchor_midpoint_distance,
    expanded_end_unclipped = anchor2_midpoint + anchor_midpoint_distance,
    expanded_start = pmax(1L, as.integer(floor(expanded_start_unclipped))),
    expanded_end = as.integer(ceiling(expanded_end_unclipped)),
    same_chromosome = chr1 == chr2
  )

# Require valid cis-loop windows before calculating positional distributions.
assert_analysis_condition(
  all(df.figure5.loop.window$same_chromosome) &&
    all(df.figure5.loop.window$anchor_midpoint_distance > 0) &&
    !any(is.na(df.figure5.loop.window$resolution)),
  "Figure 5 loop windows require ordered cis loops at 5K, 10K, or 25K."
)

gr.figure5.loop.window <- GRanges(
  seqnames = df.figure5.loop.window$chr,
  ranges = IRanges(start = df.figure5.loop.window$expanded_start, end = df.figure5.loop.window$expanded_end)
)

# Map predicted CTCF motif intervals to expanded loop windows chromosome by
# chromosome. Relative positions are aggregated into narrow bins for plotting,
# avoiding materialization of tens of millions of overlap rows in memory.
figure5.ctcf.chromosomes <- intersect(
  unique(as.character(seqnames(gr.ctcf.motif))),
  unique(df.figure5.loop.window$chr)
)
figure5.ctcf.chromosome.results <- map(
  figure5.ctcf.chromosomes,
  function(chr.i) {
    motif.index <- which(as.character(seqnames(gr.ctcf.motif)) == chr.i)
    loop.index <- which(df.figure5.loop.window$chr == chr.i)
    hit.i <- findOverlaps(
      gr.ctcf.motif[motif.index],
      gr.figure5.loop.window[loop.index],
      type = "any",
      select = "all"
    )
    if (length(hit.i) == 0L) {
      return(list(density = tibble(), summary = tibble()))
    }

    motif.index.hit <- motif.index[queryHits(hit.i)]
    loop.index.hit <- loop.index[subjectHits(hit.i)]
    relative.position <- (
      (
        start(gr.ctcf.motif)[motif.index.hit] +
          end(gr.ctcf.motif)[motif.index.hit]
      ) / 2 - df.figure5.loop.window$anchor1_midpoint[loop.index.hit]
    ) / df.figure5.loop.window$anchor_midpoint_distance[loop.index.hit]
    keep <- dplyr::between(relative.position, -1, 2)

    density.i <- tibble(
      resolution = df.figure5.loop.window$resolution[loop.index.hit[keep]],
      relative_position = relative.position[keep]
    ) %>%
      mutate(relative_position = round(relative_position / 0.0025) * 0.0025) %>%
      count(resolution, relative_position, name = "n_overlap")

    summary.i <- tibble(
      resolution = df.figure5.loop.window$resolution[loop.index.hit[keep]],
      feature_index = motif.index.hit[keep],
      loop_id = df.figure5.loop.window$loop_id[loop.index.hit[keep]]
    ) %>%
      group_by(resolution) %>%
      summarise(
        feature_type = "predicted_CTCF_motif_interval",
        n_loop_feature_overlaps = n(),
        n_unique_features = n_distinct(feature_index),
        n_unique_loops = n_distinct(loop_id),
        .groups = "drop"
      )
    list(density = density.i, summary = summary.i)
  }
)
df.figure5.ctcf.relative.position <- map_dfr(
  figure5.ctcf.chromosome.results,
  "density"
) %>%
  group_by(resolution, relative_position) %>%
  summarise(n_overlap = sum(n_overlap), .groups = "drop")
df.figure5.ctcf.feature.summary <- map_dfr(
  figure5.ctcf.chromosome.results,
  "summary"
) %>%
  group_by(feature_type, resolution) %>%
  summarise(
    n_loop_feature_overlaps = sum(n_loop_feature_overlaps),
    n_unique_features = sum(n_unique_features),
    n_unique_loops = sum(n_unique_loops),
    .groups = "drop"
  )

# Collapse transcripts sharing the same strand-aware TSS into one genomic TSS
# site while preserving transcript and gene multiplicity as descriptive fields.
df.figure5.true.tss.site <- df.true.tss.transcript %>%
  group_by(chr, true_tss_start, true_tss_end, strand) %>%
  summarise(n_transcripts = n_distinct(transcript_id), n_genes = n_distinct(gene_id), .groups = "drop") %>%
  mutate(feature_id = str_c(chr, true_tss_start, strand, sep = ":"), feature_position = as.numeric(true_tss_start), feature_type = "strand_aware_true_TSS")

gr.figure5.true.tss.site <- GRanges(
  seqnames = df.figure5.true.tss.site$chr,
  ranges = IRanges(start = df.figure5.true.tss.site$true_tss_start, end = df.figure5.true.tss.site$true_tss_end)
)

# Map every unique true-TSS site within the expanded loop windows and calculate
# its location in the original Figure 5 coordinate frame (-1 to 2).
figure5.true.tss.hit <- findOverlaps(gr.figure5.true.tss.site, gr.figure5.loop.window, type = "any", select = "all")
df.figure5.true.tss.relative.position <- tibble(
  loop_index = subjectHits(figure5.true.tss.hit),
  feature_index = queryHits(figure5.true.tss.hit)
) %>%
  transmute(
    loop_id = df.figure5.loop.window$loop_id[loop_index],
    resolution = df.figure5.loop.window$resolution[loop_index],
    feature_id = df.figure5.true.tss.site$feature_id[feature_index],
    feature_chr = df.figure5.true.tss.site$chr[feature_index],
    feature_position = df.figure5.true.tss.site$feature_position[feature_index],
    feature_strand = df.figure5.true.tss.site$strand[feature_index],
    anchor1_midpoint = df.figure5.loop.window$anchor1_midpoint[loop_index],
    anchor2_midpoint = df.figure5.loop.window$anchor2_midpoint[loop_index],
    anchor_midpoint_distance = df.figure5.loop.window$anchor_midpoint_distance[loop_index],
    relative_position = (feature_position - anchor1_midpoint) / anchor_midpoint_distance,
    feature_type = "strand_aware_true_TSS"
  ) %>%
  filter(dplyr::between(relative_position, -1, 2))

# Retain each coordinate-normalized EPD interval once and use its interval
# midpoint, matching the promoter-position definition in the original analysis.
df.figure5.promoter.site <- df.promoter.epd.rn7.1based %>%
  distinct(chr, promoter_start, promoter_end, strand, .keep_all = TRUE) %>%
  transmute(
    feature_id = promoter_annotation_id, chr, promoter_start, promoter_end, strand, gene_id, gene_name,
    feature_position = (promoter_start + promoter_end) / 2,
    feature_type = "coordinate_normalized_EPD_promoter"
  )

gr.figure5.promoter.site <- GRanges(
  seqnames = df.figure5.promoter.site$chr,
  ranges = IRanges(start = df.figure5.promoter.site$promoter_start, end = df.figure5.promoter.site$promoter_end)
)

# Map every normalized EPD promoter within the same expanded loop windows and
# apply the identical anchor-centred relative-position transformation.
figure5.promoter.hit <- findOverlaps(gr.figure5.promoter.site, gr.figure5.loop.window, type = "any", select = "all")
df.figure5.promoter.relative.position <- tibble(
  loop_index = subjectHits(figure5.promoter.hit),
  feature_index = queryHits(figure5.promoter.hit)
) %>%
  transmute(
    loop_id = df.figure5.loop.window$loop_id[loop_index],
    resolution = df.figure5.loop.window$resolution[loop_index],
    feature_id = df.figure5.promoter.site$feature_id[feature_index],
    feature_chr = df.figure5.promoter.site$chr[feature_index],
    feature_position = df.figure5.promoter.site$feature_position[feature_index],
    feature_strand = df.figure5.promoter.site$strand[feature_index],
    anchor1_midpoint = df.figure5.loop.window$anchor1_midpoint[loop_index],
    anchor2_midpoint = df.figure5.loop.window$anchor2_midpoint[loop_index],
    anchor_midpoint_distance = df.figure5.loop.window$anchor_midpoint_distance[loop_index],
    relative_position = (feature_position - anchor1_midpoint) / anchor_midpoint_distance,
    feature_type = "coordinate_normalized_EPD_promoter"
  ) %>%
  filter(dplyr::between(relative_position, -1, 2))

# Summarize the number of unique features and loop-feature overlaps contributing
# to each resolution-specific density curve.
df.figure5.revised.feature.summary <- bind_rows(
  df.figure5.ctcf.feature.summary,
  bind_rows(df.figure5.true.tss.relative.position, df.figure5.promoter.relative.position) %>%
    group_by(feature_type, resolution) %>%
    summarise(
      n_loop_feature_overlaps = n(),
      n_unique_features = n_distinct(feature_id),
      n_unique_loops = n_distinct(loop_id),
      .groups = "drop"
    )
)

# Reuse the original resolution colours and density geometry for both revised
# panels while showing the two loop anchors explicitly at x = 0 and x = 1.
create.figure5.revised.density.plot <- function(
  df.relative.position,
  panel.tag,
  panel.title
) {
  if (!"n_overlap" %in% colnames(df.relative.position)) {
    df.relative.position <- df.relative.position %>%
      mutate(n_overlap = 1)
  }
  ggplot(
    df.relative.position,
    aes(
      x = relative_position,
      weight = n_overlap,
      color = resolution,
      fill = resolution
    )
  ) +
    geom_density(alpha = 0.3, linewidth = 0.55) +
    geom_vline(
      xintercept = c(0, 1),
      color = "grey70",
      linewidth = 0.3
    ) +
    scale_color_manual(
      values = figure5.resolution.colors,
      drop = FALSE
    ) +
    scale_fill_manual(
      values = figure5.resolution.colors,
      drop = FALSE
    ) +
    coord_cartesian(xlim = c(-1, 2), ylim = c(0, 0.8)) +
    scale_x_continuous(breaks = c(-1, 0, 1, 2)) +
    labs(
      tag = panel.tag,
      title = panel.title,
      x = "Relative Position to Loop",
      y = "Density",
      color = "Resolution",
      fill = "Resolution"
    ) +
    theme_bw(base_size = 9) +
    theme(
      plot.tag = element_text(face = "bold"),
      plot.tag.position = c(0.02, 0.98),
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      panel.grid.minor = element_blank()
    )
}

plot.figure5a.ctcf.density <- create.figure5.revised.density.plot(
  df.figure5.ctcf.relative.position,
  panel.tag = "a",
  panel.title = "Predicted CTCF motif intervals"
)
plot.figure5b.true.tss.density <- create.figure5.revised.density.plot(
  df.figure5.true.tss.relative.position,
  panel.tag = "b",
  panel.title = "Strand-aware Ensembl TSSs"
)
plot.figure5c.promoter.density <- create.figure5.revised.density.plot(
  df.figure5.promoter.relative.position,
  panel.tag = "c",
  panel.title = "EPD promoters"
)
plot.figure5abc.revised.density <- patchwork::wrap_plots(
  plot.figure5a.ctcf.density,
  plot.figure5b.true.tss.density,
  plot.figure5c.promoter.density,
  nrow = 1,
  guides = "collect"
) &
  theme(legend.position = "bottom")

# Save the final a-c panel in vector PDF and 300-dpi PNG formats.
figure5.output.files <- basename(saving_plot_dual(
  plot.figure5abc.revised.density,
  filename_base = "figure5abc_revised_density_by_resolution",
  output_dir = output.dir,
  width_in = 10.5,
  height_in = 3.3
))

################################################################################
# 12-3. Legacy comparison is intentionally excluded from production
#
# Historical old-vs-new analyses remain in the separate comparison script and
# are not sourced by this production workflow.
################################################################################

