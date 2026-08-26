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
# SENSITIVITY ANALYSIS 01_2: Proximal Promoter/TSS Allocation Based on Distance
#
# Description: Evaluates non-overlapping candidate promoter/TSS annotations
#              across symmetric search windows from 1 bp to 200 kb.
#              Demonstrates the loss of specificity as proximity distance expands.
# Input:       df.loop.distinct.2mb, df.true.tss.transcript, gr.true.tss, gr.epd.promoter
# Output:      df.proximal.promoter.tss.anchor.pair,
#              df.secondary.inward.proximal.gene.assignment,
#              df.proximal.promoter.tss.cumulative.threshold.summary,
#              df.proximal.promoter.tss.summary
################################################################################

# 1. Load coordinate-normalized cache
coord.cache.object.names <- coordinate_cache_object_names()
load_coordinate_cache_objects(cache.dir = coord.cache.dir, object.names = coord.cache.object.names, envir = environment())

# 2. Prepare strand-aware TSS and EPD ranges
df.true.tss.transcript <- build_true_tss_annotation(df.transcript.ensembl.rn7.1based)
gr.true.tss <- create_true_tss_granges(df.true.tss.transcript)
list.gr.loop.anchor.by.side <- create_loop_anchor_granges_by_side(df.loop.distinct.2mb)
gr.epd.promoter <- create_epd_tss_granges(df.promoter.epd.rn7.1based)

# 3. Prepare Primary Direct Tier for candidate exclusion baseline
promoter.window.flank.bp <- 1000L
gr.true.tss.promoter.window.1kb <- expand_tss_to_promoter_windows(gr.true.tss, flank.bp = promoter.window.flank.bp)
gr.epd.tss.promoter.window.1kb <- expand_tss_to_promoter_windows(gr.epd.promoter, flank.bp = promoter.window.flank.bp)

list.primary.direct.tier <- build_direct_promoter_tss_tier(
  gr.loop.anchor.by.side = list.gr.loop.anchor.by.side,
  gr.true.tss.annotation = gr.true.tss.promoter.window.1kb,
  gr.epd.annotation = gr.epd.tss.promoter.window.1kb,
  df.true.tss = df.true.tss.transcript,
  df.epd.promoter = df.promoter.epd.rn7.1based,
  df.loop.distinct.2mb = df.loop.distinct.2mb,
  evidence.definition = "primary_1kb",
  promoter.window.flank.bp = promoter.window.flank.bp
)

df.true.tss.lookup <- list.primary.direct.tier$true_tss_lookup
df.epd.promoter.lookup <- list.primary.direct.tier$epd_promoter_lookup
list.primary.direct.summary <- list.primary.direct.tier$summary
df.direct.promoter.tss.gene.assignment <- list.primary.direct.summary$gene_assignment
df.direct.promoter.tss.loop.summary <- list.primary.direct.summary$loop_summary
df.loop.anchor.index <- list.primary.direct.summary$anchor_index

################################################################################
# 4. Generate Proximal Promoter/TSS-Assignment Tier (1 bp - 200 kb)
################################################################################

proximal.max.distance.bp <- 200000L # 200kb

# Retrieve all non-overlapping true TSS records within the symmetric search window around either anchor.
df.proximal.true.tss.anchor.hit <- map_loop_anchors_dfr(
  list.gr.loop.anchor.by.side,
  find_proximal_anchor_annotation_pairs,
  gr.annotation = gr.true.tss,
  max.distance = proximal.max.distance.bp
)

df.proximal.true.tss.anchor.pair <- df.proximal.true.tss.anchor.hit %>%
  left_join(df.true.tss.lookup, by = "annotation_index") %>%
  transmute(
    loop_id,
    resolution,
    anchor_side,
    opposite_anchor_side,
    anchor_chr,
    anchor_start,
    anchor_end,
    anchor_width_bp,
    annotation_class = "true_TSS",
    annotation_source = str_c(source, "_GTF_transcript"),
    annotation_id = true_tss_id,
    annotation_chr,
    annotation_start,
    annotation_end,
    annotation_width_bp,
    annotation_strand,
    annotation_relative_to_anchor,
    anchor_annotation_edge_distance_bp,
    proximal_max_distance_bp,
    proximal_distance_tier = classify_proximal_distance_tier(anchor_annotation_edge_distance_bp),
    gene_id,
    gene_id_versioned,
    gene_name,
    gene_biotype,
    transcript_id,
    transcript_id_versioned,
    transcript_name,
    transcript_biotype,
    transcript_start,
    transcript_end,
    promoter_annotation_id = NA_character_,
    epd_promoter_name = NA_character_,
    epd_tss_start = NA_integer_,
    epd_tss_end = NA_integer_,
    annotation_source_coordinate_system = source_coordinate_system,
    analysis_coordinate_system
  ) %>%
  mutate(proximal_evidence_id = str_c(loop_id, anchor_side, annotation_source, annotation_id, sep = "|"), .before = 1)

# Repeat the proximity search for EPD promoter intervals
df.proximal.epd.promoter.anchor.hit <- map_loop_anchors_dfr(
  list.gr.loop.anchor.by.side,
  find_proximal_anchor_annotation_pairs,
  gr.annotation = gr.epd.promoter,
  max.distance = proximal.max.distance.bp
)

df.proximal.epd.promoter.anchor.pair <- df.proximal.epd.promoter.anchor.hit %>%
  left_join(df.epd.promoter.lookup, by = "annotation_index") %>%
  transmute(
    loop_id,
    resolution,
    anchor_side,
    opposite_anchor_side,
    anchor_chr,
    anchor_start,
    anchor_end,
    anchor_width_bp,
    annotation_class = "EPD_promoter",
    annotation_source = "EPDnew_promoter",
    annotation_id = promoter_annotation_id,
    annotation_chr,
    annotation_start,
    annotation_end,
    annotation_width_bp,
    annotation_strand,
    annotation_relative_to_anchor,
    anchor_annotation_edge_distance_bp,
    proximal_max_distance_bp,
    proximal_distance_tier = classify_proximal_distance_tier(anchor_annotation_edge_distance_bp),
    gene_id,
    gene_id_versioned = gene_id,
    gene_name,
    gene_biotype = NA_character_,
    transcript_id = NA_character_,
    transcript_id_versioned = NA_character_,
    transcript_name = NA_character_,
    transcript_biotype = NA_character_,
    transcript_start = NA_integer_,
    transcript_end = NA_integer_,
    promoter_annotation_id,
    epd_promoter_name,
    epd_tss_start,
    epd_tss_end,
    annotation_source_coordinate_system = source_coordinate_system,
    analysis_coordinate_system
  ) %>%
  mutate(proximal_evidence_id = str_c(loop_id, anchor_side, annotation_source, annotation_id, sep = "|"), .before = 1)

# Combine Ensembl and EPD proximal evidence
df.proximal.promoter.tss.anchor.pair <- bind_rows(
  df.proximal.true.tss.anchor.pair,
  df.proximal.epd.promoter.anchor.pair
) %>%
  left_join(
    df.loop.distinct.2mb %>%
      dplyr::select(loop_id, loop_chr1 = chr1, loop_anchor1_start = start1, loop_anchor1_end = end1, loop_chr2 = chr2, loop_anchor2_start = start2, loop_anchor2_end = end2),
    by = "loop_id"
  ) %>%
  mutate(
    inter_anchor_start = loop_anchor1_end + 1L,
    inter_anchor_end = loop_anchor2_start - 1L,
    loop_partition_midpoint = (loop_anchor1_end + loop_anchor2_start) / 2,
    annotation_midpoint = (annotation_start + annotation_end) / 2,
    has_nonempty_inter_anchor_interval = (loop_chr1 == loop_chr2 & inter_anchor_start <= inter_anchor_end),
    annotation_fully_within_inter_anchor = (has_nonempty_inter_anchor_interval & annotation_chr == loop_chr1 & annotation_start >= inter_anchor_start & annotation_end <= inter_anchor_end),
    annotation_points_toward_loop_interior = case_when(
      anchor_side == "anchor1" ~ annotation_relative_to_anchor == "right_of_anchor",
      anchor_side == "anchor2" ~ annotation_relative_to_anchor == "left_of_anchor",
      TRUE ~ FALSE
    ),
    distance_to_anchor1_edge_bp = if_else(annotation_fully_within_inter_anchor, annotation_start - loop_anchor1_end, NA_integer_),
    distance_to_anchor2_edge_bp = if_else(annotation_fully_within_inter_anchor, loop_anchor2_start - annotation_end, NA_integer_),
    assigned_anchor_edge_distance_bp = case_when(
      anchor_side == "anchor1" ~ distance_to_anchor1_edge_bp,
      anchor_side == "anchor2" ~ distance_to_anchor2_edge_bp,
      TRUE ~ NA_integer_
    ),
    opposite_anchor_edge_distance_bp = case_when(
      anchor_side == "anchor1" ~ distance_to_anchor2_edge_bp,
      anchor_side == "anchor2" ~ distance_to_anchor1_edge_bp,
      TRUE ~ NA_integer_
    ),
    assigned_anchor_closer_than_opposite = case_when(
      !annotation_fully_within_inter_anchor ~ FALSE,
      is.na(assigned_anchor_edge_distance_bp) | is.na(opposite_anchor_edge_distance_bp) ~ FALSE,
      TRUE ~ assigned_anchor_edge_distance_bp < opposite_anchor_edge_distance_bp
    ),
    anchor_edge_distance_tie = (annotation_fully_within_inter_anchor & !is.na(assigned_anchor_edge_distance_bp) & assigned_anchor_edge_distance_bp == opposite_anchor_edge_distance_bp),
    annotation_fully_on_assigned_midpoint_side = case_when(
      !annotation_fully_within_inter_anchor ~ FALSE,
      anchor_side == "anchor1" ~ annotation_end < loop_partition_midpoint,
      anchor_side == "anchor2" ~ annotation_start > loop_partition_midpoint,
      TRUE ~ FALSE
    ),
    annotation_spans_loop_midpoint = (annotation_fully_within_inter_anchor & annotation_start <= loop_partition_midpoint & annotation_end >= loop_partition_midpoint),
    assigned_anchor_is_strictly_nearest = (annotation_fully_within_inter_anchor & annotation_fully_on_assigned_midpoint_side & assigned_anchor_closer_than_opposite),
    is_inward_proximal_evidence = (annotation_fully_within_inter_anchor & annotation_points_toward_loop_interior),
    is_outward_proximal_evidence = !is_inward_proximal_evidence,
    is_inward_proximal_evidence_10kb = (is_inward_proximal_evidence & anchor_annotation_edge_distance_bp <= 10000L),
    is_midpoint_nearest_inward_proximal_evidence = (is_inward_proximal_evidence & annotation_fully_on_assigned_midpoint_side & assigned_anchor_is_strictly_nearest),
    is_midpoint_nearest_inward_proximal_evidence_10kb = (is_midpoint_nearest_inward_proximal_evidence & anchor_annotation_edge_distance_bp <= 10000L),
    is_ambiguous_inward_proximal_evidence = (is_inward_proximal_evidence & (annotation_spans_loop_midpoint | anchor_edge_distance_tie))
  )

assert_analysis_unique_key(df.proximal.promoter.tss.anchor.pair, "proximal_evidence_id", "Proximal promoter/TSS evidence identifiers are not unique.")

# Generate gene assignment and secondary inward candidate tier
df.proximal.promoter.tss.gene.assignment <- df.proximal.promoter.tss.anchor.pair %>%
  group_by(loop_id, resolution, anchor_side, opposite_anchor_side, gene_id) %>%
  summarise(
    gene_name = {
      gene.names <- sort(unique(na.omit(gene_name)))
      if (length(gene.names) == 0L) NA_character_ else gene.names[1]
    },
    min_proximal_distance_bp = min(anchor_annotation_edge_distance_bp),
    n_proximal_evidence_records = n(),
    has_proximal_true_tss = any(annotation_class == "true_TSS"),
    has_proximal_epd_promoter = any(annotation_class == "EPD_promoter"),
    has_inward_proximal_evidence = any(is_inward_proximal_evidence),
    has_inward_proximal_evidence_10kb = any(is_inward_proximal_evidence_10kb),
    has_midpoint_nearest_inward_proximal_evidence = any(is_midpoint_nearest_inward_proximal_evidence),
    has_midpoint_nearest_inward_proximal_evidence_10kb = any(is_midpoint_nearest_inward_proximal_evidence_10kb),
    .groups = "drop"
  ) %>%
  left_join(
    df.direct.promoter.tss.gene.assignment %>%
      mutate(has_direct_same_anchor_gene = TRUE) %>%
      dplyr::select(loop_id, anchor_side, gene_id, has_direct_same_anchor_gene),
    by = c("loop_id", "anchor_side", "gene_id")
  ) %>%
  mutate(has_direct_same_anchor_gene = coalesce(has_direct_same_anchor_gene, FALSE))

df.secondary.inward.proximal.gene.assignment <- df.proximal.promoter.tss.gene.assignment %>%
  filter(has_midpoint_nearest_inward_proximal_evidence_10kb & !has_direct_same_anchor_gene)

# Cumulative distance threshold summary (Sensitivity Quantification)
proximal.distance.thresholds.bp <- c(1000L, 5000L, 10000L, 25000L, 50000L, 100000L, 200000L)
direct.supported.loop.ids <- df.direct.promoter.tss.loop.summary %>%
  filter(has_any_direct_promoter_tss) %>%
  pull(loop_id)

df.proximal.promoter.tss.cumulative.threshold.summary <- map_dfr(proximal.distance.thresholds.bp, function(distance.threshold) {
  df.threshold <- df.proximal.promoter.tss.anchor.pair %>% filter(anchor_annotation_edge_distance_bp <= distance.threshold)
  df.threshold.inward <- df.threshold %>% filter(is_inward_proximal_evidence)
  df.threshold.inward.resolved <- df.threshold %>% filter(is_midpoint_nearest_inward_proximal_evidence)
  proximal.loop.ids <- unique(df.threshold$loop_id)
  inward.loop.ids <- unique(df.threshold.inward$loop_id)
  inward.resolved.loop.ids <- unique(df.threshold.inward.resolved$loop_id)

  tibble(
    max_anchor_edge_distance_bp = distance.threshold,
    n_annotation_records = nrow(df.threshold),
    n_loop_anchor_pairs = n_distinct(df.threshold$loop_id, df.threshold$anchor_side),
    n_loop_anchor_gene_assignments = n_distinct(df.threshold$loop_id, df.threshold$anchor_side, df.threshold$gene_id),
    n_loops_with_proximity = length(proximal.loop.ids),
    n_inward_annotation_records = nrow(df.threshold.inward),
    n_inward_loop_anchor_pairs = n_distinct(df.threshold.inward$loop_id, df.threshold.inward$anchor_side),
    n_loops_with_inward_proximity = length(inward.loop.ids),
    n_midpoint_nearest_inward_annotation_records = nrow(df.threshold.inward.resolved),
    n_midpoint_nearest_inward_loop_anchor_pairs = n_distinct(df.threshold.inward.resolved$loop_id, df.threshold.inward.resolved$anchor_side),
    n_loops_with_midpoint_nearest_inward_proximity = length(inward.resolved.loop.ids),
    pct_pooled_loops_with_proximity = round(100 * length(proximal.loop.ids) / nrow(df.loop.distinct.2mb), 1),
    n_direct_and_proximity_loops = length(intersect(proximal.loop.ids, direct.supported.loop.ids)),
    n_proximity_without_direct_loops = length(setdiff(proximal.loop.ids, direct.supported.loop.ids)),
    n_direct_without_proximity_loops = length(setdiff(direct.supported.loop.ids, proximal.loop.ids)),
    n_without_direct_or_proximity_loops = nrow(df.loop.distinct.2mb) - length(union(proximal.loop.ids, direct.supported.loop.ids))
  )
})

cat("\n=================================================================\n")
cat("SENSITIVITY ANALYSIS 01_2: Cumulative Distance Threshold Summary\n")
cat("=================================================================\n")
print(df.proximal.promoter.tss.cumulative.threshold.summary)

cat("\n=================================================================\n")
cat("Proximal Evidence Retention Summary:\n")
cat("  - Total 1-200kb exploratory gene assignments :", scales::comma(nrow(df.proximal.promoter.tss.gene.assignment)), "\n")
cat("  - Secondary inward <=10-kb gene assignments  :", scales::comma(nrow(df.secondary.inward.proximal.gene.assignment)), "\n")
cat("=================================================================\n")
