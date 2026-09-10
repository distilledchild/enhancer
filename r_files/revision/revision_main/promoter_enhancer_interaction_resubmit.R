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
# Resubmission analysis: pooled rat frontal cortex chromatin-loop annotation
#
# Revised implementation order
# 1. Load rn7, chr-prefixed, one-based inclusive coordinate inputs.
# 2. Generate strand-aware true TSS coordinates from Ensembl transcripts.
# 3. Preserve every direct promoter/TSS-anchor overlap.
# 4. Generate a separate proximal-assignment tier.
# 5. Recalculate ATAC support after excluding true TSS regions.
# 6. Add transcript-contained and related positional flags.
# 7. Add predicted CTCF motif annotations without using motif thresholds for
#    loop selection or regulatory classification.
# 8. Define layered loop-evidence categories and downstream summaries.
#
# CTCF motifs are treated as structural evidence, whereas ATAC overlap supports
# open chromatin but does not by itself demonstrate enhancer function. Strain
# sharing and sequencing-depth analyses remain exploratory and appear last.
#
# This script follows the object naming, numbered section structure, and
# explicit intermediate checks used in:
# - enhancer_promoter_interaction_figures.R
# - enhancer_promoter_interaction.R
################################################################################

########################
# 0. Directories and files
########################

dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

# Coordinate preprocessing is isolated from the downstream analysis and cached.
coord.prep.script <- file.path(analysis.dir, "promoter_enhancer_interaction_resubmit_coord_prep.R")

# Require the coordinate-preparation script before using its cached outputs.
if (!file.exists(coord.prep.script)) {
  stop("Cannot locate the coordinate-preparation script: ", coord.prep.script, call. = FALSE)
} else {
  message("Located coordinate-preparation script: ", coord.prep.script)
}

################################################################################
# 1. Load coordinate-normalized cache
#
# The coordinate-preparation script is run only when one or more expected cache
# files are absent. Set REBUILD_COORD_CACHE=1 to force reconstruction after a
# source-data or coordinate-processing change.
################################################################################

coord.cache.object.names <- coordinate_cache_object_names()
expected.cache.files <- file.path(coord.cache.dir, paste0(coord.cache.object.names, ".rds"))

# Verify that all expected cache .rds files exist; stop and list any missing files.
missing.cache.files <- expected.cache.files[!file.exists(expected.cache.files)]
if (length(missing.cache.files) > 0L) {
  stop("Missing required cache .rds file(s):\n", paste(missing.cache.files, collapse = "\n"), call. = FALSE)
}

# Load all cached objects into environment
load_coordinate_cache_objects(cache.dir = coord.cache.dir, object.names = coord.cache.object.names, envir = environment())

# Extract downstream input file paths
input.file.paths <- setNames(df.analysis.input.files$input_path, df.analysis.input.files$input_name)
library.complexity.file <- input.file.paths[["library_complexity"]]
genetic.distance.file <- Sys.getenv("HRDP_GENETIC_DISTANCE_FILE", unset = input.file.paths[["genetic_distance"]])

message("Successfully verified and loaded ", length(coord.cache.object.names), " cache objects from: ", coord.cache.dir)

################################################################################
# 2. Strand-aware true TSS generation
#
# Description: Generate 1-bp true TSS coordinates (+ strand: start, - strand: end)
# Input:       df.transcript.ensembl.rn7.1based
# Output:      df.true.tss.transcript, gr.true.tss, df.true.tss.summary
################################################################################

# Extract strand-aware 1-bp TSS coordinates (+: start, -: end) and metadata per transcript
df.true.tss.transcript <- build_true_tss_annotation(df.transcript.ensembl.rn7.1based)
gr.true.tss <- create_true_tss_granges(df.true.tss.transcript)

df.true.tss.summary <- tibble(
  metric = c(
    "normalized_Ensembl_transcripts", "transcript_level_true_TSS_records", "unique_Ensembl_genes",
    "unique_strand_aware_TSS_sites", "unique_genomic_TSS_positions", "Ensembl_canonical_transcript_TSS_records",
    "genes_with_Ensembl_canonical_transcript_TSS", "plus_strand_transcript_TSS", "minus_strand_transcript_TSS",
    "transcripts_with_start_le_end", "plus_strand_TSS_at_transcript_start", "minus_strand_TSS_at_transcript_end", "EPD_promoter_TSS_records"
  ),
  n = c(
    nrow(df.transcript.ensembl.rn7.1based), nrow(df.true.tss.transcript), n_distinct(df.true.tss.transcript$gene_id),
    n_distinct(df.true.tss.transcript$chr, df.true.tss.transcript$true_tss_start, df.true.tss.transcript$strand),
    n_distinct(df.true.tss.transcript$chr, df.true.tss.transcript$true_tss_start),
    sum(df.true.tss.transcript$is_ensembl_canonical), n_distinct(df.true.tss.transcript$gene_id[df.true.tss.transcript$is_ensembl_canonical]),
    sum(df.true.tss.transcript$strand == "+"), sum(df.true.tss.transcript$strand == "-"),
    sum(df.true.tss.transcript$transcript_start <= df.true.tss.transcript$transcript_end),
    sum(df.true.tss.transcript$strand == "+" & df.true.tss.transcript$true_tss_start == df.true.tss.transcript$transcript_start),
    sum(df.true.tss.transcript$strand == "-" & df.true.tss.transcript$true_tss_start == df.true.tss.transcript$transcript_end),
    nrow(df.promoter.epd.rn7.1based)
  )
)

message("Generated ", nrow(df.true.tss.transcript), " strand-aware transcript TSS records across ", n_distinct(df.true.tss.transcript$gene_id), " Ensembl genes.")

################################################################################
# 3. Define strict and primary direct promoter/TSS-anchor overlap tiers
#
# Description: Map loop anchors to promoter/TSS using strict (exact TSS/81-bp EPD)
#              and primary (TSS +/-1-kb) direct overlap tiers.
# Input:       df.loop.universe, df.true.tss.transcript, gr.true.tss, df.promoter.epd.rn7.1based
# Output:      strict.direct.tier, primary.direct.tier, df.direct.promoter.tss.anchor.overlap
################################################################################

list.gr.loop.anchor.by.side <- create_loop_anchor_granges_by_side(df.loop.universe)
gr.epd.promoter <- create_epd_promoter_granges(df.promoter.epd.rn7.1based)

####################################################################################
# Build the strict exact-coordinate promoter/TSS tier once for both Ensembl
# transcript TSS points and EPD promoter intervals.
####################################################################################
strict.direct.tier <- build_direct_promoter_tss_tier(
  gr.loop.anchor.by.side = list.gr.loop.anchor.by.side,
  gr.true.tss.annotation = gr.true.tss,
  gr.epd.annotation = gr.epd.promoter,
  df.true.tss = df.true.tss.transcript,
  df.epd.promoter = df.promoter.epd.rn7.1based,
  df.loop.universe = df.loop.universe,
  evidence.definition = "strict"
)

# Unpacking data from strict.direct.tier
df.true.tss.lookup <- strict.direct.tier$true_tss_lookup
df.epd.promoter.lookup <- strict.direct.tier$epd_promoter_lookup
df.strict.direct.true.tss.anchor.overlap <- strict.direct.tier$true_tss_overlap
df.strict.direct.epd.promoter.anchor.overlap <- strict.direct.tier$epd_promoter_overlap
df.strict.direct.promoter.tss.anchor.overlap <- strict.direct.tier$combined_overlap
strict.direct.summary <- strict.direct.tier$summary

# Confirm that exact Ensembl TSS overlaps remain one-base intervals.
assert_analysis_condition(!any(df.strict.direct.true.tss.anchor.overlap$direct_overlap_bp != 1L), "A one-base true TSS has an unexpected direct-overlap width.")

df.strict.direct.promoter.tss.gene.assignment <- strict.direct.summary$gene_assignment
df.strict.direct.promoter.tss.anchor.count <- strict.direct.summary$anchor_count
df.strict.direct.promoter.tss.anchor.summary <- strict.direct.summary$anchor_summary
df.strict.direct.promoter.tss.anchor.wide <- strict.direct.summary$anchor_wide
df.strict.direct.promoter.tss.loop.gene.count <- strict.direct.summary$loop_gene_count
df.strict.direct.promoter.tss.loop.summary <- strict.direct.summary$loop_summary
df.strict.direct.promoter.tss.summary <- strict.direct.summary$summary

####################################################################################
# Build the primary tier from +/-1-kb TSS windows for both annotation sources.
####################################################################################
promoter.window.flank.bp <- 1000L
gr.epd.tss <- create_epd_tss_granges(df.promoter.epd.rn7.1based)
gr.true.tss.promoter.window.1kb <- expand_tss_to_promoter_windows(gr.true.tss, flank.bp = promoter.window.flank.bp)
gr.epd.tss.promoter.window.1kb <- expand_tss_to_promoter_windows(gr.epd.tss, flank.bp = promoter.window.flank.bp)

primary.direct.tier <- build_direct_promoter_tss_tier(
  gr.loop.anchor.by.side = list.gr.loop.anchor.by.side,
  gr.true.tss.annotation = gr.true.tss.promoter.window.1kb,
  gr.epd.annotation = gr.epd.tss.promoter.window.1kb,
  df.true.tss = df.true.tss.transcript,
  df.epd.promoter = df.promoter.epd.rn7.1based,
  df.loop.universe = df.loop.universe,
  evidence.definition = "primary_1kb",
  promoter.window.flank.bp = promoter.window.flank.bp
)

df.primary.direct.true.tss.anchor.overlap <- primary.direct.tier$true_tss_overlap
df.primary.direct.epd.promoter.anchor.overlap <- primary.direct.tier$epd_promoter_overlap
df.primary.direct.promoter.tss.anchor.overlap <- primary.direct.tier$combined_overlap
primary.direct.summary <- primary.direct.tier$summary

# From this point onward, df.direct.* means the primary TSS +/-1-kb direct tier.
df.direct.true.tss.anchor.overlap <- df.primary.direct.true.tss.anchor.overlap
df.direct.epd.promoter.anchor.overlap <- df.primary.direct.epd.promoter.anchor.overlap
df.direct.promoter.tss.anchor.overlap <- df.primary.direct.promoter.tss.anchor.overlap
df.direct.promoter.tss.gene.assignment <- primary.direct.summary$gene_assignment
df.loop.anchor.index <- primary.direct.summary$anchor_index
df.direct.promoter.tss.anchor.count <- primary.direct.summary$anchor_count
df.direct.promoter.tss.anchor.summary <- primary.direct.summary$anchor_summary
df.direct.promoter.tss.anchor.wide <- primary.direct.summary$anchor_wide
df.direct.promoter.tss.loop.gene.count <- primary.direct.summary$loop_gene_count
df.direct.promoter.tss.loop.summary <- primary.direct.summary$loop_summary
df.direct.promoter.tss.summary <- primary.direct.summary$summary

# Compare loop-level membership between strict and primary direct promoter definitions
df.strict.vs.primary.promoter.window.loop.comparison <- df.loop.universe %>%
  dplyr::select(loop_id, resolution) %>%
  left_join(df.strict.direct.promoter.tss.loop.summary %>% transmute(loop_id, strict_exact_direct = has_any_direct_promoter_tss), by = "loop_id") %>%
  left_join(df.direct.promoter.tss.loop.summary %>% transmute(loop_id, primary_1kb_direct = has_any_direct_promoter_tss), by = "loop_id") %>%
  mutate(
    direct_definition_membership = case_when(
      strict_exact_direct & primary_1kb_direct ~ "strict_and_primary",
      strict_exact_direct ~ "strict_only",
      primary_1kb_direct ~ "primary_1kb_only",
      TRUE ~ "neither"
    )
  )

df.strict.vs.primary.promoter.window.loop.comparison %>% count(direct_definition_membership)
# 1 neither                      14678
# 2 primary_1kb_only              1253
# 3 strict_and_primary           15090

# Extract unique assignment keys (loop_id | anchor_side | gene_id) for strict and primary tiers
strict.assignment.keys <- df.strict.direct.promoter.tss.gene.assignment %>%
  transmute(key = str_c(loop_id, anchor_side, gene_id, sep = "|")) %>%
  pull(key)
primary.assignment.keys <- df.direct.promoter.tss.gene.assignment %>%
  transmute(key = str_c(loop_id, anchor_side, gene_id, sep = "|")) %>%
  pull(key)

# Construct comparative summary table of strict vs primary promoter-anchor assignments
df.strict.vs.primary.promoter.window.summary <- bind_rows(
  df.strict.direct.promoter.tss.summary %>% mutate(evidence_definition = "strict_exact_TSS_or_81bp_EPD", .before = 1),
  df.direct.promoter.tss.summary %>% mutate(evidence_definition = "primary_TSS_plus_minus_1kb", .before = 1)
) %>%
  bind_rows(
    tibble(
      evidence_definition = "strict_vs_primary_assignment_overlap",
      metric = c("strict_loop_anchor_gene_assignments", "primary_loop_anchor_gene_assignments", "shared_loop_anchor_gene_assignments", "primary_only_loop_anchor_gene_assignments", "strict_only_loop_anchor_gene_assignments"),
      n = c(length(unique(strict.assignment.keys)), length(unique(primary.assignment.keys)), length(intersect(strict.assignment.keys, primary.assignment.keys)), length(setdiff(primary.assignment.keys, strict.assignment.keys)), length(setdiff(strict.assignment.keys, primary.assignment.keys)))
    )
  )

# Define coordinate rules and interpretations for promoter-anchor evidence tiers
df.promoter.anchor.assignment.definitions <- tribble(
  ~evidence_tier, ~coordinate_rule, ~interpretation,
  "strict_sensitivity", paste0(
    "Direct anchor overlap with a one-base Ensembl TSS or the original 81-bp ",
    "EPD promoter interval."
  ), paste0(
    "Strict coordinate sensitivity analysis; retained for robustness checks, ",
    "not as a functional-validation claim."
  ),
  "primary_direct", paste0(
    "Direct anchor overlap with a promoter window defined as TSS +/-1 kb from ",
    "either Ensembl transcript TSS or EPD TSS."
  ), paste0(
    "Primary promoter-associated anchor evidence used by downstream revised ",
    "ATAC, category, gene-resource, and GO analyses."
  ),
  "secondary_inward_proximal", paste0(
    "Non-overlapping TSS or EPD promoter interval within 10 kb of an anchor and ",
    "fully located in the inter-anchor interval on that anchor's side of the ",
    "loop midpoint. The assigned anchor must be strictly closer than the ",
    "opposite anchor; same-anchor/gene primary direct assignments are excluded."
  ), paste0(
    "Secondary proximity-supported candidate tier; weaker than direct overlap ",
    "and not treated as validated target-gene evidence."
  ),
  "exploratory_proximal_catalog", paste0(
    "Any non-overlapping TSS or EPD promoter interval 1-200 kb from an anchor in ",
    "either direction."
  ), paste0(
    "Exploratory sensitivity catalog only; unsuitable for strong regulatory or ",
    "target-gene claims."
  )
)

# Print direct promoter/TSS assignment summary message
message(
  "Direct promoter/TSS overlap retained ",
  nrow(df.direct.promoter.tss.gene.assignment),
  " loop-anchor-gene assignments across ",
  sum(
    df.direct.promoter.tss.loop.summary$has_any_direct_promoter_tss
  ),
  " pooled loops."
)

################################################################################
# 4. Generate a separate proximal promoter/TSS-assignment tier
#
# Description: Build secondary inward proximal (10-kb) and exploratory (200-kb)
#              promoter/TSS tiers for non-overlapping candidate annotations.
# Input:       list.gr.loop.anchor.by.side, gr.true.tss, gr.epd.promoter, df.loop.universe
# Output:      df.proximal.promoter.tss.anchor.pair, df.secondary.inward.proximal.*, df.proximal.promoter.tss.summary
################################################################################

proximal.max.distance.bp <- 200000L

# Retrieve all non-overlapping true TSS records within the symmetric search
# window around either anchor. Distances are measured from the nearest anchor
# boundary, with an immediately adjacent annotation assigned a distance of 1 bp.
df.proximal.true.tss.anchor.hit <- map_loop_anchors_dfr(
  list.gr.loop.anchor.by.side,
  find_proximal_anchor_annotation_pairs,
  gr.annotation = gr.true.tss,
  max.distance = proximal.max.distance.bp
)
df.proximal.true.tss.anchor.hit %>% head(3)

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

# Repeat the proximity search for EPD promoter intervals while keeping EPD and
# Ensembl evidence separate in the long table.
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

# Preserve every proximal annotation record. Loop geometry distinguishes inward
# candidates located between the two anchors from outward candidates located
# outside the loop span. Rank fields expose nearest records without deleting more
# distant candidates or tied annotations from the exploratory 200-kb catalog.
df.proximal.promoter.tss.anchor.pair <- bind_rows(
  df.proximal.true.tss.anchor.pair,
  df.proximal.epd.promoter.anchor.pair
) %>%
  left_join(
    df.loop.universe %>%
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
    annotation_fully_on_assigned_midpoint_side = case_when(
      anchor_side == "anchor1" ~ annotation_end < loop_partition_midpoint,
      anchor_side == "anchor2" ~ annotation_start > loop_partition_midpoint,
      TRUE ~ FALSE
    ),
    annotation_spans_loop_midpoint = (annotation_start <= loop_partition_midpoint & annotation_end >= loop_partition_midpoint),
    assigned_anchor_is_strictly_nearest = (!is.na(assigned_anchor_edge_distance_bp) & !is.na(opposite_anchor_edge_distance_bp) & assigned_anchor_edge_distance_bp < opposite_anchor_edge_distance_bp),
    anchor_edge_distance_tie = (!is.na(assigned_anchor_edge_distance_bp) & !is.na(opposite_anchor_edge_distance_bp) & assigned_anchor_edge_distance_bp == opposite_anchor_edge_distance_bp),
    is_inward_proximal_evidence = (annotation_fully_within_inter_anchor & annotation_points_toward_loop_interior),
    is_outward_proximal_evidence = !is_inward_proximal_evidence,
    is_inward_proximal_evidence_10kb = (is_inward_proximal_evidence & anchor_annotation_edge_distance_bp <= 10000L),
    is_midpoint_nearest_inward_proximal_evidence = (is_inward_proximal_evidence & annotation_fully_on_assigned_midpoint_side & assigned_anchor_is_strictly_nearest),
    is_midpoint_nearest_inward_proximal_evidence_10kb = (is_midpoint_nearest_inward_proximal_evidence & anchor_annotation_edge_distance_bp <= 10000L),
    is_ambiguous_inward_proximal_evidence = (is_inward_proximal_evidence & (annotation_spans_loop_midpoint | anchor_edge_distance_tie))
  ) %>%
  group_by(loop_id, anchor_side, annotation_class) %>%
  mutate(
    distance_rank_within_anchor_annotation_class = dense_rank(anchor_annotation_edge_distance_bp),
    is_nearest_within_anchor_annotation_class = distance_rank_within_anchor_annotation_class == 1L
  ) %>%
  ungroup() %>%
  arrange(loop_id, anchor_side, anchor_annotation_edge_distance_bp, gene_id, annotation_class, annotation_id)

# Enforce a unique identifier for every proximal promoter/TSS evidence record.
assert_analysis_unique_key(df.proximal.promoter.tss.anchor.pair, "proximal_evidence_id", "Proximal promoter/TSS evidence identifiers are not unique.")
# Reject missing, out-of-range, or internally inconsistent proximal distances.
assert_analysis_condition(
  !any(is.na(df.proximal.promoter.tss.anchor.pair$proximal_distance_tier)) &&
    !any(df.proximal.promoter.tss.anchor.pair$anchor_annotation_edge_distance_bp < 1L) &&
    !any(df.proximal.promoter.tss.anchor.pair$anchor_annotation_edge_distance_bp > proximal.max.distance.bp) &&
    !any(df.proximal.promoter.tss.anchor.pair$is_inward_proximal_evidence & df.proximal.promoter.tss.anchor.pair$assigned_anchor_edge_distance_bp != df.proximal.promoter.tss.anchor.pair$anchor_annotation_edge_distance_bp),
  "The proximal evidence table contains an invalid distance or tier."
)

# Collapse only duplicate transcript/promoter evidence for the same gene. Every
# loop-anchor-gene candidate remains, and genes with a direct assignment at the
# same anchor are flagged rather than removed from the proximal table.
df.proximal.promoter.tss.gene.assignment <- df.proximal.promoter.tss.anchor.pair %>%
  group_by(
    loop_id,
    resolution,
    anchor_side,
    opposite_anchor_side,
    gene_id
  ) %>%
  summarise(
    gene_name = sort(unique(gene_name))[1],
    min_proximal_distance_bp = min(anchor_annotation_edge_distance_bp),
    n_proximal_evidence_records = n(),
    has_proximal_true_tss = any(annotation_class == "true_TSS"),
    has_proximal_epd_promoter = any(annotation_class == "EPD_promoter"),
    n_proximal_true_tss_transcripts = n_distinct(transcript_id_versioned, na.rm = TRUE),
    n_proximal_epd_promoters = n_distinct(promoter_annotation_id, na.rm = TRUE),
    has_inward_proximal_evidence = any(is_inward_proximal_evidence),
    has_outward_proximal_evidence = any(is_outward_proximal_evidence),
    has_inward_proximal_evidence_10kb = any(is_inward_proximal_evidence_10kb),
    has_midpoint_nearest_inward_proximal_evidence = any(is_midpoint_nearest_inward_proximal_evidence),
    has_midpoint_nearest_inward_proximal_evidence_10kb = any(is_midpoint_nearest_inward_proximal_evidence_10kb),
    has_ambiguous_inward_proximal_evidence = any(is_ambiguous_inward_proximal_evidence),
    min_inward_proximal_distance_bp = if (any(is_inward_proximal_evidence)) min(anchor_annotation_edge_distance_bp[is_inward_proximal_evidence]) else NA_integer_,
    min_midpoint_nearest_inward_proximal_distance_bp = if (any(is_midpoint_nearest_inward_proximal_evidence)) min(anchor_annotation_edge_distance_bp[is_midpoint_nearest_inward_proximal_evidence]) else NA_integer_,
    .groups = "drop"
  ) %>%
  mutate(proximal_distance_tier = classify_proximal_distance_tier(min_proximal_distance_bp)) %>%
  left_join(
    df.direct.promoter.tss.gene.assignment %>% transmute(loop_id, anchor_side, gene_id, has_direct_same_anchor_gene = TRUE),
    by = c("loop_id", "anchor_side", "gene_id")
  ) %>%
  mutate(
    has_direct_same_anchor_gene = coalesce(has_direct_same_anchor_gene, FALSE),
    is_secondary_inward_proximal_candidate_10kb = (has_midpoint_nearest_inward_proximal_evidence_10kb & !has_direct_same_anchor_gene)
  ) %>%
  group_by(loop_id, anchor_side) %>%
  mutate(
    proximal_gene_distance_rank_at_anchor = dense_rank(min_proximal_distance_bp),
    is_nearest_proximal_gene_at_anchor = proximal_gene_distance_rank_at_anchor == 1L
  ) %>%
  ungroup() %>%
  arrange(loop_id, anchor_side, min_proximal_distance_bp, gene_id)

# The secondary tier is limited to midpoint-resolved, nearest-anchor inward
# candidates within 10 kb and excludes loop-anchor-gene assignments already
# represented by the primary direct tier.
df.secondary.inward.proximal.gene.assignment <- df.proximal.promoter.tss.gene.assignment %>%
  filter(is_secondary_inward_proximal_candidate_10kb) %>%
  arrange(loop_id, anchor_side, min_midpoint_nearest_inward_proximal_distance_bp, gene_id)

df.proximal.promoter.tss.anchor.count <- df.proximal.promoter.tss.anchor.pair %>%
  group_by(loop_id, resolution, anchor_side, opposite_anchor_side) %>%
  summarise(
    n_proximal_evidence_records = n(),
    n_proximal_genes = n_distinct(gene_id),
    min_proximal_distance_bp = min(anchor_annotation_edge_distance_bp),
    n_proximal_true_tss_records = sum(annotation_class == "true_TSS"),
    n_proximal_epd_promoter_records = sum(annotation_class == "EPD_promoter"),
    n_proximal_1bp_to_10kb = sum(proximal_distance_tier == "01_1bp_to_10kb"),
    n_proximal_gt10kb_to_50kb = sum(proximal_distance_tier == "02_gt10kb_to_50kb"),
    n_proximal_gt50kb_to_200kb = sum(proximal_distance_tier == "03_gt50kb_to_200kb"),
    n_inward_proximal_evidence_records = sum(is_inward_proximal_evidence),
    n_outward_proximal_evidence_records = sum(is_outward_proximal_evidence),
    n_inward_proximal_evidence_records_10kb = sum(is_inward_proximal_evidence_10kb),
    n_midpoint_nearest_inward_proximal_evidence_records = sum(is_midpoint_nearest_inward_proximal_evidence),
    n_midpoint_nearest_inward_proximal_evidence_records_10kb = sum(is_midpoint_nearest_inward_proximal_evidence_10kb),
    n_ambiguous_inward_proximal_evidence_records = sum(is_ambiguous_inward_proximal_evidence),
    .groups = "drop"
  )

df.secondary.inward.proximal.anchor.count <- df.secondary.inward.proximal.gene.assignment %>%
  group_by(loop_id, resolution, anchor_side, opposite_anchor_side) %>%
  summarise(
    n_secondary_inward_proximal_genes_10kb = n_distinct(gene_id),
    min_secondary_inward_proximal_distance_bp = min(min_midpoint_nearest_inward_proximal_distance_bp),
    .groups = "drop"
  )

df.proximal.promoter.tss.anchor.summary <- df.loop.anchor.index %>%
  left_join(df.proximal.promoter.tss.anchor.count, by = c("loop_id", "resolution", "anchor_side", "opposite_anchor_side")) %>%
  left_join(df.secondary.inward.proximal.anchor.count, by = c("loop_id", "resolution", "anchor_side", "opposite_anchor_side")) %>%
  mutate(
    across(c(n_proximal_evidence_records, n_proximal_genes, n_proximal_true_tss_records, n_proximal_epd_promoter_records, n_proximal_1bp_to_10kb, n_proximal_gt10kb_to_50kb, n_proximal_gt50kb_to_200kb, n_inward_proximal_evidence_records, n_outward_proximal_evidence_records, n_inward_proximal_evidence_records_10kb, n_midpoint_nearest_inward_proximal_evidence_records, n_midpoint_nearest_inward_proximal_evidence_records_10kb, n_ambiguous_inward_proximal_evidence_records, n_secondary_inward_proximal_genes_10kb), ~ coalesce(.x, 0L)),
    has_any_proximal_promoter_tss = n_proximal_evidence_records > 0L,
    has_any_inward_proximal_promoter_tss = n_inward_proximal_evidence_records > 0L,
    has_secondary_inward_proximal_10kb = n_secondary_inward_proximal_genes_10kb > 0L,
    min_proximal_distance_bp = if_else(has_any_proximal_promoter_tss, min_proximal_distance_bp, NA_integer_),
    min_secondary_inward_proximal_distance_bp = if_else(has_secondary_inward_proximal_10kb, min_secondary_inward_proximal_distance_bp, NA_integer_)
  ) %>%
  arrange(loop_id, anchor_side)

# Verify that the proximal summary retains both anchors for every pooled loop.
assert_analysis_row_count(
  df.proximal.promoter.tss.anchor.summary,
  2L * nrow(df.loop.universe),
  "Proximal anchor summary does not retain both sides of every loop."
)

df.proximal.promoter.tss.anchor.wide <- df.proximal.promoter.tss.anchor.summary %>%
  dplyr::select(
    loop_id, anchor_side, n_proximal_evidence_records, n_proximal_genes, min_proximal_distance_bp,
    has_any_proximal_promoter_tss, n_secondary_inward_proximal_genes_10kb, min_secondary_inward_proximal_distance_bp, has_secondary_inward_proximal_10kb
  ) %>%
  pivot_wider(
    names_from = anchor_side,
    values_from = -loop_id,
    names_glue = "{.value}_{anchor_side}"
  )

df.proximal.promoter.tss.loop.gene.count <- df.proximal.promoter.tss.gene.assignment %>%
  group_by(loop_id) %>%
  summarise(
    n_proximal_genes_across_anchors = n_distinct(gene_id),
    .groups = "drop"
  )

df.secondary.inward.proximal.loop.gene.count <- df.secondary.inward.proximal.gene.assignment %>%
  group_by(loop_id) %>%
  summarise(
    n_secondary_inward_proximal_genes_10kb_across_anchors = n_distinct(gene_id),
    .groups = "drop"
  )

# Combine loop-level direct and proximal flags without promoting proximity to
# direct evidence. This table will feed the later ATAC and category revisions.
df.proximal.promoter.tss.loop.summary <- df.loop.universe %>%
  left_join(df.proximal.promoter.tss.anchor.wide, by = "loop_id") %>%
  left_join(df.proximal.promoter.tss.loop.gene.count, by = "loop_id") %>%
  left_join(df.secondary.inward.proximal.loop.gene.count, by = "loop_id") %>%
  left_join(
    df.direct.promoter.tss.loop.summary %>%
      dplyr::select(
        loop_id,
        n_direct_anchor_sides,
        has_any_direct_promoter_tss,
        has_direct_promoter_tss_both_anchors
      ),
    by = "loop_id"
  ) %>%
  mutate(
    n_proximal_genes_across_anchors = coalesce(n_proximal_genes_across_anchors, 0L),
    n_secondary_inward_proximal_genes_10kb_across_anchors = coalesce(n_secondary_inward_proximal_genes_10kb_across_anchors, 0L),
    n_proximal_anchor_sides = as.integer(has_any_proximal_promoter_tss_anchor1) + as.integer(has_any_proximal_promoter_tss_anchor2),
    has_any_proximal_promoter_tss = n_proximal_anchor_sides > 0L,
    has_proximal_promoter_tss_both_anchors = n_proximal_anchor_sides == 2L,
    has_any_direct_or_proximal_promoter_tss = (has_any_direct_promoter_tss | has_any_proximal_promoter_tss),
    n_secondary_inward_proximal_anchor_sides_10kb = as.integer(has_secondary_inward_proximal_10kb_anchor1) + as.integer(has_secondary_inward_proximal_10kb_anchor2),
    has_any_secondary_inward_proximal_10kb = n_secondary_inward_proximal_anchor_sides_10kb > 0L,
    has_secondary_inward_proximal_both_anchors_10kb = n_secondary_inward_proximal_anchor_sides_10kb == 2L,
    promoter_tss_assignment_tier = case_when(
      has_any_direct_promoter_tss & has_any_secondary_inward_proximal_10kb ~ "primary_direct_and_secondary_inward_10kb",
      has_any_direct_promoter_tss ~ "direct_only",
      has_any_secondary_inward_proximal_10kb ~ "secondary_inward_10kb_only",
      has_any_proximal_promoter_tss ~ "exploratory_proximal_200kb_only",
      TRUE ~ "no_direct_or_proximal_evidence"
    )
  )

# Count distinct promoter sites rather than raw transcript rows. A site is a
# unique gene/TSS-coordinate combination, so Ensembl isoforms that share a TSS
# and matching EPD evidence do not automatically inflate anchor multiplicity.
df.primary.direct.promoter.site.count.by.anchor <- df.direct.promoter.tss.anchor.overlap %>%
  mutate(
    assigned_tss_coordinate = coalesce(tss_start, epd_tss_start),
    promoter_site_key = str_c(gene_id, annotation_chr, assigned_tss_coordinate, sep = "|")
  ) %>%
  group_by(loop_id, anchor_side) %>%
  summarise(
    n_primary_direct_promoter_sites = n_distinct(promoter_site_key),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = anchor_side,
    values_from = n_primary_direct_promoter_sites,
    names_glue = "n_primary_direct_promoter_sites_{anchor_side}",
    values_fill = 0L
  )

df.secondary.promoter.site.count.by.anchor <- df.proximal.promoter.tss.anchor.pair %>%
  filter(is_midpoint_nearest_inward_proximal_evidence_10kb) %>%
  mutate(
    assigned_tss_coordinate = coalesce(epd_tss_start, annotation_start),
    promoter_site_key = str_c(gene_id, annotation_chr, assigned_tss_coordinate, sep = "|")
  ) %>%
  group_by(loop_id, anchor_side) %>%
  summarise(
    n_secondary_promoter_sites = n_distinct(promoter_site_key),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = anchor_side,
    values_from = n_secondary_promoter_sites,
    names_glue = "n_secondary_promoter_sites_{anchor_side}",
    values_fill = 0L
  )

df.exploratory.promoter.site.count.by.anchor <- df.proximal.promoter.tss.anchor.pair %>%
  mutate(
    assigned_tss_coordinate = coalesce(epd_tss_start, annotation_start),
    promoter_site_key = str_c(gene_id, annotation_chr, assigned_tss_coordinate, sep = "|")
  ) %>%
  group_by(loop_id, anchor_side) %>%
  summarise(
    n_exploratory_promoter_sites = n_distinct(promoter_site_key),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = anchor_side,
    values_from = n_exploratory_promoter_sites,
    names_glue = "n_exploratory_promoter_sites_{anchor_side}",
    values_fill = 0L
  )

# Report the strict exact-overlap and primary +/-1-kb definitions independently.
# These rows intentionally overlap because strict exact support is nested within
# the primary promoter-window definition.
df.promoter.tss.anchor.evidence.definition.summary <- bind_rows(
  df.strict.direct.promoter.tss.loop.summary %>%
    summarise(
      evidence_definition = "strict_exact_direct_overlap",
      description = paste0(
        "A one-base Ensembl TSS or original 81-bp EPD promoter interval ",
        "overlaps a loop anchor."
      ),
      n_loops = sum(has_any_direct_promoter_tss),
      n_one_anchor_loops = sum(n_direct_anchor_sides == 1L),
      n_both_anchor_loops = sum(n_direct_anchor_sides == 2L)
    ),
  df.direct.promoter.tss.loop.summary %>%
    summarise(
      evidence_definition = "primary_TSS_plus_minus_1kb_direct_overlap",
      description = paste0(
        "An Ensembl or EPD TSS +/-1-kb promoter window overlaps a loop ",
        "anchor; strict exact-overlap loops are included."
      ),
      n_loops = sum(has_any_direct_promoter_tss),
      n_one_anchor_loops = sum(n_direct_anchor_sides == 1L),
      n_both_anchor_loops = sum(n_direct_anchor_sides == 2L)
    )
)

# Assign each pooled loop to exactly one promoter/TSS evidence category. Strict
# exact evidence takes precedence over promoter-window-only evidence, followed
# by midpoint-resolved secondary proximity and the exploratory 200-kb catalog.
# Only single-anchor primary-direct categories are eligible to enter the final
# putative P-E evaluation; opposite-anchor non-TSS ATAC support is still required.
df.promoter.tss.exclusive.loop.category <- df.loop.universe %>%
  left_join(df.strict.direct.promoter.tss.loop.summary %>% transmute(loop_id, n_strict_exact_anchor_sides = n_direct_anchor_sides, has_strict_exact_direct_overlap = has_any_direct_promoter_tss), by = "loop_id") %>%
  left_join(df.direct.promoter.tss.loop.summary %>% transmute(loop_id, n_primary_1kb_anchor_sides = n_direct_anchor_sides, has_primary_1kb_direct_overlap = has_any_direct_promoter_tss), by = "loop_id") %>%
  left_join(df.proximal.promoter.tss.loop.summary %>% dplyr::select(loop_id, n_proximal_anchor_sides, has_any_proximal_promoter_tss, n_secondary_inward_proximal_anchor_sides_10kb, has_any_secondary_inward_proximal_10kb), by = "loop_id") %>%
  left_join(df.primary.direct.promoter.site.count.by.anchor, by = "loop_id") %>%
  left_join(df.secondary.promoter.site.count.by.anchor, by = "loop_id") %>%
  left_join(df.exploratory.promoter.site.count.by.anchor, by = "loop_id") %>%
  mutate(
    across(matches("^n_.*promoter_sites_anchor[12]$"), ~ coalesce(.x, 0L)),
    mutually_exclusive_category = case_when(
      has_strict_exact_direct_overlap & n_primary_1kb_anchor_sides == 1L ~ "01_strict_exact_direct_primary_single_anchor",
      has_strict_exact_direct_overlap & n_primary_1kb_anchor_sides == 2L ~ "02_strict_exact_direct_primary_both_anchors",
      !has_strict_exact_direct_overlap & n_primary_1kb_anchor_sides == 1L ~ "03_primary_1kb_window_only_single_anchor",
      !has_strict_exact_direct_overlap & n_primary_1kb_anchor_sides == 2L ~ "04_primary_1kb_window_only_both_anchors",
      !has_primary_1kb_direct_overlap & n_secondary_inward_proximal_anchor_sides_10kb == 1L ~ "05_secondary_midpoint_nearest_10kb_single_anchor",
      !has_primary_1kb_direct_overlap & n_secondary_inward_proximal_anchor_sides_10kb == 2L ~ "06_secondary_midpoint_nearest_10kb_both_anchors",
      !has_primary_1kb_direct_overlap & !has_any_secondary_inward_proximal_10kb & has_any_proximal_promoter_tss ~ "07_exploratory_other_proximity_1bp_to_200kb",
      TRUE ~ "08_no_promoter_TSS_anchor_or_proximity_support"
    ),
    category_description = case_when(
      mutually_exclusive_category == "01_strict_exact_direct_primary_single_anchor" ~ "Strict exact TSS/EPD overlap is present and the primary +/-1-kb definition supports exactly one anchor.",
      mutually_exclusive_category == "02_strict_exact_direct_primary_both_anchors" ~ "Strict exact TSS/EPD overlap is present and the primary +/-1-kb definition supports both anchors; promoter-promoter compatible.",
      mutually_exclusive_category == "03_primary_1kb_window_only_single_anchor" ~ "No strict exact overlap; a TSS +/-1-kb promoter window supports exactly one anchor.",
      mutually_exclusive_category == "04_primary_1kb_window_only_both_anchors" ~ "No strict exact overlap; TSS +/-1-kb promoter windows support both anchors; promoter-promoter compatible.",
      mutually_exclusive_category == "05_secondary_midpoint_nearest_10kb_single_anchor" ~ "No primary direct overlap; a non-overlapping TSS/EPD promoter is within 10 kb, on the assigned midpoint side, and strictly closer to one anchor.",
      mutually_exclusive_category == "06_secondary_midpoint_nearest_10kb_both_anchors" ~ "No primary direct overlap; distinct midpoint-resolved nearest TSS/EPD candidates support both anchors within 10 kb.",
      mutually_exclusive_category == "07_exploratory_other_proximity_1bp_to_200kb" ~ "No primary or resolved secondary support; other non-overlapping TSS/EPD proximity is present within 200 kb.",
      TRUE ~ "No direct promoter/TSS anchor overlap and no non-overlapping TSS/EPD promoter within the exploratory 200-kb search range."
    ),
    n_category_promoter_sites_anchor1 = case_when(
      str_detect(mutually_exclusive_category, "^(01|02|03|04)_") ~ n_primary_direct_promoter_sites_anchor1,
      str_detect(mutually_exclusive_category, "^(05|06)_") ~ n_secondary_promoter_sites_anchor1,
      mutually_exclusive_category == "07_exploratory_other_proximity_1bp_to_200kb" ~ n_exploratory_promoter_sites_anchor1,
      TRUE ~ 0L
    ),
    n_category_promoter_sites_anchor2 = case_when(
      str_detect(mutually_exclusive_category, "^(01|02|03|04)_") ~ n_primary_direct_promoter_sites_anchor2,
      str_detect(mutually_exclusive_category, "^(05|06)_") ~ n_secondary_promoter_sites_anchor2,
      mutually_exclusive_category == "07_exploratory_other_proximity_1bp_to_200kb" ~ n_exploratory_promoter_sites_anchor2,
      TRUE ~ 0L
    ),
    anchor_promoter_site_multiplicity = case_when(
      pmax(n_category_promoter_sites_anchor1, n_category_promoter_sites_anchor2) == 0L ~ "no_assigned_promoter_site",
      pmax(n_category_promoter_sites_anchor1, n_category_promoter_sites_anchor2) == 1L ~ "one_promoter_site_at_each_supported_anchor",
      TRUE ~ "multiple_promoter_sites_at_one_or_more_anchors"
    ),
    final_P_E_input_OK = if_else(
      mutually_exclusive_category %in% c("01_strict_exact_direct_primary_single_anchor", "03_primary_1kb_window_only_single_anchor"),
      "OK",
      ""
    )
  ) %>%
  arrange(mutually_exclusive_category, loop_id)

df.promoter.tss.exclusive.loop.category.summary <- df.promoter.tss.exclusive.loop.category %>%
  group_by(mutually_exclusive_category, category_description, anchor_promoter_site_multiplicity, final_P_E_input_OK) %>%
  summarise(
    n_loops = n(),
    pct_pooled_loops = round(100 * n() / nrow(df.loop.universe), 1),
    .groups = "drop"
  ) %>%
  dplyr::select(mutually_exclusive_category, category_description, anchor_promoter_site_multiplicity, n_loops, pct_pooled_loops, final_P_E_input_OK) %>%
  arrange(mutually_exclusive_category)

# Validate complete, nested, and mutually exclusive promoter/TSS categories.
assert_analysis_condition(
  nrow(df.promoter.tss.exclusive.loop.category) == nrow(df.loop.universe) &&
    sum(df.promoter.tss.exclusive.loop.category.summary$n_loops) == nrow(df.loop.universe) &&
    !any(df.promoter.tss.exclusive.loop.category$has_strict_exact_direct_overlap & !df.promoter.tss.exclusive.loop.category$has_primary_1kb_direct_overlap),
  "The mutually exclusive promoter/TSS categories failed validation."
)

# Preserve dual-primary-anchor loops as promoter-promoter-compatible contacts.
# For gene-centric sensitivity analyses only, choose a representative promoter
# anchor by prioritizing strict exact overlap and then the shortest TSS-to-anchor
# midpoint distance. This does not turn the opposite promoter-associated anchor
# into an enhancer and does not alter the final putative P-E category.
# Preserve dual-primary-anchor loops as promoter-promoter-compatible contacts.
# For gene-centric sensitivity analyses only, choose a representative promoter
# anchor by prioritizing strict exact overlap and then the shortest TSS-to-anchor
# midpoint distance. This does not turn the opposite promoter-associated anchor
# into an enhancer and does not alter the final putative P-E category.
df.dual.primary.anchor.representative.score <- df.direct.promoter.tss.anchor.overlap %>%
  mutate(
    assigned_tss_coordinate = coalesce(tss_start, epd_tss_start),
    anchor_midpoint = (anchor_start + anchor_end) / 2,
    tss_to_anchor_midpoint_distance_bp = abs(assigned_tss_coordinate - anchor_midpoint),
    promoter_site_key = str_c(gene_id, annotation_chr, assigned_tss_coordinate, sep = "|")
  ) %>%
  group_by(loop_id, resolution, anchor_side) %>%
  summarise(
    n_unique_promoter_sites = n_distinct(promoter_site_key),
    n_unique_genes = n_distinct(gene_id),
    min_tss_to_anchor_midpoint_distance_bp = min(tss_to_anchor_midpoint_distance_bp),
    closest_promoter_site_keys = str_c(sort(unique(promoter_site_key[tss_to_anchor_midpoint_distance_bp == min(tss_to_anchor_midpoint_distance_bp)])), collapse = ";"),
    .groups = "drop"
  ) %>%
  left_join(
    df.strict.direct.promoter.tss.anchor.summary %>% transmute(loop_id, anchor_side, has_strict_exact_anchor_support = has_any_direct_promoter_tss),
    by = c("loop_id", "anchor_side")
  ) %>%
  inner_join(
    df.direct.promoter.tss.loop.summary %>% filter(n_direct_anchor_sides == 2L) %>% dplyr::select(loop_id),
    by = "loop_id"
  ) %>%
  pivot_wider(
    names_from = anchor_side,
    values_from = c(has_strict_exact_anchor_support, n_unique_promoter_sites, n_unique_genes, min_tss_to_anchor_midpoint_distance_bp, closest_promoter_site_keys),
    names_glue = "{.value}_{anchor_side}"
  ) %>%
  left_join(df.promoter.tss.exclusive.loop.category %>% dplyr::select(loop_id, mutually_exclusive_category), by = "loop_id") %>%
  mutate(
    representative_promoter_anchor = case_when(
      has_strict_exact_anchor_support_anchor1 & !has_strict_exact_anchor_support_anchor2 ~ "anchor1",
      !has_strict_exact_anchor_support_anchor1 & has_strict_exact_anchor_support_anchor2 ~ "anchor2",
      min_tss_to_anchor_midpoint_distance_bp_anchor1 < min_tss_to_anchor_midpoint_distance_bp_anchor2 ~ "anchor1",
      min_tss_to_anchor_midpoint_distance_bp_anchor2 < min_tss_to_anchor_midpoint_distance_bp_anchor1 ~ "anchor2",
      TRUE ~ NA_character_
    ),
    representative_anchor_selection_reason = case_when(
      has_strict_exact_anchor_support_anchor1 & !has_strict_exact_anchor_support_anchor2 ~ "anchor1_has_strict_exact_support_only",
      !has_strict_exact_anchor_support_anchor1 & has_strict_exact_anchor_support_anchor2 ~ "anchor2_has_strict_exact_support_only",
      !is.na(representative_promoter_anchor) ~ "shorter_TSS_to_anchor_midpoint_distance",
      TRUE ~ "unresolved_equal_evidence_and_distance"
    ),
    representative_distance_advantage_bp = abs(min_tss_to_anchor_midpoint_distance_bp_anchor1 - min_tss_to_anchor_midpoint_distance_bp_anchor2),
    representative_anchor_selection_status = if_else(is.na(representative_promoter_anchor), "unresolved_tie", "selected_for_gene_centric_sensitivity_only"),
    opposite_anchor_remains_promoter_associated = TRUE,
    eligible_for_final_P_E_reclassification = "NO"
  ) %>%
  arrange(loop_id)

df.dual.primary.anchor.representative.summary <- df.dual.primary.anchor.representative.score %>%
  count(representative_anchor_selection_reason, representative_anchor_selection_status, name = "n_loops") %>%
  mutate(pct_dual_primary_loops = round(100 * n_loops / nrow(df.dual.primary.anchor.representative.score), 1)) %>%
  arrange(desc(n_loops))

df.dual.primary.representative.gene.assignment <- df.direct.promoter.tss.gene.assignment %>%
  inner_join(
    df.dual.primary.anchor.representative.score %>%
      filter(!is.na(representative_promoter_anchor)) %>%
      transmute(loop_id, anchor_side = representative_promoter_anchor, representative_anchor_selection_reason, representative_distance_advantage_bp),
    by = c("loop_id", "anchor_side")
  ) %>%
  distinct(loop_id, anchor_side, gene_id, .keep_all = TRUE) %>%
  arrange(loop_id, anchor_side, gene_id)

# Confirm that representative-anchor scoring retains every dual-primary loop.
assert_analysis_row_count(
  df.dual.primary.anchor.representative.score,
  sum(df.direct.promoter.tss.loop.summary$n_direct_anchor_sides == 2L),
  "Dual-primary representative-anchor table lost one or more loops."
)

# Retain a compact loop-level view of the secondary inward <=10-kb tier while
# preserving all pooled loops and the primary-direct relationship for auditing.
df.secondary.inward.proximal.loop.summary <- df.proximal.promoter.tss.loop.summary %>%
  dplyr::select(
    loop_id, resolution, chr1, start1, end1, chr2, start2, end2, n_direct_anchor_sides,
    has_any_direct_promoter_tss, n_secondary_inward_proximal_anchor_sides_10kb, has_any_secondary_inward_proximal_10kb,
    has_secondary_inward_proximal_both_anchors_10kb, n_secondary_inward_proximal_genes_10kb_across_anchors, promoter_tss_assignment_tier
  )

df.secondary.inward.proximal.by.resolution <- bind_rows(
  df.secondary.inward.proximal.gene.assignment,
  df.secondary.inward.proximal.gene.assignment %>% mutate(resolution = "ALL")
) %>%
  group_by(resolution) %>%
  summarise(
    n_loop_anchor_gene_assignments = n(),
    n_loop_anchor_pairs = n_distinct(loop_id, anchor_side),
    n_loops = n_distinct(loop_id),
    n_genes = n_distinct(gene_id),
    median_midpoint_nearest_inward_edge_distance_bp = median(min_midpoint_nearest_inward_proximal_distance_bp),
    .groups = "drop"
  )

# Verify that the proximal loop summary preserves the full pooled universe.
assert_analysis_condition(
  nrow(df.proximal.promoter.tss.loop.summary) == nrow(df.loop.universe) &&
    !any(is.na(df.proximal.promoter.tss.loop.summary$n_proximal_anchor_sides)),
  "Proximal loop summary failed to preserve the pooled loop universe."
)

df.proximal.promoter.tss.by.resolution.distance <- df.proximal.promoter.tss.anchor.pair %>%
  group_by(resolution, proximal_distance_tier, annotation_class) %>%
  summarise(
    n_annotation_records = n(),
    n_loop_anchor_pairs = n_distinct(loop_id, anchor_side),
    n_loops = n_distinct(loop_id),
    n_genes = n_distinct(gene_id),
    .groups = "drop"
  ) %>%
  bind_rows(
    df.proximal.promoter.tss.anchor.pair %>%
      group_by(proximal_distance_tier, annotation_class) %>%
      summarise(
        n_annotation_records = n(),
        n_loop_anchor_pairs = n_distinct(loop_id, anchor_side),
        n_loops = n_distinct(loop_id),
        n_genes = n_distinct(gene_id),
        .groups = "drop"
      ) %>%
      mutate(resolution = "ALL", .before = 1)
  ) %>%
  arrange(resolution, proximal_distance_tier, annotation_class)

# Quantify how quickly proximity becomes non-specific as the window expands.
# These thresholds are sensitivity summaries, not validated biological cutoffs.
proximal.distance.thresholds.bp <- c(1000L, 5000L, 10000L, 25000L, 50000L, 100000L, 200000L)
direct.supported.loop.ids <- df.direct.promoter.tss.loop.summary %>% filter(has_any_direct_promoter_tss) %>% pull(loop_id)

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
    pct_pooled_loops_with_proximity = round(100 * length(proximal.loop.ids) / nrow(df.loop.universe), 1),
    n_direct_and_proximity_loops = length(intersect(proximal.loop.ids, direct.supported.loop.ids)),
    n_proximity_without_direct_loops = length(setdiff(proximal.loop.ids, direct.supported.loop.ids)),
    n_direct_without_proximity_loops = length(setdiff(direct.supported.loop.ids, proximal.loop.ids)),
    n_without_direct_or_proximity_loops = nrow(df.loop.universe) - length(union(proximal.loop.ids, direct.supported.loop.ids))
  )
})

df.proximal.promoter.tss.analysis.definition <- tribble(
  ~analysis_item, ~definition,
  "direct_vs_proximal", paste0(
    "Direct evidence overlaps an anchor; proximal evidence does not overlap ",
    "and is 1-200 kb from the nearest anchor boundary."
  ),
  "distance_measure", paste0(
    "Minimum coordinate separation between the annotation and anchor edge; ",
    "an immediately adjacent annotation is 1 bp away."
  ),
  "inward_direction", paste0(
    "An inward annotation is fully contained between the end of anchor 1 and ",
    "the start of anchor 2: right of anchor 1 or left of anchor 2."
  ),
  "midpoint_direction", paste0(
    "A secondary anchor-1 candidate must lie fully on the anchor-1 side of ",
    "the inter-anchor midpoint; an anchor-2 candidate must lie fully on the ",
    "anchor-2 side. Midpoint-spanning annotations are retained only as ",
    "ambiguous exploratory evidence."
  ),
  "nearest_anchor_requirement", paste0(
    "The candidate must be strictly closer by edge distance to its assigned ",
    "anchor than to the opposite anchor. Equal-distance candidates are not ",
    "eligible for the secondary assignment tier."
  ),
  "secondary_inward_10kb", paste0(
    "Secondary evidence requires inward location, the correct midpoint side, ",
    "strictly shorter distance to the assigned than opposite anchor, edge ",
    "distance <=10 kb, and absence of a primary direct assignment for the ",
    "same loop-anchor-gene."
  ),
  "full_200kb_catalog", paste0(
    "Exploratory sensitivity catalog only; it must not be interpreted as a ",
    "validated target-gene or promoter-enhancer assignment."
  ),
  "nearest_flags", paste0(
    "Nearest-gene ranks remain descriptive within each anchor. Midpoint or ",
    "opposite-anchor distance ambiguities are retained in the exploratory ",
    "catalog but excluded from the secondary assignment tier."
  ),
  "downstream_evidence", paste0(
    "ATAC support and transcript-containment flags are added in Sections 5 ",
    "and 6; containment remains descriptive rather than a filter."
  )
)

df.proximal.promoter.tss.summary <- tibble(
  metric = c(
    "pooled_loops",
    "proximal_true_TSS_records_1bp_to_200kb",
    "proximal_EPD_promoter_records_1bp_to_200kb",
    "unique_proximal_loop_anchor_gene_assignments",
    "loops_with_any_proximal_promoter_or_TSS",
    "loops_with_proximal_promoter_or_TSS_at_one_anchor",
    "loops_with_proximal_promoter_or_TSS_at_both_anchors",
    "secondary_inward_10kb_loop_anchor_gene_assignments",
    "loops_with_any_secondary_inward_10kb_assignment",
    "loops_with_primary_direct_and_secondary_inward_10kb",
    "loops_with_primary_direct_only",
    "loops_with_secondary_inward_10kb_only",
    "loops_with_exploratory_200kb_proximity_only",
    "loops_without_direct_or_proximal_evidence"
  ),
  n = c(
    nrow(df.loop.universe),
    nrow(df.proximal.true.tss.anchor.pair),
    nrow(df.proximal.epd.promoter.anchor.pair),
    nrow(df.proximal.promoter.tss.gene.assignment),
    sum(df.proximal.promoter.tss.loop.summary$has_any_proximal_promoter_tss),
    sum(df.proximal.promoter.tss.loop.summary$n_proximal_anchor_sides == 1L),
    sum(df.proximal.promoter.tss.loop.summary$has_proximal_promoter_tss_both_anchors),
    nrow(df.secondary.inward.proximal.gene.assignment),
    sum(df.proximal.promoter.tss.loop.summary$has_any_secondary_inward_proximal_10kb),
    sum(df.proximal.promoter.tss.loop.summary$promoter_tss_assignment_tier == "primary_direct_and_secondary_inward_10kb"),
    sum(df.proximal.promoter.tss.loop.summary$promoter_tss_assignment_tier == "direct_only"),
    sum(df.proximal.promoter.tss.loop.summary$promoter_tss_assignment_tier == "secondary_inward_10kb_only"),
    sum(df.proximal.promoter.tss.loop.summary$promoter_tss_assignment_tier == "exploratory_proximal_200kb_only"),
    sum(df.proximal.promoter.tss.loop.summary$promoter_tss_assignment_tier == "no_direct_or_proximal_evidence")
  )
)

message(
  "Proximal promoter/TSS analysis retained ",
  nrow(df.proximal.promoter.tss.gene.assignment),
  " loop-anchor-gene assignments across ",
  sum(df.proximal.promoter.tss.loop.summary$has_any_proximal_promoter_tss),
  " pooled loops."
)
message(
  "Secondary inward <=10-kb tier retained ",
  nrow(df.secondary.inward.proximal.gene.assignment),
  " loop-anchor-gene assignments across ",
  sum(df.proximal.promoter.tss.loop.summary$has_any_secondary_inward_proximal_10kb),
  " pooled loops."
)
################################################################################
# 5. Recalculate ATAC support using true TSS exclusion regions
#
# ATAC-seq overlap is interpreted only as open-chromatin support. It does not
# validate enhancer activity or a functional promoter-enhancer interaction.
# Ensembl transcript TSS and EPD TSS positions are combined, expanded by +/-1 kb,
# and removed from the ATAC union before non-TSS accessibility is evaluated.
#
# The complete two-anchor evidence table is retained for later categories. The
# primary directional summary is restricted to loops with direct promoter/TSS
# evidence at exactly one anchor; loops with direct evidence at both anchors are
# preserved separately rather than forcing one side to be an enhancer candidate.
################################################################################

atac.minimum.overlap.bp <- 50L
tss.exclusion.flank.bp <- 1000L

# Combine every transcript-level Ensembl TSS with every EPD TSS. Source records
# are retained long enough to audit the input counts, then identical genomic TSS
# positions are collapsed before construction of the exclusion intervals.
df.known.tss.source.record <- bind_rows(
  df.true.tss.transcript %>%
    transmute(
      tss_source = "Ensembl_transcript_true_TSS",
      tss_source_id = true_tss_id,
      chr,
      tss_start = true_tss_start,
      tss_end = true_tss_end
    ),
  df.promoter.epd.rn7.1based %>%
    transmute(
      tss_source = "EPDnew_TSS",
      tss_source_id = promoter_annotation_id,
      chr,
      tss_start = epd_tss_start,
      tss_end = epd_tss_end
    )
)

# Confirm that every known-TSS exclusion source is represented by one base.
assert_analysis_condition(
  !any(
    df.known.tss.source.record$tss_start !=
      df.known.tss.source.record$tss_end
  ),
  "A known TSS record is not a one-base interval."
)

df.known.tss.point <- df.known.tss.source.record %>%
  group_by(chr, tss_start, tss_end) %>%
  summarise(
    n_tss_source_records = n(),
    tss_sources = str_c(sort(unique(tss_source)), collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(chr, tss_start)

gr.known.tss.point <- GRanges(
  seqnames = df.known.tss.point$chr,
  ranges = IRanges(
    start = df.known.tss.point$tss_start,
    end = df.known.tss.point$tss_end
  ),
  n_tss_source_records = df.known.tss.point$n_tss_source_records,
  tss_sources = df.known.tss.point$tss_sources
)

# Use a symmetric 2,001-bp interval centered on each one-base TSS. This avoids
# strand-dependent promoter() width conventions because only genomic proximity
# to a known TSS, not transcriptional direction, is being excluded from ATAC.
gr.known.tss.exclusion <- GRanges(
  seqnames = seqnames(gr.known.tss.point),
  ranges = IRanges(
    start = pmax(1L, start(gr.known.tss.point) - tss.exclusion.flank.bp),
    end = end(gr.known.tss.point) + tss.exclusion.flank.bp
  )
) %>%
  GenomicRanges::reduce(ignore.strand = TRUE)

# Reduce ATAC peaks before measuring covered bases so overlapping peaks are not
# counted twice. Subtract the true-TSS exclusion union from this reduced signal.
gr.atac.union.revised <- GenomicRanges::reduce(
  gr.atac.rn7.1based,
  ignore.strand = TRUE
)
common.seqlevels.revised.atac.tss <- intersect(
  seqlevels(gr.atac.union.revised),
  seqlevels(gr.known.tss.exclusion)
)
gr.atac.union.revised.common <- keepSeqlevels(
  gr.atac.union.revised,
  common.seqlevels.revised.atac.tss,
  pruning.mode = "coarse"
)
gr.known.tss.exclusion.for.atac <- keepSeqlevels(
  gr.known.tss.exclusion,
  common.seqlevels.revised.atac.tss,
  pruning.mode = "coarse"
)
gr.atac.non.tss.revised <- GenomicRanges::setdiff(
  gr.atac.union.revised.common,
  gr.known.tss.exclusion.for.atac,
  ignore.strand = TRUE
)

df.revised.atac.true.tss.exclusion.region <- tibble(
  chr = as.character(seqnames(gr.known.tss.exclusion)),
  start = start(gr.known.tss.exclusion),
  end = end(gr.known.tss.exclusion),
  width_bp = width(gr.known.tss.exclusion),
  exclusion_definition = "combined_Ensembl_and_EPD_true_TSS_plus_minus_1kb"
)

df.revised.atac.exclusion.summary <- tibble(
  metric = c(
    "input_ATAC_peak_records",
    "reduced_raw_ATAC_intervals",
    "raw_ATAC_union_bp",
    "Ensembl_transcript_true_TSS_records",
    "EPD_true_TSS_records",
    "unique_combined_known_TSS_positions",
    "reduced_true_TSS_plus_minus_1kb_exclusion_intervals",
    "true_TSS_exclusion_union_bp",
    "non_TSS_ATAC_intervals",
    "non_TSS_ATAC_union_bp",
    "ATAC_bp_removed_by_true_TSS_exclusion"
  ),
  n = c(
    length(gr.atac.rn7.1based),
    length(gr.atac.union.revised.common),
    sum(width(gr.atac.union.revised.common)),
    nrow(df.true.tss.transcript),
    nrow(df.promoter.epd.rn7.1based),
    nrow(df.known.tss.point),
    length(gr.known.tss.exclusion),
    sum(width(gr.known.tss.exclusion)),
    length(gr.atac.non.tss.revised),
    sum(width(gr.atac.non.tss.revised)),
    sum(width(gr.atac.union.revised.common)) -
      sum(width(gr.atac.non.tss.revised))
  )
)

# Measure ATAC support at each directly assigned TSS +/-1-kb promoter window.
# This differs from raw whole-anchor ATAC: a promoter is called accessible only
# when at least one assigned promoter window itself overlaps ATAC by >=50 bp.
# Measure ATAC support at each directly assigned TSS +/-1-kb promoter window.
# This differs from raw whole-anchor ATAC: a promoter is called accessible only
# when at least one assigned promoter window itself overlaps ATAC by >=50 bp.
df.revised.promoter.window.atac.record <- df.direct.promoter.tss.anchor.overlap %>%
  distinct(loop_id, resolution, anchor_side, annotation_chr, annotation_start, annotation_end) %>%
  arrange(loop_id, anchor_side, annotation_chr, annotation_start, annotation_end)

gr.revised.promoter.window.atac.record <- GRanges(
  seqnames = df.revised.promoter.window.atac.record$annotation_chr,
  ranges = IRanges(start = df.revised.promoter.window.atac.record$annotation_start, end = df.revised.promoter.window.atac.record$annotation_end),
  loop_id = df.revised.promoter.window.atac.record$loop_id,
  resolution = df.revised.promoter.window.atac.record$resolution,
  anchor_side = df.revised.promoter.window.atac.record$anchor_side
)

df.revised.promoter.window.atac.overlap <- summarise_anchor_interval_overlap(
  gr.revised.promoter.window.atac.record,
  gr.atac.union.revised.common,
  minimum.overlap.bp = atac.minimum.overlap.bp
)

df.revised.promoter.window.atac.anchor.summary <- df.revised.promoter.window.atac.overlap %>%
  group_by(loop_id, resolution, anchor_side) %>%
  summarise(
    n_direct_promoter_windows = n(),
    n_atac_positive_promoter_windows = sum(has_feature_overlap_ge_minimum),
    promoter_window_atac_overlap_any = any(has_feature_overlap_any),
    promoter_window_atac_overlap_ge50 = any(has_feature_overlap_ge_minimum),
    promoter_window_atac_max_overlap_bp = max(feature_overlap_bp),
    .groups = "drop"
  ) %>%
  arrange(loop_id, anchor_side)

# Confirm that promoter-window ATAC evidence covers every directly annotated
# anchor once and does not create records for anchors lacking promoter evidence.
assert_analysis_condition(
  nrow(df.revised.promoter.window.atac.anchor.summary) == sum(df.direct.promoter.tss.anchor.summary$has_any_direct_promoter_tss) &&
    !anyDuplicated(str_c(df.revised.promoter.window.atac.anchor.summary$loop_id, df.revised.promoter.window.atac.anchor.summary$anchor_side, sep = "|")),
  "Promoter-window ATAC summary is not one row per direct promoter anchor."
)

# Calculate raw ATAC, TSS-exclusion, and residual non-TSS ATAC evidence once for
# all 62,042 loop anchors. Later directional and category tables join these data
# rather than rerunning overlap logic on selected subsets.
gr.loop.anchor.all <- combine_loop_anchor_granges_by_side(list.gr.loop.anchor.by.side)

df.revised.raw.atac.anchor.overlap <- summarise_anchor_interval_overlap(
  gr.loop.anchor.all,
  gr.atac.union.revised.common,
  minimum.overlap.bp = atac.minimum.overlap.bp
)
df.revised.non.tss.atac.anchor.overlap <- summarise_anchor_interval_overlap(
  gr.loop.anchor.all,
  gr.atac.non.tss.revised,
  minimum.overlap.bp = atac.minimum.overlap.bp
)
df.revised.tss.exclusion.anchor.overlap <- summarise_anchor_interval_overlap(
  gr.loop.anchor.all,
  gr.known.tss.exclusion,
  minimum.overlap.bp = atac.minimum.overlap.bp
)

df.revised.atac.anchor.evidence <- df.revised.raw.atac.anchor.overlap %>%
  transmute(
    anchor_index,
    loop_id,
    resolution,
    anchor_side,
    anchor_chr,
    anchor_start,
    anchor_end,
    anchor_width_bp,
    raw_atac_n_intervals_any = n_overlapping_feature_intervals,
    raw_atac_n_intervals_ge50 = n_overlapping_feature_intervals_ge_minimum,
    raw_atac_overlap_bp = feature_overlap_bp,
    raw_atac_overlap_fraction_of_anchor = feature_overlap_fraction_of_anchor,
    raw_atac_overlap_any = has_feature_overlap_any,
    raw_atac_overlap_ge50 = has_feature_overlap_ge_minimum,
    raw_atac_total_overlap_ge50 = has_feature_total_overlap_ge_minimum
  ) %>%
  left_join(
    df.revised.non.tss.atac.anchor.overlap %>%
      transmute(
        anchor_index,
        non_tss_atac_n_intervals_any = n_overlapping_feature_intervals,
        non_tss_atac_n_intervals_ge50 = n_overlapping_feature_intervals_ge_minimum,
        non_tss_atac_overlap_bp = feature_overlap_bp,
        non_tss_atac_overlap_fraction_of_anchor = feature_overlap_fraction_of_anchor,
        non_tss_atac_overlap_any = has_feature_overlap_any,
        non_tss_atac_overlap_ge50 = has_feature_overlap_ge_minimum,
        non_tss_atac_total_overlap_ge50 = has_feature_total_overlap_ge_minimum
      ),
    by = "anchor_index"
  ) %>%
  left_join(
    df.revised.tss.exclusion.anchor.overlap %>%
      transmute(
        anchor_index,
        known_tss_exclusion_n_intervals = n_overlapping_feature_intervals,
        known_tss_exclusion_overlap_bp = feature_overlap_bp,
        known_tss_exclusion_overlap_fraction_of_anchor = feature_overlap_fraction_of_anchor,
        overlaps_known_tss_exclusion = has_feature_overlap_any
      ),
    by = "anchor_index"
  ) %>%
  mutate(
    non_tss_anchor_bp = anchor_width_bp - known_tss_exclusion_overlap_bp,
    has_non_tss_anchor_fragment = non_tss_anchor_bp > 0L,
    non_tss_atac_overlap_fraction_of_non_tss_anchor = if_else(has_non_tss_anchor_fragment, non_tss_atac_overlap_bp / non_tss_anchor_bp, NA_real_)
  ) %>%
  left_join(
    df.direct.promoter.tss.anchor.summary %>%
      dplyr::select(loop_id, anchor_side, n_direct_evidence_records, n_direct_genes, n_true_tss_records, n_epd_promoter_records, has_direct_true_tss, has_direct_epd_promoter, has_any_direct_promoter_tss),
    by = c("loop_id", "anchor_side")
  ) %>%
  left_join(
    df.revised.promoter.window.atac.anchor.summary,
    by = c("loop_id", "resolution", "anchor_side")
  ) %>%
  mutate(
    n_direct_promoter_windows = coalesce(n_direct_promoter_windows, 0L),
    n_atac_positive_promoter_windows = coalesce(n_atac_positive_promoter_windows, 0L),
    promoter_window_atac_overlap_any = coalesce(promoter_window_atac_overlap_any, FALSE),
    promoter_window_atac_overlap_ge50 = coalesce(promoter_window_atac_overlap_ge50, FALSE),
    promoter_window_atac_max_overlap_bp = coalesce(promoter_window_atac_max_overlap_bp, 0L)
  ) %>%
  left_join(
    df.proximal.promoter.tss.anchor.summary %>%
      dplyr::select(loop_id, anchor_side, n_proximal_evidence_records, n_proximal_genes, min_proximal_distance_bp, has_any_proximal_promoter_tss),
    by = c("loop_id", "anchor_side")
  ) %>%
  arrange(loop_id, anchor_side)

# Validate two unique ATAC-evidence anchor rows and nonnegative residual bases.
assert_analysis_condition(
  nrow(df.revised.atac.anchor.evidence) == 2L * nrow(df.loop.universe) &&
    !anyDuplicated(str_c(df.revised.atac.anchor.evidence$loop_id, df.revised.atac.anchor.evidence$anchor_side, sep = "|")) &&
    !any(df.revised.atac.anchor.evidence$non_tss_anchor_bp < 0L),
  "Revised ATAC anchor evidence did not preserve two valid anchors per loop."
)

# Build one directional record for every directly annotated anchor. A loop with
# both anchors directly annotated therefore contributes two explicitly ambiguous
# orientations; only one-direct-anchor loops enter the primary candidate summary.
df.revised.atac.promoter.anchor.metric <- df.revised.atac.anchor.evidence %>%
  transmute(
    loop_id,
    promoter_anchor_side = anchor_side,
    promoter_anchor_chr = anchor_chr,
    promoter_anchor_start = anchor_start,
    promoter_anchor_end = anchor_end,
    promoter_anchor_width_bp = anchor_width_bp,
    promoter_raw_atac_overlap_any = raw_atac_overlap_any,
    promoter_raw_atac_overlap_ge50 = raw_atac_overlap_ge50,
    promoter_raw_atac_overlap_bp = raw_atac_overlap_bp,
    promoter_window_atac_overlap_any,
    promoter_window_atac_overlap_ge50,
    promoter_window_atac_max_overlap_bp,
    promoter_n_atac_positive_windows = n_atac_positive_promoter_windows,
    promoter_non_tss_atac_overlap_any = non_tss_atac_overlap_any,
    promoter_non_tss_atac_overlap_ge50 = non_tss_atac_overlap_ge50,
    promoter_non_tss_atac_overlap_bp = non_tss_atac_overlap_bp,
    promoter_known_tss_exclusion_overlap_bp = known_tss_exclusion_overlap_bp,
    promoter_non_tss_anchor_bp = non_tss_anchor_bp
  )

df.revised.atac.candidate.anchor.metric <- df.revised.atac.anchor.evidence %>%
  transmute(
    loop_id,
    candidate_anchor_side = anchor_side,
    candidate_anchor_chr = anchor_chr,
    candidate_anchor_start = anchor_start,
    candidate_anchor_end = anchor_end,
    candidate_anchor_width_bp = anchor_width_bp,
    candidate_anchor_has_direct_promoter_tss = has_any_direct_promoter_tss,
    candidate_promoter_window_atac_overlap_any = promoter_window_atac_overlap_any,
    candidate_promoter_window_atac_overlap_ge50 = promoter_window_atac_overlap_ge50,
    candidate_promoter_window_atac_max_overlap_bp = promoter_window_atac_max_overlap_bp,
    candidate_raw_atac_overlap_any = raw_atac_overlap_any,
    candidate_raw_atac_overlap_ge50 = raw_atac_overlap_ge50,
    candidate_raw_atac_overlap_bp = raw_atac_overlap_bp,
    candidate_non_tss_atac_overlap_any = non_tss_atac_overlap_any,
    candidate_non_tss_atac_overlap_ge50 = non_tss_atac_overlap_ge50,
    candidate_non_tss_atac_total_overlap_ge50 = non_tss_atac_total_overlap_ge50,
    candidate_non_tss_atac_overlap_bp = non_tss_atac_overlap_bp,
    candidate_known_tss_exclusion_overlap_bp = known_tss_exclusion_overlap_bp,
    candidate_non_tss_anchor_bp = non_tss_anchor_bp,
    candidate_has_non_tss_anchor_fragment = has_non_tss_anchor_fragment
  )

df.revised.atac.direct.orientation <- df.direct.promoter.tss.anchor.summary %>%
  filter(has_any_direct_promoter_tss) %>%
  transmute(
    loop_id,
    resolution,
    promoter_anchor_side = anchor_side,
    candidate_anchor_side = opposite_anchor_side,
    promoter_n_direct_evidence_records = n_direct_evidence_records,
    promoter_n_direct_genes = n_direct_genes,
    promoter_has_direct_true_tss = has_direct_true_tss,
    promoter_has_direct_epd_promoter = has_direct_epd_promoter
  ) %>%
  left_join(
    df.direct.promoter.tss.loop.summary %>% dplyr::select(loop_id, n_direct_anchor_sides, direct_anchor_assignment_class),
    by = "loop_id"
  ) %>%
  left_join(df.revised.atac.promoter.anchor.metric, by = c("loop_id", "promoter_anchor_side")) %>%
  left_join(df.revised.atac.candidate.anchor.metric, by = c("loop_id", "candidate_anchor_side")) %>%
  mutate(
    is_unambiguous_candidate_regulatory_orientation = n_direct_anchor_sides == 1L,
    directional_mixed_regulatory_support = (promoter_window_atac_overlap_ge50 & candidate_non_tss_atac_overlap_ge50),
    candidate_open_chromatin_support_class = case_when(
      !is_unambiguous_candidate_regulatory_orientation ~ "both_anchors_direct_promoter_TSS_no_unique_candidate",
      candidate_non_tss_atac_overlap_ge50 ~ "opposite_anchor_nonTSS_ATAC_ge50",
      candidate_non_tss_atac_overlap_any ~ "opposite_anchor_nonTSS_ATAC_lt50_only",
      !candidate_has_non_tss_anchor_fragment ~ "opposite_anchor_fully_within_known_TSS_exclusion",
      TRUE ~ "opposite_anchor_no_nonTSS_ATAC"
    )
  ) %>%
  arrange(loop_id, promoter_anchor_side)

# Subclassify loops with direct promoter/TSS evidence at both anchors. A
# direction is supported when its promoter window is ATAC-positive and the
# opposite anchor retains >=50 bp of ATAC after all known TSS windows are
# removed. These classes annotate dual-promoter contacts; they do not redefine
# the opposite promoter-containing anchor as a validated enhancer.
df.dual.promoter.atac.directional.classification <- df.revised.atac.direct.orientation %>%
  filter(n_direct_anchor_sides == 2L) %>%
  group_by(loop_id, resolution) %>%
  summarise(
    promoter_window_atac_ge50_anchor1 = any(promoter_anchor_side == "anchor1" & promoter_window_atac_overlap_ge50),
    promoter_window_atac_ge50_anchor2 = any(promoter_anchor_side == "anchor2" & promoter_window_atac_overlap_ge50),
    residual_non_tss_atac_ge50_anchor1 = any(candidate_anchor_side == "anchor1" & candidate_non_tss_atac_overlap_ge50),
    residual_non_tss_atac_ge50_anchor2 = any(candidate_anchor_side == "anchor2" & candidate_non_tss_atac_overlap_ge50),
    supports_anchor1_promoter_to_anchor2_regulatory = any(promoter_anchor_side == "anchor1" & directional_mixed_regulatory_support),
    supports_anchor2_promoter_to_anchor1_regulatory = any(promoter_anchor_side == "anchor2" & directional_mixed_regulatory_support),
    .groups = "drop"
  ) %>%
  mutate(
    n_accessible_promoter_anchors = as.integer(promoter_window_atac_ge50_anchor1) + as.integer(promoter_window_atac_ge50_anchor2),
    n_supported_regulatory_directions = as.integer(supports_anchor1_promoter_to_anchor2_regulatory) + as.integer(supports_anchor2_promoter_to_anchor1_regulatory),
    dual_promoter_directional_category = case_when(
      n_supported_regulatory_directions == 1L & n_accessible_promoter_anchors == 1L ~ "one_direction_stringent_asymmetric_mixed_regulatory_evidence",
      n_supported_regulatory_directions == 1L & n_accessible_promoter_anchors == 2L ~ "one_direction_broader_mixed_promoter_regulatory_evidence",
      n_supported_regulatory_directions == 2L ~ "both_directions_supported_bidirectional_or_ambiguous",
      TRUE ~ "neither_direction_supported_promoter_promoter_or_lower_support"
    ),
    supported_regulatory_direction = case_when(
      supports_anchor1_promoter_to_anchor2_regulatory & !supports_anchor2_promoter_to_anchor1_regulatory ~ "anchor1_promoter_to_anchor2_regulatory",
      !supports_anchor1_promoter_to_anchor2_regulatory & supports_anchor2_promoter_to_anchor1_regulatory ~ "anchor2_promoter_to_anchor1_regulatory",
      supports_anchor1_promoter_to_anchor2_regulatory & supports_anchor2_promoter_to_anchor1_regulatory ~ "both_directions_supported",
      TRUE ~ "no_supported_direction"
    )
  ) %>%
  arrange(dual_promoter_directional_category, loop_id)

# Record the four mutually exclusive dual-promoter definitions explicitly so
# the reported categories can be audited without reconstructing the case logic.
df.dual.promoter.atac.directional.category.definition <- tribble(
  ~dual_promoter_directional_category, ~promoter_ATAC_definition,
  ~residual_non_TSS_ATAC_definition, ~recommended_interpretation,
  "one_direction_stringent_asymmetric_mixed_regulatory_evidence",
  "Exactly one promoter/TSS window has >=50 bp ATAC overlap",
  "The opposite anchor has >=50 bp residual non-TSS ATAC overlap",
  "One-direction stringent mixed-regulatory evidence",
  "one_direction_broader_mixed_promoter_regulatory_evidence",
  "Both promoter/TSS windows have >=50 bp ATAC overlap",
  "Exactly one anchor has >=50 bp residual non-TSS ATAC overlap",
  "One-direction broader mixed promoter/regulatory evidence",
  "both_directions_supported_bidirectional_or_ambiguous",
  "Both promoter/TSS windows have >=50 bp ATAC overlap",
  "Both anchors have >=50 bp residual non-TSS ATAC overlap",
  "Bidirectional or directionally ambiguous support",
  "neither_direction_supported_promoter_promoter_or_lower_support",
  "Any promoter/TSS ATAC pattern not completing a supported direction",
  "No promoter-to-opposite-residual-ATAC direction is completed",
  "Promoter-promoter-compatible or lower-support contact"
)

df.dual.promoter.atac.directional.summary <- df.dual.promoter.atac.directional.classification %>%
  count(dual_promoter_directional_category, name = "n_loops") %>%
  mutate(pct_dual_promoter_loops = round(100 * n_loops / sum(n_loops), 1)) %>%
  arrange(
    factor(
      dual_promoter_directional_category,
      levels = c(
        "one_direction_stringent_asymmetric_mixed_regulatory_evidence",
        "one_direction_broader_mixed_promoter_regulatory_evidence",
        "both_directions_supported_bidirectional_or_ambiguous",
        "neither_direction_supported_promoter_promoter_or_lower_support"
      )
    )
  )

# Require one mutually exclusive directional category for all dual-promoter
# loops and preserve the expected dual-promoter universe of the current inputs.
assert_analysis_condition(
  nrow(df.dual.promoter.atac.directional.classification) == sum(df.direct.promoter.tss.loop.summary$n_direct_anchor_sides == 2L) &&
    sum(df.dual.promoter.atac.directional.summary$n_loops) == nrow(df.dual.promoter.atac.directional.classification) &&
    !anyDuplicated(df.dual.promoter.atac.directional.classification$loop_id),
  "Dual-promoter ATAC classification did not preserve all eligible loops."
)

df.revised.atac.unambiguous.orientation <- df.revised.atac.direct.orientation %>%
  filter(is_unambiguous_candidate_regulatory_orientation)

# Require one unique ATAC orientation for every single-promoter-anchor loop.
assert_analysis_condition(
  nrow(df.revised.atac.unambiguous.orientation) == sum(df.direct.promoter.tss.loop.summary$n_direct_anchor_sides == 1L) &&
    !anyDuplicated(df.revised.atac.unambiguous.orientation$loop_id),
  "Single-direct-anchor ATAC orientations are not one row per eligible loop."
)

# Reshape anchor metrics to one row per pooled loop. Candidate-regulatory fields
# are populated only when exactly one anchor has direct promoter/TSS evidence.
df.revised.atac.anchor.wide <- df.revised.atac.anchor.evidence %>%
  dplyr::select(
    loop_id, anchor_side, raw_atac_overlap_any, raw_atac_overlap_ge50, raw_atac_overlap_bp,
    non_tss_atac_overlap_any, non_tss_atac_overlap_ge50, non_tss_atac_overlap_bp,
    known_tss_exclusion_overlap_bp, non_tss_anchor_bp
  ) %>%
  pivot_wider(
    names_from = anchor_side,
    values_from = -loop_id,
    names_glue = "{.value}_{anchor_side}"
  )

df.revised.atac.loop.summary <- df.direct.promoter.tss.loop.summary %>%
  dplyr::select(loop_id, resolution, n_direct_anchor_sides, has_any_direct_promoter_tss, has_direct_promoter_tss_both_anchors, direct_anchor_assignment_class) %>%
  left_join(df.revised.atac.anchor.wide, by = "loop_id") %>%
  mutate(
    promoter_anchor_raw_atac_ge50 = case_when(
      direct_anchor_assignment_class == "direct_anchor1_only" ~ raw_atac_overlap_ge50_anchor1,
      direct_anchor_assignment_class == "direct_anchor2_only" ~ raw_atac_overlap_ge50_anchor2,
      TRUE ~ NA
    ),
    candidate_anchor_raw_atac_ge50 = case_when(
      direct_anchor_assignment_class == "direct_anchor1_only" ~ raw_atac_overlap_ge50_anchor2,
      direct_anchor_assignment_class == "direct_anchor2_only" ~ raw_atac_overlap_ge50_anchor1,
      TRUE ~ NA
    ),
    candidate_anchor_non_tss_atac_any = case_when(
      direct_anchor_assignment_class == "direct_anchor1_only" ~ non_tss_atac_overlap_any_anchor2,
      direct_anchor_assignment_class == "direct_anchor2_only" ~ non_tss_atac_overlap_any_anchor1,
      TRUE ~ NA
    ),
    candidate_anchor_non_tss_atac_ge50 = case_when(
      direct_anchor_assignment_class == "direct_anchor1_only" ~ non_tss_atac_overlap_ge50_anchor2,
      direct_anchor_assignment_class == "direct_anchor2_only" ~ non_tss_atac_overlap_ge50_anchor1,
      TRUE ~ NA
    ),
    revised_atac_evidence_class = case_when(
      n_direct_anchor_sides == 2L ~ "direct_promoter_TSS_both_anchors",
      n_direct_anchor_sides == 1L & candidate_anchor_non_tss_atac_ge50 ~ "single_direct_promoter_TSS_opposite_nonTSS_ATAC_ge50",
      n_direct_anchor_sides == 1L & candidate_anchor_non_tss_atac_any ~ "single_direct_promoter_TSS_opposite_nonTSS_ATAC_lt50_only",
      n_direct_anchor_sides == 1L ~ "single_direct_promoter_TSS_no_opposite_nonTSS_ATAC",
      TRUE ~ "no_direct_promoter_TSS_orientation"
    )
  )

# Confirm that the revised ATAC summary retains every pooled loop.
assert_analysis_row_count(df.revised.atac.loop.summary, nrow(df.loop.universe), "Revised ATAC loop summary lost pooled loops.")

# Report raw and non-TSS accessibility among unambiguous direct orientations at
# each HiCCUPS resolution and across all resolutions.
df.revised.atac.support.by.resolution <- bind_rows(
  df.revised.atac.unambiguous.orientation,
  df.revised.atac.unambiguous.orientation %>% mutate(resolution = "ALL")
) %>%
  group_by(resolution) %>%
  summarise(
    n_single_direct_promoter_tss_loops = n(),
    n_promoter_raw_atac_any = sum(promoter_raw_atac_overlap_any),
    pct_promoter_raw_atac_any = round(100 * mean(promoter_raw_atac_overlap_any), 1),
    n_promoter_raw_atac_ge50 = sum(promoter_raw_atac_overlap_ge50),
    pct_promoter_raw_atac_ge50 = round(100 * mean(promoter_raw_atac_overlap_ge50), 1),
    n_opposite_raw_atac_any = sum(candidate_raw_atac_overlap_any),
    pct_opposite_raw_atac_any = round(100 * mean(candidate_raw_atac_overlap_any), 1),
    n_opposite_raw_atac_ge50 = sum(candidate_raw_atac_overlap_ge50),
    pct_opposite_raw_atac_ge50 = round(100 * mean(candidate_raw_atac_overlap_ge50), 1),
    n_opposite_non_tss_atac_any = sum(candidate_non_tss_atac_overlap_any),
    pct_opposite_non_tss_atac_any = round(100 * mean(candidate_non_tss_atac_overlap_any), 1),
    n_opposite_non_tss_atac_ge50 = sum(candidate_non_tss_atac_overlap_ge50),
    pct_opposite_non_tss_atac_ge50 = round(100 * mean(candidate_non_tss_atac_overlap_ge50), 1),
    .groups = "drop"
  ) %>%
  arrange(resolution)

# McNemar tests use paired anchors from the same loop. Raw ATAC applies the same
# peak definition to both sides; the second comparison applies true-TSS-excluded
# ATAC to both sides and is therefore also symmetric.
df.revised.atac.paired.anchor.mcnemar <- map_dfr(
  c(sort(unique(df.revised.atac.unambiguous.orientation$resolution)), "ALL"),
  function(resolution.i) {
    df.resolution <- if (resolution.i == "ALL") {
      df.revised.atac.unambiguous.orientation
    } else {
      df.revised.atac.unambiguous.orientation %>%
        filter(resolution == resolution.i)
    }

    bind_rows(
      summarise_paired_binary_overlap(
        df.resolution,
        promoter.column = "promoter_raw_atac_overlap_ge50",
        candidate.column = "candidate_raw_atac_overlap_ge50",
        resolution.label = resolution.i,
        comparison.label = "raw_ATAC_ge50"
      ),
      summarise_paired_binary_overlap(
        df.resolution,
        promoter.column = "promoter_non_tss_atac_overlap_ge50",
        candidate.column = "candidate_non_tss_atac_overlap_ge50",
        resolution.label = resolution.i,
        comparison.label = "true_TSS_excluded_ATAC_ge50"
      )
    )
  }
) %>%
  arrange(resolution, comparison)

# Repeat candidate-anchor testing after physically subtracting TSS regions from
# each anchor. This fragment-level sensitivity analysis preserves anchor identity
# and can be checked against the primary whole-anchor/non-TSS-signal method.
df.revised.atac.unambiguous.candidate.anchor <- df.revised.atac.unambiguous.orientation %>%
  transmute(
    loop_id,
    resolution,
    anchor_side = candidate_anchor_side,
    anchor_chr = candidate_anchor_chr,
    anchor_start = candidate_anchor_start,
    anchor_end = candidate_anchor_end,
    anchor_width_bp = candidate_anchor_width_bp,
    per_anchor_non_tss_bp = candidate_non_tss_anchor_bp,
    per_anchor_non_tss_atac_any = candidate_non_tss_atac_overlap_any,
    per_anchor_non_tss_atac_ge50 = candidate_non_tss_atac_overlap_ge50,
    per_anchor_non_tss_atac_overlap_bp = candidate_non_tss_atac_overlap_bp
  )

gr.revised.atac.unambiguous.candidate.anchor <- GRanges(
  seqnames = df.revised.atac.unambiguous.candidate.anchor$anchor_chr,
  ranges = IRanges(start = df.revised.atac.unambiguous.candidate.anchor$anchor_start, end = df.revised.atac.unambiguous.candidate.anchor$anchor_end),
  loop_id = df.revised.atac.unambiguous.candidate.anchor$loop_id,
  resolution = df.revised.atac.unambiguous.candidate.anchor$resolution,
  anchor_side = df.revised.atac.unambiguous.candidate.anchor$anchor_side
)

df.revised.atac.candidate.non.tss.fragment <- subtract_exclusion_from_anchor_ranges(gr.revised.atac.unambiguous.candidate.anchor, gr.known.tss.exclusion) %>%
  mutate(fragment_row_id = row_number(), .before = 1)

gr.revised.atac.candidate.non.tss.fragment <- GRanges(
  seqnames = df.revised.atac.candidate.non.tss.fragment$fragment_chr,
  ranges = IRanges(start = df.revised.atac.candidate.non.tss.fragment$fragment_start, end = df.revised.atac.candidate.non.tss.fragment$fragment_end),
  loop_id = df.revised.atac.candidate.non.tss.fragment$loop_id,
  resolution = df.revised.atac.candidate.non.tss.fragment$resolution,
  anchor_side = df.revised.atac.candidate.non.tss.fragment$anchor_side
)

df.revised.atac.fragment.overlap <- summarise_anchor_interval_overlap(
  gr.revised.atac.candidate.non.tss.fragment,
  gr.atac.non.tss.revised,
  minimum.overlap.bp = atac.minimum.overlap.bp
)

df.revised.atac.candidate.non.tss.fragment.evidence <- df.revised.atac.candidate.non.tss.fragment %>%
  left_join(
    df.revised.atac.fragment.overlap %>%
      transmute(
        fragment_row_id = anchor_index,
        fragment_atac_n_intervals_any = n_overlapping_feature_intervals,
        fragment_atac_n_intervals_ge50 = n_overlapping_feature_intervals_ge_minimum,
        fragment_atac_overlap_bp = feature_overlap_bp,
        fragment_atac_overlap_fraction = feature_overlap_fraction_of_anchor,
        fragment_atac_overlap_any = has_feature_overlap_any,
        fragment_atac_overlap_ge50 = has_feature_overlap_ge_minimum
      ),
    by = "fragment_row_id"
  )

df.revised.atac.fragment.anchor.observed <- df.revised.atac.candidate.non.tss.fragment.evidence %>%
  group_by(loop_id, resolution, anchor_side) %>%
  summarise(
    n_non_tss_anchor_fragments = n(),
    fragment_non_tss_anchor_bp = sum(fragment_width_bp),
    fragment_non_tss_atac_overlap_bp = sum(fragment_atac_overlap_bp),
    fragment_non_tss_atac_any = any(fragment_atac_overlap_any),
    fragment_non_tss_atac_ge50 = any(fragment_atac_overlap_ge50),
    .groups = "drop"
  )

df.revised.atac.fragment.anchor.summary <- df.revised.atac.unambiguous.candidate.anchor %>%
  left_join(df.revised.atac.fragment.anchor.observed, by = c("loop_id", "resolution", "anchor_side")) %>%
  mutate(
    n_non_tss_anchor_fragments = coalesce(n_non_tss_anchor_fragments, 0L),
    fragment_non_tss_anchor_bp = coalesce(fragment_non_tss_anchor_bp, 0L),
    fragment_non_tss_atac_overlap_bp = coalesce(fragment_non_tss_atac_overlap_bp, 0L),
    fragment_non_tss_atac_any = coalesce(fragment_non_tss_atac_any, FALSE),
    fragment_non_tss_atac_ge50 = coalesce(fragment_non_tss_atac_ge50, FALSE),
    fragment_vs_anchor_ge50_status = case_when(
      per_anchor_non_tss_atac_ge50 & fragment_non_tss_atac_ge50 ~ "positive_by_both_methods",
      per_anchor_non_tss_atac_ge50 ~ "per_anchor_only",
      fragment_non_tss_atac_ge50 ~ "fragment_only",
      TRUE ~ "negative_by_both_methods"
    )
  )

# Verify that fragment subtraction and per-anchor residual lengths agree.
assert_analysis_condition(
  !any(df.revised.atac.fragment.anchor.summary$per_anchor_non_tss_bp != df.revised.atac.fragment.anchor.summary$fragment_non_tss_anchor_bp),
  "Fragment-level TSS subtraction disagrees with per-anchor residual bases."
)

df.revised.atac.method.comparison.summary <- bind_rows(
  df.revised.atac.fragment.anchor.summary,
  df.revised.atac.fragment.anchor.summary %>% mutate(resolution = "ALL")
) %>%
  count(resolution, fragment_vs_anchor_ge50_status, name = "n_anchors") %>%
  complete(
    resolution,
    fragment_vs_anchor_ge50_status = c("positive_by_both_methods", "per_anchor_only", "fragment_only", "negative_by_both_methods"),
    fill = list(n_anchors = 0L)
  ) %>%
  group_by(resolution) %>%
  mutate(pct_anchors = round(100 * n_anchors / sum(n_anchors), 1)) %>%
  ungroup() %>%
  arrange(resolution, fragment_vs_anchor_ge50_status)

df.revised.atac.analysis.definition <- tribble(
  ~analysis_item, ~definition,
  "ATAC_interpretation", paste0(
    "ATAC overlap supports open chromatin only; it does not validate enhancer ",
    "function or a promoter-enhancer interaction."
  ),
  "known_TSS_exclusion", paste0(
    "Union of strand-aware Ensembl transcript TSS and lifted EPD TSS, expanded ",
    "symmetrically by +/-1 kb in rn7 coordinates."
  ),
  "primary_directional_set", paste0(
    "Loops with primary TSS +/-1-kb promoter-window evidence at exactly one ",
    "anchor; the opposite anchor is evaluated as a candidate regulatory anchor."
  ),
  "both_direct_anchors", paste0(
    "Loops with primary TSS +/-1-kb promoter-window evidence at both anchors ",
    "are retained as a promoter-promoter-compatible class without forced ",
    "direction."
  ),
  "minimum_overlap", paste0(
    "Any-bp overlap is reported for sensitivity; >=50 bp overlap with one ",
    "reduced ATAC interval is the primary robust-support flag."
  ),
  "fragment_sensitivity", paste0(
    "TSS regions are subtracted independently from each candidate anchor and ",
    "the residual fragments are tested without losing original anchor identity."
  ),
  "paired_test_interpretation", paste0(
    "McNemar tests compare paired accessibility states conditional on selection ",
    "of loops with one direct promoter/TSS anchor; they do not validate the ",
    "opposite anchor as a functional enhancer."
  ),
  "proximal_catalog", paste0(
    "ATAC is available in the complete anchor evidence table, but the full ",
    "1-200 kb proximal catalog is not promoted to primary directional evidence."
  )
)

message(
  "Revised non-TSS ATAC support (>=50 bp) was found at the opposite anchor in ",
  sum(
    df.revised.atac.unambiguous.orientation$
      candidate_non_tss_atac_overlap_ge50
  ),
  " of ",
  nrow(df.revised.atac.unambiguous.orientation),
  " single-direct-promoter/TSS loops."
)

################################################################################
# 6. Add transcript-contained and related positional flags
#
# Transcript position is retained as annotation rather than used as a filter.
# Every Ensembl isoform linked to a direct or proximal loop-anchor-gene
# assignment is evaluated against both the complete loop span (including the
# anchors) and the inter-anchor interval. This avoids selecting one transcript,
# requiring a last exon, or discarding a gene solely because one isoform crosses
# a loop boundary.
################################################################################

df.direct.assignment.for.position <- df.direct.promoter.tss.gene.assignment %>%
  mutate(
    assignment_tier = "primary_direct_TSS_plus_minus_1kb",
    .before = 1
  )

df.proximal.assignment.for.position <- df.proximal.promoter.tss.gene.assignment %>%
  mutate(
    assignment_tier = "exploratory_proximal_1bp_to_200kb",
    .before = 1
  )

# Expand each assignment to all matching Ensembl transcripts. EPD-only genes
# without a matching Ensembl transcript remain in the output with an explicit
# missing-annotation class.
df.direct.transcript.position.detail <- annotate_assignment_transcript_positions(
  df.assignment = df.direct.promoter.tss.gene.assignment,
  df.transcript = df.true.tss.transcript,
  df.loop = df.loop.universe,
  assignment.tier = "primary_direct_TSS_plus_minus_1kb",
  promoter.window.flank.bp = promoter.window.flank.bp
)

df.proximal.transcript.position.detail <- annotate_assignment_transcript_positions(
  df.assignment = df.proximal.promoter.tss.gene.assignment,
  df.transcript = df.true.tss.transcript,
  df.loop = df.loop.universe,
  assignment.tier = "exploratory_proximal_1bp_to_200kb",
  promoter.window.flank.bp = promoter.window.flank.bp
)

# Summarise isoform-level flags at the unique loop-anchor-gene assignment level.
# Both "any" and "all" fields are reported because transcript choice can change
# the apparent containment state of a gene.
df.direct.gene.assignment.position.flags <- summarise_assignment_transcript_positions(
  df.position.detail = df.direct.transcript.position.detail,
  df.assignment = df.direct.assignment.for.position
)

df.proximal.gene.assignment.position.flags <- summarise_assignment_transcript_positions(
  df.position.detail = df.proximal.transcript.position.detail,
  df.assignment = df.proximal.assignment.for.position
)

# Confirm that transcript-position summaries preserve all gene assignments.
assert_analysis_condition(
  nrow(df.direct.gene.assignment.position.flags) ==
    nrow(df.direct.promoter.tss.gene.assignment) &&
    nrow(df.proximal.gene.assignment.position.flags) ==
      nrow(df.proximal.promoter.tss.gene.assignment),
  "Transcript-position summaries did not preserve all gene assignments."
)

# A primary Ensembl assignment must recover at least one transcript whose TSS
# +/-1-kb promoter window overlaps the assigned anchor. Exact TSS-point overlap
# remains a descriptive strict-sensitivity flag. EPD-only assignments are not
# subjected to this Ensembl-transcript consistency check.
assert_analysis_condition(
  !any(
    df.direct.gene.assignment.position.flags$has_direct_true_tss &
      !df.direct.gene.assignment.position.flags$
        any_transcript_promoter_window_overlaps_assigned_anchor
  ),
  paste0(
    "A primary Ensembl assignment lacks a matching anchor-overlapping TSS ",
    "+/-1-kb promoter window."
  )
)

df.transcript.position.summary.input <- bind_rows(
  df.direct.gene.assignment.position.flags,
  df.proximal.gene.assignment.position.flags
)

df.transcript.position.summary <- bind_rows(
  df.transcript.position.summary.input,
  df.transcript.position.summary.input %>% mutate(resolution = "ALL")
) %>%
  group_by(assignment_tier, resolution) %>%
  summarise(
    n_loop_anchor_gene_assignments = n(),
    n_with_ensembl_transcript_annotation = sum(has_ensembl_transcript_annotation),
    pct_with_ensembl_transcript_annotation = round(100 * mean(has_ensembl_transcript_annotation), 1),
    n_any_transcript_fully_within_loop_span = sum(any_transcript_fully_within_loop_span),
    pct_any_transcript_fully_within_loop_span = round(100 * mean(any_transcript_fully_within_loop_span), 1),
    n_all_transcripts_fully_within_loop_span = sum(all_transcripts_fully_within_loop_span),
    pct_all_transcripts_fully_within_loop_span = round(100 * mean(all_transcripts_fully_within_loop_span), 1),
    n_any_transcript_fully_within_inter_anchor_interval = sum(any_transcript_fully_within_inter_anchor_interval),
    pct_any_transcript_fully_within_inter_anchor_interval = round(100 * mean(any_transcript_fully_within_inter_anchor_interval), 1),
    n_overlap_without_full_transcript_containment = sum(transcript_containment_summary == "overlap_without_full_transcript_containment"),
    n_all_annotated_transcripts_outside_loop_span = sum(transcript_containment_summary == "all_annotated_transcripts_outside_loop_span"),
    .groups = "drop"
  ) %>%
  arrange(assignment_tier, resolution)

df.transcript.position.definition <- tribble(
  ~field_or_rule, ~definition,
  "loop_span", paste0(
    "Inclusive genomic span from the outer start to the outer end of the two ",
    "anchors."
  ),
  "inter_anchor_interval", paste0(
    "Open interval between the left anchor end and right anchor start; anchor ",
    "bases are excluded."
  ),
  "transcript_fully_within_loop_span", paste0(
    "The complete Ensembl transcript body lies inside the inclusive loop span."
  ),
  "all_vs_any_transcript_flags", paste0(
    "All annotated isoforms are retained; no canonical or nearest transcript ",
    "is selected for containment."
  ),
  "analysis_role", paste0(
    "Transcript containment is descriptive evidence only and is not a loop or ",
    "gene-assignment retention filter."
  )
)

message(
  "Added transcript-position flags to ",
  nrow(df.direct.gene.assignment.position.flags),
  " direct and ",
  nrow(df.proximal.gene.assignment.position.flags),
  " proximal loop-anchor-gene assignments."
)

################################################################################
# 7. Predicted CTCF motif annotation
#
# CTCF sequence predictions are retained as an independent structural annotation.
# Motif counts do not select loops, define regulatory categories, determine
# direction, or imply in-vivo CTCF occupancy.
################################################################################

df.ctcf.anchor.count <- count_loop_anchor_feature_overlaps(
  list.gr.loop.anchor.by.side,
  gr.ctcf.motif,
  count.prefix = "ctcf_count",
  ignore.strand = TRUE
)

df.ctcf.evidence <- df.loop.universe %>%
  dplyr::select(loop_id, resolution) %>%
  left_join(df.ctcf.anchor.count, by = "loop_id") %>%
  transmute(
    loop_id,
    resolution,
    predicted_ctcf_motif_interval_count_anchor1 = coalesce(ctcf_count_anchor1, 0L),
    predicted_ctcf_motif_interval_count_anchor2 = coalesce(ctcf_count_anchor2, 0L)
  ) %>%
  mutate(
    predicted_ctcf_motif_intervals_any_anchor1 = predicted_ctcf_motif_interval_count_anchor1 > 0L,
    predicted_ctcf_motif_intervals_any_anchor2 = predicted_ctcf_motif_interval_count_anchor2 > 0L,
    predicted_ctcf_motif_intervals_both_anchors = (predicted_ctcf_motif_intervals_any_anchor1 & predicted_ctcf_motif_intervals_any_anchor2),
    predicted_ctcf_motif_annotation_class = case_when(
      predicted_ctcf_motif_intervals_both_anchors ~ "predicted_motif_intervals_at_both_anchors",
      predicted_ctcf_motif_intervals_any_anchor1 | predicted_ctcf_motif_intervals_any_anchor2 ~ "predicted_motif_intervals_at_one_anchor",
      TRUE ~ "no_predicted_motif_interval_overlap"
    )
  )

################################################################################
# 8. Layered loop evidence and mutually exclusive categories
#
# Regulatory categories use direct promoter/TSS and TSS-excluded ATAC evidence.
# Predicted CTCF motif fields are joined only as independent annotations.
################################################################################

df.revised.loop.evidence <- df.loop.universe %>%
  left_join(
    df.direct.promoter.tss.loop.summary %>%
      dplyr::select(loop_id, n_direct_anchor_sides, n_direct_genes_across_anchors, has_any_direct_promoter_tss, has_direct_promoter_tss_both_anchors, direct_anchor_assignment_class),
    by = "loop_id"
  ) %>%
  left_join(
    df.proximal.promoter.tss.loop.summary %>%
      dplyr::select(loop_id, n_proximal_anchor_sides, n_proximal_genes_across_anchors, has_any_proximal_promoter_tss, n_secondary_inward_proximal_anchor_sides_10kb, n_secondary_inward_proximal_genes_10kb_across_anchors, has_any_secondary_inward_proximal_10kb, has_any_direct_or_proximal_promoter_tss),
    by = "loop_id"
  ) %>%
  left_join(
    df.revised.atac.loop.summary %>%
      dplyr::select(loop_id, promoter_anchor_raw_atac_ge50, candidate_anchor_raw_atac_ge50, candidate_anchor_non_tss_atac_any, candidate_anchor_non_tss_atac_ge50, revised_atac_evidence_class),
    by = "loop_id"
  ) %>%
  left_join(
    df.dual.promoter.atac.directional.classification %>%
      dplyr::select(loop_id, promoter_window_atac_ge50_anchor1, promoter_window_atac_ge50_anchor2, residual_non_tss_atac_ge50_anchor1, residual_non_tss_atac_ge50_anchor2, supports_anchor1_promoter_to_anchor2_regulatory, supports_anchor2_promoter_to_anchor1_regulatory, n_accessible_promoter_anchors, n_supported_regulatory_directions, dual_promoter_directional_category, supported_regulatory_direction),
    by = "loop_id"
  ) %>%
  left_join(df.ctcf.evidence %>% dplyr::select(-resolution), by = "loop_id") %>%
  mutate(
    n_direct_anchor_sides = coalesce(n_direct_anchor_sides, 0L),
    dual_promoter_directional_category = replace_na(dual_promoter_directional_category, "not_dual_promoter_loop"),
    supported_regulatory_direction = replace_na(supported_regulatory_direction, "not_dual_promoter_loop"),
    revised_putative_regulatory_support = (n_direct_anchor_sides == 1L & coalesce(candidate_anchor_non_tss_atac_ge50, FALSE)),
    revised_promoter_promoter_compatible = n_direct_anchor_sides == 2L,
    revised_single_promoter_without_opposite_atac = (n_direct_anchor_sides == 1L & !coalesce(candidate_anchor_non_tss_atac_ge50, FALSE)),
    revised_no_direct_promoter_tss = n_direct_anchor_sides == 0L,
    revised_major_category = case_when(
      revised_putative_regulatory_support ~ "putative_regulatory_direct_promoter_TSS_opposite_nonTSS_ATAC",
      revised_promoter_promoter_compatible ~ "promoter_promoter_compatible_both_direct_anchors",
      revised_single_promoter_without_opposite_atac ~ "single_direct_promoter_TSS_without_opposite_nonTSS_ATAC",
      TRUE ~ "no_direct_promoter_TSS"
    ),
    revised_detailed_category = case_when(
      revised_putative_regulatory_support ~ "single_direct_promoter_TSS_with_opposite_nonTSS_ATAC",
      revised_promoter_promoter_compatible ~ str_c("promoter_promoter_compatible__", dual_promoter_directional_category),
      revised_single_promoter_without_opposite_atac ~ "single_direct_promoter_TSS_without_opposite_nonTSS_ATAC",
      has_any_secondary_inward_proximal_10kb ~ "no_direct_promoter_TSS_secondary_inward_10kb",
      has_any_proximal_promoter_tss ~ "no_direct_promoter_TSS_exploratory_proximal_200kb_only",
      TRUE ~ "no_direct_or_proximal_promoter_TSS_evidence"
    )
  ) %>%
  arrange(chr1, start1, end1, chr2, start2, end2)

assert_analysis_condition(
  nrow(df.revised.loop.evidence) == nrow(df.loop.universe) &&
    n_distinct(df.revised.loop.evidence$loop_id) == nrow(df.loop.universe),
  "Revised loop evidence is not exactly one row per pooled loop."
)

df.revised.loop.category.definition <- tribble(
  ~revised_major_category, ~required_evidence, ~interpretation,
  "putative_regulatory_direct_promoter_TSS_opposite_nonTSS_ATAC",
  paste0(
    "Exactly one anchor directly overlaps a promoter/TSS +/-1-kb window; ",
    "the opposite anchor has >=50 bp true-TSS-excluded ATAC overlap."
  ),
  paste0(
    "Open-chromatin-supported putative promoter-to-distal-element contact; ",
    "not a validated enhancer or functional P-E interaction."
  ),
  "promoter_promoter_compatible_both_direct_anchors",
  "Both anchors directly overlap promoter/TSS +/-1-kb windows.",
  paste0(
    "Promoter-promoter-compatible contact retained without forcing either ",
    "anchor to be an enhancer."
  ),
  "single_direct_promoter_TSS_without_opposite_nonTSS_ATAC",
  paste0(
    "Exactly one anchor directly overlaps a promoter/TSS +/-1-kb window; ",
    "the opposite anchor lacks >=50 bp true-TSS-excluded ATAC overlap."
  ),
  "Directionally assignable promoter contact without distal open-chromatin support.",
  "no_direct_promoter_TSS",
  "Neither anchor directly overlaps a promoter/TSS +/-1-kb window.",
  "Pooled HiCCUPS call retained without direct promoter/TSS assignment."
)

df.revised.loop.category.summary <- df.revised.loop.evidence %>%
  count(revised_major_category, name = "n_loops") %>%
  mutate(pct_pooled_loops = round(100 * n_loops / sum(n_loops), 1)) %>%
  arrange(desc(n_loops), revised_major_category)

df.revised.loop.category.by.resolution <- df.revised.loop.evidence %>%
  count(resolution, revised_major_category, name = "n_loops") %>%
  group_by(resolution) %>%
  mutate(pct_within_resolution = round(100 * n_loops / sum(n_loops), 1)) %>%
  ungroup() %>%
  arrange(resolution, desc(n_loops), revised_major_category)

message("Revised loop categories cover all ", nrow(df.loop.universe), " calls.")

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
# 13. Exploratory strain-level loop sharing
#
# These results describe exact loop-coordinate sharing among the 10 libraries.
# They are not used to claim strain-specific chromatin biology because each
# strain has one Hi-C library and loop recovery is depth-sensitive.
################################################################################

df.sample.loop.presence <- df.sample.loop.1based %>%
  filter(passes_lt2mb) %>%
  distinct(loop_id, resolution, sample, strain)

# Keep all 10 sequenced strains in descriptive loop-sharing summaries.
excluded.strains <- character(0)

df.sample.loop.presence.pairwise <- df.sample.loop.presence %>% filter(!strain %in% excluded.strains)

df.sample.strain.key <- df.sample.loop.presence %>%
  distinct(sample, strain) %>%
  mutate(
    included_in_pairwise_loop_overlap = !strain %in% excluded.strains,
    exclusion_reason = if_else(included_in_pairwise_loop_overlap, NA_character_, "Excluded because a matching genotype sample was not available for an optional genetic-distance analysis.")
  ) %>%
  arrange(strain, sample)

strain.levels <- df.sample.loop.presence.pairwise %>%
  distinct(strain) %>%
  arrange(strain) %>%
  pull(strain)

df.strain.loop.count.by.resolution <- bind_rows(
  df.sample.loop.presence.pairwise %>% count(strain, sample, resolution, name = "n_loops"),
  df.sample.loop.presence.pairwise %>% count(strain, sample, name = "n_loops") %>% mutate(resolution = "ALL")
) %>%
  arrange(strain, sample, resolution)

df.strain.loop.presence.matrix <- df.sample.loop.presence.pairwise %>%
  distinct(loop_id, resolution, strain) %>%
  mutate(present = 1L) %>%
  pivot_wider(names_from = strain, values_from = present, values_fill = 0) %>%
  left_join(df.loop.universe %>% dplyr::select(loop_id, chr1, start1, end1, chr2, start2, end2, loop_distance), by = "loop_id") %>%
  relocate(chr1, start1, end1, chr2, start2, end2, loop_distance, .after = resolution) %>%
  arrange(resolution, chr1, start1, end1, chr2, start2, end2)

df.shared.loop.member <- df.sample.loop.presence.pairwise %>%
  distinct(loop_id, resolution, strain) %>%
  group_by(loop_id, resolution) %>%
  summarise(strains = list(sort(unique(strain))), n_strains_detected = n_distinct(strain), .groups = "drop") %>%
  filter(n_strains_detected > 1)

# Enumerate strain pairs only for loops observed in more than one strain.
if (nrow(df.shared.loop.member) > 0) {
  df.shared.loop.pair.by.strain <- map_dfr(
    seq_len(nrow(df.shared.loop.member)),
    function(i) {
      strain.pairs <- t(combn(df.shared.loop.member$strains[[i]], 2))
      tibble(loop_id = df.shared.loop.member$loop_id[[i]], resolution = df.shared.loop.member$resolution[[i]], strain1 = strain.pairs[, 1], strain2 = strain.pairs[, 2])
    }
  ) %>%
    count(resolution, strain1, strain2, name = "n_shared_loops") %>%
    arrange(resolution, strain1, strain2)
} else {
  df.shared.loop.pair.by.strain <- tibble(resolution = character(), strain1 = character(), strain2 = character(), n_shared_loops = integer())
}

df.pairwise.loop.overlap.by.strain <- bind_rows(
  create_pairwise_loop_overlap(df.sample.loop.presence.pairwise, strain.levels, "ALL"),
  map_dfr(sort(unique(df.sample.loop.presence.pairwise$resolution)), function(resolution.i) {
    create_pairwise_loop_overlap(df.sample.loop.presence.pairwise %>% filter(resolution == resolution.i), strain.levels, resolution.i)
  })
) %>%
  arrange(resolution, strain1, strain2)

df.loop.sharing.distribution.by.strain <- bind_rows(
  df.sample.loop.presence.pairwise %>%
    distinct(loop_id, resolution, strain) %>%
    group_by(loop_id, resolution) %>%
    summarise(n_strains_detected = n_distinct(strain), .groups = "drop") %>%
    count(resolution, n_strains_detected, name = "n_loops"),
  df.sample.loop.presence.pairwise %>%
    distinct(loop_id, strain) %>%
    group_by(loop_id) %>%
    summarise(n_strains_detected = n_distinct(strain), .groups = "drop") %>%
    count(n_strains_detected, name = "n_loops") %>%
    mutate(resolution = "ALL")
) %>%
  arrange(resolution, n_strains_detected)

################################################################################
# 14. Sequencing-depth sensitivity/QC
#
# Contact-normalized ratios and regression residuals are reviewer-facing QC.
# They do not replace read downsampling followed by mapping and loop re-calling,
# and they do not make one-library-per-strain comparisons biological replicates.
################################################################################

df.library.complexity <- readr::read_tsv(library.complexity.file, show_col_types = FALSE) %>% mutate(Strain = as.character(Strain))

required.depth.columns <- c("Strain", "Sequenced_RP", "Unique_Reads", "Alignable_Normal_N_Chimeric", "Hi-C_Contacts", "Long_Range_20Kb", "PCR_Duplicates", "Optical_Duplicates")
check_required_columns(df.library.complexity, required.depth.columns, "df.library.complexity")

df.depth.qc.by.strain <- df.strain.loop.count.by.resolution %>%
  left_join(df.library.complexity, by = c("strain" = "Strain")) %>%
  mutate(
    sequenced_read_pairs_millions = Sequenced_RP / 1e6,
    unique_reads_millions = Unique_Reads / 1e6,
    hic_contacts_millions = `Hi-C_Contacts` / 1e6,
    loops_per_100m_sequenced_read_pairs = safe_ratio(n_loops, Sequenced_RP / 1e8),
    loops_per_100m_unique_reads = safe_ratio(n_loops, Unique_Reads / 1e8),
    loops_per_100m_hic_contacts = safe_ratio(n_loops, `Hi-C_Contacts` / 1e8),
    unique_read_pct = 100 * safe_ratio(Unique_Reads, Sequenced_RP),
    duplicate_pct = 100 * safe_ratio(PCR_Duplicates + Optical_Duplicates, Sequenced_RP)
  ) %>%
  group_by(resolution) %>%
  group_modify(~ add_depth_residuals(.x)) %>%
  ungroup() %>%
  arrange(resolution, desc(n_loops))

# Warn when sequencing-depth metadata fail to match one or more strains.
if (any(is.na(df.depth.qc.by.strain$`Hi-C_Contacts`))) {
  warning("One or more strains did not match the library-complexity table.", call. = FALSE)
}

depth.metrics <- c("Sequenced_RP", "Unique_Reads", "Alignable_Normal_N_Chimeric", "Hi-C_Contacts", "Long_Range_20Kb")

df.depth.loop.correlation.by.resolution <- map_dfr(
  sort(unique(df.depth.qc.by.strain$resolution)),
  function(resolution.i) {
    df.depth.i <- df.depth.qc.by.strain %>% filter(resolution == resolution.i)
    map_dfr(depth.metrics, function(metric.i) {
      pearson.result <- safe_correlation_summary(df.depth.i$n_loops, df.depth.i[[metric.i]], method = "pearson")
      spearman.result <- safe_correlation_summary(df.depth.i$n_loops, df.depth.i[[metric.i]], method = "spearman")
      tibble(
        resolution = resolution.i, depth_metric = metric.i, n_strains = pearson.result$n,
        pearson_r = pearson.result$estimate, pearson_p = pearson.result$p_value,
        spearman_rho = spearman.result$estimate, spearman_p = spearman.result$p_value
      )
    })
  }
) %>%
  arrange(resolution, depth_metric)

df.depth.qc.summary <- df.depth.qc.by.strain %>%
  group_by(resolution) %>%
  summarise(
    n_strains = n_distinct(strain),
    mean_loops = round(mean(n_loops), 1),
    sd_loops = round(stats::sd(n_loops), 1),
    min_loops = min(n_loops),
    max_loops = max(n_loops),
    mean_sequenced_read_pairs_millions = round(mean(sequenced_read_pairs_millions), 1),
    mean_unique_reads_millions = round(mean(unique_reads_millions), 1),
    mean_hic_contacts_millions = round(mean(hic_contacts_millions), 1),
    mean_loops_per_100m_hic_contacts = round(mean(loops_per_100m_hic_contacts), 1),
    .groups = "drop"
  ) %>%
  arrange(resolution)

df.depth.analysis.interpretation <- tibble(
  analysis = c("loop_count_vs_depth", "loops_per_100m_hic_contacts", "depth_residual_loop_count"),
  role = "exploratory_QC",
  valid_interpretation = c(
    "Quantifies sensitivity of loop recovery to sequencing/contact depth.",
    "Descriptive contact-normalized loop-recovery rate.",
    "Descriptive residual after regressing loop count on log10 Hi-C contacts."
  ),
  invalid_interpretation = c(
    "Does not establish strain-specific biological differences.",
    "Does not replace read downsampling and loop re-calling.",
    "Does not create biological replication for one library per strain."
  )
)

################################################################################
# 15. Optional loop-sharing vs. true genetic-distance analysis
#
# A panel-type label is not a genetic-distance matrix. This analysis runs only
# when a real long-format SNP/VCF-derived distance table is available.
################################################################################

df.genetic.distance.status <- tibble(
  requested_file = genetic.distance.file,
  file_exists = file.exists(genetic.distance.file),
  analysis_status = if_else(file.exists(genetic.distance.file), "ready_to_run", "skipped_no_true_genetic_distance_matrix"),
  interpretation = if_else(
    file.exists(genetic.distance.file),
    "Pairwise loop Jaccard distance will be compared with the supplied SNP/VCF-derived genetic distance.",
    "No proxy was substituted. A real SNP/VCF-derived pairwise distance matrix is required for this analysis."
  )
)

df.pairwise.loop.genetic.distance <- tibble(
  resolution = character(), strain1 = character(), strain2 = character(), pair_strain_a = character(), pair_strain_b = character(),
  n_shared_loops = integer(), n_union_loops = integer(), jaccard_similarity = double(), loop_jaccard_distance = double(), genetic_distance = double()
)

df.loop.genetic.distance.correlation <- tibble(
  resolution = character(), n_pairs = integer(), spearman_rho = double(), spearman_p = double(), analysis_note = character()
)

# Run the genetic-distance comparison only when a true matrix is supplied.
if (file.exists(genetic.distance.file)) {
  df.genetic.distance <- read_genetic_distance_long(genetic.distance.file)

  df.pairwise.loop.genetic.distance <- df.pairwise.loop.overlap.by.strain %>%
    filter(strain1 != strain2, !is.na(jaccard_similarity)) %>%
    mutate(pair_strain_a = pmin(strain1, strain2), pair_strain_b = pmax(strain1, strain2)) %>%
    distinct(resolution, pair_strain_a, pair_strain_b, .keep_all = TRUE) %>%
    left_join(df.genetic.distance %>% dplyr::select(pair_strain_a, pair_strain_b, genetic_distance), by = c("pair_strain_a", "pair_strain_b")) %>%
    mutate(loop_jaccard_distance = jaccard_distance) %>%
    dplyr::select(resolution, strain1, strain2, pair_strain_a, pair_strain_b, n_shared_loops, n_union_loops, jaccard_similarity, loop_jaccard_distance, genetic_distance)

  df.loop.genetic.distance.correlation <- df.pairwise.loop.genetic.distance %>%
    filter(!is.na(loop_jaccard_distance), !is.na(genetic_distance)) %>%
    group_by(resolution) %>%
    group_modify(function(.x, .y) {
      correlation.result <- safe_correlation_summary(.x$loop_jaccard_distance, .x$genetic_distance, method = "spearman")
      tibble(
        n_pairs = correlation.result$n,
        spearman_rho = correlation.result$estimate,
        spearman_p = correlation.result$p_value,
        analysis_note = "Exploratory exact-loop Jaccard distance vs. supplied SNP/VCF-derived genetic distance."
      )
    }) %>%
    ungroup()

  df.genetic.distance.status <- df.genetic.distance.status %>%
    mutate(
      analysis_status = "completed_with_supplied_true_distance_matrix",
      n_distance_pairs_supplied = nrow(df.genetic.distance),
      n_loop_pairs_with_distance = sum(!is.na(df.pairwise.loop.genetic.distance$genetic_distance))
    )
} else {
  df.genetic.distance.status <- df.genetic.distance.status %>%
    mutate(n_distance_pairs_supplied = 0L, n_loop_pairs_with_distance = 0L)
}

# Overwrite legacy proxy output names with an explicit deprecation record so
# stale numeric proxy results cannot be mistaken for a valid genetic analysis.
df.genetic.proxy.deprecation <- tibble(
  analysis_status = "deprecated_not_performed",
  reason = paste0(
    "HRDP panel membership is not a quantitative genetic-distance matrix. ",
    "Provide SNP/VCF-derived pairwise distances instead."
  )
)

# Record the final production state once, after every revised component has run.
df.revised.assignment.pipeline.status <- tribble(
  ~analysis_component, ~current_status, ~replacement_section,
  ~uses_revised_direct_assignment,
  "direct_promoter_TSS_anchor_overlap",
  "revised_complete", "Section 3", TRUE,
  "proximal_promoter_TSS_assignment",
  "revised_complete_secondary_inward_10kb_and_exploratory_200kb_catalog",
  "Section 4", TRUE,
  "ATAC_support",
  "revised_complete_true_TSS_excluded_direct_orientation",
  "Section 5", TRUE,
  "transcript_containment_flags",
  "revised_complete_descriptive_not_filtering",
  "Section 6", TRUE,
  "loop_categories",
  "revised_complete_direct_promoter_TSS_and_ATAC_categories",
  "Section 8", TRUE,
  "downstream_gene_resource_GO",
  "revised_complete_multi_gene_direct_assignment",
  "Section 12", TRUE
)

################################################################################
# 15-1. Revised manuscript figure suite
#
# These figures replace legacy plots with summaries derived from the revised
# coordinate-normalized, evidence-layered analysis. Library support, predicted
# CTCF motifs, ATAC overlap, and GO enrichment are described without treating
# them as strain specificity, experimental CTCF occupancy, enhancer validation,
# or validation of individual regulatory contacts.
################################################################################

resolution.colors <- c(
  "5K" = "#4D9ACB",
  "10K" = "#2878B5",
  "25K" = "#173B7A"
)
category.colors <- c(
  "putative_regulatory_direct_promoter_TSS_opposite_nonTSS_ATAC" = "#00BA38",
  "promoter_promoter_compatible_both_direct_anchors" = "#C77CFF",
  "single_direct_promoter_TSS_without_opposite_nonTSS_ATAC" = "#619CFF",
  "no_direct_promoter_TSS" = "#F8766D"
)
category.labels <- c(
  "putative_regulatory_direct_promoter_TSS_opposite_nonTSS_ATAC" =
    "Promoter/TSS + opposite\nTSS-excluded ATAC",
  "promoter_promoter_compatible_both_direct_anchors" =
    "Promoter/TSS annotations\nat both anchors",
  "single_direct_promoter_TSS_without_opposite_nonTSS_ATAC" =
    "Single promoter/TSS;\nno opposite ATAC support",
  "no_direct_promoter_TSS" = "No direct\npromoter/TSS"
)

theme.revised.figure <- theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.subtitle = element_text(size = 9, color = "grey30", hjust = 0.5),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

# Figure 1: show the valid contact depth and resolution-specific loop-call
# yield for every library, plus their descriptive full-depth relationship.
df.figure1.library.depth <- df.depth.qc.by.strain %>%
  filter(resolution == "ALL") %>%
  distinct(strain, sample, hic_contacts_millions, n_loops) %>%
  arrange(hic_contacts_millions) %>%
  mutate(strain = factor(strain, levels = strain))
df.figure1.loop.yield <- df.strain.loop.count.by.resolution %>%
  filter(resolution != "ALL") %>%
  mutate(
    strain = factor(strain, levels = levels(df.figure1.library.depth$strain)),
    resolution = factor(resolution, levels = names(resolution.colors))
  )

plot.figure1a.library.depth <- ggplot(
  df.figure1.library.depth,
  aes(strain, hic_contacts_millions)
) +
  geom_col(fill = "#3B6B8C", width = 0.72) +
  coord_flip() +
  labs(
    tag = "a",
    title = "Valid Hi-C contacts",
    x = NULL,
    y = "Valid MAPQ >=30 contacts (millions)"
  ) +
  theme.revised.figure

plot.figure1b.loop.yield <- ggplot(
  df.figure1.loop.yield,
  aes(strain, n_loops, fill = resolution)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  scale_fill_manual(values = resolution.colors, drop = FALSE) +
  coord_flip() +
  labs(
    tag = "b",
    title = "HiCCUPS calls by resolution",
    x = NULL,
    y = "Number of exact calls",
    fill = "Resolution"
  ) +
  theme.revised.figure

plot.figure1c.depth.loop <- ggplot(
  df.figure1.library.depth,
  aes(hic_contacts_millions, n_loops)
) +
  geom_smooth(
    method = "lm", formula = y ~ x, se = FALSE,
    color = "grey55", linewidth = 0.6
  ) +
  geom_point(color = "#B23A48", size = 2.5) +
  geom_text(aes(label = sample),
    nudge_y = 170, size = 2.7,
    check_overlap = TRUE
  ) +
  labs(
    tag = "c",
    title = "Full-depth depth-call relationship",
    subtitle = "Descriptive sequencing-depth QC",
    x = "Valid MAPQ >=30 contacts (millions)",
    y = "Exact calls"
  ) +
  theme.revised.figure

# Figure 1 presentation outputs: pair contact depth with its call-yield
# relationship, and retain the resolution-specific call panel separately.
plot.figure1.depth.relationship <- patchwork::wrap_plots(
  plot.figure1a.library.depth + labs(tag = "a"),
  plot.figure1c.depth.loop + labs(tag = "b"),
  nrow = 1,
  widths = c(1, 1)
)
plot.figure1.calls.by.resolution <- plot.figure1b.loop.yield + labs(tag = NULL)

figure1.output.files <- basename(c(
  saving_plot_dual(
    plot.figure1.depth.relationship,
    "figure1_revised_contact_depth_and_call_relationship",
    output.dir,
    width_in = 9,
    height_in = 4.8
  ),
  saving_plot_dual(
    plot.figure1.calls.by.resolution,
    "figure1_revised_hiccups_calls_by_resolution",
    output.dir,
    width_in = 6.5,
    height_in = 4.8
  )
))

# Figure 2: summarize exact, resolution-specific pooled calls by chromosome.
df.figure2.chromosome.resolution <- df.loop.universe %>%
  count(chr1, resolution, name = "n_loops") %>%
  mutate(
    chromosome_label = str_remove(chr1, "^chr"),
    chromosome_order = case_when(
      chromosome_label == "X" ~ 100,
      chromosome_label == "Y" ~ 101,
      chromosome_label %in% c("M", "MT") ~ 102,
      TRUE ~ readr::parse_number(chromosome_label)
    )
  ) %>%
  arrange(chromosome_order) %>%
  mutate(
    chromosome_label = factor(
      chromosome_label,
      levels = unique(chromosome_label)
    ),
    resolution = factor(resolution, levels = names(resolution.colors))
  )

plot.figure2.revised <- ggplot(
  df.figure2.chromosome.resolution,
  aes(chromosome_label, n_loops, fill = resolution)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  scale_fill_manual(values = resolution.colors, drop = FALSE) +
  labs(
    title = "Exact pooled HiCCUPS calls by chromosome and resolution",
    subtitle = "Resolution-specific call records shorter than 2 Mb",
    x = "Chromosome",
    y = "Number of exact calls",
    fill = "Resolution"
  ) +
  theme.revised.figure +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
figure2.output.files <- basename(saving_plot_dual(
  plot.figure2.revised,
  "figure2_revised_pooled_calls_by_chromosome_resolution",
  output.dir,
  width_in = 10,
  height_in = 5.2
))

# Figure 3: report library support and pairwise exact-call similarity without
# interpreting single-library records as strain-specific biological loops.
df.figure3.library.support <- df.loop.sharing.distribution.by.strain %>%
  filter(resolution == "ALL")
df.figure3.pairwise.jaccard <- df.pairwise.loop.overlap.by.strain %>%
  filter(resolution == "ALL") %>%
  mutate(
    strain1 = factor(strain1, levels = strain.levels),
    strain2 = factor(strain2, levels = rev(strain.levels))
  )

plot.figure3a.library.support <- ggplot(
  df.figure3.library.support,
  aes(factor(n_strains_detected), n_loops)
) +
  geom_col(fill = "#4472A6", width = 0.72) +
  scale_y_log10(labels = scales::comma) +
  labs(
    tag = "a",
    title = "Exact-call library support",
    subtitle = "One Hi-C library per sampled strain",
    x = "Number of supporting libraries",
    y = "Number of exact calls (log10 scale)"
  ) +
  theme.revised.figure

plot.figure3b.pairwise.jaccard <- ggplot(
  df.figure3.pairwise.jaccard,
  aes(strain1, strain2, fill = jaccard_similarity)
) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient(
    low = "#F2F2F2",
    high = "#1F5A91",
    limits = c(0, 1),
    na.value = "white"
  ) +
  coord_fixed() +
  labs(
    tag = "b",
    title = "Pairwise exact-call similarity",
    x = NULL,
    y = NULL,
    fill = "Jaccard"
  ) +
  theme.revised.figure +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7)
  )

plot.figure3.revised <- patchwork::wrap_plots(
  plot.figure3a.library.support,
  plot.figure3b.pairwise.jaccard,
  nrow = 1,
  widths = c(0.85, 1.2)
)
figure3.output.files <- basename(saving_plot_dual(
  plot.figure3.revised,
  "figure3_revised_library_support_and_pairwise_overlap",
  output.dir,
  width_in = 11,
  height_in = 5.2
))

# Figure 4: display mutually exclusive promoter/TSS and ATAC evidence classes.
df.figure4.category.summary <- df.revised.loop.category.summary %>%
  mutate(
    category_label = recode(revised_major_category, !!!category.labels),
    category_label = forcats::fct_reorder(category_label, n_loops),
    count_label = str_c(scales::comma(n_loops), " (", pct_pooled_loops, "%)")
  )
df.figure4.category.by.resolution <- df.revised.loop.category.by.resolution %>%
  mutate(
    category_label = recode(revised_major_category, !!!category.labels),
    resolution = factor(resolution, levels = names(resolution.colors))
  )

plot.figure4a.category.summary <- ggplot(
  df.figure4.category.summary,
  aes(category_label, n_loops, fill = revised_major_category)
) +
  geom_col(width = 0.72, show.legend = FALSE) +
  geom_text(aes(label = count_label), hjust = -0.08, size = 3) +
  scale_fill_manual(values = category.colors) +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.2))
  ) +
  coord_flip() +
  labs(
    tag = "a",
    title = "Evidence-layered loop-call categories",
    x = NULL,
    y = "Number of exact calls"
  ) +
  theme.revised.figure

plot.figure4b.category.resolution <- ggplot(
  df.figure4.category.by.resolution,
  aes(resolution, pct_within_resolution / 100,
    fill = revised_major_category
  )
) +
  geom_col(width = 0.68) +
  scale_fill_manual(
    values = category.colors,
    labels = category.labels
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    tag = "b",
    title = "Category composition by resolution",
    x = "Resolution",
    y = "Fraction of exact calls",
    fill = "Evidence category"
  ) +
  theme.revised.figure +
  theme(legend.position = "bottom")

plot.figure4.revised <- patchwork::wrap_plots(
  plot.figure4a.category.summary,
  plot.figure4b.category.resolution,
  nrow = 1,
  widths = c(1.1, 1),
  guides = "collect"
) & theme(legend.position = "bottom")
figure4.output.files <- basename(saving_plot_dual(
  plot.figure4.revised,
  "figure4_revised_loop_evidence_categories",
  output.dir,
  width_in = 12,
  height_in = 5.8
))

# Figure 6: retain GO Biological Process enrichment only as exploratory
# functional context for the revised putative-regulatory gene set.
figure6.output.files <- character()
if (nrow(df.revised.go.result) > 0L) {
  df.figure6.go <- df.revised.go.result %>%
    filter(gene_set == "revised_putative_all", !is.na(p.adjust)) %>%
    arrange(p.adjust, desc(FoldEnrichment)) %>%
    slice_head(n = 15L) %>%
    mutate(
      Description = str_wrap(Description, width = 42),
      Description = factor(Description, levels = rev(Description)),
      minus_log10_adjusted_p = -log10(pmax(p.adjust, .Machine$double.xmin))
    )

  plot.figure6.revised <- ggplot(
    df.figure6.go,
    aes(minus_log10_adjusted_p, Description)
  ) +
    geom_point(aes(size = Count, color = FoldEnrichment), alpha = 0.85) +
    scale_color_gradient(low = "#2A9D8F", high = "#B23A48") +
    labs(
      title = "Exploratory GO Biological Process enrichment",
      subtitle = "Interpretive functional context; not validation of individual loop calls",
      x = expression(-log[10](adjusted ~ italic(P))),
      y = NULL,
      size = "Genes",
      color = "Fold enrichment"
    ) +
    theme.revised.figure +
    theme(legend.position = "right")
  figure6.output.files <- basename(saving_plot_dual(
    plot.figure6.revised,
    "figure6_revised_exploratory_GO_BP",
    output.dir,
    width_in = 9.5,
    height_in = 6.5
  ))
}

# Supplementary Figure S1: expose resolution-specific exact-call Jaccard
# similarity as technical sensitivity information.
df.figureS1.pairwise.jaccard <- df.pairwise.loop.overlap.by.strain %>%
  filter(resolution != "ALL") %>%
  mutate(
    resolution = factor(resolution, levels = names(resolution.colors)),
    strain1 = factor(strain1, levels = strain.levels),
    strain2 = factor(strain2, levels = rev(strain.levels))
  )
plot.figureS1.revised <- ggplot(
  df.figureS1.pairwise.jaccard,
  aes(strain1, strain2, fill = jaccard_similarity)
) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_wrap(~resolution, nrow = 1) +
  scale_fill_gradient(
    low = "#F2F2F2",
    high = "#1F5A91",
    limits = c(0, 1),
    na.value = "white"
  ) +
  coord_fixed() +
  labs(
    title = "Pairwise exact-call similarity by resolution",
    subtitle = "Technical comparison among one-library-per-strain datasets",
    x = NULL,
    y = NULL,
    fill = "Jaccard"
  ) +
  theme.revised.figure +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6)
  )
figureS1.output.files <- basename(saving_plot_dual(
  plot.figureS1.revised,
  "figureS1_revised_pairwise_exact_call_jaccard_by_resolution",
  output.dir,
  width_in = 12,
  height_in = 4.8
))

revised.figure.output.files <- c(
  figure1.output.files,
  figure2.output.files,
  figure3.output.files,
  figure4.output.files,
  figure5.output.files,
  figure6.output.files,
  figureS1.output.files
)

# Compute reproducible SHA-256 checksums for source inputs and analysis scripts.
sha256_file <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop(
      "The cross-platform 'digest' package is required for SHA-256 checksums.",
      call. = FALSE
    )
  }
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

analysis.script.files <- tibble(
  input_name = c(
    "analysis_main_script",
    "coordinate_preparation_script",
    "shared_functions_script",
    "atac_matched_null_script"
  ),
  input_path = c(
    file.path(analysis.dir, "promoter_enhancer_interaction_resubmit.R"),
    coord.prep.script,
    file.path(r.files.dir, "funcs.R"),
    file.path(revision.dir, "ATAC_validation", "atac_validation.R")
  ),
  input_group = "analysis_code",
  required = TRUE
)

df.source.data.versions <- bind_rows(df.analysis.input.files, analysis.script.files) %>%
  filter(required) %>%
  distinct(input_name, input_path, .keep_all = TRUE) %>%
  mutate(
    file_exists = file.exists(input_path),
    file_size_bytes = if_else(file_exists, as.numeric(file.info(input_path)$size), NA_real_),
    modified_time = if_else(file_exists, format(file.info(input_path)$mtime, "%Y-%m-%d %H:%M:%S %Z"), NA_character_),
    sha256 = map_chr(input_path, sha256_file),
    genome_assembly = case_when(
      input_group == "analysis_code" ~ NA_character_,
      input_name %in% c("epd_rn6_bed", "epd_coordinate") ~ "rn6 source; normalized to rn7",
      TRUE ~ "mRatBN7.2/rn7"
    ),
    provenance_note = case_when(
      str_starts(input_name, "hiccups_loop_") ~ "Juicer 1.6; Juicer Tools/HiCCUPS 1.22.01 source BEDPE",
      input_name == "ctcf_motif" ~ "Full FIMO output from 100 PWM models; predicted motif intervals",
      input_name == "atac_peak" ~ "Duttke et al. 2022 rat PFC snATAC peaks normalized to rn7",
      input_name == "ensembl_gtf" ~ "Ensembl Rattus_norvegicus mRatBN7.2 release 113 GTF",
      str_starts(input_name, "epd_") ~ "EPDnew promoter source used for strand-aware rn6-to-rn7 reconstruction",
      input_name == "library_complexity" ~ "Hi-C usable-contact depth QC input",
      input_group == "analysis_code" ~ "Versioned analysis source code",
      TRUE ~ input_group
    ),
    analysis_release = "resubmission-2026-07-28"
  ) %>%
  arrange(input_group, input_name)

assert_analysis_condition(
  all(df.source.data.versions$file_exists) && !any(is.na(df.source.data.versions$sha256)),
  "Source-data version table contains a missing file or checksum."
)

################################################################################
# 16. Write resubmission outputs
################################################################################

# Register only the curated production tables instead of registering every
# intermediate object and filtering the registry afterward.
revised.output.table.registry <- tribble(
  ~output_file, ~object_name,
  "source_data_versions.tsv", "df.source.data.versions",
  "coordinate_system_audit.tsv", "df.coordinate.system.audit",
  "hiccups_sample_loop_quality.tsv.gz", "df.hiccups.sample.loop.quality",
  "hiccups_quality_field_definitions.tsv", "df.hiccups.quality.field.definition",
  "true_tss_generation_summary.tsv", "df.true.tss.summary",
  "promoter_tss_anchor_evidence_definition_summary.tsv", "df.promoter.tss.anchor.evidence.definition.summary",
  "direct_promoter_tss_summary.tsv", "df.direct.promoter.tss.summary",
  "dual_promoter_atac_directional_category_definitions.tsv", "df.dual.promoter.atac.directional.category.definition",
  "dual_promoter_atac_directional_summary.tsv", "df.dual.promoter.atac.directional.summary",
  "revised_atac_support_by_resolution.tsv", "df.revised.atac.support.by.resolution",
  "revised_atac_paired_anchor_mcnemar.tsv", "df.revised.atac.paired.anchor.mcnemar",
  "revised_atac_analysis_definition.tsv", "df.revised.atac.analysis.definition",
  "transcript_position_summary.tsv", "df.transcript.position.summary",
  "transcript_position_definitions.tsv", "df.transcript.position.definition",
  "revised_assignment_pipeline_status.tsv", "df.revised.assignment.pipeline.status",
  "revised_loop_category_definitions.tsv", "df.revised.loop.category.definition",
  "revised_loop_category_summary.tsv", "df.revised.loop.category.summary",
  "revised_loop_category_by_resolution.tsv", "df.revised.loop.category.by.resolution",
  "revised_go_input_summary.tsv", "df.revised.go.input.summary",
  "revised_go_significance_summary.tsv", "df.revised.go.significance.summary",
  "revised_downstream_analysis_definitions.tsv", "df.revised.downstream.analysis.definition",
  "approximate_loop_locus_map.tsv", "df.approximate.loop.locus.map",
  "approximate_loop_locus_summary.tsv", "df.approximate.loop.locus.summary",
  "approximate_loop_locus_method.tsv", "df.approximate.loop.locus.method",
  "approximate_loop_locus_threshold_sensitivity.tsv", "df.approximate.loop.locus.threshold.sensitivity",
  "gene_rank_sensitivity_detail.tsv", "df.gene.rank.sensitivity.detail",
  "gene_rank_sensitivity_correlations.tsv", "df.gene.rank.sensitivity.correlation",
  "figure5_revised_feature_summary_by_resolution.tsv", "df.figure5.revised.feature.summary",
  "loop_sharing_distribution_by_library.tsv", "df.loop.sharing.distribution.by.strain",
  "pairwise_loop_overlap_by_library.tsv", "df.pairwise.loop.overlap.by.strain",
  "hrdp_sample_strain_key.tsv", "df.sample.strain.key",
  "depth_qc_by_strain.tsv", "df.depth.qc.by.strain",
  "depth_qc_summary.tsv", "df.depth.qc.summary",
  "depth_loop_correlation_by_resolution.tsv", "df.depth.loop.correlation.by.resolution",
  "depth_analysis_interpretation.tsv", "df.depth.analysis.interpretation"
)

output.tables <- resolve_output_table_registry(revised.output.table.registry, envir = environment())

# Add the selected final resource and GO-input tables with stable filenames.
selected.resource.table.names <- c("revised_pooled_loop_annotation_resource", "revised_direct_loop_gene_assignments", "revised_putative_regulatory_loops")
output.tables <- c(
  output.tables,
  set_names(revised.resource.tables[selected.resource.table.names], paste0(selected.resource.table.names, ".tsv")),
  set_names(revised.go.gene.sets, paste0("go_input_", names(revised.go.gene.sets), "_genes.tsv"))
)

# GO output is optional when no significant or reportable terms are returned.
if (nrow(df.revised.go.result) > 0L) {
  output.tables[["revised_go_enrichment_BP_all_sets.tsv"]] <- df.revised.go.result
}

assert_analysis_condition(!anyDuplicated(names(output.tables)), "The curated output table registry contains duplicate filenames.")
write_named_tsv_tables(output.tables, output.dir)

# Keep separately generated ATAC matched-null deliverables across production
# reruns and include them in the release manifest whenever they are present.
atac.matched.null.output.files <- c(
  "revised_atac_threshold_sensitivity_by_resolution.tsv", "revised_atac_matched_null_summary.tsv",
  "revised_atac_matched_null_permutation.tsv", "revised_atac_matched_null_method.tsv",
  "revised_atac_matched_control_quality.tsv", "revised_atac_matched_control_strata.tsv",
  "revised_atac_matched_null_run_metadata.tsv", "revised_atac_random_relocation_null_dumbbell_by_resolution.pdf",
  "revised_atac_random_relocation_null_dumbbell_by_resolution.png", "revised_atac_matched_null_session_info.txt"
)
available.atac.matched.null.output.files <- intersect(atac.matched.null.output.files, list.files(output.dir, all.files = FALSE, no.. = TRUE))

writeLines(str_replace(capture.output(sessionInfo()), "\\s+$", ""), con = file.path(output.dir, "resubmit_session_info.txt"))

generated.output.files <- c(names(output.tables), revised.figure.output.files, available.atac.matched.null.output.files, "resubmit_session_info.txt")
assert_analysis_condition(!anyDuplicated(generated.output.files), "The resubmission output manifest contains duplicate filenames.")

# Remove stale legacy, duplicate, and intermediate files from previous runs so
# the directory contains only the curated deliverables and provenance records.
protected.output.files <- c(generated.output.files, "resubmit_output_manifest.tsv")
stale.output.files <- setdiff(list.files(output.dir, all.files = FALSE, no.. = TRUE), protected.output.files)
if (length(stale.output.files) > 0L) {
  stale.output.removed <- file.remove(file.path(output.dir, stale.output.files))
  assert_analysis_condition(all(stale.output.removed), "One or more stale resubmission output files could not be removed.")
}

df.output.manifest <- tibble(
  output_file = generated.output.files,
  output_path = file.path(output.dir, generated.output.files),
  file_exists = file.exists(output_path),
  file_size_bytes = as.numeric(file.info(output_path)$size),
  sha256 = map_chr(output_path, sha256_file),
  generated_by = if_else(
    str_starts(output_file, "revised_atac_matched") | str_starts(output_file, "revised_atac_threshold") | str_starts(output_file, "revised_atac_random"),
    "atac_validation.R",
    "promoter_enhancer_interaction_resubmit.R"
  ),
  analysis_release = "resubmission-2026-07-28",
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
)

assert_analysis_condition(all(df.output.manifest$file_exists) && !any(is.na(df.output.manifest$sha256)), "Output manifest contains a missing file or checksum.")
readr::write_tsv(df.output.manifest, file.path(output.dir, "resubmit_output_manifest.tsv"))

################################################################################
# 17. Final checks and console summary
################################################################################

message("\nWrote outputs to: ", output.dir)
message("Pooled loop resource (<2 Mb): ", nrow(df.loop.universe))
message("Revised direct promoter/TSS-supported loops: ", sum(df.direct.promoter.tss.loop.summary$has_any_direct_promoter_tss))
message("Revised direct loop-anchor-gene assignments: ", nrow(df.direct.promoter.tss.gene.assignment))
message("Revised proximal promoter/TSS-supported loops (1-200 kb): ", sum(df.proximal.promoter.tss.loop.summary$has_any_proximal_promoter_tss))
message("Revised proximal loop-anchor-gene assignments: ", nrow(df.proximal.promoter.tss.gene.assignment))
message("Revised secondary inward <=10-kb supported loops: ", sum(df.proximal.promoter.tss.loop.summary$has_any_secondary_inward_proximal_10kb))
message("Revised secondary inward <=10-kb loop-anchor-gene assignments: ", nrow(df.secondary.inward.proximal.gene.assignment))
message("Revised direct/proximal assignment, ATAC, transcript-position, and loop-category sections are complete. Legacy analyses are excluded from this production workflow.")

message("\nRevised true-TSS-excluded ATAC support by resolution:")
print(df.revised.atac.support.by.resolution)
message("\nRevised paired-anchor ATAC comparisons:")
print(df.revised.atac.paired.anchor.mcnemar)
message("\nRevised per-anchor vs. fragment-level ATAC agreement:")
print(df.revised.atac.method.comparison.summary)
message("\nRevised transcript-position summary:")
print(df.transcript.position.summary)
message("\nRevised mutually exclusive loop categories:")
print(df.revised.loop.category.summary)
message("\nPromoter/TSS categories split by anchor-site multiplicity:")
print(df.promoter.tss.exclusive.loop.category.summary)
message("\nDual-primary representative-anchor sensitivity summary:")
print(df.dual.primary.anchor.representative.summary)
message("\nRevised gene-loop threshold summary:")
print(df.revised.gene.count.threshold.summary)
message("\nApproach-2 multiple-interaction gene summary:")
print(df.approach2.multiple.interaction.summary)
message("\nRevised GO input and significance summary:")
print(df.revised.go.significance.summary)

message("\nDepth sensitivity: loop count vs. Hi-C contacts")
print(df.depth.loop.correlation.by.resolution %>% filter(depth_metric == "Hi-C_Contacts") %>% dplyr::select(resolution, n_strains, pearson_r, pearson_p, spearman_rho, spearman_p))

message("\nGenetic-distance analysis status:")
print(df.genetic.distance.status)
