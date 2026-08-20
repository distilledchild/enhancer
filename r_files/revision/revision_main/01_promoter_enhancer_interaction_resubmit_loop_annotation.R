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
coord.prep.script <- file.path(analysis.dir, "00_promoter_enhancer_interaction_resubmit_coord_prep.R")

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

