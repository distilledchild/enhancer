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
library("RIdeogram")
library("rsvg")
library("magick")
library("cowplot")

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
# To rebuild the cache, run 00_promoter_enhancer_interaction_resubmit_coord_prep.R directly.
# (Or uncomment the rebuild block below and set Sys.setenv(REBUILD_COORD_CACHE = "1"))
################################################################################

# Optional cache reconstruction trigger (commented out for execution safety):
# rebuild.coord.cache <- identical(Sys.getenv("REBUILD_COORD_CACHE", unset = "0"), "1")
# if (rebuild.coord.cache) {
#   message("REBUILD_COORD_CACHE=1 detected. Reconstructing coordinate cache...")
#   source(coord.prep.script, local = FALSE)
# }

coord.cache.object.names <- coordinate_cache_object_names()
expected.cache.files <- file.path(coord.cache.dir, paste0(coord.cache.object.names, ".rds"))

# Verify that all expected cache .rds files exist
missing.cache.files <- expected.cache.files[!file.exists(expected.cache.files)]
if (length(missing.cache.files) > 0L) {
  stop("Missing required cache .rds file(s):\n", paste(missing.cache.files, collapse = "\n"), call. = FALSE)
}

# Load all cached objects into environment
load_coordinate_cache_objects(cache.dir = coord.cache.dir, object.names = coord.cache.object.names, envir = environment())

# Extract downstream input file paths
input.file.paths <- setNames(df.analysis.input.files$input_path, df.analysis.input.files$input_name)
library.complexity.file <- input.file.paths[["library_complexity"]]
message("Successfully verified and loaded ", length(coord.cache.object.names), " cache objects from: ", coord.cache.dir)
for (cache.obj.name in coord.cache.object.names) {
  obj <- get(cache.obj.name, envir = environment())
  obj.desc <- if (is(obj, "GRanges")) {
    paste0("GRanges (", scales::comma(length(obj)), " ranges)")
  } else if (is.data.frame(obj)) {
    paste0("tibble (", scales::comma(nrow(obj)), " rows, ", ncol(obj), " cols)")
  } else {
    class(obj)[1]
  }
  cat(sprintf("  - %-35s : %s\n", cache.obj.name, obj.desc))
}

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
# 3. Define primary DIRECT promoter/TSS-anchor overlap tier (+/-1-kb TSS window)
#
# Description: Map loop anchors to promoter/TSS using the primary promoter window
#              defined as TSS +/-1 kb (2,001 bp) from Ensembl or EPD.
# Input:       df.loop.distinct.2mb, df.true.tss.transcript, gr.true.tss, df.promoter.epd.rn7.1based
# Output:      list.primary.direct.tier, df.direct.promoter.tss.anchor.overlap, df.direct.promoter.tss.gene.assignment
################################################################################

list.gr.loop.anchor.by.side <- create_loop_anchor_granges_by_side(df.loop.distinct.2mb)
names(list.gr.loop.anchor.by.side) # [1] "anchor1" "anchor2"
gr.epd.tss <- create_epd_tss_granges(df.promoter.epd.rn7.1based)

# func to manupulate coordinate for window (here, 1000L) in GR objects; 1. reuse gr.true.tss, 2 for tss around at the end of chr, 3. integrity for metadata & index (id), 4. easy for sensitivity test
promoter.window.flank.bp <- 1000L
gr.true.tss.promoter.window.1kb <- expand_tss_to_promoter_windows(gr.true.tss, flank.bp = promoter.window.flank.bp)
gr.epd.tss.promoter.window.1kb <- expand_tss_to_promoter_windows(gr.epd.tss, flank.bp = promoter.window.flank.bp)

# Maps loop anchors to Ensembl/EPD +/-1-kb promoter windows and returns a comprehensive list of overlap details and summaries
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
names(list.primary.direct.tier)

# Unpack lookup and primary tier outputs
df.true.tss.lookup <- list.primary.direct.tier$true_tss_lookup # tibble (54,993 rows) - Ensembl transcript metadata lookup table (gene/transcript IDs, biotypes)
df.epd.promoter.lookup <- list.primary.direct.tier$epd_promoter_lookup # tibble (12,524 rows) - EPD promoter metadata lookup table (promoter IDs, gene names)
df.direct.true.tss.anchor.overlap <- list.primary.direct.tier$true_tss_overlap # tibble (43,267 rows) - Detailed 1:1 overlap records between loop anchors and Ensembl TSS +/-1-kb windows
df.direct.epd.promoter.anchor.overlap <- list.primary.direct.tier$epd_promoter_overlap # tibble (15,920 rows) - Detailed 1:1 overlap records between loop anchors and EPD promoter +/-1-kb windows
df.direct.promoter.tss.anchor.overlap <- list.primary.direct.tier$combined_overlap # tibble (59,187 rows) - Unified Ensembl + EPD direct promoter overlap evidence table
list.primary.direct.summary <- list.primary.direct.tier$summary # list (5 tibbles)     - Multi-level summary tables: gene_assignment, anchor_summary, loop_summary, etc.

df.direct.promoter.tss.gene.assignment <- list.primary.direct.summary$gene_assignment # tibble (26,920 rows) - Loop-anchor to unique gene assignment table (deduplicated across isoforms)
df.loop.anchor.index <- list.primary.direct.summary$anchor_index # tibble (62,042 rows) - Base index of all 62,042 loop anchors (anchor1 & anchor2 across 31,021 loops)
df.direct.promoter.tss.anchor.count <- list.primary.direct.summary$anchor_count # tibble (20,391 rows) - Promoter overlap counts per anchor (for anchors with >=1 promoter)
df.direct.promoter.tss.anchor.summary <- list.primary.direct.summary$anchor_summary # tibble (62,042 rows) - Complete anchor-level promoter summary for all 62,042 anchors
df.direct.promoter.tss.anchor.wide <- list.primary.direct.summary$anchor_wide # tibble (31,021 rows) - Anchor1 vs anchor2 promoter status pivoted to wide format per loop
df.direct.promoter.tss.loop.gene.count <- list.primary.direct.summary$loop_gene_count # tibble (16,343 rows) - Count of unique promoter genes assigned per loop (for loops with >=1 gene)
df.direct.promoter.tss.loop.summary <- list.primary.direct.summary$loop_summary # tibble (31,021 rows) - Comprehensive loop-level summary with promoter counts/flags for all 31,021 loops
df.direct.promoter.tss.summary <- list.primary.direct.summary$summary # tibble (8 rows)      - High-level QC summary metric counts for the primary direct tier

# Define coordinate rules and interpretations for promoter-anchor evidence tiers
df.promoter.anchor.assignment.definitions <- tribble(
  ~evidence_tier, ~coordinate_rule, ~interpretation,
  "primary_direct",
  paste0("Direct anchor overlap with a promoter window defined as TSS +/-1 kb from ", "either Ensembl transcript TSS or EPD TSS."),
  paste0("Primary promoter-associated anchor evidence used by downstream revised ", "ATAC, category, gene-resource, and GO analyses."),
  "secondary_inward_proximal",
  paste0("Non-overlapping TSS or EPD promoter interval within 10 kb of an anchor and ", "fully located in the inter-anchor interval on that anchor's side of the ", "loop midpoint. The assigned anchor must be strictly closer than the ", "opposite anchor; same-anchor/gene primary direct assignments are excluded."),
  paste0("Secondary proximity-supported candidate tier; weaker than direct overlap ", "and not treated as validated target-gene evidence."),
  "exploratory_proximal_catalog",
  paste0("Any non-overlapping TSS or EPD promoter interval 1-200 kb from an anchor in ", "either direction."),
  paste0("Exploratory sensitivity catalog only; unsuitable for strong regulatory or ", "target-gene claims.")
)
df.promoter.anchor.assignment.definitions
# Print direct promoter/TSS assignment summary message
message(
  "Direct promoter/TSS overlap retained ",
  nrow(df.direct.promoter.tss.gene.assignment),
  " loop-anchor-gene assignments across ",
  sum(df.direct.promoter.tss.loop.summary$has_any_direct_promoter_tss),
  " pooled loops."
)
# anchor-based data object (only assigned: df.direct.promoter.tss.anchor.overlap, all anchors: df.direct.promoter.tss.anchor.summary)
df.direct.promoter.tss.anchor.overlap %>%
  count() # assigned anchors ONLY: 55755
# head(2)
df.direct.promoter.tss.anchor.summary %>%
  head(2)
# count() # all anchors: 62042

# ------------------------------------------------------------------------------
# Loop-level promoter assignment summary derived from anchor summary:
# 1. loops not assigned with TSS/promoter
# 2. loops assigned with TSS/promoter
#    2-1. both anchors assigned with TSS/promoter
#    2-2. one anchor only assigned with TSS/promoter
# ------------------------------------------------------------------------------
df.direct.promoter.tss.anchor.summary.stats <- df.direct.promoter.tss.anchor.summary %>%
  group_by(loop_id) %>%
  summarise(
    n_assigned_anchors = sum(has_any_direct_promoter_tss),
    .groups = "drop"
  ) %>%
  mutate(
    promoter_assignment_category = case_when(
      n_assigned_anchors == 2L ~ "2-1. direct_both_anchors (both anchors direct)",
      n_assigned_anchors == 1L ~ "2-2. direct_one_anchor_only (one anchor only direct)",
      TRUE ~ "1. no_direct_promoter (no direct promoter/TSS assignment)"
    )
  ) %>%
  count(promoter_assignment_category, name = "n_loops") %>%
  mutate(pct_pooled_loops = round(100 * n_loops / sum(n_loops), 1))

df.direct.promoter.tss.anchor.summary.stats
#  promoter_assignment_category                               n_loops pct_pooled_loops
# 1 1. no_direct_promoter (no direct promoter/TSS assignment)   14678            47.3
# 2 2-1. direct_both_anchors (both anchors direct)               4048            13
# 3 2-2. direct_one_anchor_only (one anchor only direct)        12295            39.6

# loop-based data object
df.direct.promoter.tss.loop.summary %>%
  # count() # 31021
  head(2)

################################################################################
# 4. Note on Proximal Promoter/TSS Allocation Sensitivity Analysis
#
# Distance-based proximal search (1 bp - 200 kb) and secondary inward candidate
# tiers are isolated in:
#   01_2_promoter_enhancer_interaction_resubmit_loop_annotation_SENSITIVITY_proximal_distance.R
#
# Downstream regulatory loop categorization strictly relies on the Primary Direct
# tier (+/-1-kb TSS overlap) established in Section 3.
################################################################################
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

#######################################################################################
# 5-1. Process Ensemble TSS & EPD TSS
#######################################################################################
atac.minimum.overlap.bp <- 50L
tss.exclusion.flank.bp <- 1000L

# 5-1-1. Combine every transcript-level Ensembl TSS with every EPD TSS.
df.known.tss.source.record <- bind_rows( # 67,517
  df.true.tss.transcript %>% # 54,993
    transmute(
      tss_source = "Ensembl_transcript_true_TSS",
      tss_source_id = true_tss_id,
      chr,
      tss_start = true_tss_start,
      tss_end = true_tss_end
    ),
  df.promoter.epd.rn7.1based %>% # 12,524
    transmute(
      tss_source = "EPDnew_TSS",
      tss_source_id = promoter_annotation_id,
      chr,
      tss_start = epd_tss_start,
      tss_end = epd_tss_end
    )
)

# Check that every known-TSS exclusion source is represented by ONE BP.
assert_analysis_condition(
  !any(
    df.known.tss.source.record$tss_start !=
      df.known.tss.source.record$tss_end
  ),
  "A known TSS record is not a one-base interval."
)

# 5-1-2. Deduplicate by TSS coordinates & Create grange MERGED object of 2000bp-WINDOWED TSS points
df.known.tss.point <- df.known.tss.source.record %>%
  group_by(chr, tss_start, tss_end) %>% # group by identical genomic coordinates
  summarise(
    n_tss_source_records = n(),
    tss_sources = str_c(sort(unique(tss_source)), collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(chr, tss_start)

df.known.tss.point %>% # count() # 56,881/67,517
  head(2)

# chr   tss_start tss_end n_tss_source_records tss_sources
# 1 chr1      31564   31564                    1 Ensembl_transcript_true_TSS
# 2 chr1      41635   41635                    1 Ensembl_transcript_true_TSS

# grange object: gr.known.tss.point
gr.known.tss.point <- GRanges(
  seqnames = df.known.tss.point$chr,
  ranges = IRanges(
    start = df.known.tss.point$tss_start,
    end = df.known.tss.point$tss_end
  ),
  n_tss_source_records = df.known.tss.point$n_tss_source_records,
  tss_sources = df.known.tss.point$tss_sources
)

# add window with symmetric 2,001-bp interval centered on each one-base TSS and MERGE (reduce())
gr.known.tss.exclusion <- GRanges(
  seqnames = seqnames(gr.known.tss.point), # chr info
  ranges = IRanges(
    start = pmax(1L, start(gr.known.tss.point) - tss.exclusion.flank.bp), # -1000bp from TSS start (pmax)
    end = end(gr.known.tss.point) + tss.exclusion.flank.bp # +1000bp from TSS end
  )
) %>%
  GenomicRanges::reduce(ignore.strand = TRUE) # STRAND IGNORED!

#######################################################################################
# 5-2. Process ATAC-seq peaks
#######################################################################################
# 5-2-1. merging ATAC-seq peaks
gr.atac.union <- GenomicRanges::reduce(
  gr.atac.rn7.1based,
  ignore.strand = TRUE
)
# 5-2-2. intersect with TSS exclusion regions
# 5-2-2-1. extracting common seqlevels between ATAC and TSS exclusion regions
common.seqlevels.atac.tss <- intersect(
  seqlevels(gr.atac.union),
  seqlevels(gr.known.tss.exclusion)
)
# 5-2-2-2. keep only common seqlevels in ATAC-seq peaks
gr.atac.union.common <- keepSeqlevels(
  gr.atac.union,
  common.seqlevels.atac.tss,
  pruning.mode = "coarse"
)
# 5-2-2-3. keep only common seqlevels in TSS exclusion regions
gr.known.tss.exclusion.for.atac <- keepSeqlevels(
  gr.known.tss.exclusion,
  common.seqlevels.atac.tss,
  pruning.mode = "coarse"
)
# 5-2-3. gr.atac.union.common - gr.known.tss.exclusion.for.atac
gr.atac.non.tss <- GenomicRanges::setdiff(
  gr.atac.union.common,
  gr.known.tss.exclusion.for.atac,
  ignore.strand = TRUE
)

# convert gr(gr.atac.non.tss) to df(df.atac.true.tss.exclusion.region)
df.atac.true.tss.exclusion.region <- tibble(
  chr = as.character(seqnames(gr.known.tss.exclusion)),
  start = start(gr.known.tss.exclusion),
  end = end(gr.known.tss.exclusion),
  width_bp = width(gr.known.tss.exclusion),
  exclusion_definition = "combined_Ensembl_and_EPD_true_TSS_plus_minus_1kb"
)

#######################################################################################
# 5-2. summary
#######################################################################################
df.atac.exclusion.summary <- tibble(
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
    length(gr.atac.union.common),
    sum(width(gr.atac.union.common)),
    nrow(df.true.tss.transcript),
    nrow(df.promoter.epd.rn7.1based),
    nrow(df.known.tss.point),
    length(gr.known.tss.exclusion),
    sum(width(gr.known.tss.exclusion)),
    length(gr.atac.non.tss),
    sum(width(gr.atac.non.tss)),
    sum(width(gr.atac.union.common)) -
      sum(width(gr.atac.non.tss))
  )
)

#######################################################################################
# 5-3. [Promoter Level] Obtain Promoters Supported by ATAC-seq Peaks
#   5-3-1. Extract distinct promoter +/-1kb windows overlapping with loop anchors
#   5-3-2. Quantify ATAC peak overlap at each promoter window (>=50 bp threshold)
#   5-3-3. Summarize promoter accessibility (ATAC support) at the loop-anchor level
#######################################################################################

# 5-3-1. Extract distinct promoter +/-1kb windows overlapping with loop anchors
df.direct.promoter.tss.anchor.overlap %>% # 55,755 : excluding non-assigned anchors & including multi-assigned anchors
  count()
head(2)

gr.promoter.window.distinct.record <- df.direct.promoter.tss.anchor.overlap %>% # 55,755
  distinct(loop_id, resolution, anchor_side, annotation_chr, annotation_start, annotation_end) %>% # - 9,489: due to isoform, the same TSS between ENSEMBL and EPD
  arrange(loop_id, anchor_side, annotation_chr, annotation_start, annotation_end)

gr.promoter.window.atac.record <- GRanges( # 46,266
  seqnames = gr.promoter.window.distinct.record$annotation_chr,
  ranges = IRanges(start = gr.promoter.window.distinct.record$annotation_start, end = gr.promoter.window.distinct.record$annotation_end),
  loop_id = gr.promoter.window.distinct.record$loop_id,
  resolution = gr.promoter.window.distinct.record$resolution,
  anchor_side = gr.promoter.window.distinct.record$anchor_side
)

# 5-3-2. Quantify ATAC peak overlap at each promoter window (>=50 bp threshold)
df.promoter.window.atac.overlap <- summarise_anchor_interval_overlap(
  gr.promoter.window.atac.record,
  gr.atac.union.common,
  minimum.overlap.bp = atac.minimum.overlap.bp
)

# 5-3-3. Summarize promoter-window ATAC accessibility at the loop-anchor level
df.promoter.window.atac.anchor.summary <- df.promoter.window.atac.overlap %>%
  group_by(loop_id, resolution, anchor_side) %>%
  summarise(
    n_direct_promoter_windows = n(), # number of direct promoter windows
    n_atac_positive_promoter_windows = sum(has_feature_overlap_ge_minimum), # number of direct promoter windows with ATAC signals
    promoter_window_atac_overlap_any = any(has_feature_overlap_any), # any overlap (>= 1bp) between promoter window and ATAC peaks
    promoter_window_atac_overlap_ge50 = any(has_feature_overlap_ge_minimum), # overlap >= 50 bp between promoter window and ATAC peaks
    promoter_window_atac_max_overlap_bp = max(feature_overlap_bp), # maximum overlap between promoter window and ATAC peaks (by bp)
    .groups = "drop"
  ) %>%
  arrange(loop_id, anchor_side)

# Confirm that promoter-window ATAC evidence covers every directly annotated
# anchor once (20,391 anchors) and does not create records for anchors lacking promoter evidence (41,651 anchors).
# [Validation interpretation when TRUE]:
# 1. Matches exactly all 20,391 direct promoter-assigned anchors (no omission / no excess).
# 2. Guarantees strictly one summary row per direct promoter anchor without duplicate entries.
assert_analysis_condition(
  condition = nrow(df.promoter.window.atac.anchor.summary) == sum(df.direct.promoter.tss.anchor.summary$has_any_direct_promoter_tss) &&
    !anyDuplicated(str_c(df.promoter.window.atac.anchor.summary$loop_id, df.promoter.window.atac.anchor.summary$anchor_side, sep = "|")),
  message = "Promoter-window ATAC summary is not one row per direct promoter anchor.",
  success.message = sprintf(
    "Verified: Promoter-window ATAC summary covers all %d direct promoter anchors (1:1 match, 0 duplicate).",
    nrow(df.promoter.window.atac.anchor.summary)
  )
)

#######################################################################################
# 5-4. [Anchor Level] Obtain Loop Anchors Supported by ATAC-seq Peaks & Integrate Evidence
#   5-4-1. Combine all 62,042 loop anchors across 31,021 loops
#   5-4-2. Quantify raw ATAC, non-TSS ATAC, and TSS-exclusion overlaps per anchor
#   5-4-3. Assemble comprehensive master anchor evidence table (df.atac.anchor.evidence)
#######################################################################################

# 5-4-1. Combine all 62,042 loop anchors across 31,021 loops
gr.loop.anchor.all <- combine_loop_anchor_granges_by_side(list.gr.loop.anchor.by.side)

# 5-4-2. Quantify raw ATAC, non-TSS ATAC, and TSS-exclusion overlaps per anchor
df.raw.atac.anchor.overlap <- summarise_anchor_interval_overlap(
  gr.loop.anchor.all,
  gr.atac.union.common,
  minimum.overlap.bp = atac.minimum.overlap.bp
)
df.non.tss.atac.anchor.overlap <- summarise_anchor_interval_overlap(
  gr.loop.anchor.all,
  gr.atac.non.tss,
  minimum.overlap.bp = atac.minimum.overlap.bp
)
df.tss.exclusion.anchor.overlap <- summarise_anchor_interval_overlap(
  gr.loop.anchor.all,
  gr.known.tss.exclusion,
  minimum.overlap.bp = atac.minimum.overlap.bp
)

# 5-4-3. Assemble comprehensive master anchor evidence table
df.atac.anchor.evidence <- df.raw.atac.anchor.overlap %>%
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
    df.non.tss.atac.anchor.overlap %>%
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
    df.tss.exclusion.anchor.overlap %>%
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
    df.promoter.window.atac.anchor.summary,
    by = c("loop_id", "resolution", "anchor_side")
  ) %>%
  mutate(
    n_direct_promoter_windows = coalesce(n_direct_promoter_windows, 0L),
    n_atac_positive_promoter_windows = coalesce(n_atac_positive_promoter_windows, 0L),
    promoter_window_atac_overlap_any = coalesce(promoter_window_atac_overlap_any, FALSE),
    promoter_window_atac_overlap_ge50 = coalesce(promoter_window_atac_overlap_ge50, FALSE),
    promoter_window_atac_max_overlap_bp = coalesce(promoter_window_atac_max_overlap_bp, 0L)
  ) %>%
  arrange(loop_id, anchor_side)

# Validate two unique ATAC-evidence anchor rows and nonnegative residual bases.
assert_analysis_condition(
  condition = nrow(df.atac.anchor.evidence) == 2L * nrow(df.loop.distinct.2mb) &&
    !anyDuplicated(str_c(df.atac.anchor.evidence$loop_id, df.atac.anchor.evidence$anchor_side, sep = "|")) &&
    !any(df.atac.anchor.evidence$non_tss_anchor_bp < 0L),
  message = "ATAC anchor evidence did not preserve two valid anchors per loop.",
  success.message = sprintf(
    "Verified: ATAC anchor evidence covers all %d loop anchors (2 anchors per loop, 0 duplicate, all residual bp >= 0).",
    nrow(df.atac.anchor.evidence)
  )
)

#######################################################################################
# 5-5. Establish Directional Loop Orientations (Promoter to Candidate Regulatory Anchor)
#   5-5-1. Format promoter-anchor and candidate-anchor metric profiles
#   5-5-2. Subclassify dual-promoter loops (contacts with direct promoters at both anchors)
#   5-5-3. Isolate unambiguous single-promoter loops (12,295 candidate regulatory loops)
#   5-5-4. Reshape anchor information to pooled loop summary (df.atac.loop.summary)
#######################################################################################

# 5-5-1. Format promoter- and candidate-anchor profiles to pair directional evidence
df.atac.promoter.anchor.metric <- df.atac.anchor.evidence %>% # 5-4-3
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

df.atac.candidate.enhancer.anchor.metric <- df.atac.anchor.evidence %>%
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

df.atac.direct.orientation <- df.direct.promoter.tss.anchor.summary %>%
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
  left_join(df.atac.promoter.anchor.metric, by = c("loop_id", "promoter_anchor_side")) %>%
  left_join(df.atac.candidate.enhancer.anchor.metric, by = c("loop_id", "candidate_anchor_side")) %>%
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

# 5-5-2. Subclassify dual-promoter loops by directional promoter and residual non-TSS ATAC support
df.dual.promoter.atac.directional.classification <- df.atac.direct.orientation %>%
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

# Record the four mutually exclusive dual-promoter definitions explicitly
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

assert_analysis_condition(
  condition = nrow(df.dual.promoter.atac.directional.classification) == sum(df.direct.promoter.tss.loop.summary$n_direct_anchor_sides == 2L) &&
    sum(df.dual.promoter.atac.directional.summary$n_loops) == nrow(df.dual.promoter.atac.directional.classification) &&
    !anyDuplicated(df.dual.promoter.atac.directional.classification$loop_id),
  message = "Dual-promoter ATAC classification did not preserve all eligible loops.",
  success.message = sprintf(
    "Verified: Dual-promoter ATAC classification covers all %d dual-promoter loops (1:1 match, 0 duplicate, 100%% categorized).",
    nrow(df.dual.promoter.atac.directional.classification)
  )
)

# 5-5-3. Isolate unambiguous single-promoter loops (12,295 candidate regulatory loops)
df.atac.unambiguous.orientation <- df.atac.direct.orientation %>%
  filter(is_unambiguous_candidate_regulatory_orientation)

df.atac.unambiguous.orientation %>% head(2)

assert_analysis_condition(
  condition = nrow(df.atac.unambiguous.orientation) == sum(df.direct.promoter.tss.loop.summary$n_direct_anchor_sides == 1L) &&
    !anyDuplicated(df.atac.unambiguous.orientation$loop_id),
  message = "Single-direct-anchor ATAC orientations are not one row per eligible loop.",
  success.message = sprintf(
    "Verified: Single-direct-anchor ATAC orientation covers all %d candidate regulatory loops (1:1 match, 0 duplicate).",
    nrow(df.atac.unambiguous.orientation)
  )
)

# 5-5-4. Reshape anchor information to pooled loop summary (df.atac.loop.summary)
df.atac.anchor.wide <- df.atac.anchor.evidence %>%
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

df.atac.loop.summary <- df.direct.promoter.tss.loop.summary %>%
  dplyr::select(loop_id, resolution, n_direct_anchor_sides, has_any_direct_promoter_tss, has_direct_promoter_tss_both_anchors, direct_anchor_assignment_class) %>%
  left_join(df.atac.anchor.wide, by = "loop_id") %>%
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

assert_analysis_row_count(
  df = df.atac.loop.summary,
  expected = nrow(df.loop.distinct.2mb),
  message = "ATAC loop summary lost pooled loops.",
  success.message = sprintf(
    "Verified: ATAC loop summary retains all %d pooled loops (1:1 match).",
    nrow(df.atac.loop.summary)
  )
)

#######################################################################################
# 5-6. Statistical & Resolution-based ATAC Support Analysis
#   5-6-1. Quantify raw and non-TSS ATAC support across HiCCUPS resolutions (5K, 10K, 25K, ALL)
#   5-6-2. Perform paired-anchor McNemar tests (Promoter vs Candidate accessibility comparison)
#######################################################################################

# 5-6-1. Quantify raw and non-TSS ATAC support across HiCCUPS resolutions (5K, 10K, 25K, ALL)
df.atac.support.by.resolution <- bind_rows(
  df.atac.unambiguous.orientation,
  df.atac.unambiguous.orientation %>% mutate(resolution = "ALL")
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

# 5-6-2. Perform paired-anchor McNemar tests (Promoter vs Candidate accessibility comparison)
# McNemar tests use paired anchors from the same loop. Raw ATAC applies the same
# peak definition to both sides; the second comparison applies true-TSS-excluded
# ATAC to both sides and is therefore also symmetric.
df.atac.paired.anchor.mcnemar <- map_dfr(
  c(sort(unique(df.atac.unambiguous.orientation$resolution)), "ALL"),
  function(resolution.i) {
    df.resolution <- if (resolution.i == "ALL") {
      df.atac.unambiguous.orientation
    } else {
      df.atac.unambiguous.orientation %>%
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

#######################################################################################
# 5-7. Sensitivity Analysis: Fragment-level TSS Subtraction Validation
#   5-7-1. Subtract TSS exclusion windows from candidate anchors to generate residual fragments
#   5-7-2. Measure non-TSS ATAC overlaps across individual fragments
#   5-7-3. Cross-validate whole-anchor vs fragment-level ATAC consistency
#######################################################################################

# 5-7-1. Subtract TSS exclusion windows from candidate anchors to generate residual fragments
# Repeat candidate-anchor testing after physically subtracting TSS regions from
# each anchor. This fragment-level sensitivity analysis preserves anchor identity
# and can be checked against the primary whole-anchor/non-TSS-signal method.
df.atac.unambiguous.candidate.anchor <- df.atac.unambiguous.orientation %>%
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

gr.atac.unambiguous.candidate.anchor <- GRanges(
  seqnames = df.atac.unambiguous.candidate.anchor$anchor_chr,
  ranges = IRanges(start = df.atac.unambiguous.candidate.anchor$anchor_start, end = df.atac.unambiguous.candidate.anchor$anchor_end),
  loop_id = df.atac.unambiguous.candidate.anchor$loop_id,
  resolution = df.atac.unambiguous.candidate.anchor$resolution,
  anchor_side = df.atac.unambiguous.candidate.anchor$anchor_side
)
# take a long time
df.atac.candidate.non.tss.fragment <- subtract_exclusion_from_anchor_ranges(gr.atac.unambiguous.candidate.anchor, gr.known.tss.exclusion) %>%
  mutate(fragment_row_id = row_number(), .before = 1)

gr.atac.candidate.non.tss.fragment <- GRanges(
  seqnames = df.atac.candidate.non.tss.fragment$fragment_chr,
  ranges = IRanges(start = df.atac.candidate.non.tss.fragment$fragment_start, end = df.atac.candidate.non.tss.fragment$fragment_end),
  loop_id = df.atac.candidate.non.tss.fragment$loop_id,
  resolution = df.atac.candidate.non.tss.fragment$resolution,
  anchor_side = df.atac.candidate.non.tss.fragment$anchor_side
)

# 5-7-2. Measure non-TSS ATAC overlaps across individual fragments
df.atac.fragment.overlap <- summarise_anchor_interval_overlap(
  gr.atac.candidate.non.tss.fragment,
  gr.atac.non.tss,
  minimum.overlap.bp = atac.minimum.overlap.bp
)

df.atac.candidate.non.tss.fragment.evidence <- df.atac.candidate.non.tss.fragment %>%
  left_join(
    df.atac.fragment.overlap %>%
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

df.atac.fragment.anchor.observed <- df.atac.candidate.non.tss.fragment.evidence %>%
  group_by(loop_id, resolution, anchor_side) %>%
  summarise(
    n_non_tss_anchor_fragments = n(),
    fragment_non_tss_anchor_bp = sum(fragment_width_bp),
    fragment_non_tss_atac_overlap_bp = sum(fragment_atac_overlap_bp),
    fragment_non_tss_atac_any = any(fragment_atac_overlap_any),
    fragment_non_tss_atac_ge50 = any(fragment_atac_overlap_ge50),
    .groups = "drop"
  )

# 5-7-3. Cross-validate whole-anchor vs fragment-level ATAC consistency
df.atac.fragment.anchor.summary <- df.atac.unambiguous.candidate.anchor %>%
  left_join(df.atac.fragment.anchor.observed, by = c("loop_id", "resolution", "anchor_side")) %>%
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

assert_analysis_condition(
  condition = !any(df.atac.fragment.anchor.summary$per_anchor_non_tss_bp != df.atac.fragment.anchor.summary$fragment_non_tss_anchor_bp),
  message = "Fragment-level TSS subtraction disagrees with per-anchor residual bases.",
  success.message = sprintf(
    "Verified: Fragment-level TSS subtraction matches per-anchor residual bases across all %d candidate anchors (0 bp mismatch).",
    nrow(df.atac.fragment.anchor.summary)
  )
)

df.atac.method.comparison.summary <- bind_rows(
  df.atac.fragment.anchor.summary,
  df.atac.fragment.anchor.summary %>% mutate(resolution = "ALL")
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

df.atac.method.comparison.summary
#######################################################################################
# 5-8. Document ATAC Analysis Operational Definitions
#######################################################################################
df.atac.analysis.definition <- tribble(
  ~analysis_item, ~definition,
  "ATAC_interpretation", paste0("ATAC overlap supports open chromatin only; it does not validate enhancer ", "function or a promoter-enhancer interaction."),
  "known_TSS_exclusion", paste0("Union of strand-aware Ensembl transcript TSS and lifted EPD TSS, expanded ", "symmetrically by +/-1 kb in rn7 coordinates."),
  "primary_directional_set", paste0("Loops with primary TSS +/-1-kb promoter-window evidence at exactly one ", "anchor; the opposite anchor is evaluated as a candidate regulatory anchor."),
  "both_direct_anchors", paste0("Loops with primary TSS +/-1-kb promoter-window evidence at both anchors ", "are retained as a promoter-promoter-compatible class without forced ", "direction."),
  "minimum_overlap", paste0("Any-bp overlap is reported for sensitivity; >=50 bp overlap with one ", "reduced ATAC interval is the primary robust-support flag."),
  "fragment_sensitivity", paste0("TSS regions are subtracted independently from each candidate anchor and ", "the residual fragments are tested without losing original anchor identity."),
  "paired_test_interpretation", paste0("McNemar tests compare paired accessibility states conditional on selection ", "of loops with one direct promoter/TSS anchor; they do not validate the ", "opposite anchor as a functional enhancer."),
  "proximal_catalog", paste0("ATAC is available in the complete anchor evidence table, but the full ", "1-200 kb proximal catalog is not promoted to primary directional evidence.")
)

message(
  "Non-TSS ATAC support (>=50 bp) was found at the opposite anchor in ",
  sum(df.atac.unambiguous.orientation$candidate_non_tss_atac_overlap_ge50),
  " of ",
  nrow(df.atac.unambiguous.orientation),
  " single-direct-promoter/TSS loops."
)


################################################################################
# 6. Annotate Gene Transcript Positions Relative to Loops (inside, span across, or extend beyond the chromatin loop boundaries)
################################################################################

# 6-1. Prepare primary direct gene assignments with tier annotation
df.direct.assignment.for.position <- df.direct.promoter.tss.gene.assignment %>% # 26,920
  mutate(
    assignment_tier = "primary_direct_TSS_plus_minus_1kb",
    .before = 1
  )

# 6-2. Evaluate loop-relative positions for all individual Ensembl transcripts
# case1: gene: transcript = 1:N
# case2: EPD-only genes = 'missing-annotation'
df.direct.transcript.position.detail <- annotate_assignment_transcript_positions(
  df.assignment = df.direct.promoter.tss.gene.assignment,
  df.transcript = df.true.tss.transcript,
  df.loop = df.loop.distinct.2mb,
  assignment.tier = "primary_direct_TSS_plus_minus_1kb",
  promoter.window.flank.bp = promoter.window.flank.bp
)

# df.direct.transcript.position.detail: 57,845 (1:N) <- 26,920

# 6-3. Summarise isoform-level flags at the unique loop-anchor-gene assignment level: isoform can have multiple flags depending on the length of transcript.
df.direct.gene.assignment.position.flags <- summarise_assignment_transcript_positions(
  df.position.detail = df.direct.transcript.position.detail,
  df.assignment = df.direct.assignment.for.position
)

df.direct.gene.assignment.position.flags %>% head(2)

# Confirm that transcript-position summaries preserve all direct gene assignments.
assert_analysis_condition(
  condition = nrow(df.direct.gene.assignment.position.flags) == nrow(df.direct.promoter.tss.gene.assignment), # df.direct.promoter.tss.gene.assignment: 26,920
  message = "Transcript-position summaries did not preserve all direct gene assignments.",
  success.message = sprintf(
    "Verified: Transcript-position summaries preserve all %d direct gene assignments (1:1 match).",
    nrow(df.direct.gene.assignment.position.flags)
  )
)

# A primary Ensembl assignment must recover at least one transcript whose TSS
# +/-1-kb promoter window overlaps the assigned anchor. Exact TSS-point overlap
# remains a descriptive strict-sensitivity flag. EPD-only assignments are not
# subjected to this Ensembl-transcript consistency check.
assert_analysis_condition(
  condition = !any(
    df.direct.gene.assignment.position.flags$has_direct_true_tss &
      !df.direct.gene.assignment.position.flags$any_transcript_promoter_window_overlaps_assigned_anchor
  ),
  message = paste0(
    "A primary Ensembl assignment lacks a matching anchor-overlapping TSS ",
    "+/-1-kb promoter window."
  ),
  success.message = "Verified: All direct Ensembl assignments recover matching anchor-overlapping promoter windows."
)

# 6-4. Quantify transcript containment rates across loop resolutions
df.transcript.position.summary.input <- df.direct.gene.assignment.position.flags

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

# 6-5. Document transcript position operational definitions

df.transcript.position.definition <- tribble(
  ~field_or_rule, ~definition,
  "loop_span", "Inclusive genomic span from the outer start to the outer end of the two anchors.",
  "inter_anchor_interval", "Open interval between the left anchor end and right anchor start; anchor bases are excluded.",
  "transcript_fully_within_loop_span", "The complete Ensembl transcript body lies inside the inclusive loop span.",
  "all_vs_any_transcript_flags", "All annotated isoforms are retained; no canonical or nearest transcript is selected for containment.",
  "analysis_role", "Transcript containment is descriptive evidence only and is not a loop or gene-assignment retention filter."
)

message(
  "Added transcript-position flags to ",
  nrow(df.direct.gene.assignment.position.flags),
  " direct loop-anchor-gene assignments."
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

df.ctcf.evidence <- df.loop.distinct.2mb %>%
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

df.ctcf.evidence %>%
  count(predicted_ctcf_motif_annotation_class) %>%
  mutate(pct = round(100 * n / sum(n), 1))

df.ctcf.evidence %>%
  group_by(resolution, predicted_ctcf_motif_annotation_class) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(pct = round(100 * n / sum(n), 1))

################################################################################
# 8. Layered loop evidence and mutually exclusive categories
#
# Regulatory categories use direct promoter/TSS and TSS-excluded ATAC evidence.
# Predicted CTCF motif fields are joined only as independent annotations.
################################################################################

# 8-1. Assemble master loop evidence table and assign mutually exclusive categories
df.loop.evidence <- df.loop.distinct.2mb %>%
  left_join(
    df.direct.promoter.tss.loop.summary %>%
      dplyr::select(loop_id, n_direct_anchor_sides, n_direct_genes_across_anchors, has_any_direct_promoter_tss, has_direct_promoter_tss_both_anchors, direct_anchor_assignment_class),
    by = "loop_id"
  ) %>%
  left_join(
    df.atac.loop.summary %>%
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
      TRUE ~ "no_direct_promoter_TSS"
    )
  ) %>%
  arrange(chr1, start1, end1, chr2, start2, end2)

# Confirm that the master loop evidence table retains exactly one row per pooled loop.
assert_analysis_condition(
  condition = nrow(df.loop.evidence) == nrow(df.loop.distinct.2mb) &&
    n_distinct(df.loop.evidence$loop_id) == nrow(df.loop.distinct.2mb),
  message = "Loop evidence master table is not exactly one row per pooled loop.",
  success.message = sprintf(
    "Verified: Loop evidence master table covers all %d pooled loops (1:1 match, 0 duplicate).",
    nrow(df.loop.evidence)
  )
)

# 8-2. Document operational definitions of major loop categories
df.loop.category.definition <- tribble(
  ~revised_major_category, ~required_evidence, ~interpretation,
  "putative_regulatory_direct_promoter_TSS_opposite_nonTSS_ATAC",
  "Exactly one anchor directly overlaps a promoter/TSS +/-1-kb window; the opposite anchor has >=50 bp true-TSS-excluded ATAC overlap.",
  "Open-chromatin-supported putative promoter-to-distal-element contact; not a validated enhancer or functional P-E interaction.",
  "promoter_promoter_compatible_both_direct_anchors",
  "Both anchors directly overlap promoter/TSS +/-1-kb windows.",
  "Promoter-promoter-compatible contact retained without forcing either anchor to be an enhancer.",
  "single_direct_promoter_TSS_without_opposite_nonTSS_ATAC",
  "Exactly one anchor directly overlaps a promoter/TSS +/-1-kb window; the opposite anchor lacks >=50 bp true-TSS-excluded ATAC overlap.",
  "Directionally assignable promoter contact without distal open-chromatin support.",
  "no_direct_promoter_TSS",
  "Neither anchor directly overlaps a promoter/TSS +/-1-kb window.",
  "Pooled HiCCUPS call retained without direct promoter/TSS assignment."
)

# 8-3. Summarise major category distribution overall and by resolution
df.loop.category.summary <- df.loop.evidence %>%
  count(revised_major_category, name = "n_loops") %>%
  mutate(pct_pooled_loops = round(100 * n_loops / sum(n_loops), 1)) %>%
  arrange(desc(n_loops), revised_major_category)

df.loop.category.by.resolution <- df.loop.evidence %>%
  count(resolution, revised_major_category, name = "n_loops") %>%
  group_by(resolution) %>%
  mutate(pct_within_resolution = round(100 * n_loops / sum(n_loops), 1)) %>%
  ungroup() %>%
  arrange(resolution, desc(n_loops), revised_major_category)

message("Loop categories cover all ", nrow(df.loop.distinct.2mb), " calls.")

################################################################################
# 9. final loop categories
################################################################################
################################################################################
# 10. figures
################################################################################

# ==============================================================================
# 10-1. Figure 1: Read Category and Correlation with Loop Counts (Juicer sb-option QC)
# ==============================================================================

# Parse live Juicer sb-option QC stats (inter_30.txt) across all 10 HRDP samples
sb.qc.dir <- path.expand(Sys.getenv(
  "JUICER_SB_OPTIONS_DIR",
  unset = "/Users/pete/Library/CloudStorage/GoogleDrive-wellclouder@gmail.com/My Drive/research/juicer-w-sb-options"
))

df.sample.qc.mapping <- tribble(
  ~sample, ~strain,
  "592BB", "SHR/OlaIpcv",
  "607", "HXB10",
  "74AA", "F344/Stm",
  "A2DB", "LE/Stm",
  "D765A", "BXH6",
  "DA08A", "HXB2",
  "DA21A", "SHR/OlaIpcvxBN/NHsdMcwi",
  "DA68A", "HXB31",
  "DBA9A", "HXB23",
  "DE8BA", "BN-Lx"
)

parse_juicer_sb_qc <- function(sample_code, strain_name, root_dir) {
  qc_file <- file.path(root_dir, sample_code, "inter_30.txt")
  lines <- tryCatch(
    readLines(qc_file, warn = FALSE),
    error = function(e) NULL
  )
  if (is.null(lines) || length(lines) == 0L) {
    # Fallback candidates for offline execution
    fallback_candidates <- c(
      file.path(coord.cache.dir, "depth_qc_by_strain.tsv"),
      file.path(output.dir, "depth_qc_by_strain.tsv"),
      file.path(data.dir, "library_complexity_592BB.tsv")
    )
    for (fb_file in fallback_candidates) {
      if (file.exists(fb_file)) {
        fb_df <- tryCatch(read.delim(fb_file), error = function(e) NULL)
        if (!is.null(fb_df)) {
          strain_col <- if ("Strain" %in% names(fb_df)) "Strain" else "strain"
          matched <- fb_df %>% filter(.data[[strain_col]] == strain_name)
          if (nrow(matched) > 0) {
            return(matched %>% mutate(Strain = strain_name, sample = sample_code))
          }
        }
      }
    }
    stop("Missing required QC file: ", qc_file, call. = FALSE)
  }
  
  get_val <- function(key) {
    l <- lines[str_detect(lines, fixed(paste0(key, ":")))]
    if (length(l) == 0) return(NA_real_)
    val_str <- str_extract(l[1], "(?<=:)[0-9, ]+")
    as.numeric(gsub("[ ,]", "", val_str))
  }
  
  tibble(
    Strain = strain_name,
    sample = sample_code,
    Sequenced_RP = get_val("Sequenced Read Pairs"),
    Normal_Paired = get_val("Normal Paired"),
    Chimeric_Paired = get_val("Chimeric Paired"),
    Chimeric_Ambiguous = get_val("Chimeric Ambiguous"),
    Unmapped = get_val("Unmapped"),
    Alignable_Normal_N_Chimeric = get_val("Alignable (Normal+Chimeric Paired)"),
    Unique_Reads = get_val("Unique Reads"),
    PCR_Duplicates = get_val("PCR Duplicates"),
    Optical_Duplicates = get_val("Optical Duplicates"),
    Below_MAPQ_Threshold = get_val("Below MAPQ Threshold"),
    `Hi-C_Contacts` = get_val("Hi-C Contacts"),
    `Inter-chromosomal` = get_val("Inter-chromosomal"),
    `Intra-chromosomal` = get_val("Intra-chromosomal"),
    Short_Range_20Kb = get_val("Short Range (<20Kb)"),
    Long_Range_20Kb = get_val("Long Range (>20Kb)")
  )
}

df.figure1.qc <- map2_dfr(
  df.sample.qc.mapping$sample,
  df.sample.qc.mapping$strain,
  parse_juicer_sb_qc,
  root_dir = sb.qc.dir
) %>%
  mutate(
    Duplicates = PCR_Duplicates + Optical_Duplicates,
    Chimeric_ambiguous_and_Unmapped = Chimeric_Ambiguous + Unmapped,
    Unique_Reads_Percentage = (Unique_Reads / Sequenced_RP) * 100,
    Duplicates_Percentage = (Duplicates / Sequenced_RP) * 100,
    Chimeric_ambiguous_and_Unmapped_Percentage = (Chimeric_ambiguous_and_Unmapped / Sequenced_RP) * 100
  )

# --- Panel 1a: Read Category Stacked Horizontal Bar Plot ---
df.fig1a.melted <- df.figure1.qc %>%
  pivot_longer(
    cols = c("Unique_Reads_Percentage", "Duplicates_Percentage", "Chimeric_ambiguous_and_Unmapped_Percentage"),
    names_to = "Category",
    values_to = "Percentage"
  ) %>%
  mutate(
    Category = factor(
      Category,
      levels = c("Unique_Reads_Percentage", "Duplicates_Percentage", "Chimeric_ambiguous_and_Unmapped_Percentage"),
      labels = c("Unique Reads", "Duplicates", "Chimeric Ambiguous + Unmapped")
    ),
    Strain = factor(Strain, levels = rev(c(
      "SHR/OlaIpcvxBN/NHsdMcwi", "SHR/OlaIpcv", "LE/Stm", "HXB31", "HXB23",
      "HXB2", "HXB10", "F344/Stm", "BXH6", "BN-Lx"
    )))
  )

plot.figure1a.read.category <- ggplot(
  df.fig1a.melted,
  aes(x = Strain, y = Percentage / 100, fill = Category)
) +
  geom_bar(stat = "identity", position = "fill", width = 0.82) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_fill_manual(
    values = c(
      "Unique Reads" = "#4169E1",
      "Duplicates" = "#FFA500",
      "Chimeric Ambiguous + Unmapped" = "#1A1A1A"
    )
  ) +
  labs(
    title = "a. Read Category",
    x = "Sample",
    y = "Percent of Total",
    fill = "Category"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 9),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8))
  )

# --- Panel 1b: Read Counts vs Loop Counts by Category (<2Mb Loops) ---
df.figure1.loop.counts <- df.sample.loop.1based %>%
  filter(passes_lt2mb) %>%
  count(strain, name = "num_loop")

df.fig1b.merged <- df.figure1.loop.counts %>%
  inner_join(df.figure1.qc, by = c("strain" = "Strain")) %>%
  transmute(
    strain,
    num_loop,
    "Total Reads" = Sequenced_RP,
    "Unique Reads" = Unique_Reads,
    "Alignable Reads" = Alignable_Normal_N_Chimeric
  ) %>%
  pivot_longer(
    cols = c("Total Reads", "Unique Reads", "Alignable Reads"),
    names_to = "Sequencing_Metric",
    values_to = "Depth"
  ) %>%
  mutate(
    Sequencing_Metric = factor(
      Sequencing_Metric,
      levels = c("Alignable Reads", "Total Reads", "Unique Reads")
    )
  )

# Calculate Pearson correlation per metric
df.fig1b.correlations <- df.fig1b.merged %>%
  group_by(Sequencing_Metric) %>%
  group_modify(~ cor.test(.x$Depth, .x$num_loop) %>% broom::tidy()) %>%
  ungroup() %>%
  mutate(
    r = format(round(estimate, 2), nsmall = 2),
    p = format(round(p.value, 4), nsmall = 4),
    label = paste0("R = ", r, ", p = ", p)
  )

annot_positions <- tibble(
  Sequencing_Metric = factor(
    c("Alignable Reads", "Total Reads", "Unique Reads"),
    levels = c("Alignable Reads", "Total Reads", "Unique Reads")
  ),
  x = rep(min(df.fig1b.merged$Depth) * 0.98, 3),
  y = c(
    max(df.fig1b.merged$num_loop) * 0.98,
    max(df.fig1b.merged$num_loop) * 0.92,
    max(df.fig1b.merged$num_loop) * 0.86
  )
)

df.fig1b.annot <- left_join(df.fig1b.correlations, annot_positions, by = "Sequencing_Metric")

fig1_metric_colors <- c(
  "Alignable Reads" = "#F8766D",
  "Total Reads" = "#00BA38",
  "Unique Reads" = "#619CFF"
)

plot.figure1b.loop.correlation <- ggplot(
  df.fig1b.merged,
  aes(x = Depth, y = num_loop, color = Sequencing_Metric)
) +
  geom_point(size = 2.2) +
  geom_smooth(
    aes(fill = Sequencing_Metric),
    method = "lm",
    se = TRUE,
    linewidth = 1,
    alpha = 0.25
  ) +
  scale_x_continuous(
    labels = scales::label_number(scale = 1e-6, suffix = "M"),
    breaks = seq(300e6, 900e6, 200e6)
  ) +
  scale_y_continuous(breaks = seq(2000, 10000, 2000)) +
  scale_color_manual(values = fig1_metric_colors) +
  scale_fill_manual(values = fig1_metric_colors) +
  labs(
    title = "b. Read Counts vs Loop Counts by Category",
    x = "Number of Reads",
    y = "Number of Loops",
    color = "Category",
    fill = "Category"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8))
  ) +
  geom_text(
    data = df.fig1b.annot,
    aes(x = x, y = y, label = label, color = Sequencing_Metric),
    hjust = 0,
    size = 4.2,
    fontface = "italic",
    show.legend = FALSE
  )

# Combine Figure 1 panels
plot.figure1.revised <- patchwork::wrap_plots(
  plot.figure1a.read.category,
  plot.figure1b.loop.correlation,
  ncol = 2,
  nrow = 1
)

# Save Figure 1 to results directory
results.dir <- file.path(getwd(), "r_files", "revision", "revision_main", "results")
if (!dir.exists(results.dir)) {
  dir.create(results.dir, recursive = TRUE)
}

figure1.png.path <- file.path(results.dir, "revision_figure1_sequencing_and_loop_depth.png")
figure1.pdf.path <- file.path(results.dir, "revision_figure1_sequencing_and_loop_depth.pdf")

ggsave(
  filename = figure1.png.path,
  plot = plot.figure1.revised,
  width = 11,
  height = 5.2,
  dpi = 300
)

ggsave(
  filename = figure1.pdf.path,
  plot = plot.figure1.revised,
  width = 11,
  height = 5.2,
  device = "pdf"
)

message("Figure 1 sequencing and loop depth plots successfully saved to:")
message("  - PNG: ", figure1.png.path)
message("  - PDF: ", figure1.pdf.path)

# ==============================================================================
# 10-2. Figure 2: Genome-wide distribution of loop counts detected at three resolutions
# ==============================================================================

chr_levels_rn7 <- c(paste0("chr", 1:20), "chrX", "chrY")

df_chr_loop_counts <- df.loop.distinct.2mb %>%
  group_by(chr1, resolution) %>%
  summarise(n_loops = n_distinct(loop_id), .groups = "drop") %>%
  mutate(
    chr1 = factor(chr1, levels = chr_levels_rn7),
    resolution = factor(resolution, levels = c("5K", "10K", "25K"))
  ) %>%
  filter(!is.na(chr1))

plot.figure2.chr.loops <- ggplot(
  df_chr_loop_counts,
  aes(x = chr1, y = n_loops, fill = resolution)
) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Chromosome",
    y = "Number of Loop",
    fill = "resolution"
  ) +
  scale_fill_manual(
    values = c("5K" = "#a6cee3", "10K" = "#1f78b4", "25K" = "#1f3a93")
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

figure2.png.path <- file.path(results.dir, "revision_figure2_loop_counts_per_chr.png")
figure2.pdf.path <- file.path(results.dir, "revision_figure2_loop_counts_per_chr.pdf")

ggsave(
  filename = figure2.png.path,
  plot = plot.figure2.chr.loops,
  width = 8.5,
  height = 5.5,
  dpi = 300
)

ggsave(
  filename = figure2.pdf.path,
  plot = plot.figure2.chr.loops,
  width = 8.5,
  height = 5.5,
  device = "pdf"
)

message("Figure 2 chromosome loop distribution plots successfully saved to:")
message("  - PNG: ", figure2.png.path)
message("  - PDF: ", figure2.pdf.path)

# ==============================================================================
# 10-3. Figure 3: Shared Loops by Resolution & Sample
# ==============================================================================

# Clean strain names and unique sequencing depth annotations
strain_clean_map <- c(
  "SHR/OlaIpcv" = "SHR/OlaIpcv",
  "HXB10" = "HXB10/Ipcv",
  "F344/Stm" = "F344/Stm",
  "LE/Stm" = "LE/Stm",
  "BXH6" = "BXH6/Cub",
  "HXB2" = "HXB2/Ipcv",
  "SHR/OlaIpcvxBN/NHsdMcwi" = "SHRxBN F1",
  "HXB31" = "HXB31/Ipcv",
  "HXB23" = "HXB23/Ipcv",
  "BN-Lx" = "BN-Lx/Cub"
)

strain_depth_map <- c(
  "HXB31/Ipcv" = "491M",
  "HXB10/Ipcv" = "528M",
  "BN-Lx/Cub" = "392M",
  "SHRxBN F1" = "486M",
  "SHR/OlaIpcv" = "422M",
  "HXB2/Ipcv" = "281M",
  "HXB23/Ipcv" = "450M",
  "BXH6/Cub" = "350M",
  "LE/Stm" = "185M",
  "F344/Stm" = "181M"
)

df.sample.2mb.figure3 <- df.sample.loop.1based %>%
  filter(passes_lt2mb) %>%
  mutate(strain_label = recode(strain, !!!strain_clean_map))

# --- Panel 3a: Shared Loops by Resolution (Mean & SD per resolution) ---
shared_loops_per_res <- df.sample.2mb.figure3 %>%
  group_by(resolution, loop_id) %>%
  summarise(n_samples = n_distinct(sample), .groups = "drop") %>%
  filter(n_samples > 1) %>%
  inner_join(df.sample.2mb.figure3, by = c("resolution", "loop_id"))

df.fig3a.summary <- shared_loops_per_res %>%
  group_by(resolution, sample) %>%
  summarise(shared_count = n_distinct(loop_id), .groups = "drop") %>%
  group_by(resolution) %>%
  summarise(
    mean_shared_loops = mean(shared_count),
    sd_shared_loops = sd(shared_count),
    .groups = "drop"
  ) %>%
  mutate(resolution = factor(resolution, levels = c("5K", "10K", "25K")))

plot.figure3a.shared.resolution <- ggplot(
  df.fig3a.summary,
  aes(x = resolution, y = mean_shared_loops, fill = resolution)
) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_errorbar(
    aes(
      ymin = pmax(0, mean_shared_loops - sd_shared_loops),
      ymax = mean_shared_loops + sd_shared_loops
    ),
    width = 0.2,
    color = "black",
    linewidth = 0.5
  ) +
  scale_fill_manual(
    values = c("5K" = "#a6cee3", "10K" = "#1f78b4", "25K" = "#1f3a93")
  ) +
  coord_cartesian(ylim = c(0, 2100)) +
  labs(
    title = "a. Shared Loop by Resolution",
    x = "Resolution",
    y = "Shared Loop"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "plain", hjust = 0, size = 12),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8))
  )

# --- Panel 3b: Shared Loops by Sample (Stacked Unique vs Shared Loops) ---
loop.sharing.overall <- df.sample.2mb.figure3 %>%
  group_by(loop_id) %>%
  summarise(n_samples = n_distinct(sample), .groups = "drop") %>%
  mutate(loop_type = factor(
    ifelse(n_samples > 1, "Shared Loops", "Unique Loops"),
    levels = c("Unique Loops", "Shared Loops")
  ))

df.fig3b.data <- df.sample.2mb.figure3 %>%
  inner_join(loop.sharing.overall, by = "loop_id") %>%
  group_by(strain_label, loop_type) %>%
  summarise(n_loops = n_distinct(loop_id), .groups = "drop")

strain_order_by_total <- df.fig3b.data %>%
  group_by(strain_label) %>%
  summarise(total = sum(n_loops), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(strain_label)

df.fig3b.data <- df.fig3b.data %>%
  mutate(strain_label = factor(strain_label, levels = strain_order_by_total))

# Create dual x-axis labels with sequencing depth
strain_depth_labels <- map_chr(strain_order_by_total, function(st) {
  depth <- strain_depth_map[st]
  if (is.na(depth)) depth <- ""
  paste0(depth, "\n", st)
})
names(strain_depth_labels) <- strain_order_by_total

plot.figure3b.shared.sample <- ggplot(
  df.fig3b.data,
  aes(x = strain_label, y = n_loops, fill = loop_type)
) +
  geom_bar(stat = "identity", position = "stack", width = 0.85) +
  scale_fill_manual(
    values = c("Unique Loops" = "#FFA300", "Shared Loops" = "#00573F"),
    breaks = c("Unique Loops", "Shared Loops")
  ) +
  scale_x_discrete(labels = strain_depth_labels) +
  scale_y_continuous(breaks = seq(0, 10000, 2500), limits = c(0, 9500)) +
  labs(
    title = "b. Shared Loops by Sample",
    x = "Strain",
    y = "Number of Loops",
    fill = "Loop Type"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "plain", hjust = 0, size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "plain"),
    legend.text = element_text(size = 9),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8))
  )

# Combine Figure 3 panels
plot.figure3.revised <- patchwork::wrap_plots(
  plot.figure3a.shared.resolution,
  plot.figure3b.shared.sample,
  widths = c(1, 2.2),
  nrow = 1
)

# Save Figure 3 to results directory
results.dir <- file.path(getwd(), "r_files", "revision", "revision_main", "results")
if (!dir.exists(results.dir)) {
  dir.create(results.dir, recursive = TRUE)
}

figure3.png.path <- file.path(results.dir, "revision_figure3_shared_loops.png")
figure3.pdf.path <- file.path(results.dir, "revision_figure3_shared_loops.pdf")

ggsave(
  filename = figure3.png.path,
  plot = plot.figure3.revised,
  width = 9.5,
  height = 4.2,
  dpi = 300
)

ggsave(
  filename = figure3.pdf.path,
  plot = plot.figure3.revised,
  width = 9.5,
  height = 4.2,
  device = "pdf"
)

message("Figure 3 shared loop plots successfully saved to:")
message("  - PNG: ", figure3.png.path)
message("  - PDF: ", figure3.pdf.path)

# ==============================================================================
# 10-4. Figure 4: Chromosomal distribution of CTCF binding site density and gene density correlation
# ==============================================================================

ncbi_gff_file <- file.path(
  data.dir,
  "GCF_015227675.2_mRatBN7.2_genomic.gff"
)
if (!file.exists(ncbi_gff_file)) {
  # Fallback to local dropbox directory
  ncbi_gff_file <- "/Users/pete/Library/CloudStorage/Dropbox/Gateway_to_Hao/enhancer/data/GCF_015227675.2_mRatBN7.2_genomic.gff"
}

chr_len_file <- file.path(
  enhancer.project.dir,
  "data",
  "rn7_chromosome_length_from_ucsc.tsv"
)

chromosome_levels_fig4 <- as.character(c(1:20, "X", "Y"))
valid_chromosomes_fig4 <- c(paste0("chr", 1:20), "chrX", "chrY")

ncbi_accession_to_chr <- c(
  "NC_051336.1" = "chr1",  "NC_051337.1" = "chr2",  "NC_051338.1" = "chr3",
  "NC_051339.1" = "chr4",  "NC_051340.1" = "chr5",  "NC_051341.1" = "chr6",
  "NC_051342.1" = "chr7",  "NC_051343.1" = "chr8",  "NC_051344.1" = "chr9",
  "NC_051345.1" = "chr10", "NC_051346.1" = "chr11", "NC_051347.1" = "chr12",
  "NC_051348.1" = "chr13", "NC_051349.1" = "chr14", "NC_051350.1" = "chr15",
  "NC_051351.1" = "chr16", "NC_051352.1" = "chr17", "NC_051353.1" = "chr18",
  "NC_051354.1" = "chr19", "NC_051355.1" = "chr20",
  "NC_051356.1" = "chrX",  "NC_051357.1" = "chrY",  "NC_001665.2" = "chrM"
)

# Chromosome length data
df_chrom_fig4 <- read.table(chr_len_file, sep = "\t", col.names = c("chr", "end")) %>%
  mutate(
    start = 0,
    Chr = str_remove(chr, "^chr"),
    End = as.numeric(end)
  ) %>%
  filter(Chr %in% chromosome_levels_fig4) %>%
  arrange(factor(Chr, levels = chromosome_levels_fig4))

# Read NCBI RefSeq genes if available
if (file.exists(ncbi_gff_file)) {
  message("Reading NCBI RefSeq gene catalog for Figure 4...")
  gff_raw_fig4 <- read_tsv(
    ncbi_gff_file,
    comment = "#",
    col_names = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"),
    col_types = cols(.default = "c")
  ) %>%
    filter(feature == "gene") %>%
    filter(chr %in% names(ncbi_accession_to_chr)) %>%
    mutate(
      chr = ncbi_accession_to_chr[chr],
      start = as.numeric(start),
      end = as.numeric(end)
    ) %>%
    filter(chr %in% valid_chromosomes_fig4)

  df_ncbi_genes <- gff_raw_fig4 %>%
    mutate(
      gene_id = str_match(attribute, "ID=gene-([^;]+)")[, 2],
      gene_name = str_match(attribute, "gene=([^;]+)")[, 2],
      gene_biotype = str_match(attribute, "gene_biotype=([^;]+)")[, 2]
    ) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    dplyr::select(chr, start, end, strand, gene_id, gene_name, gene_biotype)
} else {
  # Fallback to Ensembl transcripts
  message("NCBI RefSeq GFF file not found. Using Ensembl genes for Figure 4 correlation.")
  df_ncbi_genes <- df.transcript.ensembl.rn7.1based %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    transmute(
      chr,
      start = transcript_start,
      end = transcript_end,
      strand,
      gene_id,
      gene_name,
      gene_biotype = transcript_biotype
    )
}

# CTCF data
df_ctcf_fig4 <- as_tibble(gr.ctcf.motif) %>%
  transmute(
    chr = as.character(seqnames),
    start = start,
    end = end
  ) %>%
  filter(chr %in% valid_chromosomes_fig4)

# --- Panel 4a: CTCF binding site density ideogram ---
karyotype_data_fig4 <- df_chrom_fig4 %>%
  transmute(Chr = as.character(Chr), Start = 0L, End = as.integer(End)) %>%
  arrange(factor(Chr, levels = chromosome_levels_fig4))

process_feature_bins_ideogram <- function(feature_df, chrom_df, bin_size = 1000000) {
  bin_list <- list()
  for (i in seq_len(nrow(chrom_df))) {
    chr_name <- chrom_df$Chr[i]
    chr_end <- chrom_df$End[i]
    starts <- seq(0L, chr_end, by = bin_size)
    ends <- pmin(starts + bin_size, chr_end)
    bin_list[[i]] <- tibble(Chr = chr_name, Start = as.integer(starts), End = as.integer(ends))
  }
  bins_df <- bind_rows(bin_list)

  feature_gr <- GRanges(seqnames = feature_df$Chr, ranges = IRanges(start = feature_df$Start, end = feature_df$End))
  bins_gr <- GRanges(seqnames = bins_df$Chr, ranges = IRanges(start = bins_df$Start, end = bins_df$End))

  overlaps <- countOverlaps(bins_gr, feature_gr)
  bins_df %>% mutate(Value = as.numeric(overlaps))
}

ctcf_density_for_ideogram <- process_feature_bins_ideogram(
  feature_df = df_ctcf_fig4 %>% transmute(Chr = str_remove(chr, "^chr"), Start = as.integer(start), End = as.integer(end)),
  chrom_df = karyotype_data_fig4,
  bin_size = 1000000
)

figure4_ideogram_png <- file.path(results.dir, "revision_figure4_panel_a_ctcf_ideogram.png")

# Generate SVG and convert to PNG
old_wd <- getwd()
setwd(results.dir)
ideogram(
  karyotype = karyotype_data_fig4,
  overlaid = ctcf_density_for_ideogram
)
if (file.exists("chromosome.svg")) {
  rsvg::rsvg_png("chromosome.svg", file = figure4_ideogram_png, width = 2400)
  image_read(figure4_ideogram_png) %>%
    image_background(color = "white") %>%
    image_trim(fuzz = 5) %>%
    image_border(color = "white", geometry = "80x80") %>%
    image_write(figure4_ideogram_png)
  file.remove("chromosome.svg")
}
setwd(old_wd)

# --- Panels 4b, 4c, 4d: Correlation scatter plots ---
build_figure4_density_summary <- function(chrom_data, df_gene_input, ctcf_data) {
  chrom_data %>%
    mutate(
      chr = str_remove(as.character(chr), "^chr"),
      chromosome_length_mb = end / 1e6
    ) %>%
    group_by(chr) %>%
    summarise(chromosome_length_mb = max(chromosome_length_mb, na.rm = TRUE), .groups = "drop") %>%
    left_join(
      df_gene_input %>%
        mutate(chr = str_remove(as.character(chr), "^chr")) %>%
        count(chr, name = "gene_count"),
      by = "chr"
    ) %>%
    left_join(
      ctcf_data %>%
        mutate(chr = str_remove(as.character(chr), "^chr")) %>%
        count(chr, name = "ctcf_count"),
      by = "chr"
    ) %>%
    mutate(
      gene_count = replace_na(gene_count, 0L),
      ctcf_count = replace_na(ctcf_count, 0L),
      genes_per_mb = gene_count / chromosome_length_mb,
      ctcf_per_mb = ctcf_count / chromosome_length_mb
    ) %>%
    arrange(factor(chr, levels = chromosome_levels_fig4))
}

make_figure4_scatter_plot <- function(gene_df, title_text) {
  sum_df <- build_figure4_density_summary(df_chrom_fig4, gene_df, df_ctcf_fig4)
  cor_res <- cor.test(sum_df$genes_per_mb, sum_df$ctcf_per_mb, method = "pearson")
  r_val <- round(cor_res$estimate, 3)
  p_val <- formatC(cor_res$p.value, format = "e", digits = 2)

  ggplot(sum_df, aes(x = genes_per_mb, y = ctcf_per_mb)) +
    geom_point(size = 2.2, color = "#2F5597") +
    geom_smooth(method = "lm", se = FALSE, color = "#C44E52", linewidth = 0.7) +
    geom_text(aes(label = chr), nudge_y = max(sum_df$ctcf_per_mb) * 0.03, size = 2.5, check_overlap = TRUE) +
    labs(
      title = title_text,
      subtitle = paste0("Pearson r = ", r_val, " (p = ", p_val, ")"),
      x = "Genes per Mb",
      y = "CTCF sites per Mb"
    ) +
    theme_bw(base_size = 9) +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8.5),
      panel.grid.minor = element_blank(),
      plot.margin = margin(7, 6, 4, 6)
    )
}

plot_fig4_b <- make_figure4_scatter_plot(df_ncbi_genes, "b. All genes")
plot_fig4_c <- make_figure4_scatter_plot(df_ncbi_genes %>% filter(gene_biotype == "protein_coding"), "c. Protein-coding genes")
plot_fig4_d <- make_figure4_scatter_plot(df_ncbi_genes %>% filter(gene_biotype == "lncRNA"), "d. lncRNA genes")

# Combine Figure 4 into 1x4 layout
panel_4a <- ggdraw() +
  draw_image(figure4_ideogram_png, scale = 0.9) +
  draw_label("a", x = 0.08, y = 0.985, hjust = 0, vjust = 1, fontface = "bold", size = 13)

fig4_grid <- cowplot::plot_grid(
  panel_4a, plot_fig4_b, plot_fig4_c, plot_fig4_d,
  ncol = 4,
  align = "hv",
  rel_widths = c(1.35, 1, 1, 1)
)

plot.figure4.revised <- ggdraw() +
  draw_plot(fig4_grid, x = -0.05, y = -0.055, width = 1.043, height = 1.11)

figure4.png.path <- file.path(results.dir, "revision_figure4_ctcf_and_gene_density.png")
figure4.pdf.path <- file.path(results.dir, "revision_figure4_ctcf_and_gene_density.pdf")

ggsave(
  filename = figure4.png.path,
  plot = plot.figure4.revised,
  width = 15.4,
  height = 4.3,
  dpi = 300
)

ggsave(
  filename = figure4.pdf.path,
  plot = plot.figure4.revised,
  width = 15.4,
  height = 4.3,
  device = "pdf"
)

message("Figure 4 CTCF ideogram and gene density correlation plots successfully saved to:")
message("  - PNG: ", figure4.png.path)
message("  - PDF: ", figure4.pdf.path)

# ==============================================================================
# 10-5. Figure 5: Density plots for positional distribution of CTCF, TSS, and promoters
# ==============================================================================
figure5.resolution.colors <- c(
  "5K" = "#a6cee3",
  "10K" = "#1f78b4",
  "25K" = "#1f3a93"
)

df.figure5.loop.window <- df.loop.distinct.2mb %>%
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

# 10-5-1. Raw relative positions for predicted CTCF motif intervals
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
      return(tibble())
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

    tibble(
      resolution = df.figure5.loop.window$resolution[loop.index.hit[keep]],
      relative_position = relative.position[keep]
    )
  }
)

df.figure5.ctcf.relative.position <- bind_rows(figure5.ctcf.chromosome.results)

# 10-5-2. Raw relative positions for strand-aware Ensembl TSSs
df.figure5.true.tss.site <- df.true.tss.transcript %>%
  group_by(chr, true_tss_start, true_tss_end, strand) %>%
  summarise(n_transcripts = n_distinct(transcript_id), n_genes = n_distinct(gene_id), .groups = "drop") %>%
  mutate(feature_id = str_c(chr, true_tss_start, strand, sep = ":"), feature_position = as.numeric(true_tss_start))

gr.figure5.true.tss.site <- GRanges(
  seqnames = df.figure5.true.tss.site$chr,
  ranges = IRanges(start = df.figure5.true.tss.site$true_tss_start, end = df.figure5.true.tss.site$true_tss_end)
)

figure5.true.tss.hit <- findOverlaps(gr.figure5.true.tss.site, gr.figure5.loop.window, type = "any", select = "all")
df.figure5.true.tss.relative.position <- tibble(
  loop_index = subjectHits(figure5.true.tss.hit),
  feature_index = queryHits(figure5.true.tss.hit)
) %>%
  transmute(
    loop_id = df.figure5.loop.window$loop_id[loop_index],
    resolution = df.figure5.loop.window$resolution[loop_index],
    feature_position = df.figure5.true.tss.site$feature_position[feature_index],
    anchor1_midpoint = df.figure5.loop.window$anchor1_midpoint[loop_index],
    anchor2_midpoint = df.figure5.loop.window$anchor2_midpoint[loop_index],
    anchor_midpoint_distance = df.figure5.loop.window$anchor_midpoint_distance[loop_index],
    relative_position = (feature_position - anchor1_midpoint) / anchor_midpoint_distance
  ) %>%
  filter(dplyr::between(relative_position, -1, 2))

# 10-5-3. Raw relative positions for EPD promoters
df.figure5.promoter.site <- df.promoter.epd.rn7.1based %>%
  distinct(chr, promoter_start, promoter_end, strand, .keep_all = TRUE) %>%
  transmute(
    feature_id = promoter_annotation_id, chr, promoter_start, promoter_end, strand, gene_id, gene_name,
    feature_position = (promoter_start + promoter_end) / 2
  )

gr.figure5.promoter.site <- GRanges(
  seqnames = df.figure5.promoter.site$chr,
  ranges = IRanges(start = df.figure5.promoter.site$promoter_start, end = df.figure5.promoter.site$promoter_end)
)

figure5.promoter.hit <- findOverlaps(gr.figure5.promoter.site, gr.figure5.loop.window, type = "any", select = "all")
df.figure5.promoter.relative.position <- tibble(
  loop_index = subjectHits(figure5.promoter.hit),
  feature_index = queryHits(figure5.promoter.hit)
) %>%
  transmute(
    loop_id = df.figure5.loop.window$loop_id[loop_index],
    resolution = df.figure5.loop.window$resolution[loop_index],
    feature_position = df.figure5.promoter.site$feature_position[feature_index],
    anchor1_midpoint = df.figure5.loop.window$anchor1_midpoint[loop_index],
    anchor2_midpoint = df.figure5.loop.window$anchor2_midpoint[loop_index],
    anchor_midpoint_distance = df.figure5.loop.window$anchor_midpoint_distance[loop_index],
    relative_position = (feature_position - anchor1_midpoint) / anchor_midpoint_distance
  ) %>%
  filter(dplyr::between(relative_position, -1, 2))

# 10-5-4. Build publication-ready unbinned density panels with sharp peaks
create.figure5.revised.density.plot <- function(
  df.relative.position,
  panel.tag,
  panel.title
) {
  ggplot(
    df.relative.position,
    aes(
      x = relative_position,
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
    coord_cartesian(xlim = c(-1, 2), ylim = c(0, 0.75)) +
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

# 10-5-5. Save Figure 5 density plots
figure5.png.path <- file.path(results.dir, "revision_figure5_density_plot.png")
figure5.pdf.path <- file.path(results.dir, "revision_figure5_density_plot.pdf")

ggsave(
  filename = figure5.png.path,
  plot = plot.figure5abc.revised.density,
  width = 10.5,
  height = 3.3,
  dpi = 300
)

ggsave(
  filename = figure5.pdf.path,
  plot = plot.figure5abc.revised.density,
  width = 10.5,
  height = 3.3,
  device = "pdf"
)

message("Figure 5 density plots successfully saved to:")
message("  - PNG: ", figure5.png.path)
message("  - PDF: ", figure5.pdf.path)
