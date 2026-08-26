# lintr: disable

# Resolve this script and the shared enhancer/r_files directory.
current_script_path <- function() {
  file.args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file.args) > 0L) {
    script.arg <- gsub("~+~", " ", sub("^--file=", "", file.args[[1]]), fixed = TRUE)
    return(normalizePath(script.arg, winslash = "/", mustWork = FALSE))
  }

  frame.files <- vapply(
    sys.frames(),
    function(frame) if (is.null(frame$ofile)) NA_character_ else frame$ofile,
    character(1)
  )
  frame.files <- frame.files[!is.na(frame.files) & nzchar(frame.files)]
  if (!length(frame.files)) return(NA_character_)
  normalizePath(tail(frame.files, 1L), winslash = "/", mustWork = FALSE)
}

resolve_enhancer_r_files_dir <- function() {
  script.path <- current_script_path()
  start.dirs <- c(
    if (is.na(script.path)) NA_character_ else dirname(script.path),
    getwd()
  )
  ancestor.dirs <- unique(unlist(lapply(start.dirs, function(path) {
    if (is.na(path) || !nzchar(path)) return(character())
    path <- normalizePath(path, winslash = "/", mustWork = FALSE)
    ancestors <- path
    for (i in seq_len(6L)) ancestors <- c(ancestors, dirname(tail(ancestors, 1L)))
    ancestors
  })))

  candidates <- unique(c(
    ancestor.dirs,
    Sys.getenv("ENHANCER_R_FILES_DIR", unset = ""),
    path.expand("~/dropbox/Gateway_to_Hao/enhancer/r_files"),
    path.expand("~/Dropbox/Gateway_to_Hao/enhancer/r_files"),
    Sys.glob(path.expand(
      "~/Library/CloudStorage/Dropbox*/K P/Gateway_to_Hao/enhancer/r_files"
    )),
    Sys.glob(path.expand(
      "~/Library/CloudStorage/Dropbox*/Gateway_to_Hao/enhancer/r_files"
    ))
  ))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  candidates <- candidates[dir.exists(candidates) & file.exists(file.path(candidates, "funcs.R"))]
  if (!length(candidates)) stop("Cannot locate enhancer/r_files with funcs.R.", call. = FALSE)
  normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE)
}

r.files.dir <- resolve_enhancer_r_files_dir()
enhancer.project.dir <- Sys.getenv(
  "ENHANCER_PROJECT_DIR",
  unset = dirname(r.files.dir)
)
enhancer.project.dir <- normalizePath(
  path.expand(enhancer.project.dir),
  winslash = "/",
  mustWork = TRUE
)

message("Using shared enhancer project: ", enhancer.project.dir)

library("tidyverse")
library("GenomicRanges")
library("GenomeInfoDb")

options(tibble.width = Inf)
options(tibble.print_max = Inf)
options(tibble.max_extra_cols = Inf)
options(scipen = 999)

################################################################################
# Resubmission ATAC-seq matched-null validation
#
# This analysis replaces the legacy 15,085-loop/WHERE workflow. It uses the
# coordinate-normalized pooled loop resource, strand-aware true TSS positions,
# direct promoter/TSS assignments, and Duttke et al. 2022 rat PFC snATAC peaks
# prepared by promoter_enhancer_interaction_resubmit.R.
#
# Primary population
# - Every pooled loop with direct promoter/TSS support at exactly one anchor.
# - The opposite anchor is evaluated as the unique candidate distal anchor.
# - Eligibility is defined before ATAC status, avoiding selection on outcome.
#
# Null models
# 1. rigid_loop_pair_relocation: relocate the complete loop geometry within the
#    same chromosome, preserving resolution, both anchor widths, and distance.
# 2. matched_HiC_anchor: sample promoter-free anchors from pooled Hi-C loops,
#    matching chromosome, resolution, anchor side, and loop-distance bin.
#
# ATAC supports open chromatin only. It does not establish enhancer activity or
# validate a functional promoter-enhancer interaction.
################################################################################

########################
# 0. Directories, inputs, and reproducibility settings
########################

revision.dir <- file.path(r.files.dir, "revision")
script.path <- current_script_path()
script.dir <- if (is.na(script.path)) getwd() else dirname(script.path)
revision.main.dir <- file.path(revision.dir, "revision_main")
cache.dir <- file.path(revision.main.dir, "cache_data")
resubmit.input.dir <- Sys.getenv(
  "RESUBMIT_INPUT_DIR",
  unset = file.path(revision.main.dir, "results")
)
output.dir <- Sys.getenv(
  "ATAC_VALIDATION_OUTPUT_DIR",
  unset = Sys.getenv(
    "RESUBMIT_OUTPUT_DIR",
    unset = file.path(script.dir, "results")
  )
)
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(r.files.dir, "funcs.R"))

loop.resource.file <- file.path(
  resubmit.input.dir,
  "revised_pooled_loop_annotation_resource.tsv"
)
direct.assignment.file <- file.path(
  resubmit.input.dir,
  "revised_direct_loop_gene_assignments.tsv"
)
ensembl.transcript.file <- file.path(
  cache.dir,
  "df.transcript.ensembl.rn7.1based.rds"
)
epd.promoter.file <- file.path(
  cache.dir,
  "df.promoter.epd.rn7.1based.rds"
)
atac.cache.file <- file.path(cache.dir, "gr.atac.rn7.1based.rds")
chrom.sizes.file <- Sys.getenv(
  "RN7_CHROM_SIZES_FILE",
  unset = c(
    file.path(revision.main.dir, "inputs", "data", "tracks", "rn7.chrom.sizes"),
    file.path(enhancer.project.dir, "data", "tracks", "rn7.chrom.sizes")
  )[file.exists(c(
    file.path(revision.main.dir, "inputs", "data", "tracks", "rn7.chrom.sizes"),
    file.path(enhancer.project.dir, "data", "tracks", "rn7.chrom.sizes")
  ))][1]
)
chrom.sizes.file <- path.expand(chrom.sizes.file)

required.files <- c(
  loop.resource.file,
  direct.assignment.file,
  ensembl.transcript.file,
  epd.promoter.file,
  atac.cache.file,
  chrom.sizes.file
)
missing.files <- required.files[!file.exists(required.files)]
if (length(missing.files) > 0L) {
  stop(
    "Missing required ATAC matched-null input file(s):\n",
    paste(missing.files, collapse = "\n"),
    call. = FALSE
  )
}

n.permutations <- as.integer(
  Sys.getenv("ATAC_NULL_PERMUTATIONS", unset = "1000")
)
random.seed <- as.integer(Sys.getenv("ATAC_NULL_SEED", unset = "20260727"))
detected.cores <- parallel::detectCores(logical = TRUE)
if (is.na(detected.cores)) detected.cores <- 1L
default.cores <- if (.Platform$OS.type == "windows") {
  1L
} else {
  max(1L, min(4L, detected.cores - 1L))
}
n.cores <- as.integer(Sys.getenv(
  "ATAC_NULL_CORES",
  unset = as.character(default.cores)
))
stopifnot(n.permutations > 0L, n.cores > 0L)
if (.Platform$OS.type == "windows" && n.cores != 1L) {
  warning("Windows uses sequential ATAC permutations; setting cores to 1.")
  n.cores <- 1L
}

atac.minimum.overlap.bp <- 50L
tss.exclusion.flank.bp <- 1000L
atac.fraction.thresholds <- c(0.005, 0.01)
distance.breaks <- c(-Inf, 50000, 100000, 250000, 500000, 1000000, Inf)
distance.labels <- c(
  "le50kb", "gt50_le100kb", "gt100_le250kb",
  "gt250_le500kb", "gt500_le1000kb", "gt1000kb"
)

########################
# 0-1. Local helpers
########################

# Return total and largest single-interval overlap for each query interval.
interval_overlap_stats <- function(gr.query, gr.feature) {
  n.query <- length(gr.query)
  result <- tibble(total_bp = integer(n.query), max_bp = integer(n.query))
  if (!n.query || !length(gr.feature)) return(result)

  hits <- findOverlaps(gr.query, gr.feature, ignore.strand = TRUE)
  if (!length(hits)) return(result)
  overlap.bp <- width(pintersect(
    gr.query[queryHits(hits)],
    gr.feature[subjectHits(hits)],
    ignore.strand = TRUE
  ))
  query.index <- queryHits(hits)
  total.by.query <- tapply(overlap.bp, query.index, sum)
  max.by.query <- tapply(overlap.bp, query.index, max)
  result$total_bp[as.integer(names(total.by.query))] <- as.integer(total.by.query)
  result$max_bp[as.integer(names(max.by.query))] <- as.integer(max.by.query)
  result
}

# Measure TSS-excluded ATAC overlap and width-adjusted support per anchor.
measure_non_tss_atac <- function(df.anchor) {
  gr.anchor <- GRanges(
    seqnames = df.anchor$anchor_chr,
    ranges = IRanges(df.anchor$anchor_start, df.anchor$anchor_end)
  )
  exclusion.bp <- interval_overlap_stats(gr.anchor, gr.known.tss.exclusion)$total_bp
  atac.overlap <- interval_overlap_stats(gr.anchor, gr.atac.non.tss)
  atac.overlap.bp <- atac.overlap$total_bp
  atac.max.overlap.bp <- atac.overlap$max_bp
  anchor.width.bp <- width(gr.anchor)
  non.tss.anchor.bp <- pmax(0L, anchor.width.bp - exclusion.bp)
  atac.fraction <- if_else(
    non.tss.anchor.bp > 0L,
    atac.overlap.bp / non.tss.anchor.bp,
    NA_real_
  )

  bind_cols(
    df.anchor,
    tibble(
      anchor_width_bp = anchor.width.bp,
      known_tss_exclusion_overlap_bp = exclusion.bp,
      non_tss_anchor_bp = non.tss.anchor.bp,
      non_tss_atac_overlap_bp = atac.overlap.bp,
      non_tss_atac_max_single_overlap_bp = atac.max.overlap.bp,
      atac_overlap_any = atac.overlap.bp >= 1L,
      atac_overlap_ge50 =
        atac.max.overlap.bp >= atac.minimum.overlap.bp,
      non_tss_atac_fraction = atac.fraction,
      atac_fraction_ge_0_005 = coalesce(atac.fraction >= 0.005, FALSE),
      atac_fraction_ge_0_01 = coalesce(atac.fraction >= 0.01, FALSE)
    )
  )
}

# Assign fixed loop-distance bins for matched Hi-C anchor sampling.
add_distance_match_bin <- function(df) {
  df %>%
    mutate(
      distance_match_bin = cut(
        loop_distance,
        breaks = distance.breaks,
        labels = distance.labels,
        right = TRUE,
        ordered_result = TRUE
      ),
      exact_match_key = str_c(
        anchor_chr,
        resolution,
        anchor_side,
        distance_match_bin,
        sep = "|"
      ),
      fallback_match_key = str_c(
        anchor_chr,
        resolution,
        distance_match_bin,
        sep = "|"
      ),
      relaxed_match_key = str_c(
        anchor_chr,
        resolution,
        sep = "|"
      )
    )
}

# Select one matched promoter-free Hi-C anchor for each observed candidate.
sample_matched_hic_controls <- function(
  df.observed,
  df.control.pool,
  exact.pool.index,
  fallback.pool.index,
  relaxed.pool.index
) {
  selected.index <- integer(nrow(df.observed))
  match.level <- rep(NA_character_, nrow(df.observed))

  observed.exact.groups <- split(
    seq_len(nrow(df.observed)),
    df.observed$exact_match_key
  )
  for (key.i in names(observed.exact.groups)) {
    observed.index <- observed.exact.groups[[key.i]]
    pool.index <- exact.pool.index[[key.i]]
    if (!is.null(pool.index) && length(pool.index) > 0L) {
      selected.index[observed.index] <- sample(
        pool.index,
        length(observed.index),
        replace = TRUE
      )
      match.level[observed.index] <- "exact_chr_resolution_side_distance_bin"
    }
  }

  fallback.observed.index <- which(selected.index == 0L)
  if (length(fallback.observed.index) > 0L) {
    fallback.groups <- split(
      fallback.observed.index,
      df.observed$fallback_match_key[fallback.observed.index]
    )
    for (key.i in names(fallback.groups)) {
      observed.index <- fallback.groups[[key.i]]
      pool.index <- fallback.pool.index[[key.i]]
      if (is.null(pool.index) || length(pool.index) == 0L) {
        next
      }
      selected.index[observed.index] <- sample(
        pool.index,
        length(observed.index),
        replace = TRUE
      )
      match.level[observed.index] <- "side_relaxed_same_distance_bin"
    }
  }

  relaxed.observed.index <- which(selected.index == 0L)
  if (length(relaxed.observed.index) > 0L) {
    for (observed.index in relaxed.observed.index) {
      key.i <- df.observed$relaxed_match_key[[observed.index]]
      pool.index <- relaxed.pool.index[[key.i]]
      if (is.null(pool.index) || length(pool.index) == 0L) {
        stop(
          "No chromosome/resolution-matched Hi-C control for stratum: ",
          key.i,
          call. = FALSE
        )
      }
      distance.difference <- abs(
        log1p(df.control.pool$loop_distance[pool.index]) -
          log1p(df.observed$loop_distance[[observed.index]])
      )
      nearest.pool.index <- pool.index[
        distance.difference == min(distance.difference)
      ]
      selected.index[[observed.index]] <- sample(nearest.pool.index, 1L)
      match.level[[observed.index]] <-
        "distance_bin_relaxed_nearest_same_chr_resolution"
    }
  }

  df.control.pool[selected.index, ] %>%
    mutate(
      observed_loop_id = df.observed$loop_id,
      match_level = match.level,
      exact_chr_resolution_side_distance_bin_match =
        match.level == "exact_chr_resolution_side_distance_bin"
    )
}

# Relocate each complete loop geometry to a valid position on the same chromosome.
relocate_loop_pairs <- function(df.observed) {
  max.new.span.start <- df.observed$chromosome_size -
    df.observed$loop_span_width_bp + 1L
  if (any(max.new.span.start < 1L)) {
    stop("A loop span exceeds its chromosome boundary.", call. = FALSE)
  }

  new.span.start <- floor(
    stats::runif(
      nrow(df.observed),
      min = 1,
      max = max.new.span.start + 1
    )
  )
  relocation.offset <- new.span.start - df.observed$loop_span_start

  df.observed %>%
    transmute(
      loop_id,
      resolution,
      anchor_side,
      anchor_chr,
      anchor_start = anchor_start + relocation.offset,
      anchor_end = anchor_end + relocation.offset,
      loop_distance,
      relocation_offset_bp = relocation.offset
    )
}

# Relocate both anchors of loops lacking direct promoter/TSS overlap by one shared offset so the
# original chromosome, resolution, anchor widths, and loop distance are kept.
relocate_no_direct_tss_loop_pairs <- function(df.loops) {
  max.new.span.start <- df.loops$chromosome_size -
    df.loops$loop_span_width_bp + 1L
  if (any(max.new.span.start < 1L)) {
    stop("A no-direct-TSS loop span exceeds its chromosome.", call. = FALSE)
  }

  new.span.start <- floor(stats::runif(
    nrow(df.loops),
    min = 1,
    max = max.new.span.start + 1
  ))
  relocation.offset <- new.span.start - df.loops$loop_span_start

  bind_rows(
    df.loops %>%
      transmute(
        loop_id, resolution, anchor_side = "anchor1", anchor_chr = chr1,
        anchor_start = start1 + relocation.offset,
        anchor_end = end1 + relocation.offset
      ),
    df.loops %>%
      transmute(
        loop_id, resolution, anchor_side = "anchor2", anchor_chr = chr2,
        anchor_start = start2 + relocation.offset,
        anchor_end = end2 + relocation.offset
      )
  )
}

# Summarise ATAC support for no-direct-TSS loops at the anchor level and
# at the loop level (at least one positive anchor or both positive anchors).
summarise_no_direct_tss_loop_atac <- function(
  df.anchor.metric,
  permutation = NA_integer_
) {
  metric.columns <- c("atac_overlap_any", "atac_overlap_ge50")
  resolution.labels <- c(sort(unique(df.anchor.metric$resolution)), "ALL")

  map_dfr(resolution.labels, function(resolution.i) {
    df.i <- if (resolution.i == "ALL") {
      df.anchor.metric
    } else {
      filter(df.anchor.metric, resolution == resolution.i)
    }

    map_dfr(metric.columns, function(metric.i) {
      df.loop.i <- df.i %>%
        group_by(loop_id) %>%
        summarise(
          n_positive_anchors = sum(.data[[metric.i]]),
          n_anchors = n(),
          .groups = "drop"
        )
      tibble(
        permutation,
        resolution = resolution.i,
        metric = metric.i,
        n_loops = nrow(df.loop.i),
        n_anchors = nrow(df.i),
        anchor_positive_rate = mean(df.i[[metric.i]]),
        any_anchor_positive_rate = mean(df.loop.i$n_positive_anchors >= 1L),
        both_anchors_positive_rate = mean(df.loop.i$n_positive_anchors == 2L)
      )
    })
  })
}

# Summarise one observed-vs-control permutation using paired binary metrics.
summarise_one_permutation <- function(
  df.observed,
  df.control,
  permutation,
  null.method
) {
  metric.columns <- c(
    "atac_overlap_any",
    "atac_overlap_ge50",
    "atac_fraction_ge_0_005",
    "atac_fraction_ge_0_01"
  )
  resolution.labels <- c(sort(unique(df.observed$resolution)), "ALL")

  map_dfr(resolution.labels, function(resolution.i) {
    row.index <- if (resolution.i == "ALL") {
      seq_len(nrow(df.observed))
    } else {
      which(df.observed$resolution == resolution.i)
    }

    map_dfr(metric.columns, function(metric.i) {
      observed.positive <- df.observed[[metric.i]][row.index]
      control.positive <- df.control[[metric.i]][row.index]
      discordant.observed.only <- sum(observed.positive & !control.positive)
      discordant.control.only <- sum(!observed.positive & control.positive)

      tibble(
        permutation,
        null_method = null.method,
        resolution = resolution.i,
        metric = metric.i,
        n_pairs = length(row.index),
        observed_rate = mean(observed.positive),
        null_rate = mean(control.positive),
        paired_odds_ratio = (
          discordant.observed.only + 0.5
        ) / (
          discordant.control.only + 0.5
        ),
        n_observed_only = discordant.observed.only,
        n_control_only = discordant.control.only
      )
    })
  })
}

################################################################################
# 1. Load the revised pooled loop universe and direct assignments
################################################################################

df.loop.resource <- read_tsv(loop.resource.file, show_col_types = FALSE)
df.direct.assignment <- read_tsv(
  direct.assignment.file,
  show_col_types = FALSE
)
df.chrom.sizes <- read_tsv(
  chrom.sizes.file,
  col_names = c("chr", "chromosome_size"),
  show_col_types = FALSE
) %>%
  distinct(chr, .keep_all = TRUE)

# Recover one direct and one opposite candidate side for every single-promoter loop.
df.single.direct.orientation <- df.direct.assignment %>%
  distinct(loop_id, resolution, anchor_side, opposite_anchor_side) %>%
  inner_join(
    df.loop.resource %>%
      filter(n_direct_anchor_sides == 1L) %>%
      dplyr::select(
        loop_id,
        resolution,
        chr1, start1, end1,
        chr2, start2, end2,
        loop_distance
      ),
    by = c("loop_id", "resolution")
  ) %>%
  mutate(
    anchor_side = opposite_anchor_side,
    anchor_chr = if_else(anchor_side == "anchor1", chr1, chr2),
    anchor_start = if_else(anchor_side == "anchor1", start1, start2),
    anchor_end = if_else(anchor_side == "anchor1", end1, end2),
    loop_span_start = pmin(start1, start2),
    loop_span_end = pmax(end1, end2),
    loop_span_width_bp = loop_span_end - loop_span_start + 1L
  ) %>%
  left_join(df.chrom.sizes, by = c("anchor_chr" = "chr")) %>%
  arrange(loop_id)
stopifnot(
  nrow(df.single.direct.orientation) == sum(df.loop.resource$n_direct_anchor_sides == 1L),
  !anyDuplicated(df.single.direct.orientation$loop_id),
  all(df.single.direct.orientation$chr1 == df.single.direct.orientation$chr2),
  !anyNA(df.single.direct.orientation$chromosome_size)
)

################################################################################
# 2. Rebuild true-TSS-excluded ATAC from the normalized resubmission cache
################################################################################

df.ensembl.transcript <- readRDS(ensembl.transcript.file)
df.epd.promoter <- readRDS(epd.promoter.file)
gr.atac.rn7.1based <- readRDS(atac.cache.file)

# Build one-base strand-aware Ensembl TSS and lifted EPD TSS positions.
df.known.tss <- bind_rows(
  df.ensembl.transcript %>%
    transmute(
      chr,
      tss = if_else(strand == "+", transcript_start, transcript_end)
    ),
  df.epd.promoter %>%
    transmute(chr, tss = epd_tss_start)
) %>%
  filter(!is.na(chr), !is.na(tss)) %>%
  distinct(chr, tss)

gr.known.tss.exclusion <- GRanges(
  seqnames = df.known.tss$chr,
  ranges = IRanges(
    start = pmax(1L, df.known.tss$tss - tss.exclusion.flank.bp),
    end = df.known.tss$tss + tss.exclusion.flank.bp
  )
) %>%
  reduce(ignore.strand = TRUE)

gr.atac.union <- reduce(gr.atac.rn7.1based, ignore.strand = TRUE)
common.seqlevels <- Reduce(
  intersect,
  list(
    seqlevels(gr.atac.union),
    seqlevels(gr.known.tss.exclusion),
    df.chrom.sizes$chr
  )
)
gr.atac.union <- keepSeqlevels(
  gr.atac.union,
  common.seqlevels,
  pruning.mode = "coarse"
)
gr.known.tss.exclusion <- keepSeqlevels(
  gr.known.tss.exclusion,
  common.seqlevels,
  pruning.mode = "coarse"
)
gr.atac.non.tss <- setdiff(
  gr.atac.union,
  gr.known.tss.exclusion,
  ignore.strand = TRUE
)

################################################################################
# 3. Measure observed candidate-anchor ATAC support and threshold sensitivity
################################################################################

df.observed.candidate <- df.single.direct.orientation %>%
  transmute(
    loop_id,
    resolution,
    anchor_side,
    anchor_chr,
    anchor_start,
    anchor_end,
    loop_distance,
    chr1, start1, end1,
    chr2, start2, end2,
    loop_span_start,
    loop_span_end,
    loop_span_width_bp,
    chromosome_size
  ) %>%
  add_distance_match_bin()

df.observed.metric <- measure_non_tss_atac(df.observed.candidate)

# Report absolute and anchor-width-adjusted criteria within every resolution.
df.atac.threshold.sensitivity.by.resolution <- bind_rows(
  df.observed.metric,
  df.observed.metric %>% mutate(resolution = "ALL")
) %>%
  group_by(resolution) %>%
  summarise(
    n_candidate_anchors = n(),
    n_atac_any = sum(atac_overlap_any),
    pct_atac_any = round(100 * mean(atac_overlap_any), 2),
    n_atac_ge50 = sum(atac_overlap_ge50),
    pct_atac_ge50 = round(100 * mean(atac_overlap_ge50), 2),
    median_non_tss_atac_overlap_bp = median(non_tss_atac_overlap_bp),
    median_non_tss_atac_fraction = median(
      non_tss_atac_fraction,
      na.rm = TRUE
    ),
    n_atac_fraction_ge_0_005 = sum(atac_fraction_ge_0_005),
    pct_atac_fraction_ge_0_005 = round(
      100 * mean(atac_fraction_ge_0_005),
      2
    ),
    n_atac_fraction_ge_0_01 = sum(atac_fraction_ge_0_01),
    pct_atac_fraction_ge_0_01 = round(
      100 * mean(atac_fraction_ge_0_01),
      2
    ),
    .groups = "drop"
  ) %>%
  arrange(factor(resolution, levels = c("5K", "10K", "25K", "ALL")))

################################################################################
# 4. Build chromosome/resolution/distance-matched Hi-C anchor controls
################################################################################

# Both anchors from loops lacking direct promoter/TSS overlap at either anchor
# form the control pool. These loops may still contain TSSs elsewhere and are
# therefore named no-direct-TSS loops, not TSS-free genomic intervals.
df.promoter.free.control.pool <- bind_rows(
  df.loop.resource %>%
    filter(n_direct_anchor_sides == 0L) %>%
    transmute(
      control_loop_id = loop_id,
      resolution,
      anchor_side = "anchor1",
      anchor_chr = chr1,
      anchor_start = start1,
      anchor_end = end1,
      loop_distance
    ),
  df.loop.resource %>%
    filter(n_direct_anchor_sides == 0L) %>%
    transmute(
      control_loop_id = loop_id,
      resolution,
      anchor_side = "anchor2",
      anchor_chr = chr2,
      anchor_start = start2,
      anchor_end = end2,
      loop_distance
    )
) %>%
  mutate(loop_id = control_loop_id) %>%
  add_distance_match_bin() %>%
  measure_non_tss_atac() %>%
  arrange(control_loop_id, anchor_side)

# Keep the complete no-direct-TSS loop geometry for the direct comparison
# requested in revision: real no-direct-TSS loops versus 1,000 random loci.
df.no.direct.tss.loops <- df.loop.resource %>%
  filter(n_direct_anchor_sides == 0L) %>%
  dplyr::select(
    loop_id, resolution,
    chr1, start1, end1,
    chr2, start2, end2
  ) %>%
  mutate(
    loop_span_start = pmin(start1, start2),
    loop_span_end = pmax(end1, end2),
    loop_span_width_bp = loop_span_end - loop_span_start + 1L
  ) %>%
  left_join(df.chrom.sizes, by = c("chr1" = "chr"))
stopifnot(
  !anyDuplicated(df.no.direct.tss.loops$loop_id),
  all(df.no.direct.tss.loops$chr1 == df.no.direct.tss.loops$chr2),
  !anyNA(df.no.direct.tss.loops$chromosome_size)
)

df.no.direct.tss.atac.actual <- summarise_no_direct_tss_loop_atac(
  df.promoter.free.control.pool
)
df.no.direct.tss.loop.status <- df.promoter.free.control.pool %>%
  group_by(loop_id, resolution) %>%
  summarise(
    n_anchors_atac_any = sum(atac_overlap_any),
    n_anchors_atac_ge50 = sum(atac_overlap_ge50),
    any_anchor_atac_ge50 = n_anchors_atac_ge50 >= 1L,
    both_anchors_atac_ge50 = n_anchors_atac_ge50 == 2L,
    .groups = "drop"
  )

exact.pool.index <- split(
  seq_len(nrow(df.promoter.free.control.pool)),
  df.promoter.free.control.pool$exact_match_key
)
fallback.pool.index <- split(
  seq_len(nrow(df.promoter.free.control.pool)),
  df.promoter.free.control.pool$fallback_match_key
)
relaxed.pool.index <- split(
  seq_len(nrow(df.promoter.free.control.pool)),
  df.promoter.free.control.pool$relaxed_match_key
)

df.match.availability <- df.observed.metric %>%
  count(
    exact_match_key,
    fallback_match_key,
    relaxed_match_key,
    anchor_chr,
    resolution,
    anchor_side,
    distance_match_bin,
    name = "n_observed_anchors"
  ) %>%
  mutate(
    n_exact_control_candidates = map_int(
      exact_match_key,
      ~ length(exact.pool.index[[.x]])
    ),
    n_fallback_control_candidates = map_int(
      fallback_match_key,
      ~ length(fallback.pool.index[[.x]])
    ),
    n_relaxed_control_candidates = map_int(
      relaxed_match_key,
      ~ length(relaxed.pool.index[[.x]])
    ),
    exact_match_available = n_exact_control_candidates > 0L,
    fallback_match_available = n_fallback_control_candidates > 0L,
    relaxed_match_available = n_relaxed_control_candidates > 0L
  ) %>%
  arrange(anchor_chr, resolution, anchor_side, distance_match_bin)

assert_analysis_condition(
  all(df.match.availability$relaxed_match_available),
  paste0(
    "At least one observed stratum lacks a chromosome/resolution-matched ",
    "promoter-free Hi-C control."
  )
)

################################################################################
# 5. Run paired matched-null permutations
################################################################################

# Each permutation uses an independent deterministic seed. On macOS/Linux,
# mclapply parallelizes permutations without changing their returned order.
run_one_atac_null_permutation <- function(permutation.i) {
  set.seed(random.seed + permutation.i)

  df.matched.control <- sample_matched_hic_controls(
    df.observed = df.observed.metric,
    df.control.pool = df.promoter.free.control.pool,
    exact.pool.index = exact.pool.index,
    fallback.pool.index = fallback.pool.index,
    relaxed.pool.index = relaxed.pool.index
  )
  df.relocated.control <- relocate_loop_pairs(df.observed.metric) %>%
    measure_non_tss_atac()

  bind_rows(
    summarise_one_permutation(
      df.observed = df.observed.metric,
      df.control = df.matched.control,
      permutation = permutation.i,
      null.method = "matched_HiC_anchor"
    ),
    summarise_one_permutation(
      df.observed = df.observed.metric,
      df.control = df.relocated.control,
      permutation = permutation.i,
      null.method = "rigid_loop_pair_relocation"
    )
  )
}

message(
  "Running ", n.permutations,
  " ATAC matched-null permutations on ", n.cores, " core(s)."
)
permutation.results <- if (.Platform$OS.type == "windows") {
  lapply(seq_len(n.permutations), run_one_atac_null_permutation)
} else {
  parallel::mclapply(
    seq_len(n.permutations),
    run_one_atac_null_permutation,
    mc.cores = n.cores,
    mc.preschedule = FALSE
  )
}
df.atac.null.permutation <- bind_rows(permutation.results)

# Relocate each no-direct-TSS loop as an intact pair. This comparison directly
# tests whether real no-direct-TSS Hi-C loops overlap open chromatin
# more often than random genomic loop placements with the same geometry.
run_one_no_direct_tss_permutation <- function(permutation.i) {
  set.seed(random.seed + 100000L + permutation.i)
  relocate_no_direct_tss_loop_pairs(df.no.direct.tss.loops) %>%
    measure_non_tss_atac() %>%
    summarise_no_direct_tss_loop_atac(permutation = permutation.i)
}

message(
  "Running ", n.permutations,
  " no-direct-TSS loop relocation permutations on ", n.cores,
  " core(s)."
)
no.direct.tss.permutation.results <- if (.Platform$OS.type == "windows") {
  lapply(seq_len(n.permutations), run_one_no_direct_tss_permutation)
} else {
  parallel::mclapply(
    seq_len(n.permutations),
    run_one_no_direct_tss_permutation,
    mc.cores = n.cores,
    mc.preschedule = FALSE
  )
}
df.no.direct.tss.atac.random <- bind_rows(
  no.direct.tss.permutation.results
)

################################################################################
# 6. Summarise enrichment, matched odds ratios, and empirical P-values
################################################################################

df.atac.matched.null.summary <- df.atac.null.permutation %>%
  group_by(null_method, resolution, metric) %>%
  summarise(
    n_permutations = n(),
    n_pairs = dplyr::first(n_pairs),
    observed_rate = dplyr::first(observed_rate),
    mean_null_rate = mean(null_rate),
    sd_null_rate = sd(null_rate),
    null_rate_q025 = quantile(null_rate, 0.025),
    null_rate_q975 = quantile(null_rate, 0.975),
    absolute_rate_difference = observed_rate - mean_null_rate,
    enrichment_ratio = observed_rate / mean_null_rate,
    median_paired_odds_ratio = median(paired_odds_ratio),
    paired_odds_ratio_q025 = quantile(paired_odds_ratio, 0.025),
    paired_odds_ratio_q975 = quantile(paired_odds_ratio, 0.975),
    empirical_p_greater_equal = (
      1 + sum(null_rate >= observed_rate)
    ) / (
      n() + 1
    ),
    .groups = "drop"
  ) %>%
  arrange(
    null_method,
    factor(resolution, levels = c("5K", "10K", "25K", "ALL")),
    metric
  )

# Contrast the observed no-direct-TSS loops with their geometry-preserving
# randomized placements for anchor-level, any-anchor, and both-anchor support.
df.no.direct.tss.atac.random.summary <- df.no.direct.tss.atac.random %>%
  group_by(resolution, metric) %>%
  summarise(
    n_permutations = n(),
    n_loops = dplyr::first(n_loops),
    n_anchors = dplyr::first(n_anchors),
    across(
      c(
        anchor_positive_rate,
        any_anchor_positive_rate,
        both_anchors_positive_rate
      ),
      list(
        mean = mean,
        q025 = ~ quantile(.x, 0.025),
        q975 = ~ quantile(.x, 0.975)
      )
    ),
    .groups = "drop"
  ) %>%
  left_join(
    df.no.direct.tss.atac.actual %>%
      dplyr::select(-permutation) %>%
      rename_with(
        ~ paste0("observed_", .x),
        c(
          "anchor_positive_rate",
          "any_anchor_positive_rate",
          "both_anchors_positive_rate"
        )
      ),
    by = c("resolution", "metric", "n_loops", "n_anchors")
  ) %>%
  mutate(
    anchor_enrichment_ratio = observed_anchor_positive_rate /
      anchor_positive_rate_mean,
    any_anchor_enrichment_ratio = observed_any_anchor_positive_rate /
      any_anchor_positive_rate_mean,
    both_anchors_enrichment_ratio = observed_both_anchors_positive_rate /
      both_anchors_positive_rate_mean
  )

df.no.direct.tss.atac.random.summary <- df.no.direct.tss.atac.random.summary %>%
  rowwise() %>%
  mutate(
    anchor_empirical_p_greater_equal = (
      1 + sum(
        df.no.direct.tss.atac.random$anchor_positive_rate[
          df.no.direct.tss.atac.random$resolution == resolution &
            df.no.direct.tss.atac.random$metric == metric
        ] >= observed_anchor_positive_rate
      )
    ) / (n_permutations + 1),
    any_anchor_empirical_p_greater_equal = (
      1 + sum(
        df.no.direct.tss.atac.random$any_anchor_positive_rate[
          df.no.direct.tss.atac.random$resolution == resolution &
            df.no.direct.tss.atac.random$metric == metric
        ] >= observed_any_anchor_positive_rate
      )
    ) / (n_permutations + 1),
    both_anchors_empirical_p_greater_equal = (
      1 + sum(
        df.no.direct.tss.atac.random$both_anchors_positive_rate[
          df.no.direct.tss.atac.random$resolution == resolution &
            df.no.direct.tss.atac.random$metric == metric
        ] >= observed_both_anchors_positive_rate
      )
    ) / (n_permutations + 1)
  ) %>%
  ungroup() %>%
  arrange(
    factor(resolution, levels = c("5K", "10K", "25K", "ALL")),
    metric
  )

# Put the three controls discussed with the PI in one resolution-stratified
# table. Anchor-level rates are directly comparable because each row evaluates
# one anchor rather than giving no-TSS loops two chances to overlap ATAC.
df.atac.three_way.anchor.comparison <- df.atac.matched.null.summary %>%
  filter(
    metric == "atac_overlap_ge50",
    null_method == "matched_HiC_anchor"
  ) %>%
  dplyr::select(
    resolution,
    promoter_supported_opposite_anchor_rate = observed_rate,
    matched_real_no_direct_tss_anchor_rate = mean_null_rate
  ) %>%
  left_join(
    df.atac.matched.null.summary %>%
      filter(
        metric == "atac_overlap_ge50",
        null_method == "rigid_loop_pair_relocation"
      ) %>%
      dplyr::select(
        resolution,
        promoter_supported_random_relocation_anchor_rate = mean_null_rate
      ),
    by = "resolution"
  ) %>%
  left_join(
    df.no.direct.tss.atac.actual %>%
      filter(metric == "atac_overlap_ge50") %>%
      dplyr::select(
        resolution,
        no_direct_tss_real_anchor_rate = anchor_positive_rate,
        no_direct_tss_real_loop_any_anchor_rate = any_anchor_positive_rate,
        no_direct_tss_real_loop_both_anchors_rate = both_anchors_positive_rate
      ),
    by = "resolution"
  ) %>%
  left_join(
    df.no.direct.tss.atac.random.summary %>%
      filter(metric == "atac_overlap_ge50") %>%
      dplyr::select(
        resolution,
        no_direct_tss_random_anchor_rate = anchor_positive_rate_mean,
        no_direct_tss_random_loop_any_anchor_rate =
          any_anchor_positive_rate_mean,
        no_direct_tss_random_loop_both_anchors_rate =
          both_anchors_positive_rate_mean
      ),
    by = "resolution"
  ) %>%
  arrange(factor(resolution, levels = c("5K", "10K", "25K", "ALL")))

# Record exact matching coverage independently from random control selection.
df.atac.matched.hic.control.quality <- df.match.availability %>%
  summarise(
    n_observed_with_exact_match = sum(
      n_observed_anchors[exact_match_available]
    ),
    n_observed_requiring_side_relaxed_fallback =
      sum(n_observed_anchors[!exact_match_available & fallback_match_available]),
    n_observed_requiring_distance_bin_relaxation = sum(
      n_observed_anchors[!fallback_match_available]
    ),
    n_match_strata = n(),
    n_exact_match_strata = sum(exact_match_available),
    n_fallback_match_strata = sum(!exact_match_available),
    n_observed_anchors = sum(n_observed_anchors)
  ) %>%
  mutate(
    pct_observed_with_exact_match = round(
      100 * n_observed_with_exact_match / n_observed_anchors,
      3
    )
  )

################################################################################
# 7. Save auditable tables and the resolution-stratified comparison figure
################################################################################

df.atac.matched.null.method <- tribble(
  ~analysis_item, ~definition,
  "analysis_population",
  paste0(
    "All pooled loops with direct promoter/TSS support at exactly one anchor; ",
    "the opposite anchor is evaluated regardless of its ATAC status."
  ),
  "ATAC_source",
  "Duttke et al. 2022 rat prefrontal-cortex snATAC peaks lifted to rn7.",
  "TSS_exclusion",
  paste0(
    "Union of strand-aware Ensembl transcript TSS and EPD TSS positions, ",
    "expanded by +/-1 kb and removed from the reduced ATAC union."
  ),
  "primary_absolute_overlap",
  ">=50 bp TSS-excluded ATAC overlap at the candidate anchor.",
  "threshold_sensitivity",
  paste0(
    ">=1 bp, >=50 bp, >=0.5%, and >=1% of available non-TSS anchor bases."
  ),
  "rigid_loop_pair_relocation",
  paste0(
    "The complete loop geometry is relocated uniformly within the same ",
    "chromosome, preserving resolution, anchor widths, and loop distance."
  ),
  "matched_HiC_anchor",
  paste0(
    "Promoter-free pooled Hi-C anchors matched on chromosome, resolution, ",
    "anchor side, and loop-distance bin; side is relaxed only if necessary."
  ),
  "empirical_p_value",
  "(1 + permutations with null rate >= observed rate) / (B + 1).",
  "interpretation",
  paste0(
    "ATAC enrichment supports non-TSS open chromatin but does not validate ",
    "enhancer function or a promoter-enhancer interaction."
  )
)

plot.atac.null <- df.atac.matched.null.summary %>%
  filter(metric == "atac_overlap_ge50") %>%
  mutate(
    resolution = factor(resolution, levels = c("5K", "10K", "25K", "ALL")),
    null_method = recode(
      null_method,
      matched_HiC_anchor = "Matched promoter-free Hi-C anchors",
      rigid_loop_pair_relocation = "Rigid within-chromosome loop relocation"
    )
  ) %>%
  ggplot(aes(x = resolution)) +
  geom_linerange(
    aes(ymin = null_rate_q025, ymax = null_rate_q975),
    color = "grey45",
    linewidth = 0.8
  ) +
  geom_point(
    aes(y = mean_null_rate, shape = "Null mean"),
    color = "grey25",
    size = 2.5
  ) +
  geom_point(
    aes(y = observed_rate, shape = "Observed"),
    color = "#B2182B",
    size = 2.8
  ) +
  facet_wrap(vars(null_method)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  scale_shape_manual(values = c("Null mean" = 1, "Observed" = 16)) +
  labs(
    x = "HiCCUPS resolution",
    y = "Candidate anchors with >=50 bp non-TSS ATAC overlap",
    shape = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

# Slide-ready summary of the primary matched-anchor null model. This plot uses
# the analysis results above directly; no values are entered manually.
plot.atac.matched.null.dumbbell <- df.atac.matched.null.summary %>%
  filter(
    null_method == "matched_HiC_anchor",
    metric == "atac_overlap_ge50"
  ) %>%
  mutate(
    resolution = factor(resolution, levels = c("5K", "10K", "25K", "ALL")),
    resolution_label = recode(
      as.character(resolution),
      `5K` = "5 kb",
      `10K` = "10 kb",
      `25K` = "25 kb",
      ALL = "All"
    ),
    resolution_label = factor(
      resolution_label,
      levels = c("All", "25 kb", "10 kb", "5 kb")
    ),
    observed_pct = 100 * observed_rate,
    null_pct = 100 * mean_null_rate,
    null_low_pct = 100 * null_rate_q025,
    null_high_pct = 100 * null_rate_q975
  ) %>%
  ggplot(aes(y = resolution_label)) +
  geom_segment(
    aes(x = null_low_pct, xend = null_high_pct, yend = resolution_label),
    color = "#8A99A6",
    linewidth = 1.1,
    lineend = "round"
  ) +
  geom_segment(
    aes(x = null_pct, xend = observed_pct, yend = resolution_label),
    color = "#D8A62A",
    linewidth = 1.5,
    lineend = "round"
  ) +
  geom_point(aes(x = null_pct), size = 3.8, color = "#8A99A6") +
  geom_point(aes(x = observed_pct), size = 4.2, color = "#1C918D") +
  geom_text(
    aes(x = null_low_pct - 0.35, label = sprintf("%.1f%%", null_pct)),
    hjust = 1,
    color = "#687783",
    fontface = "bold",
    size = 3.8
  ) +
  geom_text(
    aes(x = observed_pct + 0.55, label = sprintf("%.1f%%", observed_pct)),
    hjust = 0,
    color = "#0B7773",
    fontface = "bold",
    size = 3.8
  ) +
  annotate(
    "point",
    x = 69.5,
    y = 4.65,
    color = "#8A99A6",
    size = 3.2
  ) +
  annotate(
    "text",
    x = 70,
    y = 4.65,
    label = "Matched null",
    hjust = 0,
    color = "#8A99A6",
    fontface = "bold",
    size = 3.7
  ) +
  annotate(
    "point",
    x = 77.5,
    y = 4.65,
    color = "#1C918D",
    size = 3.5
  ) +
  annotate(
    "text",
    x = 78,
    y = 4.65,
    label = "Observed",
    hjust = 0,
    color = "#1C918D",
    fontface = "bold",
    size = 3.7
  ) +
  scale_x_continuous(
    limits = c(68, 93),
    breaks = c(70, 75, 80, 85, 90),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#DFE6EC", linewidth = 0.55),
    axis.text.x = element_text(color = "#657584", size = 10.5),
    axis.text.y = element_text(color = "#123E6A", face = "bold", size = 11.5),
    plot.margin = margin(16, 32, 8, 8),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

output.tables <- list(
  revised_atac_threshold_sensitivity_by_resolution =
    df.atac.threshold.sensitivity.by.resolution,
  revised_atac_matched_null_summary = df.atac.matched.null.summary,
  revised_atac_matched_null_permutation = df.atac.null.permutation,
  revised_atac_matched_null_method = df.atac.matched.null.method,
  revised_atac_matched_control_quality = df.atac.matched.hic.control.quality,
  revised_atac_matched_control_strata = df.match.availability,
  revised_atac_no_direct_tss_loop_status = df.no.direct.tss.loop.status,
  revised_atac_no_direct_tss_observed = df.no.direct.tss.atac.actual,
  revised_atac_no_direct_tss_random_permutation =
    df.no.direct.tss.atac.random,
  revised_atac_no_direct_tss_random_summary =
    df.no.direct.tss.atac.random.summary,
  revised_atac_three_way_anchor_comparison =
    df.atac.three_way.anchor.comparison,
  revised_atac_matched_null_run_metadata = tibble(
    analysis_release = "resubmission-2026-07-28",
    n_permutations = n.permutations,
    random_seed = random.seed,
    n_cores = n.cores,
    primary_overlap_threshold_bp = atac.minimum.overlap.bp,
    tss_exclusion_flank_bp = tss.exclusion.flank.bp,
    input_loop_resource = loop.resource.file,
    input_direct_assignment = direct.assignment.file,
    input_atac_cache = atac.cache.file,
    completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
)

walk2(names(output.tables), output.tables, function(file.stem, table) {
  write_tsv(table, file.path(output.dir, paste0(file.stem, ".tsv")))
})

ggsave(
  file.path(output.dir, "revised_atac_observed_vs_matched_null.pdf"),
  plot.atac.null,
  width = 8,
  height = 4.2,
  device = "pdf"
)
ggsave(
  file.path(output.dir, "revised_atac_observed_vs_matched_null.png"),
  plot.atac.null,
  width = 8,
  height = 4.2,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(
    output.dir,
    "revised_atac_matched_null_dumbbell_by_resolution.pdf"
  ),
  plot.atac.matched.null.dumbbell,
  width = 6.35,
  height = 3.05,
  device = "pdf",
  bg = "white"
)
ggsave(
  file.path(
    output.dir,
    "revised_atac_matched_null_dumbbell_by_resolution.png"
  ),
  plot.atac.matched.null.dumbbell,
  width = 6.35,
  height = 3.05,
  dpi = 300,
  bg = "white"
)

writeLines(
  sub("[[:blank:]]+$", "", capture.output(sessionInfo())),
  file.path(output.dir, "revised_atac_matched_null_session_info.txt")
)

# Record every official output with its size and SHA-256 checksum.
release.output.files <- setdiff(
  list.files(output.dir, all.files = FALSE, no.. = TRUE),
  "resubmit_output_manifest.tsv"
)
release.output.paths <- file.path(output.dir, release.output.files)
df.output.manifest <- tibble(
  output_file = release.output.files,
  output_path = file.path("results", release.output.files),
  file_exists = file.exists(release.output.paths),
  file_size_bytes = as.numeric(file.info(release.output.paths)$size),
  sha256 = map_chr(
    release.output.paths,
    ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)
  ),
  generated_by = "atac_validation.R",
  analysis_release = "resubmission-2026-07-28",
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
)
write_tsv(
  df.output.manifest,
  file.path(output.dir, "resubmit_output_manifest.tsv")
)

message("ATAC overlap threshold sensitivity by resolution:")
print(df.atac.threshold.sensitivity.by.resolution)
message("ATAC matched-null enrichment summary:")
print(df.atac.matched.null.summary)
message("Matched Hi-C control availability:")
print(df.atac.matched.hic.control.quality)
message("No-direct-TSS loops versus randomized loop placements:")
print(df.no.direct.tss.atac.random.summary)
message("Three-way ATAC anchor comparison:")
print(df.atac.three_way.anchor.comparison)
message("Resubmission ATAC matched-null analysis completed.")
