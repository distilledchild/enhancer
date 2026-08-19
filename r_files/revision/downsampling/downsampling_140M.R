# lintr: disable

# Compare the original full-depth HiCCUPS calls with the MAPQ >=30, 140M-contact
# downsampled calls. The comparison is restricted to samples with complete
# downsampling QC, .hic provenance, and 5/10/25-kb HiCCUPS outputs.

current_script_path <- function() {
  file.args <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )
  if (length(file.args) > 0L) {
    file.path <- gsub(
      "~\\+~", " ", sub("^--file=", "", file.args[[1]])
    )
    return(normalizePath(
      file.path,
      winslash = "/",
      mustWork = FALSE
    ))
  }

  frame.files <- vapply(
    sys.frames(),
    function(frame) {
      if (is.null(frame$ofile)) NA_character_ else as.character(frame$ofile)
    },
    character(1L)
  )
  frame.files <- frame.files[!is.na(frame.files) & nzchar(frame.files)]
  if (length(frame.files) > 0L) {
    return(normalizePath(
      tail(frame.files, 1L),
      winslash = "/",
      mustWork = FALSE
    ))
  }
  NA_character_
}

# Locate the shared analysis directory from the script, an environment variable,
# or the same local/Dropbox layouts used by the resubmission pipeline.
resolve_enhancer_r_files_dir <- function() {
  script.path <- current_script_path()
  start.dirs <- c(
    if (is.na(script.path)) NA_character_ else dirname(script.path),
    getwd()
  )
  ancestor.dirs <- unique(unlist(lapply(start.dirs, function(path) {
    if (is.na(path) || !nzchar(path)) return(character())
    path <- normalizePath(path, winslash = "/", mustWork = FALSE)
    paths <- path
    for (i in seq_len(6L)) {
      parent <- dirname(tail(paths, 1L))
      if (identical(parent, tail(paths, 1L))) break
      paths <- c(paths, parent)
    }
    paths
  })))
  candidates <- unique(c(
    Sys.getenv("ENHANCER_R_FILES_DIR", unset = ""),
    ancestor.dirs,
    path.expand("~/dropbox/Gateway_to_Hao/enhancer/r_files"),
    path.expand("~/Dropbox/Gateway_to_Hao/enhancer/r_files"),
    Sys.glob(path.expand(
      "~/Library/CloudStorage/Dropbox*/K P/Gateway_to_Hao/enhancer/r_files"
    )),
    Sys.glob(path.expand(
      "~/Library/CloudStorage/Dropbox*/Gateway_to_Hao/enhancer/r_files"
    ))
  ))
  candidates <- candidates[
    !is.na(candidates) & nzchar(candidates) & dir.exists(candidates)
  ]
  candidates <- candidates[file.exists(file.path(candidates, "funcs.R"))]
  if (length(candidates) == 0L) {
    stop(
      "Cannot locate enhancer/r_files with funcs.R. Set ENHANCER_R_FILES_DIR.",
      call. = FALSE
    )
  }
  normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE)
}

r.files.dir <- resolve_enhancer_r_files_dir()
enhancer.project.dir <- dirname(r.files.dir)
source(file.path(r.files.dir, "funcs.R"))
script.path <- current_script_path()
analysis.dir <- if (is.na(script.path)) {
  file.path(r.files.dir, "revision", "downsampling")
} else {
  dirname(script.path)
}

library("tidyverse")
library("GenomicRanges")
library("GenomeInfoDb")

options(tibble.width = Inf)
options(tibble.print_max = Inf)
options(tibble.max_extra_cols = Inf)
options(scipen = 999)

################################################################################
# Downsampling sensitivity analysis
#
# 1. Resolve complete full-depth and 140M HiCCUPS inputs.
# 2. Rebuild exact pooled calls from the same completed sample set.
# 3. Quantify loop-count dispersion and exact loop recovery by resolution.
# 4. Reapply promoter/TSS +/-1-kb and true-TSS-excluded ATAC annotations.
# 5. Compare category composition and putative-regulatory gene-loop counts.
# 6. Test whether pooled calls are disproportionately supported by high-depth
#    source libraries.
################################################################################

########################
# 0. Directories, sample metadata, and constants
########################

revision.dir <- file.path(r.files.dir, "revision")
input.dir <- file.path(analysis.dir, "inputs")
target.contacts <- 140000000L
max.loop.distance.bp <- 2000000L
promoter.window.flank.bp <- 1000L
atac.minimum.overlap.bp <- 50L

# Keep the exact sample/strain/contact/seed provenance used for downsampling.
df.sample.metadata <- tribble(
  ~sample, ~strain, ~source_contacts, ~seed,
  "592BB", "SHR/OlaIpcv", 259559222, 20260731,
  "607", "HXB10", 419977942, 20260737,
  "74AA", "F344/Stm", 141993330, 20260728,
  "A2DB", "LE/Stm", 145939027, 20260729,
  "D765A", "BXH6", 284296675, 20260732,
  "DA08A", "HXB2", 213727006, 20260730,
  "DA21A", "SHR/OlaIpcvxBN/NHsdMcwi", 382917950, 20260736,
  "DA68A", "HXB31", 378608379, 20260734,
  "DBA9A", "HXB23", 380788024, 20260735,
  "DE8BA", "BN-Lx", 302314076, 20260733
) %>%
  mutate(
    target_contacts = target.contacts,
    sampling_fraction = target_contacts / source_contacts
  )

# Select the first existing directory from an ordered candidate list.
resolve_existing_directory <- function(candidates, label) {
  candidates <- unique(path.expand(candidates[nzchar(candidates)]))
  candidates <- candidates[dir.exists(candidates)]
  if (length(candidates) == 0L) {
    stop("Cannot locate ", label, ".", call. = FALSE)
  }
  normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE)
}

# Prefer the candidate containing the most nonempty merged HiCCUPS outputs.
resolve_downsample_root <- function(candidates) {
  candidates <- unique(path.expand(candidates[nzchar(candidates)]))
  candidates <- candidates[dir.exists(candidates)]
  if (length(candidates) == 0L) {
    stop(
      "Cannot locate downsampled HiCCUPS results. Set DOWNSAMPLED_HICCUPS_ROOT.",
      call. = FALSE
    )
  }
  n.merged <- vapply(
    candidates,
    function(path) {
      files <- list.files(
        path,
        pattern = "merged_loops[.]bedpe$",
        recursive = TRUE,
        full.names = TRUE
      )
      sum(file.exists(files) & file.info(files)$size > 0L)
    },
    integer(1L)
  )
  normalizePath(
    candidates[[which.max(n.merged)]],
    winslash = "/",
    mustWork = TRUE
  )
}

full.depth.root <- resolve_existing_directory(
  c(
    Sys.getenv("FULL_DEPTH_HICCUPS_ROOT", unset = ""),
    file.path(input.dir, "hic", "2023A", "hic30_w_sb_options"),
    file.path(
      enhancer.project.dir, "hic", "2023A", "hic30_w_sb_options"
    ),
    path.expand(
      "~/dropbox/Gateway_to_Hao/hic/2023A/hic30_w_sb_options"
    ),
    Sys.glob(path.expand(
      paste0(
        "~/Library/CloudStorage/Dropbox*/K P/Gateway_to_Hao/",
        "hic/2023A/hic30_w_sb_options"
      )
    )),
    Sys.glob(path.expand(
      paste0(
        "~/Library/CloudStorage/Dropbox*/Gateway_to_Hao/",
        "hic/2023A/hic30_w_sb_options"
      )
    ))
  ),
  "full-depth sb-option HiCCUPS directory"
)

downsample.root <- resolve_downsample_root(c(
  Sys.getenv("DOWNSAMPLED_HICCUPS_ROOT", unset = ""),
  file.path(input.dir, "juicer_downsample_q30_140M_250M", "140M"),
  file.path(
    enhancer.project.dir,
    "data", "juicer_downsample_q30_140M_250M", "140M"
  ),
  path.expand(
    "~/dropbox/Gateway_to_Hao/enhancer/data/juicer_downsample_q30_140M"
  ),
  file.path(dirname(r.files.dir), "data", "juicer_downsample_q30_140M")
))

required.cache.files <- c(
  "df.transcript.ensembl.rn7.1based.rds",
  "df.promoter.epd.rn7.1based.rds",
  "gr.atac.rn7.1based.rds"
)
coord.cache.candidates <- unique(path.expand(c(
  Sys.getenv("DOWNSAMPLING_COORD_CACHE_DIR", unset = ""),
  file.path(revision.dir, "revision_main", "cache_data"),
  "~/dropbox/Gateway_to_Hao/enhancer/r_files/revision/revision_main/cache_data"
)))
coord.cache.candidates <- coord.cache.candidates[
  nzchar(coord.cache.candidates) & dir.exists(coord.cache.candidates)
]
coord.cache.complete <- vapply(
  coord.cache.candidates,
  function(path) all(file.exists(file.path(path, required.cache.files))),
  logical(1L)
)
if (!any(coord.cache.complete)) {
  stop(
    "Cannot locate a complete coordinate-normalization cache. Set ",
    "DOWNSAMPLING_COORD_CACHE_DIR.",
    call. = FALSE
  )
}
coord.cache.dir <- normalizePath(
  coord.cache.candidates[coord.cache.complete][[1]],
  winslash = "/",
  mustWork = TRUE
)

output.dir <- path.expand(Sys.getenv(
  "DOWNSAMPLING_140M_OUTPUT_DIR",
  unset = file.path(analysis.dir, "results", "full_140M")
))
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

message("Full-depth HiCCUPS root: ", full.depth.root)
message("Downsampled HiCCUPS root: ", downsample.root)
message("Coordinate cache: ", coord.cache.dir)
message("Output directory: ", output.dir)

########################
# 0-1. Local analysis helpers
########################

# Return the expected downsampled result paths and require all three resolution
# outputs before treating a sample as complete.
build_sample_file_status <- function(df.metadata, full.root, down.root) {
  df.metadata %>%
    mutate(
      full_file = file.path(
        full.root,
        sample,
        "hiccups_5k10k25k",
        "merged_loops.bedpe"
      ),
      downsample_result_dir = file.path(
        down.root,
        sample,
        str_c(
          sample,
          "_inter_30_140M_seed",
          seed,
          ".hiccups.5k10k25k"
        )
      ),
      downsample_file = file.path(
        downsample_result_dir,
        "merged_loops.bedpe"
      ),
      downsampling_qc_file = file.path(
        down.root, sample, "downsampling_qc.tsv"
      ),
      hic_creation_qc_file = file.path(
        down.root, sample, "hic_creation_qc.tsv"
      ),
      hic_file = file.path(
        down.root,
        sample,
        str_c(sample, "_inter_30_140M_seed", seed, ".hic")
      ),
      downsample_postprocessed_5k = file.path(
        downsample_result_dir,
        "postprocessed_pixels_5000.bedpe"
      ),
      downsample_postprocessed_10k = file.path(
        downsample_result_dir,
        "postprocessed_pixels_10000.bedpe"
      ),
      downsample_postprocessed_25k = file.path(
        downsample_result_dir,
        "postprocessed_pixels_25000.bedpe"
      ),
      full_file_ready = file.exists(full_file) & file.info(full_file)$size > 0L,
      downsample_merged_ready = file.exists(downsample_file) &
        file.info(downsample_file)$size > 0L,
      downsample_5k_ready = file.exists(downsample_postprocessed_5k) &
        file.info(downsample_postprocessed_5k)$size > 0L,
      downsample_10k_ready = file.exists(downsample_postprocessed_10k) &
        file.info(downsample_postprocessed_10k)$size > 0L,
      downsample_25k_ready = file.exists(downsample_postprocessed_25k) &
        file.info(downsample_postprocessed_25k)$size > 0L,
      provenance_ready = file.exists(downsampling_qc_file) &
        file.info(downsampling_qc_file)$size > 0L &
        file.exists(hic_creation_qc_file) &
        file.info(hic_creation_qc_file)$size > 0L &
        file.exists(hic_file) & file.info(hic_file)$size > 0L,
      analysis_ready = full_file_ready & downsample_merged_ready &
        downsample_5k_ready & downsample_10k_ready & downsample_25k_ready,
      status = if_else(
        analysis_ready & provenance_ready,
        "included_complete_5k10k25k",
        if_else(
          analysis_ready,
          "analysis_ready_but_provenance_incomplete",
          "excluded_until_complete_5k10k25k"
        )
      )
    )
}

# Add an ALL-resolution row without changing the original resolution records.
add_all_resolution <- function(df) {
  bind_rows(df, df %>% mutate(resolution = "ALL"))
}

# Calculate a coefficient of variation only when its denominator is positive.
coefficient_of_variation <- function(x) {
  if (length(x) < 2L || mean(x) <= 0) NA_real_ else sd(x) / mean(x)
}

# Return a safe Spearman estimate without failing on constant or short vectors.
safe_spearman <- function(x, y) {
  valid <- is.finite(x) & is.finite(y)
  if (sum(valid) < 3L || length(unique(x[valid])) < 2L ||
      length(unique(y[valid])) < 2L) {
    return(NA_real_)
  }
  suppressWarnings(cor(x[valid], y[valid], method = "spearman"))
}

# Return the sample size and asymptotic P value for a descriptive Spearman test.
safe_spearman_test <- function(x, y) {
  valid <- is.finite(x) & is.finite(y)
  n.valid <- sum(valid)
  if (n.valid < 3L || length(unique(x[valid])) < 2L ||
      length(unique(y[valid])) < 2L) {
    return(tibble(n_spearman = n.valid, spearman_p_value = NA_real_))
  }
  result <- suppressWarnings(cor.test(
    x[valid], y[valid], method = "spearman", exact = FALSE
  ))
  tibble(n_spearman = n.valid, spearman_p_value = result$p.value)
}

# Measure position-tolerant recovery without changing the exact-call resource.
# The 20/20/50-kb thresholds are the HiCCUPS radii configured for 5/10/25 kb.
summarise_hiccups_radius_recovery <- function(df.full, df.downsample) {
  radius.by.resolution <- c("5K" = 20000L, "10K" = 20000L, "25K" = 50000L)
  comparison.groups <- df.full %>%
    filter(passes_lt2mb) %>%
    distinct(sample, resolution) %>%
    arrange(sample, resolution)

  df.by.resolution <- pmap_dfr(
    comparison.groups,
    function(sample, resolution) {
      df.full.i <- df.full %>%
        filter(
          passes_lt2mb,
          .data$sample == .env$sample,
          .data$resolution == .env$resolution
        ) %>%
        mutate(call_index = row_number())
      df.down.i <- df.downsample %>%
        filter(
          passes_lt2mb,
          .data$sample == .env$sample,
          .data$resolution == .env$resolution
        ) %>%
        mutate(call_index = row_number())

      match.radius.bp <- unname(radius.by.resolution[[resolution]])
      if (nrow(df.full.i) == 0L || nrow(df.down.i) == 0L) {
        return(tibble(
          sample = sample,
          resolution = resolution,
          hiccups_match_radius_bp = match.radius.bp,
          n_full_depth = nrow(df.full.i),
          n_downsample_140M = nrow(df.down.i),
          n_full_depth_radius_recovered = 0L,
          n_downsample_radius_matched = 0L
        ))
      }

      common.seqinfo <- Seqinfo(seqnames = union(
        unique(df.full.i$chr1), unique(df.down.i$chr1)
      ))
      gr.full.anchor1 <- GRanges(
        seqnames = df.full.i$chr1,
        ranges = IRanges(
          start = df.full.i$hiccups_centroid1_1based,
          end = df.full.i$hiccups_centroid1_1based
        ),
        seqinfo = common.seqinfo
      )
      gr.down.anchor1.window <- GRanges(
        seqnames = df.down.i$chr1,
        ranges = IRanges(
          start = pmax(
            1L,
            df.down.i$hiccups_centroid1_1based - match.radius.bp
          ),
          end = df.down.i$hiccups_centroid1_1based + match.radius.bp
        ),
        seqinfo = common.seqinfo
      )
      hits <- findOverlaps(
        gr.full.anchor1,
        gr.down.anchor1.window,
        ignore.strand = TRUE
      )

      if (length(hits) == 0L) {
        n.full.recovered <- 0L
        n.down.matched <- 0L
      } else {
        full.index <- queryHits(hits)
        down.index <- subjectHits(hits)
        valid.anchor2 <- (
          df.full.i$chr2[full.index] == df.down.i$chr2[down.index] &
            abs(
              df.full.i$hiccups_centroid2_1based[full.index] -
                df.down.i$hiccups_centroid2_1based[down.index]
            ) <= match.radius.bp
        )
        n.full.recovered <- n_distinct(full.index[valid.anchor2])
        n.down.matched <- n_distinct(down.index[valid.anchor2])
      }

      tibble(
        sample = sample,
        resolution = resolution,
        hiccups_match_radius_bp = match.radius.bp,
        n_full_depth = nrow(df.full.i),
        n_downsample_140M = nrow(df.down.i),
        n_full_depth_radius_recovered = n.full.recovered,
        n_downsample_radius_matched = n.down.matched
      )
    }
  ) %>%
    mutate(
      pct_full_depth_radius_recovered = 100 *
        n_full_depth_radius_recovered / n_full_depth,
      pct_downsample_radius_matched = 100 *
        n_downsample_radius_matched / n_downsample_140M
    )

  df.all <- df.by.resolution %>%
    group_by(sample) %>%
    summarise(
      resolution = "ALL",
      hiccups_match_radius_bp = NA_integer_,
      n_full_depth = sum(n_full_depth),
      n_downsample_140M = sum(n_downsample_140M),
      n_full_depth_radius_recovered = sum(n_full_depth_radius_recovered),
      n_downsample_radius_matched = sum(n_downsample_radius_matched),
      pct_full_depth_radius_recovered = 100 *
        n_full_depth_radius_recovered / n_full_depth,
      pct_downsample_radius_matched = 100 *
        n_downsample_radius_matched / n_downsample_140M,
      .groups = "drop"
    )

  bind_rows(df.by.resolution, df.all) %>%
    arrange(resolution, sample)
}

# Build the normalized sample-level and exact pooled resources for one depth.
build_depth_loop_resource <- function(file.metadata, condition) {
  df.sample.raw <- read_hiccups_loop_files(file.metadata)
  df.sample.normalized <- normalize_hiccups_loop_coordinates(
    df.sample.raw,
    max.loop.distance = max.loop.distance.bp
  ) %>%
    mutate(condition = condition, .before = 1)

  pooled <- build_pooled_hiccups_loop_resource(df.sample.normalized)
  list(
    sample = df.sample.normalized,
    support = pooled$support,
    distinct = pooled$distinct,
    universe = pooled$universe %>% mutate(condition = condition, .before = 1)
  )
}

# Build promoter/TSS and residual non-TSS ATAC features once for both depths.
build_regulatory_feature_context <- function(
  df.transcript,
  df.promoter,
  gr.atac.rn7.1based.input,
  promoter.flank.bp = 1000L
) {
  df.true.tss <- build_true_tss_annotation(df.transcript)
  gr.true.tss <- create_true_tss_granges(df.true.tss)
  gr.epd.tss <- create_epd_tss_granges(df.promoter)

  gr.true.window <- expand_tss_to_promoter_windows(
    gr.true.tss,
    flank.bp = promoter.flank.bp
  )
  gr.epd.window <- expand_tss_to_promoter_windows(
    gr.epd.tss,
    flank.bp = promoter.flank.bp
  )

  gr.known.tss <- c(gr.true.tss, gr.epd.tss)
  gr.tss.exclusion <- GRanges(
    seqnames = seqnames(gr.known.tss),
    ranges = IRanges(
      start = pmax(1L, start(gr.known.tss) - promoter.flank.bp),
      end = end(gr.known.tss) + promoter.flank.bp
    )
  ) %>%
    reduce(ignore.strand = TRUE)

  gr.atac.union <- reduce(gr.atac.rn7.1based.input, ignore.strand = TRUE)
  common.seqlevels <- intersect(seqlevels(gr.atac.union), seqlevels(gr.tss.exclusion))
  gr.atac.common <- keepSeqlevels(
    gr.atac.union,
    common.seqlevels,
    pruning.mode = "coarse"
  )
  gr.exclusion.common <- keepSeqlevels(
    gr.tss.exclusion,
    common.seqlevels,
    pruning.mode = "coarse"
  )
  gr.atac.non.tss <- setdiff(
    gr.atac.common,
    gr.exclusion.common,
    ignore.strand = TRUE
  )

  list(
    df_true_tss = df.true.tss,
    gr_true_window = gr.true.window,
    gr_epd_window = gr.epd.window,
    gr_atac_union = gr.atac.common,
    gr_atac_non_tss = gr.atac.non.tss
  )
}

# Convert promoter-window hits to the lossless evidence columns required by the
# shared direct promoter/TSS summarizer.
build_direct_promoter_evidence <- function(
  gr.anchor.by.side,
  feature.context,
  df.promoter
) {
  df.true.lookup <- feature.context$df_true_tss %>%
    mutate(annotation_index = row_number()) %>%
    dplyr::select(
      annotation_index,
      true_tss_id,
      gene_id,
      gene_name,
      transcript_id_versioned
    )

  df.epd.lookup <- df.promoter %>%
    mutate(annotation_index = row_number()) %>%
    dplyr::select(
      annotation_index,
      promoter_annotation_id,
      gene_id,
      gene_name
    )

  df.true <- map_loop_anchors_dfr(
    gr.anchor.by.side,
    find_direct_anchor_annotation_overlaps,
    gr.annotation = feature.context$gr_true_window
  ) %>%
    left_join(df.true.lookup, by = "annotation_index") %>%
    transmute(
      direct_evidence_id = str_c(
        loop_id,
        anchor_side,
        "true_TSS",
        true_tss_id,
        sep = "|"
      ),
      loop_id,
      resolution,
      anchor_side,
      opposite_anchor_side,
      annotation_class = "true_TSS",
      annotation_chr,
      annotation_start,
      annotation_end,
      gene_id,
      gene_name,
      transcript_id_versioned,
      promoter_annotation_id = NA_character_
    )

  df.epd <- map_loop_anchors_dfr(
    gr.anchor.by.side,
    find_direct_anchor_annotation_overlaps,
    gr.annotation = feature.context$gr_epd_window
  ) %>%
    left_join(df.epd.lookup, by = "annotation_index") %>%
    transmute(
      direct_evidence_id = str_c(
        loop_id,
        anchor_side,
        "EPD_promoter",
        promoter_annotation_id,
        sep = "|"
      ),
      loop_id,
      resolution,
      anchor_side,
      opposite_anchor_side,
      annotation_class = "EPD_promoter",
      annotation_chr,
      annotation_start,
      annotation_end,
      gene_id,
      gene_name,
      transcript_id_versioned = NA_character_,
      promoter_annotation_id
    )

  bind_rows(df.true, df.epd) %>%
    distinct(direct_evidence_id, .keep_all = TRUE) %>%
    arrange(loop_id, anchor_side, gene_id, annotation_class)
}

# Apply the same primary promoter/TSS +/-1-kb and TSS-excluded ATAC >=50-bp
# definitions used by the revised resubmission analysis to an arbitrary loop set.
annotate_loop_resource <- function(
  df.loop,
  condition,
  feature.context,
  df.promoter
) {
  gr.anchor.by.side <- create_loop_anchor_granges_by_side(df.loop)
  df.direct.evidence <- build_direct_promoter_evidence(
    gr.anchor.by.side,
    feature.context,
    df.promoter
  )
  direct.summary <- summarise_direct_promoter_tss_evidence(
    df.direct.evidence,
    df.loop
  )

  # Measure promoter accessibility on each assigned +/-1-kb promoter window.
  df.promoter.window <- df.direct.evidence %>%
    distinct(
      loop_id,
      resolution,
      anchor_side,
      annotation_chr,
      annotation_start,
      annotation_end
    )
  gr.promoter.window <- GRanges(
    seqnames = df.promoter.window$annotation_chr,
    ranges = IRanges(
      start = df.promoter.window$annotation_start,
      end = df.promoter.window$annotation_end
    ),
    loop_id = df.promoter.window$loop_id,
    resolution = df.promoter.window$resolution,
    anchor_side = df.promoter.window$anchor_side
  )
  df.promoter.atac <- summarise_anchor_interval_overlap(
    gr.promoter.window,
    feature.context$gr_atac_union,
    minimum.overlap.bp = atac.minimum.overlap.bp
  ) %>%
    group_by(loop_id, resolution, anchor_side) %>%
    summarise(
      promoter_window_atac_ge50 = any(has_feature_overlap_ge_minimum),
      .groups = "drop"
    )

  # Measure residual ATAC after removing every known TSS +/-1-kb interval.
  gr.anchor.all <- combine_loop_anchor_granges_by_side(gr.anchor.by.side)
  df.non.tss.atac <- summarise_anchor_interval_overlap(
    gr.anchor.all,
    feature.context$gr_atac_non_tss,
    minimum.overlap.bp = atac.minimum.overlap.bp
  ) %>%
    transmute(
      loop_id,
      resolution,
      anchor_side,
      residual_non_tss_atac_any = has_feature_overlap_any,
      residual_non_tss_atac_ge50 = has_feature_overlap_ge_minimum,
      residual_non_tss_atac_bp = feature_overlap_bp
    )

  df.anchor.evidence <- direct.summary$anchor_summary %>%
    left_join(
      df.promoter.atac,
      by = c("loop_id", "resolution", "anchor_side")
    ) %>%
    left_join(
      df.non.tss.atac,
      by = c("loop_id", "resolution", "anchor_side")
    ) %>%
    mutate(
      promoter_window_atac_ge50 = coalesce(
        promoter_window_atac_ge50,
        FALSE
      ),
      residual_non_tss_atac_any = coalesce(
        residual_non_tss_atac_any,
        FALSE
      ),
      residual_non_tss_atac_ge50 = coalesce(
        residual_non_tss_atac_ge50,
        FALSE
      ),
      residual_non_tss_atac_bp = coalesce(residual_non_tss_atac_bp, 0L)
    )

  df.anchor.wide <- df.anchor.evidence %>%
    dplyr::select(
      loop_id,
      anchor_side,
      promoter_window_atac_ge50,
      residual_non_tss_atac_any,
      residual_non_tss_atac_ge50,
      residual_non_tss_atac_bp
    ) %>%
    pivot_wider(
      names_from = anchor_side,
      values_from = -loop_id,
      names_glue = "{.value}_{anchor_side}"
    )

  df.loop.evidence <- direct.summary$loop_summary %>%
    left_join(df.anchor.wide, by = "loop_id") %>%
    mutate(
      condition = condition,
      candidate_anchor_non_tss_atac_ge50 = case_when(
        direct_anchor_assignment_class == "direct_anchor1_only" ~
          residual_non_tss_atac_ge50_anchor2,
        direct_anchor_assignment_class == "direct_anchor2_only" ~
          residual_non_tss_atac_ge50_anchor1,
        TRUE ~ NA
      ),
      supports_anchor1_promoter_to_anchor2_regulatory = (
        n_direct_anchor_sides == 2L &
          promoter_window_atac_ge50_anchor1 &
          residual_non_tss_atac_ge50_anchor2
      ),
      supports_anchor2_promoter_to_anchor1_regulatory = (
        n_direct_anchor_sides == 2L &
          promoter_window_atac_ge50_anchor2 &
          residual_non_tss_atac_ge50_anchor1
      ),
      n_accessible_promoter_anchors = if_else(
        n_direct_anchor_sides == 2L,
        as.integer(promoter_window_atac_ge50_anchor1) +
          as.integer(promoter_window_atac_ge50_anchor2),
        NA_integer_
      ),
      n_supported_regulatory_directions = if_else(
        n_direct_anchor_sides == 2L,
        as.integer(supports_anchor1_promoter_to_anchor2_regulatory) +
          as.integer(supports_anchor2_promoter_to_anchor1_regulatory),
        NA_integer_
      ),
      dual_promoter_directional_category = case_when(
        n_direct_anchor_sides != 2L ~ "not_dual_promoter_loop",
        n_supported_regulatory_directions == 1L &
          n_accessible_promoter_anchors == 1L ~
          "one_direction_stringent_asymmetric_mixed_regulatory_evidence",
        n_supported_regulatory_directions == 1L &
          n_accessible_promoter_anchors == 2L ~
          "one_direction_broader_mixed_promoter_regulatory_evidence",
        n_supported_regulatory_directions == 2L ~
          "both_directions_supported_bidirectional_or_ambiguous",
        TRUE ~
          "neither_direction_supported_promoter_promoter_or_lower_support"
      ),
      revised_putative_regulatory_support = (
        n_direct_anchor_sides == 1L &
          coalesce(candidate_anchor_non_tss_atac_ge50, FALSE)
      ),
      revised_major_category = case_when(
        revised_putative_regulatory_support ~
          "putative_regulatory_direct_promoter_TSS_opposite_nonTSS_ATAC",
        n_direct_anchor_sides == 2L ~
          "promoter_promoter_compatible_both_direct_anchors",
        n_direct_anchor_sides == 1L ~
          "single_direct_promoter_TSS_without_opposite_nonTSS_ATAC",
        TRUE ~ "no_direct_promoter_TSS"
      )
    ) %>%
    arrange(chr1, start1, end1, chr2, start2, end2)

  df.gene.assignment <- direct.summary$gene_assignment %>%
    left_join(
      df.loop.evidence %>%
        dplyr::select(
          loop_id,
          condition,
          revised_putative_regulatory_support,
          revised_major_category
        ),
      by = "loop_id"
    ) %>%
    mutate(
      ensembl_gene_id = gene_id,
      gene_symbol = gene_name
    )

  df.gene.count <- df.gene.assignment %>%
    filter(revised_putative_regulatory_support, !is.na(ensembl_gene_id)) %>%
    distinct(loop_id, ensembl_gene_id, .keep_all = TRUE) %>%
    group_by(condition, ensembl_gene_id) %>%
    summarise(
      gene_symbol = {
        values <- sort(unique(na.omit(gene_symbol)))
        if (length(values) == 0L) NA_character_ else values[[1]]
      },
      n_exact_loop_calls = n_distinct(loop_id),
      .groups = "drop"
    ) %>%
    arrange(desc(n_exact_loop_calls), gene_symbol, ensembl_gene_id)

  list(
    loop = df.loop.evidence,
    anchor = df.anchor.evidence,
    gene_assignment = df.gene.assignment,
    gene_count = df.gene.count
  )
}

# Write a named table list without keeping obsolete duplicate files.
write_downsampling_tables <- function(tables, output.dir) {
  walk2(
    tables,
    names(tables),
    ~ readr::write_tsv(.x, file.path(output.dir, str_c(.y, ".tsv")))
  )
}

################################################################################
# 1. Detect complete samples and build matched full/downsampled loop resources
################################################################################

df.sample.status <- build_sample_file_status(
  df.sample.metadata,
  full.depth.root,
  downsample.root
)

df.completed.sample <- df.sample.status %>%
  filter(analysis_ready, provenance_ready)
df.pending.sample <- df.sample.status %>%
  filter(!analysis_ready | !provenance_ready)

if (nrow(df.completed.sample) == 0L) {
  stop("No complete 5/10/25-kb downsampled HiCCUPS result was found.", call. = FALSE)
}

message(
  "Analyzing ",
  nrow(df.completed.sample),
  " matched samples: ",
  str_c(df.completed.sample$sample, collapse = ", ")
)
if (nrow(df.pending.sample) > 0L) {
  message(
    "Pending samples excluded until complete: ",
    str_c(df.pending.sample$sample, collapse = ", ")
  )
}

# Confirm that every included QC record matches the sample, seed, source depth,
# and 140M target encoded in the analysis metadata.
df.downsampling.qc <- map_dfr(seq_len(nrow(df.completed.sample)), function(i) {
  row <- df.completed.sample[i, ]
  qc <- readr::read_tsv(
    row$downsampling_qc_file,
    col_types = cols(
      sample = col_character(), seed = col_double(),
      source_contacts = col_double(), downsampled_contacts = col_double(),
      output = col_character()
    ),
    show_col_types = FALSE
  )
  if (nrow(qc) != 1L) {
    stop("Expected one QC row in ", row$downsampling_qc_file, call. = FALSE)
  }
  qc %>%
    transmute(
      sample,
      seed = as.integer(seed),
      source_contacts,
      target_contacts = downsampled_contacts,
      source_output_on_hpc = output,
      qc_matches_metadata = (
        sample == row$sample & seed == row$seed &
          source_contacts == row$source_contacts &
          downsampled_contacts == target.contacts
      )
    )
})
if (any(!df.downsampling.qc$qc_matches_metadata)) {
  stop("A 140M downsampling QC row does not match the metadata.", call. = FALSE)
}

full.file.metadata <- df.completed.sample %>%
  transmute(sample, strain, file = full_file)
downsample.file.metadata <- df.completed.sample %>%
  transmute(sample, strain, file = downsample_file)

full.resource <- build_depth_loop_resource(
  full.file.metadata,
  "full_depth"
)
downsample.resource <- build_depth_loop_resource(
  downsample.file.metadata,
  "downsample_140M"
)

################################################################################
# 2. Sample loop-count dispersion before and after depth normalization
################################################################################

df.sample.loop.count <- bind_rows(
  full.resource$sample,
  downsample.resource$sample
) %>%
  filter(passes_lt2mb) %>%
  count(condition, sample, resolution, name = "n_loop_calls")

# Sum the three resolution-specific rows to one ALL row per sample and depth.
df.sample.loop.count <- bind_rows(
  df.sample.loop.count,
  df.sample.loop.count %>%
    group_by(condition, sample) %>%
    summarise(
      resolution = "ALL",
      n_loop_calls = sum(n_loop_calls),
      .groups = "drop"
    )
) %>%
  left_join(
    df.completed.sample %>%
      dplyr::select(
        sample,
        strain,
        source_contacts,
        target_contacts,
        sampling_fraction
      ),
    by = "sample"
  ) %>%
  arrange(condition, resolution, sample)

df.sample.loop.count.paired <- df.sample.loop.count %>%
  dplyr::select(sample, strain, resolution, condition, n_loop_calls) %>%
  pivot_wider(names_from = condition, values_from = n_loop_calls) %>%
  mutate(
    absolute_change = downsample_140M - full_depth,
    pct_change = 100 * absolute_change / full_depth
  ) %>%
  left_join(
    df.completed.sample %>%
      dplyr::select(sample, source_contacts, target_contacts, sampling_fraction),
    by = "sample"
  ) %>%
  arrange(resolution, source_contacts, sample)

df.loop.count.dispersion <- df.sample.loop.count %>%
  group_by(condition, resolution) %>%
  summarise(
    n_samples = n(),
    mean_loop_calls = mean(n_loop_calls),
    sd_loop_calls = sd(n_loop_calls),
    cv_loop_calls = coefficient_of_variation(n_loop_calls),
    min_loop_calls = min(n_loop_calls),
    max_loop_calls = max(n_loop_calls),
    range_loop_calls = max_loop_calls - min_loop_calls,
    .groups = "drop"
  ) %>%
  arrange(resolution, condition)

df.depth.loop.correlation <- df.sample.loop.count %>%
  group_by(condition, resolution) %>%
  summarise(
    n_samples = n(),
    spearman_source_contacts_vs_loop_calls = safe_spearman(
      source_contacts,
      n_loop_calls
    ),
    spearman_p_value = safe_spearman_test(
      source_contacts,
      n_loop_calls
    )$spearman_p_value,
    .groups = "drop"
  ) %>%
  arrange(resolution, condition)

################################################################################
# 3. Exact loop recovery and resolution-specific stability
################################################################################

df.sample.loop.membership <- bind_rows(
  full.resource$sample,
  downsample.resource$sample
) %>%
  filter(passes_lt2mb) %>%
  distinct(condition, sample, resolution, loop_id)

df.exact.recovery.by.sample.resolution <- df.sample.loop.membership %>%
  add_all_resolution() %>%
  group_by(sample, resolution) %>%
  summarise(
    n_full_depth = n_distinct(loop_id[condition == "full_depth"]),
    n_downsample_140M = n_distinct(loop_id[condition == "downsample_140M"]),
    n_exact_shared = length(intersect(
      loop_id[condition == "full_depth"],
      loop_id[condition == "downsample_140M"]
    )),
    n_full_depth_not_recovered = n_full_depth - n_exact_shared,
    n_downsample_only = n_downsample_140M - n_exact_shared,
    pct_full_depth_exact_recovered = 100 * n_exact_shared / n_full_depth,
    exact_jaccard = n_exact_shared /
      (n_full_depth + n_downsample_140M - n_exact_shared),
    .groups = "drop"
  ) %>%
  left_join(
    df.completed.sample %>%
      dplyr::select(sample, strain, source_contacts, sampling_fraction),
    by = "sample"
  ) %>%
  arrange(resolution, source_contacts, sample)

df.exact.recovery.summary <- df.exact.recovery.by.sample.resolution %>%
  group_by(resolution) %>%
  summarise(
    n_samples = n(),
    median_pct_full_depth_exact_recovered = median(
      pct_full_depth_exact_recovered
    ),
    min_pct_full_depth_exact_recovered = min(
      pct_full_depth_exact_recovered
    ),
    max_pct_full_depth_exact_recovered = max(
      pct_full_depth_exact_recovered
    ),
    median_exact_jaccard = median(exact_jaccard),
    .groups = "drop"
  ) %>%
  arrange(resolution)

# Add a location-tolerant sensitivity based on the configured HiCCUPS radius;
# exact calls remain the primary resource and are never collapsed here.
df.hiccups.radius.recovery <- summarise_hiccups_radius_recovery(
  full.resource$sample,
  downsample.resource$sample
) %>%
  left_join(
    df.completed.sample %>%
      dplyr::select(sample, strain, source_contacts, sampling_fraction),
    by = "sample"
  ) %>%
  arrange(resolution, source_contacts, sample)

df.hiccups.radius.recovery.summary <- df.hiccups.radius.recovery %>%
  group_by(resolution) %>%
  summarise(
    n_samples = n(),
    median_pct_full_depth_radius_recovered = median(
      pct_full_depth_radius_recovered
    ),
    min_pct_full_depth_radius_recovered = min(
      pct_full_depth_radius_recovered
    ),
    max_pct_full_depth_radius_recovered = max(
      pct_full_depth_radius_recovered
    ),
    .groups = "drop"
  ) %>%
  arrange(resolution)

df.pooled.loop.membership <- bind_rows(
  full.resource$universe %>% dplyr::select(condition, resolution, loop_id),
  downsample.resource$universe %>% dplyr::select(condition, resolution, loop_id)
) %>%
  distinct()

df.pooled.exact.recovery <- df.pooled.loop.membership %>%
  add_all_resolution() %>%
  group_by(resolution) %>%
  summarise(
    n_full_depth_pooled = n_distinct(loop_id[condition == "full_depth"]),
    n_downsample_140M_pooled = n_distinct(
      loop_id[condition == "downsample_140M"]
    ),
    n_exact_shared_pooled = length(intersect(
      loop_id[condition == "full_depth"],
      loop_id[condition == "downsample_140M"]
    )),
    pct_full_depth_pooled_exact_recovered =
      100 * n_exact_shared_pooled / n_full_depth_pooled,
    pooled_exact_jaccard = n_exact_shared_pooled /
      (n_full_depth_pooled + n_downsample_140M_pooled -
        n_exact_shared_pooled),
    .groups = "drop"
  ) %>%
  arrange(resolution)

################################################################################
# 4. Reapply promoter/TSS and ATAC categories to both pooled resources
################################################################################

required.feature.objects <- c(
  "df.transcript.ensembl.rn7.1based",
  "df.promoter.epd.rn7.1based",
  "gr.atac.rn7.1based"
)
feature.paths <- coordinate_cache_paths(
  coord.cache.dir,
  required.feature.objects
)
missing.feature.paths <- feature.paths[!file.exists(feature.paths)]
if (length(missing.feature.paths) > 0L) {
  stop(
    "Missing coordinate-cache feature files:\n",
    str_c(missing.feature.paths, collapse = "\n"),
    call. = FALSE
  )
}

df.transcript.ensembl.rn7.1based <- readRDS(
  feature.paths[["df.transcript.ensembl.rn7.1based"]]
)
df.promoter.epd.rn7.1based <- readRDS(
  feature.paths[["df.promoter.epd.rn7.1based"]]
)
gr.atac.rn7.1based <- readRDS(feature.paths[["gr.atac.rn7.1based"]])

feature.context <- build_regulatory_feature_context(
  df.transcript.ensembl.rn7.1based,
  df.promoter.epd.rn7.1based,
  gr.atac.rn7.1based,
  promoter.flank.bp = promoter.window.flank.bp
)

message("Annotating matched full-depth pooled calls.")
full.annotation <- annotate_loop_resource(
  full.resource$universe,
  "full_depth",
  feature.context,
  df.promoter.epd.rn7.1based
)
message("Annotating 140M downsampled pooled calls.")
downsample.annotation <- annotate_loop_resource(
  downsample.resource$universe,
  "downsample_140M",
  feature.context,
  df.promoter.epd.rn7.1based
)

df.category.count <- bind_rows(
  full.annotation$loop,
  downsample.annotation$loop
) %>%
  add_all_resolution() %>%
  count(
    condition,
    resolution,
    revised_major_category,
    name = "n_loop_calls"
  ) %>%
  group_by(condition, resolution) %>%
  mutate(pct_within_condition_resolution = 100 * n_loop_calls / sum(n_loop_calls)) %>%
  ungroup() %>%
  arrange(resolution, revised_major_category, condition)

df.category.composition.change <- df.category.count %>%
  dplyr::select(
    condition,
    resolution,
    revised_major_category,
    n_loop_calls,
    pct_within_condition_resolution
  ) %>%
  pivot_wider(
    names_from = condition,
    values_from = c(n_loop_calls, pct_within_condition_resolution),
    values_fill = 0
  ) %>%
  mutate(
    absolute_loop_change = n_loop_calls_downsample_140M -
      n_loop_calls_full_depth,
    percentage_point_change = pct_within_condition_resolution_downsample_140M -
      pct_within_condition_resolution_full_depth
  ) %>%
  arrange(resolution, revised_major_category)

df.category.exact.recovery <- full.annotation$loop %>%
  dplyr::select(loop_id, resolution, revised_major_category) %>%
  mutate(exact_recovered_at_140M = loop_id %in% downsample.annotation$loop$loop_id) %>%
  add_all_resolution() %>%
  group_by(resolution, revised_major_category) %>%
  summarise(
    n_full_depth_loop_calls = n(),
    n_exact_recovered_at_140M = sum(exact_recovered_at_140M),
    pct_exact_recovered_at_140M = 100 * mean(exact_recovered_at_140M),
    .groups = "drop"
  ) %>%
  arrange(resolution, revised_major_category)

# Exact shared coordinates should retain the same coordinate-derived category;
# report any disagreement rather than silently assuming perfect concordance.
df.shared.loop.category.concordance <- full.annotation$loop %>%
  dplyr::select(
    loop_id,
    resolution,
    full_depth_category = revised_major_category
  ) %>%
  inner_join(
    downsample.annotation$loop %>%
      dplyr::select(
        loop_id,
        downsample_140M_category = revised_major_category
      ),
    by = "loop_id"
  ) %>%
  count(
    resolution,
    full_depth_category,
    downsample_140M_category,
    name = "n_exact_shared_loop_calls"
  ) %>%
  arrange(resolution, full_depth_category, downsample_140M_category)

################################################################################
# 5. Putative-regulatory gene-loop count stability
################################################################################

df.gene.count.stability <- full.annotation$gene_count %>%
  transmute(
    ensembl_gene_id,
    gene_symbol_full_depth = gene_symbol,
    n_exact_loop_calls_full_depth = n_exact_loop_calls
  ) %>%
  full_join(
    downsample.annotation$gene_count %>%
      transmute(
        ensembl_gene_id,
        gene_symbol_downsample_140M = gene_symbol,
        n_exact_loop_calls_downsample_140M = n_exact_loop_calls
      ),
    by = "ensembl_gene_id"
  ) %>%
  mutate(
    gene_symbol = coalesce(
      gene_symbol_full_depth,
      gene_symbol_downsample_140M
    ),
    n_exact_loop_calls_full_depth = coalesce(
      n_exact_loop_calls_full_depth,
      0L
    ),
    n_exact_loop_calls_downsample_140M = coalesce(
      n_exact_loop_calls_downsample_140M,
      0L
    ),
    exact_loop_call_change = n_exact_loop_calls_downsample_140M -
      n_exact_loop_calls_full_depth,
    present_full_depth = n_exact_loop_calls_full_depth > 0L,
    present_downsample_140M = n_exact_loop_calls_downsample_140M > 0L,
    rank_full_depth = rank(
      -n_exact_loop_calls_full_depth,
      ties.method = "average"
    ),
    rank_downsample_140M = rank(
      -n_exact_loop_calls_downsample_140M,
      ties.method = "average"
    )
  ) %>%
  dplyr::select(
    ensembl_gene_id,
    gene_symbol,
    n_exact_loop_calls_full_depth,
    n_exact_loop_calls_downsample_140M,
    exact_loop_call_change,
    present_full_depth,
    present_downsample_140M,
    rank_full_depth,
    rank_downsample_140M
  ) %>%
  arrange(
    rank_full_depth,
    rank_downsample_140M,
    gene_symbol,
    ensembl_gene_id
  )

gene.intersection <- sum(
  df.gene.count.stability$present_full_depth &
    df.gene.count.stability$present_downsample_140M
)
gene.union <- nrow(df.gene.count.stability)

df.gene.count.stability.summary <- tibble(
  n_genes_full_depth = sum(df.gene.count.stability$present_full_depth),
  n_genes_downsample_140M = sum(
    df.gene.count.stability$present_downsample_140M
  ),
  n_genes_shared = gene.intersection,
  gene_presence_jaccard = gene.intersection / gene.union,
  spearman_loop_counts_union_zero_filled = safe_spearman(
    df.gene.count.stability$n_exact_loop_calls_full_depth,
    df.gene.count.stability$n_exact_loop_calls_downsample_140M
  ),
  spearman_gene_ranks_union = safe_spearman(
    df.gene.count.stability$rank_full_depth,
    df.gene.count.stability$rank_downsample_140M
  )
)

################################################################################
# 6. Pooled-resource dependence on source-library sequencing depth
################################################################################

df.pool.sample.support <- bind_rows(
  full.resource$sample,
  downsample.resource$sample
) %>%
  filter(passes_lt2mb) %>%
  distinct(condition, loop_id, sample, resolution) %>%
  group_by(condition, loop_id) %>%
  mutate(n_supporting_samples = n_distinct(sample)) %>%
  ungroup()

df.pool.size <- df.pool.sample.support %>%
  distinct(condition, loop_id) %>%
  count(condition, name = "n_pooled_exact_loop_calls")

df.sample.pool.influence <- df.pool.sample.support %>%
  group_by(condition, sample) %>%
  summarise(
    n_pooled_calls_supported = n_distinct(loop_id),
    n_pooled_calls_unique_to_sample = n_distinct(
      loop_id[n_supporting_samples == 1L]
    ),
    .groups = "drop"
  ) %>%
  left_join(df.pool.size, by = "condition") %>%
  mutate(
    pct_pooled_calls_supported = 100 * n_pooled_calls_supported /
      n_pooled_exact_loop_calls,
    pct_pooled_calls_unique_to_sample = 100 * n_pooled_calls_unique_to_sample /
      n_pooled_exact_loop_calls
  ) %>%
  left_join(
    df.completed.sample %>%
      dplyr::select(sample, strain, source_contacts, sampling_fraction),
    by = "sample"
  ) %>%
  arrange(condition, source_contacts, sample)

df.pool.influence.correlation <- df.sample.pool.influence %>%
  group_by(condition) %>%
  summarise(
    n_samples = n(),
    spearman_source_contacts_vs_supported_pooled_calls = safe_spearman(
      source_contacts,
      n_pooled_calls_supported
    ),
    spearman_source_contacts_vs_unique_pooled_calls = safe_spearman(
      source_contacts,
      n_pooled_calls_unique_to_sample
    ),
    spearman_supported_p_value = safe_spearman_test(
      source_contacts,
      n_pooled_calls_supported
    )$spearman_p_value,
    spearman_unique_p_value = safe_spearman_test(
      source_contacts,
      n_pooled_calls_unique_to_sample
    )$spearman_p_value,
    .groups = "drop"
  )

df.pool.support.distribution <- df.pool.sample.support %>%
  distinct(condition, loop_id, n_supporting_samples) %>%
  count(condition, n_supporting_samples, name = "n_pooled_exact_loop_calls") %>%
  group_by(condition) %>%
  mutate(pct_pooled_exact_loop_calls = 100 * n_pooled_exact_loop_calls /
    sum(n_pooled_exact_loop_calls)) %>%
  ungroup() %>%
  arrange(condition, n_supporting_samples)

median.source.contacts <- median(df.completed.sample$source_contacts)
df.full.pool.depth.support <- full.resource$sample %>%
  filter(passes_lt2mb) %>%
  distinct(loop_id, sample) %>%
  left_join(
    df.completed.sample %>% dplyr::select(sample, source_contacts),
    by = "sample"
  ) %>%
  group_by(loop_id) %>%
  summarise(
    n_supporting_samples_full_depth = n_distinct(sample),
    mean_supporting_source_contacts = mean(source_contacts),
    max_supporting_source_contacts = max(source_contacts),
    supported_only_by_above_median_depth_libraries = all(
      source_contacts > median.source.contacts
    ),
    .groups = "drop"
  ) %>%
  mutate(
    exact_recovered_in_downsampled_pool = loop_id %in%
      downsample.resource$universe$loop_id
  )

df.full.pool.depth.recovery <- df.full.pool.depth.support %>%
  group_by(supported_only_by_above_median_depth_libraries) %>%
  summarise(
    n_full_depth_pooled_calls = n(),
    n_exact_recovered_in_downsampled_pool = sum(
      exact_recovered_in_downsampled_pool
    ),
    pct_exact_recovered_in_downsampled_pool = 100 * mean(
      exact_recovered_in_downsampled_pool
    ),
    median_supporting_source_contacts = median(
      mean_supporting_source_contacts
    ),
    .groups = "drop"
  )

################################################################################
# 7. Save compact tables, figures, and reproducibility metadata
################################################################################

df.analysis.caveat <- tribble(
  ~analysis_issue, ~interpretive_requirement,
  "one_library_per_strain", paste0(
    "Equalizing valid contacts does not create biological replication or ",
    "separate strain effects from specimen- and library-specific effects."
  ),
  "single_random_seed", paste0(
    "One downsampling realization was analyzed, so Monte Carlo variability ",
    "from contact subsampling was not estimated."
  ),
  "valid_contacts_not_all_library_properties", paste0(
    "The analysis controls the number of duplicate-removed MAPQ >=30 valid ",
    "contacts, but not contact-distance composition or other library-quality ",
    "differences."
  ),
  "five_kb_ignore_sparsity", paste0(
    "HiCCUPS used --ignore-sparsity at 5/10/25 kb. The particularly low 5-kb ",
    "call yield at 140M must be reported by resolution and interpreted as a ",
    "depth-sensitivity result rather than evidence that 5-kb calls are absent."
  )
)

output.tables <- list(
  downsampling_sample_status = df.sample.status,
  downsampling_qc = df.downsampling.qc,
  downsampling_sample_loop_counts = df.sample.loop.count,
  downsampling_sample_loop_count_paired = df.sample.loop.count.paired,
  downsampling_loop_count_dispersion = df.loop.count.dispersion,
  downsampling_depth_loop_correlations = df.depth.loop.correlation,
  downsampling_exact_recovery_by_sample_resolution =
    df.exact.recovery.by.sample.resolution,
  downsampling_exact_recovery_summary = df.exact.recovery.summary,
  downsampling_hiccups_radius_recovery = df.hiccups.radius.recovery,
  downsampling_hiccups_radius_recovery_summary =
    df.hiccups.radius.recovery.summary,
  downsampling_pooled_exact_recovery = df.pooled.exact.recovery,
  downsampling_category_counts = df.category.count,
  downsampling_category_composition_change = df.category.composition.change,
  downsampling_category_exact_recovery = df.category.exact.recovery,
  downsampling_shared_loop_category_concordance =
    df.shared.loop.category.concordance,
  downsampling_gene_count_stability = df.gene.count.stability,
  downsampling_gene_count_stability_summary =
    df.gene.count.stability.summary,
  downsampling_sample_pool_influence = df.sample.pool.influence,
  downsampling_pool_influence_correlations = df.pool.influence.correlation,
  downsampling_pool_support_distribution = df.pool.support.distribution,
  downsampling_full_pool_depth_recovery = df.full.pool.depth.recovery,
  downsampling_analysis_caveats = df.analysis.caveat
)
write_downsampling_tables(output.tables, output.dir)

df.method.definition <- tribble(
  ~analysis_component, ~definition,
  "matched_sample_design", paste0(
    "Full-depth and 140M analyses use only samples with complete 5/10/25-kb ",
    "downsampled HiCCUPS outputs. Pending samples are excluded from both sides."
  ),
  "loop_universe", paste0(
    "Exact resolution-specific HiCCUPS calls with loop distance <2 Mb; nearby ",
    "or cross-resolution calls are not collapsed in this analysis."
  ),
  "exact_recovery", paste0(
    "A full-depth loop is recovered only when sample, anchor coordinates, and ",
    "resolution exactly match a 140M call."
  ),
  "hiccups_radius_recovery_sensitivity", paste0(
    "A sensitivity-only recovery allows both HiCCUPS centroids to differ by ",
    "no more than the configured 5/10/25-kb radii of 20/20/50 kb within the ",
    "same sample, resolution, and chromosome pair; it does not merge calls."
  ),
  "promoter_TSS", paste0(
    "Direct anchor overlap with an Ensembl or EPD TSS +/-1-kb promoter window; ",
    "all overlapping genes and transcripts are retained."
  ),
  "ATAC_support", paste0(
    "Duttke 2022 snATAC intervals are reduced, all known Ensembl/EPD TSS +/-1 ",
    "kb regions are removed, and >=50 bp residual overlap is required."
  ),
  "gene_stability", paste0(
    "Counts are distinct exact putative-regulatory loop calls per Ensembl gene; ",
    "missing genes are assigned zero before rank/count correlation."
  ),
  "pooled_depth_influence", paste0(
    "Per-sample pooled-call support, sample-unique calls, support-count ",
    "distribution, and correlations with original contact depth are reported."
  ),
  "interpretation_limit", paste0(
    "The all-ten 140M series is a depth-sensitivity analysis with one seed and ",
    "one library per strain. Correlations are descriptive and do not establish ",
    "strain-specific chromatin architecture."
  )
)
readr::write_tsv(
  df.method.definition,
  file.path(output.dir, "downsampling_method_definitions.tsv")
)

plot.loop.count <- ggplot(
  df.sample.loop.count.paired,
  aes(x = full_depth, y = downsample_140M, color = resolution, label = sample)
) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey55") +
  geom_point(size = 2.5) +
  geom_text(nudge_y = 80, size = 3, check_overlap = TRUE) +
  facet_wrap(~ resolution, scales = "free") +
  labs(
    x = "Full-depth exact loop calls (<2 Mb)",
    y = "140M exact loop calls (<2 Mb)",
    color = "Resolution"
  ) +
  theme_bw()

ggsave(
  file.path(output.dir, "downsampling_loop_counts_full_vs_140M.png"),
  plot.loop.count,
  width = 10,
  height = 7,
  dpi = 300
)
ggsave(
  file.path(output.dir, "downsampling_loop_counts_full_vs_140M.pdf"),
  plot.loop.count,
  width = 10,
  height = 7
)

plot.recovery <- ggplot(
  df.exact.recovery.by.sample.resolution,
  aes(
    x = reorder(sample, source_contacts),
    y = pct_full_depth_exact_recovered,
    fill = resolution
  )
) +
  geom_col(position = position_dodge(width = 0.8)) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(
    x = "Sample ordered by original contact depth",
    y = "Full-depth exact calls recovered at 140M (%)",
    fill = "Resolution"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(output.dir, "downsampling_exact_recovery_by_sample.png"),
  plot.recovery,
  width = 10,
  height = 6,
  dpi = 300
)
ggsave(
  file.path(output.dir, "downsampling_exact_recovery_by_sample.pdf"),
  plot.recovery,
  width = 10,
  height = 6
)

run.metadata <- tibble(
  run_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  n_completed_samples = nrow(df.completed.sample),
  completed_samples = str_c(df.completed.sample$sample, collapse = ";"),
  pending_samples = if_else(
    nrow(df.pending.sample) == 0L,
    "none",
    str_c(df.pending.sample$sample, collapse = ";")
  ),
  full_depth_root = full.depth.root,
  downsample_root = downsample.root,
  coordinate_cache = coord.cache.dir,
  target_contacts = target.contacts,
  maximum_loop_distance_bp = max.loop.distance.bp,
  promoter_window_flank_bp = promoter.window.flank.bp,
  atac_minimum_overlap_bp = atac.minimum.overlap.bp
)
readr::write_tsv(
  run.metadata,
  file.path(output.dir, "downsampling_run_metadata.tsv")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output.dir, "downsampling_session_info.txt")
)

message("Downsampling sensitivity analysis complete.")
message("Results written to: ", output.dir)
print(df.loop.count.dispersion)
print(df.exact.recovery.summary)
print(df.gene.count.stability.summary)

################################################################################
# 8. Archived 210M prototype, superseded by the matched 250M analysis
#
# This block is retained only as analysis history and is never executed. The
# production depth series is implemented in downsampling_140M_250M_comparison.R
# with the same seven libraries at full depth, 140M, and 250M.
# The 140M analysis above remains the all-library depth-normalization analysis.
# This second analysis excludes 74AA and A2DB because their full-depth usable-
# contact counts are below 210M. Full depth and 140M are rebuilt with the same
# remaining eight libraries before any full/140M/210M comparison is made.
#
# Requested comparisons:
# 1. Per-library loop counts and coefficients of variation.
# 2. Original usable-contact depth versus loop-count Spearman correlations.
# 3. Exact and HiCCUPS-radius/tolerance-based recovery.
# 4. Exact and approximate-locus pairwise Jaccard/overlap coefficients.
# 5. Exact-call and approximate-locus library-support-count distributions.
# 6. Recurrent/shared cores reproduced at both 140M and 210M.
################################################################################

if (FALSE) {

# Find an optional input root without interrupting the completed 140M analysis.
resolve_optional_downsample_root <- function(candidates) {
  candidates <- unique(path.expand(candidates[nzchar(candidates)]))
  candidates <- candidates[dir.exists(candidates)]
  if (length(candidates) == 0L) return(NA_character_)

  n.merged <- vapply(
    candidates,
    function(path) {
      files <- list.files(
        path,
        pattern = "merged_loops[.]bedpe$",
        recursive = TRUE,
        full.names = TRUE
      )
      sum(file.exists(files) & file.info(files)$size > 0L)
    },
    integer(1L)
  )
  normalizePath(candidates[[which.max(n.merged)]], winslash = "/")
}

# Resolve the target-specific HiCCUPS files using the same naming convention as
# the 140M Google Drive/Dropbox results.
build_target_file_status <- function(
  df.metadata,
  full.root,
  down.root,
  target.label
) {
  df.metadata %>%
    mutate(
      full_file = file.path(
        full.root, sample, "hiccups_5k10k25k", "merged_loops.bedpe"
      ),
      downsample_result_dir = file.path(
        down.root,
        sample,
        str_c(
          sample, "_inter_30_", target.label, "_seed", seed,
          ".hiccups.5k10k25k"
        )
      ),
      downsample_file = file.path(
        downsample_result_dir, "merged_loops.bedpe"
      ),
      downsample_postprocessed_5k = file.path(
        downsample_result_dir, "postprocessed_pixels_5000.bedpe"
      ),
      downsample_postprocessed_10k = file.path(
        downsample_result_dir, "postprocessed_pixels_10000.bedpe"
      ),
      downsample_postprocessed_25k = file.path(
        downsample_result_dir, "postprocessed_pixels_25000.bedpe"
      ),
      full_file_ready = file.exists(full_file) &
        file.info(full_file)$size > 0L,
      downsample_merged_ready = file.exists(downsample_file) &
        file.info(downsample_file)$size > 0L,
      downsample_5k_ready = file.exists(downsample_postprocessed_5k) &
        file.info(downsample_postprocessed_5k)$size > 0L,
      downsample_10k_ready = file.exists(downsample_postprocessed_10k) &
        file.info(downsample_postprocessed_10k)$size > 0L,
      downsample_25k_ready = file.exists(downsample_postprocessed_25k) &
        file.info(downsample_postprocessed_25k)$size > 0L,
      analysis_ready = full_file_ready & downsample_merged_ready &
        downsample_5k_ready & downsample_10k_ready &
        downsample_25k_ready,
      status = if_else(
        analysis_ready,
        "included_complete_5k10k25k",
        "pending_complete_5k10k25k"
      )
    )
}

# Calculate exact recovery for each ordered depth comparison.
summarise_depth_pair_exact_recovery <- function(
  df.membership,
  reference.condition,
  query.condition
) {
  df.membership %>%
    filter(condition %in% c(reference.condition, query.condition)) %>%
    add_all_resolution() %>%
    group_by(sample, resolution) %>%
    summarise(
      reference_condition = reference.condition,
      query_condition = query.condition,
      n_reference = n_distinct(
        loop_id[condition == reference.condition]
      ),
      n_query = n_distinct(loop_id[condition == query.condition]),
      n_exact_shared = length(intersect(
        loop_id[condition == reference.condition],
        loop_id[condition == query.condition]
      )),
      pct_reference_recovered = 100 * n_exact_shared / n_reference,
      exact_jaccard = n_exact_shared /
        (n_reference + n_query - n_exact_shared),
      exact_overlap_coefficient = n_exact_shared / pmin(
        n_reference, n_query
      ),
      .groups = "drop"
    )
}

# Reuse the HiCCUPS 20/20/50-kb radius sensitivity for any depth pair.
summarise_depth_pair_radius_recovery <- function(
  reference.resource,
  query.resource,
  reference.condition,
  query.condition
) {
  summarise_hiccups_radius_recovery(
    reference.resource$sample,
    query.resource$sample
  ) %>%
    transmute(
      sample,
      resolution,
      reference_condition = reference.condition,
      query_condition = query.condition,
      hiccups_match_radius_bp,
      n_reference = n_full_depth,
      n_query = n_downsample_140M,
      n_reference_radius_recovered = n_full_depth_radius_recovered,
      n_query_radius_matched = n_downsample_radius_matched,
      pct_reference_radius_recovered = pct_full_depth_radius_recovered,
      pct_query_radius_matched = pct_downsample_radius_matched
    )
}

# Compute unordered pairwise sharing for exact calls or approximate loci.
summarise_pairwise_library_similarity <- function(
  df.membership,
  feature.column,
  feature.type
) {
  comparison.strata <- df.membership %>%
    distinct(condition, resolution) %>%
    arrange(condition, resolution)

  pmap_dfr(comparison.strata, function(condition, resolution) {
    df.i <- df.membership %>%
      filter(
        .data$condition == .env$condition,
        .data$resolution == .env$resolution
      )
    samples <- sort(unique(df.i$sample))
    pairs <- combn(samples, 2L, simplify = FALSE)

    map_dfr(pairs, function(pair) {
      set.1 <- unique(df.i[[feature.column]][df.i$sample == pair[[1]]])
      set.2 <- unique(df.i[[feature.column]][df.i$sample == pair[[2]]])
      n.shared <- length(intersect(set.1, set.2))
      n.union <- length(union(set.1, set.2))
      n.minimum <- min(length(set.1), length(set.2))

      tibble(
        condition = condition,
        resolution = resolution,
        feature_type = feature.type,
        sample1 = pair[[1]],
        sample2 = pair[[2]],
        n_sample1 = length(set.1),
        n_sample2 = length(set.2),
        n_shared = n.shared,
        n_union = n.union,
        jaccard_similarity = if_else(
          n.union > 0L, n.shared / n.union, NA_real_
        ),
        overlap_coefficient = if_else(
          n.minimum > 0L, n.shared / n.minimum, NA_real_
        )
      )
    })
  })
}

# Summarize how many libraries support each exact call or approximate locus.
summarise_library_support_distribution <- function(
  df.membership,
  feature.column,
  feature.type
) {
  df.membership %>%
    distinct(condition, sample, .data[[feature.column]]) %>%
    count(
      condition,
      .data[[feature.column]],
      name = "n_supporting_libraries"
    ) %>%
    count(
      condition,
      n_supporting_libraries,
      name = "n_features"
    ) %>%
    group_by(condition) %>%
    mutate(
      feature_type = feature.type,
      pct_features = 100 * n_features / sum(n_features)
    ) %>%
    ungroup() %>%
    dplyr::select(
      condition, feature_type, n_supporting_libraries,
      n_features, pct_features
    )
}

# Test whether recurrent features (support >=k libraries) are reproduced at
# both 140M and 210M; k=2 is the primary recurrent/shared-core definition.
summarise_recurrent_core_stability <- function(
  df.membership,
  feature.column,
  feature.type,
  maximum.support = 8L
) {
  df.support <- df.membership %>%
    filter(condition %in% c("downsample_140M", "downsample_210M")) %>%
    distinct(condition, sample, .data[[feature.column]]) %>%
    count(
      condition,
      .data[[feature.column]],
      name = "n_supporting_libraries"
    )

  map_dfr(seq_len(maximum.support), function(minimum.support) {
    set.140 <- df.support %>%
      filter(
        condition == "downsample_140M",
        n_supporting_libraries >= minimum.support
      ) %>%
      pull(.data[[feature.column]]) %>%
      unique()
    set.210 <- df.support %>%
      filter(
        condition == "downsample_210M",
        n_supporting_libraries >= minimum.support
      ) %>%
      pull(.data[[feature.column]]) %>%
      unique()
    n.shared <- length(intersect(set.140, set.210))
    n.union <- length(union(set.140, set.210))
    n.minimum <- min(length(set.140), length(set.210))

    tibble(
      feature_type = feature.type,
      minimum_library_support = minimum.support,
      n_140M = length(set.140),
      n_210M = length(set.210),
      n_shared_140M_210M = n.shared,
      jaccard_140M_210M = if_else(
        n.union > 0L, n.shared / n.union, NA_real_
      ),
      overlap_coefficient_140M_210M = if_else(
        n.minimum > 0L, n.shared / n.minimum, NA_real_
      )
    )
  })
}

# Run the complete 8-library depth series after all 210M HiCCUPS outputs exist.
run_210m_depth_series_analysis <- function(
  df.metadata.210m,
  downsample.210m.root,
  depth.series.output.dir
) {
  df.status.210m <- build_target_file_status(
    df.metadata.210m,
    full.depth.root,
    downsample.210m.root,
    "210M"
  )
  readr::write_tsv(
    df.status.210m,
    file.path(depth.series.output.dir, "depth_series_210M_sample_status.tsv")
  )
  if (!all(df.status.210m$analysis_ready)) {
    pending <- df.status.210m$sample[!df.status.210m$analysis_ready]
    message(
      "210M analysis pending complete HiCCUPS results for: ",
      str_c(pending, collapse = ", ")
    )
    return(invisible(NULL))
  }

  # Rebuild all three resources with exactly the same eight libraries.
  df.status.140m.8 <- df.sample.status %>%
    filter(sample %in% df.metadata.210m$sample)
  stopifnot(all(df.status.140m.8$analysis_ready))

  resource.list <- list(
    full_depth = build_depth_loop_resource(
      df.status.210m %>% transmute(sample, strain, file = full_file),
      "full_depth"
    ),
    downsample_140M = build_depth_loop_resource(
      df.status.140m.8 %>% transmute(sample, strain, file = downsample_file),
      "downsample_140M"
    ),
    downsample_210M = build_depth_loop_resource(
      df.status.210m %>% transmute(sample, strain, file = downsample_file),
      "downsample_210M"
    )
  )

  df.sample.membership <- imap_dfr(
    resource.list,
    ~ .x$sample %>%
      filter(passes_lt2mb) %>%
      transmute(
        condition = .y, sample, strain, resolution, loop_id,
        chr1, start1, end1, chr2, start2, end2
      )
  ) %>%
    distinct()

  # Per-library loop counts, CVs, and residual source-depth correlations.
  df.loop.count.by.resolution <- df.sample.membership %>%
    count(condition, sample, strain, resolution, name = "n_loop_calls") %>%
    arrange(condition, sample, resolution)
  df.loop.count <- bind_rows(
    df.loop.count.by.resolution,
    df.loop.count.by.resolution %>%
      group_by(condition, sample, strain) %>%
      summarise(
        resolution = "ALL",
        n_loop_calls = sum(n_loop_calls),
        .groups = "drop"
      )
  ) %>%
    left_join(
      df.metadata.210m %>%
        dplyr::select(sample, source_contacts),
      by = "sample"
    ) %>%
    mutate(
      analysis_contacts = case_when(
        condition == "full_depth" ~ source_contacts,
        condition == "downsample_140M" ~ 140000000,
        condition == "downsample_210M" ~ 210000000,
        TRUE ~ NA_real_
      )
    ) %>%
    arrange(resolution, condition, source_contacts)

  df.loop.count.summary <- df.loop.count %>%
    group_by(condition, resolution) %>%
    summarise(
      n_libraries = n(),
      mean_loop_calls = mean(n_loop_calls),
      sd_loop_calls = sd(n_loop_calls),
      cv_loop_calls = coefficient_of_variation(n_loop_calls),
      min_loop_calls = min(n_loop_calls),
      max_loop_calls = max(n_loop_calls),
      spearman_original_contacts_vs_loop_calls = safe_spearman(
        source_contacts, n_loop_calls
      ),
      .groups = "drop"
    )

  # Exact recovery is reported for full-vs-140M, full-vs-210M, and 140M-vs-210M.
  exact.comparisons <- tribble(
    ~reference_condition, ~query_condition,
    "full_depth", "downsample_140M",
    "full_depth", "downsample_210M",
    "downsample_140M", "downsample_210M"
  )
  df.exact.recovery <- pmap_dfr(
    exact.comparisons,
    ~ summarise_depth_pair_exact_recovery(df.sample.membership, ..1, ..2)
  ) %>%
    left_join(
      df.metadata.210m %>% dplyr::select(sample, source_contacts),
      by = "sample"
    )

  df.radius.recovery <- bind_rows(
    summarise_depth_pair_radius_recovery(
      resource.list$full_depth,
      resource.list$downsample_140M,
      "full_depth",
      "downsample_140M"
    ),
    summarise_depth_pair_radius_recovery(
      resource.list$full_depth,
      resource.list$downsample_210M,
      "full_depth",
      "downsample_210M"
    ),
    summarise_depth_pair_radius_recovery(
      resource.list$downsample_140M,
      resource.list$downsample_210M,
      "downsample_140M",
      "downsample_210M"
    )
  ) %>%
    left_join(
      df.metadata.210m %>% dplyr::select(sample, source_contacts),
      by = "sample"
    )

  # Construct one joint approximate-locus map so locus IDs are directly
  # comparable across full depth, 140M, and 210M.
  df.joint.sample <- imap_dfr(
    resource.list,
    ~ .x$sample %>%
      filter(passes_lt2mb) %>%
      mutate(
        original_sample = sample,
        sample = str_c(.y, sample, sep = "::"),
        strain = str_c(.y, strain, sep = "::"),
        sample_loop_id = str_c(.y, sample_loop_id, sep = "::")
      )
  )
  joint.pooled <- build_pooled_hiccups_loop_resource(df.joint.sample)
  joint.approximate <- build_approximate_loop_loci(joint.pooled$universe)

  df.approximate.membership <- df.sample.membership %>%
    left_join(
      joint.approximate$map %>%
        dplyr::select(loop_id, approximate_loop_locus_id),
      by = "loop_id"
    )
  stopifnot(!any(is.na(df.approximate.membership$approximate_loop_locus_id)))

  # Exact pairwise metrics retain resolution strata; approximate loci are
  # compared across all 5/10/25-kb calls after joint canonicalization.
  df.exact.membership.by.resolution <- df.sample.membership %>%
    dplyr::select(condition, sample, resolution, loop_id) %>%
    distinct()
  df.exact.pairwise <- bind_rows(
    df.exact.membership.by.resolution,
    df.exact.membership.by.resolution %>% mutate(resolution = "ALL")
  ) %>%
    distinct() %>%
    summarise_pairwise_library_similarity("loop_id", "exact_loop_call")

  df.approximate.pairwise <- df.approximate.membership %>%
    transmute(
      condition, sample, resolution = "ALL_cross_resolution",
      approximate_loop_locus_id
    ) %>%
    distinct() %>%
    summarise_pairwise_library_similarity(
      "approximate_loop_locus_id", "approximate_loop_locus"
    )

  df.exact.support <- df.sample.membership %>%
    dplyr::select(condition, sample, loop_id) %>%
    summarise_library_support_distribution("loop_id", "exact_loop_call")
  df.approximate.support <- df.approximate.membership %>%
    dplyr::select(condition, sample, approximate_loop_locus_id) %>%
    summarise_library_support_distribution(
      "approximate_loop_locus_id", "approximate_loop_locus"
    )

  df.recurrent.core.stability <- bind_rows(
    summarise_recurrent_core_stability(
      df.sample.membership,
      "loop_id",
      "exact_loop_call"
    ),
    summarise_recurrent_core_stability(
      df.approximate.membership,
      "approximate_loop_locus_id",
      "approximate_loop_locus"
    )
  )

  # Preserve the actual recurrent features at the primary >=2-library cutoff.
  df.recurrent.exact.detail <- df.sample.membership %>%
    filter(condition %in% c("downsample_140M", "downsample_210M")) %>%
    distinct(condition, sample, loop_id) %>%
    count(condition, loop_id, name = "n_supporting_libraries") %>%
    filter(n_supporting_libraries >= 2L) %>%
    pivot_wider(
      names_from = condition,
      values_from = n_supporting_libraries,
      values_fill = 0L
    ) %>%
    filter(downsample_140M >= 2L, downsample_210M >= 2L)

  df.recurrent.approximate.detail <- df.approximate.membership %>%
    filter(condition %in% c("downsample_140M", "downsample_210M")) %>%
    distinct(condition, sample, approximate_loop_locus_id) %>%
    count(
      condition,
      approximate_loop_locus_id,
      name = "n_supporting_libraries"
    ) %>%
    filter(n_supporting_libraries >= 2L) %>%
    pivot_wider(
      names_from = condition,
      values_from = n_supporting_libraries,
      values_fill = 0L
    ) %>%
    filter(downsample_140M >= 2L, downsample_210M >= 2L) %>%
    left_join(
      joint.approximate$summary,
      by = "approximate_loop_locus_id"
    )

  output.tables.210m <- list(
    depth_series_loop_counts = df.loop.count,
    depth_series_loop_count_summary = df.loop.count.summary,
    depth_series_exact_recovery = df.exact.recovery,
    depth_series_radius_recovery = df.radius.recovery,
    depth_series_pairwise_exact = df.exact.pairwise,
    depth_series_pairwise_approximate_locus = df.approximate.pairwise,
    depth_series_library_support_exact = df.exact.support,
    depth_series_library_support_approximate_locus = df.approximate.support,
    depth_series_recurrent_core_stability = df.recurrent.core.stability,
    depth_series_recurrent_exact_core_detail = df.recurrent.exact.detail,
    depth_series_recurrent_approximate_core_detail =
      df.recurrent.approximate.detail,
    depth_series_joint_approximate_locus_method = joint.approximate$method
  )
  write_downsampling_tables(output.tables.210m, depth.series.output.dir)

  df.method.210m <- tribble(
    ~analysis_component, ~definition,
    "matched_eight_library_design",
    paste0(
      "Full-depth, 140M, and 210M resources use the same eight libraries. ",
      "74AA and A2DB are excluded because their full-depth usable-contact ",
      "counts are below 210M."
    ),
    "exact_recovery",
    "Sample, chromosome, anchor coordinates, and resolution must match.",
    "tolerance_recovery",
    paste0(
      "Both HiCCUPS centroids must be within 20/20/50 kb at 5/10/25 kb ",
      "in the same sample, resolution, and chromosome pair."
    ),
    "approximate_loop_locus",
    paste0(
      "Exact calls are jointly grouped across depth conditions using the ",
      "all-pairs-constrained 20/20/50-kb two-centroid rule; exact calls ",
      "remain unchanged."
    ),
    "recurrent_shared_core",
    paste0(
      "Primary recurrent features have support from at least two of eight ",
      "libraries independently at both 140M and 210M."
    )
  )
  readr::write_tsv(
    df.method.210m,
    file.path(depth.series.output.dir, "depth_series_method_definitions.tsv")
  )

  plot.depth.series <- ggplot(
    df.loop.count %>% filter(resolution == "ALL"),
    aes(x = reorder(sample, source_contacts), y = n_loop_calls, fill = condition)
  ) +
    geom_col(position = position_dodge(width = 0.8)) +
    labs(
      x = "Library ordered by original usable-contact depth",
      y = "Exact HiCCUPS calls <2 Mb",
      fill = "Depth condition"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(
    file.path(depth.series.output.dir, "depth_series_loop_counts.png"),
    plot.depth.series,
    width = 10,
    height = 6,
    dpi = 300
  )

  run.metadata.210m <- tibble(
    run_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    n_libraries = nrow(df.metadata.210m),
    included_libraries = str_c(df.metadata.210m$sample, collapse = ";"),
    excluded_below_210M = "74AA;A2DB",
    target_140M_contacts = 140000000,
    target_210M_contacts = 210000000,
    full_depth_root = full.depth.root,
    downsample_140M_root = downsample.root,
    downsample_210M_root = downsample.210m.root
  )
  readr::write_tsv(
    run.metadata.210m,
    file.path(depth.series.output.dir, "depth_series_run_metadata.tsv")
  )

  message("Full/140M/210M eight-library sensitivity analysis complete.")
  message("Results written to: ", depth.series.output.dir)
  print(df.loop.count.summary %>% filter(resolution == "ALL"))
  print(
    df.recurrent.core.stability %>%
      filter(minimum_library_support == 2L)
  )

  invisible(output.tables.210m)
}

target.contacts.210m <- 210000000
eligible.210m.samples <- c(
  "592BB", "607", "D765A", "DA08A",
  "DA21A", "DA68A", "DBA9A", "DE8BA"
)
df.sample.metadata.210m <- df.sample.metadata %>%
  filter(sample %in% eligible.210m.samples) %>%
  mutate(
    target_contacts = target.contacts.210m,
    sampling_fraction = target_contacts / source_contacts
  ) %>%
  arrange(match(sample, eligible.210m.samples))

depth.series.output.dir <- path.expand(Sys.getenv(
  "DOWNSAMPLING_DEPTH_SERIES_OUTPUT_DIR",
  unset = file.path(analysis.dir, "results", "full_140M_210M")
))
dir.create(depth.series.output.dir, recursive = TRUE, showWarnings = FALSE)

downsample.210m.root <- resolve_optional_downsample_root(c(
  Sys.getenv("DOWNSAMPLED_210M_HICCUPS_ROOT", unset = ""),
  file.path(input.dir, "juicer_downsample_q30_210M"),
  path.expand(
    "~/dropbox/Gateway_to_Hao/enhancer/data/juicer_downsample_q30_210M"
  ),
  file.path(dirname(r.files.dir), "data", "juicer_downsample_q30_210M")
))

if (is.na(downsample.210m.root)) {
  readr::write_tsv(
    df.sample.metadata.210m %>%
      transmute(
        sample,
        strain,
        source_contacts,
        target_contacts,
        status = "pending_210M_HiCCUPS_input_root"
      ),
    file.path(depth.series.output.dir, "depth_series_210M_sample_status.tsv")
  )
  message(
    "210M depth-series analysis is pending. Set ",
    "DOWNSAMPLED_210M_HICCUPS_ROOT after the eight complete HiCCUPS result ",
    "directories have been copied locally."
  )
} else {
  run_210m_depth_series_analysis(
    df.sample.metadata.210m,
    downsample.210m.root,
    depth.series.output.dir
  )
}
}
