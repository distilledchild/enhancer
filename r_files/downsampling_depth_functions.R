# Shared helpers for matched full-depth, 140M, and 250M Hi-C sensitivity analyses.
# Load tidyverse, GenomicRanges, and r_files/funcs.R before sourcing this file.

# Select the first existing directory from an ordered candidate list.
resolve_depth_directory <- function(candidates, label, required = TRUE) {
  candidates <- unique(path.expand(candidates[nzchar(candidates)]))
  candidates <- candidates[dir.exists(candidates)]
  if (length(candidates) == 0L) {
    if (!required) return(NA_character_)
    stop("Cannot locate ", label, ".", call. = FALSE)
  }
  normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE)
}

# Add an ALL-resolution stratum while retaining the original records.
add_all_resolution <- function(df) {
  bind_rows(df, df %>% mutate(resolution = "ALL"))
}

# Calculate the coefficient of variation when a meaningful denominator exists.
coefficient_of_variation <- function(x) {
  if (length(x) < 2L || mean(x) <= 0) NA_real_ else sd(x) / mean(x)
}

# Return a stable Spearman estimate for finite, nonconstant vectors.
safe_spearman <- function(x, y) {
  valid <- is.finite(x) & is.finite(y)
  if (sum(valid) < 3L || length(unique(x[valid])) < 2L ||
      length(unique(y[valid])) < 2L) {
    return(NA_real_)
  }
  suppressWarnings(cor(x[valid], y[valid], method = "spearman"))
}

# Return the sample size and asymptotic P value accompanying a Spearman estimate.
safe_spearman_test <- function(x, y) {
  valid <- is.finite(x) & is.finite(y)
  n.valid <- sum(valid)
  if (n.valid < 3L || length(unique(x[valid])) < 2L ||
      length(unique(y[valid])) < 2L) {
    return(tibble(
      n_spearman = n.valid,
      spearman_rho = NA_real_,
      spearman_p_value = NA_real_
    ))
  }
  result <- suppressWarnings(cor.test(
    x[valid], y[valid], method = "spearman", exact = FALSE
  ))
  tibble(
    n_spearman = n.valid,
    spearman_rho = unname(result$estimate),
    spearman_p_value = result$p.value
  )
}

# Resolve full-depth and downsampled HiCCUPS files for one target depth.
build_depth_file_status <- function(
  df.metadata,
  full.depth.root,
  downsample.root,
  target.label
) {
  df.metadata %>%
    mutate(
      target_label = target.label,
      full_file = file.path(
        full.depth.root, sample, "hiccups_5k10k25k", "merged_loops.bedpe"
      ),
      downsample_result_dir = file.path(
        downsample.root,
        sample,
        str_c(
          sample, "_inter_30_", target.label, "_seed", seed,
          ".hiccups.5k10k25k"
        )
      ),
      downsample_file = file.path(downsample_result_dir, "merged_loops.bedpe"),
      downsample_qc_file = file.path(
        downsample.root, sample, "downsampling_qc.tsv"
      ),
      hic_creation_qc_file = file.path(
        downsample.root, sample, "hic_creation_qc.tsv"
      ),
      hic_file = file.path(
        downsample.root,
        sample,
        str_c(sample, "_inter_30_", target.label, "_seed", seed, ".hic")
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
      full_file_ready = file.exists(full_file) & file.info(full_file)$size > 0L,
      downsample_qc_ready = file.exists(downsample_qc_file) &
        file.info(downsample_qc_file)$size > 0L,
      hic_creation_qc_ready = file.exists(hic_creation_qc_file) &
        file.info(hic_creation_qc_file)$size > 0L,
      hic_file_ready = file.exists(hic_file) & file.info(hic_file)$size > 0L,
      downsample_merged_ready = file.exists(downsample_file) &
        file.info(downsample_file)$size > 0L,
      downsample_5k_ready = file.exists(downsample_postprocessed_5k) &
        file.info(downsample_postprocessed_5k)$size > 0L,
      downsample_10k_ready = file.exists(downsample_postprocessed_10k) &
        file.info(downsample_postprocessed_10k)$size > 0L,
      downsample_25k_ready = file.exists(downsample_postprocessed_25k) &
        file.info(downsample_postprocessed_25k)$size > 0L,
      analysis_ready = full_file_ready & downsample_merged_ready &
        downsample_5k_ready & downsample_10k_ready & downsample_25k_ready,
      provenance_ready = downsample_qc_ready & hic_creation_qc_ready &
        hic_file_ready,
      status = if_else(
        analysis_ready & provenance_ready,
        "included_complete_5k10k25k",
        if_else(
          analysis_ready,
          "analysis_ready_but_provenance_incomplete",
          "excluded_incomplete_5k10k25k"
        )
      )
    )
}

# Read and validate the one-row downsampling QC files for one target depth.
read_depth_downsampling_qc <- function(df.status, expected.target) {
  map_dfr(seq_len(nrow(df.status)), function(i) {
    row <- df.status[i, ]
    qc <- readr::read_tsv(
      row$downsample_qc_file,
      col_types = cols(
        sample = col_character(), seed = col_double(),
        source_contacts = col_double(), downsampled_contacts = col_double(),
        output = col_character()
      ),
      show_col_types = FALSE
    )
    if (nrow(qc) != 1L) {
      stop("Expected one QC row in ", row$downsample_qc_file, call. = FALSE)
    }
    qc %>%
      transmute(
        target_label = row$target_label,
        sample,
        seed = as.integer(seed),
        source_contacts = as.numeric(source_contacts),
        target_contacts = as.numeric(downsampled_contacts),
        expected_seed = row$seed,
        expected_source_contacts = row$source_contacts,
        expected_target_contacts = expected.target,
        source_output_on_hpc = output,
        qc_matches_metadata = (
          sample == row$sample & seed == row$seed &
            source_contacts == row$source_contacts &
            target_contacts == expected.target
        )
      )
  })
}

# Build sample-level calls and the exact pooled resource for one depth condition.
build_depth_loop_resource <- function(
  file.metadata,
  condition,
  max.loop.distance.bp = 2000000L
) {
  df.sample <- read_hiccups_loop_files(file.metadata) %>%
    normalize_hiccups_loop_coordinates(max.loop.distance = max.loop.distance.bp) %>%
    mutate(condition = condition, .before = 1)
  pooled <- build_pooled_hiccups_loop_resource(df.sample)
  list(
    sample = df.sample,
    support = pooled$support,
    distinct = pooled$distinct,
    universe = pooled$universe %>% mutate(condition = condition, .before = 1)
  )
}

# Summarize exact-coordinate recovery for an ordered pair of depth conditions.
summarise_exact_depth_recovery <- function(
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
      n_reference = n_distinct(loop_id[condition == reference.condition]),
      n_query = n_distinct(loop_id[condition == query.condition]),
      n_exact_shared = length(intersect(
        loop_id[condition == reference.condition],
        loop_id[condition == query.condition]
      )),
      pct_reference_recovered = if_else(
        n_reference > 0L, 100 * n_exact_shared / n_reference, NA_real_
      ),
      exact_jaccard = if_else(
        n_reference + n_query - n_exact_shared > 0L,
        n_exact_shared / (n_reference + n_query - n_exact_shared),
        NA_real_
      ),
      exact_overlap_coefficient = if_else(
        pmin(n_reference, n_query) > 0L,
        n_exact_shared / pmin(n_reference, n_query),
        NA_real_
      ),
      .groups = "drop"
    )
}

# Measure position-tolerant recovery using the HiCCUPS 20/20/50-kb radii.
summarise_radius_depth_recovery <- function(
  reference.resource,
  query.resource,
  reference.condition,
  query.condition
) {
  radius.by.resolution <- c("5K" = 20000L, "10K" = 20000L, "25K" = 50000L)
  groups <- bind_rows(
    reference.resource$sample %>% distinct(sample, resolution),
    query.resource$sample %>% distinct(sample, resolution)
  ) %>%
    distinct() %>%
    arrange(sample, resolution)

  df.by.resolution <- pmap_dfr(groups, function(sample, resolution) {
    df.reference <- reference.resource$sample %>%
      filter(
        passes_lt2mb,
        .data$sample == .env$sample,
        .data$resolution == .env$resolution
      ) %>%
      mutate(call_index = row_number())
    df.query <- query.resource$sample %>%
      filter(
        passes_lt2mb,
        .data$sample == .env$sample,
        .data$resolution == .env$resolution
      ) %>%
      mutate(call_index = row_number())
    radius <- unname(radius.by.resolution[[resolution]])

    if (nrow(df.reference) == 0L || nrow(df.query) == 0L) {
      return(tibble(
        sample, resolution, match_radius_bp = radius,
        n_reference = nrow(df.reference), n_query = nrow(df.query),
        n_reference_radius_recovered = 0L, n_query_radius_matched = 0L
      ))
    }

    common.seqinfo <- Seqinfo(seqnames = union(
      unique(df.reference$chr1), unique(df.query$chr1)
    ))
    gr.reference <- GRanges(
      seqnames = df.reference$chr1,
      ranges = IRanges(
        start = df.reference$hiccups_centroid1_1based,
        end = df.reference$hiccups_centroid1_1based
      ),
      seqinfo = common.seqinfo
    )
    gr.query.window <- GRanges(
      seqnames = df.query$chr1,
      ranges = IRanges(
        start = pmax(1L, df.query$hiccups_centroid1_1based - radius),
        end = df.query$hiccups_centroid1_1based + radius
      ),
      seqinfo = common.seqinfo
    )
    hits <- findOverlaps(gr.reference, gr.query.window, ignore.strand = TRUE)
    reference.index <- queryHits(hits)
    query.index <- subjectHits(hits)
    valid <- if (length(hits) == 0L) logical() else (
      df.reference$chr2[reference.index] == df.query$chr2[query.index] &
        abs(
          df.reference$hiccups_centroid2_1based[reference.index] -
            df.query$hiccups_centroid2_1based[query.index]
        ) <= radius
    )

    tibble(
      sample,
      resolution,
      match_radius_bp = radius,
      n_reference = nrow(df.reference),
      n_query = nrow(df.query),
      n_reference_radius_recovered = n_distinct(reference.index[valid]),
      n_query_radius_matched = n_distinct(query.index[valid])
    )
  }) %>%
    mutate(
      reference_condition = reference.condition,
      query_condition = query.condition,
      pct_reference_radius_recovered = if_else(
        n_reference > 0L,
        100 * n_reference_radius_recovered / n_reference,
        NA_real_
      ),
      pct_query_radius_matched = if_else(
        n_query > 0L, 100 * n_query_radius_matched / n_query, NA_real_
      )
    )

  df.all <- df.by.resolution %>%
    group_by(sample, reference_condition, query_condition) %>%
    summarise(
      resolution = "ALL",
      match_radius_bp = NA_integer_,
      n_reference = sum(n_reference),
      n_query = sum(n_query),
      n_reference_radius_recovered = sum(n_reference_radius_recovered),
      n_query_radius_matched = sum(n_query_radius_matched),
      pct_reference_radius_recovered = 100 *
        n_reference_radius_recovered / n_reference,
      pct_query_radius_matched = 100 * n_query_radius_matched / n_query,
      .groups = "drop"
    )

  bind_rows(df.by.resolution, df.all) %>%
    arrange(reference_condition, query_condition, resolution, sample)
}

# Compare unordered library pairs within each condition and resolution stratum.
summarise_pairwise_library_similarity <- function(
  df.membership,
  feature.column,
  feature.type
) {
  pmap_dfr(
    df.membership %>% distinct(condition, resolution),
    function(condition, resolution) {
      df.i <- df.membership %>%
        filter(
          .data$condition == .env$condition,
          .data$resolution == .env$resolution
        )
      pairs <- combn(sort(unique(df.i$sample)), 2L, simplify = FALSE)
      map_dfr(pairs, function(pair) {
        set.1 <- unique(df.i[[feature.column]][df.i$sample == pair[[1]]])
        set.2 <- unique(df.i[[feature.column]][df.i$sample == pair[[2]]])
        n.shared <- length(intersect(set.1, set.2))
        n.union <- length(union(set.1, set.2))
        n.minimum <- min(length(set.1), length(set.2))
        tibble(
          condition, resolution, feature_type = feature.type,
          sample1 = pair[[1]], sample2 = pair[[2]],
          n_sample1 = length(set.1), n_sample2 = length(set.2),
          n_shared = n.shared, n_union = n.union,
          jaccard_similarity = if_else(
            n.union > 0L, n.shared / n.union, NA_real_
          ),
          overlap_coefficient = if_else(
            n.minimum > 0L, n.shared / n.minimum, NA_real_
          )
        )
      })
    }
  )
}

# Count features according to the number of supporting libraries per condition.
summarise_library_support_distribution <- function(
  df.membership,
  feature.column,
  feature.type
) {
  df.membership %>%
    distinct(condition, sample, .data[[feature.column]]) %>%
    count(condition, .data[[feature.column]], name = "n_supporting_libraries") %>%
    count(condition, n_supporting_libraries, name = "n_features") %>%
    group_by(condition) %>%
    mutate(
      feature_type = feature.type,
      pct_features = 100 * n_features / sum(n_features)
    ) %>%
    ungroup() %>%
    dplyr::select(
      condition, feature_type, n_supporting_libraries, n_features, pct_features
    )
}

# Compare recurrent features that meet the same library-support cutoff twice.
summarise_recurrent_core_stability <- function(
  df.membership,
  feature.column,
  feature.type,
  condition.1 = "downsample_140M",
  condition.2 = "downsample_250M",
  maximum.support = 7L
) {
  df.support <- df.membership %>%
    filter(condition %in% c(condition.1, condition.2)) %>%
    distinct(condition, sample, .data[[feature.column]]) %>%
    count(condition, .data[[feature.column]], name = "n_supporting_libraries")

  map_dfr(seq_len(maximum.support), function(minimum.support) {
    set.1 <- df.support %>%
      filter(
        condition == condition.1,
        n_supporting_libraries >= minimum.support
      ) %>%
      pull(.data[[feature.column]]) %>% unique()
    set.2 <- df.support %>%
      filter(
        condition == condition.2,
        n_supporting_libraries >= minimum.support
      ) %>%
      pull(.data[[feature.column]]) %>% unique()
    n.shared <- length(intersect(set.1, set.2))
    n.union <- length(union(set.1, set.2))
    n.minimum <- min(length(set.1), length(set.2))
    tibble(
      feature_type = feature.type,
      condition_1 = condition.1,
      condition_2 = condition.2,
      minimum_library_support = minimum.support,
      n_condition_1 = length(set.1),
      n_condition_2 = length(set.2),
      n_shared = n.shared,
      jaccard = if_else(n.union > 0L, n.shared / n.union, NA_real_),
      overlap_coefficient = if_else(
        n.minimum > 0L, n.shared / n.minimum, NA_real_
      )
    )
  })
}

# Build promoter windows and the TSS-excluded ATAC resource once per run.
build_depth_regulatory_context <- function(
  df.transcript,
  df.promoter,
  gr.atac,
  promoter.flank.bp = 1000L
) {
  df.true.tss <- build_true_tss_annotation(df.transcript)
  gr.true.tss <- create_true_tss_granges(df.true.tss)
  gr.epd.tss <- create_epd_tss_granges(df.promoter)
  gr.true.window <- expand_tss_to_promoter_windows(
    gr.true.tss, flank.bp = promoter.flank.bp
  )
  gr.epd.window <- expand_tss_to_promoter_windows(
    gr.epd.tss, flank.bp = promoter.flank.bp
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
  gr.atac.union <- reduce(gr.atac, ignore.strand = TRUE)
  common.seqlevels <- intersect(seqlevels(gr.atac.union), seqlevels(gr.tss.exclusion))
  gr.atac.common <- keepSeqlevels(
    gr.atac.union, common.seqlevels, pruning.mode = "coarse"
  )
  gr.exclusion.common <- keepSeqlevels(
    gr.tss.exclusion, common.seqlevels, pruning.mode = "coarse"
  )
  list(
    df_true_tss = df.true.tss,
    gr_true_window = gr.true.window,
    gr_epd_window = gr.epd.window,
    gr_atac_union = gr.atac.common,
    gr_atac_non_tss = setdiff(
      gr.atac.common, gr.exclusion.common, ignore.strand = TRUE
    )
  )
}

# Reapply the revised promoter/TSS and ATAC evidence definitions to one pool.
annotate_depth_loop_resource <- function(
  df.loop,
  condition,
  context,
  df.promoter,
  promoter.flank.bp = 1000L,
  minimum.atac.overlap.bp = 50L
) {
  gr.anchor.by.side <- create_loop_anchor_granges_by_side(df.loop)
  direct.tier <- build_direct_promoter_tss_tier(
    gr.loop.anchor.by.side = gr.anchor.by.side,
    gr.true.tss.annotation = context$gr_true_window,
    gr.epd.annotation = context$gr_epd_window,
    df.true.tss = context$df_true_tss,
    df.epd.promoter = df.promoter,
    df.loop.universe = df.loop,
    evidence.definition = "primary_1kb",
    promoter.window.flank.bp = promoter.flank.bp
  )

  df.promoter.window <- direct.tier$combined_overlap %>%
    distinct(
      loop_id, resolution, anchor_side,
      annotation_chr, annotation_start, annotation_end
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
    context$gr_atac_union,
    minimum.overlap.bp = minimum.atac.overlap.bp
  ) %>%
    group_by(loop_id, resolution, anchor_side) %>%
    summarise(
      promoter_window_atac_ge50 = any(has_feature_overlap_ge_minimum),
      .groups = "drop"
    )

  df.non.tss.atac <- summarise_anchor_interval_overlap(
    combine_loop_anchor_granges_by_side(gr.anchor.by.side),
    context$gr_atac_non_tss,
    minimum.overlap.bp = minimum.atac.overlap.bp
  ) %>%
    transmute(
      loop_id, resolution, anchor_side,
      residual_non_tss_atac_ge50 = has_feature_overlap_ge_minimum,
      residual_non_tss_atac_bp = feature_overlap_bp
    )

  df.anchor <- direct.tier$summary$anchor_summary %>%
    left_join(df.promoter.atac, by = c("loop_id", "resolution", "anchor_side")) %>%
    left_join(df.non.tss.atac, by = c("loop_id", "resolution", "anchor_side")) %>%
    mutate(
      promoter_window_atac_ge50 = coalesce(promoter_window_atac_ge50, FALSE),
      residual_non_tss_atac_ge50 = coalesce(residual_non_tss_atac_ge50, FALSE),
      residual_non_tss_atac_bp = coalesce(residual_non_tss_atac_bp, 0L)
    )
  df.anchor.wide <- df.anchor %>%
    dplyr::select(
      loop_id, anchor_side, promoter_window_atac_ge50,
      residual_non_tss_atac_ge50, residual_non_tss_atac_bp
    ) %>%
    pivot_wider(
      names_from = anchor_side,
      values_from = -loop_id,
      names_glue = "{.value}_{anchor_side}"
    )

  df.loop.evidence <- direct.tier$summary$loop_summary %>%
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
        n_direct_anchor_sides == 2L & promoter_window_atac_ge50_anchor1 &
          residual_non_tss_atac_ge50_anchor2
      ),
      supports_anchor2_promoter_to_anchor1_regulatory = (
        n_direct_anchor_sides == 2L & promoter_window_atac_ge50_anchor2 &
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

  df.gene.assignment <- direct.tier$summary$gene_assignment %>%
    left_join(
      df.loop.evidence %>%
        dplyr::select(
          loop_id, condition, revised_putative_regulatory_support,
          revised_major_category
        ),
      by = "loop_id"
    ) %>%
    mutate(ensembl_gene_id = gene_id, gene_symbol = gene_name)
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
    anchor = df.anchor,
    gene_assignment = df.gene.assignment,
    gene_count = df.gene.count
  )
}

# Compare gene presence, loop counts, and ranks for one ordered depth pair.
summarise_gene_depth_stability <- function(
  df.gene.count,
  reference.condition,
  query.condition,
  count.column = "n_exact_loop_calls",
  gene.count.unit = "exact_resolution_specific_call"
) {
  df.pair <- full_join(
    df.gene.count %>%
      filter(condition == reference.condition) %>%
      transmute(
        ensembl_gene_id,
        gene_symbol_reference = gene_symbol,
        n_reference = .data[[count.column]]
      ),
    df.gene.count %>%
      filter(condition == query.condition) %>%
      transmute(
        ensembl_gene_id,
        gene_symbol_query = gene_symbol,
        n_query = .data[[count.column]]
      ),
    by = "ensembl_gene_id"
  ) %>%
    mutate(
      gene_symbol = coalesce(gene_symbol_reference, gene_symbol_query),
      n_reference = coalesce(n_reference, 0L),
      n_query = coalesce(n_query, 0L),
      present_reference = n_reference > 0L,
      present_query = n_query > 0L,
      rank_reference = rank(-n_reference, ties.method = "average"),
      rank_query = rank(-n_query, ties.method = "average")
    )
  n.shared <- sum(df.pair$present_reference & df.pair$present_query)
  list(
    detail = df.pair %>%
      mutate(
        reference_condition = reference.condition,
        query_condition = query.condition,
        gene_count_unit = gene.count.unit,
        .before = 1
      ) %>%
      dplyr::select(
        reference_condition, query_condition, gene_count_unit,
        ensembl_gene_id, gene_symbol, n_reference, n_query,
        present_reference, present_query, rank_reference, rank_query
      ),
    summary = tibble(
      reference_condition = reference.condition,
      query_condition = query.condition,
      gene_count_unit = gene.count.unit,
      n_genes_reference = sum(df.pair$present_reference),
      n_genes_query = sum(df.pair$present_query),
      n_genes_shared = n.shared,
      gene_presence_jaccard = n.shared / nrow(df.pair),
      spearman_loop_counts_union_zero_filled = safe_spearman(
        df.pair$n_reference, df.pair$n_query
      ),
      spearman_gene_ranks_union = safe_spearman(
        df.pair$rank_reference, df.pair$rank_query
      )
    )
  )
}

# Write a named list of data frames as compact TSV outputs.
write_depth_tables <- function(tables, output.dir) {
  walk2(
    tables,
    names(tables),
    ~ readr::write_tsv(.x, file.path(output.dir, str_c(.y, ".tsv")))
  )
}
