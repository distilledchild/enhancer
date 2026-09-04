# Integrate pooled rat frontal-cortex Hi-C promoter-distal contacts with
# external PFC snATAC, pooled-PFC csRNA TSR, and naive-PFC bulk RNA evidence.

suppressPackageStartupMessages({
  library(dplyr)
  library(GenomicRanges)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
})

current_script_path <- function() {
  file.args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(file.args)) return(normalizePath(getwd(), mustWork = TRUE))
  normalizePath(sub("^--file=", "", file.args[[1]]), mustWork = TRUE)
}

script.dir <- dirname(current_script_path())
revision.dir <- dirname(script.dir)
revision.main.dir <- file.path(revision.dir, "revision_main")
public.input.dir <- file.path(script.dir, "inputs", "public")
output.dir <- file.path(script.dir, "results")
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

paths <- c(
  loop_resource = file.path(
    revision.main.dir, "results", "revised_pooled_loop_annotation_resource.tsv"
  ),
  direct_assignment = file.path(
    revision.main.dir, "results", "revised_direct_loop_gene_assignments.tsv"
  ),
  ensembl_transcript = file.path(
    revision.main.dir, "cache_data", "df.transcript.ensembl.rn7.1based.rds"
  ),
  epd_promoter = file.path(
    revision.main.dir, "cache_data", "df.promoter.epd.rn7.1based.rds"
  ),
  atac = file.path(
    revision.main.dir, "cache_data", "gr.atac.rn7.1based.rds"
  ),
  pfc_tsr = file.path(public.input.dir, "GSE193757_tsr.pfc.bed.gz"),
  gene_counts = file.path(public.input.dir, "GSE193757_gene_counts.txt.gz")
)
missing.paths <- paths[!file.exists(paths)]
if (length(missing.paths)) {
  stop("Missing required input(s):\n", paste(missing.paths, collapse = "\n"))
}

options(scipen = 999)
n.resamples <- as.integer(Sys.getenv("DISTAL_ANCHOR_RESAMPLES", "1000"))
random.seed <- as.integer(Sys.getenv("DISTAL_ANCHOR_SEED", "20260903"))
detected.cores <- parallel::detectCores(logical = TRUE)
if (is.na(detected.cores)) detected.cores <- 1L
n.cores <- as.integer(Sys.getenv(
  "DISTAL_ANCHOR_CORES",
  as.character(max(1L, min(4L, detected.cores - 1L)))
))
atac.minimum.overlap.bp <- 50L
tss.promoter.flank.bp <- 1000L
distal.tsr.minimum.distance.bp <- 3000L
expression.cpm.threshold <- 1

distance.breaks <- c(-Inf, 50000, 100000, 250000, 500000, 1000000, Inf)
distance.labels <- c(
  "le50kb", "gt50_le100kb", "gt100_le250kb",
  "gt250_le500kb", "gt500_le1000kb", "gt1000kb"
)

interval_overlap_stats <- function(gr.query, gr.feature) {
  result <- tibble(total_bp = integer(length(gr.query)), max_bp = integer(length(gr.query)))
  if (!length(gr.query) || !length(gr.feature)) return(result)
  hits <- findOverlaps(gr.query, gr.feature, ignore.strand = TRUE)
  if (!length(hits)) return(result)
  overlap.bp <- width(pintersect(
    gr.query[queryHits(hits)], gr.feature[subjectHits(hits)], ignore.strand = TRUE
  ))
  total.by.query <- tapply(overlap.bp, queryHits(hits), sum)
  max.by.query <- tapply(overlap.bp, queryHits(hits), max)
  result$total_bp[as.integer(names(total.by.query))] <- as.integer(total.by.query)
  result$max_bp[as.integer(names(max.by.query))] <- as.integer(max.by.query)
  result
}

count_points_in_ranges <- function(gr.ranges, gr.points) {
  answer <- integer(length(gr.ranges))
  if (!length(gr.ranges) || !length(gr.points)) return(answer)
  hits <- findOverlaps(gr.ranges, gr.points, ignore.strand = TRUE)
  if (!length(hits)) return(answer)
  counts <- table(queryHits(hits))
  answer[as.integer(names(counts))] <- as.integer(counts)
  answer
}

count_points_in_supported_feature_overlaps <- function(
  gr.ranges, gr.features, gr.points, minimum.overlap.bp
) {
  answer <- integer(length(gr.ranges))
  if (!length(gr.ranges) || !length(gr.features) || !length(gr.points)) {
    return(answer)
  }
  range.feature.hits <- findOverlaps(
    gr.ranges, gr.features, ignore.strand = TRUE
  )
  if (!length(range.feature.hits)) return(answer)
  overlap.ranges <- pintersect(
    gr.ranges[queryHits(range.feature.hits)],
    gr.features[subjectHits(range.feature.hits)],
    ignore.strand = TRUE
  )
  eligible <- width(overlap.ranges) >= minimum.overlap.bp
  if (!any(eligible)) return(answer)
  overlap.ranges <- overlap.ranges[eligible]
  mcols(overlap.ranges)$anchor_index <- queryHits(range.feature.hits)[eligible]
  point.hits <- findOverlaps(overlap.ranges, gr.points, ignore.strand = TRUE)
  if (!length(point.hits)) return(answer)
  anchor.point.pairs <- unique(data.frame(
    anchor_index = mcols(overlap.ranges)$anchor_index[queryHits(point.hits)],
    point_index = subjectHits(point.hits)
  ))
  counts <- table(anchor.point.pairs$anchor_index)
  answer[as.integer(names(counts))] <- as.integer(counts)
  answer
}

add_match_keys <- function(df) {
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
        anchor_chr, resolution, anchor_side, distance_match_bin, sep = "|"
      ),
      fallback_match_key = str_c(
        anchor_chr, resolution, distance_match_bin, sep = "|"
      ),
      relaxed_match_key = str_c(anchor_chr, resolution, sep = "|")
    )
}

read_homer_tsr_bed <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  lines <- lines[!grepl("^(#|track\\s)", lines) & nzchar(lines)]
  x <- read_tsv(
    I(paste(lines, collapse = "\n")),
    col_names = c("chr", "start0", "end0", "tsr_id", "score", "strand"),
    col_types = "ciicdc",
    show_col_types = FALSE
  )
  x %>%
    mutate(
      start = start0 + 1L,
      end = end0,
      center = as.integer(floor((start + end) / 2))
    )
}

df.loop.resource <- read_tsv(paths[["loop_resource"]], show_col_types = FALSE)
df.direct.assignment <- read_tsv(paths[["direct_assignment"]], show_col_types = FALSE)

df.single.orientation <- df.direct.assignment %>%
  distinct(loop_id, resolution, anchor_side, opposite_anchor_side) %>%
  inner_join(
    df.loop.resource %>%
      filter(n_direct_anchor_sides == 1L) %>%
      select(
        loop_id, resolution, chr1, start1, end1, chr2, start2, end2,
        loop_distance, n_supporting_libraries, n_supporting_strains,
        supporting_samples, supporting_strains,
        existing_atac_ge50 = candidate_anchor_non_tss_atac_ge50
      ),
    by = c("loop_id", "resolution")
  ) %>%
  transmute(
    loop_id,
    resolution,
    promoter_anchor_side = anchor_side,
    distal_anchor_side = opposite_anchor_side,
    anchor_side = opposite_anchor_side,
    anchor_chr = if_else(opposite_anchor_side == "anchor1", chr1, chr2),
    anchor_start = if_else(opposite_anchor_side == "anchor1", start1, start2),
    anchor_end = if_else(opposite_anchor_side == "anchor1", end1, end2),
    loop_distance,
    n_supporting_libraries,
    n_supporting_strains,
    supporting_samples,
    supporting_strains,
    existing_atac_ge50
  ) %>%
  arrange(loop_id)

stopifnot(
  nrow(df.single.orientation) == 12295L,
  !anyDuplicated(df.single.orientation$loop_id),
  all(df.single.orientation$anchor_start <= df.single.orientation$anchor_end)
)

df.loop.genes <- df.direct.assignment %>%
  semi_join(df.single.orientation, by = "loop_id") %>%
  transmute(
    loop_id,
    resolution,
    promoter_anchor_side = anchor_side,
    gene_id = ensembl_gene_id,
    assigned_gene_symbol = gene_symbol
  ) %>%
  filter(!is.na(gene_id), nzchar(gene_id)) %>%
  distinct()

df.hic.support.balance <- df.loop.resource %>%
  filter(n_direct_anchor_sides %in% c(0L, 1L)) %>%
  mutate(
    loop_group = if_else(
      n_direct_anchor_sides == 1L,
      "one_sided_promoter_anchored",
      "no_direct_tss_control_pool"
    )
  ) %>%
  group_by(loop_group, resolution) %>%
  summarise(
    n_loops = n(),
    median_supporting_libraries = median(n_supporting_libraries),
    mean_supporting_libraries = mean(n_supporting_libraries),
    median_supporting_strains = median(n_supporting_strains),
    mean_supporting_strains = mean(n_supporting_strains),
    .groups = "drop"
  )

df.ensembl.transcript <- readRDS(paths[["ensembl_transcript"]])
df.epd.promoter <- readRDS(paths[["epd_promoter"]])
gr.atac <- readRDS(paths[["atac"]])

df.known.tss <- bind_rows(
  df.ensembl.transcript %>%
    transmute(chr, tss = if_else(strand == "+", transcript_start, transcript_end)),
  df.epd.promoter %>% transmute(chr, tss = epd_tss_start)
) %>%
  filter(!is.na(chr), !is.na(tss)) %>%
  distinct(chr, tss)

gr.known.tss <- GRanges(
  seqnames = df.known.tss$chr,
  ranges = IRanges(df.known.tss$tss, df.known.tss$tss)
)
gr.tss.exclusion <- GenomicRanges::reduce(GRanges(
  seqnames = df.known.tss$chr,
  ranges = IRanges(
    pmax(1L, df.known.tss$tss - tss.promoter.flank.bp),
    df.known.tss$tss + tss.promoter.flank.bp
  )
), ignore.strand = TRUE)

common.seqlevels <- base::Reduce(
  intersect,
  list(seqlevels(gr.atac), seqlevels(gr.tss.exclusion), seqlevels(gr.known.tss))
)
gr.atac <- keepSeqlevels(gr.atac, common.seqlevels, pruning.mode = "coarse")
gr.known.tss <- keepSeqlevels(gr.known.tss, common.seqlevels, pruning.mode = "coarse")
gr.tss.exclusion <- keepSeqlevels(
  gr.tss.exclusion, common.seqlevels, pruning.mode = "coarse"
)
gr.atac.non.tss <- setdiff(
  GenomicRanges::reduce(gr.atac, ignore.strand = TRUE),
  gr.tss.exclusion,
  ignore.strand = TRUE
)

df.tsr.source <- read_homer_tsr_bed(paths[["pfc_tsr"]])
df.tsr.excluded.seqlevels <- df.tsr.source %>%
  filter(!chr %in% common.seqlevels)
df.tsr <- df.tsr.source %>%
  filter(chr %in% common.seqlevels)
gr.tsr.center <- GRanges(
  seqnames = df.tsr$chr,
  ranges = IRanges(df.tsr$center, df.tsr$center)
)
nearest.hits <- distanceToNearest(gr.tsr.center, gr.known.tss, ignore.strand = TRUE)
tsr.distance <- rep(NA_integer_, length(gr.tsr.center))
tsr.distance[queryHits(nearest.hits)] <- abs(
  start(gr.tsr.center)[queryHits(nearest.hits)] -
    start(gr.known.tss)[subjectHits(nearest.hits)]
)
df.tsr <- df.tsr %>%
  mutate(
    nearest_unified_tss_distance_bp = tsr.distance,
    is_distal_gt3kb = nearest_unified_tss_distance_bp > distal.tsr.minimum.distance.bp,
    is_outside_tss_plus_minus_1kb =
      nearest_unified_tss_distance_bp > tss.promoter.flank.bp,
    center_in_non_tss_atac = overlapsAny(
      gr.tsr.center, gr.atac.non.tss, ignore.strand = TRUE
    ),
    is_distal_gt3kb_and_atac = is_distal_gt3kb & center_in_non_tss_atac
  )
gr.tsr.all <- gr.tsr.center
gr.tsr.distal <- gr.tsr.center[df.tsr$is_distal_gt3kb]
gr.tsr.distal.atac <- gr.tsr.center[df.tsr$is_distal_gt3kb_and_atac]

measure_anchor_evidence <- function(df.anchor) {
  gr.anchor <- GRanges(
    seqnames = df.anchor$anchor_chr,
    ranges = IRanges(df.anchor$anchor_start, df.anchor$anchor_end)
  )
  atac.stats <- interval_overlap_stats(gr.anchor, gr.atac.non.tss)
  n.tsr.all <- count_points_in_ranges(gr.anchor, gr.tsr.all)
  n.tsr.distal <- count_points_in_ranges(gr.anchor, gr.tsr.distal)
  n.tsr.distal.atac <- count_points_in_supported_feature_overlaps(
    gr.anchor,
    gr.atac.non.tss,
    gr.tsr.distal,
    atac.minimum.overlap.bp
  )

  bind_cols(
    df.anchor,
    tibble(
      non_tss_atac_total_overlap_bp = atac.stats$total_bp,
      non_tss_atac_max_overlap_bp = atac.stats$max_bp,
      atac_overlap_ge50 = atac.stats$max_bp >= atac.minimum.overlap.bp,
      n_pooled_pfc_tsr_centers = n.tsr.all,
      n_pooled_pfc_distal_tsr_gt3kb_centers = n.tsr.distal,
      n_pooled_pfc_distal_tsr_gt3kb_centers_in_atac = n.tsr.distal.atac,
      pooled_pfc_tsr_present = n.tsr.all > 0L,
      pooled_pfc_distal_tsr_gt3kb_present = n.tsr.distal > 0L,
      atac_and_distal_tsr_same_anchor =
        atac_overlap_ge50 & pooled_pfc_distal_tsr_gt3kb_present,
      colocated_atac_distal_tsr =
        atac_overlap_ge50 & n.tsr.distal.atac > 0L
    )
  )
}

df.observed <- df.single.orientation %>%
  add_match_keys() %>%
  measure_anchor_evidence()

atac.discordance <- sum(
  df.observed$atac_overlap_ge50 != df.observed$existing_atac_ge50,
  na.rm = TRUE
)
stopifnot(
  sum(df.observed$atac_overlap_ge50) == 10469L,
  atac.discordance == 0L
)

df.control.pool <- bind_rows(
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
  add_match_keys() %>%
  measure_anchor_evidence()

exact.pool.index <- split(seq_len(nrow(df.control.pool)), df.control.pool$exact_match_key)
fallback.pool.index <- split(
  seq_len(nrow(df.control.pool)), df.control.pool$fallback_match_key
)
relaxed.pool.index <- split(
  seq_len(nrow(df.control.pool)), df.control.pool$relaxed_match_key
)

pool_is_available <- function(index, key) {
  pool <- index[[key]]
  !is.null(pool) && length(pool) > 0L
}

df.match.availability <- df.observed %>%
  mutate(
    match_level = case_when(
      vapply(exact_match_key, function(key) {
        pool_is_available(exact.pool.index, key)
      }, logical(1)) ~ "exact_chr_resolution_side_distance_bin",
      vapply(fallback_match_key, function(key) {
        pool_is_available(fallback.pool.index, key)
      }, logical(1)) ~ "side_relaxed_same_distance_bin",
      vapply(relaxed_match_key, function(key) {
        pool_is_available(relaxed.pool.index, key)
      }, logical(1)) ~ "nearest_same_chr_resolution",
      TRUE ~ "no_control_available"
    )
  )
stopifnot(!any(df.match.availability$match_level == "no_control_available"))

df.match.availability.summary <- bind_rows(
  df.match.availability %>%
    count(resolution, match_level, name = "n_observed_anchors"),
  df.match.availability %>%
    count(match_level, name = "n_observed_anchors") %>%
    mutate(resolution = "ALL", .before = 1)
) %>%
  group_by(resolution) %>%
  mutate(
    percent_observed_anchors = 100 * n_observed_anchors / sum(n_observed_anchors)
  ) %>%
  ungroup() %>%
  arrange(factor(resolution, c("5K", "10K", "25K", "ALL")), match_level)

sample_matched_controls <- function() {
  selected <- integer(nrow(df.observed))
  match.level <- rep(NA_character_, nrow(df.observed))
  groups <- split(seq_len(nrow(df.observed)), df.observed$exact_match_key)
  for (key in names(groups)) {
    idx <- groups[[key]]
    pool <- exact.pool.index[[key]]
    if (!is.null(pool) && length(pool)) {
      selected[idx] <- sample(pool, length(idx), replace = TRUE)
      match.level[idx] <- "exact_chr_resolution_side_distance_bin"
    }
  }
  idx.remaining <- which(selected == 0L)
  if (length(idx.remaining)) {
    groups <- split(idx.remaining, df.observed$fallback_match_key[idx.remaining])
    for (key in names(groups)) {
      idx <- groups[[key]]
      pool <- fallback.pool.index[[key]]
      if (!is.null(pool) && length(pool)) {
        selected[idx] <- sample(pool, length(idx), replace = TRUE)
        match.level[idx] <- "side_relaxed_same_distance_bin"
      }
    }
  }
  idx.remaining <- which(selected == 0L)
  if (length(idx.remaining)) {
    for (idx in idx.remaining) {
      pool <- relaxed.pool.index[[df.observed$relaxed_match_key[[idx]]]]
      if (is.null(pool) || !length(pool)) stop("No matched control available.")
      distance.delta <- abs(
        log1p(df.control.pool$loop_distance[pool]) -
          log1p(df.observed$loop_distance[[idx]])
      )
      nearest <- pool[distance.delta == min(distance.delta)]
      selected[[idx]] <- sample(nearest, 1L)
      match.level[[idx]] <- "nearest_same_chr_resolution"
    }
  }
  bind_cols(
    df.control.pool[selected, ],
    tibble(observed_loop_id = df.observed$loop_id, match_level = match.level)
  )
}

metric.columns <- c(
  "atac_overlap_ge50",
  "pooled_pfc_tsr_present",
  "pooled_pfc_distal_tsr_gt3kb_present",
  "atac_and_distal_tsr_same_anchor",
  "colocated_atac_distal_tsr"
)
resolution.levels <- c("5K", "10K", "25K", "ALL")

summarise_resample <- function(control, resample) {
  map_dfr(resolution.levels, function(resolution.i) {
    idx <- if (resolution.i == "ALL") {
      seq_len(nrow(df.observed))
    } else {
      which(df.observed$resolution == resolution.i)
    }
    map_dfr(metric.columns, function(metric.i) {
      observed.positive <- df.observed[[metric.i]][idx]
      control.positive <- control[[metric.i]][idx]
      observed.only <- sum(observed.positive & !control.positive)
      control.only <- sum(!observed.positive & control.positive)
      tibble(
        resample,
        resolution = resolution.i,
        metric = metric.i,
        n_pairs = length(idx),
        observed_rate = mean(observed.positive),
        null_rate = mean(control.positive),
        paired_odds_ratio = (observed.only + 0.5) / (control.only + 0.5),
        n_observed_only = observed.only,
        n_control_only = control.only
      )
    })
  })
}

run_one_resample <- function(i) {
  set.seed(random.seed + i)
  summarise_resample(sample_matched_controls(), i)
}

message(
  "Running ", n.resamples, " matched-control resamples on ",
  n.cores, " core(s)."
)
resample.list <- if (.Platform$OS.type == "windows" || n.cores == 1L) {
  lapply(seq_len(n.resamples), run_one_resample)
} else {
  parallel::mclapply(
    seq_len(n.resamples), run_one_resample, mc.cores = n.cores,
    mc.preschedule = TRUE
  )
}
df.resample <- bind_rows(resample.list)

df.null.summary <- df.resample %>%
  group_by(resolution, metric) %>%
  summarise(
    n_pairs = dplyr::first(n_pairs),
    n_resamples = n(),
    observed_rate = dplyr::first(observed_rate),
    mean_null_rate = mean(null_rate),
    null_rate_sd = sd(null_rate),
    null_rate_q025 = quantile(null_rate, 0.025),
    null_rate_q975 = quantile(null_rate, 0.975),
    absolute_difference = observed_rate - mean_null_rate,
    enrichment_ratio = observed_rate / mean_null_rate,
    mean_paired_odds_ratio = mean(paired_odds_ratio),
    empirical_p = (1 + sum(null_rate >= observed_rate)) / (n() + 1),
    .groups = "drop"
  ) %>%
  mutate(fdr_bh = p.adjust(empirical_p, method = "BH")) %>%
  arrange(factor(resolution, resolution.levels), match(metric, metric.columns))

df.expression <- read_tsv(paths[["gene_counts"]], show_col_types = FALSE)
required.expression.columns <- c("Geneid", "PFC-Naive-r1", "PFC-Naive-r2")
stopifnot(
  all(required.expression.columns %in% names(df.expression)),
  !anyDuplicated(df.expression$Geneid)
)
library.size.r1 <- sum(df.expression[["PFC-Naive-r1"]])
library.size.r2 <- sum(df.expression[["PFC-Naive-r2"]])
df.expression <- df.expression %>%
  transmute(
    expression_gene_symbol = Geneid,
    pfc_naive_r1_count = .data[["PFC-Naive-r1"]],
    pfc_naive_r2_count = .data[["PFC-Naive-r2"]],
    pfc_naive_r1_cpm = 1e6 * pfc_naive_r1_count / library.size.r1,
    pfc_naive_r2_cpm = 1e6 * pfc_naive_r2_count / library.size.r2,
    pfc_naive_mean_cpm = (pfc_naive_r1_cpm + pfc_naive_r2_cpm) / 2,
    expressed_cpm1_both_naive_replicates =
      pfc_naive_r1_cpm >= expression.cpm.threshold &
      pfc_naive_r2_cpm >= expression.cpm.threshold
  )

df.current.gene.annotation <- df.ensembl.transcript %>%
  transmute(
    gene_id = str_match(attribute, 'gene_id "([^"]+)"')[, 2],
    current_gene_symbol = str_match(attribute, 'gene_name "([^"]+)"')[, 2]
  ) %>%
  filter(!is.na(gene_id), nzchar(gene_id)) %>%
  group_by(gene_id) %>%
  summarise(
    n_current_symbols_for_gene_id = n_distinct(
      current_gene_symbol[!is.na(current_gene_symbol)]
    ),
    current_gene_symbol = if_else(
      n_current_symbols_for_gene_id == 1L,
      dplyr::first(current_gene_symbol[!is.na(current_gene_symbol)]),
      NA_character_
    ),
    .groups = "drop"
  )

df.current.symbol.multiplicity <- df.current.gene.annotation %>%
  filter(!is.na(current_gene_symbol)) %>%
  distinct(gene_id, current_gene_symbol) %>%
  count(
    current_gene_symbol,
    name = "n_current_ensembl_gene_ids_for_symbol"
  )

df.loop.gene.evidence <- df.loop.genes %>%
  left_join(
    df.observed %>%
      select(
        loop_id, distal_anchor_side, anchor_chr, anchor_start, anchor_end,
        n_supporting_libraries, n_supporting_strains,
        supporting_samples, supporting_strains,
        non_tss_atac_total_overlap_bp, non_tss_atac_max_overlap_bp,
        atac_overlap_ge50, pooled_pfc_distal_tsr_gt3kb_present,
        atac_and_distal_tsr_same_anchor, colocated_atac_distal_tsr,
        n_pooled_pfc_distal_tsr_gt3kb_centers,
        n_pooled_pfc_distal_tsr_gt3kb_centers_in_atac
    ),
    by = "loop_id"
  ) %>%
  left_join(df.current.gene.annotation, by = "gene_id") %>%
  left_join(df.current.symbol.multiplicity, by = "current_gene_symbol") %>%
  left_join(
    df.expression,
    by = c("current_gene_symbol" = "expression_gene_symbol")
  ) %>%
  mutate(
    expression_mapping_status = case_when(
      is.na(current_gene_symbol) ~
        "gene_id_or_symbol_not_confirmed_in_current_gtf",
      n_current_symbols_for_gene_id != 1L ~
        "ambiguous_current_symbol_for_gene_id",
      n_current_ensembl_gene_ids_for_symbol > 1L ~
        "ambiguous_symbol_in_current_gtf",
      is.na(pfc_naive_r1_count) ~ "current_symbol_unmapped_in_expression_table",
      TRUE ~ "exact_current_gene_id_symbol_expression_match"
    ),
    expression_support_evaluable =
      expression_mapping_status ==
        "exact_current_gene_id_symbol_expression_match",
    naive_pfc_expression_support =
      expression_support_evaluable &
        expressed_cpm1_both_naive_replicates %in% TRUE,
    evidence_tier = case_when(
      colocated_atac_distal_tsr &
        naive_pfc_expression_support ~
        "Hi-C+ATAC+pooled-PFC-distal-TSR+naive-PFC-expression",
      colocated_atac_distal_tsr & expression_support_evaluable ~
        "Hi-C+ATAC+pooled-PFC-distal-TSR;gene-not-detected",
      colocated_atac_distal_tsr ~
        "Hi-C+ATAC+pooled-PFC-distal-TSR;expression-unmapped",
      atac_and_distal_tsr_same_anchor ~
        "Hi-C+ATAC+pooled-PFC-distal-TSR;not-colocated",
      atac_overlap_ge50 ~ "Hi-C+ATAC",
      pooled_pfc_distal_tsr_gt3kb_present ~ "Hi-C+pooled-PFC-distal-TSR",
      TRUE ~ "Hi-C-only"
    )
  )

df.loop.expression <- df.loop.gene.evidence %>%
  group_by(loop_id) %>%
  summarise(
    n_connected_genes = n_distinct(gene_id),
    n_expression_mapped_genes = n_distinct(
      gene_id[expression_support_evaluable]
    ),
    n_expressed_connected_genes = n_distinct(
      gene_id[naive_pfc_expression_support]
    ),
    any_connected_gene_expressed = n_expressed_connected_genes > 0L,
    .groups = "drop"
  )

df.loop.evidence <- df.observed %>%
  left_join(df.loop.expression, by = "loop_id") %>%
  mutate(
    any_connected_gene_expressed = coalesce(any_connected_gene_expressed, FALSE),
    complete_four_layer_support =
      colocated_atac_distal_tsr & any_connected_gene_expressed
  )

df.evidence.summary <- bind_rows(
  df.loop.evidence,
  df.loop.evidence %>% mutate(resolution = "ALL")
) %>%
  group_by(resolution) %>%
  summarise(
    n_one_sided_promoter_anchored_loops = n(),
    n_atac_ge50 = sum(atac_overlap_ge50),
    pct_atac_ge50 = 100 * mean(atac_overlap_ge50),
    n_pooled_pfc_distal_tsr = sum(pooled_pfc_distal_tsr_gt3kb_present),
    pct_pooled_pfc_distal_tsr = 100 * mean(pooled_pfc_distal_tsr_gt3kb_present),
    n_atac_and_tsr_same_anchor = sum(atac_and_distal_tsr_same_anchor),
    pct_atac_and_tsr_same_anchor = 100 * mean(atac_and_distal_tsr_same_anchor),
    n_colocated_atac_distal_tsr = sum(colocated_atac_distal_tsr),
    pct_colocated_atac_distal_tsr = 100 * mean(colocated_atac_distal_tsr),
    n_complete_four_layer_support = sum(complete_four_layer_support),
    pct_complete_four_layer_support = 100 * mean(complete_four_layer_support),
    .groups = "drop"
  ) %>%
  arrange(factor(resolution, resolution.levels))

df.gene.evidence <- df.loop.gene.evidence %>%
  group_by(
    gene_id, current_gene_symbol, expression_mapping_status,
    pfc_naive_r1_count, pfc_naive_r2_count,
    pfc_naive_r1_cpm, pfc_naive_r2_cpm, pfc_naive_mean_cpm,
    expressed_cpm1_both_naive_replicates,
    expression_support_evaluable,
    naive_pfc_expression_support
  ) %>%
  summarise(
    n_linked_loops = n_distinct(loop_id),
    any_linked_loop_atac = any(atac_overlap_ge50),
    any_linked_loop_pooled_pfc_distal_tsr =
      any(pooled_pfc_distal_tsr_gt3kb_present),
    any_linked_loop_colocated_atac_distal_tsr = any(colocated_atac_distal_tsr),
    .groups = "drop"
  )

df.expression.mapping.summary <- df.gene.evidence %>%
  summarise(
    n_unique_ensembl_gene_ids = n(),
    n_unique_current_gene_symbols = n_distinct(
      current_gene_symbol[!is.na(current_gene_symbol)]
    ),
    n_exact_current_gene_id_symbol_expression_matches = sum(
      expression_mapping_status ==
        "exact_current_gene_id_symbol_expression_match"
    ),
    n_current_gtf_ambiguous_or_missing = sum(
      expression_mapping_status %in% c(
        "gene_id_or_symbol_not_confirmed_in_current_gtf",
        "ambiguous_current_symbol_for_gene_id",
        "ambiguous_symbol_in_current_gtf"
      )
    ),
    n_unmapped_expression_symbols = sum(
      expression_mapping_status ==
        "current_symbol_unmapped_in_expression_table"
    ),
    pct_exact_current_gene_id_symbol_expression_matches = 100 * mean(
      expression_mapping_status ==
        "exact_current_gene_id_symbol_expression_match"
    ),
    n_expressed_cpm1_both = sum(
      naive_pfc_expression_support
    )
  )

df.expression.by.support <- df.gene.evidence %>%
  filter(expression_support_evaluable) %>%
  mutate(
    distal_evidence_group = case_when(
      any_linked_loop_colocated_atac_distal_tsr ~ "colocated_ATAC_distal_TSR",
      any_linked_loop_atac & any_linked_loop_pooled_pfc_distal_tsr ~
        "ATAC_and_distal_TSR_on_different_or_unknown_sites",
      any_linked_loop_atac ~ "ATAC_only",
      any_linked_loop_pooled_pfc_distal_tsr ~ "distal_TSR_only",
      TRUE ~ "neither"
    )
  ) %>%
  group_by(distal_evidence_group) %>%
  summarise(
    n_mapped_genes = n(),
    n_expressed_cpm1_both = sum(naive_pfc_expression_support),
    pct_expressed_cpm1_both = 100 * mean(naive_pfc_expression_support),
    median_pfc_naive_mean_cpm = median(pfc_naive_mean_cpm),
    .groups = "drop"
  ) %>%
  arrange(desc(n_mapped_genes))

df.tsr.summary <- tibble(
  metric = c(
    "pooled_pfc_tsr_source_total",
    "pooled_pfc_tsr_excluded_nonanalysis_seqlevels",
    "pooled_pfc_tsr_analyzed",
    "pooled_pfc_tsr_outside_unified_tss_plus_minus_1kb",
    "pooled_pfc_distal_tsr_gt3kb",
    "pooled_pfc_distal_tsr_gt3kb_center_in_non_tss_atac"
  ),
  n = c(
    nrow(df.tsr.source),
    nrow(df.tsr.excluded.seqlevels),
    nrow(df.tsr),
    sum(df.tsr$is_outside_tss_plus_minus_1kb),
    sum(df.tsr$is_distal_gt3kb),
    sum(df.tsr$is_distal_gt3kb_and_atac)
  )
)

df.expression.test.input <- df.gene.evidence %>%
  filter(expression_support_evaluable) %>%
  mutate(
    has_colocated_distal_evidence = any_linked_loop_colocated_atac_distal_tsr,
    expression_positive = naive_pfc_expression_support
  )
expression.table <- table(
  df.expression.test.input$has_colocated_distal_evidence,
  df.expression.test.input$expression_positive
)
expression.fisher <- fisher.test(expression.table)
expression.glm <- glm(
  expression_positive ~ has_colocated_distal_evidence + log1p(n_linked_loops),
  data = df.expression.test.input,
  family = binomial()
)
expression.coef <- summary(expression.glm)$coefficients[
  "has_colocated_distal_evidenceTRUE", , drop = FALSE
]
expression.log.or <- unname(expression.coef[1, "Estimate"])
expression.se <- unname(expression.coef[1, "Std. Error"])
df.expression.association <- tibble(
  analysis_unit = "unique Ensembl gene with exact current-symbol expression mapping",
  n_genes = nrow(df.expression.test.input),
  exposed_group = "at least one linked loop with colocated ATAC-distal TSR",
  outcome = "CPM >=1 in both naive PFC replicates",
  fisher_odds_ratio = unname(expression.fisher$estimate),
  fisher_p = expression.fisher$p.value,
  adjusted_logistic_odds_ratio = exp(expression.log.or),
  adjusted_logistic_ci_low = exp(expression.log.or - 1.96 * expression.se),
  adjusted_logistic_ci_high = exp(expression.log.or + 1.96 * expression.se),
  adjusted_logistic_p = unname(expression.coef[1, "Pr(>|z|)"]),
  logistic_covariate = "log1p(number of linked loops)"
)

df.complete.candidates <- df.loop.gene.evidence %>%
  filter(colocated_atac_distal_tsr, naive_pfc_expression_support) %>%
  arrange(
    factor(resolution, c("5K", "10K", "25K")),
    desc(n_supporting_strains),
    desc(n_supporting_libraries),
    desc(n_pooled_pfc_distal_tsr_gt3kb_centers_in_atac),
    desc(non_tss_atac_max_overlap_bp),
    desc(pfc_naive_mean_cpm),
    loop_id,
    gene_id
  )

df.parameters <- tribble(
  ~parameter, ~value, ~interpretation,
  "analysis_unit", "one-sided promoter-anchored pooled Hi-C call", "Direct promoter/TSS support at exactly one anchor; that anchor may contain more than one assigned gene",
  "n_primary_loops", as.character(nrow(df.observed)), "Eligibility fixed before ATAC or TSR annotation",
  "genome_assembly", "mRatBN7.2/rn7", "Final coordinate system for overlap analysis",
  "promoter_definition", "unified Ensembl/EPD TSS +/-1000 bp", "Inherited from the revised primary loop analysis",
  "ATAC_support", ">=50 bp non-TSS ATAC overlap", "Open-chromatin support only",
  "TSR_source", "pooled PFC META-experiment", "Not naive-specific; includes all PFC treatment conditions",
  "TSR_anchor_overlap", "TSR center inside distal anchor", "Avoids counting a TSR that only clips an anchor edge",
  "distal_TSR_definition", ">3000 bp from every unified TSS", "Conservative enhancer-like transcription-initiation evidence",
  "primary_joint_support", "distal TSR center inside non-TSS ATAC and distal anchor", "Requires feature-level ATAC/TSR colocalization",
  "expression_support", "CPM >=1 in both naive PFC replicates", "Tissue expression support; not regulatory causality",
  "null_model", "matched no-direct-TSS Hi-C anchors", "Controls sampled with replacement after matching by chromosome, resolution, anchor side, and loop-distance bin",
  "n_resamples", as.character(n.resamples), "Monte Carlo matched-control resampling for a one-sided empirical tail probability"
)

write_tsv(df.loop.evidence, file.path(output.dir, "distal_anchor_evidence.tsv.gz"))
write_tsv(df.loop.gene.evidence, file.path(output.dir, "loop_gene_evidence.tsv.gz"))
write_tsv(
  df.complete.candidates,
  file.path(output.dir, "complete_four_layer_loop_gene_candidates.tsv.gz")
)
write_tsv(df.gene.evidence, file.path(output.dir, "gene_evidence.tsv.gz"))
write_tsv(df.tsr, file.path(output.dir, "pooled_pfc_tsr_classification.tsv.gz"))
write_tsv(df.control.pool, file.path(output.dir, "matched_control_anchor_pool.tsv.gz"))
write_tsv(df.resample, file.path(output.dir, "matched_null_resamples.tsv.gz"))
write_tsv(df.null.summary, file.path(output.dir, "matched_null_summary.tsv"))
write_tsv(
  df.match.availability.summary,
  file.path(output.dir, "matched_control_availability.tsv")
)
write_tsv(df.evidence.summary, file.path(output.dir, "evidence_summary_by_resolution.tsv"))
write_tsv(df.expression.mapping.summary, file.path(output.dir, "expression_mapping_summary.tsv"))
write_tsv(
  df.gene.evidence %>% filter(!expression_support_evaluable),
  file.path(output.dir, "expression_mapping_exclusions.tsv")
)
write_tsv(df.expression.by.support, file.path(output.dir, "expression_by_distal_support.tsv"))
write_tsv(
  df.expression.association,
  file.path(output.dir, "expression_association_test.tsv")
)
write_tsv(df.tsr.summary, file.path(output.dir, "pooled_pfc_tsr_summary.tsv"))
write_tsv(df.hic.support.balance, file.path(output.dir, "hic_support_balance.tsv"))
write_tsv(df.parameters, file.path(output.dir, "analysis_parameters.tsv"))

plot.summary <- df.evidence.summary %>%
  select(
    resolution,
    `Hi-C + ATAC` = pct_atac_ge50,
    `Hi-C + ATAC + colocated distal TSR` = pct_colocated_atac_distal_tsr,
    `Four-layer support` = pct_complete_four_layer_support
  ) %>%
  pivot_longer(-resolution, names_to = "evidence", values_to = "percent") %>%
  mutate(
    resolution = factor(resolution, resolution.levels),
    evidence = factor(
      evidence,
      c("Hi-C + ATAC", "Hi-C + ATAC + colocated distal TSR", "Four-layer support")
    )
  )

p1 <- ggplot(plot.summary, aes(resolution, percent, fill = evidence)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72) +
  scale_fill_manual(values = c("#3677A8", "#D28E35", "#4A8B63")) +
  labs(
    x = "Hi-C loop resolution",
    y = "One-sided promoter-anchored contacts (%)",
    fill = NULL,
    title = "External PFC evidence at Hi-C distal anchors"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

metric.labels <- c(
  atac_overlap_ge50 = "ATAC >=50 bp",
  pooled_pfc_tsr_present = "Any pooled-PFC TSR",
  pooled_pfc_distal_tsr_gt3kb_present = "Distal TSR (>3 kb from TSS)",
  atac_and_distal_tsr_same_anchor = "ATAC and distal TSR in anchor",
  colocated_atac_distal_tsr = "Colocated ATAC-distal TSR"
)
plot.null <- df.null.summary %>%
  mutate(
    resolution = factor(resolution, resolution.levels),
    metric_label = factor(metric.labels[metric], levels = unname(metric.labels))
  ) %>%
  select(resolution, metric_label, observed_rate, mean_null_rate) %>%
  pivot_longer(
    c(observed_rate, mean_null_rate), names_to = "group", values_to = "rate"
  )

p2 <- ggplot(plot.null, aes(rate * 100, metric_label, color = group)) +
  geom_line(aes(group = interaction(resolution, metric_label)), color = "grey75") +
  geom_point(size = 2.2) +
  facet_wrap(~resolution, ncol = 2) +
  scale_color_manual(
    values = c(observed_rate = "#1E5B85", mean_null_rate = "#B55B3C"),
    labels = c(observed_rate = "Observed distal anchors", mean_null_rate = "Matched Hi-C anchors")
  ) +
  labs(x = "Positive anchors (%)", y = NULL, color = NULL, title = "Matched-anchor enrichment") +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom")

ggsave(
  file.path(output.dir, "evidence_tiers_by_resolution.png"), p1,
  width = 8, height = 5, dpi = 180
)
ggsave(
  file.path(output.dir, "evidence_tiers_by_resolution.pdf"), p1,
  width = 8, height = 5
)
ggsave(
  file.path(output.dir, "matched_anchor_enrichment.png"), p2,
  width = 10, height = 7, dpi = 180
)
ggsave(
  file.path(output.dir, "matched_anchor_enrichment.pdf"), p2,
  width = 10, height = 7
)

writeLines(capture.output(sessionInfo()), file.path(output.dir, "session_info.txt"))

message("Completed distal-anchor evidence integration.")
print(df.evidence.summary)
print(df.null.summary)
print(df.expression.mapping.summary)
