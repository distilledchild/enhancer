# lintr: disable

# Compare full-depth, 140M-contact, and 250M-contact Hi-C loop calls using the
# same seven libraries that have enough valid contacts for both downsamplings.
# The original ten-library 140M analysis remains in downsampling_140M.R.

current_script_path <- function() {
  file.args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file.args) > 0L) {
    return(normalizePath(
      sub("^--file=", "", file.args[[1]]),
      winslash = "/", mustWork = FALSE
    ))
  }
  frame.files <- vapply(
    sys.frames(),
    function(frame) if (is.null(frame$ofile)) NA_character_ else frame$ofile,
    character(1L)
  )
  frame.files <- frame.files[!is.na(frame.files) & nzchar(frame.files)]
  if (length(frame.files) > 0L) {
    return(normalizePath(tail(frame.files, 1L), winslash = "/", mustWork = FALSE))
  }
  NA_character_
}

# Locate the shared R directory from the script or the established local and
# Dropbox layouts so collaborators can run the same file after synchronization.
resolve_enhancer_r_files_dir <- function() {
  script.path <- current_script_path()
  candidates <- unique(c(
    Sys.getenv("ENHANCER_R_FILES_DIR", unset = ""),
    if (is.na(script.path)) NA_character_ else dirname(script.path),
    path.expand("~/dropbox/Gateway_to_Hao/enhancer/r_files"),
    path.expand("~/Dropbox/Gateway_to_Hao/enhancer/r_files"),
    Sys.glob(path.expand(
      "~/Library/CloudStorage/Dropbox*/K P/Gateway_to_Hao/enhancer/r_files"
    )),
    path.expand("~/Desktop/playground/enhancer/r_files")
  ))
  candidates <- candidates[
    !is.na(candidates) & nzchar(candidates) & dir.exists(candidates)
  ]
  candidates <- candidates[
    file.exists(file.path(candidates, "funcs.R")) &
      file.exists(file.path(candidates, "downsampling_depth_functions.R"))
  ]
  if (length(candidates) == 0L) {
    stop(
      "Cannot locate funcs.R and downsampling_depth_functions.R. ",
      "Set ENHANCER_R_FILES_DIR.",
      call. = FALSE
    )
  }
  normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE)
}

r.files.dir <- resolve_enhancer_r_files_dir()
source(file.path(r.files.dir, "funcs.R"))

library("tidyverse")
library("GenomicRanges")
library("GenomeInfoDb")

source(file.path(r.files.dir, "downsampling_depth_functions.R"))

options(tibble.width = Inf, tibble.print_max = Inf, scipen = 999)

################################################################################
# Full-depth versus 140M versus 250M sensitivity analysis
#
# 1. Restrict every condition to the same seven libraries.
# 2. Rebuild exact, resolution-specific pooled HiCCUPS calls shorter than 2 Mb.
# 3. Compare loop counts, depth dependence, exact recovery, and radius recovery.
# 4. Build approximate cross-resolution loci only as a sensitivity analysis.
# 5. Compare library sharing, regulatory categories, and gene-loop counts.
################################################################################

########################
# 0. Directories, sample metadata, and constants
########################

revision.dir <- file.path(r.files.dir, "revision")
max.loop.distance.bp <- 2000000L
promoter.window.flank.bp <- 1000L
atac.minimum.overlap.bp <- 50L
common.library.count <- 7L

# These are the seven libraries with at least 250M valid MAPQ >=30 contacts.
df.sample.metadata <- tribble(
  ~sample, ~strain, ~source_contacts, ~seed,
  "592BB", "SHR/OlaIpcv", 259559222, 20260731,
  "607", "HXB10", 419977942, 20260737,
  "D765A", "BXH6", 284296675, 20260732,
  "DA21A", "SHR/OlaIpcvxBN/NHsdMcwi", 382917950, 20260736,
  "DA68A", "HXB31", 378608379, 20260734,
  "DBA9A", "HXB23", 380788024, 20260735,
  "DE8BA", "BN-Lx", 302314076, 20260733
)

full.depth.root <- resolve_depth_directory(
  c(
    Sys.getenv("FULL_DEPTH_HICCUPS_ROOT", unset = ""),
    path.expand("~/dropbox/Gateway_to_Hao/hic/2023A/hic30_w_sb_options"),
    Sys.glob(path.expand(paste0(
      "~/Library/CloudStorage/Dropbox*/K P/Gateway_to_Hao/",
      "hic/2023A/hic30_w_sb_options"
    )))
  ),
  "full-depth sb-option HiCCUPS directory"
)

combined.downsample.root <- resolve_depth_directory(
  c(
    Sys.getenv("DOWNSAMPLED_140M_250M_HICCUPS_ROOT", unset = ""),
    path.expand(paste0(
      "~/Library/CloudStorage/GoogleDrive-wellclouder@gmail.com/My Drive/",
      "juicer_downsample_q30_140M_250M"
    ))
  ),
  "combined 140M/250M HiCCUPS directory"
)

downsample.140m.root <- resolve_depth_directory(
  c(
    Sys.getenv("DOWNSAMPLED_140M_HICCUPS_ROOT", unset = ""),
    file.path(combined.downsample.root, "140M")
  ),
  "140M HiCCUPS directory"
)
downsample.250m.root <- resolve_depth_directory(
  c(
    Sys.getenv("DOWNSAMPLED_250M_HICCUPS_ROOT", unset = ""),
    file.path(combined.downsample.root, "250M")
  ),
  "250M HiCCUPS directory"
)

required.cache.files <- c(
  transcript = "df.ensembl.transcript.coordinate.normalized.rds",
  promoter = "df.promoter.annotation.coordinate.normalized.rds",
  atac = "gr.atac.rds"
)
cache.candidates <- unique(path.expand(c(
  Sys.getenv("DOWNSAMPLING_COORD_CACHE_DIR", unset = ""),
  path.expand("~/dropbox/Gateway_to_Hao/enhancer/r_files/revision/cache_data"),
  file.path(revision.dir, "cache_data"),
  Sys.glob(path.expand(paste0(
    "~/Library/CloudStorage/Dropbox*/K P/Gateway_to_Hao/enhancer/",
    "r_files/revision/cache_data"
  )))
)))
cache.candidates <- cache.candidates[
  nzchar(cache.candidates) & dir.exists(cache.candidates)
]
cache.complete <- vapply(
  cache.candidates,
  function(path) all(file.exists(file.path(path, required.cache.files))),
  logical(1L)
)
if (!any(cache.complete)) {
  stop(
    "Cannot locate a coordinate cache containing transcript, promoter, and ",
    "ATAC objects. Set DOWNSAMPLING_COORD_CACHE_DIR.",
    call. = FALSE
  )
}
coord.cache.dir <- normalizePath(
  cache.candidates[cache.complete][[1]], winslash = "/", mustWork = TRUE
)

default.output.parent <- c(
  path.expand("~/dropbox/Gateway_to_Hao/enhancer/r_files"),
  r.files.dir
)
default.output.parent <- default.output.parent[
  dir.exists(default.output.parent)
][[1]]
output.dir <- path.expand(Sys.getenv(
  "DOWNSAMPLING_140M_250M_OUTPUT_DIR",
  unset = file.path(default.output.parent, "downsampling_140M_250M_outputs")
))
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

message("Full-depth HiCCUPS root: ", full.depth.root)
message("140M HiCCUPS root: ", downsample.140m.root)
message("250M HiCCUPS root: ", downsample.250m.root)
message("Coordinate cache: ", coord.cache.dir)
message("Output directory: ", output.dir)

########################
# 1. Input completeness and matched seven-library design
########################

df.status.140m <- build_depth_file_status(
  df.sample.metadata, full.depth.root, downsample.140m.root, "140M"
)
df.status.250m <- build_depth_file_status(
  df.sample.metadata, full.depth.root, downsample.250m.root, "250M"
)
df.input.status <- bind_rows(df.status.140m, df.status.250m)

if (nrow(df.sample.metadata) != common.library.count ||
    any(!df.input.status$analysis_ready) ||
    any(!df.input.status$provenance_ready)) {
  incomplete <- df.input.status %>%
    filter(!analysis_ready | !provenance_ready) %>%
    dplyr::select(target_label, sample, status, everything())
  print(incomplete)
  stop(
    "The common-seven full-depth/140M/250M design is incomplete.",
    call. = FALSE
  )
}

# Verify that both targets came from the same source-contact universe and seed.
df.downsampling.qc <- bind_rows(
  read_depth_downsampling_qc(df.status.140m, 140000000),
  read_depth_downsampling_qc(df.status.250m, 250000000)
)
if (any(!df.downsampling.qc$qc_matches_metadata)) {
  stop("Downsampling QC does not match the comparison metadata.", call. = FALSE)
}

df.nested.downsampling.design <- df.downsampling.qc %>%
  dplyr::select(
    sample, target_label, seed, source_contacts, target_contacts,
    qc_matches_metadata
  ) %>%
  pivot_wider(
    names_from = target_label,
    values_from = c(seed, source_contacts, target_contacts, qc_matches_metadata)
  ) %>%
  mutate(
    same_source_contact_universe = source_contacts_140M == source_contacts_250M,
    same_random_seed = seed_140M == seed_250M,
    intended_nested_140M_within_250M = (
      same_source_contact_universe & same_random_seed &
        target_contacts_140M < target_contacts_250M
    ),
    contact_level_subset_verified_locally = FALSE,
    design_note = paste0(
      "The same sequential sampler, source record order, and sample-specific ",
      "seed were used at both targets. This is an intended nested ",
      "common-random-number design; raw contact-level subset membership was ",
      "not rechecked on the local computer."
    )
  ) %>%
  arrange(sample)
if (any(!df.nested.downsampling.design$intended_nested_140M_within_250M)) {
  stop("The intended nested 140M/250M design is inconsistent.", call. = FALSE)
}

df.full.files <- df.status.140m %>%
  transmute(sample, strain, file = full_file)
df.140m.files <- df.status.140m %>%
  transmute(sample, strain, file = downsample_file)
df.250m.files <- df.status.250m %>%
  transmute(sample, strain, file = downsample_file)

########################
# 2. Exact pooled loop resources under the same seven-library design
########################

message("Building matched exact loop resources...")
resources <- list(
  full_depth = build_depth_loop_resource(
    df.full.files, "full_depth", max.loop.distance.bp
  ),
  downsample_140M = build_depth_loop_resource(
    df.140m.files, "downsample_140M", max.loop.distance.bp
  ),
  downsample_250M = build_depth_loop_resource(
    df.250m.files, "downsample_250M", max.loop.distance.bp
  )
)

df.sample.loop.membership <- imap_dfr(
  resources,
  ~ .x$sample %>%
    filter(passes_lt2mb) %>%
    transmute(
      condition = .y, sample, strain, loop_id, resolution,
      chr1, start1, end1, chr2, start2, end2, loop_distance
    )
)
df.exact.pool.summary <- imap_dfr(
  resources,
  ~ tibble(
    condition = .y,
    n_sample_level_calls = sum(.x$sample$passes_lt2mb),
    n_exact_pooled_calls = nrow(.x$universe),
    n_private_exact_pooled_calls = sum(.x$universe$n_supporting_libraries == 1L),
    n_recurrent_exact_pooled_calls = sum(.x$universe$n_supporting_libraries >= 2L)
  )
)

########################
# 3. Loop-count variation and residual source-depth association
########################

message("Summarizing loop-count variation...")

df.loop.count.by.resolution <- df.sample.loop.membership %>%
  count(condition, sample, strain, resolution, name = "n_loop_calls")
df.loop.count <- bind_rows(
  df.loop.count.by.resolution,
  df.loop.count.by.resolution %>%
    group_by(condition, sample, strain) %>%
    summarise(resolution = "ALL", n_loop_calls = sum(n_loop_calls), .groups = "drop")
) %>%
  left_join(df.sample.metadata, by = c("sample", "strain")) %>%
  mutate(
    analysis_contacts = case_when(
      condition == "full_depth" ~ source_contacts,
      condition == "downsample_140M" ~ 140000000,
      condition == "downsample_250M" ~ 250000000
    ),
    sampling_fraction = analysis_contacts / source_contacts
  ) %>%
  arrange(resolution, sample, condition)

df.loop.count.summary <- df.loop.count %>%
  group_by(condition, resolution) %>%
  summarise(
    n_libraries = n(),
    mean_loop_calls = mean(n_loop_calls),
    sd_loop_calls = sd(n_loop_calls),
    cv_loop_calls = coefficient_of_variation(n_loop_calls),
    minimum_loop_calls = min(n_loop_calls),
    maximum_loop_calls = max(n_loop_calls),
    n_spearman = safe_spearman_test(
      source_contacts, n_loop_calls
    )$n_spearman,
    spearman_original_source_contacts_vs_loop_calls = safe_spearman(
      source_contacts, n_loop_calls
    ),
    spearman_p_value = safe_spearman_test(
      source_contacts, n_loop_calls
    )$spearman_p_value,
    correlation_interpretation = paste0(
      "Descriptive association across seven libraries; not evidence of an ",
      "independent strain effect."
    ),
    .groups = "drop"
  ) %>%
  arrange(resolution, condition)

########################
# 4. Exact and HiCCUPS-radius recovery between depth conditions
########################

message("Calculating exact and radius-based recovery...")

condition.comparisons <- tribble(
  ~reference_condition, ~query_condition,
  "full_depth", "downsample_140M",
  "full_depth", "downsample_250M",
  "downsample_140M", "downsample_250M"
)

df.exact.recovery <- pmap_dfr(
  condition.comparisons,
  ~ summarise_exact_depth_recovery(
    df.sample.loop.membership, ..1, ..2
  )
)
df.radius.recovery <- pmap_dfr(
  condition.comparisons,
  function(reference_condition, query_condition) {
    summarise_radius_depth_recovery(
      resources[[reference_condition]], resources[[query_condition]],
      reference_condition, query_condition
    )
  }
)

# Compare each pooled exact-call set directly, independent of sample-level
# recovery, so resource-level retention is explicit for every depth pair.
df.pooled.exact.recovery <- pmap_dfr(
  condition.comparisons,
  function(reference_condition, query_condition) {
    reference.ids <- resources[[reference_condition]]$universe$loop_id
    query.ids <- resources[[query_condition]]$universe$loop_id
    shared.ids <- intersect(reference.ids, query.ids)
    tibble(
      reference_condition,
      query_condition,
      n_reference_pooled_exact_calls = length(reference.ids),
      n_query_pooled_exact_calls = length(query.ids),
      n_shared_pooled_exact_calls = length(shared.ids),
      pct_reference_pooled_exact_recovered = 100 * length(shared.ids) /
        length(reference.ids),
      pooled_exact_jaccard = length(shared.ids) /
        length(union(reference.ids, query.ids)),
      pooled_exact_overlap_coefficient = length(shared.ids) /
        min(length(reference.ids), length(query.ids))
    )
  }
)

# Quantify the expected dependence of full-depth recovery on the fraction of
# source contacts retained, rather than treating recovery as replication.
df.sampling.fraction.recovery <- df.exact.recovery %>%
  filter(reference_condition == "full_depth", resolution == "ALL") %>%
  left_join(
    df.radius.recovery %>%
      filter(reference_condition == "full_depth", resolution == "ALL") %>%
      dplyr::select(
        sample, query_condition, pct_reference_radius_recovered
      ),
    by = c("sample", "query_condition")
  ) %>%
  left_join(df.sample.metadata, by = "sample") %>%
  mutate(
    target_contacts = if_else(
      query_condition == "downsample_140M", 140000000, 250000000
    ),
    sampling_fraction = target_contacts / source_contacts
  ) %>%
  arrange(query_condition, sampling_fraction, sample)

df.sampling.fraction.recovery.summary <- df.sampling.fraction.recovery %>%
  group_by(query_condition) %>%
  summarise(
    n_libraries = n(),
    spearman_sampling_fraction_vs_exact_recovery = safe_spearman(
      sampling_fraction, pct_reference_recovered
    ),
    spearman_exact_recovery_p_value = safe_spearman_test(
      sampling_fraction, pct_reference_recovered
    )$spearman_p_value,
    spearman_sampling_fraction_vs_radius_recovery = safe_spearman(
      sampling_fraction, pct_reference_radius_recovered
    ),
    spearman_radius_recovery_p_value = safe_spearman_test(
      sampling_fraction, pct_reference_radius_recovered
    )$spearman_p_value,
    .groups = "drop"
  )

########################
# 5. Approximate cross-resolution loop loci as sensitivity analysis
########################

message("Building joint approximate loci for sensitivity only...")

# Prefix condition-specific source labels so the locus builder may align the
# same biological library across depth conditions without merging two distinct
# calls from the same library within one condition.
prefix_supporting_samples <- function(x, condition) {
  map_chr(str_split(x, fixed(";")), function(samples) {
    str_c(str_c(condition, samples, sep = "::"), collapse = ";")
  })
}

df.approximate.input <- imap_dfr(
  resources,
  function(resource, condition) {
    resource$universe %>%
      mutate(
        condition = .env$condition,
        original_loop_id = loop_id,
        loop_id = str_c(.env$condition, original_loop_id, sep = "::"),
        supporting_samples = prefix_supporting_samples(
          supporting_samples, .env$condition
        )
      )
  }
)
approximate.resource <- build_approximate_loop_loci(df.approximate.input)
df.approximate.locus.summary <- approximate.resource$summary %>%
  dplyr::rename(
    n_supporting_condition_library_instances = n_supporting_libraries,
    supporting_condition_library_instances = supporting_samples
  )

df.approximate.key <- df.approximate.input %>%
  transmute(
    condition,
    original_loop_id,
    synthetic_loop_id = loop_id
  ) %>%
  left_join(
    approximate.resource$map %>%
      transmute(
        synthetic_loop_id = loop_id,
        approximate_loop_locus_id,
        n_exact_loop_calls_in_locus,
        n_resolutions_in_locus,
        is_multi_call_locus
      ),
    by = "synthetic_loop_id"
  )

if (any(is.na(df.approximate.key$approximate_loop_locus_id))) {
  stop("An exact call lacks an approximate-locus assignment.", call. = FALSE)
}

# Report position-tolerant pooled-set overlap separately from the exact-call
# comparison; these joint loci remain a sensitivity analysis only.
df.pooled.approximate.recovery <- pmap_dfr(
  condition.comparisons,
  function(reference_condition, query_condition) {
    reference.ids <- df.approximate.key %>%
      filter(condition == reference_condition) %>%
      pull(approximate_loop_locus_id) %>%
      unique()
    query.ids <- df.approximate.key %>%
      filter(condition == query_condition) %>%
      pull(approximate_loop_locus_id) %>%
      unique()
    shared.ids <- intersect(reference.ids, query.ids)
    tibble(
      reference_condition,
      query_condition,
      n_reference_approximate_loci = length(reference.ids),
      n_query_approximate_loci = length(query.ids),
      n_shared_approximate_loci = length(shared.ids),
      pct_reference_approximate_recovered = 100 * length(shared.ids) /
        length(reference.ids),
      approximate_locus_jaccard = length(shared.ids) /
        length(union(reference.ids, query.ids)),
      approximate_locus_overlap_coefficient = length(shared.ids) /
        min(length(reference.ids), length(query.ids))
    )
  }
)

df.approximate.membership <- df.sample.loop.membership %>%
  left_join(
    df.approximate.key %>%
      dplyr::select(condition, original_loop_id, approximate_loop_locus_id),
    by = c("condition", "loop_id" = "original_loop_id")
  ) %>%
  transmute(
    condition, sample, strain,
    resolution = "ALL_CROSS_RESOLUTION",
    approximate_loop_locus_id
  ) %>%
  distinct()

########################
# 6. Pairwise library sharing, support distributions, and recurrent cores
########################

message("Summarizing library sharing and recurrent support...")

df.exact.pairwise <- summarise_pairwise_library_similarity(
  add_all_resolution(df.sample.loop.membership),
  "loop_id", "exact_resolution_specific_call"
)
df.approximate.pairwise <- summarise_pairwise_library_similarity(
  df.approximate.membership,
  "approximate_loop_locus_id", "approximate_cross_resolution_locus"
)

df.exact.support <- df.sample.loop.membership %>%
  add_all_resolution() %>%
  distinct(condition, resolution, sample, loop_id) %>%
  count(condition, resolution, loop_id, name = "n_supporting_libraries") %>%
  count(condition, resolution, n_supporting_libraries, name = "n_features") %>%
  group_by(condition, resolution) %>%
  mutate(
    feature_type = "exact_resolution_specific_call",
    pct_features = 100 * n_features / sum(n_features)
  ) %>%
  ungroup()
df.approximate.support <- df.approximate.membership %>%
  distinct(condition, sample, approximate_loop_locus_id) %>%
  count(
    condition, approximate_loop_locus_id,
    name = "n_supporting_libraries"
  ) %>%
  count(condition, n_supporting_libraries, name = "n_features") %>%
  group_by(condition) %>%
  mutate(
    resolution = "ALL_CROSS_RESOLUTION",
    feature_type = "approximate_cross_resolution_locus",
    pct_features = 100 * n_features / sum(n_features)
  ) %>%
  ungroup()

df.recurrent.core.stability <- bind_rows(
  summarise_recurrent_core_stability(
    df.sample.loop.membership,
    "loop_id", "exact_resolution_specific_call",
    maximum.support = common.library.count
  ),
  summarise_recurrent_core_stability(
    df.approximate.membership,
    "approximate_loop_locus_id", "approximate_cross_resolution_locus",
    maximum.support = common.library.count
  )
)

# Preserve every feature and both support counts for the primary >=2-library
# recurrent-core sensitivity table rather than exporting only summary metrics.
build_two_depth_core_membership <- function(
  df.membership, feature.column, feature.type, minimum.support = 2L
) {
  support <- df.membership %>%
    filter(condition %in% c("downsample_140M", "downsample_250M")) %>%
    distinct(condition, sample, .data[[feature.column]]) %>%
    count(condition, .data[[feature.column]], name = "n_supporting_libraries") %>%
    pivot_wider(
      names_from = condition,
      values_from = n_supporting_libraries,
      values_fill = 0L,
      names_prefix = "n_libraries_"
    )
  support %>%
    mutate(
      feature_type = feature.type,
      minimum_library_support = minimum.support,
      recurrent_140M = n_libraries_downsample_140M >= minimum.support,
      recurrent_250M = n_libraries_downsample_250M >= minimum.support,
      recurrent_at_both_depths = recurrent_140M & recurrent_250M
    )
}

df.recurrent.core.membership <- bind_rows(
  build_two_depth_core_membership(
    df.sample.loop.membership,
    "loop_id", "exact_resolution_specific_call"
  ),
  build_two_depth_core_membership(
    df.approximate.membership,
    "approximate_loop_locus_id", "approximate_cross_resolution_locus"
  )
)

########################
# 7. Influence of individual libraries on each pooled resource
########################

message("Summarizing per-library pooled-resource influence...")

df.pool.library.influence <- df.sample.loop.membership %>%
  distinct(condition, sample, loop_id) %>%
  add_count(condition, loop_id, name = "n_supporting_libraries") %>%
  group_by(condition, sample) %>%
  summarise(
    n_sample_loop_calls = n_distinct(loop_id),
    n_sample_private_calls = n_distinct(loop_id[n_supporting_libraries == 1L]),
    .groups = "drop"
  ) %>%
  left_join(
    df.sample.loop.membership %>%
      distinct(condition, loop_id) %>%
      count(condition, name = "n_pooled_exact_calls"),
    by = "condition"
  ) %>%
  left_join(df.sample.metadata, by = "sample") %>%
  mutate(
    pct_pool_called_in_library = 100 * n_sample_loop_calls / n_pooled_exact_calls,
    pct_pool_private_to_library = 100 * n_sample_private_calls /
      n_pooled_exact_calls
  ) %>%
  arrange(condition, desc(source_contacts))

########################
# 8. Regulatory-category and gene-loop-count stability
########################

required.cache.paths <- setNames(
  file.path(coord.cache.dir, required.cache.files),
  names(required.cache.files)
)
if (any(!file.exists(required.cache.paths))) {
  stop(
    "Missing regulatory coordinate cache(s):\n",
    paste(required.cache.paths[!file.exists(required.cache.paths)], collapse = "\n"),
    call. = FALSE
  )
}

df.transcript <- readRDS(required.cache.paths[["transcript"]])
df.promoter <- readRDS(required.cache.paths[["promoter"]])
gr.atac <- readRDS(required.cache.paths[["atac"]])
regulatory.context <- build_depth_regulatory_context(
  df.transcript, df.promoter, gr.atac, promoter.window.flank.bp
)

message("Annotating regulatory evidence at all three depths...")
annotations <- imap(
  resources,
  ~ annotate_depth_loop_resource(
    .x$universe,
    .y,
    regulatory.context,
    df.promoter,
    promoter.window.flank.bp,
    atac.minimum.overlap.bp
  )
)

df.loop.annotation <- map_dfr(annotations, "loop")
df.gene.count <- map_dfr(annotations, "gene_count")
df.category.count <- df.loop.annotation %>%
  add_all_resolution() %>%
  count(
    condition, resolution, revised_major_category,
    dual_promoter_directional_category,
    name = "n_exact_loop_calls"
  ) %>%
  group_by(condition, resolution) %>%
  mutate(pct_exact_loop_calls = 100 * n_exact_loop_calls / sum(n_exact_loop_calls)) %>%
  ungroup()

df.category.concordance <- pmap_dfr(
  condition.comparisons,
  function(reference_condition, query_condition) {
    inner_join(
      df.loop.annotation %>%
        filter(condition == reference_condition) %>%
        dplyr::select(
          loop_id,
          category_reference = revised_major_category,
          dual_category_reference = dual_promoter_directional_category
        ),
      df.loop.annotation %>%
        filter(condition == query_condition) %>%
        dplyr::select(
          loop_id,
          category_query = revised_major_category,
          dual_category_query = dual_promoter_directional_category
        ),
      by = "loop_id"
    ) %>%
      count(
        category_reference, category_query,
        dual_category_reference, dual_category_query,
        name = "n_exact_shared_calls"
      ) %>%
      mutate(
        reference_condition = reference_condition,
        query_condition = query_condition,
        .before = 1
      )
  }
)

# Compare category composition directly, including calls that are not shared
# exactly between the two depth conditions.
category.composition <- pmap(
  condition.comparisons,
  function(reference_condition, query_condition) {
    detail <- full_join(
      df.category.count %>%
        filter(condition == reference_condition, resolution == "ALL") %>%
        transmute(
          revised_major_category,
          dual_promoter_directional_category,
          n_reference = n_exact_loop_calls,
          pct_reference = pct_exact_loop_calls
        ),
      df.category.count %>%
        filter(condition == query_condition, resolution == "ALL") %>%
        transmute(
          revised_major_category,
          dual_promoter_directional_category,
          n_query = n_exact_loop_calls,
          pct_query = pct_exact_loop_calls
        ),
      by = c(
        "revised_major_category",
        "dual_promoter_directional_category"
      )
    ) %>%
      mutate(
        reference_condition = reference_condition,
        query_condition = query_condition,
        n_reference = coalesce(n_reference, 0L),
        n_query = coalesce(n_query, 0L),
        pct_reference = coalesce(pct_reference, 0),
        pct_query = coalesce(pct_query, 0),
        percentage_point_difference = pct_query - pct_reference,
        absolute_percentage_point_difference = abs(
          percentage_point_difference
        ),
        .before = 1
      )
    list(
      detail = detail,
      summary = tibble(
        reference_condition,
        query_condition,
        total_variation_percentage_points =
          0.5 * sum(detail$absolute_percentage_point_difference),
        maximum_category_percentage_point_difference =
          max(detail$absolute_percentage_point_difference)
      )
    )
  }
)
df.category.composition.stability <- map_dfr(category.composition, "detail")
df.category.composition.summary <- map_dfr(category.composition, "summary")

df.gene.assignment <- map_dfr(annotations, "gene_assignment")
df.gene.approximate.count <- df.gene.assignment %>%
  filter(revised_putative_regulatory_support, !is.na(ensembl_gene_id)) %>%
  left_join(
    df.approximate.key %>%
      dplyr::select(
        condition,
        original_loop_id,
        approximate_loop_locus_id
      ),
    by = c("condition", "loop_id" = "original_loop_id")
  )
if (any(is.na(df.gene.approximate.count$approximate_loop_locus_id))) {
  stop(
    "A regulatory gene assignment lacks an approximate-locus mapping.",
    call. = FALSE
  )
}
df.gene.approximate.count <- df.gene.approximate.count %>%
  distinct(
    condition, approximate_loop_locus_id, ensembl_gene_id,
    .keep_all = TRUE
  ) %>%
  group_by(condition, ensembl_gene_id) %>%
  summarise(
    gene_symbol = {
      values <- sort(unique(na.omit(gene_symbol)))
      if (length(values) == 0L) NA_character_ else values[[1]]
    },
    n_approximate_loop_loci = n_distinct(approximate_loop_locus_id),
    .groups = "drop"
  ) %>%
  arrange(desc(n_approximate_loop_loci), gene_symbol, ensembl_gene_id)

gene.stability.exact <- pmap(
  condition.comparisons,
  ~ summarise_gene_depth_stability(
    df.gene.count, ..1, ..2,
    count.column = "n_exact_loop_calls",
    gene.count.unit = "exact_resolution_specific_call"
  )
)
gene.stability.approximate <- pmap(
  condition.comparisons,
  ~ summarise_gene_depth_stability(
    df.gene.approximate.count, ..1, ..2,
    count.column = "n_approximate_loop_loci",
    gene.count.unit = "approximate_cross_resolution_locus"
  )
)
df.gene.stability.summary <- bind_rows(
  map_dfr(gene.stability.exact, "summary"),
  map_dfr(gene.stability.approximate, "summary")
)
df.gene.stability.detail <- bind_rows(
  map_dfr(gene.stability.exact, "detail"),
  map_dfr(gene.stability.approximate, "detail")
)

########################
# 9. Outputs, figures, and reproducibility metadata
########################

df.method.definition <- tribble(
  ~analysis_component, ~definition,
  "comparison_cohort",
  paste0(
    "The same seven libraries with >=250M valid MAPQ >=30 contacts were used ",
    "at full depth, 140M, and 250M."
  ),
  "downsampling_unit",
  paste0(
    "Duplicate-removed Juicer merged_nodups contacts with MAPQ >=30 at both ",
    "ends were sampled after excluding intra-fragment pairs."
  ),
  "nested_target_design",
  paste0(
    "For each common library, 140M and 250M used the same source file, ",
    "sample-specific seed, and sequential sampler. The lower target is an ",
    "intended nested subset of the higher target; only one seed realization ",
    "was analyzed."
  ),
  "main_loop_unit",
  paste0(
    "Exact coordinate- and resolution-specific pooled HiCCUPS call records ",
    "shorter than 2 Mb; not independent biological loops."
  ),
  "approximate_locus_sensitivity",
  paste0(
    "Cross-resolution calls were grouped only for sensitivity using ",
    "all-pairs-constrained two-anchor centroid distances of 20/20/50 kb ",
    "for 5/10/25-kb source calls."
  ),
  "recurrent_core",
  paste0(
    "A recurrent feature is supported by at least two of the same seven ",
    "libraries at a given downsampling depth; support thresholds 1-7 are ",
    "reported."
  ),
  "regulatory_annotation",
  paste0(
    "Strand-aware Ensembl/EPD TSS +/-1-kb promoter evidence and >=50-bp ",
    "TSS-excluded Duttke ATAC overlap were reapplied independently at each depth."
  ),
  "gene_analysis",
  paste0(
    "Gene-loop counts are descriptive sensitivity metrics for revised ",
    "putative-regulatory exact calls. Unique approximate loci per gene are ",
    "reported separately as a position-tolerant sensitivity analysis; GO is ",
    "not run or treated as validation."
  ),
  "category_composition",
  paste0(
    "Regulatory-category percentages are compared on all calls and by ",
    "resolution. Total variation summarizes the ALL-resolution composition ",
    "difference between each ordered depth pair."
  ),
  "interpretation_limit",
  paste0(
    "The 250M series excludes 74AA, A2DB, and DA08A and is therefore a ",
    "matched seven-library depth sensitivity analysis, not a replacement for ",
    "the full ten-library pooled resource. Correlations use n=7 and are ",
    "descriptive."
  )
)

df.analysis.caveat <- tribble(
  ~analysis_issue, ~status, ~interpretive_requirement,
  "unequal_library_cohorts",
  "controlled_in_depth_series",
  paste0(
    "Full depth, 140M, and 250M are compared only within the same seven ",
    "eligible libraries; the ten-library 140M analysis is reported separately."
  ),
  "small_correlation_sample",
  "remaining_design_limit",
  paste0(
    "Depth-loop correlations use seven libraries. Report rho, P value, and n ",
    "as descriptive sensitivity statistics, not strain-level inference."
  ),
  "single_random_seed",
  "remaining_monte_carlo_limit",
  paste0(
    "One nested seed realization was analyzed. Random downsampling uncertainty ",
    "is not estimated by repeated seeds."
  ),
  "one_library_per_strain",
  "remaining_biological_limit",
  paste0(
    "Depth normalization cannot distinguish strain effects from library or ",
    "specimen effects and does not create biological replication."
  ),
  "valid_contacts_not_all_library_properties",
  "remaining_technical_limit",
  paste0(
    "The analysis equalizes duplicate-removed MAPQ >=30 valid contacts but ",
    "does not equalize contact-distance composition or all library-quality ",
    "features."
  ),
  "five_kb_ignore_sparsity",
  "resolution_specific_sensitivity",
  paste0(
    "HiCCUPS used --ignore-sparsity at all three resolutions. The low 5-kb ",
    "yield at 140M is therefore reported separately and interpreted as depth ",
    "sensitivity rather than absence of 5-kb loops."
  ),
  "full_depth_recovery_depends_on_retained_fraction",
  "expected_sampling_property",
  paste0(
    "Recovery of a full-depth call set is expected to be higher when a larger ",
    "fraction of the source contacts is retained. Recovery is therefore a ",
    "depth-sensitivity metric, not technical or biological replication."
  ),
  "approximate_locus_definition",
  "sensitivity_only",
  paste0(
    "Approximate loci are jointly grouped across the three conditions and may ",
    "produce optimistic positional concordance. Exact calls and same-resolution ",
    "HiCCUPS-radius recovery remain the primary comparisons."
  ),
  "nested_contact_subset_not_locally_verified",
  "provenance_recorded_not_contact_rechecked",
  paste0(
    "Matching source counts, seeds, and the sequential sampling algorithm ",
    "support the intended nested 140M-within-250M design, but the large raw ",
    "contact files were not compared record-by-record on this computer."
  ),
  "all_resolution_category_mix",
  "report_resolution_strata",
  paste0(
    "ALL-resolution category composition depends partly on the changing mix ",
    "of 5-, 10-, and 25-kb calls. Resolution-stratified category tables must ",
    "accompany the pooled summary."
  ),
  "recurrent_library_support",
  "technical_recurrence_only",
  paste0(
    "Support in two or more libraries is technical recurrence in this dataset, ",
    "not biological replication, conservation, or a consensus loop claim."
  )
)

df.run.metadata <- tibble(
  field = c(
    "run_timestamp", "r_files_dir", "full_depth_root", "downsample_140M_root",
    "downsample_250M_root", "coordinate_cache", "output_dir",
    "n_common_libraries", "max_loop_distance_bp", "promoter_window_flank_bp",
    "minimum_atac_overlap_bp", "downsampling_seed_realizations",
    "contact_level_subset_membership_verified_locally"
  ),
  value = c(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    r.files.dir,
    full.depth.root,
    downsample.140m.root,
    downsample.250m.root,
    coord.cache.dir,
    output.dir,
    as.character(common.library.count),
    as.character(max.loop.distance.bp),
    as.character(promoter.window.flank.bp),
    as.character(atac.minimum.overlap.bp),
    "1",
    "FALSE"
  )
)

tables <- list(
  comparison_input_status = df.input.status,
  comparison_sample_metadata = df.sample.metadata,
  downsampling_qc = df.downsampling.qc,
  nested_downsampling_design = df.nested.downsampling.design,
  exact_pool_summary = df.exact.pool.summary,
  loop_counts_by_library_resolution = df.loop.count,
  loop_count_cv_and_depth_correlation = df.loop.count.summary,
  exact_call_recovery = df.exact.recovery,
  hiccups_radius_recovery = df.radius.recovery,
  pooled_exact_recovery = df.pooled.exact.recovery,
  pooled_approximate_locus_recovery = df.pooled.approximate.recovery,
  sampling_fraction_recovery = df.sampling.fraction.recovery,
  sampling_fraction_recovery_summary =
    df.sampling.fraction.recovery.summary,
  pairwise_library_similarity_exact = df.exact.pairwise,
  pairwise_library_similarity_approximate = df.approximate.pairwise,
  library_support_distribution_exact = df.exact.support,
  library_support_distribution_approximate = df.approximate.support,
  recurrent_core_stability_140M_vs_250M = df.recurrent.core.stability,
  recurrent_core_membership_min2 = df.recurrent.core.membership,
  pooled_resource_library_influence = df.pool.library.influence,
  approximate_locus_map = df.approximate.key,
  approximate_locus_summary = df.approximate.locus.summary,
  approximate_locus_method = approximate.resource$method,
  regulatory_category_counts = df.category.count,
  regulatory_category_concordance = df.category.concordance,
  regulatory_category_composition_stability =
    df.category.composition.stability,
  regulatory_category_composition_summary =
    df.category.composition.summary,
  regulatory_gene_loop_counts = df.gene.count,
  regulatory_gene_approximate_locus_counts = df.gene.approximate.count,
  regulatory_gene_stability_summary = df.gene.stability.summary,
  regulatory_gene_stability_detail = df.gene.stability.detail,
  method_definitions = df.method.definition,
  analysis_caveats = df.analysis.caveat,
  run_metadata = df.run.metadata
)
write_depth_tables(tables, output.dir)

p.loop.count <- df.loop.count %>%
  filter(resolution == "ALL") %>%
  mutate(
    condition = factor(
      condition,
      levels = c("full_depth", "downsample_250M", "downsample_140M")
    )
  ) %>%
  ggplot(aes(condition, n_loop_calls, group = sample, color = sample)) +
  geom_line(linewidth = 0.5, alpha = 0.7) +
  geom_point(size = 2) +
  labs(
    x = NULL,
    y = "HiCCUPS calls shorter than 2 Mb",
    title = "Matched seven-library loop-call sensitivity",
    subtitle = "Full depth, 250M, and 140M valid contacts"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")
ggsave(
  file.path(output.dir, "loop_counts_full_140M_250M.png"),
  p.loop.count,
  width = 8,
  height = 5,
  dpi = 300
)

writeLines(capture.output(sessionInfo()), file.path(output.dir, "session_info.txt"))

message("Completed full-depth/140M/250M common-seven comparison.")
message("Outputs: ", output.dir)
print(df.loop.count.summary %>% filter(resolution == "ALL"))
print(df.recurrent.core.stability %>% filter(minimum_library_support == 2L))
print(df.gene.stability.summary)
