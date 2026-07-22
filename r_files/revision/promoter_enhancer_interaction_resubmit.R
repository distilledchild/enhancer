# lintr: disable
library("tidyverse")
library("GenomicRanges")
library("GenomeInfoDb")

options(tibble.width = Inf)
options(tibble.print_max = Inf)
options(tibble.max_extra_cols = Inf)
options(scipen = 999)

################################################################################
# Resubmission analysis: pooled rat frontal cortex chromatin-loop annotation
#
# Main revision framework
# 1. Use all distinct HiCCUPS loops shorter than 2 Mb as the pooled resource.
# 2. Treat CTCF motif support as structural evidence.
# 3. Treat promoter/TSS assignment as gene-annotation evidence.
# 4. Treat non-TSS ATAC overlap as open-chromatin support, not direct proof of
#    enhancer activity.
# 5. Divide the pooled resource into three mutually exclusive evidence groups:
#    - putative regulatory loops with promoter/TSS + non-TSS ATAC support
#    - structural CTCF-supported loops without that ATAC regulatory support
#    - lower-support or uncertain loops
# 6. Highlight the promoter/TSS + non-TSS ATAC + CTCF >=6 subset as a more
#    conservative subset, not as the complete resource.
# 7. Keep strain-level loop sharing and sequencing-depth analyses exploratory.
#
# This script follows the object naming, numbered section structure, and
# explicit intermediate checks used in:
# - enhancer_promoter_interaction_figures.R
# - enhancer_promoter_interaction.R
################################################################################

########################
# 0. Directories and files
########################

r.files.dir <- path.expand("~/Desktop/playground/enhancer/r_files")
if (!dir.exists(r.files.dir)) {
  r.files.dir <- path.expand("~/dropbox/Gateway_to_Hao/enhancer/r_files")
}
if (!dir.exists(r.files.dir)) {
  stop("Cannot locate the enhancer/r_files directory.", call. = FALSE)
}

setwd(r.files.dir)

revision.dir <- file.path(r.files.dir, "revision")
output.dir <- Sys.getenv(
  "RESUBMIT_OUTPUT_DIR",
  unset = file.path(revision.dir, "resubmit_outputs")
)
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

dropbox.root <- path.expand("~/dropbox/Gateway_to_Hao/enhancer")
cloud.dropbox.root <- "/Users/pete/Library/CloudStorage/Dropbox-UTHSCGGI/K P/Gateway_to_Hao/enhancer"

distinct.loop.file <- file.path(
  dropbox.root,
  "r_files/rds/df.DISTINCT.loop.deep.sample.all.rds"
)
sample.loop.file <- file.path(
  dropbox.root,
  "r_files/rds/df.loop.deep.sample.all.rds"
)
ctcf.overlap.file <- file.path(
  dropbox.root,
  "r_files/rds/df.overlapping.CTCF.w.BOTH.result.rds"
)
promoter.tss.candidate.file <- file.path(
  dropbox.root,
  "r_files/rds/df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3_final_200kb.rds"
)
current.final.loop.file <- file.path(
  dropbox.root,
  "r_files/figures/submission/lt2mb/df_final_loop_sub.4.any.lt2mb.ENSEMBL.mid.mid.final.filter.200kb.csv"
)
atac.peak.file <- file.path(
  cloud.dropbox.root,
  "data/Duttke2022_snATAC_peaks_rn7.narrowPeak"
)
tss.file <- file.path(
  dropbox.root,
  "r_files/rds/df.tss.ensembl.rds"
)
library.complexity.file <- file.path(
  dropbox.root,
  "data/library_complexity_592BB.tsv"
)

# Optional long-format table with strain1, strain2, and genetic_distance.
# The environment variable takes priority over the default revision path.
genetic.distance.file <- Sys.getenv(
  "HRDP_GENETIC_DISTANCE_FILE",
  unset = file.path(revision.dir, "hrdp_genetic_distance.tsv")
)

required.files <- c(
  distinct.loop.file,
  sample.loop.file,
  ctcf.overlap.file,
  promoter.tss.candidate.file,
  current.final.loop.file,
  atac.peak.file,
  tss.file,
  library.complexity.file
)

missing.files <- required.files[!file.exists(required.files)]
if (length(missing.files) > 0) {
  stop(
    "Missing required input file(s):\n",
    paste(missing.files, collapse = "\n"),
    call. = FALSE
  )
}

########################
# 0-1. Utility functions
########################

check_required_columns <- function(df, required.columns, object.name) {
  missing.columns <- setdiff(required.columns, colnames(df))
  if (length(missing.columns) > 0) {
    stop(
      object.name,
      " is missing required column(s): ",
      paste(missing.columns, collapse = ", "),
      call. = FALSE
    )
  }
}

normalise_resolution <- function(x) {
  case_when(
    as.character(x) %in% c("5000", "5K", "5k") ~ "5K",
    as.character(x) %in% c("10000", "10K", "10k") ~ "10K",
    as.character(x) %in% c("25000", "25K", "25k") ~ "25K",
    TRUE ~ as.character(x)
  )
}

ctcf_resolution_adjusted_threshold <- function(x) {
  case_when(
    normalise_resolution(x) == "5K" ~ 6L,
    normalise_resolution(x) == "10K" ~ 12L,
    normalise_resolution(x) == "25K" ~ 30L,
    TRUE ~ NA_integer_
  )
}

count_true <- function(x) {
  sum(x %in% TRUE, na.rm = TRUE)
}

percent_true <- function(x, digits = 1) {
  if (length(x) == 0) {
    return(NA_real_)
  }
  round(100 * mean(x %in% TRUE, na.rm = TRUE), digits)
}

safe_ratio <- function(numerator, denominator) {
  if_else(
    !is.na(denominator) & denominator > 0,
    numerator / denominator,
    NA_real_
  )
}

safe_correlation_summary <- function(x, y, method = "pearson") {
  valid.rows <- complete.cases(x, y) & is.finite(x) & is.finite(y)
  x <- x[valid.rows]
  y <- y[valid.rows]

  if (
    length(x) < 3 ||
      isTRUE(all.equal(stats::sd(x), 0)) ||
      isTRUE(all.equal(stats::sd(y), 0))
  ) {
    return(tibble(
      n = length(x),
      estimate = NA_real_,
      p_value = NA_real_
    ))
  }

  test.result <- suppressWarnings(
    stats::cor.test(x, y, method = method, exact = FALSE)
  )

  tibble(
    n = length(x),
    estimate = unname(test.result$estimate),
    p_value = test.result$p.value
  )
}

add_depth_residuals <- function(df) {
  df <- df %>%
    mutate(log10_hic_contacts = log10(`Hi-C_Contacts`))

  valid.rows <- complete.cases(df$n_loops, df$log10_hic_contacts) &
    is.finite(df$n_loops) &
    is.finite(df$log10_hic_contacts)

  if (
    sum(valid.rows) < 3 ||
      isTRUE(all.equal(stats::sd(df$n_loops[valid.rows]), 0)) ||
      isTRUE(all.equal(stats::sd(df$log10_hic_contacts[valid.rows]), 0))
  ) {
    return(
      df %>%
        mutate(
          fitted_loop_count_by_hic_contacts = NA_real_,
          depth_residual_loop_count = NA_real_
        )
    )
  }

  depth.model <- stats::lm(
    n_loops ~ log10_hic_contacts,
    data = df[valid.rows, , drop = FALSE]
  )

  fitted.values <- rep(NA_real_, nrow(df))
  fitted.values[valid.rows] <- as.numeric(stats::predict(
    depth.model,
    newdata = df[valid.rows, , drop = FALSE]
  ))

  df %>%
    mutate(
      fitted_loop_count_by_hic_contacts = fitted.values,
      depth_residual_loop_count = n_loops - fitted_loop_count_by_hic_contacts
    )
}

read_atac_narrowpeak <- function(file, sample.label = "Duttke2022_snATAC_PFC") {
  df.atac <- readr::read_tsv(
    file,
    col_names = c(
      "chr", "start", "end", "name", "score", "strand",
      "fc", "neglog10p", "neglog10q", "summit"
    ),
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  ) %>%
    mutate(
      start = as.integer(start),
      end = as.integer(end),
      score = suppressWarnings(as.numeric(score)),
      fc = suppressWarnings(as.numeric(fc)),
      neglog10q = suppressWarnings(as.numeric(neglog10q))
    )

  GRanges(
    seqnames = df.atac$chr,
    ranges = IRanges(
      start = df.atac$start + 1L,
      end = df.atac$end
    ),
    score = df.atac$score,
    fc = df.atac$fc,
    neglog10q = df.atac$neglog10q,
    sample = sample.label
  )
}

create_anchor_granges <- function(df, anchor.role) {
  if (anchor.role == "promoter") {
    df.anchor <- bind_rows(
      df %>%
        filter(WHERE == "UP") %>%
        transmute(
          loop_id,
          chr = chr1,
          start = start1,
          end = end1,
          resolution,
          component,
          gene_id,
          gene_name,
          WHERE,
          anchor_side = "anchor1"
        ),
      df %>%
        filter(WHERE == "DOWN") %>%
        transmute(
          loop_id,
          chr = chr2,
          start = start2,
          end = end2,
          resolution,
          component,
          gene_id,
          gene_name,
          WHERE,
          anchor_side = "anchor2"
        )
    )
  } else if (anchor.role == "candidate_regulatory") {
    df.anchor <- bind_rows(
      df %>%
        filter(WHERE == "UP") %>%
        transmute(
          loop_id,
          chr = chr2,
          start = start2,
          end = end2,
          resolution,
          component,
          gene_id,
          gene_name,
          WHERE,
          anchor_side = "anchor2"
        ),
      df %>%
        filter(WHERE == "DOWN") %>%
        transmute(
          loop_id,
          chr = chr1,
          start = start1,
          end = end1,
          resolution,
          component,
          gene_id,
          gene_name,
          WHERE,
          anchor_side = "anchor1"
        )
    )
  } else {
    stop("Unknown anchor role: ", anchor.role, call. = FALSE)
  }

  gr.anchor <- GRanges(
    seqnames = df.anchor$chr,
    ranges = IRanges(
      start = df.anchor$start,
      end = df.anchor$end
    ),
    loop_id = df.anchor$loop_id,
    resolution = df.anchor$resolution,
    component = df.anchor$component,
    gene_id = df.anchor$gene_id,
    gene_name = df.anchor$gene_name,
    WHERE = df.anchor$WHERE,
    anchor_role = anchor.role,
    anchor_side = df.anchor$anchor_side
  )

  gr.anchor
}

create_pairwise_loop_overlap <- function(
    df.loop.presence,
    strain.levels,
    resolution.label = "ALL") {
  loop.sets <- split(df.loop.presence$loop_id, df.loop.presence$strain)
  pair.rows <- vector("list", length(strain.levels)^2)
  row.index <- 1L

  for (strain.1 in strain.levels) {
    loops.1 <- unique(loop.sets[[strain.1]])
    if (is.null(loops.1)) {
      loops.1 <- character(0)
    }

    for (strain.2 in strain.levels) {
      loops.2 <- unique(loop.sets[[strain.2]])
      if (is.null(loops.2)) {
        loops.2 <- character(0)
      }

      n.shared <- length(intersect(loops.1, loops.2))
      n.union <- length(union(loops.1, loops.2))

      pair.rows[[row.index]] <- tibble(
        resolution = resolution.label,
        strain1 = strain.1,
        strain2 = strain.2,
        n_loops_strain1 = length(loops.1),
        n_loops_strain2 = length(loops.2),
        n_shared_loops = n.shared,
        n_union_loops = n.union,
        jaccard_similarity = if (n.union > 0) n.shared / n.union else NA_real_,
        jaccard_distance = if (n.union > 0) 1 - (n.shared / n.union) else NA_real_
      )

      row.index <- row.index + 1L
    }
  }

  bind_rows(pair.rows)
}

summarise_atac_support <- function(df, resolution.label) {
  tibble(
    resolution = resolution.label,
    n_promoter_tss_candidates = nrow(df),
    n_promoter_anchor_atac = count_true(df$atac_promoter_anchor),
    pct_promoter_anchor_atac = percent_true(df$atac_promoter_anchor),
    n_candidate_regulatory_anchor_atac = count_true(
      df$atac_candidate_regulatory_anchor
    ),
    pct_candidate_regulatory_anchor_atac = percent_true(
      df$atac_candidate_regulatory_anchor
    ),
    n_candidate_regulatory_anchor_with_non_tss_fragment = count_true(
      df$candidate_regulatory_anchor_has_non_tss_fragment
    ),
    n_candidate_regulatory_anchor_non_tss_atac = count_true(
      df$atac_candidate_regulatory_anchor_non_tss
    ),
    pct_candidate_regulatory_anchor_non_tss_atac_all = percent_true(
      df$atac_candidate_regulatory_anchor_non_tss
    ),
    pct_candidate_regulatory_anchor_non_tss_atac_eligible = {
      eligible <- df$candidate_regulatory_anchor_has_non_tss_fragment %in% TRUE
      if (sum(eligible) == 0) {
        NA_real_
      } else {
        round(
          100 * mean(
            df$atac_candidate_regulatory_anchor_non_tss[eligible] %in% TRUE
          ),
          1
        )
      }
    }
  )
}

summarise_mcnemar <- function(df, resolution.label) {
  paired.table <- table(
    promoter_ATAC = factor(
      df$atac_promoter_anchor,
      levels = c(FALSE, TRUE)
    ),
    candidate_regulatory_ATAC = factor(
      df$atac_candidate_regulatory_anchor,
      levels = c(FALSE, TRUE)
    )
  )

  promoter.only <- unname(paired.table["TRUE", "FALSE"])
  candidate.only <- unname(paired.table["FALSE", "TRUE"])

  mcnemar.p <- if ((promoter.only + candidate.only) == 0) {
    NA_real_
  } else {
    suppressWarnings(stats::mcnemar.test(
      paired.table,
      correct = TRUE
    )$p.value)
  }

  tibble(
    resolution = resolution.label,
    n_pairs = nrow(df),
    n_both_atac = unname(paired.table["TRUE", "TRUE"]),
    n_promoter_only_atac = promoter.only,
    n_candidate_regulatory_only_atac = candidate.only,
    n_neither_atac = unname(paired.table["FALSE", "FALSE"]),
    mcnemar_p = mcnemar.p,
    test_note = paste0(
      "Paired promoter and candidate-regulatory anchors from the same loop; ",
      "continuity-corrected McNemar test."
    )
  )
}

read_genetic_distance_long <- function(file) {
  df.distance <- readr::read_tsv(file, show_col_types = FALSE)
  required.pair.columns <- c("strain1", "strain2")
  check_required_columns(
    df.distance,
    required.pair.columns,
    "Genetic-distance table"
  )

  distance.column <- intersect(
    c("genetic_distance", "ibs_distance", "distance"),
    colnames(df.distance)
  )

  if (length(distance.column) == 0) {
    stop(
      "The genetic-distance table must contain genetic_distance, ",
      "ibs_distance, or distance.",
      call. = FALSE
    )
  }

  df.distance %>%
    transmute(
      strain1 = as.character(strain1),
      strain2 = as.character(strain2),
      genetic_distance = as.numeric(.data[[distance.column[1]]]),
      pair_strain_a = pmin(strain1, strain2),
      pair_strain_b = pmax(strain1, strain2)
    ) %>%
    filter(
      strain1 != strain2,
      !is.na(genetic_distance),
      is.finite(genetic_distance)
    ) %>%
    distinct(pair_strain_a, pair_strain_b, .keep_all = TRUE)
}

prepare_go_gene_set <- function(df.gene.count, set.name, minimum.loops) {
  df.gene.count %>%
    filter(
      gene_set == set.name,
      n >= minimum.loops,
      !is.na(ensembl_gene_id)
    ) %>%
    distinct(ensembl_gene_id, gene_symbol, n) %>%
    arrange(desc(n), gene_symbol)
}

run_go_enrichment <- function(df.gene, universe.ensembl, set.name, out.dir) {
  required.packages <- c(
    "clusterProfiler",
    "org.Rn.eg.db",
    "AnnotationDbi"
  )
  missing.packages <- required.packages[
    !vapply(required.packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing.packages) > 0) {
    message(
      "Skipping GO for ",
      set.name,
      ": missing ",
      paste(missing.packages, collapse = ", "),
      "."
    )
    return(tibble())
  }

  ensembl.genes <- unique(na.omit(df.gene$ensembl_gene_id))
  if (length(ensembl.genes) < 5) {
    message(
      "Skipping GO for ",
      set.name,
      ": fewer than 5 Ensembl genes."
    )
    return(tibble())
  }

  df.gene.map <- AnnotationDbi::select(
    org.Rn.eg.db::org.Rn.eg.db,
    keys = ensembl.genes,
    keytype = "ENSEMBL",
    columns = c("ENTREZID", "SYMBOL")
  ) %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(ENTREZID, .keep_all = TRUE)

  df.universe.map <- AnnotationDbi::select(
    org.Rn.eg.db::org.Rn.eg.db,
    keys = unique(na.omit(universe.ensembl)),
    keytype = "ENSEMBL",
    columns = c("ENTREZID", "SYMBOL")
  ) %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(ENTREZID, .keep_all = TRUE)

  if (nrow(df.gene.map) < 5) {
    message(
      "Skipping GO for ",
      set.name,
      ": fewer than 5 mapped Entrez genes."
    )
    return(tibble())
  }

  go.result <- clusterProfiler::enrichGO(
    gene = df.gene.map$ENTREZID,
    universe = df.universe.map$ENTREZID,
    OrgDb = org.Rn.eg.db::org.Rn.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  )

  df.go.result <- as_tibble(go.result@result) %>%
    mutate(
      gene_set = set.name,
      input_ensembl_genes = length(ensembl.genes),
      mapped_entrez_genes = nrow(df.gene.map),
      universe_mapped_entrez_genes = nrow(df.universe.map)
    ) %>%
    arrange(p.adjust, pvalue)

  readr::write_tsv(
    df.go.result,
    file.path(out.dir, paste0("go_enrichment_BP_", set.name, ".tsv"))
  )

  df.go.result
}

################################################################################
# 1. Pooled HiCCUPS loop resource
################################################################################

df.loop.distinct.all <- readRDS(distinct.loop.file)
check_required_columns(
  df.loop.distinct.all,
  c(
    "loop.id", "chr1", "x1", "x2", "chr2", "y1", "y2",
    "resolution", "distance"
  ),
  "df.loop.distinct.all"
)

df.loop.distinct.all <- df.loop.distinct.all %>%
  mutate(
    loop_id = loop.id,
    resolution = normalise_resolution(resolution),
    loop_distance = as.integer(distance),
    passes_lt2mb = loop_distance < 2000000L
  )

df.loop.universe <- df.loop.distinct.all %>%
  filter(passes_lt2mb) %>%
  transmute(
    loop_id,
    chr1,
    start1 = as.integer(x1),
    end1 = as.integer(x2),
    chr2,
    start2 = as.integer(y1),
    end2 = as.integer(y2),
    resolution,
    loop_distance,
    passes_lt2mb
  ) %>%
  distinct(loop_id, .keep_all = TRUE) %>%
  arrange(chr1, start1, end1, chr2, start2, end2)

if (nrow(df.loop.universe) != n_distinct(df.loop.universe$loop_id)) {
  stop("The pooled loop resource contains duplicate loop IDs.", call. = FALSE)
}

message("All distinct HiCCUPS loops: ", n_distinct(df.loop.distinct.all$loop_id))
message("Pooled loop resource (<2 Mb): ", nrow(df.loop.universe))

################################################################################
# 2. Exploratory strain-level loop sharing
#
# These results describe exact loop-coordinate sharing among the 10 libraries.
# They are not used to claim strain-specific chromatin biology because each
# strain has one Hi-C library and loop recovery is depth-sensitive.
################################################################################

df.sample.loop.presence <- readRDS(sample.loop.file)
check_required_columns(
  df.sample.loop.presence,
  c("sample", "strain", "loop.id", "resolution", "distance"),
  "df.sample.loop.presence"
)

df.sample.loop.presence <- df.sample.loop.presence %>%
  mutate(
    sample = as.character(sample),
    strain = as.character(strain),
    loop_id = loop.id,
    resolution = normalise_resolution(resolution),
    loop_distance = as.integer(distance),
    passes_lt2mb = loop_distance < 2000000L
  ) %>%
  filter(passes_lt2mb) %>%
  distinct(loop_id, resolution, sample, strain)

# Keep all 10 sequenced strains in descriptive loop-sharing summaries.
excluded.strains <- character(0)

df.sample.loop.presence.pairwise <- df.sample.loop.presence %>%
  filter(!strain %in% excluded.strains)

df.sample.strain.key <- df.sample.loop.presence %>%
  distinct(sample, strain) %>%
  mutate(
    included_in_pairwise_loop_overlap = !strain %in% excluded.strains,
    exclusion_reason = if_else(
      included_in_pairwise_loop_overlap,
      NA_character_,
      paste0(
        "Excluded because a matching genotype sample was not available ",
        "for an optional genetic-distance analysis."
      )
    )
  ) %>%
  arrange(strain, sample)

strain.levels <- df.sample.loop.presence.pairwise %>%
  distinct(strain) %>%
  arrange(strain) %>%
  pull(strain)

df.strain.loop.count.by.resolution <- bind_rows(
  df.sample.loop.presence.pairwise %>%
    count(strain, sample, resolution, name = "n_loops"),
  df.sample.loop.presence.pairwise %>%
    count(strain, sample, name = "n_loops") %>%
    mutate(resolution = "ALL")
) %>%
  arrange(strain, sample, resolution)

df.strain.loop.presence.matrix <- df.sample.loop.presence.pairwise %>%
  distinct(loop_id, resolution, strain) %>%
  mutate(present = 1L) %>%
  pivot_wider(
    names_from = strain,
    values_from = present,
    values_fill = 0
  ) %>%
  left_join(
    df.loop.universe %>%
      dplyr::select(
        loop_id,
        chr1, start1, end1,
        chr2, start2, end2,
        loop_distance
      ),
    by = "loop_id"
  ) %>%
  relocate(
    chr1, start1, end1,
    chr2, start2, end2,
    loop_distance,
    .after = resolution
  ) %>%
  arrange(resolution, chr1, start1, end1, chr2, start2, end2)

df.shared.loop.member <- df.sample.loop.presence.pairwise %>%
  distinct(loop_id, resolution, strain) %>%
  group_by(loop_id, resolution) %>%
  summarise(
    strains = list(sort(unique(strain))),
    n_strains_detected = n_distinct(strain),
    .groups = "drop"
  ) %>%
  filter(n_strains_detected > 1)

if (nrow(df.shared.loop.member) > 0) {
  df.shared.loop.pair.by.strain <- map_dfr(
    seq_len(nrow(df.shared.loop.member)),
    function(i) {
      strain.pairs <- t(combn(df.shared.loop.member$strains[[i]], 2))
      tibble(
        loop_id = df.shared.loop.member$loop_id[[i]],
        resolution = df.shared.loop.member$resolution[[i]],
        strain1 = strain.pairs[, 1],
        strain2 = strain.pairs[, 2]
      )
    }
  ) %>%
    count(resolution, strain1, strain2, name = "n_shared_loops") %>%
    arrange(resolution, strain1, strain2)
} else {
  df.shared.loop.pair.by.strain <- tibble(
    resolution = character(),
    strain1 = character(),
    strain2 = character(),
    n_shared_loops = integer()
  )
}

df.pairwise.loop.overlap.by.strain <- bind_rows(
  create_pairwise_loop_overlap(
    df.sample.loop.presence.pairwise,
    strain.levels,
    "ALL"
  ),
  map_dfr(
    sort(unique(df.sample.loop.presence.pairwise$resolution)),
    function(resolution.i) {
      create_pairwise_loop_overlap(
        df.sample.loop.presence.pairwise %>%
          filter(resolution == resolution.i),
        strain.levels,
        resolution.i
      )
    }
  )
) %>%
  arrange(resolution, strain1, strain2)

df.loop.sharing.distribution.by.strain <- bind_rows(
  df.sample.loop.presence.pairwise %>%
    distinct(loop_id, resolution, strain) %>%
    group_by(loop_id, resolution) %>%
    summarise(
      n_strains_detected = n_distinct(strain),
      .groups = "drop"
    ) %>%
    count(resolution, n_strains_detected, name = "n_loops"),
  df.sample.loop.presence.pairwise %>%
    distinct(loop_id, strain) %>%
    group_by(loop_id) %>%
    summarise(
      n_strains_detected = n_distinct(strain),
      .groups = "drop"
    ) %>%
    count(n_strains_detected, name = "n_loops") %>%
    mutate(resolution = "ALL")
) %>%
  arrange(resolution, n_strains_detected)

################################################################################
# 3. Sequencing-depth sensitivity/QC
#
# Contact-normalized ratios and regression residuals are reviewer-facing QC.
# They do not replace read downsampling followed by mapping and loop re-calling,
# and they do not make one-library-per-strain comparisons biological replicates.
################################################################################

df.library.complexity <- readr::read_tsv(
  library.complexity.file,
  show_col_types = FALSE
) %>%
  mutate(Strain = as.character(Strain))

required.depth.columns <- c(
  "Strain",
  "Sequenced_RP",
  "Unique_Reads",
  "Alignable_Normal_N_Chimeric",
  "Hi-C_Contacts",
  "Long_Range_20Kb",
  "PCR_Duplicates",
  "Optical_Duplicates"
)
check_required_columns(
  df.library.complexity,
  required.depth.columns,
  "df.library.complexity"
)

df.depth.qc.by.strain <- df.strain.loop.count.by.resolution %>%
  left_join(
    df.library.complexity,
    by = c("strain" = "Strain")
  ) %>%
  mutate(
    sequenced_read_pairs_millions = Sequenced_RP / 1e6,
    unique_reads_millions = Unique_Reads / 1e6,
    hic_contacts_millions = `Hi-C_Contacts` / 1e6,
    loops_per_100m_sequenced_read_pairs = safe_ratio(
      n_loops,
      Sequenced_RP / 1e8
    ),
    loops_per_100m_unique_reads = safe_ratio(
      n_loops,
      Unique_Reads / 1e8
    ),
    loops_per_100m_hic_contacts = safe_ratio(
      n_loops,
      `Hi-C_Contacts` / 1e8
    ),
    unique_read_pct = 100 * safe_ratio(Unique_Reads, Sequenced_RP),
    duplicate_pct = 100 * safe_ratio(
      PCR_Duplicates + Optical_Duplicates,
      Sequenced_RP
    )
  ) %>%
  group_by(resolution) %>%
  group_modify(~ add_depth_residuals(.x)) %>%
  ungroup() %>%
  arrange(resolution, desc(n_loops))

if (any(is.na(df.depth.qc.by.strain$`Hi-C_Contacts`))) {
  warning(
    "One or more strains did not match the library-complexity table.",
    call. = FALSE
  )
}

depth.metrics <- c(
  "Sequenced_RP",
  "Unique_Reads",
  "Alignable_Normal_N_Chimeric",
  "Hi-C_Contacts",
  "Long_Range_20Kb"
)

df.depth.loop.correlation.by.resolution <- map_dfr(
  sort(unique(df.depth.qc.by.strain$resolution)),
  function(resolution.i) {
    df.depth.i <- df.depth.qc.by.strain %>%
      filter(resolution == resolution.i)

    map_dfr(depth.metrics, function(metric.i) {
      pearson.result <- safe_correlation_summary(
        df.depth.i$n_loops,
        df.depth.i[[metric.i]],
        method = "pearson"
      )
      spearman.result <- safe_correlation_summary(
        df.depth.i$n_loops,
        df.depth.i[[metric.i]],
        method = "spearman"
      )

      tibble(
        resolution = resolution.i,
        depth_metric = metric.i,
        n_strains = pearson.result$n,
        pearson_r = pearson.result$estimate,
        pearson_p = pearson.result$p_value,
        spearman_rho = spearman.result$estimate,
        spearman_p = spearman.result$p_value
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
    mean_sequenced_read_pairs_millions = round(
      mean(sequenced_read_pairs_millions),
      1
    ),
    mean_unique_reads_millions = round(mean(unique_reads_millions), 1),
    mean_hic_contacts_millions = round(mean(hic_contacts_millions), 1),
    mean_loops_per_100m_hic_contacts = round(
      mean(loops_per_100m_hic_contacts),
      1
    ),
    .groups = "drop"
  ) %>%
  arrange(resolution)

df.depth.analysis.interpretation <- tibble(
  analysis = c(
    "loop_count_vs_depth",
    "loops_per_100m_hic_contacts",
    "depth_residual_loop_count"
  ),
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
# 4. CTCF structural evidence
################################################################################

df.ctcf.overlap <- readRDS(ctcf.overlap.file)
check_required_columns(
  df.ctcf.overlap,
  c("loop.id", "WHERE", "ctcf.id", "resolution"),
  "df.ctcf.overlap"
)

df.ctcf.overlap <- df.ctcf.overlap %>%
  mutate(
    loop_id = loop.id,
    resolution = normalise_resolution(resolution),
    WHERE = as.character(WHERE)
  )

df.ctcf.anchor.count <- df.ctcf.overlap %>%
  group_by(loop_id, WHERE) %>%
  summarise(
    ctcf_count = n_distinct(ctcf.id),
    .groups = "drop"
  ) %>%
  mutate(
    WHERE = case_when(
      WHERE == "UP" ~ "ctcf_count_anchor1",
      WHERE == "DOWN" ~ "ctcf_count_anchor2",
      TRUE ~ paste0("ctcf_count_", WHERE)
    )
  ) %>%
  pivot_wider(
    names_from = WHERE,
    values_from = ctcf_count,
    values_fill = 0
  )

df.ctcf.evidence <- df.loop.universe %>%
  dplyr::select(loop_id, resolution) %>%
  left_join(df.ctcf.anchor.count, by = "loop_id") %>%
  mutate(
    ctcf_count_anchor1 = coalesce(ctcf_count_anchor1, 0L),
    ctcf_count_anchor2 = coalesce(ctcf_count_anchor2, 0L),
    ctcf_min_anchor_count = pmin(
      ctcf_count_anchor1,
      ctcf_count_anchor2
    ),
    ctcf_both_anchors_any = (
      ctcf_count_anchor1 > 0 &
        ctcf_count_anchor2 > 0
    ),
    passes_ctcf_ge6_both_anchors = (
      ctcf_count_anchor1 >= 6 &
        ctcf_count_anchor2 >= 6
    ),
    ctcf_resolution_adjusted_min_threshold =
      ctcf_resolution_adjusted_threshold(resolution),
    passes_ctcf_resolution_adjusted_both_anchors = (
      !is.na(ctcf_resolution_adjusted_min_threshold) &
        ctcf_count_anchor1 >= ctcf_resolution_adjusted_min_threshold &
        ctcf_count_anchor2 >= ctcf_resolution_adjusted_min_threshold
    ),
    ctcf_support_level = case_when(
      passes_ctcf_ge6_both_anchors ~ "strong_both_anchors_ge6",
      ctcf_both_anchors_any ~ "both_anchors_any_motif",
      ctcf_count_anchor1 > 0 | ctcf_count_anchor2 > 0 ~
        "one_anchor_any_motif",
      TRUE ~ "no_ctcf_motif"
    )
  )

################################################################################
# 5. Promoter/TSS gene-annotation evidence
################################################################################

df.promoter.tss.candidate <- readRDS(promoter.tss.candidate.file)
check_required_columns(
  df.promoter.tss.candidate,
  c(
    "loop.id", "chr1", "x1", "x2", "chr2", "y1", "y2",
    "resolution", "component", "distance", "gene_id", "gene_name",
    "WHERE", "classification"
  ),
  "df.promoter.tss.candidate"
)

df.promoter.tss.candidate <- df.promoter.tss.candidate %>%
  dplyr::rename(loop_id = loop.id) %>%
  mutate(
    resolution = normalise_resolution(resolution),
    start1 = as.integer(x1) + 1L,
    end1 = as.integer(x2),
    start2 = as.integer(y1) + 1L,
    end2 = as.integer(y2),
    promoter_support = component == "pro",
    tss_support = component == "tss",
    promoter_or_tss_support = promoter_support | tss_support,
    promoter_anchor_side = case_when(
      WHERE == "UP" ~ "anchor1",
      WHERE == "DOWN" ~ "anchor2",
      TRUE ~ NA_character_
    ),
    candidate_regulatory_anchor_side = case_when(
      WHERE == "UP" ~ "anchor2",
      WHERE == "DOWN" ~ "anchor1",
      TRUE ~ NA_character_
    )
  )

if (
  nrow(df.promoter.tss.candidate) !=
    n_distinct(df.promoter.tss.candidate$loop_id)
) {
  stop(
    "Promoter/TSS candidate input must contain one selected annotation per loop.",
    call. = FALSE
  )
}

df.promoter.tss.evidence <- df.promoter.tss.candidate %>%
  transmute(
    loop_id,
    promoter_or_tss_support,
    promoter_support,
    tss_support,
    promoter_tss_component = component,
    promoter_tss_distance = as.integer(distance),
    promoter_tss_gene_id = gene_id,
    promoter_tss_gene_name = gene_name,
    promoter_anchor_side,
    candidate_regulatory_anchor_side,
    promoter_tss_where = WHERE,
    promoter_tss_classification = classification,
    passes_promoter_tss_200kb = (
      promoter_or_tss_support &
        promoter_tss_distance < 200000L
    )
  )

################################################################################
# 6. ATAC open-chromatin evidence
#
# The candidate-regulatory anchor is the loop anchor opposite the selected
# promoter/TSS anchor. ATAC overlap is supportive evidence of accessibility.
# It is not direct evidence that the anchor functions as an enhancer.
################################################################################

gr.atac <- read_atac_narrowpeak(atac.peak.file)
gr.atac.union <- GenomicRanges::reduce(gr.atac)

gr.promoter.anchor <- create_anchor_granges(
  df.promoter.tss.candidate,
  "promoter"
)
gr.candidate.regulatory.anchor <- create_anchor_granges(
  df.promoter.tss.candidate,
  "candidate_regulatory"
)

if (
  length(gr.promoter.anchor) != nrow(df.promoter.tss.candidate) ||
    length(gr.candidate.regulatory.anchor) !=
      nrow(df.promoter.tss.candidate)
) {
  stop(
    "Anchor construction did not preserve one promoter and one ",
    "candidate-regulatory anchor per candidate loop.",
    call. = FALSE
  )
}

promoter.atac.hit <- countOverlaps(
  gr.promoter.anchor,
  gr.atac.union,
  minoverlap = 50
) > 0

candidate.regulatory.atac.hit <- countOverlaps(
  gr.candidate.regulatory.anchor,
  gr.atac.union,
  minoverlap = 50
) > 0

df.tss <- readRDS(tss.file)
check_required_columns(
  df.tss,
  c("chr", "start", "end", "strand"),
  "df.tss"
)

gr.tss <- GRanges(
  seqnames = df.tss$chr,
  ranges = IRanges(
    start = df.tss$start,
    end = df.tss$end
  ),
  strand = df.tss$strand
)

gr.tss.region <- promoters(
  gr.tss,
  upstream = 1000,
  downstream = 1000
) %>%
  GenomicRanges::reduce(ignore.strand = TRUE)

common.seqlevels.atac.tss <- intersect(
  seqlevels(gr.atac.union),
  seqlevels(gr.tss.region)
)

gr.atac.for.tss <- keepSeqlevels(
  gr.atac.union,
  common.seqlevels.atac.tss,
  pruning.mode = "coarse"
)
gr.tss.region.for.atac <- keepSeqlevels(
  gr.tss.region,
  common.seqlevels.atac.tss,
  pruning.mode = "coarse"
)

gr.atac.non.tss <- GenomicRanges::setdiff(
  gr.atac.for.tss,
  gr.tss.region.for.atac,
  ignore.strand = TRUE
)

tss.anchor.hit <- findOverlaps(
  gr.candidate.regulatory.anchor,
  gr.tss.region,
  ignore.strand = TRUE
)

tss.bp.by.candidate.anchor <- numeric(
  length(gr.candidate.regulatory.anchor)
)

if (length(tss.anchor.hit) > 0) {
  gr.tss.anchor.overlap <- pintersect(
    gr.candidate.regulatory.anchor[queryHits(tss.anchor.hit)],
    gr.tss.region[subjectHits(tss.anchor.hit)],
    ignore.strand = TRUE
  )

  tss.bp.sum <- rowsum(
    width(gr.tss.anchor.overlap),
    group = queryHits(tss.anchor.hit),
    reorder = FALSE
  )

  tss.bp.by.candidate.anchor[
    as.integer(rownames(tss.bp.sum))
  ] <- tss.bp.sum[, 1]
}

candidate.anchor.has.non.tss.fragment <- (
  tss.bp.by.candidate.anchor <
    width(gr.candidate.regulatory.anchor)
)

candidate.regulatory.non.tss.atac.hit <- (
  candidate.anchor.has.non.tss.fragment &
    (
      countOverlaps(
        gr.candidate.regulatory.anchor,
        gr.atac.non.tss,
        minoverlap = 50
      ) > 0
    )
)

df.atac.evidence <- tibble(
  loop_id = gr.promoter.anchor$loop_id,
  atac_promoter_anchor = promoter.atac.hit,
  atac_candidate_regulatory_anchor = candidate.regulatory.atac.hit,
  candidate_regulatory_anchor_has_non_tss_fragment =
    candidate.anchor.has.non.tss.fragment,
  candidate_regulatory_anchor_non_tss_bp =
    width(gr.candidate.regulatory.anchor) -
    tss.bp.by.candidate.anchor,
  atac_candidate_regulatory_anchor_non_tss =
    candidate.regulatory.non.tss.atac.hit
)

################################################################################
# 7. Previous final set and loop-level evidence table
################################################################################

df.current.final.loop <- readr::read_csv(
  current.final.loop.file,
  show_col_types = FALSE
)
check_required_columns(
  df.current.final.loop,
  c("loop.id", "category"),
  "df.current.final.loop"
)

df.current.final.loop <- df.current.final.loop %>%
  dplyr::rename(loop_id = loop.id) %>%
  transmute(
    loop_id,
    current_final_set = TRUE,
    current_final_category = category
  )

df.loop.evidence <- df.loop.universe %>%
  left_join(
    df.ctcf.evidence %>% dplyr::select(-resolution),
    by = "loop_id"
  ) %>%
  left_join(df.promoter.tss.evidence, by = "loop_id") %>%
  left_join(df.atac.evidence, by = "loop_id") %>%
  left_join(df.current.final.loop, by = "loop_id") %>%
  mutate(
    promoter_or_tss_support = coalesce(
      promoter_or_tss_support,
      FALSE
    ),
    promoter_support = coalesce(promoter_support, FALSE),
    tss_support = coalesce(tss_support, FALSE),
    passes_promoter_tss_200kb = coalesce(
      passes_promoter_tss_200kb,
      FALSE
    ),
    atac_promoter_anchor = coalesce(
      atac_promoter_anchor,
      FALSE
    ),
    atac_candidate_regulatory_anchor = coalesce(
      atac_candidate_regulatory_anchor,
      FALSE
    ),
    candidate_regulatory_anchor_has_non_tss_fragment = coalesce(
      candidate_regulatory_anchor_has_non_tss_fragment,
      FALSE
    ),
    candidate_regulatory_anchor_non_tss_bp = coalesce(
      candidate_regulatory_anchor_non_tss_bp,
      0
    ),
    atac_candidate_regulatory_anchor_non_tss = coalesce(
      atac_candidate_regulatory_anchor_non_tss,
      FALSE
    ),
    current_final_set = coalesce(current_final_set, FALSE),
    current_final_category = replace_na(
      current_final_category,
      "not_current_final"
    ),
    ensembl_gene_id = str_extract(
      promoter_tss_gene_id,
      "ENSRNOG[0-9]+"
    ),
    gene_symbol = promoter_tss_gene_name,
    putative_regulatory_support = (
      passes_promoter_tss_200kb &
        atac_candidate_regulatory_anchor_non_tss
    ),
    putative_regulatory_support_any_atac = (
      passes_promoter_tss_200kb &
        atac_candidate_regulatory_anchor
    ),
    strongest_candidate_support = (
      putative_regulatory_support &
        passes_ctcf_ge6_both_anchors
    ),
    proposed_major_category = case_when(
      putative_regulatory_support ~
        "putative_regulatory_promoter_TSS_nonTSS_ATAC_supported",
      passes_ctcf_ge6_both_anchors ~
        "structural_CTCF_supported_not_ATAC_regulatory",
      TRUE ~ "lower_support_or_uncertain"
    ),
    proposed_detailed_category = case_when(
      strongest_candidate_support ~
        "strongest_candidate_PE_promoter_TSS_nonTSS_ATAC_CTCF",
      putative_regulatory_support &
        !passes_ctcf_ge6_both_anchors ~
        "putative_regulatory_promoter_TSS_nonTSS_ATAC_no_strong_CTCF",
      promoter_or_tss_support &
        atac_candidate_regulatory_anchor &
        !atac_candidate_regulatory_anchor_non_tss &
        passes_ctcf_ge6_both_anchors ~
        "CTCF_promoter_TSS_ATAC_signal_TSS_sensitive",
      promoter_or_tss_support &
        atac_candidate_regulatory_anchor &
        !atac_candidate_regulatory_anchor_non_tss ~
        "promoter_TSS_ATAC_signal_TSS_sensitive",
      promoter_or_tss_support &
        passes_ctcf_ge6_both_anchors ~
        "structural_CTCF_promoter_TSS_no_nonTSS_ATAC",
      passes_ctcf_ge6_both_anchors ~
        "structural_CTCF_only_no_promoter_TSS_nonTSS_ATAC",
      promoter_or_tss_support ~
        "promoter_TSS_only_no_nonTSS_ATAC_no_strong_CTCF",
      ctcf_both_anchors_any ~
        "weak_structural_CTCF_both_anchors_below_ge6",
      TRUE ~ "low_support_or_uncertain"
    )
  ) %>%
  arrange(chr1, start1, end1, chr2, start2, end2)

if (nrow(df.loop.evidence) != nrow(df.loop.universe)) {
  stop(
    "Loop-level joins changed the pooled loop-resource row count.",
    call. = FALSE
  )
}

if (n_distinct(df.loop.evidence$loop_id) != nrow(df.loop.evidence)) {
  stop(
    "The loop-level evidence table is not one row per loop.",
    call. = FALSE
  )
}

if (
  any(
    df.loop.evidence$strongest_candidate_support &
      df.loop.evidence$proposed_major_category !=
        "putative_regulatory_promoter_TSS_nonTSS_ATAC_supported"
  )
) {
  stop(
    "The conservative highlighted subset is not contained in ",
    "the putative-regulatory category.",
    call. = FALSE
  )
}

################################################################################
# 8. Filtering audit and mutually exclusive category summaries
################################################################################

n.all.distinct <- n_distinct(df.loop.distinct.all$loop_id)
n.lt2mb <- nrow(df.loop.universe)
n.ctcf.any.both <- count_true(df.loop.evidence$ctcf_both_anchors_any)
n.ctcf.ge6.both <- count_true(
  df.loop.evidence$passes_ctcf_ge6_both_anchors
)
n.promoter.tss <- count_true(
  df.loop.evidence$promoter_or_tss_support
)
n.promoter.tss.200kb <- count_true(
  df.loop.evidence$passes_promoter_tss_200kb
)
n.current.final <- count_true(df.loop.evidence$current_final_set)
n.candidate.any.atac <- count_true(
  df.loop.evidence$putative_regulatory_support_any_atac
)
n.candidate.non.tss.atac <- count_true(
  df.loop.evidence$putative_regulatory_support
)
n.strongest.candidate <- count_true(
  df.loop.evidence$strongest_candidate_support
)

df.filtering.logic.audit <- tribble(
  ~filter_step,
  ~code_source,
  ~evidence_type,
  ~input_n,
  ~output_n,
  ~filter_condition,
  ~interpretation,
  ~revision_action,
  "raw_distinct_hiccups_loops",
  "enhancer_promoter_interaction_figures.R: loop preprocessing",
  "technical",
  n.all.distinct,
  n.all.distinct,
  "Load distinct HiCCUPS loops across all 10 libraries.",
  "Pooled physical loop universe before revision filters.",
  "Retain as the starting resource universe.",
  "loop_distance_lt_2mb",
  "enhancer_promoter_interaction.R: Section 1",
  "technical",
  n.all.distinct,
  n.lt2mb,
  "loop distance < 2,000,000 bp",
  "Technical scope filter; not regulatory evidence.",
  "Retain and describe as a loop-length scope threshold.",
  "ctcf_any_motif_both_anchors",
  "enhancer_promoter_interaction.R: Section 2",
  "CTCF-driven",
  n.lt2mb,
  n.ctcf.any.both,
  "At least one predicted CTCF motif at both loop anchors.",
  "Structural CTCF support; not direct enhancer evidence.",
  "Use as structural annotation.",
  "ctcf_ge6_both_anchors_previous_filter",
  "enhancer_promoter_interaction.R: CTCF threshold filtering",
  "CTCF-driven",
  n.lt2mb,
  n.ctcf.ge6.both,
  "At least 6 predicted CTCF motifs at both loop anchors.",
  "Strong structural motif support under the previous threshold.",
  "Show threshold sensitivity; do not use alone as regulatory proof.",
  "promoter_or_tss_directional_annotation_200kb",
  "enhancer_promoter_interaction.R: directional gene assignment",
  "promoter/TSS-driven",
  n.lt2mb,
  n.promoter.tss.200kb,
  "Selected promoter/TSS is within 200 kb of its assigned anchor.",
  "Target-gene annotation support; not enhancer activity evidence.",
  "Retain as promoter/TSS evidence.",
  "current_final_loop_set",
  "enhancer_promoter_interaction.R: previous final intersection",
  "combined_CTCF_promoter_TSS",
  n.lt2mb,
  n.current.final,
  "CTCF >=6 at both anchors plus promoter/TSS annotation.",
  "Structurally strong but overclaimed if called validated P-E interactions.",
  "Reframe as the previous CTCF-supported promoter/TSS set.",
  "candidate_regulatory_anchor_ATAC_any",
  "atac_validation.R: paired anchor overlap",
  "ATAC-driven",
  n.promoter.tss,
  n.candidate.any.atac,
  "Opposite anchor overlaps an ATAC peak by at least 50 bp.",
  "Open-chromatin support at the candidate-regulatory anchor.",
  "Use as supportive, not functional, evidence.",
  "candidate_regulatory_anchor_nonTSS_ATAC",
  "atac_validation.R: TSS-exclusion sensitivity analysis",
  "ATAC-driven",
  n.promoter.tss,
  n.candidate.non.tss.atac,
  "Opposite anchor overlaps ATAC signal remaining after TSS +/-1 kb removal.",
  "Conservative accessibility support not explained only by known TSS signal.",
  "Use to define the putative-regulatory evidence category."
)

df.loop.evidence.summary <- bind_rows(
  df.loop.evidence %>%
    count(proposed_major_category, name = "n_loops") %>%
    mutate(
      summary_type = "proposed_major_category",
      category = proposed_major_category
    ) %>%
    dplyr::select(summary_type, category, n_loops),
  df.loop.evidence %>%
    count(proposed_detailed_category, name = "n_loops") %>%
    mutate(
      summary_type = "proposed_detailed_category",
      category = proposed_detailed_category
    ) %>%
    dplyr::select(summary_type, category, n_loops),
  df.loop.evidence %>%
    count(current_final_set, name = "n_loops") %>%
    mutate(
      summary_type = "current_final_set",
      category = as.character(current_final_set)
    ) %>%
    dplyr::select(summary_type, category, n_loops),
  df.loop.evidence %>%
    count(
      current_final_category,
      proposed_major_category,
      name = "n_loops"
    ) %>%
    mutate(
      summary_type = "current_category_by_new_major_category",
      category = paste(
        current_final_category,
        proposed_major_category,
        sep = " | "
      )
    ) %>%
    dplyr::select(summary_type, category, n_loops)
) %>%
  arrange(summary_type, desc(n_loops), category)

df.loop.category.by.resolution <- bind_rows(
  df.loop.evidence %>%
    count(
      resolution,
      proposed_major_category,
      name = "n_loops"
    ) %>%
    group_by(resolution) %>%
    mutate(
      pct_loops = round(100 * n_loops / sum(n_loops), 1)
    ) %>%
    ungroup(),
  df.loop.evidence %>%
    count(proposed_major_category, name = "n_loops") %>%
    mutate(
      resolution = "ALL",
      pct_loops = round(100 * n_loops / sum(n_loops), 1)
    )
) %>%
  arrange(resolution, desc(n_loops))

category.count.check <- df.loop.evidence %>%
  count(proposed_major_category, name = "n_loops")

if (
  sum(category.count.check$n_loops) != n.lt2mb ||
    nrow(category.count.check) != 3
) {
  stop(
    "The three major categories do not form a complete, mutually ",
    "exclusive partition of the pooled resource.",
    call. = FALSE
  )
}

df.category.definition <- tribble(
  ~resource_group,
  ~required_evidence,
  ~excluded_or_not_required,
  ~allowed_interpretation,
  "pooled_loop_annotation_resource",
  "Distinct HiCCUPS loop with distance <2 Mb.",
  "No CTCF, promoter/TSS, or ATAC evidence is required.",
  "Pooled rat frontal cortex chromatin-loop annotation resource.",
  "putative_regulatory_promoter_TSS_nonTSS_ATAC_supported",
  "Promoter/TSS assignment within 200 kb and non-TSS ATAC support at the opposite anchor.",
  "Strong CTCF support is not required.",
  "Putative regulatory loop with layered annotation and accessibility support.",
  "structural_CTCF_supported_not_ATAC_regulatory",
  "At least 6 predicted CTCF motifs at both anchors.",
  "Does not satisfy promoter/TSS plus non-TSS ATAC criteria.",
  "Structurally supported loop; not direct evidence of enhancer function.",
  "lower_support_or_uncertain",
  "Member of pooled loop resource.",
  "Neither of the two higher-evidence category definitions is satisfied.",
  "Lower-support or uncertain loop retained for resource completeness.",
  "highlighted_conservative_subset",
  "Promoter/TSS assignment, opposite-anchor non-TSS ATAC, and CTCF >=6 at both anchors.",
  "This is nested within the putative-regulatory category.",
  "More conservative highlighted subset; not a separate fourth category."
)

################################################################################
# 9. CTCF threshold sensitivity
################################################################################

# Fixed thresholds are retained as a sensitivity analysis. Thresholds 12 and
# 30 are included because they are the anchor-width-scaled equivalents of the
# legacy >=6 threshold at 10-kb and 25-kb resolution, respectively.
ctcf.thresholds <- c(1L, 3L, 6L, 9L, 12L, 30L)

create_ctcf_sensitivity_row <- function(df, threshold, resolution.label) {
  passes.threshold <- (
    df$ctcf_count_anchor1 >= threshold &
      df$ctcf_count_anchor2 >= threshold
  )
  n.passes <- sum(passes.threshold)

  tibble(
    resolution = resolution.label,
    ctcf_min_threshold = threshold,
    n_loops_universe = nrow(df),
    n_passes_ctcf_thr = n.passes,
    pct_passes_ctcf_thr = round(100 * mean(passes.threshold), 1),
    n_ctcf_and_promoter_tss = sum(
      passes.threshold &
        df$promoter_or_tss_support
    ),
    pct_ctcf_and_promoter_tss_of_ctcf = round(
      100 * sum(
        passes.threshold &
          df$promoter_or_tss_support
      ) / max(n.passes, 1),
      1
    ),
    n_ctcf_promoter_tss_nonTSS_ATAC = sum(
      passes.threshold &
        df$putative_regulatory_support
    ),
    pct_ctcf_promoter_tss_nonTSS_ATAC_of_ctcf = round(
      100 * sum(
        passes.threshold &
          df$putative_regulatory_support
      ) / max(n.passes, 1),
      1
    ),
    n_putative_regulatory = sum(df$putative_regulatory_support)
  )
}

df.ctcf.sensitivity.by.threshold <- map_dfr(
  ctcf.thresholds,
  function(threshold.i) {
    bind_rows(
      map_dfr(
        sort(unique(df.loop.evidence$resolution)),
        function(resolution.i) {
          create_ctcf_sensitivity_row(
            df.loop.evidence %>%
              filter(resolution == resolution.i),
            threshold.i,
            resolution.i
          )
        }
      ),
      create_ctcf_sensitivity_row(
        df.loop.evidence,
        threshold.i,
        "ALL"
      )
    )
  }
) %>%
  arrange(ctcf_min_threshold, resolution)

df.ctcf.threshold.summary <- bind_rows(
  df.loop.evidence %>%
    group_by(resolution) %>%
    summarise(
      n_loops = n(),
      n_ctcf_any_both = sum(ctcf_both_anchors_any),
      pct_ctcf_any_both = percent_true(ctcf_both_anchors_any),
      n_ctcf_ge6_both = sum(passes_ctcf_ge6_both_anchors),
      pct_ctcf_ge6_both = percent_true(
        passes_ctcf_ge6_both_anchors
      ),
      median_min_ctcf_count = stats::median(
        ctcf_min_anchor_count
      ),
      q1_min_ctcf_count = stats::quantile(
        ctcf_min_anchor_count,
        0.25
      ),
      q3_min_ctcf_count = stats::quantile(
        ctcf_min_anchor_count,
        0.75
      ),
      .groups = "drop"
    ),
  df.loop.evidence %>%
    summarise(
      resolution = "ALL",
      n_loops = n(),
      n_ctcf_any_both = sum(ctcf_both_anchors_any),
      pct_ctcf_any_both = percent_true(ctcf_both_anchors_any),
      n_ctcf_ge6_both = sum(passes_ctcf_ge6_both_anchors),
      pct_ctcf_ge6_both = percent_true(
        passes_ctcf_ge6_both_anchors
      ),
      median_min_ctcf_count = stats::median(
        ctcf_min_anchor_count
      ),
      q1_min_ctcf_count = stats::quantile(
        ctcf_min_anchor_count,
        0.25
      ),
      q3_min_ctcf_count = stats::quantile(
        ctcf_min_anchor_count,
        0.75
      )
    )
) %>%
  arrange(resolution)

# The resolution-adjusted rule holds the motif density constant relative to
# the legacy 5-kb >=6 rule: 5 kb >=6, 10 kb >=12, and 25 kb >=30 motifs at
# both anchors. This is a technical sensitivity analysis of motif burden, not
# a biological threshold for CTCF binding or enhancer activity.
summarise_ctcf_resolution_adjustment <- function(df, resolution.label) {
  fixed.pass <- df$passes_ctcf_ge6_both_anchors
  adjusted.pass <- df$passes_ctcf_resolution_adjusted_both_anchors
  putative.pass <- df$putative_regulatory_support
  adjusted.rule <- if (resolution.label == "ALL") {
    "5K>=6; 10K>=12; 25K>=30"
  } else {
    paste0(
      resolution.label,
      ">=",
      unique(df$ctcf_resolution_adjusted_min_threshold)
    )
  }

  tibble(
    resolution = resolution.label,
    resolution_adjusted_rule = adjusted.rule,
    n_loops_universe = nrow(df),
    n_fixed_ge6_both = sum(fixed.pass),
    pct_fixed_ge6_both = round(100 * mean(fixed.pass), 1),
    n_resolution_adjusted_both = sum(adjusted.pass),
    pct_resolution_adjusted_both = round(
      100 * mean(adjusted.pass),
      1
    ),
    n_pass_both_rules = sum(fixed.pass & adjusted.pass),
    n_fixed_ge6_only = sum(fixed.pass & !adjusted.pass),
    n_resolution_adjusted_only = sum(!fixed.pass & adjusted.pass),
    n_fail_both_rules = sum(!fixed.pass & !adjusted.pass),
    n_putative_regulatory = sum(putative.pass),
    n_putative_regulatory_fixed_ge6 = sum(
      putative.pass & fixed.pass
    ),
    pct_putative_regulatory_fixed_ge6 = round(
      100 * sum(putative.pass & fixed.pass) /
        max(sum(putative.pass), 1),
      1
    ),
    n_putative_regulatory_resolution_adjusted = sum(
      putative.pass & adjusted.pass
    ),
    pct_putative_regulatory_resolution_adjusted = round(
      100 * sum(putative.pass & adjusted.pass) /
        max(sum(putative.pass), 1),
      1
    )
  )
}

df.ctcf.resolution.adjusted.summary <- bind_rows(
  map_dfr(
    sort(unique(df.loop.evidence$resolution)),
    function(resolution.i) {
      summarise_ctcf_resolution_adjustment(
        df.loop.evidence %>%
          filter(resolution == resolution.i),
        resolution.i
      )
    }
  ),
  summarise_ctcf_resolution_adjustment(
    df.loop.evidence,
    "ALL"
  )
) %>%
  arrange(resolution)

df.ctcf.resolution.adjusted.loop.comparison <- df.loop.evidence %>%
  transmute(
    loop_id,
    resolution,
    ctcf_count_anchor1,
    ctcf_count_anchor2,
    ctcf_min_anchor_count,
    fixed_min_threshold = 6L,
    passes_ctcf_ge6_both_anchors,
    resolution_adjusted_min_threshold =
      ctcf_resolution_adjusted_min_threshold,
    passes_ctcf_resolution_adjusted_both_anchors,
    putative_regulatory_support,
    fixed_vs_adjusted_status = case_when(
      passes_ctcf_ge6_both_anchors &
        passes_ctcf_resolution_adjusted_both_anchors ~
        "passes_both_rules",
      passes_ctcf_ge6_both_anchors &
        !passes_ctcf_resolution_adjusted_both_anchors ~
        "passes_fixed_ge6_only",
      !passes_ctcf_ge6_both_anchors &
        passes_ctcf_resolution_adjusted_both_anchors ~
        "passes_resolution_adjusted_only",
      TRUE ~ "fails_both_rules"
    )
  ) %>%
  arrange(resolution, loop_id)

################################################################################
# 10. ATAC support by resolution and paired-anchor comparison
################################################################################

df.atac.candidate <- df.loop.evidence %>%
  filter(passes_promoter_tss_200kb)

df.atac.support.by.resolution <- bind_rows(
  map_dfr(
    sort(unique(df.atac.candidate$resolution)),
    function(resolution.i) {
      summarise_atac_support(
        df.atac.candidate %>%
          filter(resolution == resolution.i),
        resolution.i
      )
    }
  ),
  summarise_atac_support(df.atac.candidate, "ALL")
) %>%
  arrange(resolution)

# Preserve the original one-row output name for downstream compatibility.
df.atac.support.summary <- df.atac.support.by.resolution %>%
  filter(resolution == "ALL") %>%
  dplyr::select(-resolution)

df.atac.paired.anchor.mcnemar <- bind_rows(
  map_dfr(
    sort(unique(df.atac.candidate$resolution)),
    function(resolution.i) {
      summarise_mcnemar(
        df.atac.candidate %>%
          filter(resolution == resolution.i),
        resolution.i
      )
    }
  ),
  summarise_mcnemar(df.atac.candidate, "ALL")
) %>%
  arrange(resolution)

################################################################################
# 11. Resource tables and downstream gene summaries
################################################################################

resource.columns <- c(
  "loop_id",
  "chr1", "start1", "end1",
  "chr2", "start2", "end2",
  "resolution", "loop_distance",
  "proposed_major_category",
  "proposed_detailed_category",
  "putative_regulatory_support",
  "strongest_candidate_support",
  "current_final_set",
  "current_final_category",
  "ctcf_count_anchor1",
  "ctcf_count_anchor2",
  "ctcf_min_anchor_count",
  "ctcf_both_anchors_any",
  "passes_ctcf_ge6_both_anchors",
  "promoter_or_tss_support",
  "promoter_support",
  "tss_support",
  "promoter_tss_component",
  "promoter_tss_distance",
  "promoter_anchor_side",
  "candidate_regulatory_anchor_side",
  "gene_symbol",
  "ensembl_gene_id",
  "atac_promoter_anchor",
  "atac_candidate_regulatory_anchor",
  "candidate_regulatory_anchor_has_non_tss_fragment",
  "candidate_regulatory_anchor_non_tss_bp",
  "atac_candidate_regulatory_anchor_non_tss"
)

resource.tables <- list(
  pooled_loop_annotation_resource = df.loop.evidence,
  highlighted_promoter_TSS_nonTSS_ATAC_CTCF_supported_subset =
    df.loop.evidence %>%
      filter(strongest_candidate_support),
  putative_regulatory_promoter_TSS_nonTSS_ATAC_loops =
    df.loop.evidence %>%
      filter(putative_regulatory_support),
  structural_CTCF_supported_loops =
    df.loop.evidence %>%
      filter(
        proposed_major_category ==
          "structural_CTCF_supported_not_ATAC_regulatory"
      ),
  lower_support_or_uncertain_loops =
    df.loop.evidence %>%
      filter(
        proposed_major_category ==
          "lower_support_or_uncertain"
      ),
  all_loop_evidence = df.loop.evidence,
  # Backward-compatible alias for the previous output name.
  main_strongest_candidate_PE_loops =
    df.loop.evidence %>%
      filter(strongest_candidate_support)
)

df.gene.count.by.set <- bind_rows(
  df.loop.evidence %>%
    filter(strongest_candidate_support) %>%
    count(
      gene_set = "main_strongest_candidate_PE_loops",
      gene_symbol,
      ensembl_gene_id,
      sort = TRUE
    ),
  df.loop.evidence %>%
    filter(putative_regulatory_support) %>%
    count(
      gene_set =
        "putative_regulatory_promoter_TSS_nonTSS_ATAC_loops",
      gene_symbol,
      ensembl_gene_id,
      sort = TRUE
    ),
  df.loop.evidence %>%
    filter(current_final_set) %>%
    count(
      gene_set = "current_final_CTCF_promoter_TSS_loops",
      gene_symbol,
      ensembl_gene_id,
      sort = TRUE
    ),
  df.loop.evidence %>%
    filter(
      proposed_major_category ==
        "structural_CTCF_supported_not_ATAC_regulatory"
    ) %>%
    count(
      gene_set = "structural_CTCF_supported_loops",
      gene_symbol,
      ensembl_gene_id,
      sort = TRUE
    )
) %>%
  filter(!is.na(gene_symbol) | !is.na(ensembl_gene_id)) %>%
  arrange(gene_set, desc(n), gene_symbol)

df.gene.count.threshold.summary <- df.gene.count.by.set %>%
  group_by(gene_set) %>%
  summarise(
    n_genes = n(),
    n_ensembl_genes = n_distinct(
      ensembl_gene_id,
      na.rm = TRUE
    ),
    max_interactions_per_gene = max(n),
    n_genes_ge_5 = sum(n >= 5),
    n_genes_ge_10 = sum(n >= 10),
    n_genes_ge_11 = sum(n >= 11),
    n_genes_ge_12 = sum(n >= 12),
    .groups = "drop"
  ) %>%
  arrange(gene_set)

go.gene.sets <- list(
  strongest_ge10 = prepare_go_gene_set(
    df.gene.count.by.set,
    "main_strongest_candidate_PE_loops",
    10
  ),
  strongest_ge11 = prepare_go_gene_set(
    df.gene.count.by.set,
    "main_strongest_candidate_PE_loops",
    11
  ),
  putative_ge10 = prepare_go_gene_set(
    df.gene.count.by.set,
    "putative_regulatory_promoter_TSS_nonTSS_ATAC_loops",
    10
  ),
  putative_ge11 = prepare_go_gene_set(
    df.gene.count.by.set,
    "putative_regulatory_promoter_TSS_nonTSS_ATAC_loops",
    11
  ),
  current_final_ge10 = prepare_go_gene_set(
    df.gene.count.by.set,
    "current_final_CTCF_promoter_TSS_loops",
    10
  ),
  current_final_ge11 = prepare_go_gene_set(
    df.gene.count.by.set,
    "current_final_CTCF_promoter_TSS_loops",
    11
  )
)

go.universe <- df.loop.evidence %>%
  filter(
    passes_promoter_tss_200kb,
    !is.na(ensembl_gene_id)
  ) %>%
  distinct(ensembl_gene_id) %>%
  pull(ensembl_gene_id)

df.go.result <- map_dfr(
  names(go.gene.sets),
  function(set.name) {
    run_go_enrichment(
      go.gene.sets[[set.name]],
      go.universe,
      set.name,
      output.dir
    )
  }
)

################################################################################
# 12. Optional loop-sharing vs. true genetic-distance analysis
#
# A panel-type label is not a genetic-distance matrix. This analysis runs only
# when a real long-format SNP/VCF-derived distance table is available.
################################################################################

df.genetic.distance.status <- tibble(
  requested_file = genetic.distance.file,
  file_exists = file.exists(genetic.distance.file),
  analysis_status = if_else(
    file.exists(genetic.distance.file),
    "ready_to_run",
    "skipped_no_true_genetic_distance_matrix"
  ),
  interpretation = if_else(
    file.exists(genetic.distance.file),
    paste0(
      "Pairwise loop Jaccard distance will be compared with the supplied ",
      "SNP/VCF-derived genetic distance."
    ),
    paste0(
      "No proxy was substituted. A real SNP/VCF-derived pairwise ",
      "distance matrix is required for this analysis."
    )
  )
)

df.pairwise.loop.genetic.distance <- tibble(
  resolution = character(),
  strain1 = character(),
  strain2 = character(),
  pair_strain_a = character(),
  pair_strain_b = character(),
  n_shared_loops = integer(),
  n_union_loops = integer(),
  jaccard_similarity = double(),
  loop_jaccard_distance = double(),
  genetic_distance = double()
)

df.loop.genetic.distance.correlation <- tibble(
  resolution = character(),
  n_pairs = integer(),
  spearman_rho = double(),
  spearman_p = double(),
  analysis_note = character()
)

if (file.exists(genetic.distance.file)) {
  df.genetic.distance <- read_genetic_distance_long(
    genetic.distance.file
  )

  df.pairwise.loop.genetic.distance <-
    df.pairwise.loop.overlap.by.strain %>%
      filter(
        strain1 != strain2,
        !is.na(jaccard_similarity)
      ) %>%
      mutate(
        pair_strain_a = pmin(strain1, strain2),
        pair_strain_b = pmax(strain1, strain2)
      ) %>%
      distinct(
        resolution,
        pair_strain_a,
        pair_strain_b,
        .keep_all = TRUE
      ) %>%
      left_join(
        df.genetic.distance %>%
          dplyr::select(
            pair_strain_a,
            pair_strain_b,
            genetic_distance
          ),
        by = c("pair_strain_a", "pair_strain_b")
      ) %>%
      mutate(loop_jaccard_distance = jaccard_distance) %>%
      dplyr::select(
        resolution,
        strain1,
        strain2,
        pair_strain_a,
        pair_strain_b,
        n_shared_loops,
        n_union_loops,
        jaccard_similarity,
        loop_jaccard_distance,
        genetic_distance
      )

  df.loop.genetic.distance.correlation <-
    df.pairwise.loop.genetic.distance %>%
      filter(
        !is.na(loop_jaccard_distance),
        !is.na(genetic_distance)
      ) %>%
      group_by(resolution) %>%
      group_modify(function(.x, .y) {
        correlation.result <- safe_correlation_summary(
          .x$loop_jaccard_distance,
          .x$genetic_distance,
          method = "spearman"
        )

        tibble(
          n_pairs = correlation.result$n,
          spearman_rho = correlation.result$estimate,
          spearman_p = correlation.result$p_value,
          analysis_note = paste0(
            "Exploratory exact-loop Jaccard distance vs. supplied ",
            "SNP/VCF-derived genetic distance."
          )
        )
      }) %>%
      ungroup()

  df.genetic.distance.status <- df.genetic.distance.status %>%
    mutate(
      analysis_status = "completed_with_supplied_true_distance_matrix",
      n_distance_pairs_supplied = nrow(df.genetic.distance),
      n_loop_pairs_with_distance = sum(
        !is.na(
          df.pairwise.loop.genetic.distance$genetic_distance
        )
      )
    )
} else {
  df.genetic.distance.status <- df.genetic.distance.status %>%
    mutate(
      n_distance_pairs_supplied = 0L,
      n_loop_pairs_with_distance = 0L
    )
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

################################################################################
# 13. Write resubmission outputs
################################################################################

readr::write_tsv(
  df.filtering.logic.audit,
  file.path(output.dir, "filtering_logic_audit.tsv")
)
readr::write_tsv(
  df.category.definition,
  file.path(output.dir, "loop_category_definitions.tsv")
)
readr::write_tsv(
  df.loop.evidence,
  file.path(output.dir, "loop_evidence_table.tsv")
)
readr::write_tsv(
  df.loop.evidence.summary,
  file.path(output.dir, "loop_evidence_summary.tsv")
)
readr::write_tsv(
  df.loop.category.by.resolution,
  file.path(output.dir, "loop_category_by_resolution.tsv")
)
readr::write_tsv(
  df.ctcf.threshold.summary,
  file.path(output.dir, "ctcf_threshold_summary.tsv")
)
readr::write_tsv(
  df.ctcf.sensitivity.by.threshold,
  file.path(output.dir, "ctcf_sensitivity_by_threshold.tsv")
)
readr::write_tsv(
  df.ctcf.resolution.adjusted.summary,
  file.path(output.dir, "ctcf_resolution_adjusted_summary.tsv")
)
readr::write_tsv(
  df.ctcf.resolution.adjusted.loop.comparison,
  file.path(output.dir, "ctcf_resolution_adjusted_loop_comparison.tsv")
)
readr::write_tsv(
  df.atac.support.summary,
  file.path(output.dir, "atac_support_summary.tsv")
)
readr::write_tsv(
  df.atac.support.by.resolution,
  file.path(output.dir, "atac_support_by_resolution.tsv")
)
readr::write_tsv(
  df.atac.paired.anchor.mcnemar,
  file.path(output.dir, "atac_paired_anchor_mcnemar.tsv")
)
readr::write_tsv(
  df.sample.strain.key,
  file.path(output.dir, "hrdp_sample_strain_key.tsv")
)
readr::write_tsv(
  df.strain.loop.count.by.resolution,
  file.path(output.dir, "strain_loop_counts_by_resolution.tsv")
)
readr::write_tsv(
  df.depth.qc.by.strain,
  file.path(output.dir, "depth_qc_by_strain.tsv")
)
readr::write_tsv(
  df.depth.qc.summary,
  file.path(output.dir, "depth_qc_summary.tsv")
)
readr::write_tsv(
  df.depth.loop.correlation.by.resolution,
  file.path(output.dir, "depth_loop_correlation_by_resolution.tsv")
)
readr::write_tsv(
  df.depth.analysis.interpretation,
  file.path(output.dir, "depth_analysis_interpretation.tsv")
)
readr::write_tsv(
  df.strain.loop.presence.matrix,
  file.path(output.dir, "strain_loop_presence_matrix.tsv")
)
readr::write_tsv(
  df.shared.loop.pair.by.strain,
  file.path(output.dir, "shared_loop_pairs_by_strain.tsv")
)
readr::write_tsv(
  df.pairwise.loop.overlap.by.strain,
  file.path(output.dir, "pairwise_loop_overlap_by_strain.tsv")
)
readr::write_tsv(
  df.loop.sharing.distribution.by.strain,
  file.path(output.dir, "loop_sharing_distribution_by_strain.tsv")
)
readr::write_tsv(
  df.genetic.distance.status,
  file.path(output.dir, "genetic_distance_analysis_status.tsv")
)
readr::write_tsv(
  df.pairwise.loop.genetic.distance,
  file.path(output.dir, "pairwise_loop_genetic_distance.tsv")
)
readr::write_tsv(
  df.loop.genetic.distance.correlation,
  file.path(output.dir, "loop_genetic_distance_correlation.tsv")
)
readr::write_tsv(
  df.genetic.proxy.deprecation,
  file.path(output.dir, "pairwise_loop_genetic_proxy.tsv")
)
readr::write_tsv(
  df.genetic.proxy.deprecation,
  file.path(output.dir, "loop_genetic_distance_panel_summary.tsv")
)
readr::write_tsv(
  df.gene.count.by.set,
  file.path(output.dir, "top_gene_summary_by_category.tsv")
)
readr::write_tsv(
  df.gene.count.threshold.summary,
  file.path(output.dir, "top_gene_threshold_summary.tsv")
)

if (nrow(df.go.result) > 0) {
  readr::write_tsv(
    df.go.result,
    file.path(output.dir, "go_enrichment_BP_all_sets.tsv")
  )
}

for (resource.name in names(resource.tables)) {
  readr::write_tsv(
    resource.tables[[resource.name]] %>%
      dplyr::select(any_of(resource.columns)),
    file.path(output.dir, paste0(resource.name, ".tsv"))
  )
}

for (set.name in names(go.gene.sets)) {
  readr::write_tsv(
    go.gene.sets[[set.name]],
    file.path(
      output.dir,
      paste0("go_input_", set.name, "_genes.tsv")
    )
  )
  readr::write_lines(
    go.gene.sets[[set.name]]$ensembl_gene_id,
    file.path(
      output.dir,
      paste0("go_input_", set.name, "_ensembl_ids.txt")
    )
  )
}

generated.output.files <- c(
  "filtering_logic_audit.tsv",
  "loop_category_definitions.tsv",
  "loop_evidence_table.tsv",
  "loop_evidence_summary.tsv",
  "loop_category_by_resolution.tsv",
  "ctcf_threshold_summary.tsv",
  "ctcf_sensitivity_by_threshold.tsv",
  "ctcf_resolution_adjusted_summary.tsv",
  "ctcf_resolution_adjusted_loop_comparison.tsv",
  "atac_support_summary.tsv",
  "atac_support_by_resolution.tsv",
  "atac_paired_anchor_mcnemar.tsv",
  "hrdp_sample_strain_key.tsv",
  "strain_loop_counts_by_resolution.tsv",
  "depth_qc_by_strain.tsv",
  "depth_qc_summary.tsv",
  "depth_loop_correlation_by_resolution.tsv",
  "depth_analysis_interpretation.tsv",
  "strain_loop_presence_matrix.tsv",
  "shared_loop_pairs_by_strain.tsv",
  "pairwise_loop_overlap_by_strain.tsv",
  "loop_sharing_distribution_by_strain.tsv",
  "genetic_distance_analysis_status.tsv",
  "pairwise_loop_genetic_distance.tsv",
  "loop_genetic_distance_correlation.tsv",
  "pairwise_loop_genetic_proxy.tsv",
  "loop_genetic_distance_panel_summary.tsv",
  "top_gene_summary_by_category.tsv",
  "top_gene_threshold_summary.tsv",
  if (nrow(df.go.result) > 0) "go_enrichment_BP_all_sets.tsv",
  paste0(names(resource.tables), ".tsv"),
  paste0("go_input_", names(go.gene.sets), "_genes.tsv"),
  paste0("go_input_", names(go.gene.sets), "_ensembl_ids.txt")
)

df.output.manifest <- tibble(
  output_file = generated.output.files,
  output_path = file.path(output.dir, generated.output.files),
  file_exists = file.exists(output_path),
  generated_by = "promoter_enhancer_interaction_resubmit.R"
)

readr::write_tsv(
  df.output.manifest,
  file.path(output.dir, "resubmit_output_manifest.tsv")
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output.dir, "resubmit_session_info.txt")
)

################################################################################
# 14. Final checks and console summary
################################################################################

message("\nWrote outputs to: ", output.dir)
message("Pooled loop resource (<2 Mb): ", nrow(df.loop.universe))
message("Previous final CTCF + promoter/TSS set: ", n.current.final)
message("Promoter/TSS candidate loops: ", n.promoter.tss)
message(
  "Putative regulatory loops with non-TSS ATAC support: ",
  n.candidate.non.tss.atac
)
message(
  "Structural CTCF-supported loops under the previous >=6 threshold: ",
  n.ctcf.ge6.both
)
message(
  "Highlighted promoter/TSS + non-TSS ATAC + CTCF >=6 subset: ",
  n.strongest.candidate
)

message("\nThree mutually exclusive major categories:")
print(
  df.loop.evidence %>%
    count(proposed_major_category, sort = TRUE)
)

message("\nCategory counts by resolution:")
print(df.loop.category.by.resolution)

message("\nCTCF threshold sensitivity, pooled across resolutions:")
print(
  df.ctcf.sensitivity.by.threshold %>%
    filter(resolution == "ALL") %>%
    dplyr::select(
      ctcf_min_threshold,
      n_passes_ctcf_thr,
      pct_passes_ctcf_thr,
      n_ctcf_promoter_tss_nonTSS_ATAC,
      n_putative_regulatory
    )
)

message("\nFixed >=6 vs. resolution-adjusted CTCF motif thresholds:")
print(df.ctcf.resolution.adjusted.summary)

message("\nATAC support by resolution:")
print(df.atac.support.by.resolution)

message("\nPaired promoter vs. candidate-regulatory ATAC comparison:")
print(df.atac.paired.anchor.mcnemar)

message("\nDepth sensitivity: loop count vs. Hi-C contacts")
print(
  df.depth.loop.correlation.by.resolution %>%
    filter(depth_metric == "Hi-C_Contacts") %>%
    dplyr::select(
      resolution,
      n_strains,
      pearson_r,
      pearson_p,
      spearman_rho,
      spearman_p
    )
)

message("\nGenetic-distance analysis status:")
print(df.genetic.distance.status)

message("\nTop-gene threshold summary:")
print(df.gene.count.threshold.summary)

if (nrow(df.go.result) > 0) {
  message("\nGO enrichment rows written: ", nrow(df.go.result))
}
