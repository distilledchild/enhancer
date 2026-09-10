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

# Run the optional loop-sharing/genetic-distance comparison only when a true
# pairwise distance table is explicitly supplied or present in the analysis folder.
genetic.distance.file <- Sys.getenv(
  "HRDP_GENETIC_DISTANCE_FILE",
  unset = file.path(analysis.dir, "hrdp_genetic_distance.tsv")
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
  left_join(df.loop.distinct.2mb %>% dplyr::select(loop_id, chr1, start1, end1, chr2, start2, end2, loop_distance), by = "loop_id") %>%
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
  "putative_regulatory" = "#00BA38",
  "promoter_associated_without_distal_ATAC_support" = "#F8766D",
  "no_direct_promoter_TSS" = "#619CFF"
)
category.labels <- c(
  "putative_regulatory" = "Putative regulatory",
  "promoter_associated_without_distal_ATAC_support" =
    "Promoter-associated without\ndistal ATAC support",
  "no_direct_promoter_TSS" = "No direct promoter/TSS"
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
df.figure2.chromosome.resolution <- df.loop.distinct.2mb %>%
  count(chr1, resolution, name = "n_loops") %>%
  mutate(
    chromosome_label = str_remove(chr1, "^chr"),
    chromosome_order = case_when(
      chromosome_label == "X" ~ 100,
      chromosome_label == "Y" ~ 101,
      chromosome_label %in% c("M", "MT") ~ 102,
      TRUE ~ suppressWarnings(as.numeric(chromosome_label))
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
df.figure4.category.summary <- df.loop.category.summary %>%
  mutate(
    category_label = recode(revised_major_category, !!!category.labels),
    category_label = forcats::fct_reorder(category_label, n_loops),
    count_label = str_c(scales::comma(n_loops), " (", pct_pooled_loops, "%)")
  )
df.figure4.category.by.resolution <- df.loop.category.by.resolution %>%
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
  "promoter_tss_anchor_evidence_definition_summary.tsv", "df.promoter.anchor.assignment.definitions",
  "direct_promoter_tss_summary.tsv", "df.direct.promoter.tss.summary",
  "dual_promoter_atac_directional_category_definitions.tsv", "df.dual.promoter.atac.directional.category.definition",
  "dual_promoter_atac_directional_summary.tsv", "df.dual.promoter.atac.directional.summary",
  "revised_atac_support_by_resolution.tsv", "df.atac.support.by.resolution",
  "revised_atac_paired_anchor_mcnemar.tsv", "df.atac.paired.anchor.mcnemar",
  "revised_atac_analysis_definition.tsv", "df.atac.analysis.definition",
  "transcript_position_summary.tsv", "df.transcript.position.summary",
  "transcript_position_definitions.tsv", "df.transcript.position.definition",
  "revised_assignment_pipeline_status.tsv", "df.revised.assignment.pipeline.status",
  "revised_loop_category_definitions.tsv", "df.loop.category.definition",
  "revised_loop_category_summary.tsv", "df.loop.category.summary",
  "revised_loop_detailed_category_summary.tsv", "df.loop.detailed.category.summary",
  "revised_loop_category_by_resolution.tsv", "df.loop.category.by.resolution",
  "table_s3_multiple_interaction_genes_putative_regulatory.tsv", "df.table.s3.putative",
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
selected.resource.table.names <- c(
  "revised_pooled_loop_annotation_resource",
  "revised_direct_loop_gene_assignments",
  "revised_putative_regulatory_loops",
  "revised_promoter_associated_without_distal_atac_loops",
  "revised_no_direct_promoter_tss_loops"
)
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
message("Pooled loop resource (<2 Mb): ", nrow(df.loop.distinct.2mb))
message("Revised direct promoter/TSS-supported loops: ", sum(df.direct.promoter.tss.loop.summary$has_any_direct_promoter_tss))
message("Revised direct loop-anchor-gene assignments: ", nrow(df.direct.promoter.tss.gene.assignment))
message("Revised proximal promoter/TSS-supported loops (1-200 kb): ", n_distinct(df.proximal.promoter.tss.gene.assignment$loop_id))
message("Revised proximal loop-anchor-gene assignments: ", nrow(df.proximal.promoter.tss.gene.assignment))
message("Revised secondary inward <=10-kb supported loops: ", n_distinct(df.secondary.inward.proximal.gene.assignment$loop_id))
message("Revised secondary inward <=10-kb loop-anchor-gene assignments: ", nrow(df.secondary.inward.proximal.gene.assignment))
message("Revised direct/proximal assignment, ATAC, transcript-position, and loop-category sections are complete. Legacy analyses are excluded from this production workflow.")

message("\nRevised true-TSS-excluded ATAC support by resolution:")
print(df.atac.support.by.resolution)
message("\nRevised paired-anchor ATAC comparisons:")
print(df.atac.paired.anchor.mcnemar)
message("\nRevised transcript-position summary:")
print(df.transcript.position.summary)
message("\nRevised mutually exclusive loop categories:")
print(df.loop.category.summary)
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
