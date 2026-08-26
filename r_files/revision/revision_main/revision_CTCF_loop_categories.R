# lintr: disable
################################################################################
# Resubmission Analysis: Predicted CTCF Motifs Across Revised Loop Categories
#
# Purpose: Compare q-value-supported predicted CTCF motif intervals across the
#          three revised loop categories without using CTCF to select or classify
#          regulatory loops.
#
# Primary unit: Exact-coordinate-distinct predicted motif intervals with FIMO
#               q-value <= 0.05. All FIMO predictions are retained only as a
#               sensitivity analysis.
################################################################################

if (basename(getwd()) != "enhancer" && dir.exists("enhancer")) {
  setwd("./enhancer")
}

funcs.file <- "./funcs_enhancer.R"
source(funcs.file)
list2env(resolve_enhancer_analysis_paths(funcs.file), envir = environment())

library("tidyverse")
library("GenomicRanges")
library("cowplot")
library("scales")

options(tibble.width = Inf, tibble.print_max = Inf, scipen = 999)

results.dir <- file.path(analysis.dir, "results")
dir.create(results.dir, recursive = TRUE, showWarnings = FALSE)

################################################################################
# 1. Load the official pooled resource and coordinate-normalized CTCF cache
################################################################################

resource.file <- file.path(results.dir, "revised_pooled_loop_annotation_resource.tsv")
ctcf.cache.file <- file.path(coord.cache.dir, "gr.ctcf.motif.rds")

if (!file.exists(resource.file) || !file.exists(ctcf.cache.file)) {
  stop("Missing pooled loop resource or CTCF coordinate cache.", call. = FALSE)
}

df.loop.resource <- read_tsv(resource.file, show_col_types = FALSE)
gr.ctcf.all <- readRDS(ctcf.cache.file)

required.loop.columns <- c(
  "loop_id", "chr1", "start1", "end1", "chr2", "start2", "end2",
  "resolution", "direct_anchor_assignment_class",
  "dual_promoter_directional_category", "revised_putative_regulatory_support",
  "revised_promoter_promoter_compatible",
  "revised_single_promoter_without_opposite_atac",
  "revised_no_direct_promoter_tss",
  "predicted_ctcf_motif_interval_count_anchor1",
  "predicted_ctcf_motif_interval_count_anchor2"
)
check_required_columns(df.loop.resource, required.loop.columns, "Pooled loop resource")

if (!"best_q_value" %in% names(mcols(gr.ctcf.all))) {
  stop("The CTCF cache does not retain FIMO q-values.", call. = FALSE)
}

assert_analysis_condition(
  nrow(df.loop.resource) == 31021L && !anyDuplicated(df.loop.resource$loop_id),
  "The pooled resource must contain 31,021 unique exact loop-call records."
)

# FIMO q-values are monotonic with p-values; the cached representative for each
# exact interval is the prediction with the smallest p-value.
gr.ctcf.q05 <- gr.ctcf.all[
  !is.na(mcols(gr.ctcf.all)$best_q_value) &
    mcols(gr.ctcf.all)$best_q_value <= 0.05
]

message("All exact-coordinate-distinct predicted intervals: ", length(gr.ctcf.all))
message("Primary q <= 0.05 predicted intervals: ", length(gr.ctcf.q05))

################################################################################
# 2. Recalculate anchor-level CTCF annotation and revised loop categories
################################################################################

gr.anchor1 <- GRanges(
  seqnames = df.loop.resource$chr1,
  ranges = IRanges(df.loop.resource$start1, df.loop.resource$end1)
)
gr.anchor2 <- GRanges(
  seqnames = df.loop.resource$chr2,
  ranges = IRanges(df.loop.resource$start2, df.loop.resource$end2)
)

ctcf.all.anchor1 <- countOverlaps(gr.anchor1, gr.ctcf.all, ignore.strand = TRUE)
ctcf.all.anchor2 <- countOverlaps(gr.anchor2, gr.ctcf.all, ignore.strand = TRUE)

assert_analysis_condition(
  identical(
    as.integer(ctcf.all.anchor1),
    as.integer(df.loop.resource$predicted_ctcf_motif_interval_count_anchor1)
  ) && identical(
    as.integer(ctcf.all.anchor2),
    as.integer(df.loop.resource$predicted_ctcf_motif_interval_count_anchor2)
  ),
  "Recomputed unfiltered CTCF counts do not match the official resource."
)

df.ctcf.master <- df.loop.resource %>%
  mutate(
    ctcf_q05_anchor1 = as.integer(countOverlaps(gr.anchor1, gr.ctcf.q05, ignore.strand = TRUE)),
    ctcf_q05_anchor2 = as.integer(countOverlaps(gr.anchor2, gr.ctcf.q05, ignore.strand = TRUE)),
    total_ctcf_q05 = ctcf_q05_anchor1 + ctcf_q05_anchor2,
    both_anchors_ctcf_q05 = ctcf_q05_anchor1 > 0L & ctcf_q05_anchor2 > 0L,
    total_ctcf_all = ctcf.all.anchor1 + ctcf.all.anchor2,
    both_anchors_ctcf_all = ctcf.all.anchor1 > 0L & ctcf.all.anchor2 > 0L,
    sub_category_id = case_when(
      revised_putative_regulatory_support ~ "single_promoter_distal_atac",
      revised_single_promoter_without_opposite_atac ~ "single_promoter_no_distal_atac",
      revised_promoter_promoter_compatible &
        dual_promoter_directional_category !=
          "neither_direction_supported_promoter_promoter_or_lower_support" ~
        "promoter_both_anchors_with_distal_atac",
      revised_promoter_promoter_compatible ~
        "promoter_both_anchors_no_directional_distal_atac",
      revised_no_direct_promoter_tss ~ "no_direct_promoter_tss",
      TRUE ~ NA_character_
    ),
    major_category_id = case_when(
      sub_category_id %in% c(
        "single_promoter_distal_atac",
        "promoter_both_anchors_with_distal_atac"
      ) ~ "putative_regulatory",
      sub_category_id %in% c(
        "single_promoter_no_distal_atac",
        "promoter_both_anchors_no_directional_distal_atac"
      ) ~ "promoter_associated_no_distal_atac",
      sub_category_id == "no_direct_promoter_tss" ~ "no_direct_promoter_tss",
      TRUE ~ NA_character_
    )
  )

assert_analysis_condition(
  !anyNA(df.ctcf.master$major_category_id) &&
    sum(df.ctcf.master$major_category_id == "putative_regulatory") == 13376L &&
    sum(df.ctcf.master$major_category_id == "promoter_associated_no_distal_atac") == 2967L &&
    sum(df.ctcf.master$major_category_id == "no_direct_promoter_tss") == 14678L,
  "Revised loop categories are incomplete or do not match the official counts."
)

major.labels <- c(
  putative_regulatory = "1. Putative regulatory loops\nwith distal non-TSS ATAC support",
  promoter_associated_no_distal_atac = "2. Promoter-associated loops\nwithout distal ATAC support",
  no_direct_promoter_tss = "3. Loops without direct\npromoter/TSS overlap"
)
sub.labels <- c(
  single_promoter_distal_atac = "Single-promoter loops\nwith distal non-TSS ATAC",
  single_promoter_no_distal_atac = "Single-promoter loops\nwithout distal ATAC",
  promoter_both_anchors_with_distal_atac = "Promoter/TSS at both anchors\nwith directional distal ATAC",
  promoter_both_anchors_no_directional_distal_atac = "Promoter/TSS at both anchors\nwithout directional distal ATAC",
  no_direct_promoter_tss = "No direct promoter/TSS overlap"
)

df.ctcf.master <- df.ctcf.master %>%
  mutate(
    major_category = factor(
      unname(major.labels[major_category_id]),
      levels = unname(major.labels)
    ),
    sub_category = factor(
      unname(sub.labels[sub_category_id]),
      levels = unname(sub.labels)
    ),
    resolution = factor(resolution, levels = c("5K", "10K", "25K"))
  )

################################################################################
# 3. Build a correctly oriented anchor-level table (N = 62,042)
################################################################################

df.anchor.level <- bind_rows(
  df.ctcf.master %>%
    transmute(
      loop_id, resolution, major_category, sub_category, sub_category_id,
      anchor_side = "anchor1", ctcf_count_q05 = ctcf_q05_anchor1,
      promoter_tss_anchor = direct_anchor_assignment_class %in%
        c("direct_anchor1_only", "direct_both_anchors")
    ),
  df.ctcf.master %>%
    transmute(
      loop_id, resolution, major_category, sub_category, sub_category_id,
      anchor_side = "anchor2", ctcf_count_q05 = ctcf_q05_anchor2,
      promoter_tss_anchor = direct_anchor_assignment_class %in%
        c("direct_anchor2_only", "direct_both_anchors")
    )
) %>%
  mutate(
    anchor_role = case_when(
      sub_category_id == "single_promoter_distal_atac" & promoter_tss_anchor ~
        "Promoter/TSS anchor\n(putative regulatory)",
      sub_category_id == "single_promoter_distal_atac" ~
        "Distal non-TSS ATAC anchor\n(putative regulatory)",
      sub_category_id == "single_promoter_no_distal_atac" & promoter_tss_anchor ~
        "Promoter/TSS anchor\n(no distal ATAC)",
      sub_category_id == "single_promoter_no_distal_atac" ~
        "Opposite anchor\n(no distal ATAC)",
      sub_category_id == "promoter_both_anchors_with_distal_atac" ~
        "Promoter/TSS anchor\n(promoter at both anchors; distal ATAC)",
      sub_category_id == "promoter_both_anchors_no_directional_distal_atac" ~
        "Promoter/TSS anchor\n(promoter at both anchors; no directional ATAC)",
      TRUE ~ "Anchor without direct\npromoter/TSS overlap"
    )
  )

assert_analysis_condition(
  nrow(df.anchor.level) == 62042L &&
    sum(df.anchor.level$promoter_tss_anchor) == 12295L + 2L * 4048L,
  "Anchor-level reconstruction did not preserve all anchors or promoter sides."
)

################################################################################
# 4. Summary tables and statistical tests
################################################################################

summarise_ctcf <- function(data, count.column, both.column) {
  count.vector <- data[[count.column]]
  both.vector <- data[[both.column]]

  tibble(
    n_loops = nrow(data),
    n_both_anchors = sum(both.vector),
    pct_both_anchors = 100 * mean(both.vector),
    mean_intervals_per_loop = mean(count.vector),
    sd_intervals_per_loop = sd(count.vector),
    median_intervals_per_loop = median(count.vector),
    q25_intervals_per_loop = quantile(count.vector, 0.25),
    q75_intervals_per_loop = quantile(count.vector, 0.75)
  )
}

df.summary.major <- df.ctcf.master %>%
  group_by(major_category_id, major_category) %>%
  group_modify(~summarise_ctcf(.x, "total_ctcf_q05", "both_anchors_ctcf_q05")) %>%
  ungroup() %>%
  mutate(pct_loops = 100 * n_loops / sum(n_loops))

df.summary.by.resolution <- df.ctcf.master %>%
  group_by(resolution, major_category_id, major_category) %>%
  group_modify(~summarise_ctcf(.x, "total_ctcf_q05", "both_anchors_ctcf_q05")) %>%
  ungroup()

df.summary.subcategory <- df.ctcf.master %>%
  group_by(sub_category_id, sub_category) %>%
  group_modify(~summarise_ctcf(.x, "total_ctcf_q05", "both_anchors_ctcf_q05")) %>%
  ungroup() %>%
  mutate(pct_loops = 100 * n_loops / sum(n_loops))

df.summary.all.predictions <- df.ctcf.master %>%
  group_by(resolution, major_category_id, major_category) %>%
  group_modify(~summarise_ctcf(.x, "total_ctcf_all", "both_anchors_ctcf_all")) %>%
  ungroup()

run_global_tests <- function(data, scope) {
  kw <- kruskal.test(total_ctcf_q05 ~ major_category_id, data = data)
  chi <- chisq.test(table(data$major_category_id, data$both_anchors_ctcf_q05))

  tibble(
    scope = scope,
    count_test = "Kruskal-Wallis",
    count_statistic = unname(kw$statistic),
    count_df = unname(kw$parameter),
    count_p_value = kw$p.value,
    both_anchor_test = "Pearson chi-squared",
    both_anchor_statistic = unname(chi$statistic),
    both_anchor_df = unname(chi$parameter),
    both_anchor_p_value = chi$p.value
  )
}

run_pairwise_count_tests <- function(data, scope) {
  pairs <- combn(unique(as.character(data$major_category_id)), 2L, simplify = FALSE)

  map_dfr(pairs, function(pair) {
    x <- data$total_ctcf_q05[data$major_category_id == pair[1]]
    y <- data$total_ctcf_q05[data$major_category_id == pair[2]]
    test <- wilcox.test(x, y, exact = FALSE)
    u <- as.numeric(test$statistic)

    tibble(
      scope = scope,
      category_1 = pair[1], category_2 = pair[2],
      n_1 = length(x), n_2 = length(y),
      median_1 = median(x), median_2 = median(y),
      median_difference = median(x) - median(y),
      rank_biserial = 2 * u / (length(x) * length(y)) - 1,
      p_value = test$p.value
    )
  }) %>%
    mutate(p_adjust_bh = p.adjust(p_value, method = "BH"))
}

run_pairwise_binary_tests <- function(data, scope) {
  pairs <- combn(unique(as.character(data$major_category_id)), 2L, simplify = FALSE)

  map_dfr(pairs, function(pair) {
    first <- data$both_anchors_ctcf_q05[data$major_category_id == pair[1]]
    second <- data$both_anchors_ctcf_q05[data$major_category_id == pair[2]]
    contingency <- matrix(
      c(sum(first), sum(!first), sum(second), sum(!second)),
      nrow = 2L,
      byrow = TRUE
    )
    test <- fisher.test(contingency)

    tibble(
      scope = scope,
      category_1 = pair[1], category_2 = pair[2],
      n_1 = length(first), n_2 = length(second),
      pct_both_1 = 100 * mean(first), pct_both_2 = 100 * mean(second),
      percentage_point_difference = 100 * (mean(first) - mean(second)),
      odds_ratio = unname(test$estimate),
      odds_ratio_ci_low = test$conf.int[1],
      odds_ratio_ci_high = test$conf.int[2],
      p_value = test$p.value
    )
  }) %>%
    mutate(p_adjust_bh = p.adjust(p_value, method = "BH"))
}

analysis.scopes <- c(list(All = df.ctcf.master), split(df.ctcf.master, df.ctcf.master$resolution))
df.global.tests <- imap_dfr(analysis.scopes, run_global_tests)
df.pairwise.count.tests <- imap_dfr(analysis.scopes, run_pairwise_count_tests)
df.pairwise.binary.tests <- imap_dfr(analysis.scopes, run_pairwise_binary_tests)

output.tables <- list(
  revision_ctcf_major_category_summary = df.summary.major,
  revision_ctcf_major_category_summary_by_resolution = df.summary.by.resolution,
  revision_ctcf_subcategory_summary = df.summary.subcategory,
  revision_ctcf_all_predictions_sensitivity_by_resolution = df.summary.all.predictions,
  revision_ctcf_global_tests = df.global.tests,
  revision_ctcf_pairwise_count_tests = df.pairwise.count.tests,
  revision_ctcf_pairwise_both_anchor_tests = df.pairwise.binary.tests
)

# Keep line breaks in plot labels, but flatten them before writing machine-readable TSVs.
write_analysis_tsv <- function(data, name) {
  data %>%
    mutate(
      across(where(is.factor), as.character),
      across(where(is.character), ~str_replace_all(.x, "[\\r\\n]+", " "))
    ) %>%
    write_tsv(file.path(results.dir, paste0(name, ".tsv")))
}

iwalk(output.tables, write_analysis_tsv)

################################################################################
# 5. Publication figures based on q <= 0.05 predicted motif intervals
################################################################################

cat.colors <- c(
  "1. Putative regulatory loops\nwith distal non-TSS ATAC support" = "#F8766D",
  "2. Promoter-associated loops\nwithout distal ATAC support" = "#00BA38",
  "3. Loops without direct\npromoter/TSS overlap" = "#619CFF"
)

theme.publication <- function(base.size = 10) {
  theme_classic(base_size = base.size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 6)),
      plot.subtitle = element_text(color = "gray25", hjust = 0.5, margin = margin(b = 8)),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.grid.major.y = element_line(color = "gray92", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.margin = margin(10, 10, 10, 10)
    )
}

p.panel.a <- ggplot(df.ctcf.master, aes(major_category, total_ctcf_q05, fill = major_category)) +
  geom_violin(alpha = 0.30, width = 0.82, trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.22, fill = "white", outlier.shape = NA, linewidth = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.8, fill = "#DC2626") +
  scale_fill_manual(values = cat.colors) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 5, 10, 25, 50, 100, 200)) +
  labs(
    title = "A. Predicted CTCF motif intervals per loop",
    subtitle = "Exact-coordinate-distinct FIMO intervals with q <= 0.05; pooled descriptive comparison",
    x = NULL,
    y = "Intervals per loop (log1p scale)"
  ) +
  theme.publication() +
  theme(axis.text.x = element_text(size = 7.5, lineheight = 0.95))

p.panel.b <- ggplot(df.ctcf.master, aes(resolution, total_ctcf_q05, fill = major_category)) +
  geom_boxplot(
    width = 0.72,
    position = position_dodge2(width = 0.78, preserve = "single"),
    outlier.shape = NA,
    alpha = 0.78,
    linewidth = 0.45
  ) +
  scale_fill_manual(values = cat.colors) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 5, 10, 25, 50, 100, 200)) +
  labs(
    title = "B. Resolution-stratified motif interval counts",
    subtitle = "Resolution is shown explicitly because anchor width changes overlap opportunity",
    x = "HiCCUPS resolution",
    y = "Intervals per loop (log1p scale)"
  ) +
  theme.publication()

df.both.plot <- df.summary.by.resolution %>%
  mutate(label = sprintf("%.1f%%", pct_both_anchors))

p.panel.c <- ggplot(df.both.plot, aes(resolution, pct_both_anchors, fill = major_category)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68, color = "black", linewidth = 0.35) +
  geom_text(aes(label = label), position = position_dodge(width = 0.75), vjust = -0.35, size = 2.7) +
  scale_fill_manual(values = cat.colors) +
  scale_y_continuous(limits = c(0, 100), labels = label_percent(scale = 1), expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "C. Motif support at both anchors",
    subtitle = "At least one q <= 0.05 predicted motif interval at each anchor",
    x = "HiCCUPS resolution",
    y = "Loop records with both-anchor support"
  ) +
  theme.publication() +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 7))

role.levels <- c(
  "Promoter/TSS anchor\n(putative regulatory)",
  "Distal non-TSS ATAC anchor\n(putative regulatory)",
  "Promoter/TSS anchor\n(no distal ATAC)",
  "Opposite anchor\n(no distal ATAC)",
  "Promoter/TSS anchor\n(promoter at both anchors; distal ATAC)",
  "Promoter/TSS anchor\n(promoter at both anchors; no directional ATAC)",
  "Anchor without direct\npromoter/TSS overlap"
)
df.anchor.level$anchor_role <- factor(df.anchor.level$anchor_role, levels = role.levels)

p.panel.d <- ggplot(df.anchor.level, aes(anchor_role, ctcf_count_q05, fill = anchor_role)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.82, linewidth = 0.45) +
  scale_fill_manual(values = setNames(hue_pal()(length(role.levels)), role.levels)) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 2, 5, 10, 25, 50, 100)) +
  labs(
    title = "D. Anchor-level motif annotation by anchor role",
    subtitle = "All 62,042 anchors are retained; promoter sides follow the observed anchor assignment",
    x = NULL,
    y = "Intervals per anchor (log1p scale)"
  ) +
  theme.publication() +
  theme(axis.text.x = element_text(angle = 28, hjust = 1, size = 6.8, lineheight = 0.95))

p.master <- plot_grid(
  p.panel.a, p.panel.b,
  p.panel.c, p.panel.d,
  ncol = 2,
  align = "hv",
  axis = "tblr"
)

output.png <- file.path(results.dir, "revision_figure_ctcf_loop_categories.png")
output.pdf <- file.path(results.dir, "revision_figure_ctcf_loop_categories.pdf")

ggsave(output.png, p.master, width = 15, height = 11.5, dpi = 300, bg = "white")
ggsave(output.pdf, p.master, width = 15, height = 11.5, bg = "white")

message("Saved q <= 0.05 CTCF category analysis and figures to: ", results.dir)
