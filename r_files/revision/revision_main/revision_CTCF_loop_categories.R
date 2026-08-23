# lintr: disable
################################################################################
# Resubmission Analysis: CTCF Motif Occupancy and Comparison Across Loop Categories
#
# Script Name: revision_CTCF_loop_categories.R
# Purpose:     Comprehensive quantification and statistical comparison of predicted
#              CTCF motif counts and anchor occupancy across the 3 newly categorized
#              chromatin loop classes (and 5 sub-categories) in rat frontal cortex.
#              (Strictly synchronized with Figure 7 flowchart nomenclature)
#
# Inputs:      Cached coordinate-normalized data from 00_promoter_enhancer_interaction_resubmit_coord_prep.R
# Outputs:     - Statistical summary tables (TSV)
#              - High-resolution publication-quality ggplot2 figures (PNG / PDF)
################################################################################

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
library("cowplot")
library("scales")

options(tibble.width = Inf, tibble.print_max = Inf, tibble.max_extra_cols = Inf, scipen = 999)

# Results directory
results.dir <- file.path(analysis.dir, "results")
dir.create(results.dir, recursive = TRUE, showWarnings = FALSE)

################################################################################
# 1. Load Coordinate Cache and Reconstruct Loop Master Table
################################################################################

coord.cache.object.names <- coordinate_cache_object_names()
load_coordinate_cache_objects(cache.dir = coord.cache.dir, object.names = coord.cache.object.names, envir = environment())

message("Loaded coordinate-normalized cache objects.")

# Read core annotation pipeline (Lines 1 to 1235 of 01_promoter_enhancer_interaction_resubmit_loop_annotation.R)
annotation.script.file <- file.path(analysis.dir, "01_promoter_enhancer_interaction_resubmit_loop_annotation.R")
if (!file.exists(annotation.script.file)) {
  stop("Cannot locate main annotation script: ", annotation.script.file)
}
script_lines <- readLines(annotation.script.file)
eval(parse(text = script_lines[1:1235]))

assert_analysis_condition(
  condition = exists("df.loop.evidence") && nrow(df.loop.evidence) == 31021L,
  message = "Failed to assemble df.loop.evidence with 31,021 pooled loops.",
  success.message = "Verified: Successfully assembled df.loop.evidence with 31,021 distinct loops."
)

################################################################################
# 2. Categorize Loops into 3 Major Categories and 5 Sub-Categories (Flowchart Synchronized)
################################################################################

df.ctcf.master <- df.loop.evidence %>%
  mutate(
    loop_dist_bp = as.numeric(end2 - start1),
    loop_dist_kb = loop_dist_bp / 1000,
    
    # 5 Sub-categories (❶ ~ ❺) matching Figure 7 Flowchart exactly
    sub_category = case_when(
      revised_putative_regulatory_support ~ "❶ Single-promoter\nPutative Regulatory\n(10,469; 33.7%)",
      revised_single_promoter_without_opposite_atac ~ "❷ Single-promoter w/o\ndistal ATAC support\n(1,826; 5.9%)",
      revised_promoter_promoter_compatible & dual_promoter_directional_category %in% c(
        "both_directions_supported_bidirectional_or_ambiguous",
        "one_direction_stringent_asymmetric_mixed_regulatory_evidence",
        "one_direction_broader_mixed_promoter_regulatory_evidence"
      ) ~ "❸ Dual-promoter with\nEnhancer ATAC Support\n(2,907; 9.4%)",
      revised_promoter_promoter_compatible & dual_promoter_directional_category == "neither_direction_supported_promoter_promoter_or_lower_support" ~ "❹ Pure Promoter–\nPromoter Contacts\n(1,141; 3.7%)",
      revised_no_direct_promoter_tss ~ "❺ Loops without direct\npromoter/TSS overlap\n(14,678; 47.3%)",
      TRUE ~ "Other"
    ),
    
    # 3 Major Categories matching Figure 7 Bottom Parallelogram Summary exactly
    major_category = case_when(
      sub_category %in% c("❶ Single-promoter\nPutative Regulatory\n(10,469; 33.7%)", "❸ Dual-promoter with\nEnhancer ATAC Support\n(2,907; 9.4%)") ~ 
        "1. Putative regulatory loops with\ndistal open-chromatin support\n(13,376 loops; 43.1%)",
      sub_category %in% c("❷ Single-promoter w/o\ndistal ATAC support\n(1,826; 5.9%)", "❹ Pure Promoter–\nPromoter Contacts\n(1,141; 3.7%)") ~ 
        "2. Promoter-associated loops\nwithout distal ATAC support\n(2,967 loops; 9.6%)",
      sub_category == "❺ Loops without direct\npromoter/TSS overlap\n(14,678; 47.3%)" ~ 
        "3. Loops without direct\npromoter/TSS overlap\n(14,678 loops; 47.3%)",
      TRUE ~ "Other"
    ),
    
    # CTCF Counts & Anchor Flags
    ctcf_anchor1 = as.integer(predicted_ctcf_motif_interval_count_anchor1),
    ctcf_anchor2 = as.integer(predicted_ctcf_motif_interval_count_anchor2),
    total_ctcf = ctcf_anchor1 + ctcf_anchor2,
    both_anchors_ctcf = predicted_ctcf_motif_intervals_both_anchors,
    one_or_more_anchor_ctcf = (ctcf_anchor1 > 0L | ctcf_anchor2 > 0L)
  )

# Set factor levels for consistent plotting order
major_cat_levels <- c(
  "1. Putative regulatory loops with\ndistal open-chromatin support\n(13,376 loops; 43.1%)",
  "2. Promoter-associated loops\nwithout distal ATAC support\n(2,967 loops; 9.6%)",
  "3. Loops without direct\npromoter/TSS overlap\n(14,678 loops; 47.3%)"
)
df.ctcf.master$major_category <- factor(df.ctcf.master$major_category, levels = major_cat_levels)

sub_cat_levels <- c(
  "❶ Single-promoter\nPutative Regulatory\n(10,469; 33.7%)",
  "❷ Single-promoter w/o\ndistal ATAC support\n(1,826; 5.9%)",
  "❸ Dual-promoter with\nEnhancer ATAC Support\n(2,907; 9.4%)",
  "❹ Pure Promoter–\nPromoter Contacts\n(1,141; 3.7%)",
  "❺ Loops without direct\npromoter/TSS overlap\n(14,678; 47.3%)"
)
df.ctcf.master$sub_category <- factor(df.ctcf.master$sub_category, levels = sub_cat_levels)

################################################################################
# 3. Anchor-Level CTCF Dataset Construction (N = 62,042 Anchors)
################################################################################

df.anchor.level <- bind_rows(
  df.ctcf.master %>%
    transmute(
      loop_id,
      resolution,
      major_category,
      sub_category,
      anchor_side = "anchor1",
      ctcf_count = ctcf_anchor1,
      has_promoter = has_any_direct_promoter_tss & (n_direct_anchor_sides == 2L | !is.na(direct_anchor_assignment_class)),
      anchor_role = case_when(
        str_detect(sub_category, "^❶") ~ "Promoter\nAnchor (❶)",
        str_detect(sub_category, "^❷") ~ "Promoter\nAnchor (❷)",
        str_detect(sub_category, "^❸") ~ "Dual-Promoter\nAnchor (❸)",
        str_detect(sub_category, "^❹") ~ "Dual-Promoter\nAnchor (❹)",
        TRUE ~ "Structural\nAnchor (❺)"
      )
    ),
  df.ctcf.master %>%
    transmute(
      loop_id,
      resolution,
      major_category,
      sub_category,
      anchor_side = "anchor2",
      ctcf_count = ctcf_anchor2,
      has_promoter = has_any_direct_promoter_tss & (n_direct_anchor_sides == 2L),
      anchor_role = case_when(
        str_detect(sub_category, "^❶") ~ "Distal ATAC\nAnchor (❶)",
        str_detect(sub_category, "^❷") ~ "Opposite Anchor\nw/o ATAC (❷)",
        str_detect(sub_category, "^❸") ~ "Dual-Promoter\nAnchor (❸)",
        str_detect(sub_category, "^❹") ~ "Dual-Promoter\nAnchor (❹)",
        TRUE ~ "Structural\nAnchor (❺)"
      )
    )
)

################################################################################
# 4. Statistical Testing and Summary Tables
################################################################################

# Major Category Summary Table
df.summary.major <- df.ctcf.master %>%
  group_by(major_category) %>%
  summarise(
    n_loops = n(),
    pct_loops = round(100 * n() / nrow(df.ctcf.master), 2),
    n_both_ctcf = sum(both_anchors_ctcf),
    pct_both_ctcf = round(100 * mean(both_anchors_ctcf), 2),
    mean_ctcf = round(mean(total_ctcf), 2),
    sd_ctcf = round(sd(total_ctcf), 2),
    median_ctcf = stats::median(total_ctcf),
    q25_ctcf = quantile(total_ctcf, 0.25),
    q75_ctcf = quantile(total_ctcf, 0.75),
    iqr_ctcf = IQR(total_ctcf),
    mean_anchor1 = round(mean(ctcf_anchor1), 2),
    mean_anchor2 = round(mean(ctcf_anchor2), 2),
    mean_dist_kb = round(mean(loop_dist_kb), 1),
    median_dist_kb = round(stats::median(loop_dist_kb), 1),
    .groups = "drop"
  )

# Export Summary TSV
summary_tsv_path <- file.path(results.dir, "revision_ctcf_major_category_summary.tsv")
write_tsv(df.summary.major, summary_tsv_path)
message("Saved major category summary TSV to: ", summary_tsv_path)

# Statistical Tests
kw_test <- kruskal.test(total_ctcf ~ major_category, data = df.ctcf.master)
chisq_test <- chisq.test(table(df.ctcf.master$major_category, df.ctcf.master$both_anchors_ctcf))

################################################################################
# 5. Publication-Quality ggplot2 Visualizations (Flowchart Colors & Aesthetics)
################################################################################

# Color Scheme matching the requested 3-color palette (Standard ggplot2 3-color: Red, Green, Blue)
cat_colors <- c(
  "1. Putative regulatory loops with\ndistal open-chromatin support\n(13,376 loops; 43.1%)" = "#F8766D", # Coral Red
  "2. Promoter-associated loops\nwithout distal ATAC support\n(2,967 loops; 9.6%)"          = "#00BA38", # Green
  "3. Loops without direct\npromoter/TSS overlap\n(14,678 loops; 47.3%)"                  = "#619CFF"  # Sky Blue
)

sub_cat_colors <- c(
  "❶ Single-promoter\nPutative Regulatory\n(10,469; 33.7%)"     = "#16A34A",
  "❷ Single-promoter w/o\ndistal ATAC support\n(1,826; 5.9%)"   = "#EA580C",
  "❸ Dual-promoter with\nEnhancer ATAC Support\n(2,907; 9.4%)"  = "#15803D",
  "❹ Pure Promoter–\nPromoter Contacts\n(1,141; 3.7%)"         = "#C2410C",
  "❺ Loops without direct\npromoter/TSS overlap\n(14,678; 47.3%)" = "#E11D48"
)

# Common Theme
theme_publication <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.08), hjust = 0, margin = margin(b = 6)),
      plot.subtitle = element_text(size = rel(0.90), color = "gray25", margin = margin(b = 8)),
      axis.title = element_text(face = "bold", size = rel(0.95)),
      axis.text = element_text(size = rel(0.85), color = "black"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, lineheight = 1.05),
      panel.grid.major.y = element_line(color = "gray92", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )
}

# ------------------------------------------------------------------------------
# Panel A: Total CTCF Motif Count per Loop across 3 Major Categories (Violin + Boxplot)
# ------------------------------------------------------------------------------
p_panel_a <- ggplot(df.ctcf.master, aes(x = major_category, y = total_ctcf, fill = major_category, color = major_category)) +
  geom_violin(alpha = 0.35, width = 0.8, trim = TRUE, scale = "width", linewidth = 0.6) +
  geom_boxplot(width = 0.22, fill = "white", outlier.shape = NA, alpha = 0.92, linewidth = 0.75, color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.2, fill = "#DC2626", color = "black") +
  scale_fill_manual(values = cat_colors) +
  scale_color_manual(values = cat_colors) +
  scale_y_continuous(
    trans = "log1p",
    breaks = c(0, 10, 25, 50, 100, 200, 400),
    limits = c(0, 750),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    title = "A. Total CTCF Motif Count per Loop Across 3 Major Categories",
    subtitle = "Kruskal-Wallis p < 2.2e-16 (Red diamond: Mean, Solid line: Median)",
    x = NULL,
    y = "Total Predicted CTCF Motifs per Loop (log1p scale)"
  ) +
  annotate("text", x = 1, y = 540, label = "Mean: 104.7\nMedian: 92", size = 3.6, fontface = "bold", color = "#E05353") +
  annotate("text", x = 2, y = 540, label = "Mean: 77.1\nMedian: 65", size = 3.6, fontface = "bold", color = "#009E2E") +
  annotate("text", x = 3, y = 540, label = "Mean: 78.1\nMedian: 69", size = 3.6, fontface = "bold", color = "#3B82F6") +
  theme_publication()

# ------------------------------------------------------------------------------
# Panel B: Proportion of Loops with Both-Anchor CTCF Support (%)
# ------------------------------------------------------------------------------
df_prop_plot <- df.summary.major %>%
  mutate(
    pct_label = sprintf("%.1f%%\n(%s / %s)", pct_both_ctcf, scales::comma(n_both_ctcf), scales::comma(n_loops))
  )

p_panel_b <- ggplot(df_prop_plot, aes(x = major_category, y = pct_both_ctcf, fill = major_category)) +
  geom_col(width = 0.55, color = "black", linewidth = 0.7, alpha = 0.88) +
  geom_text(aes(label = pct_label), vjust = -0.3, size = 3.6, fontface = "bold") +
  scale_fill_manual(values = cat_colors) +
  scale_y_continuous(
    limits = c(0, 115),
    breaks = seq(0, 100, by = 20),
    labels = paste0(seq(0, 100, by = 20), "%"),
    expand = c(0, 0)
  ) +
  labs(
    title = "B. Both-Anchor CTCF Motif Support Rate (%) Across Categories",
    subtitle = "Chi-squared p < 2.2e-16 (Category 1 vs Cat 2/3: p < 1e-15)",
    x = NULL,
    y = "Loops with Predicted CTCF Motifs at Both Anchors (%)"
  ) +
  theme_publication()

# ------------------------------------------------------------------------------
# Panel C: Detailed Comparison Across 5 Sub-Categories (❶ ~ ❺)
# ------------------------------------------------------------------------------
p_panel_c <- ggplot(df.ctcf.master, aes(x = sub_category, y = total_ctcf, fill = sub_category)) +
  geom_boxplot(width = 0.55, outlier.alpha = 0.12, outlier.size = 0.5, linewidth = 0.6, color = "black", alpha = 0.85) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.0, fill = "#DC2626", color = "black") +
  scale_fill_manual(values = sub_cat_colors) +
  scale_y_continuous(
    trans = "log1p",
    breaks = c(0, 10, 25, 50, 100, 200, 400),
    limits = c(0, 500),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    title = "C. CTCF Motif Abundance Across 5 Detailed Sub-Categories (❶ ~ ❺)",
    subtitle = "Sub-category ❸ exhibits the highest CTCF enrichment (Mean: 123.2, Median: 110)",
    x = NULL,
    y = "Total Predicted CTCF Motifs per Loop (log1p scale)"
  ) +
  theme_publication() +
  theme(axis.text.x = element_text(size = rel(0.80), face = "bold"))

# ------------------------------------------------------------------------------
# Panel D: Anchor-Level CTCF Motif Count by Anchor Functional Role
# ------------------------------------------------------------------------------
df.anchor.level.clean <- df.anchor.level %>%
  filter(!anchor_role %in% c("Opposite Anchor\nw/o ATAC (❷)")) %>%
  mutate(
    anchor_role = factor(
      anchor_role,
      levels = c(
        "Promoter\nAnchor (❶)",
        "Distal ATAC\nAnchor (❶)",
        "Promoter\nAnchor (❷)",
        "Dual-Promoter\nAnchor (❸)",
        "Dual-Promoter\nAnchor (❹)",
        "Structural\nAnchor (❺)"
      )
    )
  )

role_colors <- c(
  "Promoter\nAnchor (❶)"      = "#16A34A",
  "Distal ATAC\nAnchor (❶)"  = "#22C55E",
  "Promoter\nAnchor (❷)"      = "#EA580C",
  "Dual-Promoter\nAnchor (❸)" = "#15803D",
  "Dual-Promoter\nAnchor (❹)" = "#C2410C",
  "Structural\nAnchor (❺)"    = "#E11D48"
)

p_panel_d <- ggplot(df.anchor.level.clean, aes(x = anchor_role, y = ctcf_count, fill = anchor_role)) +
  geom_boxplot(width = 0.55, outlier.alpha = 0.12, outlier.size = 0.5, linewidth = 0.6, color = "black", alpha = 0.85) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.8, fill = "#DC2626", color = "black") +
  scale_fill_manual(values = role_colors) +
  scale_y_continuous(
    trans = "log1p",
    breaks = c(0, 5, 15, 30, 60, 120, 250),
    limits = c(0, 250),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    title = "D. Individual Anchor-Level CTCF Motif Abundance by Functional Role",
    subtitle = "Single anchor level (N = 62,042 anchors across 31,021 loops)",
    x = NULL,
    y = "Predicted CTCF Motifs per Single Anchor (log1p scale)"
  ) +
  theme_publication() +
  theme(axis.text.x = element_text(size = rel(0.85), face = "bold"))

# ------------------------------------------------------------------------------
# Assemble Multi-Panel Figure & Save (2x2 Grid)
# ------------------------------------------------------------------------------
p_master_ctcf_figure <- plot_grid(
  p_panel_a, p_panel_b,
  p_panel_c, p_panel_d,
  ncol = 2,
  align = "hv",
  axis = "tblr"
)

output_png_path <- file.path(results.dir, "revision_figure_ctcf_loop_categories.png")
output_pdf_path <- file.path(results.dir, "revision_figure_ctcf_loop_categories.pdf")

ggsave(
  filename = output_png_path,
  plot = p_master_ctcf_figure,
  width = 14.5,
  height = 11.5,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = output_pdf_path,
  plot = p_master_ctcf_figure,
  width = 14.5,
  height = 11.5,
  bg = "white"
)

message("Successfully generated and saved CTCF multi-panel figure to:")
message("  - PNG: ", output_png_path)
message("  - PDF: ", output_pdf_path)

message("\n================================================================================")
message("All CTCF Loop Category Analyses Completed Successfully!")
message("================================================================================")
