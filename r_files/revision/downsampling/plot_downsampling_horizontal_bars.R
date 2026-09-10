# Plot the two horizontal bar charts used in the sequencing-depth slide.
#
# Usage (from this directory after copying the script into downsampling/):
#   Rscript plot_downsampling_horizontal_bars.R
#
# Optional:
#   Rscript plot_downsampling_horizontal_bars.R /path/to/downsampling
#   DOWNSAMPLING_BAR_OUTPUT_DIR=/path/to/output Rscript ...

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

current_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) == 0L) {
    return(getwd())
  }
  dirname(normalizePath(
    sub("^--file=", "", file_arg[[1]]),
    winslash = "/",
    mustWork = FALSE
  ))
}

args <- commandArgs(trailingOnly = TRUE)
analysis_dir <- if (length(args) >= 1L) {
  normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
} else {
  normalizePath(
    Sys.getenv(
      "DOWNSAMPLING_ANALYSIS_DIR",
      unset = current_script_dir()
    ),
    winslash = "/",
    mustWork = TRUE
  )
}

matched_result_dir <- file.path(
  analysis_dir, "results", "full_140M_250M"
)
all_ten_result_file <- file.path(
  analysis_dir,
  "results",
  "full_140M",
  "downsampling_pooled_exact_recovery.tsv"
)
summary_file <- file.path(
  analysis_dir, "downsampling_depth_sensitivity_summary.md"
)
output_dir <- normalizePath(
  Sys.getenv(
    "DOWNSAMPLING_BAR_OUTPUT_DIR",
    unset = file.path(analysis_dir, "figures")
  ),
  winslash = "/",
  mustWork = FALSE
)
resubmit_2nd_dir <- normalizePath(
  file.path(dirname(analysis_dir), "revision_main", "results", "2nd_resubmission"),
  winslash = "/",
  mustWork = FALSE
)
target_dirs <- unique(c(output_dir, resubmit_2nd_dir))
for (d in target_dirs) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

required_matched_files <- c(
  exact = file.path(matched_result_dir, "pooled_exact_recovery.tsv"),
  approximate = file.path(
    matched_result_dir, "pooled_approximate_locus_recovery.tsv"
  )
)
missing_matched <- required_matched_files[!file.exists(required_matched_files)]
if (length(missing_matched) > 0L) {
  stop(
    "Missing matched-depth result file(s): ",
    paste(missing_matched, collapse = ", "),
    call. = FALSE
  )
}

# Prefer the rebuildable all-ten output table. The current checkout does not
# contain results/full_140M, so the committed audit summary is a transparent
# fallback until downsampling_140M.R is rerun locally.
read_all_ten_resolution_recovery <- function() {
  if (file.exists(all_ten_result_file)) {
    message("All-ten source: ", all_ten_result_file)
    return(
      read_tsv(all_ten_result_file, show_col_types = FALSE) %>%
        filter(resolution %in% c("5K", "10K", "25K")) %>%
        transmute(
          resolution = str_replace(resolution, "K$", " kb"),
          recovered_n = n_exact_shared_pooled,
          reference_n = n_full_depth_pooled,
          recovery_pct = pct_full_depth_pooled_exact_recovered,
          source = "results/full_140M/downsampling_pooled_exact_recovery.tsv"
        )
    )
  }

  if (!file.exists(summary_file)) {
    stop(
      "Neither the all-ten result table nor the analysis summary exists.",
      call. = FALSE
    )
  }

  message(
    "All-ten result table is absent; reading the audited resolution table in ",
    summary_file
  )
  summary_lines <- read_lines(summary_file)
  table_rows <- summary_lines[
    str_detect(summary_lines, "^\\|\\s*(5|10|25)\\s*kb\\s*\\|")
  ]
  if (length(table_rows) != 3L) {
    stop(
      "Could not identify all three all-ten resolution rows in the summary.",
      call. = FALSE
    )
  }

  map_dfr(table_rows, function(line) {
    fields <- str_split(line, fixed("|"), simplify = TRUE)[1, ]
    fields <- str_trim(fields[nzchar(str_trim(fields))])
    tibble(
      resolution = fields[[1]],
      recovered_n = parse_number(fields[[4]]),
      reference_n = parse_number(fields[[2]]),
      recovery_pct = parse_number(fields[[5]]),
      source = "downsampling_depth_sensitivity_summary.md"
    )
  })
}

all_ten_recovery <- read_all_ten_resolution_recovery() %>%
  mutate(
    resolution = factor(
      resolution,
      levels = c("25 kb", "10 kb", "5 kb")
    ),
    label = sprintf(
      "%.1f%% (%s/%s)",
      recovery_pct,
      format(recovered_n, big.mark = ",", scientific = FALSE, trim = TRUE),
      format(reference_n, big.mark = ",", scientific = FALSE, trim = TRUE)
    )
  )

if (nrow(all_ten_recovery) != 3L || anyNA(all_ten_recovery$recovery_pct)) {
  stop("All-ten recovery data are incomplete.", call. = FALSE)
}

exact_recovery <- read_tsv(
  required_matched_files[["exact"]], show_col_types = FALSE
) %>%
  filter(
    reference_condition == "full_depth",
    query_condition %in% c("downsample_140M", "downsample_250M")
  ) %>%
  transmute(
    metric = "Exact calls",
    depth = recode(
      query_condition,
      downsample_140M = "140M",
      downsample_250M = "250M"
    ),
    recovered_n = n_shared_pooled_exact_calls,
    reference_n = n_reference_pooled_exact_calls,
    recovery_pct = pct_reference_pooled_exact_recovered,
    source = basename(required_matched_files[["exact"]])
  )

approximate_recovery <- read_tsv(
  required_matched_files[["approximate"]], show_col_types = FALSE
) %>%
  filter(
    reference_condition == "full_depth",
    query_condition %in% c("downsample_140M", "downsample_250M")
  ) %>%
  transmute(
    metric = "Approximate loci",
    depth = recode(
      query_condition,
      downsample_140M = "140M",
      downsample_250M = "250M"
    ),
    recovered_n = n_shared_approximate_loci,
    reference_n = n_reference_approximate_loci,
    recovery_pct = pct_reference_approximate_recovered,
    source = basename(required_matched_files[["approximate"]])
  )

matched_recovery <- bind_rows(exact_recovery, approximate_recovery) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("Approximate loci", "Exact calls")
    ),
    depth = factor(depth, levels = c("140M", "250M")),
    label = sprintf(
      "%.1f%% (%s/%s)",
      recovery_pct,
      format(recovered_n, big.mark = ",", scientific = FALSE, trim = TRUE),
      format(reference_n, big.mark = ",", scientific = FALSE, trim = TRUE)
    )
  )

if (nrow(matched_recovery) != 4L || anyNA(matched_recovery$recovery_pct)) {
  stop("Matched seven-library recovery data are incomplete.", call. = FALSE)
}

navy <- "#123E6A"
orange <- "#E68A43"
blue <- "#2F7FC1"
teal <- "#188F8B"
gray <- "#607080"
grid <- "#DCE4EA"

base_theme <- theme_minimal(base_size = 13, base_family = "sans") +
  theme(
    plot.title = element_text(
      color = navy, face = "bold", size = 15, hjust = 0
    ),
    plot.subtitle = element_text(color = gray, size = 10.5),
    axis.title = element_blank(),
    axis.text.x = element_text(color = gray, size = 10),
    axis.text.y = element_text(
      color = "#26384A", face = "bold", size = 11
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = grid, linewidth = 0.5),
    plot.margin = margin(10, 28, 10, 6),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

p_resolution <- ggplot(
  all_ten_recovery,
  aes(x = recovery_pct, y = resolution, fill = resolution)
) +
  geom_col(width = 0.48, show.legend = FALSE) +
  geom_text(
    aes(label = label, color = resolution),
    hjust = -0.20,
    size = 3.3,
    fontface = "bold",
    lineheight = 0.95,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c(
    "5 kb" = orange,
    "10 kb" = blue,
    "25 kb" = teal
  )) +
  scale_color_manual(values = c(
    "5 kb" = "#B85C1B",
    "10 kb" = "#1C5F9D",
    "25 kb" = "#08736F"
  )) +
  scale_x_continuous(
    position = "top",
    limits = c(0, 85),
    breaks = c(0, 20, 40, 60, 80),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Exact recovery at 140M is resolution dependent",
    subtitle = paste0(
      "All 10 libraries; labels show recovered / full-depth reference calls"
    )
  ) +
  base_theme

p_matched <- ggplot(
  matched_recovery,
  aes(x = recovery_pct, y = metric, fill = depth)
) +
  geom_col(
    width = 0.56,
    position = position_dodge2(
      width = 0.70, preserve = "single", reverse = TRUE
    ),
    orientation = "y"
  ) +
  geom_text(
    aes(label = label, color = depth),
    position = position_dodge2(
      width = 0.70, preserve = "single", reverse = TRUE
    ),
    hjust = -0.18,
    size = 3.3,
    fontface = "bold",
    lineheight = 0.95,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("140M" = orange, "250M" = blue)) +
  scale_color_manual(values = c("140M" = "#B85C1B", "250M" = "#1C5F9D")) +
  scale_x_continuous(
    position = "top",
    limits = c(0, 110),
    breaks = c(0, 20, 40, 60, 80, 100),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Full-depth features recovered",
    subtitle = paste0(
      "Matched 7 libraries; labels show recovered / full-depth reference"
    ),
    fill = NULL
  ) +
  base_theme +
  theme(
    legend.position = "bottom",
    legend.justification = "left",
    legend.text = element_text(color = gray, size = 10),
    legend.key.size = unit(0.55, "lines")
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

combined_plot <- p_resolution + p_matched +
  plot_layout(widths = c(0.92, 1.08))

for (out_d in target_dirs) {
  ggsave(
    file.path(out_d, "exact_recovery_140M_by_resolution.png"),
    p_resolution,
    width = 6.2,
    height = 4.1,
    dpi = 300,
    bg = "white"
  )
  ggsave(
    file.path(out_d, "exact_recovery_140M_by_resolution.pdf"),
    p_resolution,
    width = 6.2,
    height = 4.1,
    bg = "white"
  )
  ggsave(
    file.path(out_d, "matched_7_full_depth_recovery.png"),
    p_matched,
    width = 7.0,
    height = 4.1,
    dpi = 300,
    bg = "white"
  )
  ggsave(
    file.path(out_d, "matched_7_full_depth_recovery.pdf"),
    p_matched,
    width = 7.0,
    height = 4.1,
    bg = "white"
  )
  ggsave(
    file.path(out_d, "downsampling_horizontal_bars_combined.png"),
    combined_plot,
    width = 13.2,
    height = 4.2,
    dpi = 300,
    bg = "white"
  )
  ggsave(
    file.path(out_d, "downsampling_horizontal_bars_combined.pdf"),
    combined_plot,
    width = 13.2,
    height = 4.2,
    bg = "white"
  )
}

plot_data <- bind_rows(
  all_ten_recovery %>%
    transmute(
      panel = "all_10_resolution",
      category = as.character(resolution),
      depth = "140M",
      recovered_n,
      reference_n,
      recovery_pct,
      source
    ),
  matched_recovery %>%
    transmute(
      panel = "matched_7_depth",
      category = as.character(metric),
      depth = as.character(depth),
      recovered_n,
      reference_n,
      recovery_pct,
      source
    )
)
for (out_d in target_dirs) {
  write_tsv(
    plot_data,
    file.path(out_d, "downsampling_horizontal_bars_plot_data.tsv")
  )
}

message("Wrote horizontal-bar figures to: ", paste(target_dirs, collapse = ", "))
print(plot_data)
