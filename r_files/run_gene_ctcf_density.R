#!/usr/bin/env Rscript
# Standalone script: Gene density vs CTCF density correlation analysis
# Extracted from enhancer_promoter_interaction.R (L1115-L1277)

library("tidyverse")
library("ggplot2")

options(scipen = 999)

# Set working directory to match the main script
setwd("~/dropbox/Gateway_to_Hao/enhancer/r_files")
source("~/Desktop/playground/enhancer/r_files/utils_functions.R")

cat("=== Working directory:", getwd(), "===\n\n")

########################
# 1. Load chromosome data
########################
cache_file_chromosome_data <- "~/dropbox/Gateway_to_Hao/enhancer/r_files/rds/df.chromosome.data.rds"

if (file.exists(cache_file_chromosome_data)) {
  message("Loading cached chromosome data from: ", cache_file_chromosome_data)
  df.chromosome.data <- readRDS(cache_file_chromosome_data)
} else {
  message("Processing and caching chromosome data...")
  df.chromosome.data <- read.table(file = "../data/rn7_chromosome_length_from_ucsc.tsv", sep = "\t") %>%
    mutate(start = 0) %>%
    dplyr::rename(chr = V1, end = V2)
  saveRDS(df.chromosome.data, cache_file_chromosome_data)
}

df.chromosome.data <- df.chromosome.data %>% mutate(CE_start = NA, CE_end = NA)
cat("Chromosome data loaded:", nrow(df.chromosome.data), "chromosomes\n")

########################
# 2. Load CTCF data
########################
cache_file_distinct_fimo_2nd_ctcf <- "~/dropbox/Gateway_to_Hao/enhancer/r_files/rds/df.DISTINCT.fimo.2nd.trial.ctcf.rds"

if (file.exists(cache_file_distinct_fimo_2nd_ctcf)) {
  message("Loading cached DISTINCT fimo 2nd trial ctcf data from: ", cache_file_distinct_fimo_2nd_ctcf)
  df.DISTINCT.fimo.2nd.trial.ctcf <- readRDS(cache_file_distinct_fimo_2nd_ctcf)
} else {
  message("Processing and caching DISTINCT fimo 2nd trial ctcf data...")
  df.init.ctcf <- read.table(file = "~/dropbox/Gateway_to_Hao/enhancer/data/ctcf/submission/E4/fimo_E4_submission_trial.txt", header = TRUE, sep = "\t") %>%
    dplyr::rename(chr = sequence_name, end = stop) %>%
    mutate(length = end - start)

  df.DISTINCT.fimo.2nd.trial.ctcf <- df.init.ctcf %>%
    distinct(chr, start, end) %>%
    mutate(start = as.numeric(start)) %>%
    mutate(end = as.numeric(end)) %>%
    mutate(ctcf_pos = as.numeric(round((start + end) / 2))) %>%
    mutate(id = str_c(chr, "_", start, "_", end, "_", ctcf_pos))

  saveRDS(df.DISTINCT.fimo.2nd.trial.ctcf, cache_file_distinct_fimo_2nd_ctcf)
}

cat("CTCF sites loaded:", nrow(df.DISTINCT.fimo.2nd.trial.ctcf), "\n")

df_ctcf_ideogram <- df.DISTINCT.fimo.2nd.trial.ctcf %>%
  dplyr::select(chr, start, end)

########################
# 3. Load gene (TSS) data
########################
df.ensembl.gtf.for.tss.DISTINCT.geneid <- readRDS("~/dropbox/Gateway_to_Hao/enhancer/r_files/rds/df_ensembl_gtf_for_tss_DISTINCT_geneid.rds")
cat("Gene TSS data loaded:", nrow(df.ensembl.gtf.for.tss.DISTINCT.geneid), "\n")

df_gene_ideogram <- df.ensembl.gtf.for.tss.DISTINCT.geneid %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  dplyr::select(chr, start, end, gene_id, gene_name)

cat("Unique genes for density:", nrow(df_gene_ideogram), "\n\n")

########################
# 4. Compute chromosome-level density summary
########################
chromosome_density_summary <- df.chromosome.data %>%
  mutate(
    chr = str_remove(as.character(chr), "chr"),
    chromosome_length_mb = end / 1e6
  ) %>%
  dplyr::select(chr, chromosome_length_mb) %>%
  left_join(
    df_gene_ideogram %>%
      mutate(chr = str_remove(as.character(chr), "chr")) %>%
      count(chr, name = "gene_count"),
    by = "chr"
  ) %>%
  left_join(
    df_ctcf_ideogram %>%
      mutate(chr = str_remove(as.character(chr), "chr")) %>%
      count(chr, name = "ctcf_count"),
    by = "chr"
  ) %>%
  mutate(
    gene_count = replace_na(gene_count, 0L),
    ctcf_count = replace_na(ctcf_count, 0L),
    genes_per_mb = gene_count / chromosome_length_mb,
    ctcf_per_mb = ctcf_count / chromosome_length_mb
  ) %>%
  arrange(desc(ctcf_per_mb))

cat("=== Chromosome Density Summary ===\n")
print(as.data.frame(chromosome_density_summary))

########################
# 5. Correlation tests
########################
gene_ctcf_density_cor_pearson <- cor.test(
  chromosome_density_summary$genes_per_mb,
  chromosome_density_summary$ctcf_per_mb,
  method = "pearson"
)

gene_ctcf_density_cor_spearman <- cor.test(
  chromosome_density_summary$genes_per_mb,
  chromosome_density_summary$ctcf_per_mb,
  method = "spearman"
)

gene_ctcf_density_stats <- tibble(
  method = c("pearson", "spearman"),
  estimate = c(
    unname(gene_ctcf_density_cor_pearson$estimate),
    unname(gene_ctcf_density_cor_spearman$estimate)
  ),
  p_value = c(
    gene_ctcf_density_cor_pearson$p.value,
    gene_ctcf_density_cor_spearman$p.value
  )
)

cat("\n=== Correlation Results ===\n")
print(as.data.frame(gene_ctcf_density_stats))

cat("\n=== Pearson Details ===\n")
print(gene_ctcf_density_cor_pearson)

cat("\n=== Spearman Details ===\n")
print(gene_ctcf_density_cor_spearman)

########################
# 6. Save outputs
########################
output_dir <- "./figures/submission/lt2mb"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  chromosome_density_summary,
  file = file.path(output_dir, "chromosome_gene_ctcf_density_summary.csv"),
  row.names = FALSE
)

write.csv(
  gene_ctcf_density_stats,
  file = file.path(output_dir, "chromosome_gene_ctcf_density_correlation.csv"),
  row.names = FALSE
)

########################
# 7. Scatter plot
########################
plot_gene_ctcf_density <- chromosome_density_summary %>%
  ggplot(aes(x = genes_per_mb, y = ctcf_per_mb, label = chr)) +
  geom_point(size = 2.5, color = "#1B4F72") +
  geom_smooth(method = "lm", se = FALSE, color = "#C0392B", linewidth = 0.8) +
  geom_text(nudge_y = 0.6, size = 3, check_overlap = TRUE) +
  labs(
    title = "Chromosome-level association between gene density and CTCF density",
    x = "Genes per Mb",
    y = "CTCF sites per Mb"
  ) +
  theme_bw(base_size = 10)

ggsave(
  filename = file.path(output_dir, "chromosome_gene_ctcf_density_scatter.pdf"),
  plot = plot_gene_ctcf_density,
  width = 6,
  height = 5
)

ggsave(
  filename = file.path(output_dir, "chromosome_gene_ctcf_density_scatter.png"),
  plot = plot_gene_ctcf_density,
  width = 6,
  height = 5,
  dpi = 300
)

cat("\n=== Output saved to:", output_dir, "===\n")
cat("  - chromosome_gene_ctcf_density_summary.csv\n")
cat("  - chromosome_gene_ctcf_density_correlation.csv\n")
cat("  - chromosome_gene_ctcf_density_scatter.pdf\n")
cat("  - chromosome_gene_ctcf_density_scatter.png\n")
cat("\nDone!\n")
