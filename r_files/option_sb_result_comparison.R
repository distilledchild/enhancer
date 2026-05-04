library(dplyr)
library(tidyr)
library(stringr)
library(readr)

################################################################
# PART 1: Compare original vs 592BB QC stats
# Compare original vs 592BB QC stats for manuscript revision
################################################################

original <- read.delim("~/dropbox/Gateway_to_Hao/enhancer/data/library_complexity.tsv", check.names = FALSE)
updated <- read.delim("~/dropbox/Gateway_to_Hao/enhancer/data/library_complexity_592BB.tsv", check.names = FALSE)

# Loop counts per sample (from bedpe files)
loop.counts <- tibble::tribble(
  ~Strain, ~num_loop,
  "SHR/OlaIpcv", 5263,
  "BN-Lx", 6535,
  "BXH6", 6568,
  "HXB2", 4656,
  "HXB10", 7336,
  "HXB23", 7676,
  "HXB31", 9131,
  "LE/Stm", 2903,
  "F344/Stm", 2992,
  "SHR/OlaIpcvxBN/NHsdMcwi", 5932
)

cat("=== 1. Average reads per sample ===\n")
cat("ORIGINAL mean Sequenced_RP:", format(mean(original$Sequenced_RP), big.mark = ","), "\n")
cat("UPDATED  mean Sequenced_RP:", format(mean(updated$Sequenced_RP), big.mark = ","), "\n")
cat("ORIGINAL mean (million):", round(mean(original$Sequenced_RP) / 1e6, 1), "\n")
cat("UPDATED  mean (million):", round(mean(updated$Sequenced_RP) / 1e6, 1), "\n\n")

cat("=== 2. Unalignable reads % ===\n")
# Unalignable = Chimeric_Ambiguous + Unmapped (chimeric paired is alignable)
original$unalignable_pct <- (original$Chimeric_Ambiguous + original$Unmapped) / original$Sequenced_RP * 100
updated$unalignable_pct <- (updated$Chimeric_Ambiguous + updated$Unmapped) / updated$Sequenced_RP * 100
cat("ORIGINAL mean unalignable %:", round(mean(original$unalignable_pct), 1), "STD:", round(sd(original$unalignable_pct), 1), "\n")
cat("UPDATED  mean unalignable %:", round(mean(updated$unalignable_pct), 1), "STD:", round(sd(updated$unalignable_pct), 1), "\n")

# Unmapped only
original$unmapped_pct <- original$Unmapped / original$Sequenced_RP * 100
updated$unmapped_pct <- updated$Unmapped / updated$Sequenced_RP * 100
cat("ORIGINAL mean unmapped %:", round(mean(original$unmapped_pct), 1), "\n")
cat("UPDATED  mean unmapped %:", round(mean(updated$unmapped_pct), 1), "\n\n")

cat("=== 3. PCR + Optical duplicates % ===\n")
original$dup_pct <- (original$PCR_Duplicates + original$Optical_Duplicates) / original$Sequenced_RP * 100
updated$dup_pct <- (updated$PCR_Duplicates + updated$Optical_Duplicates) / updated$Sequenced_RP * 100
cat("ORIGINAL mean dup %:", round(mean(original$dup_pct), 1), "STD:", round(sd(original$dup_pct), 1), "\n")
cat("UPDATED  mean dup %:", round(mean(updated$dup_pct), 1), "STD:", round(sd(updated$dup_pct), 1), "\n")

# SHR/OlaIpcv specific
cat("ORIGINAL SHR dup %:", round(original$dup_pct[original$Strain == "SHR/OlaIpcv"], 1), "\n")
cat("UPDATED  SHR dup %:", round(updated$dup_pct[updated$Strain == "SHR/OlaIpcv"], 1), "\n\n")

cat("=== 4. Unique reads % ===\n")
original$unique_pct <- original$Unique_Reads / original$Sequenced_RP * 100
updated$unique_pct <- updated$Unique_Reads / updated$Sequenced_RP * 100
cat("ORIGINAL mean unique %:", round(mean(original$unique_pct), 0), "± STD:", round(sd(original$unique_pct), 1), "\n")
cat("UPDATED  mean unique %:", round(mean(updated$unique_pct), 0), "± STD:", round(sd(updated$unique_pct), 1), "\n\n")

cat("=== 5. Correlation: loops vs depth ===\n")
# Original
merged_orig <- left_join(loop.counts, original, by = "Strain")
cor_total_orig <- cor.test(merged_orig$num_loop, merged_orig$Sequenced_RP)
cor_unique_orig <- cor.test(merged_orig$num_loop, merged_orig$Unique_Reads)
cor_alignable_orig <- cor.test(merged_orig$num_loop, merged_orig$Alignable_Normal_N_Chimeric)

# Updated
merged_upd <- left_join(loop.counts, updated, by = "Strain")
cor_total_upd <- cor.test(merged_upd$num_loop, merged_upd$Sequenced_RP)
cor_unique_upd <- cor.test(merged_upd$num_loop, merged_upd$Unique_Reads)
cor_alignable_upd <- cor.test(merged_upd$num_loop, merged_upd$Alignable_Normal_N_Chimeric)

cat("ORIGINAL correlations:\n")
cat("  Total reads:     r =", round(cor_total_orig$estimate, 3), "p =", formatC(cor_total_orig$p.value, format = "e", digits = 2), "\n")
cat("  Unique reads:    r =", round(cor_unique_orig$estimate, 3), "p =", formatC(cor_unique_orig$p.value, format = "e", digits = 2), "\n")
cat("  Alignable reads: r =", round(cor_alignable_orig$estimate, 3), "p =", formatC(cor_alignable_orig$p.value, format = "e", digits = 2), "\n")

cat("UPDATED correlations:\n")
cat("  Total reads:     r =", round(cor_total_upd$estimate, 3), "p =", formatC(cor_total_upd$p.value, format = "e", digits = 2), "\n")
cat("  Unique reads:    r =", round(cor_unique_upd$estimate, 3), "p =", formatC(cor_unique_upd$p.value, format = "e", digits = 2), "\n")
cat("  Alignable reads: r =", round(cor_alignable_upd$estimate, 3), "p =", formatC(cor_alignable_upd$p.value, format = "e", digits = 2), "\n\n")

cat("=== 6. Loop count stats ===\n")
cat("Mean loops:", round(mean(loop.counts$num_loop), 0), "± STD:", round(sd(loop.counts$num_loop), 0), "\n")
cat("Range:", min(loop.counts$num_loop), "-", max(loop.counts$num_loop), "\n")
cat("Total:", sum(loop.counts$num_loop), "\n\n")

cat("=== 7. SHR/OlaIpcv specific changes ===\n")
shr_orig <- original %>% filter(Strain == "SHR/OlaIpcv")
shr_upd <- updated %>% filter(Strain == "SHR/OlaIpcv")
cat("ORIGINAL SHR Sequenced_RP:", format(shr_orig$Sequenced_RP, big.mark = ","), "\n")
cat("UPDATED  SHR Sequenced_RP:", format(shr_upd$Sequenced_RP, big.mark = ","), "\n")
cat("ORIGINAL SHR Unique_Reads:", format(shr_orig$Unique_Reads, big.mark = ","), "\n")
cat("UPDATED  SHR Unique_Reads:", format(shr_upd$Unique_Reads, big.mark = ","), "\n")
cat("ORIGINAL SHR PCR dup rate:", round(shr_orig$PCR_Duplicates / shr_orig$Sequenced_RP * 100, 1), "%\n")
cat("UPDATED  SHR PCR dup rate:", round(shr_upd$PCR_Duplicates / shr_upd$Sequenced_RP * 100, 1), "%\n\n")


################################################################
# PART 2: Compare old vs new (sb-option) loops
################################################################

# Read raw HICCUPS BEDPE rows
read_hiccups_bedpe <- function(file) {
  header_line <- readLines(file, n = 1)
  column_names <- str_split(str_remove(header_line, "^#"), "\t", simplify = TRUE) %>% as.character()

  readr::read_tsv(
    file,
    comment = "#",
    col_names = column_names,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  )
}

samples_to_compare <- c("592BB", "74AA", "A2DB", "D765A", "DA08A", "DA21A", "DA68A", "DBA9A", "DE8BA", "607")
key_cols <- c("chr1", "x1", "x2", "chr2", "y1", "y2")

for (sample_name in samples_to_compare) {
  cat("====================================================\n")
  cat(sprintf("=== Compare %s loops (old vs new)     ===\n", sample_name))
  cat("====================================================\n\n")

  if (sample_name == "74AA") {
    old_file <- sprintf("~/dropbox/Gateway_to_Hao/enhancer/data/loops/%s_intact_merged_loops5k10k25k.bedpe", sample_name)
  } else if (sample_name %in% c("A2DB", "592BB")) {
    old_file <- sprintf("~/dropbox/Gateway_to_Hao/enhancer/data/loops/%s_merged_loops_5k10k25k.bedpe", sample_name)
  } else {
    old_file <- sprintf("~/dropbox/Gateway_to_Hao/enhancer/data/loops/%s_intact_merged_loops_5k10k25k.bedpe", sample_name)
  }
  new_file <- sprintf("/Users/pete/Library/CloudStorage/GoogleDrive-wellclouder@gmail.com/My Drive/research/juicer-w-sb-options/%s/hiccups_5k10k25k/merged_loops.bedpe", sample_name)
  if (!file.exists(path.expand(old_file))) {
    cat(sprintf("Old file not found for %s: %s\n\n", sample_name, old_file))
    next
  }
  if (!file.exists(path.expand(new_file))) {
    cat(sprintf("New file not found for %s: %s\n\n", sample_name, new_file))
    next
  }

  # load data
  old_loops <- read_hiccups_bedpe(path.expand(old_file))
  new_loops <- read_hiccups_bedpe(path.expand(new_file))

  # extract keys
  old_keys <- old_loops %>%
    dplyr::select(all_of(key_cols)) %>%
    dplyr::distinct()
  new_keys <- new_loops %>%
    dplyr::select(all_of(key_cols)) %>%
    dplyr::distinct()

  # overlap loops
  overlap <- dplyr::inner_join(old_keys, new_keys, by = key_cols)

  # results
  n_old <- nrow(old_loops)
  n_new <- nrow(new_loops)
  n_overlap <- nrow(overlap)
  diff_count <- n_new - n_old
  diff_percent <- (diff_count / n_old) * 100
  overlap_percent <- (n_overlap / n_old) * 100

  cat("Total loops in old file:", n_old, "\n")
  cat("Total loops in new file:", n_new, "\n")
  cat(sprintf("Overlapping loops: %d (%.4f%% of old)\n", n_overlap, overlap_percent))
  cat(sprintf("Difference (new - old): %d (%.4f%% of old)\n\n", diff_count, diff_percent))
}

###################### DA08A loops (old vs new)
# Total loops in old file: 4656
# Total loops in new file: 4653
# Overlapping loops: 4650 (4650/4656 = 0.9987113)
###################### DA21A loops (old vs new)
# Total loops in old file: 5932
# Total loops in new file: 5930
# Overlapping loops: 5926 (5926/5932 = 0.9989885)
###################### DA68A loops (old vs new)
# Total loops in old file: 9131
# Total loops in new file: 9136
# Overlapping loops: 9121 (9121/9131 = 0.9989048)
###################### D765A loops (old vs new)
# Total loops in old file: 6568
# Total loops in new file: 6566
# Overlapping loops: 6563 (6563/6568 = 0.9992387)
###################### DBA9A loops (old vs new)
# Total loops in old file: 7676
# Total loops in new file: 7682
# Overlapping loops: 7668 (7668/7676 = 0.9989578)
###################### 74AA loops (old vs new)
# Total loops in old file: 2903
# Total loops in new file: 2902
# Overlapping loops: 2900 (2900/2903 = 0.9989666)
###################### A2DB loops (old vs new)
# Total loops in old file: 2992
# Total loops in new file: 2996
# Overlapping loops: 2978 (2978/2992 = 0.9953209)
###################### 592BB loops (old vs new)
# Total loops in old file: 5263
# Total loops in new file: 5264
# Overlapping loops: 5258 (5258/5263 = 0.9990500)
###################### 607 loops (old vs new)
# Total loops in old file: 7336 
# Total loops in new file: 4094 
# Overlapping loops: 1824 (24.8637% of old)
# Difference (new - old): -3242 (-44.1930% of old)
