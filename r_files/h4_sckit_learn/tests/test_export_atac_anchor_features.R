####################################################
# Test ATAC enhancer anchor feature export
# Verifies the scikit-learn input table schema
# Confirms one enhancer-anchor row per loop
# Checks binary ATAC overlap labels
####################################################
suppressPackageStartupMessages({
  library(tidyverse)
})

options(scipen = 999) # prevent scientific notation in base R
options(pillar.sigfig = 10) # display up to 10 significant digits in tibbles

setwd("/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn") # Mac
getwd() # /Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn

path.csv.features <- "/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/data/atac_enhancer_anchor_features.csv"

####################################################
# 1. Validate feature table output
####################################################
if (!file.exists(path.csv.features)) {
  stop(str_glue("Missing feature table: {path.csv.features}"))
}

df.features <- read_csv(path.csv.features, show_col_types = FALSE)

vec.required.columns <- c(
  "loop_id",
  "label_atac_overlap",
  "atac_category",
  "chr1",
  "x1",
  "x2",
  "chr2",
  "y1",
  "y2",
  "resolution",
  "distance",
  "anchor1_width",
  "anchor2_width",
  "loop_span",
  "same_chr",
  "n_assigned_genes",
  "n_components",
  "has_tss_component",
  "has_promoter_component",
  "WHERE",
  "classification"
)

vec.missing.columns <- setdiff(vec.required.columns, names(df.features))
if (length(vec.missing.columns) > 0) {
  stop(str_glue("Missing columns: {str_c(vec.missing.columns, collapse = ', ')}"))
}

stopifnot(nrow(df.features) == 15085)
stopifnot(anyDuplicated(df.features$loop_id) == 0)
stopifnot(all(df.features$label_atac_overlap %in% c(0, 1)))
stopifnot(!any(is.na(df.features$distance)))
stopifnot(all(df.features$anchor1_width > 0))
stopifnot(all(df.features$anchor2_width > 0))

print(str_glue("feature rows: {nrow(df.features)}")) # 15085
print(str_glue("positive labels: {sum(df.features$label_atac_overlap == 1)}")) # 13428
