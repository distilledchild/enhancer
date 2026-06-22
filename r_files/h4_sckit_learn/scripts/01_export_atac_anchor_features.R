####################################################
# Export ATAC enhancer anchor features
# Build the first scikit-learn input table
# Label enhancer anchors by ATAC overlap status
# Join loop-level regulatory genomics features
####################################################
library(tidyverse)

options(scipen = 999) # prevent scientific notation in base R
options(pillar.sigfig = 10) # display up to 10 significant digits in tibbles

setwd("/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn") # Mac
getwd() # /Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn

path.csv.atac.detail <- "/Users/pete/Desktop/playground/enhancer/r_files/atac_validation/atac_loop_anchor_overlap_detail.csv"
path.rds.loop.features <- "/Users/pete/Library/CloudStorage/Dropbox-UTHSCGGI/K P/Gateway_to_Hao/enhancer/r_files/rds/df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3_final_200kb.rds"
path.dir.output <- "/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/data"
path.csv.output <- file.path(path.dir.output, "atac_enhancer_anchor_features.csv")

stopifnot(file.exists(path.csv.atac.detail))
stopifnot(file.exists(path.rds.loop.features))
dir.create(path.dir.output, recursive = TRUE, showWarnings = FALSE)

# func1: Collapse unique non-missing values into one semicolon-separated string
collapse_unique <- function(vec.input) {
  vec.value <- vec.input[!is.na(vec.input) & vec.input != ""]
  vec.value <- sort(unique(vec.value))

  if (length(vec.value) == 0) {
    return(NA_character_)
  }

  str_c(vec.value, collapse = ";")
}

####################################################
# 1. Load ATAC anchor overlap labels
####################################################
df.atac.detail <- read_csv(path.csv.atac.detail, show_col_types = FALSE)

df.atac.enhancer <- df.atac.detail %>%
  filter(anchor_type == "enhancer") %>%
  mutate(
    atac_overlap_text = str_to_lower(as.character(atac_overlap)),
    label_atac_overlap = case_when(
      atac_overlap_text %in% c("true", "1") ~ 1L,
      atac_overlap_text %in% c("false", "0") ~ 0L,
      TRUE ~ NA_integer_
    )
  ) %>%
  transmute(
    loop_id,
    label_atac_overlap,
    atac_category = category
  ) %>%
  arrange(loop_id)

stopifnot(!any(is.na(df.atac.enhancer$label_atac_overlap)))
stopifnot(anyDuplicated(df.atac.enhancer$loop_id) == 0)

print(str_glue("enhancer anchors: {nrow(df.atac.enhancer)}")) # 15085
print(str_glue("positive labels: {sum(df.atac.enhancer$label_atac_overlap == 1)}")) # 13428

####################################################
# 2. Collapse RDS loop annotations to one row per loop
####################################################
df.loop.raw <- read_rds(path.rds.loop.features)

df.loop.features <- df.loop.raw %>%
  group_by(loop.id) %>%
  summarise(
    chr1 = first(chr1),
    x1 = first(x1),
    x2 = first(x2),
    chr2 = first(chr2),
    y1 = first(y1),
    y2 = first(y2),
    resolution = first(resolution),
    distance = first(distance),
    loop_chr = first(loop_chr),
    loop_start = first(loop_start),
    loop_end = first(loop_end),
    WHERE = collapse_unique(WHERE),
    classification = collapse_unique(classification),
    n_assigned_genes = n_distinct(gene_name, na.rm = TRUE),
    n_components = n_distinct(component, na.rm = TRUE),
    has_tss_component = as.integer(any(component == "tss", na.rm = TRUE)),
    has_promoter_component = as.integer(any(component == "pro", na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(across(c(x1, x2, y1, y2, resolution), parse_integer)) %>%
  mutate(
    anchor1_width = x2 - x1,
    anchor2_width = y2 - y1,
    same_chr = as.integer(chr1 == chr2),
    loop_span = if_else(same_chr == 1L, y2 - x1, NA_real_)
  ) %>%
  arrange(loop.id)

stopifnot(anyDuplicated(df.loop.features$loop.id) == 0)
stopifnot(all(df.loop.features$anchor1_width > 0))
stopifnot(all(df.loop.features$anchor2_width > 0))

print(str_glue("loop feature rows: {nrow(df.loop.features)}")) # 17648

####################################################
# 3. Join labels with loop-level features
####################################################
df.features <- df.atac.enhancer %>%
  left_join(df.loop.features, by = c("loop_id" = "loop.id")) %>%
  select(
    loop_id,
    label_atac_overlap,
    atac_category,
    chr1,
    x1,
    x2,
    chr2,
    y1,
    y2,
    resolution,
    distance,
    anchor1_width,
    anchor2_width,
    loop_span,
    same_chr,
    n_assigned_genes,
    n_components,
    has_tss_component,
    has_promoter_component,
    WHERE,
    classification
  ) %>%
  arrange(loop_id)

num.missing.features <- sum(is.na(df.features$distance))
stopifnot(num.missing.features == 0)
stopifnot(nrow(df.features) == nrow(df.atac.enhancer))
stopifnot(anyDuplicated(df.features$loop_id) == 0)

print(str_glue("missing feature rows: {num.missing.features}")) # 0

####################################################
# 4. Save scikit-learn feature table
####################################################
write_csv(df.features, path.csv.output)

print(str_glue("wrote file: {path.csv.output}")) # /Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/data/atac_enhancer_anchor_features.csv
print(str_glue("wrote rows: {nrow(df.features)}")) # 15085
