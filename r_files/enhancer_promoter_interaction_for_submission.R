suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

##############################################################
# pre-requisites for submission files
##############################################################

# Check merged-vs-replicate HICCUPS loop files for 592 and 607
# This block is for manual provenance/QC checks before deciding whether
# to use a merged sample or an individual replicate.
medium_resolution_loop_dir <- "/Users/pete/Library/CloudStorage/GoogleDrive-wellclouder@gmail.com/My Drive/medium_resolution_5k10k25k"

replicate_check_files <- tibble::tribble(
  ~sample_family, ~sample_code, ~sample_type, ~loop_file, ~qc_file,
  "592", "592", "merged", file.path(medium_resolution_loop_dir, "592_inter_30.hiccups.5k10k25k", "merged_loops.bedpe"),
  "/Users/pete/Desktop/playground/enhancer/data/QC/592_intact_inter_30.txt",
  "592", "592AA", "replicate", file.path(medium_resolution_loop_dir, "592AA_inter_30.hiccups.5k10k25k", "merged_loops.bedpe"),
  "/Users/pete/Desktop/playground/enhancer/data/QC/592AA_intact_inter_30.txt",
  "592", "592BB", "replicate_used", file.path(medium_resolution_loop_dir, "592BB_inter_30.hiccups.5k10k25k", "merged_loops.bedpe"),
  "/Users/pete/Desktop/playground/enhancer/data/QC/592BB_intact_inter_30.txt",
  "607", "607", "merged_used", file.path(medium_resolution_loop_dir, "607_inter_30.hiccups.5k10k25k", "merged_loops.bedpe"),
  "/Users/pete/Desktop/playground/enhancer/data/QC/607_intact_inter_30.txt",
  "607", "607BB", "replicate", file.path(medium_resolution_loop_dir, "607BB_inter_30.hiccups.5k10k25k", "merged_loops.bedpe"),
  "/Users/pete/Desktop/playground/enhancer/data/QC/607BB_intact_inter_30.txt",
  "607", "607CC", "replicate", file.path(medium_resolution_loop_dir, "607CC_inter_30.hiccups.5k10k25k", "merged_loops.bedpe"),
  "/Users/pete/Desktop/playground/enhancer/data/QC/607CC_intact_inter_30.txt"
)

analysis_loop_files_to_check <- tibble::tribble(
  ~sample_family, ~sample_code, ~analysis_file,
  "592", "592BB", "/Users/pete/dropbox/Gateway_to_Hao/enhancer/data/loops/592BB_merged_loops_5k10k25k.bedpe",
  "607", "607", "/Users/pete/dropbox/Gateway_to_Hao/enhancer/data/loops/607_intact_merged_loops_5k10k25k.bedpe"
)

# Calculate a SHA-256 checksum so source and analysis loop files can be compared byte-for-byte.
sha256_file <- function(file) {
  if (!file.exists(file)) {
    return(NA_character_)
  }

  str_split(system2("shasum", c("-a", "256", shQuote(file)), stdout = TRUE), "\\s+", simplify = TRUE)[1]
}

# Read one HICCUPS BEDPE file for replicate checks and add resolution plus a stable loop key.
read_hiccups_loop_check_bedpe <- function(file) {
  header_line <- readLines(file, n = 1)
  column_names <- str_split(str_remove(header_line, "^#"), "\t", simplify = TRUE) %>% as.character()

  readr::read_tsv(
    file,
    comment = "#",
    col_names = column_names,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  ) %>%
    mutate(
      x1 = as.integer(x1),
      x2 = as.integer(x2),
      y1 = as.integer(y1),
      y2 = as.integer(y2),
      resolution = x2 - x1,
      loop_key = str_c(chr1, x1, x2, chr2, y1, y2, sep = "_")
    )
}

# Pull one requested metric line, such as PCR duplicates, from a Juicer QC text report.
parse_qc_stat_line <- function(lines, pattern) {
  line <- lines[str_detect(lines, fixed(pattern))]
  if (length(line) == 0) {
    return(NA_character_)
  }

  str_squish(str_remove(line[1], paste0("^\\s*", pattern, ":\\s*")))
}

# Summarize the Juicer inter_30 QC metrics needed to compare merged and replicate samples.
read_juicer_qc_summary <- function(file) {
  if (!file.exists(file)) {
    return(tibble(
      qc_file_exists = FALSE,
      sequenced_read_pairs = NA_character_,
      normal_paired = NA_character_,
      pcr_duplicates = NA_character_,
      library_complexity_estimate = NA_character_,
      hic_contacts = NA_character_
    ))
  }

  lines <- readLines(file, warn = FALSE)
  tibble(
    qc_file_exists = TRUE,
    sequenced_read_pairs = parse_qc_stat_line(lines, "Sequenced Read Pairs"),
    normal_paired = parse_qc_stat_line(lines, "Normal Paired"),
    pcr_duplicates = parse_qc_stat_line(lines, "PCR Duplicates"),
    library_complexity_estimate = parse_qc_stat_line(lines, "Library Complexity Estimate"),
    hic_contacts = parse_qc_stat_line(lines, "Hi-C Contacts")
  )
}

replicate_loop_data <- replicate_check_files %>%
  mutate(
    loop_file_exists = file.exists(loop_file),
    qc_file_exists = file.exists(qc_file),
    sha256 = map_chr(loop_file, sha256_file),
    data = map(loop_file, read_hiccups_loop_check_bedpe)
  )

loop_qc_summary <- replicate_loop_data %>%
  transmute(
    sample_family,
    sample_code,
    sample_type,
    loop_file,
    loop_file_exists,
    sha256,
    n_loops = map_int(data, nrow),
    n_unique_loop_keys = map_int(data, ~ n_distinct(.x$loop_key)),
    resolution_distribution = map_chr(data, ~ .x %>%
      count(resolution, name = "n") %>%
      arrange(resolution) %>%
      mutate(label = str_c(resolution, "=", n)) %>%
      pull(label) %>%
      str_c(collapse = "; ")),
    qc_file
  ) %>%
  bind_cols(map_dfr(replicate_check_files$qc_file, read_juicer_qc_summary))

loop_overlap_summary <- replicate_loop_data %>%
  dplyr::select(sample_family, sample_code, data) %>%
  group_by(sample_family) %>%
  group_modify(~ {
    pair_grid <- t(combn(.x$sample_code, 2)) %>% as_tibble(.name_repair = "minimal")
    names(pair_grid) <- c("sample_a", "sample_b")

    pair_grid %>%
      rowwise() %>%
      mutate(
        loops_a = list(.x$data[[match(sample_a, .x$sample_code)]]$loop_key),
        loops_b = list(.x$data[[match(sample_b, .x$sample_code)]]$loop_key),
        n_a = length(loops_a),
        n_b = length(loops_b),
        n_intersect = length(intersect(loops_a, loops_b)),
        n_a_only = length(setdiff(loops_a, loops_b)),
        n_b_only = length(setdiff(loops_b, loops_a)),
        identical_loop_sets = setequal(loops_a, loops_b)
      ) %>%
      ungroup() %>%
      dplyr::select(-loops_a, -loops_b)
  }) %>%
  ungroup()

analysis_file_identity_summary <- analysis_loop_files_to_check %>%
  mutate(
    analysis_file_exists = file.exists(analysis_file),
    analysis_sha256 = map_chr(analysis_file, sha256_file)
  ) %>%
  left_join(
    loop_qc_summary %>%
      dplyr::select(sample_family, sample_code, source_sha256 = sha256, source_loop_file = loop_file),
    by = c("sample_family", "sample_code")
  ) %>%
  mutate(analysis_file_identical_to_source = analysis_sha256 == source_sha256)

print(loop_qc_summary, n = Inf, width = Inf)
print(loop_overlap_summary, n = Inf, width = Inf)
print(analysis_file_identity_summary, n = Inf, width = Inf)

##############################################################
# submission files
##############################################################

final_loop_rda <- path.expand("~/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/lt2mb/df_final_loop_sub.4.any.lt2mb.ENSEMBL.mid.mid.final.filter.200kb.rda")

if (!exists("df.final.loop")) {
  load(final_loop_rda)
}

# 1. loops
write.csv(df.final.loop, file = "./figures/submission/lt2mb/df_final_loop_sub.4.any.lt2mb.ENSEMBL.mid.mid.final.filter.200kb.csv", row.names = FALSE)

# 1-1. Supplementary Data 1: Raw HICCUPS loop results with functional-loop annotation
# Combine all 10 sample BEDPE files into one table with sample/strain metadata
# and flag each loop as functional (TRUE) if it appears in the final curated loop set.
loop_file_metadata <- tibble::tribble(
  ~sample, ~strain, ~file,
  "592BB", "SHR/OlaIpcv", "~/dropbox/Gateway_to_Hao/enhancer/data/loops/592BB_merged_loops_5k10k25k.bedpe",
  "607", "HXB10", "~/dropbox/Gateway_to_Hao/enhancer/data/loops/607_intact_merged_loops_5k10k25k.bedpe",
  "74AA", "F344/Stm", "~/dropbox/Gateway_to_Hao/enhancer/data/loops/74AA_intact_merged_loops5k10k25k.bedpe",
  "A2DB", "LE/Stm", "~/dropbox/Gateway_to_Hao/enhancer/data/loops/A2DB_merged_loops_5k10k25k.bedpe",
  "D765A", "BXH6", "~/dropbox/Gateway_to_Hao/enhancer/data/loops/D765A_intact_merged_loops_5k10k25k.bedpe",
  "DA08A", "HXB2", "~/dropbox/Gateway_to_Hao/enhancer/data/loops/DA08A_intact_merged_loops_5k10k25k.bedpe",
  "DA21A", "SHR/OlaIpcvxBN/NHsdMcwi", "~/dropbox/Gateway_to_Hao/enhancer/data/loops/DA21A_intact_merged_loops_5k10k25k.bedpe",
  "DA68A", "HXB31", "~/dropbox/Gateway_to_Hao/enhancer/data/loops/DA68A_intact_merged_loops_5k10k25k.bedpe",
  "DBA9A", "HXB23", "~/dropbox/Gateway_to_Hao/enhancer/data/loops/DBA9A_intact_merged_loops_5k10k25k.bedpe",
  "DE8BA", "BN-Lx", "~/dropbox/Gateway_to_Hao/enhancer/data/loops/DE8BA_intact_merged_loops_5k10k25k.bedpe"
) %>%
  mutate(file = path.expand(file))

# Read raw HICCUPS BEDPE rows preserving all original columns.
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

final_loop_ids <- df.final.loop %>%
  distinct(loop.id) %>%
  pull(loop.id)

supplementary.data1.loops <- loop_file_metadata %>%
  mutate(data = purrr::map(file, read_hiccups_bedpe)) %>%
  dplyr::select(-file) %>%
  tidyr::unnest(data) %>%
  mutate(
    loop.id = str_c(chr1, x1, x2, chr2, y1, y2, as.numeric(x2) - as.numeric(x1), sep = "_"),
    functional = loop.id %in% final_loop_ids
  ) %>%
  dplyr::select(sample, strain, functional, dplyr::everything(), -loop.id)

readr::write_tsv(
  supplementary.data1.loops,
  file = "/Users/pete/UTHSC GGI Dropbox/K P/Gateway_to_Hao/publication/PE-interaction/submission/journal/data/supplementary_data1_loops.tsv"
)

# 2. CTCF
ctcf_file_path <- "~/dropbox/Gateway_to_Hao/enhancer/data/ctcf/submission/E4/fimo_E4_submission_trial.txt"
supplementary.data2.ctcf <- read.table(file = path.expand(ctcf_file_path), header = TRUE, sep = "\t")

readr::write_tsv(
  supplementary.data2.ctcf,
  file = "/Users/pete/UTHSC GGI Dropbox/K P/Gateway_to_Hao/publication/PE-interaction/submission/journal/data/supplementary_data2_CTCF.tsv.gz"
)

########################################################
# 3. motif tables for submission
########################################################
submission_ctcf_dir <- "/Users/pete/UTHSC GGI Dropbox/K P/Gateway_to_Hao/enhancer/data/ctcf/submission"
submission_for_table_dir <- file.path(submission_ctcf_dir, "for_table")
dir.create(submission_for_table_dir, recursive = TRUE, showWarnings = FALSE)

# Read motif metadata from a MEME file and summarize one row per motif for table export.
read_meme_motif_table <- function(file, db_name) {
  lines <- readLines(file, warn = FALSE)
  motif_idx <- which(str_detect(lines, "^MOTIF\\s+"))

  map_dfr(motif_idx, function(idx) {
    motif_name <- str_split(str_squish(lines[idx]), "\\s+", simplify = TRUE)[2]

    width_idx <- which(seq_along(lines) > idx & str_detect(lines, "^letter-probability matrix:"))[1]
    width <- as.integer(str_match(lines[width_idx], "w=\\s*([0-9]+)")[, 2])

    matrix_lines <- lines[seq.int(width_idx + 1, width_idx + width)]
    consensus <- map_chr(matrix_lines, function(line) {
      probs <- as.numeric(str_split(str_squish(line), "\\s+", simplify = TRUE))
      c("A", "C", "G", "T")[which.max(probs)]
    }) %>%
      str_c(collapse = "")

    url_idx <- which(seq_along(lines) > (width_idx + width) & str_detect(lines, "^URL\\s+"))[1]
    url <- str_remove(lines[url_idx], "^URL\\s+")

    tibble(
      DB = db_name,
      MOTIF = motif_name,
      WIDTH = width,
      `BEST POSSIBLE MATCH` = consensus,
      URL = url
    )
  })
}

ctcfbsdb_table <- read_meme_motif_table(
  file.path(submission_ctcf_dir, "ctcfbsdb", "source_checked_ctcf_uthsc.meme"),
  db_name = "ctcfbsdb"
)

readr::write_csv(
  ctcfbsdb_table,
  file.path(submission_for_table_dir, "ctcfbsdb_motif_table.csv")
)
