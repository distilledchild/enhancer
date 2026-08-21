##############################################################
# Recreate library_complexity.tsv with 592BB QC for SHR/OlaIpcv
##############################################################

data_dir <- "/Users/pete/UTHSC GGI Dropbox/K P/Gateway_to_Hao/enhancer/data"
qc_dir <- "/Users/pete/Desktop/playground/enhancer/data/QC"

output_file <- file.path(data_dir, "library_complexity_592BB.tsv")
original_file <- file.path(data_dir, "library_complexity.tsv")
comparison_file <- file.path(data_dir, "library_complexity_vs_592BB_comparison.tsv")

sample_qc_files <- data.frame(
  Strain = c(
    "SHR/OlaIpcv",
    "BN-Lx",
    "BXH6",
    "HXB2",
    "HXB10",
    "HXB23",
    "HXB31",
    "LE/Stm",
    "F344/Stm",
    "SHR/OlaIpcvxBN/NHsdMcwi"
  ),
  sample_code = c(
    "592BB",
    "DE8BA",
    "D765A",
    "DA08A",
    "607",
    "DBA9A",
    "DA68A",
    "A2DB",
    "74AA",
    "DA21A"
  ),
  stringsAsFactors = FALSE
)

metric_patterns <- c(
  Sequenced_RP = "Sequenced Read Pairs",
  Normal_Paired = "Normal Paired",
  Chimeric_Paired = "Chimeric Paired",
  Chimeric_Ambiguous = "Chimeric Ambiguous",
  Unmapped = "Unmapped",
  Alignable_Normal_N_Chimeric = "Alignable \\(Normal\\+Chimeric Paired\\)",
  Unique_Reads = "Unique Reads",
  PCR_Duplicates = "PCR Duplicates",
  Optical_Duplicates = "Optical Duplicates",
  Below_MAPQ_Threshold = "Below MAPQ Threshold",
  `Hi-C_Contacts` = "Hi-C Contacts",
  `Inter-chromosomal` = "Inter-chromosomal",
  `Intra-chromosomal` = "Intra-chromosomal",
  Short_Range_20Kb = "Short Range \\(<20Kb\\)",
  Long_Range_20Kb = "Long Range \\(>20Kb\\)"
)

parse_juicer_count <- function(lines, metric_pattern) {
  metric_line <- grep(paste0("^\\s*", metric_pattern, ":"), lines, value = TRUE)

  if (length(metric_line) == 0) {
    stop("Missing metric: ", metric_pattern, call. = FALSE)
  }

  count_text <- sub("^\\s*[^:]+:\\s*([0-9,]+).*", "\\1", metric_line[1])
  as.integer(gsub(",", "", count_text))
}

read_library_complexity_row <- function(strain, sample_code) {
  qc_file <- file.path(qc_dir, paste0(sample_code, "_intact_inter_30.txt"))

  if (!file.exists(qc_file)) {
    stop("Missing QC file: ", qc_file, call. = FALSE)
  }

  lines <- readLines(qc_file, warn = FALSE)
  metric_values <- vapply(metric_patterns, parse_juicer_count, integer(1), lines = lines)

  data.frame(
    Strain = strain,
    as.data.frame(as.list(metric_values), check.names = FALSE),
    check.names = FALSE
  )
}

library_complexity_592BB <- do.call(
  rbind,
  Map(read_library_complexity_row, sample_qc_files$Strain, sample_qc_files$sample_code)
)

write.table(
  library_complexity_592BB,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

library_complexity_original <- read.delim(original_file, check.names = FALSE)

comparison <- merge(
  reshape(
    library_complexity_original,
    varying = names(library_complexity_original)[names(library_complexity_original) != "Strain"],
    v.names = "original_value",
    timevar = "metric",
    times = names(library_complexity_original)[names(library_complexity_original) != "Strain"],
    idvar = "Strain",
    direction = "long"
  ),
  reshape(
    library_complexity_592BB,
    varying = names(library_complexity_592BB)[names(library_complexity_592BB) != "Strain"],
    v.names = "library_complexity_592BB_value",
    timevar = "metric",
    times = names(library_complexity_592BB)[names(library_complexity_592BB) != "Strain"],
    idvar = "Strain",
    direction = "long"
  ),
  by = c("Strain", "metric"),
  sort = FALSE
)

comparison$difference <- comparison$library_complexity_592BB_value - comparison$original_value
comparison$changed <- comparison$difference != 0
comparison <- comparison[order(match(comparison$Strain, sample_qc_files$Strain)), ]

write.table(
  comparison,
  file = comparison_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(library_complexity_592BB)
cat("\nChanged cells by strain:\n")
print(table(comparison$Strain[comparison$changed]))
cat("\nWrote:\n")
cat(output_file, "\n")
cat(comparison_file, "\n")
