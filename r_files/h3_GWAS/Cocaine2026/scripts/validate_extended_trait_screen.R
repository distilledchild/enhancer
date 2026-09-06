#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
external <- Sys.getenv("ENHANCER_EXTERNAL", "/Volumes/external_1000GB_all/playground/enhancer")
root <- file.path(external, "r_files/h3_GWAS/Cocaine2026")
out <- file.path(root, "extended_trait_screen_2001bp")
zip <- file.path(root, "data/bb2334903t_4_1.zip")
result <- fread(file.path(out, "all_29_trait_peak_screen.tsv"))
inventory <- fread(file.path(out, "archive_MLMA_inventory.tsv"))
published <- fread(file.path(out, "source_report_lead_SNPs.tsv"))
stopifnot(nrow(result) == 29L, nrow(inventory) == 638L, !anyDuplicated(result$trait))
all_summaries <- rbindlist(lapply(result$trait, function(t)
  fread(file.path(out, "chromosome_screens", paste0(t, ".tsv")))))
stopifnot(nrow(all_summaries) == 638L, all(all_summaries$headers_read == 22L),
  all(all_summaries$n_invalid_P == 0L), all(all_summaries$n_zero_P == 0L),
  sum(all_summaries$n_rows) == sum(result$SNP_rows))

repo <- Sys.getenv("ENHANCER_REPO", "/Users/pete/Desktop/playground/enhancer")
awk <- file.path(repo, "r_files/h3_GWAS/Cocaine2026/scripts/scan_mlma_peaks.awk")
test_scanner <- function() {
  fixture <- tempfile(fileext = ".mlma")
  output <- tempfile(fileext = ".tsv")
  on.exit(unlink(c(fixture, output)))
  header <- "Chr SNP bp A1 A2 Freq b se p"
  writeLines(c(header,
    "1 1:10 10 A G 0.2 1 1 1e-7", "1 1:11 11 A G 0.2 1 1 1.00001e-7",
    "1 1:12 12 A G 0.2 1 1 NA", "1 1:13 13 A G 0.2 1 1 0",
    "1 1:14 14 A G 0.2 1 1 -0.1", "1 1:15 15 A G 0.2 1 1 1.1",
    "1 1:16 16 A G 0.2 1 1 invalid", header,
    "2 2:10 10 A G 0.2 1 1 1e-8", "2 2:11 11 A G 0.2 1 1 1e-8"), fixture)
  status <- system2("awk", c("-v", "trait=synthetic", "-f", shQuote(awk), shQuote(fixture)), stdout = output)
  x <- fread(file = output)
  a <- x[Chr == 1L]; b <- x[Chr == 2L]
  stopifnot(status == 0L, nrow(x) == 2L, all(x$headers_read == 2L),
    a$n_rows == 7L, a$n_valid == 2L, a$n_missing_P == 1L, a$n_zero_P == 1L,
    a$n_invalid_P == 3L, a$n_P_le_1e_minus7 == 1L, a$min_P == 1e-7,
    b$n_rows == 2L, b$lead_ties == 2L, b$n_P_le_1e_minus7 == 2L, b$min_P == 1e-8)
}
test_scanner()

# Independently parse complete files with R, without the awk screening program.
read_member <- function(member) {
  tmp <- tempfile(fileext = ".mlma")
  on.exit(unlink(tmp))
  status <- system2("unzip", c("-p", shQuote(zip), shQuote(member)), stdout = tmp)
  stopifnot(status == 0L)
  d <- fread(file = tmp)
  stopifnot(identical(names(d), c("Chr", "SNP", "bp", "A1", "A2", "Freq", "b", "se", "p")))
  stopifnot(!anyDuplicated(d$SNP), all(is.finite(d$p)), all(d$p > 0 & d$p <= 1))
  chr <- as.character(d$Chr)
  chr[chr == "21"] <- "X"
  chr[chr == "24"] <- "MT"
  stopifnot(all(d$SNP == paste0(chr, ":", d$bp)))
  d
}

check_traits <- unique(c(result$trait[which.min(result$min_P_nuclear)],
                         "regressedlr_pr_01_active", "regressedlr_pr_max"))
checks <- list()
for (t in check_traits) {
  message("Independent full-trait check: ", t)
  rows <- inventory[trait == t]
  source_scan <- fread(file.path(out, "chromosome_screens", paste0(t, ".tsv")))
  for (i in seq_len(nrow(rows))) {
    d <- read_member(rows$Name[i])
    stopifnot(uniqueN(d$Chr) == 1L)
    s <- source_scan[Chr == d$Chr[1]]
    stopifnot(nrow(s) == 1L, s$n_rows == nrow(d), s$n_valid == nrow(d),
              s$min_P == min(d$p), s$n_P_le_1e_minus7 == sum(d$p <= 1e-7),
              s$lead_SNP %chin% d[p == min(p), SNP])
    checks[[length(checks) + 1L]] <- data.table(trait = t, member = rows$Name[i],
      rows = nrow(d), min_P = min(d$p), rows_P_le_1e_minus7 = sum(d$p <= 1e-7),
      independent_minimum_matches = TRUE, unique_SNP_IDs = TRUE, SNP_ID_coordinates_match = TRUE)
  }
}
fwrite(rbindlist(checks), file.path(out, "independent_R_peak_checks.tsv"), sep = "\t")

lead_checks <- list()
for (i in seq_len(nrow(published))) {
  t <- paste0("regressedlr_", published$trait[i])
  id <- published$TopSNP[i]
  chr <- tolower(sub(":.*$", "", id))
  member <- inventory[trait == t & chromosome_file == chr, Name]
  stopifnot(length(member) == 1L)
  d <- read_member(member)[SNP == id]
  stopifnot(nrow(d) == 1L)
  discrepancy <- abs(-log10(d$p) - published[["-Log10(p)"]][i])
  stopifnot(discrepancy <= .00051)
  lead_checks[[i]] <- data.table(trait = t, SNP = id, raw_P = d$p,
    raw_log10P = -log10(d$p), report_log10P = published[["-Log10(p)"]][i],
    absolute_difference = discrepancy, matches_report_rounding = TRUE)
}
fwrite(rbindlist(lead_checks), file.path(out, "source_report_lead_concordance.tsv"), sep = "\t")
q <- data.table(check = c("traits_screened", "MLMA_files_screened", "trait_SNP_rows_screened",
  "traits_with_matching_source_N", "traits_without_source_N", "invalid_or_zero_P_detected",
  "missing_P_rows", "scanner_synthetic_boundary_test_passed",
  "independently_reparsed_full_traits", "independently_reparsed_chromosome_files",
  "source_report_leads_confirmed", "threshold7_nuclear_candidates", "threshold7_all_chromosome_candidates"),
  value = c(nrow(result), sum(result$files_scanned), sum(result$SNP_rows),
    sum(result$N_verified), sum(!result$N_verified),
    sum(all_summaries$n_invalid_P + all_summaries$n_zero_P),
    sum(all_summaries$n_missing_P), 1L, length(check_traits), length(checks),
    length(lead_checks), sum(result$passes_nuclear_7), sum(result$passes_all_chromosomes_7)))
fwrite(q, file.path(out, "validation_summary.tsv"), sep = "\t")
scripts <- file.path(repo, "r_files/h3_GWAS/Cocaine2026/scripts",
  c("scan_mlma_peaks.awk", "screen_extended_traits.R", "validate_extended_trait_screen.R"))
fwrite(data.table(path = scripts, md5 = unname(tools::md5sum(scripts))),
  file.path(out, "script_manifest.tsv"), sep = "\t")
print(q)
