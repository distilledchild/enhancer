#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(xml2); library(jsonlite)})
setDTthreads(1)
repo <- Sys.getenv("ENHANCER_REPO", "/Users/pete/Desktop/playground/enhancer")
external <- Sys.getenv("ENHANCER_EXTERNAL", "/Volumes/external_1000GB_all/playground/enhancer")
source_dir <- file.path(repo, "r_files/h3_GWAS/Cocaine2026/scripts")
root <- file.path(external, "r_files/h3_GWAS/Cocaine2026")
out <- file.path(root, "extended_trait_screen_2001bp")
dir.create(file.path(out, "chromosome_screens"), recursive = TRUE, showWarnings = FALSE)
zip <- file.path(root, "data/bb2334903t_4_1.zip")
report <- file.path(root, "data/bb2334903t_2_1.html")
dictionary_file <- file.path(root, "data/bb2334903t_1_1.csv")
phenotype_file <- file.path(root, "data/bb2334903t_3_1.csv")
awk <- file.path(source_dir, "scan_mlma_peaks.awk")
inputs <- c(zip, report, dictionary_file, phenotype_file, awk)
stopifnot(all(file.exists(inputs)))
fwrite(data.table(path = inputs, bytes = file.info(inputs)$size, md5 = unname(tools::md5sum(inputs))),
       file.path(out, "input_manifest.tsv"), sep = "\t")
inventory <- as.data.table(unzip(zip, list = TRUE))
m <- inventory[grepl("^gwas/regressedlr_.*_chrgwas(.*)[.]mlma$", Name)]
m[, trait := sub("_chrgwas.*$", "", basename(Name))]
m[, chromosome_file := sub("[.]mlma$", "", sub("^.*_chrgwas", "", Name))]
stopifnot(!anyDuplicated(m$Name), nrow(m) == sum(endsWith(inventory$Name, ".mlma")))
fwrite(m, file.path(out, "archive_MLMA_inventory.tsv"), sep = "\t")
traits <- sort(unique(m$trait))
for (t in traits) stopifnot(setequal(m[trait == t, chromosome_file], c(as.character(1:20), "x", "mt")))
dictionary <- fread(dictionary_file)
pheno <- fread(phenotype_file)
stopifnot(!anyDuplicated(pheno$rfid), all(traits %in% names(pheno)))

# Parse embedded Bokeh tables rather than execute the report's scripts.
decode_array <- function(v) {
  if (!is.list(v) || !identical(v$type, "ndarray")) return(unlist(v))
  if (v$dtype == "object") return(unlist(v$array))
  stopifnot(identical(v$array$type, "bytes"), v$dtype %in% c("float64", "int32"))
  readBin(base64_dec(v$array$data), if (v$dtype == "float64") "double" else "integer",
          n = prod(unlist(v$shape)), size = if (v$dtype == "float64") 8 else 4, endian = v$order)
}
heritability <- list(); lead_tables <- list(); contexts <- character()
walk <- function(v) {
  if (!is.list(v)) return()
  if (identical(v$name, "panel.models.markup.HTML")) {
    txt <- v$attributes$text
    if (!is.null(txt) && nchar(txt) < 40000 && grepl("general-information|summary-of-qtls|snp-heritability-estimates", txt)) {
      for (pass in 1:3) txt <- xml_text(read_html(paste0("<div>", txt, "</div>")))
      contexts <<- c(contexts, txt)
    }
    return()
  }
  if (identical(v$name, "ColumnDataSource")) {
    e <- v$attributes$data$entries
    keys <- vapply(e, function(z) z[[1]], "")
    if ("trait" %in% keys && any(c("n", "TopSNP") %in% keys)) {
      d <- as.data.table(setNames(lapply(e, function(z) decode_array(z[[2]])), keys))
      if ("n" %in% keys) heritability[[length(heritability) + 1L]] <<- d
      if ("TopSNP" %in% keys) lead_tables[[length(lead_tables) + 1L]] <<- d
    }
    return()
  }
  for (a in v) if (is.list(a)) walk(a)
}
doc <- read_html(report)
for (s in xml_find_all(doc, "//script[@type='application/json']")) walk(fromJSON(xml_text(s), simplifyVector = FALSE))
h <- unique(rbindlist(heritability, fill = TRUE))
stopifnot(nrow(h) > 0L)
h[, trait := paste0("regressedlr_", sub("^regressedlr_", "", trimws(trait)))]
stopifnot(!anyDuplicated(h$trait))
fwrite(h, file.path(out, "source_heritability_table.tsv"), sep = "\t")
fwrite(unique(rbindlist(lead_tables, fill = TRUE)), file.path(out, "source_report_lead_SNPs.tsv"), sep = "\t")
writeLines(unique(contexts), file.path(out, "source_report_metadata.txt"))

meta <- data.table(trait = traits)
meta[, N_phenotype := vapply(trait, function(t) sum(is.finite(pheno[[t]])), integer(1))]
meta[, original_trait := sub("^regressedlr_", "", trait)]
meta <- merge(meta, dictionary[datatype == "trait" & nzchar(alternative_name),
  .(original_trait = alternative_name, trait_label = variable, description)], by = "original_trait", all.x = TRUE)
meta <- merge(meta, h[, .(trait, N_report = as.integer(n), heritability, heritability_se)], by = "trait", all.x = TRUE)
meta[, N_verified := !is.na(N_report) & N_report == N_phenotype]
meta[, previous_four := trait %chin% c("regressedlr_sha_mean_to_01_03", "regressedlr_pc1_lga",
                                      "regressedlr_lga_total_intake", "regressedlr_shock_03_calculated")]
meta[, dictionary_documented := !is.na(trait_label)]
fwrite(meta, file.path(out, "trait_metadata.tsv"), sep = "\t")

screens <- list()
for (i in seq_along(traits)) {
  t <- traits[i]
  message(format(Sys.time(), "%H:%M:%S"), " SCREEN ", i, "/", length(traits), " ", t)
  members <- m[trait == t][order(match(chromosome_file, c(as.character(1:20), "x", "mt"))), Name]
  target <- file.path(out, "chromosome_screens", paste0(t, ".tsv"))
  # Pipefail also checks unzip's CRC/decompression status, not just awk's exit code.
  cmd <- paste("unzip -p", shQuote(zip), paste(shQuote(members), collapse = " "),
               "| LC_ALL=C awk -v", shQuote(paste0("trait=", t)), "-f", shQuote(awk))
  status <- system2("/bin/bash", c("-o", "pipefail", "-c", shQuote(cmd)),
                    stdout = target, stderr = paste0(target, ".stderr.log"))
  stopifnot(status == 0L)
  x <- fread(target)
  stopifnot(nrow(x) == 22L, all(x$headers_read == length(members)), all(x$n_invalid_P == 0L),
            all(x$n_zero_P == 0L), sum(x$n_valid) > 1000000L)
  x[, chromosome := as.character(Chr)]
  x[chromosome == "21", chromosome := "X"]
  x[chromosome == "24", chromosome := "MT"]
  stopifnot(setequal(x$chromosome, c(as.character(1:20), "X", "MT")))
  stopifnot(all(x[n_valid > 0, chromosome == sub(":.*$", "", lead_SNP)]))
  x[, nuclear := chromosome != "MT"]
  nuclear <- x[nuclear == TRUE & n_valid > 0][order(min_P)]
  best <- x[n_valid > 0][order(min_P)][1]
  s <- data.table(trait = t, files_scanned = length(members), SNP_rows = sum(x$n_rows),
    valid_P_rows = sum(x$n_valid), missing_P_rows = sum(x$n_missing_P),
    nuclear_valid_P_rows = sum(x[nuclear == TRUE, n_valid]),
    lead_SNP_all = best$lead_SNP, min_P_all = best$min_P, peak_log10P_all = -log10(best$min_P),
    lead_SNP_nuclear = nuclear$lead_SNP[1], min_P_nuclear = nuclear$min_P[1],
    peak_log10P_nuclear = -log10(nuclear$min_P[1]),
    nuclear_SNPs_P_le_1e_minus7 = sum(x[nuclear == TRUE, n_P_le_1e_minus7]),
    passes_all_chromosomes_7 = best$min_P <= 1e-7,
    passes_nuclear_7 = nuclear$min_P[1] <= 1e-7,
    lead_Wald_log10P_error = abs(-pchisq((nuclear$lead_beta[1] / nuclear$lead_se[1])^2,
      df = 1, lower.tail = FALSE, log.p = TRUE) / log(10) + log10(nuclear$min_P[1])))
  stopifnot(s$lead_Wald_log10P_error < .001)
  screens[[i]] <- s
  fwrite(rbindlist(screens), file.path(out, "screen_progress.tsv"), sep = "\t")
  message("  peak=", round(s$peak_log10P_nuclear, 4), "; lead=", s$lead_SNP_nuclear)
}
result <- merge(rbindlist(screens), meta, by = "trait")
setorder(result, min_P_nuclear)
result[, screen_rank := seq_len(.N)]
result[, selection_status := fcase(!passes_nuclear_7, "below_requested_threshold",
  !N_verified, "threshold_pass_but_N_unverified", !dictionary_documented, "threshold_pass_intermediate_trait_needs_definition",
  default = "threshold_pass_ready_for_input_QC")]
fwrite(result, file.path(out, "all_29_trait_peak_screen.tsv"), sep = "\t")
fwrite(result[passes_nuclear_7 == TRUE], file.path(out, "threshold7_candidates.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(out, "R_sessionInfo.txt"))
print(result[, .(trait, N_phenotype, peak_log10P_nuclear, lead_SNP_nuclear, selection_status)])
message("Screen finished: ", out)
