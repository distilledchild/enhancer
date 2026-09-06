#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(xml2)
  library(jsonlite)
})
setDTthreads(1)
repo <- Sys.getenv("ENHANCER_REPO", "/Users/pete/Desktop/playground/enhancer")
external <- Sys.getenv("ENHANCER_EXTERNAL", "/Volumes/external_1000GB_all/playground/enhancer")
root <- file.path(external, "r_files/h3_GWAS")
shared <- file.path(root, "revision_2001bp_loop_only")
nicotine <- file.path(root, "revision_2001bp_nicotine_loop_only")
out <- file.path(nicotine, "input_audit")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
report_dir <- "/Users/pete/Library/CloudStorage/Dropbox/Gateway_to_Hao/enhancer/data/hs_data"
reports <- list.files(report_dir, "^gwas_report.*html$", full.names = TRUE)
phenotype_url <- "https://palmerlab.s3.sdsc.edu/tsanches_dash_genotypes/gwas_results/p50_hao_chen_nicsa/processed_data_ready.csv"
phenotype_file <- file.path(out, "nicotine_source_phenotypes.csv")
if (!file.exists(phenotype_file)) download.file(phenotype_url, phenotype_file, mode = "wb", method = "libcurl")
phenotypes <- fread(file = phenotype_file)
fwrite(data.table(url = phenotype_url, local_file = phenotype_file,
  md5 = unname(tools::md5sum(phenotype_file)), cohort_rows = nrow(phenotypes),
  trait = "regressedlr_nicsa_day5_infusion", N_nonmissing = sum(!is.na(phenotypes$regressedlr_nicsa_day5_infusion))),
  file.path(out, "nicotine_phenotype_N_verification.tsv"), sep = "\t")

# Recover tables from the two locally archived report formats, without executing HTML.
decode_array <- function(v) {
  if (!is.list(v) || !identical(v$type, "ndarray")) return(unlist(v))
  if (v$dtype == "object") return(unlist(v$array))
  stopifnot(identical(v$array$type, "bytes"), v$dtype %in% c("float64", "int32"))
  readBin(base64_dec(v$array$data), if (v$dtype == "float64") "double" else "integer",
          n = prod(unlist(v$shape)), size = if (v$dtype == "float64") 8 else 4,
          endian = v$order)
}
extract_tables <- function(path) {
  doc <- read_html(path)
  scripts <- xml_find_all(doc, "//script[@type='application/json']")
  found <- list()
  context <- character()
  add <- function(d) {
    if (!"trait" %in% names(d)) return()
    d <- d[sub("^regressedlr_", "", trait) == "nicsa_day5_infusion"]
    if (nrow(d)) found[[length(found) + 1L]] <<- d
  }
  walk <- function(v) {
    if (!is.list(v)) return()
    if (identical(v$name, "panel.models.markup.HTML")) {
      t <- v$attributes$text
      if (!is.null(t) && nchar(t) < 40000L && grepl("general-information|summary-of-qtls|snp-heritability-estimates", t)) {
        decoded <- t
        for (pass in 1:3) decoded <- xml_text(read_html(paste0("<div>", decoded, "</div>")))
        context <<- c(context, decoded)
      }
      return()
    }
    if (identical(v$name, "ColumnDataSource")) {
      entries <- v$attributes$data$entries
      keys <- vapply(entries, function(e) e[[1]], "")
      if ("trait" %in% keys && any(c("n", "TopSNP") %in% keys)) {
        d <- as.data.table(setNames(lapply(entries, function(e) decode_array(e[[2]])), keys))
        add(d)
      }
      return()
    }
    for (a in v) if (is.list(a)) walk(a)
  }
  for (s in scripts) {
    j <- fromJSON(xml_text(s), simplifyVector = FALSE)
    if (!is.null(j$x$container) && !is.null(j$x$data)) {
      keys <- xml_text(xml_find_all(read_html(j$x$container), "//th"))
      if ("trait" %in% keys && any(c("n", "TopSNP", "SNP") %in% keys)) {
        d <- as.data.table(lapply(j$x$data, unlist))
        stopifnot(ncol(d) == length(keys))
        setnames(d, keys)
        add(d)
      }
    } else walk(j)
  }
  attr(found, "report_context") <- unique(context)
  found
}
ns <- list()
leads <- list()
report_context <- character()
for (f in reports) {
  message("Reading archived report: ", basename(f))
  tables <- extract_tables(f)
  if (length(attr(tables, "report_context"))) {
    report_context <- c(report_context, paste("SOURCE:", f), attr(tables, "report_context"), "")
  }
  for (d in tables) {
    if ("n" %in% names(d)) {
      ns[[length(ns) + 1L]] <- d[, .(report = f, trait, N = n, heritability, heritability_se)]
    }
    if (any(c("TopSNP", "SNP") %in% names(d)) && "significance_level" %in% names(d)) {
      snp_col <- if ("TopSNP" %in% names(d)) "TopSNP" else "SNP"
      se_col <- if ("betase" %in% names(d)) "betase" else "se"
      af_col <- if ("Freq" %in% names(d)) "Freq" else "af"
      leads[[length(leads) + 1L]] <- d[, .(report = f, trait, SNP = get(snp_col),
        report_beta = beta, report_se = get(se_col), report_af = get(af_col),
        report_logp = get("-Log10(p)"), significance_level)]
    }
  }
}
ns <- unique(rbindlist(ns))
leads <- unique(rbindlist(leads))
fwrite(ns, file.path(out, "nicotine_report_sample_sizes.tsv"), sep = "\t")
writeLines(report_context, file.path(out, "nicotine_report_metadata.txt"))
fwrite(data.table(path = reports, md5 = unname(tools::md5sum(reports))),
       file.path(out, "report_checksums.tsv"), sep = "\t")
print(ns[, .(report = basename(report), N, heritability, heritability_se)])

man <- rbindlist(lapply(c(shared, nicotine), function(d) {
  x <- fread(file.path(d, "dataset_manifest.tsv"))
  x[, analysis_root := d]
  x
}), fill = TRUE)
bfile <- file.path(external, "data/hs_data/v4/HS_genotypes_v4")
bim <- fread(paste0(bfile, ".bim"), col.names = c("chr", "SNP", "cm", "bp", "a1", "a2"))
stopifnot(!anyDuplicated(bim$SNP))
read_pairs <- function(path) {
  x <- strsplit(readLines(path), "\t", fixed = TRUE)
  unique(rbindlist(lapply(x, function(z) data.table(gene_id = z[1], SNP = z[-c(1, 2)]))))
}
new_pairs <- read_pairs(file.path(shared, "annotation/hs.revised2001bp.loop_only.hmagma.annot"))
old_pairs <- read_pairs(man$legacy_annotation[1])
direct <- read_pairs(file.path(repo, "r_files/h3_GWAS/nicotine/hs.exonpro.ONLY.annot"))
summaries <- list()
correction_checks <- list()
for (i in seq_len(nrow(man))) {
  r <- man[i]
  message("Checking all SNP P-values: ", r$cohort, " / ", r$trait)
  p <- fread(file = r$source)
  headers <- p$SNP == "SNP" & as.character(p$p) == "p"
  n_headers <- sum(headers)
  p <- p[!headers]
  p[, p := as.numeric(p)]
  stopifnot(!anyDuplicated(p$SNP), all(is.finite(p$p) & p$p > 0 & p$p <= 1),
            identical(unname(tools::md5sum(r$source)), r$source_md5))
  idx <- match(p$SNP, bim$SNP)
  g <- fread(paste0(r$output_prefix, ".genes.out"))
  cmp <- fread(file.path(r$analysis_root, "results", paste0(r$cohort, "__", r$trait, ".gene_results.tsv")))
  cmp <- cmp[match(g$GENE, gene_id)]
  stopifnot(identical(g$GENE, cmp$gene_id))
  # Independent implementation, not p.adjust, checks every reported BH value.
  ord <- order(g$P)
  bh <- numeric(nrow(g))
  bh[ord] <- pmin(1, rev(cummin(rev(g$P[ord] * nrow(g) / seq_len(nrow(g))))))
  correction_checks[[i]] <- data.table(cohort = r$cohort, trait = r$trait,
    max_BH_error = max(abs(bh - cmp$revised_BH_q)),
    max_Bonferroni_error = max(abs(pmin(1, g$P * nrow(g)) - cmp$revised_Bonferroni_p)))
  summaries[[i]] <- data.table(cohort = r$cohort, trait = r$trait, N_configured = r$N,
    SNPs = nrow(p), duplicate_header_rows = n_headers, LD_matched_SNPs = sum(!is.na(idx)),
    min_SNP_P = min(p$p), SNP_P_lt_5e_minus8 = sum(p$p < 5e-8), SNP_P_lt_1e_minus5 = sum(p$p < 1e-5),
    SNP_P_lt_0_01 = sum(p$p < .01), median_SNP_P = median(p$p),
    lambda_GC = median(qchisq(p$p, 1, lower.tail = FALSE)) / qchisq(.5, 1),
    genes_tested = nrow(g), min_gene_P = min(g$P), min_BH_q = min(bh),
    BH05_genes = sum(bh < .05), Bonf05_genes = sum(g$P < .05 / nrow(g)))
  if (r$cohort != "Nicotine") next
  for (col in c("b", "se", "Freq", "bp")) set(p, j = col, value = as.numeric(p[[col]]))
  stopifnot(all(is.finite(p$b)), all(is.finite(p$se) & p$se > 0),
            all(is.finite(p$Freq) & p$Freq >= 0 & p$Freq <= 1), !anyNA(idx))
  p[, wald_logp := -pchisq((b / se)^2, 1, lower.tail = FALSE, log.p = TRUE) / log(10)]
  delta <- abs(p$wald_logp + log10(p$p))
  expected_chr <- as.character(p$Chr)
  expected_chr[expected_chr == "21"] <- "X"
  ref_chr <- as.character(bim$chr[idx])
  ref_chr[ref_chr == "23"] <- "X"
  allele_bad <- !((p$A1 == bim$a1[idx] & p$A2 == bim$a2[idx]) |
                  (p$A1 == bim$a2[idx] & p$A2 == bim$a1[idx]))
  fwrite(data.table(SNPs = nrow(p), max_abs_Wald_logp_error = max(delta),
    n_logp_errors_above_001 = sum(delta > .001), median_abs_Wald_logp_error = median(delta),
    n_position_mismatches = sum(p$bp != bim$bp[idx]), n_chr_mismatches = sum(expected_chr != ref_chr),
    n_allele_set_mismatches = sum(allele_bad), min_MAF = min(pmin(p$Freq, 1 - p$Freq)),
    SNPs_below_filename_threshold_5_3591 = sum(p$p < 10^-5.3591),
    SNPs_below_report_5pct_threshold_5_58 = sum(p$p < 10^-5.58)),
    file.path(out, "nicotine_raw_QC.tsv"), sep = "\t")
  lead_check <- merge(leads, p[, .(SNP, raw_beta = b, raw_se = se, raw_af = Freq, raw_P = p)],
                      by = "SNP", all.x = TRUE)
  lead_check[, raw_logp := -log10(raw_P)]
  lead_check[, abs_logp_difference := abs(raw_logp - report_logp)]
  fwrite(lead_check, file.path(out, "nicotine_report_lead_SNP_comparison.tsv"), sep = "\t")
  strong <- p[p < 1e-5][order(p)]
  strong[, direct_genes := vapply(SNP, function(s) paste(direct[SNP == s, gene_id], collapse = ";"), "")]
  strong[, revised_genes := vapply(SNP, function(s) paste(new_pairs[SNP == s, gene_id], collapse = ";"), "")]
  strong[, legacy_genes := vapply(SNP, function(s) paste(old_pairs[SNP == s, gene_id], collapse = ";"), "")]
  fwrite(strong, file.path(out, "nicotine_strong_SNP_annotation.tsv"), sep = "\t")
  gs <- merge(new_pairs, p[, .(SNP, p)], by = "SNP")
  detail <- gs[, .(SNPs_with_P = .N, best_SNP_P = min(p), SNPs_P_lt_001 = sum(p < .01),
                   median_SNP_P = median(p)), by = gene_id]
  detail <- merge(detail, g, by.x = "gene_id", by.y = "GENE")
  detail[, BH_q := bh[match(gene_id, g$GENE)]]
  fwrite(detail[order(P)], file.path(out, "nicotine_gene_SNP_signal_summary.tsv"), sep = "\t")
  # Small diagnostic gene panel; original full-file data remain untouched.
  selected <- unique(c(head(g[order(P), GENE], 10), new_pairs[SNP %chin% strong$SNP, gene_id]))
  annot_lines <- readLines(file.path(shared, "annotation/hs.revised2001bp.loop_only.hmagma.annot"))
  keep <- sub("\t.*", "", annot_lines) %chin% selected
  writeLines(annot_lines[keep], file.path(out, "diagnostic_genes.annot"))
  p2 <- p[SNP %chin% new_pairs[gene_id %chin% selected, SNP], .(SNP, p)]
  fwrite(p2, file.path(out, "diagnostic_clean_SNP_P.tsv"), sep = "\t")
  fwrite(g[GENE %chin% selected], file.path(out, "diagnostic_original_gene_results.tsv"), sep = "\t")
}
fwrite(rbindlist(summaries), file.path(out, "eight_trait_signal_audit.tsv"), sep = "\t")
fwrite(rbindlist(correction_checks), file.path(out, "multiple_testing_verification.tsv"), sep = "\t")
stopifnot(all(rbindlist(correction_checks)$max_BH_error < 1e-12),
          all(rbindlist(correction_checks)$max_Bonferroni_error < 1e-12))

message("Rechecking original Lara/King MLMA conversion and sample-size sources.")
pheno <- fread(file.path(root, "Lara2024_delaydiscounting/data/ucsd_3_1.csv"))
king_data <- file.path(external, "data/King_2025")
king_n <- fread(file.path(king_data, "ucsd_bb5030313v/s3_downloads/results__heritability__heritability.tsv"))
setnames(king_n, 1L, "regressed_trait")
map <- fread(file.path(king_data, "gwas_summary_stats/king_round8_to_rn7_snp_map.tsv"))
stopifnot(!anyDuplicated(map$original_snp))
source_checks <- list()
for (i in which(man$cohort != "Nicotine")) {
  r <- man[i]
  message("Original MLMA: ", r$trait)
  if (r$cohort == "Lara2024_delaydiscounting") {
    paths <- file.path(root, "Lara2024_delaydiscounting/data/gwas_results/gwas_dd_results", paste0(r$trait, ".loco.mlma"))
    expected_n <- sum(!is.na(pheno[[r$trait]]))
    n_source <- file.path(root, "Lara2024_delaydiscounting/data/ucsd_3_1.csv")
  } else {
    paths <- file.path(king_data, "gwas_summary_stats/raw_chrgwas_mlma", r$trait,
      paste0("regressedlr_", r$trait, "_chrgwas", 1:20, ".mlma"))
    expected_n <- king_n[regressed_trait == paste0("regressedlr_", r$trait), n]
    n_source <- file.path(king_data, "ucsd_bb5030313v/s3_downloads/results__heritability__heritability.tsv")
  }
  stopifnot(all(file.exists(paths)), length(expected_n) == 1L, r$N == expected_n)
  raw <- rbindlist(lapply(paths, function(f) fread(file = f, select = c("SNP", "b", "se", "p"))))
  stopifnot(!anyDuplicated(raw$SNP), all(is.finite(raw$se) & raw$se > 0), all(is.finite(raw$b)),
            all(is.finite(raw$p) & raw$p > 0 & raw$p <= 1))
  error <- abs(-pchisq((raw$b / raw$se)^2, 1, lower.tail = FALSE, log.p = TRUE) / log(10) + log10(raw$p))
  transformed <- raw[, .(SNP, p)]
  if (r$cohort == "King2025") {
    transformed[, target := map$rn7_snp[match(SNP, map$original_snp)]]
    transformed <- transformed[!is.na(target), .(p = min(p)), by = .(SNP = target)]
  }
  used <- fread(file = r$source)
  check <- merge(transformed, used, by = "SNP", all = TRUE, suffixes = c("_expected", "_used"))
  stopifnot(!anyNA(check), nrow(check) == nrow(used))
  source_checks[[length(source_checks) + 1L]] <- data.table(cohort = r$cohort, trait = r$trait,
    N_configured = r$N, N_from_source = expected_n, N_source = n_source,
    original_paths = paste(paths, collapse = ";"), original_SNPs = nrow(raw),
    converted_SNPs = nrow(transformed), max_abs_conversion_P_error = max(abs(check$p_expected - check$p_used)),
    max_abs_Wald_logp_error = max(error), Wald_logp_errors_above_001 = sum(error > .001))
  stopifnot(all(check$p_expected == check$p_used), all(error < .001))
}
fwrite(rbindlist(source_checks), file.path(out, "Lara_King_original_source_verification.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(out, "R_sessionInfo.txt"))
message("Audit complete: ", out)
