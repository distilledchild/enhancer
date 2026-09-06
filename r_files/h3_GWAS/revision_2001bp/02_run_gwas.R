#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
repo <- Sys.getenv("ENHANCER_REPO", "/Users/pete/Desktop/playground/enhancer")
external <- Sys.getenv("ENHANCER_EXTERNAL", "/Volumes/external_1000GB_all/playground/enhancer")
root <- file.path(external, "r_files/h3_GWAS")
scope <- Sys.getenv("GWAS_SCOPE", "Lara_King")
stopifnot(scope %in% c("Lara_King", "nicotine"))
out <- file.path(root, if (scope == "nicotine") "revision_2001bp_nicotine_loop_only" else "revision_2001bp_loop_only")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
ap <- file.path(root, "revision_2001bp_loop_only/annotation")
stopifnot(file.exists(file.path(ap, "annotation_complete.rds")))
for (d in c("runs", "results")) dir.create(file.path(out, d), showWarnings = FALSE)
workers <- as.integer(Sys.getenv("GWAS_WORKERS", "3"))
stopifnot(workers >= 1L, workers <= 4L)
mode <- commandArgs(trailingOnly = TRUE)
if (!length(mode)) mode <- "run"
stopifnot(mode %in% c("run", "summarize"))
magma <- file.path(external, "tools/magma")
bfile <- file.path(external, "data/hs_data/v4/HS_genotypes_v4")
annot <- file.path(ap, "hs.revised2001bp.loop_only.hmagma.annot")

spec <- rbindlist(list(
  data.table(cohort = "Lara2024_delaydiscounting", trait = paste0("regressedlr_", c("auc", "exponential_k", "dd_indiff_2")),
             result_dir = file.path(root, "Lara2024_delaydiscounting/results")),
  data.table(cohort = "King2025", trait = c("crf_ny_incentive_value_index", "crf_ny_lever_presses", "pavca_ny_d5_response_bias", "pavca_ny_d5_index"),
             result_dir = file.path(root, "King2025"))
))
if (scope == "nicotine") {
  spec <- data.table(cohort = "Nicotine", trait = "nicsa_day5_infusion",
                     result_dir = file.path(repo, "r_files/h3_GWAS/nicotine"))
}
# Read the original commands, rather than silently choosing new inputs or N.
manifest <- rbindlist(lapply(seq_len(nrow(spec)), function(i) {
  r <- spec[i]
  prefix <- file.path(r$result_dir, paste0(r$trait, "_hmagma_strict_duttke_telese"))
  if (r$cohort == "Nicotine") {
    prefix <- file.path(r$result_dir, "nicsa_day5_hmagma_strict_duttke_telese_20260617")
  }
  log <- readLines(paste0(prefix, ".log"), warn = FALSE)
  p_line <- trimws(grep("^\\s*--pval ", log, value = TRUE))
  b_line <- trimws(grep("^\\s*--bfile ", log, value = TRUE))
  a_line <- trimws(grep("^\\s*--gene-annot ", log, value = TRUE))
  stopifnot(length(p_line) == 1L, length(b_line) == 1L, length(a_line) == 1L,
            identical(sub("^--bfile ", "", b_line), bfile),
            grepl("use=SNP,p N=[0-9]+$", p_line),
            any(grepl("MAGMA v1.08", log, fixed = TRUE)),
            any(grepl("^End time is ", log)))
  source <- sub(" use=SNP,p N=[0-9]+$", "", sub("^--pval ", "", p_line))
  source_logged <- source
  if (r$cohort == "Nicotine" && !file.exists(source)) {
    source <- "/Users/pete/Library/CloudStorage/Dropbox/Gateway_to_Hao/enhancer/sep_HMAGMA/infusion_combined1/regressedlr_nicsa_day5_infusion_chrgwas_combined.mlma"
  }
  data.table(cohort = r$cohort, trait = r$trait, source, source_logged,
    N = as.integer(sub("^.* N=", "", p_line)),
    legacy_annotation = sub("^--gene-annot ", "", a_line),
    legacy_result = paste0(prefix, ".genes.out"), legacy_log = paste0(prefix, ".log"),
    output_prefix = file.path(out, "runs", paste0(r$cohort, "__", r$trait)))
}))
stopifnot(all(file.exists(manifest$source)), all(file.exists(manifest$legacy_result)),
          all(file.exists(manifest$legacy_annotation)), nrow(manifest) == if (scope == "nicotine") 1L else 7L)
expected_old_md5 <- unname(tools::md5sum(file.path(repo, "r_files/h3_GWAS/nicotine/hs.hic.duttke_telese_atac.non_tss_promoter_filtered.annot")))
stopifnot(all(unname(tools::md5sum(manifest$legacy_annotation)) == expected_old_md5))
manifest[, source_md5 := unname(tools::md5sum(source))]
manifest[, new_annotation_md5 := unname(tools::md5sum(annot))]
manifest[, bfile_path := bfile]
manifest[, loop_window_bp := 2001L]
manifest[, source_kind := if (scope == "nicotine") "existing_local_GWAS_statistics" else "deposited_GWAS_statistics"]
manifest[, note := "Same pre-existing GWAS file, N, LD reference, MAGMA version/model and SNP rules retained; loop input only replaced"]
fwrite(manifest, file.path(out, "dataset_manifest.tsv"), sep = "\t")
reference <- data.table(path = paste0(bfile, c(".bed", ".bim", ".fam")))
reference[, size := file.info(path)$size]
reference[, mtime := as.numeric(file.info(path)$mtime)]
fwrite(reference, file.path(out, "LD_reference_metadata.tsv"), sep = "\t")

# data.table's process-local self-reference is not part of input identity.
fingerprints_match <- function(saved, current) {
  saved$reference <- as.data.frame(saved$reference)
  current$reference <- as.data.frame(current$reference)
  identical(saved, current)
}

if (mode == "run") {
  evidence <- fread(file.path(ap, "revised_HMAGMA_SNP_loop_gene_evidence.tsv.gz"))
  run_one <- function(i) {
    row <- manifest[i]
    receipt_file <- paste0(row$output_prefix, ".complete.rds")
    fp <- list(source_md5 = row$source_md5, annotation_md5 = row$new_annotation_md5,
               N = row$N, reference = as.data.frame(reference), magma_md5 = unname(tools::md5sum(magma)))
    if (file.exists(receipt_file)) {
      receipt <- readRDS(receipt_file)
      stopifnot(fingerprints_match(receipt$fingerprint, fp),
                identical(receipt$output_md5, unname(tools::md5sum(paste0(row$output_prefix, ".genes.out")))))
      return(receipt$status)
    }
    message(format(Sys.time(), "%H:%M:%S"), " START ", row$cohort, " / ", row$trait)
    started <- Sys.time()
    p <- fread(file = row$source, select = c("SNP", "p"))
    n_raw <- nrow(p)
    repeated_headers <- p$SNP == "SNP" & as.character(p$p) == "p"
    if (any(repeated_headers)) {
      stopifnot(row$cohort == "Nicotine", sum(repeated_headers) == 20L)
      p <- p[!repeated_headers]
    }
    p[, p := as.numeric(p)]
    stopifnot(!anyDuplicated(p$SNP), !anyNA(p$SNP), all(is.finite(p$p) & p$p > 0 & p$p <= 1))
    ep <- merge(evidence, p[SNP %chin% evidence$SNP], by = "SNP", sort = FALSE)
    fwrite(ep, file.path(out, "results", paste0(row$cohort, "__", row$trait, ".SNP_loop_gene_evidence.tsv.gz")), sep = "\t")
    qc <- data.table(cohort = row$cohort, trait = row$trait, n_GWAS_SNPs = nrow(p),
                     n_input_rows = n_raw, repeated_header_rows_ignored_for_evidence = sum(repeated_headers),
                     n_HiC_SNPs_with_p = uniqueN(ep$SNP), n_loops_with_HiC_SNPs = uniqueN(ep$loop_id),
                     n_genes_with_HiC_SNPs = uniqueN(ep$gene_id))
    fwrite(qc, paste0(row$output_prefix, ".input_QC.tsv"), sep = "\t")
    rm(p, ep)
    gc()
    args <- c("--bfile", shQuote(bfile), "--pval", shQuote(row$source), "use=SNP,p", paste0("N=", row$N),
              "--gene-annot", shQuote(annot), "--out", shQuote(row$output_prefix))
    status <- system2(magma, args, stdout = paste0(row$output_prefix, ".stdout.log"), stderr = paste0(row$output_prefix, ".stderr.log"))
    log <- readLines(paste0(row$output_prefix, ".log"), warn = FALSE)
    if (status != 0L || !any(grepl("^End time is ", log))) stop("MAGMA failed: ", row$output_prefix)
    result <- fread(paste0(row$output_prefix, ".genes.out"))
    stopifnot(nrow(result) > 1000L, !anyDuplicated(result$GENE), all(is.finite(result$P) & result$P > 0 & result$P <= 1))
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    status_row <- data.table(cohort = row$cohort, trait = row$trait, status = "completed", elapsed_seconds = elapsed)
    saveRDS(list(fingerprint = fp, output_md5 = unname(tools::md5sum(paste0(row$output_prefix, ".genes.out"))),
                 command = c(magma, args), status = status_row), receipt_file)
    message(format(Sys.time(), "%H:%M:%S"), " DONE ", row$cohort, " / ", row$trait, " (", round(elapsed), " seconds)")
    status_row
  }
  status <- parallel::mclapply(seq_len(nrow(manifest)), run_one, mc.cores = workers, mc.preschedule = FALSE)
  bad <- vapply(status, inherits, logical(1), "try-error")
  if (any(bad)) stop("Failed runs: ", paste(which(bad), collapse = ", "))
  fwrite(rbindlist(status), file.path(out, "run_status.tsv"), sep = "\t")
}

message("Comparing revised vs legacy H-MAGMA; no cMAGMA rerun.")
changes <- fread(file.path(ap, "gene_annotation_changes.tsv"))
read_pairs <- function(path) {
  parts <- strsplit(readLines(path), "\t", fixed = TRUE)
  unique(rbindlist(lapply(parts, function(x) data.table(gene_id = x[1], SNP = x[-c(1, 2)]))))
}
old_pairs <- read_pairs(manifest$legacy_annotation[1])
new_pairs <- read_pairs(annot)
new_unique_snps <- setdiff(new_pairs$SNP, old_pairs$SNP)
lost_unique_snps <- setdiff(old_pairs$SNP, new_pairs$SNP)
gained_pairs <- fsetdiff(new_pairs, old_pairs)
gained_pairs[, SNP_previously_unassigned := SNP %chin% new_unique_snps]
fwrite(data.table(SNP = new_unique_snps), file.path(ap, "newly_included_SNPs.tsv"), sep = "\t")
fwrite(data.table(SNP = lost_unique_snps), file.path(ap, "lost_unique_SNPs.tsv"), sep = "\t")
fwrite(gained_pairs, file.path(ap, "gained_SNP_gene_pairs.tsv"), sep = "\t")
novel_qc <- vector("list", nrow(manifest))
verification_qc <- vector("list", nrow(manifest))
tables <- vector("list", nrow(manifest))
verified_status <- vector("list", nrow(manifest))
magma_warning_qc <- vector("list", nrow(manifest))
for (i in seq_len(nrow(manifest))) {
  row <- manifest[i]
  stopifnot(file.exists(paste0(row$output_prefix, ".complete.rds")))
  receipt <- readRDS(paste0(row$output_prefix, ".complete.rds"))
  current_fp <- list(source_md5 = row$source_md5, annotation_md5 = row$new_annotation_md5,
                     N = row$N, reference = as.data.frame(reference), magma_md5 = unname(tools::md5sum(magma)))
  stopifnot(fingerprints_match(receipt$fingerprint, current_fp),
            identical(receipt$output_md5, unname(tools::md5sum(paste0(row$output_prefix, ".genes.out")))))
  verified_status[[i]] <- receipt$status
  old_warnings <- trimws(grep("WARNING:", readLines(row$legacy_log, warn = FALSE), value = TRUE))
  new_log <- readLines(paste0(row$output_prefix, ".log"), warn = FALSE)
  stopifnot(any(grepl("^End time is ", new_log)))
  new_warnings <- trimws(grep("WARNING:", new_log, value = TRUE))
  magma_warning_qc[[i]] <- data.table(cohort = row$cohort, trait = row$trait,
    n_legacy_warnings = length(old_warnings), n_revised_warnings = length(new_warnings),
    identical_warning_sets = setequal(old_warnings, new_warnings),
    revised_warnings = paste(new_warnings, collapse = " | "))
  new <- fread(paste0(row$output_prefix, ".genes.out"))
  old <- fread(row$legacy_result)
  new[, BH_q := p.adjust(P, "BH")]
  new[, Bonferroni_p := p.adjust(P, "bonferroni")]
  old[, BH_q := p.adjust(P, "BH")]
  old[, Bonferroni_p := p.adjust(P, "bonferroni")]
  setnames(new, c("GENE", "P", "NSNPS", "BH_q", "Bonferroni_p"),
                c("gene_id", "revised_P", "revised_NSNPS", "revised_BH_q", "revised_Bonferroni_p"))
  setnames(old, c("GENE", "P", "NSNPS", "BH_q", "Bonferroni_p"),
                c("gene_id", "legacy_P", "legacy_NSNPS", "legacy_BH_q", "legacy_Bonferroni_p"))
  x <- merge(new, old[, .(gene_id, legacy_P, legacy_NSNPS, legacy_BH_q, legacy_Bonferroni_p)], by = "gene_id", all = TRUE)
  x <- merge(x, changes, by = "gene_id", all.x = TRUE)
  x[, cohort := row$cohort]
  x[, trait := row$trait]
  x[, source_kind := row$source_kind]
  x[, annotation_changed := n_gained_SNPs > 0L | n_lost_SNPs > 0L]
  ep <- fread(file.path(out, "results", paste0(row$cohort, "__", row$trait, ".SNP_loop_gene_evidence.tsv.gz")))
  support <- ep[, .(n_trait_HiC_SNPs = uniqueN(SNP), n_trait_HiC_loops = uniqueN(loop_id)), by = gene_id]
  x <- merge(x, support, by = "gene_id", all.x = TRUE)
  for (col in c("n_trait_HiC_SNPs", "n_trait_HiC_loops")) set(x, which(is.na(x[[col]])), col, 0L)
  novel <- merge(ep, gained_pairs, by = c("gene_id", "SNP"))
  novel_gene <- unique(novel[, .(gene_id, SNP, p, SNP_previously_unassigned)])[order(p),
    .(n_trait_new_SNP_gene_pairs = .N,
      n_trait_previously_unassigned_SNPs = sum(SNP_previously_unassigned),
      best_new_SNP = SNP[1], best_new_SNP_P = p[1]), by = gene_id]
  x <- merge(x, novel_gene, by = "gene_id", all.x = TRUE)
  for (col in c("n_trait_new_SNP_gene_pairs", "n_trait_previously_unassigned_SNPs")) {
    set(x, which(is.na(x[[col]])), col, 0L)
  }
  fwrite(novel[order(p)], file.path(out, "results",
    paste0(row$cohort, "__", row$trait, ".new_SNP_gene_assignments.tsv.gz")), sep = "\t")
  novel_qc[[i]] <- data.table(cohort = row$cohort, trait = row$trait,
    n_previously_unassigned_SNPs_with_GWAS_p = uniqueN(novel[SNP_previously_unassigned == TRUE, SNP]),
    n_new_SNP_gene_pairs_with_GWAS_p = nrow(unique(novel[, .(gene_id, SNP)])),
    n_previously_annotated_SNPs_with_new_gene_links = uniqueN(novel[SNP_previously_unassigned == FALSE, SNP]),
    n_genes_with_new_SNP_links = uniqueN(novel$gene_id))
  x[, result_status := fifelse(is.na(revised_P), "legacy_only_tested", fifelse(is.na(legacy_P), "revised_only_tested", "both_tested"))]
  # Same reference-level SNP sets should preserve the same marginal statistic.
  unchanged <- x[!annotation_changed & !is.na(revised_P) & !is.na(legacy_P)]
  stopifnot(all(abs(unchanged$revised_P - unchanged$legacy_P) <= 1e-5),
            all(unchanged$revised_NSNPS == unchanged$legacy_NSNPS))
  verification_qc[[i]] <- data.table(cohort = row$cohort, trait = row$trait,
    n_unchanged_annotation_genes_compared = nrow(unchanged),
    maximum_absolute_P_difference = max(abs(unchanged$revised_P - unchanged$legacy_P)),
    n_unchanged_genes_with_SNP_count_difference = sum(unchanged$revised_NSNPS != unchanged$legacy_NSNPS))
  tables[[i]] <- x
}
all <- rbindlist(tables, fill = TRUE)
all[, revised_global_BH_q := NA_real_]
all[, legacy_global_BH_q := NA_real_]
all[!is.na(revised_P), revised_global_BH_q := p.adjust(revised_P, "BH")]
all[!is.na(legacy_P), legacy_global_BH_q := p.adjust(legacy_P, "BH")]
all[, newly_BH05_basis := NA_character_]
all[revised_BH_q < 0.05 & (is.na(legacy_BH_q) | legacy_BH_q >= 0.05),
    newly_BH05_basis := fcase(is.na(legacy_P), "newly_testable_gene",
      revised_P == legacy_P, "BH_threshold_only_unchanged_gene_P",
      default = "changed_gene_P")]
setorder(all, cohort, trait, revised_P, na.last = TRUE)
fwrite(all, file.path(out, "results/all_gene_results_revised_vs_legacy.tsv.gz"), sep = "\t")
for (i in seq_len(nrow(manifest))) {
  row <- manifest[i]
  fwrite(all[cohort == row$cohort & trait == row$trait],
         file.path(out, "results", paste0(row$cohort, "__", row$trait, ".gene_results.tsv")), sep = "\t")
}
summary <- all[order(revised_P, na.last = TRUE), .(
  n_legacy_tested = sum(!is.na(legacy_P)), n_revised_tested = sum(!is.na(revised_P)),
  legacy_BH05 = sum(legacy_BH_q < 0.05, na.rm = TRUE), revised_BH05 = sum(revised_BH_q < 0.05, na.rm = TRUE),
  newly_BH05 = sum(revised_BH_q < 0.05 & (is.na(legacy_BH_q) | legacy_BH_q >= 0.05), na.rm = TRUE),
  lost_BH05 = sum(legacy_BH_q < 0.05 & (is.na(revised_BH_q) | revised_BH_q >= 0.05), na.rm = TRUE),
  legacy_Bonf05 = sum(legacy_Bonferroni_p < 0.05, na.rm = TRUE), revised_Bonf05 = sum(revised_Bonferroni_p < 0.05, na.rm = TRUE),
  revised_global_BH05 = sum(revised_global_BH_q < 0.05, na.rm = TRUE),
  top_gene = gene_id[1], top_symbol = gene_symbol[1], top_P = revised_P[1], top_q = revised_BH_q[1]),
  by = .(cohort, trait, source_kind)]
setorder(summary, cohort, trait)
fwrite(summary, file.path(out, "results/trait_summary.tsv"), sep = "\t")
fwrite(all[revised_BH_q < 0.05], file.path(out, "results/revised_genes_BH05_per_trait.tsv"), sep = "\t")
fwrite(all[revised_BH_q < 0.05 & (is.na(legacy_BH_q) | legacy_BH_q >= 0.05)],
       file.path(out, "results/newly_BH05_genes_vs_legacy.tsv"), sep = "\t")
fwrite(all[revised_BH_q < 0.05 & annotation_changed], file.path(out, "results/revised_BH05_genes_with_changed_annotation.tsv"), sep = "\t")
fwrite(all[revised_BH_q < 0.05 & n_trait_HiC_SNPs > 0L],
       file.path(out, "results/revised_BH05_genes_with_actual_HiC_SNP_support.tsv"), sep = "\t")
if (nrow(manifest) > 1L) {
  fwrite(all[revised_global_BH_q < 0.05], file.path(out, "results", paste0("revised_genes_BH05_across_", nrow(manifest), "_traits.tsv")), sep = "\t")
}
fwrite(all[!is.na(revised_P)][order(revised_P), head(.SD, 20L), by = .(cohort, trait)],
       file.path(out, "results/revised_top20_per_trait.tsv"), sep = "\t")
qcs <- lapply(manifest$output_prefix, function(p) fread(paste0(p, ".input_QC.tsv")))
fwrite(rbindlist(qcs), file.path(out, "input_QC.tsv"), sep = "\t")
fwrite(rbindlist(novel_qc), file.path(out, "results/new_SNP_assignment_summary.tsv"), sep = "\t")
fwrite(rbindlist(verification_qc), file.path(out, "results/unchanged_annotation_verification.tsv"), sep = "\t")
fwrite(rbindlist(magma_warning_qc), file.path(out, "results/MAGMA_warning_QC.tsv"), sep = "\t")
fwrite(rbindlist(verified_status), file.path(out, "run_status.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "R_session_info.txt"))
print(summary)
message("Complete: ", out)
