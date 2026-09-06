#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(GenomicRanges); library(stringr)})
setDTthreads(1)
repo <- Sys.getenv("ENHANCER_REPO", "/Users/pete/Desktop/playground/enhancer")
external <- Sys.getenv("ENHANCER_EXTERNAL", "/Volumes/external_1000GB_all/playground/enhancer")
root <- file.path(external, "r_files/h3_GWAS")
out <- file.path(root, "revision_2001bp_nicotine_loop_only/input_audit")
revision <- file.path(repo, "r_files/revision/revision_main")
strong <- fread(file.path(out, "nicotine_strong_SNP_annotation.tsv"))[p < 10^-5.58]
s <- GRanges(as.character(strong$Chr), IRanges(strong$bp, strong$bp))
tt <- as.data.table(readRDS(file.path(revision, "cache_data/df.transcript.ensembl.rn7.1based.rds")))
tt[, chr := sub("^chr", "", chr)]
tt[, gene_id := str_match(attribute, 'gene_id "([^"]+)"')[, 2]]
tt[, biotype := str_match(attribute, 'transcript_biotype "([^"]+)"')[, 2]]
tt <- tt[biotype == "protein_coding" & chr %chin% c(as.character(1:20), "X")]
tt[, ps := pmax(1, ifelse(strand == "+", transcript_start - 2000, transcript_end - 500))]
tt[, pe := ifelse(strand == "+", transcript_start + 500, transcript_end + 2000)]
pro <- GRanges(tt$chr, IRanges(tt$ps, tt$pe), strand = tt$strand, gene_id = tt$gene_id)
tss <- ifelse(tt$strand == "+", tt$transcript_start, tt$transcript_end)
tss <- GRanges(tt$chr, IRanges(tss, tss), strand = tt$strand)
exc <- reduce(c(promoters(tss, upstream = 1000, downstream = 1000), pro), ignore.strand = TRUE)
atac <- readRDS(file.path(revision, "cache_data/gr.atac.rn7.1based.rds"))
seqlevels(atac) <- sub("^chr", "", seqlevels(atac))
distal <- GenomicRanges::setdiff(reduce(atac), exc, ignore.strand = TRUE)
loops <- fread(file.path(revision, "results/revised_putative_regulatory_loops.tsv"))
a1 <- GRanges(sub("^chr", "", loops$chr1), IRanges(loops$start1, loops$end1))
a2 <- GRanges(sub("^chr", "", loops$chr2), IRanges(loops$start2, loops$end2))
regulatory <- c(a2[overlapsAny(a1, pro)], a1[overlapsAny(a2, pro)])
strong[, overlaps_ATAC := overlapsAny(s, atac)]
strong[, overlaps_distal_ATAC := overlapsAny(s, distal)]
strong[, overlaps_any_revised_anchor := overlapsAny(s, c(a1, a2))]
strong[, overlaps_anchor_with_opposite_promoter := overlapsAny(s, regulatory)]
strong[, directly_assigned := nzchar(direct_genes)]
strong[, predicted_included := directly_assigned | (overlaps_distal_ATAC & overlaps_anchor_with_opposite_promoter)]
strong[, observed_included := nzchar(revised_genes)]
stopifnot(all(strong$predicted_included == strong$observed_included))
fwrite(strong, file.path(out, "nicotine_significant_SNP_filter_trace.tsv"), sep = "\t")
coverage <- strong[, .(SNPs = .N, direct = sum(directly_assigned), any_ATAC = sum(overlaps_ATAC),
  distal_ATAC = sum(overlaps_distal_ATAC), any_loop_anchor = sum(overlaps_any_revised_anchor),
  opposite_promoter_anchor = sum(overlaps_anchor_with_opposite_promoter),
  distal_ATAC_and_opposite_promoter = sum(overlaps_distal_ATAC & overlaps_anchor_with_opposite_promoter),
  final_assigned = sum(observed_included)), by = Chr]
fwrite(coverage, file.path(out, "nicotine_significant_SNP_coverage_summary.tsv"), sep = "\t")
print(coverage)

# Test input formatting and N separately on the same ten diagnostic genes.
# No multiple-testing correction on this selected panel; it is not a new discovery analysis.
magma <- file.path(external, "tools/magma")
bfile <- file.path(external, "data/hs_data/v4/HS_genotypes_v4")
old <- fread(file.path(out, "diagnostic_original_gene_results.tsv"))
stopifnot(nrow(old) == 10L)
checks <- list()
for (n in c(1317L, 1987L)) {
  prefix <- file.path(out, paste0("clean_input_N", n))
  args <- c("--bfile", shQuote(bfile), "--pval", shQuote(file.path(out, "diagnostic_clean_SNP_P.tsv")),
    "use=SNP,p", paste0("N=", n), "--gene-annot", shQuote(file.path(out, "diagnostic_genes.annot")),
    "--out", shQuote(prefix))
  message("Diagnostic MAGMA: N=", n, "; ten genes; clean SNP/P input")
  status <- system2(magma, args, stdout = paste0(prefix, ".stdout.log"), stderr = paste0(prefix, ".stderr.log"))
  stopifnot(status == 0L, any(grepl("^End time is ", readLines(paste0(prefix, ".log")))))
  new <- fread(paste0(prefix, ".genes.out"))
  x <- merge(old[, .(GENE, original_P = P, original_NSNPS = NSNPS)], new, by = "GENE")
  stopifnot(nrow(x) == nrow(old))
  x[, diagnostic_N := n]
  x[, abs_P_difference := abs(P - original_P)]
  checks[[as.character(n)]] <- x
}
checks <- rbindlist(checks)
fwrite(checks, file.path(out, "diagnostic_N_and_formatting_gene_comparison.tsv"), sep = "\t")
print(checks[, .(genes = .N, max_abs_P_difference = max(abs_P_difference),
                SNP_count_differences = sum(NSNPS != original_NSNPS)), by = diagnostic_N])
stopifnot(all(checks[diagnostic_N == 1317, abs_P_difference] == 0),
          all(checks[diagnostic_N == 1317, NSNPS == original_NSNPS]))
message("Diagnostic checks complete; main result files were not overwritten.")
