#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(stringr)
})
setDTthreads(2)
repo <- Sys.getenv("ENHANCER_REPO", "/Users/pete/Desktop/playground/enhancer")
external <- Sys.getenv("ENHANCER_EXTERNAL", "/Volumes/external_1000GB_all/playground/enhancer")
revision <- file.path(repo, "r_files/revision/revision_main")
legacy_dir <- file.path(repo, "r_files/h3_GWAS/nicotine")
out <- file.path(external, "r_files/h3_GWAS/revision_2001bp_loop_only")
ap <- file.path(out, "annotation")
dir.create(ap, recursive = TRUE, showWarnings = FALSE)
inputs <- c(
  loops = file.path(revision, "results/revised_putative_regulatory_loops.tsv"),
  legacy_loops = file.path(repo, "figures/submission/lt2mb/df_final_loop_sub.4.any.lt2mb.ENSEMBL.mid.mid.final.filter.200kb.csv"),
  atac = file.path(revision, "cache_data/gr.atac.rn7.1based.rds"),
  transcripts = file.path(revision, "cache_data/df.transcript.ensembl.rn7.1based.rds"),
  bim = file.path(external, "data/hs_data/v4/HS_genotypes_v4.bim"),
  fam = file.path(external, "data/hs_data/v4/HS_genotypes_v4.fam"),
  baseline_annotation = file.path(legacy_dir, "hs.exonpro.ONLY.annot"),
  legacy_annotation = file.path(legacy_dir, "hs.hic.duttke_telese_atac.non_tss_promoter_filtered.annot"),
  genedef = file.path(legacy_dir, "genedef_ensembl.tsv")
)
stopifnot(all(file.exists(inputs)))
manifest <- data.table(input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(inputs)))
if (file.exists(file.path(ap, "annotation_complete.rds"))) {
  stopifnot(identical(fread(file.path(ap, "input_manifest.tsv")), manifest))
  message("Verified loop-only annotation already complete.")
  quit(status = 0)
}
parse_annot <- function(path) {
  parts <- strsplit(readLines(path), "\t", fixed = TRUE)
  stopifnot(all(lengths(parts) >= 3L))
  unique(rbindlist(lapply(parts, function(x) data.table(gene_id = x[1], location = x[2], SNP = x[-c(1, 2)]))))
}
legacy <- parse_annot(inputs[["legacy_annotation"]])
positional <- parse_annot(inputs[["baseline_annotation"]])
direct_snps <- unique(positional$SNP)
# The original H-MAGMA excludes all direct SNPs from the Hi-C branch.
# This recovers its gene-specific baseline without changing the mapping.
baseline <- legacy[SNP %chin% direct_snps]
stopifnot(setequal(baseline$SNP, direct_snps))
genedef <- fread(inputs[["genedef"]], header = FALSE,
                 col.names = c("location", "chr", "start", "end", "strand", "gene_id"))
genedef[, chr := as.character(chr)]
stopifnot(!anyDuplicated(genedef$gene_id))
tt <- as.data.table(readRDS(inputs[["transcripts"]]))
tt[, chr := sub("^chr", "", chr)]
tt[, gene_id := str_match(attribute, 'gene_id "([^"]+)"')[, 2]]
tt[, gene_symbol := str_match(attribute, 'gene_name "([^"]+)"')[, 2]]
tt[, transcript_biotype := str_match(attribute, 'transcript_biotype "([^"]+)"')[, 2]]
tt <- tt[transcript_biotype == "protein_coding" & chr %chin% c(as.character(1:20), "X")]
stopifnot(setequal(tt$gene_id, genedef$gene_id))
tt[, pro_start := pmax(1L, ifelse(strand == "+", transcript_start - 2000L, transcript_end - 500L))]
tt[, pro_end := ifelse(strand == "+", transcript_start + 500L, transcript_end + 2000L)]
promoters_legacy <- GRanges(tt$chr, IRanges(tt$pro_start, tt$pro_end), strand = tt$strand, gene_id = tt$gene_id)
tss_pos <- ifelse(tt$strand == "+", tt$transcript_start, tt$transcript_end)
tss_points <- GRanges(tt$chr, IRanges(tss_pos, tss_pos), strand = tt$strand)
# Retain the original strand-aware 1000/1000 and 2000/500 exclusions.
exclusion <- reduce(c(promoters(tss_points, upstream = 1000, downstream = 1000), promoters_legacy), ignore.strand = TRUE)
atac <- readRDS(inputs[["atac"]])
seqlevels(atac) <- sub("^chr", "", seqlevels(atac))
distal_atac <- GenomicRanges::setdiff(reduce(atac), exclusion, ignore.strand = TRUE)
bim <- fread(inputs[["bim"]], header = FALSE, select = c(1, 2, 4), col.names = c("chr", "SNP", "position"))
bim[, chr := as.character(chr)]
bim <- bim[chr %chin% c(as.character(1:20), "X")]
stopifnot(!anyDuplicated(bim$SNP))
snps <- GRanges(bim$chr, IRanges(bim$position, bim$position), SNP = bim$SNP)
eligible <- snps[!snps$SNP %chin% direct_snps]
eligible <- eligible[overlapsAny(eligible, distal_atac)]
stopifnot(!any(eligible$SNP %chin% direct_snps))

map_loop_snps <- function(l) {
  orientations <- rbindlist(list(
    l[, .(loop_id, promoter_side = "anchor1", regulatory_side = "anchor2",
           promoter_chr = chr1, promoter_start = start1, promoter_end = end1,
           regulatory_chr = chr2, regulatory_start = start2, regulatory_end = end2)],
    l[, .(loop_id, promoter_side = "anchor2", regulatory_side = "anchor1",
           promoter_chr = chr2, promoter_start = start2, promoter_end = end2,
           regulatory_chr = chr1, regulatory_start = start1, regulatory_end = end1)]
  ))
  promoter_anchors <- GRanges(orientations$promoter_chr, IRanges(orientations$promoter_start, orientations$promoter_end))
  ph <- findOverlaps(promoter_anchors, promoters_legacy)
  linked <- unique(data.table(orientations[queryHits(ph)], gene_id = promoters_legacy$gene_id[subjectHits(ph)]))
  regulatory <- GRanges(linked$regulatory_chr, IRanges(linked$regulatory_start, linked$regulatory_end))
  sh <- findOverlaps(eligible, regulatory)
  evidence <- unique(data.table(linked[subjectHits(sh)], SNP = eligible$SNP[queryHits(sh)],
                                snp_position = start(eligible)[queryHits(sh)]))
  stopifnot(all(evidence$snp_position >= evidence$regulatory_start),
            all(evidence$snp_position <= evidence$regulatory_end))
  pairs <- unique(evidence[, .(gene_id, SNP)])
  pairs <- merge(pairs, genedef[, .(gene_id, location)], by = "gene_id")
  setcolorder(pairs, c("gene_id", "location", "SNP"))
  list(evidence = evidence, pairs = pairs)
}

message("Control: reproducing the existing legacy H-MAGMA SNP-gene annotation.")
old <- fread(inputs[["legacy_loops"]])
parts <- tstrsplit(old$loop.id, "_", fixed = TRUE)
old_loop <- unique(data.table(loop_id = old$loop.id, chr1 = sub("^chr", "", parts[[1]]),
  start1 = as.integer(parts[[2]]) + 1L, end1 = as.integer(parts[[3]]),
  chr2 = sub("^chr", "", parts[[4]]), start2 = as.integer(parts[[5]]) + 1L, end2 = as.integer(parts[[6]])))
old_mapped <- map_loop_snps(old_loop)
reproduced <- unique(rbind(baseline, old_mapped$pairs))
extra <- fsetdiff(reproduced, legacy)
missing <- fsetdiff(legacy, reproduced)
fwrite(extra, file.path(ap, "legacy_control_unexpected_pairs.tsv"), sep = "\t")
fwrite(missing, file.path(ap, "legacy_control_missing_pairs.tsv"), sep = "\t")
stopifnot(nrow(extra) == 0L, nrow(missing) == 0L)
message("Legacy control exact match: ", nrow(legacy), " SNP-gene pairs.")

message("Replacing only the loop input with the revised 2001-bp loop set.")
loops <- fread(inputs[["loops"]])
stopifnot(nrow(loops) == 13376L, !anyDuplicated(loops$loop_id),
          all(loops$revised_putative_regulatory_support),
          all(loops$chr1 == loops$chr2),
          all(loops$end1 - loops$start1 + 1L == loops$resolution_bp),
          all(loops$end2 - loops$start2 + 1L == loops$resolution_bp))
new_loop <- loops[, .(loop_id, chr1 = sub("^chr", "", chr1), start1, end1, chr2 = sub("^chr", "", chr2), start2, end2)]
new_mapped <- map_loop_snps(new_loop)
new <- unique(rbind(baseline, new_mapped$pairs))
stopifnot(nrow(fsetdiff(baseline, new)) == 0L, !any(new_mapped$pairs$SNP %chin% direct_snps))
sets <- new[order(SNP), .(snps = paste(SNP, collapse = "\t"), n_SNPs = .N), by = .(gene_id, location)]
stopifnot(!anyDuplicated(sets$gene_id))
writeLines(paste(sets$gene_id, sets$location, sets$snps, sep = "\t"), file.path(ap, "hs.revised2001bp.loop_only.hmagma.annot"))
fwrite(new_mapped$evidence, file.path(ap, "revised_HMAGMA_SNP_loop_gene_evidence.tsv.gz"), sep = "\t")
fwrite(unique(new_loop), file.path(ap, "revised_loop_input_1based_closed.tsv"), sep = "\t")
symbols <- tt[, .(gene_symbol = paste(sort(unique(na.omit(gene_symbol))), collapse = ";")), by = gene_id]
changes <- merge(legacy[, .(n_legacy_SNPs = .N), by = gene_id], new[, .(n_revised_SNPs = .N), by = gene_id], all = TRUE)
changes <- merge(changes, baseline[, .(n_baseline_SNPs = .N), by = gene_id], all = TRUE)
changes <- merge(changes, symbols, all.x = TRUE)
for (col in c("n_legacy_SNPs", "n_revised_SNPs", "n_baseline_SNPs")) set(changes, which(is.na(changes[[col]])), col, 0L)
changes[, n_new_HiC_assigned_SNPs := n_revised_SNPs - n_baseline_SNPs]
changes[, n_gained_SNPs := 0L]
changes[, n_lost_SNPs := 0L]
changes[fsetdiff(new, legacy)[, .(n = .N), by = gene_id], on = "gene_id", n_gained_SNPs := i.n]
changes[fsetdiff(legacy, new)[, .(n = .N), by = gene_id], on = "gene_id", n_lost_SNPs := i.n]
fwrite(changes, file.path(ap, "gene_annotation_changes.tsv"), sep = "\t")
qc <- data.table(metric = c("legacy_loops", "revised_loops", "legacy_control_extra_pairs", "legacy_control_missing_pairs",
  "baseline_SNP_gene_pairs_unchanged", "legacy_HMAGMA_SNP_gene_pairs", "revised_HMAGMA_SNP_gene_pairs",
  "revised_HiC_SNP_gene_pairs", "revised_HiC_unique_SNPs", "revised_HiC_genes", "revised_loops_with_HMAGMA_SNP_evidence",
  "legacy_genes", "revised_genes"),
  value = c(nrow(old_loop), nrow(new_loop), nrow(extra), nrow(missing), nrow(baseline), nrow(legacy), nrow(new),
            nrow(new_mapped$pairs), uniqueN(new_mapped$pairs$SNP), uniqueN(new_mapped$pairs$gene_id),
            uniqueN(new_mapped$evidence$loop_id), uniqueN(legacy$gene_id), uniqueN(new$gene_id)))
fwrite(qc, file.path(ap, "annotation_QC.tsv"), sep = "\t")
fwrite(manifest, file.path(ap, "input_manifest.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(ap, "R_session_info.txt"))
saveRDS(qc, file.path(ap, "annotation_complete.rds"))
print(qc)
