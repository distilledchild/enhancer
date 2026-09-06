#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
setDTthreads(2)
repo <- "/Users/pete/Desktop/playground/enhancer"
external <- "/Volumes/external_1000GB_all/playground/enhancer"
root <- file.path(external, "r_files/h3_GWAS/King2025")
data_root <- file.path(external, "data/King_2025/gwas_summary_stats")
revision <- file.path(external, "r_files/h3_GWAS/revision_2001bp_loop_only")
out <- file.path(root, "review_audit")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
traits <- c("crf_ny_incentive_value_index", "crf_ny_lever_presses", "pavca_ny_d5_response_bias", "pavca_ny_d5_index")
read_pairs <- function(path) {
  x <- strsplit(readLines(path), "\t", fixed = TRUE)
  unique(rbindlist(lapply(x, function(z) data.table(gene_id = z[1], SNP = z[-c(1, 2)]))))
}
old <- read_pairs(file.path(repo, "r_files/h3_GWAS/nicotine/hs.hic.duttke_telese_atac.non_tss_promoter_filtered.annot"))
new <- read_pairs(file.path(revision, "annotation/hs.revised2001bp.loop_only.hmagma.annot"))
direct <- read_pairs(file.path(repo, "r_files/h3_GWAS/nicotine/hs.exonpro.ONLY.annot"))
baseline <- old[SNP %chin% direct$SNP]
genes <- fread(file.path(revision, "annotation/gene_annotation_changes.tsv"))
selected <- genes[gene_symbol %chin% c("Lamtor1", "Klhl35", "Lrtomt", "Rnf169")]
stopifnot(nrow(selected) == 4L)
bim <- fread(file.path(external, "data/hs_data/v4/HS_genotypes_v4.bim"), header = FALSE,
             select = 2, col.names = "SNP")
mapping <- fread(file.path(data_root, "king_round8_to_rn7_snp_map.tsv"))
stopifnot(!anyDuplicated(mapping$original_snp))
duplicate_targets <- mapping[, .N, by = rn7_snp][N > 1L, rn7_snp]
collisions <- mapping[rn7_snp %chin% duplicate_targets]
collisions[, in_LD_reference := rn7_snp %chin% bim$SNP]
collisions[, n_direct_annotation_pairs := vapply(rn7_snp, function(s) sum(direct$SNP == s), integer(1))]
collisions[, n_legacy_annotation_pairs := vapply(rn7_snp, function(s) sum(old$SNP == s), integer(1))]
collisions[, n_revised_annotation_pairs := vapply(rn7_snp, function(s) sum(new$SNP == s), integer(1))]
fwrite(collisions, file.path(out, "many_to_one_targets.tsv"), sep = "\t")
summaries <- hits <- candidates <- collision_p <- vector("list", length(traits))
for (i in seq_along(traits)) {
  trait <- traits[i]
  message("Auditing ", trait)
  cm <- fread(file.path(root, paste0(trait, "_cmagma.genes.out")))
  hm <- fread(file.path(root, paste0(trait, "_hmagma_strict_duttke_telese.genes.out")))
  cm[, BH_q := p.adjust(P, "BH")]
  hm[, BH_q := p.adjust(P, "BH")]
  cm[, Bonf_p := p.adjust(P, "bonferroni")]
  hm[, Bonf_p := p.adjust(P, "bonferroni")]
  for (x in list(cm, hm)) x[, coord_key := paste(CHR, START, STOP, sep = ":")]
  stopifnot(!anyDuplicated(cm$coord_key))
  x <- merge(hm, cm[, .(coord_key, cmagma_P = P, cmagma_BH_q = BH_q,
                        cmagma_Bonf_p = Bonf_p, cmagma_NSNPS = NSNPS)], by = "coord_key", all.x = TRUE)
  x <- merge(x, genes[, .(GENE = gene_id, gene_symbol)], by = "GENE", all.x = TRUE)
  x[, hmagma_only_BH05 := BH_q < 0.05 & (is.na(cmagma_BH_q) | cmagma_BH_q >= 0.05)]
  x[, unchanged_gene_P := !is.na(cmagma_P) & P == cmagma_P]
  hits[[i]] <- x[hmagma_only_BH05 == TRUE, .(trait, gene_id = GENE, gene_symbol, CHR, START, STOP,
    hmagma_P = P, cmagma_P, hmagma_BH_q = BH_q, cmagma_BH_q, hmagma_NSNPS = NSNPS,
    cmagma_NSNPS, unchanged_gene_P)]
  p <- fread(file.path(root, paste0(trait, "_rn7.magma.tsv")))
  stopifnot(!anyDuplicated(p$SNP), all(is.finite(p$p) & p$p > 0 & p$p <= 1))
  summaries[[i]] <- data.table(trait,
    cMAGMA_tested_coordinate_IDs = nrow(cm), legacy_HMAGMA_tested_gene_IDs = nrow(hm),
    cMAGMA_BH05 = sum(cm$BH_q < 0.05), legacy_HMAGMA_BH05 = sum(hm$BH_q < 0.05),
    difference_in_significant_row_counts = sum(hm$BH_q < 0.05) - sum(cm$BH_q < 0.05),
    legacy_HMAGMA_only_BH05_by_coordinate_match = sum(x$hmagma_only_BH05),
    only_BH05_with_identical_gene_P = sum(x$hmagma_only_BH05 & x$unchanged_gene_P),
    cMAGMA_Bonf05 = sum(cm$Bonf_p < 0.05), legacy_HMAGMA_Bonf05 = sum(hm$Bonf_p < 0.05),
    legacy_HMAGMA_only_Bonf05 = sum(x$Bonf_p < 0.05 & (is.na(x$cmagma_Bonf_p) | x$cmagma_Bonf_p >= 0.05)),
    GWAS_SNPs = nrow(p), GWAS_SNPs_in_LD = sum(p$SNP %chin% bim$SNP),
    fraction_unmatched_to_LD = mean(!p$SNP %chin% bim$SNP),
    lambda_GC = median(qchisq(p$p, df = 1, lower.tail = FALSE)) / qchisq(0.5, df = 1),
    lambda_GC_without_collision_targets = median(qchisq(p[!SNP %chin% duplicate_targets, p], df = 1, lower.tail = FALSE)) / qchisq(0.5, df = 1))
  z <- x[GENE %chin% selected$gene_id]
  z[, n_direct_SNPs_in_GWAS := vapply(GENE, function(g) sum(baseline[gene_id == g, SNP] %chin% p$SNP), integer(1))]
  z[, n_reference_direct_SNPs := vapply(GENE, function(g) nrow(baseline[gene_id == g]), integer(1))]
  candidates[[i]] <- z[, .(trait, gene_id = GENE, gene_symbol, n_reference_direct_SNPs,
    n_direct_SNPs_in_GWAS, cmagma_P, hmagma_P = P, hmagma_BH_q = BH_q,
    hmagma_Bonf_p = Bonf_p, hmagma_NSNPS = NSNPS)]
  raw_collision <- rbindlist(lapply(unique(sub(":.*", "", collisions$original_snp)), function(chr) {
    raw <- fread(file.path(data_root, "raw_chrgwas_mlma", trait, paste0("regressedlr_", trait, "_chrgwas", chr, ".mlma")), select = c("SNP", "p"))
    raw[SNP %chin% collisions$original_snp]
  }))
  setnames(raw_collision, c("SNP", "p"), c("original_snp", "original_P"))
  cp <- merge(collisions, raw_collision, by = "original_snp", all.x = TRUE)
  cp <- merge(cp, p[, .(rn7_snp = SNP, retained_P = p)], by = "rn7_snp", all.x = TRUE)
  cp[, trait := trait]
  cp[, minimum_original_P := min(original_P), by = rn7_snp]
  stopifnot(all(abs(cp$retained_P - cp$minimum_original_P) < 1e-10))
  collision_p[[i]] <- cp
}
fwrite(rbindlist(summaries), file.path(out, "legacy_count_and_LD_audit.tsv"), sep = "\t")
fwrite(rbindlist(hits), file.path(out, "legacy_HMAGMA_only_hits_checked.tsv"), sep = "\t")
fwrite(rbindlist(candidates), file.path(out, "four_candidate_mapping_audit.tsv"), sep = "\t")
fwrite(rbindlist(collision_p), file.path(out, "many_to_one_P_value_audit.tsv"), sep = "\t")
print(rbindlist(summaries))
print(collisions)
message("Review audit complete: ", out)
