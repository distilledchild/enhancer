options(stringsAsFactors = FALSE)

outdir <- "/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/BMI"
genedef_path <- "/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/nicotine/genedef_ensembl.tsv"
gene_name_path <- "/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/nicotine/cmagma_vs_hmagma_duttke_telese_non_tss_promoter_snp_counts.csv"

phenotypes <- c(
  "bmi_w_tail",
  "bmi_wo_tail",
  "body_weight_g",
  "length_w_tail_cm",
  "length_wo_tail_cm"
)

genedef <- read.table(genedef_path, sep = "\t", stringsAsFactors = FALSE)
names(genedef) <- c("index", "chr", "start", "end", "strand", "ensg")

gene_names <- read.csv(gene_name_path, stringsAsFactors = FALSE)
gene_names <- gene_names[!is.na(gene_names$gene_name) & gene_names$gene_name != "", c("ensg", "gene_name")]
gene_names <- gene_names[!duplicated(gene_names$ensg), ]

read_result <- function(phenotype, method) {
  tag <- if (method == "H-MAGMA") "hmagma" else "cmagma"
  path <- file.path(outdir, sprintf("results.%s_rn7_%s.genes.out", phenotype, tag))
  x <- read.table(path, header = TRUE, stringsAsFactors = FALSE)
  x$phenotype <- phenotype
  x$method <- method
  x$source_file <- basename(path)
  x$FDR_BH <- p.adjust(x$P, method = "BH")
  x$P_Bonferroni <- p.adjust(x$P, method = "bonferroni")
  x$Bonferroni_0_05 <- x$P_Bonferroni < 0.05
  x$FDR_0_05 <- x$FDR_BH < 0.05
  x$FDR_0_10 <- x$FDR_BH < 0.10
  if (method == "cMAGMA") {
    x$ensembl_id <- genedef$ensg[match(x$GENE, genedef$index)]
  } else {
    x$ensembl_id <- x$GENE
  }
  x$gene_name <- gene_names$gene_name[match(x$ensembl_id, gene_names$ensg)]
  x
}

all_results <- do.call(
  rbind,
  unlist(
    lapply(phenotypes, function(pheno) {
      list(read_result(pheno, "H-MAGMA"), read_result(pheno, "cMAGMA"))
    }),
    recursive = FALSE
  )
)

summary_rows <- do.call(rbind, lapply(split(all_results, list(all_results$phenotype, all_results$method), drop = TRUE), function(x) {
  i <- which.min(x$P)
  data.frame(
    phenotype = x$phenotype[1],
    method = x$method[1],
    n_genes = nrow(x),
    n_fdr_0_05 = sum(x$FDR_0_05),
    n_fdr_0_10 = sum(x$FDR_0_10),
    n_bonferroni_0_05 = sum(x$Bonferroni_0_05),
    top_gene_raw = x$GENE[i],
    top_ensembl_id = x$ensembl_id[i],
    top_gene_name = x$gene_name[i],
    top_p = x$P[i],
    top_fdr_bh = x$FDR_BH[i],
    top_p_bonferroni = x$P_Bonferroni[i]
  )
}))

comparison_rows <- do.call(rbind, lapply(phenotypes, function(pheno) {
  h <- all_results[all_results$phenotype == pheno & all_results$method == "H-MAGMA", ]
  c <- all_results[all_results$phenotype == pheno & all_results$method == "cMAGMA", ]
  h_fdr <- unique(na.omit(h$ensembl_id[h$FDR_0_10]))
  c_fdr <- unique(na.omit(c$ensembl_id[c$FDR_0_10]))
  h_bonf <- unique(na.omit(h$ensembl_id[h$Bonferroni_0_05]))
  c_bonf <- unique(na.omit(c$ensembl_id[c$Bonferroni_0_05]))
  data.frame(
    phenotype = pheno,
    h_fdr_0_10 = length(h_fdr),
    c_fdr_0_10 = length(c_fdr),
    h_only_fdr_0_10 = length(setdiff(h_fdr, c_fdr)),
    c_only_fdr_0_10 = length(setdiff(c_fdr, h_fdr)),
    overlap_fdr_0_10 = length(intersect(h_fdr, c_fdr)),
    h_bonferroni_0_05 = length(h_bonf),
    c_bonferroni_0_05 = length(c_bonf),
    h_only_bonferroni_0_05 = length(setdiff(h_bonf, c_bonf)),
    c_only_bonferroni_0_05 = length(setdiff(c_bonf, h_bonf)),
    overlap_bonferroni_0_05 = length(intersect(h_bonf, c_bonf))
  )
}))

significant <- all_results[all_results$FDR_0_10 | all_results$Bonferroni_0_05, ]
significant <- significant[order(significant$phenotype, significant$method, significant$P), ]
significant <- significant[, c(
  "phenotype", "method", "GENE", "ensembl_id", "gene_name", "CHR", "START", "STOP",
  "NSNPS", "NPARAM", "N", "ZSTAT", "P", "FDR_BH", "P_Bonferroni",
  "FDR_0_05", "FDR_0_10", "Bonferroni_0_05", "source_file"
)]

write.csv(summary_rows[order(summary_rows$phenotype, summary_rows$method), ],
          file.path(outdir, "BMI_magma_rn7_summary.csv"), row.names = FALSE)
write.csv(comparison_rows,
          file.path(outdir, "BMI_magma_rn7_method_comparison.csv"), row.names = FALSE)
write.csv(significant,
          file.path(outdir, "BMI_magma_rn7_significant_genes.csv"), row.names = FALSE)

print(summary_rows[order(summary_rows$phenotype, summary_rows$method), ], row.names = FALSE)
print(comparison_rows, row.names = FALSE)
