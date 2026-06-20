suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

base_dir <- "/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/King2025"
report_dir <- file.path(base_dir, "reports")
figure_dir <- file.path(report_dir, "repeated_gene_locus_figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

hmagma_annot <- "/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/nicotine/hs.hic.duttke_telese_atac.non_tss_promoter_filtered.annot"
cmagma_annot <- "/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/nicotine/hs.exonpro.ONLY.annot"
snp_count_file <- "/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/nicotine/cmagma_vs_hmagma_duttke_telese_non_tss_promoter_snp_counts.csv"
atac_file <- "/Users/pete/Library/CloudStorage/Dropbox-UTHSCGGI/K P/Gateway_to_Hao/enhancer/data/Duttke2022_snATAC_peaks_rn7.narrowPeak"
where_rds <- "/Users/pete/Library/CloudStorage/Dropbox-UTHSCGGI/K P/Gateway_to_Hao/enhancer/data/df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3_final_200kb.rds"
key_gene_file <- file.path(report_dir, "selected_1_3_traits_hmagma_only_key_genes.tsv")

traits <- data.table(
  trait = c("crf_ny_incentive_value_index", "crf_ny_lever_presses", "pavca_ny_d5_response_bias"),
  label = c("CRF incentive value index", "CRF lever presses", "PavCA day 5 response bias"),
  n = c(1594L, 1594L, 1585L)
)

target_genes <- c("Klhl35", "Lamtor1", "Lrtomt", "Rnf169")

parse_snp <- function(x) {
  y <- tstrsplit(x, ":", fixed = TRUE)
  data.table(chr = y[[1]], pos = as.integer(y[[2]]))
}

read_selected_annotation <- function(path, target_keys, key_mode = c("gene_id", "region")) {
  key_mode <- match.arg(key_mode)
  lines <- readLines(path, warn = FALSE)
  out <- vector("list", length(target_keys))
  i <- 0L

  for (line in lines) {
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(parts) < 3L) next

    key <- if (key_mode == "gene_id") parts[[1]] else parts[[2]]
    if (!key %in% target_keys) next

    snps <- parts[-c(1, 2)]
    parsed <- parse_snp(snps)
    i <- i + 1L
    out[[i]] <- data.table(ann_key = key, region = parts[[2]], snp = snps, chr = parsed$chr, pos = parsed$pos)
  }

  if (i == 0L) return(data.table(ann_key = character(), region = character(), snp = character(), chr = character(), pos = integer()))
  rbindlist(out[seq_len(i)])
}

key <- fread(key_gene_file)
key_bh <- unique(key[criterion == "bh_fdr05" & gene_name %in% target_genes,
  .(trait, hmagma_gene_id, gene_name, CHR, START, STOP, hmagma_NSNP, cmagma_NSNP,
    hmagma_ZSTAT, cmagma_ZSTAT, hmagma_P, cmagma_P, hmagma_BH_FDR, cmagma_BH_FDR,
    n_snp_cmagma, n_snp_hmagma_filtered, added_snps_filtered, noise_snps_removed)
])

gene_meta <- unique(key_bh[, .(
  gene_id = hmagma_gene_id,
  gene_name,
  chr = as.character(CHR),
  start = as.integer(START),
  stop = as.integer(STOP),
  region = paste(CHR, START, STOP, sep = ":")
)])

h_ann <- read_selected_annotation(hmagma_annot, gene_meta$gene_id, "gene_id")
h_ann <- merge(h_ann, gene_meta[, .(gene_id, gene_name)], by.x = "ann_key", by.y = "gene_id", all.x = TRUE)
h_ann[, assignment_source := "H-MAGMA_all"]

c_ann <- read_selected_annotation(cmagma_annot, gene_meta$region, "region")
c_ann <- merge(c_ann, gene_meta[, .(region, gene_name, gene_id)], by = "region", all.x = TRUE)
c_ann[, assignment_source := "cMAGMA"]

assigned <- unique(rbindlist(list(
  h_ann[, .(gene_id = ann_key, gene_name, snp, chr, pos, in_hmagma = TRUE, in_cmagma = FALSE)],
  c_ann[, .(gene_id, gene_name, snp, chr, pos, in_hmagma = FALSE, in_cmagma = TRUE)]
), fill = TRUE))

assigned <- assigned[, .(
  in_hmagma = any(in_hmagma),
  in_cmagma = any(in_cmagma)
), by = .(gene_id, gene_name, snp, chr, pos)]
assigned[, assignment := fifelse(in_cmagma & in_hmagma, "cMAGMA local annotation", "H-MAGMA added")]

atac <- fread(atac_file, header = FALSE, select = 1:3)
setnames(atac, c("chr", "start0", "end"))
atac[, `:=`(chr = sub("^chr", "", chr), start = as.integer(start0) + 1L, end = as.integer(end))]
atac <- atac[chr %in% gene_meta$chr]

assigned[, atac_overlap := FALSE]
for (g in target_genes) {
  idx <- which(assigned$gene_name == g)
  if (!length(idx)) next
  snp_dt <- assigned[idx]
  peak_dt <- atac[chr %in% snp_dt$chr]
  overlaps <- vapply(seq_len(nrow(snp_dt)), function(i) {
    any(peak_dt$chr == snp_dt$chr[i] & peak_dt$start <= snp_dt$pos[i] & peak_dt$end >= snp_dt$pos[i])
  }, logical(1))
  assigned$atac_overlap[idx] <- overlaps
}

pval_list <- list()
for (i in seq_len(nrow(traits))) {
  tr <- traits$trait[i]
  pval_path <- file.path(base_dir, paste0(tr, "_rn7.magma.tsv"))
  p <- fread(pval_path)
  parsed <- parse_snp(p$SNP)
  p[, `:=`(chr = parsed$chr, pos = parsed$pos, trait = tr, trait_label = traits$label[i], N = traits$n[i])]
  pval_list[[tr]] <- p[chr %in% gene_meta$chr]
}
pvals <- rbindlist(pval_list)
pvals[, neglog10p := -log10(pmax(p, .Machine$double.xmin))]

where <- as.data.table(readRDS(where_rds))
loop_dt <- where[grepl(paste(gene_meta$gene_id, collapse = "|"), gene_id)]
if (nrow(loop_dt)) {
  loop_dt[, `:=`(
    x1n = as.integer(x1), x2n = as.integer(x2), y1n = as.integer(y1), y2n = as.integer(y2),
    chr_clean = sub("^chr", "", chr1)
  )]
  loop_dt[, `:=`(
    enhancer_start = fifelse(WHERE == "DOWN", x1n, y1n),
    enhancer_end = fifelse(WHERE == "DOWN", x2n, y2n),
    promoter_start = fifelse(WHERE == "DOWN", y1n, x1n),
    promoter_end = fifelse(WHERE == "DOWN", y2n, x2n)
  )]
  loop_dt[, `:=`(
    enhancer_mid = (enhancer_start + enhancer_end) / 2,
    promoter_mid = (promoter_start + promoter_end) / 2
  )]
  loop_dt[, gene_id_short := sub(".*:(ENSRNOG[0-9]+):.*", "\\1", gene_id)]
  loop_dt <- merge(loop_dt, gene_meta[, .(gene_id, gene_name)], by.x = "gene_id_short", by.y = "gene_id", all.x = TRUE)
  if ("gene_name.y" %in% names(loop_dt)) {
    loop_dt[, gene_name := gene_name.y]
  }
} else {
  loop_dt <- data.table()
}

bio <- data.table(
  gene_name = c("Klhl35", "Lamtor1", "Lrtomt", "Rnf169"),
  biological_priority = c("Exploratory", "Moderate", "Exploratory", "Exploratory"),
  core_biology = c(
    "Kelch-like protein; predicted CUL3-RING ubiquitin ligase substrate adaptor/proteostasis biology.",
    "Lysosomal Ragulator anchor/adaptor upstream of mTORC1; links lysosomal nutrient signaling to MAPK/MTOR activity.",
    "Readthrough/fusion-like leucine-rich transmembrane/O-methyltransferase locus; TOMT ortholog biology is linked to mechanotransduction and catecholamine methyltransferase-like function.",
    "Ubiquitin/chromatin DNA-damage-response regulator; limits/competes in RNF168/53BP1-type repair-factor recruitment."
  ),
  neuro_addiction_interpretation = c(
    "No strong direct addiction evidence found yet; repeated H-MAGMA signal makes it a regulatory-contact candidate, not a mature biological lead.",
    "Best biological lead among the four because mTORC1 signaling is repeatedly implicated in neuronal plasticity and drug reward/cue learning; the gene-level claim still needs expression and locus-level support.",
    "Indirect candidate only. The catecholamine-methyltransferase-like naming is tempting, but published functional work is mainly auditory hair-cell mechanotransduction/deafness, not addiction.",
    "No direct reward/addiction evidence found yet; could be relevant through chromatin/ubiquitin stress-response biology, but this remains speculative."
  ),
  key_references = c(
    "NCBI Gene 283212; Human Protein Atlas KLHL35 brain single-cell expression",
    "PMID:31001086; PMID:29039413; PMID:24666346",
    "PMID:18953341; PMID:28504928",
    "PMID:22733822; NCBI Gene 254225"
  )
)

counts <- fread(snp_count_file)
summary_tbl <- key_bh[, .(
  n_traits = uniqueN(trait),
  traits = paste(sort(unique(trait)), collapse = "; "),
  max_hmagma_Z = max(hmagma_ZSTAT, na.rm = TRUE),
  min_hmagma_P = min(hmagma_P, na.rm = TRUE),
  min_hmagma_BH_FDR = min(hmagma_BH_FDR, na.rm = TRUE),
  hmagma_NSNP_range = paste(range(hmagma_NSNP, na.rm = TRUE), collapse = "-"),
  cmagma_NSNP_observed = paste(sort(unique(na.omit(cmagma_NSNP))), collapse = ";")
), by = .(hmagma_gene_id, gene_name, CHR, START, STOP)]

summary_tbl <- merge(summary_tbl, counts[, .(
  hmagma_gene_id = ensg,
  n_snp_cmagma_annotation = n_snp_cmagma,
  n_snp_hmagma_filtered_annotation = n_snp_hmagma_filtered,
  added_snps_filtered_annotation = added_snps_filtered,
  noise_snps_removed_annotation = noise_snps_removed
)], by = "hmagma_gene_id", all.x = TRUE)
summary_tbl <- merge(summary_tbl, bio, by = "gene_name", all.x = TRUE)
setorder(summary_tbl, biological_priority, gene_name)

trait_stats <- key_bh[gene_name %in% target_genes]
setorder(trait_stats, gene_name, trait)

assigned_long <- merge(assigned, pvals[, .(snp = SNP, trait, p, neglog10p)], by = "snp", all.x = TRUE)
assigned_status <- assigned_long[, .(has_king_gwas_p_any = any(!is.na(p))), by = .(gene_id, gene_name, chr, pos, snp, assignment, atac_overlap)]
assigned_wide <- dcast(
  assigned_long,
  gene_id + gene_name + chr + pos + snp + assignment + atac_overlap ~ trait,
  value.var = c("p", "neglog10p")
)
assigned_wide <- merge(assigned_wide, assigned_status, by = c("gene_id", "gene_name", "chr", "pos", "snp", "assignment", "atac_overlap"), all.x = TRUE)
setorder(assigned_wide, gene_name, pos)

fwrite(summary_tbl, file.path(report_dir, "king_1_3_repeated_genes_biology_annotation.tsv"), sep = "\t")
fwrite(trait_stats, file.path(report_dir, "king_1_3_repeated_genes_trait_stats.tsv"), sep = "\t")
fwrite(assigned_wide, file.path(report_dir, "king_1_3_repeated_genes_snp_assignments.tsv"), sep = "\t")

palette_assign <- c("local SNPs" = "grey75", "cMAGMA local annotation" = "#2f6f9f", "H-MAGMA added" = "#c07a2c")

plot_one_gene <- function(g) {
  meta <- gene_meta[gene_name == g][1]
  snps <- assigned[gene_name == g]
  x_min <- min(c(meta$start, meta$stop, snps$pos), na.rm = TRUE) - 100000L
  x_max <- max(c(meta$start, meta$stop, snps$pos), na.rm = TRUE) + 100000L

  local_p <- pvals[chr == meta$chr & pos >= x_min & pos <= x_max]
  local_p <- merge(local_p, snps[, .(SNP = snp, assignment)], by = "SNP", all.x = TRUE)
  local_p[is.na(assignment), assignment := "local SNPs"]

  local_atac <- atac[chr == meta$chr & end >= x_min & start <= x_max]
  local_loop <- loop_dt[gene_name == g & chr_clean == meta$chr &
    ((enhancer_start <= x_max & enhancer_end >= x_min) | (promoter_start <= x_max & promoter_end >= x_min))]

  p_gwas <- ggplot(local_p, aes(x = pos / 1e6, y = neglog10p)) +
    geom_point(data = local_p[assignment == "local SNPs"], color = "grey78", size = 0.45, alpha = 0.7) +
    geom_point(data = local_p[assignment != "local SNPs"], aes(color = assignment), size = 1.35, alpha = 0.95) +
    geom_vline(xintercept = c(meta$start, meta$stop) / 1e6, linetype = "dashed", linewidth = 0.25, color = "grey40") +
    facet_wrap(~trait_label, ncol = 1, scales = "free_y") +
    scale_color_manual(values = palette_assign[c("cMAGMA local annotation", "H-MAGMA added")], drop = FALSE) +
    labs(
      title = paste0(g, " locus: King 2025 selected phenotypes"),
      subtitle = "Strict H-MAGMA uses frontal-cortex Hi-C contacts plus Duttke/Telese distal ATAC peaks with TSS/promoter-overlapping peaks excluded",
      x = NULL,
      y = expression(-log[10](GWAS~P)),
      color = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey92", color = "grey70"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )

  gene_track <- data.table(
    xmin = meta$start / 1e6,
    xmax = meta$stop / 1e6,
    ymin = 0.2,
    ymax = 0.55,
    label = g
  )
  snp_track <- copy(snps)
  snp_track[, x := pos / 1e6]

  p_reg <- ggplot() +
    geom_rect(data = local_atac, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 1.0, ymax = 1.18),
      fill = "grey60", color = NA, alpha = 0.6) +
    geom_segment(data = snp_track[assignment == "cMAGMA local annotation"], aes(x = x, xend = x, y = 1.45, yend = 1.85),
      color = palette_assign[["cMAGMA local annotation"]], linewidth = 0.35) +
    geom_segment(data = snp_track[assignment == "H-MAGMA added"], aes(x = x, xend = x, y = 1.88, yend = 2.28),
      color = palette_assign[["H-MAGMA added"]], linewidth = 0.35) +
    geom_rect(data = gene_track, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "#1f2d3a", color = "#1f2d3a") +
    geom_text(data = gene_track, aes(x = (xmin + xmax) / 2, y = 0.0, label = label),
      size = 3.2, fontface = "bold") +
    scale_y_continuous(
      breaks = c(0.38, 1.09, 1.65, 2.08, 2.75),
      labels = c("gene", "ATAC", "cMAGMA local SNP", "H-MAGMA-added SNP", "Hi-C loop"),
      limits = c(-0.15, 3.2)
    ) +
    labs(x = paste0("chr", meta$chr, " position (Mb)"), y = NULL) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(color = "grey20")
    )

  if (nrow(local_loop)) {
    p_reg <- p_reg +
      geom_curve(data = local_loop, aes(x = enhancer_mid / 1e6, xend = promoter_mid / 1e6, y = 2.72, yend = 2.72),
        curvature = 0.35, color = "#8b3f88", linewidth = 0.65, alpha = 0.8)
  }

  fig <- p_gwas / p_reg + plot_layout(heights = c(3.2, 1.15))
  png_path <- file.path(figure_dir, paste0("king_", g, "_locus_hmagma_browser_style.png"))
  pdf_path <- file.path(figure_dir, paste0("king_", g, "_locus_hmagma_browser_style.pdf"))
  ggsave(png_path, fig, width = 10.5, height = 8.2, dpi = 300)
  ggsave(pdf_path, fig, width = 10.5, height = 8.2)

  data.table(gene_name = g, png = png_path, pdf = pdf_path, x_min = x_min, x_max = x_max,
    n_local_snps = nrow(local_p), n_atac_peaks = nrow(local_atac), n_loop_arcs = nrow(local_loop))
}

plot_index <- rbindlist(lapply(target_genes, plot_one_gene))
fwrite(plot_index, file.path(report_dir, "king_1_3_repeated_genes_locus_figure_index.tsv"), sep = "\t")

md <- c(
  "# King 2025 Repeated H-MAGMA-only Gene Annotation",
  "",
  "## Repeated Genes",
  "",
  paste(capture.output(print(summary_tbl)), collapse = "\n"),
  "",
  "## Interpretation",
  "",
  "- `Lamtor1` is the strongest biological lead among the repeated genes because it sits in the lysosomal Ragulator/mTORC1 axis, which has published links to neuronal plasticity and drug reward/cue learning.",
  "- `Klhl35`, `Lrtomt`, and `Rnf169` should be treated as exploratory regulatory-contact candidates until brain expression, pathway context, and independent phenotype/dataset recurrence are confirmed.",
  "- All four genes are repeated H-MAGMA-only hits across the three selected King phenotypes at BH-FDR < 0.05; several also pass Bonferroni in CRF lever pressing or PavCA response bias.",
  "- The local cMAGMA annotation contains a small number of SNPs for these genes, but those cMAGMA-local SNPs have no matching King GWAS p-values in the rn7 p-value files and no cMAGMA `genes.out` rows were produced. The observed gene-level signal therefore comes from strict H-MAGMA added SNPs.",
  "",
  "## Output Files",
  "",
  paste0("- `", file.path(report_dir, "king_1_3_repeated_genes_biology_annotation.tsv"), "`"),
  paste0("- `", file.path(report_dir, "king_1_3_repeated_genes_trait_stats.tsv"), "`"),
  paste0("- `", file.path(report_dir, "king_1_3_repeated_genes_snp_assignments.tsv"), "`"),
  paste0("- `", file.path(report_dir, "king_1_3_repeated_genes_locus_figure_index.tsv"), "`"),
  paste0("- Figures: `", figure_dir, "`")
)
writeLines(md, file.path(report_dir, "king_1_3_repeated_genes_annotation_summary.md"))

print(summary_tbl)
print(plot_index)
