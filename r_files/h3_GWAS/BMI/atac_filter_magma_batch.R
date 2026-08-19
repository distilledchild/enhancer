#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
library(GenomicRanges)
library(dplyr)
library(readr)
library(stringr)

args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0) {
  stop("Usage: Rscript atac_filter_magma_batch.R <phenotype_base>\nExample: Rscript atac_filter_magma_batch.R results.bmi_wo_tail")
}
base_name <- args[1]

# Path Definitions
dir_path <- "/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/BMI"
path.out.hmagma <- file.path(dir_path, str_glue("{base_name}_rn7_hmagma.genes.out"))
path.out.cmagma <- file.path(dir_path, str_glue("{base_name}_rn7_cmagma.genes.out"))
path.loop <- "/Users/pete/Desktop/playground/enhancer/figures/submission/lt2mb/df_final_loop_sub.4.any.lt2mb.ENSEMBL.mid.mid.final.filter.200kb.csv"
path.narrowpeak.atac <- "/Users/pete/Library/CloudStorage/Dropbox-UTHSCGGI/K P/Gateway_to_Hao/enhancer/data/Duttke2022_snATAC_peaks_rn7.narrowPeak"

df.hmagma <- read.table(path.out.hmagma, header=TRUE, stringsAsFactors=FALSE)
df.hmagma$FDR <- p.adjust(df.hmagma$P, method="BH")

df.cmagma <- read.table(path.out.cmagma, header=TRUE, stringsAsFactors=FALSE)
df.cmagma$FDR <- p.adjust(df.cmagma$P, method="BH")

df.hmagma_sig <- df.hmagma %>% filter(FDR < 0.05)
df.cmagma_sig <- df.cmagma %>% filter(FDR < 0.05)
df.magma_sig <- bind_rows(df.hmagma_sig, df.cmagma_sig) %>% distinct(GENE, .keep_all=TRUE)

print(str_glue("Significant MAGMA genes (FDR < 0.05) for {base_name}: {nrow(df.magma_sig)}"))

if(nrow(df.magma_sig) == 0) {
  print("No significant genes found. Exiting.")
  quit(save="no")
}

df.loop <- read_csv(path.loop, show_col_types=FALSE)
df.loop_sig <- bind_rows(
  df.loop %>% filter(ensembl_id %in% df.magma_sig$GENE, is_promoter1) %>% mutate(WHERE="DOWN"),
  df.loop %>% filter(ensembl_id %in% df.magma_sig$GENE, is_promoter2) %>% mutate(WHERE="UP")
) %>% distinct(loop.id, ensembl_id, .keep_all=TRUE)

if(nrow(df.loop_sig) == 0) {
  print("No loops matched! Exiting.")
  quit(save="no")
}

df.enhancer <- bind_rows(
  df.loop_sig %>% filter(WHERE == "UP") %>% select(chr=chr2, start=y1, end=y2, loop.id, ensembl_id, gene_name),
  df.loop_sig %>% filter(WHERE == "DOWN") %>% select(chr=chr1, start=x1, end=x2, loop.id, ensembl_id, gene_name)
) %>% mutate(start = as.numeric(start), end = as.numeric(end))

gr.enhancer <- GRanges(seqnames = df.enhancer$chr, ranges = IRanges(df.enhancer$start + 1, df.enhancer$end))
df.atac <- read_tsv(path.narrowpeak.atac, col_names = c("chr", "start", "end", "name", "score", "strand", "fc", "neglog10p", "neglog10q", "summit"), show_col_types=FALSE)
gr.atac.rn7.1based <- GRanges(df.atac$chr, IRanges(df.atac$start + 1, df.atac$end))
gr.atac_union <- GenomicRanges::reduce(gr.atac.rn7.1based)

enhancer_hits <- countOverlaps(gr.enhancer, gr.atac_union, minoverlap = 50) > 0
df.atac_valid <- df.enhancer[enhancer_hits, ]
valid_genes <- unique(df.atac_valid$ensembl_id)

print(str_glue("ATAC-validated unique genes for {base_name}: {length(valid_genes)}"))

df.final_res <- df.magma_sig %>% 
  filter(GENE %in% valid_genes) %>%
  left_join(df.atac_valid %>% select(GENE=ensembl_id, gene_name) %>% distinct(), by="GENE") %>%
  arrange(FDR) %>%
  select(GENE, gene_name, CHR, START, STOP, ZSTAT, P, FDR)

output_file_name <- str_replace(base_name, "^results\\.", "")
path.csv.atac_validated <- file.path(dir_path, str_glue("ATAC_validated_{output_file_name}_genes.csv"))
write_csv(df.final_res, path.csv.atac_validated)

print(str_glue("Results saved to: {path.csv.atac_validated}"))
