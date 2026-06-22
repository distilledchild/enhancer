####################################################
# ATAC-seq Validation of MAGMA genes
# Body Weight GWAS (rn7) x Hi-C loop anchors
# Filter H-MAGMA / cMAGMA results using ATAC-seq overlap
####################################################

library(tidyverse)
library(GenomicRanges) # for genomic interval overlap and ranges operation

options(scipen = 999) # prevent scientific notation in base R
options(pillar.sigfig = 10) # display up to 10 significant digits in tibbles

setwd("/Users/pete/Desktop/playground/enhancer") # Mac
getwd() # /Users/pete/Desktop/playground/enhancer

####################################################
# 1. Loading H-MAGMA & cMAGMA results
####################################################
path.out.hmagma <- "/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/BMI/results.body_weight_g_rn7_hmagma.genes.out"
path.out.cmagma <- "/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/BMI/results.body_weight_g_rn7_cmagma.genes.out"

df.hmagma <- read.table(path.out.hmagma, header=TRUE, stringsAsFactors=FALSE)
df.hmagma$FDR <- p.adjust(df.hmagma$P, method="BH")

df.hmagma_sig_05 <- df.hmagma %>% 
  filter(FDR < 0.05) %>% 
  arrange(FDR)

df.hmagma_sig_10 <- df.hmagma %>% 
  filter(FDR < 0.1 & FDR >= 0.05) %>% 
  arrange(FDR)

print("=== H-MAGMA Results for Body Weight (3.5M unpruned) ===") # === H-MAGMA Results for Body Weight (3.5M unpruned) ===
print(str_glue("Total genes with FDR < 0.05: {nrow(df.hmagma_sig_05)}")) # Total genes with FDR < 0.05: 1762
print(str_glue("Total genes with 0.05 <= FDR < 0.1: {nrow(df.hmagma_sig_10)}")) # Total genes with 0.05 <= FDR < 0.1: 1006

if (file.exists(path.out.cmagma)) {
  df.cmagma <- read.table(path.out.cmagma, header=TRUE, stringsAsFactors=FALSE)
  df.cmagma$FDR <- p.adjust(df.cmagma$P, method="BH")
  
  df.cmagma_sig_05 <- df.cmagma %>% 
    filter(FDR < 0.05) %>% 
    arrange(FDR)
    
  df.cmagma_sig_10 <- df.cmagma %>% 
    filter(FDR < 0.1 & FDR >= 0.05) %>% 
    arrange(FDR)
    
  print("=== cMAGMA Results ===") # === cMAGMA Results ===
  print(str_glue("Total genes with FDR < 0.05: {nrow(df.cmagma_sig_05)}")) # (number)
  print(str_glue("Total genes with 0.05 <= FDR < 0.1: {nrow(df.cmagma_sig_10)}")) # (number)
} else {
  print("=== cMAGMA Results ===") # === cMAGMA Results ===
  print("cMAGMA not finished yet.") # cMAGMA not finished yet.
}

df.magma_sig <- df.hmagma_sig_05
print(str_glue("Proceeding ATAC validation with H-MAGMA FDR < 0.05 genes: {nrow(df.magma_sig)}")) # Proceeding ATAC validation with H-MAGMA FDR < 0.05 genes: 1762

####################################################
# 2. Loading Hi-C loops
####################################################
path.rds.loop <- "/Users/pete/Library/CloudStorage/Dropbox-UTHSCGGI/K P/Gateway_to_Hao/enhancer/r_files/rds/df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3_final_200kb.rds"
df.loop <- readRDS(path.rds.loop)

df.loop <- df.loop %>%
  mutate(ensembl_id = str_extract(gene_id, "ENSRNOG\\d+")) %>%
  filter(!is.na(ensembl_id), !is.na(WHERE)) %>%
  mutate(x1 = as.numeric(x1), x2 = as.numeric(x2), y1 = as.numeric(y1), y2 = as.numeric(y2))

print(str_glue("Loops with Ensembl IDs: {nrow(df.loop)}")) # Loops with Ensembl IDs: 38780

df.loop_sig <- df.loop %>% 
  filter(ensembl_id %in% df.magma_sig$GENE)

print(str_glue("Loops connected to significant MAGMA genes: {nrow(df.loop_sig)}")) # Loops connected to significant MAGMA genes: 4894

if(nrow(df.loop_sig) == 0) {
  stop("No loops matched!")
}

####################################################
# 3. Extracting Enhancer anchors
####################################################
df.enhancer <- bind_rows(
  df.loop_sig %>% filter(WHERE == "UP") %>% select(chr=chr2, start=y1, end=y2, loop.id, ensembl_id, gene_name),
  df.loop_sig %>% filter(WHERE == "DOWN") %>% select(chr=chr1, start=x1, end=x2, loop.id, ensembl_id, gene_name)
) %>% 
  mutate(start = as.numeric(start), end = as.numeric(end))

gr.enhancer <- GRanges(
  seqnames = df.enhancer$chr,
  ranges = IRanges(df.enhancer$start + 1, df.enhancer$end)
)

####################################################
# 4. Loading ATAC-seq peaks
####################################################
path.narrowpeak.atac <- "/Users/pete/Library/CloudStorage/Dropbox-UTHSCGGI/K P/Gateway_to_Hao/enhancer/data/Duttke2022_snATAC_peaks_rn7.narrowPeak"
df.atac <- read_tsv(path.narrowpeak.atac, col_names = c("chr", "start", "end", "name", "score", "strand", "fc", "neglog10p", "neglog10q", "summit"), show_col_types=FALSE)
gr.atac <- GRanges(df.atac$chr, IRanges(df.atac$start + 1, df.atac$end))
gr.atac_union <- GenomicRanges::reduce(gr.atac)

####################################################
# 5. Overlapping with ATAC peaks
####################################################
enhancer_hits <- countOverlaps(gr.enhancer, gr.atac_union, minoverlap = 50) > 0
print(str_glue("Enhancers with ATAC peak overlap: {sum(enhancer_hits)} out of {length(gr.enhancer)}")) # Enhancers with ATAC peak overlap: 1334 out of 4894

df.atac_valid <- df.enhancer[enhancer_hits, ]
valid_genes <- unique(df.atac_valid$ensembl_id)

print(str_glue("ATAC-validated unique genes: {length(valid_genes)}")) # ATAC-validated unique genes: 522

####################################################
# 6. Save results
####################################################
df.final_res <- df.magma_sig %>% 
  filter(GENE %in% valid_genes) %>%
  left_join(df.atac_valid %>% select(GENE=ensembl_id, gene_name) %>% distinct(), by="GENE") %>%
  arrange(FDR) %>%
  select(GENE, gene_name, CHR, START, STOP, ZSTAT, P, FDR)

path.csv.atac_validated <- "/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/BMI/ATAC_validated_BodyWeight_genes.csv"
write_csv(df.final_res, path.csv.atac_validated)

print("=== TOP 20 ATAC-VALIDATED GENES ===") # === TOP 20 ATAC-VALIDATED GENES ===
print(head(df.final_res, 20)) # dataframe
print(str_glue("Results saved to: {path.csv.atac_validated}")) # Results saved to: /Users/pete/Desktop/playground/enhancer/data/BMI/ATAC_validated_BodyWeight_genes.csv
