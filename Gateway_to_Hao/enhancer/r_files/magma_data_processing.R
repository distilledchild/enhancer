library(tidyverse)

# Linux
#setwd('C:\\Users\\panju\\Dropbox (UTHSC GGI)\\Gateway_to_Hao\\workshop\\2023_NIH_meeting\\loop_N_tss')
# setwd('./Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss')
# setwd('~/Desktop/temp/enhancer/dropbox_enhancer_doosan/data/enhancer_atlas2.0/all_species/neuron')
setwd('~/Desktop/temp/enhancer/dropbox_enhancer_doosan')
setwd('/home/pkim/dropbox/Gateway_to_Hao/enhancer/r_files')
source(file.path('/home/pkim/dropbox/Gateway_to_Hao/project_common_code/', 'variables.R'))
source(file.path('/home/pkim/dropbox/Gateway_to_Hao/project_common_code/', 'funcs.R'))

getwd()

# Mac
setwd('/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files')
getwd()
source(file.path('/Users/PanjunKim/dropbox/Gateway_to_Hao/project_common_code/', 'variables.R'))
source(file.path('/Users/PanjunKim/dropbox/Gateway_to_Hao/project_common_code/', 'funcs.R'))

# file reader function
file_data_reader <- function(file) {
  file_ext <- tools::file_ext(file)
  
  if (file_ext == "out" || file_ext == "raw" || file_ext == "tsv") {
    lines <- readLines(file)
    clean_lines <- lines[!grepl("^#", lines)]  # Remove lines starting with #
    
    # Read the file as tab-separated
    read.table(text = clean_lines, header = TRUE, sep = "", stringsAsFactors = FALSE)
  } else {
    stop("Unsupported file type")
  }
}

# Updated function to merge files and add filenames as identifiers
init.file.df <- function(file.list, deli_dir, delim) {
  file.list %>%
    map_dfr(file_data_reader, .id = "file") %>%
    # Extract the significant part of the filename to label the data
    # mutate(file = str_split_n(file, str_c(deli_dir, '\\/'), 2)) %>% 
    mutate(file = str_split_n(file, str_c(deli_dir, '\\/'), 2)) %>%
    mutate(file = str_split_n(file, delim, 1)) 
}

magma.out.file.list = fs::dir_ls("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/sep_HMAGMA/infusion_combined_output", regexp = ".out$")
magma.out.file.list

df.magma.out.task.init <- init.file.df(magma.out.file.list, "infusion_combined_output", '\\.')
df.magma.out.task.init %>% dim()

process_magma_data <- function(df, filter_str, suffix) {
  df %>%
    filter(str_detect(file, filter_str)) %>%
    mutate(file = str_replace(file, filter_str, '')) %>% 
    rename_with(~ paste0(., suffix), -c(file, GENE))
}

df.magma.out.task.init.w.hic <- process_magma_data(df.magma.out.task.init, '_w_hic', '_w_hic')
df.magma.out.task.init.wo.hic <- process_magma_data(df.magma.out.task.init, '_wo_hic', '_wo_hic')

df.magma.out.task.init.w.hic
df.magma.out.task.init.wo.hic


df.magma.out.full.join <- full_join(df.magma.out.task.init.wo.hic, df.magma.out.task.init.w.hic, by = c("file", "GENE"))

# abbrev.
# cm: conventional magma
# hm: hmagma
# na: not available
# sig: significant
# nonsig: non-significant

# 9 categories
# cm\hm            sig.                                         non_sig                               na
# sig.         1. cm_sig_hm_sig. (both_sig)            2. cm_sig_hm_nonsig.                     3. cm_sig_hm_na
# non_sig.     4. cm_nonsig_hm_sig.                    5. cm_nonsig_hm_nonsig. (both_nonsig)    6. cm_nonsig_hm_na
# na.          7. cm_na_hm_sig                         8. cm_na_hm_nonsig                       9. cm_na_hm_na. (both_na)

df.magma.out.full.join.category <- df.magma.out.full.join %>%
  mutate(categ_uniq = case_when(
    is.na(CHR_wo_hic) & !is.na(CHR_w_hic) ~ "uniq_in_w_hic",
    !is.na(CHR_wo_hic) & is.na(CHR_w_hic) ~ "uniq_in_wo_hic",
    TRUE ~ "common"
  )) %>% 
  mutate(category = case_when(
    P_wo_hic < 0.05 & P_w_hic < 0.05 ~ "both_sig", # case 1
    
    P_wo_hic < 0.05 & P_w_hic >= 0.05 ~ "cm_sig_hm_nonsig", # case 2
    P_wo_hic < 0.05 & is.na(P_w_hic) ~ "cm_sig_hm_na", # case 3
    
    P_wo_hic >= 0.05 & P_w_hic < 0.05 ~ "cm_nonsig_hm_sig", # case 4
    P_wo_hic >= 0.05 & P_w_hic >= 0.05 ~ "both_nonsig", # case 5
    P_wo_hic >= 0.05 & is.na(P_w_hic) ~ "cm_nonsig_hm_na", # case 6
    is.na(P_wo_hic) & P_w_hic < 0.05 ~ "cm_na_hm_sig", # case 7
    is.na(P_wo_hic) & P_w_hic >= 0.05 ~ "cm_na_hm_nonsig", # case 8
    # is.na(P_wo_hic) & is.na(P_w_hic) ~ "both_na", # case 9 THIS CANNOT HAPPEN
    
    TRUE ~ NA_character_
  )) 
  
# 114606
df.magma.out.full.join.category
df.magma.out.full.join.category %>% 
  # count(categ_uniq, categ_meaning)
  count(category)

# category     n
#         both_sig  6403 ***********
# cm_sig_hm_nonsig  1435 ?
#     cm_sig_hm_na   458 ?

# cm_nonsig_hm_sig  1434 ***********
#      both_nonsig 85498
#  cm_nonsig_hm_na  4942 ??

#  cm_na_hm_nonsig 13403 ?            : SNPs are added by HiC but the genes are non-significant
#     cm_na_hm_sig  1033 ***********  : SNPs are added by HiC but the genes are SIGNIFICANT


genedef.ncbi.bestrefseq <- read.table("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_hmagma_enhancer/0911/gene_NCBI.txt", 
                           header = FALSE, 
                           sep = "\t", 
                           col.names = c("chr", "length", "genbank", "refseq", "start", "end", "strand", "gene_id", "index"),
                           stringsAsFactors = FALSE)
genedef.ncbi.bestrefseq

df.magma.out.full.join.category.join.genedef.ncbi.bestrefseq <- left_join(df.magma.out.full.join.category, genedef.ncbi.bestrefseq, by = c("GENE" = "index" )) 
df.magma.out.full.join.category.join.genedef.ncbi.bestrefseq

# figure
df.magma.out.full.join.category %>% 
  ggplot(aes(x=-log10(P_w_hic), y=-log10(P_wo_hic)))+geom_point()

#######################
# in-depth analysis
#######################
# 1. ? cases

# 1-1. case 2?: sig -> non significant: 1435
df.magma.out.full.join.category.join.genedef.ncbi.bestrefseq %>% 
  filter(category == 'cm_sig_hm_nonsig') %>% 
  dplyr::select(file, GENE, CHR_wo_hic, gene_id, NSNPS_wo_hic, NSNPS_w_hic) %>% 
  mutate(cnt = NSNPS_w_hic - NSNPS_wo_hic)
# %>% 
  # count(CHR_wo_hic, gene_id, cnt)
  # count(CHR_wo_hic) %>% arrange(CHR_wo_hic)
  # count(file)

#           CHR   n                                                   file   n                 
# 1           1 165                        regressedlr_nicsa_day10_infusion  69                 
# 2          10 181                        regressedlr_nicsa_day11_infusion  71                 
# 3          11  46                         regressedlr_nicsa_day1_infusion  66                 
# 4          12  30                        regressedlr_nicsa_day20_infusion  70                 
# 5          13  25                         regressedlr_nicsa_day2_infusion  91                 
# 6          14  64                         regressedlr_nicsa_day3_infusion  83                 
# 7          15  26                         regressedlr_nicsa_day4_infusion  80                 
# 8          16  54                         regressedlr_nicsa_day5_infusion  90                 
# 9          17  18                         regressedlr_nicsa_day6_infusion  72                 
# 10         18  17                         regressedlr_nicsa_day7_infusion  71                 
# 11         19  32                         regressedlr_nicsa_day8_infusion  74                 
# 12          2 105                         regressedlr_nicsa_day9_infusion  81                 
# 13         20  62      regressedlr_nicsa_first_three_days_infusion_median  84                 
# 14          3  98       regressedlr_nicsa_first_three_days_infusion_total  84                 
# 15          4 117       regressedlr_nicsa_last_three_days_infusion_median  87                 
# 16          5 108        regressedlr_nicsa_last_three_days_infusion_total  88                 
# 17          6  64               regressedlr_nicsa_slope_nicotine_infusion  65                 
# 18          7  62                 regressedlr_nicsa_total_infusion_10days 109                 
# 19          8  83                       
# 20          9  78                       

# 1-2. case 3?: sig -> NA : 458
df.magma.out.full.join.category.join.genedef.ncbi.bestrefseq %>% 
  filter(category == 'cm_sig_hm_na') %>% 
  dplyr::select(file, GENE, CHR_wo_hic, gene_id, NSNPS_wo_hic, NSNPS_w_hic) %>% 
  # count(CHR_wo_hic) %>% arrange(CHR_wo_hic)
  count(file)
  # count(NSNPS_w_hic) # NA 458

#    CHR_wo_hic   n                                                            file  n  
# 1           1  57                              regressedlr_nicsa_day10_infusion 17  
# 2          10  34                              regressedlr_nicsa_day11_infusion 30  
# 3          11   2                               regressedlr_nicsa_day1_infusion 41  
# 4          12  13                              regressedlr_nicsa_day20_infusion 26  
# 5          13   6                               regressedlr_nicsa_day2_infusion 26  
# 6          14   6                               regressedlr_nicsa_day3_infusion 18  
# 7          15   3                               regressedlr_nicsa_day4_infusion 21  
# 8          16   5                               regressedlr_nicsa_day5_infusion 17  
# 9          17   8                               regressedlr_nicsa_day6_infusion 25  
# 10         18   2                               regressedlr_nicsa_day7_infusion 32  
# 11         19  10                               regressedlr_nicsa_day8_infusion 15  
# 12          2   8                               regressedlr_nicsa_day9_infusion 18  
# 13          3  11            regressedlr_nicsa_first_three_days_infusion_median 18  
# 14          4  23             regressedlr_nicsa_first_three_days_infusion_total 33  
# 15          5  24             regressedlr_nicsa_last_three_days_infusion_median 28  
# 16          6  27              regressedlr_nicsa_last_three_days_infusion_total 35  
# 17          7  35                     regressedlr_nicsa_slope_nicotine_infusion 29  
# 18          8  19                       regressedlr_nicsa_total_infusion_10days 29  
# 19          9  12             
# 20          X 153             

# 1-3. case 6?: non sig -> NA : 4942
df.magma.out.full.join.category.join.genedef.ncbi.bestrefseq %>% 
  filter(category == 'cm_nonsig_hm_na') %>% 
  dplyr::select(file, GENE, CHR_wo_hic, gene_id, NSNPS_wo_hic, NSNPS_w_hic) %>% 
  # count(CHR_wo_hic) %>% arrange(CHR_wo_hic)    ##################################### DEPENDING ON CHR
  count(file)
  # count(NSNPS_w_hic) # NA 4942

#    CHR_wo_hic    n.                                                          file   n
# 1           1  429.          1                    regressedlr_nicsa_day10_infusion 283
# 2          10  200.          2                    regressedlr_nicsa_day11_infusion 270
# 3          11   34.          3                     regressedlr_nicsa_day1_infusion 259
# 4          12   95.          4                    regressedlr_nicsa_day20_infusion 274
# 5          13  102.          5                     regressedlr_nicsa_day2_infusion 274
# 6          14   48.          6                     regressedlr_nicsa_day3_infusion 282
# 7          15   69.          7                     regressedlr_nicsa_day4_infusion 279
# 8          16  121.          8                     regressedlr_nicsa_day5_infusion 283
# 9          17   82.          9                     regressedlr_nicsa_day6_infusion 275
# 10         18   52.          10                    regressedlr_nicsa_day7_infusion 268
# 11         19  116.          11                    regressedlr_nicsa_day8_infusion 285
# 12          2  244.          12                    regressedlr_nicsa_day9_infusion 282
# 13         20   18.          13 regressedlr_nicsa_first_three_days_infusion_median 282
# 14          3  277.          14  regressedlr_nicsa_first_three_days_infusion_total 267
# 15          4  157.          15  regressedlr_nicsa_last_three_days_infusion_median 272
# 16          5  174.          16   regressedlr_nicsa_last_three_days_infusion_total 265
# 17          6  207.          17          regressedlr_nicsa_slope_nicotine_infusion 271
# 18          7  253.          18            regressedlr_nicsa_total_infusion_10days 271
# 19          8  143.          
# 20          9   96.          
# 21          X 2025.          

# 2. Extracting genes significant from HMAGMA
# 2-1. geneset from both_sig
df.extracting.geneset.from.both.sig <- df.magma.out.full.join.category.join.genedef.ncbi.bestrefseq %>% 
  filter(category == 'both_sig') %>%
  dplyr::select(file, GENE, CHR_wo_hic, gene_id, NSNPS_wo_hic, NSNPS_w_hic) %>%
  mutate(cnt = NSNPS_w_hic - NSNPS_wo_hic)

list.genes.both.sig <- split(df.extracting.geneset.from.both.sig$gene_id, df.extracting.geneset.from.both.sig$file)
list.genes.both.sig


common.genes.both.sig <- Reduce(intersect, list.genes.both.sig)
common.genes.both.sig # NO GENE common in ALL over 18 phenotypes

unique_genes <- lapply(list.genes.both.sig, function(genes) setdiff(genes, common.genes.both.sig))

# result for genesets common in phenotypes
cat("genes in common:", common.genes.both.sig, "\n")
cat("genes specific depending on phenotypes:\n")
print(unique_genes)

# input for gProfiler
gene.id.both.sig <- df.extracting.geneset.from.both.sig %>% 
  distinct(gene_id) %>% 
  pull(gene_id) 

# %>% 
#   str_c(collapse = " ") %>% 
#   cat()

gostres.both.sig <- gost(query = gene.id.both.sig, organism = "rnorvegicus")

head(gostres.both.sig$result)
p.both.sig <- gostplot(gostres.both.sig, capped = FALSE, interactive = FALSE)
p.both.sig

ggsave("./magma_data_processing/0918/p.both.sig.png", 
       plot = p.both.sig, 
       width = 10, height = 8, dpi = 300)


# 2-2. geneset from cm_nonsig_hm_sig

df.extracting.geneset.from.cm.nonsig.hm.sig <- df.magma.out.full.join.category.join.genedef.ncbi.bestrefseq %>% 
  filter(category == 'cm_nonsig_hm_sig') %>%
  dplyr::select(file, GENE, CHR_wo_hic, gene_id, NSNPS_wo_hic, NSNPS_w_hic) %>%
  mutate(cnt = NSNPS_w_hic - NSNPS_wo_hic)

list.genes.cm.nonsig.hm.sig <- split(df.extracting.geneset.from.cm.nonsig.hm.sig$gene_id, df.extracting.geneset.from.cm.nonsig.hm.sig$file)
list.genes.cm.nonsig.hm.sig


common.genes.cm.nonsig.hm.sig <- Reduce(intersect, list.genes.cm.nonsig.hm.sig)
common.genes.cm.nonsig.hm.sig # NO GENE common in ALL over 18 phenotypes

unique.genes.cm.nonsig.hm.sig <- lapply(list.genes.cm.nonsig.hm.sig, function(genes) setdiff(genes, common.genes.cm.nonsig.hm.sig))

# result for genesets common in phenotypes
cat("genes in common:", common.genes.cm.nonsig.hm.sig, "\n")
cat("genes specific depending on phenotypes:\n")
print(unique.genes.cm.nonsig.hm.sig)

# input for gProfiler
gene.id.cm.nonsig.hm.sig <- df.extracting.geneset.from.cm.nonsig.hm.sig %>% 
  distinct(gene_id) %>% 
  pull(gene_id) 
# %>% 
#   str_c(collapse = " ") %>% 
#   cat()

gostres.gene.id.cm.nonsig.hm.sig <- gost(query = gene.id.cm.nonsig.hm.sig, organism = "rnorvegicus")

head(gostres.gene.id.cm.nonsig.hm.sig$result)
p.gene.id.cm.nonsig.hm.sig <- gostplot(gostres.gene.id.cm.nonsig.hm.sig, capped = FALSE, interactive = FALSE)
p.gene.id.cm.nonsig.hm.sig

ggsave("./magma_data_processing/0918/p.gene.id.cm.nonsig.hm.sig.png", 
       plot = p.gene.id.cm.nonsig.hm.sig, 
       width = 10, height = 8, dpi = 300)

# 2-3. geneset from cm_na_hm_sig

df.extracting.geneset.from.cm.na.hm.sig <- df.magma.out.full.join.category.join.genedef.ncbi.bestrefseq %>% 
  filter(category == 'cm_na_hm_sig') %>%
  dplyr::select(file, GENE, CHR_wo_hic, gene_id, NSNPS_wo_hic, NSNPS_w_hic) %>%
  mutate(cnt = NSNPS_w_hic - NSNPS_wo_hic)

list.genes.cm.na.hm.sig <- split(df.extracting.geneset.from.cm.na.hm.sig$gene_id, df.extracting.geneset.from.cm.na.hm.sig$file)
list.genes.cm.na.hm.sig


common.genes.cm.na.hm.sig <- Reduce(intersect, list.genes.cm.na.hm.sig)
common.genes.cm.na.hm.sig # NO GENE common in ALL over 18 phenotypes

unique.genes.cm.na.hm.sig <- lapply(list.genes.cm.na.hm.sig, function(genes) setdiff(genes, common.genes.cm.na.hm.sig))

# result for genesets common in phenotypes
cat("genes in common:", common.genes.cm.na.hm.sig, "\n")
cat("genes specific depending on phenotypes:\n")
print(unique.genes.cm.na.hm.sig)

# input for gProfiler
gene.id.cm.na.hm.sig <- df.extracting.geneset.from.cm.na.hm.sig %>% 
  distinct(gene_id) %>% 
  pull(gene_id) 
# %>%
#   str_c(collapse = " ") %>% 
#   cat()

dfasdfafasfas
library(gprofiler2)

gostres.cm.na.hm.sig <- gost(query = gene.id.cm.na.hm.sig, organism = "rnorvegicus")

head(gostres$result)
p.cm.na.hm.sig <- gostplot(gostres.cm.na.hm.sig, capped = FALSE, interactive = FALSE)
p.cm.na.hm.sig

ggsave("./magma_data_processing/0918/p.gene.cm.na.hm.sig.png", 
       plot = p.cm.na.hm.sig, 
       width = 10, height = 8, dpi = 300)
