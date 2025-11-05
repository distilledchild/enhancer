library("ComplexUpset")
library("tidyverse")
library("GenomicRanges")
library("ggplot2")
library("patchwork")
library("devtools")
library("remotes")
library("rvest")
library("xml2")
library("openxlsx")
library("Vennerable")
library("scales") # figure4
options(stingAsFactors=F)

# Windows
# setwd('C:\\Users\\panju\\Dropbox (UTHSC GGI)\\Gateway_to_Hao\\workshop\\2023_NIH_meeting\\loop_N_tss')
# source(file.path('C:\\Users\\panju\\Dropbox (UTHSC GGI)\\Gateway_to_Hao\\project_common_code', 'variables.R'))
# source(file.path('C:\\Users\\panju\\Dropbox (UTHSC GGI)\\Gateway_to_Hao\\project_common_code', 'funcs.R'))

# Linux
#setwd('C:\\Users\\panju\\Dropbox (UTHSC GGI)\\Gateway_to_Hao\\workshop\\2023_NIH_meeting\\loop_N_tss')
# setwd('./Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss')
# setwd('~/Desktop/temp/enhancer/dropbox_enhancer_doosan/data/enhancer_atlas2.0/all_species/neuron')
setwd('~/Desktop/temp/enhancer/dropbox_enhancer_doosan')
source(file.path('/home/panjun/Desktop/temp/enhancer/dropbox_enhancer_doosan/', 'variables.R'))
source(file.path('/home/panjun/Desktop/temp/enhancer/dropbox_enhancer_doosan/', 'funcs.R'))

getwd()

# Linux (local of local)
# setwd('/home/panjun/dropbox/Gateway_to_Hao/enhancer/')
source(file.path('/home/panjun/dropbox/Gateway_to_Hao/project_common_code/', 'variables.R'))
source(file.path('/home/panjun/dropbox/Gateway_to_Hao/project_common_code/', 'funcs.R'))

# Mac
setwd('/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files')
source(file.path('/Users/PanjunKim/dropbox/Gateway_to_Hao/project_common_code/', 'variables.R'))
source(file.path('/Users/PanjunKim/dropbox/Gateway_to_Hao/project_common_code/', 'funcs.R'))
getwd()

#################################################
# 1. Exploratory Data analysis (EDA)
# 1.1. chr meta data preprocessing
#################################################
### 1. meta data preprocessing

url <- "https://www.ncbi.nlm.nih.gov/grc/rat/data?asm=mRatBN7.2"
page <- read_html(url)
tables <- page %>% html_table(fill = TRUE)
chromosome.table.1 <- tables[[1]]
chromosome.length <- chromosome.table.1 %>% 
  # %>% head() 
  mutate(`Total length (bp)` = as.numeric(gsub(",", "", `Total length (bp)`))) %>% 
  dplyr::rename(chr = Chromosome, length = `Total length (bp)`, genbank = `GenBank accession`, refseq = `RefSeq accession`)  

chromosome.length %>% colnames() #chr      length genbank    refseq    
chromosome.length
### 2. NCBI refseq data preprocessing
## 2-1. BestRefSeq: NCBI의 RefSeq 프로젝트에서 제공하는, 가장 신뢰할 수 있는 어노테이션을 가진 서열을 나타냅니다. 이는 엄격한 검증 과정을 거쳐 선택된 서열로, 고품질의 유전체 어노테이션을 제공합니다.
## 2-2. BestRefSeq,Gnomon: 이 값은 BestRefSeq 어노테이션과 Gnomon 어노테이션 방법을 모두 사용하여 어노테이션된 서열을 나타냅니다. Gnomon은 NCBI에서 개발한 예측 기반 어노테이션 도구입니다.
## 2-3. Curated Genomic: 전문가들에 의해 수동으로 큐레이션된 유전체 서열을 의미합니다. 이 데이터는 공개된 연구 결과나 실험적 검증을 통해 얻어진 정보를 기반으로 합니다.
## 2-4. Gnomon: NCBI의 자체 개발 어노테이션 도구로, 주로 유전자 예측에 사용됩니다. Gnomon은 컴퓨터 알고리즘을 사용하여 유전자의 위치, 구조, 기능을 예측합니다.
## 2-5. RefSeq: NCBI의 RefSeq 프로젝트에서 제공하는 참조 서열 데이터베이스를 의미합니다. RefSeq는 유전자, 단백질, 유전체 서열 등에 대한 포괄적인 어노테이션 정보를 제공합니다.
## 2-6. cmsearch: 이는 "covariance model search"의 약자로, 주로 리보핵산(RNA) 구조를 기반으로 한 서열 검색에 사용되는 도구입니다. cmsearch는 특정 RNA 구조에 대한 서열의 존재 여부를 확인하는 데 사용됩니다.
## 2-7. tRNAscan-SE: tRNA 유전자를 찾기 위해 특별히 설계된 컴퓨터 프로그램입니다. 이 도구는 유전체 데이터에서 tRNA 유전자를 식별하고 어노테이션하는 데 널리 사용됩니다.
## BestRefSeq (chosen) cf. RefSeq

#################################################
# 1. Exploratory Data analysis (EDA)
# 1.2. GTF data preprocessing
#################################################
# GTF file path
file_path <- "/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/data/gene_data/ncbi_dataset/ncbi_dataset/data/GCF_015227675.2/genomic_BestRefSeq.gtf"

df.genomic.BestRefSeq.origin <- read_delim(file_path, 
                                           delim = "\t", 
                                           col_names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"),
                                           quote = "",
                                           col_types = cols(
                                             start = col_character(),   
                                             end = col_character(),     
                                           ),
                                           trim_ws = TRUE)

df.genomic.BestRefSeq.origin %>% count(feature)
df.genomic.BestRefSeq.origin # 468,661 × 9 (PASS)

# from GCF_015227675.2
#1 CDS         192169
#2 exon        205520 *************
#3 gene          7539 *************
#4 start_codon  20613 *************
#5 stop_codon   20590 *************
#6 transcript   22230 *************

df.genomic.BestRefSeq <- df.genomic.BestRefSeq.origin %>% 
  filter(str_detect(feature, "exon|transcript")) %>%
  left_join(chromosome.length, by = c("seqname" = "refseq")) %>% 
  mutate(index = str_c(chr, start, end, sep = ':')) %>% 
  filter(!is.na(chr))

df.genomic.BestRefSeq # 227,677
df.genomic.BestRefSeq %>% count(feature)
# feature         n
# <chr>       <int>
# 1 exon       205460
# 2 transcript  22217

##################################################
# checking dups genes: START
##################################################
# Cers1 and GDF1 are produced from the same bicistronic gene, but from non-overlapping reading frames.
# same index, but different gene_id
df.exon.transcript.genomic.BestRefSeq <- df.genomic.BestRefSeq %>% 
  group_by(index) %>% 
  filter(n() > 1) %>% 
  ungroup() 

df.exon.transcript.genomic.BestRefSeq %>% 
  filter(str_detect(attributes, 'Cers1|Gdf1')) %>% 
  select(-c(source, seqname, score, frame, attributes))
##################################################
# checking dups genes: END
##################################################
# FINAL DATASET
df.genomic.BestRefSeq <- df.genomic.BestRefSeq %>%
  filter(!str_detect(attributes, 'Cers1'))
df.genomic.BestRefSeq # 227,668
df.genomic.BestRefSeq %>% count(feature)
# feature         n
# <chr>       <int>
# 1 exon       205452
# 2 transcript  22216

#################################################
# 1. Exploratory Data analysis (EDA)
# 1.3. exon data preprocessing
#################################################
#########
## EXAMPLE exons
#########
# Required format from Won
the_won_exon <- read.table("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/actual_HMAGMA/Input_Files/Gencode26_exon.bed")
the_won_exon

#########
# EXON : extracted from features exon, CDS, exon, gene, start, stop, transcript
#########
# line 104
df.genomic.BestRefSeq
df.genomic.BestRefSeq %>% head()

df.exon.genomic.BestRefSeq.extended <- df.genomic.BestRefSeq %>%
  filter(feature == "exon") %>%
  mutate(
    exon_number   = str_remove_all(str_split_n(str_split_n(attributes, 'exon_number', 2), ';', 1), '\\\\|\\"') %>% str_trim(),
    gene_id       = str_remove_all(str_split_n(str_split_n(attributes, 'gene_id', 2), ';', 1), '\\\\|\\"') %>% str_trim(),
    tag           = str_match(attributes, 'tag \\"(.*?)\\"')[, 2] %>% str_trim(),
    transcript_id = str_remove_all(str_split_n(str_split_n(attributes, 'transcript_id', 2), ';', 1), '\\\\|\\"') %>% str_trim()
  ) %>%
  dplyr::select(feature, chr, start, end, score, strand, exon_number, gene_id, tag, transcript_id, index)

df.exon.genomic.BestRefSeq.extended # 205,452 × 10
df.exon.genomic.BestRefSeq.extended %>% distinct(gene_id) # 18,621

##################################################
# checking dups exon: START
##################################################
# df.duplicated.exon.genomic.BestRefSeq.extended <- df.exon.genomic.BestRefSeq.extended %>%
#   group_by(chr, start, end, strand, gene_id, exon_number) %>%
#   filter(n() > 1) %>%
#   ungroup()
# 
# df.duplicated.exon.genomic.BestRefSeq.extended
# 
# df.exon.genomic.BestRefSeq.extended %>% 
#   distinct(chr, start, end, strand, gene_id, exon_number)
# 
# exon.shared <- df.exon.genomic.BestRefSeq.extended %>%
#   group_by(chr, start, end) %>%
#   summarise(n_genes = n_distinct(gene_id), .groups = "drop") %>%
#   filter(n_genes > 1)
# 
# shared.exons.with.genes <- df.exon.genomic.BestRefSeq.extended %>%
#   inner_join(exon.shared, by = c("chr", "start", "end"))
# 
# shared.exons.with.genes %>% print(n = Inf)
##################################################
# checking dups exon: END
##################################################
exon <- df.exon.genomic.BestRefSeq.extended
exon # 205,452

exonranges <- GRanges(
  seqnames = exon$chr, 
  ranges = IRanges(start = as.numeric(exon$start), end = as.numeric(exon$end)), 
  gene = exon$gene_id
)
exonranges # 205452

#################################################
# 1. Exploratory Data analysis (EDA)
# 1.4. promoter data preprocessing
#################################################
#########
# PROMOTER is extracting from transcript among CDS, exon, gene, start, stop, transcript
#########
df.transcript.genomic.BestRefSeq.extended <- df.genomic.BestRefSeq %>%
  filter(feature == "transcript") %>%
  mutate(length = as.numeric(end) - as.numeric(start) + 1) %>% 
  mutate(
    gene_id       = str_remove_all(str_split_n(str_split_n(attributes, 'gene_id', 2), ';', 1), '\\\\|\\"') %>% str_trim(),
    tag           = str_match(attributes, 'tag \\"(.*?)\\"')[, 2] %>% str_trim(),
    transcript_id = str_remove_all(str_split_n(str_split_n(attributes, 'transcript_id', 2), ';', 1), '\\\\|\\"') %>% str_trim()
  ) %>% 
  dplyr::select(feature, chr, start, end, score, strand, gene_id, tag, transcript_id, length, index) 

df.transcript.genomic.BestRefSeq.extended # 22,216
df.transcript.genomic.BestRefSeq.extended %>% distinct(gene_id) # 18,621 (promoter) / 18,621 (exon)

#########################################
# preprocessing clearance with tag column : START
#########################################
# STEP 1: unique transcript for a single gene_id
df.transcript.unique.by.gene_id.1 <- df.transcript.genomic.BestRefSeq.extended %>%
  group_by(gene_id, .drop = FALSE) %>%
  filter(n() == 1) %>%
  ungroup() %>% 
  mutate(tag = "RefSeq Select")
df.transcript.unique.by.gene_id.1 # 16,354
df.transcript.unique.by.gene_id.1 %>% count(gene_id) # 16,354

df.transcript.genomic.BestRefSeq.extended %>% 
  count(tag)

# STEP2: dups transcripts for a single gene_id with tag = "RefSeq Select"
df.transcript.unique.by.gene_id.2 <- df.transcript.genomic.BestRefSeq.extended %>%
  group_by(gene_id, .drop = FALSE) %>%
  filter(n() > 1) %>% 
  ungroup() %>% 
  filter(!is.na(tag))

df.transcript.unique.by.gene_id.2 # A tibble: 1,785 × 11
df.transcript.unique.by.gene_id.2 %>% distinct(gene_id) # 1,785

transcript_gene_ids2 <- df.transcript.unique.by.gene_id.2 %>% pull(gene_id) %>% unique()
transcript_gene_ids2 # 1,785

# STEP3: dups transcripts for a single gene_id with tag = NA
df.transcript.unique.by.gene_id.3 <- df.transcript.genomic.BestRefSeq.extended %>%
  group_by(gene_id, .drop = FALSE) %>%
  filter(n() > 1) %>% 
  ungroup() %>% 
  filter(!(gene_id %in% transcript_gene_ids2)) %>%
  group_by(gene_id) %>% 
  slice_max(order_by = length, n = 1, with_ties = FALSE) %>%
  ungroup()

df.transcript.unique.by.gene_id.3 # 482
df.transcript.unique.by.gene_id.3 %>% distinct(gene_id) # 482 : PASS
df.transcript.unique.by.gene_id.3 %>% count(tag)

df.transcript.genomic.BestRefSeq.extended  %>% distinct(gene_id) # 18,621 genes

df.transcript.distinct.gene_id.genomic.BestRefSeq.extended <- bind_rows(
  df.transcript.unique.by.gene_id.1, # 16354 genes
  df.transcript.unique.by.gene_id.2, # 1785 genes
  df.transcript.unique.by.gene_id.3 # 482 genes
)
#########################################
# preprocessing clearance with tag column : END
#########################################
df.transcript.distinct.gene_id.genomic.BestRefSeq.extended %>% count(gene_id) # 18,621
# > 16355 + 1785 + 482 = 18622 : PASS

########################################
# 1.4.1. checking gene set comparison: START
########################################
# 1. GENES: 18622 from transcript
df.transcript.genomic.BestRefSeq.extended # A tibble: 22,216 × 11
df.transcript.genomic.BestRefSeq.extended %>% distinct(gene_id) # 18,621

# by Refseq_Select
df.transcript.distinct.gene_id.genomic.BestRefSeq.extended # A tibble: 18,621 × 18
df.transcript.distinct.gene_id.genomic.BestRefSeq.extended %>% distinct(gene_id) # 18,621

# 3. comparison
gene_set_exon <- df.exon.genomic.BestRefSeq.extended %>% distinct(gene_id) %>% pull(gene_id)
gene_set_transcript <- df.transcript.genomic.BestRefSeq.extended %>% distinct(gene_id) %>% pull(gene_id)

gene_set_exon # 18621
gene_set_transcript # 18621

# genes from exon vs transcript
length(intersect(gene_set_exon, gene_set_transcript)) # 18621, 100% overlap

setequal(gene_set_transcript, intersect(gene_set_transcript, gene_set_exon)) # TRUE  # promoter가 exon에 포함되는가
all(gene_set_transcript %in% gene_set_exon)  # TRUE                               # TRUE면 promoter ⊆ exon

# 교집합
intersect.genes <- intersect(gene_set_exon, gene_set_transcript) # TRUE
length(intersect.genes) # TRUE

# exon에만 있는 gene
only.in.exon <- setdiff(gene_set_exon, gene_set_transcript)
only.in.exon # PASS: 0
length(only.in.exon) # PASS: 0
########################################
# 1.4.1. checking gene set comparison: END
########################################
##############
# promoter case: START
##############
df.transcript.genomic.BestRefSeq.extended # 22,216
df.transcript.distinct.gene_id.genomic.BestRefSeq.extended # 18,621

# promoter case1
promoter <- df.transcript.genomic.BestRefSeq.extended %>%
  mutate(
    start_new = ifelse(strand == "+", as.numeric(start) - 2000, end), # promoter start
    end_new   = ifelse(strand == "+", start, as.numeric(end) + 2000), # promoter end
    promoter_length = abs(as.numeric(start_new) - as.numeric(end_new)) + 1,
  )

# promoter case2
# promoter <- df.transcript.distinct.gene_id.genomic.BestRefSeq.extended %>%
#   mutate(
#     start_new = ifelse(strand == "+", as.numeric(start) - 2000, end), # promoter start
#     end_new   = ifelse(strand == "+", start, as.numeric(end) + 2000), # promoter end
#     promoter_length = abs(as.numeric(start_new) - as.numeric(end_new)) + 1,
#   )

promoter # 22,216|18,621
promoter %>% filter(start_new < 0) # PASS: 0
promoter %>% filter(end_new < 0) # PASS: 0
promoter %>% count(promoter_length) # 2001 22,216|18621
promoter %>% count(strand)
##############
# promoter case: END
##############
#########
## EXAMPLE promoter
#########
# Required format from Won
the_won_promoter <- read.table("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/actual_HMAGMA/Input_Files/Gencode26_promoter.bed")
the_won_promoter

promoterranges <- GRanges(
  seqnames = promoter$chr,
  ranges = IRanges(start = as.numeric(promoter$start_new), end = as.numeric(promoter$end_new)),
  # strand = promoter$strand,
  gene = promoter$gene_id
)
mcols(promoterranges)$promoter_id <- paste0(
  promoter$gene_id, "_",
  promoter$chr, ":", promoter$start_new, "-", promoter$end_new
)

promoterranges
exonranges

# save(exonranges, promoterranges, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_hmagma_enhancer/0911/exon_promoranges.rda")
# save(exonranges, promoterranges, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/exon_promoranges.rda")
# load ("from_hmagma_enhancer/0911/exon_promoranges.rda")

#################################################
# 2. Exploratory Data analysis (EDA)
# 2.1. SNP data preprocessing
# Generate a GenomicRanges object for the SNP annotation
#################################################
# required format comparison from Won
# EUR.bim
# 1	rs367896724	0	10177	AC	A
# 1	rs555500075	0	10352	TA	T

# HS_genotypes_v4.bim
# 1	1:22585	0	22585	A	C
# 1	1:62340	0	62340	T	G
# 1	1:120520	0	120520	G	A

snps <- read.table("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/data/hs_data/v4/HS_genotypes_v4.bim")
dim(snps) # 7358643
snps %>% head(10)
snps <- snps[,c(1,2,4)] # chr, ID, base coordinate
snps %>% head(10)
colnames(snps) <- c("chr","rsid","Position")
snps %>% head()
# snps$chr <- sub("^", "chr", snps$chr)
snps

snps.g <- GRanges(snps$chr, IRanges(snps$Position, snps$Position), rsid=snps$rsid)
snps.g # 7358643
#################################################
# 3. Asigning SNP
# 3.1. Overlap exons with SNPs
#################################################
# SNP & exon
olap.snp.exon <- findOverlaps(snps.g, exonranges)

olap.snp.exon # 87118
snpexon <- snps.g[queryHits(olap.snp.exon)]
snpexon # 87118

mcols(exonranges[subjectHits(olap.snp.exon)]) # 87118
mcols(snpexon) <- cbind(mcols(snpexon), mcols(exonranges[subjectHits(olap.snp.exon)]))
snpexon  # 87118: SNPs in exons of genes
snpexon %>% as.data.frame() %>% distinct(gene) %>% nrow() ################# 14460 genes

# Data in X chromosome 
#snpexon <- snpexon[seqnames(snpexon)!="X"] # 86239
snpexon # 87118
snpexon %>% as.data.frame() %>% distinct(gene) %>% nrow() # gene: 14460(Cers1 only)

#################################################
# 3. Asigning SNP
# 3.2. Overlap promoters with SNPs
#################################################
# SNP & promoter
olap.snp.pro <- findOverlaps(snps.g, promoterranges)

olap.snp.pro # 101325 case1/87210 case2
snpro <- snps.g[queryHits(olap.snp.pro)];
snpro # 101325 case1/87210 case2

mcols(promoterranges[subjectHits(olap.snp.pro)])
mcols(snpro) <- cbind(mcols(snpro), mcols(promoterranges[subjectHits(olap.snp.pro)]))

# Data in X chromosome 
#snpro <- snpro[seqnames(snpro)!="X"]
snpro # 101325 case1 | 86204 case2 + excluding X
snpro %>% as.data.frame() %>% distinct(gene) %>% nrow() # gene: 16121(case1 + INCLUDING X) | 15727/16089 15736/16093

snpro # 101325 case1| 86204 case2

save(snpro, snpexon, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_hmagma_enhancer/0911/snp_locating_in_exon_promoter_transcript_level.rda")
save(snpro, snpexon, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/snp_locating_in_exon_promoter_transcript_level.rda")

#load("snp_locating_in_exon_promoter_transcript_level.rda")

#################################################
# 3. Annotation file for cMAGMA
#################################################
# loading SNP in exons and promoter
# load("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_hmagma_enhancer/0911/snp_locating_in_exon_promoter_transcript_level.rda")

# snpexon snpro GR converting to df.
df.snpexon <- data.frame(rsid = snpexon$rsid, ensg = snpexon$gene)
df.snpexon # 87118/86239(incl X/excl X)
df.snpro <- data.frame(rsid = snpro$rsid, ensg = snpro$gene)
df.snpro # 101325/86204(case1 incl X/ case2 excl X)

# bind rows
bind_rows(df.snpexon %>% mutate(cate = "exon"), df.snpro %>% mutate(cate = "pro")) %>% 
  group_by(rsid, ensg) %>%
  filter(n() > 1) 

bind_rows(df.snpexon, df.snpro) %>% 
  dim() # 188443
df.distinct.snp <- bind_rows(df.snpexon, df.snpro) %>%
  distinct() 

df.distinct.snp # 164629 (case1 promoter, incl X in exon & promoter) | sub.4: 161029 (no chrX) rsid     ensg
df.distinct.snp %>% count(ensg) # 17239/16808(case1 promoter, incl X in exon & promoter / case2 promoter, excl X in exon & promoter)

###############################################
# gene info from genedef
###############################################
genes.all <- df.transcript.distinct.gene_id.genomic.BestRefSeq.extended %>% 
  dplyr::select(gene_id, chr, start, end, strand, index)
genes.all

# must assing chr:start:end values into gene column to make format of annotation
df.distinct.snp # 164629 (case1 promoter, incl X in exon & promoter) | sub.4: 161029 (no chrX), rsid     ensg
df.distinct.snp.merged <- df.distinct.snp %>%
  left_join(genes.all, by = c("ensg" = "gene_id")) %>%
  filter(!is.na(index))
df.distinct.snp.merged # 164629/161029 (case1 promoter, incl X in exon & promoter) | sub.4 (no chrX)
df.distinct.snp.merged %>% dim()
df.distinct.snp.merged %>%
  count(ensg) # 16739 + 500 = 17239 | 16808

# 3. snpaggconv
snpaggconv.cmagma.annot <- df.distinct.snp.merged %>%
  group_by(index) %>%
  summarise(rsid = paste(unique(rsid), collapse = ", ")) %>%
  mutate(ensg = index) %>%
  dplyr::select(ensg, index, rsid)

snpaggconv.cmagma.annot #17,239/16,808 (case1 promoter, incl X in exon & promoter)/sub.4 (case2 no chrX)
snpaggconv.cmagma.annot %>% head()

write.table(snpaggconv.cmagma.annot, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/snpagg_exon_promoter_only.txt", 
            quote=F, row.names=F, col.names=FALSE, sep="\t")
system("sed -e 's/, /\t/g' < /Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/snpagg_exon_promoter_only.txt > /Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/hs.exonpro.ONLY.annot.submission", wait = TRUE)

###################################################################################################
#####################END:::::::::::::::::::::::::annotation for cMAGMA ############################
###################################################################################################

#################################################
# Assign SNPs to genes using Hi-C.
#################################################

snps # line 229, data read from .bim,  7358643
snpexon # 87118
snpro # 101325
snpexon$rsid
snps$rsid

# SNPs NOT in exon area
snpranges <- snps[!(snps$rsid %in% snpexon$rsid), ]  
snpranges # 7282868 +334 = 7283202

snpranges <- snpranges[!(snpranges$rsid %in% snpro$rsid), ] # granges of snps NOT in promoter area
snpranges # 7196971 + 335 = 7197306
#snpranges: snps NOT in exon & promo -> snpranges: non_exonic_promoter_snp.rda
save(snpranges,file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/non_exonic_promoter_snp.rda")

# load("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/non_exonic_promoter_snp.rda") 
# load("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/exon_promoranges.rda")

#####################
# required Hi-C format from WON
# chrom1 start1 end1 chrom2 start2 end2
# chr10 100080000 100090000 chr10 100920000 100940000
# chr10 100080000 100090000 chr10 101020000 101040000

# load loops filtered
load("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/df_final_loop.rda")

df.final.loop # 11922
df.final.loop %>% count(str_split_n(loop.id, '_', 7))

df.formatted.final.loop.anchor1 <- df.final.loop %>%
  separate(loop.id, into = c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "resolution"), sep = "_", remove = FALSE) %>%
  mutate(chrom1 = str_remove(chrom1, "chr"),
         chrom2 = str_remove(chrom2, "chr"),
         loop.id = str_remove_all(loop.id, "chr"))
df.formatted.final.loop.anchor1

df.formatted.final.loop.anchor2 <- df.formatted.final.loop.anchor1 %>% 
  dplyr::select("loop.id", "chrom2", "start2", "end2", "chrom1", "start1", "end1", "resolution") %>% 
  dplyr::rename(chrom1 = chrom2, start1 = start2, end1 = end2, 
                chrom2 = chrom1, start2 = start1, end2 = end1) 
# %>% 
#   mutate(loop.id = str_c(chrom1, '_', start1, '_', end1, '_', chrom2, '_', start2, '_', end2, '_', resolution))
df.formatted.final.loop.anchor2

df.formatted.final.loop <- bind_rows(df.formatted.final.loop.anchor1, df.formatted.final.loop.anchor2) %>% 
  mutate(across(c(start1, end1, start2, end2), as.numeric))
df.formatted.final.loop
df.formatted.final.loop %>% count()

# 5k
df.formatted.final.loop.5k <- bind_rows(df.formatted.final.loop.anchor1, df.formatted.final.loop.anchor2) %>%
  filter(resolution == 5000) # 5k: 4442
df.formatted.final.loop <- df.formatted.final.loop.5k
df.formatted.final.loop

# 10k
df.formatted.final.loop.10k <- bind_rows(df.formatted.final.loop.anchor1, df.formatted.final.loop.anchor2) %>%
  filter(resolution == 10000) # 10k: 8844
df.formatted.final.loop <- df.formatted.final.loop.10k
df.formatted.final.loop # 8719 + 125 = 8844

# 25k
df.formatted.final.loop.25k <- bind_rows(df.formatted.final.loop.anchor1, df.formatted.final.loop.anchor2) %>%
  filter(resolution == 25000) # 5k10k: 13286
df.formatted.final.loop <- df.formatted.final.loop.25k
df.formatted.final.loop # 10433 + 125 = 10558

# 5k10k
df.formatted.final.loop.5k.10k <- bind_rows(df.formatted.final.loop.anchor1, df.formatted.final.loop.anchor2) %>%
  filter(resolution != 25000) # 5k10k: 13286
df.formatted.final.loop <- df.formatted.final.loop.5k.10k


df.formatted.final.loop
df.formatted.final.loop %>% count(loop.id) %>% filter(n > 1)

hicranges <- GRanges(df.formatted.final.loop$chrom1, 
                     IRanges(df.formatted.final.loop$start1, df.formatted.final.loop$end1), 
                     int1=df.formatted.final.loop$start2,
                     int2=df.formatted.final.loop$end2)

mcols(hicranges)$loop_id <- df.formatted.final.loop$loop.id

hicranges # bidrection: 23844
promoterranges # 22216/18621

# step 17
olap.hic.pro <- findOverlaps(hicranges,promoterranges)

generanges <- hicranges[queryHits(olap.hic.pro)]
generanges
mcols(generanges) <- cbind(mcols(hicranges[queryHits(olap.hic.pro)]), mcols(promoterranges[subjectHits(olap.hic.pro)]))

generanges # 10368

# distinct()# reverse for the other with SNPs!
genebed <- data.frame(chr=seqnames(generanges), snp.start=generanges$int1, snp.end=generanges$int2, 
                      gene.start=start(generanges), 
                      gene.end=start(generanges)+width(generanges)-1, # generanges$ranges$end
                      loop_id = generanges$loop_id, # custom col1
                      promoter_id = generanges$promoter_id, # custom col1
                      ensg=generanges$gene
)


promoterranges

genebed # 10368 = 10243 + 125 or = 10202 + 166 / 8788 # SAME DEDUP regardless of loop)id
genebed %>% dim() # 10368

# genebed %>% distinct(loop_id, ensg) %>% dim()

# > genebed
# chr snp.start   snp.end gene.start  gene.end         ensg                                      loop_id                       promoter_id
# 1    10  15000000  15025000   14900000  14925000      Wfikkn1    10:14900000-14925000_to_15000000-15025000      Wfikkn1_10:14900000-14925000
# 2    10  16875000  16900000   16525000  16550000      Rpl26l1    10:16525000-16550000_to_16875000-16900000      Rpl26l1_10:16525000-16550000
# removing dups: same gene promoter & same loop
genebed <- unique(genebed)

genebed %>% dim() # 8881: w/o custom cols (loop.id: 8881    6) // 8881: w/custom cols: 8881    7
genesnpranges <- GRanges(genebed$chr, 
                         IRanges(genebed$snp.start, genebed$snp.end), 
                         loop_id=genebed$loop_id,
                         promoter_id=genebed$promoter_id,
                         ensg=genebed$ensg
)

# snpranges: SNPs NOT in exon and promoter
snpranges.gr <- GRanges(snpranges$chr, IRanges(snpranges$Position, snpranges$Position), rsid=snpranges$rsid)

olap.snprg.genesnprg <- findOverlaps(snpranges.gr, genesnpranges) # finding overlapping SNP & gene range (snp in exon & promoter)
snpint <- snpranges.gr[queryHits(olap.snprg.genesnprg)] # getting SNP information
mcols(snpint) <- cbind(mcols(snpint), 
                       mcols(genesnpranges[subjectHits(olap.snprg.genesnprg)])) # integrating returned SNP & gene metadata

snpint # bidirection: 359871//356848 (promoter case1 + incl X + promoter contact only//promoter case2 + excl X + promoter contact only)
# 359871 w/o custom cols // 378442 w/custom cols

# load("Hi-C_transcript_interacting_snp.rda")
# load("snp_locating_in_exon_promoter_transcript_level.rda")

snpdat <- data.frame(chr=seqnames(snpint), bp=start(snpint), rsid=snpint$rsid, ensg=snpint$ensg, loop_id=snpint$loop_id)
snpdat
snpdat %>% dim() # loop_id: [1] 359871      5
###########################################
# test: START
###########################################
# Case 1: rsid + ensg 기준 중복 제거
snpdat_unique_rsid_ensg <- snpdat %>%
  dplyr::distinct(rsid, ensg)

# Case 2: rsid + ensg + loop_id 기준 중복 제거
snpdat_unique_rsid_ensg_loop <- snpdat %>%
  dplyr::distinct(rsid, ensg, loop_id)

# 비교
nrow(snpdat)  # 원래 전체 행 수
nrow(snpdat_unique_rsid_ensg)        # 중복 제거 후 (2컬럼 기준)
nrow(snpdat_unique_rsid_ensg_loop)   # 중복 제거 후 (3컬럼 기준)

# 차이 확인
nrow(snpdat_unique_rsid_ensg) - nrow(snpdat_unique_rsid_ensg_loop)

dups <- snpdat %>%
  count(rsid, ensg) %>%
  filter(n > 1)
dups

# snpdat: chr, bp, rsid, ensg, loop_id
snpdat %>% dim() # 359871
snpdat %>% distinct() %>% dim() # 359871
snpdat %>% distinct(rsid, ensg) %>% dim() # 330535********************** : same SNP + same gene
snpdat %>% distinct(rsid, ensg, loop_id) %>% dim() # 359871 ************************ : : same SNP + same gene + diff loop

snpdat.rsid.ensg <- snpdat %>% distinct(rsid, ensg) # 330535 *********************************
snpdat.rsid.ensg %>% dim()
snpdat %>% filter(rsid == '10:101554648')
# chr        bp         rsid   ensg                                             loop_id
# 1  10 101554648 10:101554648 Qrich2 10_101675000_101700000_10_101550000_101575000_25000
# 2  10 101554648 10:101554648 Qrich2 10_101690000_101700000_10_101550000_101560000_10000
###########################################
# test: END
###########################################

snpromat <- unique(data.frame(rsid=snpro$rsid, ensg=snpro$gene))
snpexonmat <- unique(data.frame(rsid=snpexon$rsid, ensg=snpexon$gene))
snpromat
snpexonmat
snpcomb.rsid.ensg <- bind_rows(snpdat.rsid.ensg, snpromat, snpexonmat) %>% unique()  # exon + pro + hic
snpcomb.rsid.ensg %>% dim()
snpcomb <- bind_rows(snpdat %>% dplyr::select(rsid, ensg), snpromat, snpexonmat) %>% unique()  # exon + pro + hic

snpcomb # 391182 + 500 = 391682 | bidirection: 522190 + 500 = 522690 | sub.4 any: 352408 + 500 = 352908
snpcomb %>% dim() # 495164
snpcomb %>% 
  count(ensg) # bidrection: 16995(exon only)/17045 | sub.4 any: 16935 
dim() # 522190 + 500 = 522690 | sub.4 any: 352408 + 500 = 352908

# making genedef
genedef <- genes.all
genedef

write.table(genedef, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_hmagma_enhancer/0911/gene_NCBI.txt", quote=F, row.names=F, col.names=F, sep="\t")
write.table(genedef, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/gene_NCBI.txt", quote=F, row.names=F, col.names=F, sep="\t")

# dedup rows grouped by gene
snpagg <- snpcomb %>%
  group_by(ensg) %>%
  summarise(rsid = paste(unique(rsid), collapse = ", ")) %>%
  ungroup()

# snpagg <- aggregate(snpcomb, list(snpcomb$ensg), unique) # dedup rows grouped by gene
# genedef <- read.table("Gencode26_gene.bed")
# colnames(genedef) <- c("chr", "start", "end", "ensg")
# genedef <- genedef[grep("chr", genedef$chr),]
# genedef$chr <- unlist(lapply(strsplit(genedef$chr, "chr"), '[[', 2))
# genedef$index <- paste(genedef$chr, genedef$start, genedef$end, sep=":")
# snpagg$index <- genedef[match(snpagg$ensg, genedef$gene_id),"index"]

matched_indexes <- match(snpagg$ensg, genedef$gene_id)
snpagg$index <- genedef$index[matched_indexes]
snpagg # 17,461 (promoter case1 + incl X + promoter contact only) | bidrection: 16,995(exon only)/17,045 | sub.4 any: 16,935

snpaggconv <- snpagg %>% 
  filter(!is.na(index)) %>% 
  mutate(ensg = index) %>% 
  dplyr::select(ensg, index, rsid)

# 17461: no unique(snpint), no custom cols // 17,395: unique(snpint)
snpaggconv #17,395 (promoter case1 + incl X + promoter contact only) | bidrection: 16873(5k10k), 16,995(all, exon only)/17,045 | sub.4 any: 16935
snpaggconv %>% head()
write.table(snpaggconv, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_hmagma_enhancer/0911/snpagg.txt", quote=F, row.names=F, col.names=FALSE, sep="\t") # change the name of the file
write.table(snpaggconv, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/snpagg.txt", quote=F, row.names=F, col.names=FALSE, sep="\t") # change the name of the file
system("sed -e 's/, /\t/g' < /Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_hmagma_enhancer/0911/snpagg.txt > /Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_hmagma_enhancer/0911/hs.hic.annot.submission", wait = TRUE)
system("sed -e 's/, /\t/g' < /Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/snpagg.txt > /Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/hs.hic.annot.submission", wait = TRUE)

# all genes
# 5k
write.table(snpaggconv, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/snpagg.all.genes_5k.txt", quote=F, row.names=F, col.names=FALSE, sep="\t") # change the name of the file
system("sed -e 's/, /\t/g' < /Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/snpagg.all.genes_5k.txt > /Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/hs.hic.annot.sub.4.any.all.genes_5k", wait = TRUE)

# 10k
write.table(snpaggconv, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/snpagg.all.genes_10k.txt", quote=F, row.names=F, col.names=FALSE, sep="\t") # change the name of the file
system("sed -e 's/, /\t/g' < /Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/snpagg.all.genes_10k.txt > /Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/hs.hic.annot.sub.4.any.all.genes_10k", wait = TRUE)

# 25k
write.table(snpaggconv, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/snpagg.all.genes_25k.txt", quote=F, row.names=F, col.names=FALSE, sep="\t") # change the name of the file
system("sed -e 's/, /\t/g' < /Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/snpagg.all.genes_25k.txt > /Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/hs.hic.annot.sub.4.any.all.genes_25k", wait = TRUE)

# all resolution
write.table(snpaggconv, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/snpagg_all_genes_all_resol.txt", quote=F, row.names=F, col.names=FALSE, sep="\t") # change the name of the file
system("sed -e 's/, /\t/g' < /Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/snpagg_all_genes_all_resol.txt > /Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/hs.hic.annot.sub.4.any.all_genes_all_resol", wait = TRUE)

# 5k10k
write.table(snpaggconv, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_hmagma_enhancer/0911/snpagg_5k10k.txt", quote=F, row.names=F, col.names=FALSE, sep="\t") # change the name of the file
write.table(snpaggconv, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/snpagg_5k10k.txt", quote=F, row.names=F, col.names=FALSE, sep="\t") # change the name of the file
system("sed -e 's/, /\t/g' < /Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_hmagma_enhancer/0911/snpagg_5k10k.txt > /Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_hmagma_enhancer/0911/hs.hic.annot.sep.2024_5k10k", wait = TRUE)
system("sed -e 's/, /\t/g' < /Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/snpagg_5k10k.txt > /Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/hs.hic.annot.sub.4.any_5k10k", wait = TRUE)

##############################################################################################################
# figure 1: Comparison for # of SNP between hic and without hic
##############################################################################################################
snpaggconv.cmagma.annot %>% head()
df.snpaggconv.cmagma.annot.without.hic <- snpaggconv.cmagma.annot %>%
  separate_rows(rsid, sep = ",\\s*") %>%  
  mutate(chr = str_extract(rsid, "^[^:]+"))

df.snpaggconv.cmagma.annot.without.hic %>% 
  count(chr)
dim()

snpaggconv %>% head()
df.snpaggconv.with_hic <- snpaggconv %>% 
  separate_rows(rsid, sep = ",\\s*") %>%  
  mutate(chr = str_extract(rsid, "^[^:]+"))

df.snpaggconv.with_hic %>% 
  count(chr)
dim()

# num of SNP per chr
snp_count_without_hic <- df.snpaggconv.cmagma.annot.without.hic %>%
  count(chr, name = "Without Hi-C")
snp_count_without_hic

snp_count_without_hic %>% 
  summarise(tot = sum(`Without Hi-C`)) # 164629

snp_count_with_hic <- df.snpaggconv.with_hic %>%
  count(chr, name = "With Hi-C")
snp_count_with_hic %>% 
  summarise(tot = sum(`With Hi-C`)) # 495164

# by Hi-C, 495164 - 164629 = 330535

# join with two df
snp_count_compare <- full_join(snp_count_with_hic, snp_count_without_hic, by = "chr") %>%
  pivot_longer(cols = c(`With Hi-C`, `Without Hi-C`), names_to = "Condition", values_to = "SNP_count") %>%
  replace_na(list(SNP_count = 0))  # NA to 0

snp_count_compare %>% view()

pdf("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/snp_count_compare_plot.pdf", width = 11*0.8, height = 8.5*0.8)

# plot
chromosome_levels <- c(as.character(1:20), "X")  # chr1 ~ chr20, X

# chr factor level
snp_count_compare <- snp_count_compare %>%
  mutate(chr = factor(chr, levels = chromosome_levels),
         Condition = factor(Condition, levels = c("Without Hi-C", "With Hi-C"))) # without_hic, and then with_hic

# 2. plot
ggplot(snp_count_compare, aes(x = chr, y = SNP_count, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +  # bar size 
  labs(
    x = "Chromosome", 
    y = "Number of SNPs"
  ) +
  theme_bw() 

dev.off()

##############################################################################################################
# figure 2: Average number of SNPs per gene (with and without Hi-C integration).
##############################################################################################################

# with Hi-C (SNP per gene)
df.snpaggconv.with_hic # 495,164
df.snpaggconv.with_hic %>% distinct(index, rsid) # 495,164 NO DUPS

gene_snp_count_with_hic <- df.snpaggconv.with_hic %>%
  group_by(ensg) %>%
  summarise(n_snp = n(), .groups = "drop") %>%
  mutate(condition = "With Hi-C")

gene_snp_count_with_hic %>% 
  full_join(genes.all, by = c('ensg' = 'index')) %>% 
  filter(!is.na(ensg)) %>% 
  dplyr::select(ensg, n_snp, condition, gene_id)

# without Hi-C (SNP per gene)
df.snpaggconv.cmagma.annot.without.hic # 164,629
df.snpaggconv.cmagma.annot.without.hic %>% distinct(ensg, rsid) # 164,629 NO DUPS

gene_snp_count_without_hic <- df.snpaggconv.cmagma.annot.without.hic %>%
  group_by(ensg) %>%
  summarise(n_snp = n(), .groups = "drop") %>%
  mutate(condition = "Without Hi-C")

gene_snp_count_without_hic %>% 
  full_join(genes.all, by = c('ensg' = 'index')) %>% 
  filter(!is.na(ensg)) %>% 
  dplyr::select(ensg, n_snp, condition, gene_id)

# bind_rows
gene_snp_counts <- bind_rows(gene_snp_count_with_hic, gene_snp_count_without_hic)
gene_snp_counts

# mean and SD
snp_stats <- gene_snp_counts %>%
  group_by(condition) %>%
  summarise(
    mean_snp = mean(n_snp),
    sd_snp = sd(n_snp),
    .groups = "drop"
  ) %>% 
  mutate(condition = factor(condition, levels = c("Without Hi-C", "With Hi-C"))) # without_hic, and then with_hic

snp_stats
# NO unique(snpint)
# condition    mean_snp sd_snp
# <fct>           <dbl>  <dbl>
#   1 With Hi-C       28.4   50.8 
# 2 Without Hi-C     9.55   8.92

# done with unique(snpint)
# condition    mean_snp sd_snp
# <fct>           <dbl>  <dbl>
#   1 With Hi-C       22.5   42.9 
# 2 Without Hi-C     9.55   8.92

pdf("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/snp_count_per_gene_compare_violin_plot_log.pdf", width = 11*0.8, height = 8.5*0.8)

gene_snp_counts_log <- gene_snp_counts %>% filter(n_snp > 0)

# LOG
label_stats_log <- gene_snp_counts_log %>%
  group_by(condition) %>%
  summarise(
    mean_val = mean(n_snp),
    median_val = median(n_snp),
    .groups = "drop"
  )

ggplot(gene_snp_counts_log, aes(x = condition, y = n_snp, fill = condition)) +
  geom_violin(trim = FALSE, scale = "width", width = 0.8, color = NA, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
  scale_y_log10() +
  geom_text(data = label_stats_log, aes(x = condition, y = mean_val + 5, label = paste0("Mean:  ", round(mean_val, 1))),
            vjust = -1, color = "black") +
  geom_text(data = label_stats_log, aes(x = condition, y = median_val - 6, label = paste0("Median: ", round(median_val, 1))),
            vjust = 2, color = "black") +
  labs(
    x = "",
    y = "Log10(Number of SNPs per Gene)"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    legend.position = "none"
  )

dev.off()

##############################################################################################################
# figure 3: AVERAGE num of SNPs in non-coding regions which are connected to promoter in the other end of loops 
##############################################################################################################
options(scipen = 999)

genes.all %>% 
  head()

pdf("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/hic_noncoding_snps_violin_by_gene_log.pdf", width = 8.5*0.8, height = 11*0.8)

snpdat
snpdat %>% filter(ensg == 'Megf8') %>% 
  count(loop_id)
snpdat %>% dim() # 418294      5
# minimum, median, mean, and maximum

snp.summary.stats.figure3<- snpdat %>%
  distinct(rsid, ensg) %>%  
  count(ensg) %>%           
  summarise(
    min_snp = min(n),
    median_snp = median(n),
    mean_snp = mean(n),
    max_snp = max(n)
  )

snp.summary.stats.figure3

min_val <- snp.summary.stats.figure3$min_snp
median_val <- snp.summary.stats.figure3$median_snp
mean_val <- snp.summary.stats.figure3$mean_snp
max_val <- snp.summary.stats.figure3$max_snp

df.snp.count.assigned.by.loops <- snpdat %>%
  distinct(rsid, ensg) %>%
  count(ensg)  

# max SNP assigned for gene 찾기
df.snp.count.assigned.by.loops %>%
  filter(n == max(n))  

snpdat %>% 
  filter(ensg == 'Notch4') %>% # dim() # 1229    5
  # distinct(rsid) # 951
  distinct(loop_id) # 4
head()

# loop and snps
df.snpdat.log10 <- snpdat %>% 
  distinct(rsid, ensg) %>% 
  count(ensg) %>% 
  mutate(log_n = log10(n),
         label = str_c("Min: ", min_val, "\nMedian: ", median_val, "\nMean: ", mean_val,"\nMax: ", max_val))

df.snpdat.log10

ggplot(df.snpdat.log10, aes(x = "", y = log_n)) +
  geom_violin(fill = "#f9756d", color = "black", scale = "width") +
  annotate(
    "text",
    x = 1,
    y = max(df.snpdat.log10$log_n) + 0.1,
    label = str_c(
      "  Max:  ", max_val,
      "\nMean: ", round(mean_val, 2),
      "\n Median: ", median_val,
      "\n     Min: ", min_val
    ),
    hjust = -1,
    vjust = 1.5
  ) +
  theme_minimal() +
  labs(
    y = "log10(Number of Non-coding SNPs)",
    x = NULL
  )

dev.off()

####################################################
# figure 4: figure for distance between SNPs and promoters by loop
####################################################
genes.all
genes.all %>% head()

snpdat
snpdat %>% dim() # 418294      5
df.transcript.distinct.gene_id.genomic.BestRefSeq.extended
df.tss <- df.transcript.distinct.gene_id.genomic.BestRefSeq.extended %>%
  mutate(
    tss = ifelse(strand == "+", as.numeric(start), as.numeric(end))
  ) %>%
  dplyr::select(gene_id, chr, strand, tss, index)

chr_levels <- as.character(c(1:20, "X"))
df.snp.tss.dist <- left_join(snpdat %>% distinct(rsid, ensg), df.tss, by = c("ensg" = "gene_id")) %>% 
  mutate(
    snp_pos = as.numeric(str_split_n(rsid, ":", 2)), 
    tss = as.numeric(tss),
    distance = abs(snp_pos - tss)
  ) %>% filter(!is.na(distance)) %>% 
  mutate(log10_distance = log10(distance)) %>% 
  mutate(chr = factor(chr, levels = chr_levels))

df.snp.tss.dist

# PDF 저장
pdf("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/snp_tss_distance_violin_by_gene_log10.pdf", width = 11*0.8, height = 8.5*0.8)

ggplot(df.snp.tss.dist, aes(x = chr, y = log10_distance)) +
  geom_violin(fill = "#66c2a5", color = "black", scale = "width") +
  # geom_jitter(width = 0.2, height = 0, alpha = 0.3) +
  theme_minimal() +
  labs(
    # title = "SNP–TSS Distance by Chromosome (log10 scale)",
    x = "Chromosome",
    y = "log10(Distance to TSS in bp)"
  )

dev.off()

####################################################
# Discussion
####################################################
library(VennDiagram)
snpaggconv.cmagma.annot # A tibble: 17,239 × 3
snpagg  # A tibble: 17,461 × 3

genes.all
df.genes.exon.pro.only <- snpaggconv.cmagma.annot %>% distinct(ensg) %>% 
  left_join(genes.all, by = c ("ensg" = "index")) %>% 
  dplyr::select(gene_id)
df.genes.exon.pro.hic <- snpagg %>% distinct(ensg)
df.genes.exon.pro.only
df.genes.exon.pro.hic

full_join(df.genes.exon.pro.only, df.genes.exon.pro.hic, by = c("gene_id" = "ensg"))

df.genes.exon.pro.only
df.genes.exon.pro.hic
# 두 벡터로 변환
genes_A <- df.genes.exon.pro.only$gene_id
genes_B <- df.genes.exon.pro.hic$ensg
genes_A %>% length()
genes_B %>% length()

# 교집합
intersect_genes <- intersect(genes_A, genes_B)

# 교집합 수 확인
length(intersect_genes)

# 순수 A-only
A_only <- setdiff(genes_A, genes_B)
length(A_only)

# 순수 B-only
B_only <- setdiff(genes_B, genes_A)
length(B_only)

###########################
ensg_conv <- snpaggconv.cmagma.annot$ensg
ensg_conv
# 2. snpagg에서 ensg_conv에 포함되지 않은 222개 유전자 추출
snpagg_extra_genes <- snpagg %>%
  filter(!(index %in% ensg_conv))
snpagg_extra_genes # 222

# parsing
snp_count_extra <- snpagg_extra_genes %>%
  mutate(num_snps = str_count(rsid, ",") + 1) %>%  # 쉼표 수 + 1 = SNP 개수
  summarise(total_snps = sum(num_snps))

snp_count_extra$total_snps

#####################
total_snps_conv <- snpaggconv.cmagma.annot %>%
  mutate(num_snps = str_count(rsid, ",") + 1) %>%
  summarise(total = sum(num_snps)) %>%
  pull(total)
total_snps_conv

total_snps_snpagg <- snpagg %>%
  mutate(num_snps = str_count(rsid, ",") + 1) %>%
  summarise(total = sum(num_snps)) %>%
  pull(total)
total_snps_snpagg
# snpaggconv.cmagma.annot의 SNP 수 총합
# 결과 출력
cat("snpagg total SNPs: ", total_snps_snpagg, "\n") # 495164
cat("snpaggconv.cmagma.annot total SNPs: ", total_snps_conv, "\n") # 164629





figure1.venn.plot <- venn.diagram(
  x = list(
    `Exon+Promoter` = df.genes.exon.pro.only$gene_id,
    `Exon+Promoter+Hi-C` = df.genes.exon.pro.hic$ensg
  ),
  filename = NULL,  
  col = "black",
  fill = c("#66c2a5", "#fc8d62"),
  alpha = 0.5,
  cex = 1.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-20, 20),
  cat.dist = 0.05
)

class(figure1.venn.plot)

pdf("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/venn_genes_exonpromoter_vs_hic.pdf", width = 11*0.8, height = 8.5*0.8)
grid.draw(figure1.venn.plot)
dev.off()


# result section 
# 1. genes assined with exon/promoter
snpaggconv.cmagma.annot

df.genes.coding.snp <- snpaggconv.cmagma.annot %>%
  left_join(genedef, by = "index") %>% 
  separate_rows(rsid, sep = ",\\s*") %>%
  group_by(gene_id) %>%
  summarise(coding_snps = n(), .groups = "drop")

df.genes.coding.snp %>% 
  summarise(tot = sum(coding_snps))



# df.genes.coding.snp %>% filter(is.na(ensg))

# 2. genes assined with non-coding SNP 
df.genes.noncoding.snp <- df.snpint.figure4 %>%
  group_by(ensg) %>%
  summarise(loop_snps = n(), .groups = "drop")

df.genes.noncoding.snp

# correct result: 3056 PASS
intersect(unique(df.genes.coding.snp$gene_id), unique(df.genes.noncoding.snp$ensg)) %>% length()

# 3. join
df.merged.genes.assinged.snp.group <- full_join(df.genes.coding.snp, df.genes.noncoding.snp, by = c("gene_id" = "ensg")) %>%
  mutate(
    group = case_when(
      !is.na(coding_snps) & is.na(loop_snps) ~ "OC",
      is.na(coding_snps) & !is.na(loop_snps) ~ "ONC",
      !is.na(coding_snps) & !is.na(loop_snps) ~ "BOTH"
    )
  )

df.snp.count.by.group <- df.merged.genes.assinged.snp.group %>%
  group_by(group) %>%
  summarise(
    total_coding_snps = sum(as.numeric(coding_snps), na.rm = TRUE),
    total_loop_snps = sum(as.numeric(loop_snps), na.rm = TRUE),
    total_snps = total_coding_snps + total_loop_snps,
    .groups = "drop"
  )

df.genes.count.by.group <- df.merged.genes.assinged.snp.group %>% count(group)
df.final.genes.group.coding.noncoding.snp <- inner_join(df.snp.count.by.group, df.genes.count.by.group, by = "group")
df.final.genes.group.coding.noncoding.snp
# group total_coding_snps total_loop_snps total_snps     n
# <chr>             <dbl>           <dbl>      <dbl> <int>
#   1 BOTH              29159          219273     248432  3056
# 2 OC               135470               0     135470 14183
# 3 ONC                   0            7780       7780   156

# result section
df.genes.hmagma.snp <- snpaggconv %>% left_join(genedef, by = "index") %>% 
  separate_rows(rsid, sep = ",\\s*") %>%
  group_by(gene_id) %>%
  summarise(hmagma_snps = n(), .groups = "drop")
df.genes.hmagma.snp

df.annot.cmagma.hmagma <- 
  full_join(df.genes.coding.snp, df.genes.hmagma.snp, by = "gene_id") %>% 
  mutate(non_coding_snps = hmagma_snps - coding_snps) %>% 
  mutate(
    group = case_when(
      non_coding_snps == 0 ~ "OC",
      non_coding_snps > 0 ~ "BOTH",
      is.na(coding_snps) ~ "ONC"
    )
  ) %>% 
  mutate(non_coding_snps = if_else(is.na(coding_snps), hmagma_snps, non_coding_snps))

df.annot.cmagma.hmagma 
df.annot.cmagma.hmagma %>% filter(is.na(coding_snps)) # 156: PASS
df.annot.cmagma.hmagma %>% filter(non_coding_snps == 0) # 14,183: PASS
df.annot.cmagma.hmagma %>% filter((non_coding_snps > 0) & (non_coding_snps != hmagma_snps)) # 3,056: PASS
df.annot.cmagma.hmagma %>% count(group)

df.annot.cmagma.hmagma %>% 
  # summarise(snp_hmagma = sum(hmagma_snps), .groups = "drop") # hmagma snps: 391682
  # summarise(snp_cmagma = sum(replace_na(coding_snps, 0)), .groups = "drop") # cmagma snps: 164629
  # filter(group == "ONC") %>% summarise(non_coding_snps_tot = sum(non_coding_snps), .groups = "drop") # 7780
  # filter(group == "BOTH") %>% summarise(non_coding_snps_tot = sum(non_coding_snps), .groups = "drop") # 219273
  summarise(whole_snps_from_hmagma = sum(hmagma_snps), .groups = "drop") # 391682
# 164629 (cmagma) + 219273(snps by loops for 3056 genes) + 7780 (snps by loops for 156 genes) = 391682: PASS

####
snp.megf8 <- snpcomb %>% filter(ensg == "Megf8")
snp.megf8 %>% dim() # 40
snp.megf8 %>% distinct() %>% dim()
snp.megf8.rsids <- snp.megf8$rsid

snp.megf8.rsids
write.table(
  snp.megf8$rsid,
  file = "/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/megf8_snps.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


