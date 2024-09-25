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
setwd('/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/')
source(file.path('/Users/PanjunKim/dropbox/Gateway_to_Hao/project_common_code/', 'variables.R'))
source(file.path('/Users/PanjunKim/dropbox/Gateway_to_Hao/project_common_code/', 'funcs.R'))
getwd()

#################################################
# data preprocessing
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

### 2. NCBI refseq data preprocessing
## 2-1. BestRefSeq: NCBI의 RefSeq 프로젝트에서 제공하는, 가장 신뢰할 수 있는 어노테이션을 가진 서열을 나타냅니다. 이는 엄격한 검증 과정을 거쳐 선택된 서열로, 고품질의 유전체 어노테이션을 제공합니다.
## 2-2. BestRefSeq,Gnomon: 이 값은 BestRefSeq 어노테이션과 Gnomon 어노테이션 방법을 모두 사용하여 어노테이션된 서열을 나타냅니다. Gnomon은 NCBI에서 개발한 예측 기반 어노테이션 도구입니다.
## 2-3. Curated Genomic: 전문가들에 의해 수동으로 큐레이션된 유전체 서열을 의미합니다. 이 데이터는 공개된 연구 결과나 실험적 검증을 통해 얻어진 정보를 기반으로 합니다.
## 2-4. Gnomon: NCBI의 자체 개발 어노테이션 도구로, 주로 유전자 예측에 사용됩니다. Gnomon은 컴퓨터 알고리즘을 사용하여 유전자의 위치, 구조, 기능을 예측합니다.
## 2-5. RefSeq: NCBI의 RefSeq 프로젝트에서 제공하는 참조 서열 데이터베이스를 의미합니다. RefSeq는 유전자, 단백질, 유전체 서열 등에 대한 포괄적인 어노테이션 정보를 제공합니다.
## 2-6. cmsearch: 이는 "covariance model search"의 약자로, 주로 리보핵산(RNA) 구조를 기반으로 한 서열 검색에 사용되는 도구입니다. cmsearch는 특정 RNA 구조에 대한 서열의 존재 여부를 확인하는 데 사용됩니다.
## 2-7. tRNAscan-SE: tRNA 유전자를 찾기 위해 특별히 설계된 컴퓨터 프로그램입니다. 이 도구는 유전체 데이터에서 tRNA 유전자를 식별하고 어노테이션하는 데 널리 사용됩니다.
## BestRefSeq (chosen) cf. RefSeq

# GTF file path
file_path <- "/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/data/gene_data/ncbi_dataset/ncbi_dataset/data/GCF_015227675.2/genomic_BestRefSeq.gtf"
# file_path2 <- "/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/data/gene_data/ncbi_dataset/hao/hao_refSeq_BestRefSeq.tsv"

df.genomic.BestRefSeq <- read_delim(file_path, 
                 delim = "\t", 
                 col_names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"),
                 quote = "",
                 trim_ws = TRUE)
df.genomic.BestRefSeq %>% count(feature)
# from GCF_015227675.2
#1 CDS         192169
#2 exon        205520 *************
#3 gene          7539 *************
#4 start_codon  20613 *************
#5 stop_codon   20590 *************
#6 transcript   22230 *************

TSS

########## LONG RUNNING TIME 
load("./from_hmagma_enhancer/0911/df_genomic_BestRefSeq_extended.rda")

df.genomic.BestRefSeq.extended <- df.genomic.BestRefSeq %>% 
  mutate(attributes = str_split(attributes, "; ")) %>%
  unnest(attributes) %>%
  separate(attributes, into = c("key", "value"), sep = " ", extra = "merge") %>%
  mutate(value = str_remove(value, '^"|"$$')) %>%
  pivot_wider(names_from = key,
              values_from = value,
              values_fn = list(value = function(x) paste(unique(x), collapse = "|")),
              values_fill = list(value = NA)) %>%
  dplyr::select(seqname, source, feature, start, end, score, strand, frame, db_xref, 
         exon_number, gbkey, gene, gene_id, protein_id, tag, transcript_id)

save(df.genomic.BestRefSeq.extended, file="./from_hmagma_enhancer/0911/df_genomic_BestRefSeq_extended.rda")

#################################################
# Make GenomicRanges objects for exons and promoter
#################################################
#########
## exons
#########

# Required format from Won
the_won_exon <- read.table("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/actual_HMAGMA/Input_Files/Gencode26_exon.bed")
the_won_exon
# V1       V2       V3              V4 (Ensembl Gene ID)
# 1  X 99883667 99884983 ENSG00000000003
# 2  X 99885756 99885863 ENSG00000000003
# 3  X 99887482 99887565 ENSG00000000003
# 239  1 209786130 209786257 ENSG00000008118
# 240  1 209786385 209787283 ENSG00000008118

##
# EXON is extracting from exon among CDS, exon, gene, start, stop, transcript
df.genomic.BestRefSeq.extended %>% head(16)

df.exon.genomic.BestRefSeq.extended <- df.genomic.BestRefSeq.extended %>%
  filter(feature == "exon") %>%
  dplyr::select(seqname, start, end, strand, gene_id) %>%
  mutate(gene_id = str_remove_all(gene_id, '\"'),
         gene_id = str_trim(gene_id)) %>% 
  arrange(gene_id, start, end)

df.genomic.BestRefSeq.extended %>% count(seqname)
df.exon.genomic.BestRefSeq.extended
df.exon.genomic.BestRefSeq.extended %>% count(seqname)
df.exon.genomic.BestRefSeq.extended %>% colnames()

chromosome.length # chr      length genbank    refseq  
df.exon.genomic.BestRefSeq.extended %>% colnames() # "seqname" "start"   "end"     "strand"  "gene_id"

exon <- right_join(chromosome.length, df.exon.genomic.BestRefSeq.extended, by = c("refseq" = "seqname")) %>%
  mutate(start = as.numeric(start),
         end = as.numeric(end)) %>% 
  filter(!is.na(chr))

exon # 177104

# WGS contig: 53
# A tibble: 6 × 2
# refseq             n
# <chr>          <int>
# 1 NW_023637717.1    12
# 2 NW_023637718.1    23
# 3 NW_023637733.1     9
# 4 NW_023637831.1     1
# 5 NW_023637849.1     3
# 6 NW_023637854.1     5

exonranges <- GRanges(seqnames = exon$chr, ranges = IRanges(start = exon$start, end = exon$end), gene = exon$gene_id)

#############
## promoters
#############
# PROMOTER is extracting from transcript among CDS, exon, gene, start, stop, transcript

df.transcript.genomic.BestRefSeq.extended <- df.genomic.BestRefSeq.extended %>%
  filter(feature == "transcript") %>%
  dplyr::select(seqname, start, end, strand, gene, gene_id) %>% 
  mutate(gene_id = str_remove_all(gene_id, '\"'),
         gene_id = str_trim(gene_id)) %>% 
  # mutate(gene = str_remove_all(gene, '\"'),
  #        gene = str_trim(gene))
  arrange(gene_id, start, end)
# all gene == gene_id

df.transcript.genomic.BestRefSeq.extended

promoter <- right_join(chromosome.length, df.transcript.genomic.BestRefSeq.extended, by = c("refseq" = "seqname")) %>%
  mutate(start_new = as.numeric(start) - 2000,
         end = as.numeric(start)) %>%
  filter(!is.na(chr))

promoter # 20448 (is.na(chr) : 12)

promoterranges <- GRanges(
  seqnames = promoter$chr,
  ranges = IRanges(start = promoter$start_new, end = promoter$end),
  gene = promoter$gene_id
)

save(exonranges, promoterranges, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_hmagma_enhancer/0911/exon_promoranges.rda")
load ("from_hmagma_enhancer/0911/exon_promoranges.rda")

wb <- createWorkbook()

addWorksheet(wb, "Sheet1")
writeData(wb, "Sheet1", exon)

addWorksheet(wb, "Sheet2")
writeData(wb, "Sheet2", promoter)

saveWorkbook(wb, "/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_hmagma_enhancer/0911/df_exon_promoter.xlsx", overwrite = TRUE)

exon.from.excel <- read.xlsx("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_hmagma_enhancer/0911/df_exon_promoter.xlsx", sheet = 1) 
promoter.from.excel <- read.xlsx("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_hmagma_enhancer/0911/df_exon_promoter.xlsx", sheet = 2)

exon.from.excel # 177104
promoter.from.excel # 20448

# line: 23~32 H-MAGMA_annotation_generation.rmd

#################################################
# Generate a GenomicRanges object for the SNP annotation
#################################################
# required format comaparison from Won
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
colnames(snps) <- c("chr","SNP","Position")
snps %>% head()
# snps$chr <- sub("^", "chr", snps$chr)
snps <- snps %>% dplyr::rename(rsid = SNP)
snps

snps.g <- GRanges(snps$chr, IRanges(snps$Position, snps$Position), rsid=snps$SNP)
snps.g # 7358643
save(snps.g, file="snps.rda")
load ("snps.rda")

# line: 23~32 H-MAGMA_annotation_generation.rmd
#################################################
# Overlap exons with SNPs
#################################################
# SNP & exon
olap <- findOverlaps(snps.g, exonranges)
olap # 76626
snpexon <- snps.g[queryHits(olap)]
snpexon # 76626

mcols(exonranges[subjectHits(olap)]) # 76626
mcols(snpexon) <- cbind(mcols(snpexon), mcols(exonranges[subjectHits(olap)]))
snpexon  # 76626

#################################################
# Overlap promoters with SNPs
#################################################
# SNP & promoter
olap <- findOverlaps(snps.g, promoterranges);
olap # 96055
snpro <- snps.g[queryHits(olap)];
snpro # 96055

mcols(promoterranges[subjectHits(olap)])
mcols(snpro) <- cbind(mcols(snpro), mcols(promoterranges[subjectHits(olap)]))

snpexon # 76626
snpro # 96055

# Excluding data in X chromosome 
snpexon <- snpexon[seqnames(snpexon)!="X"]
snpexon # 75832/76626
snpro <- snpro[seqnames(snpro)!="X"]
snpro # 94852/96055
save(snpro, snpexon, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_hmagma_enhancer/0911/snp_locating_in_exon_promoter_transcript_level.rda")

#load("snp_locating_in_exon_promoter_transcript_level.rda")

# line: 58~72 H-MAGMA_annotation_generation.rmd
################################################################################################################################################ annotation for cMAGMA ############################
###################################################################################################################################################################################################

#################################################
# Assign SNPs to genes using Hi-C.
#################################################

snps # line 229, data read from .bim,  7358643
snpexon # 75832
snpro # 94852
snpexon$rsid
snps$rsid

snpranges <- snps[!(snps$rsid %in% snpexon$rsid), ] # snps NOT in exon area 
snpranges # 7283977 != 7358643 (snps) - 75832 (snpexon)

snpranges <- snpranges[!(snpranges$rsid %in% snpro$rsid), ] # snps NOT in promoter area
snpranges # 7197182 != 7283977 (snpranges) - 94852 (snpro)
#snpranges: snps NOT in exon & promo -> snpranges: non_exonic_promoter_snp.rda
save(snpranges,file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_hmagma_enhancer/0911/non_exonic_promoter_snp.rda")

# load("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/non_exonic_promoter_snp.rda") 
# load("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/exon_promoranges.rda")

#####################
# required Hi-C format from WON
# chrom1 start1 end1 chrom2 start2 end2
# chr10 100080000 100090000 chr10 100920000 100940000
# chr10 100080000 100090000 chr10 101020000 101040000
# chr10 100180000 100190000 chr10 99530000 99550000
# chr10 100180000 100190000 chr10 100960000 100980000

# load loops filtered
load("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/0909/df_final_loop.rda")
df.final.loop

df.formatted.final.loop <- df.final.loop %>%
  separate(loop.id, into = c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "resolution"), sep = "_") %>%
  mutate(across(c(start1, end1, start2, end2), as.numeric)) %>% 
  mutate(chrom1 = str_remove(chrom1, "chr"),
         chrom2 = str_remove(chrom2, "chr"))

hicranges <- GRanges(df.formatted.final.loop$chrom1, IRanges(df.formatted.final.loop$start1, df.formatted.final.loop$end1), 
                     int1=df.formatted.final.loop$start2,int2=df.formatted.final.loop$end2)

olap <- findOverlaps(hicranges,exonranges)
exonint <- hicranges[queryHits(olap)]
mcols(exonint) <- cbind(mcols(hicranges[queryHits(olap)]), mcols(exonranges[subjectHits(olap)]))

olap <- findOverlaps(hicranges,promoterranges)
proint <- hicranges[queryHits(olap)]
mcols(proint) <- cbind(mcols(hicranges[queryHits(olap)]), mcols(promoterranges[subjectHits(olap)]))

generanges <- c(exonint, proint)
genebed <- data.frame(chr=seqnames(generanges), snp.start=generanges$int1, snp.end=generanges$int2, 
                      gene.start=start(generanges), gene.end=start(generanges)+width(generanges)-1, ensg=generanges$gene)
genebed
# 37691 + 166 = 37857
# > genebed
# chr snp.start   snp.end gene.start  gene.end     ensg
# 1     1 121900000 121925000  121250000 121275000    Ttc23
# 2     1 121900000 121925000  121250000 121275000    Ttc23

genebed <- unique(genebed) 
genebed
# 8973 + 499 = 9472/37857 = 0.25
genesnpranges <- GRanges(genebed$chr, IRanges(genebed$snp.start, genebed$snp.end), ensg=genebed$ensg)
snpranges.gr <- GRanges(snpranges$chr, IRanges(snpranges$Position, snpranges$Position), rsid=snpranges$rsid)

olap <- findOverlaps(snpranges.gr, genesnpranges) # finding overlapping SNP & gene range (snp in exon & promoter)
snpint <- snpranges.gr[queryHits(olap)] # getting SNP information
mcols(snpint) <- cbind(mcols(snpint), mcols(genesnpranges[subjectHits(olap)])) # integrating returned SNP & gene metadata
snpint # 368880
snpint <- unique(snpint)
snpint # 195013/368880

save(snpint, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_hmagma_enhancer/0911/Hi-C_transcript_interacting_snp.rda")

## Generate a H-MAGMA variant-gene annotation file -> gene annotation compatible for MAGMA

# load("Hi-C_transcript_interacting_snp.rda")
# load("snp_locating_in_exon_promoter_transcript_level.rda")

snpdat <- data.frame(chr=seqnames(snpint), bp=start(snpint), rsid=snpint$rsid, ensg=snpint$ensg)

snpromat <- unique(data.frame(rsid=snpro$rsid, ensg=snpro$gene))
snpexonmat <- unique(data.frame(rsid=snpexon$rsid, ensg=snpexon$gene))
snpcomb <- unique(rbind(snpdat[,3:4], snpromat, snpexonmat)) 
snpcomb # 359421
save(snpcomb, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_hmagma_enhancer/0911/SNP_to_transcript_comb.rda")

# making genedef
df.gene.genomic.BestRefSeq.extended <- df.genomic.BestRefSeq.extended %>%
  filter(feature == "gene") %>%
  dplyr::select(seqname, start, end, strand, gene_id) %>%
  mutate(gene_id = str_remove_all(gene_id, '\"'),
         gene_id = str_trim(gene_id)) %>% 
  arrange(gene_id, start, end)

df.gene.genomic.BestRefSeq.extended

genedef <- right_join(chromosome.length, df.gene.genomic.BestRefSeq.extended, by = c("refseq" = "seqname")) %>%
  filter(!is.na(chr)) %>% 
  mutate(index = str_c(chr, ":", start, ":", end)) %>% 
  mutate(start = as.numeric(start),
         end = as.numeric(end))
  # dplyr::select(ID, chr, start, end, strand, gene_id) %>%
genedef
write.table(genedef, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_hmagma_enhancer/0911/gene_NCBI.txt", quote=F, row.names=F, col.names=F, sep="\t")

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
snpagg # 16916
snpaggconv <- snpagg %>% filter(!is.na(index))
snpaggconv #6710
snpaggconv <- snpaggconv %>% mutate(ensg = index) %>% dplyr::select(ensg, index, rsid)
snpaggconv
write.table(snpaggconv, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_hmagma_enhancer/0911/snpagg.txt", quote=F, row.names=F, col.names=FALSE, sep="\t") # change the name of the file

system("sed -e 's/, /\t/g' < /Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_hmagma_enhancer/0911/snpagg.txt > /Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_hmagma_enhancer/0911/hs.hic.annot.sep.2024", wait = TRUE)

################################# work done on Sep 2024
# TODO 1: gene dataset comparison (ucsc_start_codon.txt vs genesets from GCF_015227675.2/genomic_BestRefSeq.gtf)
# TODO 2: promoter dataset comparison (epdnew/001/Rn_EPDnew_001_rn7_1.bed vs genesets from GCF_015227675.2/genomic_BestRefSeq.gtf)
# TODO 3: TSS

# Post-process after magma running
magam.file.list <- fs::dir_ls("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/sep_HMAGMA", regexp = ".out$")

print(magam.file.list)

# file reader function
file_data_reader <- function(file) {
  file_ext <- tools::file_ext(file)
  
  if (file_ext == "out" || file_ext == "raw" || file_ext == "tsv") {
    lines <- readLines(file)
    clean_lines <- lines[!grepl("^#", lines)]  # Remove lines starting with #
    
    # Read the file as tab-separated
    read.table(text = clean_lines, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
  } else {
    stop("Unsupported file type")
  }
}

# Updated function to merge files and add filenames as identifiers
init.file.df <- function(file.list, deli_dir) {
  file.list %>%
    map_dfr(file_data_reader, .id = "sample") %>%
    # Extract the meaningful part of the filename to label the data
    # mutate(sample = str_split_n(sample, str_c(deli_dir, '\\/'), 2)) %>% 
    mutate(sample = str_split_n(str_split_n(sample, str_c(deli_dir, '\\/'), 2), str_c(deli_dir, '\\/'), 1)) # Keep only part of the filename before the first underscore
}

df.magma.out.task.init <- init.file.df(magam.file.list, "/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/sep_HMAGMA")
df.magma.out.task.init %>% head()
names(df.magma.out.task.init)


Limitations

An important limitation of H-MAGMA that should be taken into consideration is that, although
H-MAGMA detects risk genes associated with a trait, it cannot determine the directionality of the
effects, such that risk genes may be upregulated or downregulated in the diseased state. 
This limitation can be remedied by incorporating gene expression datasets such as expression quantitative
trait loci (eQTL). For instance, when an eQTL is detected for an H-MAGMA–associated gene in a
matching tissue, the eQTL can provide an indication of whether the risk allele of the variant is
associated with upregulation or downregulation of the corresponding gene.

gene-level : a gene - qlt (SNP there are in exon, promoter, by Hi-C tran-qtl, )

all infunsion













































# loop.file.list = fs::dir_ls("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/data/loops", regexp = ".bedpe$")
# # BED fild for loops
# df.init.loop.bed <- init.bedpe.df(loop.file.list, 'loops') %>% 
#   filter(!str_detect(X.chr1, "^#")) 
# 
# df.init.loop.bed %>% head()
# 
# hic <- df.init.loop.bed %>% 
#   select(1:7) %>% 
#   mutate(chr1 = X.chr1) %>% 
#   select(-2) %>% 
#   relocate(chr1, .before=2) %>% 
#   mutate(chr1 = gsub("chr", "", chr1), chr2 = gsub("chr", "", chr2))
# 
# hic %>% head()

###################################################################
####### from protocol paper: start ################################
###################################################################
hic.int1 <- hic
hic.int2 <- hic[,c(1, 5:7,2:4)]
colnames(hic.int1) = colnames(hic.int2) = c("sample", "chrom1", "start1", "end1", "chrom2", "start2", "end2")
hic.comb <- rbind(hic.int1, hic.int2)
hicranges <- GRanges(hic.comb$chrom1, IRanges(as.numeric(hic.comb$start1), as.numeric(hic.comb$end1)), int1=hic.comb$start2, int2=hic.comb$end2, sample = hic.comb$sample)
hicranges
###################################################################
####### from protocol paper: end ##################################
###################################################################

####################################################################
# recalculation for SNPs not falling in exon and promoter - START ##
####################################################################
#################################
# Getting SNPs in NON-coding region (not in exon, promoter, but in ENDS)
#################################

# setdiff.snps.exon <- setdiff(snps.g, exonranges) # snps not in exon
# setdiff.snps.exon.pro <- setdiff(setdiff.snps.exon, promoterranges) # snps not in promotor & exon
# 
# setdiff.snps.exon.pro # snps not in promotor & exon: 7,066,472
# 
# hicranges.whole <- GRanges(hic$chr1, IRanges(as.numeric(hic$x1), as.numeric(hic$y2)))
# index.non.exon.pro.whole.hic <- findOverlaps(setdiff.snps.exon.pro, hicranges.whole)
# 
# 
# setdiff.snps.exon.pro.whole.hic <- setdiff.snps.exon.pro[queryHits(index.non.exon.pro.whole.hic)];
# setdiff.snps.exon.pro.whole.hic # 99,322,422 (one snp to multiple loops)
# 
# hicranges.nonend <- GRanges(hic$chr1, IRanges(as.numeric(hic$x2), as.numeric(hic$y1)))
# 
# snp.in.ends.of.loops <- setdiff(setdiff.snps.exon.pro.whole.hic, hicranges.nonend)
# snp.in.ends.of.loops # snps in each end: 122615
###################################
####################################################################
# recalculation for SNPs not falling in exon and promoter - END ####
####################################################################

# PAN
# chr1        x1        x2  chr2        y1        y2 sample
# 1   chr10  64550000  64575000 chr10  64650000  64675000  592BB


hicranges <- GRanges(hic$chr1, IRanges(as.numeric(hic$x1), as.numeric(hic$x2)), int1=hic$y1,int2=hic$y2)

####################################################################
####################################################################

hic %>% head()

hicranges <- GRanges(hic$chr1, IRanges(as.numeric(hic$x1), as.numeric(hic$x2)), int1=hic$y1, int2=hic$y2, sample = hic$sample)
hicranges # 58992

# OVERLAPPING hic + exon RANGES
olap <- findOverlaps(hicranges, exonranges);
exonint <- hicranges[queryHits(olap)]; # upstream: x1 & x2
exonint # hic upstream, overlapping with exons 
mcols(exonint) <- cbind(mcols(hicranges[queryHits(olap)]), mcols(exonranges[subjectHits(olap)]))
exonint

# OVERLAPPING hic + promo RANGES
olap <- findOverlaps(hicranges, promoterranges);
proint <- hicranges[queryHits(olap)]; # hic overlapping with promos
mcols(proint) <- cbind(mcols(hicranges[queryHits(olap)]), mcols(promoterranges[subjectHits(olap)]))
proint
generanges <- c(exonint, proint)
#############################
generanges %>% head()
generanges # 121482

# ranges(x1, x2) strand |      int1      int2        gene
genebed <- data.frame(chr=seqnames(generanges), 
                      snp.start=generanges$int1, # y1
                      snp.end=generanges$int2, # y2
                      gene.start=start(generanges), # x1
                      gene.end=start(generanges)+width(generanges)-1, # x2
                      ensg=generanges$gene)
genebed # 121482
genebed %>% head()
genebed <- unique(genebed) 
genebed # 17587 (14.7%)

# range with snp start and snp end
genesnpranges <- GRanges(genebed$chr, IRanges(genebed$snp.start, genebed$snp.end), ensg=genebed$ensg)

snpranges <- makeGRangesFromDataFrame(snpranges,
                                      seqnames.field = "chr",
                                      start.field = "Position",
                                      end.field = "Position",
                                      keep.extra.columns = TRUE)

olap <- findOverlaps(snpranges, genesnpranges)
snpint <- snpranges[queryHits(olap)]
mcols(snpint) <- cbind(mcols(snpranges[queryHits(olap)]), mcols(genesnpranges[subjectHits(olap)]))
snpint <- unique(snpint)

snpint

save(snpint, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/Hi-C_transcript_interacting_snp.rda")
load("Hi-C_transcript_interacting_snp.rda")

#################################################
# Generate a H-MAGMA variant-gene annotation file
#################################################

load("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/Hi-C_transcript_interacting_snp.rda")
load("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/snp_locating_in_exon_promoter_transcript_level.rda")

# snpdat from Hi-C
snpdat <- data.frame(chr=seqnames(snpint), bp=start(snpint), rsid=snpint$rsid, ensg=snpint$ensg)
snpdat

snpro$rsid
snpro$gene

snpromat <- unique(data.frame(rsid=snpro$rsid, ensg=snpro$gene))
snpexonmat <- unique(data.frame(rsid=snpexon$rsid, ensg=snpexon$gene))
snpcomb <- unique(rbind(snpdat[,3:4], snpromat, snpexonmat)) # snpcomb = all two dataset: cMAGMA(exonic, promoter) + Hi-C

snpromat
snpexonmat
snpcomb

save(snpcomb, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/SNP_to_transcript_comb.rda")
snpcomb$ensg

snpagg <- aggregate(snpcomb, list(snpcomb$ensg), unique)

?aggregate
options(dplyr.width = Inf)
print(snpagg, n = Inf)

snpagg <- snpcomb %>%
  group_by(ensg) %>%
  summarise(rsid = paste(unique(rsid), collapse = ", ")) %>%
  ungroup()

write.table(snpagg, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/snpagg.txt", quote=F, row.names=F, col.names=TRUE, sep="\t") # change the name of the file
snpagg %>% head

genedef <- read.table("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/H-MAGMA/Input_Files/Gencode26_gene.bed")
genedef 



# df.gene.genomic.BestRefSeq.extended <- df.genomic.BestRefSeq.extended %>%
#   filter(feature == "gene") %>%
#   select(seqname, start, end, strand, gene_id) %>%
#   mutate(gene_id = str_remove_all(gene_id, '\"'),
#          gene_id = str_trim(gene_id)) %>% 
#   arrange(gene_id, start, end)

df.gene.genomic.BestRefSeq.extended %>% head()

right_join(meta_chr_NC, df.gene.genomic.BestRefSeq.extended, by = c("acc" = "seqname")) %>%
  filter(!is.na(acc) & !is.na(Chromosome)) %>% 
  select(-acc) %>%
  mutate(chr = Chromosome) %>% 
  mutate(ID = str_c("hs:", chr, ":", start, ":", end)) %>% 
  # select(gene_id, chr, start, end, strand) %>% 
  select(ID, chr, start, end, strand, gene_id) %>%
  write.table(file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/gene_NCBI.txt", quote=F, row.names=F, col.names=F, sep="\t")


df.gene.genomic.BestRefSeq.extended
df.gene.genomic.BestRefSeq.extended %>% head()
gene.info <- right_join(meta_chr_NC, df.gene.genomic.BestRefSeq.extended, by = c("acc" = "seqname")) %>%
  # na.exon <- right_join(meta_chr_NC, df.exon.genomic.BestRefSeq.extended, by = c("acc" = "seqname")) %>% 
  mutate(chr = Chromosome,
         start = as.numeric(start),
         end = as.numeric(end)) %>%
  filter(!is.na(chr) & !is.na(Chromosome)) %>% 
  # filter(is.na(chr) | is.na(Chromosome)) %>%
  select(-Chromosome) %>% 
  select(chr, start, end, gene_id)



## WON
# V1        V2        V3              V4
# 1    chrX  99882106  99892101 ENSG00000000003
# 2    chrX  99839799  99854882 ENSG00000000005

# PAN
# chr       start       end gene_id
# <chr>     <dbl>     <dbl> <chr>  
# 1 1     242664657 242723239 Abcc2  
# 2 1     133217403 133298564 Abhd2 

gene.info # 7531
genedef # 59997

genedef <- gene.info
genedef

# colnames(genedef) <- c("chr", "start", "end", "gene_name")
# no need to do
# genedef <- genedef[grep("chr", genedef$chr),]
# genedef
# genedef$chr <- unlist(lapply(strsplit(genedef$chr, "chr"), '[[', 2))

genedef$index <- paste(genedef$chr, genedef$start, genedef$end, sep=":")
genedef$gene_id <- trimws(genedef$gene_id)
genedef

snpagg$ensg <- trimws(snpagg$ensg)
snpagg

# matched_indexes <- match(snpagg$upper_ensg, genedef$upper_gene_id)
matched_indexes <- match(snpagg$ensg, genedef$gene_id)
snpagg$index <- genedef$index[matched_indexes]
snpagg
write.table(snpagg, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/snpagg222.txt", quote=F, row.names=F, col.names=TRUE, sep="\t") # change the name of the file

# snpagg$index <- genedef[match(snpagg$ensg, genedef$gene_id),"index"]

colnames(genedef)
snpagg
snpagg <- snpagg[!is.na(snpagg$index),]
snpagg 
dim(snpagg) # 5513
colnames(snpagg)
write.table(snpagg, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/snpagg333.txt", quote=F, row.names=F, col.names=TRUE, sep="\t") # change the name of the file

snpaggconv <- snpagg[,c("ensg", "index", "rsid")]
write.table(snpaggconv, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/SNP_aggregate_transcript.txt", quote=F, row.names=F, col.names=F, sep="\t") # change the name of the file

###### NO!!!!
snpaggconv
writable <- format(snpaggconv)
writable
write.table(writable, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/SNP_aggregate_transcript_afterformat.txt", quote=F, row.names=F, col.names = F, sep="\t") # change the name of the file
######

system("sed -e 's/, /\t/g' < SNP_aggregate_transcript.txt > hs.hic.annot", wait = TRUE)


getwd()

####################################################################
## GWAS summary statistics
####################################################################

# PD sample
# chr	bp	A1	A2	freq	b	se	p	N_cases	N_controls	rs	N
# 11	88249377	T	C	0.9931	0.1575	0.1783	0.3771	7161	5356	rs11020170	12517
# 1	60320992	A	G	0.9336	0.0605	0.0456	0.1846	26421	442271	rs116406626	468692

column_types <- cols(
  `TopSNP` = col_character(),  # Specify TopSNP as character
  `af` = col_double(),
  `beta` = col_double(),
  `betase` = col_double(),
  `-Log10p` = col_double(),
  `significance_level` = col_double(),
  `trait` = col_character()
)

# TopSNP,af,beta,betase,-Log10p,significance_level,trait,ACI,BN,BUF,F344,M520,MR,WKY,WN
# 1:256613793,0.883,0.255,0.051,6.306,5,nicsa_day1_activelick,G G,A A,G G,G G,G G,G G,G G,G G
# 1:243673316,0.262,0.185,0.038,6.02,5,nicsa_day1_inactivelick,T T,C C,C C,C C,C C,C C,C C,C C

df.p50.nicsa <- read_csv("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/data/hs_data/phenotype/p50_hao_chen_nicsa.csv", col_types = column_types)
df.p50.nicsa

magma.input.df.p50.nicsa <- df.p50.nicsa %>%
  select(1:7) %>% 
  mutate(chr = str_split_n(TopSNP, ':', 1)) %>% 
  mutate(bp = str_split_n(TopSNP, ':', 2)) %>% 
  select(
    chr,
    bp,
    rs = TopSNP,
    P = `-Log10p`,
  ) %>% 
  mutate(P = 10^(-P))
magma.input.df.p50.nicsa
write_tsv(magma.input.df.p50.nicsa, "/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/data/hs_data/phenotype/magma_input_p50_hao_chen_nicsa.csv")







