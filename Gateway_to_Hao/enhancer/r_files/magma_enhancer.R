options(stingAsFactors=F)
library(GenomicRanges)
library(biomaRt)
library(dplyr)
library(tidyverse)
library(remotes)
library(xml2)
library(openxlsx)
library(rvest)

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
# setwd('/home/panjun/dropbox/Gateway_to_Hao/enhancer/r_files')
source(file.path('/home/panjun/dropbox/Gateway_to_Hao/project_common_code/', 'variables.R'))
source(file.path('/home/panjun/dropbox/Gateway_to_Hao/project_common_code/', 'funcs.R'))

# Mac
setwd('/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files')
source(file.path('/Users/PanjunKim/dropbox/Gateway_to_Hao/project_common_code/', 'variables.R'))
source(file.path('/Users/PanjunKim/dropbox/Gateway_to_Hao/project_common_code/', 'funcs.R'))
getwd()

#####################################
# HMAGMA annotation file (2nd) - biomart version // 3rd - biomart & protein-coding genes ONLY
#####################################
# checking ensembl dataset ver.
mart <- useEnsembl("ensembl")
datasets <- listDatasets(mart)
print(datasets) # rnorvegicus_gene_ensembl Rat genes (mRatBN7.2) mRatBN7.2

# step 1: preparing data
### 1-1. read in exonic and promoter coordinate data
# exon data from Biomart
mart <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
exon.from.biomart <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 
                                  'ensembl_exon_id', 'ensembl_transcript_id', 'transcript_biotype',
                                  'chromosome_name', 'exon_chrom_start', 'exon_chrom_end', 'strand'),
                   mart = mart) %>% 
  filter(transcript_biotype == 'protein_coding') %>% 
  filter(chromosome_name %in% c(1:20, "X"))

exon.from.biomart # 506271 rows protein-coding genes only / 525823 rows/ 30387 distinct gene_id
exon.from.biomart %>% distinct(ensembl_gene_id) # 23019 protein-coding genes only / 30387 distinct gene_id

# promoter data from Biomart (transcript)
promoter.from.biomart <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 
                                                'ensembl_transcript_id', 'chromosome_name', 'transcript_start', 'transcript_end', 'strand', 'transcript_biotype'),
                         mart = mart) %>% 
  filter(chromosome_name %in% c(1:20, "X")) %>% 
  filter(transcript_biotype == 'protein_coding') %>% 
  mutate(promoter_start = ifelse(strand == 1, transcript_start - 2000, transcript_end - 500),
         promoter_end = ifelse(strand == 1, transcript_start + 500, transcript_end + 2000))

promoter.from.biomart %>% count(transcript_biotype)

# transcript_biotype     n
# 1             IG_V_gene    38
# 2             TR_C_gene     2
# 3             TR_J_gene     2
# 4             TR_V_gene     1
# 5                 Y_RNA    18
# 6                lncRNA  4077
# 7                 miRNA   444
# 8              misc_RNA    27
# 9  processed_pseudogene   189
# 10       protein_coding 45792
# 11           pseudogene   722
# 12                 rRNA   177
# 13             ribozyme    36
# 14               scaRNA    37
# 15                snRNA  1508
# 16               snoRNA  1706
# 17            vault_RNA     1

promoter.from.biomart # 45792 rows (protein-coding genes) / 54777 rows/ 30387 distinct gene_id
promoter.from.biomart %>% distinct(ensembl_gene_id) # 23019/ 30387

promoter.from.biomart %>% 
  filter(transcript_biotype == 'protein_coding') %>% 
  distinct(ensembl_gene_id) # 23019



promoter.from.biomart %>% 
  filter(transcript_biotype == 'protein_coding') %>% 
  filter(ensembl_gene_id %in% c('ENSRNOG00000052674', # the only gene_id which doesn't code 
                                'ENSRNOG00000008994', 
                                'ENSRNOG00000012934', 
                                'ENSRNOG00000048161', 
                                'ENSRNOG00000016281', 
                                'ENSRNOG00000069004')) %>% 
  distinct(ensembl_gene_id)
# 1 ENSRNOG00000048161
# 2 ENSRNOG00000008994
# 3 ENSRNOG00000012934
# 4 ENSRNOG00000016281
# 5 ENSRNOG00000069004

# test 1: ALL three passed
promoter.from.biomart %>% filter(promoter_start < 0 | promoter_end < 0) # 0 row
promoter.from.biomart %>% filter(is.na(promoter_start) | is.na(promoter_end)) # 0 row
promoter.from.biomart %>% filter(promoter_start > promoter_end) # 0 row, data is alwasy 5' -> 3', so start coord should be alway smaller than end

# gene list from exon and promoter
gene.list.from.biomart.exon <- exon.from.biomart %>% distinct(ensembl_gene_id) %>% pull(ensembl_gene_id)
gene.list.from.biomart.promoter<- promoter.from.biomart %>% distinct(ensembl_gene_id) %>% pull(ensembl_gene_id)

gene.list.from.biomart.exon # 23019 (protein-coding) / 30387
gene.list.from.biomart.promoter # 23019 (protein-coding) / 30387

common.genes.list.from.biomart.exon.promoter <- intersect(gene.list.from.biomart.exon, gene.list.from.biomart.promoter)
common.genes.list.from.biomart.exon.promoter # 23019(protein-coding) /30387

# 1-2. Create a GRanges object for exon and promoter definitions
exonranges <- GRanges(exon.from.biomart$chromosome_name, 
                      IRanges(exon.from.biomart$exon_chrom_start, exon.from.biomart$exon_chrom_end), 
                      strand = exon.from.biomart$strand,
                      gene = exon.from.biomart$ensembl_gene_id)
promoterranges <- GRanges(promoter.from.biomart$chromosome_name, 
                          IRanges(promoter.from.biomart$promoter_start, promoter.from.biomart$promoter_end), 
                          strand = promoter.from.biomart$strand,
                          gene = promoter.from.biomart$ensembl_gene_id)

exonranges
promoterranges

save(exonranges, promoterranges, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_magma_enhancer/0924/exon_promoranges.rda")
save(exonranges, promoterranges, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_magma_enhancer/0926/exon_promoranges.rda") # protein-coding genes only
load("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_magma_enhancer/0924/exon_promoranges.rda")
load("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_magma_enhancer/0926/exon_promoranges.rda") # protein-coding genes only

# 1-3. reading SNP data and creating GRanges for the data
df.init.snps <- read.table("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/data/hs_data/v4/HS_genotypes_v4.bim")
df.init.snps # 7358643
dim(df.init.snps) 
snps <- df.init.snps %>% 
  dplyr::select(1, 2, 4) %>% 
  dplyr::rename(chr = V1, rsid = V2, Position = V4) %>% 
  filter(chr %in% c(1:20, "X"))
snps # 7358388 (1-20, X only) / 7358643 (all) 
# rsid = chr:position
snps <- GRanges(snps$chr, IRanges(snps$Position, snps$Position), rsid=snps$rsid)

save(snps, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_magma_enhancer/0924/snps.rda")
save(snps, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_magma_enhancer/0926/snps.rda") # same with that in 0924. Just added into the folder of 0926
load ("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_magma_enhancer/0924/snps.rda")

# step 2: assigning SNPs to genes by overlapping SNPs with exons and promoters
# 2-1. snps and exonranges
olap <- findOverlaps(snps, exonranges)
snpexon <- snps[queryHits(olap)]
mcols(snpexon) <- cbind(mcols(snpexon), mcols(exonranges[subjectHits(olap)]))
# snpexon <- snpexon[seqnames(snpexon)!="chrX"]
snpexon # 212799 protein-coding genes only / 234045 (exon snp) / 7358643 (1-20, X only)
# > snpexon
# GRanges object with 234045 ranges and 2 metadata columns:
#   seqnames    ranges strand |        rsid               gene
# <Rle> <IRanges>  <Rle> | <character>        <character>
#   [1]        1   1222578      * |   1:1222578 ENSRNOG00000069767
#   [2]        1   1222578      * |   1:1222578 ENSRNOG00000069767

# 2-2. snps and promoterranges
olap <- findOverlaps(snps, promoterranges)
snpro <- snps[queryHits(olap)]
mcols(snpro) <- cbind(mcols(snpro), mcols(promoterranges[subjectHits(olap)]))
# snpro <- snpro[seqnames(snpro)!="chrX"]
snpro # 261212 (protein-coding genes ONLY) / 316978 (promoter snp) / 7358643 (1-20, X only)
# > snpro
# GRanges object with 316978 ranges and 2 metadata columns:
#   seqnames    ranges strand |        rsid               gene
# <Rle> <IRanges>  <Rle> | <character>        <character>
# [1]        1   1212385      * |   1:1212385 ENSRNOG00000069767
# [2]        1   1306242      * |   1:1306242 ENSRNOG00000050129
# [3]        1   2046717      * |   1:2046717 ENSRNOG00000055877

save(snpro, snpexon, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_magma_enhancer/0924/snp_locating_in_exon_promoter_transcript_level.rda")
save(snpro, snpexon, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_magma_enhancer/0926/snp_locating_in_exon_promoter_transcript_level.rda") # protein-coding genes only 
load("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_magma_enhancer/0924/snp_locating_in_exon_promoter_transcript_level.rda")

# step3: assigning UNMAPPED SNPs to genes on the basis of Hi-C interaction data
# 3-1. identifying UNMAPPED SNPs
snpranges <- snps[!(snps$rsid %in% snpexon$rsid),] # 7240036 protein-coding genes only <> 7224824 <> 133564
snpranges
snpranges <- snpranges[!(snpranges$rsid %in% snpro$rsid),]
snpranges # snprages: snps UNMAPPED with exon and promoter, 7090277 protein-coding genes only / 7036577

save(snpranges, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_magma_enhancer/0924/non_exonic_promoter_snp.rda")
save(snpranges, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_magma_enhancer/0926/non_exonic_promoter_snp.rda")
load("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_magma_enhancer/0924/non_exonic_promoter_snp.rda")

# 3-2. reading Hi-C data and creating GRanges for the data
# load loops filtered
load("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/figures/0909/df_final_loop.rda")
df.final.loop

hic.int1 <- df.final.loop %>%
  separate(loop.id, into = c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "resolution"), sep = "_") %>%
  mutate(across(c(start1, end1, start2, end2), as.numeric)) %>% 
  mutate(chrom1 = str_remove(chrom1, "chr"),
         chrom2 = str_remove(chrom2, "chr"))

hic.int2 <- hic.int1 %>% 
  dplyr::select(c(4:6), c(1:3), 7, 8) %>% 
  dplyr::rename(chrom1 = chrom2, start1 = start2, end1 = end2, chrom2 = chrom1, start2 = start1, end2 = end1)
hic.int2 # 12493

hic.comb <- rbind(hic.int1, hic.int2)
hic.comb # 24986 (12493 * 2)

# creating GRange for hic
hicranges <- GRanges(hic.comb$chrom1, IRanges(as.numeric(hic.comb$start1), as.numeric(hic.comb$end1)),
                     int1 = hic.comb$start2, int2 = hic.comb$end2)

# identifying promoter-anchored interactions (overlapping loop anchor1 with promoters)
## |--E1--|--------------------------|--E2--|
## |-P-| : E1 means gene area // E2 far away snp area with the gene of P

olap <- findOverlaps(hicranges, promoterranges)
generanges <- hicranges[queryHits(olap)]
mcols(generanges) <- cbind(mcols(hicranges[queryHits(olap)]), mcols(promoterranges[subjectHits(olap)]))
generanges # 24172 # seqnames              ranges strand       int1      int2               gene

genebed <- data.frame(chr=seqnames(generanges), 
                      snp.start = generanges$int1, snp.end = generanges$int2,
                      gene.start = start(generanges), gene.end = start(generanges) + width(generanges) - 1, ensg = generanges$gene) # gene for the promoter
genebed # 24172
genebed <- unique(genebed)
genebed # 14280

genesnpranges <- GRanges(genebed$chr, IRanges(genebed$snp.start, genebed$snp.end), ensg = genebed$ensg) # E2 area

olap <- findOverlaps(snpranges, genesnpranges) # snpranges NOT MAPPED in exon or promoter
snpint <- snpranges[queryHits(olap)]
mcols(snpint) <- cbind(mcols(snpranges[queryHits(olap)]), mcols(genesnpranges[subjectHits(olap)]))
snpint # snps in E2 with interaction with E1 gene area, 545754
save(snpint, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_magma_enhancer/0924/Hi-C_transcript_interacting_snp.rda")
save(snpint, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_magma_enhancer/0926/Hi-C_transcript_interacting_snp.rda") # protein-coding genes only

# integrating SNP-gene relationships derived from exons, promoters and Hi-C interaction data.
load("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_magma_enhancer/0924/Hi-C_transcript_interacting_snp.rda")
load("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_magma_enhancer/0924/snp_locating_in_exon_promoter_transcript_level.rda")
snpdat <- data.frame(chr = seqnames(snpint), bp = start(snpint), rsid = snpint$rsid, ensg = snpint$ensg) # genes(ensg) in E1 interacting with snps in E2
snpdat
# > snpdat
# chr       bp       rsid               ensg
# 1     1  1419605  1:1419605 ENSRNOG00000064115
# 2     1  1419605  1:1419605 ENSRNOG00000040300
# 3     1  1804735  1:1804735 ENSRNOG00000063217
# 4     1  1804735  1:1804735 ENSRNOG00000070200
snpexonmat <- unique(data.frame(rsid = snpexon$rsid, ensg = snpexon$gene))
snpexonmat # snpexonmat: 120624 protein-coding genes only / 36469/snpexon: 234045
# > snpexonmat
# rsid               ensg
# 1    1:1222578 ENSRNOG00000069767
# 4    1:1223973 ENSRNOG00000069767
# 7    1:1353406 ENSRNOG00000067237
snpromat <- unique(data.frame(rsid = snpro$rsid, ensg = snpro$gene))
snpromat # snpromat: 174696 protein-coding genes only / 224940/snpro: 316978
# > snpexonmat 
# rsid               ensg
# 1    1:1222578 ENSRNOG00000069767
# 4    1:1223973 ENSRNOG00000069767
# 7    1:1353406 ENSRNOG00000067237
# 8    1:1966885 ENSRNOG00000062996
snpcomb <- unique(rbind(snpdat[, 3:4], snpromat, snpexonmat))
snpcomb # 690474 protein-coding genes only / 831959
save(snpcomb, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_magma_enhancer/0924/SNP_to_transcript_comb.rda")
save(snpcomb, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_magma_enhancer/0926/SNP_to_transcript_comb.rda") # protein-coding genes only

# step4: creating annotation file
# 4-1.  Aggregate SNP-gene relationships to generate the variant-gene annotation file compatible with MAGMA
# snpagg <- aggregate(snpcomb, list(snpcomb$ensg), unique)
# snpagg

snpagg <- snpcomb %>%
  group_by(ensg) %>%
  summarise(rsid = paste(unique(rsid), collapse = ", ")) %>%
  ungroup()
snpagg

# 4-2. gene definition file from Biomart
# mart <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
genedef <- getBM(attributes = c('chromosome_name', 'start_position', 'end_position', 'strand', 'ensembl_gene_id'), # external_gene_name
                              filters = 'ensembl_gene_id',
                              values = common.genes.list.from.biomart.exon.promoter,
                              mart = mart) %>% 
  dplyr::rename(chr = chromosome_name, start = start_position, end = end_position, ensg = ensembl_gene_id) %>% 
  mutate(index = str_c(chr, start, end, sep=":")) %>% 
  mutate(strand = ifelse(strand == 1, "+", ifelse(strand == -1, "-", strand))) %>% 
  dplyr::select(index, everything()) %>% 
  dplyr::select(-ensg, everything(), ensg)
  
genedef # 23019 (values = common.genes.list.from.biomart.exon.promoter) protein-coding genes only / 30387
# index chr     start       end strand               ensg
# 1   3:100064979:100083289   3 100064979 100083289      + ENSRNOG00000000009
# 2     4:28276909:28287479   4  28276909  28287479      + ENSRNOG00000000017
# 3   2:102549724:102609327   2 102549724 102609327      + ENSRNOG00000000082
# 4   2:138986471:139006307   2 138986471 139006307      + ENSRNOG00000000091
# 5   2:184401438:184445584   2 184401438 184445584      - ENSRNOG00000000095

# dedup by higher number in ensg digit part
# awk 'BEGIN {FS=OFS="\t"} {split($NF,a,"ENSRNOG"); print $0, a[2]}' /Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_magma_enhancer/0924/genedef_ensembl.tsv | sort -k1,1 -k6,6nr | awk '!seen[$1]++ {print $0}' | cut -f1-6 | cut -f1 | sort | uniq -d
genedef.for.MAGMA.annot.input <- genedef %>% dplyr::select(ensg, chr, start, end, strand)

write.table(genedef, file = "/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_magma_enhancer/0924/genedef_ensembl.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(genedef.for.MAGMA.annot.input , file = "/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_magma_enhancer/0926/genedef_ensembl.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

genedef <- read.delim("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_magma_enhancer/0924/genedef_ensembl.tsv", 
                      sep = "\t", header = FALSE, 
                      col.names = c("index", "chr", "start", "end", "strand", "ensg"))
genedef # 23019 / 30387

# genedef %>% filter(index %in% c('4:146386956:146429990', '4:146510246:146521590', '4:146511607:146511797'))
'4:146511607:146511797' : non-protein-coding snp

# 4-3. attaching the index column from Step 24 to the variant-gene annotation file (snpagg).
snpagg$index <- genedef[match(snpagg$ensg, genedef$ensg), "index"]
snpagg # 21442 protein-coding genes only / 27922
# snpagg <- snpagg[!is.na(snpagg$index),]
# snpagg # 27922

snpaggconv <- snpagg %>% 
  filter(!is.na(index)) %>% 
  dplyr::select(ensg, index, rsid)
  # mutate(ensg = index) %>% 
snpaggconv # 21442 protein-coding genes only / 27922
# > snpaggconv
# # A tibble: 27,922 × 3
# ensg               index                 rsid
write.table(snpaggconv, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_magma_enhancer/0924/SNP_aggregate_transcript.txt", quote=F, row.names=F, col.names=FALSE, sep="\t") # change the name of the file
write.table(snpaggconv, file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_magma_enhancer/0926/SNP_aggregate_transcript.txt", quote=F, row.names=F, col.names=FALSE, sep="\t") # change the name of the file
system("sed -e 's/, /\t/g' < /Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_magma_enhancer/0924/SNP_aggregate_transcript.txt > /Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_magma_enhancer/0924/hs.hic.annot.sep.2024.2nd", wait = TRUE)
system("sed -e 's/, /\t/g' < /Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_magma_enhancer/0926/SNP_aggregate_transcript.txt > /Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/from_magma_enhancer/0926/hs.hic.annot.sep.2024.3rd", wait = TRUE)

