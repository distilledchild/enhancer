library("ComplexUpset")
library("tidyverse")
library("GenomicRanges")
library("ggplot2")
library("patchwork")
library("devtools")
library("remotes")
library("CTCF")
library("AnnotationHub")
library("plyranges")
library("tracktables")
library("universalmotif")

# Linux
# setwd('~/Desktop/temp/enhancer/dropbox_enhancer_doosan')
setwd('/home/pkim/dropbox/Gateway_to_Hao/enhancer/r_files')
setwd('/home/pkim/playground/research_uthsc/enhancer/Gateway_to_Hao/enhancer/r_files')
source(file.path('/home/pkim/Dropbox/Gateway_to_Hao/project_common_code/', 'variables.R'))
source(file.path('/home/pkim/Dropbox/Gateway_to_Hao/project_common_code/', 'funcs.R'))

# Mac
setwd('/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/')
source(file.path('/Users/PanjunKim/dropbox/Gateway_to_Hao/project_common_code/', 'variables.R'))
source(file.path('/Users/PanjunKim/dropbox/Gateway_to_Hao/project_common_code/', 'funcs.R'))

getwd()

#######################################################
# convert START
#######################################################
# jaspar2022- meme provided
#######################################################
#######################################################
# hocomoco-by universalmotif
#######################################################
hocomoco.CTCF_HUMAN.H11MO.0.A.pwm <- read_matrix("/home/pkim/Desktop/temp/enhancer/dropbox_enhancer_doosan/data/ctcf/pwms/hocomoco/CTCF_HUMAN.H11MO.0.A.pwm",
                                                 type = 'PWM',
                                                 alphabet = 'DNA',
                                                 headers = '>',
                                                 positions = 'rows',
                                                 rownames = FALSE)
hocomoco.CTCF_MOUSE.H11MO.0.A.pwm <- read_matrix("/home/pkim/Desktop/temp/enhancer/dropbox_enhancer_doosan/data/ctcf/pwms/hocomoco/CTCF_MOUSE.H11MO.0.A.pwm",
                                                 type = 'PWM',
                                                 alphabet = 'DNA',
                                                 headers = '>',
                                                 positions = 'rows',
                                                 rownames = FALSE)
hocomoco.CTCF_HUMAN.H11MO.0.A.pwm
hocomoco.CTCF_MOUSE.H11MO.0.A.pwm


write_meme(hocomoco.CTCF_HUMAN.H11MO.0.A.pwm,
           "./data/ctcf/pwms/hocomoco/hocomoco.CTCF_HUMAN.H11MO.0.A.universalmotif.meme",
           version = 5, overwrite = FALSE,
           append = FALSE)
write_meme(hocomoco.CTCF_MOUSE.H11MO.0.A.pwm,
           "./data/ctcf/pwms/hocomoco/hocomoco.CTCF_MOUSE.H11MO.0.A..universalmotif.meme",
           version = 5, overwrite = FALSE,
           append = FALSE)

#######################################################
# swissregulon
#######################################################
#######################################################
# jolma-by transfact-like
#######################################################
#######################################################
# CTCFBSDB-by transfact-like
#######################################################
#######################################################
# CIS-BP-by universalmotif
#######################################################
CisBP.homo.sapiens <- read_cisbp("./data/ctcf/pwms/utoronto/CisBP_2023_10_16_11_40_am_homo_sapiens/CisBP_2023_10_16_11_40_am_homo_sapiens.txt")
CisBP.mouse <- read_cisbp("./data/ctcf/pwms/utoronto/CisBP_2023_10_16_11_41_am_mouse/CisBP_2023_10_16_11_41_am_mouse.txt")


write_meme(CisBP.homo.sapiens,
           "./data/ctcf/pwms/utoronto/CisBP_2023_10_16_11_40_am_homo_sapiens/CisBP_2023_10_16_11_40_am_homo_sapiens.universalmotif.meme",
           version = 5, overwrite = TRUE,
           append = FALSE)
write_meme(CisBP.mouse,
           "./data/ctcf/pwms/utoronto/CisBP_2023_10_16_11_41_am_mouse/CisBP_2023_10_16_11_41_am_mouse.universalmotif.meme",
           version = 5, overwrite = TRUE,
           append = FALSE)

#######################################################
# convert DONE
#######################################################




# 1. collecting position weight matrices of CTCF binding motifs and defining strand-oriented CTCF binding sites in human and mouse genome

suppressMessages(library(AnnotationHub))
ah <- AnnotationHub()
query_data <- subset(ah, preparerclass == "CTCF")

# Explore the AnnotationHub object
query_data

# Get the list of data providers
query_data$dataprovider %>% table()
# CIS-BP     CTCFBSDB 2.0 ENCODE SCREEN v3     HOCOMOCO v11 
# 6               12                2                6 
# JASPAR 2022       Jolma 2013     SwissRegulon 
# 13                6                6 

query_data$species

subset(query_data, 
       species == "Homo sapiens" & 
         genome == "hg38" & 
         dataprovider == "JASPAR 2022")

subset(query_data, dataprovider == "SwissRegulon")
query_data[["AH104722"]]

CTCF_hg38_all <- query_data[["AH104727"]]
CTCF_hg38_all

# hg38.MA0139.1
CTCF_hg38 <- query_data[["AH104729"]]
#> loading from cache
CTCF_hg38

# To sort GRanges objects and keep standard chromsomes
suppressMessages(library(plyranges))
CTCF_hg38_all <- CTCF_hg38_all %>% keepStandardChromosomes() %>% sort()
CTCF_hg38 <- CTCF_hg38 %>% keepStandardChromosomes() %>% sort()

# Note that rtracklayer::import and rtracklayer::export perform unexplained
# start coordinate conversion, likely related to 0- and 1-based coordinate
# system. We recommend converting GRanges to a data frame and save tab-separated
write.table(CTCF_hg38_all %>% sort() %>% as.data.frame(), 
            file = "CTCF_hg38_all.bed",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(CTCF_hg38 %>% sort() %>% as.data.frame(), 
            file = "CTCF_hg38.bed",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

##############################################################################
# Create an IGV XML session file out of the saved BED files using the tracktables package
##############################################################################
# Sample sheet metadata
SampleSheet <- data.frame(SampleName = c("CTCF all", "CTCF MA0139.1"),
                          Description = c("All CTCF matrices from JASPAR2022",
                                          "MA0139.1 CTCF matrix from JASPAR2022"))
# File sheet linking files with sample names
FileSheet <- data.frame(SampleName = c("CTCF all", "CTCF MA0139.1"),
                        bigwig = c(NA, NA),
                        interval = c("CTCF_hg38_all.bed", "CTCF_hg38.bed"),
                        bam = c(NA, NA))
# Creating an IGV session XML file
MakeIGVSession(SampleSheet, FileSheet, 
               igvdirectory = getwd(), "CTCF_from_JASPAR2022", "hg38")

##############################################################################
# To filter the GRanges object and keep high-confidence CTCF sites
##############################################################################

# Check length before filtering
print(paste("Number of CTCF motifs at the default 1e-4 threshold:", length(CTCF_hg38)))
#> [1] "Number of CTCF motifs at the default 1e-4 threshold: 887980"
# Filter and check length after filtering
CTCF_hg38_filtered <- CTCF_hg38 %>% plyranges::filter(pvalue < 1e-6)
print(paste("Number of CTCF motifs at the 1e-6 threshold:", length(CTCF_hg38_filtered)))
#> [1] "Number of CTCF motifs at the 1e-6 threshold: 21671"
# Similarly, filter
CTCF_hg38_all_filtered <- CTCF_hg38_all %>% plyranges::filter(pvalue < 1e-6)

# Given some databases provide multiple CTCF PWMs, 
# one CTCF site may be detected multiple times resulting in overlapping CTCF sites. 
# For example, the proportion of overlapping CTCF sites in the CTCF_hg38_all_filtered object containing CTCF sites detected by three matrices nearly 40%
# Proportion of overlapping enrtries
tmp <- findOverlaps(CTCF_hg38_all, CTCF_hg38_all)
prop_overlap <- sort(table(queryHits(tmp)) %>% table(), decreasing = TRUE)
sum(prop_overlap[which(names(prop_overlap) != "1")]) / length(CTCF_hg38_all)

# The proportion of overlapping CTCF sites in the CTCF_hg38_filtered object containing CTCF sites detected by the MA0139.1 matrix 
# is less than 2.5%
tmp <- findOverlaps(CTCF_hg38, CTCF_hg38)
prop_overlap <- sort(table(queryHits(tmp)) %>% table(), decreasing = TRUE)
sum(prop_overlap[which(names(prop_overlap) != "1")]) / length(CTCF_hg38)
#> [1] 0.0233485

# Reducing them (merging overlapping CTCF sites), combined with 1E-6 cutoff filtering, 
# yields the number of CTCF sites comparable to previously reported.
print(paste("Number of CTCF_hg38 motifs at the 1e-6 threshold AND reduced:", length(CTCF_hg38_filtered %>% reduce())))
#> [1] "Number of CTCF_hg38 motifs at the 1e-6 threshold AND reduced: 21652"
print(paste("Number of CTCF_hg38_all motifs at the 1e-6 threshold AND reduced:", length(CTCF_hg38_all_filtered %>% reduce())))
#> [1] "Number of CTCF_hg38_all motifs at the 1e-6 threshold AND reduced: 63572"

# However, regulatory elements with CTCF proteins co-occupying adjacent/overlapping CTCF binding motifs were shown to be functionally and structurally different from those with single CTCF motifs. 
# We provide non-reduced CTCF data and advise considering overlap of CTCF sites depending on the study’s goal.

### liftOver of CTCF coordinates
# As genome assemblies for model organisms continue to improve, CTCF sites for previous genome assemblies become obsolete. Typically, the actual genome sequence changes little, leading to changes in genomic coordinates. 
# The liftOver method allows for conversion of genomic coordinates between genome assemblies.
# 
# Some carefully curated CTCF sites are available only for older genome assemblies. Examples include the data from CTCFBSDB, available for hg18 and mm8 genome assemblies.
# 
# To investigate whether liftOver of CTCF sites from older genome assemblies is a viable option, we tested for overlap between CTCF sites directly detected in specific genome assemblies with those lifted over. 
# We detected CTCF sites using the MA0139.1 PWM from JASPAR 2022 database in hg18, hg19, hg38, and T2T genome assemblies and converted their genomic coordinates using the corresponding liftOver chains (download_liftOver.sh and convert_liftOver.sh scripts). 
# We observed high Jaccard overlap among CTCF sites detected in the original genome assemblies or lifted over.

# Jaccard overlaps among CTCF binding sites detected in the original and liftOver human genome assemblies. 
# CTCF sites were detected using JASPAR 2022 MA0139.1 PWM. The correlogram was clustered using Euclidean distance and Ward.D clustering. 
# White-red gradient indicate low-to-high Jaccard overlaps. Jaccard values are shown in the corresponding cells.

# Our results suggest that liftOver is a viable alternative to obtain CTCF genomic annotations for different genome assemblies. 
# We provide CTCFBSDB data converted to hg19 and hg38 genome assemblies.

##################################################################################
# CTCF Position Weight Matrices
##################################################################################
# Jaspar2022		
# MA0139.1	19	https://jaspar2022.genereg.net/matrix/MA0139.1
# MA1929.1	34	http://jaspar2022.genereg.net/matrix/MA1929.1
# MA1930.1	35	http://jaspar2022.genereg.net/matrix/MA1930.1

# HOCOMOCO v11		
# CTCF_HUMAN.H11MO.0.A	19	http://hocomoco.autosome.ru/motif/CTCF_HUMAN.H11MO.0.A
# CTCF_MOUSE.H11MO.0.A	20	http://hocomoco.autosome.ru/motif/CTCF_MOUSE.H11MO.0.A

# SwissRegulon	
# CTCF.p2	20	https://swissregulon.unibas.ch/wm/?wm=CTCF.p2&org=hg18

# Jolma 2013		

# CTCF_full	17	http://floresta.eead.csic.es/footprintdb/index.php?db=HumanTF:1.0&motif=CTCF_full

# CTCFBSDB		
# EMBL_M1, EMBL_M2, MIT_LM2, MIT_LM7, MIT_LM23, REN_20	9-20	https://insulatordb.uthsc.edu/download/CTCFBSDB_PWM.mat

# CIS-BP		
# 83 CTCF (Homo sapiens) C2H2 ZF	11-21	http://cisbp.ccbr.utoronto.ca/TFreport.php?searchTF=T094831_2.00
# 2 Ctcf (Mus musculus) C2H2 ZF	15, 20	http://cisbp.ccbr.utoronto.ca/TFreport.php?searchTF=T100985_2.00


##################################################################################
# CTCF predicted and experimental data
# Predefined CTCF binding data. “Database” - source of data; “Number” - number of binding sites; “Assembly” - genome assembly; “URL” - direct link to data download.
##################################################################################
# Database	Number	Assembly	URL

# CTCFBSDB 2.0	NA		
# Predicted human CTCF binding sites	13401	hg18	https://insulatordb.uthsc.edu/download/allcomp.txt.gz
# Predicted mouse CTCF binding sites	5504	mm8	https://insulatordb.uthsc.edu/download/allcomp.txt.gz
# SCREEN ENCODE	NA		
# Human CTCF-bound cCREs	450641	hg38	https://api.wenglab.org/screen_v13/fdownloads/cCREs/GRCh38-CTCF.bed
# Mouse CTCF-bound cCREs	82777	mm10	https://api.wenglab.org/screen_v13/fdownloads/cCREs/mm10-CTCF.bed