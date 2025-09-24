# variable for project directory
# root.dir <- "C:/Users/panju/Dropbox (UTHSC GGI)/WLIWMIB6Research" # windows dropbox
root.windows.dir <- "C:/Users/panju/Insync/pkim11@uthsc.edu/Dropbox" # windows insync
root.ubuntu.dir <- "/home/panjun/Dropbox" # linux insync
root.pp.dir <- "C:/Users/panju/Dropbox (UTHSC GGI)/WLIWMIB6ResearchPP"
root.new.dir <- "/media/panjun/fa7f5b4f-c57c-47d9-816d-eca52685f2f2/Dropbox (UTHSC GGI)/WLIWMIB6Research"
root.hao.dir <- "C:\\Users\\panju\\Insync\\pkim11@uthsc.edu\\Dropbox\\Gateway_to_Hao"

# path
path.mine.dir <- "C:/Users/panju" # mine
root.dropbox.dir <- "Dropbox (UTHSC GGI)" # dropbox root dir

# RCC oxycodone
rcc.oxy <- "RCC_oxycodone"
oxy.project.ref.dir <- "1_oxycodone/reference_files"
oxy.project.h.dir <- "1_oxycodone/2021/working files/output_files/2021_H" # 2nd dataset with BCFtools, HIGH, LOW, PASS, PART, FAIL system
oxy.project.i.dir <- "1_oxycodone/2021/working files/output_files/2021_I" # 1st dataset with BCFtools, HIGH, LOW, PASS, PART, FAIL system
oxy.project.j.dir <- "1_oxycodone/2021/working files/output_files/2021_J" # 
oxy.project.k.dir <- "1_oxycodone/2021/working files/output_files/2021_K" # 1st dataset wtih DeepVariants
oxy.project.l.dir <- "1_oxycodone/2021/working files/output_files/2021_L" # Data analysis with MERGED DATASET(1st and 2nd), Vaild, Fail, Polymorphic, non-polymorphic system
oxy.project.m.dir <- "1_oxycodone/2021/working files/output_files/2021_M" # ideogram and setting markers with their coordinates across chromosome
oxy.project.n.dir <- "1_oxycodone/2021/working files/output_files/2021_N" # based on 2021_M, picking additional markers using ensembl API
oxy.project.o.dir <- "1_oxycodone/2021/working files/output_files/2021_O" # variant calling merged files with deepvariant // and comparison with results from BCFtools
oxy.project.p.dir <- "1_oxycodone/2021/working files/output_files/2021_P" # After GATK RealignerTargetCreator, the comparison between BCFTools vs Deepvariant
oxy.project.q.dir <- "1_oxycodone/2021/working files/output_files/2021_Q" # making query for the breeding table, Vaild, Fail, Polymorphic, non-polymorphic system
oxy.project.r.dir <- "1_oxycodone/2021/working files/output_files/2021_R" # Separating picking new markers from 2021_N. 2021_N has two parts: checking correlation and pick markers using ensembl API && check the validity for the primers using primer3 and MFEprimer
oxy.project.s.dir <- "1_oxycodone/2021/working files/output_files/2021_S" # Sequencing data processing under new additional marker, Vaild, Fail, Polymorphic, non-polymorphic system
oxy.project.t.dir <- "1_oxycodone/2021/working files/output_files/2021_T" # tail vs spleen data

oxy.project.2022.a.dir <- "1_oxycodone/2022/2021_A"
oxy.project.2022.c.dir <- "Gateway_to_Hao/oxycodone/2022_C"

# RCC stroke
rcc.mice <- "RCC_mice"
str.project.ref.dir <- "2_stroke/reference_files"
str.project.20d.dir <- "2_stroke/2020/working files/output files/2020_D" # checking genotyping on targed resequencing of P & F1 for marker selection
str.project.21d.dir <- "2_stroke/2021/working files/output files/2021_D" # rechecking using new codes for 2020_D
str.project.f.dir <- "2_stroke/2021/working files/output files/2021_F" # genotyping ByJ, CRL, and SLC
str.project.g.dir <- "2_stroke/2021/working files/output files/2021_G" # working on the first dataset of 75 samples in Jan 2021 (additional work for 1st, repeat, and merged dataset)
str.project.h.dir <- "2_stroke/2021/working files/output files/2021_H" # for making geno file using merged dataset. and data QC : almost deprecated
str.project.i.dir <- "2_stroke/2021/working files/output files/2021_I" # making geno file and QC step (new phenotypes added) (mainly used more than 2021_H)
str.project.j.dir <- "2_stroke/2021/working files/output files/2021_J" # for making geno file using MiniMUGA dataset (TEMP DONE WITH seasonal project task4)
str.project.k.dir <- "2_stroke/2021/working files/output files/2021_K" # rqtl2 conversion 
str.project.l.dir <- "2_stroke/2021/working files/output files/2021_L" # joining annotation data and vcf files
str.project.m.dir <- "2_stroke/2021/working files/output files/2021_M" # data processing on Dr.Palmer's data
str.project.n.dir <- "2_stroke/2021/working files/output files/2021_N" # QTL analysis: from QC to analysis

# hic
hic.project.c.dir <- "5_hic/2021/output/2021_C" # MAGMA 
hic.project.2021.c.dir <- "/media/panjun/fa7f5b4f-c57c-47d9-816d-eca52685f2f2/Dropbox (UTHSC GGI)/Gateway_to_Hao/hic/2021_C/whole_genome_alignment"
hic.project.2022.a.dir <- "/media/panjun/fa7f5b4f-c57c-47d9-816d-eca52685f2f2/Dropbox (UTHSC GGI)/Gateway_to_Hao/hic/2021_C/matrix_eqtl/data"
hic.project.2022.f.dir <- "Gateway_to_Hao/hic/2022_F_CTC"
hic.project.2022.g.dir <- "Gateway_to_Hao/hic/2022_G"
hic.project.2023.a.dir <- "Gateway_to_Hao/hic/2023A"

# pangenome
pange.project.2022.a.dir <- "Gateway_to_Hao/pangenome.windows.10142022/2022_A" # primers
pange.project.2022.b.dir <- "Gateway_to_Hao/pangenome.windows.10142022/2022_B" # control set for primers in 2022_A
pange.project.2022.c.dir <- "Gateway_to_Hao/pangenome.windows.10142022/2022_C" # primers in chr1
pange.project.2022.d.dir <- "Gateway_to_Hao/pangenome.windows.10142022/2022_D" # calc for markers

# crispr
crispr.project.2022.a.dir <- "Gateway_to_Hao/crispr/2022_A" # deepvariant with 12 samples

# vacation
vc.project.2022.a.dir <- "Gateway_to_Hao/vacation/2022_A"

# seasonal task
season.project.ref.dir <- "references"
season.project.ref.minimuga.control.dir <- "GenotypingTargetedResequencing_MarkerPanelAndSamples/MiniMuga/ControlSamples"
season.project.ref.minimuga.sample.dir <- "MiniMUGA_Sample_Data/sample_data"
season.project.21.summer.dir <- "0_seasonal_project/21_summer_vactaion" # 4 sub-tasks

# variable for initial file
vcf.dir.name <- 'vcf'

# variable for getting sequence region
species.mouse = "mouse"
species.rat = "rat"
