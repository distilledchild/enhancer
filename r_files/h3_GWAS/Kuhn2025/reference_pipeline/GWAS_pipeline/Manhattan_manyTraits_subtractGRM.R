args <- commandArgs(trailingOnly = TRUE)

traits <- args[1]
gwas_threshold<-args[2]

gwas_threshold<-as.numeric(gwas_threshold)

exp_dir <- getwd()
project_name <- basename(exp_dir)


exp_dir<-paste0(exp_dir,"/")



dir_gwas<-paste0(exp_dir,"results/gwas")




heritability_summary=read.csv("results/snp_h2_with_N.csv",header=T,stringsAsFactors=F)



files<-list.files(dir_gwas,full.names=T)



files<-files[grep(paste0(traits,".mlma"),files)]

library("ggrepel", lib.loc="/home/aschitre/R_libs/")
library("ggman", lib.loc="/home/aschitre/R_libs/")





setwd(exp_dir)



readdata <- function(x)
{
  mydata <- read.table(x,stringsAsFactors = F,header=T)
  chr<-strsplit(basename(x)[1],split=".",fixed=T)[[1]][1]
  chr<-gsub("[^0-9.-]","",chr)
  mydata<-mydata[which(mydata$Chr==chr),]
  return(mydata)
}


gwasResults <- lapply(files, readdata)

gwasResults<-do.call("rbind",gwasResults)


  gwasResults$Chr<-as.numeric(gwasResults$Chr)


  snp_heri<-heritability_summary$snp_heritability[heritability_summary$trait==traits]
  SE<-heritability_summary$SE[heritability_summary$trait==traits]
  sample_size<-heritability_summary$N[heritability_summary$trait==traits]
  snp_heri<-round(snp_heri,digits = 3)
  SE<-round(SE, digits = 3)
  sample_size<-prettyNum(sample_size,big.mark = ",")

  plot_title<-paste0(traits,"\nN= ",sample_size,"\tSNP_h2= ",snp_heri," SE= ",SE,"\nSNPs= ",prettyNum(nrow(gwasResults),big.mark = ","))


  png(paste0("results/manhattan_plots/",traits,".png"),width=2000, height=700)
  p1=ggman(gwasResults, snp = "SNP", bp = "bp", chrom = "Chr", pvalue = "p",relative.positions = T,sigLine = gwas_threshold,pointSize = 1.0)
  print(p1 + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),title =element_text(size=20),axis.text.x = element_text(color='black',size=20),axis.text.y = element_text(color='black',size=20),plot.title = element_text(hjust = 0.5))+ ggtitle(plot_title) + scale_colour_manual(values = c("cornflowerblue", "darkblue")))
  dev.off()
