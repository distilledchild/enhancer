args <- commandArgs(trailingOnly = TRUE)

exp_dir <- args[1]

expname <- basename(exp_dir)


exp_dir<-paste0(exp_dir,"/")

setwd(exp_dir)

load("data/residuals/residuals.RData")
load("data/residuals/traits.RData")

project_name=basename(tscc_dir)

traitlist=traits

n_traits<-length(traitlist)

base_dir<-exp_dir

all=gwas


dir_heri<-paste0(base_dir,"results/heritability/")
traitlist<-traits

for (i in seq_along(traitlist)){
  command1 <- "/oasis/tscc/scratch/aschitre/gcta/gcta64 --reml --thread-num 8 "

  command2<- paste0("--pheno ",base_dir,"data/pheno/",traitlist[i],".txt ")

  command3<- paste0("--grm ",base_dir,"grm/genotypes ")

  command4<-paste0("--out ",dir_heri,traitlist[i])

  line = paste0(command1,command2,command3,command4)
  system(line)

}





library(data.table)

files<- list.files(path = paste0(dir_heri), pattern = ".hsq", recursive = TRUE,full.names=T)
filenames<-basename(files)
to_gsub<-paste0(expname,"_")

filenames<-gsub(to_gsub,"",filenames)
filenames<-gsub(".hsq","",filenames)


#grep -e 'V(G)/Vp' -e 'Pval'
#"sweet\|lemon"
readdata <- function(x)
{
  h2 <- fread(input=paste0("grep 'V(G)/Vp' ",x),stringsAsFactors = F,header=F)
  h2[,V1:=NULL]
  column_names<-c("snp_heritability","SE")
  setnames(h2, column_names)
  Pval <- fread(input=paste0("grep 'Pval' ",x),stringsAsFactors = F,header=F)
  Pval[,V1:=NULL]
  setnames(Pval, "Pval")
  mydata<-cbind(h2,Pval)
  return(mydata)
}

heritability<- lapply(files, readdata)
names(heritability) <- filenames

for(i in seq_along(filenames))
{

  heritability[[i]]$trait<- filenames[i]
}

heritability_summary <- do.call(rbind,heritability)
write.table(heritability_summary,file=paste0(exp_dir,"results/","heritability.txt"),row.names=F,quote=F)




all_heri=read.table(paste0(exp_dir,"results/heritability.txt"),header=T,stringsAsFactors = F)

h2=all_heri


h2$snp_heritability<-as.numeric(h2$snp_heritability)
h2$SE<-as.numeric(h2$SE)

resid=gwas


df=resid[,h2$trait]

total_N<-colSums(!is.na(df))
total_N=stack(total_N)

colnames(total_N)<-c("N","trait")


final=merge(total_N,h2,by="trait",all=F)
final$trait<-as.character(final$trait)
final$Pval<-NULL


write.csv(final,paste0(exp_dir,"results/snp_h2_with_N.csv"))
