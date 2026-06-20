
#get PBS arg for PBS directory

setwd("/projects/ps-palmer/apurva/p50_david_dietz_2020/")

load("data/residuals/residuals.RData")
load("data/raw_data/traits.RData")


#write GCTA pheno files
gwas$famid<-gwas$rfid

for(i in 1:length(traits)){

  write.table(gwas[,c("rfid","famid",traits[i])],file=paste0("data/pheno/",traits[i],".txt"),row.names = F,col.names = F,quote=F)


}



#also write a fam file to subset the genotype data

write.table(gwas[,c("rfid","famid")],"genotypes/rfids_n2490.txt",row.names=F,col.names=F,quote=F,sep="\t")


rm(list = ls())
