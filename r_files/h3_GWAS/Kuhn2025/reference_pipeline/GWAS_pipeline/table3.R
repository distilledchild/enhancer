library(data.table)
exp_dir <- getwd()
expname <- basename(exp_dir)


exp_dir<-paste0(exp_dir,"/")



project_name=basename(exp_dir)


topsnps=read.csv(paste0(exp_dir,"results/qtls/final_QTL_list.csv"),header=T,stringsAsFactors = F)

topsnps$chr<-sapply(strsplit(topsnps$topsnp,split=":"),`[`,1)
topsnps$pos<-sapply(strsplit(topsnps$topsnp,split=":"),`[`,2)

column_names<-c("snp1","snp2","rsquare")
LD_start_stop<-function(x,trait,chr,pos){
  data<-fread(x,header=T,stringsAsFactors =F,select=c(3,6,7))
  setnames(data, column_names)
  data$dprime<-rep(NA,nrow(data))
  data<-data[which(data$rsquare>=0.6),]

  #split snp1
  data[, c("chr", "pos") := tstrsplit(snp2, ":", fixed=TRUE)]

  data$pos <- as.numeric(data$pos)
  data<-data[order(pos),]
  LD_interval_start<-data$pos[1]
  LD_interval_stop<-data$pos[nrow(data)]
  output<-data.frame(LD_interval_start,LD_interval_stop)
  return(output)
}





out<-data.frame()


for(i in 1:nrow(topsnps)){

  trait<-topsnps$trait[i]
  topsnp<-topsnps$topsnp[i]
  chr<-topsnps$chr[i]
  pos<-topsnps$pos[i]

  system(paste0("/oasis/tscc/scratch/aschitre/software/plink-1.90/plink --bfile ",exp_dir,"genotypes/genotypes --chr ",gsub("chr","",chr)," --nonfounders --r2  --ld-snp ",topsnp," --ld-window 100000 -ld-window-kb 6000 --ld-window-r2 0 --out ",exp_dir,"temp/r2/LZ_temp"),ignore.stdout = T,ignore.stderr = T,wait = T)



  test<-LD_start_stop(paste0(exp_dir,"temp/r2/LZ_temp.ld"),trait = trait,chr=chr,pos=pos)
  out=rbind(out,test)

}






final<-cbind(topsnps,out)

final$LD_interval_size_bp <- final$LD_interval_stop - final$LD_interval_start


column_names<-c("Chr","SNP","A1","A2","Freq","b","se","p")

getassocAfBetaSe<-function(trait,topsnp,chr){
  assoc_name<-paste0("results/gwas/","chr",chr,".",trait,".mlma")
  data  =read.table(assoc_name,header=T,stringsAsFactors=F)
  #V1  V2 V3 V4 V5   V6 V7 V8 V9
  #Chr SNP bp A1 A2 Freq  b se  p

  data=data[,c(1,2,4,5,6,7,8,9)]

  colnames(data)<-column_names

  data<-data[which(data$Chr==chr),]
  rs<-data$SNP[which(data$SNP == topsnp)]
  af<-data$Freq[which(data$SNP == topsnp)]
  beta<-data$b[which(data$SNP == topsnp)]
  se<-data$se[which(data$SNP == topsnp)]
  allele1<-data$A1[which(data$SNP == topsnp)]
  allele2<-data$A2[which(data$SNP == topsnp)]
  p_score<-data$p[which(data$SNP == topsnp)]
  output<-data.frame(trait,rs,af,beta,se,allele1,allele2,p_score)
  return(output)

}

out<-NULL



for(i in 1:nrow(final)){

  trait<-final$trait[i]
  topsnp<-final$topsnp[i]
  chr<-gsub("chr","",final$chr[i])
  test<-getassocAfBetaSe(trait=trait,topsnp=topsnp,chr=chr)
  out=rbind(out,test)
}


out$trait_snp<-NA



colnames(out)[which(colnames(out) %in% "rs")]<-"topsnp"
out$trait_snp = paste0(out$trait,"_",out$topsnp)

final$trait_snp <- paste0(final$trait,"_",final$topsnp)


merged=merge(final,out,by=c("trait","topsnp","trait_snp"),all=F)


topsnps<-merged



founder_GT<-read.table("/projects/ps-palmer/apurva/round2/physiological/effect_plot_founder_GTs/founder_GT_transposedPED.txt",header=F,stringsAsFactors = F,sep="\t")
founder_GT[,c(2:3)] <- NULL
colnames(founder_GT)<-c('chr','pos','ACI','BN','BUF','F344','M520','MR','WN','WKY')
founder_GT$chr <- paste0("chr",founder_GT$chr)
founder_GT$topsnp<- paste0(founder_GT$chr,":",founder_GT$pos)


founder<-founder_GT[founder_GT$chr %in% paste0("chr",seq(1:20)),]
#founder<-founder_GT[founder_GT$chr %in% seq(1:20),]
founder$chr<-NULL
founder$pos<-NULL
#merge founder df with final_summary_HS_GT_pheno_mean by ="snp"
#founder_HS<-merge(final_summary_HS_GT_pheno_mean,founder,by="topsnp",all=F)



#topsnps<-read.csv("/oasis/tscc/scratch/aschitre/round5/pruned/physiological/physiological_kidney_missing_table3_final_missing_kidney.csv",header=T,stingsAsFactors=F)


#merge founder df with final_summary_HS_GT_pheno_mean by ="snp"
founder_HS<-merge(topsnps,founder,by="topsnp",all=F)

founder_HS[,c('ACI','BN','BUF','F344','M520','MR','WN','WKY')]<-as.data.frame(lapply(founder_HS[c('ACI','BN','BUF','F344','M520','MR','WN','WKY')], function(y) gsub("[[:space:]]", "", y)))


#write.csv(founder_HS,paste0(exp_dir,"QTL_locuszoom/",expname,"_table3_withFounder.csv"),quote=T,row.names=F)




raw=founder_HS
#pulling genes from RGD

library(httr)
library(RJSONIO,lib.loc="/home/aschitre/R_libs/")
library(openxlsx,lib.loc="/home/aschitre/R_libs/")

library(purrr,lib.loc="/home/aschitre/R_libs/")



raw$genes<-paste0('https://rest.rgd.mcw.edu/rgdws/genes/',gsub("chr","",raw$chr),'/',raw$LD_interval_start,'/',raw$LD_interval_stop,"/360")
raw$num_of_genes<-NA


for(i in 1:nrow(raw)){


  raw_genes<-fromJSON(raw$genes[i])
  raw$num_of_genes[i]<-length(raw_genes)
  temp<-paste(map_chr(raw_genes,"symbol"),collapse = ",")
  raw$gene_symbol[i]<-temp




}








raw$genes<-NULL




write.xlsx(raw, paste0(exp_dir,"results/qtls/table3_final.xlsx"))
write.csv(raw, paste0(exp_dir,"results/qtls/table3_final.csv"),row.names=F,quote=T)
