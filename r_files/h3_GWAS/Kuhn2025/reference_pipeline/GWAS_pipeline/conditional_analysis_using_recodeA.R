args <- commandArgs(trailingOnly = TRUE)

trait <- args[1]
chr <- args[2]
gwas_threshold<-args[3]

gwas_threshold<-as.numeric(gwas_threshold)

exp_dir <- getwd()


exp_dir<-paste0(exp_dir,"/")



#dosages_file="genotypes/dosages/dosages.txt"


exp_dir<-paste0(exp_dir,"/")

header<-read.table("/projects/ps-palmer/apurva/round8/gcta_assoc_header.txt",header=T,stringsAsFactors=F)



project_name=basename(exp_dir)

library(data.table)


chr <- as.numeric(gsub("chr","",chr))
cov_file=paste0(exp_dir,"temp/",chr,"_",trait,"_","temp.txt")
assoc_gcta=paste0(exp_dir,"temp/",chr,"_",trait,"_","temp_gcta_cond.assoc.txt.mlma")
fam=read.table(paste0(exp_dir,"genotypes/genotypes.fam"),header=F,stringsAsFactors = F)
fam=fam[,c("V1","V2")]
fam$V1=as.character(fam$V1)
fam$V2=as.character(fam$V2)
if (file.exists(cov_file)) {
  #Delete file if it exists
  file.remove(cov_file)
}

if (file.exists(assoc_gcta)) {
  #Delete file if it exists
  file.remove(assoc_gcta)
}



pheno<-paste0(exp_dir,"data/pheno/",trait,".txt")
assoc<-paste0(exp_dir,"results/gwas/chr",chr,".",trait,".mlma")







assoc_for_pval<-list()

topSNP_part2_list<- list()
topSNP_part2<-data.frame()
cov_counter<- 0
loop_run <- TRUE
assoc_counter<-1


while(loop_run){

  if(cov_counter==0){


    CHR<-fread(paste0('cat ',assoc,' | cut -f 1,2,3,7,9 '),header=T,stringsAsFactors =F)
    column_names<-c('Chr','SNP','bp','b','p')
    setnames(CHR, column_names)
    CHR$log10P = -log10(CHR$p)

    CHR<-CHR[Chr==chr]
    assoc_for_pval[[assoc_counter]]<-CHR


  }else{
    CHR=fread(paste0('cat ',exp_dir,'temp/',chr,'_',trait,'_','temp_gcta_cond.assoc.txt.mlma | cut -f 1,2,3,7,9 '),header=T,stringsAsFactors = F)

    column_names<-c('Chr','SNP','bp','b','p')
    setnames(CHR, column_names)
    CHR<-CHR[Chr==chr]
    CHR$log10P = -log10(CHR$p)
    assoc_for_pval[[assoc_counter]]<-CHR
  }



  if(!cov_counter==0){
   # system(paste0("rm ",exp_dir,"/",chr,'_',trait,'_','temp_gcta_cond.assoc.txt.mlma'),ignore.stderr = T,ignore.stdout = T)
    file.remove(assoc_gcta)
  }


  topsnp<-CHR$SNP[which(CHR$log10P==max(CHR$log10P,na.rm=T))]

  if(length(topsnp)>1){
    dups<-CHR[which(CHR$SNP %in% topsnp),]
    dups<-dups[order(dups$bp),]
    topsnp<-dups$SNP[which.max((abs(dups$b)))]
  }

  topsnp_log10P<-CHR$log10P[which(CHR$SNP== topsnp)]

  #update topSNP_part2 table
  trait=trait

  ifelse(cov_counter==0,covariate <- NA,covariate <- cond_topsnp)
  if(cov_counter>0){

    p_val_covar <- vector(mode="character", length=cov_counter)
    for(m in 1:length(p_val_covar)){
      p_val_covar[m]<-assoc_for_pval[[m]][which(assoc_for_pval[[m]]$SNP==covariate),"log10P"]

    }

    covar_pval<-paste(p_val_covar, sep="",collapse="_")


  }else{
    covar_pval<-NA
  }



  topsnp=topsnp
  topsnp_log10P=topsnp_log10P




  df<-data.frame()
  df <- data.frame(trait,topsnp,topsnp_log10P,covariate,covar_pval,stringsAsFactors = F)

  topSNP_part2  <- rbind (topSNP_part2,df)



  covariate <- topsnp


  pos=sapply(strsplit(topsnp, ":"), `[`, 2, simplify=T)
  #ds_command=paste0("sed -e '/chr",chr,"/!d' -e '/",pos,"/!d' ",dosages_file)
  recodeA_file <- paste0(exp_dir,"genotypes/",chr,"_",pos,".raw")
  recodeA_command <- paste0("/oasis/tscc/scratch/aschitre/software/plink-1.90/plink --bfile genotypes/genotypes --snp ",covariate," --recodeA --out ",exp_dir,"genotypes/",chr,"_",pos)
  system(recodeA_command)
  if(cov_counter == 0){

    #cov=fread(cmd=ds_command,header=F,stringsAsFactors = F)
   cov <- fread(recodeA_file,header=T,stringsAsFactors=F)
    #colsToDelete <-  colnames(cov)[c(1,2)]
   colsToDelete <-  colnames(cov)[c(1,2,3,4,5,6)]
    cov[, (colsToDelete) := NULL]
   # cov<-t(cov)
    cov <- as.data.frame(cov)
    #cov$int<-rep(1,nrow(cov))
    cov<-cbind(fam,cov)
    #cov<-cov[,c("int","V1")]
    cov_file=paste0(exp_dir,"temp/",chr,"_",trait,"_","temp.txt")
    write.table(cov,file=cov_file,row.names=F,col.names=F,quote=F)
  }else{
    c <- fread(recodeA_file,header=T,stringsAsFactors=F)
    colsToDelete <-  colnames(c)[c(1,2,3,4,5,6)]


    c[, (colsToDelete) := NULL]

    c <- as.data.frame(c)
    cov<-cbind(cov,c)
    cov_file=paste0(exp_dir,"temp/",chr,"_",trait,"_","temp.txt")
    write.table(cov,file=cov_file,row.names=F,col.names=F,quote=F)
  }

  if(cov_counter == 0){
    f<-paste0(exp_dir,"temp/code/",trait,"_",chr,".sh")
    if (file.exists(f)) file.remove(f)
    d1<-paste0("#!/bin/bash")
    d2<-paste0("#PBS -N ",trait,"_",chr)  ##This needs to be changed for each pheno file name
    d3<-paste0("#PBS -S /bin/bash")
    d4<-paste0("#PBS -l walltime=6:00:00")
    d5<-paste0("#PBS -l nodes=1:ppn=5")
    d6<-paste0("#PBS -j oe")
    d7<-paste0("#PBS -o ",exp_dir,"temp/pbs_log/",trait,"_",chr,".out")
    d8<-paste0("#PBS -q condo")
    part1 <- paste0("/projects/ps-palmer/software/local/src/gcta_1.93.2beta/gcta64 --mlma --grm ",exp_dir,"grm/genotypes ")
    part2 <- paste0("--mlma-subtract-grm ",exp_dir,"grm/chr",chr,".","genotypes --bfile ",exp_dir,"genotypes/genotypes --chr ",chr," ")
    part3<- paste0("--pheno ",pheno," --qcovar ",cov_file," ")
    part4<- paste0("--out ",exp_dir,"temp/",chr,"_",trait,"_","temp_gcta_cond.assoc.txt ") #chr10.round5_pruned.snpinfo
    part5<- paste0("--thread-num 5")
    command <- paste0(part1,part2,part3,part4,part5)
    line=c(d1,d2,d3,d4,d5,d6,d7,d8,command)
    write(line,file=paste0(exp_dir,"temp/code/",trait,"_",chr,".sh"),append=TRUE)



    system(paste0("qsub ",exp_dir,"temp/code/",trait,"_",chr,".sh"),ignore.stdout = T,ignore.stderr = T,wait = T)
    print("submitting first conditional analysis job")
  }else{

    system(paste0("qsub ",exp_dir,"temp/code/",trait,"_",chr,".sh"),ignore.stdout = T,ignore.stderr = T,wait = T)
  }



#1_u01_huda_akil_lateral_loco_score_temp_gcta_cond.assoc.txt.mlma
#1_u01_huda_akil_lateral_loco_score_temp_gcta_cond.assoc.txt.mlma

  #check if pbs_log file exists

  while (!file.exists(paste0(exp_dir,"temp/",chr,'_',trait,'_','temp_gcta_cond.assoc.txt.mlma'))) {
    Sys.sleep(100)
  }


  CHR=fread(paste0('cat ',exp_dir,'temp/',chr,'_',trait,'_','temp_gcta_cond.assoc.txt.mlma | cut -f 1,2,3,7,9 '),header=T,stringsAsFactors = F)
  column_names<-c('Chr','SNP','bp','b','p')
  setnames(CHR, column_names)
  CHR$log10P = -log10(CHR$p)

  nrow_next<-nrow(CHR[which(CHR$log10P>gwas_threshold),])


  if(nrow_next==0){
    loop_run<-FALSE
    rm(cov,recodeA_command)

  }else{
    cond_topsnp<-topsnp
    cov_counter=cov_counter+1
    assoc_counter=assoc_counter+1
  }


}

save(topSNP_part2,file=paste0(exp_dir,"temp/conditional_analysis/chr",chr,"_",trait,"_QTL_conditional_analysis.RData"))
write.csv(topSNP_part2,paste0(exp_dir,"temp/conditional_analysis/chr",chr,"_",trait,"_QTL_conditional_analysis.csv"),row.names=F,quote=F)
print("Wrote csv file and saved RData file")
topSNP_part2_list[[1]]<- topSNP_part2
#write.table(topSNP_part2,paste0(exp_dir,trait,"_",chr,".txt"),row.names=F,quote=F)


#no_covar_list<-topSNP_part2_list
no_covar<-do.call("rbind",topSNP_part2_list)
#write.csv(final,paste0(exp_dir,"temp/",expname,"_QTL_pruned_with_no.csv"),row.names=F,quote=F)
write.csv(no_covar,paste0(exp_dir,"temp/conditional_analysis/chr",chr,"_",trait,"_QTL_conditional_analysis.csv"),row.names=F,quote=F)


rm(list=ls())
