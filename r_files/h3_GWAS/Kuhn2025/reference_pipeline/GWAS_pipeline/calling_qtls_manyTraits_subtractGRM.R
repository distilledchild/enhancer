args <- commandArgs(trailingOnly = TRUE)
traits<- args[1]
gwas_threshold<-args[2]

gwas_threshold<-as.numeric(gwas_threshold)

exp_dir <- getwd()


exp_dir<-paste0(exp_dir,"/")









project_name=basename(exp_dir)



library(data.table)

setwd(paste0(exp_dir))

system(paste0("mkdir ",exp_dir,"temp/assoc"))
system(paste0("mkdir ",exp_dir,"temp/pheno"))
system(paste0("mkdir ",exp_dir,"temp/code"))
system(paste0("mkdir ",exp_dir,"temp/pbs_log"))
setwd(paste0(exp_dir,"temp/assoc"))



all_traits_chr<-list.files(path=paste0(exp_dir,"results/gwas"),full.names = T)
all_traits_chr<-all_traits_chr[grep(paste0(traits,".mlma"),all_traits_chr)]


filenames<-basename(all_traits_chr)
filenames<-gsub("allChr_","",filenames)
filenames<-gsub(".loco.mlma","",filenames)

column_names<-c('Chr','SNP','bp','b','p')

#V1  V2 V3 V4 V5   V6 V7 V8 V9
#Chr SNP bp A1 A2 Freq  b se  p
readdata<-function(x){

  data<-read.table(x,header=T,stringsAsFactors =F)
  chr<-strsplit(basename(x)[1],split=".",fixed=T)[[1]][1]
  chr<-gsub("[^0-9.-]","",chr)
  data<-data[,c(1,2,3,7,9)]
  colnames(data)<-c('Chr','SNP','bp','b','p')
  data=data[which(!is.nan(data$p)),]
  data<-data[which(data$Chr==chr),]
  data$log10P = -log10(data$p)
  return(data)

}



all_data<-lapply(all_traits_chr,readdata)


filenames<-gsub(".mlma","",filenames)

names(all_data)<-filenames

all_traits_chr<-filenames


topSNP_part1<-data.frame()

for(i in seq_along(all_traits_chr)){


  if(length(which(all_data[[all_traits_chr[i]]]$log10P > gwas_threshold  )) > 0){

    trait=all_traits_chr[i]
    above_threshold="yes"

    df <- data.frame(trait,above_threshold,stringsAsFactors = F)

  }else{
    trait=all_traits_chr[i]
    above_threshold="no"
    df <- data.frame(trait,above_threshold,stringsAsFactors = F)
  }

  topSNP_part1  <- rbind (topSNP_part1,df)

}



signi_traits_n<-length(which(topSNP_part1$above_threshold=="yes"))
if(signi_traits_n>0){



signi_traits<-topSNP_part1$trait[which(topSNP_part1$above_threshold=="yes")]

#this will contain chrs with signi results
signi_chrs <- setNames(replicate(signi_traits_n,data.frame()),signi_traits)

for(i in seq_along(signi_traits)){

  chrs<-unlist(unique(all_data[[signi_traits[i]]][which(all_data[[signi_traits[i]]]$log10P > gwas_threshold  ) ,"Chr" ]))
  chrs<-unname(chrs)

  signi_chrs[[signi_traits[i]]]<-chrs


}


summary_QTL <- list()
topSNP_part2_list <- list ()


system(paste0("mkdir ",exp_dir,"temp/r2"))
setwd(paste0(exp_dir,"temp/r2"))


#define parameters that can be soft coded
#This is log10P - 1.5 (to see if there are any supporting SNPs)
#sup_P_term<-1.5
sup_P_term<-2
#This is 0.5 MB flanking dist- to see if there are any supporting snps for top snp or if it is a rogue snp
sup_dist<-500000
#LD command parameters for plink
ld_window_kb=11000
ld_r2=0.4
#QTL boundary distance
qtl_dist<-1000000



for(i in seq_along(signi_traits)){
  topSNP_part2_list <- list ()
  chrs <-signi_chrs[[signi_traits[i]]]
  for(j in seq_along(chrs)){

    #master_list[[signi_traits[1]]][1]
    CHR <- all_data[[signi_traits[i]]][which(all_data[[signi_traits[i]]]$Chr == chrs[j]),]


    topSNP_part2<-data.frame()

    loop_run <- TRUE
    while(loop_run){
      #run until there are no topsnps
      #find top snp

      all_topsnps<-CHR[which((CHR$log10P > gwas_threshold)),]
      #topsnp<-CHR$SNP[which((CHR$log10P > 5.34) & (CHR$log10P==max(CHR$log10P)))]
      topsnp<-all_topsnps$SNP[which(all_topsnps$log10P==max(all_topsnps$log10P))]


      if(length(topsnp)>1){
        dups<-CHR[which(CHR$SNP %in% topsnp),]
        dups<-dups[order(dups$bp),]
        topsnp<-dups$SNP[which.max((abs(dups$b)))]
      }

      topsnp_log10P<-CHR$log10P[which(CHR$SNP== topsnp)]
      sup_P<-(topsnp_log10P - sup_P_term)


      chr_topsnp<-gsub("chr","",strsplit(topsnp,split=":")[[1]][1])
      ps_topsnp<-strsplit(topsnp,split=":")[[1]][2]
      ps_topsnp<-as.numeric(ps_topsnp)

      start<- ps_topsnp -sup_dist
      stop<- ps_topsnp +sup_dist


      ##This evaluates if it is not a rogue SNP
      if(nrow(CHR[which((CHR$bp %in% seq(start,stop)) & (CHR$log10P > sup_P)),])>2){
        trait=signi_traits[i]
        topsnp=topsnp
        QTL="yes"
        topsnp_log10P=topsnp_log10P
        df<-data.frame()
        df <- data.frame(trait,topsnp,QTL,topsnp_log10P,stringsAsFactors = F)

      }else{
        trait=signi_traits[i]
        topsnp=topsnp
        QTL="no"
        topsnp_log10P=topsnp_log10P
        df<-data.frame()
        df <- data.frame(trait,topsnp,QTL,topsnp_log10P,stringsAsFactors = F)
      }

      topSNP_part2  <- rbind (topSNP_part2,df)


      system(paste0("/oasis/tscc/scratch/aschitre/software/plink-1.90/plink --bfile ",exp_dir,"/genotypes/genotypes --chr ",chr_topsnp," --nonfounders --r2  --ld-snp ",topsnp," --ld-window 1000000 -ld-window-kb ",ld_window_kb," --ld-window-r2 ",ld_r2," --out ",exp_dir,"temp/r2/","temp_qtl_n_",traits),ignore.stdout = T,ignore.stderr = T,wait = T)



      f<-paste0(exp_dir,"temp/r2/temp_qtl_n_",traits,".ld")
      if (file.exists(f)){


        qtl_snps<-fread(paste0("cat ",f," | grep -v R2 | awk '{ print $6 }'"),sep=" ",header = F,stringsAsFactors = F)

        #0.5
        #1MB
        sub_start<-ps_topsnp-qtl_dist
        sub_stop<-ps_topsnp+qtl_dist

        nrow_next<-nrow(CHR[which((!CHR$SNP %in% c(qtl_snps$V1,topsnp)) & (CHR$log10P>5.34) & (!CHR$bp %in% seq(sub_start,sub_stop))),])

        if(!(nrow_next ==0)){
          CHR<-CHR[which((!CHR$SNP %in% c(qtl_snps$V1,topsnp)) & (!CHR$bp %in% seq(sub_start,sub_stop))),]
        }else{

          loop_run<-FALSE
        }

      }
      else{
        nrow_next<-nrow(CHR[which((!CHR$SNP %in% topsnp) & (CHR$log10P>5.34) ),])

        if(!(nrow_next ==0)){
          CHR<-CHR[which((!CHR$SNP %in% topsnp)),]
        }else{
          loop_run<-FALSE
        }


      }


    }

    topSNP_part2_list[[j]]<-topSNP_part2

  }

  summary_QTL[[i]] <- topSNP_part2_list

}





names(summary_QTL) <- signi_traits

temp <- unlist(summary_QTL, recursive = FALSE)
temp<-do.call("rbind",temp)


write.csv(temp,paste0(exp_dir,"results/qtls/QTL_",traits,".csv"),row.names=F,quote=F)

}else{
  file.create(paste0(exp_dir,"results/qtls/QTL_",traits,".csv"))
}