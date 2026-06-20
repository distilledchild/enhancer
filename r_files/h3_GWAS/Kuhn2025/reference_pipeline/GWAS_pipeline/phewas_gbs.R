exp_dir <- getwd()
project_name <- basename(exp_dir)


exp_dir<-paste0(exp_dir,"/")




library(data.table)

#load(paste0(exp_dir,"residuals/",expname,"_traits.RData"))



#filenames<-list.files("/oasis/tscc/scratch/aschitre/round8/unpruned",pattern="GWAS_combined_",full.names = T)

#This part reads in the most recent assoc file DB. Important to extract the date from the filenames because Apurva might refresh the decay clock on the assoc file on scratch

#dates<-gsub("[[:alpha:]]", "", basename(filenames))
#dates<-gsub("__","",dates)
#dates<-gsub(".","",dates,fixed = T)
#"%m_%d_%Y
#dates<-as.Date(dates,"%m_%d_%Y")
#dates<-format(dates, format="%m_%d_%Y")
#keep<-dates[order(dates,decreasing=T)][1]

#read in assoc and assign column names.
#also gsub from trait column to remove filename metadata
#all_assoc=fread(file=grep(paste0("*",keep,"*"),filenames,value=T),header=F,stringsAsFactors = F)

#/oasis/tscc/scratch/aschitre/round8/unpruned/GWAS_06July2020.txt
all_assoc=fread(file="/projects/ps-palmer/apurva/round8/unpruned/GWAS_29Mar2021.txt",header=F,stringsAsFactors = F)

setnames(all_assoc,c("snp","p_score","trait"))

all_assoc$trait<-gsub(".*allChr_","",all_assoc$trait)
all_assoc$trait<-gsub(".assoc.txt","",all_assoc$trait,fixed=T)

#read in query qtl table
#expname #exp_dir
query_qtl<-read.csv(paste0(exp_dir,"results/qtls/table3_final.csv"),header=T,stringsAsFactors = F)
#colnames(query_qtl)[2]<-"trait"
query_qtl<-query_qtl<-query_qtl[,c("trait","topsnp")]

query_qtl$chr<-sapply(strsplit(query_qtl$topsnp,split=":"),`[`,1)
query_qtl$pos<-sapply(strsplit(query_qtl$topsnp,split=":"),`[`,2)

##########
###for each topsnp in query qtl table, get the p_val for all other traits at that topsnp AND exclude expname from the query



#first remove all_assoc containing expname in trait column
all_assoc<-all_assoc[all_assoc$p_score<=0.0001]


all_assoc$trait = basename(all_assoc$trait)
all_assoc$trait=gsub(".loco.mlma","",all_assoc$trait)

#all_assoc<-all_assoc[!grepl(expname,all_assoc$trait)]

#table1_versionA







#https://www.dropbox.com/s/caca0nvmzixhf2v/P50.genotypes.ID.order.csv?dl=0
#https://www.dropbox.com/s/caca0nvmzixhf2v/P50.genotypes.ID.order.csv?dl=0

#https://www.dropbox.com/s/i960s0cgfhin7p9/pheWAS.R?dl=0


#table2




all_assoc$chr<-sapply(strsplit(all_assoc$snp,split=":"),`[`,1)
all_assoc$pos<-sapply(strsplit(all_assoc$snp,split=":"),`[`,2)

assoc_part2<-all_assoc

split.data.table <- function(x, f, drop = FALSE, by, flatten = FALSE, ...){
  if(missing(by) && !missing(f)) by = f
  stopifnot(!missing(by), is.character(by), is.logical(drop), is.logical(flatten), !".ll" %in% names(x), by %in% names(x))
  if(!flatten){
    .by = by[1L]
    tmp = x[, list(.ll=list(.SD)), by = .by, .SDcols = if(drop) setdiff(names(x), .by) else names(x)]
    setattr(ll <- tmp$.ll, "names", tmp[[.by]])
    if(length(by) > 1L) return(lapply(ll, split.data.table, drop = drop, by = by[-1L])) else return(ll)
  } else {
    tmp = x[, list(.ll=list(.SD)), by=by, .SDcols = if(drop) setdiff(names(x), by) else names(x)]
    setattr(ll <- tmp$.ll, 'names', tmp[, .(nm = paste(.SD, collapse = ".")), by = by, .SDcols = by]$nm)
    return(ll)
  }
}



#0.0001

#assoc_by_region<-split.data.table(all_assoc, by="trait")




table2<-data.frame(stringsAsFactors = f)
table2_list<-list()
k=0
for(j in 1:nrow(query_qtl)){
  query_trait<-query_qtl$trait[j]

  pos_query<-as.numeric(query_qtl$pos[j])
  chr_query<-query_qtl$chr[j]
  query_topsnp<-paste0(chr_query,":",pos_query)

  #chr_query=paste0("chr",chr_query)

  #assoc_by_region[[i]][chr==chr_query & pos %in% seq(start,stop)]

  start<-pos_query + 3000000
  stop<-pos_query - 3000000
  #first remove query trait from consideration
  #all_assoc<-all_assoc[!grepl(expname,all_assoc$trait)]
  query_trait_assoc<-all_assoc[!grepl(query_trait,all_assoc$trait)]
  query_trait_assoc=all_assoc[chr==chr_query & pos %in% seq(start,stop)]
  assoc_by_region<-split.data.table(query_trait_assoc, by="trait")
  #query_topsnp<-paste0("chr",query_topsnp)






  if(length(assoc_by_region)>0){


  for(i in 1:length(assoc_by_region)){

    #chr_query<-paste0("chr",chr_query)
    temp<-assoc_by_region[[i]][chr==chr_query & pos %in% seq(start,stop)]

    if(nrow(temp)>0){
      k=k+1
      table2_list[[k]]<-temp
      trait2_topsnp<-temp[which.min(temp$p_score)]$snp
      trait2_p_score<-temp[which.min(temp$p_score)]$p_score
      trait2_name<-temp[which.min(temp$p_score)]$trait
     #
      chr<-gsub("chr","",chr_query)

      if(query_topsnp!=trait2_topsnp){


        command<-paste0("/oasis/tscc/scratch/aschitre/software/plink-1.90/plink --bfile ",exp_dir,"/genotypes/genotypes --chr ",gsub("chr","",chr)," --nonfounders --r2  --ld ",paste0("","",query_topsnp)," ",paste0("","",trait2_topsnp)," --out ",exp_dir,"temp/r2/phewas_temp")
        #paste0("/oasis/tscc/scratch/aschitre/software/plink-1.90/plink --bfile tscc_dir/genotypes/u01_suzanne_mitchell_genotypes --chr ",gsub("chr","",chr)," --nonfounders --r2  --ld-snp ",topsnp," --ld-window 100000 -ld-window-kb 6000 --ld-window-r2 0 --out tscc_dir/r2/",project_name,"_LZ_","temp_qtl_n")
        r2_log<-system(command,wait=T)
        command<-paste0("if grep -q Error ",exp_dir,"temp/r2/phewas_temp.log; then echo 'No SNP'; else cat ",exp_dir,"temp/r2/phewas_temp.log | grep R-sq; fi")

        r2_log<-fread(command,header=F,stringsAsFactors = F)

        if(ncol(r2_log)>2){


          r2<-r2_log$V3
          dprime<-r2_log$V6

        }
      }else{
        r2<-1
        dprime<-NA
      }

      if(exists("r2")){
        df<-data.frame(query_trait=query_trait,query_topsnp=query_topsnp,trait2_name=trait2_name,trait2_topsnp=trait2_topsnp,trait2_p_score=trait2_p_score,r2=r2,dprime=dprime)
        table2<-rbind(df,table2)
        rm(r2)
      }

    }


  }

  }
}




write.csv(table2,paste0(exp_dir,"results/phewas/table2_3MB_window_phewas_gbs.csv"),row.names=F,quote=F)


all_assoc<-assoc_part2
#all_assoc=fread("/oasis/tscc/scratch/aschitre/round8/pruned/phewas/input/all_assoc_14May2019.txt",sep="\t",header=T,stringsAsFactors = F)
#all_assoc$chr<-sapply(strsplit(all_assoc$snp,split=":"),`[`,1)
#all_assoc$pos<-sapply(strsplit(all_assoc$snp,split=":"),`[`,2)



table1_unfiltered<-data.frame(stringsAsFactors = F)
table1_filtered<-data.frame(stringsAsFactors = F)

for(i in 1:nrow(query_qtl)){

  #remove query_trait
  query_trait<-query_qtl$trait[i]


  table1_versionA<-all_assoc[snp==query_qtl$topsnp[i]]
  table1_versionA<-table1_versionA[!grepl(query_trait,table1_versionA$trait)]
  if(nrow(table1_versionA)>0){

    query_trait<-query_qtl$trait[i]
    topSNP<-query_qtl$topsnp[i]
    df <-cbind(query_trait=rep(query_trait,times=nrow(table1_versionA)),query_topSNP=rep(topSNP,times=nrow(table1_versionA)),table1_versionA)


    table1_unfiltered <-rbind(df,table1_unfiltered)

    table1_versionB<-table1_versionA[-log10(p_score)>5.6]
    if(nrow(table1_versionB)>0){

      df <-cbind(query_trait=rep(query_trait,times=nrow(table1_versionB)),query_topSNP=rep(topSNP,nrow(table1_versionB)),table1_versionB)
      table1_filtered <-rbind(df,table1_filtered)
    }

  }



}





write.csv(table1_unfiltered,paste0(exp_dir,"results/phewas/table1_unfiltered_phewas_gbs.csv"),row.names=F,quote=F)
write.csv(table1_filtered,paste0(exp_dir,"results/phewas/table1_filtered_phewas_gbs.csv"),row.names=F,quote=F)


q(save="no")
