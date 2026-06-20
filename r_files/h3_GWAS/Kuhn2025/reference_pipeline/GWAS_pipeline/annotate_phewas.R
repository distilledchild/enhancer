load("/projects/ps-palmer/apurva/trait_descriptions/trait_descriptions_22Oct2020.RData")

all_desc$dataset="P50 2014"


library(data.table)
exp_dir <- getwd()
project_name <- basename(exp_dir)


exp_dir<-paste0(exp_dir,"/")







files=list.files(paste0("results/phewas"),full.names = T)
library(data.table)

not_empty<-list()



wc_l_table1=fread(paste0("wc -l results/phewas/table1_unfiltered_phewas.csv"))

if(wc_l_table1$V1[1] > 1){
  table1_unfiltered<-read.csv(paste0("results/phewas/table1_unfiltered_phewas.csv"),header=T,stringsAsFactors = F)
}else{
  table1_unfiltered=data.frame(dummy=character(0))
}

wc_l_table1_filt=fread(paste0("wc -l results/phewas/table1_filtered_phewas.csv"))

if(wc_l_table1_filt$V1[1] > 1){
  table1_filtered<-read.csv(paste0("results/phewas/table1_filtered_phewas.csv"),header=T,stringsAsFactors = F)
  table1_filtered<-unique(table1_filtered)
  table1_filtered$query_trait<-basename(table1_filtered$query_trait)
}else{
  table1_filtered=data.frame(dummy=character(0))
}








wc_l_table2=fread(paste0("wc -l results/phewas/table2_3MB_window_phewas.csv"))
if(wc_l_table2$V1[1] > 1){
  table2<-read.csv(paste0("results/phewas/table2_3MB_window_phewas.csv"),header=T,stringsAsFactors = F)
  table2$query_trait<-basename(table2$query_trait)
}else{
  table2=data.frame(dummy=character(0))
}









if(nrow(table1_unfiltered)>0){

  table1_unfiltered$trait_ref_name<-paste0(gsub("_"," ",table1_unfiltered$query_trait)," ",table1_unfiltered$query_topSNP)
  table1_unfiltered$trait_topsnp<-paste0(table1_unfiltered$query_trait,"_",table1_unfiltered$query_topSNP)
  table1_unfiltered$log10P<- -log10(table1_unfiltered$p_score)
  table1_unfiltered$log10P<-format(round(table1_unfiltered$log10P,3),nsmall=3)
  table1_unfiltered$p_score<-format(table1_unfiltered$p_score,scientific = T)
  temp_trait<-paste0(all_desc$experiment_name,"_",all_desc$trait)
  all_desc$trait<-temp_trait
  table1_unfiltered<-merge(table1_unfiltered,all_desc[,c("trait","trait_description","dataset")],by="trait",all.x=T)
  #table1_unfiltered<-merge(table1_unfiltered,desc[,c("trait","trait_description","dataset")],by="trait",all.x=T)
  table1_unfiltered[is.na(table1_unfiltered$dataset),"dataset"]<-"p50_david_dietz_2020"

  table1_filtered<-table1_unfiltered[which(table1_unfiltered$log10P > 5.6),]
  table1_filtered=unique(table1_filtered)
  table1_unfiltered=unique(table1_unfiltered)
}

#temp_trait<-paste0(all_desc$experiment_name,"_",all_desc$trait)
#all_desc$trait<-temp_trait

all_desc=all_desc[,c("trait","trait_description","dataset")]

if(nrow(table2)>0){

  table2$trait_ref_name<-paste0(gsub("_"," ",table2$query_trait)," ",table2$query_topsnp)
  table2$trait_topsnp<-paste0(table2$query_trait,"_",table2$query_topsnp)
  table2$log10P<- -log10(table2$trait2_p_score)
  #format(round(x, 2), nsmall = 2)
  table2$r2<-format(round(table2$r2,3),nsmall=3)
  table2$dprime<-format(round(table2$dprime,3),nsmall=3)
  table2$log10P<-format(round(table2$log10P,3),nsmall=3)
  table2$trait<-table2$trait2_name

  #all_desc$trait<-paste0(all_desc$experiment_name,"_",all_desc$trait)

  table2<-merge(table2,all_desc[,c("trait","trait_description","dataset")],by="trait",all.x=T)
  table2[is.na(table2$dataset),"dataset"]<-"p50_david_dietz_2020"
  table2=unique(table2)


}

save.image("results/phewas/annotated_phewas_tables.RData")


rm(list=ls())
