exp_dir <- getwd()
expname <- basename(exp_dir)


exp_dir<-paste0(exp_dir,"/")



project_name=basename(exp_dir)




files=list.files(paste0(exp_dir,"results/qtls"),full.names = T)
files=files[sapply(files, file.size) > 0]
read_files<-function(x){
  data=read.csv(x,header = T,stringsAsFactors = F)
  return(data)
}


all=lapply(files,read_files)


all=do.call("rbind",all)



topsnps<-all

topsnps$chr<-sapply(strsplit(topsnps$topsnp,split=":"),`[`,1)
topsnps$pos<-sapply(strsplit(topsnps$topsnp,split=":"),`[`,2)

colnames(topsnps)[1]<-"chr.trait"

topsnps$trait<-sapply(strsplit(topsnps$chr.trait,split=".",fixed=T),`[`,2)




###Traits with more than one QTL on same chromosome
QTL_same_chr_counts<-with(topsnps,table(trait,chr))
QTL_same_chr_counts<-as.data.frame(QTL_same_chr_counts,stringsAsFactors=F)


traits_chr<-QTL_same_chr_counts[which(QTL_same_chr_counts$Freq > 1),]


if(nrow(traits_chr)>0){

  cond_analysis<-TRUE
  cat("Do conditional analysis")
  write.table(traits_chr,file="results/qtls/conditional_analysis.csv",sep=",",col.names = F,row.names = F,quote=F)
  write.csv(topsnps,"results/qtls/all_including_conditional_analysis.csv")
}else{
  cond_analysis<-FALSE
  cat("No QTLs on same chromosomes for same traits")
  write.csv(topsnps,paste0(exp_dir,"results/qtls/final_QTL_list.csv"),row.names=F,quote=F)
}
