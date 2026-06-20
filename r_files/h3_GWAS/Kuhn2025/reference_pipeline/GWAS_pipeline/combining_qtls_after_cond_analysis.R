exp_dir <- getwd()
project_name <- basename(exp_dir)


exp_dir<-paste0(exp_dir,"/")




read_data<-function(x){
  data=read.csv(x,header=T,stringsAsFactors = F)
  return(data)
}


final=read.csv("results/qtls/all_including_conditional_analysis.csv",header=T,stringsAsFactors=F)





final$trait_snp <- paste0(final$trait,"_",final$topsnp)

final=unique(final)



topsnps<-final

topsnps$chr<-sapply(strsplit(topsnps$topsnp,split=":"),`[`,1)
topsnps$pos<-sapply(strsplit(topsnps$topsnp,split=":"),`[`,2)








files=list.files(paste0(exp_dir,"/temp/conditional_analysis"),pattern=".csv",full.names = T)

all=lapply(files,function(x){data=read.csv(x,header=T,stringsAsFactors = F)
return(data)
})

cond_topsnps=do.call("rbind",all)




QTL_same_chr_counts<-with(topsnps,table(trait,chr))
QTL_same_chr_counts<-as.data.frame(QTL_same_chr_counts,stringsAsFactors=F)


traits_chr<-QTL_same_chr_counts[which(QTL_same_chr_counts$Freq > 1),]



#remove QTL_same_chr_counts


for(i in 1:nrow(traits_chr)){


  topsnps<- topsnps[which(!((topsnps$trait %in% traits_chr$trait[i]) & (topsnps$chr %in% traits_chr$chr[i]))),]



}




cond_topsnps<-cond_topsnps[,c("trait","topsnp","topsnp_log10P")]


topsnps<-topsnps[,c("trait","topsnp","topsnp_log10P")]

final_QTL<-rbind(cond_topsnps,topsnps)


final_QTL$chr<-sapply(strsplit(final_QTL$topsnp,split=":"),`[`,1)
final_QTL$pos<-sapply(strsplit(final_QTL$topsnp,split=":"),`[`,2)

final_QTL=unique(final_QTL)

write.csv(final_QTL,paste0(exp_dir,"results/qtls/final_QTL_list.csv"),row.names = F,quote=F)
