

files=list.files(paste0("results/phewas"),full.names = T)
library(data.table)

not_empty<-list()

for(i in 1:length(files)){


  lines=fread(paste0("wc -l ",files[i]))
  if(lines$V1[1] > 1){

  not_empty[[i]]<-files[i]
  }
}


files=do.call("rbind",not_empty)



file_n=grep("table1_unfiltered",files,value=T)


if(length(file_n)>0){
  table1_unfiltered=list()
  for(i in 1:length(file_n)){
    table1_unfiltered[[i]] =read.csv(file_n[i],header=T,stringsAsFactors = F)
  }
  table1_unfiltered=do.call("rbind",table1_unfiltered)
  write.csv(table1_unfiltered,paste0("results/phewas/table1_unfiltered_phewas.csv"),row.names = F,quote = T)
}

file_n=grep("table1_filtered",files,value=T)


if(length(file_n)>0){
  table1_filtered=list()
  for(i in 1:length(file_n)){
    table1_filtered[[i]] =read.csv(file_n[i],header=T,stringsAsFactors = F)
  }
  table1_filtered=do.call("rbind",table1_filtered)
  write.csv(table1_filtered,paste0("results/phewas/table1_filtered_phewas.csv"),row.names = F,quote = T)

}


file_n=grep("table2_3MB",files,value=T)


if(length(file_n)>0){
  table2_3MB=list()
  for(i in 1:length(file_n)){
    table2_3MB[[i]] =read.csv(file_n[i],header=T,stringsAsFactors = F)
  }
  table2_3MB=do.call("rbind",table2_3MB)
  write.csv(table2_3MB,paste0("results/phewas/table2_3MB_window_phewas.csv"),row.names = F,quote = T)

}




rm(list=ls())
