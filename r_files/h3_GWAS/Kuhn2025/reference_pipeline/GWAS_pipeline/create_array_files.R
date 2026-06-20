exp_dir <- getwd()


exp_dir<-paste0(exp_dir,"/")

load("data/residuals/traits.RData")


write.table(unique(traits),"data/phenotypes.txt",row.names = F,col.names = F,quote=F)
