setwd("/projects/ps-palmer/apurva/p50_david_dietz_2020")

load("data/residuals/residuals.RData")
load("data/residuals/traits.RData")



total_gwas_jobs_files <- (20 * length(traits))/200
total_gwas_jobs_files<-ceiling(total_gwas_jobs_files)


chroms=seq(from=1,to=20)
chroms=paste0("chr",chroms)
nrows<-20*length(traits)



pheno_chrom<-list()
chrom_reps<-list()


for(i in 1:length(traits)){


  pheno_chrom[[i]]<-rep(traits[i],length(chroms))
  chrom_reps[[i]]<-chroms
  #write.table(data.frame(rep(traits[i],length(chroms)),chroms),paste0("data/gwas_jobs/job.",i,))

}


total_gwas_jobs_df<-data.frame(chroms=unlist(chrom_reps),pheno=unlist(pheno_chrom))


max <- 20
x <- seq_along(traits)
d1 <- split(traits, ceiling(x/max))



for(k in 1:length(d1)){


  write.table(total_gwas_jobs_df[which(total_gwas_jobs_df$pheno %in% d1[[k]]),],paste0("data/gwas_jobs/gwas.jobs.",k,".txt"),quote=F,row.names = F,col.names = F)
}







#for(k in 2:4){
#  array_index<-nrow(total_gwas_jobs_df[which(total_gwas_jobs_df$pheno %in% d1[[k]]),])
#  system(paste0("qsub -q hotel -t 1-",array_index," -l nodes=1:ppn=4 -l walltime=3:00:00 -N gwas_",k," -o gwas_",k,".log -v pheno_file=",k," /projects/ps-palmer/apurva/genetic_analysis/code/gwas_array_jobs.sh"))
#}

#for(k in 6:8){
#  array_index<-nrow(total_gwas_jobs_df[which(total_gwas_jobs_df$pheno %in% d1[[k]]),])
#  system(paste0("qsub -q condo -t 1-",array_index," -l nodes=1:ppn=4 -l walltime=3:00:00 -N gwas_",k," -o gwas_",k,".log -v pheno_file=",k," /projects/ps-palmer/apurva/genetic_analysis/code/gwas_array_jobs.sh"))
#}




#job monitoring
#qstat -t | grep aschitre
