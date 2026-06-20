exp_dir <- getwd()



exp_dir<-paste0(exp_dir,"/")



project_name=basename(exp_dir)



raw=read.csv(paste0(exp_dir,"results/qtls/table3_final.csv"),header=T,stringsAsFactors = F)



qtls=raw



library(data.table)



#the file I gave has the top variant for each gene
#total rows/ total genes= 15,231
#with a qvalue < 0.05 filter there should be 200 eQTLs

raw=fread(file="/projects/ps-palmer/apurva/eqtl_all_tissues/ratgtex_top_assoc.txt",header=T,stringsAsFactors=F)
colnames(raw)[grep("tss_distance",colnames(raw))]<-"tss"

all_tissues_top_assoc=raw[qval < 0.05]
#raw$variant_id <- gsub("chr","",raw$variant_id)

#total rows= 82649
#after qvalue threshold = 18678










results_eqtl<-data.table()
for(i in 1:nrow(qtls)){

  trait<-qtls$trait[i]
  trait_topsnp<-qtls$topsnp[i]
  chr<-qtls$chr[i]
  pos<-qtls$pos[i]



  eqtl=all_tissues_top_assoc
  eqtl$chrom=paste0("chr",eqtl$chrom)
  temp<-eqtl[chrom==qtls$chr[i] & pos %between% c(qtls$LD_interval_start[i],qtls$LD_interval_stop[i])]


  #calculate LD r2 with variant_id

  if(nrow(temp)>0){


    system(paste0("/oasis/tscc/scratch/aschitre/software/plink-1.90/plink --bfile ",exp_dir,"genotypes/genotypes --chr ",gsub("chr","",chr)," --nonfounders --r2  --ld-snp ",trait_topsnp," --ld-window 100000 -ld-window-kb 6000 --ld-window-r2 0.6 --out ",exp_dir,"temp/r2/eqtl_temp"),ignore.stdout = T,ignore.stderr = T,wait = T)


    column_names<-c("trait_topsnp","variant_id","rsquare")
    readLDdata<-function(x,trait,chr,pos){

      data<-fread(x,header=T,stringsAsFactors =F,select=c(3,6,7))
      setnames(data, column_names)
      return(data)
    }


    LD=readLDdata(paste0(exp_dir,"temp/r2/eqtl_temp.ld"),trait = trait,chr=chr,pos=pos)

    if(nrow(LD)>0){

      merged=merge(LD,temp,by="variant_id",all=F)



      if(nrow(merged)>0){



        colnames(merged)[1]<-"eqtl_topsnp"


        merged$trait<-trait

       merged$comments<-NA

       colnames(merged)[grep("pval_nominal",colnames(merged))][1]<-"eqtl_topsnp_pval"
       results_eqtl<- rbind(merged,results_eqtl)
      }
    }}else{


      no_coloc=setNames(data.table(matrix(nrow = 1, ncol = 22)), c("eqtl_topsnp","trait_topsnp","rsquare","tissue","gene_id","gene_name","tss","num_var","chrom","pos","ref","alt","af","eqtl_topsnp_pval","slope","slope_se","pval_beta","qval","pval_nominal_threshold","log2_aFC","trait","comments"))



      no_coloc$trait_topsnp[1] <- trait_topsnp
      no_coloc$trait[1] <- trait
      no_coloc$comments[1] <- "No EQTL info found for this locus"


      results_eqtl<- rbind(no_coloc,results_eqtl)
    }


}


write.csv(results_eqtl,paste0(exp_dir,"results/eqtl/eQTL_results.csv"),row.names = F,quote=T)

#annotate eqtl results



.libPaths(c("/home/aschitre/R_libs", .libPaths()))
library(data.table)
library(plyr)

library(annotables, lib.loc="/home/aschitre/R_libs/")




all_regions_df<-results_eqtl



colnames(all_regions_df)[5]<-"ensembl_gene"

ensembl_genes<-unique(all_regions_df$ensembl_gene)
ensembl_genes<-ensembl_genes[!is.na(ensembl_genes)]



all_regions_df$RGD_link<-NA
all_regions_df$gene_long_name<-NA
all_regions_df$gene_name<-NA

if(nrow(rnor6[which(rnor6$ensgene %in% ensembl_genes),])>0){

  for(i in 1:length(ensembl_genes))  {
    if(ensembl_genes[i] %in% rnor6$ensgene){
      all_regions_df$gene_name[which(all_regions_df$ensembl_gene %in% ensembl_genes[i])]<-as.character(rnor6[which(rnor6$ensgene %in% ensembl_genes[i]),"symbol"] )
      #rnor6[which(rnor6$ensgene %in% ensembl_genes[i]),"description"]
      acc<-gsub("]","",strsplit(as.character(rnor6[which(rnor6$ensgene %in% ensembl_genes[i]),"description"]),split="Acc:")[[1]][2])
      if(is.na(acc)){
        all_regions_df$RGD_link[which(all_regions_df$ensembl_gene %in% ensembl_genes[i])]<-NA
        all_regions_df$gene_long_name[which(all_regions_df$ensembl_gene %in% ensembl_genes[i])]<-NA

      }else{
        all_regions_df$RGD_link[which(all_regions_df$ensembl_gene %in% ensembl_genes[i])]<-paste0("https://rgd.mcw.edu/rgdweb/report/gene/main.html?id=",acc)
        all_regions_df$gene_long_name[which(all_regions_df$ensembl_gene %in% ensembl_genes[i])]<-strsplit(as.character(rnor6[which(rnor6$ensgene %in% ensembl_genes[i]),"description"]),split=" \\[")[[1]][1]

      }


    }




  }




}



library(openxlsx)
#write.xlsx(all_regions_df,paste0(base_dir,"/",expname,"/eqtl/brain_eQTL_results_annotated.csv"))

raw=all_regions_df



#Pubmed links




gene_name<-unique(raw$gene_name)
gene_name<-gene_name[!is.na(gene_name)]
gene_name<-toupper(gene_name)


raw$human_entrez_id<-NA
raw$pubmed_url<-NA



if(nrow(grch38[which(grch38$symbol %in% gene_name),])>0){

  for(i in 1:length(gene_name))  {
    if(gene_name[i] %in% grch38$symbol){
      raw$human_entrez_id[which(toupper(raw$gene_name) %in% gene_name[i])]<-as.character(grch38[which(grch38$symbol %in% gene_name[i]),"entrez"])
      #grch38[which(grch38$symbol %in% gene_name[i]),"description"]
      raw[grep("(",raw$human_entrez_id,fixed=T),"human_entrez_id"]<-as.numeric(gsub(".*?([0-9]+).*", "\\1", grep("(",raw$human_entrez_id,fixed=T,value=T)))

      raw$pubmed_url[which(toupper(raw$gene_name) %in% gene_name[i])]<-paste0("https://www.ncbi.nlm.nih.gov/kis/ortholog/",raw$human_entrez_id[which(toupper(raw$gene_name) %in% gene_name[i])],"/?scope=7776#literature-tab")

    }




  }




}



#GWAS catalog links
library(data.table)
gene_catalog=fread(file="/projects/ps-palmer/apurva/gwas_catalog/GWAS_catalog_genes_lookup_23Nov2020.csv",header=T,stringsAsFactors = F)

#https://www.ebi.ac.uk/gwas/genes/CSMD1


raw$gwas_catalog_url<-NA
for(i in 1:nrow(raw)){

  if(toupper(raw$gene_name[i]) %in% toupper(gene_catalog$gene_names)) {

    raw$gwas_catalog_url[i] <- paste0("https://www.ebi.ac.uk/gwas/genes/",raw$gene_name[i])

  }



}

write.csv(raw,paste0(exp_dir,"results/eqtl/eqtl_results_annotated.csv"),quote=T,row.names=F)
