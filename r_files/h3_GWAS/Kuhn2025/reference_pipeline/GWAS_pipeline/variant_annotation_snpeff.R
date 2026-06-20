library(data.table)
exp_dir <- getwd()
project_name <- basename(exp_dir)


exp_dir<-paste0(exp_dir,"/")



project_name=basename(exp_dir)

.libPaths(c("/home/aschitre/R_libs", .libPaths()))

library("data.table", lib.loc="/home/aschitre/R_libs/")

library("BSgenome.Rnorvegicus.UCSC.rn6",lib.loc="/home/aschitre/R_libs/")


topsnps=read.csv(paste0(exp_dir,"results/qtls/table3_final.csv"),header=T,stringsAsFactors=F)

topsnps$chr<-sapply(strsplit(topsnps$topsnp,split=":"),`[`,1)
topsnps$pos<-sapply(strsplit(topsnps$topsnp,split=":"),`[`,2)



qtls<-topsnps







#create directory
#dir.create(paste0("/oasis/tscc/scratch/aschitre/round8/unpruned/snp_annotations/",project_name), showWarnings = FALSE)
dir.create(paste0("/oasis/tscc/scratch/aschitre/round8/unpruned/snp_annotations/snpeff_input/",project_name), showWarnings = FALSE)

#LD file

LD_ref_snps<-function(x,trait,chr,pos){
  data<-fread(x,header=T,stringsAsFactors =F,select=c(6,7))
  column_names_ld<-c("snp","rsquare")
  setnames(data, column_names_ld)
  data<-data[which(data$rsquare>=0.6),]
  #get reference allele
  data$reference_allele<-NA
  for(j in 1:nrow(data)){
    #strsplit(data$snp[j],split=":")[[1]][1]  ##chr
    #as.numeric(strsplit(data$snp[j],split=":")[[1]][2]) ##pos
    chr_ref<-strsplit(data$snp[j],split=":")[[1]][1]
    #chr_ref=paste0("chr",chr_ref)
    pos_ref<-as.numeric(strsplit(data$snp[j],split=":")[[1]][2])
    data$reference_allele[j]<-as.character(getSeq(BSgenome.Rnorvegicus.UCSC.rn6, chr_ref, start=pos_ref, end=pos_ref) )
  }

  #getSeq(BSgenome.Rnorvegicus.UCSC.rn6, 'chr13', start=75076290, end=75076290)
  return(data)
}



results_dir = paste0(exp_dir,"results/")
#allele1 = minor allele
#allele0= major allele
#V1  V2 V3 V4 V5   V6 V7 V8 V9
#Chr SNP bp A1 A2 Freq  b se  p
#Referenceallele=A1         Otherallele=A2
readAssoc<-function(x){

  data<-fread(paste0('cat ',x,' | cut -f 2,4,5 '),header=T,stringsAsFactors =F)
  column_names<-c("snp","allele1","allele2")
  setnames(data, column_names)
  return(data)

}





qtls$assoc<- paste0(exp_dir,"results/gwas/",qtls$chr,".",qtls$trait,".mlma")
for(i in 1:nrow(qtls)){

  trait<-qtls$trait[i]
  topsnp<-qtls$topsnp[i]
  chr<-qtls$chr[i]
  pos<-qtls$pos[i]
  assoc_name = qtls$assoc[i]

  system(paste0("/oasis/tscc/scratch/aschitre/software/plink-1.90/plink --bfile ",exp_dir,"genotypes/genotypes --chr ",gsub("chr","",chr)," --nonfounders --r2  --ld-snp ",topsnp," --ld-window 100000 -ld-window-kb 6000 --ld-window-r2 0 --out ",exp_dir,"temp/r2/var_anno"),ignore.stdout = T,ignore.stderr = T,wait = T)


  LD_df<-LD_ref_snps(paste0(exp_dir,"temp/r2/var_anno.ld"),trait = trait,chr=chr,pos=pos)


  assoc_df<- readAssoc(assoc_name)

  merged<-merge(LD_df,assoc_df,by="snp",all=F)


  merged$vep_allele<-NA

  for(k in 1:nrow(merged)){

    chr_ref<-strsplit(merged$snp[k],split=":")[[1]][1]
    chr_ref=paste0("",chr_ref)
    pos_ref<-as.numeric(strsplit(merged$snp[k],split=":")[[1]][2])

    if(merged$reference_allele[k] %in% merged$allele1[k]){
      merged$vep_allele[k]<-merged$allele2[k]

    }else{
      merged$vep_allele[k]<-merged$allele1[k]

    }






  }

  ##CHROM  POS ID      REF ALT QUAL    FILTER  INFO
  setwd(paste0("/oasis/tscc/scratch/aschitre/round8/unpruned/snp_annotations/snpeff_input/",project_name))
  colnames(merged)[grepl("vep_allele",colnames(merged))]<-"ALT"

  merged$CHROM<-NA
  merged$POS<-NA
  merged$CHROM<-sapply(strsplit(merged$snp,split=":"),`[`,1)
  merged$POS<-sapply(strsplit(merged$snp,split=":"),`[`,2)
  merged$ID<-paste0(merged$CHR,"_",merged$POS)
  #/home/apurva/ensembl-vep

  colnames(merged)[grepl("reference_allele",colnames(merged))]<-"REF"


  colnames(merged)[grepl("CHROM",colnames(merged))]<-"##CHROM"

  #addr = paste(name,  address, cityState, sep="\n")
  #vcf_header=paste(name,  address, cityState, sep="\n")

  cat("##fileformat=VCFv4.0 \n",file=paste0(trait,"_",chr,"_",pos,".vcf"))


  #a = matrix(1:9, nrow = 3, ncol = 3, dimnames = list(LETTERS[1:3], LETTERS[1:3]))
  #cat("World \n",file="a.csv")
  #write.table(a, 'a.csv',sep=",",append=TRUE, col.names=NA)

  merged$QUAL<-"."
  merged$FILTER<-"."
  merged$INFO<-"."

  write.table(merged[,c("##CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")],paste0(trait,"_",chr,"_",pos,".vcf"),col.names=T,append=T,row.names=F,quote=F,sep="\t",na = ".")


  #fConn <- file(paste0(trait,"_",chr,"_",pos,".txt"), 'r+')
  #Lines <- readLines(fConn)
  #writeLines(paste("Text at beginning of file", Lines, sep = "\n"),con = fConn)



}




#run snpeff for one file
dir.create(paste0("/oasis/tscc/scratch/aschitre/round8/unpruned/snp_annotations/snpeff_output/",project_name), showWarnings = FALSE)



#/usr/lib/jvm/jre-1.8.0/bin/java -Xms40g -Xmx40g -jar snpEff.jar Rnor_6.0.99  -no-intergenic -no-intron -noStats /oasis/tscc/scratch/aschitre/round8/unpruned/snp_annotations/snpeff_input/ccp/ccp_post_test_time_on_cs_chr13_75349353.vcf > ccp_post_test_time_on_cs_chr13_75349353.ann.vcf


#for file in list of files, create these commands
#paste0(trait,"_",chr,"_",pos,".vcf")
qtls$vcf_snpeff<-paste0("/oasis/tscc/scratch/aschitre/round8/unpruned/snp_annotations/snpeff_input/",project_name,"/",qtls$trait,"_",qtls$chr,"_",qtls$pos,".vcf")

if (file.exists(paste0("/projects/ps-palmer/software/local/src/snpEff/",project_name,".sh"))) {
  file.remove(paste0("/projects/ps-palmer/software/local/src/snpEff/",project_name,".sh"))
}



for(i in 1:nrow(qtls)){
  input<-paste0("/oasis/tscc/scratch/aschitre/round8/unpruned/snp_annotations/snpeff_input/",project_name,"/",qtls$trait[i],"_",qtls$chr[i],"_",qtls$pos[i],".vcf")
  part1<-"/usr/lib/jvm/jre-1.8.0/bin/java -Xms40g -Xmx40g -jar /projects/ps-palmer/software/local/src/snpEff/snpEff.jar Rnor_6.0.99  -no-intergenic -no-intron -noStats "
  part2<-" > "
  output<-paste0("/oasis/tscc/scratch/aschitre/round8/unpruned/snp_annotations/snpeff_output/",project_name,"/",qtls$trait[i],"_",qtls$chr[i],"_",qtls$pos[i],"_anno.vcf")

  command<- paste0(part1,input,part2,output)

  #snpeff<-paste(dir_snpeff,command,sep="\n")
  write(command,file=paste0("/projects/ps-palmer/software/local/src/snpEff/",project_name,".sh"),append=TRUE)

  #addr = paste(name,  address, cityState, sep="\n")

}


#fConn <- file(paste0(trait,"_",chr,"_",pos,".txt"), 'r+')
#Lines <- readLines(fConn)
#writeLines(paste("Text at beginning of file", Lines, sep = "\n"),con = fConn)

#system(paste0("cp ",project_name,".sh /projects/ps-palmer/software/local/src/snpEff/"))
system(paste0("chmod +x /projects/ps-palmer/software/local/src/snpEff/",project_name,".sh"))
system(paste0("/projects/ps-palmer/software/local/src/snpEff/",project_name,".sh"),wait=T)

#run .sh command on an interactive node

if (file.exists(paste0("/projects/ps-palmer/software/local/src/snpEff/",project_name,"_filter.sh"))) {
  file.remove(paste0("/projects/ps-palmer/software/local/src/snpEff/",project_name,"_filter.sh"))
}

for(i in 1:nrow(qtls)){

  part1<-paste0("cat /oasis/tscc/scratch/aschitre/round8/unpruned/snp_annotations/snpeff_output/",project_name,"/",qtls$trait[i],"_",qtls$chr[i],"_",qtls$pos[i],"_anno.vcf \\")
  part2<-paste0("    | /projects/ps-palmer/software/local/src/snpEff/scripts/vcfEffOnePerLine.pl \\")
  part3<-paste0("| /usr/lib/jvm/jre-1.8.0/bin/java -Xms512m -Xmx512m -jar /projects/ps-palmer/software/local/src/snpEff/SnpSift.jar extractFields - CHROM POS REF ALT \"ANN[*].ALLELE\" \"ANN[*].ANNOTATION\" \"ANN[*].IMPACT\" \"ANN[*].GENE\" \"ANN[*].GENEID\" \"ANN[*].FEATURE\" \"ANN[*].FEATUREID\" \"ANN[*].BIOTYPE\" \"ANN[*].RANK\" \"ANN[*].HGVS_C\" \"ANN[*].HGVS_P\" \"ANN[*].CDNA_POS\" \"ANN[*].CDNA_LEN\" \"ANN[*].CDS_POS\" \"ANN[*].CDS_LEN\" \"ANN[*].AA_POS\" \"ANN[*].AA_LEN\" \"ANN[*].DISTANCE\" \"ANN[*].ERRORS\" >")
  output<- paste0("/oasis/tscc/scratch/aschitre/round8/unpruned/snp_annotations/snpeff_output/",project_name,"/",qtls$trait[i],"_",qtls$chr[i],"_",qtls$pos[i],"_fields.txt")
  output<-paste0(part3,output)
  command<-paste(part1,part2,output,sep = "\n")
  write(command,file=paste0("/projects/ps-palmer/software/local/src/snpEff/",project_name,"_filter.sh"),append=TRUE)




  #addr = paste(name,  address, cityState, sep="\n")

}



#system(paste0("cp ",project_name,"_filter.sh /projects/ps-palmer/software/local/src/snpEff/"))
system(paste0("chmod +x /projects/ps-palmer/software/local/src/snpEff/",project_name,"_filter.sh"))
system(paste0("/projects/ps-palmer/software/local/src/snpEff/",project_name,"_filter.sh"),wait=T)
#setwd(paste0("/projects/ps-palmer/software/local/src/snpEff/"))
#system(paste0("./",project_name,"_filter.sh"))

.libPaths(c("/home/aschitre/R_libs", .libPaths()))
library("data.table", lib.loc="/home/aschitre/R_libs/")

setwd(paste0("/oasis/tscc/scratch/aschitre/round8/unpruned/snp_annotations/snpeff_output/",project_name))

#header<-read.table("/oasis/tscc/scratch/aschitre/round8/assoc_header.txt",header=F,stringsAsFactors = F)

anno_files<-list.files(path=".",full.names = F,pattern=paste0("*_fields.txt_*"))


readdata<-function(x){

  data<-fread(file=x,header=T,stringsAsFactors =F,na.strings = "")
  column_names=c('chr','pos','ref','alt','snp_allele','snp_effect','snp_impact','snp_gene','snp_geneid','snp_feature','snp_featureid','snp_biotype','snp_rank','snp_hgvs_c','snp_hgvs_p','snp_cdna_pos','snp_cdna_length','snp_cds_pos','snp_cds_len','snp_aa_pos','snp_aa_length','snp_dist_to_feature','snp_errors')
  setnames(data, column_names)
  data=data[snp_impact %in% c("HIGH","MODERATE")]
  trait_snp<-gsub("_fields.txt","",x)
  data$trait_snp<-trait_snp
  return(data)

}

anno_snpeff<-lapply(anno_files,readdata)
anno_files<-gsub("_fields.txt","",anno_files)
names(anno_snpeff)<-anno_files


annotations<-do.call("rbind",anno_snpeff)


write.csv(annotations,paste0(exp_dir,"results/snp_annotations_final_withWarnings.csv"),row.names=F,quote=T)



annotations<-annotations[is.na(annotations$snp_errors),]



annotations$gene_name<-NA
annotations$RGD_link<-NA
annotations$gene_long_name<-NA

#genedesc <- getBM(attributes=c('external_gene_name','description'), filters = 'ensembl_gene_id', values = annotations$snp_geneid[i], mart =ensembl)

ensembl_genes<-unique(annotations$snp_geneid)
ensembl_genes<-ensembl_genes[!is.na(ensembl_genes)]



.libPaths(c("/home/aschitre/R_libs", .libPaths()))

library(vctrs, lib.loc="/home/aschitre/R_libs/")
library(annotables, lib.loc="/home/aschitre/R_libs/")


if(nrow(rnor6[which(rnor6$ensgene %in% ensembl_genes),])>0){

  for(i in 1:length(ensembl_genes))  {
    if(ensembl_genes[i] %in% rnor6$ensgene){
      annotations$gene_name[which(annotations$snp_geneid %in% ensembl_genes[i])]<-as.character(rnor6[which(rnor6$ensgene %in% ensembl_genes[i]),"symbol"] )
      #rnor6[which(rnor6$ensgene %in% ensembl_genes[i]),"description"]
      acc<-gsub("]","",strsplit(as.character(rnor6[which(rnor6$ensgene %in% ensembl_genes[i]),"description"]),split="Acc:")[[1]][2])
      if(is.na(acc)){
        annotations$RGD_link[which(annotations$snp_geneid %in% ensembl_genes[i])]<-NA
        annotations$gene_long_name[which(annotations$snp_geneid %in% ensembl_genes[i])]<-NA


      }else{
        annotations$RGD_link[which(annotations$snp_geneid %in% ensembl_genes[i])]<-paste0("https://rgd.mcw.edu/rgdweb/report/gene/main.html?id=",acc)
        annotations$gene_long_name[which(annotations$snp_geneid %in% ensembl_genes[i])]<-strsplit(as.character(rnor6[which(rnor6$ensgene %in% ensembl_genes[i]),"description"]),split=" \\[")[[1]][1]

      }


    }




  }




}

















#genedesc <- getBM(attributes=c('external_gene_name','description','ensembl_gene_id'), filters = 'ensembl_gene_id', values = ensembl_gene_id, mart =ensembl)

#colnames(genedesc)[which(colnames(genedesc) %in% "ensembl_gene_id")]<-"snp_geneid"

#annotations<-merge(annotations,genedesc,by="snp_geneid",all=T)

#for(i in seq_along(annotations$snp_geneid)){

#  if(is.na(annotations$description[i])){

#    annotations$RGD_link[i]<-"No RGD entry"

#  }else{
#    acc<-gsub("]","",strsplit(annotations$description[i],split="Acc:")[[1]][2])
#    annotations$RGD_link[i]<-paste0("https://rgd.mcw.edu/rgdweb/report/gene/main.html?id=",acc)


#  }



#}


annotations$r2_with_trait_topsnp<-NA
annotations$dprime_with_trait_topsnp<-NA


for(i in 1:nrow(annotations)){
  topsnp<-gsub(".*_chr","",annotations$trait_snp[i])
  chr<-strsplit(topsnp,split="_")[[1]][1]
  trait_topsnp<-paste0("chr",chr,":",strsplit(topsnp,split="_")[[1]][2])
  snp<-paste0(annotations$chr[i],":",annotations$pos[i])
  command<-paste0("/oasis/tscc/scratch/aschitre/software/plink-1.90/plink --bfile ",exp_dir,"genotypes/genotypes --chr ",chr," --nonfounders --r2  --ld ",trait_topsnp," ",snp," --out ",exp_dir,"temp/r2/temp_snpeff")
  system(command,wait=T)
  r2_log=fread(paste0("cat ",exp_dir,"temp/r2/temp_snpeff.log | grep R-sq"),header=F,stringsAsFactors = F)
  if(nrow(r2_log)<2){


    annotations$r2_with_trait_topsnp[i]<-r2_log$V3
    annotations$dprime_with_trait_topsnp[i]<-r2_log$V6
  }else{
    annotations$r2_with_trait_topsnp[i]<-as.numeric(r2_log[which.max(r2_log$V3),"V3"])
    annotations$dprime_with_trait_topsnp[i]<-as.numeric(r2_log[which.max(r2_log$V3),"V6"])
  }

}

annotations$r2_with_trait_topsnp<-unlist(annotations$r2_with_trait_topsnp)
annotations$dprime_with_trait_topsnp<-unlist(annotations$dprime_with_trait_topsnp)
write.csv(annotations,paste0(exp_dir,"results/snp_annotations_final.csv"),row.names=F,quote=T)
