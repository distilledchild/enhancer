#step1 prepare locus zoom snp set input file
cd genotypes

module load plink
plink --bfile genotypes --write-snplist --out genotypes


#cp to locus zoom bin directory
cp genotypes.snplist /projects/ps-palmer/software/local/src/locuszoom/bin/genotypes.snplist

#edit above file
#start R session
library(data.table)
snps=fread(file="/projects/ps-palmer/software/local/src/locuszoom/bin/genotypes.snplist",header=F,stringsAsFactors=F,sep=" ")

snps[,c("chr","pos") := tstrsplit(V1, ":", fixed=TRUE)]
colnames(snps)[1]<-"snp"
snps$chr<-gsub("chr","",snps$chr)
write.table(snps,"/projects/ps-palmer/software/local/src/locuszoom/bin/genotypes_snp_pos_file.txt",row.names=F,sep = "\t",quote=F)
q()

#step 2 Prepare refFlat file
http://genome.ucsc.edu/cgi-bin/hgTables?db=rn6&hgta_group=genes&hgta_track=refSeqComposite&hgta_table=ncbiRefSeqCurated

#use same P50 RP4 locuszoom refflat file
ori<-fread(file="/projects/ps-palmer/software/local/src/locuszoom/bin/rn6_refFlat_13Sep2020_raw.txt",stringsAsFactors = F,header=T)

head(ori)
#name2=geneName
#name=name
#chrom=chrom
#strand=strand
#txStart=txStart
#txEnd=txEnd
#cdsStart=cdsStart
#cdsEnd=cdsEnd
#exonCount=exonCount
#exonStarts=exonStarts
#exonEnds=exonEnds



colnames(ori)[which(colnames(ori)=="name2")]<-"geneName"
colnames(ori)
write.table(ori[,c("geneName","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds")],"/projects/ps-palmer/software/local/src/locuszoom/bin/rn6_refFlat_13Sep2020.txt",row.names=F,sep = "\t",quote=F)
q()

#step3 create custom locus zoom database
CWD=$(pwd)
project_name=basename $CWD
module load python

cd /projects/ps-palmer/software/local/src/locuszoom/bin

./dbmeister.py --db $project_name.db --snp_pos genotypes_snp_pos_file.txt

./dbmeister.py --db $project_name.db --refflat rn6_refFlat_13Sep2020.txt


#step4 modifying LZ conf file
cd /projects/ps-palmer/software/local/src/locuszoom/conf
#edit the conf file using vi. Add the following line to SQLITE_DB portion of the conf file
#vi commands: esc + i = insert a line
#esc :wq! = save the edits
vi m2zfast.conf

 'cfw_zou' : "/projects/ps-palmer/software/local/src/locuszoom/bin/cfw_zou.db"
'p50_david_dietz_2020' : "/projects/ps-palmer/software/local/src/locuszoom/bin/p50_david_dietz_2020.db"
