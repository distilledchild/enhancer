
module load plink

infile=/projects/ps-palmer/hs_rats/round9_1/Heterogenous-stock_n14780_10182022_QC_Sex_Het_pass_n13548.vcf.gz
outfile=/projects/ps-palmer/apurva/riptide/genotypes/Heterogenous-stock_n14780_10182022_QC_Sex_Het_pass_n13548



plink --vcf $infile --chr 1-20 --double-id --set-missing-var-ids @:# --make-bed --out $outfile


plink --bfile $outfile --geno 0.1 --maf 0.005 --hwe 1e-10 --make-bed --out /projects/ps-palmer/apurva/riptide/genotypes/round9_1



#subset to project samples
plink --bfile /projects/ps-palmer/apurva/riptide/genotypes/round9_1 --keep /projects/ps-palmer/apurva/p50_david_dietz_2020/genotypes/rfids_n2490.txt --make-bed --out /projects/ps-palmer/apurva/p50_david_dietz_2020/genotypes/genotypes


#calculate GRM for all chromosomes and for subtract GRM

for i in {1..20}
do
/projects/ps-palmer/software/local/src/gcta/gcta64 --thread-num 4 --bfile /projects/ps-palmer/apurva/p50_david_dietz_2020/genotypes/genotypes --chr $i --make-grm-bin --out /projects/ps-palmer/apurva/p50_david_dietz_2020/grm/chr${i}.genotypes
done

/projects/ps-palmer/software/local/src/gcta/gcta64 --thread-num 4 --bfile /projects/ps-palmer/apurva/p50_david_dietz_2020/genotypes/genotypes --autosome-num 20 --autosome --thread-num 1 --make-grm-bin --out /projects/ps-palmer/apurva/p50_david_dietz_2020/grm/genotypes
