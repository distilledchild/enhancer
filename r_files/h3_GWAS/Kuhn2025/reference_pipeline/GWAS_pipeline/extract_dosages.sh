module load bcftools

infile=/projects/ps-palmer/hs_rats/round9_1/Heterogenous-stock_n14780_10182022_QC_Sex_Het_pass_n13548.vcf.gz

CWD=$(pwd)



bcftools query -f "%CHROM\t%POS[\t%DS]\n" $infile > $CWD/genotypes/dosages/dosages.txt

bcftools query -l $infile > $CWD/genotypes/dosages/sample_ids.txt
