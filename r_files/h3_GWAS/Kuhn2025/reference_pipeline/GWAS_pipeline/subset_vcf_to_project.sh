module load bcftools

infile=/projects/ps-palmer/hs_rats/round9_1/Heterogenous-stock_n14780_10182022_QC_Sex_Het_pass_n13548.vcf.gz


cd $PBS_O_WORKDIR

project_name=basename "$PWD"



bcftools view --force-samples -S ./genotypes/rfids_bcftools.txt $infile > ./genotypes/filtered.vcf
