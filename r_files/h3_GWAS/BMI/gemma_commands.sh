#Constructing GRM
#Since we use the LOCO (leave-one-chromosome-out) method when we do GWAS, we use a relatedness matrix (-k) that has been formed using all the SNPs EXCLUDING those on the chromosome currently being tested.
#We need 20 dosage files [bimbam file format], each missing a chromosome of SNPs, in order to make the matrix. As well as 20 individual dosage files for when you do the testing chromosome by chromosome.
#I found this to be very easily accomplished with a simple grep command. If you use grep -v "${chrom}\." it will give you the original file, minus the chromosome you want to leave out.

gemma -g /oasis/tscc/scratch/aschitre/round2_unpruned/dosages/allExcept.${chrom}.P50_round2_3473_unpruned.bimbam -p bodyweight.txt -gk 1 -o allExcept.${chrom}.round2_unpruned_3473.dosages

#GWAS command
#-g = genotype file
#-p = phenotype file and should just be 1 columns per phenotype in the same sample order as the genotype file, NO HEADER
#-k is the relatedness matrix that you make in the first step
#-lmm 4 tells it to run all the tests and give you a beta estimate for the SNP effect.
#-a is a SNPINFO file that looks like this
#-o is just the output prefix


gemma -g ${chrom}.round2_impute2_3473.bimbam -p round2_physiological_bmi_bodylength_w_tail.txt -n 1 -k allExcept.${chrom}.round2_impute2_LDpruned0.95_3473.dosages.cXX.txt  -a ${chrom}.round2_128447.snpinfo -lmm 4 -o ${chrom}.round2_impute2_LDpruned0.95.gwas.physiological_bmi_bodylength_w_tail
