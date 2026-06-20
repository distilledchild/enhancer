#cd to project directory before executing all the steps
#if there is no qsub command in front of the script, it was not a computationally intensive step and I did not queue it on the cluster.
module load R #Rscript requires this

#Step 1
#Process phenotypes

Rscript process_phenotypes.R


#########################################################################################


#check the PVE explained by covariates. These RData frames are saved in pheno_processing_summary/covs folder
#if everything looks ok, proceed to the next steps

#########################################################################################


#Step 2
#Prep genotype data, extract dosages for conditional analysis, estimate GRMs, prepare locuszoom standalone database for regional association plots
#These can be run as bash job dependencies:

Rscript create_famid.R

geno_filters=`qsub -q condo -l nodes=1:ppn=6 -l walltime=3:00:00  /projects/ps-palmer/apurva/genetic_analysis/code/genotype_filters_on_vcf.sh`
qsub -q condo -W depend=afterok:$geno_filters -l nodes=1:ppn=6 -l walltime=8:00:00  /projects/ps-palmer/apurva/genetic_analysis/code/extract_dosages.sh

qsub -q hotel -l nodes=1:ppn=6 -l walltime=20:00:00  /projects/ps-palmer/apurva/genetic_analysis/code/subset_vcf_to_project.sh
qsub -q condo -l nodes=1:ppn=6 -l walltime=8:00:00  /projects/ps-palmer/apurva/genetic_analysis/code/extract_dosages.sh
qsub -q condo -W depend=afterok:$geno_filters -l nodes=1:ppn=6 -l walltime=3:00:00  /projects/ps-palmer/apurva/genetic_analysis/code//subtract_grm.sh
all_grm=`qsub -q condo -W depend=afterok:$geno_filters -l nodes=1:ppn=6 -l walltime=3:00:00  /projects/ps-palmer/apurva/genetic_analysis/code/all_chromosomes_grm.sh`
qsub -q condo -W depend=afterok:$all_grm -l nodes=1:ppn=2 -l walltime=1:00:00 /projects/ps-palmer/apurva/genetic_analysis/code/snp_heritability.sh


#locuszoom standalone database
#This is a series of bash and R commands. It also requires editing conf/m2zfast.conf, and changing the SQLITE_DB variable. For more information check out the locuszoom standalone instructions PDF.
LZ_prep.sh


#########################################################################################


#check the SNP heritability estimates. If these are close to zero, it indicates some issue such as data being scrambled or some issue with how the traits were calculated.
#if everything looks ok, proceed to the next steps

#########################################################################################


#Step 3

Rscript /projects/ps-palmer/apurva/genetic_analysis/code/gwas_jobs_files.R

#depending on the number of traits you're analyzing, you'll need to wait before you submit all the array jobs because you'll hit the max jobs in the queue limit. if less than 75 traits
#Option 1
#Submits GWAS jobs as individual PBS Batch jobs
gwas_jobs=$(Rscript /projects/ps-palmer/apurva/genetic_analysis/code/qsub_gwas_jobs.R)


#loco-mlma command
#Option 2
#Submits GWAS jobs as array jobs

#This will create phenotypes.txt file required to submit an array job
Rscript create_array_files.R

num_phenos=`wc -l < data/phenotypes.txt`
gwas_jobs=`qsub -q hotel -t 1-$num_phenos -l nodes=1:ppn=5 -l walltime=10:00:00 /projects/ps-palmer/apurva/genetic_analysis/code/gwas_loco_array.sh`



#create phewas db for riptide data. This step is dependent on the GWAS jobs.
qsub -q hotel -l nodes=1:ppn=5 -l walltime=10:00:00 /projects/ps-palmer/apurva/genetic_analysis/code/phewas_prep.sh


#manhattan
num_phenos=`wc -l < data/phenotypes.txt`
qsub -q condo -t 1-$num_phenos -l nodes=1:ppn=3 -W depend=afterok:$gwas_jobs -l walltime=3:00:00 -v gwas_threshold=5.73 /projects/ps-palmer/apurva/genetic_analysis/code/Manhattan_manyTraits_subtractGRM.sh

#calling qtls
num_phenos=`wc -l < data/phenotypes.txt`
calling_qtls=`qsub -q condo -t 1-$num_phenos -l nodes=1:ppn=3 -W depend=afterok:$gwas_jobs -l walltime=3:00:00 -v gwas_threshold=5.73 /projects/ps-palmer/apurva/genetic_analysis/code/calling_qtls_manyTraits_subtractGRM.sh`


#Check if conditional analysis required and summarize the QTLs
Rscript /projects/ps-palmer/apurva/genetic_analysis/code/combine_calling_QTL_results.R

#if final QTL doesn't exist then
#do conditional analysis
if [ -f "results/qtls/conditional_analysis.csv" ];
then

# Multiple QTLs exist on the same chromosome for a trait
echo "Submitting conditional analysis jobs"

num_jobs=`wc -l < results/qtls/conditional_analysis.csv`
qsub -q condo -t 1-$num_jobs -l nodes=1:ppn=3 -l walltime=3:00:00 -v gwas_threshold=5.73 /projects/ps-palmer/apurva/genetic_analysis/code/conditional_analysis_using_recodeA.sh

#create final list of QTLs
Rscript /projects/ps-palmer/apurva/genetic_analysis/code/combining_qtls_after_cond_analysis.R

else
#submit all other jobs [eQTL, PheWAS etc]
echo "Conditional analysis not required"
fi


#proceed to LZ, PheWAS, eQTL coloc and variant annotation

#table3
#This code annotates the QTLs. (Defines LD r2 based intervals, gets genes from Rat Genome Database, adds Strain Distribution Pattern (SDP) for HS founders
qsub -q condo -l nodes=1:ppn=3 -l walltime=3:00:00 /projects/ps-palmer/apurva/genetic_analysis/code/table3.sh

#phewas
Rscript phewas_riptide.R

qsub -q hotel -l nodes=1:ppn=2 -l walltime=10:00:00 /projects/ps-palmer/apurva/genetic_analysis/code/phewas_gbs.sh

#combines GBS and Riptide PheWAS results files
Rscript combine_phewas.R

#annotates PheWAS results [annotate = add trait descriptions]
annotate_phewas.R



#LZ
qsub -q hotel -l nodes=1:ppn=2 -l walltime=3:00:00 /projects/ps-palmer/apurva/genetic_analysis/code/locuszoom.sh


#eqtl
#check with Daniel Munro if there's an update to the eQTL results for HS rats
qsub -q hotel -l nodes=1:ppn=2 -l walltime=3:00:00 /projects/ps-palmer/apurva/genetic_analysis/code/eqtl_rat_gtex.sh


#variant annotation
qsub -q hotel -l nodes=1:ppn=4 -l walltime=6:00:00 /projects/ps-palmer/apurva/genetic_analysis/code/variant_annotation_snpeff.sh


#compress images - this might be computationally intensive, I did it within an interactive job
bash convert_locuszoom_to_png.sh ./results/locuszoom_plots ./results/locuszoom_plots/locuszoom_pngs

#compress manhattan plots Done locally on my Desktop
compress_manhattan_plots.sh

#knit report on the cluster
#currently knitting the Rmd locally
p50_david_dietz_2020.Rmd

#clean up intermediate files from project directory to free up storage
#pending!
