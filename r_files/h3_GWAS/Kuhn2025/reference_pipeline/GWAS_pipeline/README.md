Within this pipeline, the following software is used:

• PLINK v1.90p 64-bit (2 Jun 2015) gcta-1.94.1

• R version 4.1.2

R packages used can be installed using conda.
Please contact Faith Okamoto for R libraries and other dependencies

------------------------------------------------------------------------

-   This document describes what each script and step of the GWAS
    analysis pipeline does. For detailed commands, check out the
    **all_commands.sh** document.
-   There are 2 main checkpoints in the pipeline:

1.  **checkpoint 1**: After processing phenotypes, look at the
    percentage of variance explained by different covariates. If
    everything looks ok, proceed to the next step.
2.  **checkpoint 2**: If the SNP heritability estimates look ok, submit
    the next steps.

------------------------------------------------------------------------

The computational requirements for each step are based on a project with
N = 2,300, unpruned set of SNPs (\~6M SNPs) and 90 traits.

------------------------------------------------------------------------

1.  **Clone the git repo**

------------------------------------------------------------------------

2.  **Make all files in the code directory executable**
    This is the directory that contains the GWAS scripts

``` bash
chmod +x /projects/ps-palmer/apurva/genetic_analysis/code

```

------------------------------------------------------------------------

3.  **Create directory structure**
    This will create the directory structure for the project
    *p50_david_dietz_2020*

``` bash
/projects/ps-palmer/apurva/genetic_analysis/code/directory_structure.sh p50_david_dietz_2020
```

------------------------------------------------------------------------

4.  **Copy csv file with raw phenotypes and trait descriptions to
    data/raw_data**

-   Ideally, we’d like to query the **gwas_phenotypes** and the
    **descriptions** tables from the database for GWAS analysis. We’ll
    need to upload the tables in the database regularly to make this
    possible. For now, we’ll copy these files to the
    **project_name/data/raw_data** folder.
-   The project metadata \[PI name, project title etc\] should also be
    queried from the database. This information is in the
    **project_metadata** table in the sample_tracking schema. I’m
    storing this information in the **traits.RData** object for now.
-   The **process_phenotypes.R** script will look for *raw_data.csv* and
    the *traits.RData* (data dictionary) in the
    **project_name/data/raw_data** folder.

Here’s what the **process_phenotypes.R** script does:
- removes white spaces in any trait or covariate values and sets them to
NA
- quantile normalize the traits within each sex and center/project_name
value. (This is similar to using sex as a covariate)

There are 2 sets of covariates:
- common covariates : coat color, cohort number. These are not specific
to an experiment.
- covariates specific to an experiment : experiment box, age at
experiment. These should have the following naming convention: append
the experiment string in front of the covariate. For example, “loco_age”
will be used as a covariate for these locomotor traits: “loco_center”,
“loco_rear”, “loco_distance”, “loco_activity”

*output: data/residuals/residuals.RData*

**Before proceeding to the next steps, check the PVE explained by the
different covariates. These RData frames are saved in
pheno_processing_summary/covs folder**

------------------------------------------------------------------------

5.  **Subset genotype data and calculate GRM, SNP h2 estimates**

-   *genotype_filters_on_vcf.sh* This will prepare the genotype data for
    GWAS. If we’re planning to get the PLINK binary file as an output of
    the genotyping pipeline, then we’ll only need to subset the file to
    the project rfids. This can be done using the *–keep* option in
    PLINK.
-   *create_famid.R* This will create the sample id list. This list is
    an input to the –keep command for *genotype_filters_on_vcf.sh*
-   **Calculate 2 sets of GRMs**:
-   *all_chromosomes_grm.sh*: GRM estimated from all SNPs on all
    chromosomes.
-   *subtract_grm.sh*: calculated from SNPs on one chromosome. This GRM
    will be used for *–mlma-subtract-grm* option designed to parallelise
    the MLMA-LOCO analysis for large data set.
    *More info*:
    <https://gcta.freeforums.net/thread/173/mixed-linear-model-association-analysis> -
    *LZ_prep.sh* detailed steps for preparing the locuszoom backend
    SQLite database are outlined in this document:
    <https://www.dropbox.com/s/2ugc710e2sp9v98/locuszoom_standalone_instructions.pdf?dl=0>

``` bash
Rscript create_famid.R

qsub -q condo -W depend=afterok:$geno_filters -l nodes=1:ppn=6 -l walltime=8:00:00  /projects/ps-palmer/apurva/genetic_analysis/code/extract_dosages.sh

geno_filters=`qsub -q condo -l nodes=1:ppn=6 -l walltime=3:00:00  /projects/ps-palmer/apurva/genetic_analysis/code/genotype_filters_on_vcf.sh`

qsub -q condo -W depend=afterok:$geno_filters -l nodes=1:ppn=6 -l walltime=3:00:00  /projects/ps-palmer/apurva/genetic_analysis/code//subtract_grm.sh

all_grm=`qsub -q condo -W depend=afterok:$geno_filters -l nodes=1:ppn=6 -l walltime=3:00:00  /projects/ps-palmer/apurva/genetic_analysis/code/all_chromosomes_grm.sh`

qsub -q condo -W depend=afterok:$all_grm -l nodes=1:ppn=2 -l walltime=1:00:00 /projects/ps-palmer/apurva/genetic_analysis/code/snp_heritability.sh
```

------------------------------------------------------------------------

6.  **Submit GWAS jobs**

-   We’re using MLMA-LOCO analysis for large data sets in GCTA to
    perform mapping.
    *More info*:
    <https://gcta.freeforums.net/thread/173/mixed-linear-model-association-analysis>
-   I’ve included 2 options for submitting GWAS jobs in the
    **all_commands.sh** file.
-   Option 2 is better since it’s easier to manage array jobs.
-   Option 2 requires phenotypes.txt file in the data directory.
    **create_array_files.R** will create this file.

``` bash
#example command
#gcta64 --mlma --grm test_all --mlma-subtract-grm test_chr1 --bfile test --chr 1 --pheno test.phen --out test_loco_chr1 --thread-num 10

##This will create phenotypes.txt file required to submit an array job
Rscript create_array_files.R

num_phenos=`wc -l < data/phenotypes.txt`
gwas_jobs=`qsub -q hotel -t 1-$num_phenos -l nodes=1:ppn=5 -l walltime=10:00:00 /projects/ps-palmer/apurva/genetic_analysis/code/gwas_loco_array.sh`
```

------------------------------------------------------------------------

7.  **Annotate and visualize GWAS results**
-   All downstream steps in the pipeline are dependent on the completion
    of the GWAS jobs.
-   These include: prepping PheWAS database, calling QTLs, submitting
    Manhattan plot jobs, creating regional association plots (Locuszoom)
    etc.
-   Algorithm for calling QTLs (calling_qtls_manyTraits_subtractGRM.sh)
    is described in Genetic mapping section of this paper:
    <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7511439/>
-   If the algorithm for calling QTLs, calls multiple QTLs on the same
    chromosome, then we perform conditional analysis to establish
    independence.
-   *table3.R* : This script annotates the QTLs. (Defines LD r2 based
    intervals, gets genes from Rat Genome Database, adds Strain
    Distribution Pattern (SDP) for HS founders
-   We used SnpEff to perform variant annotation.

------------------------------------------------------------------------

8.  **Compress images, knit report and create minio public link**
-   I convert the Locuszoom PDFs to png and compress the Manhattan
    plots using command-line tools.
-   After knitting the report, copy the html file to the **/projects/ps-palmer/s3/data/** directory.
    p50_reports = directory for all P50 projects. u01_reports = U01 projects.
-   Please follow the pinned document on TSCC channel for minio access key.

``` bash
convert_locuszoom_to_png.sh
#I converted LZ PDFs locally using pdftoppm. Found the following package that looks like a wrapper around pdftoppm.
#can be installed using conda
#try https://anaconda.org/conda-forge/pdf2image


#compress manhattan plots Done locally on my Desktop
compress_manhattan_plots.sh

#knit report on the cluster
#currently knitting the Rmd locally
p50_david_dietz_2020.Rmd


#copy the html file to the /projects/ps-palmer/s3/data/ directory
#This will generate web-link for the report
./mc anonymous set public minio/p50_reports --recursive
./mc anonymous links minio/p50_reports --recursive
```
