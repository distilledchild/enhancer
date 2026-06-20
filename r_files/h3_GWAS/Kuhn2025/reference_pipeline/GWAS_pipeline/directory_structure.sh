project_name="$1"


cd "$PWD"



if [ -d $project_name ]
then
    echo "Error: Directory already exists. Please use a different name"
else
    mkdir -p  "$project_name"

    cd $project_name

	mkdir -p ./results/{eqtl,phewas,gwas,manhattan_plots,qtls,locuszoom_plots,heritability,genetic_correlation}
	mkdir -p ./data/{raw_data,pheno,residuals,gwas_jobs}
	mkdir -p ./pheno_processing_summary/{plots/{age,traits},age,covs,misc,N}
	mkdir -p ./code
	mkdir -p ./log
	mkdir -p ./temp/{assoc,code,pbs_log,pheno,r2,conditional_analysis}
	mkdir -p ./genotypes/{grm,dosages}
	mkdir -p ./grm


echo "Directory structure created"



fi
