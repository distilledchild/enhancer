
cd $PBS_O_WORKDIR



project_name="$(basename $PWD)"



find ./results/gwas -name "*.mlma" -exec awk ' $9 <=0.0001 {print $2,$9,FILENAME}' {} \; > /projects/ps-palmer/phewas_db/$project_name.txt

cat /projects/ps-palmer/apurva/u01_phewas_db/u01_phewas_db.txt /projects/ps-palmer/phewas_db/$project_name.txt > /projects/ps-palmer/phewas_db/riptide_gwas.txt
