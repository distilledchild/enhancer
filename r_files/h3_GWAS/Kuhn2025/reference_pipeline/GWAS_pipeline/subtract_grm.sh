#usage :  qsub -q condo -l nodes=1:ppn=6 -l walltime=3:00:00  subtract_grm.sh


for i in {1..20}
do
/projects/ps-palmer/software/local/src/gcta/gcta64 --thread-num 4 --bfile /projects/ps-palmer/apurva/p50_david_dietz_2020/genotypes/genotypes --chr $i --make-grm-bin --out /projects/ps-palmer/apurva/p50_david_dietz_2020/grm/chr${i}.genotypes
done
