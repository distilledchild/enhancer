# HRDP genotype input

Place `hrdp_genotype_by_chr.tgz` in this directory. The analysis extracts the
ten mapped strain genotypes, retains complete autosomal biallelic variable
SNPs, and uses PLINK 2.0 for LD pruning, IBS distance, and PCA.

This large archive is staged once in Google Drive and excluded from
post-commit deletion.
