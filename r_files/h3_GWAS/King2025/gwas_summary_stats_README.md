# King et al. 2025 GWAS Summary Statistics

Source dataset: UCSD Digital Collections object `bb5030313v`, "Data from: Genetic Loci Influencing Cue-Reactivity in Heterogeneous Stock Rats".

Primary public report URL:
https://palmerlab.s3.sdsc.edu/tsanches_dash_genotypes/gwas_results/p50_paul_meyer_2014_rn6mfsplit/results/gwas_report.html

Downloaded files:
- `raw_chrgwas_mlma/`: chromosome-wise GCTA MLMA files for 59 QTL-positive traits, autosomes 1-20.
- `king_qtl_positive_autosome_mlma_manifest.tsv`: download manifest for the 1,180 `.mlma` files.
- `trait_n_from_heritability.tsv`: MAGMA-ready sample size lookup for each downloaded trait.

Related extracted/public files are in:
- `../ucsd_bb5030313v/extracted_tables/`
- `../ucsd_bb5030313v/s3_downloads/`

Notes:
- The `.mlma` files contain `Chr`, `SNP`, `bp`, `A1`, `A2`, `Freq`, `b`, `se`, and `p` columns.
- For MAGMA, use `SNP` and `p` from the merged chromosome-wise MLMA file for a trait, and use the corresponding `n` value in `trait_n_from_heritability.tsv`.
- The public report also includes phenotype/covariate tables and top-QTL summaries. The full genome-wide SNP-level summary statistics are available as chromosome-wise MLMA files on S3 rather than as one combined file per trait.
