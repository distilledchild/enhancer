# HRDP genotype diversity analysis

This analysis addresses the reviewer's request to define the genetic diversity
of the ten strains represented by the Hi-C libraries. It is not designed to
test strain-specific loop differences, genotype-loop associations, or whether
these ten strains capture the full genetic diversity of the HRDP.

The analysis uses autosomal, biallelic, variable SNPs with complete genotypes
in all ten libraries. Native PLINK 2.0 performs LD pruning, pairwise IBS-distance
calculation, and PCA. Since only ten samples are analyzed, LD pruning and PCA
use PLINK 2.0's `--bad-ld` and `--bad-freqs` small-sample overrides.
`--indep-order 1` fixes the pruning order for reproducibility.

Open `hrdp_genotype_diversity.R` and execute it from top to bottom. Set the R
working directory to either the `enhancer` repository root or this analysis
directory before running the file. The complete script can also be run with:

```bash
Rscript r_files/revision/hrdp_genotype_diversity/hrdp_genotype_diversity.R
```

On Apple Silicon, the R file downloads the official native ARM64 PLINK 2.0
alpha 7.3 binary to `.tools/plink2` when it is not already available; Rosetta
and a Conda environment are not required. The local tool directory is ignored
by Git, and primary outputs are written to `results/plink2_primary`. On other
platforms, set `PLINK2_BIN` to the platform-specific PLINK 2.0 executable
before running the R file.
