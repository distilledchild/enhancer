# Gunturkun 2022 HS Rat GWAS Pipeline

This folder contains the local GWAS pipeline for generating MAGMA-ready
summary statistics from the GeneNetwork `HSNIH-PalmerPublish` phenotype table
and the local HS rat genotype PLINK files.

## Current Status

- Work root:
  `/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022`
- Batch session:
  `screen` session named `gunturkun_gwas`
- Batch log:
  `/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022/logs/full_all_screen_20260617_235230.log`
- Full-all batch start:
  `2026-06-17 23:52:30 CDT`
- Completed validation trait:
  `10424`
- Validated full LOCO output:
  `/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022/results/trait_10424.loco.mlma`
- Validated MAGMA p-value input:
  `/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022/results/trait_10424.magma.tsv`
- `trait_10424` measured runtime:
  `20 minutes 24 seconds`
- `trait_10424.magma.tsv` valid SNP p-values:
  `5,478,632`

The running `full-all` job skips any trait whose
`results/trait_<trait>.magma.tsv` already exists. This makes the pipeline
restartable.

## Expected Runtime

The runtime estimate is based on the measured full LOCO run for trait `10424`
using 10 GCTA threads:

- One trait: about `20 minutes`
- Remaining full-all work after `10424` validation: `22 traits`
- Total expected batch runtime: about `7.5 hours`
- Expected completion window from the `2026-06-17 23:52:30 CDT` batch start:
  about `2026-06-18 07:00-08:30 CDT`

The estimate can drift if the external SSD is busy, the machine sleeps, or
larger-N traits run slower.

## Trait Set

The pipeline prepares GeneNetwork trait IDs `10418-10440`:

- `10418-10423`: Novel object interaction test traits
- `10424-10429`: Open field test traits
- `10430-10440`: Social interaction test traits

Prepared phenotype metadata:

`/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022/data/prepared/trait_metadata.tsv`

Union keep file:

`/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022/data/prepared/target_traits_union.keep`

## Input Data

Phenotype source:

`https://genenetwork.org/api/v_pre1/sample_data/HSNIH-PalmerPublish.csv`

Local phenotype CSV:

`/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022/data/HSNIH-PalmerPublish.csv`

Local genotype prefix:

`/Volumes/external_1000GB_all/playground/enhancer/data/hs_data/v4/HS_genotypes_v4`

Genotype files:

- `/Volumes/external_1000GB_all/playground/enhancer/data/hs_data/v4/HS_genotypes_v4.bed`
- `/Volumes/external_1000GB_all/playground/enhancer/data/hs_data/v4/HS_genotypes_v4.bim`
- `/Volumes/external_1000GB_all/playground/enhancer/data/hs_data/v4/HS_genotypes_v4.fam`

The genotype data are analyzed as rat autosomes with:

```bash
--autosome-num 20 --autosome
```

## Tools

GCTA wrapper:

`/Volumes/external_1000GB_all/playground/enhancer/tools/gcta-1.95.2-macOS-arm64/gcta-1.95.2-macOS-arm64/run.sh`

GCTA binary:

`/Volumes/external_1000GB_all/playground/enhancer/tools/gcta-1.95.2-macOS-arm64/gcta-1.95.2-macOS-arm64/bin/gcta64`

MAGMA binary:

`/Volumes/external_1000GB_all/playground/enhancer/tools/magma`

Python:

`/Users/pete/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3`

macOS detached session tool:

`/usr/bin/screen`

## Scripts

Prepare phenotypes:

`/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022/scripts/prepare_gunturkun2022_traits.py`

Run GCTA GWAS:

`/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022/scripts/run_gcta_gwas.sh`

Convert GCTA MLMA output to MAGMA p-value input:

`/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022/scripts/mlma_to_magma.py`

## Phenotype Processing

The phenotype preparation script:

1. Reads the GeneNetwork `HSNIH-PalmerPublish` CSV.
2. Matches sample IDs to the local PLINK `.fam`.
3. Prepares one phenotype file per trait.
4. Applies sex-stratified rank inverse-normal transformation when sex metadata
   are available.
5. Residualizes categorical covariates
   `batch_within_center`, `coat_color`, and `center` when combined R2 is above
   `0.02`.
6. Applies a final rank inverse-normal transformation.

Important limitation:

The paper mentions age adjustment, but age was not available in the
GeneNetwork phenotype table used here. Therefore this local pipeline is a
practical MAGMA/H-MAGMA summary-stat generation pipeline, not a perfect
byte-for-byte reproduction of the paper preprocessing.

## GRM

Union GRM output prefix:

`/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022/results/target_traits_union_grm`

GRM was generated with:

```bash
THREADS=10 /Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022/scripts/run_gcta_gwas.sh full-grm
```

Observed GRM details:

- Samples: `1,346`
- Autosomal SNPs included before MAF filtering: `7,182,588`
- Valid SNPs used for GRM: `5,412,081`
- Runtime: about `1 minute 19 seconds`

## Commands

Prepare inputs only:

```bash
/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022/scripts/run_gcta_gwas.sh prepare
```

Run a smoke test on trait `10424`:

```bash
THREADS=10 /Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022/scripts/run_gcta_gwas.sh smoke 10424
```

Run one full LOCO trait:

```bash
THREADS=10 /Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022/scripts/run_gcta_gwas.sh full-trait 10424
```

Run all traits, skipping completed MAGMA inputs:

```bash
THREADS=10 /Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022/scripts/run_gcta_gwas.sh full-all
```

Run all traits in a detached `screen` session:

```bash
screen -dmS gunturkun_gwas /bin/bash -lc "THREADS=10 /Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022/scripts/run_gcta_gwas.sh full-all > /Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022/logs/full_all_screen_YYYYMMDD_HHMMSS.log 2>&1"
```

Check running session:

```bash
screen -ls
ps aux | grep -E 'SCREEN|run_gcta_gwas|gcta64' | grep -v grep
```

Attach to the running session:

```bash
screen -r gunturkun_gwas
```

Detach without stopping the job:

```text
Ctrl-a, then d
```

Stop the detached batch:

```bash
screen -S gunturkun_gwas -X quit
```

## Output Files

For each trait, GCTA writes:

```text
results/trait_<trait>.loco.mlma
```

The MAGMA-ready p-value input is:

```text
results/trait_<trait>.magma.tsv
```

MAGMA input format:

```text
SNP    p
```

Rows with invalid, non-finite, zero, or `nan` p-values are filtered out during
conversion.

## Notes For MAGMA/H-MAGMA

The `.magma.tsv` files are p-value inputs. A typical next step is to combine
these with a rat SNP-to-gene annotation or H-MAGMA annotation compatible with
the SNP IDs in the genotype `.bim` file. The current SNP IDs are mostly
chromosome-position style IDs such as:

```text
1:22585
1:198606
```

