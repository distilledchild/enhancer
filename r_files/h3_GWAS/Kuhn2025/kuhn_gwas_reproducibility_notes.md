# Kuhn et al. 2025 GWAS Reproducibility Notes

Status: input and method reconstruction only. GWAS was not executed.

## Local data

Base directory:

`/Volumes/external_1000GB_all/playground/enhancer/data/Khun_2025`

Genotype prefix prepared for PLINK/GCTA:

`/Volumes/external_1000GB_all/playground/enhancer/data/Khun_2025/hs_rats_round10_2`

This prefix is a symlink set to the downloaded PLINK files:

- `hs_rats_round10_2.bed -> bb29129987_2_1.bed`
- `hs_rats_round10_2.bim -> bb29129987_3_1.bim`
- `hs_rats_round10_2.fam -> bb29129987_4_1.fam`

Current checks:

- `.bim`: 7,358,643 variants
- `.fam`: 17,812 rats
- `processed_data_ready.csv`: 874 rats, 241 columns
- `data_dict_u01_peter_kalivas.csv`: 132 rows
- All 874 phenotype RFIDs match IIDs in the `.fam` file.

Run-ready project directory:

`/Volumes/external_1000GB_all/playground/enhancer/data/Khun_2025/gwas_inputs/u01_peter_kalivas`

This directory contains symlinks to:

- `processed_data_ready.csv`
- `data_dict_u01_peter_kalivas.csv`
- `data_distributions.html`
- `locuszoom/`

and the required project/output folders expected by the pipeline:

- `data/pheno/`
- `genotypes/`
- `grm/`
- `results/gwas/`
- `results/heritability/`
- `results/qtls/`

As of this check, `data/pheno/` has also been pre-populated from `processed_data_ready.csv` with 96 `regressedlr_*` phenotype files, each with 874 rows and three columns:

`FID IID PHENO`

The sample keep lists were also generated:

- `genotypes/keep_rfids.txt`: 874 rats
- `genotypes/keep_rfids_males.txt`
- `genotypes/keep_rfids_females.txt`

## Public supplement contents

The downloaded Nature/Molecular Psychiatry supplement archive

`supplementary/41380_2025_2922_MOESM3_ESM_GWAS_dataset.zip`

contains only:

- `gwas_report_u01_peter_kalivas_round10.2.1_threshold5.3591_n874_date2024-01-26_gwasversion_v0.2.0-11-g990844d.html`

The full SNP-level GWAS summary statistics are not present in that zip. The report points to the phenotype and data dictionary files, which were downloaded separately into `gwas_inputs/`.

## Pipeline code

Exact Python pipeline commit matching the report filename:

`/Volumes/external_1000GB_all/playground/enhancer/data/Khun_2025/gwas_inputs/reference_pipeline/sanchestm_GWAS-pipeline_v0.2.0-11-g990844d`

Git commit:

`990844d79a7e69a6619d5c394d98b06f81ab2cb2`

Git describe:

`v0.2.0-11-g990844d`

This resolves the first open issue: the actual GWAS pipeline version referenced by the Kuhn GWAS report is available locally.

## GCTA/PLINK version evidence

The exact report-linked pipeline commit does not pin executable versions. Its `environment.yml` lists:

- `plink`
- `gcta`

without version constraints, and `gwas_class_auto.py` uses:

- `plink` on `PATH`
- `gcta64` on `PATH`

The report HTML also does not include `plink --version` or `gcta64 --version` output.

Best current reconstruction:

- GCTA was probably the Bioconda `gcta` package available before the report date, i.e. `1.94.1`.
- PLINK was probably the Bioconda `plink` package available before the report date, i.e. `1.90b6.21`.

This is a strong environment-based inference, not an exact recorded fact, because no conda lock file or runtime version log is present in the report-linked commit.

The later/current GitHub README says the pipeline is centered around PLINK data from "PLINK 1.96", but that wording is not present in the exact Kuhn report-linked commit and should not override the report-linked environment evidence.

## Phenotype/covariate handling

This resolves the second open issue.

The distributed `processed_data_ready.csv` already contains raw traits and preprocessed `regressedlr_*` traits. In normal no-`regressout` mode, `gwas_cli.py` reads:

`processed_data_ready.csv`

and selects all columns beginning with:

`regressedlr_`

The `gwas_pipe.__init__()` method then writes one GCTA phenotype file per trait:

`data/pheno/<regressedlr_trait>.txt`

with three columns:

`FID IID PHENO`

where both FID and IID are the rat `rfid`.

Therefore, for the downloaded Kuhn data, the GWAS phenotype passed to GCTA is the already processed `regressedlr_*` column. Raw phenotype plus covariates are used only if the CLI is run with the `regressout` option, which regenerates `processed_data_ready.csv` and the `data/pheno/*.txt` files.

Main traits currently present and complete:

- `regressedlr_total_heroin_consumption`: 874 non-missing
- `regressedlr_bp`: 874 non-missing
- `regressedlr_escalation_of_heroin_intake_12h`: 874 non-missing
- `regressedlr_tail_flick_bl1`: 874 non-missing

## Genotype filtering implemented by the exact pipeline

The exact pipeline's `SubsetAndFilter()` method:

1. Writes the 874 matched RFIDs to `genotypes/keep_rfids.txt`.
2. Computes PLINK missingness, MAF, and HWE stats in the 874-rat subset.
3. Keeps SNPs passing:
   - missingness: `F_MISS < 0.1`
   - MAF: `MAF >= 0.005`
   - HWE: `P > 1e-10`, with chromosome-specific exceptions in code
4. Writes accepted SNP IDs to `genotypes/accepted_snps.txt`.
5. Creates the filtered genotype prefix:

`u01_peter_kalivas/genotypes/genotypes`

using PLINK `--extract`, `--keep`, `--make-bed`, `--keep-allele-order`, `--set-hh-missing`, and `--make-founders`.

The report states 7,358,643 SNPs before filtering and 5,424,892 after filtering, matching the local `.bim` before-filter variant count.

## GRM and GWAS model

The exact pipeline's `generateGRM()` method:

1. Creates per-autosome GRMs for chromosomes 1-20:

`grm/<chr>chrGRM`

2. Writes:

`grm/listofchrgrms.txt`

3. Merges the autosomal chromosome GRMs into:

`grm/AllchrGRM`

The exact pipeline's `fastGWAS()` method then runs per-trait, per-chromosome GCTA MLMA commands of this form:

```bash
gcta64 --thread-num 1 \
  --pheno /path/to/u01_peter_kalivas/data/pheno/<trait>.txt \
  --bfile /path/to/u01_peter_kalivas/genotypes/genotypes \
  --grm /path/to/u01_peter_kalivas/grm/AllchrGRM \
  --autosome-num 20 \
  --chr <chromosome> \
  --mlma-subtract-grm /path/to/u01_peter_kalivas/grm/<chromosome>chrGRM \
  --mlma \
  --out /path/to/u01_peter_kalivas/results/gwas/<trait>_chrgwas<chromosome>
```

For X and Y, the code omits `--mlma-subtract-grm`.

## Dry-run command template

Do not run this until GCTA/PLINK availability and Python environment are checked.

```bash
cd /Volumes/external_1000GB_all/playground/enhancer/data/Khun_2025/gwas_inputs/reference_pipeline/sanchestm_GWAS-pipeline_v0.2.0-11-g990844d

python3 gwas_cli.py \
  path=/Volumes/external_1000GB_all/playground/enhancer/data/Khun_2025/gwas_inputs \
  project=u01_peter_kalivas \
  genotypes=/Volumes/external_1000GB_all/playground/enhancer/data/Khun_2025/hs_rats_round10_2 \
  n_autosome=20 \
  genome=rn7 \
  round=10.2.1 \
  gwas_version=v0.2.0-11-g990844d \
  threshold=5.3591 \
  subset=1 \
  grm=1 \
  gwas=1
```

Optional downstream report/QTL steps in the same CLI:

- `h2=1`
- `qtl=1`
- `manhattanplot=1`
- `locuszoom=1`
- `report=1`

## Remaining caveats before execution

1. The paper text and report/genotype files appear to differ in total/final SNP counts. The local genotype and report agree on 7,358,643 before filtering; the paper text may refer to a different genotype release or reporting layer.
2. The public supplement does not include full SNP-level summary statistics. If exact published summary stats are needed for direct H-MAGMA/MAGMA input without rerunning GWAS, they still need to be located elsewhere or requested.
3. Reproducing the report should use the exact checked-out pipeline commit above, not the current HEAD of the GitHub repository.
4. The dry-run command should be executed only after confirming the Python package environment plus `gcta64` and `plink` behavior.
