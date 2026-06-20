# h3_GWAS Storage Policy

This directory is the canonical location for GWAS project code, scripts,
analysis helpers, and lightweight documentation.

Large inputs, downloaded data, generated GWAS outputs, MAGMA outputs, logs,
figures, archives, and other heavy artifacts stay on the external SSD.

## Canonical Roots

Local code/documentation root:

```text
/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS
```

External data/result root:

```text
/Volumes/external_1000GB_all/playground/enhancer
```

## Project Layout

| Project | Local code/docs | External heavy data/results |
|---|---|---|
| BMI | `r_files/h3_GWAS/BMI/` | `data/BMI/`, `data/GWAS_3173_Outbred_Rats_BodyWeight_Adiposity_and_Fasting_Glucose/`, `data/obesity/`, `r_files/h3_GWAS/BMI/` |
| Cocaine2026 | `r_files/h3_GWAS/Cocaine2026/` | `r_files/h3_GWAS/Cocaine2026/` |
| Gunturkun2022 | `r_files/h3_GWAS/Gunturkun2022/` | `r_files/h3_GWAS/Gunturkun2022/` |
| King2025 | `r_files/h3_GWAS/King2025/` | `data/King_2025/`, `r_files/h3_GWAS/King2025/` |
| Kuhn2025 | `r_files/h3_GWAS/Kuhn2025/` | `data/Khun_2025/`, `r_files/h3_GWAS/Kuhn2025/` |
| Lara2024_delaydiscounting | `r_files/h3_GWAS/Lara2024_delaydiscounting/` | `r_files/h3_GWAS/Lara2024_delaydiscounting/` |

## What Belongs In Git

- `.sh`, `.py`, `.R`, `.Rmd`, `.ipynb`
- `.md`
- lightweight `.yml` / `.yaml` environment or configuration files
- small text manifests needed to understand how the pipeline was run

## What Belongs On The External SSD

- raw genotype, phenotype, summary-stat, annotation, archive, and PDF inputs
- generated `.mlma`, MAGMA p-value inputs, `.genes.out`, `.genes.raw`, and logs
- generated result `.csv` / `.tsv` files
- figures, downloaded reports, and large recovery folders
- tool binaries and large third-party databases

## Current Move Check

As of the cleanup commit that introduced this policy, the external
`data/` and `r_files/h3_GWAS/` trees contain no remaining project code or
documentation files with these extensions:

```text
.sh .py .R .Rmd .md .yml .yaml .ipynb
```

External `.git` metadata under GWAS project folders was also removed to avoid
confusion about which location is canonical.
