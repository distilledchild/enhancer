# h3_GWAS Project Files

This directory contains lightweight code and documentation for HS rat GWAS,
cMAGMA, and H-MAGMA analyses. Large generated GWAS/MAGMA outputs are stored on
the external SSD and ignored by git.

## Dataset Folders

| Folder | Role |
|---|---|
| `BMI/` | Wright BMI proof-of-principle GEMMA summary-stat and cMAGMA/H-MAGMA workflow. |
| `Cocaine2026/` | Public cocaine self-administration `.mlma` summary-stat conversion and cMAGMA/H-MAGMA analysis. |
| `Gunturkun2022/` | Practical local GCTA MLMA/MLMA-LOCO pipeline for GeneNetwork HS rat behavioral traits. |
| `King2025/` | Cue-reactivity/sign-tracking summary-stat validation and cMAGMA/H-MAGMA discovery analysis. |
| `Kuhn2025/` | Near-exact heroin/nociception GWAS reproduction/validation plus downstream cMAGMA/H-MAGMA. |
| `Lara2024_delaydiscounting/` | Deposited delay-discounting `.mlma` validation and selected cMAGMA/H-MAGMA analysis. |
| `nicotine/` | Historical nicotine SA MAGMA/H-MAGMA scripts and strict annotation files used by the project. |
| `docs/` | General GCTA GWAS pipeline and result interpretation guides. |

## Storage Policy

Project code, scripts, analysis helpers, and lightweight documentation live in
this local directory. Heavy input/output files stay on the external SSD. See
`EXTERNAL_STORAGE_POLICY.md` for the six-project layout and move rules.

## Git Policy

Track:

- `.md` documentation
- `.R`, `.py`, `.sh` scripts
- lightweight reproducibility notes
- reference pipeline source snapshots needed to understand reproduction commands

Ignore:

- `.mlma`, `.genes.out`, `.genes.raw`, `.log`
- MAGMA p-value inputs
- large CSV/TSV outputs
- downloaded data folders
- PDFs and compressed archives

The corresponding external work root is:

```text
/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS
```
