# Cocaine2026 Expanded Trait Screen

Date: 2026-09-06. Scope: public summary-statistic screening for a possible revised-loop H-MAGMA follow-up.

## Decision

**Zero of 29 archived traits passed the requested lead-SNP threshold of -log10(P) >= 7 (P <= 1e-7).** The strongest signal was `regressedlr_lga_mean_11_14`, SNP `6:34284564`, P = 5.62593e-7 (-log10(P) = 6.249806). No new H-MAGMA run was launched; the threshold was not relaxed. Existing four-trait MAGMA results were not modified.

This is a genome-wide input screen, not an analysis filtered by our loops. `2001bp` in the output-directory name identifies the intended follow-up context; the revised 13,376-loop annotation was not applied in this step.

## Actual Available Data

- The complete local archive contains 29 traits x 22 chromosome files = **638 MLMA files**, representing 27 dictionary-documented traits and two additional variables (`shock_03_pre`, `shock_03_avg1h`). There are 25 traits beyond the previously tested four, but these are correlated measurements from the same cohort, not independent datasets.
- Available behavioral categories are short access (ShA), long access (LgA), progressive ratio (PR), footshock-related behavior, and total intake. **No locomotor-sensitization, extinction, or reinstatement GWAS files were found in this archive.**
- `pr_01_active` measures total active lever presses in PR01; `pr_max` is the maximum total active lever presses across PR02/PR03. These are not reported breakpoint variables. Neither these measurements nor their GWAS associations establish exclusive dependence on the PFC.
- The source specifies mRatBN7.2 coordinates. Each trait has 5,446,333 SNP rows on chromosomes 1-20, X, and MT. Raw MLMA chromosome codes 21 and 24 denote X and MT, respectively; SNP identifiers were used for cross-checking.
- Trait-specific phenotype N ranges from 611 to 836. N agrees with the original report for **28/29 traits**. `shock_03_avg1h` has 625 nonmissing regressed phenotype values but no corresponding report N; it remains unverified for downstream MAGMA use.

## Key Results

| Trait | Definition | Verified N | Lead SNP | Minimum P | Peak -log10(P) |
|---|---|---:|---|---:|---:|
| `lga_mean_11_14` | Mean infusions, LgA days 11-14 | 810 | 6:34284564 | 5.62593e-7 | 6.2498 |
| `sha_mean_to_01_03` | Mean active presses minus infusions, ShA days 1-3 | 835 | 9:57520529 | 5.77650e-7 | 6.2383 |
| `lga_titration_intake_01_03` | Last-hour infusion sum, LgA days 1-3 | 783 | X:65656263 | 6.84495e-7 | 6.1646 |
| `pr_max` | Maximum PR02/PR03 total active presses | 805 | 4:114813329 | 2.10303e-5 | 4.6772 |
| `pr_01_active` | PR01 total active presses | 831 | 6:20688058 | 4.76762e-5 | 4.3217 |

All archive trait identifiers have the prefix `regressedlr_`. The complete ranked table includes all 29, their definitions, N provenance, previous-four status, and nuclear/all-chromosome peak values.

## Workflow And Verification

1. Inventory the ZIP, parse the local dictionary/phenotypes and embedded source-report tables, and record input hashes. No dataset was downloaded or GWAS regenerated.
2. Stream every chromosome file through awk, checking decompression status and repeated headers; count valid, missing, invalid, zero, and threshold-passing P values and record each chromosome's minimum.
3. Screen nuclear chromosomes including X, while retaining an all-chromosome result including MT. Both definitions yield zero threshold-7 candidates. All **157,943,657 trait-SNP rows** were read; no missing, invalid, or zero P values were detected. All 29 nuclear lead P values agree with beta/SE-derived Wald P values within 0.000042 log10 units.
4. Independently reparse all 66 chromosome files for the strongest trait and both PR traits with R. Minimum P values, row counts, threshold counts, SNP-ID uniqueness within each chromosome file, and ID-coordinate consistency all pass. Separately confirm all seven source-report lead SNP P values within the report's rounding precision. Synthetic tests check the exact threshold boundary, ties, repeated headers, and invalid/missing/zero P handling.

## Interpretation Limits

The paper's permutation-derived 5% threshold is -log10(P) = **5.58**, less stringent than the requested 7. Eight archived traits have at least one peak above 5.58; this is **not eight independent loci** and is not a correction for selecting among all traits. The paper reports six genome-wide significant associations. [Published study](https://www.nature.com/articles/s41467-026-73694-w).

Failure to pass 7 is not a formal power calculation, evidence that all associations are null, or proof that MAGMA gene-level aggregation cannot yield a significant result. Selection on observed GWAS peaks would be post-hoc exploration, not independent validation of our regulatory loops. No downstream LD/reference harmonization or new MAGMA inference was performed because no trait met the requested gate.

## Files And Reproduction

Input root: `/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Cocaine2026/data/`.
- `bb2334903t_4_1.zip`: chromosome-split GWAS summary statistics.
- `bb2334903t_1_1.csv`: phenotype dictionary.
- `bb2334903t_3_1.csv`: phenotype table, 836 unique animals.
- `bb2334903t_2_1.html`: original report generated 2025-07-23, including N, thresholds, build, and lead SNPs.
- Source: [UCSD dataset, DOI 10.6075/J0QN675X](https://doi.org/10.6075/J0QN675X).

Outputs:
- [All 29 traits, ranked TSV](/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Cocaine2026/extended_trait_screen_2001bp/all_29_trait_peak_screen.tsv)
- [Threshold-7 candidate TSV, header only](/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Cocaine2026/extended_trait_screen_2001bp/threshold7_candidates.tsv)
- [Validation summary](/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Cocaine2026/extended_trait_screen_2001bp/validation_summary.tsv)

The same output directory contains per-chromosome summaries, independent R checks, source-report concordance, trait metadata, archive inventory, input/script MD5 manifests, and R session information.

Run in order:
```bash
Rscript /Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/Cocaine2026/scripts/screen_extended_traits.R
Rscript /Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/Cocaine2026/scripts/validate_extended_trait_screen.R
```
The first script calls `/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/Cocaine2026/scripts/scan_mlma_peaks.awk`. Dependencies: R packages data.table, xml2, jsonlite; system unzip, awk, and bash. Re-running overwrites only this screen's outputs, not previous MAGMA analyses.
