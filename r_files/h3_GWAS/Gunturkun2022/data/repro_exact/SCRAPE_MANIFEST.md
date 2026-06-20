# Gunturkun 2022 Exact-Reproduction Data Scrape Manifest

Date: 2026-06-18

Root:

`/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022/data/repro_exact`

## What Was Collected

### Frontiers Article

Saved official article sources:

- `frontiers/frontiers_full.html`
- `frontiers/frontiers_article.xml`
- `frontiers/frontiers_article.pdf`

Extracted from XML:

- `frontiers/T1.tsv`: heritability table
- `frontiers/T2.tsv`: Table 2 QTL table
- `frontiers/T3.tsv`: candidate genes table
- `frontiers/methods_key_points.txt`: methods sections needed for GWAS reproduction
- `frontiers/supplementary_materials.tsv`: XML supplementary-material declaration

Important result:

- Table 2 extraction produced 30 QTL rows.
- XML declares one supplement: `Data_Sheet_1.pdf`.
- The direct supplement file was not recovered. Probe evidence is in:
  - `frontiers/frontiers_supplement_probe_manifest.tsv`
  - `frontiers/frontiers_supplement_pattern_probe.tsv`

### RatGenes / C-GORD / DOI Evidence

Saved official pages:

- `ratgenes/ratgenes_gunturkun_page.html`
- `ratgenes/ratgenes_related_dois.html`
- `ratgenes/ratgenes_fair_statement.html`
- `ratgenes/scicrunch_SCR_021866.html`

DOI checks saved:

- `doi/10.48810_P44W2.summary.txt`
- `doi/10.48810_P44W2Q.summary.txt`

Important result:

- RatGenes official data page and related DOI page point to `doi.org/10.48810/P44W2Q`.
- The paper XML/PDF text points to `10.48810/P44W2`.
- Both DOI resolver checks currently return DOI-not-found style pages, so C-GORD DOI access is blocked from this machine as of 2026-06-18.

### GeneNetwork

Saved:

- `genenetwork/HSNIH-PalmerPublish.csv`: full public phenotype matrix
- `genenetwork/HSNIH-PalmerPublish_trait_metadata.tsv`: metadata for all 508 HSR traits
- `genenetwork/sample_data_traits/`: sample-level JSON for 220 selected behavior/covariate candidate traits
- `genenetwork/paper_trait_metadata_10418_10446.tsv`: paper behavior traits and covariates
- `genenetwork/behavior_trait_nonmissing_counts.tsv`: per-trait non-missing N and sex/batch/center counts
- `genenetwork/paper_cohort_candidate_sets.tsv`: candidate cohort set sizes
- `genenetwork/cohort_oft_totaldistance_10424_n1246.samples.txt`: strongest reconstructed paper cohort candidate
- `genenetwork/cohort_oft_totaldistance_10424_n1246.tsv`: reconstructed cohort with 23 behavior traits plus covariates

Important result:

- Paper trait IDs are strongly mapped as:
  - NOIT: `10418-10423`
  - OFT: `10424-10429`
  - SIT: `10430-10440`
  - covariates: `10443 sex`, `10444 batchnumber`, `10445 color`, `10446 center`
- `open_field_totaldistance` / trait `10424` has exactly `1,246` non-missing samples.
- That trait's cohort has exactly `620 female` and `626 male`, matching the paper.
- All 1,246 reconstructed cohort samples have `center=3.0` (University of Tennessee).

### Local Genotype Audit

Audited prefix:

`/Volumes/external_1000GB_all/playground/enhancer/data/King_2025/bb15123938_2_1/round8_unpruned`

Saved:

- `local_audit/round8_file_manifest.tsv`: file sizes, mtimes, SHA256 hashes
- `local_audit/round8_chr_distribution.tsv`: chromosome SNP counts
- `local_audit/table2_snps_in_round8.tsv`: Table 2 peak SNP presence check
- `local_audit/round8_cohort_overlap_summary.tsv`: reconstructed cohort vs `.fam`

Important result:

- `.bim` SNP count: `3,513,494`
- `.fam` sample count: `6,147`
- chromosomes present: `1-20`
- reconstructed 1,246-sample cohort is present in `.fam`: `1,246 / 1,246`
- Table 2 top SNPs present in `.bim`: `30 / 30`

Key SHA256 values:

- `round8_unpruned.bed`: `e6437b05f71d06be3efbf7ec49859879d1a38516e2e833996f057416f6a1d12c`
- `round8_unpruned.bim`: `843bcb77b4537f3d16c39799f383121cd41513d3a283a1ddad149a9ba40cf327`
- `round8_unpruned.fam`: `1ebe31637d2bd4c83bef2d0e511d69e9e262416be7a3f6403d58e2530fea126c`
- `round8_unpruned.sample`: `953b128d8d2dca301a988540af7addd4e65001b3462df8f5a6b9ccecac87584f`

## Reproduction Implications

Use this genotype for exact-reproduction attempts:

`/Volumes/external_1000GB_all/playground/enhancer/data/King_2025/bb15123938_2_1/round8_unpruned`

Use this sample list first:

`genenetwork/cohort_oft_totaldistance_10424_n1246.samples.txt`

Use this target sanity check first:

- Trait: `10424 open_field_totaldistance`
- Expected Table 2 loci:
  - `chr10:94549701`, `-log10P = 7.286`
  - `chr11:33359859`, `-log10P = 8.268`

## Remaining Gaps

1. The exact C-GORD dataset DOI is referenced but not currently resolvable.
2. `Data_Sheet_1.pdf` is declared in Frontiers XML but was not recoverable through direct/current Frontiers file endpoints.
3. The exact `age` covariate column was not found in GeneNetwork metadata. The public trait descriptions encode nominal test age, but not per-sample exact age.
4. The paper's covariate-selection implementation is still not public in the collected files: exact test, alpha threshold, R2 cutoff implementation, and retained covariates per trait remain to be reconstructed or obtained from authors.
5. GeneNetwork SIT traits have 1,328 non-missing samples, while the paper reports 1,246 total animals. The OFT trait `10424` exactly matches paper N and sex counts, so it is the strongest cohort anchor, but this should be validated against C-GORD or author code.

## Re-run Scripts

Scripts are in `scripts/`:

- `extract_frontiers_xml.py`
- `scrape_genenetwork.py`
- `analyze_gn_paper_cohort.py`
- `audit_round8_genotype.py`
- `probe_frontiers_supplement.py`
