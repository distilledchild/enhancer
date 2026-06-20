# Exact Reproduction Audit: Gunturkun et al. 2022

This note separates the current exploratory MAGMA-summary pipeline from the
requirements for an exact reproduction of Gunturkun et al. 2022.

## Bottom Line

The current local run is not an exact reproduction. It is an exploratory local
GWAS summary-stat pipeline.

The most important fixes before rerunning are:

1. Use the paper-era 3.5M GBS genotype set, not the current 7.3M v4 set.
2. Restrict to the paper's 1,246 adolescent HS rats.
3. Recover the exact age covariate and paper sample list from C-GORD/RatGenes
   or author data.
4. Reimplement trait-level covariate selection using significance and >2%
   variance explained.
5. Validate against paper invariants before running all traits.

## Paper Invariants To Match

These are the non-negotiable checks.

- Study cohort:
  `626 males + 620 females = 1,246 rats`
- Tests:
  `OFT`, `NOIT`, `SIT`, one per day in that order
- Trait count:
  `23`
- Phenotype preprocessing:
  - quantile-normalize each trait separately within males and females
  - identify age, batch number, and coat color covariates per trait
  - regress covariates out only if significant and explaining >2% variance
  - quantile-normalize residuals again
  - pool males and females
- Genotype data:
  - GBS
  - approximately `3.5 million` SNPs
  - X and Y variants not called
- GWAS:
  - GCTA linear mixed model
  - GRM for relatedness
  - LOCO
- Threshold:
  `-log10(P) > 5.609` / about `5.6`
- QTL calling:
  - at least one SNP beyond permutation threshold
  - support SNP within 0.5 Mb and within 2 log10 units of top SNP
  - conditional chromosome-specific GWAS to split multiple loci
  - LD interval with r2 = 0.6

## Current Local Pipeline Problems

Current exploratory pipeline root:

`/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022`

Current exploratory genotype prefix:

`/Volumes/external_1000GB_all/playground/enhancer/data/hs_data/v4/HS_genotypes_v4`

Observed current genotype:

- `.bim` SNPs: `7,358,643`
- autosomal SNPs used by GCTA: `7,182,588`
- `.fam` samples: `17,812`

This does not match the paper's approximately 3.5M SNP GBS set.

Current prepared phenotype sample counts:

- NOIT union: `1,237`
- OFT union: `1,238`
- SIT union: `1,320`
- all trait union: `1,346`
- all trait intersection: `1,121`

These do not match the paper's `1,246` rats. In particular, SIT N = `1,320`
is impossible under the paper cohort, so the GeneNetwork `HSNIH-PalmerPublish`
matrix includes rows beyond, or not identical to, the exact paper analysis set.

Current completed trait:

`trait_10424`

Current top result:

`11:32550603`, `p = 1.7827e-07`, `-log10(p) ~= 6.75`

Paper Table 2 for OFT total travel distance includes:

- `chr10:94549701`, `-log10P = 7.286`
- `chr11:33359859`, `-log10P = 8.268`

The paper top SNP `chr11:33359859` is not present in the current v4 result
file, and the nearby current signal around chr11:33.36Mb is weaker. This is
consistent with a marker-set/build mismatch.

## Local Data That Looks Closer To The Paper

Likely paper-era genotype candidate:

`/Volumes/external_1000GB_all/playground/enhancer/data/King_2025/bb15123938_2_1/round8_unpruned`

Observed:

- `.bim` SNPs: `3,513,494`
- `.fam` samples: `6,147`
- chromosomes: `1-20`
- file timestamp: `2022-04-18`

This matches the paper statement of approximately 3.5M SNPs much better than
the v4 genotype.

Additional check completed:

- All 24 distinct Table 2 top SNP coordinates tested were present exactly in
  `round8_unpruned.bim`.
- Examples:
  - `chr10:94549701`
  - `chr11:33359859`
  - `chr4:112234344`
  - `chr19:55339863`

This makes `round8_unpruned` the best local candidate genotype for an exact
paper reproduction.

Do not assume it is exact until these checks pass:

1. Confirm paper sample IDs can be selected from its `.fam`.
2. Confirm exact allele coding/frequency direction if beta effects are compared.
3. Confirm the original authors used this exact unpruned file, or an equivalent
   filtered derivative.

## Required Verification Steps Before Rerunning

### 1. Resolve The Dataset DOI

The paper text says:

`10.48810/P44W2`

RatGenes lists:

`10.48810/P44W2Q`

Action:

- Try both DOI strings.
- If neither resolves to a downloadable C-GORD object, contact:
  `GWAS@ratgenes.org`
- Ask specifically for:
  - exact analyzed sample list of 1,246 rats
  - raw phenotype table for the 23 traits
  - age, sex, batch number, coat color covariates
  - exact genotype PLINK/VCF prefix or SNP inclusion list
  - exact GCTA command/script if available
  - permutation procedure or seed/count

### 2. Verify Trait ID Mapping

Known GeneNetwork old trait page example:

`https://gn1.genenetwork.org/webqtl/main.py?FormID=showDatabase&RISet=HSNIH-Palmer&database=HSNIH-PalmerPublish&ProbeSetID=10418`

Observed 10418 description:

`Distance traveled (locomotion) in an open field (1 x 1m) that includes a novel object in arena center ... [cm/hr]`

Action:

- Scrape/record descriptions for `10418-10440`.
- Map each GeneNetwork trait ID to the exact paper trait name:
  - OFT: six traits
  - NOIT: six traits
  - SIT: eleven traits
- Save mapping in `data/repro_exact/trait_map.tsv`.

### 3. Recover The Exact 1,246 Sample List

Required output:

`data/repro_exact/paper_1246.keep`

Validation:

- line count must be `1246`
- sex count must be `626 male`, `620 female`
- all selected IDs must exist in the chosen genotype `.fam`
- selected IDs must have expected raw trait/covariate coverage

If the exact paper list is unavailable, do not call the run exact. Label it
`paper-like reproduction`.

### 4. Build Exact Phenotype Processor

The current script is approximate:

`scripts/prepare_gunturkun2022_traits.py`

For exact reproduction, write a separate script:

`scripts/prepare_gunturkun2022_traits_exact.py`

Required behavior:

1. Load only `paper_1246.keep`.
2. For each trait:
   - split by sex
   - quantile-normalize within sex
   - test age, batch number, and coat color covariates
   - retain covariates only if significant and explaining >2% variance
   - regress retained covariates
   - quantile-normalize residuals
   - pool sexes
3. Save:
   - per-trait phenotype files
   - retained covariates
   - covariate p-values
   - covariate R2
   - N per trait

Important uncertainty:

The paper does not specify the exact significance test or alpha threshold for
covariate selection. This must be obtained from author code or inferred and
documented.

### 5. Use Paper-Era Genotypes

Candidate prefix:

`/Volumes/external_1000GB_all/playground/enhancer/data/King_2025/bb15123938_2_1/round8_unpruned`

Required checks:

```bash
wc -l round8_unpruned.bim
wc -l round8_unpruned.fam
awk '{print $1}' round8_unpruned.bim | sort -V | uniq -c
```

Expected:

- about `3,513,494` SNPs
- chromosomes `1-20`

Check whether Table 2 SNPs exist:

```bash
awk '$2=="chr11:33359859" || $2=="11:33359859" {print}' round8_unpruned.bim
```

Repeat for all Table 2 top SNPs.

### 6. Match GCTA Output Against Paper Before Full Run

Run only OFT total travel distance first. This is the strongest sanity check
because Table 2 gives two clear loci.

Expected top signals:

- `chr10:94549701`, `-log10P = 7.286`
- `chr11:33359859`, `-log10P = 8.268`

Pass criteria:

- same top SNP IDs, or exact proxy SNPs if marker IDs differ but coordinates
  align
- -log10(P) within small numerical tolerance
- same genome-wide significant loci
- threshold `5.609`

If this fails, do not run all 23 traits.

### 7. Reproduce Downstream QTL Calls

After GWAS p-values match:

1. Implement threshold filter.
2. Implement 0.5Mb support-SNP rule.
3. Implement conditional chromosome reruns with top SNPs as covariates.
4. Compute LD intervals at r2 = 0.6.
5. Compare Table 2 locus list.

Only after Table 2 is reproduced should MAGMA/H-MAGMA summary stats be treated
as paper-reproduced.

## Evidence Links

- Frontiers article:
  `https://www.frontiersin.org/journals/psychiatry/articles/10.3389/fpsyt.2022.790566/full`
- RatGenes data page:
  `https://ratgenes.org/genome-wide-association-study-on-three-behaviors-tested-in-an-open-field-in-heterogeneous-stock-rats-identifies-multiple-loci-implicated-in-psychiatric-disorders/`
- GeneNetwork 10418 page:
  `https://gn1.genenetwork.org/webqtl/main.py?FormID=showDatabase&RISet=HSNIH-Palmer&database=HSNIH-PalmerPublish&ProbeSetID=10418`

## 2026-06-18 Data Scrape Update

Detailed scrape manifest:

`/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022/data/repro_exact/SCRAPE_MANIFEST.md`

Key new evidence:

- Frontiers XML/PDF/HTML were saved locally and Tables 1-3 were extracted to TSV.
- RatGenes official pages point to dataset DOI `10.48810/P44W2Q`; paper text points to `10.48810/P44W2`. Both DOI resolver checks currently fail and are saved under `data/repro_exact/doi/`.
- GeneNetwork metadata for all 508 `HSNIH-PalmerPublish` traits was scraped successfully.
- Paper trait IDs are mapped as NOIT `10418-10423`, OFT `10424-10429`, SIT `10430-10440`, covariates `10443-10446`.
- `10424 open_field_totaldistance` has exactly `1,246` non-missing samples, with `620` female and `626` male, matching the paper. This sample list is saved at:
  `data/repro_exact/genenetwork/cohort_oft_totaldistance_10424_n1246.samples.txt`
- `round8_unpruned` has `3,513,494` SNPs, chromosomes `1-20`, and all `30/30` Table 2 peak SNP coordinates are present.
- The reconstructed 1,246-sample cohort is fully present in `round8_unpruned.fam`.

Remaining blocker:

The exact per-sample age covariate and exact covariate-selection implementation were not recovered from public GeneNetwork/Frontiers/RatGenes pages. These likely require C-GORD access or author code/contact.
