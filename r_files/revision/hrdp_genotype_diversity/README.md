# HRDP common-variant coverage analysis

`hrdp_variant_coverage.R` and `hrdp_variant_coverage.awk` analyze the unpruned
autosomal gVCFs in `228_gvcf`. The denominator is fixed to 94 HRDP strains:
30 classic inbred strains, 30 HXB/BXH recombinant-inbred strains, and 34
FXLE/LEXF recombinant-inbred strains. One named gVCF specimen is selected for
each strain; the remaining specimen columns are not included in this
denominator.

The ten Hi-C libraries represent nine unique strain backgrounds plus one SHR x
BN F1. Variant coverage is calculated from the nine parental/inbred
backgrounds because the F1 cannot introduce an allele absent from both of its
parents. Autosomal records with `FILTER` equal to `.` or `PASS` are analyzed.
`QUAL >= 30`, `MAF > 0.10` across the 94 strains, and site-level coverage are
primary; `MAF > 0.20`, `QUAL >= 40`, no-MAF-filter, and allele-level coverage
are comparison or sensitivity results. SNPs, indels, and other variant types
are reported separately.

At Q30 and MAF > 0.10, the 94-strain panel contains 7,566,911 common SNP sites
and 1,952,961 common indel sites. The nine study backgrounds represent
6,890,890 SNP sites (91.1%) and 1,752,026 indel sites (89.7%). At MAF > 0.20,
the corresponding estimates are 93.7% and 92.8%; without the MAF filter, they
are 62.6% and 66.9%. These values describe common-variant representation; they
do not test genotype-loop associations or strain-specific loop biology.

Outputs are written to `results/variant_coverage`:

- `hrdp_panel_sample_mapping.tsv`
- `hrdp_variant_coverage_summary.tsv`
- `hrdp_variant_coverage_by_chromosome.tsv`
- `hrdp_variant_coverage_run_metadata.tsv`
- `hrdp_variant_coverage.pdf` and `.png`

Set `HRDP_GVCF_DIR` only when the `228_gvcf` directory is stored elsewhere.

Run the R file from top to bottom, or execute:

```bash
Rscript hrdp_variant_coverage.R
```
