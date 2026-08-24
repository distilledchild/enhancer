# lintr: disable
library("data.table")
library("tidyverse")

options(tibble.width = Inf)
options(tibble.print_max = Inf)
options(scipen = 999)

########################
# 1. Input files and output directory
########################

script.argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script.dir <- if (length(script.argument)) {
  dirname(normalizePath(sub("^--file=", "", script.argument[1])))
} else {
  normalizePath(getwd())
}

gvcf.candidates <- c(
  Sys.getenv("HRDP_GVCF_DIR", unset = ""),
  "/Volumes/external_1000GB_all/playground/enhancer/data/228_gvcf"
)
gvcf.candidates <- gvcf.candidates[dir.exists(gvcf.candidates)]
if (!length(gvcf.candidates)) {
  stop("Set HRDP_GVCF_DIR to the directory containing the 228-specimen gVCFs.")
}

gvcf.dir <- gvcf.candidates[1]
awk.file <- file.path(script.dir, "hrdp_variant_coverage.awk")
output.dir <- file.path(script.dir, "results/variant_coverage")
chromosomes <- paste0("chr", 1:20)
gvcf.files <- setNames(
  file.path(
    gvcf.dir,
    paste0(
      "deepvariant140_228_rats_", chromosomes,
      "_Qual30_hiHet_removed.gvcf.gz"
    )
  ),
  chromosomes
)

if (!file.exists(awk.file)) stop("Missing AWK helper: ", awk.file)
if (!all(file.exists(gvcf.files))) stop("At least one autosomal gVCF is missing.")
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

########################
# 2. One gVCF specimen for each of 94 official HRDP strains
########################

# The supplied 228-specimen gVCFs include repeated specimens and related
# substrains. This table fixes one representative specimen for each of the
# 30 classic, 30 HXB/BXH RI, and 34 FXLE/LEXF RI strains used as the denominator.
df.hrdp.map <- tribble(
  ~panel_group, ~official_strain, ~gvcf_sample,
  "classic", "ACI/EurMcwi", "ACI_EurMcwi_200534",
  "classic", "BDIX/NemOdaMcwi", "BDIX_NemOda_440",
  "classic", "BN/NHsdMcwi", "BN_NHsdMcwi_199632",
  "classic", "BN-Lx/CubMcwi", "BN_Lx_Cub_200538",
  "classic", "BUF/MnaMcwi", "BUF_Mna_441",
  "classic", "DA/OlaHsd", "DA_OlaHsd_201571",
  "classic", "F344/DuCrl", "F344_DuCrl_201577",
  "classic", "F344/NCrl", "F344_NCrl_199613",
  "classic", "F344/StmMcwi", "F344-stm",
  "classic", "FHH/EurMcwi", "FHH_EurMcwi_198627",
  "classic", "GK/FarMcwi", "GK_FarMcwi_199107",
  "classic", "LE/StmMcwi", "LE-stm",
  "classic", "LEW/Crl", "LEW_Crl_201570",
  "classic", "LH/MavRrrcAek", "LH_MavRrrc_593_GSPMC",
  "classic", "LL/MavRrrcAek", "LL_MavRrrc_594_GSPMC",
  "classic", "LN/MavRrrcAek", "LN_MavRrrc_595_GSPMC",
  "classic", "M520/NRrrcMcwi", "M520_NRrrcMcwi_NOV",
  "classic", "MNS/Gib", "ERR224459_MNS_Gib_TA.ILM",
  "classic", "MR/NRrrc", "MR_N_439",
  "classic", "MWF/Hsd", "MWF_Hsd_198844",
  "classic", "PVG/SeacMcwi", "PVG_Seac_1",
  "classic", "RCS/LavRrrcMcwi", "RCS_LavRrrc_673",
  "classic", "SBH/Ygl", "ERR224460_SBH_Ygl_TA.ILM",
  "classic", "SBN/Ygl", "ERR224461_SBN_Ygl_TA.ILM",
  "classic", "SHR/OlaIpcvMcwi", "SHR_Olalpcv_199265",
  "classic", "SHRSP/A3NCrl", "SHRSP_A3NCrl_201576",
  "classic", "SR/JrHsd", "SR_JrHsd_605_GSPMC",
  "classic", "SS/JrHsdMcwi", "SS_JrHsdMcwi_199478",
  "classic", "WAG/RijCrl", "WAG_RijCrl_604_GSPMC",
  "classic", "WKY/NCrl", "WKY_NCrl_201573",
  "HXB_BXH_RI", "BXH10", "BXH10_mRatNor1",
  "HXB_BXH_RI", "BXH11", "BXH11_mRatNor1",
  "HXB_BXH_RI", "BXH12", "BXH12_mRatNor1",
  "HXB_BXH_RI", "BXH13", "BXH13_mRatNor1",
  "HXB_BXH_RI", "BXH2", "BXH2_606_GSPMC",
  "HXB_BXH_RI", "BXH3", "BXH3_607_GSPMC",
  "HXB_BXH_RI", "BXH5", "BXH5_mRatNor1",
  "HXB_BXH_RI", "BXH6", "BXH6_Cub_678",
  "HXB_BXH_RI", "BXH8", "BXH8_mRatNor1",
  "HXB_BXH_RI", "BXH9", "BXH9_mRatNor1",
  "HXB_BXH_RI", "HXB1", "HXB1_mRatNor1",
  "HXB_BXH_RI", "HXB10", "HXB10_200103",
  "HXB_BXH_RI", "HXB13", "HXB13_mRatNor1",
  "HXB_BXH_RI", "HXB15", "HXB15_mRatNor1",
  "HXB_BXH_RI", "HXB17", "HXB17_Ipcv_674",
  "HXB_BXH_RI", "HXB18", "HXB18_Ipcv_675",
  "HXB_BXH_RI", "HXB2", "HXB2_198893",
  "HXB_BXH_RI", "HXB20", "HXB20_609_GSPMC",
  "HXB_BXH_RI", "HXB21", "HXB21_mRatNor1",
  "HXB_BXH_RI", "HXB22", "HXB22_mRatNor1",
  "HXB_BXH_RI", "HXB23", "HXB23_Ipcv_672",
  "HXB_BXH_RI", "HXB24", "HXB24_mRatNor1",
  "HXB_BXH_RI", "HXB25", "HXB25_mRatNor1",
  "HXB_BXH_RI", "HXB27", "HXB27",
  "HXB_BXH_RI", "HXB29", "HXB29_mRatNor1",
  "HXB_BXH_RI", "HXB3", "HXB3_mRatNor1",
  "HXB_BXH_RI", "HXB31", "HXB31_200315",
  "HXB_BXH_RI", "HXB4", "HXB4_608_GSPMC",
  "HXB_BXH_RI", "HXB5", "HXB5_mRatNor1",
  "HXB_BXH_RI", "HXB7", "HXB7_mRatNor1",
  "FXLE_LEXF_RI", "FXLE12", "FXLE12",
  "FXLE_LEXF_RI", "FXLE13", "FXLE13",
  "FXLE_LEXF_RI", "FXLE14", "FXLE14",
  "FXLE_LEXF_RI", "FXLE15", "FXLE15",
  "FXLE_LEXF_RI", "FXLE16", "FXLE16_612_GSPMC",
  "FXLE_LEXF_RI", "FXLE17", "FXLE17",
  "FXLE_LEXF_RI", "FXLE18", "FXLE18_613_GSPMC",
  "FXLE_LEXF_RI", "FXLE19", "FXLE19_Stm_442",
  "FXLE_LEXF_RI", "FXLE20", "FXLE20",
  "FXLE_LEXF_RI", "FXLE21", "FXLE21",
  "FXLE_LEXF_RI", "FXLE22", "FXLE22",
  "FXLE_LEXF_RI", "FXLE23", "FXLE23",
  "FXLE_LEXF_RI", "FXLE24", "FXLE24",
  "FXLE_LEXF_RI", "FXLE25", "FXLE25",
  "FXLE_LEXF_RI", "FXLE26", "FXLE26",
  "FXLE_LEXF_RI", "LEXF10A", "LXF10A_610_GSPMC",
  "FXLE_LEXF_RI", "LEXF10B", "LEXF10B_Stm_447",
  "FXLE_LEXF_RI", "LEXF10C", "LEXF10C",
  "FXLE_LEXF_RI", "LEXF11", "LEXF11",
  "FXLE_LEXF_RI", "LEXF1A", "LEXF1A",
  "FXLE_LEXF_RI", "LEXF1C", "LEXF1C",
  "FXLE_LEXF_RI", "LEXF2A", "LEXF2A",
  "FXLE_LEXF_RI", "LEXF2B", "LEXF2B_9",
  "FXLE_LEXF_RI", "LEXF2C", "LEXF2C_Stm_443",
  "FXLE_LEXF_RI", "LEXF3", "LEXF3_615_GSPMC",
  "FXLE_LEXF_RI", "LEXF4", "LEXF4",
  "FXLE_LEXF_RI", "LEXF5", "LEXF5",
  "FXLE_LEXF_RI", "LEXF6B", "LEXF6B",
  "FXLE_LEXF_RI", "LEXF7A", "LEXF7A",
  "FXLE_LEXF_RI", "LEXF7B", "LEXF7B",
  "FXLE_LEXF_RI", "LEXF7C", "LEXF7C",
  "FXLE_LEXF_RI", "LEXF8A", "LEXF8A",
  "FXLE_LEXF_RI", "LEXF8D", "LEXF8D_Stm_449",
  "FXLE_LEXF_RI", "LEXF9", "LEXF9"
) %>%
  mutate(
    is_study_background = official_strain %in% c(
      "SHR/OlaIpcvMcwi", "HXB10", "F344/StmMcwi", "LE/StmMcwi",
      "BXH6", "HXB2", "HXB31", "HXB23", "BN-Lx/CubMcwi"
    )
  )

header.command <- sprintf(
  "gzip -dc %s | awk '/^#CHROM/{print; exit}'",
  shQuote(gvcf.files[[1]])
)
header.line <- system(header.command, intern = TRUE)
gvcf.sample.names <- strsplit(header.line, "\t", fixed = TRUE)[[1]][-seq_len(9)]

expected.group.counts <- c(classic = 30L, HXB_BXH_RI = 30L, FXLE_LEXF_RI = 34L)
observed.group.counts <- table(df.hrdp.map$panel_group)
if (
  nrow(df.hrdp.map) != 94L ||
  n_distinct(df.hrdp.map$official_strain) != 94L ||
  n_distinct(df.hrdp.map$gvcf_sample) != 94L ||
  !all(
    as.integer(observed.group.counts[names(expected.group.counts)]) ==
      as.integer(expected.group.counts)
  )
) {
  stop("The HRDP mapping must contain 30 classic, 30 HXB/BXH, and 34 FXLE/LEXF strains.")
}
if (!all(df.hrdp.map$gvcf_sample %in% gvcf.sample.names)) {
  stop("At least one selected HRDP specimen is absent from the gVCF header.")
}
if (sum(df.hrdp.map$is_study_background) != 9L) {
  stop("The study panel must map to exactly nine unique inbred backgrounds.")
}

panel.names <- df.hrdp.map$gvcf_sample
study.names <- df.hrdp.map$gvcf_sample[df.hrdp.map$is_study_background]
write_tsv(df.hrdp.map, file.path(output.dir, "hrdp_panel_sample_mapping.tsv"))

########################
# 3. Autosomal SNP and indel coverage
########################

run.chromosome <- function(chromosome) {
  command <- sprintf(
    "gzip -dc %s | awk -v panel_names=%s -v study_names=%s -f %s",
    shQuote(gvcf.files[[chromosome]]),
    shQuote(paste(panel.names, collapse = ",")),
    shQuote(paste(study.names, collapse = ",")),
    shQuote(awk.file)
  )
  output <- system(command, intern = TRUE)
  if (!length(output)) stop("No coverage output for ", chromosome)
  fread(text = paste(output, collapse = "\n"), data.table = FALSE) %>%
    mutate(chromosome = chromosome, .before = 1)
}

workers <- min(4L, parallel::detectCores(logical = FALSE), length(chromosomes))
coverage.parts <- if (.Platform$OS.type == "unix") {
  parallel::mclapply(chromosomes, run.chromosome, mc.cores = workers)
} else {
  lapply(chromosomes, run.chromosome)
}

df.coverage.by.chromosome <- bind_rows(coverage.parts) %>%
  mutate(
    coverage_percent = 100 * n_study_background_variants /
      n_hrdp_panel_variants
  )

df.coverage.summary <- df.coverage.by.chromosome %>%
  group_by(measurement_unit, variant_type, qual_threshold, maf_threshold) %>%
  summarize(
    n_hrdp_panel_variants = sum(n_hrdp_panel_variants),
    n_study_background_variants = sum(n_study_background_variants),
    .groups = "drop"
  ) %>%
  mutate(
    coverage_percent = 100 * n_study_background_variants /
      n_hrdp_panel_variants
  )

write_tsv(
  df.coverage.by.chromosome,
  file.path(output.dir, "hrdp_variant_coverage_by_chromosome.tsv")
)
write_tsv(
  df.coverage.summary,
  file.path(output.dir, "hrdp_variant_coverage_summary.tsv")
)

########################
# 4. Figure and run metadata
########################

df.plot <- df.coverage.summary %>%
  filter(
    measurement_unit == "site",
    variant_type %in% c("SNP", "INDEL"),
    maf_threshold == 0.1
  ) %>%
  mutate(quality = factor(qual_threshold, levels = c(30, 40)))

p.coverage <- ggplot(
  df.plot,
  aes(variant_type, coverage_percent, fill = quality)
) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  geom_text(
    aes(label = sprintf("%.1f%%", coverage_percent)),
    position = position_dodge(width = 0.75),
    vjust = -0.35,
    size = 3.6
  ) +
  scale_fill_manual(values = c("#4E79A7", "#F28E2B")) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.06))
  ) +
  labs(
    x = NULL,
    y = "Official HRDP variant sites represented (%)",
    fill = "Minimum QUAL",
    title = "Common-variant coverage of the nine study strain backgrounds",
    subtitle = "94 HRDP strains; minor allele frequency > 0.10"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

for (extension in c("pdf", "png")) {
  ggsave(
    file.path(output.dir, paste0("hrdp_variant_coverage.", extension)),
    p.coverage,
    width = 7.5,
    height = 5.2,
    dpi = 300
  )
}

df.run.metadata <- tribble(
  ~field, ~value,
  "input_directory", gvcf.dir,
  "n_gvcf_specimens", as.character(length(gvcf.sample.names)),
  "denominator_definition", paste(
    "alternate alleles observed in one representative specimen for each of",
    "94 HRDP strains: 30 classic, 30 HXB/BXH RI, and 34 FXLE/LEXF RI"
  ),
  "study_definition", paste(
    "nine unique inbred backgrounds represented by the ten Hi-C libraries;",
    "the SHR x BN F1 was not counted as an additional allele source"
  ),
  "chromosomes", "chr1-chr20",
  "primary_quality_threshold", "QUAL >= 30",
  "sensitivity_quality_threshold", "QUAL >= 40",
  "primary_maf_threshold", "MAF > 0.10 across 94 HRDP strains",
  "sensitivity_maf_threshold", "MAF > 0.20 across 94 HRDP strains",
  "unfiltered_maf_threshold", "MAF > 0 retained for comparison",
  "accepted_filter_values", ". or PASS",
  "primary_measurement", "variant sites",
  "sensitivity_measurement", "alternate alleles at multiallelic records"
)
write_tsv(
  df.run.metadata,
  file.path(output.dir, "hrdp_variant_coverage_run_metadata.tsv")
)

print(df.coverage.summary)
