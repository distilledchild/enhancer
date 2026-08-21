# lintr: disable

getwd()
setwd("./enhancer") # Please set root directory: enhancer

funcs.file <- "./funcs_enhancer.R"
source(funcs.file)
list2env(resolve_enhancer_analysis_paths(funcs.file), envir = environment())

library("tidyverse")
library("GenomicRanges")
library("GenomeInfoDb")

options(tibble.width = Inf)
options(tibble.print_max = Inf)
options(tibble.max_extra_cols = Inf)
options(scipen = 999)

# BED / BEDPE / narrowPeak (ATAC, Hi-C) : 0-based -> start + 1 coord shift required
# GTF / FIMO (Ensembl, CTCF) : 1-based -> coord shif NOT required

################################################################################
# Resubmission Step 1: coordinate-system audit and normalization
#
# Target analysis coordinates for every genomic dataset:
# - genome assembly: mRatBN7.2 / rn7
# - chromosome names: chr-prefixed (for example, chr1 and chrX)
# - interval convention: 1-based, end-inclusive, as used by GenomicRanges
#
# This script reads the original coordinate sources, validates their documented
# conventions, creates normalized analysis objects, and stores those objects in
# revision/revision_main/cache_data. The main script loads this cache instead of
# rebuilding the coordinate objects on every run.
################################################################################

########################
# 0. Directories and source files
########################

dir.create(coord.cache.dir, recursive = TRUE, showWarnings = FALSE)

# Inputs are read from revision_main/inputs in the shared Google Drive tree.
if (!dir.exists(data.dir) || !dir.exists(hiccups.loop.root)) {
  stop(
    paste0(
      "The Google Drive input bundle is incomplete. Expected data under: ",
      input.root
    ),
    call. = FALSE
  )
}

# Record the upstream Juicer and HiCCUPS provenance associated with the loops.
df.hic.source.provenance <- tibble(
  juicer_pipeline_version = "1.6",
  bwa_version = "0.7.17-r1188",
  juicer_stage = "chimeric",
  genome_assembly = "rn7",
  restriction_enzyme_option = "Arima",
  restriction_site_file = "rn7chr_juicer_arima4.txt",
  threads = 30L,
  java_max_heap = "64g",
  hiccups_juicer_tools_version = "1.22.01"
)

loop.file.metadata <- tribble(
  ~sample, ~strain,
  "592BB", "SHR/OlaIpcv",
  "607", "HXB10",
  "74AA", "F344/Stm",
  "A2DB", "LE/Stm",
  "D765A", "BXH6",
  "DA08A", "HXB2",
  "DA21A", "SHR/OlaIpcvxBN/NHsdMcwi",
  "DA68A", "HXB31",
  "DBA9A", "HXB23",
  "DE8BA", "BN-Lx"
) %>%
  mutate(
    file = file.path(
      hiccups.loop.root,
      sample,
      "hiccups_5k10k25k",
      "merged_loops.bedpe"
    )
  )
########################################
# input file paths
# ctcf: fimo.4.tsv
# atac: Duttke2022_snATAC_peaks_rn7.narrowPeak
# ensembl: Rattus_norvegicus.mRatBN7.2.113.gtf
# epd:Rn_EPDnew_001_rn6.bed, promoter_coordinate.txt, rn6ToRn7.over.chain, promoter_ensembl.txt
# library complexity: library_complexity.tsv
########################################

ctcf.motif.file <- file.path(
  data.dir,
  "ctcf/submission/E4/fimo.4.tsv"
)
atac.peak.file <- file.path(
  data.dir,
  "Duttke2022_snATAC_peaks_rn7.narrowPeak"
)
ensembl.gtf.file <- file.path(
  data.dir,
  "Rattus_norvegicus.mRatBN7.2.113.gtf"
)
epd.rn6.bed.file <- file.path(
  data.dir,
  "epdnew/001/Rn_EPDnew_001_rn6.bed"
)
epd.coordinate.file <- file.path(
  data.dir,
  "epdnew/001/db/promoter_coordinate.txt"
)
epd.rn6.to.rn7.chain.file <- file.path(
  data.dir,
  "epdnew/rn6ToRn7.over.chain"
)
epd.promoter.mapping.file <- file.path(
  data.dir,
  "epdnew/001/db/promoter_ensembl.txt"
)
library.complexity.file <- file.path(
  data.dir,
  "library_complexity.tsv"
)

# Keep the optional true genetic-distance input configurable without requiring
# it for cache construction or the main pooled-resource analysis.
genetic.distance.file <- Sys.getenv(
  "HRDP_GENETIC_DISTANCE_FILE",
  unset = file.path(
    revision.dir,
    "hrdp_genotype_diversity",
    "results",
    "plink2_primary",
    "hrdp_plink2_ibs_distance.tsv"
  )
)

# Store all source and downstream input paths in one cached provenance table so
# the main script does not need to redefine the same paths.
df.analysis.input.files <- bind_rows(
  loop.file.metadata %>%
    transmute(
      input_name = str_c("hiccups_loop_", sample),
      input_path = file,
      input_group = "coordinate_source",
      required = TRUE
    ),
  tribble(
    ~input_name, ~input_path, ~input_group, ~required,
    "ctcf_motif", ctcf.motif.file, "coordinate_source", TRUE,
    "atac_peak", atac.peak.file, "coordinate_source", TRUE,
    "ensembl_gtf", ensembl.gtf.file, "coordinate_source", TRUE,
    "epd_rn6_bed", epd.rn6.bed.file, "coordinate_source", TRUE,
    "epd_coordinate", epd.coordinate.file, "coordinate_source", TRUE,
    "epd_rn6_to_rn7_chain", epd.rn6.to.rn7.chain.file,
    "coordinate_source", TRUE,
    "epd_promoter_mapping", epd.promoter.mapping.file,
    "coordinate_source", TRUE,
    "library_complexity", library.complexity.file,
    "downstream_depth_QC", TRUE,
    "genetic_distance", genetic.distance.file,
    "optional_downstream_analysis", FALSE
  )
)

required.files <- df.analysis.input.files %>%
  filter(required) %>%
  pull(input_path)

missing.files <- required.files[!file.exists(required.files)]
if (length(missing.files) > 0L) {
  stop(
    "Missing required input file(s):\n",
    paste(missing.files, collapse = "\n"),
    "\nMake the required Google Drive files available offline before rebuilding the cache.",
    call. = FALSE
  )
} else {
  message(
    "All required input files verified (",
    length(required.files),
    " files present)."
  )
}

################################################################################
# 1. Coordinate-system audit and normalization
################################################################################

########################
# 1-1. HiCCUPS loop anchors
#
# Source: BEDPE, 0-based and half-open.
# Target: 1-based and end-inclusive.
# Conversion: add 1 to x1/y1; leave x2/y2 unchanged.
########################

df.sample.loop.raw <- read_hiccups_loop_files(loop.file.metadata)
df.sample.loop.raw %>% head(4)

# 1. converting coordinate system of HiCCUPS loop anchors from 0-based and half-open to 1-based and end-inclusive.
# 2. filtering out loops whose length is longer than 2 Mb.
df.sample.loop.1based <- normalize_hiccups_loop_coordinates(
  df.sample.loop.raw,
  max.loop.distance = 2000000L
)
df.sample.loop.1based %>% count(passes_lt2mb)
#   passes_lt2mb     n
# 1 FALSE         1463
# 2 TRUE         57537
df.sample.loop.1based %>% head()

# three columns: n_supporting_libraries n_supporting_strains supporting_samples
list.loop.resource <- build_pooled_hiccups_loop_resource(df.sample.loop.1based)
names(list.loop.resource)
purrr::map(list.loop.resource, dim)
head(list.loop.resource$distinct_2mb)

df.loop.pooled.support.summary <- list.loop.resource$support
df.loop.distinct <- list.loop.resource$distinct
df.loop.distinct.2mb <- list.loop.resource$distinct_2mb

df.loop.pooled.support.summary %>% head(2)
df.loop.distinct %>% head(2)
df.loop.distinct.2mb %>% head(2)

# These counts were verified from the current original merged_loops.bedpe files on 2026-07-22 (from sb option processing)
# the older copied inputs contained 58,992/31,773/31,019.
df.loop.source.count.check <- check_hiccups_loop_counts(
  df.sample.loop.1based = df.sample.loop.1based,
  df.loop.distinct = df.loop.distinct,
  df.loop.distinct.2mb = df.loop.distinct.2mb,
  expected.n = c(59000L, 31778L, 31021L)
)

message("Original sample-level HiCCUPS rows: ", nrow(df.sample.loop.raw))
message("All distinct HiCCUPS loops: ", n_distinct(df.loop.distinct$loop_id))
message("Pooled loop resource (<2 Mb): ", nrow(df.loop.distinct.2mb))

########################
# 1-2. CTCF FIMO motif coordinates
#
# Source and target: 1-based and end-inclusive.
# Conversion: NO coordinate shift; validate and deduplicate intervals.
########################

gr.ctcf.motif <- read_ctcf_fimo(ctcf.motif.file)
df.ctcf.fimo.summary <- S4Vectors::metadata(
  gr.ctcf.motif
)$fimo_summary

# Verify that full FIMO provenance and representative statistical fields remain
# attached after exact-coordinate interval consolidation.
required.ctcf.metadata.columns <- c(
  "best_motif_id", "best_score", "best_p_value", "best_q_value",
  "n_fimo_predictions", "n_distinct_motif_ids"
)
if (
  is.null(df.ctcf.fimo.summary) ||
    !all(required.ctcf.metadata.columns %in% names(mcols(gr.ctcf.motif)))
) {
  stop(
    "CTCF cache preparation did not preserve FIMO motif/statistical metadata.",
    call. = FALSE
  )
}
message("Full CTCF FIMO cache summary:")
print(df.ctcf.fimo.summary)

# PASS!!

########################
# 1-3. Ensembl GTF transcript coordinates
#
# Source: GTF, 1-based and end-inclusive, chromosomes without a chr prefix.
# Target: 1-based and end-inclusive with chr-prefixed chromosome names.
########################

# feature == "transcript"
# chr
# checking coord & strand (+/-): source_coordinate_system, analysis_coordinate_system
# keeping 'attribute' column as is
df.transcript.ensembl.rn7.1based <- read_ensembl_gtf_transcripts(ensembl.gtf.file)
df.transcript.ensembl.rn7.1based %>% head(3)

########################
# 1-4. EPD promoter coordinates
#
# Source: original rn6 EPD BED, coordinate table, and rn6-to-rn7 chain.
# Target: strand-aware one-base rn7 TSS and an 81-bp promoter interval.
########################

# source_coordinate_system, analysis_coordinate_system
df.promoter.epd.rn7.1based <- read_epd_rn6_liftover_promoters(
  bed.file = epd.rn6.bed.file,
  chain.file = epd.rn6.to.rn7.chain.file,
  coordinate.file = epd.coordinate.file,
  mapping.file = epd.promoter.mapping.file
)
df.promoter.epd.rn7.1based %>% head(3)

# Verify the production EPD reconstruction directly from the original rn6
# source and liftover result without requiring a copied legacy rn7 annotation.
epd.liftover.counts <- c(
  attr(df.promoter.epd.rn7.1based, "n_rn6_input"),
  attr(df.promoter.epd.rn7.1based, "n_rn7_lifted"),
  attr(df.promoter.epd.rn7.1based, "n_failed_liftover"),
  nrow(df.promoter.epd.rn7.1based)
)
if (!identical(epd.liftover.counts, c(12601L, 12528L, 73L, 12524L))) {
  stop(
    "EPD source counts do not match the verified coordinate reconstruction.",
    call. = FALSE
  )
}

########################
# 1-5. Frontal-cortex ATAC-seq peak coordinates
#
# Source: narrowPeak/BED, 0-based and half-open.
# Target: 1-based and end-inclusive.
########################

# coord shift to 1-based and end-inclusive
gr.atac.rn7.1based <- read_atac_narrowpeak(atac.peak.file)

########################
# 1-6. Coordinate audit summary
########################
# ==============================================================================
# Coordinate System Audit & Normalization Summary
# ------------------------------------------------------------------------------
# Input Dataset (input_name)   | Source Coordinate System | Coordinate Conversion                | Final Object (normalized_object) | Record Count | Status
# Hi-C Sample Loops (10 Libs)  | BEDPE: 0-based, half-open| start + 1; end unchanged             | df.sample.loop.1based            | 59,000       | passed
# Hi-C Pooled Distinct Loops   | 1-based inclusive        | deduplicate by loop_id (no shift)    | df.loop.distinct                 | 31,778       | passed
# CTCF FIMO Motifs             | FIMO: 1-based, inclusive | exact coordinates retained (no shift)| gr.ctcf.motif                    | 462,810      | passed
# Ensembl Transcripts (mRatBN7)| GTF: 1-based, inclusive  | chr prefix normalized (no shift)     | df.transcript.ensembl.rn7.1based | 54,993       | passed
# EPD Promoters                | rn6 BED (1-based TSS)    | strand-aware TSS liftOver to rn7     | df.promoter.epd.rn7.1based       | 12,524       | passed
# Frontal-Cortex ATAC Peaks    | narrowPeak: 0-based      | start + 1; end unchanged             | gr.atac.rn7.1based               | 154,007      | passed
# ==============================================================================

df.coordinate.system.audit <- tribble(
  ~input_name,
  ~input_path,
  ~source_coordinate_system,
  ~coordinate_system_basis,
  ~target_coordinate_system,
  ~conversion,
  ~normalized_object,
  ~n_records,
  ~validation_result,
  "HiCCUPS loop anchors from 10 libraries",
  str_c(loop.file.metadata$file, collapse = "; "),
  "BEDPE: 0-based, half-open",
  "HiCCUPS BEDPE format; normalized widths equal 5K/10K/25K bins",
  "rn7, chr-prefixed, 1-based inclusive",
  "start + 1; end unchanged",
  "df.sample.loop.1based",
  nrow(df.sample.loop.1based),
  "passed",
  "Pooled distinct HiCCUPS loop resource",
  "derived from df.sample.loop.1based",
  "Derived from normalized source intervals",
  "Distinct genomic BEDPE interval and resolution across 10 libraries",
  "rn7, chr-prefixed, 1-based inclusive",
  "deduplicate by stable loop_id; no coordinate shift",
  "df.loop.distinct",
  nrow(df.loop.distinct),
  "passed",
  "CTCF FIMO motifs",
  ctcf.motif.file,
  "FIMO genomic coordinates: 1-based, inclusive",
  "FIMO genomic output provenance and valid inclusive intervals",
  "rn7, chr-prefixed, 1-based inclusive",
  "no coordinate shift",
  "gr.ctcf.motif",
  length(gr.ctcf.motif),
  "passed_full_FIMO_metadata_preserved",
  "Ensembl transcript annotation",
  ensembl.gtf.file,
  "GTF: 1-based, inclusive; chromosomes not chr-prefixed",
  "GTF specification, transcript intervals, and strand validation",
  "rn7, chr-prefixed, 1-based inclusive",
  "chromosome prefix normalized; coordinates unchanged",
  "df.transcript.ensembl.rn7.1based",
  nrow(df.transcript.ensembl.rn7.1based),
  "passed",
  "EPD promoter annotation",
  str_c(
    epd.rn6.bed.file,
    epd.coordinate.file,
    epd.rn6.to.rn7.chain.file,
    epd.promoter.mapping.file,
    sep = "; "
  ),
  "rn6 BED plus EPD 1-based TSS coordinate table",
  "BED fields, exact coordinate-table agreement, and liftOver validation",
  "rn7, chr-prefixed, 1-based inclusive",
  "strand-aware rn6 TSS; one-base liftOver; rebuild +/-40-bp promoter",
  "df.promoter.epd.rn7.1based",
  nrow(df.promoter.epd.rn7.1based),
  "passed",
  "Frontal-cortex ATAC peaks",
  atac.peak.file,
  "narrowPeak/BED: 0-based, half-open",
  "narrowPeak format and valid BED intervals",
  "rn7, chr-prefixed, 1-based inclusive",
  "start + 1; end unchanged",
  "gr.atac.rn7.1based",
  length(gr.atac.rn7.1based),
  "passed"
)

########################
# 1-7. Save normalized coordinate cache
########################

df.coord.cache.manifest <- save_coordinate_cache_objects(
  cache.dir = coord.cache.dir,
  envir = environment()
)

message("Coordinate cache written to: ", coord.cache.dir)
print(df.coord.cache.manifest)
