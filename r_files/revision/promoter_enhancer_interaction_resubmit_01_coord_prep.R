# lintr: disable
library("tidyverse")
library("GenomicRanges")
library("GenomeInfoDb")

options(tibble.width = Inf)
options(tibble.print_max = Inf)
options(tibble.max_extra_cols = Inf)
options(scipen = 999)

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
# revision/cache_data. The main resubmission script loads this cache instead of
# rebuilding the coordinate objects on every run.
################################################################################

########################
# 0. Directories and source files
########################

r.files.dir <- path.expand("~/Desktop/playground/enhancer/r_files")
if (!dir.exists(r.files.dir)) {
  r.files.dir <- path.expand("~/dropbox/Gateway_to_Hao/enhancer/r_files")
}
if (!dir.exists(r.files.dir)) {
  stop("Cannot locate the enhancer/r_files directory.", call. = FALSE)
}

setwd(r.files.dir)

revision.dir <- file.path(r.files.dir, "revision")
coord.cache.dir <- file.path(revision.dir, "cache_data")
dir.create(coord.cache.dir, recursive = TRUE, showWarnings = FALSE)

dropbox.root <- path.expand("~/dropbox/Gateway_to_Hao/enhancer")
hiccups.loop.root <- path.expand(
  paste0(
    "~/Library/CloudStorage/Dropbox-UTHSCGGI/K P/Gateway_to_Hao/",
    "hic/2023A/hic30_w_sb_options"
  )
)

source(file.path(r.files.dir, "funcs.R"))

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

ctcf.motif.file <- file.path(
  dropbox.root,
  "data/ctcf/submission/E4/fimo.4.tsv"
)
atac.peak.file <- file.path(
  dropbox.root,
  "data/Duttke2022_snATAC_peaks_rn7.narrowPeak"
)
ensembl.gtf.file <- file.path(
  dropbox.root,
  "data/Rattus_norvegicus.mRatBN7.2.113.gtf"
)
epd.rn6.bed.file <- file.path(
  dropbox.root,
  "data/epdnew/001/Rn_EPDnew_001_rn6.bed"
)
epd.coordinate.file <- file.path(
  dropbox.root,
  "data/epdnew/001/db/promoter_coordinate.txt"
)
epd.rn6.to.rn7.chain.file <- file.path(
  dropbox.root,
  "data/epdnew/rn6ToRn7.over.chain"
)
epd.promoter.mapping.file <- file.path(
  dropbox.root,
  "data/epdnew/001/db/promoter_ensembl.txt"
)
library.complexity.file <- file.path(
  dropbox.root,
  "data/library_complexity_592BB.tsv"
)

# Keep the optional true genetic-distance input configurable without requiring
# it for cache construction or the main pooled-resource analysis.
genetic.distance.file <- Sys.getenv(
  "HRDP_GENETIC_DISTANCE_FILE",
  unset = file.path(revision.dir, "hrdp_genetic_distance.tsv")
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
    "\nMake required Dropbox files available offline before rebuilding the cache.",
    call. = FALSE
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

df.sample.loop.coordinate.normalized <- normalize_hiccups_loop_coordinates(
  df.sample.loop.raw,
  max.loop.distance = 2000000L
)

list.loop.resource <- build_pooled_hiccups_loop_resource(
  df.sample.loop.coordinate.normalized
)
df.loop.pooled.support.summary <- list.loop.resource$support
df.loop.distinct.coordinate.normalized <- list.loop.resource$distinct
df.loop.universe <- list.loop.resource$universe

# These counts were verified from the current original merged_loops.bedpe files
# on 2026-07-22; the older copied inputs contained 58,992/31,773/31,019.
df.loop.source.count.check <- check_hiccups_loop_counts(
  df.sample.loop.coordinate.normalized =
    df.sample.loop.coordinate.normalized,
  df.loop.distinct.coordinate.normalized =
    df.loop.distinct.coordinate.normalized,
  df.loop.universe = df.loop.universe,
  expected.n = c(59000L, 31778L, 31021L)
)

message("Original sample-level HiCCUPS rows: ", nrow(df.sample.loop.raw))
message(
  "All distinct HiCCUPS loops: ",
  n_distinct(df.loop.distinct.coordinate.normalized$loop_id)
)
message("Pooled loop resource (<2 Mb): ", nrow(df.loop.universe))

########################
# 1-2. CTCF FIMO motif coordinates
#
# Source and target: 1-based and end-inclusive.
# Conversion: no coordinate shift; validate and deduplicate intervals.
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

########################
# 1-3. Ensembl GTF transcript coordinates
#
# Source: GTF, 1-based and end-inclusive, chromosomes without a chr prefix.
# Target: 1-based and end-inclusive with chr-prefixed chromosome names.
########################

df.ensembl.transcript.coordinate.normalized <-
  read_ensembl_gtf_transcripts(
  ensembl.gtf.file
)

########################
# 1-4. EPD promoter coordinates
#
# Source: original rn6 EPD BED, coordinate table, and rn6-to-rn7 chain.
# Target: strand-aware one-base rn7 TSS and an 81-bp promoter interval.
########################

df.promoter.annotation.coordinate.normalized <-
  read_epd_rn6_liftover_promoters(
  bed.file = epd.rn6.bed.file,
  chain.file = epd.rn6.to.rn7.chain.file,
  coordinate.file = epd.coordinate.file,
  mapping.file = epd.promoter.mapping.file
)

# Verify the production EPD reconstruction directly from the original rn6
# source and liftover result without requiring a copied legacy rn7 annotation.
epd.normalization.counts <- c(
  attr(df.promoter.annotation.coordinate.normalized, "n_rn6_input"),
  attr(df.promoter.annotation.coordinate.normalized, "n_rn7_lifted"),
  attr(df.promoter.annotation.coordinate.normalized, "n_failed_liftover"),
  nrow(df.promoter.annotation.coordinate.normalized)
)
if (!identical(epd.normalization.counts, c(12601L, 12528L, 73L, 12524L))) {
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

gr.atac <- read_atac_narrowpeak(atac.peak.file)

########################
# 1-6. Coordinate audit summary
########################

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
  "df.sample.loop.coordinate.normalized",
  nrow(df.sample.loop.coordinate.normalized),
  "passed",
  "Pooled distinct HiCCUPS loop resource",
  "derived from df.sample.loop.coordinate.normalized",
  "Derived from normalized source intervals",
  "Distinct genomic BEDPE interval and resolution across 10 libraries",
  "rn7, chr-prefixed, 1-based inclusive",
  "deduplicate by stable loop_id; no coordinate shift",
  "df.loop.distinct.coordinate.normalized",
  nrow(df.loop.distinct.coordinate.normalized),
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
  "df.ensembl.transcript.coordinate.normalized",
  nrow(df.ensembl.transcript.coordinate.normalized),
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
  "df.promoter.annotation.coordinate.normalized",
  nrow(df.promoter.annotation.coordinate.normalized),
  "passed",
  "Frontal-cortex ATAC peaks",
  atac.peak.file,
  "narrowPeak/BED: 0-based, half-open",
  "narrowPeak format and valid BED intervals",
  "rn7, chr-prefixed, 1-based inclusive",
  "start + 1; end unchanged",
  "gr.atac",
  length(gr.atac),
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
