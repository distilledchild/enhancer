options(stringsAsFactors = FALSE)
library(GenomicRanges)
library(dplyr)
library(tidyverse)

######################################################
# 0. Path Definitions and Setup
######################################################

setwd("/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/nicotine")
# Working directory set to: getwd()

# Input File Paths
snp_bim_file <- "/Users/pete/Library/CloudStorage/OneDrive-UniversityofTennessee/v4/HS_genotypes_v4.bim"
final_loop_file <- "/Users/pete/Desktop/playground/enhancer/figures/submission/lt2mb/df_final_loop_sub.4.any.lt2mb.ENSEMBL.mid.mid.final.filter.200kb.csv"
ensembl_gtf_file <- "/Users/pete/Library/CloudStorage/Dropbox-UTHSCGGI/K P/Gateway_to_Hao/enhancer/data/Rattus_norvegicus.mRatBN7.2.113.gtf"

# Francesca Telese / Duttke et al. 2022 snATAC-seq peak file.
# H-MAGMA uses only distal ATAC signal after removing known TSS/promoter regions.
atac_file <- "/Users/pete/Library/CloudStorage/Dropbox-UTHSCGGI/K P/Gateway_to_Hao/enhancer/data/Duttke2022_snATAC_peaks_rn7.narrowPeak"
hmagma_annot_name <- "hs.hic.duttke_telese_atac.non_tss_promoter_filtered.annot"
hmagma_counts_name <- "cmagma_vs_hmagma_duttke_telese_non_tss_promoter_snp_counts.csv"

# Output Directory Setup
output_dir <- "/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/nicotine"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Verify Input Files
if (!file.exists(snp_bim_file)) {
  stop("SNP BIM file not found at: ", snp_bim_file)
}
if (!file.exists(final_loop_file)) {
  stop("Hi-C Loop file not found at: ", final_loop_file)
}
if (!file.exists(ensembl_gtf_file)) {
  stop("Ensembl GTF file not found at: ", ensembl_gtf_file)
}
if (!file.exists(atac_file)) {
  stop("ATAC-seq peak file not found at: ", atac_file)
}

# Helper function to parse attributes efficiently from GTF
parse_attribute <- function(attributes, key) {
  pattern <- paste0(key, " \"(.*?)\"")
  str_match(attributes, pattern)[, 2]
}

# Helper function to read narrowPeak safely without rtracklayer dependency
read_peak_file <- function(filepath) {
  df <- read.table(filepath, header = FALSE, stringsAsFactors = FALSE)
  chroms <- gsub("^chr", "", as.character(df$V1))
  GRanges(
    seqnames = chroms,
    ranges = IRanges(start = as.numeric(df$V2) + 1, end = as.numeric(df$V3))
  )
}

######################################################
# 1. Load and Parse Ensembl Reference (Local GTF)
######################################################
# Reading and parsing local Ensembl GTF: ensembl_gtf_file

# Read raw GTF lines for exons and transcripts
gtf_raw <- read_tsv(
  ensembl_gtf_file,
  comment = "#",
  col_names = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"),
  col_types = cols(
    chr = col_character(),
    start = col_double(),
    end = col_double(),
    strand = col_character(),
    feature = col_character(),
    attribute = col_character()
  )
) %>%
  filter(feature %in% c("exon", "transcript")) %>%
  filter(chr %in% c(as.character(1:20), "X"))
gtf_raw # 580,600
gtf_raw %>% distinct(feature)

# Extracting attributes for transcripts...
transcripts <- gtf_raw %>%
  filter(feature == "transcript") %>%
  mutate(
    gene_id = parse_attribute(attribute, "gene_id"),
    gene_name = parse_attribute(attribute, "gene_name"),
    transcript_id = parse_attribute(attribute, "transcript_id"),
    transcript_biotype = parse_attribute(attribute, "transcript_biotype")
  ) %>%
  filter(transcript_biotype == "protein_coding")

# Extracting attributes for exons...
exons <- gtf_raw %>%
  filter(feature == "exon") %>%
  mutate(
    gene_id = parse_attribute(attribute, "gene_id"),
    gene_name = parse_attribute(attribute, "gene_name"),
    transcript_id = parse_attribute(attribute, "transcript_id"),
    transcript_biotype = parse_attribute(attribute, "transcript_biotype")
  ) %>%
  filter(transcript_biotype == "protein_coding")

# Extract common protein-coding genes
common_genes <- intersect(exons$gene_id, transcripts$gene_id)
# Total protein-coding genes identified: length(common_genes)

# 1-1. Define Exon coordinates
exon.from.gtf <- exons %>%
  filter(gene_id %in% common_genes)

# 1-2. Define Promoter coordinates
# promoter_start = transcript_start - 2000, promoter_end = transcript_start + 500 (strand == "+")
# promoter_start = transcript_end - 500, promoter_end = transcript_end + 2000 (strand == "-")
promoter.from.gtf <- transcripts %>%
  filter(gene_id %in% common_genes) %>%
  mutate(
    promoter_start = ifelse(strand == "+", start - 2000, end - 500),
    promoter_end = ifelse(strand == "+", start + 500, end + 2000)
  )

# Ensure coordinates are valid (cannot be negative)
promoter.from.gtf <- promoter.from.gtf %>%
  mutate(
    promoter_start = ifelse(promoter_start < 1, 1, promoter_start),
    promoter_end = ifelse(promoter_end < 1, 1, promoter_end)
  )

# 1-3. Create GenomicRanges Objects for Exon & Promoter
exonranges <- GRanges(
  seqnames = as.character(exon.from.gtf$chr),
  ranges = IRanges(start = exon.from.gtf$start, end = exon.from.gtf$end),
  strand = exon.from.gtf$strand,
  gene = as.character(exon.from.gtf$gene_id)
)

promoterranges <- GRanges(
  seqnames = as.character(promoter.from.gtf$chr),
  ranges = IRanges(start = promoter.from.gtf$promoter_start, end = promoter.from.gtf$promoter_end),
  strand = promoter.from.gtf$strand,
  gene = as.character(promoter.from.gtf$gene_id)
)

######################################################
# 2. Load SNP Data (.bim)
######################################################
# Loading SNP data from: snp_bim_file
df.init.snps <- read.table(snp_bim_file)
snps <- df.init.snps %>%
  dplyr::select(V1, V2, V4) %>%
  dplyr::rename(chr = V1, rsid = V2, Position = V4) %>%
  filter(chr %in% c(1:20, "X"))

snps.g <- GRanges(
  seqnames = as.character(snps$chr),
  ranges = IRanges(start = snps$Position, end = snps$Position),
  rsid = as.character(snps$rsid)
)
snps.g # 7358388
# Total SNPs loaded: length(snps.g)

######################################################
# 3. Create Gene Definition File (for MAGMA)
######################################################
# Generating gene definition file...
# Select unique genes with their min start and max end from GTF
genedef <- transcripts %>%
  filter(gene_id %in% common_genes) %>%
  group_by(gene_id, chr, strand) %>%
  summarise(
    start_pos = min(start),
    end_pos = max(end),
    .groups = "drop"
  ) %>%
  dplyr::rename(ensg = gene_id) %>%
  mutate(
    index = str_c(chr, start_pos, end_pos, sep = ":")
  ) %>%
  dplyr::select(index, chr, start = start_pos, end = end_pos, strand, ensg)

stopifnot(!any(duplicated(genedef$ensg)))
genedef_out_path <- file.path(output_dir, "genedef_ensembl.tsv")
write.table(genedef, file = genedef_out_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
# Gene definition saved to: genedef_out_path

######################################################
# 4. cMAGMA Annotation (Exons and Promoters ONLY)
######################################################
# Calculating cMAGMA (Exon + Promoter only) overlaps...

# Exon overlap
olap.snp.exon <- findOverlaps(snps.g, exonranges)
df.snpexon <- data.frame(
  rsid = as.character(mcols(snps.g[queryHits(olap.snp.exon)])$rsid),
  ensg = as.character(mcols(exonranges[subjectHits(olap.snp.exon)])$gene)
)

# Promoter overlap
olap.snp.pro <- findOverlaps(snps.g, promoterranges)
df.snpro <- data.frame(
  rsid = as.character(mcols(snps.g[queryHits(olap.snp.pro)])$rsid),
  ensg = as.character(mcols(promoterranges[subjectHits(olap.snp.pro)])$gene)
)

# Merge directly mapped SNPs
df.distinct.snp.cmagma <- bind_rows(df.snpexon, df.snpro) %>%
  distinct()

# Map back to index (chr:start:end)
df.distinct.snp.cmagma.merged <- df.distinct.snp.cmagma %>%
  left_join(genedef, by = c("ensg" = "ensg")) %>%
  filter(!is.na(index))

# Aggregate by gene (cMAGMA format)
snpagg.cmagma <- df.distinct.snp.cmagma.merged %>%
  group_by(index) %>%
  summarise(rsid = paste(unique(rsid), collapse = ", ")) %>%
  mutate(ensg = index) %>%
  dplyr::select(ensg, index, rsid)

snpagg.cmagma %>% count()

# Save cMAGMA annotation
cmagma_txt_path <- file.path(output_dir, "snpagg_exon_promoter_only.txt")
cmagma_final_path <- file.path(output_dir, "hs.exonpro.ONLY.annot")

write.table(snpagg.cmagma, file = cmagma_txt_path, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# Use sed to replace ', ' with '\t' to meet MAGMA formatting
system(sprintf("sed -e 's/, /\\t/g' < %s > %s", shQuote(cmagma_txt_path), shQuote(cmagma_final_path)), wait = TRUE)
# cMAGMA annotation successfully created at: cmagma_final_path

######################################################
# 5. H-MAGMA Annotation (Exons + Promoters + Hi-C Interaction + ATAC filtering)
######################################################
# Calculating H-MAGMA (Exon + Promoter + Hi-C Loops) overlaps...

# 5-1. Load ATAC-seq peaks and keep only distal open chromatin signal.
combined_atac_ranges <- GenomicRanges::reduce(read_peak_file(atac_file))

tss_points <- GRanges(
  seqnames = as.character(transcripts$chr),
  ranges = IRanges(
    start = ifelse(transcripts$strand == "+", transcripts$start, transcripts$end),
    end = ifelse(transcripts$strand == "+", transcripts$start, transcripts$end)
  ),
  strand = transcripts$strand
)
tss_regions <- promoters(tss_points, upstream = 1000, downstream = 1000)
tss_promoter_ranges <- GenomicRanges::reduce(c(tss_regions, promoterranges), ignore.strand = TRUE)
distal_atac_ranges <- GenomicRanges::setdiff(combined_atac_ranges, tss_promoter_ranges, ignore.strand = TRUE)

# Identify SNPs not mapped to exons or promoters
snpranges <- snps.g[!(mcols(snps.g)$rsid %in% df.distinct.snp.cmagma$rsid), ]
# Unmapped SNPs left for Hi-C analysis: length(snpranges)

# Filter SNPs by distal ATAC-seq peaks.
olap.snps.atac <- findOverlaps(snpranges, distal_atac_ranges)
snpranges_atac_filtered <- snpranges[unique(queryHits(olap.snps.atac)), ]
# SNPs remaining after distal ATAC filtering: length(snpranges_atac_filtered)

# Load Hi-C loops
df.final.loop <- read.csv(final_loop_file, stringsAsFactors = FALSE, check.names = FALSE)

hic.int1 <- df.final.loop %>%
  separate(loop.id, into = c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "resolution"), sep = "_") %>%
  mutate(across(c(start1, end1, start2, end2), as.numeric)) %>%
  mutate(
    chrom1 = str_remove(chrom1, "chr"),
    chrom2 = str_remove(chrom2, "chr")
  ) %>%
  dplyr::select(chrom1, start1, end1, chrom2, start2, end2, resolution)

hic.int2 <- hic.int1 %>%
  dplyr::select(chrom2, start2, end2, chrom1, start1, end1, resolution) %>%
  dplyr::rename(
    chrom1 = chrom2, start1 = start2, end1 = end2,
    chrom2 = chrom1, start2 = start1, end2 = end1
  )

hic.comb <- rbind(hic.int1, hic.int2)

hicranges <- GRanges(
  seqnames = hic.comb$chrom1,
  ranges = IRanges(start = hic.comb$start1 + 1, end = hic.comb$end1),
  int1 = hic.comb$start2,
  int2 = hic.comb$end2
)

# Overlap loop anchors with promoter regions
olap.hic.pro <- findOverlaps(hicranges, promoterranges)
generanges <- hicranges[queryHits(olap.hic.pro)]

# Assign metadata columns explicitly to avoid nested DataFrame structures
mcols(generanges)$int1 <- hicranges$int1[queryHits(olap.hic.pro)]
mcols(generanges)$int2 <- hicranges$int2[queryHits(olap.hic.pro)]
mcols(generanges)$gene <- promoterranges$gene[subjectHits(olap.hic.pro)]

# Create target gene-end anchor ranges (E2)
genebed <- data.frame(
  chr = seqnames(generanges),
  snp.start = generanges$int1,
  snp.end = generanges$int2,
  gene.start = start(generanges),
  gene.end = start(generanges) + width(generanges) - 1,
  ensg = generanges$gene
) %>%
  distinct()

genesnpranges <- GRanges(
  seqnames = genebed$chr,
  ranges = IRanges(start = genebed$snp.start + 1, end = genebed$snp.end),
  ensg = genebed$ensg
)

# Overlap SNPs with target anchor ranges:
# A) Unfiltered (Original H-MAGMA)
olap.snps.hic.unfiltered <- findOverlaps(snpranges, genesnpranges)
df.snpint.unfiltered <- data.frame(
  rsid = mcols(snpranges)$rsid[queryHits(olap.snps.hic.unfiltered)],
  ensg = mcols(genesnpranges)$ensg[subjectHits(olap.snps.hic.unfiltered)]
) %>%
  distinct()

# B) Filtered (ATAC-seq filtered H-MAGMA)
olap.snps.hic.filtered <- findOverlaps(snpranges_atac_filtered, genesnpranges)
df.snpint.filtered <- data.frame(
  rsid = mcols(snpranges_atac_filtered)$rsid[queryHits(olap.snps.hic.filtered)],
  ensg = mcols(genesnpranges)$ensg[subjectHits(olap.snps.hic.filtered)]
) %>%
  distinct()

# Combine all SNP relationships: Exon + Promoter + Hi-C Loops
snpcomb.hmagma.unfiltered <- unique(rbind(df.snpint.unfiltered, df.snpro, df.snpexon))
snpcomb.hmagma.filtered <- unique(rbind(df.snpint.filtered, df.snpro, df.snpexon))

# Aggregate H-MAGMA (We will output the ATAC-filtered version as our final H-MAGMA annot)
snpagg.hmagma <- snpcomb.hmagma.filtered %>%
  group_by(ensg) %>%
  summarise(rsid = paste(unique(as.character(rsid)), collapse = ", ")) %>%
  ungroup()

snpagg.hmagma$index <- genedef$index[match(snpagg.hmagma$ensg, genedef$ensg)]

snpaggconv.hmagma <- snpagg.hmagma %>%
  filter(!is.na(index)) %>%
  dplyr::select(ensg, index, rsid)

# Save H-MAGMA annotation (ATAC filtered)
hmagma_txt_path <- file.path(output_dir, "SNP_aggregate_transcript.txt")
hmagma_final_path <- file.path(output_dir, hmagma_annot_name)

write.table(snpaggconv.hmagma, file = hmagma_txt_path, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# Use sed to replace ', ' with '\t' to meet MAGMA formatting
system(sprintf("sed -e 's/, /\\t/g' < %s > %s", shQuote(hmagma_txt_path), shQuote(hmagma_final_path)), wait = TRUE)
# H-MAGMA annotation successfully created at: hmagma_final_path

# Preprocessing complete. Output files are in: output_dir

######################################################
# 6. Comparison: cMAGMA vs H-MAGMA (Unfiltered & ATAC Filtered)
######################################################
# Count SNPs per gene for all 3 cases
cmagma_counts <- df.distinct.snp.cmagma %>%
  group_by(ensg) %>%
  summarise(n_snp_cmagma = n_distinct(rsid), .groups = "drop")

hmagma_unfiltered_counts <- snpcomb.hmagma.unfiltered %>%
  group_by(ensg) %>%
  summarise(n_snp_hmagma_unfiltered = n_distinct(rsid), .groups = "drop")

hmagma_filtered_counts <- snpcomb.hmagma.filtered %>%
  group_by(ensg) %>%
  summarise(n_snp_hmagma_filtered = n_distinct(rsid), .groups = "drop")

# Merge counts and calculate differences
comparison_df <- full_join(cmagma_counts, hmagma_unfiltered_counts, by = "ensg") %>%
  full_join(hmagma_filtered_counts, by = "ensg") %>%
  mutate(
    n_snp_cmagma = coalesce(n_snp_cmagma, 0L),
    n_snp_hmagma_unfiltered = coalesce(n_snp_hmagma_unfiltered, 0L),
    n_snp_hmagma_filtered = coalesce(n_snp_hmagma_filtered, 0L),
    added_snps_unfiltered = n_snp_hmagma_unfiltered - n_snp_cmagma,
    added_snps_filtered = n_snp_hmagma_filtered - n_snp_cmagma,
    noise_snps_removed = n_snp_hmagma_unfiltered - n_snp_hmagma_filtered
  )

# Add Gene Name mapping
gene_mapping <- transcripts %>%
  dplyr::select(gene_id, gene_name) %>%
  distinct()

comparison_df <- comparison_df %>%
  left_join(gene_mapping, by = c("ensg" = "gene_id")) %>%
  dplyr::select(
    ensg, gene_name,
    n_snp_cmagma,
    n_snp_hmagma_unfiltered,
    n_snp_hmagma_filtered,
    added_snps_unfiltered,
    added_snps_filtered,
    noise_snps_removed
  ) %>%
  arrange(desc(added_snps_filtered))

# Save the detailed 3-way comparison to CSV
comparison_out_path <- file.path(output_dir, hmagma_counts_name)
write.csv(comparison_df, file = comparison_out_path, row.names = FALSE)

# Print Summary Statistics
total_added_unfiltered <- sum(comparison_df$added_snps_unfiltered)
total_added_filtered <- sum(comparison_df$added_snps_filtered)
total_removed_noise <- sum(comparison_df$noise_snps_removed)
pct_noise_removed <- (total_removed_noise / total_added_unfiltered) * 100

cat("\n=== 3-Way SNP Assignment Comparison Summary ===\n")
cat("Total SNPs added by original Hi-C (unfiltered): ", total_added_unfiltered, "\n")
cat("Total SNPs added by Hi-C + ATAC-seq (filtered): ", total_added_filtered, "\n")
cat("Noise SNPs filtered out by ATAC-seq peaks:       ", total_removed_noise,
  " (", round(pct_noise_removed, 2), "% reduction)\n",
  sep = ""
)
cat("Mean SNPs added per gene (ATAC filtered):        ", round(mean(comparison_df$added_snps_filtered), 2), "\n")
cat("Max SNPs added to a single gene (ATAC filtered): ", max(comparison_df$added_snps_filtered), "\n")
cat("Comparison table saved to:                        ", comparison_out_path, "\n\n")
