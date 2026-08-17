# lintr: disable
library("data.table")
library("tidyverse")

options(tibble.width = Inf)
options(tibble.print_max = Inf)
options(tibble.max_extra_cols = Inf)
options(scipen = 999)

########################
# 1. Input files and output directories
########################

# Locate this script when it is run line by line or with Rscript.
script.argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script.dir <- if (length(script.argument)) {
  dirname(normalizePath(sub("^--file=", "", script.argument[1])))
} else {
  candidate.dirs <- c(
    getwd(),
    file.path(getwd(), "r_files/revision/hrdp_genotype_diversity")
  )
  candidate.dirs <- candidate.dirs[
    file.exists(file.path(candidate.dirs, "hrdp_genotype_diversity.R"))
  ]
  if (!length(candidate.dirs)) {
    stop(
      "Set the working directory to the enhancer repository root or ",
      "the hrdp_genotype_diversity directory before running this file."
    )
  }
  normalizePath(candidate.dirs[1])
}
nearby.archive <- file.path(
  normalizePath(file.path(script.dir, "../../.."), mustWork = FALSE),
  "genotype_hdrp/hrdp_genotype_by_chr.tgz"
)
drive.archive <- path.expand(paste0(
  "~/Library/CloudStorage/GoogleDrive-wellclouder@gmail.com/My Drive/",
  "research/enhancer/genotype_hdrp/hrdp_genotype_by_chr.tgz"
))
archive.file <- Sys.getenv(
  "HRDP_GENOTYPE_ARCHIVE",
  unset = if (file.exists(nearby.archive)) nearby.archive else drive.archive
)
local.plink2 <- file.path(script.dir, ".tools/plink2/bin/plink2")
plink2.bin <- Sys.getenv(
  "PLINK2_BIN",
  unset = if (file.exists(local.plink2)) local.plink2 else Sys.which("plink2")
)
output.dir <- file.path(script.dir, "results/plink2_primary")

if (!file.exists(archive.file)) stop("Missing HRDP genotype archive: ", archive.file)
if (!nzchar(plink2.bin) || !file.exists(plink2.bin)) {
  machine <- tolower(Sys.info()[["machine"]])
  if (Sys.info()[["sysname"]] != "Darwin" || !machine %in% c("arm64", "aarch64")) {
    stop("PLINK 2.0 was not found. Set PLINK2_BIN to an executable for this platform.")
  }

  plink2.url <- paste0(
    "https://s3.amazonaws.com/plink2-assets/alpha7/",
    "plink2_mac_arm64_20260808.zip"
  )
  plink2.archive <- tempfile(fileext = ".zip")
  dir.create(dirname(local.plink2), recursive = TRUE, showWarnings = FALSE)
  download.file(plink2.url, plink2.archive, mode = "wb")
  unzip(plink2.archive, exdir = dirname(local.plink2))
  unlink(plink2.archive)
  Sys.chmod(local.plink2, mode = "0755")
  plink2.bin <- local.plink2
}
if (file.access(plink2.bin, mode = 1L) != 0L) stop("PLINK 2.0 is not executable: ", plink2.bin)
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

########################
# 2. Sample mapping and complete autosomal SNPs
########################

df.sample.map <- tribble(
  ~sample_id, ~strain, ~vcf_sample, ~mapping_status,
  "592BB", "SHR/OlaIpcv", "SHR_OlaIpcv", "exact strain-name match",
  "607", "HXB10", "HXB10_Ipcv", "exact strain-name match",
  "74AA", "F344/Stm", "F344_StmMcwi", "exact strain-name match",
  "A2DB", "LE/Stm", "LE_StmMcwi", "exact strain-name match",
  "D765A", "BXH6", "BXH6_Cub", "exact strain-name match",
  "DA08A", "HXB2", "HXB2_Ipcv", "exact strain-name match",
  "DA21A", "SHR/OlaIpcvxBN/NHsdMcwi", "SHR_BN_F1",
  "confirmed sample-to-VCF mapping",
  "DA68A", "HXB31", "HXB31_Ipcv", "exact strain-name match",
  "DBA9A", "HXB23", "HXB23_Ipcv", "exact strain-name match",
  "DE8BA", "BN-Lx", "BN-Lx_Cub", "exact strain-name match"
)

# Convert a biallelic GT field to alternate-allele dosage (0, 1, or 2).
gt.to.dosage <- function(x) {
  gt <- chartr("|", "/", sub(":.*$", "", x))
  recode(
    gt,
    "0/0" = 0,
    "0/1" = 1,
    "1/0" = 1,
    "1/1" = 2,
    .default = NA_real_
  )
}

variant.parts <- vector("list", 20L)
dosage.parts <- vector("list", 20L)

for (chromosome.index in seq_len(20L)) {
  chromosome <- paste0("chr", chromosome.index)
  archive.member <- paste0(
    "hrdp_118strains_plus_F1_genotype_", chromosome, ".gvcf.gz"
  )
  command <- sprintf(
    "tar -xOzf %s %s | gzip -dc | awk '!/^##/'",
    shQuote(archive.file),
    shQuote(archive.member)
  )
  df.vcf <- fread(
    cmd = command,
    select = c("#CHROM", "POS", "REF", "ALT", df.sample.map$vcf_sample),
    data.table = FALSE,
    showProgress = FALSE
  ) %>%
    filter(nchar(REF) == 1L, nchar(ALT) == 1L, !str_detect(ALT, fixed(",")))

  dosage <- vapply(
    df.vcf[df.sample.map$vcf_sample],
    gt.to.dosage,
    numeric(nrow(df.vcf))
  )
  informative <- rowSums(!is.na(dosage)) == ncol(dosage) &
    apply(dosage, 1L, function(x) length(unique(x)) > 1L)

  variant.parts[[chromosome.index]] <- df.vcf[informative, ] %>%
    transmute(
      chromosome_number = chromosome.index,
      position = POS,
      ref = REF,
      alt = ALT,
      variant_id = paste(chromosome, POS, REF, ALT, sep = ":")
    )
  dosage.parts[[chromosome.index]] <- dosage[informative, , drop = FALSE]
}

df.variants <- bind_rows(variant.parts)
dosage.matrix <- do.call(rbind, dosage.parts)
colnames(dosage.matrix) <- df.sample.map$sample_id
if (anyDuplicated(df.variants$variant_id)) {
  df.variants$variant_id <- make.unique(df.variants$variant_id)
}

########################
# 3. PLINK PED/MAP input
########################

# Export complete genotypes as PED/MAP so PLINK uses the same primary SNP set.
plink.input <- file.path(output.dir, "hrdp_ten_libraries_complete")
fwrite(
  df.variants %>%
    transmute(chromosome_number, variant_id, cm = 0, position),
  paste0(plink.input, ".map"),
  sep = "\t",
  col.names = FALSE
)
ped.rows <- map(seq_len(nrow(df.sample.map)), function(sample.index) {
  dosage <- dosage.matrix[, sample.index]
  allele.one <- ifelse(dosage == 2, df.variants$alt, df.variants$ref)
  allele.two <- ifelse(dosage == 0, df.variants$ref, df.variants$alt)
  c(
    df.sample.map$sample_id[sample.index],
    df.sample.map$sample_id[sample.index],
    0, 0, 1, -9,
    as.vector(rbind(allele.one, allele.two))
  )
})
fwrite(
  as.data.table(do.call(rbind, ped.rows)),
  paste0(plink.input, ".ped"),
  sep = " ",
  col.names = FALSE
)

########################
# 4. LD pruning, IBS distance, and PCA with PLINK 2.0
########################

# Run one PLINK 2.0 command and stop if it fails.
run.plink2 <- function(arguments) {
  status <- system2(plink2.bin, arguments)
  if (status != 0L) stop("PLINK 2.0 failed: ", paste(arguments, collapse = " "))
}

binary.prefix <- file.path(output.dir, "hrdp_complete")
prune.prefix <- file.path(output.dir, "hrdp_complete_ld_pruned")
distance.prefix <- file.path(output.dir, "hrdp_complete_ld_pruned_ibs")
pca.prefix <- file.path(output.dir, "hrdp_complete_ld_pruned_pca")
chromosome.arguments <- c("--chr-set", "20", "no-x", "no-y", "no-xy", "no-mt")

run.plink2(c(
  "--pedmap", plink.input, chromosome.arguments,
  "--make-pgen", "--out", binary.prefix
))
run.plink2(c(
  "--pfile", binary.prefix, chromosome.arguments,
  "--indep-pairwise", "50", "5", "0.2",
  "--indep-order", "1", "--bad-ld", "--out", prune.prefix
))
run.plink2(c(
  "--pfile", binary.prefix, chromosome.arguments,
  "--extract", paste0(prune.prefix, ".prune.in"),
  "--make-king-table", "cols=fid,id,nsnp,ibs",
  "--out", distance.prefix
))
run.plink2(c(
  "--pfile", binary.prefix, chromosome.arguments,
  "--extract", paste0(prune.prefix, ".prune.in"),
  "--pca", "5", "--bad-freqs", "--out", pca.prefix
))

########################
# 5. PLINK distance and PCA tables
########################

# Reconstruct a symmetric distance matrix from PLINK 2.0's pairwise IBS table.
df.plink.pairs <- fread(paste0(distance.prefix, ".kin0"), data.table = FALSE)
plink.ids <- df.sample.map$sample_id
plink.distance <- matrix(
  0,
  nrow = length(plink.ids),
  ncol = length(plink.ids),
  dimnames = list(plink.ids, plink.ids)
)
for (pair.index in seq_len(nrow(df.plink.pairs))) {
  sample.one <- df.plink.pairs$IID1[pair.index]
  sample.two <- df.plink.pairs$IID2[pair.index]
  plink.distance[sample.one, sample.two] <- df.plink.pairs$IBS[pair.index]
  plink.distance[sample.two, sample.one] <- df.plink.pairs$IBS[pair.index]
}

df.distance.long <- as.data.frame(as.table(plink.distance)) %>%
  as_tibble() %>%
  rename(sample1 = Var1, sample2 = Var2, plink2_ibs_distance = Freq) %>%
  left_join(df.sample.map %>% select(sample1 = sample_id, strain1 = strain), by = "sample1") %>%
  left_join(df.sample.map %>% select(sample2 = sample_id, strain2 = strain), by = "sample2")

df.pca.eigenvec <- fread(
  paste0(pca.prefix, ".eigenvec"),
  header = TRUE,
  data.table = FALSE
) %>%
  rename(fid = `#FID`, sample_id = IID) %>%
  left_join(df.sample.map, by = "sample_id")
pca.eigenval <- scan(paste0(pca.prefix, ".eigenval"), quiet = TRUE)

df.summary <- tibble(
  plink_version = system2(plink2.bin, "--version", stdout = TRUE)[1],
  n_libraries = nrow(df.sample.map),
  n_complete_variable_snps = nrow(df.variants),
  n_ld_pruned_snps = length(readLines(paste0(prune.prefix, ".prune.in"))),
  minimum_pairwise_ibs_distance = min(plink.distance[upper.tri(plink.distance)]),
  maximum_pairwise_ibs_distance = max(plink.distance[upper.tri(plink.distance)]),
  pruning_parameters = paste(
    "PLINK 2.0 --indep-pairwise 50 5 0.2",
    "--indep-order 1 --bad-ld"
  ),
  distance_method = "PLINK 2.0 --make-king-table cols=fid,id,nsnp,ibs",
  pca_method = "PLINK 2.0 --pca 5 --bad-freqs"
)

########################
# 6. Output tables and figures
########################

write_tsv(df.summary, file.path(output.dir, "hrdp_plink2_primary_summary.tsv"))
write_tsv(df.sample.map, file.path(output.dir, "hrdp_plink2_sample_mapping.tsv"))
write_tsv(df.distance.long, file.path(output.dir, "hrdp_plink2_ibs_distance.tsv"))
write_tsv(df.pca.eigenvec, file.path(output.dir, "hrdp_plink2_ld_pruned_pca.tsv"))
write_tsv(
  tibble(
    principal_component = paste0("PC", seq_along(pca.eigenval)),
    eigenvalue = pca.eigenval
  ),
  file.path(output.dir, "hrdp_plink2_ld_pruned_pca_eigenvalues.tsv")
)
write.table(
  plink.distance,
  file.path(output.dir, "hrdp_plink2_ibs_distance_matrix.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

cluster.fit <- hclust(as.dist(plink.distance))
strain.order <- df.sample.map$strain[match(
  cluster.fit$labels[cluster.fit$order],
  df.sample.map$sample_id
)]
plot.heatmap <- df.distance.long %>%
  mutate(
    strain1 = factor(strain1, levels = strain.order),
    strain2 = factor(strain2, levels = rev(strain.order))
  ) %>%
  ggplot(aes(strain1, strain2, fill = plink2_ibs_distance)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_viridis_c(option = "C", direction = -1) +
  coord_equal() +
  labs(x = NULL, y = NULL, fill = "LD-pruned\nPLINK 2 IBS") +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

plot.pca <- ggplot(df.pca.eigenvec, aes(PC1, PC2, label = strain)) +
  geom_point(size = 2.5, color = "#1F5A94") +
  ggrepel::geom_text_repel(size = 3, max.overlaps = Inf) +
  labs(x = "PLINK 2 PC1", y = "PLINK 2 PC2") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank())

for (extension in c("pdf", "png")) {
  ggsave(
    file.path(output.dir, paste0("hrdp_plink2_ibs_heatmap.", extension)),
    plot.heatmap,
    width = 6.5,
    height = 5.5,
    dpi = 300
  )
  ggsave(
    file.path(output.dir, paste0("hrdp_plink2_ld_pruned_pca.", extension)),
    plot.pca,
    width = 6.5,
    height = 4.8,
    dpi = 300
  )
}

writeLines(capture.output(sessionInfo()), file.path(output.dir, "hrdp_plink2_session_info.txt"))
message("Completed primary PLINK 2.0 analysis: ", output.dir)
