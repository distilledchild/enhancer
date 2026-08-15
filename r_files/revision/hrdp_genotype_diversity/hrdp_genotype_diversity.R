# Quantify genome-wide SNP relationships among the ten Hi-C libraries.

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
})

options(scipen = 999)

# Keep results beside this script and use the matching project archive first.
script.argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script.dir <- if (length(script.argument)) {
  dirname(normalizePath(sub("^--file=", "", script.argument[1])))
} else getwd()
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
output.dir <- Sys.getenv(
  "HRDP_GENOTYPE_OUTPUT_DIR",
  unset = file.path(script.dir, "results")
)

if (!file.exists(archive.file)) stop("Missing HRDP genotype archive: ", archive.file)
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

sample.map <- tribble(
  ~sample_id, ~strain, ~vcf_sample, ~mapping_status,
  "592BB", "SHR/OlaIpcv", "SHR_OlaIpcv", "exact strain-name match",
  "607", "HXB10", "HXB10_Ipcv", "exact strain-name match",
  "74AA", "F344/Stm", "F344_StmMcwi", "exact strain-name match",
  "A2DB", "LE/Stm", "LE_StmMcwi", "exact strain-name match",
  "D765A", "BXH6", "BXH6_Cub", "exact strain-name match",
  "DA08A", "HXB2", "HXB2_Ipcv", "exact strain-name match",
  "DA21A", "SHR/OlaIpcvxBN/NHsdMcwi", "SHR_BN_F1",
  "inferred F1 match; confirm with genotype-data provider",
  "DA68A", "HXB31", "HXB31_Ipcv", "exact strain-name match",
  "DBA9A", "HXB23", "HXB23_Ipcv", "exact strain-name match",
  "DE8BA", "BN-Lx", "BN-Lx_Cub", "exact strain-name match"
)

# Convert the biallelic GT subfield to alternate-allele dosage (0, 1, or 2).
gt_to_dosage <- function(x) {
  gt <- gsub("\\|", "/", sub(":.*$", "", x))
  recode(
    gt,
    "0/0" = 0,
    "0/1" = 1,
    "1/0" = 1,
    "1/1" = 2,
    .default = NA_real_
  )
}

chromosomes <- paste0("chr", 1:20)
genotype.parts <- vector("list", length(chromosomes))
variant.qc <- vector("list", length(chromosomes))

for (list.index in seq_along(chromosomes)) {
  chromosome.i <- chromosomes[[list.index]]
  archive.member <- paste0(
    "hrdp_118strains_plus_F1_genotype_", chromosome.i, ".gvcf.gz"
  )
  command <- sprintf(
    "tar -xOzf %s %s | gzip -dc | awk '!/^##/'",
    shQuote(archive.file),
    shQuote(archive.member)
  )
  df.vcf <- data.table::fread(
    cmd = command,
    select = c("#CHROM", "POS", "REF", "ALT", sample.map$vcf_sample),
    data.table = FALSE,
    showProgress = FALSE
  ) %>%
    filter(nchar(REF) == 1L, nchar(ALT) == 1L, !str_detect(ALT, fixed(",")))

  dosage <- vapply(
    df.vcf[sample.map$vcf_sample],
    gt_to_dosage,
    numeric(nrow(df.vcf))
  )
  colnames(dosage) <- sample.map$strain
  informative <- rowSums(!is.na(dosage)) >= 8L &
    apply(dosage, 1L, function(x) length(unique(x[!is.na(x)])) > 1L)

  genotype.parts[[list.index]] <- dosage[informative, , drop = FALSE]
  variant.qc[[list.index]] <- tibble(
    chromosome = chromosome.i,
    n_biallelic_snp_records = nrow(dosage),
    n_informative_for_ten_libraries = sum(informative)
  )
  message(chromosome.i, ": ", sum(informative), " informative SNPs retained.")
}

genotype.matrix <- do.call(rbind, genotype.parts)
variant.qc <- bind_rows(variant.qc)
if (!nrow(genotype.matrix)) stop("No informative SNPs remained for the ten libraries.")

# Mean absolute dosage difference / 2 is an allele-sharing distance: zero for
# identical genotypes and one for opposite homozygotes at every SNP.
n.strains <- ncol(genotype.matrix)
distance.matrix <- matrix(
  0,
  nrow = n.strains,
  ncol = n.strains,
  dimnames = list(colnames(genotype.matrix), colnames(genotype.matrix))
)
n.compared.matrix <- distance.matrix
for (i in seq_len(n.strains)) {
  for (j in i:n.strains) {
    comparable <- !is.na(genotype.matrix[, i]) & !is.na(genotype.matrix[, j])
    n.compared.matrix[i, j] <- n.compared.matrix[j, i] <- sum(comparable)
    distance.matrix[i, j] <- distance.matrix[j, i] <- mean(
      abs(genotype.matrix[comparable, i] - genotype.matrix[comparable, j]) / 2
    )
  }
}

# Mean-impute the few missing dosages and perform sample-level PCA.
pca.matrix <- genotype.matrix
missing.index <- which(is.na(pca.matrix), arr.ind = TRUE)
if (nrow(missing.index))
  pca.matrix[missing.index] <- rowMeans(pca.matrix, na.rm = TRUE)[missing.index[, 1]]
pca.fit <- prcomp(t(pca.matrix), center = TRUE, scale. = FALSE)
pca.variance <- 100 * pca.fit$sdev^2 / sum(pca.fit$sdev^2)
df.pca <- as_tibble(pca.fit$x[, 1:3, drop = FALSE], rownames = "strain") %>%
  left_join(sample.map, by = "strain")

df.distance.long <- as.data.frame(as.table(distance.matrix)) %>%
  as_tibble() %>%
  rename(strain1 = Var1, strain2 = Var2, allele_sharing_distance = Freq) %>%
  mutate(
    n_snps_compared = n.compared.matrix[cbind(
      match(strain1, rownames(n.compared.matrix)),
      match(strain2, colnames(n.compared.matrix))
    )]
  )

strain.clustering <- hclust(as.dist(distance.matrix))
strain.order <- strain.clustering$labels[strain.clustering$order]
plot.heatmap <- df.distance.long %>%
  mutate(
    strain1 = factor(strain1, levels = strain.order),
    strain2 = factor(strain2, levels = rev(strain.order))
  ) %>%
  ggplot(aes(strain1, strain2, fill = allele_sharing_distance)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_viridis_c(option = "C", direction = -1) +
  coord_equal() +
  labs(x = NULL, y = NULL, fill = "SNP distance") +
  theme_bw(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

plot.pca <- ggplot(df.pca, aes(PC1, PC2, label = strain)) +
  geom_point(size = 2.5, color = "#1F5A94") +
  ggrepel::geom_text_repel(size = 3, max.overlaps = Inf) +
  labs(
    x = sprintf("PC1 (%.1f%%)", pca.variance[1]),
    y = sprintf("PC2 (%.1f%%)", pca.variance[2])
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank())

output.tables <- list(
  hrdp_ten_library_sample_mapping = sample.map,
  hrdp_genotype_variant_qc = variant.qc,
  hrdp_pairwise_snp_distance = df.distance.long,
  hrdp_genotype_pca_coordinates = df.pca,
  hrdp_genotype_pca_variance = tibble(
    principal_component = paste0("PC", seq_along(pca.variance)),
    variance_explained_percent = pca.variance
  )
)
walk2(names(output.tables), output.tables, ~ write_tsv(
  .y, file.path(output.dir, paste0(.x, ".tsv"))
))
write.table(
  distance.matrix,
  file.path(output.dir, "hrdp_pairwise_snp_distance_matrix.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
ggsave(
  file.path(output.dir, "hrdp_pairwise_snp_distance_heatmap.pdf"),
  plot.heatmap,
  width = 6.5,
  height = 5.5
)
ggsave(
  file.path(output.dir, "hrdp_genotype_pca.pdf"),
  plot.pca,
  width = 6.5,
  height = 4.8
)

write_tsv(
  tibble(
    analysis = "Genome-wide SNP relationships among ten Hi-C libraries",
    reference_build = "mRatBN7.2/rn7",
    n_libraries = nrow(sample.map),
    n_autosomes = length(chromosomes),
    n_informative_snps = nrow(genotype.matrix),
    distance_definition = paste0(
      "Mean absolute diploid alternate-allele dosage difference divided by 2; ",
      "computed over pairwise nonmissing informative biallelic SNPs."
    ),
    source_note = paste0(
      "The supplied HRDP archive was previously filtered and merged with ",
      "missing-to-reference handling; results describe relationships in this ",
      "provided panel rather than de novo variant calling."
    )
  ),
  file.path(output.dir, "hrdp_genotype_analysis_metadata.tsv")
)
writeLines(
  sub("[[:blank:]]+$", "", capture.output(sessionInfo())),
  file.path(output.dir, "hrdp_genotype_session_info.txt")
)

message(
  "Completed genotype analysis using ", nrow(genotype.matrix),
  " informative SNPs. Results: ", output.dir
)
