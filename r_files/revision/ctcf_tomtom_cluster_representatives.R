#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3L) {
  stop(
    paste(
      "Usage: Rscript ctcf_tomtom_cluster_representatives.R",
      "<tomtom.tsv> <input.meme> <output_dir>"
    ),
    call. = FALSE
  )
}

tomtom.file <- args[[1]]
motif.file <- args[[2]]
output.dir <- args[[3]]

dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

# Use the strong-similarity threshold used by XSTREME motif grouping.
tomtom.evalue.threshold <- 0.05

# Read all pairwise Tomtom matches and preserve the complete motif universe.
df.tomtom <- read_tsv(
  tomtom.file,
  comment = "#",
  show_col_types = FALSE
)

motif.ids <- sort(unique(c(df.tomtom$Query_ID, df.tomtom$Target_ID)))

if (length(motif.ids) != 100L) {
  warning("Expected 100 historical motifs but found ", length(motif.ids), ".")
}

# Build the directed strong-similarity table used for deterministic seed groups.
df.strong.similarity <- df.tomtom %>%
  filter(`E-value` <= tomtom.evalue.threshold) %>%
  select(
    query_id = Query_ID,
    target_id = Target_ID,
    p_value = `p-value`,
    e_value = `E-value`,
    q_value = `q-value`,
    overlap = Overlap,
    orientation = Orientation
  )

# Greedily select the seed covering the most unassigned motifs; use aggregate
# Tomtom significance and motif ID as deterministic tie breakers.
unassigned.ids <- motif.ids
cluster.members <- list()
cluster.representatives <- character()
cluster.index <- 0L

while (length(unassigned.ids) > 0L) {
  df.candidate.seed <- map_dfr(unassigned.ids, function(seed.id) {
    df.seed.hits <- df.strong.similarity %>%
      filter(query_id == seed.id, target_id %in% unassigned.ids)

    tibble(
      representative_motif_id = seed.id,
      n_unassigned_motifs_covered = n_distinct(df.seed.hits$target_id),
      aggregate_similarity_strength = sum(
        -log10(pmax(df.seed.hits$e_value, 1e-300))
      )
    )
  }) %>%
    arrange(
      desc(n_unassigned_motifs_covered),
      desc(aggregate_similarity_strength),
      representative_motif_id
    )

  representative.id <- df.candidate.seed$representative_motif_id[[1]]
  member.ids <- df.strong.similarity %>%
    filter(
      query_id == representative.id,
      target_id %in% unassigned.ids
    ) %>%
    pull(target_id) %>%
    unique()

  member.ids <- sort(unique(c(representative.id, member.ids)))
  cluster.index <- cluster.index + 1L
  cluster.representatives[[cluster.index]] <- representative.id
  cluster.members[[cluster.index]] <- member.ids
  unassigned.ids <- setdiff(unassigned.ids, member.ids)
}

# Export one row per original motif with its cluster and representative motif.
df.cluster.mapping <- map2_dfr(
  cluster.members,
  seq_along(cluster.members),
  function(member.ids, cluster.id) {
    tibble(
      cluster_id = cluster.id,
      representative_motif_id = cluster.representatives[[cluster.id]],
      motif_id = member.ids,
      cluster_size = length(member.ids),
      tomtom_evalue_threshold = tomtom.evalue.threshold,
      clustering_method = "deterministic_greedy_seed_maximum_coverage"
    )
  }
) %>%
  arrange(cluster_id, motif_id)

df.representatives <- df.cluster.mapping %>%
  distinct(
    cluster_id,
    representative_motif_id,
    cluster_size,
    tomtom_evalue_threshold,
    clustering_method
  ) %>%
  arrange(cluster_id)

# Extract complete MEME blocks for the selected representative motifs.
motif.lines <- readLines(motif.file, warn = FALSE)
motif.starts <- which(str_detect(motif.lines, "^MOTIF\\s+"))

if (length(motif.starts) != length(motif.ids)) {
  stop(
    "The MEME motif count does not match the Tomtom motif universe.",
    call. = FALSE
  )
}

motif.ends <- c(motif.starts[-1L] - 1L, length(motif.lines))
motif.block.ids <- str_match(
  motif.lines[motif.starts],
  "^MOTIF\\s+(\\S+)"
)[, 2]

representative.blocks <- map(
  df.representatives$representative_motif_id,
  function(representative.id) {
    block.index <- match(representative.id, motif.block.ids)

    if (is.na(block.index)) {
      stop("Representative motif missing from MEME file: ", representative.id)
    }

    motif.lines[motif.starts[[block.index]]:motif.ends[[block.index]]]
  }
)

representative.meme.lines <- c(
  motif.lines[seq_len(motif.starts[[1]] - 1L)],
  unlist(map(representative.blocks, ~ c(.x, "")), use.names = FALSE)
)

write_tsv(
  df.cluster.mapping,
  file.path(output.dir, "historical_100_tomtom_cluster_mapping.tsv")
)
write_tsv(
  df.representatives,
  file.path(output.dir, "historical_100_tomtom_representatives.tsv")
)
writeLines(
  representative.meme.lines,
  file.path(output.dir, "historical_100_tomtom_5_representatives.meme")
)

print(df.representatives)
