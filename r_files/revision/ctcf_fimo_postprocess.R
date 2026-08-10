# lintr: disable
library("tidyverse")
library("GenomicRanges")

options(scipen = 999)

arguments <- commandArgs(trailingOnly = TRUE)

if (length(arguments) != 2L) {
  stop(
    "Usage: Rscript ctcf_fimo_postprocess.R <fimo.tsv> <output_directory>",
    call. = FALSE
  )
}

fimo.file <- normalizePath(arguments[[1L]], mustWork = TRUE)
output.dir <- arguments[[2L]]
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

# Read the complete FIMO table while retaining motif statistics and strand.
df.fimo <- read_tsv(
  fimo.file,
  comment = "#",
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
) %>%
  filter(!is.na(motif_id), !is.na(sequence_name)) %>%
  transmute(
    motif_id,
    motif_alt_id,
    chr = sequence_name,
    start = as.integer(start),
    end = as.integer(stop),
    strand,
    score = as.numeric(score),
    p_value = as.numeric(`p-value`),
    q_value = as.numeric(`q-value`),
    matched_sequence
  )

if (
  anyNA(df.fimo[, c("chr", "start", "end", "strand", "p_value", "q_value")]) ||
    any(df.fimo$start < 1L) ||
    any(df.fimo$end < df.fimo$start) ||
    any(!df.fimo$strand %in% c("+", "-"))
) {
  stop("FIMO output contains invalid coordinates, statistics, or strand.", call. = FALSE)
}

# Collapse overlapping shifted FIMO windows into loci and retain the strongest
# window as the representative hit without discarding orientation ambiguity.
collapse_fimo_hits <- function(df, q.threshold, output.prefix) {
  df.filtered <- df %>%
    filter(q_value <= q.threshold) %>%
    arrange(chr, start, end, strand, q_value, p_value)

  gr.hit <- GRanges(
    seqnames = df.filtered$chr,
    ranges = IRanges(start = df.filtered$start, end = df.filtered$end),
    strand = df.filtered$strand,
    hit_index = seq_len(nrow(df.filtered))
  )

  # min.gapwidth=0 merges overlapping windows but does not merge merely
  # adjacent motif occurrences.
  gr.locus <- reduce(
    gr.hit,
    ignore.strand = TRUE,
    min.gapwidth = 0L
  )

  df.hit.to.locus <- as_tibble(as.data.frame(findOverlaps(
    gr.hit,
    gr.locus,
    ignore.strand = TRUE
  ))) %>%
    transmute(
      hit_index = queryHits,
      locus_index = subjectHits
    )

  if (
    nrow(df.hit.to.locus) != nrow(df.filtered) ||
      anyDuplicated(df.hit.to.locus$hit_index)
  ) {
    stop("Each filtered FIMO hit must map to exactly one collapsed locus.", call. = FALSE)
  }

  df.locus.interval <- tibble(
    locus_index = seq_along(gr.locus),
    chr = as.character(seqnames(gr.locus)),
    locus_start = start(gr.locus),
    locus_end = end(gr.locus)
  )

  df.locus <- df.filtered %>%
    mutate(hit_index = row_number(), .before = 1) %>%
    inner_join(df.hit.to.locus, by = "hit_index") %>%
    group_by(locus_index) %>%
    arrange(q_value, p_value, desc(score), start, end, strand, .by_group = TRUE) %>%
    mutate(
      n_overlapping_windows = n(),
      n_plus_windows = sum(strand == "+"),
      n_minus_windows = sum(strand == "-"),
      strand_support_class = case_when(
        n_plus_windows > 0L & n_minus_windows > 0L ~ "both_strands",
        n_plus_windows > 0L ~ "plus_only",
        TRUE ~ "minus_only"
      )
    ) %>%
    dplyr::slice(1L) %>%
    ungroup() %>%
    left_join(df.locus.interval, by = c("locus_index", "chr")) %>%
    transmute(
      ctcf_predicted_locus_id = str_c(
        "CTCF_", chr, "_", locus_start, "_", locus_end
      ),
      chr,
      locus_start,
      locus_end,
      representative_hit_start = start,
      representative_hit_end = end,
      representative_strand = strand,
      motif_id,
      motif_alt_id,
      representative_score = score,
      representative_p_value = p_value,
      representative_q_value = q_value,
      n_overlapping_windows,
      n_plus_windows,
      n_minus_windows,
      strand_support_class,
      representative_matched_sequence = matched_sequence
    ) %>%
    arrange(chr, locus_start, locus_end)

  df.bed <- df.locus %>%
    transmute(
      chrom = chr,
      chromStart = locus_start - 1L,
      chromEnd = locus_end,
      name = ctcf_predicted_locus_id,
      score = pmin(
        1000L,
        as.integer(round(-10 * log10(representative_q_value)))
      ),
      strand = representative_strand
    )

  write_tsv(
    df.locus,
    file.path(output.dir, str_c(output.prefix, "_overlap_collapsed_loci.tsv"))
  )
  write_tsv(
    df.bed,
    file.path(output.dir, str_c(output.prefix, "_overlap_collapsed_loci.bed")),
    col_names = FALSE
  )

  tibble(
    q_value_threshold = q.threshold,
    n_fimo_windows = nrow(df.filtered),
    n_overlap_collapsed_loci = nrow(df.locus),
    n_loci_plus_only = sum(df.locus$strand_support_class == "plus_only"),
    n_loci_minus_only = sum(df.locus$strand_support_class == "minus_only"),
    n_loci_both_strands = sum(df.locus$strand_support_class == "both_strands"),
    median_locus_width = median(df.locus$locus_end - df.locus$locus_start + 1L),
    max_locus_width = max(df.locus$locus_end - df.locus$locus_start + 1L)
  )
}

df.summary <- bind_rows(
  collapse_fimo_hits(df.fimo, 0.05, "fimo_q0.05"),
  collapse_fimo_hits(df.fimo, 0.01, "fimo_q0.01")
)

write_tsv(
  df.summary,
  file.path(output.dir, "fimo_overlap_collapse_summary.tsv")
)

print(df.summary)
