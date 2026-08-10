# lintr: disable
library("tidyverse")
library("GenomicRanges")

options(scipen = 999)

arguments <- commandArgs(trailingOnly = TRUE)

if (length(arguments) != 4L) {
  stop(
    paste0(
      "Usage: Rscript ctcf_fimo_loop_annotation.R <loop_resource.tsv> ",
      "<predicted_loci.tsv> <output_directory> <analysis_label>"
    ),
    call. = FALSE
  )
}

loop.file <- normalizePath(arguments[[1L]], mustWork = TRUE)
ctcf.locus.file <- normalizePath(arguments[[2L]], mustWork = TRUE)
output.dir <- arguments[[3L]]
analysis.label <- arguments[[4L]]
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

# Read the pooled loop resource and one nonredundant predicted-motif locus set.
df.loop <- read_tsv(loop.file, show_col_types = FALSE)
df.ctcf <- read_tsv(ctcf.locus.file, show_col_types = FALSE)

required.loop.columns <- c(
  "loop_id", "chr1", "start1", "end1", "chr2", "start2", "end2",
  "resolution", "revised_putative_regulatory_support"
)
required.ctcf.columns <- c(
  "ctcf_predicted_locus_id", "chr", "locus_start", "locus_end",
  "representative_strand", "representative_p_value", "representative_q_value"
)

if (
  !all(required.loop.columns %in% names(df.loop)) ||
    !all(required.ctcf.columns %in% names(df.ctcf))
) {
  stop("Loop or predicted CTCF locus table lacks required columns.", call. = FALSE)
}

if (
  any(df.loop$chr1 != df.loop$chr2) ||
    any(df.loop$start1 > df.loop$start2) ||
    any(!df.ctcf$representative_strand %in% c("+", "-"))
) {
  stop("Expected ordered intrachromosomal loops and valid motif strands.", call. = FALSE)
}

# Represent both loop anchors as one table while retaining their genomic order.
df.anchor <- bind_rows(
  df.loop %>%
    transmute(
      loop_id,
      resolution,
      anchor_side = "anchor1",
      chr = chr1,
      start = start1,
      end = end1
    ),
  df.loop %>%
    transmute(
      loop_id,
      resolution,
      anchor_side = "anchor2",
      chr = chr2,
      start = start2,
      end = end2
    )
) %>%
  mutate(anchor_index = row_number(), .before = 1)

gr.anchor <- GRanges(
  seqnames = df.anchor$chr,
  ranges = IRanges(start = df.anchor$start, end = df.anchor$end),
  anchor_index = df.anchor$anchor_index
)
gr.ctcf <- GRanges(
  seqnames = df.ctcf$chr,
  ranges = IRanges(start = df.ctcf$locus_start, end = df.ctcf$locus_end),
  strand = df.ctcf$representative_strand,
  ctcf_index = seq_len(nrow(df.ctcf))
)

# Count predicted loci by strand at each anchor; strand is then used only for a
# sequence-based convergent-orientation flag, not as evidence of CTCF binding.
df.overlap <- as_tibble(as.data.frame(findOverlaps(
  gr.anchor,
  gr.ctcf,
  ignore.strand = TRUE
))) %>%
  transmute(
    anchor_index = queryHits,
    ctcf_index = subjectHits
  ) %>%
  left_join(
    df.ctcf %>%
      mutate(ctcf_index = row_number(), .before = 1) %>%
      select(ctcf_index, representative_strand),
    by = "ctcf_index"
  )

df.anchor.count <- df.anchor %>%
  left_join(
    df.overlap %>%
      group_by(anchor_index) %>%
      summarise(
        n_predicted_ctcf_loci = n_distinct(ctcf_index),
        n_predicted_ctcf_plus = n_distinct(
          ctcf_index[representative_strand == "+"]
        ),
        n_predicted_ctcf_minus = n_distinct(
          ctcf_index[representative_strand == "-"]
        ),
        .groups = "drop"
      ),
    by = "anchor_index"
  ) %>%
  mutate(
    across(
      c(
        n_predicted_ctcf_loci,
        n_predicted_ctcf_plus,
        n_predicted_ctcf_minus
      ),
      ~replace_na(.x, 0L)
    )
  ) %>%
  select(-anchor_index, -chr, -start, -end) %>%
  pivot_wider(
    names_from = anchor_side,
    values_from = c(
      n_predicted_ctcf_loci,
      n_predicted_ctcf_plus,
      n_predicted_ctcf_minus
    ),
    names_glue = "{.value}_{anchor_side}"
  )

df.loop.annotation <- df.loop %>%
  select(
    loop_id,
    resolution,
    revised_putative_regulatory_support
  ) %>%
  left_join(df.anchor.count, by = c("loop_id", "resolution")) %>%
  mutate(
    predicted_ctcf_any_anchor1 = n_predicted_ctcf_loci_anchor1 > 0L,
    predicted_ctcf_any_anchor2 = n_predicted_ctcf_loci_anchor2 > 0L,
    predicted_ctcf_both_anchors = (
      predicted_ctcf_any_anchor1 & predicted_ctcf_any_anchor2
    ),
    predicted_ctcf_convergent_pair = (
      n_predicted_ctcf_plus_anchor1 > 0L &
        n_predicted_ctcf_minus_anchor2 > 0L
    )
  ) %>%
  arrange(resolution, loop_id)

# Summarise the full pooled resource and the putative-regulatory subset without
# using predicted CTCF motifs as a loop-retention criterion.
df.summary.input <- bind_rows(
  df.loop.annotation %>% mutate(loop_subset = "all_pooled_loops"),
  df.loop.annotation %>%
    filter(revised_putative_regulatory_support) %>%
    mutate(loop_subset = "putative_regulatory_loops")
)

df.summary.input <- bind_rows(
  df.summary.input,
  df.summary.input %>% mutate(resolution = "ALL")
)

df.summary <- df.summary.input %>%
  group_by(loop_subset, resolution) %>%
  summarise(
    n_loops = n(),
    n_any_anchor1 = sum(predicted_ctcf_any_anchor1),
    pct_any_anchor1 = round(100 * mean(predicted_ctcf_any_anchor1), 1),
    n_any_anchor2 = sum(predicted_ctcf_any_anchor2),
    pct_any_anchor2 = round(100 * mean(predicted_ctcf_any_anchor2), 1),
    n_both_anchors = sum(predicted_ctcf_both_anchors),
    pct_both_anchors = round(100 * mean(predicted_ctcf_both_anchors), 1),
    n_convergent_pairs = sum(predicted_ctcf_convergent_pair),
    pct_convergent_pairs = round(100 * mean(predicted_ctcf_convergent_pair), 1),
    .groups = "drop"
  ) %>%
  mutate(analysis_label = analysis.label, .before = 1) %>%
  arrange(loop_subset, resolution)

write_tsv(
  df.loop.annotation,
  file.path(output.dir, str_c(analysis.label, "_loop_annotation.tsv"))
)
write_tsv(
  df.summary,
  file.path(output.dir, str_c(analysis.label, "_loop_summary.tsv"))
)

print(df.summary)
