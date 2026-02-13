#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(patchwork)
  library(rtracklayer)
  library(GenomicRanges)
  library(strawr)
  library(grid) # unit()
})

# Generates publication-style panels:
# - triangular Hi-C heatmap from .hic (via strawr)
# - CTCF binned density (bedGraph)
# - HiCCUPS loops arcs (bedpe)
# - Gene track (GTF)
#
# Output: PDF files under `out_dir`

hic_path <- "/Users/pete/UTHSC GGI Dropbox/K P/Gateway_to_Hao/hic/2023A/hic_analysis/juicer/DA68A/intact/DA68A_intact_inter_30.hic"
ctcf_bedgraph_path <- path.expand("~/dropbox/Gateway_to_Hao/enhancer/data/tracks/ctcf_density_5kb.sorted.bedGraph")
loops_bedpe_path <- path.expand("~/dropbox/Gateway_to_Hao/enhancer/data/loops/DA68A_intact_merged_loops_5k10k25k.bedpe")
gtf_path <- path.expand("~/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/ucsc_refGene.gtf")
chrom_sizes_path <- path.expand("~/dropbox/Gateway_to_Hao/enhancer/data/tracks/rn7.chrom.sizes")

out_dir <- path.expand("~/dropbox/Gateway_to_Hao/enhancer/figures/hic_panels")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(hic_path))
stopifnot(file.exists(ctcf_bedgraph_path))
stopifnot(file.exists(loops_bedpe_path))
stopifnot(file.exists(gtf_path))

if (!file.exists(chrom_sizes_path)) {
  stop("chrom.sizes not found at: ", chrom_sizes_path, "\nCreate it from df.chromosome.data or rn7 TSV first.")
}

chrom_sizes <- read_tsv(chrom_sizes_path, col_names = c("chr", "end"), show_col_types = FALSE) %>%
  mutate(chr = as.character(chr), end = as.integer(end)) %>%
  filter(!is.na(chr), !is.na(end), end > 0)

genes_of_interest <- tibble::tribble(
  ~gene,  ~chr,    ~gene_start,   ~gene_end, ~strand,
  "Foxo1","chr2",  136312168L,    136390603L, "+",
  "Gja1", "chr20",  35756007L,     35768481L, "+",
  "Spry2","chr15",  82692291L,     82697408L, "-",
  "Tgfb2","chr13",  98160075L,     98261771L, "-"
)

binsize <- 5000L
flank_bp <- 2000000L  # +/- 2Mb around gene midpoint
norm_try <- c("KR", "NONE")  # try KR first, fall back to NONE
hic_gamma <- 0.22  # smaller => darker (boosts low intensities after normalization)

read_loops <- function(path) {
  cols <- c(
    "chr1","x1","x2","chr2","y1","y2","name","score","strand1","strand2","color",
    "observed","expectedBL","expectedDonut","expectedH","expectedV",
    "fdrBL","fdrDonut","fdrH","fdrV","numCollapsed","centroid1","centroid2","radius"
  )
  read_tsv(path, comment = "#", col_names = cols, show_col_types = FALSE) %>%
    mutate(
      chr1 = as.character(chr1), chr2 = as.character(chr2),
      x1 = as.integer(x1), x2 = as.integer(x2), y1 = as.integer(y1), y2 = as.integer(y2),
      centroid1 = as.integer(centroid1), centroid2 = as.integer(centroid2)
    )
}

loops_all <- read_loops(loops_bedpe_path)
ctcf_all <- read_tsv(ctcf_bedgraph_path, col_names = c("chr","start","end","score"), show_col_types = FALSE) %>%
  mutate(chr = as.character(chr), start = as.integer(start), end = as.integer(end), score = as.numeric(score))

pack_rows <- function(df, start_col = "start", end_col = "end") {
  # Greedy interval packing into rows (minimize overlaps).
  if (nrow(df) == 0) return(df %>% mutate(row = integer()))

  df <- df %>% arrange(.data[[start_col]], .data[[end_col]])
  row_ends <- integer()
  rows <- integer(nrow(df))

  for (i in seq_len(nrow(df))) {
    s <- df[[start_col]][i]
    e <- df[[end_col]][i]
    placed <- FALSE
    if (length(row_ends) > 0) {
      for (r in seq_along(row_ends)) {
        if (s > row_ends[r]) {
          row_ends[r] <- e
          rows[i] <- r
          placed <- TRUE
          break
        }
      }
    }
    if (!placed) {
      row_ends <- c(row_ends, e)
      rows[i] <- length(row_ends)
    }
  }
  df %>% mutate(row = rows)
}

compute_triangle_coords <- function(pos1, pos2, region_start, binsize) {
  i <- (pos1 - region_start) / binsize
  j <- (pos2 - region_start) / binsize
  px <- (i + j) / 2
  py <- (j - i) / 2
  list(px = px, py = py)
}

make_triangle_hic <- function(chr, start, end, binsize, loops_df = NULL) {
  loc <- sprintf("%s:%d:%d", chr, start, end)

  last_err <- NULL
  mat <- NULL
  for (norm in norm_try) {
    mat <- tryCatch(
      strawr::straw(norm, hic_path, loc, loc, "BP", binsize),
      error = function(e) { last_err <<- e; NULL }
    )
    if (!is.null(mat) && nrow(mat) > 0) break
  }
  if (is.null(mat) || nrow(mat) == 0) {
    stop("Failed to fetch Hi-C data for ", loc, " (last error: ", if (!is.null(last_err)) last_err$message else "unknown", ")")
  }

  tri <- mat %>%
    mutate(
      i = (x - start) / binsize,
      j = (y - start) / binsize
    ) %>%
    filter(i >= 0, j >= 0, j >= i) %>%
    mutate(
      px = (i + j) / 2,
      py = (j - i) / 2,
      # Match the paper-like "0-4" scale feel using log10.
      z = log10(counts + 1)
    )

  # Paper-like contrast: normalize by a high quantile, then apply gamma to boost low values.
  # This makes the heatmap much darker without losing dynamic range.
  zcap <- as.numeric(stats::quantile(tri$z, probs = 0.99, na.rm = TRUE))
  if (!is.finite(zcap) || zcap <= 0) zcap <- max(tri$z, na.rm = TRUE)
  tri <- tri %>%
    mutate(
      z = pmin(z, zcap),
      z = (z / zcap),
      z = (z ^ hic_gamma) * 4,
      z = pmin(z, 4)
    )

  p <- ggplot(tri, aes(px, py, fill = z)) +
    # Use tiles (px/py are on a regular 0.5-grid; raster warns about uneven spacing).
    geom_tile(width = 1, height = 1) +
    scale_fill_gradientn(
      # Faster ramp to strong reds (matches typical Capture-C/Hi-C figure feel)
      colors = c("#ffffff", "#fee5d9", "#fcae91", "#fb6a4a", "#de2d26", "#a50f15"),
      limits = c(0, 4),
      breaks = c(0, 2, 4),
      oob = scales::squish,
      name = NULL
    ) +
    coord_fixed(expand = FALSE) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.key.height = unit(1.6, "in"),
      legend.key.width = unit(0.12, "in"),
      legend.text = element_text(size = 9, color = "black"),
      plot.margin = margin(2, 2, 0, 2)
    )

  # Overlay loop markers (arrow + label) on the heatmap for a paper-like callout style.
  if (!is.null(loops_df) && nrow(loops_df) > 0) {
    # Pick a modest number to avoid clutter.
    loops_pick <- loops_df %>%
      mutate(span_bp = abs(centroid2 - centroid1), span_mb = span_bp / 1e6) %>%
      arrange(desc(observed)) %>%
      slice_head(n = 80) %>%
      mutate(
        c1 = pmin(centroid1, centroid2),
        c2 = pmax(centroid1, centroid2)
      )

    coords <- compute_triangle_coords(loops_pick$c1, loops_pick$c2, start, binsize)
    loops_pick <- loops_pick %>%
      mutate(px = coords$px, py = coords$py) %>%
      filter(is.finite(px), is.finite(py), py >= 0)

    # Professor preference: show loop location with black dots only (no arrows/labels).
    p <- p +
      geom_point(
        data = loops_pick,
        aes(x = px, y = py),
        inherit.aes = FALSE,
        color = "black",
        size = 0.9,
        alpha = 0.9
      )
  }

  p
}

make_ctcf_track <- function(chr, start, end) {
  ctcf <- ctcf_all %>%
    filter(chr == !!chr, end >= start, start <= end) %>%
    mutate(mid = (start + end) / 2)

  ggplot(ctcf, aes(x = mid, y = score)) +
    geom_col(width = binsize, fill = "#3b0f70") +
    coord_cartesian(xlim = c(start, end), expand = FALSE) +
    theme_minimal(base_size = 9) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(0, 2, 0, 2)
    ) +
  ylab("CTCF")
}

make_gene_track <- function(chr, start, end, highlight_gene) {
  gr <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end))
  gtf <- rtracklayer::import(gtf_path, which = gr)

  tx <- gtf[gtf$type == "transcript"]
  ex <- gtf[gtf$type == "exon"]
  if (length(tx) == 0 || length(ex) == 0) {
    return(ggplot() + theme_void())
  }

  tx_df <- tibble(
    transcript_id = as.character(mcols(tx)$transcript_id),
    gene = if (!is.null(mcols(tx)$gene_name)) as.character(mcols(tx)$gene_name) else as.character(mcols(tx)$gene_id),
    start = start(tx),
    end = end(tx),
    strand = as.character(strand(tx))
  ) %>%
    filter(!is.na(transcript_id), !is.na(gene), end >= start, start <= end) %>%
    mutate(mid = (start + end) / 2, is_hi = gene == highlight_gene)

  tx_df <- pack_rows(tx_df, "start", "end")

  ex_df <- tibble(
    transcript_id = as.character(mcols(ex)$transcript_id),
    gene = if (!is.null(mcols(ex)$gene_name)) as.character(mcols(ex)$gene_name) else as.character(mcols(ex)$gene_id),
    start = start(ex),
    end = end(ex),
    strand = as.character(strand(ex))
  ) %>%
    filter(!is.na(transcript_id), !is.na(gene), end >= start, start <= end)

  ex_df <- ex_df %>%
    inner_join(tx_df %>% select(transcript_id, row, is_hi), by = "transcript_id") %>%
    mutate(
      ymin = row - 0.25,
      ymax = row + 0.25
    )

  gene_df <- tx_df %>%
    group_by(gene) %>%
    summarise(
      gene_start = min(start),
      gene_end = max(end),
      mid = as.integer((gene_start + gene_end) / 2),
      row = min(row),
      is_hi = any(is_hi),
      .groups = "drop"
    ) %>%
    arrange(desc(is_hi), mid)

  # Limit labels to reduce clutter: always keep highlight gene, then keep others spaced apart.
  min_sep_bp <- 250000L
  keep <- rep(FALSE, nrow(gene_df))
  last_mid <- NA_integer_
  for (i in seq_len(nrow(gene_df))) {
    if (gene_df$is_hi[i]) {
      keep[i] <- TRUE
      last_mid <- gene_df$mid[i]
      next
    }
    if (is.na(last_mid) || abs(gene_df$mid[i] - last_mid) >= min_sep_bp) {
      keep[i] <- TRUE
      last_mid <- gene_df$mid[i]
    }
    if (sum(keep) >= 12) break
  }
  label_df <- gene_df[keep, , drop = FALSE]

  ggplot() +
    # transcript lines
    geom_segment(
      data = tx_df,
      aes(x = start, xend = end, y = row, yend = row, color = is_hi),
      linewidth = 0.35
    ) +
    # exon boxes (BED12-like)
    geom_rect(
      data = ex_df,
      aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = is_hi),
      color = NA,
      alpha = 0.95
    ) +
    {
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        ggrepel::geom_text_repel(
          data = label_df,
          aes(x = mid, y = row + 0.65, label = gene, color = is_hi),
          size = 2.8,
          fontface = ifelse(label_df$is_hi, "bold", "plain"),
          min.segment.length = 0,
          segment.color = "grey70",
          box.padding = 0.25,
          point.padding = 0.1,
          direction = "x",
          max.overlaps = Inf,
          seed = 1,
          inherit.aes = FALSE
        )
      } else {
        geom_text(
          data = label_df,
          aes(x = mid, y = row + 0.65, label = gene, color = is_hi),
          size = 2.8,
          fontface = ifelse(label_df$is_hi, "bold", "plain"),
          inherit.aes = FALSE
        )
      }
    } +
    scale_color_manual(values = c(`TRUE` = "#d7191c", `FALSE` = "black"), guide = "none") +
    scale_fill_manual(values = c(`TRUE` = "#d7191c", `FALSE` = "black"), guide = "none") +
    coord_cartesian(xlim = c(start, end), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(0, 2, 2, 2))
}

make_distance_track <- function(start, end, origin, step_bp = 500000L) {
  # Simple "distance from origin" track (e.g., -2000kb ... 0 ... +2000kb)
  # that sits above the gene annotation.
  start <- as.integer(start)
  end <- as.integer(end)
  origin <- as.integer(origin)

  offsets <- seq(from = start - origin, to = end - origin, by = step_bp)
  ticks <- tibble(
    x = origin + offsets,
    label = ifelse(offsets == 0, "0", sprintf("%+dkb", as.integer(offsets / 1000)))
  ) %>%
    filter(x >= start, x <= end)

  ggplot() +
    geom_segment(aes(x = start, xend = end, y = 0, yend = 0), linewidth = 0.35, color = "black") +
    geom_segment(data = ticks, aes(x = x, xend = x, y = 0, yend = 0.18), linewidth = 0.35, color = "black") +
    geom_text(data = ticks, aes(x = x, y = 0.35, label = label), size = 2.6, color = "black") +
    coord_cartesian(xlim = c(start, end), ylim = c(-0.1, 0.55), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(0, 2, 0, 2))
}

clamp_region <- function(chr, start, end) {
  chr_end <- chrom_sizes$end[match(chr, chrom_sizes$chr)]
  if (is.na(chr_end)) stop("Unknown chr in chrom.sizes: ", chr)
  s <- max(0L, as.integer(start))
  e <- min(as.integer(chr_end), as.integer(end))
  if (e <= s) stop("Invalid clamped region for ", chr, ": ", s, "-", e)
  list(start = s, end = e)
}

for (i in seq_len(nrow(genes_of_interest))) {
  g <- genes_of_interest[i, ]
  mid <- as.integer((g$gene_start + g$gene_end) / 2)
  reg <- clamp_region(g$chr, mid - flank_bp, mid + flank_bp)

  message("Rendering ", g$gene, " ", g$chr, ":", reg$start, "-", reg$end)

  loops_region <- loops_all %>%
    filter(chr1 == !!g$chr, chr2 == !!g$chr) %>%
    filter(pmin(centroid1, centroid2) >= reg$start, pmax(centroid1, centroid2) <= reg$end)

  p_hic <- make_triangle_hic(g$chr, reg$start, reg$end, binsize, loops_df = loops_region)
  p_ctcf <- make_ctcf_track(g$chr, reg$start, reg$end)
  p_dist <- make_distance_track(reg$start, reg$end, origin = mid, step_bp = 500000L)
  p_genes <- make_gene_track(g$chr, reg$start, reg$end, g$gene)

  panel <- p_hic / p_ctcf / p_dist / p_genes +
    plot_layout(heights = c(7.0, 1.6, 0.7, 2.2))

  out_pdf <- file.path(out_dir, sprintf("DA68A_%s_%s_%d_%d.pdf", g$gene, g$chr, reg$start, reg$end))
  ggsave(out_pdf, panel, width = 9.5, height = 7.2, units = "in", dpi = 300)
  message("Wrote: ", out_pdf)

  out_png <- file.path(out_dir, sprintf("DA68A_%s_%s_%d_%d.png", g$gene, g$chr, reg$start, reg$end))
  ggsave(out_png, panel, width = 9.5, height = 7.2, units = "in", dpi = 300, bg = "white")
  message("Wrote: ", out_png)
}
