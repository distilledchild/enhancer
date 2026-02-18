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
tad_bedpe_path <- path.expand("~/UTHSC GGI Dropbox/K P/Gateway_to_Hao/hic/2023A/hic_analysis/juicer/DA68A/intact/DA68A_intact_arrowhead/DA68A_intact_arrowhead_25000/25000_blocks.bedpe")
tad_bedpe_path_50kb <- path.expand("~/UTHSC GGI Dropbox/K P/Gateway_to_Hao/hic/2023A/hic_analysis/juicer/DA68A/intact/DA68A_intact_arrowhead/DA68A_intact_arrowhead_50000/50000_blocks.bedpe")
gtf_path <- path.expand("~/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/ucsc_refGene.gtf")
chrom_sizes_path <- path.expand("~/dropbox/Gateway_to_Hao/enhancer/data/tracks/rn7.chrom.sizes")
gene_loop_map_rds_path <- path.expand("~/dropbox/Gateway_to_Hao/enhancer/data/df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3_final_200kb.rds")

out_dir <- path.expand("~/dropbox/Gateway_to_Hao/enhancer/figures/hic_panels")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(hic_path))
stopifnot(file.exists(ctcf_bedgraph_path))
stopifnot(file.exists(loops_bedpe_path))
stopifnot(file.exists(tad_bedpe_path))
stopifnot(file.exists(tad_bedpe_path_50kb))
stopifnot(file.exists(gtf_path))
stopifnot(file.exists(gene_loop_map_rds_path))

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
flank_bp <- 1000000L  # zoom-in: keep middle 50% of previous window (4Mb -> 2Mb)
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
      centroid1 = as.integer(centroid1), centroid2 = as.integer(centroid2),
      end.distance = as.integer(x2 - x1),
      loop.id = str_c(chr1, "_", x1, "_", x2, "_", chr2, "_", y1, "_", y2, "_", end.distance)
    )
}

read_tads <- function(path) {
  # Arrowhead blocks.bedpe has duplicated "score" column name in header comments.
  cols <- c(
    "chr1", "x1", "x2", "chr2", "y1", "y2", "name", "score1", "strand1", "strand2",
    "color", "score2", "uVarScore", "lVarScore", "upSign", "loSign"
  )
  read_tsv(path, comment = "#", col_names = cols, show_col_types = FALSE) %>%
    mutate(
      chr1 = as.character(chr1), chr2 = as.character(chr2),
      x1 = as.integer(x1), x2 = as.integer(x2), y1 = as.integer(y1), y2 = as.integer(y2)
    )
}

loops_all <- read_loops(loops_bedpe_path)
tads_all <- read_tads(tad_bedpe_path)
tads_all_50kb <- read_tads(tad_bedpe_path_50kb)
ctcf_all <- read_tsv(ctcf_bedgraph_path, col_names = c("chr","start","end","score"), show_col_types = FALSE) %>%
  mutate(chr = as.character(chr), start = as.integer(start), end = as.integer(end), score = as.numeric(score))

# Gene-loop mapping from enhancer_promoter_interaction pipeline (strict filtered result).
gene_loop_map_raw <- readRDS(gene_loop_map_rds_path) %>%
  as_tibble() %>%
  filter(!is.na(gene_name), !is.na(loop.id)) %>%
  mutate(
    gene_key = tolower(gene_name),
    chr1 = as.character(chr1),
    chr2 = as.character(chr2),
    x1 = as.integer(x1),
    x2 = as.integer(x2),
    y1 = as.integer(y1),
    y2 = as.integer(y2),
    centroid1 = as.integer((x1 + x2) / 2),
    centroid2 = as.integer((y1 + y2) / 2),
    loop_left = pmin(centroid1, centroid2),
    loop_right = pmax(centroid1, centroid2)
  ) %>%
  distinct(gene_key, loop.id, .keep_all = TRUE)

gene_to_loop_ids <- gene_loop_map_raw %>%
  group_by(gene_key) %>%
  summarise(loop_ids = list(unique(loop.id)), n_related_total = n_distinct(loop.id), .groups = "drop")

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

# Mark loop as target-related using loop.id mapping from enhancer_promoter_interaction.R results.
annotate_target_related_loops <- function(df, related_loop_ids = character()) {
  if (is.null(df) || nrow(df) == 0) return(df)
  df %>%
    mutate(
      is_target_related = loop.id %in% related_loop_ids
    )
}

make_triangle_hic <- function(chr, start, end, binsize, loops_df = NULL, tads_df = NULL) {
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

  # Overlay TAD domains as thin triangular outlines (Arrowhead blocks).
  if (!is.null(tads_df) && nrow(tads_df) > 0) {
    tads_tri <- tads_df %>%
      mutate(
        s = pmin(x1, y1),
        e = pmax(x2, y2),
        sbin = (s - start) / binsize,
        ebin = (e - start) / binsize
      ) %>%
      filter(is.finite(sbin), is.finite(ebin), ebin > sbin) %>%
      transmute(
        x_left = sbin,
        y_left = 0,
        x_apex = (sbin + ebin) / 2,
        y_apex = (ebin - sbin) / 2,
        x_right = ebin,
        y_right = 0
      )

    if (nrow(tads_tri) > 0) {
      p <- p +
        geom_segment(
          data = tads_tri,
          aes(x = x_left, y = y_left, xend = x_apex, yend = y_apex),
          inherit.aes = FALSE,
          color = "#2ca25f",
          linewidth = 0.25,
          alpha = 0.45
        ) +
        geom_segment(
          data = tads_tri,
          aes(x = x_apex, y = y_apex, xend = x_right, yend = y_right),
          inherit.aes = FALSE,
          color = "#2ca25f",
          linewidth = 0.25,
          alpha = 0.45
        )
    }
  }

  # Overlay loop markers (arrow + label) on the heatmap for a paper-like callout style.
  if (!is.null(loops_df) && nrow(loops_df) > 0) {
    # Show all loop dots (per preference). This can be dense, so use small size + alpha.
    loops_pick <- loops_df %>%
      mutate(span_bp = abs(centroid2 - centroid1), span_mb = span_bp / 1e6) %>%
      arrange(desc(observed)) %>%
      mutate(
        c1 = pmin(centroid1, centroid2),
        c2 = pmax(centroid1, centroid2)
      )

    coords <- compute_triangle_coords(loops_pick$c1, loops_pick$c2, start, binsize)
    loops_pick <- loops_pick %>%
      mutate(px = coords$px, py = coords$py) %>%
      filter(is.finite(px), is.finite(py), py >= 0)

    # Professor preference: show loop location with black dots only (no arrows/labels).
    if (!("loop_color" %in% colnames(loops_pick))) {
      loops_pick <- loops_pick %>%
        mutate(loop_color = ifelse(is_target_related, "#2c7fb8", "black"))
    }

    p <- p +
      geom_point(
        data = loops_pick,
        aes(x = px, y = py),
        inherit.aes = FALSE,
        color = loops_pick$loop_color,
        size = 0.78, # +20%
        alpha = 0.55
      )
  }

  p
}

add_aux_loop_legend <- function(p, start, end, binsize, n_blue = NA_integer_, n_purple = NA_integer_) {
  x_max <- (end - start) / binsize
  y_max <- x_max / 2

  x0 <- x_max * 0.03
  y0 <- y_max * 0.92
  dy <- y_max * 0.05
  blue_label <- if (is.na(n_blue)) {
    "Related in main panel"
  } else {
    sprintf("Related in main panel (%d)", as.integer(n_blue))
  }
  purple_label <- if (is.na(n_purple)) {
    "Added related"
  } else {
    sprintf("Added related (%d)", as.integer(n_purple))
  }

  p +
    annotate("point", x = x0, y = y0, color = "#2c7fb8", size = 2.3, alpha = 0.9) +
    annotate("text", x = x0 + x_max * 0.02, y = y0, label = blue_label, hjust = 0, vjust = 0.5, size = 2.8, color = "black") +
    annotate("point", x = x0, y = y0 - dy, color = "#c7a0ff", size = 2.3, alpha = 0.95) +
    annotate("text", x = x0 + x_max * 0.02, y = y0 - dy, label = purple_label, hjust = 0, vjust = 0.5, size = 2.8, color = "black")
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

make_coord_track <- function(start, end, step_bp = 250000L) {
  # Genomic coordinate track (real coordinates) that sits above the gene annotation.
  start <- as.integer(start)
  end <- as.integer(end)
  step_bp <- as.integer(step_bp)

  first_tick <- as.integer(ceiling(start / step_bp) * step_bp)
  xs <- seq(from = first_tick, to = end, by = step_bp)
  ticks <- tibble(
    x = xs,
    label = sprintf("%.2fMb", xs / 1e6)
  )

  ggplot() +
    geom_segment(aes(x = start, xend = end, y = 0, yend = 0), linewidth = 0.35, color = "black") +
    geom_segment(data = ticks, aes(x = x, xend = x, y = 0, yend = 0.18), linewidth = 0.35, color = "black") +
    geom_text(data = ticks, aes(x = x, y = 0.35, label = label), size = 2.6, color = "black") +
    coord_cartesian(xlim = c(start, end), ylim = c(-0.1, 0.55), expand = FALSE) +
    theme_void() +
    # Keep vertical whitespace tight so it sits close to the gene track below.
    theme(plot.margin = margin(0, 2, -2, 2))
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

  gene_key_target <- tolower(g$gene)
  related_tbl <- gene_to_loop_ids %>% filter(gene_key == !!gene_key_target)
  related_loop_ids <- if (nrow(related_tbl) == 0) character() else related_tbl$loop_ids[[1]]
  related_loops_gene <- gene_loop_map_raw %>%
    filter(gene_key == !!gene_key_target, chr1 == !!g$chr, chr2 == !!g$chr)

  loops_region <- annotate_target_related_loops(
    loops_region,
    related_loop_ids = related_loop_ids
  )
  loops_region <- loops_region %>%
    mutate(loop_color = ifelse(is_target_related, "#2c7fb8", "black"))
  n_related_in_window <- loops_region %>%
    filter(is_target_related) %>%
    summarise(n = n_distinct(loop.id)) %>%
    pull(n)
  n_related_total <- if (nrow(related_tbl) == 0) 0L else as.integer(related_tbl$n_related_total[[1]])
  visible_related_ids_main <- loops_region %>%
    filter(is_target_related) %>%
    distinct(loop.id) %>%
    pull(loop.id)
  message(
    "[loop-map] ", g$gene,
    " related loops total=", n_related_total,
    ", in-window=", n_related_in_window
  )

  # TADs from Arrowhead blocks (25kb), intrachromosomal domains in region.
  tads_region <- tads_all %>%
    filter(chr1 == !!g$chr, chr2 == !!g$chr) %>%
    filter(x1 == y1, x2 == y2) %>%
    filter(x2 >= reg$start, x1 <= reg$end)

  p_hic <- make_triangle_hic(g$chr, reg$start, reg$end, binsize, loops_df = loops_region, tads_df = tads_region)
  p_ctcf <- make_ctcf_track(g$chr, reg$start, reg$end)
  p_coord <- make_coord_track(reg$start, reg$end, step_bp = 250000L)
  p_genes <- make_gene_track(g$chr, reg$start, reg$end, g$gene)

  panel <- p_hic / p_ctcf / p_coord / p_genes +
    plot_layout(heights = c(7.0, 1.6, 0.45, 2.2))

  out_pdf <- file.path(out_dir, sprintf("DA68A_%s_%s_%d_%d.pdf", g$gene, g$chr, reg$start, reg$end))
  ggsave(out_pdf, panel, width = 9.5, height = 7.2, units = "in", dpi = 300)
  message("Wrote: ", out_pdf)

  out_png <- file.path(out_dir, sprintf("DA68A_%s_%s_%d_%d.png", g$gene, g$chr, reg$start, reg$end))
  ggsave(out_png, panel, width = 9.5, height = 7.2, units = "in", dpi = 300, bg = "white")
  message("Wrote: ", out_png)

  # AUX figure: expand window to include all target-related loops.
  if (length(related_loop_ids) > 0) {
    if (nrow(related_loops_gene) > 0) {
      aux_start_raw <- min(c(reg$start, related_loops_gene$loop_left), na.rm = TRUE)
      aux_end_raw <- max(c(reg$end, related_loops_gene$loop_right), na.rm = TRUE)
      reg_aux <- clamp_region(g$chr, aux_start_raw, aux_end_raw)

      loops_aux <- loops_all %>%
        filter(chr1 == !!g$chr, chr2 == !!g$chr) %>%
        filter(pmin(centroid1, centroid2) >= reg_aux$start, pmax(centroid1, centroid2) <= reg_aux$end) %>%
        annotate_target_related_loops(related_loop_ids = related_loop_ids) %>%
        mutate(
          in_main_window = (pmin(centroid1, centroid2) >= reg$start) & (pmax(centroid1, centroid2) <= reg$end),
          added_for_aux = loop.id %in% setdiff(related_loop_ids, visible_related_ids_main),
          loop_color = case_when(
            is_target_related & added_for_aux ~ "#c7a0ff",   # light purple: loops not shown in main panel
            is_target_related & !added_for_aux ~ "#2c7fb8",  # blue: loops already shown in main panel
            TRUE ~ "black"
          )
        )
      added_related_ids <- setdiff(related_loop_ids, visible_related_ids_main)
      related_aux_map_only <- related_loops_gene %>%
        filter(loop.id %in% added_related_ids) %>%
        filter(loop_left >= reg_aux$start, loop_right <= reg_aux$end) %>%
        filter(!(loop.id %in% loops_aux$loop.id)) %>%
        mutate(
          is_target_related = TRUE,
          in_main_window = (loop_left >= reg$start) & (loop_right <= reg$end),
          added_for_aux = TRUE,
          loop_color = "#c7a0ff"
        )
      if (nrow(related_aux_map_only) > 0) {
        loops_aux <- bind_rows(loops_aux, related_aux_map_only)
      }

      n_related_aux <- length(unique(related_loop_ids))
      n_related_purple <- length(unique(added_related_ids))
      message(
        "[aux-loop-map] ", g$gene,
        " related in-aux=", n_related_aux,
        ", added(light-purple)=", n_related_purple
      )

      tads_aux <- tads_all %>%
        filter(chr1 == !!g$chr, chr2 == !!g$chr) %>%
        filter(x1 == y1, x2 == y2) %>%
        filter(x2 >= reg_aux$start, x1 <= reg_aux$end)

      p_hic_aux <- make_triangle_hic(g$chr, reg_aux$start, reg_aux$end, binsize, loops_df = loops_aux, tads_df = tads_aux)
      n_blue <- length(unique(visible_related_ids_main))
      n_purple <- length(unique(added_related_ids))
      p_hic_aux <- add_aux_loop_legend(
        p_hic_aux, reg_aux$start, reg_aux$end, binsize,
        n_blue = n_blue,
        n_purple = n_purple
      )
      p_ctcf_aux <- make_ctcf_track(g$chr, reg_aux$start, reg_aux$end)
      p_coord_aux <- make_coord_track(reg_aux$start, reg_aux$end, step_bp = 250000L)
      p_genes_aux <- make_gene_track(g$chr, reg_aux$start, reg_aux$end, g$gene)

      panel_aux <- p_hic_aux / p_ctcf_aux / p_coord_aux / p_genes_aux +
        plot_layout(heights = c(7.0, 1.6, 0.45, 2.2))

      out_pdf_aux <- file.path(out_dir, sprintf("aux_DA68A_%s_%s_%d_%d.pdf", g$gene, g$chr, reg_aux$start, reg_aux$end))
      ggsave(out_pdf_aux, panel_aux, width = 9.5, height = 7.2, units = "in", dpi = 300)
      message("Wrote: ", out_pdf_aux)

      out_png_aux <- file.path(out_dir, sprintf("aux_DA68A_%s_%s_%d_%d.png", g$gene, g$chr, reg_aux$start, reg_aux$end))
      ggsave(out_png_aux, panel_aux, width = 9.5, height = 7.2, units = "in", dpi = 300, bg = "white")
      message("Wrote: ", out_png_aux)
    }
  }

  # Add one comparison panel with 50kb TAD annotation (requested: one extra panel).
  if (g$gene == "Tgfb2") {
    tads_region_50kb <- tads_all_50kb %>%
      filter(chr1 == !!g$chr, chr2 == !!g$chr) %>%
      filter(x1 == y1, x2 == y2) %>%
      filter(x2 >= reg$start, x1 <= reg$end)

    p_hic_50kb <- make_triangle_hic(g$chr, reg$start, reg$end, binsize, loops_df = loops_region, tads_df = tads_region_50kb)
    panel_50kb <- p_hic_50kb / p_ctcf / p_coord / p_genes +
      plot_layout(heights = c(7.0, 1.6, 0.45, 2.2))

    out_pdf_50kb <- file.path(out_dir, sprintf("DA68A_%s_%s_%d_%d_TAD50kb.pdf", g$gene, g$chr, reg$start, reg$end))
    ggsave(out_pdf_50kb, panel_50kb, width = 9.5, height = 7.2, units = "in", dpi = 300)
    message("Wrote: ", out_pdf_50kb)

    out_png_50kb <- file.path(out_dir, sprintf("DA68A_%s_%s_%d_%d_TAD50kb.png", g$gene, g$chr, reg$start, reg$end))
    ggsave(out_png_50kb, panel_50kb, width = 9.5, height = 7.2, units = "in", dpi = 300, bg = "white")
    message("Wrote: ", out_png_50kb)
  }
}
