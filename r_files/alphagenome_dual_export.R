################################################################################################
################################################################################################
# AlphaGenome dual-anchor input generation (enhancer + promoter)
# Flow summary:
# - choose top genes from filtered loops
# - keep single-sided loops
# - define enhancer anchor (opposite side) and promoter anchor (annotated side)
# - build 1,048,576bp windows centered on each anchor midpoint (role-specific windows)
# - keep enhancer/promoter windows separate (no shared-window merge)
# - fetch rn7 DNA via UCSC API for unique windows
# - export A_dual.csv for dual-role AlphaGenome analysis
################################################################################################
################################################################################################

ag.script.path <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
ag.script.dir <- ifelse(length(ag.script.path) > 0, dirname(normalizePath(ag.script.path)), getwd())
ag.project.root <- normalizePath(file.path(ag.script.dir, ".."))

if (!exists("parse_loop_id") || !exists("str_split_n")) {
  source(file.path(ag.script.dir, "utils_functions.R"))
}

if (!exists("df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3")) {
  ag.cache.candidates <- c(
    file.path(ag.project.root, "data", "df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3_final_200kb.rds"),
    path.expand("~/dropbox/Gateway_to_Hao/enhancer/data/df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3_final_200kb.rds"),
    path.expand("~/Dropbox/Gateway_to_Hao/enhancer/data/df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3_final_200kb.rds")
  )
  ag.cache.rds <- ag.cache.candidates[file.exists(ag.cache.candidates)][1]
  if (!is.na(ag.cache.rds)) {
    df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3 <- readRDS(ag.cache.rds)
    message("[AlphaGenome dual] Loaded cached object: ", ag.cache.rds)
  } else {
    stop(
      "Missing object: df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3 and cache not found in:\n  - ",
      paste(ag.cache.candidates, collapse = "\n  - ")
    )
  }
}

ag.threshold.distance <- 200000
ag.top.n.genes <- 54
ag.output.dir <- file.path(ag.project.root, "data", "alphagenome", "top54_genes_loop_anchors_200kb")
ag.genome <- "rn7"
ag.api.sleep.sec <- 0.02
ag.log.every <- 25
ag.window.bp <- 1048576L
ag.window.half.size <- as.integer(ag.window.bp / 2L)
ag.chrom.size.path <- file.path(ag.project.root, "data", "rn7_chromosome_length_from_ucsc.tsv")

dir.create(ag.output.dir, recursive = TRUE, showWarnings = FALSE)

########################
# STEP 1. Top genes + single-sided loops
########################
df.ag.input <- df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3 %>%
  dplyr::rename(gene_id_id = gene_id) %>%
  mutate(gene_symbol = str_split_n(str_split_n(component_id, ":", 6), "\\|", 1))

df.ag.top.genes <- df.ag.input %>%
  filter(distance <= ag.threshold.distance) %>%
  count(gene_symbol, sort = TRUE) %>%
  slice_max(n, n = ag.top.n.genes, with_ties = FALSE) %>%
  mutate(top_rank = row_number()) %>%
  dplyr::rename(n_loops = n)

df.ag.gene.loop <- df.ag.input %>%
  filter(distance <= ag.threshold.distance, gene_symbol %in% df.ag.top.genes$gene_symbol) %>%
  distinct(gene_symbol, loop.id, .keep_all = TRUE) %>%
  dplyr::select(gene_symbol, loop.id, distance, component, WHERE)

df.ag.single.side <- df.ag.gene.loop %>%
  group_by(loop.id) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  filter(WHERE %in% c("UP", "DOWN"))

df.ag.loop.coords <- parse_loop_id(unique(df.ag.single.side$loop.id)) %>% # utils_functions.R
  dplyr::select(loop.id, chr1, x1, x2, chr2, y1, y2, end.distance)

df.ag.single.side.wide <- df.ag.single.side %>%
  left_join(df.ag.loop.coords, by = "loop.id") %>%
  left_join(df.ag.top.genes, by = "gene_symbol") %>%
  relocate(top_rank, gene_symbol, n_loops, loop.id, chr1, x1, x2, chr2, y1, y2, end.distance, .before = distance)

########################
# STEP 2. Build enhancer/promoter anchor intervals per loop
########################
# WHERE == UP   -> annotated side is anchor1 (promoter), enhancer is anchor2
# WHERE == DOWN -> annotated side is anchor2 (promoter), enhancer is anchor1
df.ag.pairs <- df.ag.single.side.wide %>%
  mutate(
    pair_uid = str_c("gene=", gene_symbol, "|loop=", loop.id),
    promoter_chr = ifelse(WHERE == "UP", chr1, chr2),
    promoter_start0 = ifelse(WHERE == "UP", as.integer(x1), as.integer(y1)),
    promoter_end = ifelse(WHERE == "UP", as.integer(x2), as.integer(y2)),
    enhancer_chr = ifelse(WHERE == "UP", chr2, chr1),
    enhancer_start0 = ifelse(WHERE == "UP", as.integer(y1), as.integer(x1)),
    enhancer_end = ifelse(WHERE == "UP", as.integer(y2), as.integer(x2)),
    promoter_mid = as.integer(floor((promoter_start0 + promoter_end) / 2)),
    enhancer_mid = as.integer(floor((enhancer_start0 + enhancer_end) / 2))
  )

if (file.exists(ag.chrom.size.path)) {
  df.chrom.sizes <- read.table(ag.chrom.size.path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(df.chrom.sizes) <- c("chr", "chr_len")
} else {
  message("[AlphaGenome dual] Chrom size file not found: ", ag.chrom.size.path, " (window capped by 0 only)")
  df.chrom.sizes <- tibble(chr = character(0), chr_len = integer(0))
}

get_chr_len <- function(chr) {
  v <- df.chrom.sizes$chr_len[df.chrom.sizes$chr == chr]
  if (length(v) == 0) {
    return(NA_integer_)
  }
  as.integer(v[[1]])
}

compute_window_1mb <- function(chr, mid, window_bp = 1048576L) {
  half <- as.integer(window_bp / 2L)
  chr_len <- get_chr_len(chr)
  ws <- as.integer(max(0L, as.integer(mid) - half))
  we <- as.integer(ws + window_bp)
  if (!is.na(chr_len)) {
    if (we > chr_len) {
      we <- as.integer(chr_len)
      ws <- as.integer(max(0L, we - window_bp))
    }
  }
  list(start0 = ws, end = we)
}

anchors_fit_in_window <- function(window_start0, window_end, s1, e1, s2, e2) {
  s1 >= window_start0 && e1 <= window_end && s2 >= window_start0 && e2 <= window_end
}

compute_pair_windows <- function(promoter_chr, promoter_start0, promoter_end, promoter_mid,
                                 enhancer_chr, enhancer_start0, enhancer_end, enhancer_mid) {
  # Use anchor-centered windows only: each role gets its own midpoint-based 1,048,576bp window.
  # Do not merge to shared windows, to keep per-anchor midpoint semantics exact.
  shared_window_used <- FALSE
  promoter_window_chr <- promoter_chr
  enhancer_window_chr <- enhancer_chr

  prom_win <- compute_window_1mb(promoter_chr, promoter_mid, ag.window.bp)
  enh_win <- compute_window_1mb(enhancer_chr, enhancer_mid, ag.window.bp)
  list(
    shared_window_used = shared_window_used,
    promoter_window_chr = promoter_window_chr,
    promoter_window_start0 = prom_win$start0,
    promoter_window_end = prom_win$end,
    enhancer_window_chr = enhancer_window_chr,
    enhancer_window_start0 = enh_win$start0,
    enhancer_window_end = enh_win$end
  )
}

df.ag.pairs <- df.ag.pairs %>%
  rowwise() %>%
  mutate(
    window_info = list(compute_pair_windows(
      promoter_chr = promoter_chr, promoter_start0 = promoter_start0, promoter_end = promoter_end, promoter_mid = promoter_mid,
      enhancer_chr = enhancer_chr, enhancer_start0 = enhancer_start0, enhancer_end = enhancer_end, enhancer_mid = enhancer_mid
    ))
  ) %>%
  ungroup() %>%
  tidyr::unnest_wider(window_info)

########################
# STEP 3. Long role table (one row per pair-role)
########################
df.ag.roles <- bind_rows(
  df.ag.pairs %>%
    transmute(
      pair_uid = pair_uid,
      gene_symbol = gene_symbol,
      loop.id = loop.id,
      top_rank = top_rank,
      n_loops = n_loops,
      WHERE = WHERE,
      role = "enhancer",
      role_box_color = "yellow",
      anchor_chr = enhancer_chr,
      anchor_start0 = as.integer(enhancer_start0),
      anchor_end = as.integer(enhancer_end),
      anchor_mid = as.integer(enhancer_mid),
      anchor_width_bp = as.integer(enhancer_end - enhancer_start0),
      window_chr = enhancer_window_chr,
      window_start0 = as.integer(enhancer_window_start0),
      window_end = as.integer(enhancer_window_end),
      shared_window_used = shared_window_used
    ),
  df.ag.pairs %>%
    transmute(
      pair_uid = pair_uid,
      gene_symbol = gene_symbol,
      loop.id = loop.id,
      top_rank = top_rank,
      n_loops = n_loops,
      WHERE = WHERE,
      role = "promoter",
      role_box_color = "pink",
      anchor_chr = promoter_chr,
      anchor_start0 = as.integer(promoter_start0),
      anchor_end = as.integer(promoter_end),
      anchor_mid = as.integer(promoter_mid),
      anchor_width_bp = as.integer(promoter_end - promoter_start0),
      window_chr = promoter_window_chr,
      window_start0 = as.integer(promoter_window_start0),
      window_end = as.integer(promoter_window_end),
      shared_window_used = shared_window_used
    )
) %>%
  mutate(
    anchor_uid = str_c(anchor_chr, ":", anchor_start0, "-", anchor_end),
    window_uid = str_c(window_chr, ":", window_start0, "-", window_end),
    role_uid = str_c(pair_uid, "|role=", role)
  ) %>%
  arrange(top_rank, gene_symbol, loop.id, role)

df.ag.roles %>% count(role)
df.ag.roles %>% count(shared_window_used)

########################
# STEP 4. Fetch sequences for unique windows (UCSC API)
########################
if (!requireNamespace("httr2", quietly = TRUE)) {
  stop("Package 'httr2' is required for UCSC API sequence fetch. Install with install.packages('httr2')")
}

fetch_ucsc_sequence <- function(chr, start0, end1, genome = "rn7") {
  tryCatch(
    {
      req <- httr2::request("https://api.genome.ucsc.edu/getData/sequence") %>%
        httr2::req_url_query(
          genome = genome,
          chrom = chr,
          start = as.integer(start0),
          end = as.integer(end1)
        )
      res <- httr2::req_perform(req)
      body <- httr2::resp_body_json(res, simplifyVector = TRUE)
      dna <- body$dna
      if (is.null(dna) || length(dna) == 0 || is.na(dna)) {
        return(NA_character_)
      }
      toupper(as.character(dna))
    },
    error = function(e) {
      message("[AlphaGenome dual][UCSC API] Failed: ", chr, ":", start0, "-", end1, " | ", conditionMessage(e))
      NA_character_
    }
  )
}

df.ag.windows.unique <- df.ag.roles %>%
  distinct(window_chr, window_start0, window_end, window_uid) %>%
  mutate(sequence = NA_character_)

n.total.windows <- nrow(df.ag.windows.unique)
n.failed.windows <- 0L
t0.fetch <- Sys.time()

if (n.total.windows > 0) {
  message("[AlphaGenome dual][UCSC API] Start fetching ", n.total.windows, " windows...")
  for (i in seq_len(n.total.windows)) {
    chr.i <- df.ag.windows.unique$window_chr[[i]]
    start.i <- df.ag.windows.unique$window_start0[[i]]
    end.i <- df.ag.windows.unique$window_end[[i]]
    Sys.sleep(ag.api.sleep.sec)
    seq.i <- fetch_ucsc_sequence(chr = chr.i, start0 = start.i, end1 = end.i, genome = ag.genome)
    df.ag.windows.unique$sequence[[i]] <- seq.i
    if (is.na(seq.i)) {
      n.failed.windows <- n.failed.windows + 1L
    }
    if (i == 1L || i %% ag.log.every == 0L || i == n.total.windows) {
      elapsed.sec <- as.numeric(difftime(Sys.time(), t0.fetch, units = "secs"))
      message(sprintf(
        "[AlphaGenome dual][UCSC API] progress %d/%d (%.1f%%) | failed=%d | elapsed=%.1fs",
        i, n.total.windows, (100 * i / n.total.windows), n.failed.windows, elapsed.sec
      ))
    }
  }
  message("[AlphaGenome dual][UCSC API] Fetch completed. failed=", n.failed.windows, "/", n.total.windows)
}

########################
# STEP 5. Export A_dual.csv
########################
df.ag.dual <- df.ag.roles %>%
  left_join(
    df.ag.windows.unique %>% dplyr::select(window_uid, sequence),
    by = "window_uid"
  ) %>%
  filter(!is.na(sequence)) %>%
  transmute(
    pair_uid = pair_uid,
    gene_symbol = gene_symbol,
    loop_id = loop.id,
    role = role,
    role_box_color = role_box_color,
    shared_window_used = shared_window_used,
    anchor_uid = anchor_uid,
    anchor_chr = anchor_chr,
    anchor_start0 = as.integer(anchor_start0),
    anchor_end = as.integer(anchor_end),
    anchor_width_bp = as.integer(anchor_width_bp),
    window_uid = window_uid,
    chr = window_chr,
    start0 = as.integer(window_start0),
    end = as.integer(window_end),
    sequence = sequence,
    sequence_length = nchar(sequence)
  )

write.csv(
  df.ag.dual,
  file = file.path(ag.output.dir, "A_dual.csv"),
  row.names = FALSE,
  quote = TRUE
)

message("[AlphaGenome dual] Completed. Output: ", file.path(ag.output.dir, "A_dual.csv"))
message(
  "[AlphaGenome dual] n_pairs=", n_distinct(df.ag.dual$pair_uid),
  ", n_role_rows=", nrow(df.ag.dual),
  ", n_unique_windows=", n_distinct(df.ag.dual$window_uid),
  ", shared_window_rows=", sum(df.ag.dual$shared_window_used, na.rm = TRUE)
)
