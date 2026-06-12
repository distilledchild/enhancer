####################################################
# ATAC-seq Validation of P-E Loop Anchors
# Duttke et al. 2022 snATAC-seq (rn6 → rn7 liftOver) x Hi-C loop anchors
# Source: GSM5820551 (rat PFC, filtered peak set)
# Q: Is enhancer anchor in open chromatin region?
####################################################
library(tidyverse)
library(GenomicRanges) # for genomic interval overlap and ranges operation

options(scipen = 999) # prevent scientific notation in base R
options(pillar.sigfig = 10) # display up to 10 significant digits in tibbles

setwd("/Users/pete/Desktop/playground/enhancer") # Mac
getwd() # /Users/pete/Desktop/playground/enhancer

path.csv.final.loop <- "/Users/pete/Desktop/playground/enhancer/figures/submission/lt2mb/df_final_loop_sub.4.any.lt2mb.ENSEMBL.mid.mid.final.filter.200kb.csv"
path.rds.final.loop <- "~/Dropbox/Gateway_to_Hao/enhancer/data/df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3_final_200kb.rds"
path.narrowpeak.atac <- "/Users/pete/Desktop/playground/enhancer/data/Duttke2022_snATAC/Duttke2022_snATAC_peaks_rn7.narrowPeak" # Duttke 2022 snATAC (rn6→rn7 liftOver)
path.dir.out <- "/Users/pete/Desktop/playground/enhancer/r_files"

####################################################
# 1. Loop load and anchor parsing
####################################################
df.loops.raw <- read_csv(path.csv.final.loop, show_col_types = FALSE) %>%
  dplyr::rename(loop_id = loop.id)
print(str_c("loops: ", nrow(df.loops.raw))) # 15085

# Split loop_id into coordinates and convert to integers (drop 7th+ fields)
df.loops <- df.loops.raw %>%
  mutate(loop_id_orig = loop_id) %>%
  separate_wider_delim(
    cols = loop_id_orig,
    delim = "_",
    names = c("chr1", "start1", "end1", "chr2", "start2", "end2"),
    too_many = "drop"
  ) %>%
  mutate(
    across(c(start1, end1, start2, end2), parse_integer),
    start1 = start1 + 1L,
    start2 = start2 + 1L
  )

df.loops %>% head(3)

####################################################
# 2. Load WHERE (promoter/enhancer anchor direction)
####################################################
rds.final.loop <- read_rds(path.expand(path.rds.final.loop))
df.where.loop <- rds.final.loop %>%
  dplyr::select(loop_id = loop.id, WHERE) %>%
  dplyr::distinct(loop_id, .keep_all = TRUE)
print(str_c("WHERE matched loops: ", nrow(df.where.loop))) # 17648

# Merge WHERE with loops
df.loops <- left_join(df.loops, df.where.loop, by = "loop_id")
print(str_c("WHERE NA: ", sum(is.na(df.loops$WHERE)))) # 0
t_dist <- table(df.loops$WHERE, useNA = "always")
print(str_c("WHERE distribution: ", str_c(str_c(coalesce(names(t_dist), "NA"), t_dist, sep = ": "), collapse = ", "))) # DOWN: 7510, UP: 7575, <NA>: 0

####################################################
# 3. Create Promoter/Enhancer anchor GRanges
####################################################
#    WHERE == "UP"   -> anchor1 is promoter, anchor2 is enhancer
#    WHERE == "DOWN" -> anchor2 is promoter, anchor1 is enhancer

# Use loops with WHERE information
df.loops.w <- df.loops %>% filter(!is.na(WHERE))
print(str_c("Loops with WHERE: ", nrow(df.loops.w))) # 15085

# promoter anchor GRanges
df.promoter <- bind_rows(
  df.loops.w %>%
    dplyr::filter(WHERE == "UP") %>%
    dplyr::select(chr = chr1, start = start1, end = end1, loop_id, category),
  df.loops.w %>%
    dplyr::filter(WHERE == "DOWN") %>%
    dplyr::select(chr = chr2, start = start2, end = end2, loop_id, category)
)
gr.promoter <- GRanges(
  seqnames = df.promoter$chr,
  ranges = IRanges(df.promoter$start, df.promoter$end),
  loop_id = df.promoter$loop_id,
  category = df.promoter$category,
  anchor_type = "promoter"
)

# enhancer anchor GRanges
df.enhancer <- bind_rows(
  df.loops.w %>%
    dplyr::filter(WHERE == "UP") %>%
    dplyr::select(chr = chr2, start = start2, end = end2, loop_id, category),
  df.loops.w %>%
    dplyr::filter(WHERE == "DOWN") %>%
    dplyr::select(chr = chr1, start = start1, end = end1, loop_id, category)
)
gr.enhancer <- GRanges(
  seqnames = df.enhancer$chr,
  ranges = IRanges(df.enhancer$start, df.enhancer$end),
  loop_id = df.enhancer$loop_id,
  category = df.enhancer$category,
  anchor_type = "enhancer"
)

print(str_c("Promoter anchors: ", length(gr.promoter))) # 15085
print(str_c("Enhancer anchors: ", length(gr.enhancer))) # 15085

# All anchors (reference)
gr.all <- c(gr.promoter, gr.enhancer)

####################################################
# 4. Load Duttke 2022 snATAC-seq peaks (rn7 liftOver)
####################################################
# func1: Load ATAC peaks and convert to GRanges
load_atac <- function(path, label) {
  df <- read_tsv(path,
    col_names = c(
      "chr", "start", "end", "name", "score", "strand",
      "fc", "neglog10p", "neglog10q", "summit"
    ),
    show_col_types = FALSE
  )
  print(str_c("peaks: ", nrow(df)))
  GRanges(df$chr, IRanges(df$start + 1L, df$end),
    score = df$score, fc = df$fc, sample = label
  )
}

gr.atac <- load_atac(path.narrowpeak.atac, "Duttke2022_snATAC_PFC")

# Merge overlapping peaks
gr.atac.union <- GenomicRanges::reduce(gr.atac)
print(str_c("Duttke2022 snATAC peaks (rn7): ", length(gr.atac)))
print(str_c("After reduce (merged): ", length(gr.atac.union)))

####################################################
# 5. Overlap analysis: promoter vs enhancer anchor x ATAC peaks
####################################################
# Check if each anchor overlaps with ATAC peak (TRUE/FALSE)
promoter.hits <- countOverlaps(gr.promoter, gr.atac.union) > 0
enhancer.hits <- countOverlaps(gr.enhancer, gr.atac.union) > 0

n.promoter <- length(gr.promoter) # 15085
n.enhancer <- length(gr.enhancer) # 15085
n.promoter.atac <- sum(promoter.hits) # 12160
n.enhancer.atac <- sum(enhancer.hits) # 11980
pct.promoter <- round(100 * n.promoter.atac / n.promoter, 1) # 12160/15085 = 80.6
pct.enhancer <- round(100 * n.enhancer.atac / n.enhancer, 1) # 11980/15085  = 79.4

### ATAC overlap results
print(str_c("Promoter anchors with ATAC peak: ", n.promoter.atac, " / ", n.promoter, " ", sprintf("(%.1f%%)", pct.promoter))) # 12160 / 15085 (80.6%)
print(str_c("Enhancer anchors with ATAC peak: ", n.enhancer.atac, " / ", n.enhancer, " ", sprintf("(%.1f%%)", pct.enhancer))) # 11980 / 15085 (79.4%)

# Fisher's exact test: is enhancer anchor more overlapped with ATAC peak?
mat <- matrix(
  c(
    n.enhancer.atac,  n.enhancer - n.enhancer.atac,
    n.promoter.atac,  n.promoter - n.promoter.atac
  ),
  nrow = 2, byrow = TRUE,
  dimnames = list(c("enhancer", "promoter"), c("ATAC_yes", "ATAC_no"))
)
### Contingency table
print(mat) # enhancer: 11980 / 3105, promoter: 12160 / 2925
#          ATAC_yes ATAC_no
# enhancer    11980    3105
# promoter    12160    2925

ft <- fisher.test(mat)
### Fisher's exact test
print(str_c("OR: ", round(ft$estimate, 3))) # 0.928
print(str_c("95% CI: ", round(ft$conf.int[1], 3), " - ", round(ft$conf.int[2], 3))) # 0.877 - 0.982
print(str_c("p-value: ", ft$p.value)) # 0.00996260331658994

####################################################
# 6. Permutation test: random genomic region vs ATAC overlap
####################################################
# set.seed for Random shift
set.seed(42)
obs.enh.pct <- n.enhancer.atac / n.enhancer

# Random shift
perm.pct <- map_dbl(1:1000, function(i) {
  g <- GenomicRanges::shift(
    gr.enhancer,
    sample(-5e6:5e6, length(gr.enhancer), replace = TRUE)
  )
  g <- g[start(g) >= 1]
  sum(countOverlaps(g, gr.atac.union) > 0) / length(g)
})

pm <- mean(perm.pct) # 0.464
ps <- sd(perm.pct) # 0.00385
zs <- (obs.enh.pct - pm) / ps # Z-score: 85.7
pp <- mean(perm.pct >= obs.enh.pct) # Empirical p: 0

print(str_c("Enhancer anchor ATAC overlap rate (obs): ", round(obs.enh.pct * 100, 1), " %")) # 79.4 %
print(str_c("Permutation mean: ", round(pm * 100, 1), " %")) # 46.4 %
print(str_c("Z-score: ", round(zs, 3))) # 85.702
print(str_c("Empirical p: ", pp)) # 0

####################################################
# 7. ATAC overlap by category (CP vs CT)
####################################################
### Promoter anchor category
for (cat_val in c("CP", "CT")) {
  idx <- gr.promoter$category == cat_val
  n <- sum(idx)
  h <- sum(promoter.hits[idx])
  print(str_c(cat_val, " : ", h, " / ", n, " (", sprintf("%.1f%%", 100 * h / n), ")")) # CP: 4094 / 4960 (82.5%), CT: 8066 / 10125 (79.7%)
}

### Enhancer anchor category
for (cat_val in c("CP", "CT")) {
  idx <- gr.enhancer$category == cat_val
  n <- sum(idx)
  h <- sum(enhancer.hits[idx])
  print(str_c(cat_val, " : ", h, " / ", n, " (", sprintf("%.1f%%", 100 * h / n), ")")) # CP: 3991 / 4960 (80.5%), CT: 7989 / 10125 (78.9%)
}

####################################################
# 8. Duttke2022 snATAC single-sample summary
####################################################
# (Single sample — no per-sample breakdown needed)
print("Duttke2022_snATAC_PFC (single sample = union)")
print(str_c("  Promoter: ", n.promoter.atac, " / ", n.promoter, " (", sprintf("%.1f%%", pct.promoter), ")"))
print(str_c("  Enhancer: ", n.enhancer.atac, " / ", n.enhancer, " (", sprintf("%.1f%%", pct.enhancer), ")"))

####################################################
# 9. Save results
####################################################
df.result.summary <- tibble(
  anchor_type      = c("promoter", "enhancer"),
  n_anchors        = c(n.promoter, n.enhancer),
  n_atac_overlap   = c(n.promoter.atac, n.enhancer.atac),
  pct_atac_overlap = c(pct.promoter, pct.enhancer),
  fisher_OR        = c(NA, round(ft$estimate, 3)),
  fisher_p         = c(NA, ft$p.value),
  perm_z           = c(NA, round(zs, 3)),
  perm_p           = c(NA, pp)
)
write_csv(df.result.summary, file.path(path.dir.out, "atac_loop_anchor_overlap_summary.csv"))

# Detailed ATAC overlap by anchor
df.detail <- tibble(
  loop_id = c(gr.promoter$loop_id, gr.enhancer$loop_id),
  anchor_type = c(
    rep("promoter", length(gr.promoter)),
    rep("enhancer", length(gr.enhancer))
  ),
  category = c(gr.promoter$category, gr.enhancer$category),
  atac_overlap = c(promoter.hits, enhancer.hits)
)
write_csv(df.detail, file.path(path.dir.out, "atac_loop_anchor_overlap_detail.csv"))

# Category (C=CTCF structural, P=Promoter functional, T=TSS functional)
#                           Promoter Anchor       Enhancer Anchor
# CP (CTCF + Promoter)     4,485/4,960 (90.4%)   4,451/4,960 (89.7%)
# CT (CTCF + TSS)          8,973/10,125 (88.6%)  8,977/10,125 (88.7%)
