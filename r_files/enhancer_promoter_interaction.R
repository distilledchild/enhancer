library("ComplexUpset")
library("tidyverse")
library("GenomicRanges")
library("ggplot2")
library("devtools")
library("remotes")
library("gridExtra")
library("patchwork")
library("cowplot")
library("biomaRt")
library("reshape2")
library("ggvenn")
library("gtools")
library("scales") # for better axis formatting
library("ggrepel")
library("circlize")
library("RColorBrewer")
library("broom")
library("igraph")
library("ggraph")
library("writexl")
library("ggrepel")
library("pdftools")
library("magick")
library("tools")
library("httpgd")
library("RIdeogram")
library("rsvg")
library("rtracklayer")

options(tibble.width = Inf)
options(tibble.print_max = Inf)
options(tibble.max_extra_cols = Inf)
options(scipen = 999)

getwd()

# Linux
# setwd('C:\\Users\\panju\\Dropbox (UTHSC GGI)\\Gateway_to_Hao\\workshop\\2023_NIH_meeting\\loop_N_tss')
# setwd('./Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss')
# setwd("/home/pkim/dropbox/Gateway_to_Hao/enhancer/r_files")
# getwd()

# Mac
setwd("~/dropbox/Gateway_to_Hao/enhancer/r_files")
getwd()
source("~/Desktop/playground/enhancer/r_files/utils_functions.R") # Load all utility functions

########################
# 1. Loop
# 1-3. Loop data preprocessing: loops less than 2mb: df.DISTINCT.loop.deep.sample.all.lt.2mb (31,019)
# 1-4. Loop data preprocessing: padding on loops (1 distance)
# -> OVERALL.df.DISTINCT.loop.deep.sample.all : OBJECT to be used for the distribution of sth over loops
# x0, y3, new_distance, new.loop.id
########################
# loading all DISTINCT loops
cache_file_distinct_loop <- "../data/df.DISTINCT.loop.deep.sample.all.rds"
df.DISTINCT.loop.deep.sample.all <- readRDS(cache_file_distinct_loop) # 31,773
df.DISTINCT.loop.deep.sample.all.lt.2mb <- df.DISTINCT.loop.deep.sample.all %>% # 31,019: less than 2mb
  filter(distance < 2000000) # 31,019/31,773, only use less than 2mb
df.DISTINCT.loop.deep.sample.all.lt.2mb %>% dim() # 31,019 (31,773 - 754 (longer than 2mb))
df.DISTINCT.loop.deep.sample.all.lt.2mb %>% head(2)

# padding 1 distance with all loops: 31,773
OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance <- df.DISTINCT.loop.deep.sample.all %>% # 31,773
  mutate(
    x0 = ifelse(mid_x - distance < 0, 0, mid_x - distance),
    y3 = ifelse(mid_y + distance > chr.end.coord, chr.end.coord, mid_y + distance)
  ) # 1 distance for padding

# padding 1 distance && w/o capping only from DISTINCT loops: 31417
OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping <- OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance %>%
  filter(!(x0 == 0 | chr.end.coord == y3)) # 31773 - 444 (282 + 162) = 31329 + 88 = 31417
# filter(x0 == 0) # 282
# filter(chr.end.coord == y3)  # 162
OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping %>% dim() # 31417
OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance %>% filter(x0 == 0 & chr.end.coord == y3) # 88 |

# padding 1 distance && w/o capping && less than 2mb: 31329
OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping.lt.2mb <- OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance %>%
  filter(distance < 2000000) %>% # 31019                          ################ distance is the original loop distance (prior to padding)
  filter(!(x0 == 0 | chr.end.coord == y3)) # 31773 - 444 = 31329, 88// # 30928
# count(loop.id)  %>% view()
# filter(x0 == 0) # 282
# filter(chr.end.coord == y3)  # 162
OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping.lt.2mb %>% dim() # 30928 Q1: 150,000, Median: 195,000, Q3: 375,000

#########################
# 2. CTCF
#########################

cache_file_overlapping_CTCF_BOTH <- "../data/df.overlapping.CTCF.w.BOTH.result.rds"
df.overlapping.CTCF.w.BOTH.result <- readRDS(cache_file_overlapping_CTCF_BOTH)

# for Q1: quartile among loops with CTCF, so min = 1
df.ctcf.counts <- df.overlapping.CTCF.w.BOTH.result %>%
  group_by(loop.id, WHERE, resolution) %>%
  summarise(ctcf_count = n_distinct(ctcf.id), .groups = "drop")

df.ctcf.counts %>% count(ctcf_count)

df.ctcf.case <- df.ctcf.counts %>%
  group_by(loop.id) %>%
  summarise(
    direction_set = list(unique(WHERE)),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(case = case_when(
    length(direction_set) == 2 ~ "BOTH",
    length(direction_set) == 1 ~ direction_set[1],
    TRUE ~ NA_character_
  )) %>%
  ungroup() %>%
  dplyr::select(loop.id, case)

df.ctcf.case
df.ctcf.case %>% count(case)

df.loop.with.ctcf.case <- df.DISTINCT.loop.deep.sample.all.lt.2mb %>% # 31773
  left_join(df.ctcf.case, by = "loop.id") %>%
  mutate(case = ifelse(is.na(case), "NONE", case))

df.loop.with.ctcf.case %>% count(case)

df.DISTINCT.loop.deep.sample.all.lt.2mb %>% dim() # [1] 31019 16
# lt2mb
# case     n
# 1 BOTH 29978
# 2 DOWN   467
# 3 NONE   111
# 4   UP   463
# 29978 + 467 + 111 + 463 = 31019 PASS

ctcf.stats.by.resolution <- calculate_grouped_stats( # utils_functions.R
  df.ctcf.counts,
  value_col = "ctcf_count",
  group_cols = c("WHERE", "resolution")
) # Using utility function
ctcf.stats.by.resolution

### lt2mb, any
# # A tibble: 6 × 9
# WHERE resolution   Min    Q1 Median    Q3  Mean    SD   Max
# <fct> <fct>      <int> <dbl>  <dbl> <dbl> <dbl> <dbl> <int>
# 1 UP    5K             1     9     29    40  29.1  22.2   155
# 2 UP    10K            1    15     34    52  37.5  28.7   491
# 3 UP    25K            1    31     52    84  60.8  43.0   635
# 4 DOWN  5K             1    10     29    41  29.7  22.8   357
# 5 DOWN  10K            1    16     34    53  38.2  28.6   282
# 6 DOWN  25K            1    31     51    83  60.9  42.6   389


########################
# 3. TSS
########################
####################################################
# retrieving TSS data from Ensembl GTF             : df.ensembl.gtf.for.tss.DISTINCT.geneid
####################################################
df.ensembl.gtf.for.tss.DISTINCT.geneid <- readRDS("../data/df_ensembl_gtf_for_tss_DISTINCT_geneid.rds")
df.ensembl.gtf.for.tss.DISTINCT.geneid %>% head(2) # chr, start, end, strand, gene_id, gene_name, tss.id: chr1:157231467:157231469:+:ENSRNOG00000070168:Or51f23c
df.ensembl.gtf.for.tss.DISTINCT.geneid %>% dim() # 21,725
####################################################
# retrieving EXON data 1 from RefSeq GTF           : df.refgene.gtf.for.exon
####################################################
df.refgene.gtf.for.exon <- readRDS("../data/df_refgene_gtf_for_exon.rds")
df.refgene.gtf.for.exon %>% head(2) # chr, start, end, strand, exon_id, exon_number, gene_id, gene_name, transcript_id refseq_exon_id: chr7:92494376:92494506:-:A1bg:A1bg:8
df.refgene.gtf.for.exon %>% dim() # 17,488
####################################################
# retrieving EXON data 2 from Ensembl GTF           : df.tss.ensembl
####################################################
df.ensembl.gtf.for.exon.attribute <- readRDS("../data/df_ensembl_gtf_for_exon_attribute.rds")
df.tss.ensembl <- readRDS("../data/df.tss.ensembl.rds")
df.tss.ensembl %>% head(2)
# chr, start, end, strand, gene_id, gene_name
# tss.id: chr1:157231467:157231469:+:ENSRNOG00000070168:Or51f23c
# ensembl_exon_id: chr1:157231467:157232417:+:ENSRNOG00000070168:Or51f23c:1
# refseq_exon_id: chr1:80126460:80131881:+:Irgq:Irgq:3
df.tss.ensembl %>% dim() # 21,725

########################
# 4. promoter from EPD
########################
df.promoter.rn7.epd <- readRDS("../data/df.promoter.rn7.epd.rds")
df.promoter.rn7.epd %>% head(2) # seqnames, start, end, width, strand, promoter_id, score, gene_id, gene_name, ensembl_exon_id, refseq_exon_id, promoter.id
df.promoter.rn7.epd.GR <- GRanges(
  seqnames = df.promoter.rn7.epd$seqnames,
  ranges = IRanges(
    start = df.promoter.rn7.epd$start,
    end = df.promoter.rn7.epd$end
  )
)
# ensembl_exon_id
# 1 chr1:1807644:1807710:-:ENSRNOG00000040300:LOC120093164:6
# refseq_exon_id                         start     end
# 2 chr1:2106343:2108110:+:Lrp11:Lrp11:8 2079583 2079663
# promoter.id
# 1 chr1:1402581:1402661:-:ENSRNOG00000040300:Raet1e:chr1:1402621:1402622

# metadata
mcols(df.promoter.rn7.epd.GR) <- df.promoter.rn7.epd[, c("promoter.id", "gene_id", "gene_name", "ensembl_exon_id", "refseq_exon_id")]

##########################################################
##########################################################
# checking hits of TSS & promoter in an anchor using distanceToNearest
##########################################################
##########################################################

##########################################################
# 5. Combining TSS and Promoter datasets with exon information
##########################################################

# Prepare TSS data & Promoter data (already has exon information)
df.tss.ensembl %>% head(2) #      chr, start, end,        strand, gene_id, gene_name, tss.id,             ensembl_exon_id, refseq_exon_id
df.promoter.rn7.epd %>% head(2) # seqnames, start, end, width, strand, gene_id, gene_name, promoter_id, score, ensembl_exon_id, refseq_exon_id

# Combine TSS and Promoter datasets
df.gene_tss_and_pro <- bind_rows(
  df.tss.ensembl %>%
    dplyr::select(chr, start, end, gene_id, gene_name, tss.id, ensembl_exon_id, refseq_exon_id) %>%
    mutate(across(c(start, end), as.numeric)) %>%
    mutate(component = "tss") %>%
    dplyr::rename(component_id = tss.id) %>%
    mutate(component_id = str_c(component_id, component, sep = "|")),
  df.promoter.rn7.epd %>%
    dplyr::rename(chr = seqnames) %>%
    # mutate(promoter.id = str_c(chr, ':', start, ':', end, ':', gene_id, ':', gene_name)) %>%
    dplyr::select(chr, start, end, gene_id, gene_name, promoter.id, ensembl_exon_id, refseq_exon_id) %>%
    mutate(across(c(start, end), as.numeric)) %>%
    mutate(component = "pro") %>%
    dplyr::rename(component_id = promoter.id) %>%
    mutate(component_id = str_c(component_id, component, sep = "|"))
)

df.gene_tss_and_pro %>% head(2)
df.gene_tss_and_pro %>% dim() # TSS: 21,725 + Promoter: 12,529 = 34,254
df.gene_tss_and_pro_prep <- df.gene_tss_and_pro %>%
  mutate(gene_id = coalesce(ensembl_exon_id, refseq_exon_id)) %>% # priority: ensembl_exon_id > refseq_exon_id
  mutate(gene_chr = str_split_n(gene_id, ":", 1), gene_start = str_split_n(gene_id, ":", 2), gene_end = str_split_n(gene_id, ":", 3))
df.gene_tss_and_pro_prep %>% head(2)
df.gene_tss_and_pro_prep %>% dim()

# Create GRanges object with exon information
df.gene_tss_and_pro.GR <- GRanges(
  seqnames = df.gene_tss_and_pro_prep$chr,
  ranges = IRanges(
    start = df.gene_tss_and_pro_prep$start,
    end = df.gene_tss_and_pro_prep$end
  ),
  gene_id = df.gene_tss_and_pro_prep$gene_id,
  gene_name = df.gene_tss_and_pro_prep$gene_name,
  component_id = df.gene_tss_and_pro_prep$component_id,
  ensembl_exon_id = df.gene_tss_and_pro_prep$ensembl_exon_id,
  refseq_exon_id = df.gene_tss_and_pro_prep$refseq_exon_id,
  component = df.gene_tss_and_pro_prep$component,
  gene_chr = df.gene_tss_and_pro_prep$gene_chr,
  gene_start = df.gene_tss_and_pro_prep$gene_start,
  gene_end = df.gene_tss_and_pro_prep$gene_end
)

df.gene_tss_and_pro.GR # 34,254

df.DISTINCT.loop.deep.sample.all.lt.2mb %>% dim() # 31,019
df.DISTINCT.loop.deep.sample.all.lt.2mb %>% head(2)
colnames(df.DISTINCT.loop.deep.sample.all.lt.2mb) # "loop.id" "chr1" "x1" "x2" "chr2" "y1" "y2" "end.distance" "distance" "resolution" "x0" "x3" "y0" "y3" "mid_x" "mid_y" "mid_loop" "chr.end.coord" "start" "case"

# Create GRanges for UP and DOWN anchors using creating_granges function (PARAM: point = TRUE)
df.DISTINCT.loop.deep.sample.all.lt.2mb.UP.point.GR <- creating_granges(df.DISTINCT.loop.deep.sample.all.lt.2mb,
  direction = "up", use_anchor = TRUE, point = TRUE
) # utils_functions.R
df.DISTINCT.loop.deep.sample.all.lt.2mb.DOWN.point.GR <- creating_granges(df.DISTINCT.loop.deep.sample.all.lt.2mb,
  direction = "down", use_anchor = TRUE, point = TRUE
) # utils_functions.R

#############
# directional filtering for distanceToNearest
#############
# CACHING: df.final.up.down.directional processing
cache_file_directional <- "../data/df_final_up_down_mid_mid_directional.rds"

if (file.exists(cache_file_directional)) {
  message("Loading cached directional filtering results from: ", cache_file_directional)
  df.final.up.down.directional <- readRDS(cache_file_directional)
} else {
  message("Processing directional filtering...")

  # For UP loops: filter genes where gene end <= loop end (x2)
  # For DOWN loops: filter genes where gene start >= loop start (y1)

  # Create filtering vectors for components
  component_chr <- as.character(seqnames(df.gene_tss_and_pro.GR))
  component_start <- start(df.gene_tss_and_pro.GR)
  component_end <- end(df.gene_tss_and_pro.GR)

  # Process UP loops with directional filtering
  up_results_list <- lapply(seq_along(df.DISTINCT.loop.deep.sample.all.lt.2mb.UP.point.GR), function(i) {
    loop_chr <- as.character(seqnames(df.DISTINCT.loop.deep.sample.all.lt.2mb.UP.point.GR[i]))
    loop_mid <- as.numeric(mcols(df.DISTINCT.loop.deep.sample.all.lt.2mb.UP.point.GR)$mid_loop[i])

    # Filter genes: same chr AND gene end <= loop x2
    valid_component_up_idx <- which(component_chr == loop_chr & component_start <= loop_mid)

    if (length(valid_component_up_idx) > 0) {
      component_up_subset <- df.gene_tss_and_pro.GR[valid_component_up_idx]
      hits <- distanceToNearest(df.DISTINCT.loop.deep.sample.all.lt.2mb.UP.point.GR[i],
        component_up_subset,
        select = "all"
      ) # multiple hits

      if (length(hits) > 0) {
        # Convert to data frame with original indices
        tibble(
          loop_idx = i,
          component_idx = valid_component_up_idx[subjectHits(hits)],
          distance = mcols(hits)$distance
        )
      } else {
        NULL
      }
    } else {
      NULL
    }
  })

  # Combine UP results
  up_results_df <- bind_rows(up_results_list)
  up_results_df # 36,490  // mid mid: 31,045

  # Process DOWN loops with directional filtering
  down_results_list <- lapply(seq_along(df.DISTINCT.loop.deep.sample.all.lt.2mb.DOWN.point.GR), function(i) {
    loop_chr <- as.character(seqnames(df.DISTINCT.loop.deep.sample.all.lt.2mb.DOWN.point.GR[i]))
    loop_mid <- as.numeric(mcols(df.DISTINCT.loop.deep.sample.all.lt.2mb.DOWN.point.GR)$mid_loop[i])

    # Filter genes: same chr AND gene start >= loop y1
    valid_component_down_idx <- which(component_chr == loop_chr & component_end >= loop_mid)

    if (length(valid_component_down_idx) > 0) {
      component_down_subset <- df.gene_tss_and_pro.GR[valid_component_down_idx]
      hits <- distanceToNearest(df.DISTINCT.loop.deep.sample.all.lt.2mb.DOWN.point.GR[i],
        component_down_subset,
        select = "all"
      ) # multiple hits

      if (length(hits) > 0) {
        # Convert to data frame with original indices
        tibble(
          loop_idx = i,
          component_idx = valid_component_down_idx[subjectHits(hits)],
          distance = mcols(hits)$distance
        )
      } else {
        NULL
      }
    } else {
      NULL
    }
  })

  # Combine DOWN results
  down_results_df <- bind_rows(down_results_list)
  down_results_df

  message("UP hits: ", nrow(up_results_df), ", DOWN hits: ", nrow(down_results_df))
  # Midpoint in anchor (select all) and midpoint in loop- UP hits: 31,045, DOWN hits: 31,066
  # Convert indices to actual data for UP results
  up_results_full <- up_results_df %>%
    mutate(
      # Loop information from UP GRanges
      loop.id = mcols(df.DISTINCT.loop.deep.sample.all.lt.2mb.UP.point.GR)$loop.id[loop_idx],
      loop.mid = mcols(df.DISTINCT.loop.deep.sample.all.lt.2mb.UP.point.GR)$mid_loop[loop_idx],
      loop_chr = as.character(seqnames(df.DISTINCT.loop.deep.sample.all.lt.2mb.UP.point.GR[loop_idx])),
      loop_start = BiocGenerics::start(df.DISTINCT.loop.deep.sample.all.lt.2mb.UP.point.GR[loop_idx]),
      loop_end = BiocGenerics::end(df.DISTINCT.loop.deep.sample.all.lt.2mb.UP.point.GR[loop_idx]),

      # Gene information from gene GRanges
      gene_id = mcols(df.gene_tss_and_pro.GR)$gene_id[component_idx],
      gene_name = mcols(df.gene_tss_and_pro.GR)$gene_name[component_idx],
      component_id = mcols(df.gene_tss_and_pro.GR)$component_id[component_idx],
      component = mcols(df.gene_tss_and_pro.GR)$component[component_idx],
      component_chr = as.character(seqnames(df.gene_tss_and_pro.GR[component_idx])),
      component_start = BiocGenerics::start(df.gene_tss_and_pro.GR[component_idx]),
      component_end = BiocGenerics::end(df.gene_tss_and_pro.GR[component_idx]),

      # Gene coordinates from metadata
      gene_chr = mcols(df.gene_tss_and_pro.GR)$gene_chr[component_idx],
      gene_start = as.numeric(mcols(df.gene_tss_and_pro.GR)$gene_start[component_idx]),
      gene_end = as.numeric(mcols(df.gene_tss_and_pro.GR)$gene_end[component_idx]),

      # Add WHERE column
      WHERE = "UP"
    ) %>%
    dplyr::select(-loop_idx, -component_idx) # Remove index columns

  # Convert indices to actual data for DOWN results
  down_results_full <- down_results_df %>%
    mutate(
      # Loop information from DOWN GRanges
      loop.id = mcols(df.DISTINCT.loop.deep.sample.all.lt.2mb.DOWN.point.GR)$loop.id[loop_idx],
      loop.mid = mcols(df.DISTINCT.loop.deep.sample.all.lt.2mb.DOWN.point.GR)$mid_loop[loop_idx],
      loop_chr = as.character(seqnames(df.DISTINCT.loop.deep.sample.all.lt.2mb.DOWN.point.GR[loop_idx])),
      loop_start = BiocGenerics::start(df.DISTINCT.loop.deep.sample.all.lt.2mb.DOWN.point.GR[loop_idx]),
      loop_end = BiocGenerics::end(df.DISTINCT.loop.deep.sample.all.lt.2mb.DOWN.point.GR[loop_idx]),

      # Gene information from gene GRanges
      gene_id = mcols(df.gene_tss_and_pro.GR)$gene_id[component_idx],
      gene_name = mcols(df.gene_tss_and_pro.GR)$gene_name[component_idx],
      component_id = mcols(df.gene_tss_and_pro.GR)$component_id[component_idx],
      component = mcols(df.gene_tss_and_pro.GR)$component[component_idx],
      component_chr = as.character(seqnames(df.gene_tss_and_pro.GR[component_idx])),
      component_start = BiocGenerics::start(df.gene_tss_and_pro.GR[component_idx]),
      component_end = BiocGenerics::end(df.gene_tss_and_pro.GR[component_idx]),

      # Gene coordinates from metadata
      gene_chr = mcols(df.gene_tss_and_pro.GR)$gene_chr[component_idx],
      gene_start = as.numeric(mcols(df.gene_tss_and_pro.GR)$gene_start[component_idx]),
      gene_end = as.numeric(mcols(df.gene_tss_and_pro.GR)$gene_end[component_idx]),

      # Add WHERE column
      WHERE = "DOWN"
    ) %>%
    dplyr::select(-loop_idx, -component_idx) # Remove index columns

  # Combine UP and DOWN results
  df.final.up.down.directional <- bind_rows(up_results_full, down_results_full)

  # Save to cache
  saveRDS(df.final.up.down.directional, cache_file_directional)
}

df.final.up.down.directional %>% dim() # 72,772// point: 62,070/ mid mid:62,111 (select = all): must filter multiple cases from an anchor
df.final.up.down.directional %>% head(2)
df.final.up.down.directional.point <- df.final.up.down.directional %>%
  group_by(loop.id, WHERE) %>%
  filter(distance == min(distance)) %>% # multiple hits with the same distance in a loop and an anchor
  add_count(name = "tie.count.loop.where") %>%
  ungroup()

df.final.up.down.directional.point %>% dim() # 61,974 // 61,936// mid mid: 62,111 with_ties = TRUE// FYI: 62038 = 31019 * 2
df.final.up.down.directional.point %>% head(2)

########################################
# Categorizing and filtering steps
########################################
# 1. only one case in a loop and an anchor
df.final.up.down.directional.point.one.in.loop.where <- df.final.up.down.directional.point %>%
  filter(tie.count.loop.where == 1) # 61,838

# 2. multiple cases in a loop and an anchor : filter(tie.count.loop.where > 1): 273
# for df.final.combined.using.labelings from df.final.up.down.directional.point.multi.in.loop.where
df.final.up.down.directional.point.multi.in.loop.where <- df.final.up.down.directional.point %>%
  filter(tie.count.loop.where > 1) # 273
df.final.up.down.directional.point.multi.in.loop.where %>%
  # count(loop.id) # 136
  count(tie.count.loop.where) # 273 (2:270 + 3(3))
#   tie.count.loop.where     n
# 1                    2   270
# 2                    3     3

# step for filtering only cases of PRO and TSS both in an anchor with the same gene_name (only removing TSS cases)
df.final.up.down.directional.point.multi.in.loop.where.two.labelings <- df.final.up.down.directional.point.multi.in.loop.where %>% # head(2)
  filter(!is.na(gene_chr)) %>% # view() # 245 (28: is.na)
  add_count(loop.id, WHERE, name = "multi.1st.filter") %>% # multi.1st.filter: is.na 처리하고 재통계 (multi.1st.filter > 1이 문제들)
  # count(multi.1st.filter) %>%
  # view()
  # count(multi.1st.filter) %>% # AFTER CLEANING NA
  group_by(loop.id, WHERE, gene_name) %>% ###### REQUIRED: gene_name because we use both if genes are different for cases in an anchor in a loop
  # filter(is.na(gene_name)) %>% view() # 0
  mutate(
    has_pro_tss = all(c("pro", "tss") %in% component), # true or false: depending on pro와 tss 둘 다 있는지, (pro/pro) 인지, (tss/tss) 인지 확인
    keep_row = case_when(
      multi.1st.filter == 1 ~ TRUE, # These are FALSE by has_pro_tss
      has_pro_tss & component == "pro" ~ TRUE, # priority: pro from pro/tss both
      has_pro_tss & component == "tss" ~ FALSE, # tss is removed
      TRUE ~ TRUE # keep all
    )
  ) %>%
  ungroup() %>%
  filter(keep_row) %>%
  add_count(loop.id, WHERE, name = "multi.2nd.filter")
#   RESULT at 461.                                  RESULT at end of code
#   multi.1st.filter     n                          multi.1st.filter     n
# 1                1    24 : by NA in gene_chr    1                1    24 # protected by line 471
# 2                2   218                        2                2   147 # just removed by this step: 71
# 3                3     3                        3                3     3

df.final.up.down.directional.point.multi.in.loop.where.two.labelings %>% count(multi.1st.filter)
df.final.up.down.directional.point.multi.in.loop.where.two.labelings %>% dim() # 174
df.final.up.down.directional.point.multi.in.loop.where.two.labelings %>% count(has_pro_tss)
df.final.up.down.directional.point.multi.in.loop.where.two.labelings %>% count(multi.2nd.filter)
df.final.up.down.directional.point.multi.in.loop.where.two.labelings %>% head(2)

df.final.up.down.directional.point.multi.in.loop.where.two.labelings.selected <- df.final.up.down.directional.point.multi.in.loop.where.two.labelings %>% # 95/174
  filter(multi.2nd.filter == 1)
df.final.up.down.directional.point.multi.in.loop.where.two.labelings.selected %>% dim() # 95
df.final.up.down.directional.point.multi.in.loop.where.two.labelings.selected %>% head(2)

# 3. AFTER PROMOTER PRIORITY filteration, 2 rows from loop.id and WHERE
df.final.up.down.directional.point.multi.in.loop.where.two.labelings.ongoing <- df.final.up.down.directional.point.multi.in.loop.where.two.labelings %>%
  filter(multi.2nd.filter > 1) # 79
df.final.up.down.directional.point.multi.in.loop.where.two.labelings.ongoing %>% dim() # 79
df.final.up.down.directional.point.multi.in.loop.where.two.labelings.ongoing %>% head(2)

df.final.up.down.directional.point.multi.in.loop.where.three.labelings <- df.final.up.down.directional.point.multi.in.loop.where.two.labelings.ongoing %>%
  arrange(loop.id, WHERE) %>%
  dplyr::select(multi.2nd.filter, gene_name, component, loop.id, everything()) %>%
  add_count(loop.id, WHERE, name = "multi.3rd.filter")

df.final.up.down.directional.point.multi.in.loop.where.three.labelings %>% dim()
df.final.up.down.directional.point.multi.in.loop.where.three.labelings %>% head(2)

# gene_name 필터링
df.final.up.down.directional.point.multi.in.loop.where.three.labelings.actual.same.gene.selected <- df.final.up.down.directional.point.multi.in.loop.where.three.labelings %>%
  filter(multi.3rd.filter > 1) %>% # 79
  mutate(gene_id_trimmed = sapply(strsplit(gene_id, ":"), function(x) paste(x[-5], collapse = ":"))) %>%
  group_by(loop.id, WHERE) %>%
  mutate(all_same_trimmed = n_distinct(gene_id_trimmed) == 1) %>%
  ungroup() %>%
  # count(all_same_trimmed) %>% view() # FALSE: 30, TRUE: 49
  filter(all_same_trimmed) %>% # # mutaul exclusive to L533 :::::::::: 49
  group_by(loop.id, WHERE) %>%
  filter(row_number() == 1) %>%
  ungroup()

df.final.up.down.directional.point.multi.in.loop.where.three.labelings.actual.same.gene.selected %>% dim() # 24/49
df.final.up.down.directional.point.multi.in.loop.where.three.labelings.actual.same.gene.selected %>% head(2)

df.final.up.down.directional.point.multi.in.loop.where.three.labelings.exon.number.and.mannually.selected <- df.final.up.down.directional.point.multi.in.loop.where.three.labelings %>%
  filter(multi.3rd.filter > 1) %>%
  mutate(gene_id_trimmed = sapply(strsplit(gene_id, ":"), function(x) paste(x[-5], collapse = ":"))) %>%
  group_by(loop.id, WHERE) %>%
  mutate(all_same_trimmed = n_distinct(gene_id_trimmed) == 1) %>%
  ungroup() %>%
  filter(!all_same_trimmed) %>% # mutaul exclusive to L519 :::::::::: 30
  # view() # 30
  mutate(
    gene_name_from_id = sapply(strsplit(gene_id_trimmed, ":"), function(x) x[length(x) - 1]),
    last_element = as.numeric(sapply(strsplit(gene_id_trimmed, ":"), function(x) x[length(x)]))
  ) %>%
  # gene_name이 같을 때 last_element가 가장 큰 것 선택
  group_by(loop.id, WHERE, gene_name_from_id) %>%
  slice_max(last_element, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(-gene_id_trimmed, -all_same_trimmed, -gene_name_from_id, -last_element) %>%
  # view() # Now:::::::::::::::::::::::::::::::::::::::::::::::: 21 (removed 9)
  # count(loop.id, WHERE) %>% count(n) %>% view() # 1: 9/ 2: 6
  filter(!str_detect(gene_name, "LOC|Gzmbl2")) %>% # removing LOC100911295, LOC120100814// Gzmbl2
  view() # now:::::::::::::::::::::::::::::::::::::::::::::::::: 18 (removed 3)
# count(loop.id, WHERE) %>% count(n) %>% view() # 1: 12/ 2: 3. (L:574)

# group_by(loop.id, WHERE) %>%
# filter(n() > 1) %>% # 그룹 내 row가 2개인 것만
# view()

df.final.up.down.directional.point.multi.in.loop.where.three.labelings.exon.number.and.mannually.selected %>% dim() # 18/30
df.final.up.down.directional.point.multi.in.loop.where.three.labelings.exon.number.and.mannually.selected %>% head(2)

df.final.combined.using.labelings <- bind_rows( # 137 = 95 + 24 + 18
  # df.final.up.down.directional.point.one.in.loop.where,                                    # 61,838
  df.final.up.down.directional.point.multi.in.loop.where.two.labelings.selected, # 95
  df.final.up.down.directional.point.multi.in.loop.where.three.labelings.actual.same.gene.selected, # 24
  df.final.up.down.directional.point.multi.in.loop.where.three.labelings.exon.number.and.mannually.selected # 18
)

df.final.up.down.directional.point.final <- bind_rows(df.final.up.down.directional.point.one.in.loop.where, df.final.combined.using.labelings)

df.final.up.down.directional.point.final %>% dim()
df.final.up.down.directional.point.final %>% head(3)
# 예상: 61,838 + 137 = 61,975

df.final.up.down.directional.point.final %>%
  add_count(loop.id, WHERE, name = "final_count") %>%
  count(final_count)
# final_count count
#1           1 61969
#2           2     6 (L:548)

df.final.up.down.directional.point.final %>% dim() # 62,070 // 61,936// mid mid: 61,975/ FYI: 62038 = 31019 * 2

# validity of the gene’s genomic position based on TSS or promoters selected by the nearestToDistance
df.final.up.down.directional.point.decision <- df.final.up.down.directional.point.final %>%
  separate(
    col = loop.id,
    into = c("chr1", "x1", "x2", "chr2", "y1", "y2", "resolution"),
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(
    UP = case_when(
      (WHERE == "UP") & (gene_end < loop.mid) & (gene_start > as.numeric(x1)) ~ "OK",
      WHERE == "UP" ~ "FAIL",
      TRUE ~ NA_character_
    ),
    DOWN = case_when(
      (WHERE == "DOWN") & (gene_start > loop.mid) & (gene_end < as.numeric(y2)) ~ "OK",
      WHERE == "DOWN" ~ "FAIL",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(loop.id) %>%
  mutate(
    # loop.id당 UP에 "OK"가 있는지 확인
    has_UP_OK = any(UP == "OK", na.rm = TRUE),
    # loop.id당 DOWN에 "OK"가 있는지 확인
    has_DOWN_OK = any(DOWN == "OK", na.rm = TRUE),
    classification = case_when(
      has_UP_OK & has_DOWN_OK ~ "Both_OK", # UP과 DOWN 둘 다 OK 존재
      has_UP_OK | has_DOWN_OK ~ "One_OK", # UP 또는 DOWN 중 하나만 OK 존재
      TRUE ~ "Both_FAIL" # 둘 다 OK 없음
    )
  ) %>%
  ungroup()

df.final.up.down.directional.point.decision %>%
  count(classification)

#   classification     n
# 1 Both_FAIL      30497| Both_FAIL      30498
# 2 Both_OK         7606| Both_OK         7606
# 3 One_OK         23872| One_OK         23871

###################################
# 1. ONE OK case in either anchor
###################################

df.final.up.down.directional.point.decision %>%
  distinct(loop.id) %>%
  count()
df.final.up.down.directional.point.decision %>%
  filter(classification == "One_OK") %>% #
  # distinct(gene_id) %>% # 7,894
  distinct(loop.id) %>% # 11,935/31,019
  count()

# stats for distance with ALL loops
stats_decision <- df.final.up.down.directional.point.decision %>%
  # filter(classification == "One_OK") %>%
  # group_by(WHERE) %>%
  summarise(
    Q1 = quantile(distance, 0.25),
    Median = median(distance),
    Q3 = quantile(distance, 0.75)
  )
stats_decision
# all loops           
#     Q1 Median    Q3 
# 1 7400.  23371 66691

df.final.up.down.directional.point.decision %>% head(2)

# all loops boxplot for distance
df.final.up.down.directional.point.all.in.one.decision <- df.final.up.down.directional.point.decision %>%
  ggplot(aes(y = distance)) +
  geom_boxplot(fill = "skyblue") +
  annotate("text", x = 1, y = stats_decision$Q1, label = paste0("Q1=", round(stats_decision$Q1, 2)), vjust = -0.5, hjust = 1, color = "blue") +
  annotate("text", x = 1, y = stats_decision$Median, label = paste0("Median=", round(stats_decision$Median, 2)), vjust = -0.5, hjust = 1, color = "red") +
  annotate("text", x = 1, y = stats_decision$Q3, label = paste0("Q3=", round(stats_decision$Q3, 2)), vjust = -0.5, hjust = 1, color = "blue") +
  labs(
    title = "Boxplot of Distance (All Data)",
    y = "Distance"
  ) +
  scale_y_log10() +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
df.final.up.down.directional.point.all.in.one.decision

# boxplot by classification
distribution.all.cases.distance <- df.final.up.down.directional.point.decision %>%
  ggplot(aes(x = classification, y = distance, fill = classification)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "Boxplot of Distance by Classification",
    x = "Classification",
    y = "Distance"
  ) +
  scale_y_log10() +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
distribution.all.cases.distance

# all loops           
#     Q1 Median    Q3 
# 1 7400.  23371 66691

df.final.up.down.directional.point.decision.one.OK.filtered <- df.final.up.down.directional.point.decision %>%
  filter(classification == "One_OK") %>% # 23,874
  group_by(loop.id) %>%
  filter(UP == "OK" | DOWN == "OK") %>% # UP 또는 DOWN 중 OK인 행만 선택
  ungroup()

df.final.up.down.directional.point.decision.one.OK.filtered %>%
  distinct(loop.id) %>%
  count() # 11,935
df.final.up.down.directional.point.decision.one.OK.filtered %>% head(2) # 11,937 // 11,935

approach_2nd_analyze_loops_by_threshold(
  df.final.up.down.directional.point.decision.one.OK.filtered %>%
    dplyr::rename(gene_id_id = gene_id) %>%
    mutate(gene_id = str_split_n(str_split_n(component_id, ":", 6), "\\|", 1)),
  threshold_distance = stats_decision$Q3,
  top_n_genes = 70,
  print_top_n = 50
) # utils_functions.R

# distribution.OK.distance.q3 %>% dim()
# distribution.OK.distance.q3 %>% distinct(loop.id) %>% count()
# distribution.OK.distance.q3 %>% head(2)

###################################
# 2. TWO OK case in either anchor
###################################

# two cases in both anchors
df.final.up.down.directional.point.decision.both.OK <- df.final.up.down.directional.point.decision %>%
  filter(classification == "Both_OK") %>% # 3,803/31,019
  # distinct(loop.id) %>%
  # count() %>% 
  view()

df.final.up.down.directional.point.decision.both.OK %>% head(2)
df.final.up.down.directional.point.decision.both.OK %>% dim() # 7,606 //// 7606/3,801 * 2

###################################
# 2-1. distance difference between UP and DOWN in BOTH OK case
###################################

dist.diff.decision.both.OK <- df.final.up.down.directional.point.decision.both.OK %>%
  filter(WHERE %in% c("UP", "DOWN")) %>%
  group_by(loop.id) %>%
  filter(n_distinct(WHERE) == 2) %>% # 3810
  summarise(
    dist_up = distance[WHERE == "UP"],
    dist_down = distance[WHERE == "DOWN"],
    dist_diff = abs(dist_up - dist_down),
    .groups = "drop"
  )

dist.diff.decision.both.OK

stats_dist_diff <- dist.diff.decision.both.OK %>%
  summarise(
    Q1 = quantile(dist_diff, 0.25, na.rm = TRUE),
    Q3 = quantile(dist_diff, 0.75, na.rm = TRUE)
  )

stats_dist_diff

# Boxplot for distance difference
ggplot(dist.diff.decision.both.OK, aes(x = "", y = dist_diff)) +
  geom_boxplot(fill = "#6CABDD", color = "black") +
  scale_y_log10() +
  labs(y = "Absolute distance difference (UP vs DOWN)", x = NULL) +
  geom_hline(yintercept = stats_dist_diff$Q1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = stats_dist_diff$Q3, linetype = "dashed", color = "blue") +
  annotate("text",
    x = 0.8, y = stats_dist_diff$Q1, label = paste0("Q1=", stats_dist_diff$Q1),
    hjust = 1, vjust = -0.5, color = "red"
  ) +
  annotate("text",
    x = 0.8, y = stats_dist_diff$Q3, label = paste0("Q3=", stats_dist_diff$Q3),
    hjust = 1, vjust = -0.5, color = "blue"
  )


df.final.up.down.directional.point.decision.both.OK %>%
  count(distance)

##################################
# 큰 distance 제거하고 하나만 남기기
##################################
df.final.up.down.directional.point.decision.both.OK.filtered <- df.final.up.down.directional.point.decision.both.OK %>%
  group_by(loop.id) %>%
  filter(distance == min(distance)) %>% # 최소 distance만 남김
  ungroup()

df.final.up.down.directional.point.decision.both.OK.filtered %>% dim() # 3,803
df.final.up.down.directional.point.decision.both.OK.filtered %>% head(2)

# 동일 distance가 존재하는 loop.id 파악하기
df.final.up.down.directional.point.decision.both.OK.tied <- df.final.up.down.directional.point.decision.both.OK %>%
  group_by(loop.id) %>%
  filter(n_distinct(distance) == 1, n() > 1) %>% # 같은 loop.id 내 distance가 모두 동일
  ungroup()

df.final.up.down.directional.point.decision.both.OK.tied # 0

################################
# WHAT? NO tied cases? 0 - 0
################################
df.final.up.down.directional.point.decision.both.OK.filtered %>%
  count(loop.id) %>%
  filter(n != 1) %>%
  count(n) # n(1) = 0

stats_decision$Q3 # 66,892
dist.diff.decision.both.OK
# 히스토그램 + Q3선
ggplot(df.final.up.down.directional.point.decision.both.OK, aes(x = distance)) +
  geom_histogram(binwidth = 1000, fill = "#6CABDD", color = "black") +
  labs(x = "Absolute distance difference (UP vs DOWN)", y = "Count") +
  geom_vline(xintercept = stats_decision$Q3, linetype = "dashed", color = "red") +
  annotate("text",
    x = stats_decision$Q3, y = 0, label = paste0("Q3=", round(stats_decision$Q3, 0)),
    angle = 90, vjust = -0.5, hjust = 0, color = "red"
  )

df.final.up.down.directional.point.decision.both.OK.filtered %>%
  summarise(
    perc_below_Q3 = mean(distance < stats_decision$Q3) * 100
  )

df.final.up.down.directional.point.decision.both.OK.filtered %>% # dim() 3801
  filter(distance < stats_decision$Q3) %>% dim() # 3730/3801
# perc_below_Q3
# 1          98.1

# Log scale histogram
ggplot(dist.diff.decision.both.OK, aes(x = dist_diff)) +
  geom_histogram(fill = "#6CABDD", color = "black") +
  scale_x_log10() +
  labs(x = "Absolute distance difference (UP vs DOWN, log scale)", y = "Count")

###################################
# 2. One OK + TWO OK case in either anchor
###################################

df.final.up.down.directional.point.decision %>% dim() # 61,975
df.final.up.down.directional.point.decision %>% head(2)
df.final.up.down.directional.point.decision %>% count(classification)
# 1 Both_FAIL      30498|||||||||||            1 Both_FAIL      30498  1 Both_FAIL      30497
# 2 Both_OK         7606|||||||||||            2 Both_OK         7602  2 Both_OK         7606
# 3 One_OK         23871|||||||||||            3 One_OK         23874  3 One_OK         23872

# CACHING: df.final.up.down.directional.point.decision.COMBINED.OK.filtered.lt.Q3
cache_file_combined_ok_filtered <- "../data/df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3.rds"

if (!file.exists(cache_file_combined_ok_filtered)) {
  message("Saving df.final.up.down.directional.point.decision.COMBINED.OK.filtered.lt.Q3 to cache: ", cache_file_combined_ok_filtered)

  df.final.up.down.directional.point.decision.COMBINED.OK.filtered.lt.Q3 <- bind_rows(
    df.final.up.down.directional.point.decision.one.OK.filtered,
    df.final.up.down.directional.point.decision.both.OK.filtered
  ) %>%
    filter(distance < stats_decision$Q3) # 3,730

  df.final.up.down.directional.point.decision.COMBINED.OK.filtered.lt.Q3 %>% dim() # 14210
  df.final.up.down.directional.point.decision.COMBINED.OK.filtered.lt.Q3 %>% head(2)
  df.final.up.down.directional.point.decision.COMBINED.OK.filtered.lt.Q3 %>% count(loop.id) # 14,210 PASS!

  df.final.up.down.directional.point.decision.COMBINED.OK.filtered.lt.Q3 %>%
    dplyr::rename(gene_id_id = gene_id) %>%
    mutate(gene_id = str_split_n(str_split_n(component_id, ":", 6), "\\|", 1)) %>%
    count(gene_id, sort = TRUE) # desc(n)
    # slice_max(n, n = top_n_genes) %>%

  saveRDS(df.final.up.down.directional.point.decision.COMBINED.OK.filtered.lt.Q3, cache_file_combined_ok_filtered)
} else {
  message("Cache file already exists: ", cache_file_combined_ok_filtered)
}

df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3 <- readRDS(cache_file_combined_ok_filtered)
df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3 %>% dim() # 14210
df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3 %>% head(2)

approach_2nd_analyze_loops_by_threshold(
  df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3 %>%
    dplyr::rename(gene_id_id = gene_id) %>%
    mutate(gene_id = str_split_n(str_split_n(component_id, ":", 6), "\\|", 1)),
  threshold_distance = stats_decision$Q3,
  top_n_genes = 30,
  print_top_n = 50
) # utils_functions.R

df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3 %>%
  count(component)
# 1 pro        4776  1 pro        4569|||||||| 1 pro        4569
# 2 tss       10891  2 tss        9641|||||||| 2 tss        9641

######################################################################
######################################################################

# final
final.loops.from.tss.step <- df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3 %>%
  filter(component == "tss")
final.loops.from.tss.step %>% dim() # 9,641

final.loops.from.promoter.step <- df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3 %>%
  filter(component == "pro")
final.loops.from.promoter.step %>% dim() # 4,569

df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3 %>% head(3)
df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3 %>% dim() # 23524// 14210
df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3 %>% count(loop.id, component) %>% # view()
  count(n)
# 1 1 9177 // 14210
# 2 2 1983

df_single_loop <- df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3 %>% # 9177
  group_by(loop.id) %>%
  filter(n() == 1) %>%
  ungroup()

# WHERE distribution
where_dist <- df_single_loop %>%
  count(WHERE)
where_dist
# 1 DOWN   4485  // 1 DOWN   6943
# 2 UP     4692  // 2 UP     7267

# component distribution
component_dist <- df_single_loop %>%
  count(component)
component_dist
# 1 pro        4035  // 1 pro        4569
# 2 tss        5142  // 2 tss        9641

################################################################################################
################################################################################################
# threshold value visualization
################################################################################################
################################################################################################
# 4-1. CTCF
df.ctcf.counts %>% head(2)
df.ctcf.counts %>% dim()
df.ctcf.counts %>% count(resolution, WHERE)
df.ctcf.counts %>%
  count(resolution, WHERE) %>%
  group_by(resolution) %>%
  summarise(total = sum(n), .groups = "drop")
df.ctcf.counts %>%
  group_by(resolution, WHERE) %>%
  summarise(max_ctcf = max(ctcf_count), .groups = "drop")

df.ctcf.filtered <- plotting_and_filtering_summary( # utils_functions.R
  df_counts = df.ctcf.counts,
  count_col = "ctcf_count",
  output_prefix = "ctcf"
)
df.ctcf.counts %>% head()
ctcf.stats.by.resolution

# 2^2.5
# [1] 5.656854
# log2(5.656854)
# [1] 2.5

df.loops.above.ctcf.threshold <- df.ctcf.counts %>%
  filter(ctcf_count >= 6) %>%
  group_by(loop.id) %>%
  filter(n_distinct(WHERE) == 2) %>%
  ungroup()

df.loops.above.ctcf.threshold %>% dim() # 50,536
df.loops.above.ctcf.threshold %>% head(2)

final.loops.from.ctcf.step <- df.loops.above.ctcf.threshold %>%
  distinct(loop.id) %>%
  mutate(end.distance = str_split_n(loop.id, "_", 7)) %>% # utils_functions.R
  mutate(resolution = case_when(
    end.distance == 5000 ~ "5K",
    end.distance == 10000 ~ "10K",
    end.distance == 25000 ~ "25K",
    TRUE ~ NA
  ))

# %>%
#   mutate(loop.id = str_remove(loop.id, "_[^_]+$"))
final.loops.from.ctcf.step %>% dim() # 18,327 + 10 | sub.4 any, threshold > 6 : 25,620| lt2mb: 25,268
final.loops.from.ctcf.step %>% head(2)
final.loops.from.ctcf.step %>% count(resolution)
# sub.4
# any:
# new threshold
# lt2mb              # resolution     n
# 1 10K         9375 # 1 10K         9420
# 2 25K        11750 # 2 25K        12049
# 3 5K          4143 # 3 5K          4151

####################################
####################################
# Venn Diagram
####################################
####################################
final.loops.from.promoter.step$loop.id # 4569
final.loops.from.tss.step$loop.id # 9641

venn_plot_submission <- create_venn_plot(final.loops.from.ctcf.step, final.loops.from.promoter.step, final.loops.from.tss.step, "CTCF") # utils_functions.R
venn_plot_submission

saving_plot_dual( # utils_functions.R
  plot_obj = venn_plot_submission,
  filename_base = "ctcf_vs_promoter_tss_venn_diagrams_approach2_lt_2mb_final",
  output_dir = "./figures/submission/lt2mb",
  scale_x = 1,
  scale_y = 1
)

# initial: 4569 + 9641 = 14210
# final including CTCF: 3973 + 8261 = 12234
####################################
# functional loop extraction
####################################

final.loops.from.ctcf.step %>% dim() # sub.4 any, either 6: 25,620// span ENSEMBL 25268
final.loops.from.tss.step %>% dim() # mid mid: 9,641
final.loops.from.promoter.step %>% dim() # mid mid: 4,569

overlapping_loops <- extract_overlapping_loops( # utils_functions.R
  ctcf_data = final.loops.from.ctcf.step,
  promoter_data = final.loops.from.promoter.step,
  tss_data = final.loops.from.tss.step
)

overlapping_loops$ctcf_promoter_only # PASS// 6648 // final: 3,973
overlapping_loops$ctcf_tss_only # PASS// 8613 // final: 8261
overlapping_loops$ctcf_promoter_tss_overlap # PASS// 3423// final: 0

df.ctcf.promoter.only.loop <- data.frame(loop.id = overlapping_loops$ctcf_promoter_only, category = "CP", stringsAsFactors = FALSE)
df.ctcf.tss.only.loop <- data.frame(loop.id = overlapping_loops$ctcf_tss_only, category = "CT", stringsAsFactors = FALSE)
df.ctcf.promoter.tss.overlap.loop <- data.frame(loop.id = overlapping_loops$ctcf_promoter_tss_overlap, category = "CPT", stringsAsFactors = FALSE)

df.ctcf.promoter.only.loop %>% dim() # 5497//3973
df.ctcf.tss.only.loop %>% dim() # 7799// 8261
df.ctcf.promoter.tss.overlap.loop %>% dim() # 2302//0

# sub.4 dedups: 11891 = 679 + 3096 + 8116
df.final.loop <- bind_rows(df.ctcf.promoter.only.loop, df.ctcf.tss.only.loop)

df.final.loop # sub.4 any, either 6// 12234/31019(lt2mb) (0.3944034)
df.final.loop %>% head()
df.final.loop %>% dim() # 15598: 5656 + 7980 + 1962// final: 12,234
df.final.loop %>% count(category) # CP 3973 CT 8261

df.final.loop %>%
  mutate(resolution = str_split_n(loop.id, "_", 7)) %>%
  count(resolution) # utils_functions.R
#  resolution    n                      resolution       n  resolution    n     resolution    n
# 1      10000 4369                      1      10000  86751      10000 5656     1      10000 6968
# 2      25000 6104                      2      25000 111442      25000 7980     2      25000 8617
# 3       5000 1761                      3       5000  38213       5000 1962     3       5000 3099

save(df.final.loop, file = "./figures/submission/lt2mb/df_final_loop_sub.4.any.lt2mb.ENSEMBL.mid.mid.final.filter.rda")
write.csv(df.final.loop, file = "./figures/submission/lt2mb/df_final_loop_sub.4.any.lt2mb.ENSEMBL.mid.mid.final.filter.csv", row.names = FALSE)

df.final.loop %>% dim() # 23640// 18681// 18684// 12234
df.final.loop %>% head(2)
df.DISTINCT.loop.deep.sample.all %>% head()

df.final.DISTINCT.loop.joined <- df.final.loop %>%
  left_join(df.DISTINCT.loop.deep.sample.all, by = "loop.id")
df.final.DISTINCT.loop.joined %>% dim()
df.final.DISTINCT.loop.joined %>% head(2)

df.final.loop %>% distinct(loop.id) # 23460// 18681// 18684// 12234

##########################################################
# ideogram                                        Figure 4
##########################################################
df.chromosome.data %>% head(3)
df.chromosome.data <- df.chromosome.data %>% mutate(CE_start = NA, CE_end = NA)
df.chromosome.data

# CTCF for ideogram: df_ctcf
df.DISTINCT.fimo.2nd.trial.ctcf %>% head()
df.DISTINCT.fimo.2nd.trial.ctcf # 3191859 (sub .4)

df_ctcf_ideogram <- df.DISTINCT.fimo.2nd.trial.ctcf %>%
  dplyr::select(chr, start, end)

df_ctcf_ideogram

bin_size <- 1000000

df_binned <- df_ctcf_ideogram %>%
  mutate(
    StartBin = floor(start / bin_size) * bin_size + 1,
    EndBin = StartBin + bin_size - 1
  )

# Step 4: counting components per bin
gene_density <- df_binned %>%
  group_by(Chr = chr, Start = StartBin, End = EndBin) %>%
  summarise(Value = n(), .groups = "drop") %>%
  arrange(Chr, Start)

head(gene_density)

getwd()
ctcf_density <- process_feature_bins(df_ctcf_ideogram, df.chromosome.data, bin_size, label = "CTCF") %>% dplyr::select(-last_col()) # color = "#E41A1C" # utils_functions.R
ctcf_density

ideogram(
  karyotype = df.chromosome.data %>% dplyr::rename(Chr = chr, Start = start, End = end) %>% mutate(Chr = str_remove(Chr, "chr")),
  overlaid = ctcf_density %>% dplyr::rename(Chr = chr, Start = start, End = end) %>% mutate(Chr = str_remove(Chr, "chr")),
)

base_dir <- getwd()

# filenames
svg_file <- file.path(base_dir, "chromosome.svg")
pdf_file <- file.path(base_dir, "chromosome.pdf")
png_file <- file.path(base_dir, "chromosome.png")

# SVG → PDF
rsvg_pdf(svg_file, pdf_file)

# SVG → PNG (dpi=300 for high resolution)
rsvg_png(svg_file, png_file)

########################
########################
# circos plot for loops
########################
########################
df.final.loop %>% head()
df.circos.final.loop <- df.final.loop %>%
  separate(loop.id, into = c("chr1", "x1", "x2", "chr2", "y1", "y2", "end.distance"), sep = "_", convert = TRUE, remove = FALSE) %>%
  relocate(chr1, x1, x2, chr2, y1, y2, end.distance, .after = loop.id) %>%
  mutate(y12 = (y2 + y1) / 2, x12 = (x2 + x1) / 2) %>%
  mutate(distance = abs(y12 - x12)) %>%
  mutate(resolution = convert_to_resolution(end.distance)) # utils_functions.R
df.circos.final.loop %>% head()

# Log transform the distance to highlight small distances
df.circos.input.final.loop <- df.circos.final.loop %>%
  mutate(
    distance_log = log10(distance + 0.01),
    midx = x12,
    midy = y12
  )

# Normalize height between 0.1 and 0.5
min_height <- 0.1
max_height <- 0.5
df.circos.input.log.final.loop <- df.circos.input.final.loop %>%
  mutate(height = min_height + (max_height - min_height) *
    (distance_log - min(distance_log)) /
    (max(distance_log) - min(distance_log)))


# Convert resolution to factor
df.circos.input.log.final.loop$distance <- as.factor(df.circos.input.log.final.loop$distance)

# Convert chromosome to factor
chromosome_order <- c(as.character(1:20), "X")

df.circos.input.log.final.loop <- df.circos.input.log.final.loop %>%
  mutate(
    chr1_clean = gsub("chr", "", chr1),
    chr1_clean = factor(chr1_clean, levels = chromosome_order)
  )

# Assign colors to each resolution level
colors <- brewer.pal(n = 3, name = "Set1")
resolution_levels <- levels(df.circos.input.log.final.loop$resolution)
# resolution_colors <- setNames(colors, resolution_levels)
resolution_colors <- c(
  # "5K" = "#377EB8",
  # "10K" = "lightcoral",
  # "25K" = "springgreen"
  # "5K" = "#a6cee3",
  # "10K" = "#1f78b4",
  # "25K" = "#1f3a93"
  "5K" = "#f8766d",
  "10K" = "#629bfe",
  "25K" = "#32ba36"
  # values = c("5K" = "#a6cee3", "10K" = "#1f78b4", "25K" = "#1f3a93")
)
# Output PDF file
pdf(file = "./figures/submission/lt2mb/circos_loops_by_resolution_either.6.lt2mb.ENSEMBL.mid.mid.pdf", height = 11 * 0.8, width = 8.5 * 0.8)

# Initialize circos with ideogram for rat genome
circos.initializeWithIdeogram(species = "rn7")

df.circos.input.log.final.loop %>% head()

# Draw global links
for (i in 1:nrow(df.circos.input.log.final.loop)) {
  circos.genomicLink(
    region1 = df.circos.input.log.final.loop[i, c("chr1", "midx", "midx")],
    region2 = df.circos.input.log.final.loop[i, c("chr2", "midy", "midy")],
    col = resolution_colors[as.character(df.circos.input.log.final.loop$resolution[i])],
    h = df.circos.input.log.final.loop$height[i],
    border = "black"
  )
}

# Add legend
legend("bottomright", legend = resolution_levels, fill = resolution_colors, title = "Resolution")
# Add title
# title("Circos Plot for Functional Loop", line = -1)
# Clear global plot
circos.clear()

# Plot per chromosome
unique_chromosomes <- levels(df.circos.input.log.final.loop$chr1_clean)
for (chr in unique_chromosomes) {
  plot_circos_for_chromosome(chr) # utils_functions.R
}

dev.off()

# 1. Load first two pages from PDF
pdf_path <- "./figures/submission/lt2mb/circos_loops_by_resolution_either.6.lt2mb.ENSEMBL.mid.mid.pdf"
img_list <- image_read_pdf(pdf_path, pages = 1:2, density = 300)

# Create labeled plots using ggdraw
plot_a <- ggdraw() +
  draw_image(img_list[[1]]) +
  draw_label("a", x = 0.02, y = 0.88, hjust = 0, vjust = 1, fontface = "bold", size = 16)

plot_b <- ggdraw() +
  draw_image(img_list[[2]]) +
  draw_label("b", x = 0.02, y = 0.88, hjust = 0, vjust = 1, fontface = "bold", size = 16)

# Combine into 1-row, 2-column layout
combined_plot <- plot_grid(plot_a, plot_b, nrow = 1)

# 5. Save as PNG
saving_plot_dual( # utils_functions.R
  output_dir = "./figures/submission/lt2mb",
  filename_base = "circos_first_two_panels_F8.ENSEMBL.mid.mid",
  plot_obj = combined_plot,
  scale_x = 1,
  scale_y = 1
)
