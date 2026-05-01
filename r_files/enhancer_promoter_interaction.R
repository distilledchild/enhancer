library("tidyverse")
library("GenomicRanges")
library("ggplot2")
library("cowplot")
library("circlize")
library("RColorBrewer")
library("magick")
library("RIdeogram")
library("rsvg")

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
source("./utils_functions.R") # Load all utility functions

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

# padding 1 distance && w/o capping only from DISTINCT loops: 31,417
OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping <- OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance %>%
  filter(!(x0 == 0 | chr.end.coord == y3)) # 31773 - 444 (282 + 162) = 31329 + 88 = 31,417
# filter(x0 == 0) # 282
# filter(chr.end.coord == y3)  # 162
OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping %>% dim() # 31,417
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

df.promoter.rn7.epd %>% dim() # 12,529
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
df.final.up.down.directional %>%
  distinct(loop.id) %>%
  nrow() # 31,019

df.final.up.down.directional.point <- df.final.up.down.directional %>%
  group_by(loop.id, WHERE) %>%
  filter(distance == min(distance)) %>% # multiple hits with the same distance in a loop and an anchor
  add_count(name = "tie.count.loop.where") %>%
  ungroup()

df.final.up.down.directional.point %>% dim() # 61,974 // 61,936// mid mid: 62,111 with_ties = TRUE// FYI: 62038 = 31019 * 2
df.final.up.down.directional.point %>% head(2)

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
# STAGE 1: How filtering muli-cases in an anchor
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

########################################
# Categorizing and filtering steps FOR MULTIPLE CASES AT AN ANCHOR IN A LOOP
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

# step for filtering only cases of PRO and TSS both in an anchor with the same gene_name (only removing TSS cases): multi케이스중 오직 pro tss 둘다 있는 케이스 처리: 두개 중 pro 우선권 (TSS제외)
df.final.up.down.directional.point.multi.in.loop.where.two.labelings <- df.final.up.down.directional.point.multi.in.loop.where %>% # head(2)
  filter(!is.na(gene_chr)) %>% # view() # 245 (28: is.na): 일단 gene이 NA인 케이스 제외
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

df.final.up.down.directional.point.multi.in.loop.where.two.labelings.selected <- df.final.up.down.directional.point.multi.in.loop.where.two.labelings %>% # 95/174 (273->174 after multiple case들 중 NA나 pro에게 우선권 준뒤...)
  filter(multi.2nd.filter == 1)
df.final.up.down.directional.point.multi.in.loop.where.two.labelings.selected %>% dim() # 95
df.final.up.down.directional.point.multi.in.loop.where.two.labelings.selected %>% head(2)

# 3. AFTER PROMOTER PRIORITY filteration, 2 rows from loop.id and WHERE
df.final.up.down.directional.point.multi.in.loop.where.two.labelings.ongoing <- df.final.up.down.directional.point.multi.in.loop.where.two.labelings %>%
  filter(multi.2nd.filter > 1) # 79 (이제 이 multi 케이스들 처리작업 해야함.)
df.final.up.down.directional.point.multi.in.loop.where.two.labelings.ongoing %>% dim() # 79
df.final.up.down.directional.point.multi.in.loop.where.two.labelings.ongoing %>% head(2)

df.final.up.down.directional.point.multi.in.loop.where.three.labelings <- df.final.up.down.directional.point.multi.in.loop.where.two.labelings.ongoing %>%
  arrange(loop.id, WHERE) %>%
  dplyr::select(multi.2nd.filter, gene_name, component, loop.id, everything()) %>%
  add_count(loop.id, WHERE, name = "multi.3rd.filter")

df.final.up.down.directional.point.multi.in.loop.where.three.labelings %>% dim()
df.final.up.down.directional.point.multi.in.loop.where.three.labelings %>% head(2)

# gene_name 필터링: 실제로 같은 gene인지 확인 (gene_id의 뒤에서 두번째 요소가 gene_name이므로, 앞에서 5번째 요소-ENSEMBL id 제거 후 비교)
df.final.up.down.directional.point.multi.in.loop.where.three.labelings.actual.same.gene.selected <- df.final.up.down.directional.point.multi.in.loop.where.three.labelings %>%
  filter(multi.3rd.filter > 1) %>%
  # print() # 79
  mutate(gene_id_trimmed = sapply(strsplit(gene_id, ":"), function(x) paste(x[-5], collapse = ":"))) %>% # dplyr::select(gene_id_trimmed) %>% view()
  group_by(loop.id, WHERE) %>%
  mutate(all_same_trimmed = n_distinct(gene_id_trimmed) == 1) %>%
  ungroup() %>%
  # count(all_same_trimmed) %>% view() # FALSE: 30, TRUE: 49
  filter(all_same_trimmed) %>% # # mutaul exclusive to L537 :::::::::: 49
  group_by(loop.id, WHERE) %>%
  filter(row_number() == 1) %>%
  ungroup()

df.final.up.down.directional.point.multi.in.loop.where.three.labelings.actual.same.gene.selected %>% dim() # 24/49
df.final.up.down.directional.point.multi.in.loop.where.three.labelings.actual.same.gene.selected %>% head(2)

# mannually filtering
df.final.up.down.directional.point.multi.in.loop.where.three.labelings.exon.number.and.mannually.selected <- df.final.up.down.directional.point.multi.in.loop.where.three.labelings %>%
  filter(multi.3rd.filter > 1) %>%
  mutate(gene_id_trimmed = sapply(strsplit(gene_id, ":"), function(x) paste(x[-5], collapse = ":"))) %>%
  group_by(loop.id, WHERE) %>%
  mutate(all_same_trimmed = n_distinct(gene_id_trimmed) == 1) %>%
  ungroup() %>%
  filter(!all_same_trimmed) %>% # mutaul exclusive to L522 :::::::::: 30
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

df.final.combined.using.labelings %>% head(2)
df.final.combined.using.labelings %>% dim() # 137
df.final.combined.using.labelings %>%
  count(loop.id, WHERE) %>%
  count(n)
# 1     1   131
# 2     2     3
df.final.up.down.directional.point.final <- bind_rows(df.final.up.down.directional.point.one.in.loop.where, df.final.combined.using.labelings) # 61,838 + 137 = 61,975

df.final.up.down.directional.point.final %>% dim() # 61,975
df.final.up.down.directional.point.final %>% head(3)
# 예상: 61,838 + 137 = 61,975

df.final.up.down.directional.point.final %>%
  distinct(loop.id) %>%
  nrow() # 31,019

df.final.up.down.directional.point.final %>%
  add_count(loop.id, WHERE, name = "final_count") %>%
  count(final_count)
# final_count count
# 1           1 61969
# 2           2     6 (L:549)

df.final.up.down.directional.point.final %>% dim() # 62,070 // 61,936// mid mid: 61,975/ FYI: 62038 = 31019 * 2

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
# STAGE 2: How getting the loops having only one PE interaction in a loop
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

#############################################################
# validity of the gene’s genomic position based on TSS or promoters selected by the nearestToDistance
#############################################################
df.final.up.down.directional.point.decision <- df.final.up.down.directional.point.final %>%
  separate(
    col = loop.id,
    into = c("chr1", "x1", "x2", "chr2", "y1", "y2", "resolution"),
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(
    UP = case_when(
      # (WHERE == "UP") & (gene_end < loop.mid) & (gene_start > as.numeric(x1)) ~ "OK",
      (WHERE == "UP") & (gene_end < as.numeric(y2)) & (gene_start > as.numeric(x1)) ~ "OK",
      WHERE == "UP" ~ "FAIL",
      TRUE ~ NA_character_
    ),
    DOWN = case_when(
      # (WHERE == "DOWN") & (gene_start > loop.mid) & (gene_end < as.numeric(y2)) ~ "OK",
      (WHERE == "DOWN") & (gene_start > as.numeric(x1)) & (gene_end < as.numeric(y2)) ~ "OK",
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

#   classification     n | latest.              | final
# 1 Both_FAIL      30497 | Both_FAIL      30498 | 1 Both_FAIL      26022
# 2 Both_OK         7606 | Both_OK         7606 | 2 Both_OK        10538
# 3 One_OK         23872 | One_OK         23871 | 3 One_OK         25415

###################################
# 1. ONE OK case in either anchor
###################################

df.final.up.down.directional.point.decision %>%
  distinct(loop.id) %>%
  count() # final: 31,019
df.final.up.down.directional.point.decision %>%
  filter(classification == "One_OK") %>% #
  # distinct(gene_id) %>% # 7,894
  distinct(loop.id) %>% # 11,935/31,019 # final: 12,707
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
# 1 7400.  23,371 66,691

# final
# 1 7400.  23,371 66,691

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
  # scale_y_log10() +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
df.final.up.down.directional.point.all.in.one.decision

df.final.up.down.directional.point.all.in.one.decision.histogram <- df.final.up.down.directional.point.decision %>%
  ggplot(aes(x = distance)) +
  geom_histogram(fill = "skyblue", color = "black", bins = 50) +
  geom_vline(xintercept = stats_decision$Q1, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = stats_decision$Median, linetype = "solid", color = "red") +
  geom_vline(xintercept = stats_decision$Q3, linetype = "dashed", color = "blue") +
  annotate("text", x = stats_decision$Q1, y = Inf, label = paste0("Q1=", round(stats_decision$Q1, 2)), vjust = 1.5, hjust = 1.1, color = "blue") +
  annotate("text", x = stats_decision$Median, y = Inf, label = paste0("Median=", round(stats_decision$Median, 2)), vjust = 3.0, hjust = 0.5, color = "red") +
  annotate("text", x = stats_decision$Q3, y = Inf, label = paste0("Q3=", round(stats_decision$Q3, 2)), vjust = 1.5, hjust = -0.1, color = "blue") +
  labs(
    title = "Histogram of Distance (All Data)",
    x = "Distance",
    y = "Count"
  ) +
  scale_x_log10() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
df.final.up.down.directional.point.all.in.one.decision.histogram

# Figure 8
df.final.up.down.directional.point.all.in.one.decision.200kb <- df.final.up.down.directional.point.decision %>%
  ggplot(aes(x = distance)) +
  # geom_density(fill = "skyblue") +
  geom_histogram(bins = 100, fill = "skyblue") +
  # annotate("text", y = 1, x = stats_decision$Q1, label = paste0("Q1=", round(stats_decision$Q1, 2)), vjust = -0.5, hjust = 1, color = "blue") +
  # annotate("text", y = 1, x = stats_decision$Median, label = paste0("Median=", round(stats_decision$Median, 2)), vjust = -0.5, hjust = 1, color = "red") +
  # annotate("text", y = 1, x = stats_decision$Q3, label = paste0("Q3=", round(stats_decision$Q3, 2)), vjust = -0.5, hjust = 1, color = "blue") +
  labs(
    # title = "Density plot of distance between TSS/Promoter and anchor",
    y = "Count",
    x = "Distance"
  ) +
  # scale_y_log10() +
  theme_minimal() +
  xlim(c(0, 5e5)) +
  theme(plot.title = element_text(hjust = 0.5))

df.final.up.down.directional.point.all.in.one.decision.200kb

saving_plot_dual( # utils_functions.R
  plot_obj = df.final.up.down.directional.point.all.in.one.decision.200kb,
  filename_base = "df.final.up.down.directional.point.all.in.one.decision.200kb",
  output_dir = "./figures/submission/lt2mb",
)

# boxplot by classification
# distribution.all.cases.distance <- df.final.up.down.directional.point.decision %>%
#   ggplot(aes(x = classification, y = distance, fill = classification)) +
#   geom_boxplot() +
#   theme_minimal() +
#   labs(
#     title = "Boxplot of Distance by Classification",
#     x = "Classification",
#     y = "Distance"
#   ) +
#   scale_y_log10() +
#   theme_minimal() +
#   theme(plot.title = element_text(hjust = 0.5))
# distribution.all.cases.distance

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
  count() # 11,935// final: 12,707
df.final.up.down.directional.point.decision.one.OK.filtered %>% head(2) # 11,937 // 11,935

# for just one OK
# approach_2nd_analyze_loops_by_threshold(
#   df.final.up.down.directional.point.decision.one.OK.filtered %>%
#     dplyr::rename(gene_id_id = gene_id) %>%
#     mutate(gene_id = str_split_n(str_split_n(component_id, ":", 6), "\\|", 1)),
#   threshold_distance = stats_decision$Q3,
#   top_n_genes = 70,
#   print_top_n = 50
# ) # utils_functions.R // final One_OK only gprofiler: https://biit.cs.ut.ee/gplink/l/aRdAXDm53Re

# distribution.OK.distance.q3 %>% dim()
# distribution.OK.distance.q3 %>% distinct(loop.id) %>% count()
# distribution.OK.distance.q3 %>% head(2)

###################################
# 2. TWO OK case in either anchor
###################################

# two cases in both anchors
df.final.up.down.directional.point.decision.both.OK <- df.final.up.down.directional.point.decision %>%
  filter(classification == "Both_OK") %>% # 3,803(7,606)/31,019
  # distinct(loop.id) %>%
  # count() %>%
  view()

# final: 5,269(10,538)

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
    min.dist = min(dist_diff, na.rm = TRUE),
    Q1 = quantile(dist_diff, 0.25, na.rm = TRUE),
    Q3 = quantile(dist_diff, 0.75, na.rm = TRUE),
    max.dist = max(dist_diff, na.rm = TRUE)
  )

stats_dist_diff

# Boxplot for distance difference
barplot.dist.diff.decision.both.OK <- ggplot(dist.diff.decision.both.OK, aes(x = "", y = dist_diff)) +
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
barplot.dist.diff.decision.both.OK

df.final.up.down.directional.point.decision.both.OK %>%
  count(distance)

##################################
# 큰 distance 제거하고 하나만 남기기
##################################
df.final.up.down.directional.point.decision.both.OK.filtered <- df.final.up.down.directional.point.decision.both.OK %>%
  group_by(loop.id) %>%
  mutate(strand = str_split_n(gene_id, ":", 4)) %>% # Extract strand (+ or -)
  # 1순위: distance(오름차순)
  # 2순위: Strand에 따른 전사 방향 (Tie-breaking)
  #        (+) Strand: 5' (Upstream) 이 좌측이므로 UP 앵커 선호
  #        (-) Strand: 5' (Upstream) 이 우측이므로 DOWN 앵커 선호
  arrange(distance, case_when(
    strand == "+" & WHERE == "UP" ~ 1,
    strand == "+" & WHERE == "DOWN" ~ 2,
    strand == "-" & WHERE == "DOWN" ~ 1,
    strand == "-" & WHERE == "UP" ~ 2,
    TRUE ~ 3
  )) %>%
  dplyr::slice(1) %>% # 최우선 순위 하나만 선택 (dplyr 충돌 방지)
  ungroup() %>%
  dplyr::select(-strand)

df.final.up.down.directional.point.decision.both.OK.filtered %>% dim() # 3,803: final: 5,269
df.final.up.down.directional.point.decision.both.OK.filtered %>% head(2)
# df.final.up.down.directional.point.decision.both.OK.filtered %>% filter(gene_name == "Trpa1")
# 동일 distance가 존재하는 loop.id 파악하기
# df.final.up.down.directional.point.decision.both.OK.tied <- df.final.up.down.directional.point.decision.both.OK %>%
#   group_by(loop.id) %>%
#   filter(n_distinct(distance) == 1, n() > 1) %>% # 같은 loop.id 내 distance가 모두 동일
#   ungroup()

# df.final.up.down.directional.point.decision.both.OK.tied # 0// final: 2
#   distance loop.id                                         chr1  x1           x2  chr2  y1      y2      resolution loop.mid loop_chr loop_start loop_end    gene_id                                     gene_name  component_id                                        component component_chr         component_start component_end gene_chr gene_start gene_end WHERE tie.count.loop.where multi.1st.filter has_pro_tss keep_row multi.2nd.filter  multi.3rd.filter gene_id_trimmed all_same_trimmed UP    DOWN  has_UP_OK  has_DOWN_OK classification
# 1    54998 chr5_4320000_4330000_chr5_4430000_4440000_10000 chr5  4320000 4330000  chr5  4430000 4440000 10000       4380000 chr5        4325000  4325000 chr5:4433021:4433570:+:ENSRNOG00000007354:Trpa1:27 Trpa1  chr5:4379999:4380001:+:ENSRNOG00000007354:Trpa1|tss tss       chr5                  4379999       4380001 chr5        4433021  4433570 UP            1               NA NA          NA                     NA                     NA NA              NA               OK    NA    TRUE                     TRUE        Both_OK
# 2    54998 chr5_4320000_4330000_chr5_4430000_4440000_10000 chr5  4320000 4330000  chr5  4430000 4440000 10000       4380000 chr5        4435000  4435000 chr5:4433021:4433570:+:ENSRNOG00000007354:Trpa1:27 Trpa1  chr5:4379999:4380001:+:ENSRNOG00000007354:Trpa1|tss tss       chr5                  4379999       4380001 chr5        4433021  4433570 DOWN          1               NA NA          NA                     NA                     NA NA              NA               NA    OK    TRUE                     TRUE        Both_OK

###################################
# 2. One OK + TWO OK case in either anchor
###################################

df.final.up.down.directional.point.decision %>% dim() # 61,975
df.final.up.down.directional.point.decision %>% head(2)
df.final.up.down.directional.point.decision %>% count(classification) #| final
# 1 Both_FAIL      30498|||||||||||            1 Both_FAIL      30498  1 Both_FAIL      30497| 1 Both_FAIL      26022
# 2 Both_OK         7606|||||||||||            2 Both_OK         7602  2 Both_OK         7606| 2 Both_OK        10538
# 3 One_OK         23871|||||||||||            3 One_OK         23874  3 One_OK         23872| 3 One_OK         25415

# CACHING: df.final.up.down.directional.point.decision.COMBINED.OK.filtered.lt.Q3
cache_file_combined_ok_filtered <- "../data/df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3_final_200kb.rds"

if (!file.exists(cache_file_combined_ok_filtered)) {
  message("Saving df.final.up.down.directional.point.decision.COMBINED.OK.filtered.lt.Q3_final to cache: ", cache_file_combined_ok_filtered)

  df.final.up.down.directional.point.decision.COMBINED.OK.filtered.lt.Q3 <- bind_rows(
    df.final.up.down.directional.point.decision.one.OK.filtered,
    df.final.up.down.directional.point.decision.both.OK.filtered
  ) %>% # view() # 15,738// final: 17,976
    # filter(distance < stats_decision$Q3) # %>% view() # 3,730
    # filter(distance < 150000) # %>% view() # 3,730
    filter(distance < 200000) # %>% view() # 3,730

  df.final.up.down.directional.point.decision.COMBINED.OK.filtered.lt.Q3 %>% dim() # 14,210// final: 16,113// 150kb: 17,417// 200kb: 17,648
  df.final.up.down.directional.point.decision.COMBINED.OK.filtered.lt.Q3 %>% head(2)
  df.final.up.down.directional.point.decision.COMBINED.OK.filtered.lt.Q3 %>% count(loop.id) # 14,210 PASS! (final: 16,113)

  saveRDS(df.final.up.down.directional.point.decision.COMBINED.OK.filtered.lt.Q3, cache_file_combined_ok_filtered)
} else {
  message("Cache file already exists: ", cache_file_combined_ok_filtered)
}

df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3 <- readRDS(cache_file_combined_ok_filtered)
df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3 %>% dim() # 14210// final: 16,113// 150kb: 17,417// 200kb: 17,648
df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3 %>% head(2)

approach_2nd_analyze_loops_by_threshold(
  df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3 %>%
    dplyr::rename(gene_id_id = gene_id) %>%
    mutate(gene_id = str_split_n(str_split_n(component_id, ":", 6), "\\|", 1)),
  # threshold_distance = stats_decision$Q3,
  threshold_distance = 200000,
  top_n_genes = 30,
  print_top_n = 80
) # utils_functions.R// final gprofiler top 30: https://biit.cs.ut.ee/gplink/l/aHWE_iUdCTy

################################################################################################
################################################################################################
# AlphaGenome input generation (STEPWISE, UCSC API only)
################################################################################################
################################################################################################
# source("./alphagenome_export.R") # Not needed per user request

df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3 %>%
  count(component) #  |     # 200kb
# 1 pro        4776  1 pro        4569|||||||| 1 pro        4569||||| 1 pro        5429| 1 pro        5729
# 2 tss       10891  2 tss        9641|||||||| 2 tss        9641||||| 2 tss       10684| 2 tss       11919

######################################################################
######################################################################

# final
final.loops.from.tss.step <- df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3 %>%
  filter(component == "tss")
final.loops.from.tss.step %>% dim() # 9,641// 10,684

final.loops.from.promoter.step <- df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3 %>%
  filter(component == "pro")
final.loops.from.promoter.step %>% dim() # 4,569// 5,429

df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3 %>% head(3)
df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3 %>% dim() # 23524// 14210// final: 16,113// 200kb: 17,648
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
component_dist # 200kb
# 1 pro        4035  // 1 pro        4569| 1 pro        5729
# 2 tss        5142  // 2 tss        9641| 2 tss       11919

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
final.loops.from.promoter.step$loop.id # 4569// final: 5,429// 200kb: 5,729
final.loops.from.tss.step$loop.id # 9641// final: 10,684// 200kb: 11,919

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
final.loops.from.tss.step %>% dim() # mid mid: 9,641// final:10,684// 200kb: 11,919
final.loops.from.promoter.step %>% dim() # mid mid: 4,569// final:5,429// 200kb: 5,729

overlapping_loops <- extract_overlapping_loops( # utils_functions.R
  ctcf_data = final.loops.from.ctcf.step,
  promoter_data = final.loops.from.promoter.step,
  tss_data = final.loops.from.tss.step
)

overlapping_loops$ctcf_promoter_only # PASS// 6648 // final: 3,973// final: 4712// 200kb: 4,960
overlapping_loops$ctcf_tss_only # PASS// 8613 // final: 8261// final: 9111// 200kb: 10,125
overlapping_loops$ctcf_promoter_tss_overlap # PASS// 3423// final: 0 // 200kb: 0

df.ctcf.promoter.only.loop <- data.frame(loop.id = overlapping_loops$ctcf_promoter_only, category = rep("CP", length(overlapping_loops$ctcf_promoter_only)), stringsAsFactors = FALSE)
df.ctcf.tss.only.loop <- data.frame(loop.id = overlapping_loops$ctcf_tss_only, category = rep("CT", length(overlapping_loops$ctcf_tss_only)), stringsAsFactors = FALSE)
df.ctcf.promoter.tss.overlap.loop <- data.frame(loop.id = overlapping_loops$ctcf_promoter_tss_overlap, category = rep("CPT", length(overlapping_loops$ctcf_promoter_tss_overlap)), stringsAsFactors = FALSE)

df.ctcf.promoter.only.loop %>% dim() # 5497//3973// final : 4712// 200kb: 4,960
df.ctcf.tss.only.loop %>% dim() # 7799// 8261// final: 9111// 200kb: 10,125
df.ctcf.promoter.tss.overlap.loop %>% dim() # 2302//0// 200kb: 0

# sub.4 dedups: 11891 = 679 + 3096 + 8116
df.final.loop <- bind_rows(df.ctcf.promoter.only.loop, df.ctcf.tss.only.loop)

df.final.loop # sub.4 any, either 6// 12234/31019(lt2mb) (0.3944034)// final: 13823/31019(0.4456301)// 200kb: 15,085/30,109 (0.501013)
df.final.loop %>% head()
df.final.loop %>% dim() # 15598: 5656 + 7980 + 1962// final(loop-mid): 12,234// final: 13823// 200kb: 15,085
df.final.loop %>% count(category) # CP 3973 CT 8261// final CP 4712 CT 9111// 200kb: 4,960CT 10,125

df.final.loop %>%
  mutate(resolution = str_split_n(loop.id, "_", 7)) %>%
  count(resolution) # utils_functions.R
#  resolution    n |||          200kb
# 1      10000 4369|||   | 1      10000 5430
# 2      25000 6104|||   | 2      25000 7416
# 3       5000 1761|||   | 3       5000 2239

save(df.final.loop, file = "./figures/submission/lt2mb/df_final_loop_sub.4.any.lt2mb.ENSEMBL.mid.mid.final.filter.200kb.rda")
write.csv(df.final.loop, file = "./figures/submission/lt2mb/df_final_loop_sub.4.any.lt2mb.ENSEMBL.mid.mid.final.filter.200kb.csv", row.names = FALSE)

df.final.loop %>% dim() # 23640// 18681// 18684// 12234// final: 13823// 200kb: 15,085
df.final.loop %>% head(2)
df.DISTINCT.loop.deep.sample.all %>% head()

df.final.DISTINCT.loop.joined <- df.final.loop %>%
  left_join(df.DISTINCT.loop.deep.sample.all, by = "loop.id")
df.final.DISTINCT.loop.joined %>% dim()
df.final.DISTINCT.loop.joined %>% head(2)

df.final.loop %>% distinct(loop.id) # 15085

##########################################################
# Chromosome-level gene/CTCF density correlations
# Supplementary six-panel figure
##########################################################

########################
# Loading df.chromosome.data (originally in enhancer_promoter_interaction_figures.R)
########################
cache_file_chromosome_data <- "../data/df.chromosome.data.rds"
cache_file_chromosome_data

if (file.exists(cache_file_chromosome_data)) {
  message("Loading cached chromosome data from: ", cache_file_chromosome_data)
  df.chromosome.data <- readRDS(cache_file_chromosome_data)
} else {
  message("Processing and caching chromosome data...")
  df.chromosome.data <- read.table(file = "../data/rn7_chromosome_length_from_ucsc.tsv", sep = "\t") %>%
    mutate(start = 0) %>%
    dplyr::rename(chr = V1, end = V2)

  saveRDS(df.chromosome.data, cache_file_chromosome_data)
}

df.chromosome.data %>% head(3)
df.chromosome.data <- df.chromosome.data %>% mutate(CE_start = NA, CE_end = NA)
df.chromosome.data

########################
# Loading df.DISTINCT.fimo.2nd.trial.ctcf (originally in enhancer_promoter_interaction_figures.R)
########################
cache_file_distinct_fimo_2nd_ctcf <- "../data/df.DISTINCT.fimo.2nd.trial.ctcf.rds"

if (file.exists(cache_file_distinct_fimo_2nd_ctcf)) {
  message("Loading cached DISTINCT fimo 2nd trial ctcf data from: ", cache_file_distinct_fimo_2nd_ctcf)
  df.DISTINCT.fimo.2nd.trial.ctcf <- readRDS(cache_file_distinct_fimo_2nd_ctcf)
} else {
  message("Processing and caching DISTINCT fimo 2nd trial ctcf data...")
  df.init.ctcf <- read.table(file = "~/dropbox/Gateway_to_Hao/enhancer/data/ctcf/submission/E4/fimo_E4_submission_trial.txt", header = TRUE, sep = "\t") %>%
    dplyr::rename(chr = sequence_name, end = stop) %>%
    mutate(length = end - start)

  df.DISTINCT.fimo.2nd.trial.ctcf <- df.init.ctcf %>%
    distinct(chr, start, end) %>% # .4:3191859 ************** NO STRAND INFO
    mutate(start = as.numeric(start)) %>%
    mutate(end = as.numeric(end)) %>%
    mutate(ctcf_pos = as.numeric(round((start + end) / 2))) %>%
    mutate(id = str_c(chr, "_", start, "_", end, "_", ctcf_pos))

  saveRDS(df.DISTINCT.fimo.2nd.trial.ctcf, cache_file_distinct_fimo_2nd_ctcf)
}

# CTCF positions for chromosome-level density comparison
df.DISTINCT.fimo.2nd.trial.ctcf %>% head()
df.DISTINCT.fimo.2nd.trial.ctcf # 3191859 (sub .4)

df_ctcf_ideogram <- df.DISTINCT.fimo.2nd.trial.ctcf %>%
  dplyr::select(chr, start, end)

df_ctcf_ideogram

# Gene positions for chromosome-level density comparison
output_dir <- "./figures/submission/lt2mb"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

chromosome_levels <- c(as.character(1:20), "X", "Y")
valid_chromosomes <- paste0("chr", chromosome_levels)
make_file_tag <- function(...) {
  str_c(...) %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "")
}
##############################################
# functions for reading gtf files (Ensembl)
##############################################
read_ensembl_gene_catalog <- function(gtf_file) {
  message("Reading Ensembl gene catalog from: ", gtf_file)

  ensembl_gene_raw <- read_tsv(
    gtf_file,
    comment = "#",
    col_names = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"),
    col_types = cols(.default = "c")
  ) %>%
    filter(feature == "gene") %>%
    filter(chr %in% c(as.character(1:20), "X", "Y")) %>%
    mutate(
      chr = str_c("chr", chr),
      start = as.numeric(start),
      end = as.numeric(end)
    )

  ensembl_gene_keys <- get_attribute_keys(ensembl_gene_raw$attribute)

  ensembl_gene_raw %>%
    bind_cols(ensembl_gene_raw$attribute %>% map_dfr(~ extracting_attributes(.x, keys = ensembl_gene_keys))) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    dplyr::select(chr, start, end, strand, gene_id, gene_name, gene_biotype)
}

##############################################
# functions for reading gtf files (NCBI)
##############################################
# NC_* accession to chr name mapping for mRatBN7.2 (GCF_015227675.2)
# NC_051336.1=chr1 .. NC_051355.1=chr20, NC_051356.1=chrX, NC_051357.1=chrY, NC_001665.2=chrM
ncbi_accession_to_chr <- c(
  "NC_051336.1" = "chr1",  "NC_051337.1" = "chr2",  "NC_051338.1" = "chr3",
  "NC_051339.1" = "chr4",  "NC_051340.1" = "chr5",  "NC_051341.1" = "chr6",
  "NC_051342.1" = "chr7",  "NC_051343.1" = "chr8",  "NC_051344.1" = "chr9",
  "NC_051345.1" = "chr10", "NC_051346.1" = "chr11", "NC_051347.1" = "chr12",
  "NC_051348.1" = "chr13", "NC_051349.1" = "chr14", "NC_051350.1" = "chr15",
  "NC_051351.1" = "chr16", "NC_051352.1" = "chr17", "NC_051353.1" = "chr18",
  "NC_051354.1" = "chr19", "NC_051355.1" = "chr20",
  "NC_051356.1" = "chrX",  "NC_051357.1" = "chrY",  "NC_001665.2" = "chrM"
)

read_ncbi_gff3_gene_catalog <- function(gff3_file) {
  message("Reading NCBI RefSeq gene catalog from GFF3: ", gff3_file)

  gff3_raw <- read_tsv(
    gff3_file,
    comment = "#",
    col_names = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"),
    col_types = cols(.default = "c")
  ) %>%
    filter(feature == "gene") %>%
    filter(chr %in% names(ncbi_accession_to_chr)) %>%
    mutate(
      chr = ncbi_accession_to_chr[chr],
      start = as.numeric(start),
      end = as.numeric(end)
    ) %>%
    filter(chr %in% valid_chromosomes)

  # GFF3 attribute: key=value;key=value format
  gff3_raw %>%
    mutate(
      gene_id = str_match(attribute, "ID=gene-([^;]+)")[, 2],
      gene_name = str_match(attribute, "gene=([^;]+)")[, 2],
      gene_biotype = str_match(attribute, "gene_biotype=([^;]+)")[, 2]
    ) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    dplyr::select(chr, start, end, strand, gene_id, gene_name, gene_biotype)
}

# function to build density summary (6 CSV files, A, B, C(ENSEMBL) and D, E, F(NCBI))
build_density_summary <- function(chrom_data, df_gene_input, ctcf_data) {
  chrom_data %>%
    mutate(
      chr = str_remove(as.character(chr), "^chr"),
      chromosome_length_mb = end / 1e6 # Unit for mb
    ) %>%
    group_by(chr) %>%
    summarise(chromosome_length_mb = max(chromosome_length_mb, na.rm = TRUE), .groups = "drop") %>%
    left_join(
      df_gene_input %>%
        mutate(chr = str_remove(as.character(chr), "^chr")) %>%
        count(chr, name = "gene_count"),
      by = "chr"
    ) %>%
    left_join(
      ctcf_data %>%
        mutate(chr = str_remove(as.character(chr), "^chr")) %>%
        count(chr, name = "ctcf_count"),
      by = "chr"
    ) %>%
    mutate(
      gene_count = replace_na(gene_count, 0L),
      ctcf_count = replace_na(ctcf_count, 0L),
      genes_per_mb = gene_count / chromosome_length_mb,
      ctcf_per_mb = ctcf_count / chromosome_length_mb
    ) %>%
    arrange(factor(chr, levels = chromosome_levels))
}

run_density_analysis <- function(df_gene_input, panel_label, source_label, gene_category, chrom_data, ctcf_data, out_dir) {
  message(paste0(">>> Running panel ", panel_label, ": ", source_label, " / ", gene_category))

  density_summary <- build_density_summary(chrom_data, df_gene_input, ctcf_data)
  cor_pearson <- cor.test(density_summary$genes_per_mb, density_summary$ctcf_per_mb, method = "pearson")

  stats_df <- tibble(
    panel = panel_label,
    annotation_source = source_label,
    gene_category = gene_category,
    method = "pearson",
    estimate = unname(cor_pearson$estimate),
    p_value = cor_pearson$p.value,
    gene_count_total = sum(density_summary$gene_count),
    ctcf_count_total = sum(density_summary$ctcf_count)
  )

  density_summary <- density_summary %>%
    mutate(
      panel = panel_label,
      annotation_source = source_label,
      gene_category = gene_category
    ) %>%
    relocate(panel, annotation_source, gene_category)

  file_tag <- make_file_tag(panel_label, source_label, gene_category)
  write.csv(density_summary, file = file.path(out_dir, paste0("chromosome_gene_ctcf_density_summary_", file_tag, ".csv")), row.names = FALSE)

  plot_obj <- density_summary %>%
    ggplot(aes(x = genes_per_mb, y = ctcf_per_mb)) +
    geom_point(size = 2.2, color = "#2F5597") +
    geom_smooth(method = "lm", se = FALSE, color = "#C44E52", linewidth = 0.7) +
    geom_text(aes(label = chr), nudge_y = max(density_summary$ctcf_per_mb, na.rm = TRUE) * 0.03, size = 2.5, check_overlap = TRUE) +
    labs(
      title = paste0(source_label, ": ", gene_category),
      subtitle = paste0("Pearson r = ", round(cor_pearson$estimate, 3), " (p = ", formatC(cor_pearson$p.value, format = "e", digits = 2), ")"),
      x = "Genes per Mb",
      y = "CTCF sites per Mb"
    ) +
    theme_bw(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 8),
      panel.grid.minor = element_blank()
    )

  return(list(summary = density_summary, stats = stats_df, plot = plot_obj))
}

# data location
ensembl_gtf_file <- "../data/Rattus_norvegicus.mRatBN7.2.113.gtf"
ncbi_refseq_gff3_file <- "../data/GCF_015227675.2_mRatBN7.2_genomic.gff" # https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/227/675/GCF_015227675.2_mRatBN7.2/
# ../data/ncbiRefSeq.gtf is from https://hgdownload.soe.ucsc.edu/goldenPath/rn7/bigZips/genes/

################################################
# Run density analysis for each category (panels A–D)
################################################
df_ensembl_gene_catalog <- read_ensembl_gene_catalog(ensembl_gtf_file)
df_ncbi_refseq_gene_catalog <- read_ncbi_gff3_gene_catalog(ncbi_refseq_gff3_file)

df_ensembl_gene_catalog %>% head()
df_ncbi_refseq_gene_catalog %>% head()

df_ensembl_gene_catalog %>%
  count(gene_biotype, name = "count") %>%
  arrange(desc(count))
#    gene_biotype         count
#    <chr>                <int>
#  1 protein_coding       23038
#  2 lncRNA                2465

df_ncbi_refseq_gene_catalog %>%
  count(gene_biotype, name = "count") %>%
  arrange(desc(count))
#    gene_biotype   count
#  1 protein_coding 21944
#  2 lncRNA          7815

gene_sets <- list(
  ensembl_total = df_ensembl_gene_catalog,
  ensembl_protein_coding = df_ensembl_gene_catalog %>% filter(gene_biotype == "protein_coding"),
  ensembl_lncRNA = df_ensembl_gene_catalog %>% filter(gene_biotype == "lncRNA"),
  ncbi_refseq_total = df_ncbi_refseq_gene_catalog,
  ncbi_refseq_protein_coding = df_ncbi_refseq_gene_catalog %>% filter(gene_biotype == "protein_coding"),
  ncbi_refseq_lncRNA = df_ncbi_refseq_gene_catalog %>% filter(gene_biotype == "lncRNA")
)

analysis_design <- tribble(
  ~panel, ~annotation_source, ~gene_category, ~gene_set_key,
  "A", "Ensembl", "All genes", "ensembl_total",
  "B", "Ensembl", "Protein-coding genes", "ensembl_protein_coding",
  "C", "Ensembl", "lncRNA genes", "ensembl_lncRNA",
  "D", "NCBI RefSeq", "All genes", "ncbi_refseq_total",
  "E", "NCBI RefSeq", "Protein-coding genes", "ncbi_refseq_protein_coding",
  "F", "NCBI RefSeq", "lncRNA genes", "ncbi_refseq_lncRNA"
)

# generating 6 csv files
density_results_list <- pmap(
  analysis_design,
  function(panel, annotation_source, gene_category, gene_set_key) {
    run_density_analysis(
      df_gene_input = gene_sets[[gene_set_key]],
      panel_label = panel,
      source_label = annotation_source,
      gene_category = gene_category,
      chrom_data = df.chromosome.data,
      ctcf_data = df_ctcf_ideogram,
      out_dir = output_dir
    )
  }
)
names(density_results_list) <- analysis_design$panel

comparison_stats <- bind_rows(map(density_results_list, "stats"))
comparison_summary <- bind_rows(map(density_results_list, "summary"))

# 6 CSV files merged into one file: chromosome_gene_ctcf_density_summary_six_panel.csv, AND correlation: chromosome_gene_ctcf_density_correlation_six_panel.csv
write.csv(comparison_stats, file = file.path(output_dir, "chromosome_gene_ctcf_density_correlation_six_panel.csv"), row.names = FALSE)
write.csv(comparison_summary, file = file.path(output_dir, "chromosome_gene_ctcf_density_summary_six_panel.csv"), row.names = FALSE)

six_panel_gene_ctcf_density <- cowplot::plot_grid(
  plotlist = map(density_results_list, "plot"),
  labels = analysis_design$panel,
  ncol = 3,
  align = "hv",
  label_size = 14,
  label_fontface = "bold"
)

# saving 6 panel plot (ENSEMBL vs NCBI RefSeq genes with CTCF)
saving_plot_dual(
  plot_obj = six_panel_gene_ctcf_density,
  filename_base = "chromosome_gene_ctcf_density_six_panel_ensembl_ncbi_refseq",
  output_dir = output_dir,
  width_in = 11,
  height_in = 8.5,
  scale_x = 0.8,
  scale_y = 0.8,
  dpi = 300
)

message("\n>>> Six-panel Pearson correlation statistics:")
print(comparison_stats)

##########################################################
# Figure 4: CTCF density ideogram + NCBI RefSeq density
##########################################################

# ── Panel a: CTCF binding site density ideogram ──────────────────────────────
# Prepare karyotype data for RIdeogram (Chr, Start, End format)
karyotype_data <- df.chromosome.data %>%
  mutate(
    Chr = str_remove(chr, "^chr"),
    Start = 0,
    End = as.numeric(end)
  ) %>%
  filter(Chr %in% chromosome_levels) %>%
  dplyr::select(Chr, Start, End) %>%
  arrange(factor(Chr, levels = chromosome_levels))

# CTCF density bins for ideogram heatmap overlay
bin_size_ideogram <- 1000000
ctcf_density_for_ideogram <- process_feature_bins(
  feature_df = df_ctcf_ideogram %>% mutate(chr = str_remove(chr, "^chr")),
  chromosome_ends = karyotype_data %>% dplyr::rename(chr = Chr, end = End),
  bin_size = bin_size_ideogram,
  label = "CTCF"
) %>%
  dplyr::rename(Chr = chr, Start = start, End = end, Value = Value) %>%
  dplyr::select(Chr, Start, End, Value)

# Generate ideogram SVG via RIdeogram
ideogram(
  karyotype = karyotype_data,
  overlaid = ctcf_density_for_ideogram
)

# Convert SVG to PNG
ideogram_svg <- "chromosome.svg"
ideogram_png <- file.path(output_dir, "figure4_panel_a_ctcf_ideogram.png")
rsvg::rsvg_png(ideogram_svg, file = ideogram_png, width = 2400)

# White background and trim excess whitespace so panel (a) fills its slot.
image_read(ideogram_png) %>%
  image_background(color = "white") %>%
  image_trim(fuzz = 5) %>%
  image_border(color = "white", geometry = "80x80") %>%
  image_write(ideogram_png)

# Clean up SVG
if (file.exists(ideogram_svg)) file.remove(ideogram_svg)

message(">>> Figure 4 panel (a) ideogram saved: ", ideogram_png)

# ── Panels b, c, d: NCBI RefSeq density scatter plots ──────────────────────
# Extract the three NCBI RefSeq plots (D, E, F from the six-panel)
format_figure4_correlation_plot <- function(plot_obj, plot_title) {
  plot_obj +
    labs(title = plot_title) +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 11),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9),
      plot.margin = margin(7, 6, 4, 6)
    )
}

plot_ncbi_all <- format_figure4_correlation_plot(density_results_list[["D"]]$plot, "b. All genes")
plot_ncbi_protein <- format_figure4_correlation_plot(density_results_list[["E"]]$plot, "c. Protein-coding genes")
plot_ncbi_lncRNA <- format_figure4_correlation_plot(density_results_list[["F"]]$plot, "d. lncRNA genes")

# ── Combine all panels into Figure 4 ────────────────────────────────────────
# Panel a: ideogram (read as ggdraw image, no label here — labels added by plot_grid)
panel_a <- ggdraw() +
  draw_image(ideogram_png, scale = 0.9) +
  draw_label("a", x = 0.08, y = 0.985, hjust = 0, vjust = 1, fontface = "bold", size = 13)

# ── Figure 4: 1 row x 4 columns layout ─────────────────────────────────────
figure_4_1x4_grid <- cowplot::plot_grid(
  panel_a, plot_ncbi_all, plot_ncbi_protein, plot_ncbi_lncRNA,
  ncol = 4,
  align = "hv",
  rel_widths = c(1.35, 1, 1, 1)
)

figure_4_1x4 <- ggdraw() +
  draw_plot(figure_4_1x4_grid, x = -0.05, y = -0.055, width = 1.043, height = 1.11)

saving_plot_dual(
  plot_obj = figure_4_1x4,
  filename_base = "figure_4_ctcf_ideogram_ncbi_refseq_density_1x4",
  output_dir = output_dir,
  width_in = 15.4,
  height_in = 4.3,
  scale_x = 0.9,
  scale_y = 0.9,
  dpi = 300
)

message(">>> Figure 4 (1x4) saved to: ", output_dir)

##########################################################
# additional for discussion
##########################################################
# Check correlation between chromosome length and CTCF density
message("\n>>> Checking correlation between Chromosome Length and CTCF Density:")
df_len_ctcf <- density_results_list[["A"]]$summary %>%
  dplyr::select(chr, chromosome_length_mb, ctcf_per_mb)

cor_len_ctcf_pearson <- cor.test(df_len_ctcf$chromosome_length_mb, df_len_ctcf$ctcf_per_mb, method = "pearson")
cor_len_ctcf_spearman <- cor.test(df_len_ctcf$chromosome_length_mb, df_len_ctcf$ctcf_per_mb, method = "spearman")

len_ctcf_correlation_stats <- tibble(
  comparison = "Chromosome length vs CTCF density",
  method = c("pearson", "spearman"),
  estimate = c(
    unname(cor_len_ctcf_pearson$estimate),
    unname(cor_len_ctcf_spearman$estimate)
  ),
  p_value = c(
    cor_len_ctcf_pearson$p.value,
    cor_len_ctcf_spearman$p.value
  )
)

message(paste0(
  "Pearson correlation (Length vs CTCF Density): r = ", round(cor_len_ctcf_pearson$estimate, 3),
  ", p-value = ", formatC(cor_len_ctcf_pearson$p.value, format = "e", digits = 2)
))
message(paste0(
  "Spearman correlation (Length vs CTCF Density): rho = ", round(cor_len_ctcf_spearman$estimate, 3),
  ", p-value = ", formatC(cor_len_ctcf_spearman$p.value, format = "e", digits = 2)
))

# Save chromosome-level summary and correlation statistics for record
write.csv(df_len_ctcf, file = file.path(output_dir, "chromosome_length_vs_ctcf_density_summary.csv"), row.names = FALSE)
write.csv(len_ctcf_correlation_stats, file = file.path(output_dir, "chromosome_length_vs_ctcf_density_correlation.csv"), row.names = FALSE)

plot_len_ctcf <- df_len_ctcf %>%
  ggplot(aes(x = chromosome_length_mb, y = ctcf_per_mb)) +
  geom_point(size = 2.6, color = "#2F5597") +
  geom_smooth(method = "lm", se = FALSE, color = "#C44E52", linewidth = 0.8) +
  geom_text(aes(label = chr), nudge_y = max(df_len_ctcf$ctcf_per_mb, na.rm = TRUE) * 0.03, size = 3, check_overlap = TRUE) +
  labs(
    title = "Chromosome length vs CTCF density",
    subtitle = paste0(
      "Pearson r = ", round(cor_len_ctcf_pearson$estimate, 3),
      " (p = ", formatC(cor_len_ctcf_pearson$p.value, format = "e", digits = 2), "); ",
      "Spearman rho = ", round(cor_len_ctcf_spearman$estimate, 3),
      " (p = ", formatC(cor_len_ctcf_spearman$p.value, format = "e", digits = 2), ")"
    ),
    x = "Chromosome length (Mb)",
    y = "CTCF sites per Mb"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

saving_plot_dual(
  plot_obj = plot_len_ctcf,
  filename_base = "chromosome_length_vs_ctcf_density_correlation",
  output_dir = output_dir,
  width_in = 5.5,
  height_in = 4.5,
  scale_x = 1,
  scale_y = 1,
  dpi = 300
)

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

# Convert chromosome to factor (include chrY)
chromosome_order <- c(as.character(1:20), "X", "Y")

df.circos.input.log.final.loop <- df.circos.input.log.final.loop %>%
  mutate(
    chr1_clean = gsub("chr", "", chr1),
    chr1_clean = factor(chr1_clean, levels = chromosome_order),
    chr2_clean = gsub("chr", "", chr2),
    chr2_clean = factor(chr2_clean, levels = chromosome_order)
  )

# Assign colors to each resolution level
colors <- brewer.pal(n = 3, name = "Set1")
resolution_levels <- levels(df.circos.input.log.final.loop$resolution)
# resolution_colors <- setNames(colors, resolution_levels)
resolution_colors <- c(
  "5K" = "#f8766d",
  "10K" = "#629bfe",
  "25K" = "#32ba36"
)
# Output PDF file
pdf(
  file = "./figures/submission/lt2mb/circos_loops_by_resolution_either.6.lt2mb.ENSEMBL.mid.mid.final.200kb.pdf",
  height = 5.2, width = 5.2
)

# Plot per chromosome in fixed order: chr1..chr20, chrX, chrY
present_chromosomes <- union(
  as.character(na.omit(unique(df.circos.input.log.final.loop$chr1_clean))),
  as.character(na.omit(unique(df.circos.input.log.final.loop$chr2_clean)))
)
unique_chromosomes <- chromosome_order[chromosome_order %in% present_chromosomes]
for (chr in unique_chromosomes) {
  plot_circos_for_chromosome(
    chr = chr,
    page_label = paste("Chromosome", chr),
    show_legend = identical(chr, "1")
  ) # utils_functions.R
}

dev.off()

# 1. Load first two pages from PDF
# Create labeled plots using ggdraw
# Build Figure 8 directly so panel a is whole-genome and panel b is chr1.
panel_a_path <- "./figures/submission/lt2mb/circos_all_chromosomes_F8_panel_a.ENSEMBL.mid.mid.final.200kb.png"
panel_b_path <- "./figures/submission/lt2mb/circos_chr1_F8_panel_b.ENSEMBL.mid.mid.final.200kb.png"

png(filename = panel_a_path, width = 2400, height = 1800, res = 300, bg = "white")
par(mar = c(1.2, 1.2, 1.2, 1.2), xpd = NA)
plot_circos_all_chromosomes(
  page_label = NULL,
  show_legend = FALSE
) # utils_functions.R
dev.off()

png(filename = panel_b_path, width = 2400, height = 1800, res = 300, bg = "white")
par(mar = c(1.2, 1.2, 1.2, 1.2))
plot_circos_for_chromosome(
  chr = "1",
  page_label = NULL,
  show_legend = FALSE
) # utils_functions.R
dev.off()

plot_a <- ggdraw() +
  draw_image(panel_a_path, scale = 1.1) +
  draw_label("a", x = 0.02, y = 0.88, hjust = 0, vjust = 1, fontface = "bold", size = 16)

plot_b <- ggdraw() +
  draw_image(panel_b_path, scale = 0.9345) +
  draw_label("b", x = 0.02, y = 0.88, hjust = 0, vjust = 1, fontface = "bold", size = 16)

legend_df <- tibble(
  resolution = factor(names(resolution_colors), levels = names(resolution_colors)),
  x = 1,
  y = 1
)

legend_plot <- ggplot(legend_df, aes(x = x, y = y, fill = resolution)) +
  geom_point(shape = 22, size = 2.7, stroke = 0.24) +
  scale_fill_manual(values = resolution_colors) +
  guides(
    fill = guide_legend(
      title = "Resolution",
      title.position = "top",
      ncol = 1,
      byrow = TRUE
    )
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.title = element_text(size = 7.2, face = "plain"),
    legend.text = element_text(size = 6.6),
    legend.key.height = grid::unit(0.108, "in"),
    legend.key.width = grid::unit(0.108, "in"),
    legend.spacing.x = grid::unit(0.048, "in"),
    legend.spacing.y = grid::unit(0.012, "in"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  )

legend_grob <- get_legend(legend_plot)

base_panels <- plot_grid(plot_a, plot_b, nrow = 1, rel_widths = c(1, 1))

combined_plot <- ggdraw() +
  draw_plot(base_panels, x = 0, y = 0.12, width = 1, height = 0.88) +
  draw_grob(legend_grob, x = 0.468, y = 0.228, width = 0.06, height = 0.132) +
  draw_label("All chromosomes", x = 0.24, y = 0.25, fontface = "bold", size = 10) +
  draw_label("Chromosome 1", x = 0.76, y = 0.25, fontface = "bold", size = 10)

# 5. Save as PNG for FIGURE 8
saving_plot_dual( # utils_functions.R
  output_dir = "./figures/submission/lt2mb",
  filename_base = "circos_first_two_panels_F8.ENSEMBL.mid.mid.final.200kb",
  plot_obj = combined_plot
)

f8_base <- "./figures/submission/lt2mb/circos_first_two_panels_F8.ENSEMBL.mid.mid.final.200kb"
f8_png_path <- paste0(f8_base, ".png")
f8_pdf_path <- paste0(f8_base, ".pdf")
f8_img <- image_read(f8_png_path)
f8_data <- image_data(f8_img, channels = "rgb")
img_width <- dim(f8_data)[2]
img_height <- dim(f8_data)[3]
non_white_rows <- apply(
  f8_data,
  3,
  function(slice) any(slice != as.raw(255))
)

if (any(non_white_rows)) {
  max_content_row <- max(which(non_white_rows))
  bottom_whitespace <- img_height - max_content_row
  crop_pixels <- floor(bottom_whitespace / 2)

  if (crop_pixels > 0) {
    cropped_height <- img_height - crop_pixels
    f8_img <- image_crop(
      f8_img,
      geometry = paste0(img_width, "x", cropped_height, "+0+0")
    )
    image_write(f8_img, path = f8_png_path, format = "png")
    image_write(f8_img, path = f8_pdf_path, format = "pdf")
  }
}

# 입력/출력 경로
infile <- "~/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/lt2mb/circos_loops_by_resolution_either.6.lt2mb.ENSEMBL.mid.mid.final.200kb.pdf"
outfile <- "~/dropbox/Gateway_to_Hao/enhancer/r_files/figures/submission/lt2mb/circos_loops_by_resolution_either.6.lt2mb.ENSEMBL.mid.mid.final.200kb.png"

# 1) PDF → 이미지 리스트로 읽기
imgs <- image_read_pdf(infile, density = 600) # density 높이면 더 선명

# Preserve the full page legend area when rasterizing the PDF pages for collage.
page_info <- image_info(imgs)
padded_width <- ceiling(max(page_info$width) * 1.12)
padded_height <- max(page_info$height)
imgs <- image_extent(
  imgs,
  geometry = paste0(padded_width, "x", padded_height),
  gravity = "west",
  color = "white"
)

# 2) 전체 페이지를 5열 row-major(좌→우, 위→아래)로 병합
# image_montage geometry 이슈를 피하기 위해 수동 병합
n_pages <- length(imgs)
n_cols <- 5
row_tiles <- list()

for (i in seq(1, n_pages, by = n_cols)) {
  row_imgs <- imgs[i:min(i + n_cols - 1, n_pages)]
  row_tiles[[length(row_tiles) + 1]] <- image_append(row_imgs, stack = FALSE)
}

collage <- image_append(do.call(c, row_tiles), stack = TRUE)

# 3) PNG로 저장
image_write(collage, path = outfile, format = "png")


##############################################################
##############################################################
# pre-requisites for submission files
##############################################################
##############################################################
##############################################################
# Check merged-vs-replicate HICCUPS loop files for 592 and 607
##############################################################
# This block is for manual provenance/QC checks before deciding whether
# to use a merged sample or an individual replicate.
medium_resolution_loop_dir <- "/Users/pete/Library/CloudStorage/GoogleDrive-wellclouder@gmail.com/My Drive/medium_resolution_5k10k25k"

replicate_check_files <- tibble::tribble(
  ~sample_family, ~sample_code, ~sample_type, ~loop_file, ~qc_file,
  "592", "592", "merged", file.path(medium_resolution_loop_dir, "592_inter_30.hiccups.5k10k25k", "merged_loops.bedpe"),
  "/Users/pete/Desktop/playground/enhancer/data/QC/592_intact_inter_30.txt",
  "592", "592AA", "replicate", file.path(medium_resolution_loop_dir, "592AA_inter_30.hiccups.5k10k25k", "merged_loops.bedpe"),
  "/Users/pete/Desktop/playground/enhancer/data/QC/592AA_intact_inter_30.txt",
  "592", "592BB", "replicate_used", file.path(medium_resolution_loop_dir, "592BB_inter_30.hiccups.5k10k25k", "merged_loops.bedpe"),
  "/Users/pete/Desktop/playground/enhancer/data/QC/592BB_intact_inter_30.txt",
  "607", "607", "merged_used", file.path(medium_resolution_loop_dir, "607_inter_30.hiccups.5k10k25k", "merged_loops.bedpe"),
  "/Users/pete/Desktop/playground/enhancer/data/QC/607_intact_inter_30.txt",
  "607", "607BB", "replicate", file.path(medium_resolution_loop_dir, "607BB_inter_30.hiccups.5k10k25k", "merged_loops.bedpe"),
  "/Users/pete/Desktop/playground/enhancer/data/QC/607BB_intact_inter_30.txt",
  "607", "607CC", "replicate", file.path(medium_resolution_loop_dir, "607CC_inter_30.hiccups.5k10k25k", "merged_loops.bedpe"),
  "/Users/pete/Desktop/playground/enhancer/data/QC/607CC_intact_inter_30.txt"
)

analysis_loop_files_to_check <- tibble::tribble(
  ~sample_family, ~sample_code, ~analysis_file,
  "592", "592BB", "/Users/pete/dropbox/Gateway_to_Hao/enhancer/data/loops/592BB_merged_loops_5k10k25k.bedpe",
  "607", "607", "/Users/pete/dropbox/Gateway_to_Hao/enhancer/data/loops/607_intact_merged_loops_5k10k25k.bedpe"
)

# sha256 for integrity check
sha256_file <- function(file) {
  if (!file.exists(file)) {
    return(NA_character_)
  }

  str_split(system2("shasum", c("-a", "256", shQuote(file)), stdout = TRUE), "\\s+", simplify = TRUE)[1]
}

read_hiccups_loop_check_bedpe <- function(file) {
  header_line <- readLines(file, n = 1)
  column_names <- str_split(str_remove(header_line, "^#"), "\t", simplify = TRUE) %>% as.character()

  readr::read_tsv(
    file,
    comment = "#",
    col_names = column_names,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  ) %>%
    mutate(
      x1 = as.integer(x1),
      x2 = as.integer(x2),
      y1 = as.integer(y1),
      y2 = as.integer(y2),
      resolution = x2 - x1,
      loop_key = str_c(chr1, x1, x2, chr2, y1, y2, sep = "_")
    )
}

parse_qc_stat_line <- function(lines, pattern) {
  line <- lines[str_detect(lines, fixed(pattern))]
  if (length(line) == 0) {
    return(NA_character_)
  }

  str_squish(str_remove(line[1], paste0("^\\s*", pattern, ":\\s*")))
}

read_juicer_qc_summary <- function(file) {
  if (!file.exists(file)) {
    return(tibble(
      qc_file_exists = FALSE,
      sequenced_read_pairs = NA_character_,
      normal_paired = NA_character_,
      pcr_duplicates = NA_character_,
      library_complexity_estimate = NA_character_,
      hic_contacts = NA_character_
    ))
  }

  lines <- readLines(file, warn = FALSE)
  tibble(
    qc_file_exists = TRUE,
    sequenced_read_pairs = parse_qc_stat_line(lines, "Sequenced Read Pairs"),
    normal_paired = parse_qc_stat_line(lines, "Normal Paired"),
    pcr_duplicates = parse_qc_stat_line(lines, "PCR Duplicates"),
    library_complexity_estimate = parse_qc_stat_line(lines, "Library Complexity Estimate"),
    hic_contacts = parse_qc_stat_line(lines, "Hi-C Contacts")
  )
}

replicate_loop_data <- replicate_check_files %>%
  mutate(
    loop_file_exists = file.exists(loop_file),
    qc_file_exists = file.exists(qc_file),
    sha256 = map_chr(loop_file, sha256_file),
    data = map(loop_file, read_hiccups_loop_check_bedpe)
  )

loop_qc_summary <- replicate_loop_data %>%
  transmute(
    sample_family,
    sample_code,
    sample_type,
    loop_file,
    loop_file_exists,
    sha256,
    n_loops = map_int(data, nrow),
    n_unique_loop_keys = map_int(data, ~ n_distinct(.x$loop_key)),
    resolution_distribution = map_chr(data, ~ .x %>%
      count(resolution, name = "n") %>%
      arrange(resolution) %>%
      mutate(label = str_c(resolution, "=", n)) %>%
      pull(label) %>%
      str_c(collapse = "; ")),
    qc_file
  ) %>%
  bind_cols(map_dfr(replicate_check_files$qc_file, read_juicer_qc_summary))

loop_overlap_summary <- replicate_loop_data %>%
  dplyr::select(sample_family, sample_code, data) %>%
  group_by(sample_family) %>%
  group_modify(~ {
    pair_grid <- t(combn(.x$sample_code, 2)) %>% as_tibble(.name_repair = "minimal")
    names(pair_grid) <- c("sample_a", "sample_b")

    pair_grid %>%
      rowwise() %>%
      mutate(
        loops_a = list(.x$data[[match(sample_a, .x$sample_code)]]$loop_key),
        loops_b = list(.x$data[[match(sample_b, .x$sample_code)]]$loop_key),
        n_a = length(loops_a),
        n_b = length(loops_b),
        n_intersect = length(intersect(loops_a, loops_b)),
        n_a_only = length(setdiff(loops_a, loops_b)),
        n_b_only = length(setdiff(loops_b, loops_a)),
        identical_loop_sets = setequal(loops_a, loops_b)
      ) %>%
      ungroup() %>%
      dplyr::select(-loops_a, -loops_b)
  }) %>%
  ungroup()

analysis_file_identity_summary <- analysis_loop_files_to_check %>%
  mutate(
    analysis_file_exists = file.exists(analysis_file),
    analysis_sha256 = map_chr(analysis_file, sha256_file)
  ) %>%
  left_join(
    loop_qc_summary %>%
      dplyr::select(sample_family, sample_code, source_sha256 = sha256, source_loop_file = loop_file),
    by = c("sample_family", "sample_code")
  ) %>%
  mutate(analysis_file_identical_to_source = analysis_sha256 == source_sha256)

print(loop_qc_summary, n = Inf, width = Inf)
print(loop_overlap_summary, n = Inf, width = Inf)
print(analysis_file_identity_summary, n = Inf, width = Inf)

##############################################################
# submission files
##############################################################
# 1. loops

write.csv(df.final.loop, file = "./figures/submission/lt2mb/df_final_loop_sub.4.any.lt2mb.ENSEMBL.mid.mid.final.filter.200kb.csv", row.names = FALSE)

# 1-1. Raw HICCUPS loop results with functional-loop annotation
# Keep every original BEDPE column and prepend sample metadata plus a functional flag.
loop_file_metadata <- tibble::tribble(
  ~sample_code, ~sample_name, ~file,
  "592BB", "SHR/OlaIpcv", "~/dropbox/Gateway_to_Hao/enhancer/data/loops/592BB_merged_loops_5k10k25k.bedpe",
  "607", "HXB10", "~/dropbox/Gateway_to_Hao/enhancer/data/loops/607_intact_merged_loops_5k10k25k.bedpe",
  "74AA", "F344/Stm", "~/dropbox/Gateway_to_Hao/enhancer/data/loops/74AA_intact_merged_loops5k10k25k.bedpe",
  "A2DB", "LE/Stm", "~/dropbox/Gateway_to_Hao/enhancer/data/loops/A2DB_merged_loops_5k10k25k.bedpe",
  "D765A", "BXH6", "~/dropbox/Gateway_to_Hao/enhancer/data/loops/D765A_intact_merged_loops_5k10k25k.bedpe",
  "DA08A", "HXB2", "~/dropbox/Gateway_to_Hao/enhancer/data/loops/DA08A_intact_merged_loops_5k10k25k.bedpe",
  "DA21A", "SHR/OlaIpcvxBN/NHsdMcwi", "~/dropbox/Gateway_to_Hao/enhancer/data/loops/DA21A_intact_merged_loops_5k10k25k.bedpe",
  "DA68A", "HXB31", "~/dropbox/Gateway_to_Hao/enhancer/data/loops/DA68A_intact_merged_loops_5k10k25k.bedpe",
  "DBA9A", "HXB23", "~/dropbox/Gateway_to_Hao/enhancer/data/loops/DBA9A_intact_merged_loops_5k10k25k.bedpe",
  "DE8BA", "BN-Lx", "~/dropbox/Gateway_to_Hao/enhancer/data/loops/DE8BA_intact_merged_loops_5k10k25k.bedpe"
) %>%
  mutate(file = path.expand(file))
loop_file_metadata

read_hiccups_bedpe <- function(file) {
  header_line <- readLines(file, n = 1)
  column_names <- str_split(str_remove(header_line, "^#"), "\t", simplify = TRUE) %>% as.character()

  readr::read_tsv(
    file,
    comment = "#",
    col_names = column_names,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  )
}

final_loop_ids <- df.final.loop %>%
  distinct(loop.id) %>%
  pull(loop.id)

df.raw.hiccups.loops.with.functional <- loop_file_metadata %>%
  mutate(data = purrr::map(file, read_hiccups_bedpe)) %>%
  dplyr::select(-file) %>%
  tidyr::unnest(data) %>%
  mutate(
    loop.id = str_c(chr1, x1, x2, chr2, y1, y2, as.numeric(x2) - as.numeric(x1), sep = "_"),
    functional = loop.id %in% final_loop_ids
  ) %>%
  dplyr::select(sample_code, sample_name, functional, dplyr::everything(), -loop.id)

df.raw.hiccups.loops.with.functional %>% head(2)

write.csv(
  df.raw.hiccups.loops.with.functional,
  file = "./figures/submission/lt2mb/raw_hiccups_loops_with_functional_annotation.csv",
  row.names = FALSE
)
