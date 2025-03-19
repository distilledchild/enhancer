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

getwd()

# Windows
# setwd('C:\\Users\\panju\\Dropbox (UTHSC GGI)\\Gateway_to_Hao\\workshop\\2023_NIH_meeting\\loop_N_tss')
# source(file.path('C:\\Users\\panju\\Dropbox (UTHSC GGI)\\Gateway_to_Hao\\project_common_code', 'variables.R'))
# source(file.path('C:\\Users\\panju\\Dropbox (UTHSC GGI)\\Gateway_to_Hao\\project_common_code', 'funcs.R'))

getwd()

# Linux
# setwd('C:\\Users\\panju\\Dropbox (UTHSC GGI)\\Gateway_to_Hao\\workshop\\2023_NIH_meeting\\loop_N_tss')
# setwd('./Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss')
# setwd('~/Desktop/temp/enhancer/dropbox_enhancer_doosan/data/enhancer_atlas2.0/all_species/neuron')
setwd('~/Desktop/temp/enhancer/dropbox_enhancer_doosan')
setwd('/home/pkim/dropbox/Gateway_to_Hao/enhancer/r_files')
source(file.path('/home/pkim/dropbox/Gateway_to_Hao/project_common_code/', 'variables.R'))
source(file.path('/home/pkim/dropbox/Gateway_to_Hao/project_common_code/', 'funcs.R'))

getwd()

# Mac
setwd('/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files')
getwd()
source(file.path('/Users/PanjunKim/dropbox/Gateway_to_Hao/project_common_code/', 'variables.R'))
source(file.path('/Users/PanjunKim/dropbox/Gateway_to_Hao/project_common_code/', 'funcs.R'))

# source(file.path('/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/enhancer_synced_data_preparation.R'))
# source(file.path('/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/data_analysis.R'))
# source(file.path('/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/data_analysis.R'))
<<<<<<< HEAD

########################
# 3. TSS
# 3-1. Exploratory Data analysis (EDA) : TSS
# 3-2. TSS data preprocessing: id and dedup GRange Obj.
########################
# tss file list
# Linux
file.tss.list = fs::dir_ls("/home/pkim/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/", regexp = "\\.txt$")
# Mac
file.tss.list = fs::dir_ls("/Users/PanjunKim/UTHSC GGI Dropbox/K P/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss", regexp = "\\.txt$")

file.tss.list
# csRNA.NuAcc.tss.txt
# csRNA.PFC.tss.txt
# ucsc_start_codon.txt

##### tss with UCSC
# get the TSS location into a GenomicRange object, note gene name is added, also there are some duplicated lines, so use uniq
# Linux
tss<-read.table(file="/home/pkim/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/ucsc_start_codon.txt", sep="\t", head=F)
# Mac
tss<-read.table(file="/Users/PanjunKim/UTHSC GGI Dropbox/K P/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/ucsc_start_codon.txt", sep="\t", head=F)
tss
head(tss)[,c(1,4,5,7,9)]
tss.select<-tss[,c(1,4,5,7,9)] # onlty take the relevant columns
tss.select
tss.select.sparate <- tss.select %>% 
  separate(V9, into = c("gene_id", "transcript_id", "exon_number", "exon_id", "gene_name"), sep = "; ") %>%
  mutate(across(everything(), ~ gsub(".+\\s", "", .))) %>% 
  mutate(gene_name = str_replace(gene_name, ';', ''))

tss.select.sparate
names(tss.select.sparate)<-c("chr", "start", "end", "strand","gene_id", "transcript_id", "exon_number", "exon_id", "gene_name")

df.tss.ucsc <- tss.select.sparate %>% 
  mutate(tss.id = paste(chr, start, end, strand, gene_id, transcript_id, sep = ":"))

df.tss.ucsc %>% head()
df.tss.ucsc %>% count() # 17849

# GRanges Obj.: df.tss.ucsc.GR
df.tss.ucsc.GR <- GRanges(
  seqnames = as.character(df.tss.ucsc$chr),
  ranges = IRanges(
    start = as.numeric(df.tss.ucsc$start),
    end = as.numeric(df.tss.ucsc$end)
  ),
  strand = df.tss.ucsc$strand
)

# meta data
mcols(df.tss.ucsc.GR)$tss_id <- df.tss.ucsc$tss.id
mcols(df.tss.ucsc.GR)$gene_id <- df.tss.ucsc$gene_id
mcols(df.tss.ucsc.GR)$transcript_id <- df.tss.ucsc$transcript_id

df.tss.ucsc.GR
##### tss with nuacc, 96563
# nuacc <- read.csv("/home/pkim/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/csRNA.NuAcc.tss.txt", header = T, sep = '\t') %>% 
# nuacc <- read.csv("/Users/PanjunKim/UTHSC GGI Dropbox/K P/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/csRNA.NuAcc.tss.txt", header = T, sep = '\t') %>% 
#   mutate(Anno3 = Detailed.Annotation) %>% 
#   dplyr::select(Chr, Start, End, Annotation, Anno3) %>% 
#   separate(Annotation, into=c("Anno1", "Anno2"), sep = ' ', remove=FALSE) %>%
#   mutate(file = "nuacc") %>% 
#   mutate(Anno2 = str_replace_all(Anno2, "\\(|\\)|,", "")) %>% 
#   # filter(Anno2 == 'NR_132639') %>%
#   # count(Anno2) %>% 
#   view()

##### tss with pfc ver1, 131647
# pfc <- read.csv("/home/pkim/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/csRNA.PFC.tss.txt", header = T, sep = '\t') %>% 
# pfc <- read.csv("/Users/PanjunKim/UTHSC GGI Dropbox/K P/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/csRNA.PFC.tss.txt", header = T, sep = '\t') %>% 
#   dplyr::select(Chr, Start, End, Annotation) %>% 
#   separate(Annotation, into=c("Anno1", "Anno2"), sep = ' ', remove=FALSE) %>%
#   mutate(file = "PFC") %>% 
#   mutate(Anno2 = str_replace_all(Anno2, "\\(|\\)|,", "")) %>% 
#   # count(Anno2) %>%
#   view()

# nuacc + pfc: 228210
# df.tss.nuacc.pfc <- bind_rows(nuacc, pfc) %>% 
#   view()

# genomic range for tss from nuacc, 96563
# df.tss.nuacc.GR <-GRanges(seqnames=nuacc$Chr, ranges=IRanges(start=nuacc$Start, end=nuacc$End), loc=nuacc$Anno1, geneid=nuacc$Anno2)
# df.tss.nuacc.GR

# genomic range for tss from pfc ver1, 131647
# df.tss.pfc.GR <- GRanges(seqnames=pfc$Chr, ranges=IRanges(start=pfc$Start, end=pfc$End), loc=pfc$Anno1, geneid = pfc$Anno2)
# df.tss.pfc.GR

# data.frame(df.tss.ucsc.GR) %>% 
# data.frame(df.tss.pfc.GR) %>% 
# data.frame(df.tss.nuacc.GR) %>% 
# count(seqnames)

# chr19_NW_023637717v1_random
# chrUn_NW_023637854v1
# chrY_NW_023637718v1_random

########################
# 3. TSS
# 3-3. overall distribution of TSS on loops
########################

index.distinct.tss.w.overall.whole.loop <- findOverlaps(
  df.tss.ucsc.GR, 
  overall.df.DISTINCT.loop.deep.sample.all.GR, 
  type = "within",
  select = "all"
)
overall.loop.for.tss.hits <- subjectHits(index.distinct.tss.w.overall.whole.loop)
overall.tss.on.loop.hits <- queryHits(index.distinct.tss.w.overall.whole.loop)

df.tss.dist.result <- tibble(
  loop.id = overall.df.DISTINCT.loop.deep.sample.all$new.loop.id[overall.loop.for.tss.hits],
  loop.start = overall.df.DISTINCT.loop.deep.sample.all$x0[overall.loop.for.tss.hits],
  loop.end = overall.df.DISTINCT.loop.deep.sample.all$y3[overall.loop.for.tss.hits],
  loop.res=overall.df.DISTINCT.loop.deep.sample.all$resolution[overall.loop.for.tss.hits],
  tss_chr = df.tss.ucsc$chr[overall.tss.on.loop.hits],
  tss_start = df.tss.ucsc$start[overall.tss.on.loop.hits],
  tss_end = df.tss.ucsc$end[overall.tss.on.loop.hits],
  tss_id = df.tss.ucsc$tss.id[overall.tss.on.loop.hits],
  tss_geneid = df.tss.ucsc$gene_id[overall.tss.on.loop.hits],
  tss_strand = df.tss.ucsc$strand[overall.tss.on.loop.hits]
) %>% 
  mutate(across(where(is.numeric), ~ format(., scientific = FALSE)))

df.tss.dist.result %>% dim() # 375567      10
df.tss.dist.result %>% head()

relative.pos.df.tss.dist.result <- df.tss.dist.result %>%
  mutate(tss_start = as.numeric(tss_start),
         tss_end = as.numeric(tss_end),
         loop.start = as.numeric(loop.start),
         loop.end = as.numeric(loop.end)) %>% 
  mutate(pos_coord = (tss_start + tss_end) / 2) %>% # coordinate for midpoint
  mutate(loop_length = (loop.end - loop.start)) %>% # overall length of the loop
  mutate(relative_pos = pos_coord - loop.start) %>% # relative pos
  # mutate(value = (relative_pos/(loop_length / 2)) - 1) %>% for .5 distance
  mutate(value = (relative_pos/(loop_length / 3)) - 1) %>% # 1 distance
  mutate(resolution = case_when(
    loop.res == 5000 ~ "5K",
    loop.res  == 10000 ~ "10K",
    loop.res  == 25000 ~ "25K",
    TRUE ~ NA_character_ 
  )) %>%
  dplyr::select(loop.id, tss_geneid, value, loop.res)

relative.pos.df.tss.dist.result %>% head()

########################
# 3. TSS
# 3-3. overall distribution of TSS on loops: figures
# 3-3-1. by CHROMOSOME
########################

chromosomes <- c(1:20, "X", "Y")

plot_tss_list <- list()

for (chr in chromosomes) {
  tryCatch({
    
    message("START: Processing chromosome: ", chr)
    
    relative.pos.df.tss.dist.result.chr <- relative.pos.df.tss.dist.result %>% 
      filter(str_detect(loop.id, paste0("_chr", chr, "_")))
    
    plot1.tss.dens <- relative.pos.df.tss.dist.result.chr %>% 
      ggplot(aes(x = value)) +
      geom_density(fill = "skyblue", color = "black", alpha = 0.7) +
      ylim(c(0,1))+
      labs(title = paste0("Density of TSS Found over Loops on Chr", chr),
           x = "Relative position to loop",
           y = "Density"
      ) + 
      theme(plot.title = element_text(hjust = 0.5))
    
    plot1.tss.hist <- relative.pos.df.tss.dist.result.chr %>% 
      ggplot(aes(x = value)) +
      geom_histogram(fill = "skyblue", color = "black", alpha = 0.7, bins=200) +
      labs(title = paste0("Histogram of CTCF found over Loops on Chr", chr),
           x = "Relative position to loop",
           y = "Counts"
      ) + 
      theme(plot.title = element_text(hjust = 0.5))
    
    combined_plot_tss <- plot1.tss.hist + plot1.tss.dens # + can be used instead of |
    plot_tss_list[[chr]] <- combined_plot_tss
    
    message("END: Processing chromosome: ", chr)
  }, error = function(e) {
    
    message("Error processing chromosome: ", chr)
    message("Error message: ", e$message)
  })
}

pdf("./figures/submission/overall_distribution_of_TSS_by_chromosome.pdf", width = 11, height = 8.5)
tss_num_plots <- length(plot_tss_list) # 22 chromosomes
tss_plots_per_page <- 4

for (i in seq(1, tss_num_plots, by = tss_plots_per_page)) {
  end_idx <- min(i + tss_plots_per_page - 1, tss_num_plots)
  page_plots <- plot_tss_list[i:end_idx]
  
  combined_page <- plot_grid(plotlist = page_plots, ncol = 1, nrow = 4)
  
  print(combined_page)
}

dev.off()

########################
# 3. TSS
# 3-3. overall distribution of TSS on loops: figures
# 3-3-2. by resolution
########################
# TSS Histogram (all resolution combined, left)
plot.tss.hist <- ggplot(relative.pos.df.tss.dist.result, aes(x = value)) +
  geom_histogram(fill = "skyblue", color = "black", alpha = 0.7, bins=200) +
  labs(title = paste0("Histogram of TSS over Loops on ALL Chromosomes (", res, ")"),
       x = "Relative position to loop",
       y = "Counts"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

plot.tss.hist

# TSS Density Plot (by resolution combined, right)
plot.tss.dens <- ggplot(relative.pos.df.tss.dist.result, aes(x = value, color = loop.res, fill = loop.res)) +
  geom_density(alpha = 0.3) +  # transparent
  ylim(c(0, 1)) +
  scale_color_manual(values = c("5K" = "#619CFF", "10K" = "#F8766D", "25K" = "#00BA38")) + 
  scale_fill_manual(values = c("5K" = "#619CFF", "10K" = "#F8766D", "25K" = "#00BA38")) +
  labs(title = "Density of TSS Found over Loops",
       x = "Relative position to loop",
       y = "Density",
       color = "Resolution",
       fill = "Resolution") + 
  theme(plot.title = element_text(hjust = 0.5))
plot.tss.dens

pdf("figures/submission/overall_distribution_of_TSS_by_resolution.pdf", width = 11, height = 5)

final_tss_plot <- plot.tss.hist | plot.tss.dens 
print(final_tss_plot)

dev.off()

# plot_tss_all_res <- create_tss_plots(relative.pos.df.tss.dist.result, "All")
# plot_tss_5K <- create_tss_plots(relative.pos.df.tss.dist.result, "5K")
# plot_tss_10K <- create_tss_plots(relative.pos.df.tss.dist.result, "10K")
# plot_tss_25K <- create_tss_plots(relative.pos.df.tss.dist.result, "25K")
# pdf("./figures/submission/overall_distribution_tss_on_all_chromosomes_1_distance.pdf", width = 11, height = 8.5)
# (plot_tss_all_res/plot_tss_5K)
# (plot_tss_10K / plot_tss_25K)
# dev.off()

########################
# 3. TSS
# 3-4. Distribution of TSS at each end in a loop: for the number of TSS used in filtering valid loops: figures
# so, the object should be used one WITH padding on df.DISTINCT.loop.deep.sample.all: 1/2 distance for OUTER & 1/4 distance for INNER
# 3-4-1. adding padding at each end, x12 & y12
########################
# df.DISTINCT.loop.deep.sample.all

df.DISTINCT.loop.deep.sample.all.padded.for.TSS <- df.DISTINCT.loop.deep.sample.all %>% 
  mutate(x12 = (x1 + x2)/2, # middle point of each end
         y12 = (y1 + y2)/2, # middle point of each end
         distance = y12 - x12,
         x0 = ifelse(x12 - (distance / 2) < 0, 0, x12 - (distance / 2)), 
         y3 = y12 + (distance/2), # 1/2 distance for OUTER PADDING
         x12 = (x12 + (distance / 4)), 
         y12 = ifelse((y12 - (distance/4)) < 0, 0, y12 - (distance/4)), # 1/4 distance for INNER PADDING
         new.loop.id = paste0(chr1, '_', x0, '_', x12, '_', chr2, '_', y12, '_', y3, '_', end.distance))

# GRange Obj. and meta data
df.DISTINCT.loop.deep.sample.all.padded.for.TSS.UP.GR <- GRanges(
  seqnames = df.DISTINCT.loop.deep.sample.all.padded.for.TSS$chr1,
  ranges = IRanges(
    start = df.DISTINCT.loop.deep.sample.all.padded.for.TSS$x0, 
    end = df.DISTINCT.loop.deep.sample.all.padded.for.TSS$x12
  )
)

mcols(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.UP.GR) <- df.DISTINCT.loop.deep.sample.all.padded.for.TSS %>%
  dplyr::select(loop.id, end.distance, resolution, new.loop.id)

########################
# 3. TSS
# 3-4. Distribution of TSS at ends in a loop: for the number of TSS used in filtering valid loops: figures
# so, the object should be used one WITH padding: df.DISTINCT.loop.deep.sample.all.padded.for.TSS
# 3-4-3. data processing for getting information of TSS at ends in loops 
########################
# 1. TSS + UPSTREAM (df.DISTINCT.loop.deep.sample.all.padded.for.TSS.UP.GR, df.tss.ucsc.GR)
index.distinct.tss.w.up.loop.each.end <- findOverlaps(
  df.tss.ucsc.GR, 
  df.DISTINCT.loop.deep.sample.all.padded.for.TSS.UP.GR, 
  type = "within",
  select = "all"
)

end.loop.up.tss.hits <- subjectHits(index.distinct.tss.w.up.loop.each.end)
end.tss.up.hits <- queryHits(index.distinct.tss.w.up.loop.each.end)
df.tss.ucsc.GR
df.overlapping.TSS.w.UPSTREAM.result.each.end <- tibble(
  up.loop.id = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.UP.GR)$loop.id[end.loop.up.tss.hits],
  end.up.distance = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.UP.GR)$end.distance[end.loop.up.tss.hits],
  resolution = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.UP.GR)$resolution[end.loop.up.tss.hits],
  tss.id = mcols(df.tss.ucsc.GR)$tss_id[end.tss.up.hits],
  WHERE = "UP"
) %>% 
  mutate(tss.loop.up.id = str_c(up.loop.id, '|', tss.id, '|', WHERE))

df.overlapping.TSS.w.UPSTREAM.result.each.end %>% dim() # 95831
df.overlapping.TSS.w.UPSTREAM.result.each.end %>% head()

df.overlapping.TSS.w.UPSTREAM.result.each.end %>% 
  count(end.up.distance)

# end.up.distance     n
# 1            5000 21393
# 2           10000 28339
# 3           25000 46099
# total: 95831

# 2. TSS + DOWNSTREAM
df.DISTINCT.loop.deep.sample.all.padded.for.TSS.DOWN.GR <- GRanges(
  seqnames=df.DISTINCT.loop.deep.sample.all.padded.for.TSS$chr1, 
  ranges=IRanges(
    start=df.DISTINCT.loop.deep.sample.all.padded.for.TSS$y12, 
    end=df.DISTINCT.loop.deep.sample.all.padded.for.TSS$y3
  )
)
mcols(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.DOWN.GR) <- df.DISTINCT.loop.deep.sample.all.padded.for.TSS %>%
  dplyr::select(loop.id, end.distance, resolution, new.loop.id)

# 1. TSS + DOWNSTREAM (df.DISTINCT.loop.deep.sample.all.padded.for.TSS.DOWN.GR, df.tss.ucsc.GR)
index.distinct.tss.w.down.loop.each.end <- findOverlaps(
  df.tss.ucsc.GR, 
  df.DISTINCT.loop.deep.sample.all.padded.for.TSS.DOWN.GR, 
  type = "within",
  select = "all"
)
end.loop.down.tss.hits <- subjectHits(index.distinct.tss.w.down.loop.each.end)
end.tss.down.hits.end <- queryHits(index.distinct.tss.w.down.loop.each.end)

# df.DISTINCT.loop.deep.sample.all.padded.for.TSS %>% 
df.tss.ucsc.GR %>% 
  head()

df.overlapping.TSS.w.DOWNSTREAM.result.each.end <- tibble(
  down.loop.id = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.UP.GR)$loop.id[end.loop.down.tss.hits],
  end.down.distance = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.UP.GR)$end.distance[end.loop.down.tss.hits],
  resolution = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.UP.GR)$resolution[end.loop.down.tss.hits],
  tss.id = mcols(df.tss.ucsc.GR)$tss_id[end.tss.down.hits.end],
  WHERE = "DOWN"
) %>%
  mutate(tss.loop.down.id = str_c(down.loop.id, '|', tss.id, '|', WHERE))

df.overlapping.TSS.w.DOWNSTREAM.result.each.end %>% dim() # 111334
df.overlapping.TSS.w.DOWNSTREAM.result.each.end %>% head()
df.overlapping.TSS.w.DOWNSTREAM.result.each.end %>% 
  count(end.down.distance)

# end.down.distance     n
# 1              5000 25473
# 2             10000 31699
# 3             25000 54162
# total: 111334

########## bind_rows(UPSTREAM & DOWNSTREAM) -> BOTH
df.overlapping.TSS.w.BOTH.result <- bind_rows(df.overlapping.TSS.w.UPSTREAM.result.each.end %>% 
                                                mutate(loop.id = up.loop.id, end.distance = end.up.distance) %>% 
                                                mutate(case.id = tss.loop.up.id) %>% 
                                                dplyr::select(-c(up.loop.id, end.up.distance, tss.loop.up.id)), 
                                              df.overlapping.TSS.w.DOWNSTREAM.result.each.end %>% 
                                                mutate(loop.id = down.loop.id, end.distance = end.down.distance) %>% 
                                                mutate(case.id = tss.loop.down.id) %>% 
                                                dplyr::select(-c(down.loop.id, end.down.distance, tss.loop.down.id))) %>% 
  mutate(chr = ifelse(WHERE == "UP", str_split_n(loop.id, '_', 1), str_split_n(loop.id, '_', 4))) %>% 
  mutate(chr = factor(chr, levels = c(paste0("chr", 1:20), "chrX", "chrY"))) %>%
  mutate(WHERE = fct_relevel(WHERE, "UP", "DOWN")) %>% 
  mutate(resolution = case_when(
    end.distance == 5000 ~ "5K",
    end.distance == 10000 ~ "10K",
    end.distance == 25000 ~ "25K",
    TRUE ~ NA
  )) %>% 
  mutate(resolution = fct_relevel(resolution, "5K", "10K", "25K"))

df.overlapping.TSS.w.BOTH.result %>% dim() # 207165
df.overlapping.TSS.w.BOTH.result %>% head()
df.overlapping.TSS.w.BOTH.result %>% count(resolution)

# resolution      n
# 1         5K 21393 + 25473 =  46866
# 2        10K 28339 + 31699 =  60038
# 3        25K 46099 + 54162 = 100261

########################
# 3. TSS
# 3-4. Distribution of TSS at ends in a loop: for the number of TSS used in filtering valid loops: figures
# so, the object should be used one WITH padding: df.DISTINCT.loop.deep.sample.all
# 3-4-4. Distribution of TSS at both ends in a loop: figure
########################
# for boxplot
df.overlapping.TSS.w.BOTH.result %>% head(2) # resolution, tss.id, WHERE, loop.id, end.distance, case.id, chr

df.overlapping.TSS.w.BOTH.result.boxplot <- df.overlapping.TSS.w.BOTH.result %>% 
  group_by(chr, loop.id, WHERE, resolution) %>%
  # group_by(loop.id, WHERE, resolution) %>%
  # group_by(chr, loop.id, resolution) %>%
  summarise(tss_count_by_loop_id_each_end = n_distinct(tss.id), .groups = 'drop')

df.overlapping.TSS.w.BOTH.result.boxplot

# for Q1
df.overlapping.TSS.w.BOTH.result %>% head()

df.tss.counts.each.end <- df.overlapping.TSS.w.BOTH.result %>%
  group_by(loop.id, WHERE, resolution) %>%
  summarise(tss_count_each_end = n_distinct(tss.id), .groups = 'drop')

df.tss.counts.each.end

tss.stats.by.resolution.each.end <- df.tss.counts.each.end %>%
  # group_by(resolution) %>%
  group_by(WHERE, resolution) %>%
  summarise(
    Min = min(tss_count_each_end),
    Q1 = quantile(tss_count_each_end, 0.25, na.rm = TRUE),
    Median = median(tss_count_each_end, na.rm = TRUE),
    Q3 = quantile(tss_count_each_end, 0.75, na.rm = TRUE),
    Mean = mean(tss_count_each_end, na.rm = TRUE),
    SD = sd(tss_count_each_end, na.rm = TRUE),
    Max = max(tss_count_each_end),
    .groups = 'drop'
  )

tss.stats.by.resolution.each.end

# resolution   Min    Q1 Median    Q3  Mean    SD   Max
# <fct>      <int> <dbl>  <dbl> <dbl> <dbl> <dbl> <int>
# 1 5K             1     1      2     3  5.33  35.8  1282
# 2 10K            1     1      2     3  3.59  18.0   993
# 3 25K            1     1      2     4  5.26  28.3  1282

# WHERE resolution   Min    Q1 Median    Q3  Mean    SD   Max
# <fct> <fct>      <int> <dbl>  <dbl> <dbl> <dbl> <dbl> <int>
# 1 UP    5K             1     1      2     3  4.80  28.7   921
# 2 UP    10K            1     1      2     3  3.37  12.4   428
# 3 UP    25K            1     1      2     4  4.84  21.5   921
# 4 DOWN  5K             1     1      2     3  5.87  41.8  1282
# 5 DOWN  10K            1     1      2     3  3.82  22.2   993
# 6 DOWN  25K            1     1      2     4  5.68  33.7  1282

# quantile for EACH END// NOT in a loop : results are loops have ctcfs more than 18 at least in each end
q3_tss_count_each_end <- quantile(df.overlapping.TSS.w.BOTH.result.boxplot$tss_count_by_loop_id_each_end, 0.75, na.rm = TRUE)
q3_tss_count_each_end # 3 :: integrity PASS

q1_tss_count_each_end <- quantile(df.overlapping.TSS.w.BOTH.result.boxplot$tss_count_by_loop_id_each_end, 0.25, na.rm = TRUE)
q1_tss_count_each_end # 1 :: integrity PASS

# drawing boxplot
df.overlapping.TSS.w.BOTH.result.boxplot
df.overlapping.TSS.w.BOTH.result.boxplot %>% head(4)


# figure by chr & res
boxplot.w.TSS.by.chr.and.res.each.end <- df.overlapping.TSS.w.BOTH.result.boxplot %>% 
  ggplot(aes(x = chr, y = tss_count_by_loop_id_each_end, fill = resolution)) +
  geom_boxplot() +
  stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 1)), 
               position = position_dodge(width = 0.75), vjust = -0.5, size = 2, color = "black") +
  stat_summary(fun.data = function(y) {
    data.frame(
      y = quantile(y, 0.25),
      label = round(quantile(y, 0.25), 1)
    )
  }, geom = "text", aes(label = after_stat(label)), position = position_dodge(width = 0.75), vjust = 2.5, size = 2, color = "blue") +
  stat_summary(fun.data = function(y) {
    data.frame(
      y = quantile(y, 0.75),
      label = round(quantile(y, 0.75), 1)
    )
  }, geom = "text", aes(label = after_stat(label)), position = position_dodge(width = 0.75), vjust = -2.5, size = 2, color = "blue") +
  stat_summary(fun.data = function(y) {
    data.frame(
      y = min(y),
      label = round(min(y), 1)
    )
  }, geom = "text", aes(label = after_stat(label)), position = position_dodge(width = 0.75), vjust = 1.5, hjust = -0.5, size = 2, color = "red") +
  stat_summary(fun.data = function(y) {
    data.frame(
      y = max(y),
      label = round(max(y), 1)
    )
  }, geom = "text", aes(label = after_stat(label)), position = position_dodge(width = 0.75), vjust = -1.5, hjust = -0.5, size = 2, color = "red") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("5K" = "skyblue", "10K" = "lightcoral", "25K" = "springgreen"), breaks = c("5K", "10K", "25K")) +
  scale_y_continuous(limits = c(0, 2 * q3_tss_count_each_end), 
                     breaks = seq(0, 2 * q3_tss_count_each_end, by = 1)) +
  labs(title = "Boxplot for Number of TSS by Chromosome and Resolution", x = "Chromosomes", y = "Number of TSS in a loop")

boxplot.w.TSS.by.chr.and.res.each.end

pdf("figures/submission/end_boxplot_num_tss_by_chr_res.pdf", width = 16.5, height = 23.5)
grid.arrange(boxplot.w.TSS.by.chr.and.res.each.end, ncol = 1)
dev.off()

df.overlapping.TSS.w.BOTH.result.boxplot %>% 
  head()

boxplot.w.TSS.by.res <- df.overlapping.TSS.w.BOTH.result.boxplot %>% 
  ggplot(aes(x = resolution, y = tss_count_by_loop_id_each_end, fill = resolution)) +
  geom_boxplot() +
  stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 1)), 
               position = position_dodge(width = 0.75), vjust = -0.5, size = 2) + 
  stat_summary(fun.data = function(y) {
    data.frame(y = quantile(y, 0.25), label = paste0("Q1: ", round(quantile(y, 0.25), 1)))
  }, geom = "text", aes(label = after_stat(label)), 
  position = position_dodge(width = 0.75), vjust = 6.5, size = 1.5) +
  stat_summary(fun.data = function(y) {
    data.frame(y = quantile(y, 0.75), label = paste0("Q3: ", round(quantile(y, 0.75), 1)))
  }, geom = "text", aes(label = after_stat(label)), 
  position = position_dodge(width = 0.75), vjust = -6.5, size = 2, color = "blue") +
  stat_summary(fun.data = function(y) {
    data.frame(y = min(y), label = paste0("Min: ", round(min(y), 1)))
  }, geom = "text", aes(label = after_stat(label)), 
  position = position_dodge(width = 0.75), vjust = 1.5, size = 2, color = "red") +
  stat_summary(fun.data = function(y) {
    data.frame(y = max(y), label = paste0("Max: ", round(max(y), 1)))
  }, geom = "text", aes(label = after_stat(label)), 
  position = position_dodge(width = 0.75), vjust = -1.5, size = 2, color = "red") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("5K" = "skyblue", "10K" = "lightcoral", "25K" = "springgreen"), 
                    breaks = c("5K", "10K", "25K")) +
  scale_y_continuous(limits = c(0, 2 * q3_tss_count_each_end), 
                     breaks = seq(0, 2 * q3_tss_count_each_end, by = 10)) +
  labs(title = "Boxplot for Number of TSS by Resolution", x = "Resolution", y = "Number of TSS in a loop")

boxplot.w.TSS.by.res

pdf("./figures/submission/end_boxplot_num_tss_by_res.pdf", width = 16.5, height = 23.5)
grid.arrange(boxplot.w.TSS.by.res, ncol = 1)
dev.off()

pdf("./figures/submission/boxplot_of_tss_count.pdf", width = 16.5, height = 23.5)
grid.arrange(boxplot.w.TSS.by.chr.and.res.each.end, boxplot.w.TSS.by.res, ncol = 1)
dev.off()

tss.stats.by.resolution.each.end %>% head()

########################
# 3. TSS
# 3-4. Distribution of TSS at ends in a loop: for the number of TSS used in filtering valid loops: figures
# so, the object should be used one without padding: df.DISTINCT.loop.deep.sample.all
# 3-4-5. Distribution of TSS at each end in a loop: figure
########################
# figure by end
boxplot.w.TSS.ALL.chr <- ggplot(df.overlapping.TSS.w.BOTH.result.boxplot, aes(x = WHERE, y = tss_count_by_loop_id_each_end, fill = WHERE)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 1)), 
               position = position_dodge(width = 0.75), vjust = -0.5, size = 2.5) +
  stat_summary(fun.data = function(y) {
    data.frame(y = median(y), label = paste0("Q1: ", round(quantile(y, 0.25), 1), "\nQ3: ", round(quantile(y, 0.75), 1)))
  }, geom = "text", aes(label = after_stat(label)), position = position_dodge(width = 0.75), vjust = 1.5, size = 2.5) +  # font-size
  theme_minimal() + 
  theme(
    plot.title = element_text(hjust = 0.5),  # title center
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("UP" = "yellow", "DOWN" = "purple")) +
  scale_y_continuous(limits = c(0, 2 * q3_tss_count_each_end),  # y축을 Q3의 두 배로 제한
                     breaks = seq(0, 2 * q3_tss_count_each_end, by = 5)) +
  labs(title = "Boxplot by UP/DOWNSTREAM End", x = "Upstream and Downstream End", y = "Number of CTCF inside each ends in a loop")

boxplot.w.TSS.ALL.chr

# figure by end & res
boxplot.w.TSS.ALL.chr.by.res <- ggplot(df.overlapping.TSS.w.BOTH.result.boxplot, aes(x = WHERE, y = tss_count_by_loop_id_each_end, fill = resolution)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 1)), 
               position = position_dodge(width = 0.75), vjust = -0.5, size = 2.5) +
  stat_summary(fun.data = function(y) {
    data.frame(y = median(y), label = paste0("Q1: ", round(quantile(y, 0.25), 1), "\nQ3: ", round(quantile(y, 0.75), 1)))
  }, geom = "text", aes(label = after_stat(label)), position = position_dodge(width = 0.75), vjust = 1.5, size = 2.5) +  # font-size
  stat_summary(fun.data = function(y) {
    data.frame(y = min(y), label = paste0("Min: ", round(min(y), 1)))
  }, geom = "text", aes(label = after_stat(label)), position = position_dodge(width = 0.75), vjust = 1.5, hjust = -0.2, size = 2.5, color = "red") +
  theme_minimal() + 
  theme(
    plot.title = element_text(hjust = 0.5),  # title center
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("5K" = "skyblue", "10K" = "lightcoral", "25K" = "springgreen"), breaks = c("5K", "10K", "25K")) +
  scale_y_continuous(limits = c(0, 2 * q3_tss_count_each_end), 
                     breaks = seq(0, 2 * q3_tss_count_each_end, by = 5)) +
  labs(title = "Boxplot for Number of TSS in UP/DOWNSTREAM End by Resolution", x = "Upstream and Downstream End", y = "Number of TSS", fill = "Resolution")

boxplot.w.TSS.ALL.chr.by.res

pdf("./figures/submission/boxplot_for_no.tss_by_end_res.pdf", width = 16.5, height = 23.5)
grid.arrange(boxplot.w.TSS.ALL.chr, boxplot.w.TSS.ALL.chr.by.res, ncol = 1)
dev.off()

###############################################
# (1/2) Histogram of TSS Counts per Loop (Exclusive)
###############################################

df.overlapping.TSS.w.BOTH.result.boxplot %>% head()

tss_count_histogram_data <- df.overlapping.TSS.w.BOTH.result.boxplot %>%
  group_by(loop.id) %>%
  summarise(tss_count_total = sum(tss_count_by_loop_id), .groups = 'drop')

total_distinct_loops <- n_distinct(tss_count_histogram_data$loop.id)

more_than_5_loops <- tss_count_histogram_data %>%
  filter(tss_count_total > 5) %>%
  nrow()

more_than_5_loops_percent <- round(more_than_5_loops / total_distinct_loops * 100, 1)

tss.count.per.loop.histogram <- ggplot(tss_count_histogram_data, aes(x = tss_count_total)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  scale_x_continuous(limits = c(0, 6.5), breaks = seq(0, 6, by = 1)) +
  labs(title = "Histogram of TSS Counts per Loop (Exclusive)", 
       x = "Total TSS Count per Loop", 
       y = "Number of Loops") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = rel(1.2))) +
  stat_count(aes(label = paste0(..count.., " (", round(..count../total_distinct_loops * 100, 1), "%)")), 
             geom = "text", 
             vjust = -0.5, size = 3) +
  annotate("text", x = 5, y = max(table(tss_count_histogram_data$tss_count_total)) * 1.1, 
           label = paste("Total Loops:", total_distinct_loops), size = 4, hjust = 1) +
  annotate("text", x = 5, y = max(table(tss_count_histogram_data$tss_count_total)) * 1.05, 
           label = paste("Loops > 5 TSS: ", more_than_5_loops, " (", more_than_5_loops_percent, "%)"), 
           size = 3, hjust = 1, color = "red")

tss.count.per.loop.histogram

df.overlapping.TSS.w.BOTH.result.boxplot %>% head()

###############################################
# (2/2) Histogram of Loops with TSS Concentrated in UP or DOWN
###############################################
tss.count.per.loop <- df.overlapping.TSS.w.BOTH.result.boxplot %>%
  group_by(loop.id, WHERE) %>%
  summarise(total_tss = sum(tss_count_by_loop_id), .groups = 'drop')

df.overlapping.TSS.w.BOTH.result.boxplot %>% head()
tss.count.per.loop %>% head()

# loops which have only in an end
only_up_or_down_loops <- tss.count.per.loop %>%
  group_by(loop.id) %>%
  filter(n() == 1) %>%
  ungroup()
only_up_or_down_loops %>% head()

# counting the loops above
histogram.data <- only_up_or_down_loops %>%
  group_by(total_tss) %>%
  summarise(one_sided_tss_loop_count = n(), .groups = 'drop')
# histogram.data %>% view()

tss.sum.per.loop <- tss.count.per.loop %>%
  group_by(loop.id) %>%
  summarise(total_tss_sum = sum(total_tss), .groups = 'drop')
tss.sum.per.loop %>% head()

tss.count.by.sum <- tss.sum.per.loop %>%
  group_by(total_tss_sum) %>%
  summarise(loop_count = n(), .groups = 'drop')
histogram.data %>% head()
tss.count.by.sum %>% head()

one.sided.loop.histogram.data.final <- right_join(histogram.data, tss.count.by.sum, by = c("total_tss" = "total_tss_sum")) %>%
  mutate(ratio = one_sided_tss_loop_count / loop_count)

one.sided.loop.histogram.data.final %>% head()

one.sided.loop.histogram <- ggplot(one.sided.loop.histogram.data.final, aes(x = total_tss, y = one_sided_tss_loop_count)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  scale_x_continuous(limits = c(0, 6.5), breaks = seq(0, 6, by = 1)) +
  labs(title = "Histogram of Loops with TSS Concentrated in UP or DOWN",
       x = "Total TSS Count per Loop",
       y = "Number of Loops") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = paste0(one_sided_tss_loop_count, "/", loop_count, " (", round(ratio * 100, 1), "%)")),
            vjust = -0.5, size = 3, color = "black")

one.sided.loop.histogram
one.sided.loop.histogram.data.final %>% view()

pdf("./figures/0909/tss_count_per_loop_histogram.pdf", width = 16.5, height = 23.5)

grid.arrange(tss.count.per.loop.histogram, one.sided.loop.histogram, ncol = 1)
dev.off()

=======
>>>>>>> d754844 (changing variables, generating figures for submission, splitting file of enhancer_synced.R into three.)
