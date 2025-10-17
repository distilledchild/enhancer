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

########################
# 4. promoter
# 4-1. Exploratory Data analysis (EDA) : promoter
# 4-2. promoter data preprocessing: id and dedup GRange Obj.
########################
file_path <- "/Users/PanjunKim/UTHSC GGI Dropbox/K P/Gateway_to_Hao/enhancer/data/epdnew/001/Rn_EPDnew_001_rn7.bed"
file_path <- "/Users/PanjunKim/UTHSC GGI Dropbox/K P/Gateway_to_Hao/enhancer/data/epdnew/001/Rn_EPDnew_001_rn7_1.bed"
df.promoter.rn7.raw <- read_tsv(file_path, col_names = c("chr", "start", "end", "gene", "score"))
df.promoter.rn7.raw %>% head()

df.promoter.rn7 <- df.promoter.rn7.raw %>% 
  mutate(length = end - start, 
         center = round((end + start)/2), 
         check.id = str_c(chr, '_', start, '_', end),
         promoter.id = str_c(check.id, '_', gene),)

# dedup
df.promoter.rn7 %>% head()
df.promoter.rn7 %>% dim() # [1] 12463     9
df.promoter.rn7 %>% distinct(check.id) %>% dim() # 12427     1
df.promoter.rn7 %>% distinct(promoter.id) %>% dim() # 12463     1
df.promoter.rn7 %>% filter(check.id %in% check.id[duplicated(check.id)])

# dedup GRange Obj.
df.promoter.rn7.GR <- GRanges(
  seqnames = df.promoter.rn7$chr,
  ranges = IRanges(
    start = df.promoter.rn7$start, 
    end = df.promoter.rn7$end)
)

# metadata
mcols(df.promoter.rn7.GR) <- df.promoter.rn7[, c("gene", "promoter.id", "center", "length")]

########################
# 4. promoter
# 4-3. overall distribution of promoter on loops
########################
index.promoter.w.overall.whole.loop <- findOverlaps(
  df.promoter.rn7.GR, 
  overall.df.DISTINCT.loop.deep.sample.all.GR, 
  type = "within",
  select = "all"
)

overall.loop.for.promoter.hits <- subjectHits(index.promoter.w.overall.whole.loop)
overall.promoter.on.loop.hits <- queryHits(index.promoter.w.overall.whole.loop)

df.promoter.dist.result <- data.frame(
  loop.id = overall.df.DISTINCT.loop.deep.sample.all$new.loop.id[overall.loop.for.promoter.hits],
  loop.start = overall.df.DISTINCT.loop.deep.sample.all$x0[overall.loop.for.promoter.hits],
  loop.end = overall.df.DISTINCT.loop.deep.sample.all$y3[overall.loop.for.promoter.hits],
  loop.res = overall.df.DISTINCT.loop.deep.sample.all$resolution[overall.loop.for.promoter.hits],
  # loop.new.distance = overall.df.DISTINCT.loop.deep.sample.all$new_distance[overall.loop.for.promoter.hits],
  promoter_id = df.promoter.rn7$promoter.id[overall.promoter.on.loop.hits],
  promoter_start = df.promoter.rn7$start[overall.promoter.on.loop.hits],
  promoter_end = df.promoter.rn7$end[overall.promoter.on.loop.hits],
  promoter_gene = df.promoter.rn7$gene[overall.promoter.on.loop.hits],
  promoter_length = df.promoter.rn7$length[overall.promoter.on.loop.hits]
) %>% 
  mutate(across(where(is.numeric), ~ format(., scientific = FALSE)))

df.promoter.dist.result %>% head()
df.promoter.dist.result %>% dim() # 265565      9(1 distance)

relative.pos.df.promoter.dist.result <- df.promoter.dist.result %>%
  mutate(promoter_start = as.numeric(promoter_start),
         promoter_end = as.numeric(promoter_end),
         loop.start = as.numeric(loop.start),
         loop.end = as.numeric(loop.end)) %>% 
  mutate(pos_coord = round((promoter_start + promoter_end) / 2)) %>%
  mutate(loop_length = (loop.end - loop.start)) %>% # x0, y3
  mutate(relative_pos = pos_coord - loop.start) %>% # relative position from loop start
  # mutate(value = (relative_pos / (loop_length / 2)) - 1) %>% # for padding .5x distance
  mutate(value = relative_pos/(loop_length/3)-1) %>% # for padding 1x distance
  mutate(resolution = case_when(
    loop.res == 5000 ~ "5K",
    loop.res  == 10000 ~ "10K",
    loop.res  == 25000 ~ "25K",
    TRUE ~ NA_character_ 
  )) %>%
  dplyr::select(loop.id, promoter_id, value, loop.res, promoter_gene)
relative.pos.df.promoter.dist.result
relative.pos.df.promoter.dist.result %>% 
  dim() # 265565      5

########################
# 4. promoter
# 4-3. overall distribution of promoter on loops: figures
# 4-3-1. by CHROMOSOME
########################
plot_promoter_list <- list()

for (chr in chromosomes) {
  tryCatch({
    
    message("START: Processing chromosome: ", chr)
    
    relative.pos.df.promoter.dist.result.chr <- relative.pos.df.promoter.dist.result %>% 
      filter(str_detect(loop.id, paste0("_chr", chr, "_")))
    
    plot1.promoter.dens <- relative.pos.df.promoter.dist.result.chr %>% 
      ggplot(aes(x = value)) +
      geom_density(fill = "skyblue", color = "black", alpha = 0.7) +
      ylim(c(0,1))+
      labs(title = paste0("Density of Promoter Found near Loops on chr", chr),
           x = "Relative position to loop",
           y = "Density"
      )
    
    plot1.promoter.hist <- relative.pos.df.promoter.dist.result.chr %>% 
      ggplot(aes(x = value)) +
      geom_histogram(fill = "skyblue", color = "black", alpha = 0.7, bins=200) +
      labs(title = paste0("Histogram of Promoter found Loops on Chr", chr),
           x = "Relative position to loop",
           y = "Counts"
      )
    
    combined_plot <- plot1.promoter.hist + plot1.promoter.dens # + can be used instead of |
    plot_promoter_list[[chr]] <- combined_plot
    
    message("END: Processing chromosome: ", chr)
  }, error = function(e) {
    
    message("Error processing chromosome: ", chr)
    message("Error message: ", e$message)
  })
}

pdf("./figures/submission/overall_distribution_of_promoter_by_chromosome.pdf", width = 11, height = 8.5)

promoter_num_plots <- length(plot_promoter_list) # 22 chromosomes
promoter_plots_per_page <- 4

for (i in seq(1, promoter_num_plots, by = promoter_plots_per_page)) {
  end_idx <- min(i + promoter_plots_per_page - 1, promoter_num_plots)
  page_plots <- plot_promoter_list[i:end_idx]
  
  combined_page <- plot_grid(plotlist = page_plots, ncol = 1, nrow = 4)
  
  print(combined_page)
}

dev.off()

########################
# 4. promoter
# 4-3. overall distribution of promoter on loops: figures
# 4-3-2. by resolution
########################
# promoter Histogram (all resolution combined, left)
plot.promoter.hist <- relative.pos.df.promoter.dist.result %>% 
  ggplot(aes(x = value)) +
  geom_histogram(fill = "skyblue", color = "black", alpha = 0.7, bins=200) +
  labs(title = paste0("Histogram of Promoter over Loops on ALL Chromosomes (", res, ")"),
       x = "Relative position to loop",
       y = "Counts"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

plot.promoter.hist

# Promoter Density Plot (by resolution combined, right)
plot.promoter.dens <- ggplot(relative.pos.df.promoter.dist.result, aes(x = value, color = loop.res, fill = loop.res)) +
  geom_density(alpha = 0.3) +  # transparent
  ylim(c(0, 1)) +
  scale_color_manual(values = c("5K" = "#619CFF", "10K" = "#F8766D", "25K" = "#00BA38")) + 
  scale_fill_manual(values = c("5K" = "#619CFF", "10K" = "#F8766D", "25K" = "#00BA38")) +
  labs(title = "Density of Promoter Found over Loops",
       x = "Relative position to loop",
       y = "Density",
       color = "Resolution",
       fill = "Resolution") + 
  theme(plot.title = element_text(hjust = 0.5))

plot.promoter.dens

return(plot.promoter.hist + plot.promoter.dens)


pdf("./figures/submission/overall_distribution_of_promoter_by_resolution.pdf", width = 11, height = 5)

final_promoter_plot <- plot.promoter.hist | plot.promoter.dens 
print(final_promoter_plot)

dev.off()

# plot_promoter_all_res <- create_promoter_plots(relative.pos.df.promoter.dist.result, "All")
# plot_promoter_5K <- create_promoter_plots(relative.pos.df.promoter.dist.result, "5K")
# plot_promoter_10K <- create_promoter_plots(relative.pos.df.promoter.dist.result, "10K")
# plot_promoter_25K <- create_promoter_plots(relative.pos.df.promoter.dist.result, "25K")
# pdf("./figures/0909/overall_distribution_promoter_on_all_chromosomes_1_distance.pdf", width = 11, height = 8.5) # pdf for a distance as a padding
# (plot_promoter_all_res/plot_promoter_5K)
# (plot_promoter_10K / plot_promoter_25K)
# dev.off()

########################
# 4. promoter
# 4-4. Distribution of promoter at each end in a loop: for the number of promoter used in filtering valid loops: figures
# so, the object should be used one without padding: df.DISTINCT.loop.deep.sample.all
# 4-4-1. inner distance analysis 1: checking inner distance between x3 and y0
########################
df.DISTINCT.loop.deep.sample.all.padded.for.promoter <- df.DISTINCT.loop.deep.sample.all %>% 
  mutate(x12 = (x1 + x2)/2, y12 = (y1 + y2)/2) %>% # middle point of each end
  mutate(distance = y12 - x12) %>% 
  mutate(x0 = ifelse(x12 - (distance / 2) < 0, 0, x12 - (distance / 2)), y3 = y12 + (distance/2)) %>% # 1/2 distance for OUTER PADDING
  mutate(x12 = (x12 + (distance / 4)), y12 = ifelse((y12 - (distance/4)) < 0, 0, y12 - (distance/4))) %>% # 1/4 distance for INNER PADDING
  mutate(new.loop.id = paste0(chr1, '_', x0, '_', x12, '_', chr2, '_', y12, '_', y3, '_', end.distance))

df.DISTINCT.loop.deep.sample.all.padded.for.promoter %>% head()

df.DISTINCT.loop.deep.sample.all.padded.for.promoter.UP.GR<- GRanges(seqnames=df.DISTINCT.loop.deep.sample.all.padded.for.promoter$chr1, 
                                                                     ranges=IRanges(start=df.DISTINCT.loop.deep.sample.all.padded.for.promoter$x0, 
                                                                                    end=df.DISTINCT.loop.deep.sample.all.padded.for.promoter$x12), 
                                                                     loop.id=df.DISTINCT.loop.deep.sample.all.padded.for.promoter$loop.id, 
                                                                     end.distance=df.DISTINCT.loop.deep.sample.all.padded.for.promoter$end.distance,
                                                                     resolution=df.DISTINCT.loop.deep.sample.all.padded.for.promoter$resolution,
                                                                     new.loop.id=df.DISTINCT.loop.deep.sample.all.padded.for.promoter$new.loop.id)

# 1. Promoter + UPSTREAM
index.distinct.promoter.w.up.loop.each.end <- findOverlaps(
  df.promoter.rn7.GR, 
  df.DISTINCT.loop.deep.sample.all.padded.for.promoter.UP.GR, 
  type = "within",
  select = "all"
)

loop.up.hits.promoter.each.end <- subjectHits(index.distinct.promoter.w.up.loop.each.end)
promoter.up.hits.each.end <- queryHits(index.distinct.promoter.w.up.loop.each.end)

df.overlapping.promoter.w.UPSTREAM.result.each.end <- data.frame(
  up.loop.id = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.promoter.UP.GR)$loop.id[loop.up.hits.promoter.each.end],
  end.up.distance = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.promoter.UP.GR)$end.distance[loop.up.hits.promoter.each.end],
  resolution = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.promoter.UP.GR)$resolution[loop.up.hits.promoter.each.end],
  promoter.id = mcols(df.promoter.rn7.GR)$promoter.id[promoter.up.hits.each.end],
  WHERE = "UP"
) %>% 
  mutate(promoter.loop.up.id = str_c(up.loop.id, '|', promoter.id, '|', WHERE))

df.overlapping.promoter.w.UPSTREAM.result.each.end %>% dim() # 69110/95831(TSS)
df.overlapping.promoter.w.UPSTREAM.result.each.end %>% head()

df.overlapping.promoter.w.UPSTREAM.result.each.end %>% 
  count(end.up.distance)

# end.up.distance     n
# 1            5000 15431
# 2           10000 20634
# 3           25000 33045
# total: 69110

# 2. Promoter + DOWNSTREAM
df.DISTINCT.loop.deep.sample.all.padded.for.promoter.DOWN.GR<- GRanges(seqnames=df.DISTINCT.loop.deep.sample.all.padded.for.promoter$chr1, 
                                                                       ranges=IRanges(start=df.DISTINCT.loop.deep.sample.all.padded.for.promoter$y12, 
                                                                                      end=df.DISTINCT.loop.deep.sample.all.padded.for.promoter$y3), 
                                                                       loop.id=df.DISTINCT.loop.deep.sample.all.padded.for.promoter$loop.id, 
                                                                       end.distance=df.DISTINCT.loop.deep.sample.all.padded.for.promoter$end.distance,
                                                                       resolution=df.DISTINCT.loop.deep.sample.all.padded.for.promoter$resolution,
                                                                       new.loop.id=df.DISTINCT.loop.deep.sample.all.padded.for.promoter$new.loop.id)

# 1. promoter + DOWNSTREAM (df.DISTINCT.loop.deep.sample.all.padded.for.promoter.DOWN.GR, df.promoter.rn7.GR)
index.distinct.promoter.w.down.loop.each.end <- findOverlaps(
  df.promoter.rn7.GR, 
  df.DISTINCT.loop.deep.sample.all.padded.for.promoter.DOWN.GR, 
  type = "within",
  select = "all"
)

end.loop.down.promoter.hits <- subjectHits(index.distinct.promoter.w.down.loop.each.end)
end.promoter.down.hits <- queryHits(index.distinct.promoter.w.down.loop.each.end)

df.overlapping.promoter.w.DOWNSTREAM.result.each.end <- tibble(
  down.loop.id = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.promoter.UP.GR)$loop.id[end.loop.down.promoter.hits],
  end.down.distance = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.promoter.UP.GR)$end.distance[end.loop.down.promoter.hits],
  resolution = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.promoter.UP.GR)$resolution[end.loop.down.promoter.hits],
  promoter.id = mcols(df.promoter.rn7.GR)$promoter.id[end.promoter.down.hits],
  WHERE = "DOWN"
) %>% 
  mutate(promoter.loop.down.id = str_c(down.loop.id, '|', promoter.id, '|', WHERE))

df.overlapping.promoter.w.DOWNSTREAM.result.each.end %>% dim() # 79241     6
df.overlapping.promoter.w.DOWNSTREAM.result.each.end %>% head()

df.overlapping.promoter.w.DOWNSTREAM.result.each.end %>% 
  count(end.down.distance)
# end.down.distance     n
# 1              5000 18525
# 2             10000 23067
# 3             25000 37649

# total: 79241 = 18525 + 23067 + 37649

df.overlapping.promoter.w.UPSTREAM.result.each.end %>% head()

########## bind_rows(UPSTREAM & DOWNSTREAM) -> BOTH
df.overlapping.promoter.w.BOTH.result <- bind_rows(df.overlapping.promoter.w.UPSTREAM.result.each.end %>% 
                                                     mutate(loop.id = up.loop.id, end.distance = end.up.distance) %>% 
                                                     mutate(case.id = promoter.loop.up.id) %>% 
                                                     dplyr::select(-c(up.loop.id, end.up.distance, promoter.loop.up.id)), 
                                                   df.overlapping.promoter.w.DOWNSTREAM.result.each.end %>% 
                                                     mutate(loop.id = down.loop.id, end.distance = end.down.distance) %>% 
                                                     mutate(case.id = promoter.loop.down.id) %>% 
                                                     dplyr::select(-c(down.loop.id, end.down.distance, promoter.loop.down.id))) %>% 
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

df.overlapping.promoter.w.BOTH.result %>% dim() # 148351
df.overlapping.promoter.w.BOTH.result %>% head()
df.overlapping.promoter.w.BOTH.result %>% count(resolution)

########################
# 4. promoter
# 2-4. Distribution of promoter at ends in a loop: for the number of promoter used in filtering valid loops: figures
# so, the object should be used one without padding: df.DISTINCT.loop.deep.sample.all
# 2-4-3. data processing for getting information of promoter at ends in loops 
########################
# for boxplot
df.overlapping.promoter.w.BOTH.result %>% head(2) # distance, resolution, ctcf.id, WHERE, loop.id, end.distance, case.id, chr

df.overlapping.promoter.w.BOTH.result.boxplot <- df.overlapping.promoter.w.BOTH.result %>% 
  group_by(chr, loop.id, WHERE, resolution) %>%
  # group_by(loop.id, WHERE, resolution) %>%
  # group_by(chr, loop.id, resolution) %>%
  summarise(promoter_count_by_loop_id_each_end = n_distinct(promoter.id), .groups = 'drop')

df.overlapping.promoter.w.BOTH.result.boxplot

# for Q1
df.overlapping.promoter.w.BOTH.result %>% head()

df.promoter.counts.each.end <- df.overlapping.promoter.w.BOTH.result %>%
  group_by(loop.id, WHERE, resolution) %>%
  summarise(promoter_count_each_end = n_distinct(promoter.id), .groups = 'drop')

df.promoter.counts.each.end

promoter.stats.by.resolution.each.end <- df.promoter.counts.each.end %>%
  # group_by(resolution) %>%
  group_by(WHERE, resolution) %>%
  summarise(
    Min = min(promoter_count_each_end),
    Q1 = quantile(promoter_count_each_end, 0.25, na.rm = TRUE),
    Median = median(promoter_count_each_end, na.rm = TRUE),
    Q3 = quantile(promoter_count_each_end, 0.75, na.rm = TRUE),
    Mean = mean(promoter_count_each_end, na.rm = TRUE),
    SD = sd(promoter_count_each_end, na.rm = TRUE),
    Max = max(promoter_count_each_end),
    .groups = 'drop'
  )

promoter.stats.by.resolution.each.end

# resolution   Min    Q1 Median    Q3  Mean    SD   Max
# 1 5K             1     1      2     3  4.46  25.9   746
# 2 10K            1     1      2     3  3.07  13.6   692
# 3 25K            1     1      2     3  4.34  19.9   746

# WHERE resolution   Min    Q1 Median    Q3  Mean    SD   Max
# 1 UP    5K             1     1      2     3  3.98 20.2    592
# 2 UP    10K            1     1      2     3  2.86  9.09   259
# 3 UP    25K            1     1      2     3  4.03 15.6    592
# 4 DOWN  5K             1     1      1     2  4.96 30.8    746
# 5 DOWN  10K            1     1      2     3  3.28 17.1    692
# 6 DOWN  25K            1     1      2     3  4.65 23.5    746

# quantile for EACH END// NOT in a loop : results are loops have ctcfs more than 18 at least in each end
q3_promoter_count_each_end <- quantile(df.overlapping.promoter.w.BOTH.result.boxplot$promoter_count_by_loop_id_each_end, 0.75, na.rm = TRUE)
q3_promoter_count_each_end # 3 :: integrity PASS

q1_promoter_count_each_end <- quantile(df.overlapping.promoter.w.BOTH.result.boxplot$promoter_count_by_loop_id_each_end, 0.25, na.rm = TRUE)
q1_promoter_count_each_end # 1 :: integrity PASS

# drawing boxplot
df.overlapping.promoter.w.BOTH.result.boxplot
df.overlapping.promoter.w.BOTH.result.boxplot %>% head(4)

# figure by chr & res
boxplot.w.promoter.by.chr.and.res.each.end <- df.overlapping.promoter.w.BOTH.result.boxplot %>% 
  ggplot(aes(x = chr, y = promoter_count_by_loop_id_each_end, fill = resolution)) +
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
  scale_y_continuous(limits = c(0, 2 * q3_promoter_count_each_end), 
                     breaks = seq(0, 2 * q3_promoter_count_each_end, by = 1)) +
  labs(title = "Boxplot for Number of Promoter by Chromosome and Resolution", x = "Chromosomes", y = "Number of Promoter in a loop")

boxplot.w.promoter.by.chr.and.res.each.end

pdf("./figures/submission/end_boxplot_num_promoter_by_chr_res.pdf", width = 16.5, height = 23.5)
grid.arrange(boxplot.w.promoter.by.chr.and.res.each.end, ncol = 1)
dev.off()

boxplot.w.promoter.by.res <- df.overlapping.promoter.w.BOTH.result.boxplot %>% 
  ggplot(aes(x = resolution, y = promoter_count_by_loop_id_each_end, fill = resolution)) +
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
  scale_y_continuous(limits = c(0, 2 * q3_promoter_count_each_end), 
                     breaks = seq(0, 2 * q3_promoter_count_each_end, by = 10)) +
  labs(title = "Boxplot for Number of Promote by Resolution", x = "Resolution", y = "Number of Promoter in a loop")

boxplot.w.promoter.by.res

pdf("./figures/submission/end_boxplot_num_promoter_by_res.pdf", width = 16.5, height = 23.5)
grid.arrange(boxplot.w.promoter.by.res, ncol = 1)
dev.off()

# figures above in a PDF
pdf("./figures/submission/boxplot_of_promoter_count.pdf", width = 16.5, height = 23.5)
grid.arrange(boxplot.w.promoter.by.chr.and.res.each.end, boxplot.w.promoter.by.res, ncol = 1)
dev.off()
promoter.stats.by.resolution.each.end %>% head()

########################
# 4. promoter
# 4-4. Distribution of promoter at ends in a loop: for the number of promoter used in filtering valid loops: figures
# so, the object should be used one without padding: df.DISTINCT.loop.deep.sample.all
# 4-4-4. Distribution of promoter at both ends in a loop: figure
########################
# figure by end
boxplot.w.promoter.ALL.chr <- ggplot(df.overlapping.promoter.w.BOTH.result.boxplot, aes(x = WHERE, y = promoter_count_by_loop_id_each_end, fill = WHERE)) +
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
  scale_y_continuous(limits = c(0, 2 * q3_promoter_count_each_end),  # y축을 Q3의 두 배로 제한
                     breaks = seq(0, 2 * q3_promoter_count_each_end, by = 5)) +
  labs(title = "Boxplot by UP/DOWNSTREAM End", x = "Upstream and Downstream End", y = "Number of CTCF inside each ends in a loop")

boxplot.w.promoter.ALL.chr

# figure by end & res
boxplot.w.promoter.ALL.chr.by.res <- ggplot(df.overlapping.promoter.w.BOTH.result.boxplot, aes(x = WHERE, y = promoter_count_by_loop_id_each_end, fill = resolution)) +
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
  scale_y_continuous(limits = c(0, 2 * q3_promoter_count_each_end), 
                     breaks = seq(0, 2 * q3_promoter_count_each_end, by = 5)) +
  labs(title = "Boxplot for Number of Promoter in UP/DOWNSTREAM End by Resolution", x = "Upstream and Downstream End", y = "Number of Promoter", fill = "Resolution")

boxplot.w.promoter.ALL.chr.by.res

pdf("./figures/submission/boxplot_for_no.promoter_by_end_res.pdf", width = 16.5, height = 23.5)
grid.arrange(boxplot.w.promoter.ALL.chr, boxplot.w.promoter.ALL.chr.by.res, ncol = 1)
dev.off()

###############################################
# (1/2) Histogram of promoter Counts per Loop (Exclusive)
###############################################

df.overlapping.promoter.w.BOTH.result.boxplot %>% head()

promoter_count_histogram_data <- df.overlapping.promoter.w.BOTH.result.boxplot %>%
  group_by(loop.id) %>%
  summarise(promoter_count_total = sum(promoter_count_by_loop_id_each_end), .groups = 'drop')

total_distinct_loops_promoter <- n_distinct(promoter_count_histogram_data$loop.id)

more_than_5_loops <- promoter_count_histogram_data %>%
  filter(promoter_count_total > 5) %>%
  nrow()

more_than_5_loops_percent <- round(more_than_5_loops / total_distinct_loops_promoter * 100, 1)

promoter.count.per.loop.histogram <- ggplot(promoter_count_histogram_data, aes(x = promoter_count_total)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  scale_x_continuous(limits = c(0, 6.5), breaks = seq(0, 6, by = 1)) +
  labs(title = "Histogram of Promoter Counts per Loop (Exclusive)", 
       x = "Total Promoter Count per Loop", 
       y = "Number of Loops") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = rel(1.2))) +
  stat_count(aes(label = paste0(..count.., " (", round(..count../total_distinct_loops_promoter * 100, 1), "%)")), 
             geom = "text", 
             vjust = -0.5, size = 3) +
  annotate("text", x = 5, y = max(table(promoter_count_histogram_data$promoter_count_total)) * 1.1, 
           label = paste("Total Loops:", total_distinct_loops_promoter), size = 4, hjust = 1) +
  annotate("text", x = 5, y = max(table(promoter_count_histogram_data$promoter_count_total)) * 1.05, 
           label = paste("Loops > 5 Promoter: ", more_than_5_loops, " (", more_than_5_loops_percent, "%)"), 
           size = 3, hjust = 1, color = "red")

promoter.count.per.loop.histogram

df.overlapping.promoter.w.BOTH.result.boxplot %>% head()

###############################################
# (2/2) Histogram of Loops with promoter Concentrated in UP or DOWN
###############################################
promoter.count.per.loop <- df.overlapping.promoter.w.BOTH.result.boxplot %>%
  group_by(loop.id, WHERE) %>%
  summarise(total_promoter = sum(promoter_count_by_loop_id_each_end), .groups = 'drop')

df.overlapping.promoter.w.BOTH.result.boxplot %>% head()
promoter.count.per.loop %>% head()

# loops which have only in an end
only_up_or_down_loops_promoter <- promoter.count.per.loop %>%
  group_by(loop.id) %>%
  filter(n() == 1) %>%
  ungroup()
only_up_or_down_loops_promoter %>% head()

# counting the loops above
histogram.data <- only_up_or_down_loops_promoter %>%
  group_by(total_promoter) %>%
  summarise(one_sided_promoter_loop_count = n(), .groups = 'drop')
# histogram.data %>% view()

promoter.sum.per.loop <- promoter.count.per.loop %>%
  group_by(loop.id) %>%
  summarise(total_promoter_sum = sum(total_promoter), .groups = 'drop')
promoter.sum.per.loop %>% head()

promoter.count.by.sum <- promoter.sum.per.loop %>%
  group_by(total_promoter_sum) %>%
  summarise(loop_count = n(), .groups = 'drop')
histogram.data %>% head()
promoter.count.by.sum %>% head()

one.sided.loop.histogram.data.promoter.final <- right_join(histogram.data, promoter.count.by.sum, by = c("total_promoter" = "total_promoter_sum")) %>%
  mutate(ratio = one_sided_promoter_loop_count / loop_count)
one.sided.loop.histogram.data.promoter.final %>% head()

one.sided.loop.histogram <- ggplot(one.sided.loop.histogram.data.promoter.final, aes(x = total_promoter, y = one_sided_promoter_loop_count)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  scale_x_continuous(limits = c(0, 6.5), breaks = seq(0, 6, by = 1)) +
  labs(title = "Histogram of Loops with Promoter Concentrated in UP or DOWN",
       x = "Total Promoter Count per Loop",
       y = "Number of Loops") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = paste0(one_sided_promoter_loop_count, "/", loop_count, " (", round(ratio * 100, 1), "%)")),
            vjust = -0.5, size = 3, color = "black")

one.sided.loop.histogram
one.sided.loop.histogram.data.promoter.final %>% view()

pdf("./figures/0909/promoter_count_per_loop_histogram.pdf", width = 16.5, height = 23.5)

grid.arrange(promoter.count.per.loop.histogram, one.sided.loop.histogram, ncol = 1)
dev.off()
