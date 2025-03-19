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
<<<<<<< HEAD
# setwd('~/Desktop/temp/enhancer/dropbox_enhancer_doosan')
# setwd('/home/pkim/dropbox/Gateway_to_Hao/enhancer/r_files')
# source(file.path('/home/pkim/dropbox/Gateway_to_Hao/project_common_code/', 'variables.R'))
# source(file.path('/home/pkim/dropbox/Gateway_to_Hao/project_common_code/', 'funcs.R'))
=======
setwd('~/Desktop/temp/enhancer/dropbox_enhancer_doosan')
setwd('/home/pkim/dropbox/Gateway_to_Hao/enhancer/r_files')
source(file.path('/home/pkim/dropbox/Gateway_to_Hao/project_common_code/', 'variables.R'))
source(file.path('/home/pkim/dropbox/Gateway_to_Hao/project_common_code/', 'funcs.R'))
>>>>>>> d754844 (changing variables, generating figures for submission, splitting file of enhancer_synced.R into three.)

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
# 2. CTCF
# 2-1. Exploratory Data analysis (EDA) : CTCF
########################
# Linux
df.init.ctcf<-read.table(file="/home/pkim/Desktop/temp/enhancer/dropbox_enhancer_doosan/data/ctcf/fimo_2nd_trial_meme/fimo_2nd_trial.txt", header=TRUE, sep="\t")
df.init.ctcf<-read.table(file="/home/pkim/dropbox/Gateway_to_Hao/enhancer/data/ctcf/fimo_2nd_trial_meme/fimo_2nd_trial.txt", header=TRUE, sep="\t")
df.init.ctcf %>% count() # 5767921

# Mac
# ctcf<-read.table(file="/home/pkim/Desktop/temp/enhancer/dropbox_enhancer_doosan/data/ctcf/fimo_2nd_trial_meme/fimo_2nd_trial.txt", sep="\t", col.names = c("chr", "start", "end", "strand", "length"), header = FALSE) %>% 
df.init.ctcf<-read.table(file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/data/ctcf/fimo_2nd_trial_meme/fimo_2nd_trial.txt", header=TRUE, sep="\t")
df.init.ctcf %>% dim() # 5767921
df.init.ctcf %>% head()

# strand checking
df.init.ctcf %>%                              # +: 2891072, -: 2876849 = 5767921
  count(strand)
# removing dups including strand
df.init.ctcf %>% 
  distinct(chr, start, end, strand)           # 2701585/5767921
# removing dups without strand
df.init.ctcf %>% 
  distinct(chr, start, end)                   # 2544216/5767921
# removing dups with all including length
df.init.ctcf %>% 
  distinct()                                  # 2701585/5767921
# verifying for length column
df.init.ctcf %>% 
  mutate(test = ifelse((end - start) == length, TRUE, FALSE)) %>% 
  count(test)

########################
# 2. CTCF
# 2-2. CTCF data preprocessing: id and dedup GRange Obj. (df.DISTINCT.fimo.2nd.trial.ctcf/ df.DISTINCT.ctcf.2nd.fimo.GR)
########################
# generating id column (long running time)
df.DISTINCT.fimo.2nd.trial.ctcf <- df.init.ctcf %>%
  distinct() %>% 
  mutate(id = str_c(chr, "_", start, "_", end, "_", strand, '_', length)) %>% 
  mutate(start = as.numeric(start)) %>% 
  mutate(end = as.numeric(end)) 

df.DISTINCT.fimo.2nd.trial.ctcf # 2701585/5767921 : 0.4683811
df.DISTINCT.fimo.2nd.trial.ctcf %>% head()  # id column done

# GRanges Obj.: df.DISTINCT.ctcf.2nd.fimo.GR
df.DISTINCT.ctcf.2nd.fimo.GR <- GRanges(
  seqnames = as.character(df.DISTINCT.fimo.2nd.trial.ctcf$chr),
  ranges = IRanges(start = df.DISTINCT.fimo.2nd.trial.ctcf$start, 
                   end = df.DISTINCT.fimo.2nd.trial.ctcf$end),
  strand = df.DISTINCT.fimo.2nd.trial.ctcf$strand  # strand 
)

# metadata: id
mcols(df.DISTINCT.ctcf.2nd.fimo.GR)$id <- df.DISTINCT.fimo.2nd.trial.ctcf$id

df.DISTINCT.ctcf.2nd.fimo.GR

# data file for Phenogen
df.ctcf.for.phenogen <- df.DISTINCT.fimo.2nd.trial.ctcf %>% 
  dplyr::select(chr, start) %>% 
  mutate(chr = str_remove(chr, "^chr")) %>% 
  mutate(phenotype = "CTCF")

df.ctcf.for.phenogen

write_tsv(df.ctcf.for.phenogen, "~/dropbox/Gateway_to_Hao/publication/PE-interaction/ideogram/phenogen_input_ctcf.tsv")
########################
# 2. CTCF
# 2-3. overall distribution of CTCF on loops
########################

# adding 1/2 distance in each end (total distance becomes 2*distance)
# overall.df.DISTINCT.loop.deep.sample.all <- overall.df.DISTINCT.loop.deep.sample.all.5.distance
# adding 1 distance in each end (total distance becomes 3*distance)
overall.df.DISTINCT.loop.deep.sample.all <- overall.df.DISTINCT.loop.deep.sample.all.1.distance

overall.df.DISTINCT.loop.deep.sample.all.GR <- GRanges(
  seqnames = as.character(overall.df.DISTINCT.loop.deep.sample.all$chr1),  
  ranges = IRanges(
    start = overall.df.DISTINCT.loop.deep.sample.all$x0, 
    end = overall.df.DISTINCT.loop.deep.sample.all$y3
  )
)

# metadata: id
mcols(overall.df.DISTINCT.loop.deep.sample.all.GR) <- data.frame(
  id = overall.df.DISTINCT.loop.deep.sample.all$loop.id,
  new.loop.id = overall.df.DISTINCT.loop.deep.sample.all$new.loop.id,
  mcols(overall.df.DISTINCT.loop.deep.sample.all.GR)$resolution <- overall.df.DISTINCT.loop.deep.sample.all$resolution
)

index.distinct.ctcf.w.overall.whole.loop <- findOverlaps(
  df.DISTINCT.ctcf.2nd.fimo.GR, 
  overall.df.DISTINCT.loop.deep.sample.all.GR, 
  type = "within",
  select = "all"
)

overall.loop.for.ctcf.hits <- subjectHits(index.distinct.ctcf.w.overall.whole.loop)
overall.ctcf.on.loop.hits <- queryHits(index.distinct.ctcf.w.overall.whole.loop)

df.ctcf.dist.result <- tibble(
  loop.id = overall.df.DISTINCT.loop.deep.sample.all$new.loop.id[overall.loop.for.ctcf.hits],
  loop.start = overall.df.DISTINCT.loop.deep.sample.all$x0[overall.loop.for.ctcf.hits],
  loop.end = overall.df.DISTINCT.loop.deep.sample.all$y3[overall.loop.for.ctcf.hits],
  loop.res = overall.df.DISTINCT.loop.deep.sample.all$resolution[overall.loop.for.ctcf.hits],
  # loop.new.distance = overall.df.DISTINCT.loop.deep.sample.all$new_distance[overall.loop.for.ctcf.hits],
  ctcf.id = df.DISTINCT.fimo.2nd.trial.ctcf$id[overall.ctcf.on.loop.hits],
  ctcf.start = df.DISTINCT.fimo.2nd.trial.ctcf$start[overall.ctcf.on.loop.hits],
  ctcf.end = df.DISTINCT.fimo.2nd.trial.ctcf$end[overall.ctcf.on.loop.hits]
) %>% 
  mutate(across(where(is.numeric), ~ format(., scientific = FALSE)))

df.ctcf.dist.result %>% head()
df.ctcf.dist.result %>% dim() # 40861258/57612728(0.5 distance/1 distance)

relative.pos.df.ctcf.dist.result <- df.ctcf.dist.result %>%
  mutate(ctcf.start = as.numeric(ctcf.start),
         ctcf.end = as.numeric(ctcf.end),
         loop.start = as.numeric(loop.start),
         loop.end = as.numeric(loop.end),
         pos_coord = round((ctcf.start + ctcf.end) / 2),
         loop_length = (loop.end - loop.start),
         relative_pos = pos_coord - loop.start,
         value = relative_pos / (loop_length / 3) - 1
  ) %>%
  # mutate(value = (relative_pos / (loop_length / 2)) - 1) %>% # for padding .5x distance
  dplyr::select(loop.id, ctcf.id, value, loop.res)

relative.pos.df.ctcf.dist.result %>% 
  head()

########################
# 2. CTCF
# 2-3. overall distribution of CTCF on loops: figures
# 2-3-1. by CHROMOSOME
########################

chromosomes <- c(1:20, "X", "Y")

plot_ctcf_list <- list()

for (chr in chromosomes) {
  tryCatch({
    
    message("START: Processing chromosome: ", chr)
    
    relative.pos.df.ctcf.dist.result.chr <- relative.pos.df.ctcf.dist.result %>% 
      filter(str_detect(loop.id, paste0("_chr", chr, "_")))
    
    plot1.ctcf.dens <- relative.pos.df.ctcf.dist.result.chr %>% 
      ggplot(aes(x = value)) +
      geom_density(fill = "skyblue", color = "black", alpha = 0.7) +
      ylim(c(0,1))+
      labs(title = paste0("Density of CTCF Found over Loops on Chr", chr),
           x = "Relative position to loop",
           y = "Density"
      )
    
    plot1.ctcf.hist <- relative.pos.df.ctcf.dist.result.chr %>% 
      ggplot(aes(x = value)) +
      geom_histogram(fill = "skyblue", color = "black", alpha = 0.7, bins=200) +
      labs(title = paste0("Histogram of CTCF Found over Loops on Chr", chr),
           x = "Relative position to loop",
           y = "Counts"
      )
    
    combined_plot_ctcf <- plot1.ctcf.hist + plot1.ctcf.dens # + can be used instead of |
    plot_ctcf_list[[chr]] <- combined_plot_ctcf
    
    message("END: Processing chromosome: ", chr)
  }, error = function(e) {
    
    message("Error processing chromosome: ", chr)
    message("Error message: ", e$message)
  })
}

#pdf("overall_distribution_of_CTCF_on_each_chr_.5_distance.pdf", width = 11, height = 8.5)
pdf("./figures/submission/overall_distribution_of_CTCF_by_chrmosome.pdf", width = 11, height = 8.5)

num_plots <- length(plot_ctcf_list) # 22 chromosomes
plots_per_page <- 4

for (i in seq(1, num_plots, by = plots_per_page)) {
  end_idx <- min(i + plots_per_page - 1, num_plots)
  page_plots <- plot_ctcf_list[i:end_idx]
  
  combined_page <- plot_grid(plotlist = page_plots, ncol = 1, nrow = 4)
  
  print(combined_page)
}

dev.off()

relative.pos.df.ctcf.dist.result %>% head()

########################
# 2. CTCF
# 2-3. overall distribution of CTCF on loops: figures
# 2-3-2. by resolution
########################
# CTCF Histogram (all resolution combined, left)
plot.ctcf.hist <- ggplot(relative.pos.df.ctcf.dist.result, aes(x = value)) +
  geom_histogram(fill = "skyblue", color = "black", alpha = 0.7, bins = 200) +
  labs(title = "Histogram of CTCF Found over Loops",
       x = "Relative position to loop",
       y = "Counts") + 
  theme(plot.title = element_text(hjust = 0.5))

# CTCF Density Plot (by resolution combined, right)
plot.ctcf.dens <- ggplot(relative.pos.df.ctcf.dist.result, aes(x = value, color = loop.res, fill = loop.res)) +
  geom_density(alpha = 0.3) +  # transparent
  ylim(c(0, 1)) +
  scale_color_manual(values = c("5K" = "#619CFF", "10K" = "#F8766D", "25K" = "#00BA38")) + 
  scale_fill_manual(values = c("5K" = "#619CFF", "10K" = "#F8766D", "25K" = "#00BA38")) +
  labs(title = "Density of CTCF Found over Loops",
       x = "Relative position to loop",
       y = "Density",
       color = "Resolution",
       fill = "Resolution") + 
  theme(plot.title = element_text(hjust = 0.5))

pdf("figures/submission/overall_distribution_of_CTCF_by_resolution.pdf", width = 11, height = 5)

final_ctcf_plot <- plot.ctcf.hist | plot.ctcf.dens 
print(final_ctcf_plot)

dev.off()

########################
# 2. CTCF
# 2-4. Distribution of CTCF at each end in a loop: for the number of CTCF used in filtering valid loops: figures
# so, the object should be used one without padding: df.DISTINCT.loop.deep.sample.all
# 2-4-1. inner distance analysis 1: checking inner distance between x3 and y0
########################

df.DISTINCT.loop.deep.sample.all %>% head(2)

apply_padding <- function(df, padding_factor) {
  padded_distance <- df$end.distance * padding_factor
  
  df_padded <- df %>%
    mutate(x0 = ifelse(x1 - padded_distance < 0, 0, x1 - padded_distance)) %>%
    mutate(x3 = ifelse(x2 + padded_distance < 0, 0, x2 + padded_distance)) %>%
    mutate(y0 = y1 - padded_distance) %>%
    mutate(y3 = y2 + padded_distance) %>%
    mutate(inner.distance = y0 - x3) %>%
    mutate(inner.distance = as.numeric(format(inner.distance, scientific = FALSE)))
  
  return(df_padded)
}

# IQR based filtering func
generate_boxplot <- function(data, title) {
  q1 <- quantile(data$inner.distance, 0.25)
  q3 <- quantile(data$inner.distance, 0.75)
  median_value <- median(data$inner.distance)
  min_value <- min(data$inner.distance)
  max_value <- max(data$inner.distance)
  iqr <- q3 - q1
  
  filtered_inner_distance <- data$inner.distance[
    data$inner.distance >= (q1 - 1.5 * iqr) &
      data$inner.distance <= (q3 + 1.5 * iqr)
  ]
  
  boxplot(filtered_inner_distance, 
          main = title, 
          ylab = "Inner Distance", 
          col = "lightblue", 
          border = "black",
          frame = FALSE)
  
  text(1.3, min(filtered_inner_distance), paste("Min:", round(min_value)), pos = 4, col = "black")
  text(1.3, max(filtered_inner_distance), paste("Max:", round(max_value)), pos = 4, col = "black")
  text(0.7, q1, paste("Q1:", round(q1)), pos = 2, col = "red")
  text(0.7, q3, paste("Q3:", round(q3)), pos = 2, col = "red")
  text(1.3, median_value, paste("Median:", round(median_value)), pos = 4, col = "black")
}

# 50% padding
df_50_padded <- apply_padding(df.DISTINCT.loop.deep.sample.all, padding_factor = 0.5)
# 100% padding
df_100_padded <- apply_padding(df.DISTINCT.loop.deep.sample.all, padding_factor = 1.0)

pdf("./figures/submission/inner_distance_distribution_boxplots_after_padding.pdf", width = 16, height = 8)
par(mfrow = c(1, 2))  

# 100% padding boxplot
generate_boxplot(df_100_padded, "Inner Distance Distribution after Loop Padding")

dev.off()

########################
# 2. CTCF
# 2-4. Distribution of CTCF at ends in a loop: for the number of CTCF used in filtering valid loops: figures
# so, the object should be used one without padding: df.DISTINCT.loop.deep.sample.all
# 2-4-3. data processing for getting information of CTCF at ends in loops 
########################

df_50_padded %>% head(2)
df_100_padded %>% head(2)
df.DISTINCT.loop.deep.sample.all %>% head(2)

df.DISTINCT.loop.deep.sample.all
df.DISTINCT.loop.deep.sample.all.up.GR

# 1. CTCF + UPSTREAM (df.DISTINCT.loop.deep.sample.all.up, df.DISTINCT.fimo.2nd.trial.ctcf)
index.distinct.ctcf.w.up.loop <- findOverlaps(df.DISTINCT.ctcf.2nd.fimo.GR, 
                                              df.DISTINCT.loop.deep.sample.all.up.GR, 
                                              type = "within",
                                              select = "all")

end.loop.up.ctcf.hits <- subjectHits(index.distinct.ctcf.w.up.loop)
end.ctcf.up.hits <- queryHits(index.distinct.ctcf.w.up.loop)

df.overlapping.CTCF.w.UPSTREAM.result <- tibble(
  up.loop.id = mcols(df.DISTINCT.loop.deep.sample.all.up.GR)$id[end.loop.up.ctcf.hits],
  end.up.distance = mcols(df.DISTINCT.loop.deep.sample.all.up.GR)$end.distance[end.loop.up.ctcf.hits],
  distance = mcols(df.DISTINCT.loop.deep.sample.all.up.GR)$distance[end.loop.up.ctcf.hits],
  resolution = mcols(df.DISTINCT.loop.deep.sample.all.up.GR)$resolution[end.loop.up.ctcf.hits],
  ctcf.id = mcols(df.DISTINCT.ctcf.2nd.fimo.GR)$id[end.ctcf.up.hits],
  WHERE = "UP"
) %>% 
  unite(ctcf.loop.up.id, up.loop.id, distance, ctcf.id, WHERE, sep = "|", remove = FALSE)

df.overlapping.CTCF.w.UPSTREAM.result %>% dim() # 1206765
df.overlapping.CTCF.w.UPSTREAM.result %>% head()

df.overlapping.CTCF.w.UPSTREAM.result %>% 
  count(end.up.distance)

# end.up.distance      n
# 1            5000 154950
# 2           10000 383562
# 3           25000 668253
# total: 1206765

# 2. CTCF + DOWNSTREAM (df.DISTINCT.loop.deep.sample.all.down, df.DISTINCT.fimo.2nd.trial.ctcf)
index.distinct.ctcf.w.down.loop <- findOverlaps(df.DISTINCT.ctcf.2nd.fimo.GR, 
                                                df.DISTINCT.loop.deep.sample.all.down.GR, 
                                                type = "within",
                                                select = "all")

end.loop.down.ctcf.hits <- subjectHits(index.distinct.ctcf.w.down.loop)
end.ctcf.down.hits <- queryHits(index.distinct.ctcf.w.down.loop)

df.overlapping.CTCF.w.DOWNSTREAM.result <- tibble(
  down.loop.id = mcols(df.DISTINCT.loop.deep.sample.all.down.GR)$id[end.loop.down.ctcf.hits],
  end.down.distance = mcols(df.DISTINCT.loop.deep.sample.all.down.GR)$end.distance[end.loop.down.ctcf.hits],
  distance = mcols(df.DISTINCT.loop.deep.sample.all.down.GR)$distance[end.loop.down.ctcf.hits],
  resolution = mcols(df.DISTINCT.loop.deep.sample.all.down.GR)$resolution[end.loop.down.ctcf.hits],
  ctcf.id = mcols(df.DISTINCT.ctcf.2nd.fimo.GR)$id[end.ctcf.down.hits],
  WHERE = "DOWN"
) %>% 
  unite(ctcf.loop.down.id, down.loop.id, distance, ctcf.id, WHERE, sep = "|", remove = FALSE)

df.overlapping.CTCF.w.DOWNSTREAM.result %>% dim() # 1213237
df.overlapping.CTCF.w.DOWNSTREAM.result %>% 
  count(end.down.distance)
# #end.down.distance      n
# 1              5000 158453
# 2             10000 390476
# 3             25000 664308
# total: 1213237

########## bind_rows(UPSTREAM & DOWNSTREAM) -> BOTH
df.overlapping.CTCF.w.BOTH.result <- bind_rows(df.overlapping.CTCF.w.UPSTREAM.result %>% 
                                                 mutate(loop.id = up.loop.id, end.distance = end.up.distance) %>% 
                                                 mutate(case.id = ctcf.loop.up.id) %>% 
                                                 dplyr::select(-c(up.loop.id, end.up.distance, ctcf.loop.up.id)), 
                                               df.overlapping.CTCF.w.DOWNSTREAM.result %>% 
                                                 mutate(loop.id = down.loop.id, end.distance = end.down.distance) %>% 
                                                 mutate(case.id = ctcf.loop.down.id) %>% 
                                                 dplyr::select(-c(down.loop.id, end.down.distance, ctcf.loop.down.id))) %>% 
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

df.overlapping.CTCF.w.BOTH.result %>% dim() # 2420002
df.overlapping.CTCF.w.BOTH.result %>% head()
df.overlapping.CTCF.w.BOTH.result %>% count(resolution)

# resolution       n
# 1         5K 154950 + 158453 =  313403: integrity PASS
# 2        10K 383562 + 390476 =  774038: integrity PASS
# 3        25K 668253 + 664308 = 1332561 : integrity PASS

# for boxplot
df.overlapping.CTCF.w.BOTH.result %>% head(2) # distance, resolution, ctcf.id, WHERE, loop.id, end.distance, case.id, chr

# for Q1
df.ctcf.counts <- df.overlapping.CTCF.w.BOTH.result %>%
  group_by(loop.id, WHERE, resolution) %>%
  summarise(ctcf_count = n_distinct(ctcf.id), .groups = 'drop')

ctcf.stats.by.resolution <- df.ctcf.counts %>%
  # group_by(resolution) %>%
  group_by(WHERE, resolution) %>%
  summarise(
    Min = min(ctcf_count),
    Q1 = quantile(ctcf_count, 0.25, na.rm = TRUE),
    Median = median(ctcf_count, na.rm = TRUE),
    Q3 = quantile(ctcf_count, 0.75, na.rm = TRUE),
    Mean = mean(ctcf_count, na.rm = TRUE),
    SD = sd(ctcf_count, na.rm = TRUE),
    Max = max(ctcf_count),
    .groups = 'drop'
  )

ctcf.stats.by.resolution

# resolution    Q1 Median    Q3  Mean    SD
# <fct>      <dbl>  <dbl> <dbl> <dbl> <dbl>
# 1 5K             7     26    35  25.6  20.1
# 2 10K           12     30    46  32.7  25.4
# 3 25K           25     43    71  51.7  37.4

# WHERE resolution    Q1 Median    Q3  Mean    SD
# <fct> <fct>      <dbl>  <dbl> <dbl> <dbl> <dbl>
# 1 UP    5K             7     26    34  25.3  19.8
# 2 UP    10K           12     30    45  32.4  25.4
# 3 UP    25K           25     44    72  51.8  37.8
# 4 DOWN  5K             8     26    35  26.0  20.3
# 5 DOWN  10K           12     30    47  32.9  25.4
# 6 DOWN  25K           25     43    71  51.5  37.1

# drawing boxplot
df.overlapping.CTCF.w.BOTH.result.boxplot <- df.overlapping.CTCF.w.BOTH.result %>% 
  group_by(chr, loop.id, WHERE, resolution) %>%
  # group_by(loop.id, WHERE, resolution) %>%
  # group_by(chr, loop.id, resolution) %>%
  summarise(ctcf_count_by_loop_id = n_distinct(ctcf.id), .groups = 'drop')

# quantile for EACH END// NOT in a loop : results are loops have ctcfs more than 18 at least in each end
q3_ctcf_count <- quantile(df.overlapping.CTCF.w.BOTH.result.boxplot$ctcf_count_by_loop_id, 0.75, na.rm = TRUE)
q3_ctcf_count # 54 :: integrity PASS

q1_ctcf_count <- quantile(df.overlapping.CTCF.w.BOTH.result.boxplot$ctcf_count_by_loop_id, 0.25, na.rm = TRUE)
q1_ctcf_count # 16 :: integrity PASS

########################
# 2. CTCF
# 2-4. Distribution of CTCF at ends in a loop: for the number of CTCF used in filtering valid loops: figures
# so, the object should be used one without padding: df.DISTINCT.loop.deep.sample.all
# 2-4-4. Distribution of CTCF at both ends in a loop: figure
########################
# figure by chr & res
boxplot.w.CTCF.by.chr.and.res <- df.overlapping.CTCF.w.BOTH.result.boxplot %>% 
  ggplot(aes(x = chr, y = ctcf_count_by_loop_id, fill = resolution)) +
  geom_boxplot() +
  stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 1)), position = position_dodge(width = 0.75), vjust = -0.5, size = 2) + 
  stat_summary(fun.data = function(y) {
    data.frame(y = median(y), label = paste0("Q1: ", round(quantile(y, 0.25), 1)))
  }, geom = "text", aes(label = after_stat(label)), position = position_dodge(width = 0.75), vjust = 6.5, size = 1.5) +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("5K" = "skyblue", "10K" = "lightcoral", "25K" = "springgreen"), breaks = c("5K", "10K", "25K")) +
  scale_y_continuous(breaks = seq(min(df.overlapping.CTCF.w.BOTH.result.boxplot$ctcf_count_by_loop_id), max(df.overlapping.CTCF.w.BOTH.result.boxplot$ctcf_count_by_loop_id), by = 10)) +
  labs(title = "Distribution for Number of CTCF by Chromosome and Resolution", x = "Chromosomes", y = "Number of CTCF in A Loop")

boxplot.w.CTCF.by.chr.and.res

pdf("figures/submission/boxplot_of_ctcf_count_by_chr_res.pdf", width = 16.5, height = 23.5)
grid.arrange(boxplot.w.CTCF.by.chr.and.res, ncol = 1)
dev.off()

boxplot.w.CTCF.by.res <- df.overlapping.CTCF.w.BOTH.result.boxplot %>% 
  ggplot(aes(x = resolution, y = ctcf_count_by_loop_id, fill = resolution)) +
  geom_boxplot() +
  stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 1)), position = position_dodge(width = 0.75), vjust = -0.5, size = 2) + 
  stat_summary(fun.data = function(y) {
    data.frame(y = median(y), label = paste0("Q1: ", round(quantile(y, 0.25), 1)))
  }, geom = "text", aes(label = after_stat(label)), position = position_dodge(width = 0.75), vjust = 6.5, size = 1.5) +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("5K" = "skyblue", "10K" = "lightcoral", "25K" = "springgreen"), breaks = c("5K", "10K", "25K")) +
  scale_y_continuous(breaks = seq(min(df.overlapping.CTCF.w.BOTH.result.boxplot$ctcf_count_by_loop_id), max(df.overlapping.CTCF.w.BOTH.result.boxplot$ctcf_count_by_loop_id), by = 10)) +
  labs(title = "Distribution for Number of CTCF by Resolution", x = "Resolution", y = "Number of CTCF in A Loop")

boxplot.w.CTCF.by.res

pdf("./figures/submission/boxplot_of_ctcf_count_by_res.pdf", width = 10.5, height = 13.5) #letter size
grid.arrange(boxplot.w.CTCF.by.res, ncol = 1)
dev.off()

# two figures above in one pdf (boxplot.w.CTCF.by.chr.and.res, boxplot.w.CTCF.by.res)
pdf("./figures/submission/boxplot_of_ctcf_count.pdf", width = 16.5, height = 23.5)
grid.arrange(boxplot.w.CTCF.by.chr.and.res, boxplot.w.CTCF.by.res, ncol = 1)
dev.off()

########################
# 2. CTCF
# 2-4. Distribution of CTCF at ends in a loop: for the number of CTCF used in filtering valid loops: figures
# so, the object should be used one without padding: df.DISTINCT.loop.deep.sample.all
# 2-4-5. Distribution of CTCF at each end in a loop: figure
########################
# deta processing
df.overlapping.CTCF.w.BOTH.result.boxplot.ALL.chr <- df.overlapping.CTCF.w.BOTH.result %>% 
  group_by(loop.id, WHERE, resolution) %>%
  summarise(ctcf_count_by_loop_id = n_distinct(ctcf.id), .groups = 'drop')

df.overlapping.CTCF.w.BOTH.result.boxplot.ALL.chr

# figure by end
boxplot.w.CTCF.ALL.chr <- ggplot(df.overlapping.CTCF.w.BOTH.result.boxplot.ALL.chr, aes(x = WHERE, y = ctcf_count_by_loop_id, fill = WHERE)) +
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +  # hiding outlier
  stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 1)), 
               position = position_dodge(width = 0.75), vjust = -0.5, size = 2.5) +
  stat_summary(fun.data = function(y) {
    data.frame(y = median(y), label = paste0("Q1: ", round(quantile(y, 0.25), 1), "\nQ3: ", round(quantile(y, 0.75), 1)))
  }, geom = "text", aes(label = after_stat(label)), position = position_dodge(width = 0.75), vjust = 1.5, size = 2.5) +
  theme_minimal() + 
  theme(
    plot.title = element_text(hjust = 0.5),  # title center
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("UP" = "yellow", "DOWN" = "purple")) +
  scale_y_continuous(breaks = seq(min(df.overlapping.CTCF.w.BOTH.result.boxplot.ALL.chr$ctcf_count_by_loop_id), max(df.overlapping.CTCF.w.BOTH.result.boxplot.ALL.chr$ctcf_count_by_loop_id), by = 5)) +  # y축 눈금 5단위
  coord_cartesian(ylim = c(quantile(df.overlapping.CTCF.w.BOTH.result.boxplot.ALL.chr$ctcf_count_by_loop_id, 0.05), 
                           quantile(df.overlapping.CTCF.w.BOTH.result.boxplot.ALL.chr$ctcf_count_by_loop_id, 0.95))) +  # y-axis range
  labs(title = "Distribution for Number of CTCF by Up/Downstream End", x = "Upstream and Downstream End", y = "Number of CTCF")

boxplot.w.CTCF.ALL.chr

# figure by end & res
boxplot.w.CTCF.ALL.chr.by.res <- ggplot(df.overlapping.CTCF.w.BOTH.result.boxplot.ALL.chr, aes(x = WHERE, y = ctcf_count_by_loop_id, fill = resolution)) +
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +  # removing outliers
  stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 1)), 
               position = position_dodge(width = 0.75), vjust = -0.5, size = 2.5) +
  stat_summary(fun.data = function(y) {
    data.frame(y = median(y), label = paste0("Q1: ", round(quantile(y, 0.25), 1), "\nQ3: ", round(quantile(y, 0.75), 1)))
  }, geom = "text", aes(label = after_stat(label)), position = position_dodge(width = 0.75), vjust = 1.5, size = 2.5) +
  theme_minimal() + 
  theme(
    plot.title = element_text(hjust = 0.5),  # title center
    legend.position = "none"  # removing legend
  ) +
  scale_fill_manual(values = c("5K" = "skyblue", "10K" = "lightcoral", "25K" = "springgreen"), breaks = c("5K", "10K", "25K")) +
  scale_y_continuous(breaks = seq(min(df.overlapping.CTCF.w.BOTH.result.boxplot.ALL.chr$ctcf_count_by_loop_id), 
                                  max(df.overlapping.CTCF.w.BOTH.result.boxplot.ALL.chr$ctcf_count_by_loop_id), by = 5)) +  # y-axis tip
  coord_cartesian(ylim = c(quantile(df.overlapping.CTCF.w.BOTH.result.boxplot.ALL.chr$ctcf_count_by_loop_id, 0.05), 
                           quantile(df.overlapping.CTCF.w.BOTH.result.boxplot.ALL.chr$ctcf_count_by_loop_id, 0.95))) + 
  labs(title = "Boxplot for Number of CTCF in UP/DOWNSTREAM End by Resolution", x = "Upstream and Downstream End", y = "Number of CTCF", fill = "Resolution")


boxplot.w.CTCF.ALL.chr.by.res

pdf("./figures/submission/boxplot_for_no.ctcf_by_end_res.pdf", width = 16.5, height = 23.5)
grid.arrange(boxplot.w.CTCF.ALL.chr, boxplot.w.CTCF.ALL.chr.by.res, ncol = 1)
dev.off()


##########################################################
# 2. CTCF
# 2-5.ideogram

library(RIdeogram)
library(rsvg)
library(magick)

getwd()

df.DISTINCT.fimo.2nd.trial.ctcf
df.DISTINCT.fimo.2nd.trial.ctcf %>% 
  count(length)

chromosome_data <- tss<-read.table(file="/Users/PanjunKim/UTHSC\ GGI\ Dropbox/K\ P/Gateway_to_Hao/enhancer/data/rn7_chromosome_length.tsv", sep="\t", head=T) %>% 
  mutate(start = as.numeric(start), end = str_remove_all(end, ",") %>% as.numeric())
chromosome_data %>% head()
chromosome_data <- chromosome_data %>% 
  mutate(start = start - 1)
chromosome_data 

bin_size <- 1000000
df_ctcf <- df.DISTINCT.fimo.2nd.trial.ctcf %>% 
  dplyr::select(chr, start, end) %>% 
  mutate(chr = str_remove(chr, "chr"))
df_ctcf

process_chromosome_bins <- function(df_ctcf, chromosome_ends, bin_size) {
  result <- list()
  
  for (chr in unique(df_ctcf$chr)) {
    # end of chromosome coordinate
    chr_data <- df_ctcf %>% filter(chr == !!chr)
    chr_end <- chromosome_ends %>% filter(chr == !!chr) %>% pull(end)
    
    if (length(chr_end) != 1) {
      stop(paste("Invalid chr_end length for", chr)) # error
    }
    
    # bin
    bins <- data.frame(
      Start = seq(1, chr_end, by = bin_size),
      End = c(seq(bin_size, chr_end, by = bin_size), chr_end)
    )
    bins <- bins %>%
      mutate(chr = chr, Value = 0)
    
    # getting # of CTCF in a bin
    for (i in 1:nrow(bins)) {
      bins$Value[i] <- sum(
        chr_data$start >= bins$Start[i] & chr_data$start <= bins$End[i]
      )
    }
    result[[chr]] <- bins
  }
  
  bind_rows(result)
}

ctcf_density <- process_chromosome_bins(df_ctcf, chromosome_data, bin_size)

ctcf_density <- ctcf_density %>%
  dplyr::select(chr, Start, End, Value) %>%
  dplyr::rename(Chr = chr) %>% 
  arrange(chr, Start)
chromosome_data <- chromosome_data %>% 
  dplyr::rename(Chr = chr, Start = start, End = end)





chromosome_data %>% head()


ideogram(karyotype = chromosome_data, overlaid = ctcf_density, colorset1 = c("#008000", "#FFA500", "#d73027"))
# SVG file save
svg_file <- "chromosome.svg"
convertSVG(svg_file, device = "png")

# converting to PNG with size
png_file <- "chromosome_high_res.png"
rsvg::rsvg_png(svg_file, file = png_file, width = 6000, height = 6000)

image <- image_read(png_file) %>%             # reading PNG
  image_background(color = "white") %>%       # white background
  image_write(png_file)     

data(liriodendron_karyotype, package="RIdeogram") #load the karyotype data
data(Fst_between_CE_and_CW, package="RIdeogram") #load the Fst data for overlaid heatmap
data(Pi_for_CE_and_CW, package="RIdeogram")
ideogram(karyotype = liriodendron_karyotype, overlaid = Fst_between_CE_and_CW, label = Pi_for_CE_and_CW, label_type = "polygon", colorset1 = c("#e5f5f9", "#99d8c9", "#2ca25f"))
convertSVG("chromosome.svg", device = "png")


liriodendron_karyotype %>% head()
Fst_between_CE_and_CW %>% head()
Pi_for_CE_and_CW %>% head()


liriodendron_karyotype

##########################################################################################
# refactoring code
##########################################################################################

df.DISTINCT.loop.deep.sample.all.50.padded <- df_50_padded %>% 
  mutate(padded.distance = 50) %>% 
  mutate(x1 = x0, x2 = x3, y1 = y0, y2 = y3) %>% 
  mutate(loop.id = str_c(loop.id, '_', 50))

df.DISTINCT.loop.deep.sample.all.100.padded <- df_100_padded %>% 
  mutate(padded.distance = 100) %>% 
  mutate(x1 = x0, x2 = x3, y1 = y0, y2 = y3) %>% 
  mutate(loop.id = str_c(loop.id, '_', 100))

# analysis func
perform_analysis <- function(loop_data, ctcf_data, padding_label) {
  
  # 1. GRanges object
  loop_up_GR <- GRanges(seqnames=loop_data$chr1, 
                        ranges=IRanges(start=(loop_data$x1), end=(loop_data$x2)), 
                        id=loop_data$loop.id, 
                        end.distance = loop_data$end.distance, 
                        resolution = loop_data$resolution, 
                        distance = loop_data$distance)
  
  loop_down_GR <- GRanges(seqnames=loop_data$chr1, 
                          ranges=IRanges(start=(loop_data$y1), end=(loop_data$y2)), 
                          id=loop_data$loop.id, 
                          end.distance = loop_data$end.distance, 
                          resolution = loop_data$resolution, 
                          distance = loop_data$distance)
  
  # 2. CTCF + UPSTREAM
  index_ctcf_w_up_loop <- findOverlaps(ctcf_data, loop_up_GR, type = "within")
  end_loop_up_ctcf_hits <- subjectHits(index_ctcf_w_up_loop)
  end_ctcf_up_hits <- queryHits(index_ctcf_w_up_loop)
  
  df_up_result <- data.frame(
    up.loop.id = mcols(loop_up_GR)$id[end_loop_up_ctcf_hits],
    end.up.distance = mcols(loop_up_GR)$end.distance[end_loop_up_ctcf_hits],
    distance = mcols(loop_up_GR)$distance[end_loop_up_ctcf_hits],
    resolution = mcols(loop_up_GR)$resolution[end_loop_up_ctcf_hits],
    ctcf.id = mcols(ctcf_data)$id[end_ctcf_up_hits],
    WHERE = "UP"
  ) %>% 
    mutate(ctcf.loop.up.id = str_c(up.loop.id, '|', distance, '|', ctcf.id, '|', WHERE))
  
  # 3. CTCF + DOWNSTREAM
  index_ctcf_w_down_loop <- findOverlaps(ctcf_data, loop_down_GR, type = "within")
  end_loop_down_ctcf_hits <- subjectHits(index_ctcf_w_down_loop)
  end_ctcf_down_hits <- queryHits(index_ctcf_w_down_loop)
  
  df_down_result <- data.frame(
    down.loop.id = mcols(loop_down_GR)$id[end_loop_down_ctcf_hits],
    end.down.distance = mcols(loop_down_GR)$end.distance[end_loop_down_ctcf_hits],
    distance = mcols(loop_down_GR)$distance[end_loop_down_ctcf_hits],
    resolution = mcols(loop_down_GR)$resolution[end_loop_down_ctcf_hits],
    ctcf.id = mcols(ctcf_data)$id[end_ctcf_down_hits],
    WHERE = "DOWN"
  ) %>% 
    mutate(ctcf.loop.down.id = str_c(down.loop.id, '|', distance, '|', ctcf.id, '|', WHERE))
  
  # 4. UPSTREAM & DOWNSTREAM bind_rows
  df_both_result <- bind_rows(
    df_up_result %>% 
      mutate(loop.id = up.loop.id, end.distance = end.up.distance) %>% 
      mutate(case.id = ctcf.loop.up.id) %>% 
      dplyr::select(-c(up.loop.id, end.up.distance, ctcf.loop.up.id)), 
    df_down_result %>% 
      mutate(loop.id = down.loop.id, end.distance = end.down.distance) %>% 
      mutate(case.id = ctcf.loop.down.id) %>% 
      dplyr::select(-c(down.loop.id, end.down.distance, ctcf.loop.down.id))
  ) %>% 
    mutate(chr = ifelse(WHERE == "UP", str_split_n(loop.id, '_', 1), str_split_n(loop.id, '_', 4))) %>% 
    mutate(chr = factor(chr, levels = c(paste0("chr", 1:20), "chrX", "chrY"))) %>%
    mutate(WHERE = fct_relevel(WHERE, "UP", "DOWN")) %>% 
    mutate(resolution = case_when(
      end.distance == 5000 ~ "5K",
      end.distance == 10000 ~ "10K",
      end.distance == 25000 ~ "25K",
      TRUE ~ NA_character_
    )) %>% 
    mutate(resolution = fct_relevel(resolution, "5K", "10K", "25K"))
  
  # 5. df for boxplot
  df_boxplot <- df_both_result %>% 
    group_by(chr, loop.id, WHERE, resolution) %>%
    summarise(ctcf_count_by_loop_id = n_distinct(ctcf.id), .groups = 'drop')
  
  return(df_boxplot)
}

# boxplot func
draw_boxplots <- function(df_boxplot, padding_label, min_value, max_value) {
  
  # by Chromosome and Resolution
  boxplot_by_chr_res <- df_boxplot %>% 
    ggplot(aes(x = chr, y = ctcf_count_by_loop_id, fill = resolution)) +
    geom_boxplot() +
    stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 1)), 
                 position = position_dodge(width = 0.75), vjust = -0.5, size = 2) + 
    stat_summary(fun.data = function(y) {
      data.frame(y = median(y), label = paste0("Q1: ", round(quantile(y, 0.25), 1)))
    }, geom = "text", aes(label = after_stat(label)), 
    position = position_dodge(width = 0.75), vjust = 6.5, size = 1.5) +
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = c("5K" = "skyblue", "10K" = "lightcoral", "25K" = "springgreen"), breaks = c("5K", "10K", "25K")) +
    scale_y_continuous(breaks = seq(min_value, max_value, by = 10)) +
    labs(title = paste("Boxplot for Number of CTCF by Chromosome and Resolution (", padding_label, " Padding)", sep = ""), 
         x = "Chromosomes", y = "Number of CTCF in a loop")
  
  # by Resolution only
  boxplot_by_res <- df_boxplot %>% 
    ggplot(aes(x = resolution, y = ctcf_count_by_loop_id, fill = resolution)) +
    geom_boxplot() +
    stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 1)), 
                 position = position_dodge(width = 0.75), vjust = -0.5, size = 2) + 
    stat_summary(fun.data = function(y) {
      data.frame(y = median(y), label = paste0("Q1: ", round(quantile(y, 0.25), 1)))
    }, geom = "text", aes(label = after_stat(label)), 
    position = position_dodge(width = 0.75), vjust = 6.5, size = 1.5) +
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = c("5K" = "skyblue", "10K" = "lightcoral", "25K" = "springgreen"), breaks = c("5K", "10K", "25K")) +
    scale_y_continuous(breaks = seq(min_value, max_value, by = 10)) +
    labs(title = paste("Boxplot for Number of CTCF by Resolution (", padding_label, " Padding)", sep = ""), 
         x = "Resolution", y = "Number of CTCF in a loop")
  
  return(list(boxplot_by_chr_res = boxplot_by_chr_res, boxplot_by_res = boxplot_by_res))
}

# 1. perform_analysis
df_ctcf_boxplot_none <- perform_analysis(df.DISTINCT.loop.deep.sample.all, df.DISTINCT.ctcf.2nd.fimo.GR, "No")
df_ctcf_boxplot_50   <- perform_analysis(df.DISTINCT.loop.deep.sample.all.50.padded, df.DISTINCT.ctcf.2nd.fimo.GR, "50%")
df_ctcf_boxplot_100  <- perform_analysis(df.DISTINCT.loop.deep.sample.all.100.padded, df.DISTINCT.ctcf.2nd.fimo.GR, "100%")

# df_boxplot_none min/max 
min_value_ctcf_none <- min(df_ctcf_boxplot_none$ctcf_count_by_loop_id, na.rm = TRUE)
max_value_ctcf_none <- max(df_ctcf_boxplot_none$ctcf_count_by_loop_id, na.rm = TRUE)

# df_boxplot_50 min/max
min_value_ctcf_50 <- min(df_ctcf_boxplot_50$ctcf_count_by_loop_id, na.rm = TRUE)
max_value_ctcf_50 <- max(df_ctcf_boxplot_50$ctcf_count_by_loop_id, na.rm = TRUE)

# df_boxplot_100 min/max
min_value_ctcf_100 <- min(df_ctcf_boxplot_100$ctcf_count_by_loop_id, na.rm = TRUE)
max_value_ctcf_100 <- max(df_ctcf_boxplot_100$ctcf_count_by_loop_id, na.rm = TRUE)

# 2. boxplot
ctcf_boxplots_none <- draw_boxplots(df_ctcf_boxplot_none, "No", min_value_ctcf_none, max_value_ctcf_none)
ctcf_boxplots_50 <- draw_boxplots(df_ctcf_boxplot_50, "50%", min_value_ctcf_50, max_value_ctcf_50)
ctcf_boxplots_100 <- draw_boxplots(df_ctcf_boxplot_100, "100%", min_value_ctcf_100, max_value_ctcf_100)

# 3. PDF
pdf("./figures/0909/end_boxplot_num_ctcf_combined_no_padding.pdf", width = 16.5, height = 23.5)
grid.arrange(ctcf_boxplots_none$boxplot_by_chr_res, ctcf_boxplots_none$boxplot_by_res, ncol = 1)
dev.off()

pdf("./figures/0909/end_boxplot_num_ctcf_combined_50_padding.pdf", width = 16.5, height = 23.5)
grid.arrange(ctcf_boxplots_50$boxplot_by_chr_res, ctcf_boxplots_50$boxplot_by_res, ncol = 1)
dev.off()

pdf("./figures/0909/end_boxplot_num_ctcf_combined_100_padding.pdf", width = 16.5, height = 23.5)
grid.arrange(ctcf_boxplots_100$boxplot_by_chr_res, ctcf_boxplots_100$boxplot_by_res, ncol = 1)
dev.off()

pdf("./figures/0909/end_boxplot_num_ctcf_combined_all.pdf", width = 16.5, height = 23.5)
grid.arrange(ctcf_boxplots_none$boxplot_by_chr_res, ctcf_boxplots_none$boxplot_by_res, 
             ctcf_boxplots_50$boxplot_by_chr_res, ctcf_boxplots_50$boxplot_by_res,
             ctcf_boxplots_100$boxplot_by_chr_res, ctcf_boxplots_100$boxplot_by_res,
             ncol = 2)
dev.off()

# stats
ctcf_boxplots_none
ctcf_boxplots_50
ctcf_boxplots_100

calculate_boxplot_stats <- function(df) {
  df %>%
    group_by(resolution) %>%
    summarize(
      Q1 = quantile(ctcf_count_by_loop_id, 0.25),
      Median = median(ctcf_count_by_loop_id),
      Q3 = quantile(ctcf_count_by_loop_id, 0.75)
    )
}

ctcf_stats_none <- calculate_boxplot_stats(df_ctcf_boxplot_none)
ctcf_stats_50   <- calculate_boxplot_stats(df_ctcf_boxplot_50)
ctcf_stats_100  <- calculate_boxplot_stats(df_ctcf_boxplot_100)
#      ctcf_stats_none                        ctcf_stats_50                   ctcf_stats_100
# resolution      Q1 Median    Q3  |  resolution    Q1 Median    Q3  | resolution    Q1 Median    Q3
# 1 5K             7     26    35  |  5K            19     34    53  |  5K            27     43    66
# 2 10K           12     30    46  |  10K           28     45    72  |  10K           37     61    94
# 3 25K           25     43    71  |  25K           47     78   120  |  25K           67    108   164
=======
>>>>>>> d754844 (changing variables, generating figures for submission, splitting file of enhancer_synced.R into three.)
