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

##############################################################################################
## LOOP processing: BASIC
##############################################################################################

# BEDPE file list (10 files)
# Linux
# loop.file.list = fs::dir_ls("/home/pkim/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/", regexp = ".bedpe$")
# loop.file.list = fs::dir_ls("/home/pkim/Desktop/temp/enhancer/dropbox_enhancer_doosan/data/loops/", regexp = ".bedpe$")
loop.file.list = fs::dir_ls("/home/pkim/dropbox/Gateway_to_Hao/enhancer/data/loops", regexp = ".bedpe$")
loop.file.list

# Mac
loop.file.list = fs::dir_ls("/Users/PanjunKim/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/", regexp = ".bedpe$")
loop.file.list = fs::dir_ls("/Users/PanjunKim/Desktop/temp/enhancer/dropbox_enhancer_doosan/data/loops/", regexp = ".bedpe$")
loop.file.list = fs::dir_ls("/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/data/loops", regexp = ".bedpe$")
loop.file.list

# BED fild for loops
df.init.loop.bed <- init.bedpe.df(loop.file.list, 'loops') %>% 
  filter(!str_detect(X.chr1, "^#")) 
# %>% #58,992
#   view()
df.init.loop.bed

df.init.loop.bed %>% 
  count(sample)

# sample    n            strain    n
# 1   592BB 5263           BXH6 6568
# 2     607 7336          Bn-Lx 6535
# 3    74AA 2903       F344/Stm 2903
# 4    A2DB 2992          HXB10 7336
# 5   D765A 6568           HXB2 4656
# 6   DA08A 4656          HXB23 7676
# 7   DA21A 5932          HXB31 9131
# 8   DA68A 9131         LE/Stm 2992
# 9   DBA9A 7676    SHR/OlaIpcv 5263
# 10  DE8BA 6535           SOBN 5932

# creating loop.id, sample.loop.id 
df.loop.deep.sample.all <- df.init.loop.bed %>% 
  mutate(strain = case_when(
    sample == '592BB' ~ "SHR/OlaIpcv",
    sample == '607' ~ "HXB10",
    sample == '74AA' ~ "F344/Stm",
    sample == 'A2DB' ~ "LE/Stm",
    sample == 'D765A' ~ "BXH6",
    sample == 'DE8BA' ~ "Bn-Lx",
    sample == 'DBA9A' ~ "HXB23",
    sample == 'DA21A' ~ "SOBN",
    sample == 'DA08A' ~ "HXB2",
    sample == 'DA68A' ~ "HXB31",
    TRUE ~ NA
  )) %>% 
  mutate(end.distance = x2 - x1) %>% 
  mutate(resolution = case_when(
    end.distance == 5000 ~ "5K",
    end.distance == 10000 ~ "10K",
    end.distance == 25000 ~ "25K",
    TRUE ~ NA
  )) %>% 
  mutate(resolution = factor(resolution, levels = c("5K", "10K", "25K"))) %>%
  mutate(loop.id = str_c(X.chr1, '_', x1, '_', x2, '_', chr2, '_', y1, '_', y2, '_', end.distance)) %>% # loop.id
  mutate(sample.loop.id = str_c(strain, '_', loop.id)) %>% #sample.loop.id
  dplyr::select(sample, strain, X.chr1, x1, x2, chr2, y1, y2, distance, end.distance, resolution, loop.id, sample.loop.id) %>% 
  dplyr::rename(chr1 = X.chr1)

df.loop.deep.sample.all %>% 
  head()
df.loop.deep.sample.all %>% 
  count(strain)

########################
# 1. Loop
# 1-1. Exploratory Data analysis (EDA) 
########################
# 1. checking duplication loops in a strain
df.loop.deep.sample.all %>% 
  #count(strain, loop.id) %>% # NO dups in a strain
  count(strain, resolution, loop.id) %>% # NO dup in the same resolution in a strain
  filter(n > 1) 

# 2. checking how much common loops are in samples
resolutions <- c("5K", "10K", "25K")
plots <- list()

for(res in resolutions) {
  df.loop.deep.sample.all.res <- df.loop.deep.sample.all %>% filter(resolution == res)
  
  loop_counts <- df.loop.deep.sample.all.res %>% group_by(strain) %>% summarise(total_loops = n_distinct(loop.id))
  
  location_list <- split(df.loop.deep.sample.all.res$loop.id, df.loop.deep.sample.all.res$strain)
  
  common_counts <- matrix(0, nrow = length(location_list), ncol = length(location_list))
  rownames(common_counts) <- colnames(common_counts) <- names(location_list)
  
  for(i in 1:length(location_list)) {
    for(j in 1:length(location_list)) {
      common <- length(intersect(location_list[[i]], location_list[[j]]))
      percentage <- (common / length(location_list[[i]])) * 100
      common_counts[i, j] <- percentage
    }
  }
  
  # long
  common_counts_long <- melt(common_counts)
  
  # heatmap
  p <- ggplot(common_counts_long, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.1f%%", value)), color = "black", size = 4) +
    scale_fill_gradient(low = "white", high = "darkred", 
                        name="Common Loop %", limits = c(0, 100)) +
    theme_minimal() +
    labs(x = "Strain", y = "Strain", title = paste("Heatmap of Common Loops Percentage -", res, "Resolution")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))  # align title in the center
  
  # total loops
  loop_counts_table <- loop_counts %>% 
    mutate(Strain = factor(strain, levels = names(location_list))) %>% 
    arrange(Strain)
  
  colnames(common_counts) <- paste0(colnames(common_counts), "\n(", loop_counts$total_loops, ")")
  rownames(common_counts) <- paste0(rownames(common_counts), "\n(", loop_counts$total_loops, ")")
  
  plots[[res]] <- p
}

# all figures by resolution into one figure
pdf("figures/submission/common_loops_heatmap_percentage_bw_strains.pdf", width = 12, height = 12)
grid.arrange(plots[["5K"]], plots[["10K"]], plots[["25K"]], nrow = 3)
dev.off()

########################
# 1. Loop
# 1-2. Loop data preprocessing: distinct loops
########################

df.loop.deep.sample.all %>% 
  #  head()
  count() # 58992

# common loops between samples picked only one: df.DISTINCT.loop.deep.sample.all
df.DISTINCT.loop.deep.sample.all <- df.loop.deep.sample.all %>%
  dplyr::select(loop.id) %>%
  distinct() %>% 
  separate(loop.id, into = c("chr1", "x1", "x2", "chr2", "y1", "y2", "end.distance"), sep = "_", remove = FALSE, convert = TRUE) %>% 
  mutate(distance = y2 - x2) %>%
  mutate(resolution = case_when(
    end.distance == 5000 ~ "5K",
    end.distance == 10000 ~ "10K",
    end.distance == 25000 ~ "25K",
    TRUE ~ NA
  )) %>% 
  mutate(resolution = factor(resolution, levels = c("5K", "10K", "25K"))) %>% 
  mutate(padded.distance = 0) %>% 
  mutate(loop.id = str_c(loop.id, '_', padded.distance))

df.DISTINCT.loop.deep.sample.all

df.DISTINCT.loop.deep.sample.all %>% 
  count(resolution) # 31773/58992
head()
# resolution     n
# 1         5K  6680
# 2        10K 12162
# 3        25K 12931

########################
# 1. Loop
# 1-3. Loop data preprocessing: Creating GRange Obj. from df.DISTINCT.loop.deep.sample.all (whole, up, down)
########################
t = 0
# function for creating GRange Obj.
create_granges <- function(start, end) {
  GRanges(
    seqnames = as.character(df.DISTINCT.loop.deep.sample.all$chr1),  
    ranges = IRanges(start = start, end = end),  
    id = df.DISTINCT.loop.deep.sample.all$loop.id,  
    end.distance = df.DISTINCT.loop.deep.sample.all$end.distance,
    resolution = df.DISTINCT.loop.deep.sample.all$resolution,
    distance = df.DISTINCT.loop.deep.sample.all$distance
  )
}

# whole (x1 - y2) 
df.DISTINCT.loop.deep.sample.all.whole.GR <- create_granges(
  start = df.DISTINCT.loop.deep.sample.all$x1 - t,
  end = df.DISTINCT.loop.deep.sample.all$y2 + t
)
# up (x1 - x2)
df.DISTINCT.loop.deep.sample.all.up.GR <- create_granges(
  start = df.DISTINCT.loop.deep.sample.all$x1 - t,
  end = df.DISTINCT.loop.deep.sample.all$x2 + t
)
# down (y1 - y2)
df.DISTINCT.loop.deep.sample.all.down.GR <- create_granges(
  start = df.DISTINCT.loop.deep.sample.all$y1 - t,
  end = df.DISTINCT.loop.deep.sample.all$y2 + t
)
# middle 
df.DISTINCT.loop.deep.sample.all.middle.GR <- create_granges(
  start = df.DISTINCT.loop.deep.sample.all$x2 - t,
  end = df.DISTINCT.loop.deep.sample.all$y1 + t
)

########################
# 1. Loop
# 1-4. Loop data preprocessing: padding on loops (1 distance, .5 distance) -> overall.df.DISTINCT.loop.deep.sample.all : OBJECT to be used for the distribution of sth over loops
# x0, y3, new_distance, new.loop.id
########################
overall.df.DISTINCT.loop.deep.sample.all.5.distance <- df.DISTINCT.loop.deep.sample.all %>% 
  mutate(x12 = (x1 + x2)/2, y12 = (y1 + y2)/2) %>% # middle point of each end
  mutate(x0 = ifelse(x12 - (distance / 2) < 0, 0, x12 - (distance / 2)), y3 = y12 + (distance/2)) %>% # 1/2 distance for padding
  # mutate(new_distance = abs(y3 - x0)) %>% 
  mutate(new.loop.id = str_c(loop.id, '|', x0, '_', y3))

overall.df.DISTINCT.loop.deep.sample.all.1.distance <- df.DISTINCT.loop.deep.sample.all %>% 
  mutate(x12 = (x1 + x2)/2, y12 = (y1 + y2)/2) %>% # middle point of each end
  mutate(x0 = ifelse(x12 - (distance) < 0, 0, x12 - (distance)), y3 = y12 + (distance)) %>% # 1 distance for padding
  # mutate(new_distance = abs(y3 - x0)) %>% 
  mutate(new.loop.id = str_c(loop.id, '|', x0, '_', y3))

# new.df.loop.deep.sample.all %>%
# overall.df.DISTINCT.loop.deep.sample.all.5.distance %>% 
# overall.df.DISTINCT.loop.deep.sample.all.1.distance %>% 
# filter(x0 < 0) %>%   
#  head()

# data integrity: PASS
# overall.df.DISTINCT.loop.deep.sample.all %>% mutate(test = ifelse(y3-x0 == 2*distance, TRUE, FALSE)) %>% 
#   filter(test == FALSE) %>% 
#   count(x0)
# count(test)

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
# 4. CTCF
# 4-4. Distribution of CTCF at each end in a loop: for the number of CTCF used in filtering valid loops: figures
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
# 4. TSS
# 2-4. Distribution of TSS at ends in a loop: for the number of TSS used in filtering valid loops: figures
# so, the object should be used one without padding: df.DISTINCT.loop.deep.sample.all
# 2-4-3. data processing for getting information of TSS at ends in loops 
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

#########################################################
# 5. diagram with loops from 3 previous steps
#########################################################
# distinct loops: 31773

df_ctcf_boxplot_none
df_ctcf_boxplot_50  
df_ctcf_boxplot_100 

ctcf_stats_none
ctcf_stats_50  
ctcf_stats_100 

#      ctcf_stats_none                        ctcf_stats_50                   ctcf_stats_100
# resolution      Q1 Median    Q3  |  resolution    Q1 Median    Q3  | resolution    Q1 Median    Q3
# 1 5K             7     26    35  |  5K            19     34    53  |  5K            27     43    66
# 2 10K           12     30    46  |  10K           28     45    72  |  10K           37     61    94
# 3 25K           25     43    71  |  25K           47     78   120  |  25K           67    108   164

# 4-1. CTCF
# important object
df.overlapping.CTCF.w.BOTH.result %>% 
  dim() # 2420002
head()
df.ctcf.counts

process_final_loops <- function(df_ctcf_counts, ctcf_stats) {
  q1_values_by_resolution <- ctcf_stats %>%
    dplyr::select(resolution, Q1)
  # resolution    Q1
  # 1 5K             7
  # 2 10K           12
  # 3 25K           25
  
  loops_with_ctcf_above_q1 <- df_ctcf_counts %>%
    left_join(q1_values_by_resolution, by = "resolution") %>%
    filter(ctcf_count_by_loop_id >= Q1)
  loops_with_ctcf_above_q1 # 46,851
  
  loops.with.ctcf.both.ends <- loops_with_ctcf_above_q1 %>%
    group_by(loop.id) %>%
    filter(all(c("UP", "DOWN") %in% WHERE)) %>%
    ungroup()
  loops.with.ctcf.both.ends # 36,664 + 10
  
  final.loops.from.ctcf.step <- loops.with.ctcf.both.ends %>% 
    distinct(loop.id) %>% 
    mutate(end.distance = as.numeric(str_split_n(loop.id, '_', 7))) %>% 
    mutate(resolution = case_when(
      end.distance == 5000 ~ "5K",
      end.distance == 10000 ~ "10K",
      end.distance == 25000 ~ "25K",
      TRUE ~ NA
    ))
  # %>%  
  #   mutate(loop.id = str_remove(loop.id, "_[^_]+$"))
  return(final.loops.from.ctcf.step)
}

ctcf_final_loops_none <- process_final_loops(df_ctcf_boxplot_none, ctcf_stats_none)
ctcf_final_loops_50 <- process_final_loops(df_ctcf_boxplot_50, ctcf_stats_50)
ctcf_final_loops_100 <- process_final_loops(df_ctcf_boxplot_100, ctcf_stats_100)

ctcf_final_loops_none # 18337                  // initial loops: 31773
ctcf_final_loops_50   # 19298 ( + 961)         // initial loops: 31773
ctcf_final_loops_100  # 19392 (+ 1055 / + 94 ) // initial loops: 31773

############################
############################
############################
# loop analysis by number of ctcf
df.DISTINCT.loop.deep.sample.all %>% # 31773 (UP/DOWN = 63546)
  head(2)
df.ctcf.counts # 61713

####### Analysis 1: These 198 loops means they DO NOT have CTCF in their ends #######
missing_loops <- anti_join(df.DISTINCT.loop.deep.sample.all,
                           df.ctcf.counts %>% dplyr::select(loop.id) %>% distinct(), by = "loop.id")
missing_loops # 198 loop.id
distribution_of_loops_without_ctcf_bindings <- missing_loops %>% 
  count(str_split_n(loop.id, '_', 1)) %>%
  dplyr::rename(chr = `str_split_n(loop.id, "_", 1)`, count = n) %>% 
  mutate(chr = factor(chr, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                                      "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", 
                                      "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                                      "chr20", "chrX", "chrY"))) %>% 
  ggplot(aes(x = chr, y = count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  labs(title = "Distribution of Loops Without CTCF Bindings on Chrnomosome", x = "Chromosome", y = "Count") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
ggsave(filename = "./figures/0909/distribution_of_loops_without_ctcf_bindings.pdf", plot = distribution_of_loops_without_ctcf_bindings, width = 8, height = 6)
missing_loops

distribution_of_loops_without_ctcf_bindings_per_resolution <- missing_loops %>%
  mutate(end.distance = str_split_n(loop.id, '_', 7)) %>% 
  mutate(resolution = case_when(
    end.distance == 5000 ~ "5K",
    end.distance == 10000 ~ "10K",
    end.distance == 25000 ~ "25K",
    TRUE ~ NA
  )) %>% 
  group_by(loop.id) %>%
  ungroup() %>% 
  mutate(chr = str_split_n(loop.id, '_', 1)) %>% 
  count(chr, resolution) %>% 
  mutate(chr = factor(chr, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                                      "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", 
                                      "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                                      "chr20", "chrX", "chrY"))) %>% 
  ggplot(aes(x = chr, y = n, fill = resolution)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = c("5K" = "skyblue", "10K" = "orange", "25K" = "lightgreen")) +
  labs(title = "Distribution of Loops without CTCF by Chromosome and Resolution",
       x = "Chromosome",
       y = "Total Count",
       fill = "Resolution") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave(filename = "./figures/0909/distribution_of_loops_without_ctcf_bindings_per_resolution.pdf", plot = distribution_of_loops_without_ctcf_bindings_per_resolution, width = 8, height = 6)

####### Analysis 2: These 1437 loops means they DO NOT have CTCF in their 'ONE ENDs ONLY': 63150 - 61713 = 1437 #######
df.ctcf.counts.paired.completed <- df.ctcf.counts %>%
  complete(loop.id, WHERE = c("UP", "DOWN"), fill = list(ctcf_count = 0)) %>% 
  mutate(resolution = case_when(
    str_split_n(loop.id, '_', 7) == 5000 ~ "5K",
    str_split_n(loop.id, '_', 7) == 10000 ~ "10K",
    str_split_n(loop.id, '_', 7) == 25000 ~ "25K",
    TRUE ~ NA
  ))

df.ctcf.counts.paired.completed # 63150 <- 61713
# 63546 - 63150 = 396 (198 loops) (no ctcf in both ends) - identical result to Analaysis 1: PASS

loops_only_one_end_per_resolution <- df.ctcf.counts %>%
  group_by(loop.id) %>%
  filter(n_distinct(WHERE) < 2) %>%
  ungroup() %>% 
  mutate(chr = str_split_n(loop.id, '_', 1)) %>% 
  count(chr, resolution) %>% 
  mutate(chr = factor(chr, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                                      "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", 
                                      "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                                      "chr20", "chrX", "chrY"))) %>% 
  ggplot(aes(x = chr, y = n, fill = resolution)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = c("5K" = "skyblue", "10K" = "orange", "25K" = "lightgreen")) +
  labs(title = "Distribution of Loops by Chromosome and Resolution",
       x = "Chromosome",
       y = "Total Count",
       fill = "Resolution") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave(filename = "./figures/0909/loops_with_ctcf_only_one_end_per_resolution.pdf", plot = loops_only_one_end_per_resolution, width = 8, height = 6)

####### Analysis 3: loops less than Q1
df.ctcf.counts # 61,713
df.loops.with.ctcf.pairing.filtered <- df.ctcf.counts %>%
  group_by(loop.id) %>%
  filter(n_distinct(WHERE) == 2) %>% # n_distinct(WHERE) < 2 loop.ids are NO CTCF in the one end
  ungroup()

df.loops.with.ctcf.pairing.filtered # 60276
# 61713 - 60276 = 1437 : PASS

ctcf_stats_none
ctcf_stats_50
ctcf_stats_100

q1_values_by_resolution <- ctcf_stats_none %>%
  dplyr::select(resolution, Q1)

q1_values_by_resolution # EMPTY

loops.with.ctcf.above.q1.complement <- df.loops.with.ctcf.pairing.filtered %>% 
  left_join(q1_values_by_resolution, by = "resolution") %>%
  filter(ctcf_count < Q1)

loops.with.ctcf.above.q1.complement # 14190

loops.with.less.ctcf.than.q1.on.both.ends <- loops.with.ctcf.above.q1.complement %>%
  group_by(loop.id) %>%
  filter(n_distinct(WHERE) == 2) %>%
  ungroup() %>% 
  mutate(ends = "BOTH_ENDS")
loops.with.less.ctcf.than.q1.on.both.ends # 4778 (/2 = 2389 loops)

loops.with.less.ctcf.than.q1.on.only.one.end <- loops.with.ctcf.above.q1.complement %>%
  group_by(loop.id) %>%
  filter(n_distinct(WHERE) < 2) %>%
  mutate(ends = "ONE_END") %>% 
  ungroup()

loops.with.ctcf.on.only.one.end # 9412 (the other ends have CTCF more than q1)

# PASS: 4778 + 9412 = 14190

df.loops.with.less.ctcf.than.q1 <- bind_rows(loops.with.less.ctcf.than.q1.on.both.ends, loops.with.less.ctcf.than.q1.on.only.one.end)
df.loops.with.less.ctcf.than.q1

q1_data <- df.loops.with.less.ctcf.than.q1 %>% distinct(resolution, Q1)

# Histogram plot
df.loops.with.less.ctcf.than.q1.hist <- ggplot(df.loops.with.less.ctcf.than.q1, aes(x = ctcf_count, fill = ends)) +
  geom_histogram(binwidth = 1, position = "dodge", color = "black") +
  scale_fill_manual(values = c("BOTH_ENDS" = "skyblue", "ONE_END" = "orange")) +
  theme_minimal() +
  labs(title = "Histogram of CTCF Count by Ends and Resolution", x = "CTCF Count", y = "Frequency", fill = "Ends") +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~ resolution, scales = "free_y", nrow = 1, labeller = label_both) +
  geom_vline(data = q1_data, aes(xintercept = Q1), color = "red", linetype = "dashed", linewidth = 1) +
  geom_text(data = q1_data, aes(x = Q1, y = 0, label = paste("Q1 =", Q1)), 
            color = "red", vjust = 1, hjust = 1, inherit.aes = FALSE)

df.loops.with.less.ctcf.than.q1 %>% distinct(resolution, Q1)

ggsave(filename = "./figures/0909/distribution_of_loops_with_less_ctcf_than_q1_hist.pdf", plot = df.loops.with.less.ctcf.than.q1.hist, width = 8, height = 6)

# Save all three plots to a single PDF file
pdf("./figures/0909/plot_loops_with_ctcf_above_q1_complement_paired_by_resolution.pdf", height = 12, width = 8.5)
grid.arrange(plot.5k.loops.with.ctcf.above.q1.complement.paired, 
             plot.10k.loops.with.ctcf.above.q1.complement.paired, 
             plot.25k.loops.with.ctcf.above.q1.complement.paired, ncol = 1)
dev.off()

########################
# loops.with.above.q1
########################
# 4-1. CTCF
loops.with.ctcf.above.q1 <- df.ctcf.counts %>% 
  left_join(ctcf.stats.by.resolution %>% 
              dplyr::select(resolution, Q1), by = "resolution") %>% 
  filter(ctcf_count >= Q1)

loops.with.ctcf.both.ends <- loops.with.ctcf.above.q1 %>%
  group_by(loop.id) %>%
  filter(all(c("UP", "DOWN") %in% WHERE)) %>%
  ungroup()

final.loops.from.ctcf.step <- loops.with.ctcf.both.ends %>% 
  distinct(loop.id) %>% 
  mutate(end.distance = str_split_n(loop.id, '_', 7)) %>% 
  mutate(resolution = case_when(
    end.distance == 5000 ~ "5K",
    end.distance == 10000 ~ "10K",
    end.distance == 25000 ~ "25K",
    TRUE ~ NA
  ))
# %>% 
#   mutate(loop.id = str_remove(loop.id, "_[^_]+$"))
final.loops.from.ctcf.step # 18,327 + 10


# 4-2. TSS
df.overlapping.TSS.w.BOTH.result %>% head()
df.tss.counts.each.end # 44,553 + 10
tss.stats.by.resolution.each.end

loops.with.tss.above.q1 <- df.tss.counts.each.end %>% 
  left_join(tss.stats.by.resolution.each.end %>% 
              dplyr::select(resolution, Q1), by = "resolution") %>% 
  filter(tss_count_each_end >= Q1)

loops.with.tss.both.ends <- loops.with.tss.above.q1 %>%
  group_by(loop.id) %>%
  filter(all(c("UP", "DOWN") %in% WHERE)) %>%
  ungroup()

final.loops.from.tss.step <- loops.with.tss.both.ends %>% 
  distinct(loop.id) %>% 
  mutate(end.distance = str_split_n(loop.id, '_', 7)) %>% 
  mutate(resolution = case_when(
    end.distance == 5000 ~ "5K",
    end.distance == 10000 ~ "10K",
    end.distance == 25000 ~ "25K",
    TRUE ~ NA
  ))
# %>% 
#   mutate(loop.id = str_remove(loop.id, "_[^_]+$"))
final.loops.from.tss.step # 17,036 + 10

# 4-3. promoter
df.overlapping.promoter.w.BOTH.result
df.promoter.counts.each.end # 38,162
promoter.stats.by.resolution.each.end

loops.with.promoter.above.q1 <- df.promoter.counts.each.end %>% 
  left_join(promoter.stats.by.resolution.each.end %>% 
              dplyr::select(resolution, Q1), by = "resolution") %>% 
  filter(promoter_count_each_end >= Q1)

loops.with.promoter.both.ends <- loops.with.promoter.above.q1 %>% # 38,152 + 10
  group_by(loop.id) %>%
  filter(all(c("UP", "DOWN") %in% WHERE)) %>%
  ungroup()

final.loops.from.promoter.step <- loops.with.promoter.both.ends %>% 
  distinct(loop.id) %>% 
  mutate(end.distance = str_split_n(loop.id, '_', 7)) %>% 
  mutate(resolution = case_when(
    end.distance == 5000 ~ "5K",
    end.distance == 10000 ~ "10K",
    end.distance == 25000 ~ "25K",
    TRUE ~ NA
  ))
# %>% 
#   mutate(loop.id = str_remove(loop.id, "_[^_]+$"))
final.loops.from.promoter.step # 13,003 + 10
############################
############################
############################


##############
# Venn Diagram
##############

ctcf_final_loops_none
ctcf_final_loops_50  
ctcf_final_loops_100 

final.loops.from.promoter.step$loop.id
final.loops.from.tss.step$loop.id

create_venn_plot <- function(ctcf_data, promoter_data, tss_data, ctcf_label) {
  
  ctcf_loops <- str_extract(ctcf_data$loop.id, "^(?:[^_]+_){6}[^_]+")
  promoter_loops <- str_extract(promoter_data$loop.id, "^(?:[^_]+_){6}[^_]+")
  tss_loops <- str_extract(tss_data$loop.id, "^(?:[^_]+_){6}[^_]+")
  
  venn_data <- list(
    CTCF = ctcf_loops,
    Promoter = promoter_loops,
    TSS = tss_loops
  )
  
  ctcf_count <- length(unique(ctcf_loops))
  promoter_count <- length(unique(promoter_loops))
  tss_count <- length(unique(tss_loops))
  
  venn_plot <- ggvenn(
    venn_data, 
    fill_color = c("#E41A1C", "#377EB8", "#4DAF4A"),
    show_elements = FALSE  
  ) + 
    ggtitle(paste(ctcf_label, "vs Promoter vs TSS")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  venn_plot <- venn_plot +
    annotate("text", x = -1.5, y = 1.5, label = paste("CTCF:", ctcf_count), size = 3.5, color = "#E41A1C") +
    annotate("text", x = 1.5, y = 1.5, label = paste("Promoter:", promoter_count), size = 3.5, color = "#377EB8") +
    annotate("text", x = 0, y = -1.5, label = paste("TSS:", tss_count), size = 3.5, color = "#4DAF4A")
  
  
  return(venn_plot)
}

create_ctcf_venn_diagrams <- function(ctcf_final_loops_none, ctcf_final_loops_50, ctcf_final_loops_100,
                                      final_loops_from_promoter_step, final_loops_from_tss_step, output_path) {
  
  # Venn Diagrams from create_venn_plot
  venn_plot_none <- create_venn_plot(ctcf_final_loops_none, final_loops_from_promoter_step, final_loops_from_tss_step, "CTCF w/o padding")
  venn_plot_50 <- create_venn_plot(ctcf_final_loops_50, final_loops_from_promoter_step, final_loops_from_tss_step, "CTCF padded with 50% RES")
  venn_plot_100 <- create_venn_plot(ctcf_final_loops_100, final_loops_from_promoter_step, final_loops_from_tss_step, "CTCF padded with 100% RES")
  
  pdf(file = output_path)
  print(venn_plot_none)
  print(venn_plot_50)
  print(venn_plot_100)
  dev.off()
  # grid.arrange(
  #   venn_plot_none, 
  #   venn_plot_50, 
  #   venn_plot_100, 
  #   ncol = 3, nrow = 1
  # )
  dev.off()
}

create_ctcf_venn_diagrams(
  ctcf_final_loops_none = ctcf_final_loops_none,
  ctcf_final_loops_50 = ctcf_final_loops_50,
  ctcf_final_loops_100 = ctcf_final_loops_100,
  final_loops_from_promoter_step = final.loops.from.promoter.step,
  final_loops_from_tss_step = final.loops.from.tss.step,
  output_path = "./figures/0909/ctcf_vs_promoter_tss_venn_diagrams.pdf"
)

####################################
####################################
####################################
# Venn Diagram SUBMISSION
####################################
####################################
####################################
venn_plot_submission <- create_venn_plot(final.loops.from.ctcf.step, final.loops.from.promoter.step, final.loops.from.tss.step, "CTCF")
ggsave(filename = "./figures/submission/ctcf_vs_promoter_tss_venn_diagrams.pdf", plot = venn_plot_submission, width = 6, height = 6, dpi = 300)

venn_plot_submission

##################
# loop extraction
##################

# func for common loop extraction
extract_overlapping_loops <- function(ctcf_data, promoter_data, tss_data) {
  
  ctcf_loops <- unique(str_extract(ctcf_data$loop.id, "^(?:[^_]+_){6}[^_]+"))
  promoter_loops <- unique(str_extract(promoter_data$loop.id, "^(?:[^_]+_){6}[^_]+"))
  tss_loops <- unique(str_extract(tss_data$loop.id, "^(?:[^_]+_){6}[^_]+"))
  
  # 1. CTCF .vs Promoter (NONE TSS)
  ctcf_promoter_overlap <- intersect(ctcf_loops, promoter_loops)
  ctcf_promoter_only <- setdiff(ctcf_promoter_overlap, tss_loops)
  
  # 2. CTCF .vs TSS (NONE Promoter)
  ctcf_tss_overlap <- intersect(ctcf_loops, tss_loops)
  ctcf_tss_only <- setdiff(ctcf_tss_overlap, promoter_loops)
  
  # 3. CTCF, Promoter, TSS
  ctcf_promoter_tss_overlap <- Reduce(intersect, list(ctcf_loops, promoter_loops, tss_loops))
  
  return(list(
    ctcf_promoter_only = ctcf_promoter_only,
    ctcf_tss_only = ctcf_tss_only,
    ctcf_promoter_tss_overlap = ctcf_promoter_tss_overlap
  ))
}

final.loops.from.ctcf.step
final.loops.from.promoter.step
final.loops.from.tss.step

overlapping_loops <- extract_overlapping_loops(
  ctcf_data = final.loops.from.ctcf.step,
  promoter_data = final.loops.from.promoter.step,
  tss_data = final.loops.from.tss.step
)

overlapping_loops$ctcf_promoter_only # 730 :PASS
overlapping_loops$ctcf_tss_only # 3266 :PASS
overlapping_loops$ctcf_promoter_tss_overlap # 8537 :PASS

df.ctcf.promoter.only.loop <- data.frame(loop.id = overlapping_loops$ctcf_promoter_only, category = "CP", stringsAsFactors = FALSE)
df.ctcf.tss.only.loop <- data.frame(loop.id = overlapping_loops$ctcf_tss_only, category = "CT", stringsAsFactors = FALSE)
df.ctcf.promoter.tss.overlap.loop <- data.frame(loop.id = overlapping_loops$ctcf_promoter_tss_overlap, category = "CPT", stringsAsFactors = FALSE)

df.final.loop <- bind_rows(df.ctcf.promoter.only.loop, df.ctcf.tss.only.loop, df.ctcf.promoter.tss.overlap.loop)

df.final.loop # 12493/31773(0.3931955) /58992
save(df.final.loop, file="./figures/submission/df_final_loop.rda")
write.csv(df.final.loop, file = "./figures/submission/df_final_loop.csv", row.names = FALSE)



