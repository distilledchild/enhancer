library("ComplexUpset")
library("tidyverse")
library("GenomicRanges")
library(ggplot2)
#library("devtools")
#library("remotes")
library(gridExtra)
library("patchwork")
#library("cowplot")
#library(biomaRt)
library(reshape2)
library(ggvenn)

#getwd()

# Windows
# setwd('C:\\Users\\panju\\Dropbox (UTHSC GGI)\\Gateway_to_Hao\\workshop\\2023_NIH_meeting\\loop_N_tss')
# source(file.path('C:\\Users\\panju\\Dropbox (UTHSC GGI)\\Gateway_to_Hao\\project_common_code', 'variables.R'))
# source(file.path('C:\\Users\\panju\\Dropbox (UTHSC GGI)\\Gateway_to_Hao\\project_common_code', 'funcs.R'))

#getwd()

# Linux
#setwd('C:\\Users\\panju\\Dropbox (UTHSC GGI)\\Gateway_to_Hao\\workshop\\2023_NIH_meeting\\loop_N_tss')
# setwd('./Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss')
# setwd('~/Desktop/temp/enhancer/dropbox_enhancer_doosan/data/enhancer_atlas2.0/all_species/neuron')
#setwd('~/Desktop/temp/enhancer/dropbox_enhancer_doosan')
#setwd('/home/pkim/dropbox/Gateway_to_Hao/enhancer/r_files')
#source(file.path('/home/pkim/dropbox/Gateway_to_Hao/project_common_code/', 'variables.R'))
#source(file.path('/home/pkim/dropbox/Gateway_to_Hao/project_common_code/', 'funcs.R'))

#getwd()

# Mac
#setwd('/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files')
#getwd()
#source(file.path('/Users/PanjunKim/dropbox/Gateway_to_Hao/project_common_code/', 'variables.R'))
#source(file.path('/Users/PanjunKim/dropbox/Gateway_to_Hao/project_common_code/', 'funcs.R'))

# source(file.path('/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/enhancer_synced_data_preparation.R'))
# source(file.path('/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/data_analysis.R'))
# source(file.path('/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/data_analysis.R'))

##############################################################################################
## LOOP processing: BASIC
##############################################################################################

# BEDPE file list (10 files)
# Linux
# loop.file.list = fs::dir_ls("/home/pkim/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/", regexp = ".bedpe$")
#loop.file.list = fs::dir_ls("/home/pkim/Desktop/temp/enhancer/dropbox_enhancer_doosan/data/loops/", regexp = ".bedpe$")
#loop.file.list = fs::dir_ls("/home/pkim/dropbox/Gateway_to_Hao/enhancer/data/loops", regexp = ".bedpe$")
#loop.file.list

# Mac
#loop.file.list = fs::dir_ls("/Users/PanjunKim/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/", regexp = ".bedpe$")
#loop.file.list = fs::dir_ls("/Users/PanjunKim/Desktop/temp/enhancer/dropbox_enhancer_doosan/data/loops/", regexp = ".bedpe$")
loop.file.list = fs::dir_ls("../data/loops", regexp = ".bedpe$")
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
  mutate(loop.size= (y2+y1)/2 - (x2+x1)/2) |> 
  mutate(resolution = case_when(
    end.distance == 5000 ~ "5K",
    end.distance == 10000 ~ "10K",
    end.distance == 25000 ~ "25K",
    TRUE ~ NA
  )) %>% 
  mutate(Resolution = factor(resolution, levels = c("5K", "10K", "25K"))) %>%
  mutate(loop.id = str_c(X.chr1, '_', x1, '_', x2, '_', chr2, '_', y1, '_', y2, '_', end.distance)) %>% # loop.id
  mutate(sample.loop.id = str_c(strain, '_', loop.id)) #sample.loop.id

names(df.loop.deep.sample.all)
#  [1] "sample"         "X.chr1"         "x1"             "x2"            
#  [5] "chr2"           "y1"             "y2"             "name"          
#  [9] "score"          "strand1"        "strand2"        "color"         
# [13] "observed"       "expectedBL"     "expectedDonut"  "expectedH"     
# [17] "expectedV"      "fdrBL"          "fdrDonut"       "fdrH"          
# [21] "fdrV"           "numCollapsed"   "centroid1"      "centroid2"     
# [25] "radius"         "distance"       "strain"         "end.distance"  
# [29] "loop.size"      "resolution"     "loop.id"        "sample.loop.id"

plot_raw_loop_length_by_rez<-ggplot(data=df.loop.deep.sample.all, aes(x=loop.size/1000000, group=Resolution, fill=Resolution))+geom_density(alpha=.5)+scale_x_log10() + xlab("Mb")
pdf(file="hao_loop_size_before_remove_duplicate_loops.pdf", width=10,height=6)
plot_raw_loop_length_by_rez
dev.off()

## maybe create an upset to find shared vs unique 
## max 10 (because we have 10 samples)
df_count_loopid<- df.loop.deep.sample.all |> count(loop.id , resolution)
head(df_count_loopid)
plot_loopid_counts<-ggplot(df_count_loopid, aes(x=n, group=resolution, fill=resolution))+geom_histogram()
#plot_loopid_counts
n10_loopids<-subset(df_count_loopid, n==10)$loop.id
subset(df.loop.deep.sample.all, loop.id %in% n10_loopids[1])

df.loop.deep.sample.all %>% 
  head()
df.loop.deep.sample.all %>% 
  count(strain)




########################
# loop analysis
########################
# TODO HOLDING
# 1. loop qualification for filtering by observed/expected
df.loop.deep.sample.all <- df.loop.deep.sample.all %>%
  mutate(Q1_value = quantile(observed, 0.25, na.rm = TRUE)) %>%
  group_by(strain, end.distance) %>%
  ungroup()
df.loop.deep.sample.all %>% count(end.distance)

# create the boxplot with grouping by strain and end.distance
plot_resolution_only <- ggplot(df.loop.deep.sample.all, aes(x = resolution, y = observed, fill = resolution)) +
  geom_boxplot(outlier.shape = NA) +  # Remove outliers
  stat_summary(fun.data = function(y) { 
    data.frame(y = quantile(y, 0.25), label = paste0("Q1: ", round(quantile(y, 0.25), 1)))
  }, geom = "text", aes(label = after_stat(label)), vjust = 2.5, size = 3) +  # Add Q1 value as text
  theme_minimal() +
  theme(plot.title = element_text(hjus <- .5)) +
  labs(title = "Observed Values by Resolution", y = "Frequency", x = "Resolution") +
  scale_y_continuous(limits = c(0, 250)) +  # Set y-axis limit to 250
  scale_fill_brewer(palette = "Set3")  # Choose color palette for resolutions

# create the boxplot with grouping by strain and end.distance
plot_resolution_and_strain <- ggplot(df.loop.deep.sample.all, aes(x = interaction(strain, end.distance), y = observed, fill = strain, linetype = resolution)) +
  geom_boxplot(outlier.shape = NA, aes(group = interaction(strain, end.distance))) +  # Remove outliers and differentiate with linetype
  stat_summary(fun.data = function(y) { 
    data.frame(y = quantile(y, 0.25), label = paste0("Q1: ", round(quantile(y, 0.25), 1)))
  }, geom = "text", aes(label = after_stat(label)), vjust = 2.5, size = 1.5) +  # Add Q1 value as text
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Observed Values by Strain and Resolution", y = "Frequency", fill = "Strain", linetype = "Resolution") +
  scale_y_continuous(limits = c(0, 250)) +  # Set y-axis limit to 250
  scale_fill_brewer(palette = "Set3") +  # Choose color palette for strains
  theme(
    axis.text.x = element_blank(),  # Remove x-axis labels
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    strip.text = element_blank(),  # Remove facet labels
    strip.background = element_blank(),
    legend.position = "right"  # Place legend on the right
  ) +
  facet_grid(~ strain + end.distance, scales = "free_x", space = "free_x") +  # Group by strain and end.distance
  theme(axis.title.x = element_blank())  # Remove x-axis title

pdf("loop_quality_control.pdf", width = 10, height = 14)  # Set the PDF output file
grid.arrange(plot_resolution_only, plot_resolution_and_strain, ncol = 1)  # Arrange plots vertically
dev.off()  # Close the PDF device
########################################### ongoing

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
      # way1 - percent: common loops / ave of two samples
      # percentage <- (common / mean(c(length(location_list[[i]]), length(location_list[[j]])))) * 100
      # way2 - percent: common loops / sample loop in X-axis 
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
          plot.title = element_text(hjust = 0.5))  # 제목을 가운데 정렬
  
  # total loops
  loop_counts_table <- loop_counts %>% 
    mutate(Strain = factor(strain, levels = names(location_list))) %>% 
    arrange(Strain)
  
  colnames(common_counts) <- paste0(colnames(common_counts), "\n(", loop_counts$total_loops, ")")
  rownames(common_counts) <- paste0(rownames(common_counts), "\n(", loop_counts$total_loops, ")")
  
  plots[[res]] <- p
}

# all figures by resolution into one figure
pdf("common_Loops_Heatmap_Percentage.pdf", width = 12, height = 12)
grid.arrange(plots[["5K"]], plots[["10K"]], plots[["25K"]], nrow = 3)
dev.off()
  
# loop.setting with GRanges: x1,x2 = the coordinates of the UPSTREAM | y1,y2 = the coordinates of the DOWNSTREAM

########### distinct loops
df.loop.deep.sample.all %>% 
  count() # 58992
#  head()
  

df.DISTINCT.loop.deep.sample.all <- df.loop.deep.sample.all %>%
  dplyr::select(loop.id) %>%
  distinct() %>% 
  separate(loop.id, into = c("chr1", "x1", "x2", "chr2", "y1", "y2", "end.distance"), sep = "_", remove = FALSE) %>% 
  mutate(
    x1 = as.numeric(x1),
    x2 = as.numeric(x2),
    y1 = as.numeric(y1),
    y2 = as.numeric(y2),
    end.distance = as.numeric(end.distance)
  ) %>% 
  mutate(distance = y2 - x2) %>%
  mutate(resolution = case_when(
    end.distance == 5000 ~ "5K",
    end.distance == 10000 ~ "10K",
    end.distance == 25000 ~ "25K",
    TRUE ~ NA
  )) %>% 
  mutate(Resolution = factor(resolution, levels = c("5K", "10K", "25K"))) %>% 
  mutate(padded.distance = 0) %>% 
  mutate(loop.id = str_c(loop.id, '_', padded.distance))

df.DISTINCT.loop.deep.sample.all$loop.id

df.DISTINCT.loop.deep.sample.all %>% 
  count(resolution) # 31773/58992
  head()
# resolution     n
# 1         5K  6680
# 2        10K 12162
# 3        25K 12931
dim(df.DISTINCT.loop.deep.sample.all)
# [1] 31773    11
  

df.DISTINCT.loop.deep.sample.all <- df.DISTINCT.loop.deep.sample.all |> 
  mutate(loop.size= (y2+y1)/2 - (x2+x1)/2)  

plot_distinct_loop_length_by_res<-ggplot(data=df.DISTINCT.loop.deep.sample.all, aes(x=loop.size/1000000, group=Resolution, fill=Resolution))+geom_density(alpha=.5)+scale_x_log10() + xlab("Mb")+ggtitle("Size of loops")+theme(legend.position="none")
plot_num_loops_by_res<-ggplot(df.DISTINCT.loop.deep.sample.all, aes(x=Resolution, fill=Resolution))+geom_histogram(stat="count")+ggtitle("Number of loops")
pdf(file="hao_loop_size_after_remove_duplicate_loops.pdf", width=8,height=4)
plot_distinct_loop_length_by_res+plot_num_loops_by_res
dev.off()



########### GRANGE with object from line 81 
# loop.setting with GRanges - loop.deep.sample.all (whole, up, down)
df.loop.deep.sample.all %>% head(2)
t = 0
loop.deep.sample.all.whole.GR<-GRanges(seqnames=df.loop.deep.sample.all$X.chr1, ranges=IRanges(start=(df.loop.deep.sample.all$x1-t), end=(df.loop.deep.sample.all$y2+t)), id=df.loop.deep.sample.all$loop.id)
loop.deep.sample.all.up.GR<-GRanges(seqnames=df.loop.deep.sample.all$X.chr1, ranges=IRanges(start=(df.loop.deep.sample.all$x1-t), end=(df.loop.deep.sample.all$x2+t)), id=df.loop.deep.sample.all$loop.id)
loop.deep.sample.all.down.GR<-GRanges(seqnames=df.loop.deep.sample.all$X.chr1, ranges=IRanges(start=(df.loop.deep.sample.all$y1-t), end=(df.loop.deep.sample.all$y2+t)), id=df.loop.deep.sample.all$loop.id)
loop.deep.sample.all.middle.GR<-GRanges(seqnames=df.loop.deep.sample.all$X.chr1, ranges=IRanges(start=(df.loop.deep.sample.all$x2-t), end=(df.loop.deep.sample.all$y1+t)), id=df.loop.deep.sample.all$loop.id)

# DISTINCT
t = 0
df.DISTINCT.loop.deep.sample.all.whole.GR<-GRanges(seqnames=df.DISTINCT.loop.deep.sample.all$chr1, ranges=IRanges(start=(df.DISTINCT.loop.deep.sample.all$x1-t), end=(df.DISTINCT.loop.deep.sample.all$y2+t)), id=df.DISTINCT.loop.deep.sample.all$loop.id, end.distance = df.DISTINCT.loop.deep.sample.all$end.distance, resolution = df.DISTINCT.loop.deep.sample.all$resolution, distance = df.DISTINCT.loop.deep.sample.all$distance)
df.DISTINCT.loop.deep.sample.all.up.GR<-GRanges(seqnames=df.DISTINCT.loop.deep.sample.all$chr1, ranges=IRanges(start=(df.DISTINCT.loop.deep.sample.all$x1-t), end=(df.DISTINCT.loop.deep.sample.all$x2+t)), id=df.DISTINCT.loop.deep.sample.all$loop.id, end.distance = df.DISTINCT.loop.deep.sample.all$end.distance, resolution = df.DISTINCT.loop.deep.sample.all$resolution, distance = df.DISTINCT.loop.deep.sample.all$distance)
df.DISTINCT.loop.deep.sample.all.down.GR<-GRanges(seqnames=df.DISTINCT.loop.deep.sample.all$chr1, ranges=IRanges(start=(df.DISTINCT.loop.deep.sample.all$y1-t), end=(df.DISTINCT.loop.deep.sample.all$y2+t)), id=df.DISTINCT.loop.deep.sample.all$loop.id, end.distance = df.DISTINCT.loop.deep.sample.all$end.distance, resolution = df.DISTINCT.loop.deep.sample.all$resolution, distance = df.DISTINCT.loop.deep.sample.all$distance)
df.DISTINCT.loop.deep.sample.all.middle.GR<-GRanges(seqnames=df.DISTINCT.loop.deep.sample.all$chr1, ranges=IRanges(start=(df.DISTINCT.loop.deep.sample.all$x2-t), end=(df.DISTINCT.loop.deep.sample.all$y1+t)), id=df.DISTINCT.loop.deep.sample.all$loop.id, end.distance = df.DISTINCT.loop.deep.sample.all$end.distance, resolution = df.DISTINCT.loop.deep.sample.all$resolution, distance = df.DISTINCT.loop.deep.sample.all$distance)


########### NEW OBJECT for OVERALL DISTRIBUTION: new.df.loop.deep.sample.all, overall.df.DISTINCT.loop.deep.sample.all
# new.loop.id, x0, y3, new_distance
# new.df.loop.deep.sample.all <- df.loop.deep.sample.all %>% 
overall.df.DISTINCT.loop.deep.sample.all.5.distance <- df.DISTINCT.loop.deep.sample.all %>% 
  mutate(x12 = (x1 + x2)/2, y12 = (y1 + y2)/2) %>% # middle point of each end
  # mutate(x0 = ifelse(x12 - (distance) < 0, 0, x12 - (distance)), y3 = y12 + (distance)) %>% # a distance for padding
  mutate(x0 = ifelse(x12 - (distance / 2) < 0, 0, x12 - (distance / 2)), y3 = y12 + (distance/2)) %>% # 1/2 distance for padding
  # mutate(new_distance = abs(y3 - x0)) %>% 
  mutate(new.loop.id = str_c(loop.id, '|', x0, '_', y3))

overall.df.DISTINCT.loop.deep.sample.all.1.distance <- df.DISTINCT.loop.deep.sample.all %>% 
  mutate(x12 = (x1 + x2)/2, y12 = (y1 + y2)/2) %>% # middle point of each end
  mutate(x0 = ifelse(x12 - (distance) < 0, 0, x12 - (distance)), y3 = y12 + (distance)) %>% # a distance for padding
  # mutate(x0 = ifelse(x12 - (distance / 2) < 0, 0, x12 - (distance / 2)), y3 = y12 + (distance/2)) %>% # 1/2 distance for padding
  # mutate(new_distance = abs(y3 - x0)) %>% 
  mutate(new.loop.id = str_c(loop.id, '|', x0, '_', y3))

# new.df.loop.deep.sample.all %>%
# overall.df.DISTINCT.loop.deep.sample.all.5.distance %>% 
overall.df.DISTINCT.loop.deep.sample.all.1.distance %>% 
  # filter(x0 < 0) %>%   
  head()

# data integrity: PASS
# overall.df.DISTINCT.loop.deep.sample.all %>% mutate(test = ifelse(y3-x0 == 2*distance, TRUE, FALSE)) %>% 
#   filter(test == FALSE) %>% 
#   count(x0)
# count(test)

#### overall.df.DISTINCT.loop.deep.sample.all : OBJECT to be used for the distribution of sth over loops

####################################################
## DISTRIBUTION: CTCF distribution
####################################################

## 1. CTCF Processing

# Linux
#df.init.ctcf<-read.table(file="/home/pkim/Desktop/temp/enhancer/dropbox_enhancer_doosan/data/ctcf/fimo_2nd_trial_meme/fimo_2nd_trial.txt", header=TRUE, sep="\t")
#df.init.ctcf<-read.table(file="/home/pkim/dropbox/Gateway_to_Hao/enhancer/data/ctcf/fimo_2nd_trial_meme/fimo_2nd_trial.txt", header=TRUE, sep="\t")
#df.init.ctcf %>% count() # 5767921

# Mac
# ctcf<-read.table(file="/home/pkim/Desktop/temp/enhancer/dropbox_enhancer_doosan/data/ctcf/fimo_2nd_trial_meme/fimo_2nd_trial.txt", sep="\t", col.names = c("chr", "start", "end", "strand", "length"), header = FALSE) %>% 
#df.init.ctcf<-read.table(file="/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/data/ctcf/fimo_2nd_trial_meme/fimo_2nd_trial.txt", header=TRUE, sep="\t")
df.init.ctcf<-read.table(file="../data/ctcf/fimo_2nd_trial_meme/fimo_2nd_trial.txt", header=TRUE, sep="\t")
df.init.ctcf %>% dim() # 5767921
df.init.ctcf %>% head()

# data exploration: CTCF
df.init.ctcf %>%                              # +: 2891072, -: 2876849 = 5767921
  count(strand)
# removing dups including strand
df.init.ctcf %>% 
  distinct(chr, start, end, strand)  |> count()         # 2701585/5767921
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

# generating id column (long running time)
df.DISTINCT.fimo.2nd.trial.ctcf <- df.init.ctcf %>%
  distinct() %>% 
  mutate(id = str_c(chr, "_", start, "_", end, "_", strand, '_', length)) %>% 
  mutate(start = as.numeric(start)) %>% 
  mutate(end = as.numeric(end)) 
dim(df.DISTINCT.fimo.2nd.trial.ctcf)
# [1] 2701585       6
df.DISTINCT.fimo.2nd.trial.ctcf # 2701585/5767921 : 0.4683811
df.DISTINCT.fimo.2nd.trial.ctcf %>% head()  # id column done


ctcf_chr4<-subset(df.DISTINCT.fimo.2nd.trial.ctcf, chr=="chr4")
ctcf_chr19<-subset(df.DISTINCT.fimo.2nd.trial.ctcf, chr=="chr19")
ctcf_chr5<-subset(df.DISTINCT.fimo.2nd.trial.ctcf, chr=="chr5")
write.csv(file="chr4_ctcf.tab", ctcf_chr4, row.names=F)
write.csv(file="chr19_ctcf.tab", ctcf_chr19, row.names=F)
write.csv(file="chr5_ctcf.tab", ctcf_chr5, row.names=F)




# GRANGE for ctcf from FIMO 2nd trial, 17849
df.DISTINCT.ctcf.2nd.fimo.GR <- GRanges(seqnames=df.DISTINCT.fimo.2nd.trial.ctcf$chr, ranges=IRanges(start=df.DISTINCT.fimo.2nd.trial.ctcf$start, end=df.DISTINCT.fimo.2nd.trial.ctcf$end), id=df.DISTINCT.fimo.2nd.trial.ctcf$id, strand=df.DISTINCT.fimo.2nd.trial.ctcf$strand)
#df.DISTINCT.ctcf.2nd.fimo.G5R

####################################################
##### 1. OVERALL on loops Overlapping CTCF & loops 
####################################################

overall.df.DISTINCT.loop.deep.sample.all <- overall.df.DISTINCT.loop.deep.sample.all.5.distance

overall.df.DISTINCT.loop.deep.sample.all.GR <- GRanges(seqnames=overall.df.DISTINCT.loop.deep.sample.all$chr1, ranges=IRanges(start=overall.df.DISTINCT.loop.deep.sample.all$x0, end=overall.df.DISTINCT.loop.deep.sample.all$y3), id=overall.df.DISTINCT.loop.deep.sample.all$loop.id, new.loop.id=overall.df.DISTINCT.loop.deep.sample.all$new.loop.id)

index.distinct.ctcf.w.overall.whole.loop <- findOverlaps(df.DISTINCT.ctcf.2nd.fimo.GR, overall.df.DISTINCT.loop.deep.sample.all.GR, type = "within")

overall.loop.for.ctcf.hits <- subjectHits(index.distinct.ctcf.w.overall.whole.loop)
overall.ctcf.on.loop.hits <- queryHits(index.distinct.ctcf.w.overall.whole.loop)

df.ctcf.dist.result <- data.frame(
  loop.id = overall.df.DISTINCT.loop.deep.sample.all$new.loop.id[overall.loop.for.ctcf.hits],
  loop.start0 = overall.df.DISTINCT.loop.deep.sample.all$x0[overall.loop.for.ctcf.hits], # modified by hao
  loop.start12 = overall.df.DISTINCT.loop.deep.sample.all$x12[overall.loop.for.ctcf.hits], # modified by hao
  loop.end3 = overall.df.DISTINCT.loop.deep.sample.all$y3[overall.loop.for.ctcf.hits], # modified by hao
  loop.end12 = overall.df.DISTINCT.loop.deep.sample.all$y12[overall.loop.for.ctcf.hits], # modified by hao
  loop.res = overall.df.DISTINCT.loop.deep.sample.all$resolution[overall.loop.for.ctcf.hits],
  # loop.new.distance = overall.df.DISTINCT.loop.deep.sample.all$new_distance[overall.loop.for.ctcf.hits],
  ctcf.id = df.DISTINCT.fimo.2nd.trial.ctcf$id[overall.ctcf.on.loop.hits],
  ctcf.start = df.DISTINCT.fimo.2nd.trial.ctcf$start[overall.ctcf.on.loop.hits],
  ctcf.end = df.DISTINCT.fimo.2nd.trial.ctcf$end[overall.ctcf.on.loop.hits]
) 

df.ctcf.dist.result %>% head()
df.ctcf.dist.result %>% dim() # 40861258(both DISTINCT, 0.5 distance), 6/57612728(both DISTINCT, 1 distance)

## Start By Hao
relative.pos.df.ctcf.dist.result <- df.ctcf.dist.result %>%
  mutate(pos_coord = round((ctcf.start + ctcf.end) / 2)) %>%
  mutate(loop_length03 = (loop.end3 - loop.start0)) %>% # x0, y3
  mutate(loop_length = (loop.end12 - loop.start12)) %>% # x12, y12
  mutate(relative_pos03 = pos_coord - loop.start0) %>% # relative position from loop start
  mutate(relative_pos = pos_coord - loop.start12) %>% # relative position from loop start
  mutate(value03 = (relative_pos03 / (loop_length03 / 2)) - 0.5) %>% # for padding .5x distance
  mutate(value = relative_pos / loop_length - 0.5 ) %>% # for padding .5x distance
  # mutate(value = relative_pos/(loop_length/3)-1) %>% # for padding 1x distance
  dplyr::select(loop.id, ctcf.id, value, value03, loop.res)

## End By Hao

relative.pos.df.ctcf.dist.result %>% 
  head()

## PDF file by CHROMOSOME

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
    labs(title = paste0("A. Density of CTCF found near Hi-C loops on chr", chr),
         x = "Relative position to loop",
         y = "Density"
    )
  
  plot1.ctcf.hist <- relative.pos.df.ctcf.dist.result.chr %>% 
    ggplot(aes(x = value)) +
    geom_histogram(fill = "skyblue", color = "black", alpha = 0.7, bins=200) +
    labs(title = paste0("A. Histogram of CTCF found near Hi-C loops on chr", chr),
         x = "Relative position to loop",
         y = "Counts"
    )

  combined_plot <- plot1.ctcf.hist + plot1.ctcf.dens # + can be used instead of |
  plot_ctcf_list[[chr]] <- combined_plot

  message("END: Processing chromosome: ", chr)
  }, error = function(e) {
    
    message("Error processing chromosome: ", chr)
    message("Error message: ", e$message)
  })
}

pdf("overall_distribution_of_CTCF_on_each_chr_.5_distance.pdf", width = 11, height = 8.5)
#pdf("overall_distribution_of_CTCF_on_each_chr_1_distance.pdf", width = 11, height = 8.5)

num_plots <- length(plot_ctcf_list) # 22 chromosomes
plots_per_page <- 4

for (i in seq(1, num_plots, by = plots_per_page)) {
  end_idx <- min(i + plots_per_page - 1, num_plots)
  page_plots <- plot_ctcf_list[i:end_idx]
  
  combined_page <- plot_grid(plotlist = page_plots, ncol = 1, nrow = 4)
  
  print(combined_page)
}

dev.off()

relative.pos.df.ctcf.dist.result <- 
relative.pos.df.ctcf.dist.result |> 
  mutate(Resolution = factor(loop.res, levels=c("5K", "10K", "25K")))



## Start By Hao: 
plot.ctcf.hist_03 <- relative.pos.df.ctcf.dist.result %>% 
  ggplot(aes(x = value03)) +
  geom_histogram(fill = "skyblue", color = "black", alpha = 0.7, bins=200) +
  labs(title = paste0("A. value03, Histogram of CTCF over Loops on ALL Chromosomes" ),
       x = "Relative position to loop",
       y = "Counts"
  )

plot.ctcf.hist_12 <- relative.pos.df.ctcf.dist.result %>% 
  ggplot(aes(x = value)) +
  geom_histogram(fill = "skyblue", color = "black", alpha = 0.7, bins=200) +
  labs(title = paste0("A. value12, Histogram of CTCF over Loops on ALL Chromosomes" ),
       x = "Relative position to loop",
       y = "Counts"
  )
                  names(relative.pos.df.ctcf.dist.result)
# [1] "loop.id"  "ctcf.id"  "value03"  "value12" 
# [5] "loop.res"
plot.ctcf.density_12 <- relative.pos.df.ctcf.dist.result %>% 
  ggplot(aes(x = value , group=Resolution, fill=Resolution)) +
  geom_density(alpha = 0.5) +
  labs(title = paste0("Density of CTCF over Loops on ALL Chromosomes" ),
       x = "Relative position to loop",
       y = "Counts"
  )
plot.ctcf.hist_12
pdf(file="hao.ctcf.hist.03.vs12.pdf", width=12, height=8) 
plot.ctcf.hist_03+plot.ctcf.hist_12
plot.ctcf.hist_12+facet_wrap(~loop.res)
dev.off()
## End By Hao 


## PDF file for ALL CHROMOSOMES by resolution
create_ctcf_plots <- function(data, res) {
  if (res == "All" ) {
    df.plot<-data
  } else { 
    df.plot<- data %>% filter(loop.res==res)
  }

  plot.ctcf.hist <- df.plot %>% 
  ggplot(aes(x = value)) +
  geom_histogram(fill = "skyblue", color = "black", alpha = 0.7, bins=200) +
  labs(title = paste0("A. Histogram of CTCF over Loops on ALL Chromosomes (", res, ")"),
       x = "Relative position to loop",
       y = "Counts"
  )
  
  plot.ctcf.dens <- df.plot %>% 
    ggplot(aes(x = value)) +
    geom_density(fill = "skyblue", color = "black", alpha = 0.7) +
    ylim(c(0,1))+
    labs(title = paste0("A. Density of CTCF found near Hi-C loops on ALL chromosomes"),
         x = "Relative position to loop",
         y = "Density"
    )

    return(plot.ctcf.hist + plot.ctcf.dens)
}

plot_ctcf_all_res <- create_ctcf_plots(relative.pos.df.ctcf.dist.result, "All")
plot_ctcf_5K <- create_ctcf_plots(relative.pos.df.ctcf.dist.result, "5K")
plot_ctcf_10K <- create_ctcf_plots(relative.pos.df.ctcf.dist.result, "10K")
plot_ctcf_25K <- create_ctcf_plots(relative.pos.df.ctcf.dist.result, "25K")
pdf("overall_distribution_ctcf_on_all_chromosomes_.5_distance.pdf", width = 11, height = 8.5) # pdf for .5 distance as a padding
# pdf("overall_distribution_ctcf_on_all_chromosomes_1_distance.pdf", width = 11, height = 8.5) # pdf for a distance as a padding
(plot_ctcf_all_res/plot_ctcf_5K)
(plot_ctcf_10K / plot_ctcf_25K)
dev.off()

################################################################
# 2. EACH END, boxplots
# CTCF at each end in a loop: for the number of CTCF used in filtering valid loops
# so, the object should be used one without padding
# df.DISTINCT.loop.deep.sample.all
################################################################

# 1. inner distance analysis 1: checking inner distance between x3 and y0
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
          border = "black")
  
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

pdf("./figures/inner_distance_boxplots.pdf", width = 16, height = 8)
par(mfrow = c(1, 2))  

# 50% padding boxplot
generate_boxplot(df_50_padded, "Boxplot of Inner Distance (50% of End Size Padding)")

# 100% padding boxplot
generate_boxplot(df_100_padded, "Boxplot of Inner Distance (100% of End Size Padding)")

dev.off()

################################################################################################

df_50_padded %>% head(2)
df_100_padded %>% head(2)
df.DISTINCT.loop.deep.sample.all %>% head(2)

df.DISTINCT.loop.deep.sample.all

df.DISTINCT.loop.deep.sample.all.whole.GR<-GRanges(seqnames=df.DISTINCT.loop.deep.sample.all$chr1, ranges=IRanges(start=(df.DISTINCT.loop.deep.sample.all$x1-t), end=(df.DISTINCT.loop.deep.sample.all$y2+t)), id=df.DISTINCT.loop.deep.sample.all$loop.id, end.distance = df.DISTINCT.loop.deep.sample.all$end.distance, resolution = df.DISTINCT.loop.deep.sample.all$resolution, distance = df.DISTINCT.loop.deep.sample.all$distance)
df.DISTINCT.loop.deep.sample.all.up.GR<-GRanges(seqnames=df.DISTINCT.loop.deep.sample.all$chr1, ranges=IRanges(start=(df.DISTINCT.loop.deep.sample.all$x1-t), end=(df.DISTINCT.loop.deep.sample.all$x2+t)), id=df.DISTINCT.loop.deep.sample.all$loop.id, end.distance = df.DISTINCT.loop.deep.sample.all$end.distance, resolution = df.DISTINCT.loop.deep.sample.all$resolution, distance = df.DISTINCT.loop.deep.sample.all$distance)
df.DISTINCT.loop.deep.sample.all.down.GR<-GRanges(seqnames=df.DISTINCT.loop.deep.sample.all$chr1, ranges=IRanges(start=(df.DISTINCT.loop.deep.sample.all$y1-t), end=(df.DISTINCT.loop.deep.sample.all$y2+t)), id=df.DISTINCT.loop.deep.sample.all$loop.id, end.distance = df.DISTINCT.loop.deep.sample.all$end.distance, resolution = df.DISTINCT.loop.deep.sample.all$resolution, distance = df.DISTINCT.loop.deep.sample.all$distance)
df.DISTINCT.loop.deep.sample.all.middle.GR<-GRanges(seqnames=df.DISTINCT.loop.deep.sample.all$chr1, ranges=IRanges(start=(df.DISTINCT.loop.deep.sample.all$x2-t), end=(df.DISTINCT.loop.deep.sample.all$y1+t)), id=df.DISTINCT.loop.deep.sample.all$loop.id, end.distance = df.DISTINCT.loop.deep.sample.all$end.distance, resolution = df.DISTINCT.loop.deep.sample.all$resolution, distance = df.DISTINCT.loop.deep.sample.all$distance)

df.DISTINCT.loop.deep.sample.all.up.GR

# 1. CTCF + UPSTREAM (df.DISTINCT.loop.deep.sample.all, df.DISTINCT.fimo.2nd.trial.ctcf)
index.distinct.ctcf.w.up.loop <- findOverlaps(df.DISTINCT.ctcf.2nd.fimo.GR, df.DISTINCT.loop.deep.sample.all.up.GR, type = "within")
index.distinct.ctcf.w.up.loop
end.loop.up.ctcf.hits <- subjectHits(index.distinct.ctcf.w.up.loop)
end.ctcf.up.hits <- queryHits(index.distinct.ctcf.w.up.loop)
end.loop.up.ctcf.hits

df.DISTINCT.loop.deep.sample.all %>% 
  head()
df.DISTINCT.fimo.2nd.trial.ctcf %>% 
  head()

df.overlapping.CTCF.w.UPSTREAM.result <- data.frame(
  up.loop.id = mcols(df.DISTINCT.loop.deep.sample.all.up.GR)$id[end.loop.up.ctcf.hits],
  end.up.distance = mcols(df.DISTINCT.loop.deep.sample.all.up.GR)$end.distance[end.loop.up.ctcf.hits],
  distance = mcols(df.DISTINCT.loop.deep.sample.all.up.GR)$distance[end.loop.up.ctcf.hits],
  resolution = mcols(df.DISTINCT.loop.deep.sample.all.up.GR)$resolution[end.loop.up.ctcf.hits],
  ctcf.id = mcols(df.DISTINCT.ctcf.2nd.fimo.GR)$id[end.ctcf.up.hits],
  WHERE = "UP"
) %>% 
  mutate(ctcf.loop.up.id = str_c(up.loop.id, '|', distance, '|', ctcf.id, '|', WHERE))

df.overlapping.CTCF.w.UPSTREAM.result %>% dim() # 1206765
df.overlapping.CTCF.w.UPSTREAM.result %>% head()

df.overlapping.CTCF.w.UPSTREAM.result %>% 
  count(end.up.distance)

# end.up.distance      n
# 1            5000 154950
# 2           10000 383562
# 3           25000 668253
# total: 1206765

# 2. CTCF + DOWNSTREAM
index.distinct.ctcf.w.down.loop <- findOverlaps(df.DISTINCT.ctcf.2nd.fimo.GR, df.DISTINCT.loop.deep.sample.all.down.GR, type = "within")
df.DISTINCT.loop.deep.sample.all.down.GR
end.loop.down.ctcf.hits <- subjectHits(index.distinct.ctcf.w.down.loop)
end.ctcf.down.hits <- queryHits(index.distinct.ctcf.w.down.loop)

df.overlapping.CTCF.w.DOWNSTREAM.result <- data.frame(
  down.loop.id = mcols(df.DISTINCT.loop.deep.sample.all.down.GR)$id[end.loop.down.ctcf.hits],
  end.down.distance = mcols(df.DISTINCT.loop.deep.sample.all.down.GR)$end.distance[end.loop.down.ctcf.hits],
  distance = mcols(df.DISTINCT.loop.deep.sample.all.down.GR)$distance[end.loop.down.ctcf.hits],
  resolution = mcols(df.DISTINCT.loop.deep.sample.all.down.GR)$resolution[end.loop.down.ctcf.hits],
  ctcf.id = mcols(df.DISTINCT.ctcf.2nd.fimo.GR)$id[end.ctcf.down.hits],
  WHERE = "DOWN"
) %>% 
  mutate(ctcf.loop.down.id = str_c(down.loop.id, '|', distance, '|', ctcf.id, '|', WHERE))

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

df.overlapping.CTCF.w.BOTH.result.boxplot <- df.overlapping.CTCF.w.BOTH.result %>% 
  group_by(chr, loop.id, WHERE, resolution) %>%
  # group_by(loop.id, WHERE, resolution) %>%
  # group_by(chr, loop.id, resolution) %>%
  summarise(ctcf_count_by_loop_id = n_distinct(ctcf.id), .groups = 'drop')
df.overlapping.CTCF.w.BOTH.result.boxplot

# for Q1
df.overlapping.CTCF.w.BOTH.result %>% head()

# no chr, the rest is the same
df.ctcf.counts <- df.overlapping.CTCF.w.BOTH.result %>%
  group_by(loop.id, WHERE, resolution) %>%
  summarise(ctcf_count = n_distinct(ctcf.id), .groups = 'drop')

#subset(df.ctcf.counts, loop.id=="chr1_101070000_101080000_chr1_101560000_101570000_10000_0")
# # A tibble: 2 × 4
#   loop.id                                                   WHERE resolution ctcf_count
#   <chr>                                                     <fct> <fct>           <int>
# 1 chr1_101070000_101080000_chr1_101560000_101570000_10000_0 UP    10K                23
# 2 chr1_101070000_101080000_chr1_101560000_101570000_10000_0 DOWN  10K                50


ctcf_stats_by_resolution <- df.ctcf.counts %>%
  group_by(resolution) %>%
  # group_by(WHERE, resolution) %>%
  summarise(
    Q1 = quantile(ctcf_count, 0.25, na.rm = TRUE),
    Median = median(ctcf_count, na.rm = TRUE),
    Q3 = quantile(ctcf_count, 0.75, na.rm = TRUE),
    Mean = mean(ctcf_count, na.rm = TRUE),
    SD = sd(ctcf_count, na.rm = TRUE),
    .groups = 'drop'
  )

ctcf_stats_by_resolution

# with resolution
# resolution    Q1 Median    Q3  Mean    SD
# 1 5K            22     40    76  55.2  50.0
# 2 10K           23     43    78  56.7  50.3
# 3 25K           26     49    92  66.5  60.7

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

# quantile for EACH END// NOT in a loop : results are loops have ctcfs more than 18 at least in each end
q3_ctcf_count <- quantile(df.overlapping.CTCF.w.BOTH.result.boxplot$ctcf_count_by_loop_id, 0.75, na.rm = TRUE)
q3_ctcf_count # 54 :: integrity PASS

q1_ctcf_count <- quantile(df.overlapping.CTCF.w.BOTH.result.boxplot$ctcf_count_by_loop_id, 0.25, na.rm = TRUE)
q1_ctcf_count # 16 :: integrity PASS

# drawing boxplot
df.overlapping.CTCF.w.BOTH.result.boxplot
df.overlapping.CTCF.w.BOTH.result.boxplot %>% head(4)

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
    labs(title = "Boxplot for Number of CTCF by Chromosome and Resolution", x = "Chromosomes", y = "Number of CTCF in a loop")

boxplot.w.CTCF.by.chr.and.res

pdf("end_boxplot_num_ctcf_by_chr_res.pdf", width = 16.5, height = 23.5)
grid.arrange(boxplot.w.CTCF.by.chr.and.res, ncol = 1)
dev.off()
pdf("hao_end_boxplot_num_ctcf_by_chr_res.pdf", width=12, height=8)
boxplot.w.CTCF.by.chr.and.res+ylim(c(0,200))
dev.off()


# this is only for one end
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
  labs(title = "", x = "Resolution", y = "Number of CTCF in a loop")

boxplot.w.CTCF.by.res

pdf("hao_end_boxplot_num_ctcf_by_res.pdf", width = 8, height = 6)
grid.arrange(boxplot.w.CTCF.by.res+ylim(c(0,200)), ncol = 1)
dev.off()

pdf("end_boxplot_num_ctcf_combined.pdf", width = 16.5, height = 23.5)
grid.arrange(boxplot.w.CTCF.by.chr.and.res, boxplot.w.CTCF.by.res, ncol = 1)
dev.off()


# figure by end
boxplot.w.CTCF.ALL.chr <- ggplot(df.overlapping.CTCF.w.BOTH.result.boxplot.ALL.chr, aes(x = WHERE, y = ctcf_count_by_loop_id, fill = WHERE)) +
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
  scale_y_continuous(breaks = seq(min(df.overlapping.CTCF.w.BOTH.result.boxplot.ALL.chr$ctcf_count_by_loop_id), max(df.overlapping.CTCF.w.BOTH.result.boxplot.ALL.chr$ctcf_count_by_loop_id), by = 5)) +
  labs(title = "Boxplot by UP/DOWNSTREAM End", x = "Upstream and Downstream End", y = "Number of CTCF inside each ends in a loop")

boxplot.w.CTCF.ALL.chr

# figure by end & res
boxplot.w.CTCF.ALL.chr.by.res <- ggplot(df.overlapping.CTCF.w.BOTH.result.boxplot.ALL.chr, aes(x = WHERE, y = ctcf_count_by_loop_id, fill = resolution)) +
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
  scale_fill_manual(values = c("5K" = "skyblue", "10K" = "lightcoral", "25K" = "springgreen"), breaks = c("5K", "10K", "25K")) +
  scale_y_continuous(breaks = seq(min(df.overlapping.CTCF.w.BOTH.result.boxplot.ALL.chr$ctcf_count_by_loop_id), max(df.overlapping.CTCF.w.BOTH.result.boxplot.ALL.chr$ctcf_count_by_loop_id), by = 5)) +
  labs(title = "Boxplot for Number of CTCF in UP/DOWNSTREAM End by Resolution", x = "Upstream and Downstream End", y = "Number of CTCF", fill = "Resolution")
  # facet_wrap(~ resolution, scales = "free_y")

boxplot.w.CTCF.ALL.chr.by.res

pdf("end_boxplot_no.ctcf_by_end_res.pdf", width = 16.5, height = 23.5)
grid.arrange(boxplot.w.CTCF.ALL.chr, boxplot.w.CTCF.ALL.chr.by.res, ncol = 1)
dev.off()

##########################################################################################
# refactoring code
##########################################################################################

df.DISTINCT.loop.deep.sample.all.50.padded <- df_50_padded %>% 
  mutate(padded.distance = 50) %>% 
  mutate(x1 = x0, x2 = x3, y1 = y0, y2 = y3) #%>% 
#  mutate(loop.id = str_c(loop.id, '_', 50))

df.DISTINCT.loop.deep.sample.all.100.padded <- df_100_padded %>% 
  mutate(padded.distance = 100) %>% 
  mutate(x1 = x0, x2 = x3, y1 = y0, y2 = y3)# %>% 
#  mutate(loop.id = str_c(loop.id, '_', 100))

# analysis func
perform_analysis <- function(loop_data, ctcf_data, padding_label) {
  
  # 1. GRanges objec
  loop_whole_GR <- GRanges(seqnames=loop_data$chr1, 
                           ranges=IRanges(start=(loop_data$x1), end=(loop_data$y2)), 
                           id=loop_data$loop.id, 
                           end.distance = loop_data$end.distance, 
                           resolution = loop_data$resolution, 
                           distance = loop_data$distance)
  
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
  
  loop_middle_GR <- GRanges(seqnames=loop_data$chr1, 
                            ranges=IRanges(start=(loop_data$x2), end=(loop_data$y1)), 
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
head(df_ctcf_boxplot_100)$loop.id
# [1] "chr1_101070000_101080000_chr1_101560000_101570000_10000_0"
# [2] "chr1_101070000_101080000_chr1_101560000_101570000_10000_0"
# [3] "chr1_101075000_101100000_chr1_101550000_101575000_25000_0"
# [4] "chr1_101075000_101100000_chr1_101550000_101575000_25000_0"
# [5] "chr1_101200000_101210000_chr1_101390000_101400000_10000_0"
# [6] "chr1_101200000_101210000_chr1_101390000_101400000_10000_0"

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
pdf("./figures/end_boxplot_num_ctcf_combined_no_padding.pdf", width = 16.5, height = 23.5)
grid.arrange(ctcf_boxplots_none$boxplot_by_chr_res, ctcf_boxplots_none$boxplot_by_res, ncol = 1)
dev.off()

pdf("./figures/end_boxplot_num_ctcf_combined_50_padding.pdf", width = 16.5, height = 23.5)
grid.arrange(ctcf_boxplots_50$boxplot_by_chr_res, ctcf_boxplots_50$boxplot_by_res, ncol = 1)
dev.off()

pdf("./figures/end_boxplot_num_ctcf_combined_100_padding.pdf", width = 16.5, height = 23.5)
grid.arrange(ctcf_boxplots_100$boxplot_by_chr_res, ctcf_boxplots_100$boxplot_by_res, ncol = 1)
dev.off()

pdf("./figures/end_boxplot_num_ctcf_combined_all.pdf", width = 16.5, height = 23.5)
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


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
# DISTRIBUTION: TSS distribution on loops
##########################################################################################
##########################################################################################

##############################################################################################
## TSS processing
##############################################################################################

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
tss<-read.table(file="../../workshop/2023_NIH_meeting/loop_N_tss/ucsc_start_codon.txt", sep="\t", head=F)
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

df.tss.ucsc.GR <-GRanges(seqnames=df.tss.ucsc$chr, ranges=IRanges(start=as.numeric(df.tss.ucsc$start), end=as.numeric(df.tss.ucsc$end)), strand=df.tss.ucsc$strand, gene_id=df.tss.ucsc$gene_id, transcript_id=df.tss.ucsc$transcript_id, tss.id = df.tss.ucsc$tss.id)

###### legacy - start
# filter(gene_id != gene_name)

# tss$geneid<-gsub(".+transcript_id ", "", tss$geneid)
# tss$geneid<-gsub(";.+", "", tss$geneid)
tss.uniq<-unique(tss)
ucsc <-tss.uniq

ucsc %>% count() # 17849
ucsc %>% head()
ucsc %>% count(geneid)
# ucsc %>% filter(str_detect(geneid, "NW")) %>% view()
###### legacy -end

##### tss with nuacc, 96563
# nuacc <- read.csv("/home/pkim/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/csRNA.NuAcc.tss.txt", header = T, sep = '\t') %>% 
nuacc <- read.csv("/Users/PanjunKim/UTHSC GGI Dropbox/K P/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/csRNA.NuAcc.tss.txt", header = T, sep = '\t') %>% 
  mutate(Anno3 = Detailed.Annotation) %>% 
  dplyr::select(Chr, Start, End, Annotation, Anno3) %>% 
  separate(Annotation, into=c("Anno1", "Anno2"), sep = ' ', remove=FALSE) %>%
  mutate(file = "nuacc") %>% 
  mutate(Anno2 = str_replace_all(Anno2, "\\(|\\)|,", "")) %>% 
  # filter(Anno2 == 'NR_132639') %>%
  # count(Anno2) %>% 
  view()

##### tss with pfc ver1, 131647
# pfc <- read.csv("/home/pkim/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/csRNA.PFC.tss.txt", header = T, sep = '\t') %>% 
pfc <- read.csv("/Users/PanjunKim/UTHSC GGI Dropbox/K P/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/csRNA.PFC.tss.txt", header = T, sep = '\t') %>% 
  dplyr::select(Chr, Start, End, Annotation) %>% 
  separate(Annotation, into=c("Anno1", "Anno2"), sep = ' ', remove=FALSE) %>%
  mutate(file = "PFC") %>% 
  mutate(Anno2 = str_replace_all(Anno2, "\\(|\\)|,", "")) %>% 
  # count(Anno2) %>%
  view()

# nuacc + pfc: 228210
df.tss.nuacc.pfc <- bind_rows(nuacc, pfc) %>% 
  view()

# legacy - start
# genomic range for tss from UCSC, 17849
# sapply(df.tss.ucsc, length)
# sum(is.na(df.tss.ucsc))
# df.tss.ucsc %>% dim()
# df.tss.ucsc %>% count(strand)
# str(df.tss.ucsc)
# colnames(df.tss.ucsc)
# legacy - end

# genomic range for tss from nuacc, 96563
df.tss.nuacc.GR <-GRanges(seqnames=nuacc$Chr, ranges=IRanges(start=nuacc$Start, end=nuacc$End), loc=nuacc$Anno1, geneid=nuacc$Anno2)
df.tss.nuacc.GR

# genomic range for tss from pfc ver1, 131647
df.tss.pfc.GR <- GRanges(seqnames=pfc$Chr, ranges=IRanges(start=pfc$Start, end=pfc$End), loc=pfc$Anno1, geneid = pfc$Anno2)
df.tss.pfc.GR

data.frame(df.tss.ucsc.GR) %>% 
  # data.frame(df.tss.pfc.GR) %>% 
  # data.frame(df.tss.nuacc.GR) %>% 
  count(seqnames)

# chr19_NW_023637717v1_random
# chrUn_NW_023637854v1
# chrY_NW_023637718v1_random

### original

overall.df.DISTINCT.loop.deep.sample.all.1.distance

overall.df.DISTINCT.loop.deep.sample.all.1.distance.GR <- GRanges(seqnames=overall.df.DISTINCT.loop.deep.sample.all.1.distance$chr1, ranges=IRanges(start=overall.df.DISTINCT.loop.deep.sample.all.1.distance$x0, end=overall.df.DISTINCT.loop.deep.sample.all.1.distance$y3), id=overall.df.DISTINCT.loop.deep.sample.all.1.distance$loop.id, new.loop.id=overall.df.DISTINCT.loop.deep.sample.all.1.distance$new.loop.id)

index.distinct.tss.w.whole.loop.1.distance <- findOverlaps(df.tss.ucsc.GR, overall.df.DISTINCT.loop.deep.sample.all.1.distance.GR, type = "within")

overall.loop.for.tss.hits <- subjectHits(index.distinct.tss.w.whole.loop.1.distance)
overall.tss.on.loop.hits <- queryHits(index.distinct.tss.w.whole.loop.1.distance)

df.tss.dist.result <- data.frame(
  loop.id = overall.df.DISTINCT.loop.deep.sample.all.1.distance$new.loop.id[overall.loop.for.tss.hits],
  loop.start = overall.df.DISTINCT.loop.deep.sample.all.1.distance$x0[overall.loop.for.tss.hits],
  loop.end = overall.df.DISTINCT.loop.deep.sample.all.1.distance$y3[overall.loop.for.tss.hits],
   loop.start.x12 = overall.df.DISTINCT.loop.deep.sample.all.1.distance$x12[overall.loop.for.tss.hits], ## Hao
  loop.end.y12 = overall.df.DISTINCT.loop.deep.sample.all.1.distance$y12[overall.loop.for.tss.hits],  ##
  loop.res=overall.df.DISTINCT.loop.deep.sample.all.1.distance$resolution[overall.loop.for.tss.hits],
  tss_chr = seqnames(df.tss.ucsc.GR)[overall.tss.on.loop.hits],
  tss_start = start(df.tss.ucsc.GR)[overall.tss.on.loop.hits],
  tss_end = end(df.tss.ucsc.GR)[overall.tss.on.loop.hits],
  tss_geneid = mcols(df.tss.ucsc.GR)$gene_id[overall.tss.on.loop.hits],
  tss_strand = strand(df.tss.ucsc.GR)[overall.tss.on.loop.hits]
)

df.tss.dist.result %>% dim() # 375567      9
df.tss.dist.result %>% head()

relative.pos.df.tss.dist.result <- df.tss.dist.result %>%
  mutate(pos_coord = (tss_start + tss_end) / 2) %>% # coordinate for midpoint
  mutate(loop_length = loop.end.y12 - loop.start.x12) %>% # overall length of the loop
  mutate(relative_pos = pos_coord - loop.start.x12) %>% # relative pos
  # mutate(value = (relative_pos/(loop_length / 2)) - 1) %>% .5 distance
  mutate(value = (relative_pos/loop_length ) - 1) %>% # 1 distance
  mutate(resolution = case_when(
    loop.res == 5000 ~ "5K",
    loop.res  == 10000 ~ "10K",
    loop.res  == 25000 ~ "25K",
    TRUE ~ NA_character_ 
  )) %>%
  dplyr::select(loop.id, tss_geneid, value, loop.res)

relative.pos.df.tss.dist.result %>% head()

## PDF file for ALL CHROMOSOMES by resolution
create_tss_plots <- function(data, res) {
  if (res == "All" ) {
    df.plot<-data
  } else { 
    df.plot<- data %>% filter(loop.res==res)
  }
# create_tss_plots <- function(data, res) {
  plot.tss.hist <- df.plot %>% 
    ggplot(aes(x = value)) +
    geom_histogram(fill = "skyblue", color = "black", alpha = 0.7, bins=200) +
    labs(title = paste0("A. Histogram of TSS over Loops on ALL Chromosomes (", res, ")"),
         x = "Relative position to loop",
         y = "Counts"
    )
  
  plot.tss.dens <- df.plot %>% 
    ggplot(aes(x = value)) +
    geom_density(fill = "skyblue", color = "black", alpha = 0.7) +
    ylim(c(0,1))+
    labs(title = paste0("A. Density of TSS over Loops on ALL chromosomes (", res, ")"),
         x = "Relative position to loop",
         y = "Density"
    )
  
  return(plot.tss.hist + plot.tss.dens)
}

plot_tss_all_res <- create_tss_plots(relative.pos.df.tss.dist.result, "All")
plot_tss_5K <- create_tss_plots(relative.pos.df.tss.dist.result, "5K")
plot_tss_10K <- create_tss_plots(relative.pos.df.tss.dist.result, "10K")
plot_tss_25K <- create_tss_plots(relative.pos.df.tss.dist.result, "25K")
pdf("./figures/overall_distribution_tss_on_all_chromosomes_1_distance.pdf", width = 11, height = 8.5)
# pdf("./figures/overall_distribution_tss_on_all_chromosomes_.5_distance.pdf", width = 11, height = 8.5)
(plot_tss_all_res/plot_tss_5K)
(plot_tss_10K / plot_tss_25K)
dev.off()

############################################
##### Overlapping TSS and loop processing
############################################
overall.df.DISTINCT.loop.deep.sample.all.1.distance

overall.df.DISTINCT.loop.deep.sample.all.1.distance.GR <- GRanges(seqnames=overall.df.DISTINCT.loop.deep.sample.all.1.distance$chr1, ranges=IRanges(start=overall.df.DISTINCT.loop.deep.sample.all.1.distance$x0, end=overall.df.DISTINCT.loop.deep.sample.all.1.distance$y3), id=overall.df.DISTINCT.loop.deep.sample.all.1.distance$loop.id, new.loop.id=overall.df.DISTINCT.loop.deep.sample.all.1.distance$new.loop.id )

index.distinct.tss.w.whole.loop.1.distance <- findOverlaps(df.tss.ucsc.GR, overall.df.DISTINCT.loop.deep.sample.all.1.distance.GR)

overall.loop.for.tss.hits <- subjectHits(index.distinct.tss.w.whole.loop.1.distance)
overall.tss.on.loop.hits <- queryHits(index.distinct.tss.w.whole.loop.1.distance)
head(overall.loop.for.tss.hits)
head(overall.tss.on.loop.hits)
head(index.distinct.tss.w.whole.loop.1.distance)
df.tss.ucsc.GR[1,]
overall.df.DISTINCT.loop.deep.sample.all.1.distance.GR[c(1822,1823,1824),]


df.tss.dist.result <- data.frame(
  loop.id = overall.df.DISTINCT.loop.deep.sample.all.1.distance$new.loop.id[overall.loop.for.tss.hits],
  loop.start = overall.df.DISTINCT.loop.deep.sample.all.1.distance$x0[overall.loop.for.tss.hits], ## Hao
  loop.end = overall.df.DISTINCT.loop.deep.sample.all.1.distance$y3[overall.loop.for.tss.hits], ## Hao
  loop.x12 = overall.df.DISTINCT.loop.deep.sample.all.1.distance$x12[overall.loop.for.tss.hits], ## Hao
  loop.y12 = overall.df.DISTINCT.loop.deep.sample.all.1.distance$y12[overall.loop.for.tss.hits],  ## Hao
  loop.res=overall.df.DISTINCT.loop.deep.sample.all.1.distance$resolution[overall.loop.for.tss.hits],
  tss_chr = seqnames(df.tss.ucsc.GR)[overall.tss.on.loop.hits],
  tss_start = start(df.tss.ucsc.GR)[overall.tss.on.loop.hits],
  tss_end = end(df.tss.ucsc.GR)[overall.tss.on.loop.hits],
  tss_geneid = mcols(df.tss.ucsc.GR)$gene_id[overall.tss.on.loop.hits],
  tss_strand = strand(df.tss.ucsc.GR)[overall.tss.on.loop.hits]
)

df.tss.dist.result %>% dim() # 375567      9
df.tss.dist.result %>% head()
df.tss.dist.result[,c("loop.x12")]
names(df.tss.dist.result)
relative.pos.df.tss.dist.result <- df.tss.dist.result %>%
  mutate(pos_coord = (tss_start + tss_end) / 2) %>% # coordinate for midpoint
  mutate(loop_length = (loop.y12 - loop.x12)*3) %>% # overall length of the loop # Hao
  mutate(relative_pos = pos_coord - loop.x12) %>% # relative pos 
  mutate(value = (relative_pos/(loop_length/3 )-1)) %>% # 1 distance
  mutate(Resolution=fct_relevel(loop.res, c("5K", "10K", "25K"))) |> 
  dplyr::select(loop.id, tss_geneid, value, Resolution, loop_length, relative_pos)

relative.pos.df.tss.dist.result %>% head()
unique(relative.pos.df.tss.dist.result$Resolution)



names(df.tss.dist.result)
## Start by Hao 
names(df.tss.dist.result)
absolute.pos.df.tss.dist.on.loop.result <- df.tss.dist.result %>%
  mutate(pos_coord = (tss_start + tss_end) / 2) %>% # coordinate for midpoint
  mutate(dist_to_x = pos_coord-loop.x12) %>% # overall length of the loop
  mutate(dist_to_y = pos_coord-loop.y12) %>% # overall length of the loop
  mutate(Resolution=fct_relevel(loop.res, c("5K", "10K", "25K")))|>
  mutate(dist_to_loop = ifelse( abs(dist_to_x) < abs(dist_to_y), dist_to_x, -1 * dist_to_y)) |> 
  mutate(loop.length=loop.y12-loop.x12) |>
  mutate(resolution = case_when(
    loop.res == 5000 ~ "5K",
    loop.res  == 10000 ~ "10K",
    loop.res  == 25000 ~ "25K",
  )) %>%
  dplyr::select(loop.id, tss_geneid, dist_to_loop, Resolution, loop.length)

head(absolute.pos.df.tss.dist.on.loop.result)

plot.absolute.distance.density.tss<-ggplot(absolute.pos.df.tss.dist.on.loop.result, aes(x=dist_to_loop/1000, group=Resolution, fill=Resolution))+geom_density(alpha=.5)+xlim(c(-1e3, 1e3))+xlab("Absolute distance, kb")

plot.relative.distance.hist.tss<-ggplot(relative.pos.df.tss.dist.result, aes(x=value, group=Resolution, fill=Resolution))+geom_histogram(bins=100) +xlab("Relative distance to loops")+theme(legend.position="none")

names(absolute.pos.df.tss.dist.on.loop.result)
plot.absolute.distance.hist.tss<- ggplot(absolute.pos.df.tss.dist.on.loop.result, aes(x=dist_to_loop/1000, group=Resolution, fill=Resolution))+geom_histogram(binwidth=250)+coord_cartesian(ylim=c(0,5000))+xlab("Absolute distance, kb")


pdf(file="hao_absolute_position_of_tss_to_loop.pdf", width=10, height=5)
plot.relative.distance.hist.tss+plot.absolute.distance.density.tss
dev.off()

## End By Hao

## PDF file by CHROMOSOME

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
      labs(title = paste0("A. Density of CTCF found near Hi-C loops on chr", chr),
           x = "Relative position to loop",
           y = "Density"
      )
    
    plot1.tss.hist <- relative.pos.df.tss.dist.result.chr %>% 
      ggplot(aes(x = value)) +
      geom_histogram(fill = "skyblue", color = "black", alpha = 0.7, bins=200) +
      labs(title = paste0("A. Histogram of CTCF found near Hi-C loops on chr", chr),
           x = "Relative position to loop",
           y = "Counts"
      )
    
    combined_plot <- plot1.tss.hist + plot1.tss.dens # + can be used instead of |
    plot_tss_list[[chr]] <- combined_plot
    
    message("END: Processing chromosome: ", chr)
  }, error = function(e) {
    
    message("Error processing chromosome: ", chr)
    message("Error message: ", e$message)
  })
}

pdf("./figures/overall_distribution_of_TSS_on_each_chr_.5_distance.pdf", width = 11, height = 8.5)
pdf("./figures/overall_distribution_of_TSS_on_each_chr_1_distance.pdf", width = 11, height = 8.5)

tss_num_plots <- length(plot_tss_list) # 22 chromosomes
tss_plots_per_page <- 4

for (i in seq(1, tss_num_plots, by = tss_plots_per_page)) {
  end_idx <- min(i + tss_plots_per_page - 1, tss_num_plots)
  page_plots <- plot_tss_list[i:end_idx]
  
  combined_page <- plot_grid(plotlist = page_plots, ncol = 1, nrow = 4)
  
  print(combined_page)
}

dev.off()

## PDF file for ALL CHROMOSOMES by resolution
create_tss_plots <- function(data, res) {
  if (res == "All" ) {
    df.plot<-data
  } else { 
    df.plot<- data %>% filter(loop.res==res)
  }
# create_tss_plots <- function(data, res) {
  plot.tss.hist <- df.plot %>% 
    ggplot(aes(x = value)) +
    geom_histogram(fill = "skyblue", color = "black", alpha = 0.7, bins=200) +
    labs(title = paste0("A. Histogram of TSS over Loops on ALL Chromosomes (", res, ")"),
         x = "Relative position to loop",
         y = "Counts"
    )
  
  plot.tss.dens <- df.plot %>% 
    ggplot(aes(x = value)) +
    geom_density(fill = "skyblue", color = "black", alpha = 0.7) +
    ylim(c(0,1))+
    labs(title = paste0("A. Density of TSS over Loops on ALL chromosomes (", res, ")"),
         x = "Relative position to loop",
         y = "Density"
    )
  
  return(plot.tss.hist + plot.tss.dens)
}

plot_tss_all_res <- create_tss_plots(relative.pos.df.tss.dist.result, "All")
plot_tss_5K <- create_tss_plots(relative.pos.df.tss.dist.result, "5K")
plot_tss_10K <- create_tss_plots(relative.pos.df.tss.dist.result, "10K")
plot_tss_25K <- create_tss_plots(relative.pos.df.tss.dist.result, "25K")
pdf("./figures/overall_distribution_tss_on_all_chromosomes_1_distance.pdf", width = 11, height = 8.5)
# pdf("./figures/overall_distribution_tss_on_all_chromosomes_.5_distance.pdf", width = 11, height = 8.5)
(plot_tss_all_res/plot_tss_5K)
(plot_tss_10K / plot_tss_25K)
dev.off()

################################################################
# 2. EACH END, boxplots
# TSS at each end in a loop: for the number of TSS in each end
# so, the object should be used one with padding
# df.DISTINCT.loop.deep.sample.all
################################################################
# df.DISTINCT.loop.deep.sample.all
# 2-1. adding padding at each end, x12 & y12
  
df.DISTINCT.loop.deep.sample.all

df.DISTINCT.loop.deep.sample.all.padded.for.TSS <- df.DISTINCT.loop.deep.sample.all %>% 
  mutate(x12 = (x1 + x2)/2, y12 = (y1 + y2)/2) %>% # middle point of each end
  mutate(distance = y12 - x12) %>% 
  mutate(x0 = ifelse(x12 - (distance / 2) < 0, 0, x12 - (distance / 2)), y3 = y12 + (distance/2)) %>% # 1/2 distance for OUTER PADDING
  mutate(x12 = (x12 + (distance / 4)), y12 = ifelse((y12 - (distance/4)) < 0, 0, y12 - (distance/4))) %>% # 1/4 distance for INNER PADDING
  mutate(new.loop.id = paste0(chr1, '_', x0, '_', x12, '_', chr2, '_', y12, '_', y3, '_', end.distance)) # new coordinates with the OUTER & INNER padding

df.DISTINCT.loop.deep.sample.all.padded.for.TSS %>% head()

df.DISTINCT.loop.deep.sample.all.padded.for.TSS.UP.GR<- GRanges(seqnames=df.DISTINCT.loop.deep.sample.all.padded.for.TSS$chr1, 
                                                                ranges=IRanges(start=df.DISTINCT.loop.deep.sample.all.padded.for.TSS$x0, 
                                                                               end=df.DISTINCT.loop.deep.sample.all.padded.for.TSS$x12), 
                                                                loop.id=df.DISTINCT.loop.deep.sample.all.padded.for.TSS$loop.id, 
                                                                end.distance=df.DISTINCT.loop.deep.sample.all.padded.for.TSS$end.distance,
                                                                resolution=df.DISTINCT.loop.deep.sample.all.padded.for.TSS$resolution,
                                                                new.loop.id=df.DISTINCT.loop.deep.sample.all.padded.for.TSS$new.loop.id)

# 1. TSS + UPSTREAM (df.DISTINCT.loop.deep.sample.all.padded.for.TSS.UP.GR, df.tss.ucsc.GR)
index.distinct.tss.w.up.loop.each.end <- findOverlaps(df.tss.ucsc.GR, df.DISTINCT.loop.deep.sample.all.padded.for.TSS.UP.GR, type = "within")
end.loop.up.tss.hits <- subjectHits(index.distinct.tss.w.up.loop.each.end)
end.tss.up.hits <- queryHits(index.distinct.tss.w.up.loop.each.end)

# df.DISTINCT.loop.deep.sample.all.padded.for.TSS %>% 
df.tss.ucsc.GR %>% 
  head()


df.overlapping.TSS.w.UPSTREAM.result.each.end <- data.frame(
  up.loop.id = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.UP.GR)$loop.id[end.loop.up.tss.hits],
  end.up.distance = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.UP.GR)$end.distance[end.loop.up.tss.hits],
  resolution = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.UP.GR)$resolution[end.loop.up.tss.hits],
  tss.id = mcols(df.tss.ucsc.GR)$tss.id[end.tss.up.hits],
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
df.DISTINCT.loop.deep.sample.all.padded.for.TSS.DOWN.GR<- GRanges(seqnames=df.DISTINCT.loop.deep.sample.all.padded.for.TSS$chr1, 
                                                                ranges=IRanges(start=df.DISTINCT.loop.deep.sample.all.padded.for.TSS$y12, 
                                                                               end=df.DISTINCT.loop.deep.sample.all.padded.for.TSS$y3), 
                                                                loop.id=df.DISTINCT.loop.deep.sample.all.padded.for.TSS$loop.id, 
                                                                end.distance=df.DISTINCT.loop.deep.sample.all.padded.for.TSS$end.distance,
                                                                resolution=df.DISTINCT.loop.deep.sample.all.padded.for.TSS$resolution,
                                                                new.loop.id=df.DISTINCT.loop.deep.sample.all.padded.for.TSS$new.loop.id)

# 1. TSS + DOWNSTREAM (df.DISTINCT.loop.deep.sample.all.padded.for.TSS.DOWN.GR, df.tss.ucsc.GR)
index.distinct.tss.w.down.loop.each.end <- findOverlaps(df.tss.ucsc.GR, df.DISTINCT.loop.deep.sample.all.padded.for.TSS.DOWN.GR, type = "within")
end.loop.down.tss.hits <- subjectHits(index.distinct.tss.w.down.loop.each.end)
end.tss.down.hits.end <- queryHits(index.distinct.tss.w.down.loop.each.end)

# df.DISTINCT.loop.deep.sample.all.padded.for.TSS %>% 
df.tss.ucsc.GR %>% 
  head()

df.overlapping.TSS.w.DOWNSTREAM.result.each.end <- data.frame(
  down.loop.id = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.UP.GR)$loop.id[end.loop.down.tss.hits],
  end.down.distance = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.UP.GR)$end.distance[end.loop.down.tss.hits],
  resolution = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.UP.GR)$resolution[end.loop.down.tss.hits],
  tss.id = mcols(df.tss.ucsc.GR)$tss.id[end.tss.down.hits.end],
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
df.overlapping.TSS.w.UPSTREAM.result.each.end %>% head()
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
head(df.overlapping.TSS.w.BOTH.result)$loop.id

# resolution      n
# 1         5K 21393 + 25473 =  46866
# 2        10K 28339 + 31699 =  60038
# 3        25K 46099 + 54162 = 100261

# for boxplot
df.overlapping.TSS.w.BOTH.result %>% head(2) # distance, resolution, ctcf.id, WHERE, loop.id, end.distance, case.id, chr

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
  group_by(resolution) %>%
  # group_by(WHERE, resolution) %>%
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


pdf("end_boxplot_num_tss_by_chr_res.pdf", width = 16.5, height = 23.5)
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

pdf("end_boxplot_num_tss_by_res.pdf", width = 16.5, height = 23.5)
grid.arrange(boxplot.w.TSS.by.res, ncol = 1)
dev.off()

pdf("end_boxplot_num_tss_combined.pdf", width = 16.5, height = 23.5)
grid.arrange(boxplot.w.TSS.by.chr.and.res.each.end, boxplot.w.TSS.by.res, ncol = 1)
dev.off()
tss.stats.by.resolution.each.end %>% head()

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

pdf("end_boxplot_no.tss_by_end_res.pdf", width = 16.5, height = 23.5)
grid.arrange(boxplot.w.TSS.ALL.chr, boxplot.w.TSS.ALL.chr.by.res, ncol = 1)
dev.off()

###############################################
# (1/2) Histogram of TSS Counts per Loop (Exclusive)
###############################################

df.overlapping.TSS.w.BOTH.result.boxplot %>% head()

tss_count_histogram_data <- df.overlapping.TSS.w.BOTH.result.boxplot %>%
  group_by(loop.id) %>%
  summarise(tss_count_total = sum(tss_count_by_loop_id_each_end), .groups = 'drop')

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
  summarise(total_tss = sum(tss_count_by_loop_id_each_end), .groups = 'drop')

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

# legacy
# one.sided.loop.histogram <- ggplot(one.sided.loop.histogram.data.final, aes(x = total_tss, y = one_sided_tss_loop_count)) +
#   geom_bar(stat = "identity", fill = "skyblue", color = "black") +
#   scale_x_continuous(limits = c(0, 6.5), breaks = seq(0, 6, by = 1)) +
#   labs(title = "Histogram of Loops with TSS Concentrated in UP or DOWN",
#        x = "Total TSS Count per Loop",
#        y = "Number of Loops") +
#   theme_minimal() +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   geom_text(aes(label = paste0(one_sided_tss_loop_count, "/", loop_count, " (", round(ratio * 100, 1), "%)")),
#             position = position_stack(vjust = 1.017), size = 3, color = "black")

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

pdf("./figures/tss_count_per_loop_histogram.pdf", width = 16.5, height = 23.5)

grid.arrange(tss.count.per.loop.histogram, one.sided.loop.histogram, ncol = 1)
dev.off()

##########################################################################################
##########################################################################################
# DISTRIBUTION: promoter distribution on loops
##########################################################################################
##########################################################################################

# promoter processing

file_path <- "/Users/PanjunKim/UTHSC GGI Dropbox/K P/Gateway_to_Hao/enhancer/data/epdnew/001/Rn_EPDnew_001_rn7.bed"
file_path <- "/Users/PanjunKim/UTHSC GGI Dropbox/K P/Gateway_to_Hao/enhancer/data/epdnew/001/Rn_EPDnew_001_rn7_1.bed"
file_path <- "../data/epdnew/001/Rn_EPDnew_001_rn7_1.bed"
df.promoter.rn7.raw <- read_tsv(file_path, col_names = c("chr", "start", "end", "gene", "score"))
df.promoter.rn7.raw %>% head()

df.promoter.rn7 <- df.promoter.rn7.raw %>% 
  mutate(length = end - start, 
         center = round((end + start)/2), 
         check.id = str_c(chr, '_', start, '_', end),
         promoter.id = str_c(check.id, '_', gene),)

df.promoter.rn7 %>% head()
df.promoter.rn7 %>% dim() # [1] 12463     9
df.promoter.rn7 %>% distinct(check.id) %>% dim() # 12427     1
df.promoter.rn7 %>% distinct(promoter.id) %>% dim() # 12463     1
df.promoter.rn7 %>% filter(check.id %in% check.id[duplicated(check.id)])

df.promoter.rn7.GR <- GRanges(seqnames=df.promoter.rn7$chr, 
                              ranges=IRanges(start=(df.promoter.rn7$start), 
                                             end=(df.promoter.rn7$end)), 
                              gene=df.promoter.rn7$gene, 
                              promoter.id=df.promoter.rn7$promoter.id, 
                              center=df.promoter.rn7$center, 
                              length=df.promoter.rn7$length)

# promoter processing done

####################################################
##### 1. OVERALL on loops Overlapping promoter & loops 
####################################################
overall.df.DISTINCT.loop.deep.sample.all.1.distance

overall.df.DISTINCT.loop.deep.sample.all <- overall.df.DISTINCT.loop.deep.sample.all.1.distance

overall.df.DISTINCT.loop.deep.sample.all.GR <- GRanges(seqnames=overall.df.DISTINCT.loop.deep.sample.all$chr1, ranges=IRanges(start=overall.df.DISTINCT.loop.deep.sample.all$x0, end=overall.df.DISTINCT.loop.deep.sample.all$y3), id=overall.df.DISTINCT.loop.deep.sample.all$loop.id, new.loop.id=overall.df.DISTINCT.loop.deep.sample.all$new.loop.id)

index.promoter.w.overall.whole.loop <- findOverlaps(df.promoter.rn7.GR, overall.df.DISTINCT.loop.deep.sample.all.GR, type = "within")

loop.for.promoter.hits <- subjectHits(index.promoter.w.overall.whole.loop)
promoter.on.loop.hits <- queryHits(index.promoter.w.overall.whole.loop)

df.promoter.dist.result <- data.frame(
  loop.id = overall.df.DISTINCT.loop.deep.sample.all$new.loop.id[loop.for.promoter.hits],
  loop.start = overall.df.DISTINCT.loop.deep.sample.all$x0[loop.for.promoter.hits],
  loop.end = overall.df.DISTINCT.loop.deep.sample.all$y3[loop.for.promoter.hits],
  loop.res = overall.df.DISTINCT.loop.deep.sample.all$resolution[loop.for.promoter.hits],
  loop.x12 = overall.df.DISTINCT.loop.deep.sample.all$x12[loop.for.promoter.hits],
  loop.y12 = overall.df.DISTINCT.loop.deep.sample.all$y12[loop.for.promoter.hits],
  # loop.new.distance = overall.df.DISTINCT.loop.deep.sample.all$new_distance[loop.for.promoter.hits],
  promoter.id = df.promoter.rn7$promoter.id[promoter.on.loop.hits],
  promoter.start = df.promoter.rn7$start[promoter.on.loop.hits],
  promoter.end = df.promoter.rn7$end[promoter.on.loop.hits]
) 

df.promoter.dist.result %>% head()
df.promoter.dist.result %>% dim() # 265565      7(1 distance)

relative.pos.df.promoter.dist.result <- df.promoter.dist.result %>%
  mutate(pos_coord = round((promoter.start + promoter.end) / 2)) %>%
  mutate(loop_length = (loop.y12 - loop.x12)*3) %>% # x0, y3
  mutate(relative_pos = pos_coord - loop.x12) %>% # relative position from loop start
  mutate(Resolution=fct_relevel(loop.res, c("5K", "10K", "25K"))) |>
  mutate(value = relative_pos/(loop_length/3)-1) %>% # for padding 1x distance
  dplyr::select(loop.id, promoter.id, value, Resolution)

relative.pos.df.promoter.dist.result %>% 
  dim() # 265565      4
  head()

  
absolute.pos.df.promoter.dist.result <- df.promoter.dist.result %>%
  mutate(pos_coord = (promoter.start + promoter.end) / 2) %>% # coordinate for midpoint
  mutate(dist_to_x = pos_coord-loop.x12) %>% # overall length of the loop
  mutate(dist_to_y = pos_coord-loop.y12) %>% # overall length of the loop
  mutate(Resolution=fct_relevel(loop.res, c("5K", "10K", "25K")))|>
  mutate(dist_to_loop = ifelse( abs(dist_to_x) < abs(dist_to_y), dist_to_x, -1 * dist_to_y)) |> 
  mutate(loop.length=loop.y12-loop.x12) |>
  dplyr::select(loop.id, dist_to_loop, Resolution, loop.length)

plot.absolute.distance.density.promoter<-ggplot(absolute.pos.df.promoter.dist.result, aes(x=dist_to_loop/1000, group=Resolution, fill=Resolution))+geom_density(alpha=0.5)+xlim(c(-1e3, 1e3))+xlab("Absolute distance, kb")

plot.relative.distance.hist.promoter<-ggplot(relative.pos.df.promoter.dist.result, aes(x=value, group=Resolution, fill=Resolution))+geom_histogram(bins=100) +xlab("Relative distance to loops")+theme(legend.position="none")
#plot.relative.distance.hist.promoter

names(absolute.pos.df.promoter.dist.on.loop.result)
plot.absolute.distance.hist.promoter<- ggplot(absolute.pos.df.promoter.dist.result, aes(x=dist_to_loop/1000, group=Resolution, fill=Resolution))+geom_histogram(binwidth=250)+coord_cartesian(ylim=c(0,5000))+xlab("Absolute distance, kb")


pdf(file="hao_position_of_promoter_to_loop.pdf", width=10, height=5)
plot.relative.distance.hist.promoter+plot.absolute.distance.density.promoter
dev.off()




## PDF file by CHROMOSOME

chromosomes <- c(1:20, "X", "Y")

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
      labs(title = paste0("A. Density of Promoter found near Hi-C loops on chr", chr),
           x = "Relative position to loop",
           y = "Density"
      )
    
    plot1.promoter.hist <- relative.pos.df.promoter.dist.result.chr %>% 
      ggplot(aes(x = value)) +
      geom_histogram(fill = "skyblue", color = "black", alpha = 0.7, bins=200) +
      labs(title = paste0("A. Histogram of Promoter found near Hi-C loops on chr", chr),
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

# pdf("./figures/overall_distribution_of_Promoter_on_each_chr_.5_distance.pdf", width = 11, height = 8.5)
pdf("./figures/overall_distribution_of_Promoter_on_each_chr_1_distance.pdf", width = 11, height = 8.5)

num_plots_promoter <- length(plot_promoter_list) # 22 chromosomes
plots_per_page_promoter <- 4

for (i in seq(1, num_plots_promoter, by = plots_per_page_promoter)) {
  end_idx <- min(i + plots_per_page_promoter - 1, num_plots_promoter)
  page_plots <- plot_promoter_list[i:end_idx]
  
  combined_page <- plot_grid(plotlist = page_plots, ncol = 1, nrow = 4)
  
  print(combined_page)
}

dev.off()

relative.pos.df.promoter.dist.result %>% head()

## PDF file for ALL CHROMOSOMES by resolution
create_promoter_plots <- function(data, res) {
  if (res == "All" ) {
    df.plot<-data
  } else { 
    df.plot<- data %>% filter(loop.res==res)
  }
  
  plot.promoter.hist <- df.plot %>% 
    ggplot(aes(x = value)) +
    geom_histogram(fill = "skyblue", color = "black", alpha = 0.7, bins=200) +
    labs(title = paste0("A. Histogram of Promoter over Loops on ALL Chromosomes (", res, ")"),
         x = "Relative position to loop",
         y = "Counts"
    )
  
  plot.promoter.dens <- df.plot %>% 
    ggplot(aes(x = value)) +
    geom_density(fill = "skyblue", color = "black", alpha = 0.7) +
    ylim(c(0,1))+
    labs(title = paste0("A. Density of Promoter found near Hi-C loops on ALL chromosomes"),
         x = "Relative position to loop",
         y = "Density"
    )
  
  return(plot.promoter.hist + plot.promoter.dens)
}

plot_promoter_all_res <- create_promoter_plots(relative.pos.df.promoter.dist.result, "All")
plot_promoter_5K <- create_promoter_plots(relative.pos.df.promoter.dist.result, "5K")
plot_promoter_10K <- create_promoter_plots(relative.pos.df.promoter.dist.result, "10K")
plot_promoter_25K <- create_promoter_plots(relative.pos.df.promoter.dist.result, "25K")
# pdf("overall_distribution_promoter_on_all_chromosomes_.5_distance.pdf", width = 11, height = 8.5) # pdf for .5 distance as a padding
pdf("./figure/overall_distribution_promoter_on_all_chromosomes_1_distance.pdf", width = 11, height = 8.5) # pdf for a distance as a padding
(plot_promoter_all_res/plot_promoter_5K)
(plot_promoter_10K / plot_promoter_25K)
dev.off()

################################################################
# 2. EACH END, boxplots
# Promoter at each end in a loop: for the number of promoter in each end
# so, the object should be used one with padding
# df.DISTINCT.loop.deep.sample.all
################################################################
# 2-1. adding padding at each end, x12 & y12
df.DISTINCT.loop.deep.sample.all

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

# 1. Promoter + UPSTREAM (df.DISTINCT.loop.deep.sample.all.padded.for.promoter.UP.GR, df.promoter.rn7.GR)
index.distinct.promoter.w.up.loop.each.end <- findOverlaps(df.promoter.rn7.GR, df.DISTINCT.loop.deep.sample.all.padded.for.promoter.UP.GR, type = "within")
loop.up.hits.promoter.each.end <- subjectHits(index.distinct.promoter.w.up.loop.each.end)
promoter.up.hits.each.end <- queryHits(index.distinct.promoter.w.up.loop.each.end)

# df.DISTINCT.loop.deep.sample.all.padded.for.promoter %>% 
df.promoter.rn7.GR %>% 
  head()

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
index.distinct.promoter.w.down.loop.each.end <- findOverlaps(df.promoter.rn7.GR, df.DISTINCT.loop.deep.sample.all.padded.for.promoter.DOWN.GR, type = "within")
end.loop.down.promoter.hits <- subjectHits(index.distinct.promoter.w.down.loop.each.end)
end.promoter.down.hits <- queryHits(index.distinct.promoter.w.down.loop.each.end)

# df.DISTINCT.loop.deep.sample.all.padded.for.promoter %>% 
df.promoter.rn7.GR %>% 
  head()

df.overlapping.promoter.w.DOWNSTREAM.result.each.end <- data.frame(
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
  group_by(resolution) %>%
  # group_by(WHERE, resolution) %>%
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


pdf("./figures/end_boxplot_num_promoter_by_chr_res.pdf", width = 16.5, height = 23.5)
grid.arrange(boxplot.w.promoter.by.chr.and.res.each.end, ncol = 1)
dev.off()

df.overlapping.promoter.w.BOTH.result.boxplot %>% 
  head()

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

pdf("./figures/end_boxplot_num_promoter_by_res.pdf", width = 16.5, height = 23.5)
grid.arrange(boxplot.w.promoter.by.res, ncol = 1)
dev.off()

pdf("end_boxplot_num_promoter_combined.pdf", width = 16.5, height = 23.5)
grid.arrange(boxplot.w.promoter.by.chr.and.res.each.end, boxplot.w.promoter.by.res, ncol = 1)
dev.off()
promoter.stats.by.resolution.each.end %>% head()

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

pdf("./figures/end_boxplot_no.promoter_by_end_res.pdf", width = 16.5, height = 23.5)
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

# legacy
# one.sided.loop.histogram <- ggplot(one.sided.loop.histogram.data.promoter.final, aes(x = total_promoter, y = one_sided_promoter_loop_count)) +
#   geom_bar(stat = "identity", fill = "skyblue", color = "black") +
#   scale_x_continuous(limits = c(0, 6.5), breaks = seq(0, 6, by = 1)) +
#   labs(title = "Histogram of Loops with Promoter Concentrated in UP or DOWN",
#        x = "Total Promoter Count per Loop",
#        y = "Number of Loops") +
#   theme_minimal() +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   geom_text(aes(label = paste0(one_sided_promoter_loop_count, "/", loop_count, " (", round(ratio * 100, 1), "%)")),
#             position = position_stack(vjust = 1.017), size = 3, color = "black")

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

pdf("./figures/promoter_count_per_loop_histogram.pdf", width = 16.5, height = 23.5)

grid.arrange(promoter.count.per.loop.histogram, one.sided.loop.histogram, ncol = 1)
dev.off()

#########################################################
# 4. diagram with loops from 3 previous steps
#########################################################
# distinct loops: 31773

df_ctcf_boxplot_none
df_ctcf_boxplot_50  
df_ctcf_boxplot_100 

ctcf_stats_none
ctcf_stats_50  
ctcf_stats_100 

# 4-1. CTCF
# important object
df.overlapping.CTCF.w.BOTH.result %>% 
  head()
  dim() # 2420002
df.ctcf.counts

head(df_ctcf_boxplot_50)$loop.id

process_final_loops <- function(df_ctcf_counts, ctcf_stats) {
  q1_values_by_resolution <- ctcf_stats %>%
    dplyr::select(resolution, Q1)
  # resolution    Q1
  # 1 5K             7
  # 2 10K           12
  # 3 25K           25
 # 
  loops_with_ctcf_above_q1 <- df_ctcf_counts %>%
    left_join(q1_values_by_resolution, by = "resolution") %>%
    filter(ctcf_count_by_loop_id >= Q1)
  loops_with_ctcf_above_q1 # 46,851
 # 
  loops.with.ctcf.both.ends <- loops_with_ctcf_above_q1 %>%
    group_by(loop.id) %>%
    filter(all(c("UP", "DOWN") %in% WHERE)) %>%
    ungroup()
  loops.with.ctcf.both.ends # 36,664 + 10
 # 
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
head(ctcf_final_loops_100)$loop.id
head(ctcf_final_loops_none)$loop.id

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
ggsave(filename = "./figures/distribution_of_loops_without_ctcf_bindings.pdf", plot = distribution_of_loops_without_ctcf_bindings, width = 8, height = 6)
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
ggsave(filename = "./figures/distribution_of_loops_without_ctcf_bindings_per_resolution.pdf", plot = distribution_of_loops_without_ctcf_bindings_per_resolution, width = 8, height = 6)

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

loops_only_one_end <- df.ctcf.counts %>%
  group_by(loop.id) %>%
  filter(n_distinct(WHERE) < 2) %>%
  ungroup() %>% 
  count(str_split_n(loop.id, '_', 1)) %>%
  dplyr::rename(chr = `str_split_n(loop.id, "_", 1)`, count = n) %>% 
  mutate(chr = factor(chr, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                                      "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", 
                                      "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                                      "chr20", "chrX", "chrY"))) 


  ggplot(aes(x = chr, y = count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  labs(title = "Distribution of Loops with CTCF Bindings in only One End per Chrnomosome", x = "Chromosome", y = "Count") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
ggsave(filename = "./figures/loops_with_ctcf_only_one_end.pdf", plot = loops_only_one_end, width = 8, height = 6)

loops_only_one_end

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
ggsave(filename = "./figures/loops_with_ctcf_only_one_end_per_resolution.pdf", plot = loops_only_one_end_per_resolution, width = 8, height = 6)

####### Analysis 3: loops less than Q1
df.ctcf.counts # 61,713
df.loops.with.ctcf.pairing.filtered <- df.ctcf.counts %>%
  group_by(loop.id) %>%
  filter(n_distinct(WHERE) == 2) %>% # n_distinct(WHERE) < 2 loop.ids are NO CTCF in the one end
  ungroup()

df.loops.with.ctcf.pairing.filtered # 60276
# 61713 - 60276 = 1437 : PASS

q1_values_by_resolution

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

ggsave(filename = "./figures/distribution_of_loops_with_less_ctcf_than_q1_hist.pdf", plot = df.loops.with.less.ctcf.than.q1.hist, width = 8, height = 6)

# Save all three plots to a single PDF file
pdf("./figures/plot_loops_with_ctcf_above_q1_complement_paired_by_resolution.pdf", height = 12, width = 8.5)
grid.arrange(plot.5k.loops.with.ctcf.above.q1.complement.paired, 
             plot.10k.loops.with.ctcf.above.q1.complement.paired, 
             plot.25k.loops.with.ctcf.above.q1.complement.paired, ncol = 1)
dev.off()

# 4-2. TSS
df.overlapping.TSS.w.BOTH.result %>% head()
df.tss.counts.each.end # 44,553 + 10
tss.stats.by.resolution.each.end
head(df.tss.counts.each.end)

loops.with.tss.above.q1 <- df.tss.counts.each.end %>% 
  left_join(tss.stats.by.resolution.each.end %>% 
              dplyr::select(resolution, Q1), by = "resolution") %>% 
  filter(tss_count_each_end >= 1)

loops.with.tss.both.ends <- loops.with.tss.above.q1 %>%
  group_by(loop.id) %>%
  filter(all(c("UP", "DOWN") %in% WHERE)) %>%
  ungroup()
loops.with.tss.both.ends
final.loops.from.tss.step <- loops.with.tss.above.q1 %>% 
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
  filter(promoter_count_each_end >= 1)

loops.with.promoter.both.ends <- loops.with.promoter.above.q1 %>% # 38,152 + 10
  group_by(loop.id) %>%
  filter(all(c("UP", "DOWN") %in% WHERE)) %>%
  ungroup()
dim(loops.with.promoter.both.ends) 

names(df.promoter.counts.each.end)
# [1] "loop.id"                
# [2] "WHERE"                  
# [3] "resolution"             
# [4] "promoter_count_each_end"
dim(df.promoter.counts.each.end)

final.loops.from.promoter.step <- df.promoter.counts.each.end %>% 
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

##############
# Venn Diagram
##############
dim(final.loops.from.promoter.step)
# [1] 25149     3
dim(final.loops.from.tss.step)
# [1] 27517     3
final.loops.from.ctcf.step<-ctcf_final_loops_100
dim(final.loops.from.ctcf.step)
# [1] 19392     3
head(final.loops.from.ctcf.step)[1,1]
# 1 chr1_101075000_101100000_chr1_101550000_101575000_25000_0
head(final.loops.from.promoter.step)[1,1]
# 1 chr10_1000000_1025000_chr10_1200000_1225000_25000_0
head(final.loops.from.tss.step)[1,1]
# 1 chr10_100320000_100330000_chr10_100390000_100400000_10000_0

loop.data.for.venn <- list(
  Promoter = final.loops.from.promoter.step$loop.id,
  TSS = final.loops.from.tss.step$loop.id,
  CTCF = final.loops.from.ctcf.step$loop.id
)

overall.loops.common.venn.plot <- ggvenn(
  loop.data.for.venn, 
  fill_color = c("#E41A1C", "#377EB8", "#4DAF4A")
) + ggtitle("Overall Loop Overlap") +
  theme(
    plot.title = element_text(hjust = 0.5),       
    legend.title = element_text(size = 0.3)         
  )

pdf("./figures/loops_common_venn_diagram.pdf")
pdf("hao_loops_common_venn_diagram.pdf", width=6, height=6)
print(overall.loops.common.venn.plot)
dev.off()



list_of_loops<-unique(c(final.loops.from.ctcf.step$loop.id, final.loops.from.tss.step$loop.id, final.loops.from.promoter.step$loop.id))
length(list_of_loops)
idx.ctcf<-list_of_loops %in%  final.loops.from.ctcf.step$loop.id  
idx.tss<-list_of_loops %in%  final.loops.from.tss.step$loop.id
idx.promoter<-list_of_loops %in%  final.loops.from.promoter.step$loop.id  
sum(idx.ctcf)
idx.final<- idx.ctcf & (idx.promoter | idx.tss) 
final.loops<-data.frame(loop.id=list_of_loops[idx.final])
dim(final.loops)
# [1] 18216

# Split the strings and separate into columns using separate() from tidyr
final.loops.df <- final.loops %>%
  separate(loop.id, 
           into = c("chr1", "start1", "end1", "chr2", "start2", "end2", "resolution", "none"), 
           sep = "_")
head(final.loops.df)
write.table(file="hao_final_loops.tab", row.names=F, final.loops.df)

# by resolution
loop.data.by.resolution <- list(
  `25K` = c(final.loops.from.promoter.step$loop.id[final.loops.from.promoter.step$end.distance == "25000"],
            final.loops.from.tss.step$loop.id[final.loops.from.tss.step$end.distance == "25000"],
            final.loops.from.ctcf.step$loop.id[final.loops.from.ctcf.step$end.distance == "25000"]),
  `10K` = c(final.loops.from.promoter.step$loop.id[final.loops.from.promoter.step$end.distance == "10000"],
            final.loops.from.tss.step$loop.id[final.loops.from.tss.step$end.distance == "10000"],
            final.loops.from.ctcf.step$loop.id[final.loops.from.ctcf.step$end.distance == "10000"]),
  `5K` = c(final.loops.from.promoter.step$loop.id[final.loops.from.promoter.step$end.distance == "5000"],
           final.loops.from.tss.step$loop.id[final.loops.from.tss.step$end.distance == "5000"],
           final.loops.from.ctcf.step$loop.id[final.loops.from.ctcf.step$end.distance == "5000"])
)

loop.data.for.venn.by.resolution <- list()
for (resolution in c("5K", "10K", "25K")) {
  loop.data.for.venn.by.resolution[[resolution]] <- ggvenn(
    list(
      Promoter = final.loops.from.promoter.step$loop.id[final.loops.from.promoter.step$resolution == resolution],
      TSS = final.loops.from.tss.step$loop.id[final.loops.from.tss.step$resolution == resolution],
      CTCF = final.loops.from.ctcf.step$loop.id[final.loops.from.ctcf.step$resolution == resolution]
    ),
    fill_color = c("#E41A1C", "#377EB8", "#4DAF4A")
  ) + ggtitle(paste("Resolution Overlap", "(", resolution, ")" )) +
    theme(
      plot.title = element_text(hjust = 0.5),       
      legend.title = element_text(size =0.3)       
    )
}

pdf("./figures/venn_diagrams_by_resolution.pdf")
for (plot in loop.data.for.venn.by.resolution) {
  print(plot)
}
dev.off()

pdf("./figures/combined_venn_diagrams.pdf")
grid.arrange(
  overall.loops.common.venn.plot, 
  loop.data.for.venn.by.resolution[["5K"]], 
  loop.data.for.venn.by.resolution[["10K"]], 
  loop.data.for.venn.by.resolution[["25K"]], 
  ncol = 2, nrow = 2
)
dev.off()



###############
# refactoring
###############

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
  output_path = "./figures/0829/ctcf_vs_promoter_tss_venn_diagrams.pdf"
)





















##########################################################################################
##########################################################################################
## TSS VS promoter check
##########################################################################################
##########################################################################################

# grep -E '^(ID|DR|SE)' Rn_EPDnew_001_rn6.dat | awk '{if ($1 == "ID") $2 = "ID; " $2; else if ($1 == "SE") $2 = "Seq; " $2; print}' | awk '{$1=""; sub(/^ /, ""); print}' | awk -F";" '{print $1, substr($0, index($0, $2))}' OFS="\t" | sed 's/Gene Symbol/GeneSymbol/g' | awk '{first_col=$1; $1=""; sub(/^ /, ""); gsub(/[\t ]+/, ":", $0); print first_col, $0}' OFS="\t" | sed 's/[;.]//g'| sed 's/:\([^:]*\)/\1/' > Rn_EPDnew_001_rn6_2col.ts

# making df from data file: input for liftover
file_path <- "/Users/PanjunKim/UTHSC GGI Dropbox/K P/Gateway_to_Hao/enhancer/data/epdnew/001/Rn_EPDnew_001_rn6_2col.tsv"
df.promoter.rn6.data <- read_tsv(file_path, col_names = c("Key", "Value")) 

list.promoter.rn6.data <- split(df.promoter.rn6.data, cumsum(df.promoter.rn6.data$Key == "ID"))

df.promoter.rn6.data # 98,536 + 10

required_keys <- c("ID", "UCSC", "Ensembl", "RefSeq", "GeneSymbol", "Seq")

processed.list.promoter.rn6.data <- map(list.promoter.rn6.data, function(df) {
  
  df <- df %>% 
    complete(Key = required_keys, fill = list(Value = NA))
  
  df <- df %>%
    group_by(Key) %>%
    summarise(Value = paste(na.omit(Value), collapse = "|"), .groups = 'drop')
  
  return(df)
})

df.processed.list.promoter.rn6.data <- map_dfr(processed.list.promoter.rn6.data, function(df) {
  df_wide <- df %>%
    pivot_wider(names_from = Key, values_from = Value) %>%
    dplyr::select(ID, UCSC, Ensembl, RefSeq, GeneSymbol, Seq) %>% 
    mutate(UCSC = str_c(UCSC, ':', (as.numeric(str_split_n(UCSC, ':', 4)) + str_length(Seq) -1))) 
  return(df_wide)
})

df.processed.list.promoter.rn6.data # 12601
df.processed.list.promoter.rn6.data.liftover.input <- df.processed.list.promoter.rn6.data %>%
  mutate(chr = str_split_n(UCSC, ':', 2), strand = str_split_n(UCSC, ':', 3), start = as.numeric(str_split_n(UCSC, ':', 4)), end = as.numeric(str_split_n(UCSC, ':', 5))) %>% 
  mutate(start = as.numeric(start),
         end = as.numeric(end)) %>% 
  mutate(new.loop.id = str_c(str_split_n(ID, ':', 1), ':', Ensembl, ':', RefSeq, ':', GeneSymbol)) %>% 
  dplyr::select(chr, start, end, new.loop.id)
  
write.table(df.processed.list.promoter.rn6.data.liftover.input, "Rn_EPDnew_001_rn6_liftover_input_1.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


df.promoter.rn7 # 12463
df.promoter.rn7.GR 

df.tss.ucsc # 17849 
# chr,start,end,strand,gene_id,transcript_id,exon_number,exon_id,gene_name

###################################
# MATCHING with gene id & RefSeq
###################################
df.promoter.rn7.join <- df.promoter.rn7 %>%
  separate(gene, into = c("gene_id_part", "Ensembl", "RefSeq", "GeneSymbol"), sep = ":", extra = "merge", remove = FALSE)
df.promoter.rn7.join 

matched_by_gene_refseq <- df.tss.ucsc %>%
  inner_join(df.promoter.rn7.join, by = c("gene_id" = "GeneSymbol", "transcript_id" = "RefSeq"), relationship = "many-to-many")

# chr,start,end,gene,gene_id_part,Ensembl,RefSeq,GeneSymbol,score

matched_by_gene_refseq.distance <- matched_by_gene_refseq %>% # 8714 + 62
  mutate(distance = ifelse(chr.x == chr.y, ifelse((start.x + 1 >= end.y), abs(start.x - end.y), ifelse((start.x + 1 <= start.y), start.x - start.y, 0)), NA)) %>% 
  filter(between(distance, -distance_threshold, distance_threshold))
# matched_by_gene_refseq.distance %>% filter(is.na(distance)) # 12 NA
#   count(distance)
# chr.x, start.x, end.x, strand, gene_id, transcript_id, exon_number, exon_id, gene_name, chr.y, start.y, end.y, gene, gene_id_part, Ensembl, score

matched_by_gene_refseq.distance.boxplot <- ggplot(matched_by_gene_refseq.distance %>% filter(!is.na(distance)) %>% count(distance), aes(y = distance)) +
  geom_boxplot(fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Boxplot of of Distance between TSS & promoter by gene & refseq",
       y = "Distance",
       x = "") +
  theme_minimal()

# PDF
ggsave("matched_by_gene_refseq.distance.boxplot.pdf", plot = matched_by_gene_refseq.distance.boxplot, width = 8, height = 6, units = "in")


matched_by_gene_refseq_distance_hist <- 
  ggplot(matched_by_gene_refseq.distance %>% filter(!is.na(distance)), aes(x = distance)) +
  geom_histogram(binwidth = 1000, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Distance between TSS & promoter by gene & refseq",
       x = "Distance",
       y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("matched_by_gene_refseq_distance_hist.pdf", plot = matched_by_gene_refseq.distance.hist, 
       width = 8, height = 6, units = "in")

############################
# MATCHING with coordinates
############################
distance_threshold <- 2000

df.tss.ucsc <- df.tss.ucsc %>%
  mutate(start = as.numeric(start), end = as.numeric(end))

df.promoter.rn7 <- df.promoter.rn7 %>%
  mutate(start = as.numeric(start), end = as.numeric(end))

df.promoter.rn7
df.tss.ucsc

df.tss.ucsc.match.GR <- GRanges(
  seqnames = df.tss.ucsc$chr,
  ranges = IRanges(start = df.tss.ucsc$start, end = df.tss.ucsc$end),
  # strand = df.tss.ucsc$strand,
  gene_id = df.tss.ucsc$gene_id,
  transcript_id = df.tss.ucsc$transcript_id
)

df.promoter.rn7.threshold.match.GR <- GRanges(
  seqnames = df.promoter.rn7$chr,
  ranges = IRanges(start = df.promoter.rn7$start - distance_threshold, end = df.promoter.rn7$end + distance_threshold),
  # strand = df.promoter.rn7$strand,
  gene = df.promoter.rn7$gene
)

# same level of chromosomes
common_seqlevels_match <- intersect(seqlevels(df.tss.ucsc.match.GR), seqlevels(df.promoter.rn7.threshold.match.GR))

df.tss.ucsc.common.match.GR <- keepSeqlevels(df.tss.ucsc.match.GR, common_seqlevels_match, pruning.mode="coarse")
df.promoter.rn7.common.threshold.match.GR <- keepSeqlevels(df.promoter.rn7.threshold.match.GR, common_seqlevels_match, pruning.mode="coarse")

tss.promoter.match.hits <- findOverlaps(df.tss.ucsc.common.match.GR, df.promoter.rn7.common.threshold.match.GR)

df.tss.ucsc.common.match.GR$idx <- 1:length(df.tss.ucsc.common.match.GR)
df.promoter.rn7.common.threshold.match.GR$idx <- 1:length(df.promoter.rn7.common.threshold.match.GR)

matched_tss_promoter <- as.data.frame(tss.promoter.match.hits) %>%
  mutate(
    tss_idx = queryHits(tss.promoter.match.hits),
    promoter_idx = subjectHits(tss.promoter.match.hits)
  ) %>%
  inner_join(as.data.frame(df.tss.ucsc.common.match.GR), by = c("tss_idx" = "idx")) %>%
  inner_join(as.data.frame(df.promoter.rn7.common.threshold.match.GR), by = c("promoter_idx" = "idx")) %>%
  # view()
  # x: tss, y: promoter
  # mutate(distance = abs(start.y - start.x)) %>%
  mutate(distance = ifelse(seqnames.x == seqnames.y, ifelse((start.x + 1 >= (end.y - distance_threshold)), start.x - (end.y - distance_threshold), ifelse((start.x + 1 <= (start.y + distance_threshold)), start.x - (start.y + distance_threshold), 0)), NA))
  dplyr::select(
    chr = seqnames.x, 
    start = start.x, 
    end = end.x, 
    gene_id, 
    tss_start = start.y, 
    tss_end = end.y, 
    tss_geneid = gene_id, 
    tss_transcript_id = transcript_id,
    distance
  )

print(matched_tss_promoter)

matched_tss_promoter %>% count(distance)

total_n <- nrow(matched_tss_promoter)

matched_tss_promoter_distance_hist <- ggplot(matched_tss_promoter, aes(x = distance)) +
  geom_histogram(binwidth = 100, fill = "skyblue", color = "black", alpha = 0.7) +
  geom_text(stat = 'bin', aes(label = after_stat(count)), vjust = -0.5, binwidth = 100, color = "black", size = 3.5) +
  labs(title = "Histogram of Distance between TSS & Promoter by Coordinates",
       x = "Distance",
       y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  
  annotate("text", x = Inf, y = Inf, label = paste("Total n =", total_n), 
           hjust = 1.1, vjust = 2, size = 5, color = "black")

ggsave("matched_tss_promoter_distance_histogram.pdf", plot = matched_tss_promoter_distance_hist, width = 8, height = 6, units = "in")

######################################
# tss vs loops overlapping
######################################

# to filter valid loop, when checking a loop, I will check the number of CTCF at an END in a loop
# for example Q1  = 22 in 5k, a loop has 21 in upstream, 23 in downstream
# it would be classified invalid loop

# step 1: deciding to padding through distribution of tss
# step 1-1: checking distribution of loop lengths
#########################################
# LENTH of loops distribution: Start
#########################################

loop.length.stat.summary.table <- overall.df.DISTINCT.loop.deep.sample.all %>%
  mutate(loop.length = y12 - x12) %>%
  filter(!is.na(loop.length) & is.finite(loop.length)) %>%
  group_by(resolution) %>%
  summarise(
    Min = min(loop.length),
    Q1 = quantile(loop.length, 0.25),
    Median = median(loop.length),
    Q3 = quantile(loop.length, 0.75),
    Max = max(loop.length)
  )

print(loop.length.stat.summary.table)
35000 - 5000*2 = 25000
# loop length stats
# resolution      Min     Q1 Median     Q3       Max
# 1 5K          35000  65000 125000 260000 211685000
# 2 10K         50000  90000 170000 330000 219090000
# 3 25K        100000 150000 275000 525000 211675000


# loops
# loops with TSS only: in one only, in both ends
# loops with promoter only: in one only, in both ends
# loops with both: in one only, in both ends

q3_values <- overall.df.DISTINCT.loop.deep.sample.all %>%
  mutate(loop.length = y12 - x12) %>%
  group_by(resolution) %>%
  summarise(Q3 = quantile(loop.length, 0.75))

y_max <- max(q3_values$Q3) * 3

loop_length_distribution_by_resolution <- overall.df.DISTINCT.loop.deep.sample.all %>%
  mutate(loop.length = y12 - x12) %>%
  # filter(is.na(loop.length)) %>% 
  # view()
  ggplot(aes(x = resolution, y = loop.length)) +
  geom_boxplot() +
  stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red") +
  stat_summary(fun.data = function(y) {
    return(data.frame(y = quantile(y, probs = c(0.25, 0.5, 0.75))))
  }, geom = "text", aes(label = round(..y.., 1)), position = position_nudge(x = 0.2)) +
  labs(title = "Loop Length Distribution by Resolution",
       x = "Resolution",
       y = "Loop Length") +
  theme_minimal() +
  ylim(NA, y_max) +  
  theme(plot.title = element_text(hjust = 0.5))

loop_length_distribution_by_resolution
ggsave("loop_length_distribution.pdf", plot = loop_length_distribution_by_resolution, width = 8, height = 6)

#########################################
# LENTH of loops distribution: step 1-1 End
#########################################
# step 1-2: deciding each end
# inner: distance*1/4
# outer: distance*1/2

# step 2: 

df.ctcf.counts
ctcf_stats_by_resolution

# valid loop
df.filtered.loops.by.ctcf <- df.ctcf.counts %>%
  inner_join(ctcf_stats_by_resolution, by = c("resolution", "WHERE")) %>% 
  filter(ctcf_count >= Q1) %>%
  group_by(loop.id) %>%
  filter(n_distinct(WHERE) == 2) %>%
  ungroup()

df.filtered.loops.by.ctcf %>% head()

df.filtered.loops.by.ctcf %>% dim()

# checking overlapping with TSS
df.filtered.loops.by.ctcf.processed <- df.filtered.loops.by.ctcf %>% 
  distinct(loop.id) %>% 
  separate(loop.id, into = c("chr1", "x1", "x2", "chr2", "y1", "y2", "end.distance"), sep = "_", remove = FALSE) %>% 
  mutate(x1 = as.numeric(x1), x2 = as.numeric(x2), y1 = as.numeric(y1), y2 = as.numeric(y2), end.distance = as.numeric(end.distance)) %>% 
  mutate(x12 = (x1 + x2)/2, y12 = (y1 + y2)/2) %>% # middle point of each end
  mutate(distance = y12 - x12) %>% 
  mutate(x0 = ifelse(x12 - (distance / 2) < 0, 0, x12 - (distance / 2)), y3 = y12 + (distance/2)) %>% # 1/2 distance for OUTER PADDING
  mutate(x12 = (x12 + (distance / 4)), y12 = ifelse((y12 - (distance/4)) < 0, 0, y12 - (distance/4))) %>% # 1/4 distance for INNER PADDING
  mutate(resolution = case_when(
    end.distance == 5000 ~ "5K",
    end.distance  == 10000 ~ "10K",
    end.distance  == 25000 ~ "25K",
    TRUE ~ NA
  ))

df.filtered.loops.by.ctcf.processed %>%
  head()
  count() #18270
  # dplyr::select(loop.id)
# loop.id                                                   chr1         x1        x2 chr2         y1        y2 end.distance resolution

# TSS with UP
df.filtered.UP.loops.by.ctcf.processed.GR <- GRanges(seqnames=df.filtered.loops.by.ctcf.processed$chr1, ranges=IRanges(start=(df.filtered.loops.by.ctcf.processed$x0), end=(df.filtered.loops.by.ctcf.processed$x12)), loop.id=df.filtered.loops.by.ctcf.processed$loop.id, end.distance = df.filtered.loops.by.ctcf.processed$end.distance, resolution = df.filtered.loops.by.ctcf.processed$resolution)

index.valid.UP.loop.with.tss.and.ctcf <- findOverlaps(df.tss.ucsc.GR, df.filtered.UP.loops.by.ctcf.processed.GR, type = "within")

loop.valid.UP.hits <- subjectHits(index.valid.UP.loop.with.tss.and.ctcf)
tss.on.loop.UP.hits <- queryHits(index.valid.UP.loop.with.tss.and.ctcf)


df.loop.UPTREAM.with.tss.result <- data.frame(
  up.loop.id = mcols(df.filtered.UP.loops.by.ctcf.processed.GR)$loop.id[loop.valid.UP.hits],
  end.up.distance = mcols(df.filtered.UP.loops.by.ctcf.processed.GR)$end.distance[loop.valid.UP.hits],
  resolution = mcols(df.filtered.UP.loops.by.ctcf.processed.GR)$resolution[loop.valid.UP.hits],
  tss.id = mcols(df.tss.ucsc.GR)$tss.id[tss.on.loop.UP.hits],
  WHERE = "UP"
) %>% 
  mutate(case.tss.up.id = str_c(up.loop.id, '|', tss.id))

df.loop.UPTREAM.with.tss.result %>% dim() # 45168
df.loop.UPTREAM.with.tss.result %>% 
  count(resolution)
# resolution     n
# 1        10K 14737
# 2        25K 23894
# 3         5K  6537

# total: 45168

df.loop.UPTREAM.with.tss.result %>% 
  head()
  count()

df.filtered.DOWN.loops.by.ctcf.processed.GR <- GRanges(seqnames=df.filtered.loops.by.ctcf.processed$chr1, ranges=IRanges(start=(df.filtered.loops.by.ctcf.processed$y12), end=(df.filtered.loops.by.ctcf.processed$y3)), loop.id=df.filtered.loops.by.ctcf.processed$loop.id, end.distance = df.filtered.loops.by.ctcf.processed$end.distance, resolution = df.filtered.loops.by.ctcf.processed$resolution)

index.valid.DOWN.loop.with.tss.and.ctcf <- findOverlaps(df.tss.ucsc.GR, df.filtered.DOWN.loops.by.ctcf.processed.GR, type = "within")

loop.valid.DOWN.hits <- subjectHits(index.valid.DOWN.loop.with.tss.and.ctcf)
tss.on.loop.DOWN.hits <- queryHits(index.valid.DOWN.loop.with.tss.and.ctcf)


df.loop.DOWNTREAM.with.tss.result <- data.frame(
  down.loop.id = mcols(df.filtered.DOWN.loops.by.ctcf.processed.GR)$loop.id[loop.valid.DOWN.hits],
  end.down.distance = mcols(df.filtered.DOWN.loops.by.ctcf.processed.GR)$end.distance[loop.valid.DOWN.hits],
  resolution = mcols(df.filtered.DOWN.loops.by.ctcf.processed.GR)$resolution[loop.valid.DOWN.hits],
  tss.id = mcols(df.tss.ucsc.GR)$tss.id[tss.on.loop.DOWN.hits],
  WHERE = "DOWN"
) %>% 
  mutate(case.tss.down.id = str_c(down.loop.id, '|', tss.id))

df.loop.DOWNTREAM.with.tss.result %>% dim() # 45639
df.loop.DOWNTREAM.with.tss.result %>% 
  count(resolution)
# resolution     n
# 1        10K 14388
# 2        25K 24957
# 3         5K  6294

########## bind_rows(UPSTREAM & DOWNSTREAM) -> BOTH
df.overlapping.TSS.w.BOTH.for.valid.loop.result <- bind_rows(df.loop.UPTREAM.with.tss.result %>% 
                                                 mutate(loop.id = up.loop.id, end.distance = end.up.distance) %>% 
                                                 mutate(case.id = case.tss.up.id) %>% 
                                                 dplyr::select(-c(up.loop.id, end.up.distance, case.tss.up.id)), 
                                               df.loop.DOWNTREAM.with.tss.result %>% 
                                                 mutate(loop.id = down.loop.id, end.distance = end.down.distance) %>% 
                                                 mutate(case.id = case.tss.down.id) %>% 
                                                 dplyr::select(-c(down.loop.id, end.down.distance, case.tss.down.id))) %>% 
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

df.overlapping.TSS.w.BOTH.for.valid.loop.result %>% dim() # 90807
df.overlapping.TSS.w.BOTH.for.valid.loop.result %>% head()
df.overlapping.TSS.w.BOTH.for.valid.loop.result %>% count(resolution)

# resolution     n
# 1         5K 12831
# 2        10K 29125
# 3        25K 48851

# loop.counts.by.chromosome.and.resolution <- 
  
  
  df.overlapping.TSS.w.BOTH.for.valid.loop.result %>% 
    head()
  # distinct(loop.id) %>% 
  separate(loop.id, into = c("chr1", "x1", "x2", "chr2", "y1", "y2", "end.distance"), sep = "_") %>% 
  mutate(resolution = case_when(
    end.distance == 5000 ~ "5K",
    end.distance == 10000 ~ "10K",
    end.distance == 25000 ~ "25K",
    TRUE ~ NA
  )) %>% 
  count(chr1, resolution)

  # count(loop.id) # 16637
valid.loop <- ggplot(loop.counts.by.chromosome.and.resolution, aes(x = chr1, y = n, fill = resolution)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_discrete(limits = paste0("chr", c(1:20, "X", "Y"))) + 
  labs(title = "Number of Loops by Chromosome and Resolution",
       x = "Chromosome",
       y = "Number of Loops") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

valid.loop

ggsave("valid_loop_by_chromosome_and_resolution.pdf", plot = valid.loop, width = 10, height = 6)





########################################################################################################################
# Found from slack that I sent this to you back in Nov 2022. 
# Please modify the code in two places, 
# 1. change the ensembl annotation to refseq annotation, 
# 2. make sure the bed file contains loop data from all four samples

library("GenomicRanges")
library("ggplot2")

# get the TSS location into a GenomicRange object
tss<-read.table(file="./ensembl_mRatBN7.2_TSS.gtf", sep="\t", head=F)
tss<-tss[,c(1,4,5)]
names(tss)<-c("chr", "start", "end")
tss$chr<-paste0("chr", tss$chr)
head(tss)

#genomic range for tss
grTSS<-GRanges(seqnames=tss$chr, ranges=IRanges(start=tss$start, end=tss$end))
head(grTSS)

# get the loop coordinates into a GR
loop<-read.table(file="./merged_loops.bedpe")
head(loop)
grLoop<-GRanges(seqnames=loop$V1, ranges=IRanges(start=loop$V2, end=loop$V3))
head(grLoop)

## find the nearest location
idxnear<-nearest(grLoop, grTSS)
# reorder grTSS to follow the nearest grLoop
reorderedTSS<-grTSS[idxnear,]

# put the loops and their nearest TSS together
out<-data.frame(grLoop, reorderedTSS)
names(out)<-c("chr","loopStart", "loopEnd", "loopWidth", "loopStrand", "chrTSS", "TSSStart", "TSSEnd", "TSSWidth", "TSSStrand")
# calcualte the middle of the Loop and TSS
out$loopMid<-(out$loopEnd-out$loopStart)/2+out$loopStart
# distance between Loop and TSS
head(out)
out$loopToTSS<-out$loopMid-out$TSSStart
head(out)


pdf(file="Distance_between_loop_and_TSS.pdf", width=6, height=5)
ggplot(data=out, aes(x=loopToTSS))+geom_histogram(fill="darkblue", color="grey")+xlim(c(-30000,30000))+xlab("Distance between loop and nearest TSS, bp")
dev.off()
########################################################################################################################




########################################################################################################################
# get the TSS location into a GenomicRange object, note gene name is added, also there are some duplicated lines, so use uniq
tss<-read.table(file="./ucsc_start_codon.txt", sep="\t", head=F)
head(tss)[,c(1,4,5,9)]
tss<-tss[,c(1,4,5,9)] # onlty take the relevant columns
names(tss)<-c("chr", "start", "end", "geneid")
tss$geneid<-gsub(".+transcript_id ", "", tss$geneid)
tss$geneid<-gsub(";.+", "", tss$geneid)
tss.uniq<-unique(tss)
tss<-tss.uniq
tss
#genomic range for tss
grTSS<-GRanges(seqnames=tss$chr, ranges=IRanges(start=tss$start, end=tss$end), geneid=tss$geneid)

# combines loop with their closest TSS.
out<-data.frame(grLoop, reorderedTSS)
# take only two columns from the dataframe above, 
dfoverlap<-out[,c("Gene_ID", "strain")]
# create a pseudo column so that I can change it to wide format
dfoverlap$val<-1 
# change to wide
dfoverlapWide<-spread(unique(dfoverlap), strain, val )

##############################
# setup for upset plot
##############################
strains<-names(dfoverlapWide)[-1]
upsetplot<-upset(dfoverlapWide,strains)
upsetplot
########################################################################################################################
# the data is not correct but the figure shows what it will look like.
########################################################################################################################


findOverlaps.loop.SHR.OlaIpcv.neuron.EAtlas2.GR <- findOverlaps(loop.SHR.OlaIpcv.GR, df.enhancer.neuron.EAtlas2.GR)
findOverlaps.loop.HXB10.BB.neuron.EAtlas2.GR <- findOverlaps(loop.HXB10.BB.GR, df.enhancer.neuron.EAtlas2.GR)
findOverlaps.loop.HXB10.CC.neuron.EAtlas2.GR <- findOverlaps(loop.HXB10.CC.GR, df.enhancer.neuron.EAtlas2.GR)
findOverlaps.loop.F344.Stm.neuron.EAtlas2.GR <- findOverlaps(loop.F344.Stm.GR, df.enhancer.neuron.EAtlas2.GR)
findOverlaps.loop.LE.Stm.neuron.EAtlas2.GR <- findOverlaps(loop.LE.Stm.GR, df.enhancer.neuron.EAtlas2.GR)
findOverlaps.loop.BXH6.neuron.EAtlas2.GR <- findOverlaps(loop.BXH6.GR, df.enhancer.neuron.EAtlas2.GR)
findOverlaps.loop.Bn.Lx.neuron.EAtlas2.GR <- findOverlaps(loop.Bn.Lx.GR, df.enhancer.neuron.EAtlas2.GR)

hits.SHR.OlaIpcv <- findOverlaps.loop.SHR.OlaIpcv.neuron.EAtlas2.GR
hits.HXB10.BB <- findOverlaps.loop.HXB10.BB.neuron.EAtlas2.GR
hits.HXB10.CC <- findOverlaps.loop.HXB10.CC.neuron.EAtlas2.GR
hits.F344.Stm <- findOverlaps.loop.F344.Stm.neuron.EAtlas2.GR
hits.LE.Stm <- findOverlaps.loop.LE.Stm.neuron.EAtlas2.GR
hits.BXH6 <- findOverlaps.loop.BXH6.neuron.EAtlas2.GR
hits.Bn.Lx <- findOverlaps.loop.Bn.Lx.neuron.EAtlas2.GR

idx.SHR.OlaIpcv <- unique(subjectHits(hits.SHR.OlaIpcv))
idx.SHR.OlaIpcv

# sugg1: unique number loop UPSET (barchart - like multiple vendiagram)
# sugg2: At the end, one + the other end in a loop (each samples -> combine to find unique)
# sugg3: indivisual as it is. we can give weight on loops over multiple samples.



data <- data.frame(
  A = c(1, 1, 0, 1, 0, 1, 0, 1, 1, 0),
  B = c(1, 0, 1, 0, 1, 1, 0, 0, 1, 0),
  C = c(1, 0, 1, 0, 0, 0, 1, 0, 1, 1)
)

upset(data, order.by = "freq", sets.bar.color = "#56B4E9", keep.order = TRUE, main.bar.color = "#D55E00")

################################################################ hold for figure ####################################################
## figure1
vec <- c(421468530,	338261639, 333215722,	236756631,	241132264)
fig.1$seq_RP <- vec

fig.1 %>% 
  mutate(loop = n) %>% 
  select(sample, loop, seq_RP ) %>% 
  view()

fig.1 %>% 
  ggplot(aes(x=seq_RP, y=n, color=sample)) +
  labs(x = "number of sequenced read pairs", y = "number of loops") +
  geom_point(size = 3)
################################################################ hold for figure #####################################################



