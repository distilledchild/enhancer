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
library("scales")  # for better axis formatting
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
options(tibble.max_extra_cols = Inf)
options(scipen = 999)

getwd()

# Linux
# setwd('C:\\Users\\panju\\Dropbox (UTHSC GGI)\\Gateway_to_Hao\\workshop\\2023_NIH_meeting\\loop_N_tss') 
# setwd('./Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss')
# setwd('~/Desktop/temp/enhancer/dropbox_enhancer_doosan/data/enhancer_atlas2.0/all_species/neuron')
setwd("~/Desktop/temp/enhancer/dropbox_enhancer_doosan")
setwd("/home/pkim/dropbox/Gateway_to_Hao/enhancer/r_files")

getwd()

# Mac
setwd("~/dropbox/Gateway_to_Hao/enhancer/r_files")
getwd()
source("utils_functions.R")  # Load all utility functions

####################################
# Sequencing stats ###### Figure 1.a
####################################
seq.data <- read.table("../data/library_complexity.tsv", header = TRUE, sep = "\t")
seq.data[, -1] <- lapply(seq.data[, -1], function(x) as.numeric(as.character(x)))
colnames.of.seq.data <- colnames(seq.data)
colnames.of.seq.data
seq.data <- seq.data %>% mutate(Duplicates = PCR_Duplicates + Optical_Duplicates,
                                Chimeric_ambiguous_and_Unmmapped = Chimeric_Ambiguous + Unmapped)
seq.data %>% head(15)
seq.data %>% dplyr::select(Strain)

seq.data$Unique_Reads_Percentage <- (seq.data$Unique_Reads / seq.data$Sequenced_RP) * 100
seq.data$Duplicates_Percentage <- (seq.data$Duplicates / seq.data$Sequenced_RP) * 100
seq.data$Chimeric_ambiguous_and_Unmapped_Percentage <- (seq.data$Chimeric_ambiguous_and_Unmmapped / seq.data$Sequenced_RP) * 100

seq.data_melted <- melt(seq.data, id.vars = "Strain",
                        measure.vars = c("Unique_Reads_Percentage", 
                                         "Duplicates_Percentage", 
                                         "Chimeric_ambiguous_and_Unmapped_Percentage"),
                        variable.name = "Category", value.name = "Percentage")

# 4. Category order
seq.data_melted$Category <- factor(seq.data_melted$Category, 
                                   levels = c("Unique_Reads_Percentage", 
                                              "Duplicates_Percentage", 
                                              "Chimeric_ambiguous_and_Unmapped_Percentage"),
                                   labels = c("Unique Reads", "Duplicates", "Chimeric Ambiguous + Unmapped"))
seq.data_melted

summary_stats <- seq.data_melted %>%
  group_by(Category) %>%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE)
  )

summary_stats

# 5. Stacked Bar Plot
sequencing_basic_stats <- ggplot(seq.data_melted, aes(x = Strain, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = "fill") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "a. Read Category",
    x = "Sample",
    y = "Percent of Total"
  ) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) + 
  scale_fill_manual(
    values = c(
      "Unique Reads" = "#4269e1",
      "Duplicates" = "#f8a607",
      "Chimeric Ambiguous + Unmapped" = "#1a1a1a"
    )
  )

sequencing_basic_stats

saving_plot_dual( # utils_functions.R
  plot_obj = sequencing_basic_stats,
  filename_base = "sequencing_basic_stats_F1a",
  output_dir = "./figures/submission/lt2mb",
)
##############################################################################################
## LOOP processing: BASIC
##############################################################################################

# BEDPE file list (10 files)
# Linux
# loop.file.list = fs::dir_ls("/home/pkim/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/", regexp = ".bedpe$")
# loop.file.list = fs::dir_ls("/home/pkim/Desktop/temp/enhancer/dropbox_enhancer_doosan/data/loops/", regexp = ".bedpe$")
loop.file.list = fs::dir_ls("/home/hao/Dropbox/Gateway_to_Hao/enhancer/data/loops", regexp = ".bedpe$")
loop.file.list

# Mac
# loop.file.list = fs::dir_ls("~/dropbox/K\ P/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/", regexp = ".bedpe$")
# loop.file.list = fs::dir_ls("/Users/PanjunKim/Desktop/temp/enhancer/dropbox_enhancer_doosan/data/loops/", regexp = ".bedpe$")
# loop.file.list = fs::dir_ls("~/Desktop/temp/enhancer/dropbox_enhancer_doosan/data/loops/", regexp = ".bedpe$")
loop.file.list = fs::dir_ls("~/dropbox/Gateway_to_Hao/enhancer/data/loops", regexp = ".bedpe$")
loop.file.list

# BED fild for loops
df.init.loop.bed <- init.bedpe.df(loop.file.list, 'loops') %>%  # utils_functions.R
  filter(!str_detect(X.chr1, "^#"))
# %>% #58,992
#   view()

df.init.loop.bed %>% 
  count(sample)

# sample    n        strain    n
# 1   592BB 5263   SHR/OlaIpcv 5263
# 2     607 7336         HXB10 7336
# 3    74AA 2903      F344/Stm 2903
# 4    A2DB 2992        LE/Stm 2992
# 5   D765A 6568          BXH6 6568
# 6   DA08A 4656          HXB2 4656
# 7   DA21A 5932          SHR/OlaIpcvxBN/NHsdMcwi 5932 
# 8   DA68A 9131         HXB31 9131
# 9   DBA9A 7676         HXB23 7676
# 10  DE8BA 6535         BN-Lx 6535

# creating loop.id, sample.loop.id 
df.loop.deep.sample.all <- df.init.loop.bed %>% 
  mutate(strain = case_when(
    sample == '592BB' ~ "SHR/OlaIpcv",
    sample == '607' ~ "HXB10",
    sample == '74AA' ~ "F344/Stm",
    sample == 'A2DB' ~ "LE/Stm",
    sample == 'D765A' ~ "BXH6",
    sample == 'DE8BA' ~ "BN-Lx",
    sample == 'DBA9A' ~ "HXB23",
    sample == 'DA21A' ~ "SHR/OlaIpcvxBN/NHsdMcwi",
    sample == 'DA08A' ~ "HXB2",
    sample == 'DA68A' ~ "HXB31",
    TRUE ~ NA
  )) %>% 
  mutate(end.distance = x2 - x1) %>% 
  mutate(resolution = convert_to_resolution(end.distance)) %>%  # Using utility function # utils_functions.R
  mutate(loop.id = str_c(X.chr1, '_', x1, '_', x2, '_', chr2, '_', y1, '_', y2, '_', end.distance)) %>% # loop.id
  mutate(sample.loop.id = str_c(strain, '_', loop.id)) %>% #sample.loop.id
  dplyr::select(sample, strain, X.chr1, x1, x2, chr2, y1, y2, distance, end.distance, resolution, loop.id, sample.loop.id) %>% 
  dplyr::rename(chr1 = X.chr1)

df.loop.deep.sample.all # 51300 + 7692 = 58992
df.loop.deep.sample.all %>% head()
df.loop.deep.sample.all %>% count(strain)

################################################
# 1. Loop
# 1-1. Exploratory Data analysis (EDA)
# 1-1.1. data processing ############# Figure S3
################################################
# 1. checking duplication loops in a strain
df.loop.deep.sample.all %>% 
  # count(strain, loop.id) %>%             # NO dups in a strain
  count(strain, resolution, loop.id) %>% # NO dup in the same resolution in a strain
  filter(n > 1) # 0

# 2. checking how much common loops are in samples
resolutions <- c("5K", "10K", "25K")
plots <- list()

df.loop.deep.sample.all %>% head(1)

for(res in resolutions) {
  message("Start processing: ", res)
  df.loop.deep.sample.all.res <- df.loop.deep.sample.all %>% filter(resolution == res)
  
  loop_counts <- df.loop.deep.sample.all.res %>% group_by(strain) %>% summarise(total_loops = n_distinct(loop.id))
  
  location_list <- split(df.loop.deep.sample.all.res$loop.id, df.loop.deep.sample.all.res$strain)
  
  message("location_list created with ", length(location_list), " strains.")
  
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
  common_counts_long <- reshape2::melt(common_counts)
  message("Melted common_counts into long format: ", nrow(common_counts_long), " rows.")
  
  # heatmap
  p <- ggplot(common_counts_long, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.1f%%", value)), color = "black") +
    scale_fill_gradient(low = "white", high = "darkred", 
                        name="Common Loop %", limits = c(0, 100)) +
    theme_minimal() +
    labs(x = "Strain", y = "Strain", title = paste("Heatmap of Common Loop Percentage -", res, "Resolution")) +
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

p

# all figures by resolution into one figure
common_combined_plot <- arrangeGrob(
  plots[["5K"]],
  plots[["10K"]],
  plots[["25K"]],
  nrow = 3
)

saving_plot_dual( # utils_functions.R
  plot_obj = common_combined_plot,
  filename_base = "common_loops_heatmap_percentage_bw_strains_S3",
  output_dir = "./figures/submission/lt2mb",
  width_in = 11*1.6,
  height_in = 8.5*3*1.6, 
  dpi = 300
)

########################################################################
# 1. Loop
# 1-1-2. figure: shared loops - bar plot ##################### Figure 3a
########################################################################
# checking shared loops
# Step 1: shared loop
df.loop.deep.sample.all %>% head()

shared_loops <- df.loop.deep.sample.all %>%
  group_by(resolution, loop.id) %>%
  summarise(n_samples = n_distinct(sample), .groups = "drop") %>%
  filter(n_samples > 1) %>%  # loops shared by more than 2 sample
  inner_join(df.loop.deep.sample.all, by = c("resolution", "loop.id"))
shared_loops
# Step 2: shared loop by resolution by sample
sample_loop_counts <- shared_loops %>%
  group_by(resolution, sample) %>%
  summarise(shared_loop_count = n_distinct(loop.id), .groups = "drop")

# Step 3: mean & SD per resolution
plot_df <- sample_loop_counts %>%
  group_by(resolution) %>%
  summarise(
    mean_shared_loops = mean(shared_loop_count),
    sd_shared_loops = sd(shared_loop_count),
    .groups = "drop"
  )

plot_df

# Step 4: bar plot
shared_loops_by_resolution <- ggplot(plot_df, aes(x = resolution, y = mean_shared_loops, fill = resolution)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_errorbar(aes(ymin = mean_shared_loops - sd_shared_loops,
                    ymax = mean_shared_loops + sd_shared_loops),
                width = 0.2, color = "black") +
  labs(
    title = "a. Shared Loop by Resolution",
    x = "Resolution",
    y = "Shared Loops (Average)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),  
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("5K" = "#a6cee3", "10K" = "#1f78b4", "25K" = "#1f3a93"))

shared_loops_by_resolution

saving_plot_dual( # utils_functions.R
  plot_obj = shared_loops_by_resolution,
  filename_base = "shared_loops_by_resolution_F3a",
  output_dir = "./figures/submission/lt2mb",
)
########################################################################
# 1. Loop
# 1-1-2. figure: shared loops - network plot ################# Figure S2
########################################################################
figure_green <- "#00573F"
figure_orange <- "#FFA300"

# Generate strain pairs based on shared loops
strain_pairs <- shared_loops %>%
  dplyr::select(loop.id, strain) %>%
  distinct() %>%
  group_by(loop.id) %>%
  summarise(strains = list(unique(strain)), .groups = "drop") %>%
  mutate(pairs = map(strains, ~ combn(.x, 2, simplify = FALSE))) %>%
  dplyr::select(pairs) %>%
  unnest(pairs) %>%
  mutate(from = map_chr(pairs, 1),
         to   = map_chr(pairs, 2)) %>%
  count(from, to, name = "weight")  # number of shared loops

strain_pairs

# Convert to igraph
network_graph <- graph_from_data_frame(strain_pairs, directed = FALSE)

# Plot
network_plot_for_shared_loops <- ggraph(network_graph, layout = "fr") +  # fr: force-directed layout
  geom_edge_link(aes(width = weight), alpha = 0.7, color = figure_green) +
  geom_node_point(size = 5, color = figure_orange) +
  geom_node_text(aes(label = name), repel = TRUE, size = 4, fontface = "bold") +
  scale_edge_width(range = c(0.5, 3)) +
  theme_void()

network_plot_for_shared_loops

saving_plot_dual( # utils_functions.R
  plot_obj = network_plot_for_shared_loops,
  filename_base = "network_plot_for_shared_loops_S2",
  output_dir = "./figures/submission/lt2mb"
)

################################################################################################
# 1. Loop
# 1-1-4. figure: loops by sequencing reads ########################################### Figure 1b
################################################################################################
# Step 1: loops by strain: df.loop.deep.sample.all
loop_counts_by_sample <- df.loop.deep.sample.all %>%
  count(strain) %>% 
  dplyr::rename(num_loop = n)
loop_counts_by_sample
seq.data
seq.data <- seq.data %>% dplyr::rename(strain = Strain)
seq.data

# Step 2: sequencing info & join 
merged_df <- loop_counts_by_sample %>%
  inner_join(seq.data, by = "strain")

merged_df %>% head(5)

# Step 3: long-format for plot
df_long_for_plot <- merged_df %>%
  mutate(
    Sequenced_RP = as.numeric(gsub(",", "", Sequenced_RP)),
    Unique_Reads = as.numeric(gsub(",", "", Unique_Reads)),
    Alignable_Normal_N_Chimeric = as.numeric(gsub(",", "", Alignable_Normal_N_Chimeric))
  ) %>%
  dplyr::select(strain, num_loop,
                "Total Reads" = Sequenced_RP,
                "Unique Reads" = Unique_Reads,
                "Alignable Reads" = Alignable_Normal_N_Chimeric) %>%
  pivot_longer(cols = c("Total Reads", "Unique Reads", "Alignable Reads"),
               names_to = "Sequencing_Metric",
               values_to = "Depth")

df_long_for_plot

# r & p-value
correlation_results <- df_long_for_plot %>%
  group_by(Sequencing_Metric) %>%
  group_modify(~ cor.test(.x$Depth, .x$num_loop) %>% broom::tidy()) %>%
  ungroup() %>%
  mutate(
    r = format(round(estimate, 2), nsmall = 2),
    p = format(round(p.value, 4), nsmall = 4),
    label = paste0("R = ", r, ", p = ", p)
  )

correlation_results
#   Sequencing_Metric estimate statistic  p.value parameter conf.low conf.high
#   <chr>                <dbl>     <dbl>    <dbl>     <int>    <dbl>     <dbl>
# 1 Alignable Reads      0.787      3.61 0.00689          8    0.312     0.947
# 2 Total Reads          0.780      3.52 0.00783          8    0.295     0.945
# 3 Unique Reads         0.884      5.35 0.000688         8    0.574     0.972

annot_positions <- data.frame(
  Sequencing_Metric = c("Alignable Reads", "Total Reads", "Unique Reads"),
  x = rep(min(df_long_for_plot$Depth) * 1.05, 3),
  y = c(
    max(df_long_for_plot$num_loop) * 0.97,
    max(df_long_for_plot$num_loop) * 0.92,
    max(df_long_for_plot$num_loop) * 0.87
  )
)

# correlation_results와 병합
annot_df <- left_join(correlation_results, annot_positions, by = "Sequencing_Metric")

# Original labels
original_labels <- levels(factor(df_long_for_plot$Sequencing_Metric))
original_labels

line_graph_for_loops_per_depth <- ggplot(df_long_for_plot, aes(x = Depth, y = num_loop, color = Sequencing_Metric)) +
  geom_point() +
  geom_smooth(aes(color = Sequencing_Metric, fill = Sequencing_Metric), method = "lm", se = TRUE, linewidth = 1, linetype = "solid") +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  labs(
    title = "b. Read Counts vs Loop Counts by Category",
    x = "Number of Reads",
    y = "Number of Loops",
    color = "Category",
    fill = "Category" # Add fill to labs
  ) +
  scale_color_discrete(labels = original_labels) + # Use the combined labels
  scale_fill_discrete(labels = original_labels) + # and for fill as well
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  ) +
  geom_text(
    data = annot_df,
    aes(x = x, y = y, label = label, color = Sequencing_Metric),
    hjust = 0,
    size = 4.5,
    fontface = "italic"
  )

line_graph_for_loops_per_depth 

saving_plot_dual( # utils_functions.R
  plot_obj = line_graph_for_loops_per_depth,
  filename_base = "line_graph_for_loops_per_depth_hao_w_new_label_F1b",
  output_dir = "./figures/submission/lt2mb"
)

figure1_combined <- sequencing_basic_stats + line_graph_for_loops_per_depth +
  plot_layout(ncol = 2, nrow = 1) 

figure1_combined

saving_plot_dual( # utils_functions.R
  plot_obj = figure1_combined,
  filename_base = "figure1_combined_F1",
  output_dir = "./figures/submission/lt2mb",
  height = 5.5,        ############ NOT 8.5
  scale_x = 1,
  scale_y = 1
)

########################################################################
# 1. Loop
# 1-1-6. figure: loops per chromosome by resolution ########### Figure 2
########################################################################
df_chr_loop_counts <- df.loop.deep.sample.all %>%
  group_by(chr1, resolution) %>%
  summarise(n_loops = n_distinct(loop.id), .groups = "drop")

# chr factor (chr1, chr2, ..., chrX)
df_chr_loop_counts$chr1 <- factor(df_chr_loop_counts$chr1, levels = mixedsort(unique(df_chr_loop_counts$chr1)))

df_chr_loop_counts_fig <- ggplot(df_chr_loop_counts, aes(x = chr1, y = n_loops, fill = resolution)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Chromosome",
    y = "Number of Loops"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_manual(values = c("5K" = "#a6cee3", "10K" = "#1f78b4", "25K" = "#1f3a93"))

df_chr_loop_counts_fig

saving_plot_dual( # utils_functions.R
  plot_obj = df_chr_loop_counts_fig,
  filename_base = "loop_counts_per_chr_F2",
  output_dir = "./figures/submission/lt2mb"
)

########################
# 1. Loop
# 1-2-1. Loop data preprocessing: distinct loops
########################

chromosome_data <- read.table(file="../data/rn7_chromosome_length_from_ucsc.tsv", sep="\t") %>%
  mutate(start = 0) %>%
  dplyr::rename(chr = V1, end = V2)
chromosome_data
chromosome_data %>% head()

df.loop.deep.sample.all %>% dim() # 58992   12
df.loop.deep.sample.all %>% head(3)

# common loops between samples picked only one: df.DISTINCT.loop.deep.sample.all
df.DISTINCT.loop.deep.sample.all <- df.loop.deep.sample.all %>%
  dplyr::select(loop.id) %>% # line 180
  distinct() %>% 
  separate(loop.id, into = c("chr1", "x1", "x2", "chr2", "y1", "y2", "end.distance"), sep = "_", remove = FALSE, convert = TRUE) %>% 
  mutate(distance = y2 - x2) %>%
  mutate(resolution = convert_to_resolution(end.distance)) %>%  # Using utility function # utils_functions.R
  mutate(x0 = x1, x3 = x2, y0 = y1, y3 = y2) %>%                                # x0, x3, y0, y3
left_join(chromosome_data, by = c('chr1' = 'chr')) %>%                          # chromosome length info
dplyr::rename(chr.end.coord = end)

df.DISTINCT.loop.deep.sample.all %>% dim() # 31773    11
df.DISTINCT.loop.deep.sample.all %>% head(3)

df.DISTINCT.loop.deep.sample.all %>% 
  count(resolution) # 31773/58992
# resolution     n
# 1         5K  6680
# 2        10K 12162
# 3        25K 12931

########################
# 1. Loop
# 1-3. Loop data preprocessing: Creating GRange Obj. from df.DISTINCT.loop.deep.sample.all (whole, up, down)
########################

########################
# 1. Loop
# 1-4. Loop data preprocessing: padding on loops (1 distance) 
# -> OVERALL.df.DISTINCT.loop.deep.sample.all : OBJECT to be used for the distribution of sth over loops
# x0, y3, new_distance, new.loop.id
########################

df.DISTINCT.loop.deep.sample.all %>% dim() # 31773
df.DISTINCT.loop.deep.sample.all %>% head(2) 

# loops only less than 2mb: 31019
df.DISTINCT.loop.deep.sample.all.lt.2mb <- df.DISTINCT.loop.deep.sample.all %>% 
  filter(distance < 2000000) # 31019, only use less than 2mb
df.DISTINCT.loop.deep.sample.all.lt.2mb %>% dim() # 31019 (31773 - 754 (longer than 2mb))
df.DISTINCT.loop.deep.sample.all.lt.2mb %>% head(3)

################################################
# padding 1 distance with all loops: 31773
################################################
OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance <- df.DISTINCT.loop.deep.sample.all %>% 
  mutate(x_mid = (x1 + x2)/2, y_mid = (y1 + y2)/2) %>% # middle point of each end
  mutate(x0 = ifelse(x_mid - distance < 0, 0, x_mid - distance), 
         y3 = ifelse(y_mid + distance > chr.end.coord, chr.end.coord, y_mid + distance))   # 1 distance for padding

OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance %>% dim() # 31773
################################################
# padding 1 distance && w/o capping only from all loops: 31417
################################################
OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping <- OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance %>% 
  filter(!(x0 == 0 | chr.end.coord == y3)) # 31773 - 444 (282 + 162) = 31329 + 88 = 31417
# filter(x0 == 0) # 282
# filter(chr.end.coord == y3)  # 162
OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping %>% dim() # 31417
OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance %>% filter(x0 == 0 & chr.end.coord == y3) # 88 | 

################################################
# padding 1 distance && w/o capping && less than 2mb: 31329
################################################
OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping.lt.2mb <- OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance %>% 
  filter(distance < 2000000)   %>% # 31019                          ################ distance is the original loop distance (prior to padding)
  filter(!(x0 == 0 | chr.end.coord == y3)) #%>% # 31773 - 444 = 31329, 88// # 30928
  count(loop.id)  %>% view()
  # filter(x0 == 0) # 282
  # filter(chr.end.coord == y3)  # 162
OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping.lt.2mb %>% dim() # 30928

df.DISTINCT.loop.deep.sample.all.lt.2mb.stats <- df.DISTINCT.loop.deep.sample.all.lt.2mb %>%
  summarise(
    Q1 = quantile(distance, 0.25), # 150000
    Median = median(distance), # 195000
    Q3 = quantile(distance, 0.75) # 375000
  )

ggplot(df.DISTINCT.loop.deep.sample.all.lt.2mb, aes(y = distance)) +
  geom_boxplot(fill = "skyblue", color = "darkblue") +
  # Q1, Median, Q3
  geom_text(data = df.DISTINCT.loop.deep.sample.all.lt.2mb.stats, aes(x = 0, y = Q1, label = paste0("Q1: ", round(Q1))), hjust = 0) +
  geom_text(data = df.DISTINCT.loop.deep.sample.all.lt.2mb.stats, aes(x = 0, y = Median, label = paste0("Median: ", round(Median))), hjust = 0) +
  geom_text(data = df.DISTINCT.loop.deep.sample.all.lt.2mb.stats, aes(x = 0, y = Q3, label = paste0("Q3: ", round(Q3))), hjust = 0) +
  labs(
    y = "Loop Length (bp)",
    title = "Distribution of Loop Length"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)  # title centering
  )

########################
# 2. CTCF
# 2-1. Exploratory Data analysis (EDA) : CTCF
########################
# Linux
# df.init.ctcf<-read.table(file="/home/pkim/Desktop/temp/enhancer/dropbox_enhancer_doosan/data/ctcf/fimo_2nd_trial_meme/fimo_2nd_trial.txt", header=TRUE, sep="\t")
df.init.ctcf<-read.table(file="../data/ctcf/fimo_2nd_trial_meme/fimo_2nd_trial.txt", header=TRUE, sep="\t")
df.init.ctcf %>% count() # 5767921

# Mac
# ctcf<-read.table(file="/home/pkim/Desktop/temp/enhancer/dropbox_enhancer_doosan/data/ctcf/fimo_2nd_trial_meme/fimo_2nd_trial.txt", sep="\t", col.names = c("chr", "start", "end", "strand", "length"), header = FALSE) %>% 
# df.init.ctcf<-read.table(file="~/dropbox/K P/Gateway_to_Hao/enhancer/data/ctcf/fimo_2nd_trial_meme/fimo_2nd_trial.txt", header=TRUE, sep="\t")
df.init.ctcf<-read.table(file="~/dropbox/Gateway_to_Hao/enhancer/data/ctcf/submission/E4/fimo_E4_submission_trial.txt", header=TRUE, sep="\t") %>% 
  dplyr::rename(chr = sequence_name, end = stop) %>%
  mutate(length = end - start)

df.init.ctcf %>% dim() # .4:6551641
df.init.ctcf %>% head(3) # chr    start      end strand length
df.init.ctcf %>% count(length) # length distribution: 9 ~ 39

# strand checking
df.init.ctcf %>%                              # +: 2891072, -: 2876849 = 5767921//+:3267352, -:3284289
  count(strand)
# removing dups including strand
df.init.ctcf %>% 
  distinct(chr, start, end, strand)           # .4: 3331233/6551641
# removing dups excluding strand
df.init.ctcf %>% 
  distinct(chr, start, end)                   # .4: 3191859/6551641
# removing dups with all including length
df.init.ctcf %>% 
  distinct()                                  # .4: 3331233/6551641

########################
# 2. CTCF
# 2-2. CTCF data preprocessing: id and dedup GRange Obj. (df.DISTINCT.fimo.2nd.trial.ctcf/ df.DISTINCT.ctcf.2nd.fimo.GR)
########################
# generating id column (long running time)
# df.init.ctcf %>% head() 6551641
# chr     start       end strand length

df.DISTINCT.fimo.2nd.trial.ctcf <- df.init.ctcf %>%
  distinct(chr, start, end) %>% # .4:3191859 ************** NO STRAND INFO
  mutate(start = as.numeric(start)) %>% 
  mutate(end = as.numeric(end)) %>% 
  mutate(ctcf_pos = as.numeric(round((start + end) / 2))) %>% 
  mutate(id = str_c(chr, "_", start, "_", end, "_", ctcf_pos))

df.DISTINCT.fimo.2nd.trial.ctcf %>% dim() # .4: 3191859
df.DISTINCT.fimo.2nd.trial.ctcf # .4 : 3191859/6551641 :0.4871847
df.DISTINCT.fimo.2nd.trial.ctcf %>% head()  # chr    start      end ctcf_pos                               id
df.DISTINCT.fimo.2nd.trial.ctcf %>% count(chr)

# GRanges Obj.: df.DISTINCT.ctcf.2nd.fimo.GR
df.DISTINCT.ctcf.2nd.fimo.GR <- GRanges(
  seqnames = as.character(df.DISTINCT.fimo.2nd.trial.ctcf$chr),
  ranges = IRanges(start = df.DISTINCT.fimo.2nd.trial.ctcf$start, 
                   end = df.DISTINCT.fimo.2nd.trial.ctcf$end)
)

# metadata: id
mcols(df.DISTINCT.ctcf.2nd.fimo.GR)$id <- df.DISTINCT.fimo.2nd.trial.ctcf$id
mcols(df.DISTINCT.ctcf.2nd.fimo.GR)$ctcf_pos <- df.DISTINCT.fimo.2nd.trial.ctcf$ctcf_pos

########################
# 2. CTCF
# 2-3. overall distribution of CTCF on loops
########################
OVERALL.df.DISTINCT.loop.deep.sample.all <- OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping.lt.2mb ###################### 30928    16

OVERALL.df.DISTINCT.loop.deep.sample.all %>% head(3)
OVERALL.df.DISTINCT.loop.deep.sample.all.GR <- creating_granges(OVERALL.df.DISTINCT.loop.deep.sample.all) # utils_functions.R

index.distinct.ctcf.w.OVERALL.whole.loop <- findOverlaps(
  df.DISTINCT.ctcf.2nd.fimo.GR, 
  OVERALL.df.DISTINCT.loop.deep.sample.all.GR, 
  type = "any",
  select = "all"
)

OVERALL.loop.for.ctcf.hits <- subjectHits(index.distinct.ctcf.w.OVERALL.whole.loop)
OVERALL.ctcf.on.loop.hits <- queryHits(index.distinct.ctcf.w.OVERALL.whole.loop)

df.ctcf.dist.result <- tibble(
  loop.id = OVERALL.df.DISTINCT.loop.deep.sample.all$loop.id[OVERALL.loop.for.ctcf.hits],
  loop.start = OVERALL.df.DISTINCT.loop.deep.sample.all$x0[OVERALL.loop.for.ctcf.hits],
  loop.end = OVERALL.df.DISTINCT.loop.deep.sample.all$y3[OVERALL.loop.for.ctcf.hits],
  loop.res = OVERALL.df.DISTINCT.loop.deep.sample.all$resolution[OVERALL.loop.for.ctcf.hits],
  ctcf.id = df.DISTINCT.fimo.2nd.trial.ctcf$id[OVERALL.ctcf.on.loop.hits],
  ctcf.pos = df.DISTINCT.fimo.2nd.trial.ctcf$ctcf_pos[OVERALL.ctcf.on.loop.hits],
) %>% 
  mutate(across(where(is.numeric), ~ format(., scientific = FALSE)))

df.ctcf.dist.result %>% dim() # 40250853(any, w/o capping + lt 2mb)

# converting data types
relative.pos.df.ctcf.dist.result <- df.ctcf.dist.result %>%
  mutate(loop.start = as.numeric(loop.start),
         loop.end = as.numeric(loop.end),
         pos_coord = as.numeric(ctcf.pos),
         loop_length = (loop.end - loop.start),
         relative_pos = pos_coord - loop.start,
         value = relative_pos / (loop_length / 3) - 1
  ) %>%
  dplyr::select(loop.id, ctcf.id, value, loop_length, loop.res)

relative.pos.df.ctcf.dist.result %>% head(3)

########################
# 2. CTCF
# 2-3. overall distribution of CTCF on loops: figures
# 2-3-1. by CHROMOSOME
########################
plot.ctcf.hist <- plot_histogram(relative.pos.df.ctcf.dist.result) # utils_functions.R
plot.ctcf.dens <- plot_density(relative.pos.df.ctcf.dist.result) # utils_functions.R

saving_combined_plot(plot.ctcf.hist, plot.ctcf.dens, # utils_functions.R
                   "figures/submission/lt2mb/overall_distribution_of_CTCF_by_resolution_wo_capping_lt2mb.pdf")

########################
# 2. CTCF
# 2-4. Distribution of CTCF at ends in a loop: for the number of CTCF used in filtering valid loops: figures
# so, the object should be used one without padding: df.DISTINCT.loop.deep.sample.all
# 2-4-3. data processing for getting information of CTCF at ends in loops 
########################
df.DISTINCT.loop.deep.sample.all.lt.2mb # 31019
df.DISTINCT.loop.deep.sample.all.up.GR   <- creating_granges(df.DISTINCT.loop.deep.sample.all.lt.2mb, direction = "up") # utils_functions.R
df.DISTINCT.loop.deep.sample.all.up.GR # 31773/ 31019

# 1. CTCF + UPSTREAM (df.DISTINCT.loop.deep.sample.all.lt.2mb, df.DISTINCT.fimo.2nd.trial.ctcf)
index.distinct.ctcf.w.up.loop <- findOverlaps(df.DISTINCT.ctcf.2nd.fimo.GR, 
                                              df.DISTINCT.loop.deep.sample.all.up.GR, 
                                              type = "any",
                                              select = "all")

end.loop.up.ctcf.hits <- subjectHits(index.distinct.ctcf.w.up.loop)
end.ctcf.up.hits <- queryHits(index.distinct.ctcf.w.up.loop)

df.overlapping.CTCF.w.UPSTREAM.result <- tibble(
  up.loop.id = mcols(df.DISTINCT.loop.deep.sample.all.up.GR)$loop.id[end.loop.up.ctcf.hits],
  end.up.distance = mcols(df.DISTINCT.loop.deep.sample.all.up.GR)$end.distance[end.loop.up.ctcf.hits],
  distance = mcols(df.DISTINCT.loop.deep.sample.all.up.GR)$distance[end.loop.up.ctcf.hits],
  resolution = mcols(df.DISTINCT.loop.deep.sample.all.up.GR)$resolution[end.loop.up.ctcf.hits],
  ctcf.id = mcols(df.DISTINCT.ctcf.2nd.fimo.GR)$id[end.ctcf.up.hits],
  WHERE = "UP"
) %>% 
  unite(ctcf.loop.up.id, up.loop.id, distance, ctcf.id, WHERE, sep = "|", remove = FALSE)

df.overlapping.CTCF.w.UPSTREAM.result %>% dim() # sub.4 any: 1402928| .4 lt2mb 1381954
df.overlapping.CTCF.w.UPSTREAM.result %>% head()

df.overlapping.CTCF.w.UPSTREAM.result %>% count(end.up.distance)

# sub.4
# ANY
# end.up.distance      n              # sub.4 any lt2mb
# 1            5000 179049  # 1            5000 178245
# 2           10000 444448  # 2           10000 442386
# 3           25000 779431  # 3           25000 761323

df.DISTINCT.loop.deep.sample.all.down.GR   <- creating_granges(df.DISTINCT.loop.deep.sample.all.lt.2mb, direction = "down") # utils_functions.R
df.DISTINCT.loop.deep.sample.all.down.GR

# 2. CTCF + DOWNSTREAM (df.DISTINCT.loop.deep.sample.all.lt.2mb, df.DISTINCT.fimo.2nd.trial.ctcf)
index.distinct.ctcf.w.down.loop <- findOverlaps(df.DISTINCT.ctcf.2nd.fimo.GR, 
                                                df.DISTINCT.loop.deep.sample.all.down.GR, 
                                                type = "any",
                                                select = "all")

end.loop.down.ctcf.hits <- subjectHits(index.distinct.ctcf.w.down.loop)
end.ctcf.down.hits <- queryHits(index.distinct.ctcf.w.down.loop)

df.overlapping.CTCF.w.DOWNSTREAM.result <- tibble(
  down.loop.id = mcols(df.DISTINCT.loop.deep.sample.all.down.GR)$loop.id[end.loop.down.ctcf.hits],
  end.down.distance = mcols(df.DISTINCT.loop.deep.sample.all.down.GR)$end.distance[end.loop.down.ctcf.hits],
  distance = mcols(df.DISTINCT.loop.deep.sample.all.down.GR)$distance[end.loop.down.ctcf.hits],
  resolution = mcols(df.DISTINCT.loop.deep.sample.all.down.GR)$resolution[end.loop.down.ctcf.hits],
  ctcf.id = mcols(df.DISTINCT.ctcf.2nd.fimo.GR)$id[end.ctcf.down.hits],
  WHERE = "DOWN"
) %>% 
  unite(ctcf.loop.down.id, down.loop.id, distance, ctcf.id, WHERE, sep = "|", remove = FALSE)

df.overlapping.CTCF.w.DOWNSTREAM.result %>% dim() # sub.4 any: 1411588 | .4lt2mb 1393938    
df.overlapping.CTCF.w.DOWNSTREAM.result %>% count(end.down.distance)
# sub.4lt2mb                              # sub.4
# end.down.distance      n                # end.down.distance      n
# any                                     # any
# <int>  <int>                            # <int>  <int>
#1              5000 182341 # 1              5000 182974
#2             10000 449962 # 2             10000 452423
#3             25000 761635 # 3             25000 776191

########## bind_rows(UPSTREAM & DOWNSTREAM) -> BOTH
df.overlapping.CTCF.w.BOTH.result <- bind_rows(df.overlapping.CTCF.w.UPSTREAM.result %>% 
                                                 mutate(loop.id = up.loop.id, end.distance = end.up.distance) %>% 
                                                 mutate(case.id = ctcf.loop.up.id) %>% 
                                                 dplyr::select(-c(up.loop.id, end.up.distance, ctcf.loop.up.id)), 
                                               df.overlapping.CTCF.w.DOWNSTREAM.result %>% 
                                                 mutate(loop.id = down.loop.id, end.distance = end.down.distance) %>% 
                                                 mutate(case.id = ctcf.loop.down.id) %>% 
                                                 dplyr::select(-c(down.loop.id, end.down.distance, ctcf.loop.down.id))) %>% 
  mutate(chr = ifelse(WHERE == "UP", str_split_n(loop.id, '_', 1), str_split_n(loop.id, '_', 4))) %>%  # utils_functions.R
  mutate(chr = factor(chr, levels = c(paste0("chr", 1:20), "chrX", "chrY"))) %>%
  mutate(WHERE = fct_relevel(WHERE, "UP", "DOWN")) %>% 
  mutate(resolution = convert_to_resolution(end.distance))  # Using utility function # utils_functions.R

df.overlapping.CTCF.w.BOTH.result %>% dim() # sub.4 any: 2814516 |.4 lt2mb 2775892 
df.overlapping.CTCF.w.BOTH.result %>% colnames() # "distance" "resolution" "ctcf.id" "WHERE" "loop.id" "end.distance" "case.id" "chr"
df.overlapping.CTCF.w.BOTH.result %>% count(resolution)
#  .4 lt2mb            # any
# 1 5K          360586 # 1 5K          362023
# 2 10K         892348 # 2 10K         896871
# 3 25K        1522958 # 3 25K        1555622

# for boxplot
df.overlapping.CTCF.w.BOTH.result %>% head(2) # distance, resolution, ctcf.id, WHERE, loop.id, end.distance, case.id, chr
df.overlapping.CTCF.w.BOTH.result %>% distinct(loop.id) # .4 any 31773 | .4 any lt2mb 30,908

# for Q1: quartile among loops with CTCF, so min = 1
df.ctcf.counts <- df.overlapping.CTCF.w.BOTH.result %>%
  group_by(loop.id, WHERE, resolution) %>%
  summarise(ctcf_count = n_distinct(ctcf.id), .groups = 'drop')

df.ctcf.counts %>% count(ctcf_count)

###############
###############
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
# case     n
# 1 both 31634
# 2   no   139
# > 31634 + 139 = 31773

# lt2mb 
# case     n
# 1 BOTH 29978
# 2 DOWN   467
# 3 NONE   111
# 4   UP   463
# 29978 + 467 + 111 + 463 = 31019 PASS

# proofreading
# check 1. DOWN
down_ctcf_ids <- df.ctcf.counts %>%
  group_by(loop.id) %>%
  filter(all(WHERE == "DOWN")) %>%
  distinct(loop.id)
down_ctcf_ids # 467

# check 2. UP
up_ctcf_ids <- df.ctcf.counts %>%
  group_by(loop.id) %>%
  filter(all(WHERE == "UP")) %>%
  distinct(loop.id)

up_ctcf_ids # 463

# 3. UP & DOWN
both_ctcf_ids <- df.ctcf.counts %>%
  group_by(loop.id) %>%
  filter(all(c("UP", "DOWN") %in% WHERE) & n_distinct(WHERE) == 2) %>%
  distinct(loop.id)

both_ctcf_ids # 29,978

# 3 cases: 30573 + 540 + 521 = 31634 + 139 (no CTCF on both ends) = TOTAL 31773
# 3 cases: 29978 + 467 + 463 = 30908 + 111 (no CTCF on both ends) = TOTAL 31019 lt2mb

###############
###############

ctcf.stats.by.resolution <- calculate_grouped_stats( # utils_functions.R
  df.ctcf.counts, 
  value_col = "ctcf_count", 
  group_cols = c("WHERE", "resolution")
)  # Using utility function
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
# 3-1. Exploratory Data analysis (EDA) : TSS
# 3-2. TSS data preprocessing: id and dedup GRange Obj.
########################
#############################
# tss resource 1: deprecated
#############################
# tss file list
# Linux
file.tss.list = fs::dir_ls("/home/pkim/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/", regexp = "\\.txt$")
# Mac
file.tss.list = fs::dir_ls("~/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss", regexp = "\\.txt$")

file.tss.list
# csRNA.NuAcc.tss.txt
# csRNA.PFC.tss.txt
# ucsc_start_codon.txt

##### tss with nuacc, 96563
# nuacc <- read.csv("/home/pkim/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/csRNA.NuAcc.tss.txt", header = T, sep = '\t')
# nuacc <- read.csv("/Users/PanjunKim/UTHSC GGI Dropbox/K P/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/csRNA.NuAcc.tss.txt", header = T, sep = '\t')

##### tss with pfc ver1, 131647
# pfc <- read.csv("/home/pkim/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/csRNA.PFC.tss.txt", header = T, sep = '\t')
# pfc <- read.csv("/Users/PanjunKim/UTHSC GGI Dropbox/K P/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/csRNA.PFC.tss.txt", header = T, sep = '\t')

##### tss with UCSC
# get the TSS location into a GenomicRange object, note gene name is added, also there are some duplicated lines, so use uniq
# Linux
tss<-read.table(file="/home/pkim/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/ucsc_start_codon.txt", sep="\t", head=F)
# Mac
tss<-read.table(file="~/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/ucsc_start_codon.txt", sep="\t", head=F)
tss # 17849
tss %>% head()
head(tss)[,c(1,4,5,7,9)]
tss.select<-tss[,c(1,4,5,7,9)] # onlty take the relevant columns
tss.select
tss.select.sparate <- tss.select %>% 
  separate(V9, into = c("gene_id", "transcript_id", "exon_number", "exon_id", "gene_name"), sep = "; ") %>%
  mutate(across(everything(), ~ gsub(".+\\s", "", .))) %>% 
  mutate(gene_name = str_replace(gene_name, ';', ''))
tss.select.sparate %>% head()
tss.select.sparate %>% dim() # 17849
# tss.select.sparate %>% filter(gene_id != gene_name) # 0 rows, so gene_id and gene_name are consistent

#############################
# tss resource 2: deprecated (EXON used)
#############################
df.refgene.gtf <- read_tsv("~/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/ucsc_refGene.gtf", # download from ucsc: https://hgdownload.soe.ucsc.edu/goldenPath/rn7/bigZips/genes/
                           comment = "#",  
                           col_names = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"),
                           col_types = cols(.default = "c"))  # char

df.refgene.gtf %>% 
  # head()
  count(feature)
# feature          n
# <chr>        <int>
# 1 3UTR         16571
# 2 5UTR         21453
# 3 CDS         164583
# 4 exon        174505
# 5 start_codon  17849 ***********
# 6 stop_codon   17806
# 7 transcript   18570

############################
# tss resource 2-1: RefSeq exon data processing
############################
#################################################### 
# retrieving exon data from RefSeq GTF
#################################################### 
df.refgene.gtf.for.exon.raw <- df.refgene.gtf %>% 
  filter(feature == "exon") %>% 
  filter(chr %in% c(paste0("chr", 1:20), "chrX", "chrY"))

# checking attribute keys
refgene.exon.attribute.keys <- get_attribute_keys(df.refgene.gtf.for.exon.raw$attribute) # utils_functions.R
refgene.exon.attribute.keys
# [1] "exon_id"       "exon_number"   "gene_id"       "gene_name"     "transcript_id"

# adding columns from attribute
df.refgene.gtf.for.exon.attribute <- df.refgene.gtf.for.exon.raw %>% 
  bind_cols(df.refgene.gtf.for.exon.raw$attribute %>% map_dfr(~extracting_attributes(.x, keys = refgene.exon.attribute.keys))) %>% # utils_functions.R
  dplyr::select(-c(source, feature, attribute, score, frame))

df.refgene.gtf.for.exon.attribute # 174,505
df.refgene.gtf.for.exon.attribute %>% # 174,505
  filter(gene_id != gene_name)

# filtering: get the last exon (largest exon_number) per gene_id
df.refgene.gtf.for.exon <- df.refgene.gtf.for.exon.attribute %>%
  mutate(exon_number = as.numeric(exon_number)) %>% # convert to numeric
  group_by(gene_id) %>% 
  slice_max(order_by = exon_number, n = 1, with_ties = FALSE) %>% # keep only the largest exon_number
  ungroup() %>% 
  mutate(refseq_exon_id = str_c(chr, ':', start, ':', end, ':', gene_name, ':', gene_id, ':', exon_number))

df.refgene.gtf.for.exon %>% dim() # 17,488
df.refgene.gtf.for.exon %>% head(3)
df.refgene.gtf.for.exon %>% add_count(gene_id) %>% filter(n > 1) # should be 0
df.refgene.gtf.for.exon


############################
# tss resource 3
############################
df.ensembl.gtf <- read_tsv("~/dropbox/Gateway_to_Hao/enhancer/data/Rattus_norvegicus.mRatBN7.2.113.gtf", # download (https://ftp.ensembl.org/pub/release-113/gtf/rattus_norvegicus/)
                comment = "#", 
                col_names = FALSE) # 1,284,446 × 9

df.ensembl.gtf %>% head()
df.ensembl.gtf %>% count(X3)
#   X3                   n
# 1 CDS             471937
# 2 Selenocysteine      25
# 3 exon            526642 ** last exon of gene
# 4 five_prime_utr   63845
# 5 gene             30562
# 6 start_codon      42990 ********
# 7 stop_codon       45157
# 8 three_prime_utr  48295
# 9 transcript       54993
df.ensembl.gtf %>% dplyr::select(X9) %>% head(3) # attribute

colnames(df.ensembl.gtf) <- c("chr", "source", "feature", "start", "end", 
                   "score", "strand", "frame", "attribute")

#################################################### 
# retrieving tss data
#################################################### 
df.ensembl.gtf.for.tss <- df.ensembl.gtf %>% 
  filter(feature == "start_codon") %>% 
  filter(chr %in% c(as.character(1:20), "X", "Y")) %>% 
  mutate(chr = str_c("chr", chr))

df.ensembl.gtf.for.tss %>% dim() # 42925
df.ensembl.gtf.for.tss %>% head(3)
df.ensembl.gtf.for.tss %>% dplyr::select(attribute) %>% head(3)

df.ensembl.gtf.for.tss %>% filter(str_detect(attribute, 'ENSRNOG00000042691'))

df.ensembl.gtf.for.tss # 42,925

# checking attribute keys
tss.attribute.keys <- get_attribute_keys(df.ensembl.gtf.for.tss$attribute) # utils_functions.R
tss.attribute.keys              
#  [1] "exon_number"                  "gene_biotype"                
#  [3] "gene_id"                      "gene_name"                   
#  [5] "gene_source"                  "gene_version"                
#  [7] "projection_parent_transcript" "tag"                         
#  [9] "transcript_biotype"           "transcript_id"               
# [11] "transcript_name"              "transcript_source"           
# [13] "transcript_version" 

# adding columns from attribute
df.ensembl.gtf.for.tss.attribute <- df.ensembl.gtf.for.tss %>% 
  bind_cols(df.ensembl.gtf.for.tss$attribute %>% map_dfr(~extracting_attributes(.x, keys = tss.attribute.keys))) %>% # utils_functions.R
  dplyr::select(-c(source, feature, attribute, score, frame, exon_number)) 

df.ensembl.gtf.for.tss.attribute

# filtering
df.ensembl.gtf.for.tss.attribute.filtered.tag.gnbt <- df.ensembl.gtf.for.tss.attribute %>% 
  filter(tag == "Ensembl_canonical") %>% # 21766
  filter(gene_biotype == "protein_coding") %>% # 21760
  mutate(
    transcript_version = as.numeric(transcript_version)
  )
df.ensembl.gtf.for.tss.attribute.filtered.tag.gnbt %>% dim() # 21760
df.ensembl.gtf.for.tss.attribute.filtered.tag.gnbt %>% head() # 21760

# step 1. non-redundant gene_id dataset
df.ensembl.gtf.for.tss.attribute.filtered.tag.gnbt.geneid.NOdup <- df.ensembl.gtf.for.tss.attribute.filtered.tag.gnbt %>% 
  add_count(gene_id) %>%
  filter(n == 1) # 21690

# step 2. duplicated dataset per gene_id & filtering
df.ensembl.gtf.for.tss.attribute.filtered.tag.gnbt.geneid.dup.1pick <- df.ensembl.gtf.for.tss.attribute.filtered.tag.gnbt %>% 
  add_count(gene_id) %>%
  filter(n > 1) %>%
  group_by(gene_id) %>%
  arrange(desc(transcript_version)) %>%
  slice_head(n = 1) %>%
  ungroup() # 35

# step 3. combining step 1 & step 2 using bind_rows
df.ensembl.gtf.for.tss.DISTINCT.geneid <- bind_rows(df.ensembl.gtf.for.tss.attribute.filtered.tag.gnbt.geneid.NOdup,
                                                                                   df.ensembl.gtf.for.tss.attribute.filtered.tag.gnbt.geneid.dup.1pick)  %>% # 21776
                                          dplyr::select(-n) %>% 
                                          mutate(tss.id = paste(chr, start, end, strand, gene_id, gene_name, sep = ":")) %>% 
                                          dplyr::select(chr, start, end, strand, gene_id, gene_name, tss.id) 

df.ensembl.gtf.for.tss.DISTINCT.geneid %>% # 21725
  # count(chr, start, end, gene_id) #%>% # 21,725
  count(chr, start, end, strand, gene_id) #%>% # 21,725

df.ensembl.gtf.for.tss.DISTINCT.geneid %>% count(gene_id) # 21,725

# tss data integrity check
df.ensembl.gtf.for.tss.DISTINCT.geneid # %>% 
  # count(chr) %>% print(n = Inf)
  # head(3) # chr     start       end strand gene_id        gene_name
  # dim() # 21,725
  # distinct(gene_id) # 21,725
  # filter(str_detect(gene_id, 'LOC|RGD')) # 0

#################################################### 
# retrieving exon data
#################################################### 
df.ensembl.gtf.for.exon.raw <- df.ensembl.gtf %>% 
  filter(feature == "exon") %>% 
  filter(chr %in% c(as.character(1:20), "X", "Y")) %>% 
  mutate(chr = str_c("chr", chr))

# checking attribute keys
exon.attribute.keys <- get_attribute_keys(df.ensembl.gtf.for.exon.raw$attribute) # utils_functions.R
exon.attribute.keys      

#  [1] "exon_id"                      "exon_number"                 
#  [3] "exon_version"                 "gene_biotype"                
#  [5] "gene_id"                      "gene_name"                   
#  [7] "gene_source"                  "gene_version"                
#  [9] "projection_parent_transcript" "tag"                         
# [11] "transcript_biotype"           "transcript_id"               
# [13] "transcript_name"              "transcript_source"           
# [15] "transcript_version"  

# adding columns from attribute
df.ensembl.gtf.for.exon.attribute <- df.ensembl.gtf.for.exon.raw %>% 
  bind_cols(df.ensembl.gtf.for.exon.raw$attribute %>% map_dfr(~extracting_attributes(.x, keys = exon.attribute.keys))) %>% # utils_functions.R
  dplyr::select(-c(source, feature, attribute, score, frame)) 

df.ensembl.gtf.for.exon.attribute # 526,204

# filtering
df.ensembl.gtf.for.exon <- df.ensembl.gtf.for.exon.attribute %>% # 526,204
  filter(tag == "Ensembl_canonical") %>%          # 1st filter
  filter(gene_biotype == "protein_coding") %>%     # 2nd filter
  mutate(exon_number = as.numeric(exon_number)) %>% # convert to numeric
  group_by(gene_id) %>% 
  slice_max(order_by = exon_number, n = 1, with_ties = FALSE) %>%      # keep only the largest exon_number
  ungroup() %>% 
  mutate(ensembl_exon_id = str_c(chr, ':', start, ':', end, ':', gene_name, ':', gene_id, ':', exon_number)) %>% 
  dplyr::select(-c(gene_biotype, tag, transcript_biotype, transcript_version))

df.ensembl.gtf.for.exon %>% dim() # 23,024
df.ensembl.gtf.for.exon %>% head(3)
df.ensembl.gtf.for.exon %>% add_count(gene_id) %>% filter(n == 1) # 23,024
df.ensembl.gtf.for.exon

df.ensembl.gtf.for.exon %>% 
  filter(gene_id == "ENSG00000157764")

df.tss.ensembl <- df.ensembl.gtf.for.tss.DISTINCT.geneid %>% 
  left_join(df.ensembl.gtf.for.exon %>% dplyr::select(gene_id, ensembl_exon_id), by = "gene_id") %>%  # 21,725
  left_join(df.refgene.gtf.for.exon %>% dplyr::select(gene_name, refseq_exon_id), by = "gene_name")

# integrity check
df.tss.ensembl %>% filter(!is.na(ensembl_exon_id) & !is.na(refseq_exon_id)) # 14,071/21,725
df.tss.ensembl %>% filter(is.na(ensembl_exon_id) & is.na(refseq_exon_id)) # 1,360/21,725
df.tss.ensembl %>% filter(is.na(ensembl_exon_id) & !is.na(refseq_exon_id)) # 0/21,725
df.tss.ensembl %>% filter(!is.na(ensembl_exon_id) & is.na(refseq_exon_id)) # 6,294/21,725

# GRanges Obj.: df.tss.ensembl.GR: chr   start     end gene_id                      tss.id
df.tss.ensembl.GR <- GRanges(
  seqnames = as.character(df.tss.ensembl$chr),
  ranges = IRanges(
    start = as.numeric(df.tss.ensembl$start),
    end = as.numeric(df.tss.ensembl$end)
  ),
  strand = df.tss.ensembl$strand
)

# meta data
mcols(df.tss.ensembl.GR)$tss_id <- df.tss.ensembl$tss.id
mcols(df.tss.ensembl.GR)$gene_id <- df.tss.ensembl$gene_id
mcols(df.tss.ensembl.GR)$gene_name <- df.tss.ensembl$gene_name
mcols(df.tss.ensembl.GR)$ensembl_exon_id <- df.tss.ensembl$ensembl_exon_id
mcols(df.tss.ensembl.GR)$refseq_exon_id <- df.tss.ensembl$refseq_exon_id

df.tss.ensembl.GR # 21,725

########################
# 3. TSS
# 3-3. overall distribution of TSS on loops
########################
OVERALL.df.DISTINCT.loop.deep.sample.all.GR # 30928

index.distinct.tss.w.OVERALL.whole.loop <- findOverlaps(
  df.tss.ensembl.GR, 
  OVERALL.df.DISTINCT.loop.deep.sample.all.GR, 
  type = "any",
  select = "all"
)

index.distinct.tss.w.OVERALL.whole.loop # any: 359130/ w/o capping lt2mb: 208296/ ENSEMBL: 259900

OVERALL.loop.for.tss.hits <- subjectHits(index.distinct.tss.w.OVERALL.whole.loop)
OVERALL.tss.on.loop.hits <- queryHits(index.distinct.tss.w.OVERALL.whole.loop)

df.tss.dist.result <- tibble(
  loop.id = OVERALL.df.DISTINCT.loop.deep.sample.all$loop.id[OVERALL.loop.for.tss.hits],
  loop.start = OVERALL.df.DISTINCT.loop.deep.sample.all$x0[OVERALL.loop.for.tss.hits],
  loop.end = OVERALL.df.DISTINCT.loop.deep.sample.all$y3[OVERALL.loop.for.tss.hits],
  loop.res=OVERALL.df.DISTINCT.loop.deep.sample.all$resolution[OVERALL.loop.for.tss.hits],
  tss_chr = df.tss.ensembl$chr[OVERALL.tss.on.loop.hits],
  tss_start = df.tss.ensembl$start[OVERALL.tss.on.loop.hits],
  tss_end = df.tss.ensembl$end[OVERALL.tss.on.loop.hits],
  tss_id = df.tss.ensembl$tss.id[OVERALL.tss.on.loop.hits],
  tss_geneid = df.tss.ensembl$gene_id[OVERALL.tss.on.loop.hits],
  tss_strand = df.tss.ensembl$strand[OVERALL.tss.on.loop.hits]
) %>% 
  mutate(across(where(is.numeric), ~ format(., scientific = FALSE)))

df.tss.dist.result %>% dim() # any: 359130  9/ w/o capping lt2mb: 208296      9/ ENSEMBL: 259900  9
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
  dplyr::select(loop.id, value, loop.res)

relative.pos.df.tss.dist.result %>% head()

########################
# 3. TSS
# 3-3. overall distribution of TSS on loops: figures
# 3-3-1. by CHROMOSOME
########################

checking_component_distribution(relative.pos.df.tss.dist.result, "tss") # utils_functions.R

relative.pos.df.tss.dist.result %>% head()

########################
# 3. TSS
# 3-3. overall distribution of TSS on loops: figures
# 3-3-2. by resolution
########################

plot.tss.hist <- plot_histogram(relative.pos.df.tss.dist.result) # utils_functions.R
plot.tss.dens <- plot_density(relative.pos.df.tss.dist.result) # utils_functions.R

saving_combined_plot(plot.tss.hist, plot.tss.dens, # utils_functions.R
                   "figures/submission/lt2mb/overall_distribution_of_TSS_by_resolution_wo_capping_lt2mb_ENSEMBL.pdf")
                   
########################
# 3. TSS
# 3-4. Distribution of TSS at each end in a loop: for the number of TSS used in filtering valid loops: figures
# so, the object should be used one WITH padding on df.DISTINCT.loop.deep.sample.all: 1/2 distance for OUTER & 1/4 distance for INNER
# 3-4-1. adding padding at each end, x12 & y12
########################
df.DISTINCT.loop.deep.sample.all # 31773
df.DISTINCT.loop.deep.sample.all.lt.2mb  %>% dim()# 31019
df.DISTINCT.loop.deep.sample.all.lt.2mb %>% head(2) # 31019

########################
# 4. promoter
# 4-1. Exploratory Data analysis (EDA) : promoter
# 4-2. promoter data preprocessing: id and dedup GRange Obj.
########################
############################
# promoter resource 1: deprecated
############################
# https://epd.expasy.org/epd/EPDnew_select.php using GUI
epd.rat.promoter.rn6.bed <- import("~/dropbox/Gateway_to_Hao/enhancer/data/epdnew/R_norvegicus_epdnew_rn6.bed", format = "BED")
epd.rat.promoter.rn6.bed # 12022

############################
# promoter resource 2
############################
# Download: https://epd.expasy.org/ftp/epdnew/R_norvegicus/
Rn_EPDnew_001_rn6.bed.raw <- read.table("~/dropbox/Gateway_to_Hao/enhancer/data/epdnew/001/Rn_EPDnew_001_rn6.bed", sep = " ") %>% dplyr::rename(
  seqnames = V1,
  start = V7,
  name = V4,
  score = V5,
  strand = V6) %>% dplyr::select(-starts_with("V")) %>% 
  mutate(start = as.numeric(start), end = start + 1, score = 1)

Rn_EPDnew_001_rn6.bed.raw # 12601
Rn_EPDnew_001_rn6.bed.gr <- GRanges(
  seqnames = Rn_EPDnew_001_rn6.bed.raw$seqnames,
  ranges = IRanges(
    start = Rn_EPDnew_001_rn6.bed.raw$start,
    end = Rn_EPDnew_001_rn6.bed.raw$end
  ),
  strand = Rn_EPDnew_001_rn6.bed.raw$strand,
  name = Rn_EPDnew_001_rn6.bed.raw$name,
  score = Rn_EPDnew_001_rn6.bed.raw$score
)
Rn_EPDnew_001_rn6.bed.gr
BiocIO::export(Rn_EPDnew_001_rn6.bed.gr, "/Users/pete/dropbox/Gateway_to_Hao/enhancer/data/epdnew/001/Rn_EPDnew_001_rn6.bed.gr", format = "BED")
Rn_EPDnew_001_rn6.bed <- import("/Users/pete/dropbox/Gateway_to_Hao/enhancer/data/epdnew/001/Rn_EPDnew_001_rn6.bed.gr", format = "BED")
Rn_EPDnew_001_rn6.bed

# chain file
chain.rn6.to.rn7 <- import.chain("~/dropbox/Gateway_to_Hao/enhancer/data/epdnew/rn6ToRn7.over.chain")

# 3. liftOver: rn6 → rn7
Rn_EPDnew_001_rn7.list <- liftOver(Rn_EPDnew_001_rn6.bed, chain.rn6.to.rn7)
Rn_EPDnew_001_rn7.list # 12601

# 4. GRangesList → GRanges
Rn_EPDnew_001_rn7 <- unlist(Rn_EPDnew_001_rn7.list)

# 5. BED export
export(Rn_EPDnew_001_rn7, "~/dropbox/Gateway_to_Hao/enhancer/data/epdnew/001/Rn_EPDnew_001_rn7.bed", format = "BED")
# export(Rn_EPDnew_001_rn7, "~/dropbox/Gateway_to_Hao/enhancer/data/epdnew/R_norvegicus_epdnew_rn7.bed", format = "BED")

df.Rn_EPDnew_001_rn7 <- as_tibble(Rn_EPDnew_001_rn7)
df.Rn_EPDnew_001_rn7
df.Rn_EPDnew_001_rn7 %>% add_count(seqnames, start, end, name) %>% filter(n > 1) # 0
df.Rn_EPDnew_001_rn7 %>% dim() # 12533/12601
df.Rn_EPDnew_001_rn7 %>% count(name) %>% filter(n > 1) # 0

# ENSEMBL ID with gene symbol
df.gene.mapping.for.promoter <- read_tsv("~/dropbox/Gateway_to_Hao/enhancer/data/epdnew/001/db/promoter_ensembl.txt", 
                        col_names = c("promoter_id", "gene_id"),
                        col_types = cols(
                          promoter_id = col_character(),
                          gene_id = col_character()
                        )) %>% # 12,793
                        mutate(gene_id = if_else(str_detect(promoter_id, 'Cfb_1'), "ENSRNOG00000051158.3", gene_id)) %>% # https://useast.ensembl.org/Rattus_norvegicus/Gene/Idhistory?g=ENSRNOG00000051158
                        distinct() # 12,600

df.gene.mapping.for.promoter #12,600
df.gene.mapping.for.promoter %>% add_count(promoter_id) %>% filter(n > 1) # 0
df.gene.mapping.for.promoter %>% add_count(gene_id) %>% filter(n > 1) # 1072
# example 
#    promoter_id gene_id                n
#  1 Sgk1_1      ENSRNOG00000011815     3
#  2 Sgk1_2      ENSRNOG00000011815     3
#  3 Sgk1_3      ENSRNOG00000011815     3

# adding gene_id (ENSEMBL)
df.promoter.rn7.raw <- left_join(df.Rn_EPDnew_001_rn7, df.gene.mapping.for.promoter, by = c("name" = "promoter_id"))
df.promoter.rn7.raw %>% dim() # 11,953| 12533
df.promoter.rn7.raw %>% head()
df.promoter.rn7.raw %>% filter(is.na(gene_id)) # 0| 1 : 1 chr3     115016225 115016226     2 -      AABR07053687_2     1 NA     

df.gene.mapping.for.promoter %>% filter(str_detect(promoter_id, "AABR07053687"))

df.ensembl.gtf.for.exon

# adding exon information from ENSEMBL and RefSeq
df.promoter.rn7.exon_id <- df.promoter.rn7.raw %>% # seqnames   start     end width strand name     score gene_id
  mutate(gene_id = if_else(str_detect(name, "AABR07053687"), "ENSRNOG00000015756", gene_id)) %>% 
  dplyr::rename(promoter_id = name) %>% 
  mutate(gene_name = str_split_n(promoter_id, '_', 1)) %>%  # utils_functions.R
  left_join(df.ensembl.gtf.for.exon %>% dplyr::select(gene_id, ensembl_exon_id), by = "gene_id") %>% 
  left_join(df.refgene.gtf.for.exon %>% dplyr::select(gene_name, refseq_exon_id), by = "gene_name")

# integrity check
df.promoter.rn7.exon_id %>% filter(!is.na(ensembl_exon_id) & !is.na(refseq_exon_id)) # both 10,068
df.promoter.rn7.exon_id %>% filter(is.na(ensembl_exon_id) & is.na(refseq_exon_id)) # none 445
df.promoter.rn7.exon_id %>% filter(is.na(ensembl_exon_id) & !is.na(refseq_exon_id)) # RefSeq only 533
df.promoter.rn7.exon_id %>% filter(!is.na(ensembl_exon_id) & is.na(refseq_exon_id)) # ENSEMBL only 1487


df.promoter.rn7.exon_id %>% dim() # 12533
df.promoter.rn7.exon_id %>% head()

df.promoter.rn7.exon_id %>% add_count(promoter_id) %>% filter(n > 1) # 0
df.promoter.rn7.exon_id %>% add_count(gene_id) %>% filter(n > 1) # 1,060
df.promoter.rn7.exon_id %>% add_count(seqnames, start, end, gene_id) %>% filter(n > 1) # 8

# Deduplication: keep only _1 promoters when seqnames, start, end, gene_id are identical
df.promoter.rn7 <- df.promoter.rn7.exon_id %>% 
  add_count(seqnames, start, end, gene_id, name = "dup_count") %>% 
  filter(dup_count == 1 | str_detect(promoter_id, "_1$")) %>% 
  dplyr::select(-dup_count)

df.promoter.rn7 # seqnames   start     end width strand promoter_id  score gene_id ensembl_exon_id refseq_exon_id gene_name
df.promoter.rn7 %>% dim() # 12,529 (removed 4 duplicates)
df.promoter.rn7 %>% add_count(seqnames, start, end, gene_id) %>% filter(n > 1) # 0

# df.promoter.rn7 (dups)
df.promoter.rn7.epd <- df.promoter.rn7 %>% 
  dplyr::rename(tss_start = start, tss_end = end) %>% 
  mutate(
      start = tss_start - 40,
      end   = tss_start + 40
    ) %>% 
  mutate(promoter.id = str_c(seqnames,':', start,':', end, ':', gene_id, ':', gene_name, ':', seqnames, ':', tss_start, ':', tss_end))
df.promoter.rn7.epd

df.promoter.rn7.epd.GR <- GRanges(
  seqnames = df.promoter.rn7.epd$seqnames,
  ranges = IRanges(
    start = df.promoter.rn7.epd$start, 
    end = df.promoter.rn7.epd$end)
)
# metadata
mcols(df.promoter.rn7.epd.GR) <- df.promoter.rn7.epd[, c("promoter.id", "gene_id", "gene_name", "ensembl_exon_id", "refseq_exon_id")]

########################
# 4. promoter
# 4-3. overall distribution of promoter on loops
########################
index.promoter.w.OVERALL.whole.loop <- findOverlaps(
  df.promoter.rn7.epd.GR, 
  OVERALL.df.DISTINCT.loop.deep.sample.all.GR, 
  type = "any",
  select = "all"
)
index.promoter.w.OVERALL.whole.loop # any: 264932 // w/o capping lt2mb: 158517// ENSMBL 159944

OVERALL.loop.for.promoter.hits <- subjectHits(index.promoter.w.OVERALL.whole.loop)
OVERALL.promoter.on.loop.hits <- queryHits(index.promoter.w.OVERALL.whole.loop)

df.promoter.dist.result <- data.frame(
  loop.id = OVERALL.df.DISTINCT.loop.deep.sample.all$loop.id[OVERALL.loop.for.promoter.hits],
  loop.start = OVERALL.df.DISTINCT.loop.deep.sample.all$x0[OVERALL.loop.for.promoter.hits],
  loop.end = OVERALL.df.DISTINCT.loop.deep.sample.all$y3[OVERALL.loop.for.promoter.hits],
  loop.res = OVERALL.df.DISTINCT.loop.deep.sample.all$resolution[OVERALL.loop.for.promoter.hits],
  promoter_id = df.promoter.rn7.epd$promoter.id[OVERALL.promoter.on.loop.hits],
  promoter_start = df.promoter.rn7.epd$start[OVERALL.promoter.on.loop.hits],
  promoter_end = df.promoter.rn7.epd$end[OVERALL.promoter.on.loop.hits],
  promoter_gene_id = df.promoter.rn7.epd$gene_id[OVERALL.promoter.on.loop.hits]
) %>% 
  mutate(across(where(is.numeric), ~ format(., scientific = FALSE)))

df.promoter.dist.result %>% head()
df.promoter.dist.result %>% dim() # any: 264932(dedups)// w/o capping lt2mb: 158517// ENSMBL: 159944

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
  dplyr::select(loop.id, promoter_id, value, loop.res, promoter_gene_id)

relative.pos.df.promoter.dist.result %>% head()
relative.pos.df.promoter.dist.result %>% dim() # any: 264932(dedups) // w/o capping lt2mb: 158517 // ENSMBL: 159944
########################
# 4. promoter
# 4-3. overall distribution of promoter on loops: figures
# 4-3-1. by CHROMOSOME
########################

checking_component_distribution(relative.pos.df.promoter.dist.result, "promoter") # utils_functions.R

relative.pos.df.promoter.dist.result %>% head()

########################
# 4. promoter
# 4-3. overall distribution of promoter on loops: figures
# 4-3-2. by resolution
########################

plot.promoter.hist <- plot_histogram(relative.pos.df.promoter.dist.result) # utils_functions.R
plot.promoter.dens <- plot_density(relative.pos.df.promoter.dist.result) # utils_functions.R

saving_combined_plot(plot.promoter.hist, plot.promoter.dens, # utils_functions.R
                   "figures/submission/lt2mb/overall_distribution_of_promoter_by_resolution_wo_capping_lt2mb_ENSEMBL.pdf")

#########################################
#########################################
# Figure for new ones
#########################################
#########################################

combining_and_save_plots(plot.ctcf.hist, plot.tss.hist, plot.promoter.hist, "histogram_combined_all_ENSEMBL.png") # utils_functions.R
combining_and_save_plots(plot.ctcf.dens, plot.tss.dens, plot.promoter.dens, "density_combined_all_ENSEMBL.png") # utils_functions.R

##########################################################
##########################################################
# checking hits of TSS & promoter in an anchor using distanceToNearest
##########################################################
##########################################################

##########################################################
# 5. Combining TSS and Promoter datasets with exon information
##########################################################

# Prepare TSS data
df.tss.ensembl %>% head(3) # chr start end strand gene_id gene_name tss.id ensembl_exon_id refseq_exon_id

# Prepare Promoter data (already has exon information)
df.promoter.rn7.epd %>% head(3) # seqnames start end width strand promoter_id score gene_id gene_name ensembl_exon_id refseq_exon_id

# Combine TSS and Promoter datasets
df.gene_tss_and_pro <- bind_rows(
  df.tss.ensembl %>% 
    dplyr::select(chr, start, end, gene_id, gene_name, tss.id, ensembl_exon_id, refseq_exon_id) %>%
    mutate(across(c(start, end), as.numeric)) %>%
    mutate(component = "tss") %>%
    dplyr::rename(component_id = tss.id) %>%
    mutate(component_id = str_c(component_id, component, sep='|')),
  df.promoter.rn7.epd %>% 
    dplyr::rename(chr = seqnames) %>%
    mutate(promoter.id = str_c(chr, ':', start, ':', end, ':', gene_id, ':', gene_name)) %>%
    dplyr::select(chr, start, end, gene_id, gene_name, promoter.id, ensembl_exon_id, refseq_exon_id) %>%
    mutate(across(c(start, end), as.numeric)) %>%
    mutate(component = "pro") %>%
    dplyr::rename(component_id = promoter.id) %>%
    mutate(component_id = str_c(component_id, component, sep='|'))
)

df.gene_tss_and_pro %>% head(3)
df.gene_tss_and_pro %>% dim() # TSS: 21,725 + Promoter: 12,529 = 34,254

# Create GRanges object with exon information
df.gene_tss_and_pro.GR <- GRanges(
  seqnames = df.gene_tss_and_pro$chr,
  ranges = IRanges(start = df.gene_tss_and_pro$start, end = df.gene_tss_and_pro$end),
  gene_id = df.gene_tss_and_pro$gene_id,
  gene_name = df.gene_tss_and_pro$gene_name,
  component_id = df.gene_tss_and_pro$component_id,
  component = df.gene_tss_and_pro$component,
  ensembl_exon_id = df.gene_tss_and_pro$ensembl_exon_id,
  refseq_exon_id = df.gene_tss_and_pro$refseq_exon_id
)

df.gene_tss_and_pro.GR # 34,254

df.DISTINCT.loop.deep.sample.all.lt.2mb %>% dim() # 31019
df.DISTINCT.loop.deep.sample.all.lt.2mb %>% head(3)

df.DISTINCT.loop.deep.sample.all.lt.2mb.prep <- df.DISTINCT.loop.deep.sample.all.lt.2mb %>% 
  mutate(mid_x = floor((x0 + x3) / 2), mid_y = floor((y0 + y3) / 2))

# Create GRanges for UP and DOWN anchors using creating_granges function
df.DISTINCT.loop.deep.sample.all.lt.2mb.prep.UP.GR <- creating_granges(df.DISTINCT.loop.deep.sample.all.lt.2mb.prep, direction = "up", use_anchor = TRUE) # utils_functions.R
df.DISTINCT.loop.deep.sample.all.lt.2mb.prep.DOWN.GR <- creating_granges(df.DISTINCT.loop.deep.sample.all.lt.2mb.prep, direction = "down", use_anchor = TRUE) # utils_functions.R

hits_up <- distanceToNearest(df.DISTINCT.loop.deep.sample.all.lt.2mb.prep.UP.GR, df.gene_tss_and_pro.GR)
hits_up
df_up <- data.frame(
  queryHits = queryHits(hits_up),
  subjectHits = subjectHits(hits_up),
  distance = mcols(hits_up)$distance
) %>%
  mutate(
    loop.id = mcols(df.DISTINCT.loop.deep.sample.all.lt.2mb.prep.UP.GR)[queryHits, "loop.id"],
    resolution = mcols(df.DISTINCT.loop.deep.sample.all.lt.2mb.prep.UP.GR)[queryHits, "resolution"],
    gene_id = mcols(df.gene_tss_and_pro.GR)[subjectHits, "gene_id"],
    gene_name = mcols(df.gene_tss_and_pro.GR)[subjectHits, "gene_name"],
    component = mcols(df.gene_tss_and_pro.GR)[subjectHits, "component"],
    component_id = mcols(df.gene_tss_and_pro.GR)[subjectHits, "component_id"],
    ensembl_exon_id = mcols(df.gene_tss_and_pro.GR)[subjectHits, "ensembl_exon_id"],
    refseq_exon_id = mcols(df.gene_tss_and_pro.GR)[subjectHits, "refseq_exon_id"],
    WHERE = "UP"
  )

hits_down <- distanceToNearest(df.DISTINCT.loop.deep.sample.all.lt.2mb.prep.DOWN.GR, df.gene_tss_and_pro.GR)
hits_down
df_down <- data.frame(
  queryHits = queryHits(hits_down),
  subjectHits = subjectHits(hits_down),
  distance = mcols(hits_down)$distance
) %>%
  mutate(
    loop.id = mcols(df.DISTINCT.loop.deep.sample.all.lt.2mb.prep.DOWN.GR)[queryHits, "loop.id"],
    resolution = mcols(df.DISTINCT.loop.deep.sample.all.lt.2mb.prep.DOWN.GR)[queryHits, "resolution"],
    gene_id = mcols(df.gene_tss_and_pro.GR)[subjectHits, "gene_id"],
    gene_name = mcols(df.gene_tss_and_pro.GR)[subjectHits, "gene_name"],
    component = mcols(df.gene_tss_and_pro.GR)[subjectHits, "component"],
    component_id = mcols(df.gene_tss_and_pro.GR)[subjectHits, "component_id"],
    ensembl_exon_id = mcols(df.gene_tss_and_pro.GR)[subjectHits, "ensembl_exon_id"],
    refseq_exon_id = mcols(df.gene_tss_and_pro.GR)[subjectHits, "refseq_exon_id"],
    WHERE = "DOWN"
  )

df_up %>% head(3)
df_down %>% head(3)

df.final.up.down.tss.pro.nearest <- bind_rows(df_up, df_down)
df.final.up.down.tss.pro.nearest %>% dim() # 62038
df.final.up.down.tss.pro.nearest %>% head(3)

df.final.up.down.tss.pro.nearest %>% count(loop.id) %>% count(n) # 2 31019 : all UP and DOWN = 2 rows per loop.id

# quantile boxplot
approach2.stats <- df.final.up.down.tss.pro.nearest %>%
  summarise(
    Q1 = quantile(distance, 0.25),
    Median = median(distance),
    Q3 = quantile(distance, 0.75)
  )
approach2.stats$Q1
approach2.stats$Median
approach2.stats$Q3

mean(df.final.up.down.tss.pro.nearest$distance) # 55118

ggplot(df.final.up.down.tss.pro.nearest, aes(y = distance)) +
  geom_boxplot(fill = "#A6CEE3", color = "#1F78B4", outlier.color = "red", outlier.shape = 16) +
  # scale_y_continuous() +
  scale_y_log10() +
  geom_hline(yintercept = approach2.stats$Q1, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = approach2.stats$Median, linetype = "dashed", color = "darkgreen") +
  geom_hline(yintercept = approach2.stats$Q3, linetype = "dashed", color = "purple") +
  annotate("text", x = 1.2, y = approach2.stats$Q1, label = paste("Q1:", round(approach2.stats$Q1)), color = "blue", vjust = -0.5) +
  annotate("text", x = 1.2, y = approach2.stats$Median, label = paste("Median:", round(approach2.stats$Median)), color = "darkgreen", vjust = -0.5) +
  annotate("text", x = 1.2, y = approach2.stats$Q3, label = paste("Q3:", round(approach2.stats$Q3)), color = "purple", vjust = -0.5) +
  labs(
    title = "Distance between TSS or Promoter and Anchor Distribution (Log Scale)",
    y = "Distance (bp, log10)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

df.final.up.down.tss.pro.nearest %>% dim() # 62038
df.final.up.down.tss.pro.nearest %>% head()

# getting cases per loop with minimum distance including 0
df.final.up.down.tss.pro.nearest.filtered <- df.final.up.down.tss.pro.nearest %>%
  group_by(loop.id) %>%
  filter(if (all(distance == 0)) {
    TRUE   
  } else {
    distance == min(distance)  
  }) %>%
  ungroup()

df.final.up.down.tss.pro.nearest.filtered %>% dim() # 33399/62038 = 31019 * 2 (2380 are 0 in both/28639 are only one in either one of 2 anchors)
df.final.up.down.tss.pro.nearest.filtered %>% head() 
df.final.up.down.tss.pro.nearest.filtered %>% count(loop.id) %>% 
  count(n)

###############################################################################################
# columns:
# loop_validity: the number of overlapping: 1 overlapping or 2 overlapping: clear or vague
# functionality: target gene location: functional, unfunctional, unclear
###############################################################################################

# Apply both functions to the filtered dataset
df.final.up.down.tss.pro.nearest.filtered.validity.functionality <- df.final.up.down.tss.pro.nearest.filtered %>%
  adding_loop_validity() %>% # utils_functions.R
  adding_functionality() # utils_functions.R

df.final.up.down.tss.pro.nearest.filtered.validity.functionality %>% head()
df.final.up.down.tss.pro.nearest.filtered.validity.functionality %>% count(loop_validity)
df.final.up.down.tss.pro.nearest.filtered.validity.functionality %>% count(functionality)

df.final.up.down.tss.pro.nearest.filtered.validity.functionality %>% count(loop_validity, functionality)

# Apply the function
df.final.up.down.tss.pro.nearest.filtered.validity.functionality.final <- df.final.up.down.tss.pro.nearest.filtered.validity.functionality %>%
  adding_final_decision() # utils_functions.R

df.final.up.down.tss.pro.nearest.filtered.validity.functionality.final

# Check the results
df.final.up.down.tss.pro.nearest.filtered.validity.functionality.final %>% count(loop_validity, final)
df.final.up.down.tss.pro.nearest.filtered.validity.functionality.final %>% count(final, functionality)
df.final.up.down.tss.pro.nearest.filtered.validity.functionality.final %>% 
  filter(loop_validity == "vague") %>% 
  count(final, functionality)

df.final.up.down.tss.pro.nearest.filtered.only.one.anchor <- df.final.up.down.tss.pro.nearest.filtered.validity.functionality.final %>% 
  filter(loop_validity == "clear") # Only keep clear cases
df.final.up.down.tss.pro.nearest.filtered.only.one.anchor

# Extract functional loops from clear and vague cases
# 1. Clear + functional + decided
df.functional.clear <- df.final.up.down.tss.pro.nearest.filtered.validity.functionality.final %>%
  filter(loop_validity == "clear" & functionality == "functional" & final == "decided")
df.functional.clear # 11,636

# 2. Vague + functional + decidable
df.functional.vague.decidable <- df.final.up.down.tss.pro.nearest.filtered.validity.functionality.final %>%
  filter(loop_validity == "vague" & functionality == "functional" & final == "decidable")
df.functional.vague.decidable # 1,086

# Combined: all functional loops (clear decided + vague decidable)
df.functional.loops <- bind_rows(
  df.functional.clear,
  df.functional.vague.decidable
)

df.functional.clear %>% dim()
df.functional.vague.decidable %>% dim()
df.functional.loops %>% dim()

# overlapping in both anchors
df.final.up.down.tss.pro.nearest.only.one.anchor.both <- df.final.up.down.tss.pro.nearest.only.one.anchor %>% 
  add_count(loop.id) %>% 
  filter(n == 2) %>% 
  arrange(loop.id)
df.final.up.down.tss.pro.nearest.only.one.anchor.both

missing_229 <- df.final.up.down.tss.pro.nearest.only.one.anchor.both  %>% # 4760
  # count(component) # pro        2909 tss        1851
  filter(component == 'pro') %>% 
  left_join(df.ensembl.gtf.for.exon, by = ('gene_name')) %>%  # 401
  filter(is.na(exon_number)) %>% 
  dplyr::select(gene_id.x, gene_name) %>% # 401
  left_join(df.refgene.gtf.parsed.exon.number, by = ('gene_name')) %>% 
  filter(is.na(exon_number)) %>% # 229
  dplyr::rename(gene_id = gene_id.x)

missing_229
missing_229 %>% distinct(gene_id.x)


df.ensembl.gtf.for.exon

left_join(df.ensembl.gtf.for.exon, by = ('gene_id')) %>% 
  filter(is.na(exon_number))
df.refgene.gtf.parsed.exon.number

df.ensembl.gtf.for.exon
df.final.up.down.tss.pro.nearest.only.one.anchor %>% 
# df.final.up.down.tss.pro.nearest %>% 
  ggplot(aes(x = log2(distance + 1))) +
  geom_histogram(binwidth = 1) +   
  labs(
    title = "Histogram of Distances",
    x = "Distance (bp)",
    y = "Count"
  ) +
  theme_minimal() +
  facet_wrap(~resolution)

df.final.up.down.tss.pro.nearest %>% dim()

df.final.up.down.tss.pro.nearest.gene.location <- df.final.up.down.tss.pro.nearest %>% 
  left_join(df.ensembl.gtf.for.exon, by = "gene_id") %>% 
  mutate(loop_parts = str_split(loop.id, "_")) %>%
  mutate(
    chr1 = map_chr(loop_parts, 1),
    x1   = as.numeric(map_chr(loop_parts, 2)),
    x2   = as.numeric(map_chr(loop_parts, 3)),
    chr2 = map_chr(loop_parts, 4),
    y1   = as.numeric(map_chr(loop_parts, 5)),
    y2   = as.numeric(map_chr(loop_parts, 6)),
    anchor_distance = as.numeric(map_chr(loop_parts, 7))
  ) %>%
  mutate(
    x_mid = (x1 + x2) / 2,
    y_mid = (y1 + y2) / 2
  ) %>%
  mutate(location = start >= x1 & end <= y2)

df.final.up.down.tss.pro.nearest.gene.location.functional <- df.final.up.down.tss.pro.nearest.gene.location %>% 
  group_by(loop.id) %>%
  mutate(
    functional = case_when(
      all(distance == 0) ~ "ambiguous",                     # both 0
      any(distance == 0) ~ "highly",                       # one 0
      TRUE              ~ "unfunctional"                   # neither 0
    )
  ) %>%
  ungroup()

df.final.up.down.tss.pro.nearest.gene.location.functional %>% 
  filter(functional == 'ambiguous') %>% 
  filter(count(location)
  filter(distance < get_distance(15)) %>%  # utils_functions.R
  count(location)
  head()
# 1    FALSE 34873
# 2     TRUE 25478
# 3       NA  1687
  
  %>% 
  filter(is.na(exon_number)) %>% # 1687 from 100% promoter
  count(component)
  dim()
  head()

get_distance(0)   # 0:       gprofiler: https://biit.cs.ut.ee/gplink/l/aLDYW6PA7QD| https://biit.cs.ut.ee/gplink/l/ajrmuiTFCQJ # utils_functions.R
get_distance(15)  # 32767:   gprofiler: https://biit.cs.ut.ee/gplink/l/ae6cRE94rQa| https://biit.cs.ut.ee/gplink/l/abVmsUcaRQ8 # utils_functions.R
get_distance(20)  # 1048575: gprofiler: https://biit.cs.ut.ee/gplink/l/awPCWBYS9RJ| https://biit.cs.ut.ee/gplink/l/alXrPBsdGRv # utils_functions.R

df.final.up.down.tss.pro.nearest.only.one.anchor %>% filter(distance == 0) %>%  count(gene_id) # 6697

df.final.up.down.tss.pro.nearest <- df.final.up.down.tss.pro.nearest.only.one.anchor %>% filter(distance != 0)
df.final.up.down.tss.pro.nearest.only.one.anchor %>% filter(distance == 0) %>% add_count(loop.id) %>% filter(n > 1) %>% 
  count(component)
  count(resolution)

df.final.up.down.tss.pro.nearest

approach_2nd_analyze_loops_by_threshold(df.final.up.down.tss.pro.nearest, threshold_distance = get_distance(0), top_n_genes = 70, print_top_n = 50) # utils_functions.R
approach_2nd_analyze_loops_by_threshold(df.final.up.down.tss.pro.nearest %>% filter(distance != 0), threshold_distance = get_distance(15), top_n_genes = 20, print_top_n = 20) # utils_functions.R
approach_2nd_analyze_loops_by_threshold(df.final.up.down.tss.pro.nearest, threshold_distance = get_distance(20), top_n_genes = 70, print_top_n = 50) # utils_functions.R


mean(df.final.up.down.tss.pro.nearest$distance) # 109537.7// midpoint, ENSEMBL61564.34// span, ENSEMBL 55191.25| 55118
threshold_distance <- approach2.stats$Q3 # 82165.75// midpoint, ENSEMBL 62222.75//  span, ENSEMBL 54824| 54779
threshold_distance <- approach2.stats$Median # 28943.5// midpoint, ENSEMBL 23377.5//  span, ENSEMBL 16545.5| 16549.5
threshold_distance <- approach2.stats$Q1 # 8807// midpoint, ENSEMBL 7422 //  span, ENSEMBL 514.25| 512.25
threshold_distance <- 10000 # 2k:46.36%, 3k:49.2%, 5k:54.33%, 10k: 63.87%
threshold_distance
approach_2nd_analyze_loops_by_threshold(df.final.up.down.tss.pro.nearest, threshold_distance = threshold_distance, top_n_genes = 70, print_top_n = 50) # utils_functions.R
# ENSEMBL point: https://biit.cs.ut.ee/gplink/l/aEqIdfCD3TD
# ENSEMBL span: https://biit.cs.ut.ee/gplink/l/aDqh-aFfbQP
df.final.up.down.tss.pro.nearest %>% head()

df.unique_loops <- df.final.up.down.tss.pro.nearest %>%
  filter(distance <= threshold_distance)


ggplot(df.final.up.down.tss.pro.nearest, aes(x = distance)) +
  geom_density(fill = "lightgreen", alpha = 0.5) +
  labs(title = "Density Plot of Loop Distances", x = "Distance (bp)", y = "Density") +
  theme_bw(base_size = 12) + 
  facet_wrap(~resolution)

df.unique_loops %>% dim() # 23524// point ENSEMBL: 24987// span ENSEMBL: 31019

# final
final.loops.from.tss.step <- df.unique_loops %>%
  filter(component == "tss")
final.loops.from.tss.step %>% dim() # point ENSEMBL:13407//span ENSEMBL: 15748// new resources: 17107

final.loops.from.promoter.step <- df.unique_loops %>%
  filter(component == "pro")
final.loops.from.promoter.step %>% dim() # point ENSEMBL: 10117// span ENSEMBL: 13209// new resources:  13912

df.unique_loops  %>% head()
df.unique_loops %>%
  add_count(loop.id) %>%                 # loop.id별로 몇 번 등장했는지 세기
  filter(n >= 2) %>%                     # 2개 이상 등장한 loop.id만 필터링
  arrange(loop.id) %>%
  view()

  distinct(loop.id, resolution) %>%     # 중복 제거 (loop.id-resolution 조합)
  count(resolution)   

library(ggplot2)

df.unique_loops %>% head(3)
df.unique_loops %>% dim() # 23524
df.unique_loops %>% count(loop.id, component) %>%# view()
  count(n) 
# 1 1 9177
# 2 2 1983

df_single_loop <- df.unique_loops %>% # 9177
  group_by(loop.id) %>%
  filter(n() == 1) %>%
  ungroup()

# WHERE distribution
where_dist <- df_single_loop %>%
  count(WHERE)
where_dist
# 1 DOWN   4485
# 2 UP     4692

# component distribution
component_dist <- df_single_loop %>%
  count(component)
component_dist
# 1 pro        4035
# 2 tss        5142

df_multi <- df.unique_loops %>%
  group_by(loop.id) %>%
  filter(n() > 1) %>%
  summarise(where_combo = paste(sort(unique(WHERE)), collapse = "_")) %>%
  count(where_combo)
df_multi # DOWN_UP      1983

df_multi_component <- df.unique_loops %>%
  group_by(loop.id) %>%
  filter(n() > 1) %>%
  summarise(component_combo = paste(sort(unique(component)), collapse = "_")) %>%
  count(component_combo)
df_multi_component
#   component_combo     n 473 + 954 + 556 = 1983
# 1 pro               473
# 2 pro_tss           954
# 3 tss               556

# TODO: comparing join VS bind_rows results
##########################################################
# approach2 for TSS and promoter
##########################################################
df.DISTINCT.loop.deep.sample.all.lt.2mb # 31773/31019
df.gene_tss_and_pro.GR # seqnames              ranges strand ,           gene_id, component_id   component
# 1. tss_and_pro + UPSTREAM (df.DISTINCT.loop.deep.sample.all.lt.2mb, df.gene_tss_and_pro.GR)
index.distinct.tss_and_pro.w.up.loop <- findOverlaps(df.gene_tss_and_pro.GR, 
                                              df.DISTINCT.loop.deep.sample.all.up.GR, 
                                              type = "any",
                                              select = "all")

end.loop.up.tss_and_pro.hits <- subjectHits(index.distinct.tss_and_pro.w.up.loop)
end.tss_and_pro.up.hits <- queryHits(index.distinct.tss_and_pro.w.up.loop)

df.overlapping.tss_and_pro.w.UPSTREAM.result <- tibble(
  up.loop.id = mcols(df.DISTINCT.loop.deep.sample.all.up.GR)$loop.id[end.loop.up.tss_and_pro.hits],
  end.up.distance = mcols(df.DISTINCT.loop.deep.sample.all.up.GR)$end.distance[end.loop.up.tss_and_pro.hits],
  distance = mcols(df.DISTINCT.loop.deep.sample.all.up.GR)$distance[end.loop.up.tss_and_pro.hits],
  resolution = mcols(df.DISTINCT.loop.deep.sample.all.up.GR)$resolution[end.loop.up.tss_and_pro.hits],
  gene.id = mcols(df.gene_tss_and_pro.GR)$gene_id[end.tss_and_pro.up.hits],
  component.id = mcols(df.gene_tss_and_pro.GR)$component_id[end.tss_and_pro.up.hits],
  component = mcols(df.gene_tss_and_pro.GR)$component[end.tss_and_pro.up.hits],
  WHERE = "UP"
) %>% 
  unite(tss_and_pro.loop.up.id, gene.id, up.loop.id, distance, component.id, component, WHERE, sep = "|", remove = FALSE)

df.overlapping.tss_and_pro.w.UPSTREAM.result %>% dim() # 12947
df.overlapping.tss_and_pro.w.UPSTREAM.result %>% head()

df.overlapping.tss_and_pro.w.UPSTREAM.result %>% count(end.up.distance)
#   end.up.distance     n
# 1            5000  1599
# 2           10000  4197
# 3           25000  7151


# 2. tss_and_pro + DOWNSTREAM (df.DISTINCT.loop.deep.sample.all.lt.2mb, df.gene_tss_and_pro.GR)
index.distinct.tss_and_pro.w.down.loop <- findOverlaps(df.gene_tss_and_pro.GR, 
                                                df.DISTINCT.loop.deep.sample.all.down.GR, 
                                                type = "any",
                                                select = "all")

end.loop.down.tss_and_pro.hits <- subjectHits(index.distinct.tss_and_pro.w.down.loop)
end.tss_and_pro.down.hits <- queryHits(index.distinct.tss_and_pro.w.down.loop)

df.overlapping.tss_and_pro.w.DOWNSTREAM.result <- tibble(
  down.loop.id = mcols(df.DISTINCT.loop.deep.sample.all.down.GR)$loop.id[end.loop.down.tss_and_pro.hits],
  end.down.distance = mcols(df.DISTINCT.loop.deep.sample.all.down.GR)$end.distance[end.loop.down.tss_and_pro.hits],
  distance = mcols(df.DISTINCT.loop.deep.sample.all.down.GR)$distance[end.loop.down.tss_and_pro.hits],
  resolution = mcols(df.DISTINCT.loop.deep.sample.all.down.GR)$resolution[end.loop.down.tss_and_pro.hits],
  gene.id = mcols(df.gene_tss_and_pro.GR)$gene_id[end.tss_and_pro.down.hits],
  component.id = mcols(df.gene_tss_and_pro.GR)$component_id[end.tss_and_pro.down.hits],
  component = mcols(df.gene_tss_and_pro.GR)$component[end.tss_and_pro.down.hits],
  WHERE = "DOWN"
) %>% 
  unite(tss_and_pro.loop.down.id, gene.id, down.loop.id, distance, component.id, component, WHERE, sep = "|", remove = FALSE)

df.overlapping.tss_and_pro.w.DOWNSTREAM.result %>% dim() # 12530    
df.overlapping.tss_and_pro.w.DOWNSTREAM.result %>% count(end.down.distance)
#   end.down.distance     n
# 1              5000  1502
# 2             10000  4069
# 3             25000  6959

########## bind_rows(UPSTREAM & DOWNSTREAM) -> BOTH
df.overlapping.tss_and_pro.w.BOTH.result <- bind_rows(df.overlapping.tss_and_pro.w.UPSTREAM.result %>% 
                                                 mutate(loop.id = up.loop.id, end.distance = end.up.distance) %>% 
                                                 mutate(case.id = tss_and_pro.loop.up.id) %>% 
                                                 dplyr::select(-c(up.loop.id, end.up.distance, tss_and_pro.loop.up.id)), 
                                               df.overlapping.tss_and_pro.w.DOWNSTREAM.result %>% 
                                                 mutate(loop.id = down.loop.id, end.distance = end.down.distance) %>% 
                                                 mutate(case.id = tss_and_pro.loop.down.id) %>% 
                                                 dplyr::select(-c(down.loop.id, end.down.distance, tss_and_pro.loop.down.id))) %>% 
  mutate(chr = ifelse(WHERE == "UP", str_split_n(loop.id, '_', 1), str_split_n(loop.id, '_', 4))) %>%  # utils_functions.R
  mutate(chr = factor(chr, levels = c(paste0("chr", 1:20), "chrX", "chrY"))) %>%
  mutate(WHERE = fct_relevel(WHERE, "UP", "DOWN")) %>% 
  mutate(resolution = convert_to_resolution(end.distance))  # Using utility function # utils_functions.R

df.overlapping.tss_and_pro.w.BOTH.result %>% dim() # 25477
df.overlapping.tss_and_pro.w.BOTH.result %>% colnames() # [1] "distance"     "resolution"   "gene.id"      "component.id" "component"  "WHERE"        "loop.id"      "end.distance" "case.id"      "chr" 
df.overlapping.tss_and_pro.w.BOTH.result %>% count(resolution)
# 1 5K          3101
# 2 10K         8266
# 3 25K        14110

# for boxplot
df.overlapping.tss_and_pro.w.BOTH.result %>% head(2) # istance resolution gene.id  component.id component WHERE loop.id end.distance
df.overlapping.tss_and_pro.w.BOTH.result %>% distinct(loop.id) # 12,349 × 1

# for Q1: quartile among loops with CTCF, so min = 1
df.tss_and_pro.counts <- df.overlapping.tss_and_pro.w.BOTH.result %>%
  group_by(loop.id, WHERE, resolution) %>%
  summarise(tss_and_pro_count = n_distinct(component.id), .groups = 'drop')

df.tss_and_pro.counts %>% count(tss_and_pro_count)

###############
###############
df.tss_and_pro.case <- df.tss_and_pro.counts %>%
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

df.tss_and_pro.case
df.tss_and_pro.case %>% count(case)

df.loop.with.tss_and_pro.case <- df.DISTINCT.loop.deep.sample.all.lt.2mb %>% # 31773
  left_join(df.tss_and_pro.case, by = "loop.id") %>%
  mutate(case = ifelse(is.na(case), "NONE", case))

df.loop.with.tss_and_pro.case %>% count(case)

tss_and_pro.stats.by.resolution <- df.tss_and_pro.counts %>%
  # group_by(resolution) %>%
  group_by(WHERE, resolution) %>%
  summarise(
    Min = min(tss_and_pro_count),
    Q1 = quantile(tss_and_pro_count, 0.25, na.rm = TRUE),
    Median = median(tss_and_pro_count, na.rm = TRUE),
    Q3 = quantile(tss_and_pro_count, 0.75, na.rm = TRUE),
    Mean = mean(tss_and_pro_count, na.rm = TRUE),
    SD = sd(tss_and_pro_count, na.rm = TRUE),
    Max = max(tss_and_pro_count),
    .groups = 'drop'
  )
ctcf.stats.by.resolution




##########################################################
##########################################################
# UP-all
df.final.loop.dataset.tss.unique <- df.final.loop.dataset.tss %>% 
  filter(str_starts(case, "u.")) %>% 
  dplyr::select(loop.id, case.id, WHERE, tss.gene_id) %>% 
  mutate(case.id = str_c(case.id, "TSS", sep = '|')) %>% 
  mutate(gene.id = str_split_n(tss.gene_id, '\\|', 1)) %>% # utils_functions.R
  dplyr::rename(case.where.id = tss.gene_id)

# view()
# colnames() %>%
# head()
df.final.loop.dataset.promoter.unique <- df.final.loop.dataset.promoter %>%
  filter(str_starts(case, "u.")) %>% 
  dplyr::select(loop.id, case.id, WHERE, promoter.gene_id) %>% 
  mutate(case.id = str_c(case.id, "promoter", sep = '|')) %>% 
  mutate(gene.id = str_split_n(promoter.gene_id, '\\|', 1)) %>% # utils_functions.R
  dplyr::rename(case.where.id = promoter.gene_id)
# colnames()
# head()
# filter(is.na(gene.id)) %>% 
# view()

df.final.loop.dataset.pro_tss.unique.all <- bind_rows(df.final.loop.dataset.tss.unique, df.final.loop.dataset.promoter.unique) 
df.final.loop.dataset.pro_tss.unique.all %>% head()
df.final.loop.dataset.pro_tss.unique.all %>% dim() # default: 15701| case6: 16732| case7: 18314| case8: 18957|case9: 19204

###################
# Loops per Gene
###################
gene_loop_counts_all <- df.final.loop.dataset.pro_tss.unique.all %>%
  group_by(gene.id) %>%
  summarise(n_loops = n_distinct(loop.id), .groups = "drop") %>%
  arrange(desc(n_loops))

gene_loop_counts_all 
gene_loop_counts_all %>% count(n_loops)

top50_genes <- gene_loop_counts_all %>% 
  slice_head(n = 50) %>% 
  dplyr::select(gene.id)

print(top50_genes, n = Inf)

gene_list_for_gprofiler <- top50_genes$gene.id
gene_list_for_gprofiler

# TOP50
ggplot(gene_loop_counts_all %>% slice_max(n_loops, n = 50, with_ties = FALSE),
       aes(x = reorder(gene.id, n_loops), y = n_loops)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = n_loops), vjust = -0.3, size = 3.5) +
  labs(title = "Number of Loops per Gene",
       x = "Gene Symbol",
       y = "Loop Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

top50_genes_for_gprofiler <- df.final.loop.dataset.pro_tss.unique.all %>%
  group_by(gene.id) %>%
  summarise(
    n_loops = n_distinct(loop.id),  # gene별 unique loop 수 계산
    .groups = "drop"
  ) %>%
  arrange(n_loops) %>%        # loop 수 내림차순 정렬
  slice_head(n = 50)                # 상위 50개 선택

# 결과 확인
print(top50_genes_for_gprofiler)

# gProfiler용 gene 리스트만 추출 (문자 벡터 형태)
gene_list_for_gprofiler <- top50_genes_for_gprofiler$gene.id
print(gene_list_for_gprofiler)

###################
# Genes per Loop
###################
loop_gene_counts_all <- df.final.loop.dataset.pro_tss.unique.all %>%
  group_by(loop.id) %>%
  summarise(n_genes = n_distinct(gene.id), .groups = "drop")

df.unique_loops %>% head()
loop_gene_counts_all <- df.unique_loops %>%
  group_by(loop.id) %>%
  summarise(n_genes = n_distinct(gene_id), .groups = "drop")

the_number_of_genes_per_loops <- ggplot(loop_gene_counts_all, aes(x = n_genes)) +
  geom_bar(fill = "darkorange") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.3, size = 3.5) +
  scale_x_continuous(breaks = seq(min(loop_gene_counts_all$n_genes),
                                  max(loop_gene_counts_all$n_genes), by = 1)) +
  labs(title = "Distribution of Number of Genes per Loop",
       x = "Number of Unique Genes per Loop",
       y = "Count of Loops") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

saving_plot_dual( # utils_functions.R
  plot_obj = the_number_of_genes_per_loops,
  output_dir = "./figures/submission/lt2mb",
  filename_base = "the_number_of_genes_per_loops")

##########################################################
##########################################################

df.final.loop.dataset.promoter %>% head(3)
df.final.loop.dataset.tss %>% head(3)

##########################################################
# ideogram                                        Figure 4
##########################################################
chromosome_data %>% head()
chromosome_data <- chromosome_data %>% mutate(CE_start = NA, CE_end = NA)
chromosome_data

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
ctcf_density <- process_feature_bins(df_ctcf_ideogram, chromosome_data, bin_size, label = "CTCF") %>% dplyr::select(-last_col()) # color = "#E41A1C" # utils_functions.R
ctcf_density

ideogram(
  karyotype = chromosome_data %>% dplyr::rename(Chr = chr, Start = start, End = end) %>% mutate(Chr = str_remove(Chr, 'chr')),
  overlaid = ctcf_density %>% dplyr::rename(Chr = chr, Start = start, End = end) %>% mutate(Chr = str_remove(Chr, 'chr')),
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
################################################################################################
################################################################################################
# threshold value visualization
################################################################################################
################################################################################################
# 4-1. CTCF
df.ctcf.counts
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

df.loops.above.ctcf.threshold

final.loops.from.ctcf.step <- df.loops.above.ctcf.threshold %>% 
  distinct(loop.id) %>% 
  mutate(end.distance = str_split_n(loop.id, '_', 7)) %>%  # utils_functions.R
  mutate(resolution = case_when(
    end.distance == 5000 ~ "5K",
    end.distance == 10000 ~ "10K",
    end.distance == 25000 ~ "25K",
    TRUE ~ NA
  ))

final.loops.from.ctcf.step
# %>% 
#   mutate(loop.id = str_remove(loop.id, "_[^_]+$"))
final.loops.from.ctcf.step # 18,327 + 10 | sub.4 any, threshold > 6 : 25,620| lt2mb: 25,268
final.loops.from.ctcf.step %>% count(resolution)
# sub.4
# any: 
# new threshold
# lt2mb              # resolution     n
# 1 10K         9375 # 1 10K         9420
# 2 25K        11750 # 2 25K        12049
# 3 5K          4143 # 3 5K          4151


# 4-2. TSS & pro
df.tss_and_pro.counts
df.tss_and_pro.counts %>%
  group_by(resolution, WHERE) %>%
  summarise(max_ctcf = max(tss_and_pro_count), .groups = "drop")


df.tss_and_pro.filtered <- plotting_and_filtering_summary( # utils_functions.R
  df_counts = df.tss_and_pro.counts,
  count_col = "tss_and_pro_count",
  output_prefix = "tss_and_pro"
)
df.tss_and_pro.counts %>% head()
tss_and_pro.stats.by.resolution

####################################
####################################
# Venn Diagram
####################################
####################################
final.loops.from.promoter.step$loop.id
final.loops.from.tss.step$loop.id

venn_plot_submission <- create_venn_plot(final.loops.from.ctcf.step, final.loops.from.promoter.step, final.loops.from.tss.step, "CTCF") # utils_functions.R
venn_plot_submission

saving_plot_dual( # utils_functions.R
  plot_obj = venn_plot_submission,
  filename_base = "ctcf_vs_promoter_tss_venn_diagrams_approach2_lt_2mb_final",
  output_dir = "./figures/submission/lt2mb",
  scale_x = 1,
  scale_y = 1)

####################################
# functional loop extraction
####################################

final.loops.from.ctcf.step %>% dim()# sub.4 any, either 6: 25,620// span ENSEMBL 25268
final.loops.from.tss.step %>% dim() # sub.4 dedups any, either: 27,518// 13407// span ENSEMBL  17810// 17107
final.loops.from.promoter.step %>% dim() # sub.4 dedups any, either: 25,156// span ENSEMBL  13209// 13912

overlapping_loops <- extract_overlapping_loops( # utils_functions.R
  ctcf_data = final.loops.from.ctcf.step,
  promoter_data = final.loops.from.promoter.step,
  tss_data = final.loops.from.tss.step
)

overlapping_loops$ctcf_promoter_only # PASS// 6648
overlapping_loops$ctcf_tss_only # PASS// 8613
overlapping_loops$ctcf_promoter_tss_overlap # PASS// 3423

df.ctcf.promoter.only.loop <- data.frame(loop.id = overlapping_loops$ctcf_promoter_only, category = "CP", stringsAsFactors = FALSE)
df.ctcf.tss.only.loop <- data.frame(loop.id = overlapping_loops$ctcf_tss_only, category = "CT", stringsAsFactors = FALSE)
df.ctcf.promoter.tss.overlap.loop <- data.frame(loop.id = overlapping_loops$ctcf_promoter_tss_overlap, category = "CPT", stringsAsFactors = FALSE)

df.ctcf.promoter.only.loop %>% dim() # 5497
df.ctcf.tss.only.loop %>% dim() # 7799
df.ctcf.promoter.tss.overlap.loop %>% dim() # 2302

# sub.4 dedups: 11891 = 679 + 3096 + 8116
df.final.loop <- bind_rows(df.ctcf.promoter.only.loop, df.ctcf.tss.only.loop, df.ctcf.promoter.tss.overlap.loop)

df.final.loop # 23640/31773(0.7440279) sub.4 any, either 6// 11526/31019(lt2mb) (0.3627608)|EBSEMBL  18684/31019(0.6023405)
df.final.loop %>% head()
df.final.loop %>% dim() # 15598: 5656 + 7980 + 1962

df.final.loop %>% mutate(resolution = str_split_n(loop.id, '_', 7)) %>% count(resolution) # utils_functions.R
# resolution       n  resolution    n     resolution    n 
# 1      10000  86751      10000 5656     1      10000 6968
# 2      25000 111442      25000 7980     2      25000 8617
# 3       5000  38213       5000 1962     3       5000 3099

save(df.final.loop, file="./figures/submission/lt2mb/df_final_loop_sub.4.any.lt2mb.ENSEMBL.rda")
write.csv(df.final.loop, file = "./figures/submission/lt2mb/df_final_loop_sub.4.dedup.any.lt2mb.ENSEMBL.csv", row.names = FALSE)

df.final.loop # 23640// 18681// 18684
df.final.loop %>% head()
df.DISTINCT.loop.deep.sample.all %>% head()

df.final.DISTINCT.loop.joined <- df.final.loop %>% 
  inner_join(df.DISTINCT.loop.deep.sample.all, by = "loop.id")
df.final.DISTINCT.loop.joined %>% dim()
df.final.DISTINCT.loop.joined


df.final.loop %>% distinct(loop.id) # 23460// 18681// 18684

########################
########################
# 1. Loop
# 1-2-2. figure:circos plot for loops
########################
########################
df.final.loop %>% head()
df.circos.final.loop <- df.final.loop %>% 
  separate(loop.id, into = c("chr1", "x1", "x2", "chr2", "y1", "y2", "end.distance"), sep = "_", convert = TRUE, remove = FALSE) %>%
  relocate(chr1, x1, x2, chr2, y1, y2, end.distance, .after = loop.id) %>% 
  mutate(y12 = (y2 + y1)/2, x12 = (x2 + x1)/2) %>% 
  mutate(distance = abs(y12 - x12)) %>% 
  resolution = convert_to_resolution(end.distance)) # utils_functions.R
df.circos.final.loop %>% head()

# Log transform the distance to highlight small distances
df.circos.input.final.loop <- df.circos.final.loop %>% 
  mutate(distance_log = log10(distance + 0.01),
         midx = x12,
         midy = y12 )

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
  mutate(chr1_clean = gsub("chr", "", chr1),
         chr1_clean = factor(chr1_clean, levels = chromosome_order))

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
pdf(file = "./figures/submission/lt2mb/circos_loops_by_resolution_either.6.lt2mb.ENSEMBL.pdf", height = 11*0.8, width = 8.5*0.8)

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
pdf_path <- "./figures/submission/lt2mb/circos_loops_by_resolution_either.6.lt2mb.ENSEMBL.pdf"
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
  filename_base = "circos_first_two_panels_F8.ENSEMBL",
  plot_obj = combined_plot,
  scale_x = 1,
  scale_y = 1
)

################################################################################
# Gene Location-based Filtering and Distance Distribution Analysis
# Goal: 1. Filter loops - remove cases where gene is OUTSIDE loop
#       2. Visualize distance distribution after filtering
#       3. Categorize by loop.id row counts (2 rows, 1 row, 0 rows)
################################################################################

# Step 1: Add gene location information
df.with.loop.coords <- df.final.up.down.tss.pro.nearest %>%
  mutate(loop_parts = str_split(loop.id, "_")) %>%
  mutate(
    chr1 = map_chr(loop_parts, 1),
    x1   = as.numeric(map_chr(loop_parts, 2)),
    x2   = as.numeric(map_chr(loop_parts, 3)),
    chr2 = map_chr(loop_parts, 4),
    y1   = as.numeric(map_chr(loop_parts, 5)),
    y2   = as.numeric(map_chr(loop_parts, 6)),
    anchor_distance = as.numeric(map_chr(loop_parts, 7))
  ) %>%
  dplyr::select(-loop_parts)

cat("Original data rows:", nrow(df.with.loop.coords), "\n") # 62038
cat("Unique loops:", n_distinct(df.with.loop.coords$loop.id), "\n\n") # 31019

# Join with gene boundaries
df.with.gene.location <- df.with.loop.coords %>%
  left_join(
    df.ensembl.gtf.for.exon %>% 
      group_by(gene_id) %>%
      summarise(
        gene_start = min(start),
        gene_end = max(end),
        gene_chr = dplyr::first(chr),
        .groups = "drop"
      ),
    by = "gene_id"
  ) %>%
  mutate(gene_inside_loop = !is.na(gene_start) & 
                            chr1 == gene_chr &
                            gene_start >= x1 & 
                            gene_end <= y2)

cat("Before filtering:\n")
cat("  Total rows:", nrow(df.with.gene.location), "\n") # 62038
cat("  Rows with gene info:", sum(!is.na(df.with.gene.location$gene_start)), "\n") # 60351
cat("  Genes INSIDE loop:", sum(df.with.gene.location$gene_inside_loop, na.rm = TRUE), "\n") # 25478
cat("  Genes OUTSIDE loop:", sum(df.with.gene.location$gene_inside_loop == FALSE, na.rm = TRUE), "\n\n") # 34873

# Step 2: Filter - Keep only genes INSIDE loop
df.gene.inside.only <- df.with.gene.location %>%
  filter(gene_inside_loop == TRUE)

cat("After filtering:\n")
cat("  Total rows:", nrow(df.gene.inside.only), "\n") # 25,478
cat("  Unique loops:", n_distinct(df.gene.inside.only$loop.id), "\n\n") # 18,796

# Step 3: Classify loops by remaining row counts
loop_row_counts <- df.gene.inside.only %>%
  count(loop.id, name = "rows_remaining")

original_loop_ids <- unique(df.with.gene.location$loop.id)

loop_classification <- tibble(loop.id = original_loop_ids) %>%
  left_join(loop_row_counts, by = "loop.id") %>%
  mutate(rows_remaining = replace_na(rows_remaining, 0)) %>%
  mutate(category = case_when(
    rows_remaining == 2 ~ "Both_UP_DOWN_inside",
    rows_remaining == 1 ~ "Only_one_inside",
    rows_remaining == 0 ~ "Both_outside_removed"
  ))

cat("\nLoop Classification Summary:\n")
loop_classification %>%
  count(category) %>%
  mutate(percentage = round(n / sum(n) * 100, 2)) %>%
  arrange(desc(n)) %>%
  print()

df.gene.inside.with.category <- df.gene.inside.only %>%
  left_join(loop_classification %>% dplyr::select(loop.id, category), by = "loop.id")

# Step 4: Distance statistics
overall_stats <- df.gene.inside.only %>%
  summarise(
    min = min(distance),
    Q1 = quantile(distance, 0.25),
    median = median(distance),
    mean = mean(distance),
    Q3 = quantile(distance, 0.75),
    max = max(distance)
  )

cat("\nOverall Distance Statistics (genes inside loop):\n")
print(overall_stats)

category_stats <- df.gene.inside.with.category %>%
  group_by(category) %>%
  summarise(
    n = n(),
    min = min(distance),
    Q1 = quantile(distance, 0.25),
    median = median(distance),
    mean = mean(distance),
    Q3 = quantile(distance, 0.75),
    max = max(distance),
    .groups = "drop"
  )

cat("\nDistance Statistics by Category:\n")
print(category_stats)

# Step 5: Visualizations
# Plot 1: Overall histogram (linear scale)
plot_overall_hist_linear <- ggplot(df.gene.inside.only, aes(x = distance/1000)) +
  geom_histogram(bins = 50, fill = "#2E86AB", color = "white", alpha = 0.8) +
  geom_vline(xintercept = overall_stats$median/1000, 
             linetype = "dashed", color = "#E63946", linewidth = 1.2) +
  geom_vline(xintercept = overall_stats$mean/1000, 
             linetype = "dotted", color = "#F77F00", linewidth = 1.2) +
  annotate("text", x = overall_stats$median/1000, y = Inf,
           label = paste0("Median: ", round(overall_stats$median/1000, 1), " kb"),
           hjust = -0.1, vjust = 1.5, color = "#E63946", fontface = "bold") +
  annotate("text", x = overall_stats$mean/1000, y = Inf,
           label = paste0("Mean: ", round(overall_stats$mean/1000, 1), " kb"),
           hjust = -0.1, vjust = 3, color = "#F77F00", fontface = "bold") +
  labs(
    title = "Distance Distribution (Genes Inside Loop Only)",
    subtitle = paste0("n = ", format(nrow(df.gene.inside.only), big.mark = ","), " rows"),
    x = "Distance from Anchor (kb)",
    y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30")
  )
plot_overall_hist_linear
# Plot 2: Overall histogram (log scale)
plot_overall_hist_log <- ggplot(df.gene.inside.only, aes(x = distance + 1)) +
  geom_histogram(bins = 50, fill = "#A23B72", color = "white", alpha = 0.8) +
  scale_x_log10(labels = label_number(scale_cut = cut_short_scale())) +
  geom_vline(xintercept = overall_stats$median, 
             linetype = "dashed", color = "#E63946", linewidth = 1.2) +
  labs(
    title = "Distance Distribution (Log Scale)",
    subtitle = "Genes inside loop only",
    x = "Distance from Anchor (bp, log10)",
    y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30")
  )
plot_overall_hist_log
# Plot 3: Boxplot
plot_overall_boxplot <- ggplot(df.gene.inside.only, aes(y = distance/1000)) +
  geom_boxplot(fill = "#50C878", color = "#1B5E20", 
               outlier.color = "#E63946", outlier.alpha = 0.5) +
  geom_hline(yintercept = overall_stats$Q1/1000, 
             linetype = "dashed", color = "blue", alpha = 0.7) +
  geom_hline(yintercept = overall_stats$median/1000, 
             linetype = "dashed", color = "darkgreen", alpha = 0.7) +
  geom_hline(yintercept = overall_stats$Q3/1000, 
             linetype = "dashed", color = "purple", alpha = 0.7) +
  annotate("text", x = 1.3, y = overall_stats$Q1/1000, 
           label = paste0("Q1: ", round(overall_stats$Q1/1000, 1), " kb"), 
           color = "blue", hjust = 0) +
  annotate("text", x = 1.3, y = overall_stats$median/1000, 
           label = paste0("Median: ", round(overall_stats$median/1000, 1), " kb"), 
           color = "darkgreen", hjust = 0) +
  annotate("text", x = 1.3, y = overall_stats$Q3/1000, 
           label = paste0("Q3: ", round(overall_stats$Q3/1000, 1), " kb"), 
           color = "purple", hjust = 0) +
  labs(
    title = "Distance Distribution (Boxplot)",
    y = "Distance from Anchor (kb)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )
plot_overall_boxplot
# Plot 4: Histogram by category (linear)
plot_category_hist_linear <- ggplot(df.gene.inside.with.category, 
                                    aes(x = distance/1000, fill = category)) +
  geom_histogram(bins = 40, color = "white", alpha = 0.8) +
  facet_wrap(~category, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(
    "Both_UP_DOWN_inside" = "#2E86AB",
    "Only_one_inside" = "#F77F00"
  )) +
  geom_vline(data = category_stats, 
             aes(xintercept = median/1000), 
             linetype = "dashed", color = "red", linewidth = 1) +
  labs(
    title = "Distance Distribution by Loop Category",
    subtitle = "After filtering genes outside loop",
    x = "Distance from Anchor (kb)",
    y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 11)
  )
plot_category_hist_linear
# Plot 5: Histogram by category (log)
plot_category_hist_log <- ggplot(df.gene.inside.with.category, 
                                 aes(x = distance + 1, fill = category)) +
  geom_histogram(bins = 40, color = "white", alpha = 0.8) +
  facet_wrap(~category, ncol = 1, scales = "free_y") +
  scale_x_log10(labels = label_number(scale_cut = cut_short_scale())) +
  scale_fill_manual(values = c(
    "Both_UP_DOWN_inside" = "#2E86AB",
    "Only_one_inside" = "#F77F00"
  )) +
  labs(
    title = "Distance Distribution by Category (Log Scale)",
    x = "Distance from Anchor (bp, log10)",
    y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 11)
  )
plot_category_hist_log
# Plot 6: Boxplot by category
plot_category_boxplot <- ggplot(df.gene.inside.with.category, 
                                aes(x = category, y = distance/1000, fill = category)) +
  geom_boxplot(outlier.alpha = 0.3, outlier.size = 1) +
  scale_fill_manual(values = c(
    "Both_UP_DOWN_inside" = "#2E86AB",
    "Only_one_inside" = "#F77F00"
  )) +
  labs(
    title = "Distance Comparison by Category",
    x = "Loop Category",
    y = "Distance from Anchor (kb)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    legend.position = "none",
    axis.text.x = element_text(angle = 15, hjust = 1)
  ) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, 
               fill = "red", color = "darkred")
plot_category_boxplot
# Display plots
print(plot_overall_hist_linear)
print(plot_overall_hist_log)
print(plot_overall_boxplot)
print(plot_category_hist_linear)
print(plot_category_hist_log)
print(plot_category_boxplot)

# Save plots
saving_plot_dual( # utils_functions.R
  plot_obj = plot_overall_hist_linear,
  filename_base = "gene_inside_distance_dist_overall_linear",
  output_dir = "./figures/submission/lt2mb"
)

saving_plot_dual( # utils_functions.R
  plot_obj = plot_overall_hist_log,
  filename_base = "gene_inside_distance_dist_overall_log",
  output_dir = "./figures/submission/lt2mb"
)

saving_plot_dual( # utils_functions.R
  plot_obj = plot_overall_boxplot,
  filename_base = "gene_inside_distance_dist_overall_boxplot",
  output_dir = "./figures/submission/lt2mb"
)

saving_plot_dual( # utils_functions.R
  plot_obj = plot_category_hist_linear,
  filename_base = "gene_inside_distance_dist_by_category_linear",
  output_dir = "./figures/submission/lt2mb",
  height_in = 10
)

saving_plot_dual( # utils_functions.R
  plot_obj = plot_category_hist_log,
  filename_base = "gene_inside_distance_dist_by_category_log",
  output_dir = "./figures/submission/lt2mb",
  height_in = 10
)

saving_plot_dual( # utils_functions.R
  plot_obj = plot_category_boxplot,
  filename_base = "gene_inside_distance_dist_by_category_boxplot",
  output_dir = "./figures/submission/lt2mb"
)

################################################################################
# Analysis of Loops with 2 Rows (Both UP and DOWN inside loop)
# Task 1: Count rows by distance pattern (0-0, 0-!0, !0-!0)
# Task 2: Visualize distance difference distribution by pattern
################################################################################
# Step 1: Extract loops with 2 rows (Both_UP_DOWN_inside)
df.both.up.down.inside <- df.gene.inside.with.category %>%
  filter(category == "Both_UP_DOWN_inside")

cat("Loops with both UP and DOWN inside:\n")
cat("  Total rows:", nrow(df.both.up.down.inside), "\n") # 13364
cat("  Unique loops:", n_distinct(df.both.up.down.inside$loop.id), "\n\n") # 6682

# Step 2: Reshape to wide format for distance comparison
df.both.wide <- df.both.up.down.inside %>%
  dplyr::select(loop.id, WHERE, distance, resolution, component, gene_id, gene_name) %>%
  pivot_wider(
    id_cols = c(loop.id, resolution),
    names_from = WHERE,
    values_from = c(distance, component, gene_id, gene_name),
    names_sep = "_"
  )

cat("Wide format data created\n")
cat("  Rows:", nrow(df.both.wide), "\n\n") # 6682

# Step 3: Classify by distance pattern
df.both.classified <- df.both.wide %>%
  mutate(
    distance_pattern = case_when(
      distance_UP == 0 & distance_DOWN == 0 ~ "0-0",
      (distance_UP == 0 & distance_DOWN != 0) | 
      (distance_UP != 0 & distance_DOWN == 0) ~ "0-!0",
      distance_UP != 0 & distance_DOWN != 0 ~ "!0-!0",
      TRUE ~ "other"
    )
  ) %>%
  mutate(
    distance_diff = abs(distance_UP - distance_DOWN),
    min_distance = pmin(distance_UP, distance_DOWN),
    max_distance = pmax(distance_UP, distance_DOWN)
  )

# Task 1: Count rows by pattern
pattern_counts <- df.both.classified %>%
  count(distance_pattern) %>%
  mutate(percentage = round(n / sum(n) * 100, 2)) %>%
  arrange(desc(n))

print(pattern_counts)

cat("\nDetailed breakdown:\n")
cat("  0-0   (both distances = 0):", 
    sum(df.both.classified$distance_pattern == "0-0"), "loops\n") # 1190
cat("  0-!0  (one distance = 0, other != 0):", 
    sum(df.both.classified$distance_pattern == "0-!0"), "loops\n") # 1070
cat("  !0-!0 (both distances != 0):", 
    sum(df.both.classified$distance_pattern == "!0-!0"), "loops\n\n") # 4422

# Additional statistics per pattern
distance_stats_by_pattern <- df.both.classified %>%
  group_by(distance_pattern) %>%
  summarise(
    n = n(),
    mean_UP = mean(distance_UP),
    median_UP = median(distance_UP),
    mean_DOWN = mean(distance_DOWN),
    median_DOWN = median(distance_DOWN),
    mean_diff = mean(distance_diff),
    median_diff = median(distance_diff),
    .groups = "drop"
  )

print(distance_stats_by_pattern)

# Task 2: Visualize distance difference distribution
# Remove 0-0 pattern for difference analysis (difference is always 0)
df.both.for.diff <- df.both.classified %>%
  filter(distance_pattern != "0-0")

cat("Analyzing distance differences (excluding 0-0 pattern):\n")
cat("  Total loops:", nrow(df.both.for.diff), "\n\n") # 5492

# Plot 1: Distance difference histogram by pattern (linear scale)
plot_diff_hist_linear <- ggplot(df.both.for.diff, 
                                aes(x = distance_diff/1000, fill = distance_pattern)) +
  geom_histogram(bins = 50, color = "white", alpha = 0.8) +
  facet_wrap(~distance_pattern, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(
    "0-!0" = "#E63946",
    "!0-!0" = "#2E86AB"
  )) +
  labs(
    title = "Distance Difference Distribution by Pattern",
    subtitle = "Absolute difference between UP and DOWN distances",
    x = "Distance Difference (kb)",
    y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 11)
  )

# Plot 2: Distance difference histogram by pattern (log scale)
plot_diff_hist_log <- ggplot(df.both.for.diff, 
                             aes(x = distance_diff + 1, fill = distance_pattern)) +
  geom_histogram(bins = 50, color = "white", alpha = 0.8) +
  facet_wrap(~distance_pattern, ncol = 1, scales = "free_y") +
  scale_x_log10(labels = label_number(scale_cut = cut_short_scale())) +
  scale_fill_manual(values = c(
    "0-!0" = "#E63946",
    "!0-!0" = "#2E86AB"
  )) +
  labs(
    title = "Distance Difference Distribution (Log Scale)",
    subtitle = "Absolute difference between UP and DOWN distances",
    x = "Distance Difference (bp, log10)",
    y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 11)
  )

# Plot 3: Boxplot comparison
plot_diff_boxplot <- ggplot(df.both.for.diff, 
                            aes(x = distance_pattern, y = distance_diff/1000, 
                                fill = distance_pattern)) +
  geom_boxplot(outlier.alpha = 0.3, outlier.size = 1) +
  scale_fill_manual(values = c(
    "0-!0" = "#E63946",
    "!0-!0" = "#2E86AB"
  )) +
  labs(
    title = "Distance Difference Comparison by Pattern",
    x = "Distance Pattern",
    y = "Distance Difference (kb)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    legend.position = "none"
  ) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
               fill = "yellow", color = "darkred")

# Plot 4: Scatter plot - UP vs DOWN distances
plot_scatter_up_down <- ggplot(df.both.classified, 
                               aes(x = distance_UP/1000, y = distance_DOWN/1000, 
                                   color = distance_pattern)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
              color = "black", linewidth = 1) +
  scale_color_manual(values = c(
    "0-0" = "#50C878",
    "0-!0" = "#E63946",
    "!0-!0" = "#2E86AB"
  )) +
  labs(
    title = "UP vs DOWN Distance Comparison",
    subtitle = "Dashed line represents equal distances",
    x = "UP Distance (kb)",
    y = "DOWN Distance (kb)",
    color = "Pattern"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
    legend.position = "bottom"
  )

# Plot 5: Scatter plot (log scale)
plot_scatter_up_down_log <- ggplot(df.both.classified, 
                                   aes(x = distance_UP + 1, y = distance_DOWN + 1, 
                                       color = distance_pattern)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
              color = "black", linewidth = 1) +
  scale_x_log10(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_log10(labels = label_number(scale_cut = cut_short_scale())) +
  scale_color_manual(values = c(
    "0-0" = "#50C878",
    "0-!0" = "#E63946",
    "!0-!0" = "#2E86AB"
  )) +
  labs(
    title = "UP vs DOWN Distance Comparison (Log Scale)",
    subtitle = "Dashed line represents equal distances",
    x = "UP Distance (bp, log10)",
    y = "DOWN Distance (bp, log10)",
    color = "Pattern"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
    legend.position = "bottom"
  )

# Plot 6: Violin plot for distance difference
plot_diff_violin <- ggplot(df.both.for.diff, 
                           aes(x = distance_pattern, y = distance_diff/1000, 
                               fill = distance_pattern)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.2, fill = "white", alpha = 0.8, outlier.alpha = 0.3) +
  scale_fill_manual(values = c(
    "0-!0" = "#E63946",
    "!0-!0" = "#2E86AB"
  )) +
  labs(
    title = "Distance Difference Distribution (Violin Plot)",
    x = "Distance Pattern",
    y = "Distance Difference (kb)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    legend.position = "none"
  )

# Display plots
print(plot_diff_hist_linear)
print(plot_diff_hist_log)
print(plot_diff_boxplot)
print(plot_scatter_up_down)
print(plot_scatter_up_down_log)
print(plot_diff_violin)

# Save plots
saving_plot_dual( # utils_functions.R
  plot_obj = plot_diff_hist_linear,
  filename_base = "distance_diff_by_pattern_linear",
  output_dir = "./figures/submission/lt2mb",
  height_in = 8
)

saving_plot_dual( # utils_functions.R
  plot_obj = plot_diff_hist_log,
  filename_base = "distance_diff_by_pattern_log",
  output_dir = "./figures/submission/lt2mb",
  height_in = 8
)

saving_plot_dual( # utils_functions.R
  plot_obj = plot_diff_boxplot,
  filename_base = "distance_diff_by_pattern_boxplot",
  output_dir = "./figures/submission/lt2mb"
)

saving_plot_dual( # utils_functions.R
  plot_obj = plot_scatter_up_down,
  filename_base = "distance_up_vs_down_scatter",
  output_dir = "./figures/submission/lt2mb"
)

saving_plot_dual( # utils_functions.R
  plot_obj = plot_scatter_up_down_log,
  filename_base = "distance_up_vs_down_scatter_log",
  output_dir = "./figures/submission/lt2mb"
)

saving_plot_dual( # utils_functions.R
  plot_obj = plot_diff_violin,
  filename_base = "distance_diff_by_pattern_violin",
  output_dir = "./figures/submission/lt2mb"
)
