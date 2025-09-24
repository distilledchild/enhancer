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
setwd('~/dropbox/K\ P/Gateway_to_Hao/enhancer/r_files')
getwd()
source(file.path('~/dropbox/K\ P/Gateway_to_Hao/project_common_code/', 'variables.R'))
source(file.path('~/dropbox/K\ P/Gateway_to_Hao/project_common_code/', 'funcs.R'))

# source(file.path('~/dropbox/K P/Gateway_to_Hao/enhancer/r_files/enhancer_synced_data_preparation.R'))
# source(file.path('~/dropbox/K P/Gateway_to_Hao/enhancer/r_files/data_analysis.R'))
# source(file.path('~/dropbox/K P/Gateway_to_Hao/enhancer/r_files/data_analysis.R'))

####################################
# Sequencing stats
####################################
seq.data <- read.table("../data/library_complexity.tsv", header = TRUE, sep = "\t")
seq.data[ , -1] <- lapply(seq.data[ , -1], function(x) as.numeric(as.character(x)))
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
    x = "Strain",
    y = "Proportion of Read Categories (%)"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(labels = scales::percent)

sequencing_basic_stats

ggsave(
  "figures/submission/sequencing_basic_stats.pdf",
  plot = sequencing_basic_stats,
  device = "pdf",
  width = 11*0.6,
  height = 8.5*0.6,
  units = "in" # letter size
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
loop.file.list = fs::dir_ls("~/dropbox/K\ P/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/", regexp = ".bedpe$")
loop.file.list = fs::dir_ls("/Users/PanjunKim/Desktop/temp/enhancer/dropbox_enhancer_doosan/data/loops/", regexp = ".bedpe$")
loop.file.list = fs::dir_ls("~/Desktop/temp/enhancer/dropbox_enhancer_doosan/data/loops/", regexp = ".bedpe$")
loop.file.list = fs::dir_ls("~/dropbox/K P/Gateway_to_Hao/enhancer/data/loops", regexp = ".bedpe$")
loop.file.list

# BED fild for loops
df.init.loop.bed <- init.bedpe.df(loop.file.list, 'loops') %>% 
  filter(!str_detect(X.chr1, "^#")) 
# %>% #58,992
#   view()

df.init.loop.bed %>% 
  count()

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

df.loop.deep.sample.all
df.loop.deep.sample.all %>% filter(chr1 == 'chrY')

df.loop.deep.sample.all %>% head()
df.loop.deep.sample.all %>% count(strain)
df.loop.deep.sample.all %>% count(resolution)
df.loop.deep.sample.all %>% count(chr1, resolution)

########################
# 1. Loop
# 1-1. Exploratory Data analysis (EDA)
# 1-1.1. data processing
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
pdf("figures/submission/common_loops_heatmap_percentage_bw_strains.pdf", width = 8.5, height = 11)
grid.arrange(plots[["5K"]], plots[["10K"]], plots[["25K"]], nrow = 3)
dev.off()


########################
# 1. Loop
# 1-1-2. figure: shared loops - bar plot
########################
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
  # geom_text(aes(label = paste0("±", round(sd_shared_loops, 1)), 
  #               y = mean_shared_loops + sd_shared_loops + 50),  
  #           vjust = 0) +
  # geom_text(aes(label = round(mean_shared_loops, 1), 
  #               y = mean_shared_loops),  
  #           vjust = -1.5, fontface = "bold") +
  labs(
    # title = "Average Number of Shared Loop by Resolution",
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

ggsave("figures/submission/shared_loops_by_resolution.pdf", 
       plot = shared_loops_by_resolution, 
       device = "pdf", 
       width = 8.5*0.6, 
       height = 11*0.6, 
       units = "in")
########################
# 1. Loop
# 1-1-2. figure: shared loops - network plot
########################
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

ggsave("figures/submission/network_plot_for_shared_loops.pdf", plot = network_plot_for_shared_loops, device = "pdf", width = 11*0.6, height = 8.5*0.6, units = "in")

########################
# 1. Loop
# 1-1-3. figure: shared loops - upset plot
########################
# df.loop.deep.sample.all %>% head()
# 
# # Step 1: loop.id - strain mapping table
# dfoverlap <- df.loop.deep.sample.all %>%
#   dplyr::select(loop.id, strain) %>%
#   distinct() %>%
#   mutate(val = 1)
# dfoverlap
# 
# # Step 2: Convert to wide format
# dfoverlap_wide <- dfoverlap %>%
#   pivot_wider(names_from = strain, values_from = val, values_fill = 0) %>%
#   mutate(
#     chr1 = str_split_n(loop.id, "_", 1),
#     x1   = as.numeric(str_split_n(loop.id, "_", 2)),
#     x2   = as.numeric(str_split_n(loop.id, "_", 3)),
#     chr2 = str_split_n(loop.id, "_", 4),
#     y1   = as.numeric(str_split_n(loop.id, "_", 5)),
#     y2   = as.numeric(str_split_n(loop.id, "_", 6)),
#     end_distance = str_split_n(loop.id, "_", 7),
#     end_distance = str_replace(end_distance, "000$", "K"),
#     distance = abs(((y1 + y2)/2) - ((x1 + x2)/2))
#   )
# 
# strain_cols <- names(dfoverlap_wide)[2:11]
# 
# dfoverlap_wide <- dfoverlap_wide %>% rowwise() %>%
#   mutate(
#     count_1s = sum(c_across(all_of(strain_cols)) == 1),
#     `5K` = ifelse(end_distance == "5K", count_1s, 0),
#     `10K` = ifelse(end_distance == "10K", count_1s, 0),
#     `25K` = ifelse(end_distance == "25K", count_1s, 0)
#   ) %>%
#   ungroup() %>%
#   dplyr::select(-count_1s)
# 
# # test
# # dfoverlap_wide %>% dplyr::select('end_distance', '5K', '10K', '25K') %>% filter(end_distance == '10K')
# 
# # The number of combination for sharing loops BASED ON MY DATASET
# # Step 1: strain presence matrix
# strain_matrix <- dfoverlap_wide %>%
#   dplyr::select(all_of(strain_list))
# 
# # Step 2: binary pattern per loop
# dfoverlap_wide <- dfoverlap_wide %>%
#   mutate(pattern = apply(strain_matrix, 1, function(row) paste0(ifelse(row == 1, names(row), ""), collapse = ",")))
# 
# # Step 3: filter shared loops (more than 1 sample)
# dfoverlap_wide <- dfoverlap_wide %>%
#   mutate(num_strains = rowSums(across(all_of(strain_list)))) %>%
#   filter(num_strains >= 2)
# 
# # Step 4: # of sample combination
# shared_pattern_counts <- dfoverlap_wide %>%
#   count(pattern, sort = TRUE)
# shared_pattern_counts
# 
# # total number
# nrow(shared_pattern_counts)
# 
# 
# 
# strain_list <- names(dfoverlap_wide)[2:11]
# dfoverlap_wide
# 
# # Step 3: Generate upset plot
# upset_plot <- upset(
#   dfoverlap_wide,
#   intersect = strain_list,
#   name = "Shared Loops",
#   width_ratio = 0.15,
#   sort_sets = "descending",
#   sort_intersections_by = "cardinality",
#   n_intersection = 42,
#   base_annotations=list(
#     'Intersection size'=intersection_size(
#       counts=TRUE,
#       mapping=aes(fill='bars_color'),
#       text = list(size = 2, fontface = "bold")
#       ) +
#     scale_fill_manual(values=c('bars_color'=figure_green), guide='none')
#   ),
#   set_sizes = upset_set_size()
#     + ylab('Number of loop')
#     + geom_bar(fill = figure_orange)  # Set (horizontal) bars: orange (UTHSC color)
#     + scale_y_continuous(
#       breaks = c(2500, 5000, 7500),
#       labels = c("2.5K", "5K", "7.5K"))
# )
# 
# upset_plot
# 
# # Step 4: Save to PDF
# ggsave(
#   filename = "./figures/submission/shared_loops_upset_plot.pdf",
#   plot = upset_plot,
#   width = 11, height = 8.5,
#   units = "in"
# )

################################################
# 1. Loop
# 1-1-4. figure: loops by sequencing reads
################################################
df.loop.deep.sample.all %>% 
  count(chr1, resolution) %>% 
  arrange(desc(n))
# Step 1: loops by strain: df.loop.deep.sample.all
loop_counts_by_sample <- df.loop.deep.sample.all %>%
  count(strain) %>% 
  dplyr::rename(num_loop = n)
loop_counts_by_sample
# strain
# SHR/OlaIpcv
# HXB10
# F344/Stm
# BXH6
# Bn-Lx
# HXB23
# SHR/OlaIpcvxBN/NHsdMcwi
# HXB2
# HXB31
seq.data
seq.data <- seq.data %>% dplyr::rename(strain = Strain)
seq.data
# strain    n
# 1         BXH6 6568
# 2        Bn-Lx 6535        
# 3     F344/Stm 2903
# 4        HXB10 7336
# 5         HXB2 4656
# 6        HXB23 7676
# 7        HXB31 9131
# 8       LE/Stm 2992
# 9  SHR/OlaIpcv 5263
# 10        SOBN 5932        


# Step 2: sequencing info & join 
merged_df <- loop_counts_by_sample %>%
  inner_join(seq.data, by = "strain")

merged_df %>% head(500)

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

# p-value and cor
correlation_results <- df_long_for_plot %>%
  group_by(Sequencing_Metric) %>%
  group_modify(~ cor.test(.x$Depth, .x$num_loop) %>% tidy()) %>%
  dplyr::select(Sequencing_Metric, estimate, p.value) %>% 
  # float
  mutate(
    r = round(estimate, 3),
    p = ifelse(p.value < 0.001, "< 0.001", format(round(p.value, 3), nsmall = 3)),
    label = paste0("r = ", r, ", p = ", p)
  )

# Sequencing_Metric estimate  p.value
# <chr>                <dbl>    <dbl>
# 1 Alignable_Reads      0.787 0.00689 
# 2 Total_Reads          0.780 0.00783 
# 3 Unique_Reads         0.884 0.000688

correlation_results

label_positions <- df_long_for_plot %>%
  group_by(Sequencing_Metric) %>%
  summarise(x = max(Depth), .groups = "drop") %>%
  group_modify(~ {
    metric <- .x$Sequencing_Metric[1]
    model <- lm(num_loop ~ Depth, data = df_long_for_plot %>% filter(Sequencing_Metric == metric))
    y_pred <- predict(model, newdata = data.frame(Depth = .x$x))
    .x$y <- y_pred
    .x
  }) %>%
  mutate(x = ifelse(Sequencing_Metric == "Unique Reads", x*0.65, x)) %>%
  mutate(y = ifelse(Sequencing_Metric == "Unique Reads", 8214, y)) %>%
  mutate(x = ifelse(Sequencing_Metric == "Alignable Reads", x*0.86, x)) %>%
  mutate(y = ifelse(Sequencing_Metric == "Alignable Reads", 8314, y)) %>%
  mutate(x = ifelse(Sequencing_Metric == "Total Reads", x*0.78, x)) %>%
  mutate(y = ifelse(Sequencing_Metric == "Total Reads", y*0.54, y))


# label_positions
# correlation_results_position <- correlation_results %>%
#   left_join(label_positions, by = "Sequencing_Metric") %>% 
#   mutate(Sequencing_Metric = str_replace(Sequencing_Metric, '_', ' '))

# correlation_results_position
# Step 4: line graph

library(ggrepel)
library(scales)
# Calculate linear model statistics
correlation_results <- df_long_for_plot %>%
  group_by(Sequencing_Metric) %>%
  do({
    model <- lm(num_loop ~ Depth, data = .)
    data.frame(
      r.squared = summary(model)$r.squared,
      p.value = summary(model)$coefficients["Depth", "Pr(>|t|)"]
    )
  }) %>%
  ungroup() %>%
  mutate(label = paste0("R2 = ", round(r.squared, 2), ", p = ", format(p.value, digits = 2)))

correlation_results

# R^2 to R
correlation_results <- df_long_for_plot %>%
  group_by(Sequencing_Metric) %>%
  group_modify(~ cor.test(.x$Depth, .x$num_loop) %>% broom::tidy()) %>%
  ungroup() %>%
  mutate(
    r = round(estimate, 3),
    p = ifelse(p.value < 0.001, "< 0.001", format(round(p.value, 3), nsmall = 3)),
    label = paste0("r = ", r, ", p = ", p)
  )

correlation_results

# Original labels
original_labels <- levels(factor(df_long_for_plot$Sequencing_Metric))

# Combine original labels with statistics
combined_labels <- paste0(original_labels, " (", correlation_results$label, ")")

line_graph_for_loops_per_depth <- ggplot(df_long_for_plot, aes(x = Depth, y = num_loop, color = Sequencing_Metric)) +
  geom_point() +
  geom_text_repel(aes(label = strain), show.legend = FALSE, max.overlaps = 50) +
  geom_smooth(aes(color = Sequencing_Metric, fill = Sequencing_Metric), method = "lm", se = TRUE, linewidth = 1, linetype = "solid") +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  labs(
    x = "Number of Reads",
    y = "Number of Loops",
    color = "Category",
    fill = "Category" # Add fill to labs
  ) +
  scale_color_discrete(labels = combined_labels) + # Use the combined labels
  scale_fill_discrete(labels = combined_labels) + # and for fill as well
  theme(
    legend.position = c(0.7,0.15)
  )

line_graph_for_loops_per_depth # by hao

ggsave("figures/submission/line_graph_for_loops_per_depth_hao_w_new_label.pdf", plot = line_graph_for_loops_per_depth, device = "pdf", width = 11*0.6, height = 8.5*0.6, units = "in")


########################
# 1. Loop
# 1-1-5. figure: loops per sample by resolution
########################
df.loop.counts <- df.loop.deep.sample.all %>%
  group_by(strain, resolution) %>%
  summarise(n_loops = n_distinct(loop.id), .groups = "drop")

loops_per_sample_by_resolution <- ggplot(df.loop.counts, aes(x = strain, y = n_loops, fill = resolution)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    # title = "Number of Loop per Strain at Resolution",
    x = "Strain",
    y = "Number of Loops"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_manual(values = c("5K" = "#a6cee3", "10K" = "#1f78b4", "25K" = "#1f3a93"))

loops_per_sample_by_resolution

ggsave("figures/submission/loops_per_sample_by_resolution.pdf", plot = loops_per_sample_by_resolution, device = "pdf", width = 11*0.6, height = 8.5*0.6 , units = "in")

########################
# 1. Loop
# 1-1-6. figure: loops per chromosome by resolution
########################
df_chr_loop_counts <- df.loop.deep.sample.all %>%
  group_by(chr1, resolution) %>%
  summarise(n_loops = n_distinct(loop.id), .groups = "drop")

# chr factor (chr1, chr2, ..., chrX)
df_chr_loop_counts$chr1 <- factor(df_chr_loop_counts$chr1, levels = mixedsort(unique(df_chr_loop_counts$chr1)))

df_chr_loop_counts_fig <- ggplot(df_chr_loop_counts, aes(x = chr1, y = n_loops, fill = resolution)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    # title = "Number of Loop per Chromosome at Resolution",
    x = "Chromosome",
    y = "Number of Loops"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_manual(values = c("5K" = "#a6cee3", "10K" = "#1f78b4", "25K" = "#1f3a93"))

df_chr_loop_counts_fig

ggsave("figures/submission/loop_counts_per_chr.pdf", plot = df_chr_loop_counts_fig, device = "pdf", width = 11*0.6, height = 8.5*0.6, units = "in")

########################
# 1. Loop
# 1-2-1. Loop data preprocessing: distinct loops
########################

df.loop.deep.sample.all %>% 
  # head()
  count() # 58992

# common loops between samples picked only one: df.DISTINCT.loop.deep.sample.all
df.DISTINCT.loop.deep.sample.all <- df.loop.deep.sample.all %>%
  dplyr::select(loop.id) %>% # line 180
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

df.DISTINCT.loop.deep.sample.all %>% 
  dim()

df.DISTINCT.loop.deep.sample.all %>% 
  group_by(resolution) %>% 
  summarise(ave = mean(distance)) %>% 
  ungroup() %>% 
  mutate(ave_bin = ave/200)

df.DISTINCT.loop.deep.sample.all %>% 
  count(resolution) # 31773/58992
#head()
# resolution     n
# 1         5K  6680
# 2        10K 12162
# 3        25K 12931

df.DISTINCT.loop.deep.sample.all

########################
# 1. Loop
# 1-2-1. relation between number of Loops and chromosome length
########################

# chromosome_data_old <- read.table(file="../data/rn7_chromosome_length.tsv", sep="\t", head=T) %>%
#   mutate(start = as.numeric(start), end = str_remove_all(end, ",") %>% as.numeric()) %>% 
#   mutate(start = start - 1)
# chromosome_data_old

chromosome_data <- read.table(file="../data/rn7_chromosome_length_from_ucsc.tsv", sep="\t") %>%
  mutate(start = 0) %>%
  dplyr::rename(chr = V1, end = V2) %>% 
  mutate(chr = str_remove(chr, '^chr'))
chromosome_data
chromosome_data %>% head()

df.DISTINCT.loop.deep.sample.all %>% head()

cor_chr_length_loop_data <- df.DISTINCT.loop.deep.sample.all %>% 
  count(chr1) %>% mutate(chr = chr1, chr = str_remove(chr, "chr")) %>% inner_join(chromosome_data, by = 'chr') %>% 
  summarise(cor_test = list(cor.test(n, end, method = "pearson"))) %>%
  pull(cor_test) %>%
  .[[1]]
cor_chr_length_loop_data

########################
# 1. Loop
# 1-3. Loop data preprocessing: Creating GRange Obj. from df.DISTINCT.loop.deep.sample.all (whole, up, down)
########################

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

# up (x1 - x2)
df.DISTINCT.loop.deep.sample.all.up.GR <- create_granges(
  start = df.DISTINCT.loop.deep.sample.all$x1,
  end = df.DISTINCT.loop.deep.sample.all$x2
)
# down (y1 - y2)
df.DISTINCT.loop.deep.sample.all.down.GR <- create_granges(
  start = df.DISTINCT.loop.deep.sample.all$y1,
  end = df.DISTINCT.loop.deep.sample.all$y2
)

########################
# 1. Loop
# 1-4. Loop data preprocessing: padding on loops (1 distance) 
# -> overall.df.DISTINCT.loop.deep.sample.all : OBJECT to be used for the distribution of sth over loops
# x0, y3, new_distance, new.loop.id
########################
df.DISTINCT.loop.deep.sample.all %>% head()

overall.df.DISTINCT.loop.deep.sample.all.1.distance <- df.DISTINCT.loop.deep.sample.all %>% 
  mutate(x12 = (x1 + x2)/2, y12 = (y1 + y2)/2) %>% # middle point of each end
  mutate(x0 = ifelse(x12 - (distance) < 0, 0, x12 - (distance)), y3 = y12 + (distance))  # 1 distance for padding


# new.df.loop.deep.sample.all %>%
# overall.df.DISTINCT.loop.deep.sample.all.1.distance %>% 
# filter(x0 < 0) %>%   
#  head()

# data integrity: PASS
# overall.df.DISTINCT.loop.deep.sample.all %>% mutate(test = ifelse(y3-x0 == 2*distance, TRUE, FALSE)) %>% 
#   filter(test == FALSE) %>% 
#   count(x0)
# count(test)

# source('~/playground/research_uthsc/enhancer/Gateway_to_Hao/enhancer/r_files/enhancer_synced_ctcf.R')
# source('~/playground/research_uthsc/enhancer/Gateway_to_Hao/enhancer/r_files/enhancer_synced_tss.R')
# source('~/playground/research_uthsc/enhancer/Gateway_to_Hao/enhancer/r_files/enhancer_synced_promoter.R')

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
df.init.ctcf<-read.table(file="~/dropbox/K P/Gateway_to_Hao/enhancer/data/ctcf/fimo_2nd_trial_meme/fimo_2nd_trial.txt", header=TRUE, sep="\t")
df.init.ctcf<-read.table(file="~/dropbox/K P/Gateway_to_Hao/enhancer/data/ctcf/submission/E4/fimo_E4_submission_trial.txt", header=TRUE, sep="\t") %>% 
  dplyr::rename(chr = sequence_name, end = stop)

df.init.ctcf %>% dim() # 5767921/.4:6551641
df.init.ctcf %>% head() # chr    start      end strand ctcf_pos
df.init.ctcf %>% count(length)

# strand checking
df.init.ctcf %>%                              # +: 2891072, -: 2876849 = 5767921//+:3267352, -:3284289
  count(strand)
# removing dups including strand
df.init.ctcf %>% 
  distinct(chr, start, end, strand)           # 2701585/5767921 | 3331233/6551641
# removing dups excluding strand
df.init.ctcf %>% 
  distinct(chr, start, end)                   # 2544216/5767921 | 3191859/6551641
# removing dups with all including length
df.init.ctcf %>% 
  distinct()                                  # 2701585/5767921 | 3331233/6551641
# verifying for length column
df.init.ctcf %>% 
  mutate(test = ifelse((end - start) == length, TRUE, FALSE)) %>% 
  count(test)


########################
# 2. CTCF
# 2-2. CTCF data preprocessing: id and dedup GRange Obj. (df.DISTINCT.fimo.2nd.trial.ctcf/ df.DISTINCT.ctcf.2nd.fimo.GR)
########################
# generating id column (long running time)
# df.init.ctcf %>% head() 6551641
# chr     start       end strand length

df.DISTINCT.fimo.2nd.trial.ctcf <- df.init.ctcf %>%
  distinct(chr, start, end) %>% # 2701585 | .4:3191859 ************** NO STRAND INFO
  mutate(start = as.numeric(start)) %>% 
  mutate(end = as.numeric(end)) %>% 
  mutate(ctcf_pos = round((start + end) / 2)) %>% 
  mutate(id = str_c(chr, "_", start, "_", end, "_", ctcf_pos))

df.DISTINCT.fimo.2nd.trial.ctcf %>% dim() # 3191859
df.DISTINCT.fimo.2nd.trial.ctcf # 2701585/5767921 : 0.4683811 | 3191859/6551641 :0.4871847
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

df.DISTINCT.ctcf.2nd.fimo.GR

########################
# 2. CTCF
# 2-3. overall distribution of CTCF on loops
########################
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
mcols(overall.df.DISTINCT.loop.deep.sample.all.GR)$id <- overall.df.DISTINCT.loop.deep.sample.all$loop.id
mcols(overall.df.DISTINCT.loop.deep.sample.all.GR)$resolution <- overall.df.DISTINCT.loop.deep.sample.all$resolution

index.distinct.ctcf.w.overall.whole.loop <- findOverlaps(
  df.DISTINCT.ctcf.2nd.fimo.GR, 
  overall.df.DISTINCT.loop.deep.sample.all.GR, 
  type = "any",
  select = "all"
)

overall.loop.for.ctcf.hits <- subjectHits(index.distinct.ctcf.w.overall.whole.loop)
overall.ctcf.on.loop.hits <- queryHits(index.distinct.ctcf.w.overall.whole.loop)
overall.loop.for.ctcf.hits
overall.ctcf.on.loop.hits
# new.loop.id
df.ctcf.dist.result <- tibble(
  loop.id = overall.df.DISTINCT.loop.deep.sample.all$loop.id[overall.loop.for.ctcf.hits],
  loop.start = overall.df.DISTINCT.loop.deep.sample.all$x0[overall.loop.for.ctcf.hits],
  loop.end = overall.df.DISTINCT.loop.deep.sample.all$y3[overall.loop.for.ctcf.hits],
  loop.res = overall.df.DISTINCT.loop.deep.sample.all$resolution[overall.loop.for.ctcf.hits],
  ctcf.id = df.DISTINCT.fimo.2nd.trial.ctcf$id[overall.ctcf.on.loop.hits],
  ctcf.pos = df.DISTINCT.fimo.2nd.trial.ctcf$ctcf_pos[overall.ctcf.on.loop.hits],
) %>% 
  mutate(across(where(is.numeric), ~ format(., scientific = FALSE)))

df.ctcf.dist.result %>% head()
df.ctcf.dist.result %>% dim() # 40861258/57612728| sub.4: 68014890(any), 68012481(within)

head(df.ctcf.dist.result)

relative.pos.df.ctcf.dist.result <- df.ctcf.dist.result %>%
  mutate(loop.start = as.numeric(loop.start),
         loop.end = as.numeric(loop.end),
         pos_coord = as.numeric(ctcf.pos),
         loop_length = (loop.end - loop.start),
         relative_pos = pos_coord - loop.start,
         value = relative_pos / (loop_length / 3) - 1
  ) %>%
  dplyr::select(loop.id, ctcf.id, value, loop_length, loop.res)

relative.pos.df.ctcf.dist.result %>% 
  head()

########################
# 2. CTCF
# 2-3. overall distribution of CTCF on loops: figures
# 2-3-1. by CHROMOSOME
########################

relative.pos.df.ctcf.dist.result


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
      ylim(c(0,1)) +
      labs(title = paste0("Density of CTCF Found over Loop on Chr", chr),
           x = "Relative Position to Loop",
           y = "Density"
      ) +
      theme(plot.title = element_text(hjust = 0.5)) 
    
    plot1.ctcf.hist <- relative.pos.df.ctcf.dist.result.chr %>% 
      mutate(loop.res.num = as.numeric(gsub("K", "000", as.character(loop.res)))) %>% 
      ggplot(aes(x = value)) +
      geom_histogram(fill = "skyblue", color = "black", alpha = 0.7, bins = 200) +
      labs(title = paste0("Histogram of CTCF Found over Loop on Chr", chr),
           x = "Relative Position to Loop",
           y = "Count"
      ) +
      theme(plot.title = element_text(hjust = 0.5)) 
    
    combined_plot_ctcf <- plot1.ctcf.hist + plot1.ctcf.dens # + can be used instead of |
    plot_ctcf_list[[chr]] <- combined_plot_ctcf
    
    message("END: Processing chromosome: ", chr)
  }, error = function(e) {
    
    message("Error processing chromosome: ", chr)
    message("Error message: ", e$message)
  })
}

pdf("./figures/submission/overall_distribution_of_CTCF_by_chromosome_latest_test.pdf", width = 11*0.8, height = 8.5*0.8)

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
  # geom_histogram(fill = "skyblue", color = "black", alpha = 0.7, bins = 200) +
  # geom_histogram(fill = "skyblue", color = NA, alpha = 0.7, bins = 200) +
  geom_histogram(fill = "skyblue", color = "grey70", alpha = 0.7, bins = 200, linewidth = 0.1) +
  labs(# title = "Histogram of CTCF Found over Loop",
    x = "Relative Position to Loop",
    y = "Count") + 
  theme(plot.title = element_text(hjust = 0.5))
plot.ctcf.hist

# CTCF Density Plot (by resolution combined, right)
plot.ctcf.dens <- ggplot(relative.pos.df.ctcf.dist.result, aes(x = value, color = loop.res, fill = loop.res)) +
  geom_density(alpha = 0.3) +  # transparent
  ylim(c(0, 0.75)) +
  scale_color_manual(values = c("5K" = "#a6cee3", "10K" = "#1f78b4", "25K" = "#1f3a93")) + 
  scale_fill_manual(values = c("5K" = "#a6cee3", "10K" = "#1f78b4", "25K" = "#1f3a93")) +
  labs(# title = "Density of CTCF Found over Loop",
    x = "Relative Position to Loop",
    y = "Density",
    color = "Resolution",
    fill = "Resolution") + 
  theme(plot.title = element_text(hjust = 0.5))
plot.ctcf.dens

pdf("figures/submission/overall_distribution_of_CTCF_by_resolution.pdf", width = 11*0.8, height = 8.5*0.4)

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

pdf("./figures/submission/inner_distance_distribution_boxplots_after_padding.pdf", width = 11*0.8, height = 8.5*0.8)
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
                                              type = "any",
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

df.overlapping.CTCF.w.UPSTREAM.result %>% dim() # 1206765| sub.4 any: 1402928 within: 1399302      
df.overlapping.CTCF.w.UPSTREAM.result %>% head()

df.overlapping.CTCF.w.UPSTREAM.result %>% 
  count(end.up.distance)

# sub.4
# ANY
# end.up.distance      n
# <int>  <int>
# 1            5000 179049
# 2           10000 444448
# 3           25000 779431

# 2. CTCF + DOWNSTREAM (df.DISTINCT.loop.deep.sample.all.down, df.DISTINCT.fimo.2nd.trial.ctcf)
index.distinct.ctcf.w.down.loop <- findOverlaps(df.DISTINCT.ctcf.2nd.fimo.GR, 
                                                df.DISTINCT.loop.deep.sample.all.down.GR, 
                                                type = "any",
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

df.overlapping.CTCF.w.DOWNSTREAM.result %>% dim() # sub.4 any: 1411588 within: 1408116       7
df.overlapping.CTCF.w.DOWNSTREAM.result %>% 
  count(end.down.distance)
# sub.4
# A tibble: 3 × 2
# end.down.distance      n
# any
# <int>  <int>
# 1              5000 182974
# 2             10000 452423
# 3             25000 776191

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

df.overlapping.CTCF.w.BOTH.result %>% dim() # sub.4 any: 2814516/ within: 2807418       
df.overlapping.CTCF.w.BOTH.result %>% names()
df.overlapping.CTCF.w.BOTH.result %>% count(resolution)

# resolution       n
# 1         5K 154950 + 158453 =  313403: integrity PASS
# 2        10K 383562 + 390476 =  774038: integrity PASS
# 3        25K 668253 + 664308 = 1332561 : integrity PASS
# any
# resolution       n
# <fct>        <int>
# 1 5K          362023
# 2 10K         896871
# 3 25K        1555622

# for boxplot
df.overlapping.CTCF.w.BOTH.result %>% head(2) # distance, resolution, ctcf.id, WHERE, loop.id, end.distance, case.id, chr

df.overlapping.CTCF.w.BOTH.result %>% distinct(loop.id) # 31634 / 31773 (df.DISTINCT.loop.deep.sample.all %>% distinct(loop.id) %>% dim())

# for Q1
df.ctcf.counts <- df.overlapping.CTCF.w.BOTH.result %>%
  group_by(loop.id, WHERE, resolution) %>%
  summarise(ctcf_count = n_distinct(ctcf.id), .groups = 'drop')

df.overlapping.CTCF.w.BOTH.result
df.ctcf.counts

df.ctcf.counts %>% 
  filter(ctcf_count == 0)

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

df.ctcf.case %>% count(case)

df.loop.with.ctcf.case <- df.DISTINCT.loop.deep.sample.all %>% # 31773
  left_join(df.ctcf.case, by = "loop.id") %>%
  mutate(case = ifelse(is.na(case), "NONE", case))

df.loop.with.ctcf.case

df.loop.with.ctcf.case %>% count(case)

df.DISTINCT.loop.deep.sample.all %>% dim() # [1] 31773    11
# case     n
# 1 both 31634
# 2   no   139
# > 31634 + 139 = 31773

# proofreading
# check 1. DOWN
down_ctcf_ids <- df.ctcf.counts %>%
  group_by(loop.id) %>%
  filter(all(WHERE == "DOWN")) %>%
  distinct(loop.id)
down_ctcf_ids # 521

# check 2. UP
up_ctcf_ids <- df.ctcf.counts %>%
  group_by(loop.id) %>%
  filter(all(WHERE == "UP")) %>%
  distinct(loop.id)

up_ctcf_ids # 540

# 3. UP & DOWN
both_ctcf_ids <- df.ctcf.counts %>%
  group_by(loop.id) %>%
  filter(all(c("UP", "DOWN") %in% WHERE) & n_distinct(WHERE) == 2) %>%
  distinct(loop.id)

both_ctcf_ids # 30573

# 3 cases: 30573 + 540 + 521 = 31634 + 139(no CTCF on both ends) = TOTAL 31773

###############
###############
# panjun, checkout this section. 
# threshould for ctcf is about 2^2.5, regardless resolution
p.hist.ctcf.loopend.count<-ggplot(df.ctcf.counts, aes(x=log2(ctcf_count)))+
  geom_histogram()+
  facet_wrap(~WHERE+resolution, scales="free_y")

pdf(file="histogram_number_of_ctcf_within_ends_of_loops_by_resolution.pdf", p.hist.ctcf.loopend.count, width=10, height=8)
p.hist.ctcf.loopend.count
dev.off()
2^2.5
# [1] 5.656854
log2(5.656854)
# [1] 2.5
log2(CTCF) = 2.5

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
### any
# # A tibble: 6 × 9
# WHERE resolution   Min    Q1 Median    Q3  Mean    SD   Max
# <fct> <fct>      <int> <dbl>  <dbl> <dbl> <dbl> <dbl> <int>
# 1 UP    5K             1     8     28    40  28.7  22.3   155
# 2 UP    10K            1    14     34    52  37.2  28.7   491
# 3 UP    25K            1    30     51    83  60.4  43.2   635
# 4 DOWN  5K             1     9     29    41  29.3  22.8   357
# 5 DOWN  10K            1    15     34    53  37.9  28.6   282
# 6 DOWN  25K            1    30     51    82  60.1  42.5   389
### within
# WHERE resolution   Min    Q1 Median    Q3  Mean    SD   Max
# <fct> <fct>      <int> <dbl>  <dbl> <dbl> <dbl> <dbl> <int>
# 1 UP    5K             1     8     28    40  28.5  22.2   155
# 2 UP    10K            1    14     34    52  37.1  28.6   491
# 3 UP    25K            1    30     51    83  60.3  43.2   635
# 4 DOWN  5K             1     9     29    40  29.2  22.7   357
# 5 DOWN  10K            1    15     34    53  37.8  28.6   282
# 6 DOWN  25K            1    30     51    82  60.0  42.5   389

# drawing boxplot
df.overlapping.CTCF.w.BOTH.result.boxplot <- df.overlapping.CTCF.w.BOTH.result %>% 
  group_by(chr, loop.id, WHERE, resolution) %>%
  # group_by(loop.id, WHERE, resolution) %>%
  # group_by(chr, loop.id, resolution) %>%
  summarise(ctcf_count_by_loop_id = n_distinct(ctcf.id), .groups = 'drop')

# quantile for EACH END// NOT in a loop : results are loops have ctcfs more than 18 at least in each end
q3_ctcf_count <- quantile(df.overlapping.CTCF.w.BOTH.result.boxplot$ctcf_count_by_loop_id, 0.75, na.rm = TRUE)
q3_ctcf_count # 54 :: integrity PASS | sub.4 62

q1_ctcf_count <- quantile(df.overlapping.CTCF.w.BOTH.result.boxplot$ctcf_count_by_loop_id, 0.25, na.rm = TRUE)
q1_ctcf_count # 16 :: integrity PASS | sub.4 19

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
  stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 1)), position = position_dodge(width = 0.75), vjust = -0.5) + 
  stat_summary(fun.data = function(y) {
    data.frame(y = median(y), label = paste0("Q1: ", round(quantile(y, 0.25), 1)))
  }, geom = "text", aes(label = after_stat(label)), position = position_dodge(width = 0.75), vjust = 6.5) +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("5K" = "#a6cee3", "10K" = "#1f78b4", "25K" = "#1f3a93"), breaks = c("5K", "10K", "25K")) +
  scale_y_continuous(breaks = seq(min(df.overlapping.CTCF.w.BOTH.result.boxplot$ctcf_count_by_loop_id), max(df.overlapping.CTCF.w.BOTH.result.boxplot$ctcf_count_by_loop_id), by = 10)) +
  labs(title = "Boxplot for Number of CTCF by Chromosome and Resolution", 
       x = "Chromosomes", y = "Number of CTCF Sites per Loop")

boxplot.w.CTCF.by.chr.and.res

pdf("figures/submission/boxplot_of_ctcf_count_by_chr_res.pdf", width = 11*0.8, height = 8.5*0.8)
grid.arrange(boxplot.w.CTCF.by.chr.and.res, ncol = 1)
dev.off()

boxplot.w.CTCF.by.res <- df.overlapping.CTCF.w.BOTH.result.boxplot %>% 
  ggplot(aes(x = resolution, y = ctcf_count_by_loop_id, fill = resolution)) +
  geom_boxplot() +
  stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 1)), position = position_dodge(width = 0.75), vjust = -0.5) + 
  stat_summary(fun.data = function(y) {
    data.frame(y = median(y), label = paste0("Q1: ", round(quantile(y, 0.25), 1)))
  }, geom = "text", aes(label = after_stat(label)), position = position_dodge(width = 0.75), vjust = 6.5) +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("5K" = "#a6cee3", "10K" = "#1f78b4", "25K" = "#1f3a93"), breaks = c("5K", "10K", "25K")) +
  scale_y_continuous(breaks = seq(min(df.overlapping.CTCF.w.BOTH.result.boxplot$ctcf_count_by_loop_id), max(df.overlapping.CTCF.w.BOTH.result.boxplot$ctcf_count_by_loop_id), by = 10)) +
  labs(title = "Boxplot for Number of CTCF by Resolution", 
       x = "Resolutions", y = "Number of CTCF Sites per Loop")

boxplot.w.CTCF.by.res

pdf("./figures/submission/boxplot_of_ctcf_count_by_res.pdf", width = 11*0.8, height = 8.5*0.8) #letter size
grid.arrange(boxplot.w.CTCF.by.res, ncol = 1)
dev.off()

# two figures above in one pdf (boxplot.w.CTCF.by.chr.and.res, boxplot.w.CTCF.by.res)
pdf("./figures/submission/boxplot_of_ctcf_count.pdf", width = 11*0.8, height = 8.5*0.8)
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
               position = position_dodge(width = 0.75), vjust = -0.5) +
  stat_summary(fun.data = function(y) {
    data.frame(y = median(y), label = paste0("Q1: ", round(quantile(y, 0.25), 1), "\nQ3: ", round(quantile(y, 0.75), 1)))
  }, geom = "text", aes(label = after_stat(label)), position = position_dodge(width = 0.75), vjust = 1.5) +
  theme_minimal() + 
  theme(
    plot.title = element_text(hjust = 0.5),  # title center
    legend.position = "none"
  ) +
  # scale_fill_manual(values = c("UP" = "yellow", "DOWN" = "purple")) +
  scale_fill_manual(values = c("UP" = "#b2df8a", "DOWN" = "#33a02c")) +
  scale_y_continuous(breaks = seq(min(df.overlapping.CTCF.w.BOTH.result.boxplot.ALL.chr$ctcf_count_by_loop_id), max(df.overlapping.CTCF.w.BOTH.result.boxplot.ALL.chr$ctcf_count_by_loop_id), by = 5)) +  # y축 눈금 5단위
  coord_cartesian(ylim = c(quantile(df.overlapping.CTCF.w.BOTH.result.boxplot.ALL.chr$ctcf_count_by_loop_id, 0.05), 
                           quantile(df.overlapping.CTCF.w.BOTH.result.boxplot.ALL.chr$ctcf_count_by_loop_id, 0.95))) +  # y-axis range
  labs(#title = "Boxplot for CTCF by Up/Downstream End", 
    x = "Upstream and Downstream Ends", y = "Number of CTCF Sites")

boxplot.w.CTCF.ALL.chr

# figure by end & res
boxplot.w.CTCF.ALL.chr.by.res <- ggplot(df.overlapping.CTCF.w.BOTH.result.boxplot.ALL.chr, aes(x = WHERE, y = ctcf_count_by_loop_id, fill = resolution)) +
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +  # removing outliers
  stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 1)), 
               position = position_dodge(width = 0.75), vjust = -0.5) +
  stat_summary(fun.data = function(y) {
    data.frame(y = median(y), label = paste0("Q1: ", round(quantile(y, 0.25), 1), "\nQ3: ", round(quantile(y, 0.75), 1)))
  }, geom = "text", aes(label = after_stat(label)), position = position_dodge(width = 0.75), vjust = 1.5) +
  theme_minimal() + 
  theme(
    plot.title = element_text(hjust = 0.5),  # title center
    legend.position = "none"  # removing legend
  ) +
  scale_fill_manual(values = c("5K" = "#a6cee3", "10K" = "#1f78b4", "25K" = "#1f3a93"), breaks = c("5K", "10K", "25K")) +
  scale_y_continuous(breaks = seq(min(df.overlapping.CTCF.w.BOTH.result.boxplot.ALL.chr$ctcf_count_by_loop_id), 
                                  max(df.overlapping.CTCF.w.BOTH.result.boxplot.ALL.chr$ctcf_count_by_loop_id), by = 5)) +  # y-axis tip
  coord_cartesian(ylim = c(quantile(df.overlapping.CTCF.w.BOTH.result.boxplot.ALL.chr$ctcf_count_by_loop_id, 0.05), 
                           quantile(df.overlapping.CTCF.w.BOTH.result.boxplot.ALL.chr$ctcf_count_by_loop_id, 0.95))) + 
  labs(#title = "Boxplot for CTCF in Up/Downstream End by Resolution", 
    x = "Upstream and Downstream Ends", y = "Number of CTCF Sites", fill = "Resolution")


boxplot.w.CTCF.ALL.chr.by.res

pdf("./figures/submission/boxplot_for_no.ctcf_by_end_res.pdf", width = 8.5*0.8, height = 11*0.8)
grid.arrange(boxplot.w.CTCF.ALL.chr, boxplot.w.CTCF.ALL.chr.by.res, ncol = 1)
dev.off()

##########################################################################################
##########################################################################################
##########################################################################################
# refactoring code
##########################################################################################
##########################################################################################
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
  index_ctcf_w_up_loop <- findOverlaps(ctcf_data, loop_up_GR, type = "any")
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
  index_ctcf_w_down_loop <- findOverlaps(ctcf_data, loop_down_GR, type = "any")
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
                 position = position_dodge(width = 0.75), vjust = -0.5) + 
    stat_summary(fun.data = function(y) {
      data.frame(y = median(y), label = paste0("Q1: ", round(quantile(y, 0.25), 1)))
    }, geom = "text", aes(label = after_stat(label)), 
    position = position_dodge(width = 0.75), vjust = 6.5) +
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = c("5K" = "#a6cee3", "10K" = "#1f78b4", "25K" = "#1f3a93"), breaks = c("5K", "10K", "25K")) +
    scale_y_continuous(breaks = seq(min_value, max_value, by = 10)) +
    labs(title = paste("Boxplot for Number of CTCF by Chromosome and Resolution (", padding_label, " Padding)", sep = ""), 
         x = "Chromosomes", y = "Number of CTCF in a loop")
  
  # by Resolution only
  boxplot_by_res <- df_boxplot %>% 
    ggplot(aes(x = resolution, y = ctcf_count_by_loop_id, fill = resolution)) +
    geom_boxplot() +
    stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 1)), 
                 position = position_dodge(width = 0.75), vjust = -0.5) + 
    stat_summary(fun.data = function(y) {
      data.frame(y = median(y), label = paste0("Q1: ", round(quantile(y, 0.25), 1)))
    }, geom = "text", aes(label = after_stat(label)), 
    position = position_dodge(width = 0.75), vjust = 6.5) +
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = c("5K" = "#a6cee3", "10K" = "#1f78b4", "25K" = "#1f3a93"), breaks = c("5K", "10K", "25K")) +
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



#########################################################
# 1. Loops
# 1-5. diagram with loops from 3 previous steps
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

########################
# 3. TSS
# 3-1. Exploratory Data analysis (EDA) : TSS
# 3-2. TSS data preprocessing: id and dedup GRange Obj.
########################
# tss file list
# Linux
file.tss.list = fs::dir_ls("/home/pkim/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/", regexp = "\\.txt$")
# Mac
file.tss.list = fs::dir_ls("~/dropbox/K P/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss", regexp = "\\.txt$")

file.tss.list
# csRNA.NuAcc.tss.txt
# csRNA.PFC.tss.txt
# ucsc_start_codon.txt

# GTF 
df.refgene.gtf <- read_tsv("~/dropbox/K P/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/ucsc_refGene.gtf",
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
# 5 start_codon  17849
# 6 stop_codon   17806
# 7 transcript   18570

df.refgene.gtf.parsed <- df.refgene.gtf %>%
  mutate(
    gene_id = gsub('.*gene_id "([^"]+)".*', '\\1', attribute),
    transcript_id = gsub('.*transcript_id "([^"]+)".*', '\\1', attribute),
    gene_name = gsub('.*gene_name "([^"]+)".*', '\\1', attribute)
  )

df.refgene.gtf.parsed.start.codon <- df.refgene.gtf.parsed %>%
  filter(feature == "start_codon") %>%
  dplyr::select(chr, start, end, strand, gene_id, transcript_id, gene_name) 
# %>%
#   head()

df.refgene.gtf.parsed.start.codon # 17849
df.refgene.gtf.parsed.start.codon %>% count(strand)


##### tss with UCSC
# get the TSS location into a GenomicRange object, note gene name is added, also there are some duplicated lines, so use uniq
# Linux
tss<-read.table(file="/home/pkim/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/ucsc_start_codon.txt", sep="\t", head=F)
# Mac
tss<-read.table(file="~/dropbox/K P/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/ucsc_start_codon.txt", sep="\t", head=F)
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
tss.select.sparate %>% dplyr::rename(chr = V1, start = V4, end = V5, strand = V7) %>% head()

df.tss.ucsc.tss.id <- tss.select.sparate %>% 
  dplyr::rename(chr = V1, start = V4, end = V5, strand = V7 ) %>% 
  dplyr::mutate(tss.id = paste(chr, start, end, strand, gene_id, transcript_id, sep = ":"))

df.tss.ucsc.tss.id %>% head()
df.tss.ucsc.tss.id %>% count() # 17849
df.tss.ucsc.tss.id.distinct <- df.tss.ucsc.tss.id %>% distinct(chr, start, end) %>% # 17080
  mutate(tss.id = paste(chr, as.character(start), as.character(end), sep = ":"))
df.tss.ucsc.tss.id.distinct.strand <- df.tss.ucsc.tss.id %>% distinct(chr, start, end, strand) %>%  # 17080
  mutate(tss.id = paste(chr, as.character(start), as.character(end), sep = ":"))

df.tss.ucsc.tss.id %>% distinct(chr, start, end, strand, transcript_id) # 17849
df.tss.ucsc.tss.id %>% distinct(chr, start, end, strand, transcript_id, exon_id) # 17849
tss.id.common <- df.tss.ucsc.tss.id.distinct %>% 
  semi_join(df.tss.ucsc.tss.id.distinct.strand, by = c("tss.id" = "tss.id"))
# tss.id.common # 17080
# df.tss.ucsc.tss.id.distinct = df.tss.ucsc.tss.id.distinct.strand 

df.tss.ucsc.tss.id.distinct
df.tss.ucsc.tss.id.distinct.strand
df.tss.ucsc.tss.id %>% distinct(chr, start, end, gene_id) %>% dim() # 17139
df.tss.ucsc.tss.id %>% distinct(chr, start, end, strand, gene_id) # 17139
df.tss.ucsc.tss.id %>% distinct(chr, start, end, strand, gene_id, gene_name) # 17139

########

# 기본 카운트
n_coord      <- df.tss.ucsc.tss.id %>% distinct(chr, start, end) %>% nrow()
n_coord_gene <- df.tss.ucsc.tss.id %>% distinct(chr, start, end, gene_id) %>% nrow()
delta <- n_coord_gene - n_coord
delta

multi_gene_per_coord <- df.tss.ucsc.tss.id %>%
  distinct(chr, start, end, gene_id) %>%
  group_by(chr, start, end) %>%
  summarise(
    n_genes = n_distinct(gene_id),
    gene_ids = paste(sort(unique(gene_id)), collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(n_genes > 1) %>%
  arrange(desc(n_genes), chr, start, end)

# 2) 위 케이스들의 총 초과 개수가 실제 차이와 일치하는지 검증
#    (각 좌표에서 n_genes-1의 합이 delta와 같아야 함)
sum_excess <- multi_gene_per_coord %>% summarise(sum(n_genes - 1)) %>% pull()
all.equal(sum_excess, delta)

# 3) 예시 몇 개만 확인
head(multi_gene_per_coord, 20)

# 4) 만약 동일 좌표·동일 gene_id가 데이터에 중복 존재했는지(집계 전 원시 중복) 확인
raw_dups_same_gene <- df.tss.ucsc.tss.id %>%
  count(chr, start, end, gene_id, name = "n") %>%
  filter(n > 1) %>%
  arrange(desc(n))
head(raw_dups_same_gene)

# 5) gene_id가 비어 있는 데이터가 있는지 참고(숫자 차이를 만들진 않지만 품질 점검용)
has_na_gene <- df.tss.ucsc.tss.id %>% filter(is.na(gene_id)) %>% nrow()
has_na_gene

########
map_gene <- df.tss.ucsc.tss.id %>%
  group_by(chr, start, end) %>%
  summarise(
    gene_id = {
      ids <- unique(gene_id)
      ids <- ids[!is.na(ids)]                 # NA 제외
      if (length(ids) == 0) NA_character_ else paste(sort(ids), collapse = "]")
    },
    .groups = "drop"
  )

map_gene %>% filter(str_detect(gene_id, "\\]") == TRUE) # 54

df.tss.ucsc <- df.tss.ucsc.tss.id %>% distinct(chr, start, end) %>% # 17080******
  left_join(map_gene, by = c("chr", "start", "end")) %>%
  mutate(tss.id = paste(chr, start, end, gene_id, sep = ":"))
df.tss.ucsc # 17080
df.tss.ucsc %>% filter(str_detect(gene_id, "_"))
df.tss.ucsc %>% dim()

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
# mcols(df.tss.ucsc.GR)$transcript_id <- df.tss.ucsc$transcript_id

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
  type = "any",
  select = "all"
)

index.distinct.tss.w.overall.whole.loop # any: 359130/ within: 359129

overall.loop.for.tss.hits <- subjectHits(index.distinct.tss.w.overall.whole.loop)
overall.tss.on.loop.hits <- queryHits(index.distinct.tss.w.overall.whole.loop)

df.tss.dist.result <- tibble(
  loop.id = overall.df.DISTINCT.loop.deep.sample.all$loop.id[overall.loop.for.tss.hits],
  loop.start = overall.df.DISTINCT.loop.deep.sample.all$x0[overall.loop.for.tss.hits],
  loop.end = overall.df.DISTINCT.loop.deep.sample.all$y3[overall.loop.for.tss.hits],
  loop.res=overall.df.DISTINCT.loop.deep.sample.all$resolution[overall.loop.for.tss.hits],
  tss_chr = df.tss.ucsc$chr[overall.tss.on.loop.hits],
  tss_start = df.tss.ucsc$start[overall.tss.on.loop.hits],
  tss_end = df.tss.ucsc$end[overall.tss.on.loop.hits],
  tss_id = df.tss.ucsc$tss.id[overall.tss.on.loop.hits],
  # tss_geneid = df.tss.ucsc$gene_id[overall.tss.on.loop.hits],
  tss_strand = df.tss.ucsc$strand[overall.tss.on.loop.hits]
) %>% 
  mutate(across(where(is.numeric), ~ format(., scientific = FALSE)))

df.tss.dist.result %>% dim() # any: 359130  9/ within: 359129      9
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
      labs(title = paste0("Density of TSS Found over Loop on Chr", chr),
           x = "Relative Position to Loop",
           y = "Density"
      ) + 
      theme(plot.title = element_text(hjust = 0.5))
    
    plot1.tss.hist <- relative.pos.df.tss.dist.result.chr %>% 
      ggplot(aes(x = value)) +
      geom_histogram(fill = "skyblue", color = "black", alpha = 0.7, bins=200) +
      labs(title = paste0("Histogram of TSS Found over Loop on Chr", chr),
           x = "Relative Position to Loop",
           y = "Count"
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

pdf("./figures/submission/overall_distribution_of_TSS_by_chromosome.pdf", width = 11*0.8, height = 8.5*0.8)
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
  # geom_histogram(fill = "skyblue", color = "black", alpha = 0.7, bins=200) +
  geom_histogram(fill = "skyblue", color = "grey70", alpha = 0.7, bins = 200, linewidth = 0.1) +
  labs(# title = "Histogram of TSS Found over Loop",
    x = "Relative Position to Loop",
    y = "Count"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

plot.tss.hist

# TSS Density Plot (by resolution combined, right)
plot.tss.dens <- ggplot(relative.pos.df.tss.dist.result, aes(x = value, color = loop.res, fill = loop.res)) +
  geom_density(alpha = 0.3) +  # transparent
  ylim(c(0, 0.75)) +
  scale_color_manual(values = c("5K" = "#a6cee3", "10K" = "#1f78b4", "25K" = "#1f3a93")) + 
  scale_fill_manual(values = c("5K" = "#a6cee3", "10K" = "#1f78b4", "25K" = "#1f3a93")) +
  labs(# title = "Density of TSS Found over Loop",
    x = "Relative Position to Loop",
    y = "Density",
    color = "Resolution",
    fill = "Resolution") + 
  theme(plot.title = element_text(hjust = 0.5))

plot.tss.dens

pdf("figures/submission/overall_distribution_of_TSS_by_resolution.pdf", width = 11*0.8, height = 8.5*0.4)

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
df.DISTINCT.loop.deep.sample.all
chromosome_data
chrom_ends2 <- chromosome_data %>%
  transmute(chr2 = paste0("chr", as.character(chr)),
            chr2_end = end)
chrom_ends2

count_each_filter <- function(df) {
  list(
    # capping
    x0_eq_0 = df %>% filter(x0 == 0) %>% nrow(),
    y3_raw_gt_chr2_end = df %>% filter(y3_raw > chr2_end) %>% nrow(),
    y3_eq_chr2_end = df %>% filter(y3 == chr2_end) %>% nrow(),
    delta_x_neq_delta_y = df %>% filter(delta_x != delta_y) %>% nrow(),
    # crossing
    x3_gt_y0_raw = df %>% filter(x3 > y0_raw) %>% nrow(),
    x3_eq_y0 = df %>% filter(x3 == y0) %>% nrow()
  )
}

############################################
# padding case draft: L 1/2 outer, L 1/4 inner
############################################
df.DISTINCT.loop.deep.sample.all.padded.for.TSS <- df.DISTINCT.loop.deep.sample.all %>% 
  left_join(chrom_ends2, by = "chr2") %>%
  mutate(x_mid = (x1 + x2)/2, # middle point of each end
         y_mid = (y1 + y2)/2, # middle point of each end
         distance = y_mid - x_mid,
         x0 = pmax(0, x_mid - (distance * 0.5)), # distance for OUTER PADDING
         x3 = (x_mid + (distance*0.25)), # distance for INNER PADDING
         y0_raw = y_mid - (distance * 0.25),
         y0 = pmax(x3, y0_raw), # distance distance for INNER PADDING
         y3_raw = y_mid + (distance * 0.5),
         y3 = pmin(y3_raw, chr2_end), # distance for OUTER PADDING
         new.loop.id = paste0(chr1, '_', x0, '_', x3, '_', chr2, '_', y0, '_', y3, '_', distance),
         delta_x  = x3 - x0,
         delta_y  = y3 - y0,
         inner.distance = y0 - x3)

# x0 = ifelse(x_mid - (distance / 2) < 0, 0, x_mid - (distance / 2)), # 1/2 distance for OUTER PADDING
# y3 = y_mid + (distance/2), # 1/2 distance for OUTER PADDING
# x3 = (x_mid + (distance / 4)), # 1/4 distance for INNER PADDING
# y0 = ifelse((y_mid - (distance/4)) < 0, 0, y_mid - (distance/4)), # 1/4 distance for INNER PADDING
# new.loop.id = paste0(chr1, '_', x0, '_', x3, '_', chr2, '_', y0, '_', y3, '_', end.distance))
df.DISTINCT.loop.deep.sample.all.padded.for.TSS %>% filter(y3 < y0)

print(count_each_filter(df.DISTINCT.loop.deep.sample.all.padded.for.TSS))

# quantile
ylim_vals_TSS <- df.DISTINCT.loop.deep.sample.all.padded.edo1_edi1.for.TSS %>%
  pull(inner.distance) %>%
  quantile(c(0.05, 0.95), na.rm = TRUE)
ylim_vals_TSS

# min & median 
stats_TSS <- df.DISTINCT.loop.deep.sample.all.padded.edo1_edi1.for.TSS %>%
  summarise(
    min_val = min(inner.distance, na.rm = TRUE),
    median_val = median(inner.distance, na.rm = TRUE)
  )
stats_TSS

############################################
# padding by resolution
############################################
df.DISTINCT.loop.deep.sample.all.padded.1trial.for.TSS <- df.DISTINCT.loop.deep.sample.all %>% 
  left_join(chrom_ends2, by = "chr2") %>%
  mutate(x_mid = (x1 + x2)/2, # middle point of each end
         y_mid = (y1 + y2)/2, # middle point of each end
         e     = end.distance,
         pad_scale = case_when(
           resolution == "25K" ~ 0.5,   # 원래 루프 (중점 ± 0.5e)
           resolution == "10K" ~ 0.75,  # half resolution 추가
           resolution == "5K"  ~ 1.5,   # 1x resolution 추가
           TRUE ~ 0.5                   # default (25K처럼)
         ),
         x0 = pmax(0, x_mid - pad_scale * e),
         x3 = x_mid + pad_scale * e,
         y0_raw = y_mid - pad_scale * e,
         y0 = pmax(0, y0_raw),
         y3_raw = y_mid + pad_scale * e,
         y3 = pmin(y3_raw, chr2_end),
         new.loop.id = paste0(chr1, '_', x0, '_', x3, '_', chr2, '_', y0, '_', y3, '_', end.distance),
         delta_x  = x3 - x0,
         delta_y  = y3 - y0,
         inner.distance = y0 - x3)

result_counts_1trial_TSS <- count_each_filter(df.DISTINCT.loop.deep.sample.all.padded.1trial.for.TSS)
print(result_counts_1trial_TSS)

# quantile
ylim_vals_1trial_TSS <- df.DISTINCT.loop.deep.sample.all.padded.1trial.for.TSS %>%
  pull(inner.distance) %>%
  quantile(c(0.05, 0.95), na.rm = TRUE)
ylim_vals_1trial_TSS

# min & median 
stats_1trial_TSS <- df.DISTINCT.loop.deep.sample.all.padded.1trial.for.TSS %>%
  summarise(
    min_val = min(inner.distance, na.rm = TRUE),
    median_val = median(inner.distance, na.rm = TRUE)
  )
stats_1trial_TSS
############################################
# padding case 2: ed(end.distance) outer, ed 1/2 inner
############################################
df.DISTINCT.loop.deep.sample.all.padded.edo1_edi.5.for.TSS <- df.DISTINCT.loop.deep.sample.all %>% 
  left_join(chrom_ends2, by = "chr2") %>%
  mutate(x_mid = (x1 + x2)/2, # middle point of each end
         y_mid = (y1 + y2)/2, # middle point of each end
         x0 = pmax(0, x_mid - (end.distance * 1.5)), # end.distance for OUTER PADDING
         x3 = (x_mid + (end.distance)), # end.distance for INNER PADDING
         y0_raw = y_mid - (end.distance),
         y0 = pmax(x3, y0_raw), # end.distance distance for INNER PADDING
         y3_raw = y_mid + (end.distance * 1.5),
         y3 = pmin(y3_raw, chr2_end), # end.distance for OUTER PADDING
         new.loop.id = paste0(chr1, '_', x0, '_', x3, '_', chr2, '_', y0, '_', y3, '_', end.distance),
         delta_x  = x3 - x0,
         delta_y  = y3 - y0,
         inner.distance = y0 - x3)

result_counts_1.5_TSS <- count_each_filter(df.DISTINCT.loop.deep.sample.all.padded.edo1_edi.5.for.TSS)
print(result_counts_1.5_TSS)

# quantile
ylim_vals_1.5_TSS <- df.DISTINCT.loop.deep.sample.all.padded.edo1_edi.5.for.TSS %>%
  pull(inner.distance) %>%
  quantile(c(0.05, 0.95), na.rm = TRUE)
ylim_vals_1.5_TSS

# min & median 
stats_1.5_TSS <- df.DISTINCT.loop.deep.sample.all.padded.edo1_edi.5.for.TSS %>%
  summarise(
    min_val = min(inner.distance, na.rm = TRUE),
    median_val = median(inner.distance, na.rm = TRUE)
  )
stats_1.5_TSS

############################################
# padding case 3: ed(end.distance) 1/2 outer, ed inner
############################################

df.DISTINCT.loop.deep.sample.all.padded.edo.5_edi1.for.TSS <- df.DISTINCT.loop.deep.sample.all %>% 
  left_join(chrom_ends2, by = "chr2") %>%
  mutate(x_mid = (x1 + x2)/2, # middle point of each end
         y_mid = (y1 + y2)/2, # middle point of each end
         x0 = pmax(0, x_mid - (end.distance)), # end.distance for OUTER PADDING
         x3 = (x_mid + (end.distance*1.5)), # end.distance for INNER PADDING
         y0_raw = y_mid - (end.distance*1.5),
         y0 = pmax(x3, y0_raw), # end.distance distance for INNER PADDING
         y3_raw = y_mid + (end.distance),
         y3 = pmin(y3_raw, chr2_end), # end.distance for OUTER PADDING
         new.loop.id = paste0(chr1, '_', x0, '_', x3, '_', chr2, '_', y0, '_', y3, '_', end.distance),
         delta_x  = x3 - x0,
         delta_y  = y3 - y0,
         inner.distance = y0 - x3)

result_counts_.5.1_TSS <- count_each_filter(df.DISTINCT.loop.deep.sample.all.padded.edo.5_edi1.for.TSS)
print(result_counts_.5.1_TSS)

# quantile
ylim_vals_.5.1_TSS <- df.DISTINCT.loop.deep.sample.all.padded.edo.5_edi1.for.TSS %>%
  pull(inner.distance) %>%
  quantile(c(0.05, 0.95), na.rm = TRUE)
ylim_vals_.5.1_TSS

# min & median 
stats_.5.1_TSS <- df.DISTINCT.loop.deep.sample.all.padded.edo.5_edi1.for.TSS %>%
  summarise(
    min_val = min(inner.distance, na.rm = TRUE),
    median_val = median(inner.distance, na.rm = TRUE)
  )
stats_.5.1_TSS

############################################
# padding case 4: ed(end.distance) 1/2 outer, ed 1/2 inner
############################################
df.DISTINCT.loop.deep.sample.all.padded.edo.5_edi.5.for.TSS <- df.DISTINCT.loop.deep.sample.all %>% 
  left_join(chrom_ends2, by = "chr2") %>%
  mutate(x_mid = (x1 + x2)/2, # middle point of each end
         y_mid = (y1 + y2)/2, # middle point of each end
         x0 = pmax(0, x_mid - (end.distance)), # end.distance for OUTER PADDING
         x3 = (x_mid + (end.distance)), # end.distance for INNER PADDING
         y0_raw = y_mid - (end.distance),
         y0 = pmax(x3, y0_raw), # end.distance distance for INNER PADDING
         y3_raw = y_mid + (end.distance),
         y3 = pmin(y3_raw, chr2_end), # end.distance for OUTER PADDING
         new.loop.id = paste0(chr1, '_', x0, '_', x3, '_', chr2, '_', y0, '_', y3, '_', end.distance),
         delta_x  = x3 - x0,
         delta_y  = y3 - y0,
         inner.distance = y0 - x3)

result_counts_.5.5_TSS <- count_each_filter(df.DISTINCT.loop.deep.sample.all.padded.edo.5_edi.5.for.TSS)
print(result_counts_.5.5_TSS)
result_counts_.5.5_TSS

# quantile
ylim_vals_.5.5_TSS <- df.DISTINCT.loop.deep.sample.all.padded.edo.5_edi.5.for.TSS %>%
  pull(inner.distance) %>%
  quantile(c(0.05, 0.95), na.rm = TRUE)
ylim_vals_.5.5_TSS

# min & median 
stats_.5.5_TSS <- df.DISTINCT.loop.deep.sample.all.padded.edo.5_edi.5.for.TSS %>%
  summarise(
    min_val = min(inner.distance, na.rm = TRUE),
    median_val = median(inner.distance, na.rm = TRUE)
  )
stats_.5.5_TSS

############################################
# padding case 5: ed(end.distance) 1/2 outer, ed 1/4 inner
############################################
df.DISTINCT.loop.deep.sample.all.padded.edo.5_edi.25.for.TSS <- df.DISTINCT.loop.deep.sample.all %>% 
  left_join(chrom_ends2, by = "chr2") %>%
  mutate(x_mid = (x1 + x2)/2, # middle point of each end
         y_mid = (y1 + y2)/2, # middle point of each end
         x0 = pmax(0, x_mid - (end.distance)), # end.distance for OUTER PADDING
         x3 = (x_mid + (end.distance*0.75)), # end.distance for INNER PADDING
         y0_raw = y_mid - (end.distance*0.75),
         y0 = pmax(x3, y0_raw), # end.distance distance for INNER PADDING
         y3_raw = y_mid + (end.distance),
         y3 = pmin(y3_raw, chr2_end), # end.distance for OUTER PADDING
         new.loop.id = paste0(chr1, '_', x0, '_', x3, '_', chr2, '_', y0, '_', y3, '_', end.distance),
         delta_x  = x3 - x0,
         delta_y  = y3 - y0,
         inner.distance = y0 - x3)

result_counts_.5.25_TSS <- count_each_filter(df.DISTINCT.loop.deep.sample.all.padded.edo.5_edi.25.for.TSS)
print(result_counts_.5.25_TSS)
result_counts_.5.25_TSS

# quantile
ylim_vals_.5.25_TSS <- df.DISTINCT.loop.deep.sample.all.padded.edo.5_edi.25.for.TSS %>%
  pull(inner.distance) %>%
  quantile(c(0.05, 0.95), na.rm = TRUE)
ylim_vals_.5.25_TSS

# min & median 
stats_.5.25_TSS <- df.DISTINCT.loop.deep.sample.all.padded.edo.5_edi.25.for.TSS %>%
  summarise(
    min_val = min(inner.distance, na.rm = TRUE),
    median_val = median(inner.distance, na.rm = TRUE)
  )
stats_.5.25_TSS

# padding case draft: L 1/2 outer, L 1/4 inner:::: df.DISTINCT.loop.deep.sample.all.padded.for.TSS
# padding case 1: ed(end.distance) outer, ed inner:::: df.DISTINCT.loop.deep.sample.all.padded.edo1_edi1.for.TSS
# padding case 2: ed(end.distance) outer, ed 1/2 inner:::: df.DISTINCT.loop.deep.sample.all.padded.edo1_edi.5.for.TSS
# padding case 3: ed(end.distance) 1/2 outer, ed inner:::: df.DISTINCT.loop.deep.sample.all.padded.edo.5_edi1.for.TSS
# padding case 4: ed(end.distance) 1/2 outer, ed 1/2 inner:::: df.DISTINCT.loop.deep.sample.all.padded.edo.5_edi.5.for.TSS
# padding case 5: ed(end.distance) 1/2 outer, ed 1/4 inner:::: df.DISTINCT.loop.deep.sample.all.padded.edo.5_edi.25.for.TSS

plot_inner_distance_boxplot <- function(df, ylim_vals, stats, title_text) {
  df %>%
    mutate(inner.distance = y0 - x3) %>%
    ggplot(aes(y = inner.distance)) +
    geom_boxplot(outlier.shape = NA, fill = "#4C9F70", color = "black") +
    geom_hline(yintercept = stats$min_val, linetype = "dashed", color = "blue", linewidth = 0.6) +
    geom_hline(yintercept = stats$median_val, linetype = "dotted", color = "red", linewidth = 0.8) +
    annotate("text", x = 0.2, y = stats$min_val, label = paste0("Min: ", stats$min_val), color = "blue", hjust = 0) +
    annotate("text", x = 0.2, y = stats$median_val, label = paste0("Median: ", stats$median_val), color = "red", hjust = 0) +
    coord_cartesian(ylim = ylim_vals) +
    labs(
      title = title_text,
      y = "Inner Distance (bp)",
      x = ""
    ) +
    theme_minimal(base_size = 14)
}

plot_inner_distance_boxplot(df.DISTINCT.loop.deep.sample.all.padded.for.TSS,
                            ylim_vals_TSS, stats_TSS,
                            "Padding Draft: L 1/2 Outer, L 1/4 Inner")

plot_inner_distance_boxplot(df.DISTINCT.loop.deep.sample.all.padded.1trial.for.TSS,
                            ylim_vals_1trial_TSS, stats_1trial_TSS,
                            "Padding Case 1: ed 1 Outer, ed 1 Inner")

plot_inner_distance_boxplot(df.DISTINCT.loop.deep.sample.all.padded.edo1_edi.5.for.TSS,
                            ylim_vals_1.5_TSS, stats_1.5_TSS,
                            "Padding Case 2: ed 1 Outer, ed 1/2 Inner")

plot_inner_distance_boxplot(df.DISTINCT.loop.deep.sample.all.padded.edo.5_edi1.for.TSS,
                            ylim_vals_.5.1_TSS, stats_.5.1_TSS,
                            "Padding Case 3: ed 1/2 Outer, ed 1 Inner")

plot_inner_distance_boxplot(df.DISTINCT.loop.deep.sample.all.padded.edo.5_edi.5.for.TSS,
                            ylim_vals_.5.5_TSS, stats_.5.5_TSS,
                            "Padding Case 4: ed 1/2 Outer, ed 1/2 Inner")

plot_inner_distance_boxplot(df.DISTINCT.loop.deep.sample.all.padded.edo.5_edi.25.for.TSS,
                            ylim_vals_.5.25_TSS, stats_.5.25_TSS,
                            "Padding Case 5: ed 1/2 Outer, ed 1/4 Inner")

filtered_names <- ls() %>%
  grep("^(stats|ylim).*_TSS$", ., value = TRUE)
filtered_objects <- mget(filtered_names)
str(filtered_objects)
########################
# 3. TSS
# 3-4. Distribution of TSS at ends in a loop: for the number of TSS used in filtering valid loops: figures
# so, the object should be used one WITH padding: df.DISTINCT.loop.deep.sample.all.padded.for.TSS
# 3-4-3. data processing for getting information of TSS at ends in loops 
########################
df.DISTINCT.loop.deep.sample.all.padded.for.TSS <- df.DISTINCT.loop.deep.sample.all.padded.edo1_edi1.for.TSS
df.DISTINCT.loop.deep.sample.all.padded.for.TSS <- df.DISTINCT.loop.deep.sample.all.padded.1trial.for.TSS ##################
df.DISTINCT.loop.deep.sample.all.padded.for.TSS <- df.DISTINCT.loop.deep.sample.all.padded.edo1_edi.5.for.TSS
df.DISTINCT.loop.deep.sample.all.padded.for.TSS <- df.DISTINCT.loop.deep.sample.all.padded.edo.5_edi1.for.TSS
df.DISTINCT.loop.deep.sample.all.padded.for.TSS <- df.DISTINCT.loop.deep.sample.all.padded.edo.5_edi.5.for.TSS
df.DISTINCT.loop.deep.sample.all.padded.for.TSS <- df.DISTINCT.loop.deep.sample.all.padded.edo.5_edi.25.for.TSS



# 1. TSS + UPSTREAM (df.DISTINCT.loop.deep.sample.all.padded.for.TSS.UP.GR, df.tss.ucsc.GR)
df.DISTINCT.loop.deep.sample.all.padded.for.TSS.UP.GR <- GRanges(
  seqnames=df.DISTINCT.loop.deep.sample.all.padded.for.TSS$chr1, 
  ranges=IRanges(
    start=df.DISTINCT.loop.deep.sample.all.padded.for.TSS$x0, 
    end=df.DISTINCT.loop.deep.sample.all.padded.for.TSS$x3
  )
)
mcols(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.UP.GR) <- df.DISTINCT.loop.deep.sample.all.padded.for.TSS %>%
  dplyr::select(loop.id, end.distance, resolution, new.loop.id)

index.distinct.tss.w.up.loop.each.end <- findOverlaps(
  df.tss.ucsc.GR, 
  df.DISTINCT.loop.deep.sample.all.padded.for.TSS.UP.GR, 
  type = "any",
  select = "all"
)

index.distinct.tss.w.up.loop.each.end # any: 91559// 16251// 14443// 13285// 11477// 10508
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

df.overlapping.TSS.w.UPSTREAM.result.each.end %>% dim() # any: 91559     6 // 16251// 14443// 13285
df.overlapping.TSS.w.UPSTREAM.result.each.end %>% head()

df.overlapping.TSS.w.UPSTREAM.result.each.end %>% 
  count(end.up.distance)

# any
# end.up.distance     n
# 1            5000 20414
# 2           10000 27145
# 3           25000 44000
# total: 91559 (dedups) //95831 (dups)

# 2. TSS + DOWNSTREAM
df.DISTINCT.loop.deep.sample.all.padded.for.TSS
df.DISTINCT.loop.deep.sample.all.padded.for.TSS.DOWN.GR <- GRanges(
  seqnames=df.DISTINCT.loop.deep.sample.all.padded.for.TSS$chr1, 
  ranges=IRanges(
    start=df.DISTINCT.loop.deep.sample.all.padded.for.TSS$y0, 
    end=df.DISTINCT.loop.deep.sample.all.padded.for.TSS$y3
  )
)

mcols(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.DOWN.GR) <- df.DISTINCT.loop.deep.sample.all.padded.for.TSS %>%
  dplyr::select(loop.id, end.distance, resolution, new.loop.id)

# 1. TSS + DOWNSTREAM (df.DISTINCT.loop.deep.sample.all.padded.for.TSS.DOWN.GR, df.tss.ucsc.GR)
index.distinct.tss.w.down.loop.each.end <- findOverlaps(
  df.tss.ucsc.GR, 
  df.DISTINCT.loop.deep.sample.all.padded.for.TSS.DOWN.GR, 
  type = "any",
  select = "all"
)
index.distinct.tss.w.down.loop.each.end # any: 106271, within: 106269 // 15731// 14002// 12815// 11086// 10097

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

df.overlapping.TSS.w.DOWNSTREAM.result.each.end %>% dim() # 106271 (dedup) // 111334 (dups) // 15731// 14002// 12815// 11086
df.overlapping.TSS.w.DOWNSTREAM.result.each.end %>% head()
df.overlapping.TSS.w.DOWNSTREAM.result.each.end %>% 
  count(end.down.distance)
# any
# end.down.distance     n
# 1              5000 24273
# 2             10000 30197
# 3             25000 51801
# total: 106271 (dedup) // 111334 (dups)

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

df.overlapping.TSS.w.BOTH.result %>% dim() # any: 197830      7 (dedup) // 31982// 28445 //22563// 20605
df.overlapping.TSS.w.BOTH.result %>% head()
df.overlapping.TSS.w.BOTH.result %>% count(resolution)
# any
# resolution      n
# 1 5K         44687
# 2 10K        57342
# 3 25K        95801


# for Q1
df.tss.counts <- df.overlapping.TSS.w.BOTH.result %>%
  group_by(loop.id, WHERE, resolution) %>%
  summarise(tss_count = n_distinct(tss.id), .groups = 'drop')

df.tss.counts
df.tss.counts %>% filter(tss_count == 0)

###############
###############
df.tss.case <- df.tss.counts %>%
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

df.tss.case %>% count(case)

df.loop.with.tss.case <- df.DISTINCT.loop.deep.sample.all %>% # 31773
  left_join(df.tss.case, by = "loop.id") %>%
  mutate(case = ifelse(is.na(case), "NONE", case))

df.loop.with.tss.case %>% count(case)
# case     n
# 1 BOTH 17046
# 2 DOWN  5131
# 3 NONE  4255
# 4   UP  5341

# 4255/31773 : 13%

df.DISTINCT.loop.deep.sample.all %>% dim() # [1] 31773    11

# proofreading
# check 1. DOWN
down_tss_ids <- df.tss.counts %>%
  group_by(loop.id) %>%
  filter(all(WHERE == "DOWN")) %>%
  distinct(loop.id)
down_tss_ids # 5131

# check 2. UP
up_tss_ids <- df.tss.counts %>%
  group_by(loop.id) %>%
  filter(all(WHERE == "UP")) %>%
  distinct(loop.id)

up_tss_ids # 5341

# 3. UP & DOWN
both_tss_ids <- df.tss.counts %>%
  group_by(loop.id) %>%
  filter(all(c("UP", "DOWN") %in% WHERE) & n_distinct(WHERE) == 2) %>%
  distinct(loop.id)

both_tss_ids # 17046

# 3 cases: 5131 + 5341 + 17046 = 27518 + 4255(no TSS on both ends) = TOTAL 31773

###############
###############
df.tss.counts
# panjun, checkout this section. 
# threshould for ctcf is about 2^2.5, regardless resolution
p.hist.tss.loopend.count<-ggplot(df.tss.counts, aes(x=log2(tss_count)+1))+
  geom_histogram()+
  facet_wrap(~WHERE+resolution, scales="free_y")

pdf(file="histogram_number_of_tss_within_ends_of_loops_by_resolution_padded.1trial.pdf", p.hist.tss.loopend.count, width=10, height=8)
p.hist.tss.loopend.count
dev.off()

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

df.overlapping.TSS.w.BOTH.result %>% distinct(loop.id) # 27,518

df.tss.counts.each.end <- df.overlapping.TSS.w.BOTH.result %>%
  group_by(loop.id, WHERE, resolution) %>%
  summarise(tss_count_each_end = n_distinct(tss.id), .groups = 'drop')

df.tss.counts.each.end %>% filter(tss_count_each_end == 0)

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

# TODO: min value check
tss.stats.by.resolution.each.end

# resolution   Min    Q1 Median    Q3  Mean    SD   Max
# <fct>      <int> <dbl>  <dbl> <dbl> <dbl> <dbl> <int>
# 1 5K             1     1      2     3  5.33  35.8  1282
# 2 10K            1     1      2     3  3.59  18.0   993
# 3 25K            1     1      2     4  5.26  28.3  1282

# any
# dedup
# WHERE resolution   Min    Q1 Median    Q3  Mean    SD   Max
# <fct> <fct>      <int> <dbl>  <dbl> <dbl> <dbl> <dbl> <int>
# 1 UP    5K             1     1      2     3  4.58  27.3   869
# 2 UP    10K            1     1      2     3  3.23  11.9   406
# 3 UP    25K            1     1      2     4  4.62  20.4   869
# 4 DOWN  5K             1     1      2     3  5.59  39.7  1219
# 5 DOWN  10K            1     1      2     3  3.64  21.0   933
# 6 DOWN  25K            1     1      2     4  5.44  32.1  1219

# quantile for EACH END// NOT in a loop : results are loops have ctcfs more than 18 at least in each end
q3_tss_count_each_end <- quantile(df.overlapping.TSS.w.BOTH.result.boxplot$tss_count_by_loop_id_each_end, 0.75, na.rm = TRUE)
q3_tss_count_each_end # 3 :: integrity PASS

q1_tss_count_each_end <- quantile(df.overlapping.TSS.w.BOTH.result.boxplot$tss_count_by_loop_id_each_end, 0.25, na.rm = TRUE)
q1_tss_count_each_end # 1 :: integrity PASS

# drawing boxplot
df.overlapping.TSS.w.BOTH.result.boxplot # any: 44564 (dedup)
df.overlapping.TSS.w.BOTH.result.boxplot %>% head(4)


# figure by chr & res
boxplot.w.TSS.by.chr.and.res.each.end <- df.overlapping.TSS.w.BOTH.result.boxplot %>% 
  ggplot(aes(x = chr, y = tss_count_by_loop_id_each_end, fill = resolution)) +
  geom_boxplot() +
  stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 1)), 
               position = position_dodge(width = 0.75), vjust = -0.5, color = "black") +
  stat_summary(fun.data = function(y) {
    data.frame(
      y = quantile(y, 0.25),
      label = round(quantile(y, 0.25), 1)
    )
  }, geom = "text", aes(label = after_stat(label)), position = position_dodge(width = 0.75), vjust = 2.5, color = "blue") +
  stat_summary(fun.data = function(y) {
    data.frame(
      y = quantile(y, 0.75),
      label = round(quantile(y, 0.75), 1)
    )
  }, geom = "text", aes(label = after_stat(label)), position = position_dodge(width = 0.75), vjust = -2.5, color = "blue") +
  stat_summary(fun.data = function(y) {
    data.frame(
      y = min(y),
      label = round(min(y), 1)
    )
  }, geom = "text", aes(label = after_stat(label)), position = position_dodge(width = 0.75), vjust = 1.5, hjust = -0.5, color = "red") +
  stat_summary(fun.data = function(y) {
    data.frame(
      y = max(y),
      label = round(max(y), 1)
    )
  }, geom = "text", aes(label = after_stat(label)), position = position_dodge(width = 0.75), vjust = -1.5, hjust = -0.5, color = "red") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("5K" = "#a6cee3", "10K" = "#1f78b4", "25K" = "#1f3a93"), breaks = c("5K", "10K", "25K")) +
  scale_y_continuous(limits = c(0, 2 * q3_tss_count_each_end), 
                     breaks = seq(0, 2 * q3_tss_count_each_end, by = 1)) +
  labs(#title = "Boxplot for Number of TSS by Chromosome and Resolution", 
    x = "Chromosome", y = "Number of TSS per Loop")

boxplot.w.TSS.by.chr.and.res.each.end

pdf("figures/submission/end_boxplot_num_tss_by_chr_res.pdf", width = 11*0.8, height = 8.5*0.8)
grid.arrange(boxplot.w.TSS.by.chr.and.res.each.end, ncol = 1)
dev.off()

df.overlapping.TSS.w.BOTH.result.boxplot %>% 
  head()

boxplot.w.TSS.by.res <- df.overlapping.TSS.w.BOTH.result.boxplot %>% 
  ggplot(aes(x = resolution, y = tss_count_by_loop_id_each_end, fill = resolution)) +
  geom_boxplot() +
  stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 1)), 
               position = position_dodge(width = 0.75), vjust = -0.5) + 
  stat_summary(fun.data = function(y) {
    data.frame(y = quantile(y, 0.25), label = paste0("Q1: ", round(quantile(y, 0.25), 1)))
  }, geom = "text", aes(label = after_stat(label)), 
  position = position_dodge(width = 0.75), vjust = 6.5) +
  stat_summary(fun.data = function(y) {
    data.frame(y = quantile(y, 0.75), label = paste0("Q3: ", round(quantile(y, 0.75), 1)))
  }, geom = "text", aes(label = after_stat(label)), 
  position = position_dodge(width = 0.75), vjust = -6.5, color = "blue") +
  stat_summary(fun.data = function(y) {
    data.frame(y = min(y), label = paste0("Min: ", round(min(y), 1)))
  }, geom = "text", aes(label = after_stat(label)), 
  position = position_dodge(width = 0.75), vjust = 1.5, color = "red") +
  stat_summary(fun.data = function(y) {
    data.frame(y = max(y), label = paste0("Max: ", round(max(y), 1)))
  }, geom = "text", aes(label = after_stat(label)), 
  position = position_dodge(width = 0.75), vjust = -1.5, color = "red") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("5K" = "#a6cee3", "10K" = "#1f78b4", "25K" = "#1f3a93"), 
                    breaks = c("5K", "10K", "25K")) +
  scale_y_continuous(limits = c(0, 2 * q3_tss_count_each_end), 
                     breaks = seq(0, 2 * q3_tss_count_each_end, by = 10)) +
  labs(# title = "Boxplot for Number of TSS by Resolution", 
    x = "Resolution", y = "Number of TSS per Loop")

boxplot.w.TSS.by.res

pdf("./figures/submission/end_boxplot_num_tss_by_res.pdf", width = 11*0.8, height = 8.5*0.8)
grid.arrange(boxplot.w.TSS.by.res, ncol = 1)
dev.off()

pdf("./figures/submission/boxplot_of_tss_count.pdf", width = 11*0.8, height = 8.5*0.8)
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
               position = position_dodge(width = 0.75), vjust = -0.5) +
  stat_summary(fun.data = function(y) {
    data.frame(y = median(y), label = paste0("Q1: ", round(quantile(y, 0.25), 1), "\nQ3: ", round(quantile(y, 0.75), 1)))
  }, geom = "text", aes(label = after_stat(label)), position = position_dodge(width = 0.75), vjust = 1.5) +  
  theme_minimal() + 
  theme(
    plot.title = element_text(hjust = 0.5),  # title center
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("UP" = "#b2df8a", "DOWN" = "#33a02c")) +
  scale_y_continuous(limits = c(0, 2 * q3_tss_count_each_end),  
                     breaks = seq(0, 2 * q3_tss_count_each_end, by = 5)) +
  labs(# title = "Boxplot for TSS by Up/Downstream End", 
    x = "Upstream and Downstream Ends", y = "Number of TSS Sites")

boxplot.w.TSS.ALL.chr

# figure by end & res
boxplot.w.TSS.ALL.chr.by.res <- ggplot(df.overlapping.TSS.w.BOTH.result.boxplot, aes(x = WHERE, y = tss_count_by_loop_id_each_end, fill = resolution)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 1)), 
               position = position_dodge(width = 0.75), vjust = -0.5) +
  stat_summary(fun.data = function(y) {
    data.frame(y = median(y), label = paste0("Q1: ", round(quantile(y, 0.25), 1), "\nQ3: ", round(quantile(y, 0.75), 1)))
  }, geom = "text", aes(label = after_stat(label)), position = position_dodge(width = 0.75), vjust = 1.5) +  
  stat_summary(fun.data = function(y) {
    data.frame(y = min(y), label = paste0("Min: ", round(min(y), 1)))
  }, geom = "text", aes(label = after_stat(label)), position = position_dodge(width = 0.75), vjust = 1.5, hjust = -0.2, color = "red") +
  theme_minimal() + 
  theme(
    plot.title = element_text(hjust = 0.5),  # title center
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("5K" = "#a6cee3", "10K" = "#1f78b4", "25K" = "#1f3a93"), breaks = c("5K", "10K", "25K")) +
  scale_y_continuous(limits = c(0, 2 * q3_tss_count_each_end), 
                     breaks = seq(0, 2 * q3_tss_count_each_end, by = 5)) +
  labs(# title = "Boxplot for TSS in Up/Downstream End by Resolution", 
    x = "Upstream and Downstream Ends", y = "Number of TSS Sites", fill = "Resolution")
boxplot.w.TSS.ALL.chr.by.res

pdf("./figures/submission/boxplot_for_no.tss_by_end_res.pdf", width = 8.5*0.8, height = 11*0.8)
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
       y = "Number of Loop") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5) +
          stat_count(aes(label = paste0(..count.., " (", round(..count../total_distinct_loops * 100, 1), "%)")), 
                     geom = "text", 
                     vjust = -0.5) +
          annotate("text", x = 5, y = max(table(tss_count_histogram_data$tss_count_total)) * 1.1, 
                   label = paste("Total Loops:", total_distinct_loops), hjust = 1) +
          annotate("text", x = 5, y = max(table(tss_count_histogram_data$tss_count_total)) * 1.05, 
                   label = paste("Loops > 5 TSS: ", more_than_5_loops, " (", more_than_5_loops_percent, "%)"), 
                   hjust = 1, color = "red")
        
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
               y = "Number of Loop") +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5)) +
          geom_text(aes(label = paste0(one_sided_tss_loop_count, "/", loop_count, " (", round(ratio * 100, 1), "%)")),
                    vjust = -0.5, color = "black")
        
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
        # file_path <- "/Users/PanjunKim/UTHSC GGI Dropbox/K P/Gateway_to_Hao/enhancer/data/epdnew/001/Rn_EPDnew_001_rn7.bed"
        file_path <- "~/dropbox/K P/Gateway_to_Hao/enhancer/data/epdnew/001/Rn_EPDnew_001_rn7_1.bed"
        df.promoter.rn7.raw <- read_tsv(file_path, col_names = c("chr", "start", "end", "gene", "score"))
        df.promoter.rn7.raw # 12,463
        
        df.promoter.rn7 <- df.promoter.rn7.raw %>% 
          distinct(chr, start, end) %>% 
          mutate(length = end - start, 
                 center = round((end + start)/2), 
                 promoter.id = str_c(chr, '_', start, '_', end))
        
        df.promoter.rn7 # 12427 (dedup)
        
        ###############
        # 1) distinct로 '첫 번째'만 남긴 버전(참고용)
        df_promoter_distinct <- df.promoter.rn7.raw %>%
          distinct(chr, start, end, .keep_all = TRUE)
        
        # 2) distinct에서 '제거된' 행(중복의 2번째 이후만)
        dropped_rows <- df.promoter.rn7.raw %>%
          group_by(chr, start, end) %>%
          mutate(ord = row_number()) %>%
          filter(ord > 1) %>%
          ungroup() %>%
          dplyr::select(-ord)
        
        nrow(dropped_rows)         # 기대: 36
        head(dropped_rows, 10)
        
        # 3) 어떤 좌표에서 몇 개씩 중복인지 요약
        dup_summary <- df.promoter.rn7.raw %>%
          count(chr, start, end, name = "n") %>%
          filter(n > 1) %>%
          arrange(desc(n), chr, start, end)
        
        dup_summary
        sum(dup_summary$n - 1)     # 기대: 36 (총 초과 개수 = 제거된 행 수)
        
        # 4) 좌표별로 중복된 gene 표기를 한 줄로 묶어 보기(진단용)
        dup_with_genes <- df.promoter.rn7.raw %>%
          semi_join(dup_summary, by = c("chr","start","end")) %>%
          group_by(chr, start, end) %>%
          summarise(
            n = n(),
            genes = str_c(sort(unique(gene)), collapse = "]"),
            .groups = "drop"
          ) %>%
          arrange(desc(n), chr, start, end)
        
        dup_with_genes
        #########################
        
        gene_map <- df.promoter.rn7.raw %>%
          group_by(chr, start, end) %>%
          summarise(
            gene_id = {
              g <- unique(gene)
              g <- g[!is.na(g)]
              if (length(g) == 0) NA_character_ else str_c(sort(g), collapse = "]")
            },
            .groups = "drop"
          )
        gene_map %>% filter(str_detect(gene_id, "\\]") == TRUE)
        df.promoter.rn7
        # 2) df.promoter.rn7에 gene 열 붙이기
        df.promoter.rn7 <- df.promoter.rn7 %>%
          left_join(gene_map, by = c("chr","start","end")) %>% 
          mutate(promoter.id = str_c(promoter.id, "_", gene_id))
        df.promoter.rn7
        df.promoter.rn7 %>% filter(str_detect(gene_id, "\\]") == TRUE)
        
        
        # dedup GRange Obj.
        df.promoter.rn7.GR <- GRanges(
          seqnames = df.promoter.rn7$chr,
          ranges = IRanges(
            start = df.promoter.rn7$start, 
            end = df.promoter.rn7$end)
        )
        df.promoter.rn7.GR
        
        # metadata
        mcols(df.promoter.rn7.GR) <- df.promoter.rn7[, c("promoter.id", "center", "length", "gene_id")]
        
        df.promoter.rn7.GR
        ########################
        # 4. promoter
        # 4-3. overall distribution of promoter on loops
        ########################
        index.promoter.w.overall.whole.loop <- findOverlaps(
          df.promoter.rn7.GR, 
          overall.df.DISTINCT.loop.deep.sample.all.GR, 
          type = "any",
          select = "all"
        )
        index.promoter.w.overall.whole.loop # any: 264932 // within: 264906
        
        overall.loop.for.promoter.hits <- subjectHits(index.promoter.w.overall.whole.loop)
        overall.promoter.on.loop.hits <- queryHits(index.promoter.w.overall.whole.loop)
        
        df.promoter.dist.result <- data.frame(
          loop.id = overall.df.DISTINCT.loop.deep.sample.all$loop.id[overall.loop.for.promoter.hits],
          loop.start = overall.df.DISTINCT.loop.deep.sample.all$x0[overall.loop.for.promoter.hits],
          loop.end = overall.df.DISTINCT.loop.deep.sample.all$y3[overall.loop.for.promoter.hits],
          loop.res = overall.df.DISTINCT.loop.deep.sample.all$resolution[overall.loop.for.promoter.hits],
          # loop.new.distance = overall.df.DISTINCT.loop.deep.sample.all$new_distance[overall.loop.for.promoter.hits],
          promoter_id = df.promoter.rn7$promoter.id[overall.promoter.on.loop.hits],
          promoter_start = df.promoter.rn7$start[overall.promoter.on.loop.hits],
          promoter_end = df.promoter.rn7$end[overall.promoter.on.loop.hits],
          promoter_length = df.promoter.rn7$length[overall.promoter.on.loop.hits]
        ) %>% 
          mutate(across(where(is.numeric), ~ format(., scientific = FALSE)))
        
        df.promoter.dist.result %>% head()
        df.promoter.dist.result %>% dim() # any: 264932(dedups)// 265565      9(1 distance) dups
        
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
          dplyr::select(loop.id, promoter_id, value, loop.res)
        relative.pos.df.promoter.dist.result %>% head()
        relative.pos.df.promoter.dist.result %>% 
          dim() # any: 264932(dedups) // 265565      5 (dups)
        relative.pos.df.promoter.dist.result
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
              labs(title = paste0("Density of Promoter Found over Loop on Chr", chr),
                   x = "Relative Position to Loop",
                   y = "Density"
              ) +
              theme(plot.title = element_text(hjust = 0.5)) 
            
            plot1.promoter.hist <- relative.pos.df.promoter.dist.result.chr %>% 
              ggplot(aes(x = value)) +
              geom_histogram(fill = "skyblue", color = "black", alpha = 0.7, bins=200) +
              labs(title = paste0("Histogram of Promoter Found over Loop on Chr", chr),
                   x = "Relative Position to Loop",
                   y = "Count"
              ) +
              theme(plot.title = element_text(hjust = 0.5)) 
            
            
            combined_plot_promoter <- plot1.promoter.hist + plot1.promoter.dens # + can be used instead of |
            plot_promoter_list[[chr]] <- combined_plot_promoter
            
            message("END: Processing chromosome: ", chr)
          }, error = function(e) {
            
            message("Error processing chromosome: ", chr)
            message("Error message: ", e$message)
          })
        }
        
        pdf("./figures/submission/overall_distribution_of_promoter_by_chromosome.pdf", width = 11*0.8, height = 8.5*0.8)
        
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
          # geom_histogram(fill = "skyblue", color = "black", alpha = 0.7, bins=200) +
          geom_histogram(fill = "skyblue", color = "grey70", alpha = 0.7, bins = 200, linewidth = 0.1) +
          labs(# title = "Histogram of Promoter Found over Loop",
            x = "Relative Position to Loop",
            y = "Count"
          ) +
          theme(plot.title = element_text(hjust = 0.5))
        
        plot.promoter.hist
        
        # Promoter Density Plot (by resolution combined, right)
        plot.promoter.dens <- ggplot(relative.pos.df.promoter.dist.result, aes(x = value, color = loop.res, fill = loop.res)) +
          geom_density(alpha = 0.3) +  # transparent
          ylim(c(0, 0.75)) +
          scale_color_manual(values = c("5K" = "#a6cee3", "10K" = "#1f78b4", "25K" = "#1f3a93")) + 
          scale_fill_manual(values = c("5K" = "#a6cee3", "10K" = "#1f78b4", "25K" = "#1f3a93")) +
          labs(# title = "Density of Promoter Found over Loop",
            x = "Relative Position to Loop",
            y = "Density",
            color = "Resolution",
            fill = "Resolution") + 
          theme(plot.title = element_text(hjust = 0.5))
        
        plot.promoter.dens
        
        pdf("./figures/submission/overall_distribution_of_promoter_by_resolution.pdf", width = 11*0.8, height = 8.5*0.4)
        
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
        df.DISTINCT.loop.deep.sample.all %>% head()
        
        ############################################
        # padding case draft: L 1/2 outer, L 1/4 inner
        ############################################
        df.DISTINCT.loop.deep.sample.all.padded.for.promoter <- df.DISTINCT.loop.deep.sample.all %>% 
          left_join(chrom_ends2, by = "chr2") %>%
          mutate(x_mid = (x1 + x2)/2, # middle point of each end
                 y_mid = (y1 + y2)/2, # middle point of each end
                 distance = y_mid - x_mid,
                 x0 = ifelse(x_mid - (distance / 2) < 0, 0, x_mid - (distance / 2)), # 1/2 distance for OUTER PADDING
                 x3 = (x_mid + (distance / 4)), # 1/4 distance for INNER PADDING
                 y0_raw = y_mid - (distance * 0.25),
                 y0 = ifelse((y_mid - (distance/4)) < 0, 0, y_mid - (distance/4)), # 1/4 distance for INNER PADDING
                 y3_raw = y_mid + (distance * 0.5),
                 y3 = y_mid + (distance/2), # 1/2 distance for OUTER PADDING
                 new.loop.id = paste0(chr1, '_', x0, '_', x3, '_', chr2, '_', y0, '_', y3, '_', end.distance),
                 delta_x  = x3 - x0,
                 delta_y  = y3 - y0,
                 inner.distance = y0 - x3)
        
        df.DISTINCT.loop.deep.sample.all.padded.for.promoter %>% head()
        
        
        print(count_each_filter(df.DISTINCT.loop.deep.sample.all.padded.for.promoter))
        
        chrom_ends2
        df.DISTINCT.loop.deep.sample.all
        ############################################
        # padding by resolution
        ############################################
        df.DISTINCT.loop.deep.sample.all.padded.1trial.for.promoter <- df.DISTINCT.loop.deep.sample.all %>% 
          left_join(chrom_ends2, by = "chr2") %>%
          mutate(x_mid = (x1 + x2)/2, # middle point of each end
                 y_mid = (y1 + y2)/2, # middle point of each end
                 e     = end.distance,
                 pad_scale = case_when(
                   resolution == "25K" ~ 0.5,   # 원래 루프 (중점 ± 0.5e)
                   resolution == "10K" ~ 0.75,  # half resolution 추가
                   resolution == "5K"  ~ 1.5,   # 1x resolution 추가
                   TRUE ~ 0.5                   # default (25K처럼)
                 ),
                 x0 = pmax(0, x_mid - pad_scale * e),
                 x3 = x_mid + pad_scale * e,
                 y0_raw = y_mid - pad_scale * e,
                 y0 = pmax(0, y0_raw),
                 y3_raw = y_mid + pad_scale * e,
                 y3 = pmin(y3_raw, chr2_end),
                 new.loop.id = paste0(chr1, '_', x0, '_', x3, '_', chr2, '_', y0, '_', y3, '_', end.distance),
                 delta_x  = x3 - x0,
                 delta_y  = y3 - y0,
                 inner.distance = y0 - x3)
        
        result_counts_1trial_promoter <- count_each_filter(df.DISTINCT.loop.deep.sample.all.padded.1trial.for.promoter)
        print(result_counts_1trial_promoter)
        
        # quantile
        ylim_vals_1trial_promoter <- df.DISTINCT.loop.deep.sample.all.padded.1trial.for.promoter %>%
          pull(inner.distance) %>%
          quantile(c(0.05, 0.95), na.rm = TRUE)
        ylim_vals_1trial_promoter
        
        # min & median 
        stats_1trial_promoter <- df.DISTINCT.loop.deep.sample.all.padded.1trial.for.promoter %>%
          summarise(
            min_val = min(inner.distance, na.rm = TRUE),
            median_val = median(inner.distance, na.rm = TRUE)
          )
        stats_1trial_promoter
        
        df.DISTINCT.loop.deep.sample.all.padded.1trial.for.promoter %>%
          ggplot(aes(y = inner.distance)) +
          geom_boxplot(outlier.shape = NA, fill = "#4C9F70", color = "black") +
          geom_hline(yintercept = stats_1trial_promoter$min_val, linetype = "dashed", color = "blue", linewidth = 0.6) +
          geom_hline(yintercept = stats_1trial_promoter$median_val, linetype = "dotted", color = "red", linewidth = 0.8) +
          annotate("text", x = 0.2, y = stats_1trial_promoter$min_val, label = paste0("Min: ", stats_1trial_promoter$min_val), color = "blue", hjust = 0) +
          annotate("text", x = 0.2, y = stats_1trial_promoter$median_val, label = paste0("Median: ", stats_1trial_promoter$median_val), color = "red", hjust = 0) +
          coord_cartesian(ylim = ylim_vals_1trial_promoter) +
          labs(
            title = "Inner Distance (y0 - x3) — Min & Median Annotated",
            y = "Inner Distance (bp)",
            x = ""
          ) +
          theme_minimal(base_size = 14)
        
        df.DISTINCT.loop.deep.sample.all.padded.for.promoter <- df.DISTINCT.loop.deep.sample.all.padded.1trial.for.promoter
        
        # 1. Promoter + UPSTREAM
        df.DISTINCT.loop.deep.sample.all.padded.for.promoter.UP.GR<- GRanges(seqnames=df.DISTINCT.loop.deep.sample.all.padded.for.promoter$chr1, 
                                                                             ranges=IRanges(start=df.DISTINCT.loop.deep.sample.all.padded.for.promoter$x0, 
                                                                                            end=df.DISTINCT.loop.deep.sample.all.padded.for.promoter$x3), 
                                                                             loop.id=df.DISTINCT.loop.deep.sample.all.padded.for.promoter$loop.id, 
                                                                             end.distance=df.DISTINCT.loop.deep.sample.all.padded.for.promoter$end.distance,
                                                                             resolution=df.DISTINCT.loop.deep.sample.all.padded.for.promoter$resolution,
                                                                             new.loop.id=df.DISTINCT.loop.deep.sample.all.padded.for.promoter$new.loop.id)
        
        index.distinct.promoter.w.up.loop.each.end <- findOverlaps(
          df.promoter.rn7.GR, 
          df.DISTINCT.loop.deep.sample.all.padded.for.promoter.UP.GR, 
          type = "any",
          select = "all"
        )
        index.distinct.promoter.w.up.loop.each.end # any: 68983 // within : 68958
        
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
        
        df.overlapping.promoter.w.UPSTREAM.result.each.end %>% dim() # 68983     6 (dedups) // dups: 69110/95831(promoter)
        df.overlapping.promoter.w.UPSTREAM.result.each.end %>% head()
        
        df.overlapping.promoter.w.UPSTREAM.result.each.end %>% 
          count(end.up.distance)
        # any
        # dedups
        # end.up.distance     n
        # 1            5000 15392
        # 2           10000 20602
        # 3           25000 32989
        
        # 2. Promoter + DOWNSTREAM
        df.DISTINCT.loop.deep.sample.all.padded.for.promoter.DOWN.GR<- GRanges(seqnames=df.DISTINCT.loop.deep.sample.all.padded.for.promoter$chr1, 
                                                                               ranges=IRanges(start=df.DISTINCT.loop.deep.sample.all.padded.for.promoter$y0, 
                                                                                              end=df.DISTINCT.loop.deep.sample.all.padded.for.promoter$y3), 
                                                                               loop.id=df.DISTINCT.loop.deep.sample.all.padded.for.promoter$loop.id, 
                                                                               end.distance=df.DISTINCT.loop.deep.sample.all.padded.for.promoter$end.distance,
                                                                               resolution=df.DISTINCT.loop.deep.sample.all.padded.for.promoter$resolution,
                                                                               new.loop.id=df.DISTINCT.loop.deep.sample.all.padded.for.promoter$new.loop.id)
        
        # 1. promoter + DOWNSTREAM (df.DISTINCT.loop.deep.sample.all.padded.for.promoter.DOWN.GR, df.promoter.rn7.GR)
        index.distinct.promoter.w.down.loop.each.end <- findOverlaps(
          df.promoter.rn7.GR, 
          df.DISTINCT.loop.deep.sample.all.padded.for.promoter.DOWN.GR, 
          type = "any",
          select = "all"
        )
        index.distinct.promoter.w.down.loop.each.end # any: 79080 // within: 79050
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
        
        df.overlapping.promoter.w.DOWNSTREAM.result.each.end %>% dim() # any: 79080     6 (dedups) // 79241     6 (dups)
        df.overlapping.promoter.w.DOWNSTREAM.result.each.end %>% head()
        
        df.overlapping.promoter.w.DOWNSTREAM.result.each.end %>% 
          count(end.down.distance)
        # any
        # dedups
        # end.down.distance     n
        # 1              5000 18483
        # 2             10000 23025
        # 3             25000 37572
        # total: 79080 = 18483 + 23025 + 37572 
        # dups
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
        
        df.overlapping.promoter.w.BOTH.result %>% dim() # any: 148063 (dedups) // 148351 (dups)
        df.overlapping.promoter.w.BOTH.result %>% head()
        df.overlapping.promoter.w.BOTH.result %>% count(resolution)
        
        # any
        # dedups
        # resolution     n
        # 1         5K 33875
        # 2        10K 43627
        # 3        25K 70561
        
        # for Q1
        df.promoter.counts <- df.overlapping.promoter.w.BOTH.result %>%
          group_by(loop.id, WHERE, resolution) %>%
          summarise(promoter_count = n_distinct(promoter.id), .groups = 'drop')
        
        df.promoter.counts
        df.promoter.counts %>% 
          filter(promoter_count == 0)
        
        ###############
        ###############
        df.promoter.case <- df.promoter.counts %>%
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
        
        df.promoter.case %>% count(case)
        # A tibble: 3 × 2
        # case      n
        # <chr> <int>
        # 1 BOTH  13018
        # 2 DOWN   5860
        # 3 UP     6278
        
        df.loop.with.promoter.case <- df.DISTINCT.loop.deep.sample.all %>% # 31773
          left_join(df.promoter.case, by = "loop.id") %>%
          mutate(case = ifelse(is.na(case), "NONE", case))
        
        df.loop.with.promoter.case
        df.loop.with.promoter.case %>% count(case)
        # case     n
        # 1 BOTH 13018
        # 2 DOWN  5860
        # 3 NONE  6617 : 21%
        # 4   UP  6278
        
        df.DISTINCT.loop.deep.sample.all %>% dim() # [1] 31773    11
        
        # proofreading
        # check 1. DOWN
        down_promoter_ids <- df.promoter.counts %>%
          group_by(loop.id) %>%
          filter(all(WHERE == "DOWN")) %>%
          distinct(loop.id)
        down_promoter_ids # 5860
        
        # check 2. UP
        up_promoter_ids <- df.promoter.counts %>%
          group_by(loop.id) %>%
          filter(all(WHERE == "UP")) %>%
          distinct(loop.id)
        
        up_promoter_ids # 6278
        
        # 3. UP & DOWN
        both_promoter_ids <- df.promoter.counts %>%
          group_by(loop.id) %>%
          filter(all(c("UP", "DOWN") %in% WHERE) & n_distinct(WHERE) == 2) %>%
          distinct(loop.id)
        
        both_promoter_ids # 13018
        
        # 3 cases: 13018 + 5860 + 6278 = 25156 + 6617 (no promoter on both ends) = TOTAL 31773
        
        df.promoter.counts %>% 
          count(resolution)
        ###############
        ###############
        # panjun, checkout this section. 
        # threshould for ctcf is about 2^2.5, regardless resolution
        p.hist.promoter.loopend.count<-ggplot(df.promoter.counts, aes(x=log2(promoter_count) + 1))+
          geom_histogram()+
          facet_wrap(~WHERE+resolution, scales="free_y")
        
        pdf(file="histogram_number_of_promoter_within_ends_of_loops_by_resolution_1trial.pdf", p.hist.promoter.loopend.count, width=10, height=8)
        p.hist.promoter.loopend.count
        
        dev.off()
        
        
        promoter.stats.by.resolution <- df.promoter.counts %>%
          # group_by(resolution) %>%
          group_by(WHERE, resolution) %>%
          summarise(
            Min = min(promoter_count),
            Q1 = quantile(promoter_count, 0.25, na.rm = TRUE),
            Median = median(promoter_count, na.rm = TRUE),
            Q3 = quantile(promoter_count, 0.75, na.rm = TRUE),
            Mean = mean(promoter_count, na.rm = TRUE),
            SD = sd(promoter_count, na.rm = TRUE),
            Max = max(promoter_count),
            .groups = 'drop'
          )
        
        
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
        
        # TODO inquiries on Min value
        promoter.stats.by.resolution.each.end
        
        # any
        # dedups
        # WHERE resolution   Min    Q1 Median    Q3  Mean    SD   Max
        # <fct> <fct>      <int> <dbl>  <dbl> <dbl> <dbl> <dbl> <int>
        # 1 UP    5K             1     1      2     3  3.97 20.1    589
        # 2 UP    10K            1     1      2     3  2.86  9.07   259
        # 3 UP    25K            1     1      2     3  4.02 15.6    589
        # 4 DOWN  5K             1     1      1     2  4.95 30.6    743
        # 5 DOWN  10K            1     1      2     3  3.27 17.0    688
        # 6 DOWN  25K            1     1      2     3  4.64 23.4    743
        
        # dups
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
                       position = position_dodge(width = 0.75), vjust = -0.5, color = "black") +
          stat_summary(fun.data = function(y) {
            data.frame(
              y = quantile(y, 0.25),
              label = round(quantile(y, 0.25), 1)
            )
          }, geom = "text", aes(label = after_stat(label)), position = position_dodge(width = 0.75), vjust = 2.5, color = "blue") +
          stat_summary(fun.data = function(y) {
            data.frame(
              y = quantile(y, 0.75),
              label = round(quantile(y, 0.75), 1)
            )
          }, geom = "text", aes(label = after_stat(label)), position = position_dodge(width = 0.75), vjust = -2.5, color = "blue") +
          stat_summary(fun.data = function(y) {
            data.frame(
              y = min(y),
              label = round(min(y), 1)
            )
          }, geom = "text", aes(label = after_stat(label)), position = position_dodge(width = 0.75), vjust = 1.5, hjust = -0.5, color = "red") +
          stat_summary(fun.data = function(y) {
            data.frame(
              y = max(y),
              label = round(max(y), 1)
            )
          }, geom = "text", aes(label = after_stat(label)), position = position_dodge(width = 0.75), vjust = -1.5, hjust = -0.5, color = "red") +
          theme_minimal() + 
          theme(plot.title = element_text(hjust = 0.5)) +
          scale_fill_manual(values = c("5K" = "#a6cee3", "10K" = "#1f78b4", "25K" = "#1f3a93"), breaks = c("5K", "10K", "25K")) +
          scale_y_continuous(limits = c(0, 2 * q3_promoter_count_each_end), 
                             breaks = seq(0, 2 * q3_promoter_count_each_end, by = 1)) +
          labs(# title = "Boxplot for Number of Promoter by Chromosome and Resolution", 
            x = "Chromosome", y = "Number of Promoters per Loop")
        
        boxplot.w.promoter.by.chr.and.res.each.end
        
        pdf("./figures/submission/end_boxplot_num_promoter_by_chr_res.pdf", width = 11*0.8, height = 8.5*0.8)
        grid.arrange(boxplot.w.promoter.by.chr.and.res.each.end, ncol = 1)
        dev.off()
        
        boxplot.w.promoter.by.res <- df.overlapping.promoter.w.BOTH.result.boxplot %>% 
          ggplot(aes(x = resolution, y = promoter_count_by_loop_id_each_end, fill = resolution)) +
          geom_boxplot() +
          stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 1)), 
                       position = position_dodge(width = 0.75), vjust = -0.5) + 
          stat_summary(fun.data = function(y) {
            data.frame(y = quantile(y, 0.25), label = paste0("Q1: ", round(quantile(y, 0.25), 1)))
          }, geom = "text", aes(label = after_stat(label)), 
          position = position_dodge(width = 0.75), vjust = 6.5) +
          stat_summary(fun.data = function(y) {
            data.frame(y = quantile(y, 0.75), label = paste0("Q3: ", round(quantile(y, 0.75), 1)))
          }, geom = "text", aes(label = after_stat(label)), 
          position = position_dodge(width = 0.75), vjust = -6.5, color = "blue") +
          stat_summary(fun.data = function(y) {
            data.frame(y = min(y), label = paste0("Min: ", round(min(y), 1)))
          }, geom = "text", aes(label = after_stat(label)), 
          position = position_dodge(width = 0.75), vjust = 1.5, color = "red") +
          stat_summary(fun.data = function(y) {
            data.frame(y = max(y), label = paste0("Max: ", round(max(y), 1)))
          }, geom = "text", aes(label = after_stat(label)), 
          position = position_dodge(width = 0.75), vjust = -1.5, color = "red") +
          theme_minimal() + 
          theme(plot.title = element_text(hjust = 0.5)) +
          scale_fill_manual(values = c("5K" = "#a6cee3", "10K" = "#1f78b4", "25K" = "#1f3a93"), 
                            breaks = c("5K", "10K", "25K")) +
          scale_y_continuous(limits = c(0, 2 * q3_promoter_count_each_end), 
                             breaks = seq(0, 2 * q3_promoter_count_each_end, by = 10)) +
          labs(# title = "Boxplot for Number of Promote by Resolution", 
            x = "Resolution", y = "Number of Promoters per Loop")
        
        boxplot.w.promoter.by.res
        
        pdf("./figures/submission/end_boxplot_num_promoter_by_res.pdf", width = 11*0.8, height = 8.5*0.8)
        grid.arrange(boxplot.w.promoter.by.res, ncol = 1)
        dev.off()
        
        # figures above in a PDF
        pdf("./figures/submission/boxplot_of_promoter_count.pdf", width = 11*0.8, height = 8.5*0.8)
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
                       position = position_dodge(width = 0.75), vjust = -0.5) +
          stat_summary(fun.data = function(y) {
            data.frame(y = median(y), label = paste0("Q1: ", round(quantile(y, 0.25), 1), "\nQ3: ", round(quantile(y, 0.75), 1)))
          }, geom = "text", aes(label = after_stat(label)), position = position_dodge(width = 0.75), vjust = 1.5) +  
          theme_minimal() + 
          theme(
            plot.title = element_text(hjust = 0.5),  # title center
            legend.position = "none"
          ) +
          scale_fill_manual(values = c("UP" = "#b2df8a", "DOWN" = "#33a02c")) +
          scale_y_continuous(limits = c(0, 2 * q3_promoter_count_each_end),  
                             breaks = seq(0, 2 * q3_promoter_count_each_end, by = 5)) +
          labs(# title = "Boxplot for Promoter by Up/Downstream End", 
            x = "Upstream and Downstream Ends", y = "Number of Promoters")
        
        boxplot.w.promoter.ALL.chr
        
        # figure by end & res
        boxplot.w.promoter.ALL.chr.by.res <- ggplot(df.overlapping.promoter.w.BOTH.result.boxplot, aes(x = WHERE, y = promoter_count_by_loop_id_each_end, fill = resolution)) +
          geom_boxplot(position = position_dodge(width = 0.75)) +
          stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 1)), 
                       position = position_dodge(width = 0.75), vjust = -0.5) +
          stat_summary(fun.data = function(y) {
            data.frame(y = median(y), label = paste0("Q1: ", round(quantile(y, 0.25), 1), "\nQ3: ", round(quantile(y, 0.75), 1)))
          }, geom = "text", aes(label = after_stat(label)), position = position_dodge(width = 0.75), vjust = 1.5) +  
          stat_summary(fun.data = function(y) {
            data.frame(y = min(y), label = paste0("Min: ", round(min(y), 1)))
          }, geom = "text", aes(label = after_stat(label)), position = position_dodge(width = 0.75), vjust = 1.5, hjust = -0.2, color = "red") +
          theme_minimal() + 
          theme(
            plot.title = element_text(hjust = 0.5),  # title center
            legend.position = "none"
          ) +
          scale_fill_manual(values = c("5K" = "#a6cee3", "10K" = "#1f78b4", "25K" = "#1f3a93"), breaks = c("5K", "10K", "25K")) +
          scale_y_continuous(limits = c(0, 2 * q3_promoter_count_each_end), 
                             breaks = seq(0, 2 * q3_promoter_count_each_end, by = 5)) +
          labs(# title = "Boxplot for Promoter in Up/Downstream End by Resolution", 
            x = "Upstream and Downstream Ends", y = "Number of Promoters", fill = "Resolution")
        
        boxplot.w.promoter.ALL.chr.by.res
        
        pdf("./figures/submission/boxplot_for_no.promoter_by_end_res.pdf", width = 8.5*0.8, height = 11*0.8)
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
               y = "Number of Loop") +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5)) +
          stat_count(aes(label = paste0(..count.., " (", round(..count../total_distinct_loops_promoter * 100, 1), "%)")), 
                     geom = "text", 
                     vjust = -0.5) +
          annotate("text", x = 5, y = max(table(promoter_count_histogram_data$promoter_count_total)) * 1.1, 
                   label = paste("Total Loops:", total_distinct_loops_promoter), hjust = 1) +
          annotate("text", x = 5, y = max(table(promoter_count_histogram_data$promoter_count_total)) * 1.05, 
                   label = paste("Loops > 5 Promoter: ", more_than_5_loops, " (", more_than_5_loops_percent, "%)"), 
                   hjust = 1, color = "red")
        
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
               y = "Number of Loop") +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5)) +
          geom_text(aes(label = paste0(one_sided_promoter_loop_count, "/", loop_count, " (", round(ratio * 100, 1), "%)")),
                    vjust = -0.5, color = "black")
        
        one.sided.loop.histogram
        one.sided.loop.histogram.data.promoter.final %>% view()
        
        pdf("./figures/0909/promoter_count_per_loop_histogram.pdf", width = 16.5, height = 23.5)
        
        grid.arrange(promoter.count.per.loop.histogram, one.sided.loop.histogram, ncol = 1)
        dev.off()
        
        ##########################################################
        ##########################################################
        ##########################################################
        # ideogram (HAO)
        ##########################################################
        ##########################################################
        ##########################################################
        library(RIdeogram)
        library(rsvg)
        library(magick)
        library(scales)
        
        getwd()
        
        chromosome_data_ideo_old <- read.table(file="~/dropbox/K\ P/Gateway_to_Hao/enhancer/data/rn7_chromosome_length.tsv", sep="\t", head=T) %>% 
          mutate(start = as.numeric(start), end = str_remove_all(end, ",") %>% as.numeric()) %>%
          mutate(start = start - 1) %>%
          dplyr::rename(Chr = chr, Start = start, End = end) %>%
          mutate(CE_start = NA, CE_end = NA) %>% 
          view()
        
        chromosome_data_ideo_old
        
        chromosome_data <- read.table(file="~/dropbox/K\ P/Gateway_to_Hao/enhancer/data/rn7_chromosome_length_from_ucsc.tsv", sep="\t") %>%
          dplyr::rename(Chr = V1, End = V2) %>% 
          mutate(Chr = str_remove(Chr, '^chr')) %>% 
          mutate(Start = 0) %>% 
          mutate(Start = as.numeric(Start), End = as.numeric(End)) %>% 
          mutate(CE_start = NA, CE_end = NA)
        
        chromosome_data
        chromosome_data %>% head()
        
        # CTCF for ideogram: df_ctcf
        df.DISTINCT.fimo.2nd.trial.ctcf %>% head()
        df.DISTINCT.fimo.2nd.trial.ctcf # 3191859 (sub .4)
        
        df_ctcf_ideogram <- df.DISTINCT.fimo.2nd.trial.ctcf %>% 
          dplyr::select(chr, start, end) %>% 
          mutate(chr = str_remove(chr, "chr"))
        df_ctcf_ideogram
        
        
        df_binned <- df %>%
          mutate(
            StartBin = floor(start / bin_size) * bin_size + 1,
            EndBin = StartBin + bin_size - 1
          )
        
        # Step 4: 각 bin마다 개수 세기
        gene_density <- df_binned %>%
          group_by(Chr = chr, Start = StartBin, End = EndBin) %>%
          summarise(Value = n(), .groups = "drop") %>%
          arrange(Chr, Start)
        
        # 결과 확인
        head(gene_density)
        
        # TSS for ideogram: df_tss_ucsc
        df.tss.ucsc %>% head()
        df_tss_ucsc_ideogram <- df.tss.ucsc %>% 
          dplyr::select(chr, start, end) %>% 
          mutate(chr = str_remove(chr, "chr"))
        df_tss_ucsc_ideogram %>% head()
        
        # promoter for ideogram: df_promoter
        df.promoter.rn7 %>% head()
        df_promoter_ideogram <- df.promoter.rn7 %>% 
          dplyr::select(chr, start, end) %>% 
          mutate(chr = str_remove(chr, "chr"))
        df_promoter_ideogram %>% head()
        
        bin_size <- 1000000
        
        process_feature_bins <- function(feature_df, chromosome_ends, bin_size = 1e6, label = "Feature") {
          options(scipen = 999)
          
          # col name
          colnames(chromosome_ends) <- c("Chr", "Start", "End", "CE_start", "CE_end")
          colnames(feature_df) <- c("Chr", "Start", "End")
          
          # col type
          chromosome_ends$Chr <- as.character(chromosome_ends$Chr)
          feature_df$Chr <- as.character(feature_df$Chr)
          
          result <- list()
          
          valid_chromosomes <- as.character(c(1:20, "X", "Y"))
          feature_df <- feature_df %>% filter(Chr %in% valid_chromosomes)
          
          for (current_chr in unique(feature_df$Chr)) {
            chr_data <- feature_df %>% filter(Chr == current_chr)
            chr_end <- chromosome_ends %>% filter(Chr == current_chr) %>% pull(End)
            
            if (length(chr_end) != 1 || !is.finite(chr_end) || chr_end < bin_size) {
              warning(paste("No valid end position for", current_chr, "- skipping"))
              next
            }
            
            # interval: Start = End+1
            bin_edges <- seq(0, chr_end, by = bin_size)
            if (length(bin_edges) < 2) {
              warning(paste("Not enough bin edges for", current_chr, "- skipping"))
              next
            }
            bins <- data.frame(
              Start = bin_edges[-length(bin_edges)],
              End = bin_edges[-1]
            )
            
            if (tail(bins$End, 1) < chr_end) {
              bins <- rbind(bins, data.frame(Start = bins$End[nrow(bins)] + 1, End = chr_end))
            }
            
            bins <- bins %>%
              mutate(Chr = current_chr,
                     Value = 0,
                     Feature = label)
            
            for (i in 1:nrow(bins)) {
              bins$Value[i] <- sum(
                chr_data$Start >= bins$Start[i] & chr_data$Start <= bins$End[i],
                na.rm = TRUE
              )
            }
            
            result[[current_chr]] <- bins
          }
          
          final <- bind_rows(result)
          final <- final %>%
            dplyr::select(Chr, Start, End, Value, Feature)
          
          return(final)
        }
        
        # HAO SLIDE!! [ # of ]
        process_sliding_feature_bins <- function(feature_df, chromosome_ends, bin_size = 2e6, step_size = 1e6, label = "Feature") {
          options(scipen = 999)
          result <- list()
          
          valid_chromosomes <- as.character(c(1:20, "X", "Y", "x", "y"))
          feature_df <- feature_df %>% filter(chr %in% valid_chromosomes)
          
          for (current_chr in unique(feature_df$chr)) {
            chr_data <- feature_df %>% filter(chr == current_chr)
            chr_end <- chromosome_ends %>% filter(chr == current_chr) %>% pull(end)
            
            if (length(chr_end) != 1) {
              stop(paste("Invalid chr_end length for", current_chr))
            }
            
            bin_starts <- seq(1, chr_end - bin_size + 1, by = step_size)
            bin_ends <- bin_starts + bin_size - 1
            bin_ends[bin_ends > chr_end] <- chr_end
            
            bins <- data.frame(
              chr = current_chr,
              start = bin_starts,
              end = bin_ends,
              value = 0,
              feature = label
            )
            
            for (i in 1:nrow(bins)) {
              bins$value[i] <- sum(
                chr_data$start >= bins$start[i] & chr_data$start <= bins$end[i]
              )
            }
            
            result[[current_chr]] <- bins
          }
          
          final <- bind_rows(result)
          final <- final %>%
            dplyr::select(chr, start, end, value, feature)
          
          return(final)
        }
        
        ctcf_density <- process_feature_bins(df_ctcf_ideogram, chromosome_data, bin_size, label = "CTCF") %>% dplyr::select(-last_col()) # color = "#E41A1C"
        tss_density  <- process_feature_bins(df_tss_ucsc_ideogram, chromosome_data, bin_size, label = "TSS") %>% dplyr::select(-last_col()) # color = "#4DAF4A"
        promoter_density <- process_feature_bins(df_promoter_ideogram, chromosome_data, bin_size, label = "Promoter") %>% dplyr::select(-last_col())# color = "#377EB8"
        
        chromosome_data
        
        ctcf_density_norm <- ctcf_density %>% 
          mutate(Normalized = (Value - min(Value)) / (max(Value) - min(Value))) %>% 
          dplyr::select(Chr, Start, End, Value = Normalized)
        
        ctcf_density_norm
        
        ctcf_density_norm_chr <- ctcf_density %>%
          group_by(Chr) %>%
          mutate(Value = (Value - min(Value)) / (max(Value) - min(Value))) %>%
          ungroup()
        
        ctcf_density_norm_chr
        
        ctcf_density_norm_log <- ctcf_density %>% 
          mutate(Normalized = log1p(Value)) %>%  # log(1 + x)
          dplyr::select(Chr, Start, End, Value = Normalized)
        
        ctcf_density_norm_log
        
        
        ideogram(
          karyotype = chromosome_data,
          # overlaid = ctcf_density,
          overlaid = ctcf_density_norm_log,
          # overlaid = ctcf_density_norm,
          # overlaid = ctcf_density_norm_chr,
          label = NULL,  # No additional labels for now
          label_type = "heatmap"
        )
        
        # svg_file <- "chromosome.svg"
        # convertSVG(svg_file, device = "png")
        
        # Convert the SVG output to PNG (or other format)
        convertSVG("chromosome.svg", device = "png")
        
        
        
        ctcf_density_s <- process_sliding_feature_bins(feature_df = df_ctcf_ideogram, chromosome_ends = chromosome_data, bin_size = 2e6, step_size = 1e6, label = "CTCF") %>% dplyr::select(-last_col())
        tss_density_s <- process_sliding_feature_bins(feature_df = df_tss_ucsc_ideogram, chromosome_ends = chromosome_data, bin_size = 2e6, step_size = 1e6, label = "TSS") %>% dplyr::select(-last_col())
        promoter_density_s <- process_sliding_feature_bins(feature_df = df_promoter_ideogram, chromosome_ends = chromosome_data, bin_size = 2e6, step_size = 1e6, label = "Promoter") %>% dplyr::select(-last_col())
        
        chromosome_data <- chromosome_data %>% mutate(end = as.numeric(end)) %>% dplyr::rename(Chr = chr, Start = start, End = end)
        chromosome_data
        
        ctcf_density_s
        tss_density_s
        promoter_density_s
        
        # fill_color = c("#f8766d", "#32ba36", "#629bfe"),
        ctcf_density_format_s <- ctcf_density_format_s <- ctcf_density_s %>% dplyr::rename(Chr = chr, Start = start, End = end, Value = value) # %>% mutate(Value = Value/(End - Start))
        tss_promoter_density_s <- bind_cols(tss_density_s %>% dplyr::rename(Value_1 = value) %>% mutate(Color_1 = "32ba36"), promoter_density_s %>% dplyr::select(value) %>% dplyr::rename(Value_2 = value) %>% mutate(Color_2 = "629bfe")) %>% 
          dplyr::rename(Chr = chr, Start = start, End = end) #%>% mutate(Value_1 = Value_1 / (End - Start), Value_2 / (End - Start))
        
        ctcf_density_format_s %>% head()
        tss_promoter_density_s %>% head()
        
        TSS "#32ba36"
        Promoter "#629bfe"
        
        # Pi_for_CE_and_CW
        # Chr     Start       End    Value_1 Color_1    Value_2 Color_2
        # 1     1         1   2000000 0.00273566  fc8d62 0.00385702  8da0cb
        # 2     1   1000001   3000000 0.00239580  fc8d62 0.00331109  8da0cb
        # 3     1   2000001   4000000 0.00319407  fc8d62 0.00374530  8da0cb
        ideogram(karyotype = liriodendron_karyotype, 
                 overlaid = Fst_between_CE_and_CW, 
                 label = Pi_for_CE_and_CW, 
                 label_type = "polygon", 
                 colorset1 = c("#e5f5f9", "#99d8c9", "#2ca25f"))
        
        # tss_promoter_density_s
        # > tss_promoter_density_s
        # Chr     Start       End Value_1 Color_1 Value_2 Color_2
        # 1     1         1   2000000     894 #32ba36       1 #629bfe
        # 2     1   1000001   3000000    1573 #32ba36       8 #629bfe
        # 3     1   2000001   4000000     695 #32ba36       7 #629bfe
        
        ideogram(karyotype = chromosome_data, 
                 overlaid = ctcf_density_format_s, 
                 label = tss_promoter_density_s, 
                 label_type = "polygon", 
                 colorset1 = c("#fde0dd", "#fa9fb5", "#E41A1C")
        )
        
        # SVG file save
        svg_file <- "chromosome.svg"
        convertSVG(svg_file, device = "png") # taking time too long
        
        # converting to PNG with size
        png_file <- "chromosome.png"
        rsvg::rsvg_png(svg_file, file = png_file, width = 3300, height = 2550)
        
        image <- image_read(png_file) %>%             # reading PNG
          image_background(color = "white") %>%       # white background
          image_write(png_file)     
        
        # liriodendron_karyotype %>% head()
        # Fst_between_CE_and_CW %>% head()
        # Pi_for_CE_and_CW %>% head()
        
        ###########################################
        ###########################################
        ###########################################
        dat <- ctcf_density %>%
          transmute(chr = as.character(Chr), start = Start, end = End, count = Value) %>%
          mutate(bin_width = end - start + 1,
                 rate_per_mb = ifelse(bin_width > 0, count / (bin_width/1e6), 0))
        
        # 극단값 완화 + 로그
        p99 <- quantile(dat$rate_per_mb, 0.99, na.rm = TRUE)
        rate_w <- pmin(dat$rate_per_mb, p99)
        overlaid_global <- dat %>%
          transmute(chr, start, end, value = log1p(rate_w))  # ideogram은 내부 min–max 매핑
        
        overlaid_z_bychr <- dat %>%
          group_by(chr) %>%
          mutate(z = as.numeric(scale(rate_per_mb))) %>%
          ungroup() %>%
          transmute(chr, start, end, value = z)
        
        # 또는 분위수(ECDF)로 0–1
        overlaid_ecdf_bychr <- dat %>%
          group_by(chr) %>%
          mutate(q = ecdf(rate_per_mb)(rate_per_mb)) %>%
          ungroup() %>%
          transmute(chr, start, end, value = q)
        
        library(stats)
        
        model_df <- dat %>%
          left_join(motif_gc_df, by = c("chr","start","end")) %>%   # motif_gc_df: (chr,start,end,motif_count,gc,map,cutsite_density ...)
          mutate(eff_len_mb = (end - start + 1)/1e6,
                 motif_off  = log(pmax(motif_count, 1)),            # 0 방지
                 offset_val = log(pmax(eff_len_mb, 1e-6)))
        
        fit <- glm(count ~ log1p(motif_count) + poly(gc, 2) + map + cutsite_density + offset(offset_val),
                   family = quasipoisson, data = model_df)
        
        model_df$resid <- residuals(fit, type = "pearson")  # 또는 deviance 잔차
        overlaid_resid <- model_df %>% transmute(chr, start, end, value = resid)
        
        
        
        ###########################################
        ###########################################
        ###########################################
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
          geom_text(aes(label = n), position = position_stack(vjust = 0.5)) +
          scale_fill_manual(values = c("5K" = "#a6cee3", "10K" = "#1f78b4", "25K" = "#1f3a93")) +
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
          geom_text(aes(label = n), position = position_stack(vjust = 0.5)) +
          scale_fill_manual(values = c("5K" = "#a6cee3", "10K" = "#1f78b4", "25K" = "#1f3a93")) +
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
        
        ################################################################################################
        ################################################################################################
        ################################################################################################
        # loops.with.above.q1
        ################################################################################################
        ################################################################################################
        ################################################################################################
        ################################################################################################
        # 4-1. CTCF
        df.ctcf.counts
        ctcf.stats.by.resolution
        
        p.hist.ctcf.loopend.count <- ggplot(df.ctcf.counts, aes(x=log2(ctcf_count)))+
          # p.hist.ctcf.loopend.count <- ggplot(df.ctcf.counts, aes(x=ctcf_count))+
          geom_histogram()+
          facet_wrap(~WHERE+resolution, scales="free_y")
        
        p.hist.ctcf.loopend.count
        
        pdf(file = "./figures/submission/histogram_number_of_ctcf_within_ends_of_loops_by_resolution_end.pdf", width = 11*0.8, height = 8.5*0.8)
        
        print(p.hist.ctcf.loopend.count)
        
        dev.off()
        
        # 2^2.5
        # [1] 5.656854
        # log2(5.656854)
        # [1] 2.5
        
        # loops.with.ctcf.above.q1 <- df.ctcf.counts %>% 
        #   left_join(ctcf.stats.by.resolution %>% 
        #               dplyr::select(resolution, WHERE, Q1), by = c("resolution", "WHERE")) %>% 
        #   # filter(ctcf_count >= Q1)
        #   filter(ctcf_count >= 6)
        
        df.ctcf.counts
        
        df.loops.above.ctcf.threshold <- df.ctcf.counts %>%
          filter(ctcf_count >= 6) %>%
          group_by(loop.id) %>%
          filter(n_distinct(WHERE) == 2) %>%
          ungroup()
        
        df.loops.above.ctcf.threshold
        
        final.loops.from.ctcf.step <- df.loops.above.ctcf.threshold %>% 
          distinct(loop.id) %>% 
          mutate(end.distance = str_split_n(loop.id, '_', 7)) %>% 
          mutate(resolution = case_when(
            end.distance == 5000 ~ "5K",
            end.distance == 10000 ~ "10K",
            end.distance == 25000 ~ "25K",
            TRUE ~ NA
          ))
        
        final.loops.from.ctcf.step
        # %>% 
        #   mutate(loop.id = str_remove(loop.id, "_[^_]+$"))
        final.loops.from.ctcf.step # 18,327 + 10 | sub.4 any, threshold > 6 : 25,620
        final.loops.from.ctcf.step %>% 
          count(resolution)
        # sub.4
        # any: 
        # new threshold
        # resolution     n
        # <chr>      <int>
        # 1 10K         9420
        # 2 25K        12049
        # 3 5K          4151
        
        # 4-2. TSS
        df.overlapping.TSS.w.BOTH.result %>% head()
        df.tss.counts.each.end # 44564 (.4 any dedup)
        
        df.tss.counts.each.end %>%  # 608 + 242 + 1980 = 2790 loops 5 cases ALL have At least more than 2790
          count(tss_count_each_end) %>% filter(tss_count_each_end == 0)
        
        # TSS histogram
        p.hist.tss.loopend.count <- ggplot(df.tss.counts.each.end, aes(x = log2(tss_count_each_end + 1))) +
          geom_histogram() +
          facet_wrap(~WHERE + resolution, scales = "free_y")
        
        p.hist.tss.loopend.count
        
        print(p.hist.tss.loopend.count)
        
        # PDF
        pdf(file = "./figures/submission/histogram_number_of_tss_within_ends_of_loops_by_resolution_end_1trial.pdf", width = 11*0.8, height = 8.5*0.8)
        print(p.hist.tss.loopend.count)
        dev.off()
        
        loops.with.tss.above.q1 <- df.tss.counts.each.end %>% 
          left_join(tss.stats.by.resolution.each.end %>% 
                      dplyr::select(resolution, WHERE, Q1), by = c("resolution", "WHERE")) %>% 
          filter(tss_count_each_end >= Q1)
        
        loops.with.tss.above.q1 # 44564
        
        loops.with.tss.both.ends <- loops.with.tss.above.q1 %>%
          group_by(loop.id) %>%
          filter(any(c("UP", "DOWN") %in% WHERE)) %>% 
          ungroup()
        
        loops.with.tss.both.ends %>% distinct(loop.id) # 27518
        
        final.loops.from.tss.step <- loops.with.tss.both.ends %>% 
          distinct(loop.id) %>% 
          mutate(end.distance = str_split_n(loop.id, '_', 7)) %>% 
          mutate(resolution = case_when(
            end.distance == 5000 ~ "5K",
            end.distance == 10000 ~ "10K",
            end.distance == 25000 ~ "25K",
            TRUE ~ NA
          ))
        
        # 1//0
        # 0//1
        # bigger number of loops qualified with Q1 of TSS
        # bigger number of loops qualified with Q1 of promoter
        
        # %>% 
        #   mutate(loop.id = str_remove(loop.id, "_[^_]+$"))
        final.loops.from.tss.step # any: 27,518 (dedups) w/ either // 11,993
        final.loops.from.tss.step %>% count(resolution)
        # dedups
        # any
        # resolution     n
        # <chr>      <int>
        # 1 10K         6303
        # 2 25K         7515
        # 3 5K          3228
        
        # A tibble: 3 × 2
        # resolution     n
        # <chr>      <int>
        #   1 10K         4279
        # 2 25K         5176
        # 3 5K          2538
        
        # 4-3. promoter
        df.overlapping.promoter.w.BOTH.result
        df.promoter.counts.each.end # 38,174
        promoter.stats.by.resolution.each.end
        
        df.promoter.counts.each.end %>% 
          count(promoter_count_each_end)
        
        # Promoter histogram
        p.hist.promoter.loopend.count <- ggplot(df.promoter.counts.each.end, aes(x = log2(promoter_count_each_end+1))) +
          geom_histogram() +
          facet_wrap(~WHERE + resolution, scales = "free_y")
        
        p.hist.promoter.loopend.count
        print(p.hist.promoter.loopend.count)
        
        # PDF
        pdf(file = "./figures/submission/histogram_number_of_promoters_within_ends_of_loops_by_resolution_end_1trial.pdf", width = 11*0.8, height = 8.5*0.8)
        print(p.hist.promoter.loopend.count)
        dev.off()
        
        
        
        loops.with.promoter.above.q1 <- df.promoter.counts.each.end %>% 
          left_join(promoter.stats.by.resolution.each.end %>% 
                      dplyr::select(resolution, WHERE, Q1), by = c("resolution", "WHERE")) %>% 
          filter(promoter_count_each_end >= Q1)
        
        loops.with.promoter.above.q1
        
        loops.with.promoter.both.ends <- loops.with.promoter.above.q1 %>% # 38,152 + 10
          group_by(loop.id) %>%
          filter(any(c("UP", "DOWN") %in% WHERE)) %>%
          ungroup()
        
        loops.with.promoter.both.ends %>% distinct(loop.id) # 25156
        
        final.loops.from.promoter.step <- loops.with.promoter.both.ends %>% 
          distinct(loop.id) %>% 
          mutate(end.distance = str_split_n(loop.id, '_', 7)) %>% 
          mutate(resolution = case_when(
            end.distance == 5000 ~ "5K",
            end.distance == 10000 ~ "10K",
            end.distance == 25000 ~ "25K",
            TRUE ~ NA
          ))
        
        final.loops.from.promoter.step
        
        # %>% 
        #   mutate(loop.id = str_remove(loop.id, "_[^_]+$"))
        final.loops.from.promoter.step # .4 any (dedups): 25,156// 9,509
        final.loops.from.promoter.step %>% count(resolution)
        # dedups
        # any
        # resolution     n
        # <chr>      <int>
        # 1 10K         9517
        # 2 25K        10511
        # 3 5K          5128
        
        ############################
        ############################
        ############################
        df.tss.counts.each.end
        tss.stats.by.resolution.each.end
        
        tss.loop.type <- df.tss.counts.each.end %>%
          left_join(tss.stats.by.resolution.each.end %>% dplyr::select(WHERE, resolution, Q1), by = c("resolution", "WHERE")) %>%
          filter(tss_count_each_end >= Q1) %>% 
          group_by(loop.id) %>%
          mutate(
            tss.type = if(n() == 2) "both" else as.character(WHERE)
          ) %>%
          ungroup() %>% 
          distinct(loop.id, tss.type)
        
        tss.loop.type %>% count(loop.id, tss.type) # 27518
        
        
        promoter.loop.type <- df.promoter.counts.each.end %>%
          left_join(promoter.stats.by.resolution.each.end %>% dplyr::select(WHERE, resolution, Q1), by = c("resolution", "WHERE")) %>%
          filter(promoter_count_each_end >= Q1) %>%
          group_by(loop.id) %>%
          mutate(
            pro.type = if(n() == 2) "both" else as.character(WHERE)
          ) %>%
          ungroup() %>% 
          distinct(loop.id, pro.type)
        
        promoter.loop.type %>% count(loop.id, pro.type) # 25,156
        
        loop.combined.categories <- inner_join(tss.loop.type, promoter.loop.type, by = "loop.id")
        loop.combined.categories %>% count(tss.type, pro.type)
        
        loop.combined.categories %>% # 24306
          distinct(loop.id)
        
        
        
        ####################################
        ####################################
        ####################################
        # Venn Diagram SUBMISSION
        ####################################
        ####################################
        ####################################
        
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
            fill_color = c("#f8766d", "#629bfe", "#32ba36"),
            show_elements = FALSE  
          ) 
          # + 
          # ggtitle(paste(ctcf_label, "vs Promoter vs TSS")) +
          # theme(plot.title = element_text(hjust = 0.5, size = 22))
          
          venn_plot <- venn_plot +
            annotate("text", x = -1.5, y = 1.5, label = paste("CTCF:", ctcf_count), color = "#f8766d") +
            annotate("text", x = 1.55, y = 1.5, label = paste("Promoter:", promoter_count), color = "#629bfe") +
            annotate("text", x = 0, y = -1.7, label = paste("TSS:", tss_count), color = "#32ba36")
          
          
          return(venn_plot)
        }
        
        venn_plot_submission <- create_venn_plot(final.loops.from.ctcf.step, final.loops.from.promoter.step, final.loops.from.tss.step, "CTCF")
        ggsave(filename = "./figures/submission/ctcf_vs_promoter_tss_venn_diagrams_latest_histo_1trial.pdf", plot = venn_plot_submission, width = 11, height = 8.5, units = "in")
        
        venn_plot_submission
        
        ##################
        # functional loop extraction
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
        
        final.loops.from.ctcf.step # sub.4 any, either 6: 25,620
        final.loops.from.promoter.step # sub.4 dedups any, either: 25,156
        final.loops.from.tss.step # sub.4 dedups any, either: 27,518
        
        overlapping_loops <- extract_overlapping_loops(
          ctcf_data = final.loops.from.ctcf.step,
          promoter_data = final.loops.from.promoter.step,
          tss_data = final.loops.from.tss.step
        )
        
        overlapping_loops$ctcf_promoter_only # PASS// 608       
        overlapping_loops$ctcf_tss_only # PASS// 2519           
        overlapping_loops$ctcf_promoter_tss_overlap # PASS// 20513 
        
        df.ctcf.promoter.only.loop <- data.frame(loop.id = overlapping_loops$ctcf_promoter_only, category = "CP", stringsAsFactors = FALSE)
        df.ctcf.tss.only.loop <- data.frame(loop.id = overlapping_loops$ctcf_tss_only, category = "CT", stringsAsFactors = FALSE)
        df.ctcf.promoter.tss.overlap.loop <- data.frame(loop.id = overlapping_loops$ctcf_promoter_tss_overlap, category = "CPT", stringsAsFactors = FALSE)
        
        # sub.4 dedups: 11891 = 679 + 3096 + 8116
        df.final.loop <- bind_rows(df.ctcf.promoter.only.loop, df.ctcf.tss.only.loop, df.ctcf.promoter.tss.overlap.loop)
        
        df.final.loop # 23640/31773(0.7440279) sub.4 any, either 6// 11526/31773 (0.3627608)
        df.final.loop %>% dim()
        df.final.loop %>% head()
        
        df.final.loop %>% mutate(resolution = str_split_n(loop.id, '_', 7)) %>% count(resolution)
        # resolution     n
        # 1      10000  8675
        # 2      25000 11144
        # 3       5000  3821
        
        save(df.final.loop, file="./figures/submission/df_final_loop_sub.4.any.either.6.rda")
        write.csv(df.final.loop, file = "./figures/submission/df_final_loop_sub.4.dedup.any.either.6.csv", row.names = FALSE)
        
        df.final.loop # 23640
        df.DISTINCT.loop.deep.sample.all
        
        df.final.DISTINCT.loop.joined <- df.final.loop %>% 
          inner_join(df.DISTINCT.loop.deep.sample.all %>% mutate(loop.id = str_remove(loop.id, "_[^_]*$")), by = "loop.id")
        df.final.DISTINCT.loop.joined %>% dim()
        df.final.DISTINCT.loop.joined
        
        
        df.final.loop %>% distinct(loop.id) # 23460
        
        
        ############################
        
        df.final.loop
        
        df.box <- df.final.loop %>% filter(category != "CP") %>% left_join(df.overlapping.TSS.w.BOTH.result %>% mutate(loop.id = str_remove(loop.id, "_[^_]*$")), by = "loop.id") %>% 
          group_by(tss.id) %>%
          summarise(unique_loop_count = n_distinct(loop.id)) %>%
          arrange(desc(unique_loop_count))
        
        stats <- df.box %>%
          summarise(
            median_val = median(unique_loop_count),
            min_val = min(unique_loop_count)
          )
        
        
        stats
        df.box
        
        
        df.overlapping.TSS.w.BOTH.result
        df.overlapping.promoter.w.BOTH.result
        
        ggplot(df.box, aes(x = unique_loop_count)) +
          geom_histogram(binwidth = 1, fill = "skyblue", color = "white") +
          labs(
            title = "Histogram of Unique Loop Counts per TSS",
            x = "Unique Loop Count",
            y = "Frequency"
          ) +
          theme_minimal()
        
        ############################
        
        
        
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
          mutate(resolution = case_when(
            end.distance == 5000 ~ "5K",
            end.distance == 10000 ~ "10K",
            end.distance == 25000 ~ "25K",
            TRUE ~ NA
          ),
          resolution = factor(resolution, levels = c("5K", "10K", "25K")))
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
        pdf(file = "./figures/submission/circos_loops_by_resolution_either.6.pdf", height = 11*0.8, width = 8.5*0.8)
        
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
        
        # Function to plot per chromosome
        plot_circos_for_chromosome <- function(chr) {
          df_filtered <- subset(df.circos.input.log.final.loop, chr1_clean == chr)
          
          if (nrow(df_filtered) > 0) {
            chromosomes_to_display <- unique(c(df_filtered$chr1, df_filtered$chr2))
            
            circos.par(gap.degree = 25)
            circos.initializeWithIdeogram(species = "rn7", chromosome.index = chromosomes_to_display)
            
            for (i in 1:nrow(df_filtered)) {
              circos.genomicLink(
                region1 = df_filtered[i, c("chr1", "midx", "midx")],
                region2 = df_filtered[i, c("chr2", "midy", "midy")],
                col = resolution_colors[as.character(df_filtered$resolution[i])],
                h = df_filtered$height[i],
                border = "black"
              )
            }
            circos.clear()
          }
        }
        
        # Plot per chromosome
        unique_chromosomes <- levels(df.circos.input.log.final.loop$chr1_clean)
        for (chr in unique_chromosomes) {
          plot_circos_for_chromosome(chr)
        }
        
        dev.off()
        
        #########################################
        # Figure for new ones
        #########################################
        
        p1 <- plot.ctcf.dens
        p2 <- plot.tss.dens
        p3 <- plot.promoter.dens
        
        p1 <- plot.ctcf.hist
        p2 <- plot.tss.hist
        p3 <- plot.promoter.hist
        
        base_theme <- theme_bw(base_size = 12) +
          theme(plot.title.position = "plot",
                panel.grid.minor = element_blank())
        
        p1 <- p1 + base_theme
        p2 <- p2 + base_theme
        p3 <- p3 + base_theme
        
        # ylim to c(0,0.75)
        
        
        # 1row 3col legend
        fig_combined <-
          (p1 | p2 | p3) +
          plot_layout(ncol = 3, guides = "collect", widths = c(1, 1, 1)) &
          theme(legend.position = "bottom",
                legend.margin = margin(2, 6, 2, 6),
                plot.tag = element_text(face = "bold", size = 12))
        
        # title/subtitle/pannel tag(a, b, c)
        fig_combined <- fig_combined +
          plot_annotation(
            # title = "Figure X. Your overall title",
            # subtitle = "Optional subtitle or notes",
            tag_levels = "a"
          )
        
        fig_combined
        
        # save (10 x 8.5 inch, 300 dpi)
        ggsave("./figures/submission/density_combined_all.png", fig_combined,
               width = 11, height = 8.5, dpi = 300, bg = "white")
        ggsave("./figures/submission/histogram_combined_all.png", fig_combined,
               width = 11, height = 8.5, dpi = 300, bg = "white")
        
        
        
        
        
        
        
        