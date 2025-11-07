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

getwd()

# Linux
# setwd('C:\\Users\\panju\\Dropbox (UTHSC GGI)\\Gateway_to_Hao\\workshop\\2023_NIH_meeting\\loop_N_tss') # nolint: commented_code_linter.
# setwd('./Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss')
# setwd('~/Desktop/temp/enhancer/dropbox_enhancer_doosan/data/enhancer_atlas2.0/all_species/neuron')
setwd("~/Desktop/temp/enhancer/dropbox_enhancer_doosan")
setwd("/home/pkim/dropbox/Gateway_to_Hao/enhancer/r_files")
source(file.path("/home/pkim/dropbox/Gateway_to_Hao/project_common_code/", "variables.R"))
source(file.path("/home/pkim/dropbox/Gateway_to_Hao/project_common_code/", "funcs.R"))

getwd()

# Mac
setwd("~/dropbox/Gateway_to_Hao/enhancer/r_files")
getwd()
source(file.path("~/dropbox/Gateway_to_Hao/project_common_code/", "variables.R"))
source(file.path("~/dropbox/Gateway_to_Hao/project_common_code/", "funcs.R"))

####################################
# Sequencing stats ###### Figure 1.a
####################################
seq.data <- read.table("../data/library_complexity.tsv", header = TRUE, sep = "\t")
seq.data[, -1] <- lapply(seq.data[, -1], function(x) as.numeric(as.character(x)))
colnames.of.seq.data <- colnames(seq.data) # nolint
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

# func 0-1. saving figures
saving_plot_dual <- function(plot_obj,
                           filename_base,
                           output_dir = "./figures",
                           width_in = 11,
                           height_in = 8.5,
                           scale_x = 0.6,
                           scale_y = 0.6,
                           dpi = 300,
                           bg = "white") {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # path
  pdf_path <- file.path(output_dir, paste0(filename_base, ".pdf"))
  png_path <- file.path(output_dir, paste0(filename_base, ".png"))

  # size calculation
  final_width <- width_in * scale_x
  final_height <- height_in * scale_y

  # PDF
  ggsave(
    filename = pdf_path,
    plot = plot_obj,
    device = "pdf",
    width = final_width,
    height = final_height,
    units = "in"
  )

  # PNG
  ggsave(
    filename = png_path,
    plot = plot_obj,
    device = "png",
    width = final_width,
    height = final_height,
    units = "in",
    dpi = dpi
  )
}
# func 0-2. END

saving_plot_dual(
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
df.init.loop.bed <- init.bedpe.df(loop.file.list, 'loops') %>% 
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
pdf("./figures/submission/lt2mb/common_loops_heatmap_percentage_bw_strains_S3.pdf", width = 8.5, height = 11)
grid.arrange(plots[["5K"]], 
             plots[["10K"]], 
             plots[["25K"]], 
             nrow = 3
)
dev.off()

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

saving_plot_dual(
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

saving_plot_dual(
  plot_obj = network_plot_for_shared_loops,
  filename_base = "network_plot_for_shared_loops_S2",
  output_dir = "./figures/submission/lt2mb")

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
#                    strain num_loop
# 1                    BN-Lx     6535
# 2                     BXH6     6568
# 3                 F344/Stm     2903
# 4                    HXB10     7336
# 5                     HXB2     4656
# 6                    HXB23     7676
# 7                    HXB31     9131
# 8                   LE/Stm     2992
# 9              SHR/OlaIpcv     5263
# 10 SHR/OlaIpcvxBN/NHsdMcwi     5932        

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

saving_plot_dual(
  plot_obj = line_graph_for_loops_per_depth,
  filename_base = "line_graph_for_loops_per_depth_hao_w_new_label_F1b",
  output_dir = "./figures/submission/lt2mb")

figure1_combined <- sequencing_basic_stats + line_graph_for_loops_per_depth +
  plot_layout(ncol = 2, nrow = 1) 

figure1_combined

saving_plot_dual(
  plot_obj = figure1_combined,
  filename_base = "figure1_combined_F1",
  output_dir = "./figures/submission/lt2mb",
  height = 5.5,        ############ NOT 8.5
  scale_x = 1,
  scale_y = 1)

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

saving_plot_dual(
  plot_obj = df_chr_loop_counts_fig,
  filename_base = "loop_counts_per_chr_F2",
  output_dir = "./figures/submission/lt2mb")

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
  mutate(resolution = case_when(
    end.distance == 5000 ~ "5K",
    end.distance == 10000 ~ "10K",
    end.distance == 25000 ~ "25K",
    TRUE ~ NA
  )) %>% 
  mutate(resolution = factor(resolution, levels = c("5K", "10K", "25K"))) %>% 
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

# func0. creating GRange Obj.
create_granges <- function(df, direction = NULL, chr.col = "chr1", metadata.cols = NULL) {
  if (is.null(direction)) {
    start.col <- "x0"
    end.col <- "y3"
  } else if (direction == "up") {
    start.col <- "x0"
    end.col <- "x3"
  } else if (direction == "down") {
    start.col <- "y0"
    end.col <- "y3"
  } else {
    stop("Invalid direction. Use 'up', 'down', or leave NULL for overall.")
  }
  
  gr <- GRanges(
    seqnames = df[[chr.col]],
    ranges = IRanges(start = as.integer(df[[start.col]]), end = as.integer(df[[end.col]]))
  )
  
  if (is.null(metadata.cols)) {
    metadata.cols <- setdiff(names(df), c(chr.col, start.col, end.col))
  }
  
  mcols(gr) <- df[, metadata.cols, drop = FALSE]
  return(gr)
}

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
OVERALL.df.DISTINCT.loop.deep.sample.all.GR <- create_granges(OVERALL.df.DISTINCT.loop.deep.sample.all)

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

df.ctcf.dist.result %>% dim() # sub.4: 68014890(any), 68012481(within) | 40250853(any, w/o capping + lt 2mb)

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

# func 1. merging figures into one for the overall distribution of CTCF on loops
merging_pdf_pages_to_single_image <- function(pdf_path, output_pdf_path, output_png_path, dpi = 150, stack = TRUE) {
  # Step 1: checking page count
  n_pages <- pdf_info(pdf_path)$pages
  # Step 2: converting each page to image
  pdf_images <- lapply(1:n_pages, function(p) {
    image_read(pdf_render_page(pdf_path, page = p, dpi = dpi))
  })
  # Step 3: merging images vertically
  merged_img <- image_append(do.call(c, pdf_images), stack = stack)
  
  # Step 4: generating output paths
  base_path <- file_path_sans_ext(pdf_path)
  output_pdf_path <- paste0(base_path, "_merged.pdf")
  output_png_path <- paste0(base_path, "_merged.png")
  
  # Step 5: saving PDF and PNG
  image_write(merged_img, path = output_pdf_path, format = "pdf")
  image_write(merged_img, path = output_png_path, format = "png")
  
  return(list(pdf = output_pdf_path, png = output_png_path))
}
# func1 END

# func2. checking component distribution by chromosome
checking_component_distribution <- function(data, component_name, output_dir = "./figures/submission/lt2mb") {
  
  chromosomes <- c(1:20, "X", "Y")
  plot_list <- list()
  
  for (chr in chromosomes) {
    tryCatch({
      message("START: Processing chromosome: ", chr)
      
      data_chr <- data %>% 
        filter(str_detect(loop.id, paste0("_chr", chr, "_")))
      
      # Density plot
      plot_dens <- data_chr %>%
        ggplot(aes(x = value)) +
        geom_density(fill = "skyblue", color = "black", alpha = 0.5) +
        ylim(c(0,1)) +
        labs(
          title = paste0("Density of ", toupper(component_name), " Found over Loop on Chr", chr),
          x = "Relative Position to Loop",
          y = "Density"
        ) +
        theme(plot.title = element_text(hjust = 0.5))
      
      # Histogram
      plot_hist <- data_chr %>%
        ggplot(aes(x = value)) +
        geom_histogram(fill = "skyblue", color = "black", alpha = 0.5, bins = 200) +
        labs(
          title = paste0("Histogram of ", toupper(component_name), " Found over Loop on Chr", chr),
          x = "Relative Position to Loop",
          y = "Count"
        ) +
        theme(plot.title = element_text(hjust = 0.5))
      
      # merging density & histogram
      combined_plot <- plot_hist + plot_dens
      plot_list[[chr]] <- combined_plot
      
      message("END: Processing chromosome: ", chr)
      
    }, error = function(e) {
      message("Error processing chromosome: ", chr)
      message("Error message: ", e$message)
    })
  }

  pdf_path <- file.path(output_dir, paste0("overall_distribution_of_", component_name, "_by_chromosome_wo_capping_lt2mb.pdf"))
  png_path <- paste0(file_path_sans_ext(pdf_path), ".png")
  
  # saving PDF
  pdf(pdf_path, width = 11*0.8, height = 8.5*0.8)
  num_plots <- length(plot_list)
  plots_per_page <- 4
  
  for (i in seq(1, num_plots, by = plots_per_page)) {
    end_idx <- min(i + plots_per_page - 1, num_plots)
    page_plots <- plot_list[i:end_idx]
    combined_page <- plot_grid(plotlist = page_plots, ncol = 1, nrow = 4)
    print(combined_page)
  }

  dev.off()
  
  # converting PDF → PNG 
  merging_pdf_pages_to_single_image(pdf_path, output_pdf_path = pdf_path, output_png_path = png_path)
  
  message("Completed processing for: ", toupper(component_name))
}
# func2 END

checking_component_distribution(relative.pos.df.ctcf.dist.result, "ctcf")

relative.pos.df.ctcf.dist.result %>% head()

########################
# 2. CTCF
# 2-3. overall distribution of CTCF on loops: figures
# 2-3-2. by resolution
########################
# func3. plot_histogram for overall distribution
plot_histogram <- function(df, xvar = value, xlab = "Relative Position to Loop", ylab = "Count") {
  ggplot(df, aes(x = {{ xvar }})) +
    geom_histogram(fill = "skyblue", color = "grey70", alpha = 0.7, bins = 200, linewidth = 0.1) +
    labs(x = xlab, y = ylab) +
    theme(plot.title = element_text(hjust = 0.5))
}
# func3 END

# func4. plot_density for overall distribution
plot_density <- function(df, xvar = value, groupvar = loop.res,
                         xlab = "Relative Position to Loop", ylab = "Density",
                         colors = c("5K" = "#a6cee3", "10K" = "#1f78b4", "25K" = "#1f3a93")) {
  ggplot(df, aes(x = {{ xvar }}, color = {{ groupvar }}, fill = {{ groupvar }})) +
    geom_density(alpha = 0.3) +
    ylim(c(0, 0.75)) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    labs(x = xlab, y = ylab, color = "Resolution", fill = "Resolution") +
    theme(plot.title = element_text(hjust = 0.5))
}
# func4 END

# func5. saving histogram & desity plot combined
saving_combined_plot <- function(plot_left, plot_right, filename, width = 11 * 0.8, height = 8.5 * 0.4) {
  pdf(file = filename, width = width, height = height)
  combined_plot <- plot_left | plot_right
  print(combined_plot)
  dev.off()
}
# func5 END

plot.ctcf.hist <- plot_histogram(relative.pos.df.ctcf.dist.result)
plot.ctcf.dens <- plot_density(relative.pos.df.ctcf.dist.result)

saving_combined_plot(plot.ctcf.hist, plot.ctcf.dens,
                   "figures/submission/lt2mb/overall_distribution_of_CTCF_by_resolution_wo_capping_lt2mb.pdf")

########################
# 2. CTCF
# 2-4. Distribution of CTCF at ends in a loop: for the number of CTCF used in filtering valid loops: figures
# so, the object should be used one without padding: df.DISTINCT.loop.deep.sample.all
# 2-4-3. data processing for getting information of CTCF at ends in loops 
########################
df.DISTINCT.loop.deep.sample.all.lt.2mb # 31019
df.DISTINCT.loop.deep.sample.all.up.GR   <- create_granges(df.DISTINCT.loop.deep.sample.all.lt.2mb, direction = "up")
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

df.DISTINCT.loop.deep.sample.all.down.GR   <- create_granges(df.DISTINCT.loop.deep.sample.all.lt.2mb, direction = "down")
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
# tss file list
# Linux
file.tss.list = fs::dir_ls("/home/pkim/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/", regexp = "\\.txt$")
# Mac
file.tss.list = fs::dir_ls("~/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss", regexp = "\\.txt$")

file.tss.list
# csRNA.NuAcc.tss.txt
# csRNA.PFC.tss.txt
# ucsc_start_codon.txt

# GTF 
df.refgene.gtf <- read_tsv("~/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/ucsc_refGene.gtf",
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

# df.refgene.gtf.parsed.start.codon == df.tss.ucsc.tss.id 

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

# Locus ID or RGD Gene ID (21) to be allocated with gene symbol (gene_name)
loc.map.tss <- c(
  "LOC100360846" = "Psmb6l1", # Locus ID duplicated: NM_001329883 picked
  "LOC100909648" = "Avpr1b",
  "LOC100911576" = "Hnrnpc",
  "LOC102553861" = "Gzmb",
  "LOC100911664" = "C17h6orf62",
  "LOC100365112" = "Ces2",
  "LOC103689947" = "Selenbp1",
  "LOC100911796" = "Adora3",  # Locus ID duplicated: SAME COORDINATES
  "LOC361914"     = "Slc7a12",
  "LOC108348108"  = "Hspa1a",
  "LOC102548514"  = "Crnkl1",
  "LOC100911068"  = "Robo4",
  "LOC689600"     = "1700020N15Rikl", # SAME COORDINATES ****
  "RGD1564599"    = "1700020N15Rikl", # SAME COORDINATES ****
  "LOC497796"     = "Klra17",
  "RGD1565355"    = "Cd36-ps1", # RGD ID duplicated: NM_001109218 picked
  "LOC259246"     = "Mup5", # SAME COORDINATES ***
  "LOC298116"     = "Mup5", # SAME COORDINATES ***
  "LOC498750"    = "C17h6orf62-ps1",
  "LOC688090"    = "RT1-Bb",
  "RGD2301395"    = "Klrb1al"
)

# Mup4l1 (Ensembl: Mup5)
# tss.select.sparate %>% filter(str_detect(gene_id, "LOC259246|LOC298109|LOC298116|LOC366379|LOC688514"))
# tss.select.sparate %>% filter(str_detect(gene_id, "Hnrnpc"))

gene_list_tss <- c(
  "Psmb6l1", "Avpr1b", "Hnrnpc", "Gzmb", "C17h6orf62", 
  "Ces2", "Selenbp1", "Adora3", "Slc7a12", "Hspa1a", 
  "Crnkl1", "Robo4", "1700020N15Rikl", "Klra17", "Cd36-ps1", 
  "Mup4l1", "Mup5", "Mup", "Mup4l2", "C17h6orf62-ps1", "RT1-Bb", "Klrb1al"
)

gene_pattern_tss <- str_c(gene_list_tss, collapse = "|")
gene_pattern_tss
loc_pattern_tsssss <- str_c(names(loc.map.tss), collapse = "|")
loc_pattern_tsssss
gene_matching <- tss.select.sparate %>% filter(str_detect(gene_id, loc_pattern_tsssss)) 
gene_matching # 24 rows (duplicated: LOC100360846, LOC100911796, RGD1565355)

tss.select.sparate %>% dim() # 17849
df.tss.ucsc.tss.id <- tss.select.sparate %>% 
  dplyr::rename(chr = V1, start = V4, end = V5, strand = V7 ) %>% 
  dplyr::mutate(tss.id = paste(chr, start, end, strand, gene_id, transcript_id, sep = ":")) %>% 
  mutate(gene_id = str_replace_all(gene_id, loc.map.tss)) %>% # 17,849
  filter(gene_id != 'LOC259245') # LOC259245: withdrawn| 17,848 ********************************

df.tss.ucsc.tss.id %>% head(20)
df.tss.ucsc.tss.id %>% dim() # dim: 17848
# df.tss.ucsc.tss.id %>% filter(str_detect(gene_id, "Hnrnpc")) %>% View()

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
df.tss.ucsc.tss.id %>% distinct(chr, start, end, gene_id) %>% dim() # 17125
df.tss.ucsc.tss.id %>% distinct(chr, start, end, strand, gene_id) # 17125
df.tss.ucsc.tss.id %>% distinct(chr, start, end, strand, gene_id, gene_name) # 17138

###################
# exploring data of TSS
###################

# checking duplicates by counting rows b/w rows and gene rows
n_coord_tss      <- df.tss.ucsc.tss.id %>% distinct(chr, start, end) %>% nrow() # 17080
n_coord_tss
n_coord_gene_tss <- df.tss.ucsc.tss.id %>% distinct(chr, start, end, gene_id) %>% nrow() # 17125 **********************
n_coord_gene_tss
delta_tss <- n_coord_gene_tss - n_coord_tss
delta_tss # 45 duplicated coord with different gene rows: max gene number lt 45

multi_gene_per_coord_tss <- df.tss.ucsc.tss.id %>% 
  distinct(chr, start, end, gene_id) %>%
  group_by(chr, start, end) %>%
  summarise(
    n_genes = n_distinct(gene_id),
    gene_ids = paste(sort(unique(gene_id)), collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(n_genes > 1) %>%
  arrange(desc(n_genes), chr, start, end)

multi_gene_per_coord_tss # mutiple genes at same coord: 40 rows
print(multi_gene_per_coord_tss, n=Inf)
# write_xlsx(multi_gene_per_coord_tss, "multi_gene_per_coord_tss.xlsx")

# 2) 위 케이스들의 총 초과 개수가 실제 차이와 일치하는지 검증
#    (각 좌표에서 n_genes-1의 합이 delta와 같아야 함)
sum_excess_tss <- multi_gene_per_coord_tss %>% summarise(sum(n_genes - 1)) %>% pull()
all.equal(sum_excess_tss, delta_tss)

# 3) 예시 몇 개만 확인
head(multi_gene_per_coord_tss, 20)
print(multi_gene_per_coord_tss, n = Inf)

# 4) 만약 동일 좌표 & 동일 gene_id가 데이터에 중복 존재했는지(집계 전 원시 중복) 확인
raw_dups_same_gene <- df.tss.ucsc.tss.id %>%
  count(chr, start, end, gene_id, name = "n") %>%
  filter(n > 1) %>%
  arrange(desc(n))
raw_dups_same_gene

# 5) gene_id가 비어 있는 데이터가 있는지 참고(숫자 차이를 만들진 않지만 품질 점검용)
has_na_gene_tss <- df.tss.ucsc.tss.id %>% filter(is.na(gene_id)) %>% nrow()
has_na_gene_tss # 0

########
map_gene_tss <- df.tss.ucsc.tss.id %>%
  group_by(chr, start, end) %>%
  summarise(
    gene_id = {
      ids <- unique(gene_id)
      ids <- ids[!is.na(ids)]                 # NA 제외
      if (length(ids) == 0) NA_character_ else paste(sort(ids), collapse = "]")
    },
    .groups = "drop"
  )
map_gene_tss
map_gene_tss %>% filter(str_detect(gene_id, "\\]") == TRUE) # 40

df.tss.ucsc <- df.tss.ucsc.tss.id %>%  # 17848
  distinct(chr, start, end) %>% # 17080******
  left_join(map_gene_tss, by = c("chr", "start", "end")) %>%
  mutate(tss.id = paste(chr, start, end, gene_id, sep = ":")) %>% 
  filter(!str_detect(gene_id, "\\]")) # - 40

df.tss.ucsc %>% head()# 17040/17080: chr   start     end gene_id                      tss.id
df.tss.ucsc %>% dim() # 17040
df.tss.ucsc %>% add_count(gene_id) %>% filter(n > 1)
df.tss.ucsc %>% filter(str_detect(gene_id, gene_pattern_tss)) # 0 rows

df.tss.ucsc %>% add_count(chr, start, end, gene_id) %>% filter(n > 1) # 0 rows

############################
############################
df.ensembl.gtf <- read_tsv("~/dropbox/Gateway_to_Hao/enhancer/data/Rattus_norvegicus.mRatBN7.2.113.gtf", 
                comment = "#", 
                col_names = FALSE) 

df.ensembl.gtf %>% head()
df.ensembl.gtf %>% count(X3)
df.ensembl.gtf %>% dplyr::select(X9) %>% head(3)

colnames(df.ensembl.gtf) <- c("chr", "source", "feature", "start", "end", 
                   "score", "strand", "frame", "attribute")

df.ensembl.gtf.for.tss <- df.ensembl.gtf %>% 
  filter(feature == "start_codon") %>% 
  filter(chr %in% c(as.character(1:20), "X", "Y")) %>% 
  mutate(chr = str_c("chr", chr))

df.ensembl.gtf.for.tss %>% dim() # 42925
df.ensembl.gtf.for.tss %>% head(3)
df.ensembl.gtf.for.tss %>% dplyr::select(attribute) %>% head(3)

# checking attribute keys
df.ensembl.gtf.for.tss %>%
  pull(attribute) %>%                           
  strsplit(";") %>%                             
  unlist() %>%                                  
  trimws() %>%                                  
  str_extract("^[^ ]+") %>%                     
  unique() %>%                                   
  sort()              

# function for attributes
extract_attributes <- function(attr_string) {
  safe_extract <- function(key, string) {
    match <- str_match(string, paste0(key, ' "([^"]*)?"'))[,2]
    if (is.na(match) || trimws(match) == "") return(NA_character_) else return(match)
  }
  
  tibble(
    gene_id        = safe_extract("gene_id", attr_string),
    gene_name      = safe_extract("gene_name", attr_string),
    gene_biotype   = safe_extract("gene_biotype", attr_string),
    tag                  = safe_extract("tag", attr_string),
    transcript_biotype   = safe_extract("transcript_biotype", attr_string),
    transcript_version   = safe_extract("transcript_version", attr_string),
    exon_number    = safe_extract("exon_number", attr_string)
  )
}

df.ensembl.gtf.for.tss # 42,925

df.ensembl.gtf.for.tss.attribute <- df.ensembl.gtf.for.tss %>% 
  bind_cols( df.ensembl.gtf.for.tss$attribute %>% lapply(extract_attributes) %>% bind_rows() ) %>%
  dplyr::select(-c(source, feature, attribute, score, frame)) 

df.ensembl.gtf.for.tss.attribute

df.ensembl.gtf.for.tss.attribute.filtered.tag.gnbt <- df.ensembl.gtf.for.tss.attribute %>% 
  filter(tag == "Ensembl_canonical") %>% # 21766
  filter(gene_biotype == "protein_coding") %>% # 21760
  mutate(
    transcript_version = as.numeric(transcript_version),
    exon_number = as.numeric(exon_number)
  )
df.ensembl.gtf.for.tss.attribute.filtered.tag.gnbt %>% dim() # 21760
df.ensembl.gtf.for.tss.attribute.filtered.tag.gnbt %>% head() # 21760

df.ensembl.gtf.for.tss.attribute.filtered.tag.gnbt.geneid.NOdup <- df.ensembl.gtf.for.tss.attribute.filtered.tag.gnbt %>% 
  add_count(gene_id) %>%
  filter(n == 1) # 21690

df.ensembl.gtf.for.tss.attribute.filtered.tag.gnbt.geneid.dup.1pick <- df.ensembl.gtf.for.tss.attribute.filtered.tag.gnbt %>% 
  add_count(gene_id) %>%
  filter(n > 1) %>%
  group_by(gene_id) %>%
  arrange(desc(transcript_version), exon_number) %>%
  slice_head(n = 1) %>%
  ungroup() # 35
  
df.ensembl.gtf.for.tss.DISTINCT.geneid <- bind_rows(df.ensembl.gtf.for.tss.attribute.filtered.tag.gnbt.geneid.NOdup,
                                                                                   df.ensembl.gtf.for.tss.attribute.filtered.tag.gnbt.geneid.dup.1pick)  %>% # 21776
                                          dplyr::select(-n) %>% 
                                          mutate(tss.id = paste(chr, start, end, strand, gene_id, sep = ":")) %>% 
                                          dplyr::select(chr, start, end, strand, gene_id, gene_name, tss.id) 

  
df.ensembl.gtf.for.tss.DISTINCT.geneid %>% # 21725
  # count(chr, start, end, gene_id) #%>% # 21,725
  count(chr, start, end, strand, gene_id) #%>% # 21,725

df.ensembl.gtf.for.tss.DISTINCT.geneid # %>% 
  # count(chr) %>% print(n = Inf)
  # head(3) # chr     start       end strand gene_id        gene_name
  # dim() # 21,725
  # distinct(gene_id) # 21,725
  # filter(str_detect(gene_id, 'LOC|RGD')) # 0

df.tss.ucsc <- df.ensembl.gtf.for.tss.DISTINCT.geneid # 21725

# GRanges Obj.: df.tss.ucsc.GR: chr   start     end gene_id                      tss.id
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

df.tss.ucsc.GR # 17080//21,725
##### tss with nuacc, 96563
# nuacc <- read.csv("/home/pkim/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/csRNA.NuAcc.tss.txt", header = T, sep = '\t')
# nuacc <- read.csv("/Users/PanjunKim/UTHSC GGI Dropbox/K P/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/csRNA.NuAcc.tss.txt", header = T, sep = '\t')

##### tss with pfc ver1, 131647
# pfc <- read.csv("/home/pkim/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/csRNA.PFC.tss.txt", header = T, sep = '\t')
# pfc <- read.csv("/Users/PanjunKim/UTHSC GGI Dropbox/K P/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/csRNA.PFC.tss.txt", header = T, sep = '\t')

########################
# 3. TSS
# 3-3. overall distribution of TSS on loops
########################
OVERALL.df.DISTINCT.loop.deep.sample.all.GR # 30928

index.distinct.tss.w.OVERALL.whole.loop <- findOverlaps(
  df.tss.ucsc.GR, 
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
  tss_chr = df.tss.ucsc$chr[OVERALL.tss.on.loop.hits],
  tss_start = df.tss.ucsc$start[OVERALL.tss.on.loop.hits],
  tss_end = df.tss.ucsc$end[OVERALL.tss.on.loop.hits],
  tss_id = df.tss.ucsc$tss.id[OVERALL.tss.on.loop.hits],
  tss_geneid = df.tss.ucsc$gene_id[OVERALL.tss.on.loop.hits],
  tss_strand = df.tss.ucsc$strand[OVERALL.tss.on.loop.hits]
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

checking_component_distribution(relative.pos.df.tss.dist.result, "tss")

relative.pos.df.tss.dist.result %>% head()

########################
# 3. TSS
# 3-3. overall distribution of TSS on loops: figures
# 3-3-2. by resolution
########################

plot.tss.hist <- plot_histogram(relative.pos.df.tss.dist.result)
plot.tss.dens <- plot_density(relative.pos.df.tss.dist.result)

saving_combined_plot(plot.tss.hist, plot.tss.dens,
                   "figures/submission/lt2mb/overall_distribution_of_TSS_by_resolution_wo_capping_lt2mb.pdf")
                   
########################
# 3. TSS
# 3-4. Distribution of TSS at each end in a loop: for the number of TSS used in filtering valid loops: figures
# so, the object should be used one WITH padding on df.DISTINCT.loop.deep.sample.all: 1/2 distance for OUTER & 1/4 distance for INNER
# 3-4-1. adding padding at each end, x12 & y12
########################
df.DISTINCT.loop.deep.sample.all # 31773
df.DISTINCT.loop.deep.sample.all.lt.2mb  %>% dim()# 31019
df.DISTINCT.loop.deep.sample.all.lt.2mb %>% head(2) # 31019

# func7. Computing inner distance statistics
computing_inner_distance_stats <- function(df) {
  df <- df %>% mutate(inner.distance = y0 - x3)
  
  # quantile
  ylim_vals <- quantile(df$inner.distance, c(0.05, 0.95), na.rm = TRUE)
  
  # min & median
  stats <- df %>%
    summarise(
      min_val = min(inner.distance, na.rm = TRUE),
      median_val = median(inner.distance, na.rm = TRUE)
    )
  
  # return: list
  list(ylim_vals = ylim_vals, stats = stats)
}
# func7 END

# func8. Plotting inner distance boxplot
plotting_inner_distance_boxplot <- function(df, ylim_vals, stats, title_text) {
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
    theme(
      axis.text.x = element_text(angle = 0, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom"
    )
}
# func8 END

########################
# 3. TSS
# 3-4. Distribution of TSS at ends in a loop: for the number of TSS used in filtering valid loops: figures
# so, the object should be used one WITH padding: df.DISTINCT.loop.deep.sample.all.padded.for.TSS
# 3-4-3. data processing for getting information of TSS at ends in loops 
########################
df.DISTINCT.loop.deep.sample.all.padded.for.TSS.promoter <- df.DISTINCT.loop.deep.sample.all.lt.2mb # 31019

df.DISTINCT.loop.deep.sample.all.padded.for.TSS.promoter %>% dim() # 31019
df.DISTINCT.loop.deep.sample.all.padded.for.TSS.promoter %>% head(3)
df.DISTINCT.loop.deep.sample.all.padded.for.TSS.promoter %>% colnames()

# 1. TSS + UPSTREAM (df.DISTINCT.loop.deep.sample.all.padded.for.TSS.promoter.UP.GR, df.tss.ucsc.GR)
df.DISTINCT.loop.deep.sample.all.padded.for.TSS.promoter.UP.GR   <- create_granges(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.promoter, direction = "up")
df.DISTINCT.loop.deep.sample.all.padded.for.TSS.promoter.UP.GR # 31019

index.distinct.tss.w.up.loop.each.end <- findOverlaps(
  df.tss.ucsc.GR, 
  df.DISTINCT.loop.deep.sample.all.padded.for.TSS.promoter.UP.GR, 
  type = "any",
  select = "all"
)

index.distinct.tss.w.up.loop.each.end # any: 91559// default trial: 8074| 6case: 9290| 7case: 12050| 8case: 14702| 9case: 19908

end.loop.up.tss.hits <- subjectHits(index.distinct.tss.w.up.loop.each.end)
end.tss.up.hits <- queryHits(index.distinct.tss.w.up.loop.each.end)

df.tss.ucsc.GR

df.overlapping.TSS.w.UPSTREAM.result.each.end <- tibble(
  up.loop.id = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.promoter.UP.GR)$loop.id[end.loop.up.tss.hits],
  end.up.distance = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.promoter.UP.GR)$end.distance[end.loop.up.tss.hits],
  resolution = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.promoter.UP.GR)$resolution[end.loop.up.tss.hits],
  tss.id = mcols(df.tss.ucsc.GR)$tss_id[end.tss.up.hits],
  tss.gene_id = mcols(df.tss.ucsc.GR)$gene_id[end.tss.up.hits],
  WHERE = "UP"
) %>% 
  mutate(tss.loop.up.id = str_c(up.loop.id, '|', tss.id, '|', WHERE))

df.overlapping.TSS.w.UPSTREAM.result.each.end %>% dim() # any: 91559     6 // default trial: 8074| 6case: 9290| 7case: 12050| 8case: 14702| 9case: 19908
df.overlapping.TSS.w.UPSTREAM.result.each.end %>% head()

df.overlapping.TSS.w.UPSTREAM.result.each.end %>% 
  count(end.up.distance)
colnames(df.overlapping.TSS.w.UPSTREAM.result.each.end)
df.overlapping.TSS.w.UPSTREAM.result.each.end %>% dplyr::select(tss.loop.up.id)

# any                         # case default            # case6                   # case7                 # case8                   # case9
# end.up.distance     n
# 1            5000 20414     1            5000  1706   1            5000  2167   1            5000  2643 1            5000  3126   # 1            5000  4056
# 2           10000 27145     2           10000  2787   2           10000  3542   2           10000  4295 2           10000  4996   # 2           10000  6483
# 3           25000 44000     3           25000  3581   3           25000  3581   3           25000  5112 3           25000  6580   # 3           25000  9369
# total: 91559 (dedups) 

# 2. TSS + DOWNSTREAM (df.DISTINCT.loop.deep.sample.all.padded.for.TSS.promoter.DOWN.GR, df.tss.ucsc.GR)
df.DISTINCT.loop.deep.sample.all.padded.for.TSS.promoter.DOWN.GR   <- create_granges(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.promoter, direction = "down")
df.DISTINCT.loop.deep.sample.all.padded.for.TSS.promoter.DOWN.GR # 31019

index.distinct.tss.w.down.loop.each.end <- findOverlaps(
  df.tss.ucsc.GR, 
  df.DISTINCT.loop.deep.sample.all.padded.for.TSS.promoter.DOWN.GR, 
  type = "any",
  select = "all"
)
index.distinct.tss.w.down.loop.each.end # any: 106271, within: 106269 // default trial: 7932| 6case: 9040| 7case: 11556| 8case: 13989| 9case: 18910

end.loop.down.tss.hits <- subjectHits(index.distinct.tss.w.down.loop.each.end)
end.tss.down.hits.end <- queryHits(index.distinct.tss.w.down.loop.each.end)

df.overlapping.TSS.w.DOWNSTREAM.result.each.end <- tibble(
  down.loop.id = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.promoter.DOWN.GR)$loop.id[end.loop.down.tss.hits],
  end.down.distance = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.promoter.DOWN.GR)$end.distance[end.loop.down.tss.hits],
  resolution = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.promoter.DOWN.GR)$resolution[end.loop.down.tss.hits],
  tss.id = mcols(df.tss.ucsc.GR)$tss_id[end.tss.down.hits.end],
  tss.gene_id = mcols(df.tss.ucsc.GR)$gene_id[end.tss.down.hits.end],
  WHERE = "DOWN"
) %>%
  mutate(tss.loop.down.id = str_c(down.loop.id, '|', tss.id, '|', WHERE))

df.overlapping.TSS.w.DOWNSTREAM.result.each.end %>% dim() # 106271 (dedup) // default trial: 7932| 6case: 9040| 7case: 11556| 8case: 13989| 9case: 18910
df.overlapping.TSS.w.DOWNSTREAM.result.each.end %>% head()
df.overlapping.TSS.w.DOWNSTREAM.result.each.end %>% 
  count(end.down.distance)
colnames(df.overlapping.TSS.w.DOWNSTREAM.result.each.end)
# any
# end.down.distance     n     # default                    # case6                     # case7                   # case8                       # case9
# 1              5000 24273   1              5000  1652    1              5000  2066   1              5000  2484 1              5000  2877     1              5000  3745
# 2             10000 30197   2             10000  2746    2             10000  3440   2             10000  4180 2             10000  4819     2             10000  6101
# 3             25000 51801   3             25000  3534    3             25000  3534   3             25000  4892 3             25000  6293     3             25000  9064
# total: 106271 (dedup) // 111334 (dups)

########## bind_rows(UPSTREAM & DOWNSTREAM) -> BOTH
df.overlapping.TSS.w.DOWNSTREAM.result.each.end %>% head(3)

df.overlapping.TSS.w.BOTH.result <- bind_rows(df.overlapping.TSS.w.UPSTREAM.result.each.end %>% 
                                                dplyr::rename(loop.id = up.loop.id, end.distance = end.up.distance, case.id = tss.loop.up.id) %>% 
                                                mutate(tss.gene_id = str_c(tss.gene_id, '|', 'UP')),
                                              df.overlapping.TSS.w.DOWNSTREAM.result.each.end %>% 
                                                dplyr::rename(loop.id = down.loop.id, end.distance = end.down.distance, case.id = tss.loop.down.id) %>% 
                                                mutate(tss.gene_id = str_c(tss.gene_id, '|', 'DOWN'))) %>% 
  mutate(tss_result_both_id = str_c(loop.id, '|', tss.gene_id)) %>% 
  mutate(case_sort_id = str_c(loop.id, '|', WHERE)) %>%
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

df.overlapping.TSS.w.BOTH.result %>% dim() # any: 197830      7 (dedup) // default: 16006| case6: 18330| case7: 23606| case8: 28691| case9: 38818
df.overlapping.TSS.w.BOTH.result %>% colnames() # "resolution"   "tss.id"       "WHERE"        "loop.id"      "end.distance" "case.id"      "chr"         
df.overlapping.TSS.w.BOTH.result %>% head()
df.overlapping.TSS.w.BOTH.result %>% count(resolution)
df.overlapping.TSS.w.BOTH.result %>% dplyr::select(case.id) %>% filter(!is.na(case.id))
# any
# resolution      n   # case default          # case6             # case7              # case8              # case9
# 1 5K         446871 1 5K          3358      5K          4233    1 5K          5127   1 5K          6003   1 5K          7801   
# 2 10K        573422 2 10K         5533      10K         6982    2 10K         8475   2 10K         9815   2 10K        12584   
# 3 25K        958013 3 25K         7115      25K         7115    3 25K        10004   3 25K        12873   3 25K        18433

df.case.from.tss <- df.overlapping.TSS.w.BOTH.result %>% 
  group_by(loop.id) %>% 
  summarise(
    WHERE_list = list(WHERE),
    WHERE_uniq = list(base::unique(WHERE)),
    .groups = "drop"
  ) %>% 
  rowwise() %>%
  mutate(case = case_when(
    length(WHERE_uniq) == 2 ~ "BOTH",
    length(WHERE_uniq) == 1 & WHERE_uniq[[1]] == "UP" & length(unlist(WHERE_list)) == 1 ~ "u.UP",
    length(WHERE_uniq) == 1 & WHERE_uniq[[1]] == "UP" & length(unlist(WHERE_list)) > 1 ~ "m.UP",
    length(WHERE_uniq) == 1 & WHERE_uniq[[1]] == "DOWN" & length(unlist(WHERE_list)) == 1 ~ "u.DOWN",
    length(WHERE_uniq) == 1 & WHERE_uniq[[1]] == "DOWN" & length(unlist(WHERE_list)) > 1 ~ "m.DOWN",
    TRUE ~ NA_character_
  )) %>% 
  ungroup()

df.case.from.tss %>% head()
df.case.from.tss %>% count(case)
# default       n      # case6               # case7               # case8             # case9           
# 1 BOTH    2139    1 BOTH    2615        1 BOTH    3719        1 BOTH    4803      1 BOTH    6852
# 2 m.DOWN   437    2 m.DOWN   508        2 m.DOWN   715        2 m.DOWN   899      2 m.DOWN  1187
# 3 m.UP     445    3 m.UP     564        3 m.UP     807        3 m.UP    1055      3 m.UP    1485
# 4 u.DOWN  4368    4 u.DOWN  4651        4 u.DOWN  5021        4 u.DOWN  5122      4 u.DOWN  5105
# 5 u.UP    4478    5 u.UP    4745        5 u.UP    5179        5 u.UP    5325      5 u.UP    5298

## debug
# df.tss.debug <- df.overlapping.TSS.w.BOTH.result %>%
#   group_by(loop.id) %>%
#   summarise(
#     WHERE_list = list(WHERE),
#     WHERE_uniq = list(base::unique(WHERE)),
#     .groups = "drop"
#   )

# View(df.tss.debug)

# checking: loop.id case exclusive 
check.exclusive.tss <- df.case.from.tss %>%
  group_by(loop.id) %>%
  summarise(
    n_case = n_distinct(case),
    case_list = paste(unique(case), collapse = ","),
    .groups = "drop"
  ) %>%
  filter(n_case > 1)  # loop.id involved in more than 1 case

check.exclusive.tss # 0

df.final.loop.dataset.tss <- df.overlapping.TSS.w.BOTH.result %>% left_join(df.case.from.tss, by = 'loop.id')

df.final.loop.dataset.tss %>% head(3)
df.final.loop.dataset.tss %>% dim() # 31019    17
df.final.loop.dataset.tss %>% count(case)

# default       n      # case6           # case7           # case8             # case9
# 1 BOTH    5277    1 BOTH    6631    1 BOTH   10055    1 BOTH   13848      1 BOTH   22215
# 2 m.DOWN   919    2 m.DOWN  1072    2 m.DOWN  1545    2 m.DOWN  1999      2 m.DOWN  2742
# 3 m.UP     964    3 m.UP    1231    3 m.UP    1806    3 m.UP    2397      3 m.UP    3458
# 4 u.DOWN  4368    4 u.DOWN  4651    4 u.DOWN  5021    4 u.DOWN  5122      4 u.DOWN  5105
# 5 u.UP    4478    5 u.UP    4745    5 u.UP    5179    5 u.UP    5325      5 u.UP    5298


###############
###############
# for Q1: quartile among loops with TSS, so min = 1
df.tss.counts <- df.overlapping.TSS.w.BOTH.result %>%
  group_by(loop.id, WHERE, resolution) %>%
  summarise(tss_count = n_distinct(tss.id), .groups = 'drop')

df.tss.counts %>% count(tss_count)

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

df.tss.case
df.tss.case %>% count(case)

df.loop.with.tss.case <- df.DISTINCT.loop.deep.sample.all.lt.2mb %>% # 31773
  left_join(df.tss.case, by = "loop.id") %>%
  mutate(case = ifelse(is.na(case), "NONE", case))

df.loop.with.tss.case %>% count(case)
df.tss.counts

########################
# 4. promoter
# 4-1. Exploratory Data analysis (EDA) : promoter
# 4-2. promoter data preprocessing: id and dedup GRange Obj.
########################
# file_path <- "/Users/PanjunKim/UTHSC GGI Dropbox/K P/Gateway_to_Hao/enhancer/data/epdnew/001/Rn_EPDnew_001_rn7.bed"
file_path <- "~/dropbox/Gateway_to_Hao/enhancer/data/epdnew/001/Rn_EPDnew_001_rn7_1.bed"
df.promoter.rn7.raw <- read_tsv(file_path, col_names = c("chr", "start", "end", "gene", "score"))
gene_pattern_promoter <- 'LOC103689961|LOC100909700|NEWGENE_1310561|LOC100911365|LOC100911881|LOC100912537|LOC103689931|LOC100909661|LOC100912534|LOC100911558|LOC100911417|LOC103691744|LOC108348108|LOC103690139|LOC103690006|RGD1566134|LOC100912026|LOC100911130|LOC103689993|LOC108348144'

############################
df.promoter.rn7.raw %>% dim()# 12,463
df.promoter.rn7.raw %>% head(3)
df.promoter.rn7.raw.ENSMBL <- df.promoter.rn7.raw %>% # 12,463
   filter(!str_detect(gene, ":::")) %>% # 1
   mutate(gene = str_replace_all(gene, "::", ":NA:")) %>%
  # mutate(colon_count = stringr::str_count(gene, ":")) %>% count(colon_count) # 3 12463
  separate(
    gene,
    into = c("iso_id", "gene_id", "refseq_id", "gene_name"),
    sep = ":",
    fill = "right",
    extra = "drop"
  ) 
  
df.promoter.rn7.raw.ENSMBL.single <- df.promoter.rn7.raw.ENSMBL %>% # 12,462
  group_by(chr, start, end, gene_id) %>%
  filter(n() == 1) %>% # 12454
  ungroup()

df.promoter.rn7.raw.ENSMBL.dup <- df.promoter.rn7.raw.ENSMBL %>% # 12,462
  group_by(chr, start, end, gene_id) %>%
  filter(n() > 1) %>% # 8
  filter(!str_detect(refseq_id, "NA")) %>%
  ungroup() # 4

df.promoter.rn7.raw.ENSMBL.dedup <- bind_rows(
  df.promoter.rn7.raw.ENSMBL.single,
  df.promoter.rn7.raw.ENSMBL.dup
) # 12,458

df.promoter.rn7.raw.ENSMBL.dedup %>% 
  # distinct(gene_id) %>% nrow() # 12,458
  # count(chr, start, end) %>% filter(n > 1) # 0
  # count(chr, start, end, gene_id) %>% filter(n > 1) # 0
  # count(chr, start, end, gene_id, gene_name) %>% filter(n > 1) # 0
#  filter(str_detect(gene_id, "NA")) %>% # 0
  head(3)
df.promoter.rn7 <- df.promoter.rn7.raw.ENSMBL.dedup %>%
  mutate(length = end - start, 
         center = round((end + start)/2), 
         promoter.id = paste(chr, start, end, gene_id, sep = ':'))

df.promoter.rn7 %>% dim() # 12,458
df.promoter.rn7 %>% head()

# checking duplicates by counting rows b/w rows and gene rows
n_coord_promoter <- df.promoter.rn7.raw %>% distinct(chr, start, end) %>% nrow() 
n_coord_promoter # 12427
n_coord_gene_promoter <- df.promoter.rn7.raw %>% 
  ######## data data cleaning
  mutate(gene = str_replace(gene, "LOC100912405", "Rupl-ps1")) %>% 
  filter(!str_detect(gene, gene_pattern_promoter)) %>% 
  ######## gene data cleaning
  mutate(gene_id = ifelse(str_detect(gene, ":"), str_extract(gene, "[^:]+(?=$)"), gene)) %>% 
  mutate(gene_id = ifelse(is.na(gene_id), str_extract(gene, "^[^:]+"), gene_id)) %>% 
  distinct(chr, start, end, gene_id) %>% nrow() 
n_coord_gene_promoter # 12438
delta_promoter <- n_coord_gene_promoter - n_coord_promoter
delta_promoter # 11 duplicated gene rows: min gene number lt 11

multi_gene_per_coord_promoter <- df.promoter.rn7.raw %>%
  ######## data data cleaning
  mutate(gene = str_replace(gene, "LOC100912405", "Rupl-ps1")) %>% 
  filter(!str_detect(gene, gene_pattern_promoter)) %>% 
  ######## gene data cleaning
  mutate(gene_id = ifelse(str_detect(gene, ":"), str_extract(gene, "[^:]+(?=$)"), gene)) %>% 
  mutate(gene_id = ifelse(is.na(gene_id), str_extract(gene, "^[^:]+"), gene_id)) %>%  
  distinct(chr, start, end, gene_id) %>%
  group_by(chr, start, end) %>%
  summarise(
    n_genes = n_distinct(gene_id),
    gene_ids = paste(sort(unique(gene_id)), collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(n_genes > 1) %>%
  arrange(desc(n_genes), chr, start, end)

multi_gene_per_coord_promoter # 11 coords for multi genes

# write_xlsx(multi_gene_per_coord_promoter, "multi_gene_per_coord_promoter.xlsx")

# 2) 위 케이스들의 총 초과 개수가 실제 차이와 일치하는지 검증
#    (각 좌표에서 n_genes-1의 합이 delta와 같아야 함)
sum_excess_promoter <- multi_gene_per_coord_promoter %>% summarise(sum(n_genes - 1)) %>% pull()
all.equal(sum_excess_promoter, delta_promoter)

# 3) 예시 몇 개만 확인
head(multi_gene_per_coord_promoter, 20)
print(multi_gene_per_coord_promoter, n = Inf)

# 4) 만약 동일 좌표·동일 gene_id가 데이터에 중복 존재했는지(집계 전 원시 중복) 확인
raw_dups_same_gene_promoter <- df.promoter.rn7.raw %>% mutate(gene_id = ifelse(str_detect(gene, ":"), str_extract(gene, "[^:]+(?=$)"), gene)) %>%
  count(chr, start, end, gene_id, name = "n") %>%
  filter(n > 1) %>%
  arrange(desc(n))
head(raw_dups_same_gene_promoter)

# 5) gene_id가 비어 있는 데이터가 있는지 참고(숫자 차이를 만들진 않지만 품질 점검용)
has_na_gene_promoter <- df.promoter.rn7.raw %>% filter(is.na(gene)) %>% nrow() # 0
has_na_gene_promoter <- df.promoter.rn7.raw %>% mutate(gene_id = ifelse(str_detect(gene, ":"), str_extract(gene, "[^:]+(?=$)"), gene)) %>% 
  nrow() # 12,463
has_na_gene_promoter

########
map_gene_promoter <- df.promoter.rn7.raw %>% 
  ######## data data cleaning
  mutate(gene = str_replace(gene, "LOC100912405", "Rupl-ps1")) %>% 
  filter(!str_detect(gene, gene_pattern_promoter)) %>% 
  ######## gene data cleaning
  mutate(gene_id = ifelse(str_detect(gene, ":"), str_extract(gene, "[^:]+(?=$)"), gene)) %>% 
  mutate(gene_id = ifelse(is.na(gene_id), str_extract(gene, "^[^:]+"), gene_id)) %>% 
  group_by(chr, start, end) %>%
  summarise(
    gene_id = {
      ids <- unique(gene_id)
      ids <- ids[!is.na(ids)]                 # NA 제외
      if (length(ids) == 0) NA_character_ else paste(sort(ids), collapse = "]")
    },
    .groups = "drop"
  )
map_gene_promoter
map_gene_promoter %>% filter(str_detect(gene_id, "\\]") == TRUE) # 11

df.promoter.rn7 <- df.promoter.rn7.raw %>%
  ######## data data cleaning
  mutate(gene = str_replace(gene, "LOC100912405", "Rupl-ps1")) %>% 
  filter(!str_detect(gene, gene_pattern_promoter)) %>% 
  ######## gene data cleaning
  distinct(chr, start, end) %>% # 12427******
  left_join(map_gene_promoter, by = c("chr", "start", "end")) %>%
  mutate(length = end - start, 
         center = round((end + start)/2), 
         promoter.id = paste(chr, start, end, gene_id, sep = ':'))

df.promoter.rn7 %>% dim() # g.cleanded: 12427
df.promoter.rn7 %>% head(3) # chr     start     end gene_id length  center promoter.id 
df.promoter.rn7 %>% filter(str_detect(gene_id, '\\]')) # 11

# checking genes with the ones from tss
df_filtered.test <- df.promoter.rn7 %>%
  filter(str_detect(gene_id, gene_pattern_tss))
df_filtered.test
gene_pattern_tss

# dedup GRange Obj.
df.promoter.rn7.GR <- GRanges(
  seqnames = df.promoter.rn7$chr,
  ranges = IRanges(
    start = df.promoter.rn7$start, 
    end = df.promoter.rn7$end)
)

# metadata
mcols(df.promoter.rn7.GR) <- df.promoter.rn7[, c("promoter.id", "center", "length", "gene_id")]

########################
# 4. promoter
# 4-3. overall distribution of promoter on loops
########################
index.promoter.w.OVERALL.whole.loop <- findOverlaps(
  df.promoter.rn7.GR, 
  OVERALL.df.DISTINCT.loop.deep.sample.all.GR, 
  type = "any",
  select = "all"
)
index.promoter.w.OVERALL.whole.loop # any: 264932 // w/o capping lt2mb: 158517// ENSMBL 158846

OVERALL.loop.for.promoter.hits <- subjectHits(index.promoter.w.OVERALL.whole.loop)
OVERALL.promoter.on.loop.hits <- queryHits(index.promoter.w.OVERALL.whole.loop)

df.promoter.dist.result <- data.frame(
  loop.id = OVERALL.df.DISTINCT.loop.deep.sample.all$loop.id[OVERALL.loop.for.promoter.hits],
  loop.start = OVERALL.df.DISTINCT.loop.deep.sample.all$x0[OVERALL.loop.for.promoter.hits],
  loop.end = OVERALL.df.DISTINCT.loop.deep.sample.all$y3[OVERALL.loop.for.promoter.hits],
  loop.res = OVERALL.df.DISTINCT.loop.deep.sample.all$resolution[OVERALL.loop.for.promoter.hits],
  # loop.new.distance = OVERALL.df.DISTINCT.loop.deep.sample.all$new_distance[OVERALL.loop.for.promoter.hits],
  promoter_id = df.promoter.rn7$promoter.id[OVERALL.promoter.on.loop.hits],
  promoter_start = df.promoter.rn7$start[OVERALL.promoter.on.loop.hits],
  promoter_end = df.promoter.rn7$end[OVERALL.promoter.on.loop.hits],
  promoter_gene_id = df.promoter.rn7$gene_id[OVERALL.promoter.on.loop.hits],
  promoter_length = df.promoter.rn7$length[OVERALL.promoter.on.loop.hits]
) %>% 
  mutate(across(where(is.numeric), ~ format(., scientific = FALSE)))

df.promoter.dist.result %>% head()
df.promoter.dist.result %>% dim() # any: 264932(dedups)// w/o capping lt2mb: 158517// ENSMBL: 158846

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
relative.pos.df.promoter.dist.result %>% dim() # any: 264932(dedups) // w/o capping lt2mb: 158517 // ENSMBL: 158846
########################
# 4. promoter
# 4-3. overall distribution of promoter on loops: figures
# 4-3-1. by CHROMOSOME
########################

checking_component_distribution(relative.pos.df.promoter.dist.result, "promoter")

relative.pos.df.promoter.dist.result %>% head()

########################
# 4. promoter
# 4-3. overall distribution of promoter on loops: figures
# 4-3-2. by resolution
########################

plot.promoter.hist <- plot_histogram(relative.pos.df.promoter.dist.result)
plot.promoter.dens <- plot_density(relative.pos.df.promoter.dist.result)

saving_combined_plot(plot.promoter.hist, plot.promoter.dens,
                   "figures/submission/lt2mb/overall_distribution_of_promoter_by_resolution_wo_capping_lt2mb.pdf")

#########################################
#########################################
# Figure for new ones
#########################################
#########################################

# func9
combining_and_save_plots <- function(p1, p2, p3,
                                   filename = "histogram_combined_all.png",
                                   output_dir = "./figures/submission/lt2mb",
                                   width = 11, height = 8.5, dpi = 300) {
  # base theme
  base_theme <- theme_bw(base_size = 12) +
    theme(plot.title.position = "plot",
          panel.grid.minor = element_blank())
  
  # theme
  p1 <- p1 + base_theme
  p2 <- p2 + base_theme
  p3 <- p3 + base_theme
  
  # combining plots
  fig_combined <- (p1 | p2 | p3) +
    plot_layout(ncol = 3, guides = "collect", widths = c(1, 1, 1)) &
    theme(legend.position = "bottom",
          legend.margin = margin(2, 6, 2, 6),
          plot.tag = element_text(face = "bold", size = 12))
  
  # adding tag
  fig_combined <- fig_combined +
    plot_annotation(tag_levels = "a")
  
  # saving
  ggsave(file.path(output_dir, filename), fig_combined, width = width, height = height, dpi = dpi, bg = "white")
  
  # return(fig_combined)
}
# func9 END

combining_and_save_plots(plot.ctcf.hist, plot.tss.hist, plot.promoter.hist, "histogram_combined_all.png")
combining_and_save_plots(plot.ctcf.dens, plot.tss.dens, plot.promoter.dens, "density_combined_all.png")

########################
# 4. promoter
# 4-4. Distribution of promoter at each end in a loop: for the number of promoter used in filtering valid loops: figures
# so, the object should be used one without padding: df.DISTINCT.loop.deep.sample.all
# 4-4-1. inner distance analysis 1: checking inner distance between x3 and y0
########################
df.DISTINCT.loop.deep.sample.all %>% head()

############################################
# padding by resolution
############################################
# 1. Promoter + UPSTREAM
index.distinct.promoter.w.up.loop.each.end <- findOverlaps(
  df.promoter.rn7.GR, 
  df.DISTINCT.loop.deep.sample.all.padded.for.TSS.promoter.UP.GR, 
  type = "any",
  select = "all"
)
index.distinct.promoter.w.up.loop.each.end # any: 68983 // lt2mb default: 6710| 6case: 7573| 7case: 9691| 8case: 11776| 9case: 15656

loop.up.hits.promoter.each.end <- subjectHits(index.distinct.promoter.w.up.loop.each.end)
promoter.up.hits.each.end <- queryHits(index.distinct.promoter.w.up.loop.each.end)

df.overlapping.promoter.w.UPSTREAM.result.each.end <- data.frame(
  up.loop.id = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.promoter.UP.GR)$loop.id[loop.up.hits.promoter.each.end],
  end.up.distance = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.promoter.UP.GR)$end.distance[loop.up.hits.promoter.each.end],
  resolution = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.promoter.UP.GR)$resolution[loop.up.hits.promoter.each.end],
  promoter.id = mcols(df.promoter.rn7.GR)$promoter.id[promoter.up.hits.each.end],
  promoter.gene_id = mcols(df.promoter.rn7.GR)$gene_id[promoter.up.hits.each.end],
  WHERE = "UP"
) %>% 
  mutate(promoter.loop.up.id = str_c(up.loop.id, '|', promoter.id, '|', WHERE))

df.overlapping.promoter.w.UPSTREAM.result.each.end %>% dim() # 68983     6 (dedups) // lt2mb default: 6710| 6case: 7573| 7case: 9691| 8case: 11776| 9case: 15656
df.overlapping.promoter.w.UPSTREAM.result.each.end %>% head()

df.overlapping.promoter.w.UPSTREAM.result.each.end %>% 
  count(end.up.distance)
# any
# default                  # caes 6                    # caes 7                     # caes 8                  # caes 9                   # end.up.distance     n
#1            5000 1485    #1            5000 1813     #1            5000 2159      # 1            5000 2549  #1            5000 3238    # 1            5000 15392
#2           10000 2367    #2           10000 2902     #2           10000 3477      # 2           10000 3999  #2           10000 5129    # 2           10000 20602
#3           25000 2858    #3           25000 2858     #3           25000 4055      # 3           25000 5228  #3           25000 7289    # 3           25000 329

# 2. promoter + DOWNSTREAM (df.DISTINCT.loop.deep.sample.all.padded.for.TSS.promoter.DOWN.GR, df.promoter.rn7.GR)
index.distinct.promoter.w.down.loop.each.end <- findOverlaps(
  df.promoter.rn7.GR, 
  df.DISTINCT.loop.deep.sample.all.padded.for.TSS.promoter.DOWN.GR, 
  type = "any",
  select = "all"
)

index.distinct.promoter.w.down.loop.each.end # .4 any: 79080 // // lt2mb default: 6342| 6case: 7178| 7case: 9100| 8case: 10945| 9case: 14884
end.loop.down.promoter.hits <- subjectHits(index.distinct.promoter.w.down.loop.each.end)
end.promoter.down.hits <- queryHits(index.distinct.promoter.w.down.loop.each.end)

df.overlapping.promoter.w.DOWNSTREAM.result.each.end <- tibble(
  down.loop.id = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.promoter.DOWN.GR)$loop.id[end.loop.down.promoter.hits],
  end.down.distance = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.promoter.DOWN.GR)$end.distance[end.loop.down.promoter.hits],
  resolution = mcols(df.DISTINCT.loop.deep.sample.all.padded.for.TSS.promoter.DOWN.GR)$resolution[end.loop.down.promoter.hits],
  promoter.id = mcols(df.promoter.rn7.GR)$promoter.id[end.promoter.down.hits],
  promoter.gene_id = mcols(df.promoter.rn7.GR)$gene_id[end.promoter.down.hits],
  WHERE = "DOWN"
) %>% 
  mutate(promoter.loop.down.id = str_c(down.loop.id, '|', promoter.id, '|', WHERE))

df.overlapping.promoter.w.DOWNSTREAM.result.each.end %>% dim() # .4 any: 79080     6 (dedups) // lt2mb default: 6342| 6case: 7178| 7case: 9100| 8case: 10945| 9case: 14884
df.overlapping.promoter.w.DOWNSTREAM.result.each.end %>% head()

df.overlapping.promoter.w.DOWNSTREAM.result.each.end %>% 
  count(end.down.distance)

colnames(df.overlapping.promoter.w.DOWNSTREAM.result.each.end)
# any
# dedups
# case default                    # case6                      # case7                      # case8                       # case9                    # end.down.distance     n
#1              5000  1358        #1              5000  1668   #1              5000  1992   # 1              5000  2283   #1              5000  3013 # 1              5000 18483
#2             10000  2248        #2             10000  2774   #2             10000  3314   # 2             10000  3792   #2             10000  4806 # 2             10000 23025
#3             25000  2736        #3             25000  2736   #3             25000  3794   # 3             25000  4870   #3             25000  7065 # 3             25000 37572

########## bind_rows(UPSTREAM & DOWNSTREAM) -> BOTH
df.overlapping.promoter.w.BOTH.result <- bind_rows(df.overlapping.promoter.w.UPSTREAM.result.each.end %>% 
                                                     mutate(loop.id = up.loop.id, end.distance = end.up.distance) %>% 
                                                     mutate(case.id = promoter.loop.up.id) %>% 
                                                     mutate(promoter.gene_id = str_c(promoter.gene_id, '|', 'UP')) %>% 
                                                     dplyr::select(-c(up.loop.id, end.up.distance, promoter.loop.up.id)), 
                                                   df.overlapping.promoter.w.DOWNSTREAM.result.each.end %>% 
                                                     mutate(loop.id = down.loop.id, end.distance = end.down.distance) %>% 
                                                     mutate(case.id = promoter.loop.down.id) %>% 
                                                     mutate(promoter.gene_id = str_c(promoter.gene_id, '|', 'DOWN')) %>% 
                                                     dplyr::select(-c(down.loop.id, end.down.distance, promoter.loop.down.id))) %>% 
  mutate(promoter_result_both_id = str_c(loop.id, '|', promoter.gene_id)) %>% 
  mutate(case_sort_id = str_c(loop.id, '|', WHERE)) %>%
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

colnames(df.overlapping.promoter.w.BOTH.result)

df.overlapping.promoter.w.BOTH.result %>% dim() # any: 148063 (dedups) // lt2mb default: 13052| 6case: 14751| 7case: 18791| 8case: 22721| 9case: 30540
df.overlapping.promoter.w.BOTH.result %>% head() 
df.overlapping.promoter.w.BOTH.result %>% count(resolution)
df.overlapping.promoter.w.BOTH.result %>% dplyr::select(promoter.gene_id) %>% filter(is.na(promoter.gene_id))

df.case.from.promoter <- df.overlapping.promoter.w.BOTH.result %>% 
  group_by(loop.id) %>% 
  summarise(
    WHERE_list = list(WHERE),
    WHERE_uniq = list(base::unique(WHERE)),
    .groups = "drop"
  ) %>% 
  rowwise() %>%
  mutate(case = case_when(
    length(WHERE_uniq) == 2 ~ "BOTH",
    length(WHERE_uniq) == 1 & WHERE_uniq[[1]] == "UP" & length(unlist(WHERE_list)) == 1 ~ "u.UP",
    length(WHERE_uniq) == 1 & WHERE_uniq[[1]] == "UP" & length(unlist(WHERE_list)) > 1 ~ "m.UP",
    length(WHERE_uniq) == 1 & WHERE_uniq[[1]] == "DOWN" & length(unlist(WHERE_list)) == 1 ~ "u.DOWN",
    length(WHERE_uniq) == 1 & WHERE_uniq[[1]] == "DOWN" & length(unlist(WHERE_list)) > 1 ~ "m.DOWN",
    TRUE ~ NA_character_
  )) %>% 
  ungroup()

df.case.from.promoter %>% count(case)
# case5                # case7              # case8           # case9           # default       n
#1 BOTH    1836        #1 BOTH    2570      # 1 BOTH    3303  #1 BOTH    4744   # 1 BOTH    1540
#2 m.DOWN   583        #2 m.DOWN   757      # 2 m.DOWN   918  #2 m.DOWN  1265   # 2 m.DOWN   506
#3 m.UP     608        #3 m.UP     846      # 3 m.UP    1090  #3 m.UP    1488   # 3 m.UP     528
#4 u.DOWN  3508        #4 u.DOWN  3887      # 4 u.DOWN  4064  #4 u.DOWN  4276   # 4 u.DOWN  3269
#5 u.UP    3828        #5 u.UP    4227      # 5 u.UP    4446  #5 u.UP    4525   # 5 u.UP    3586

## debug
# df.promoter.debug <- df.overlapping.promoter.w.BOTH.result %>%
#   group_by(loop.id) %>%
#   summarise(
#     WHERE_list = list(WHERE),
#     WHERE_uniq = list(base::unique(WHERE)),
#     .groups = "drop"
#   )

# View(df.promoter.debug)

# checking: loop.id case exclusive 
check.exclusive.promoter <- df.case.from.promoter %>%
  group_by(loop.id) %>%
  summarise(
    n_case = n_distinct(case),
    case_list = paste(unique(case), collapse = ","),
    .groups = "drop"
  ) %>%
  filter(n_case > 1)  # loop.id involved in more than 1 case

check.exclusive.promoter # 0

df.final.loop.dataset.promoter <- df.overlapping.promoter.w.BOTH.result %>% left_join(df.case.from.promoter, by = 'loop.id')

df.final.loop.dataset.promoter %>% head(3)
df.final.loop.dataset.promoter %>% dim() # 22721    13
df.final.loop.dataset.promoter %>% count(case)


# case6               # case7           # case8           # case9           # default
#    case    n       # case    n        # case    n       # case    n       # case    n
#1   BOTH 4865       #1   BOTH 7177     # 1   BOTH 9753   # 1   BOTH 15430  # 1   BOTH 4000
#2 m.DOWN 1249       #2 m.DOWN 1668     # 2 m.DOWN 2059   # 2 m.DOWN  2943  # 2 m.DOWN 1071
#3   m.UP 1301       #3   m.UP 1832     # 3   m.UP 2399   # 3   m.UP  3366  # 3   m.UP 1126
#4 u.DOWN 3508       #4 u.DOWN 3887     # 4 u.DOWN 4064   # 4 u.DOWN  4276  # 4 u.DOWN 3269
#5   u.UP 3828       #5   u.UP 4227     # 5   u.UP 4446   # 5   u.UP  4525  # 5   u.UP 3586

###############
###############
# for Q1: quartile among loops with TSS, so min = 1
df.promoter.counts <- df.overlapping.promoter.w.BOTH.result %>%
  group_by(loop.id, WHERE, resolution) %>%
  summarise(promoter_count = n_distinct(promoter.id), .groups = 'drop')

df.promoter.counts %>% count(promoter_count)

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

df.promoter.case
df.promoter.case %>% count(case)

df.loop.with.promoter.case <- df.DISTINCT.loop.deep.sample.all.lt.2mb %>% # 31773
  left_join(df.promoter.case, by = "loop.id") %>%
  mutate(case = ifelse(is.na(case), "NONE", case))

df.loop.with.promoter.case %>% count(case)

df.promoter.counts

##########################################################
##########################################################
# checking hits of TSS & promoter in an anchor using distanceToNearest
##########################################################
##########################################################

df.tss.ucsc %>% head(3) # chr  start    end        gene_id                          tss.id
df.promoter.rn7 %>% head(3) # chr     start     end gene_id length  center promoter.id 
df.gene_tss_and_pro <- bind_rows(df.tss.ucsc %>% 
                                  dplyr::select(chr, start, end, gene_id, tss.id) %>%
                                  mutate(across(c(start, end), as.numeric)) %>%
                                  mutate(component = "tss") %>%
                                  dplyr::rename(component_id = tss.id) %>%
                                  mutate(component_id = str_c(component_id, component, sep='|')),
                                df.promoter.rn7 %>% 
                                  dplyr::select(chr, start, end, gene_id, promoter.id) %>%
                                  mutate(across(c(start, end), as.numeric)) %>%
                                  mutate(component = "pro") %>%
                                  dplyr::rename(component_id = promoter.id) %>%
                                  mutate(component_id = str_c(component_id, component, sep='|')))

df.gene_tss_and_pro %>% head(3)
df.gene_tss_and_pro %>% dim() # 29507  // 34183

df.gene_tss_and_pro.GR <- GRanges(seqnames=df.gene_tss_and_pro$chr,
                                  ranges=IRanges(start=df.gene_tss_and_pro$start, end=df.gene_tss_and_pro$end),
                                  gene_id=df.gene_tss_and_pro$gene_id,
                                  component_id=df.gene_tss_and_pro$component_id,
                                  component=df.gene_tss_and_pro$component)

df.gene_tss_and_pro.GR # 29507// 34183

df.DISTINCT.loop.deep.sample.all.lt.2mb %>% dim() # 31019
df.DISTINCT.loop.deep.sample.all.lt.2mb %>% head(3)

df.DISTINCT.loop.deep.sample.all.lt.2mb.prep <- df.DISTINCT.loop.deep.sample.all.lt.2mb %>% 
  mutate(mid_x = floor((x0 + x3) / 2), mid_y = floor((y0 + y3) / 2))

df.DISTINCT.loop.deep.sample.all.lt.2mb.prep.UP.GR <- GRanges(
  seqnames = df.DISTINCT.loop.deep.sample.all.lt.2mb.prep$chr1,
  # ranges = IRanges(start = df.DISTINCT.loop.deep.sample.all.lt.2mb.prep$mid_x, end = df.DISTINCT.loop.deep.sample.all.lt.2mb.prep$mid_x + 1),
  ranges = IRanges(start = df.DISTINCT.loop.deep.sample.all.lt.2mb.prep$x1, end = df.DISTINCT.loop.deep.sample.all.lt.2mb.prep$x2),
  strand = "*",
  loop.id = df.DISTINCT.loop.deep.sample.all.lt.2mb.prep$loop.id,
  distance = df.DISTINCT.loop.deep.sample.all.lt.2mb.prep$distance,
  resolution = df.DISTINCT.loop.deep.sample.all.lt.2mb.prep$resolution
)
df.DISTINCT.loop.deep.sample.all.lt.2mb.prep.DOWN.GR <- GRanges(
  seqnames = df.DISTINCT.loop.deep.sample.all.lt.2mb.prep$chr1,
  # ranges = IRanges(start = df.DISTINCT.loop.deep.sample.all.lt.2mb.prep$mid_y, end = df.DISTINCT.loop.deep.sample.all.lt.2mb.prep$mid_y + 1),
  ranges = IRanges(start = df.DISTINCT.loop.deep.sample.all.lt.2mb.prep$y1, end = df.DISTINCT.loop.deep.sample.all.lt.2mb.prep$y2),
  strand = "*",
  loop.id = df.DISTINCT.loop.deep.sample.all.lt.2mb.prep$loop.id,
  distance = df.DISTINCT.loop.deep.sample.all.lt.2mb.prep$distance,
  resolution = df.DISTINCT.loop.deep.sample.all.lt.2mb.prep$resolution
)

hits_up <- distanceToNearest(df.DISTINCT.loop.deep.sample.all.lt.2mb.prep.UP.GR, df.gene_tss_and_pro.GR)

df_up <- data.frame(
  queryHits = queryHits(hits_up),
  subjectHits = subjectHits(hits_up),
  distance = mcols(hits_up)$distance
) %>%
  mutate(
    loop.id = mcols(df.DISTINCT.loop.deep.sample.all.lt.2mb.prep.UP.GR)[queryHits, "loop.id"],
    resolution = mcols(df.DISTINCT.loop.deep.sample.all.lt.2mb.prep.UP.GR)[queryHits, "resolution"],
    gene_id = mcols(df.gene_tss_and_pro.GR)[subjectHits, "gene_id"],
    component = mcols(df.gene_tss_and_pro.GR)[subjectHits, "component"],
    WHERE = "UP"
  )

hits_down <- distanceToNearest(df.DISTINCT.loop.deep.sample.all.lt.2mb.prep.DOWN.GR, df.gene_tss_and_pro.GR)

df_down <- data.frame(
  queryHits = queryHits(hits_down),
  subjectHits = subjectHits(hits_down),
  distance = mcols(hits_down)$distance
) %>%
  mutate(
    loop.id = mcols(df.DISTINCT.loop.deep.sample.all.lt.2mb.prep.DOWN.GR)[queryHits, "loop.id"],
    resolution = mcols(df.DISTINCT.loop.deep.sample.all.lt.2mb.prep.DOWN.GR)[queryHits, "resolution"],
    gene_id = mcols(df.gene_tss_and_pro.GR)[subjectHits, "gene_id"],
    component = mcols(df.gene_tss_and_pro.GR)[subjectHits, "component"],
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
mean(df.final.up.down.tss.pro.nearest$distance)

ggplot(df.final.up.down.tss.pro.nearest, aes(y = distance)) +
  geom_boxplot(fill = "#A6CEE3", color = "#1F78B4", outlier.color = "red", outlier.shape = 16) +
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

# func.
approach_2nd_analyze_loops_by_threshold <- function(df, threshold_distance = 2e5, top_n_genes = 70, print_top_n = 50) {

  # Task 1: Filtering and top N genes by loop count
  df.filtered <- df %>% filter(distance < threshold_distance)
# dim(df.filtered) # 31019
  df.gene_loop_count <- df.filtered %>%
    count(gene_id, sort = TRUE)

  df.top_genes <- df.gene_loop_count %>%
    slice_max(n, n = top_n_genes) %>%
    arrange(n)

  plot <- ggplot(df.top_genes, aes(x = reorder(gene_id, n), y = n)) +
    geom_bar(stat = "identity", fill = "#1F78B4") +
    geom_text(aes(label = n), vjust = -0.3, color = "red", size = 3) +
    coord_flip() +
    labs(
      title = paste("Top", top_n_genes, "Genes by Number of Loops (distance <", threshold_distance, ")"),
      x = "Gene",
      y = "Number of Loops"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

  print(plot)

  # Task 2: Percentage of loops retained
  total_loops <- df %>% distinct(loop.id)
  filtered_loops <- df.filtered %>% distinct(loop.id)

  percent_retained <- round((nrow(filtered_loops) / nrow(total_loops)) * 100, 2)
  message("total loops: ", nrow(total_loops), " ")
  message("filtered_loops: ", nrow(filtered_loops), " ")
  message("Filtered loops retain ", percent_retained, "% of total loops")

  # Task 3: Print top genes
  top_genes <- df.gene_loop_count %>%
    slice_max(n, n = print_top_n) %>%
    pull(gene_id)

  cat("Top", print_top_n, "genes:\n")
  cat(top_genes, sep = "\n")
}

mean(df.final.up.down.tss.pro.nearest$distance) # 109537.7// midpoint, ENSEMBL61564.34// span, ENSEMBL 55191.25
threshold_distance <- approach2.stats$Q3 # 82165.75// midpoint, ENSEMBL 62222.75//  span, ENSEMBL 54824
threshold_distance <- approach2.stats$Median # 28943.5// midpoint, ENSEMBL 23377.5//  span, ENSEMBL 16545.5
threshold_distance <- approach2.stats$Q1 # 8807// midpoint, ENSEMBL 7422 //  span, ENSEMBL 514.25
threshold_distance
approach_2nd_analyze_loops_by_threshold(df.final.up.down.tss.pro.nearest, threshold_distance = threshold_distance, top_n_genes = 70, print_top_n = 50)
# ENSEMBL point: https://biit.cs.ut.ee/gplink/l/aEqIdfCD3TD
# ENSEMBL span: https://biit.cs.ut.ee/gplink/l/aDqh-aFfbQP
df.final.up.down.tss.pro.nearest %>% head()

df.unique_loops <- df.final.up.down.tss.pro.nearest %>%
  # filter(
  #   case_when(
  #     resolution == "5K"  ~ distance <= 5000*1.5,
  #     resolution == "10K" ~ distance <= 10000*1.5,
  #     resolution == "25K" ~ distance <= 25000*1.5,
  #     TRUE                ~ FALSE  
  #   )
  # ) %>% 
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
final.loops.from.tss.step %>% dim() # point ENSEMBL:13407//span ENSEMBL: 15748

final.loops.from.promoter.step <- df.unique_loops %>%
  filter(component == "pro")
final.loops.from.promoter.step %>% dim() # point ENSEMBL: 10117// span ENSEMBL: 13209

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
##########################################################
# UP-all
df.final.loop.dataset.tss.unique <- df.final.loop.dataset.tss %>% 
  filter(str_starts(case, "u.")) %>% 
  dplyr::select(loop.id, case.id, WHERE, tss.gene_id) %>% 
  mutate(case.id = str_c(case.id, "TSS", sep = '|')) %>% 
  mutate(gene.id = str_split_n(tss.gene_id, '\\|', 1)) %>%
  dplyr::rename(case.where.id = tss.gene_id)

# view()
# colnames() %>%
# head()
df.final.loop.dataset.promoter.unique <- df.final.loop.dataset.promoter %>%
  filter(str_starts(case, "u.")) %>% 
  dplyr::select(loop.id, case.id, WHERE, promoter.gene_id) %>% 
  mutate(case.id = str_c(case.id, "promoter", sep = '|')) %>% 
  mutate(gene.id = str_split_n(promoter.gene_id, '\\|', 1)) %>%
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

saving_plot_dual(
  plot_obj = the_number_of_genes_per_loops,
  output_dir = "./figures/submission/lt2mb",
  filename_base = "the_number_of_genes_per_loops")
# case8
# 13177 loops with 1 gene
# 1780 loops with 2 genes

###################
# UP-unique
# df.final.loop.dataset.tss %>% dim() # 16006    13
# df.final.loop.dataset.promoter %>% dim() # 13052    13

# df.final.loop.dataset.tss %>% head(3)
# df.final.loop.dataset.promoter %>% head(3)

# df.final.loop.dataset.tss %>% colnames()
# df.final.loop.dataset.promoter %>% colnames()

# df.uUP.final.loop.tss <- df.final.loop.dataset.tss %>%
#   filter(case == "u.UP") %>% # 4478
#   mutate(gene_id = str_c('tss|', tss.gene_id)) %>% 
#   dplyr::select(loop.id, gene_id, case)

# df.uUP.final.loop.promoter <- df.final.loop.dataset.promoter %>% 
#   filter(case == "u.UP") %>% # 3586
#   mutate(gene_id = str_c('promoter|', promoter.gene_id)) %>% 
#   dplyr::select(loop.id, gene_id, case)

# # full join
# df.uUP.final.tss.promoter <- full_join(df.uUP.final.loop.tss, df.uUP.final.loop.promoter, by = 'loop.id') %>% 
#   mutate(tss_gene = str_split_n(gene_id.x, '\\|', 2), promoter_gene = str_split_n(gene_id.y, '\\|', 2))

# # DOWN
# df.uDOWN.final.loop.tss <- df.final.loop.dataset.tss %>%
#   filter(case == "u.DOWN") %>% # 4368
#   mutate(gene_id = str_c('tss|', tss.gene_id)) %>% 
#   dplyr::select(loop.id, gene_id, case)


# df.uDOWN.final.loop.promoter <- df.final.loop.dataset.promoter %>% 
#   filter(case == "u.DOWN") %>% # 3,269
#   mutate(gene_id = str_c('promoter|', promoter.gene_id)) %>% 
#   dplyr::select(loop.id, gene_id, case)

# # full join
# df.uDOWN.final.tss.promoter <- df.uDOWN.final.loop.tss %>% 
#   full_join(df.uDOWN.final.loop.promoter, by = 'loop.id') %>% 
#   mutate(tss_gene = str_split_n(gene_id.x, '\\|', 2), promoter_gene = str_split_n(gene_id.y, '\\|', 2))

# ##########################################
# df.uUP.final.tss.promoter # 5752
# df.uDOWN.final.tss.promoter # 5523

# df.unique.final <- bind_rows(df.uUP.final.tss.promoter, df.uDOWN.final.tss.promoter)
# all.genes.from.u.final <- c(df.unique.final$tss_gene, df.unique.final$promoter_gene)

# df.unique.final
# all.genes.from.u.final

# # deleteing NA & dups
# df.all.genes.from.u.final.cleaned <- data.frame(gene = unique(na.omit(all.genes.from.u.final)))

# df.all.genes.from.u.final.cleaned # 4859

# df.gene.loop.map <- df.unique.final %>% 
#   dplyr::select(loop.id, tss_gene, promoter_gene) %>% # view()
#   pivot_longer(cols = c(tss_gene, promoter_gene), values_to = "gene") %>% # view()
#   filter(!is.na(gene)) %>%
#   group_by(gene) %>% # view()
#   summarise(loop_ids = list(unique(loop.id)),
#             n_loops = length(unique(loop.id)),
#             .groups = "drop")

# df.gene.loop.map %>% head(3)

# df.gene.loop.map.for.export <- df.gene.loop.map %>%
#   mutate(loop_ids = sapply(loop_ids, function(x) paste(x, collapse = ",")))

# df.gene.loop.map.for.export %>% head(3)

# ggplot(df.gene.loop.map.for.export, aes(x = n_loops)) +
#   geom_histogram(binwidth = 1, fill = "#1f78b4", color = "white") +
#   stat_bin(binwidth = 1, geom = "text", aes(label = ..count..), vjust = -0.5, color = "black") +
#   labs(
#     title = "Distribution of Loop Counts per Gene",
#     x = "Number of Loops (n_loops)",
#     y = "Number of Genes"
#   ) +
#   theme(
#     axis.text.x = element_text(angle = 0, hjust = 1),
#     plot.title = element_text(hjust = 0.5),
#     legend.position = "bottom"
#   )


# # TSV
# library("vroom")
# write_tsv(df.gene.loop.map.for.export, "./gene_loop_map.tsv")

# getwd()

########## 
# Step 1: tss_gene & promoter_gene: long-format
# df.long.unique.final <- df.unique.final %>%
#   dplyr::select(loop.id, tss_gene, promoter_gene) %>%
#   pivot_longer(cols = c(tss_gene, promoter_gene), names_to = "source", values_to = "gene") %>%
#   filter(!is.na(gene))  # NA deletion

# # Step 2: genes per loop
# genes_per_loop <- df.long.unique.final %>%
#   group_by(loop.id) %>%
#   summarise(n_genes = n_distinct(gene), .groups = "drop")

# # Step 3: ggplot
# plot_data <- genes_per_loop %>%
#   count(n_genes)

# ggplot(plot_data, aes(x = n_genes, y = n)) +
#   geom_bar(stat = "identity", fill = "steelblue") +
#   geom_text(aes(label = n), vjust = -0.3, size = 3.5) +
#   scale_x_continuous(breaks = seq(min(plot_data$n_genes), max(plot_data$n_genes), by = 1)) +
#   labs(title = "Distribution of Number of Genes per Loop",
#        x = "Number of Unique Genes per Loop",
#        y = "Count of Loops") +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 0, hjust = 1),
#     plot.title = element_text(hjust = 0.5),
#     legend.position = "none"
#   )

# top.genes.50 <- df.long.unique.final %>%
#   group_by(gene) %>%
#   summarise(n_loops = n_distinct(loop.id)) %>%
#   arrange(desc(n_loops)) %>%
#   slice_head(n = 50)

# writeLines(top.genes.50$gene)

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

process_feature_bins <- function(feature_df, chromosome_ends, bin_size = 1e6, label = "Feature") {
  options(scipen = 999)
  
  # col type
  chromosome_ends <- chromosome_ends %>% mutate(chr = str_remove(as.character(chr), 'chr'))
  feature_df$chr <- as.character(feature_df$chr)
  
  result <- list()
  
  valid_chromosomes <- as.character(c(1:20, "X", "Y"))
  feature_df <- feature_df %>% mutate(chr = str_remove(chr, 'chr')) %>% filter(chr %in% valid_chromosomes)
  
  for (current_chr in unique(feature_df$chr)) {
    chr_data <- feature_df %>% filter(chr == current_chr)
    chr_end <- chromosome_ends %>% filter(chr == current_chr) %>% pull(end)
    
    if (length(chr_end) != 1 || !is.finite(chr_end) || chr_end < bin_size) {
      warning(paste("No valid end position for", current_chr, "- skipping"))
      next
    }
    
    # interval: Start = End+1
    bin_edges <- seq(1, chr_end, by = bin_size)
    if (length(bin_edges) < 2) {
      warning(paste("Not enough bin edges for", current_chr, "- skipping"))
      next
    }
    bins <- data.frame(
      start = bin_edges[-length(bin_edges)],
      end = bin_edges[-1] -1
    )
    
    if (tail(bins$end, 1) < chr_end) {
      bins <- rbind(bins, data.frame(start = bins$end[nrow(bins)] + 1, end = chr_end))
    }
    
    bins <- bins %>%
      mutate(chr = current_chr,
             Value = 0,
             Feature = label)
    
    for (i in 1:nrow(bins)) {
      bins$Value[i] <- sum(
        chr_data$start >= bins$start[i] & chr_data$start <= bins$end[i],
        na.rm = TRUE
      )
    }
    
    result[[current_chr]] <- bins
  }
  
  final <- bind_rows(result)
  final <- final %>%
    dplyr::select(chr, start, end, Value, Feature)
  
  return(final)
}

getwd()
ctcf_density <- process_feature_bins(df_ctcf_ideogram, chromosome_data, bin_size, label = "CTCF") %>% dplyr::select(-last_col()) # color = "#E41A1C"
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
# func10. filtering loops by threshold
plotting_and_filtering_summary <- function(df_counts, count_col, output_prefix) {

  # 1. hitogram
  p.hist <- ggplot(df_counts, aes(x = log2(.data[[count_col]]))) +
    geom_histogram() +
    facet_wrap(~WHERE + resolution, scales = "free_y") +
    labs(
      y = "Count"
    )

  #   geom_histogram() +
  #   facet_wrap(~WHERE + resolution, scales = "free_y")

  # 2. PDF
  pdf_path <- paste0("./figures/submission/lt2mb/histogram_number_of_", output_prefix, "_within_ends_of_loops_by_resolution_end.pdf")
  pdf(file = pdf_path, width = 11 * 0.8, height = 8.5 * 0.8)
  print(p.hist)
  dev.off()

  # 3. png
  png_path <- paste0("./figures/submission/lt2mb/histogram_number_of_", output_prefix, "_within_ends_of_loops_by_resolution_end.png")
  ggsave(filename = png_path, plot = p.hist, width = 11 * 0.8, height = 8.5 * 0.8, dpi = 300)
}
# func10 END

# 4-1. CTCF
df.ctcf.counts
df.ctcf.counts %>%
  group_by(resolution, WHERE) %>%
  summarise(max_ctcf = max(ctcf_count), .groups = "drop")


df.ctcf.filtered <- plotting_and_filtering_summary(
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
final.loops.from.ctcf.step # 18,327 + 10 | sub.4 any, threshold > 6 : 25,620| lt2mb: 25,268
final.loops.from.ctcf.step %>% count(resolution)
# sub.4
# any: 
# new threshold
# lt2mb              # resolution     n
# 1 10K         9375 # 1 10K         9420
# 2 25K        11750 # 2 25K        12049
# 3 5K          4143 # 3 5K          4151

# 4-2. TSS
df.overlapping.TSS.w.BOTH.result %>% head()
df.tss.counts.each.end # .4 any: 44564| lt2mb: 14,006

df.tss.filtered <- plotting_and_filtering_summary(
  df_counts = df.tss.counts.each.end,
  count_col = "tss_count_each_end",
  output_prefix = "tss"
)

loops.with.tss.above.threshold <- df.tss.counts.each.end %>% 
  left_join(tss.stats.by.resolution.each.end %>% 
              dplyr::select(resolution, WHERE, Q1), by = c("resolution", "WHERE")) %>% 
  filter(tss_count_each_end >= Q1)

loops.with.tss.above.threshold # 44564

loops.with.tss.both.ends <- loops.with.tss.above.threshold %>%
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

final.loops.from.tss.step <- df.final.loop.dataset.tss %>% 
  filter(str_starts(case, "u"))

final.loops.from.tss.step
final.loops.from.tss.step %>% head()
final.loops.from.tss.step %>% dim()# any: 27,518 (dedups) w/ either // 11,993| lt2mb+up default: 8846, case6: 9396, case7: 10200, case8: 10447, case9: 10403
final.loops.from.tss.step %>% count(resolution)

# default                # case6                # case7                # case8                 # case9
#1 5K          1896      #1 5K          2066    #1 5K          2214    # 1 5K          2191    #1 5K          2140
#2 10K         3283      #2 10K         3663    #2 10K         3887    # 2 10K         4029    #2 10K         4103
#3 25K         3667      #3 25K         3667    #3 25K         4099    # 3 25K         4227    #3 25K         4160

# 4-3. promoter
df.overlapping.promoter.w.BOTH.result
df.promoter.counts.each.end # 38,174
promoter.stats.by.resolution.each.end

df.promoter.filtered <- plotting_and_filtering_summary(
  df_counts = df.promoter.counts.each.end,
  count_col = "promoter_count_each_end",
  output_prefix = "promoter"
)

threshold.promoter = 1

loops.with.promoter.above.threshold <- df.promoter.counts.each.end %>% 
  left_join(promoter.stats.by.resolution.each.end %>% 
              dplyr::select(resolution, WHERE, Q1), by = c("resolution", "WHERE")) %>% 
  filter(promoter_count_each_end >= threshold.promoter)

loops.with.promoter.above.threshold

loops.with.promoter.both.ends <- loops.with.promoter.above.threshold %>% # 38,152 + 10
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

final.loops.from.promoter.step <- df.final.loop.dataset.promoter %>% 
  filter(str_starts(case, "u"))

final.loops.from.promoter.step
final.loops.from.promoter.step %>% head()
final.loops.from.promoter.step %>% dim()
final.loops.from.promoter.step %>% count(resolution)

# default               # case6                # case7              # case8               # case9
# 1         5K 1522     #1         5K 1660     #1         5K 1800   # 1         5K 1836   # 1         5K 1892
# 2        10K 2531     #2        10K 2874     #2        10K 3113   # 2        10K 3292   # 2        10K 3495
# 3        25K 2802     #3        25K 2802     #3        25K 3201   # 3        25K 3382   # 3        25K 3414

####################################
####################################
# Venn Diagram
####################################
####################################

final.loops.from.promoter.step$loop.id
final.loops.from.tss.step$loop.id

# func11. venn plot
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
  
  venn_plot <- venn_plot +
    annotate("text", x = -1.5, y = 1.5, label = paste("CTCF:", ctcf_count), color = "#f8766d") +
    annotate("text", x = 1.55, y = 1.5, label = paste("Promoter:", promoter_count), color = "#629bfe") +
    annotate("text", x = 0, y = -1.7, label = paste("TSS:", tss_count), color = "#32ba36")
  
  
  return(venn_plot)
}
# func11. END

venn_plot_submission <- create_venn_plot(final.loops.from.ctcf.step, final.loops.from.promoter.step, final.loops.from.tss.step, "CTCF")
venn_plot_submission

saving_plot_dual(
  plot_obj = venn_plot_submission,
  filename_base = "ctcf_vs_promoter_tss_venn_diagrams_approach2_lt_2mb",
  output_dir = "./figures/submission/lt2mb",
  scale_x = 1,
  scale_y = 1)
####################################
# functional loop extraction
####################################

# func12. extracting common loops
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
# func12. END

final.loops.from.ctcf.step %>% dim()# sub.4 any, either 6: 25,620// span ENSEMBL 25268
final.loops.from.tss.step %>% dim() # sub.4 dedups any, either: 27,518// 13407// span ENSEMBL  17810
final.loops.from.promoter.step %>% dim() # sub.4 dedups any, either: 25,156// span ENSEMBL  13209

overlapping_loops <- extract_overlapping_loops(
  ctcf_data = final.loops.from.ctcf.step,
  promoter_data = final.loops.from.promoter.step,
  tss_data = final.loops.from.tss.step
)

overlapping_loops$ctcf_promoter_only # PASS// 608 // 5497      // 6240
overlapping_loops$ctcf_tss_only # PASS// 2519 // 7799 // 9015
overlapping_loops$ctcf_promoter_tss_overlap # PASS// 20513 // 2302 // 3426

df.ctcf.promoter.only.loop <- data.frame(loop.id = overlapping_loops$ctcf_promoter_only, category = "CP", stringsAsFactors = FALSE)
df.ctcf.tss.only.loop <- data.frame(loop.id = overlapping_loops$ctcf_tss_only, category = "CT", stringsAsFactors = FALSE)
df.ctcf.promoter.tss.overlap.loop <- data.frame(loop.id = overlapping_loops$ctcf_promoter_tss_overlap, category = "CPT", stringsAsFactors = FALSE)

df.ctcf.promoter.only.loop %>% dim() # 5497
df.ctcf.tss.only.loop %>% dim() # 7799
df.ctcf.promoter.tss.overlap.loop %>% dim() # 2302

# sub.4 dedups: 11891 = 679 + 3096 + 8116
df.final.loop <- bind_rows(df.ctcf.promoter.only.loop, df.ctcf.tss.only.loop, df.ctcf.promoter.tss.overlap.loop)

df.final.loop # 23640/31773(0.7440279) sub.4 any, either 6// 11526/31019(lt2mb) (0.3627608)
df.final.loop %>% head()
df.final.loop %>% dim() # 15598: 5656 + 7980 + 1962

df.final.loop %>% mutate(resolution = str_split_n(loop.id, '_', 7)) %>% count(resolution)
# resolution       n  resolution    n
# 1      10000  86751      10000 5656
# 2      25000 111442      25000 7980
# 3       5000  38213       5000 1962

save(df.final.loop, file="./figures/submission/lt2mb/df_final_loop_sub.4.any.lt2mb.rda")
write.csv(df.final.loop, file = "./figures/submission/lt2mb/df_final_loop_sub.4.dedup.any.lt2mb.csv", row.names = FALSE)

df.final.loop # 23640// 18681
df.final.loop %>% head()
df.DISTINCT.loop.deep.sample.all %>% head()

df.final.DISTINCT.loop.joined <- df.final.loop %>% 
  inner_join(df.DISTINCT.loop.deep.sample.all, by = "loop.id")
df.final.DISTINCT.loop.joined %>% dim()
df.final.DISTINCT.loop.joined


df.final.loop %>% distinct(loop.id) # 23460// 18681

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
pdf(file = "./figures/submission/lt2mb/circos_loops_by_resolution_either.6.lt2mb.pdf", height = 11*0.8, width = 8.5*0.8)

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

# 1. Load first two pages from PDF
pdf_path <- "./figures/submission/lt2mb/circos_loops_by_resolution_either.6.lt2mb.pdf"
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
saving_plot_dual(
  output_dir = "./figures/submission/lt2mb",
  filename_base = "circos_first_two_panels_F8",
  plot_obj = combined_plot,
  scale_x = 1,
  scale_y = 1
)

