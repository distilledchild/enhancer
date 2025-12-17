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
options(tibble.max_extra_cols = Inf)
options(scipen = 999)

getwd()

# Linux
# setwd('C:\\Users\\panju\\Dropbox (UTHSC GGI)\\Gateway_to_Hao\\workshop\\2023_NIH_meeting\\loop_N_tss')
# setwd('./Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss')
setwd("/home/pkim/dropbox/Gateway_to_Hao/enhancer/r_files")
getwd()

# Mac
setwd("~/dropbox/Gateway_to_Hao/enhancer/r_files")
getwd()
source("~/Desktop/playground/enhancer/r_files/utils_functions.R") # Load all utility functions

####################################
# Sequencing stats ###### Figure 1.a
####################################
seq.data <- read.table("../data/library_complexity.tsv", header = TRUE, sep = "\t")
seq.data[, -1] <- lapply(seq.data[, -1], function(x) as.numeric(as.character(x)))
colnames.of.seq.data <- colnames(seq.data)
colnames.of.seq.data
seq.data <- seq.data %>% mutate(
    Duplicates = PCR_Duplicates + Optical_Duplicates,
    Chimeric_ambiguous_and_Unmmapped = Chimeric_Ambiguous + Unmapped
)
seq.data %>% head(15)
seq.data %>% dplyr::select(Strain)

seq.data$Unique_Reads_Percentage <- (seq.data$Unique_Reads / seq.data$Sequenced_RP) * 100
seq.data$Duplicates_Percentage <- (seq.data$Duplicates / seq.data$Sequenced_RP) * 100
seq.data$Chimeric_ambiguous_and_Unmapped_Percentage <- (seq.data$Chimeric_ambiguous_and_Unmmapped / seq.data$Sequenced_RP) * 100

seq.data_melted <- melt(seq.data,
    id.vars = "Strain",
    measure.vars = c(
        "Unique_Reads_Percentage",
        "Duplicates_Percentage",
        "Chimeric_ambiguous_and_Unmapped_Percentage"
    ),
    variable.name = "Category", value.name = "Percentage"
)

# 4. Category order
seq.data_melted$Category <- factor(seq.data_melted$Category,
    levels = c(
        "Unique_Reads_Percentage",
        "Duplicates_Percentage",
        "Chimeric_ambiguous_and_Unmapped_Percentage"
    ),
    labels = c("Unique Reads", "Duplicates", "Chimeric Ambiguous + Unmapped")
)
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
    theme(
        axis.text.x = element_text(angle = 0, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom"
    ) +
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

# CACHING for df.loop.deep.sample.all
cache_file_loop_deep_sample_all <- "../data/df.loop.deep.sample.all.rds"

if (file.exists(cache_file_loop_deep_sample_all)) {
    message("Loading cached loop data from: ", cache_file_loop_deep_sample_all)
    df.loop.deep.sample.all <- readRDS(cache_file_loop_deep_sample_all)
} else {
    message("Processing and caching loop data...")

    # BEDPE file list (10 files)
    # Linux
    # loop.file.list = fs::dir_ls("/home/pkim/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/", regexp = ".bedpe$")
    # loop.file.list = fs::dir_ls("/home/pkim/Desktop/temp/enhancer/dropbox_enhancer_doosan/data/loops/", regexp = ".bedpe$")
    # loop.file.list = fs::dir_ls("/home/hao/Dropbox/Gateway_to_Hao/enhancer/data/loops", regexp = ".bedpe$")

    # Mac
    # loop.file.list = fs::dir_ls("~/dropbox/K\ P/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/", regexp = ".bedpe$")
    # loop.file.list = fs::dir_ls("/Users/PanjunKim/Desktop/temp/enhancer/dropbox_enhancer_doosan/data/loops/", regexp = ".bedpe$")
    # loop.file.list = fs::dir_ls("~/Desktop/temp/enhancer/dropbox_enhancer_doosan/data/loops/", regexp = ".bedpe$")

    # Detect OS and set path
    if (Sys.info()["sysname"] == "Linux") {
        loop.file.list <- fs::dir_ls("/home/hao/Dropbox/Gateway_to_Hao/enhancer/data/loops", regexp = ".bedpe$")
    } else {
        loop.file.list <- fs::dir_ls("~/dropbox/Gateway_to_Hao/enhancer/data/loops", regexp = ".bedpe$")
    }

    print(loop.file.list)

    # BED fild for loops
    df.init.loop.bed <- init.bedpe.df(loop.file.list, "loops") %>% # utils_functions.R
        filter(!str_detect(X.chr1, "^#"))
    # %>% #58,992
    #   view()

    print(count(df.init.loop.bed, sample))

    # creating loop.id, sample.loop.id
    df.loop.deep.sample.all <- df.init.loop.bed %>%
        mutate(strain = case_when(
            sample == "592BB" ~ "SHR/OlaIpcv",
            sample == "607" ~ "HXB10",
            sample == "74AA" ~ "F344/Stm",
            sample == "A2DB" ~ "LE/Stm",
            sample == "D765A" ~ "BXH6",
            sample == "DE8BA" ~ "BN-Lx",
            sample == "DBA9A" ~ "HXB23",
            sample == "DA21A" ~ "SHR/OlaIpcvxBN/NHsdMcwi",
            sample == "DA08A" ~ "HXB2",
            sample == "DA68A" ~ "HXB31",
            TRUE ~ NA
        )) %>%
        mutate(end.distance = x2 - x1) %>%
        mutate(resolution = convert_to_resolution(end.distance)) %>% # Using utility function # utils_functions.R
        mutate(loop.id = str_c(X.chr1, "_", x1, "_", x2, "_", chr2, "_", y1, "_", y2, "_", end.distance)) %>% # loop.id
        mutate(sample.loop.id = str_c(strain, "_", loop.id)) %>% # sample.loop.id
        dplyr::select(sample, strain, X.chr1, x1, x2, chr2, y1, y2, distance, end.distance, resolution, loop.id, sample.loop.id) %>%
        dplyr::rename(chr1 = X.chr1)

    saveRDS(df.loop.deep.sample.all, cache_file_loop_deep_sample_all)
}

df.loop.deep.sample.all %>% dim() # 51300 + 7692 = 58992
df.loop.deep.sample.all %>% head(3)
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

for (res in resolutions) {
    message("Start processing: ", res)
    df.loop.deep.sample.all.res <- df.loop.deep.sample.all %>% filter(resolution == res)

    loop_counts <- df.loop.deep.sample.all.res %>%
        group_by(strain) %>%
        summarise(total_loops = n_distinct(loop.id))

    location_list <- split(df.loop.deep.sample.all.res$loop.id, df.loop.deep.sample.all.res$strain)

    message("location_list created with ", length(location_list), " strains.")

    common_counts <- matrix(0, nrow = length(location_list), ncol = length(location_list))
    rownames(common_counts) <- colnames(common_counts) <- names(location_list)

    for (i in 1:length(location_list)) {
        for (j in 1:length(location_list)) {
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
        scale_fill_gradient(
            low = "white", high = "darkred",
            name = "Common Loop %", limits = c(0, 100)
        ) +
        theme_minimal() +
        labs(x = "Strain", y = "Strain", title = paste("Heatmap of Common Loop Percentage -", res, "Resolution")) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5)
        ) # align title in the center

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
    width_in = 11 * 1.6,
    height_in = 8.5 * 3 * 1.6,
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
    filter(n_samples > 1) %>% # loops shared by more than 2 sample
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
    geom_errorbar(
        aes(
            ymin = mean_shared_loops - sd_shared_loops,
            ymax = mean_shared_loops + sd_shared_loops
        ),
        width = 0.2, color = "black"
    ) +
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
    mutate(
        from = map_chr(pairs, 1),
        to = map_chr(pairs, 2)
    ) %>%
    count(from, to, name = "weight") # number of shared loops

strain_pairs

# Convert to igraph
network_graph <- graph_from_data_frame(strain_pairs, directed = FALSE)

# Plot
network_plot_for_shared_loops <- ggraph(network_graph, layout = "fr") + # fr: force-directed layout
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
        "Alignable Reads" = Alignable_Normal_N_Chimeric
    ) %>%
    pivot_longer(
        cols = c("Total Reads", "Unique Reads", "Alignable Reads"),
        names_to = "Sequencing_Metric",
        values_to = "Depth"
    )

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
#   Sequencing_Metric estimate statistic  p.value parameter conf.low conf.high      label
#   <chr>                <dbl>     <dbl>    <dbl>     <int>    <dbl>     <dbl>      <chr>
# 1 Alignable Reads      0.787      3.61 0.00689          8    0.312     0.9471      R = 0.79, p = 0.0069
# 2 Total Reads          0.780      3.52 0.00783          8    0.295     0.9452      R = 0.78, p = 0.0078
# 3 Unique Reads         0.884      5.35 0.000688         8    0.574     0.9723      R = 0.88, p = 0.0007

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
    height = 5.5, ############ NOT 8.5
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

# CACHING for df.chromosome.data
cache_file_chromosome_data <- "../data/df.chromosome.data.rds"

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

df.chromosome.data
df.chromosome.data %>% head()

df.loop.deep.sample.all %>% dim() # 58992   12
df.loop.deep.sample.all %>% head(3)

# CACHING for df.DISTINCT.loop.deep.sample.all
cache_file_distinct_loop <- "../data/df.DISTINCT.loop.deep.sample.all.rds"

if (file.exists(cache_file_distinct_loop)) {
    message("Loading cached distinct loop data from: ", cache_file_distinct_loop)
    df.DISTINCT.loop.deep.sample.all <- readRDS(cache_file_distinct_loop)
} else {
    message("Processing and caching distinct loop data...")
    # common loops between samples picked only one: df.DISTINCT.loop.deep.sample.all
    df.DISTINCT.loop.deep.sample.all <- df.loop.deep.sample.all %>%
        dplyr::select(loop.id) %>% # line 180
        distinct() %>%
        separate(loop.id, into = c("chr1", "x1", "x2", "chr2", "y1", "y2", "end.distance"), sep = "_", remove = FALSE, convert = TRUE) %>%
        mutate(distance = y2 - x2) %>%
        mutate(resolution = convert_to_resolution(end.distance)) %>% # Using utility function # utils_functions.R
        mutate(x0 = x1, x3 = x2, y0 = y1, y3 = y2) %>% # x0, x3, y0, y3
        mutate(mid_x = round((x1 + x2) / 2), mid_y = round((y1 + y2) / 2)) %>%
        mutate(mid_loop = round((mid_x + mid_y) / 2)) %>%
        left_join(df.chromosome.data, by = c("chr1" = "chr")) %>% # chromosome length info
        dplyr::rename(chr.end.coord = end)

    saveRDS(df.DISTINCT.loop.deep.sample.all, cache_file_distinct_loop)
}

df.DISTINCT.loop.deep.sample.all %>% dim() # 31773    11
df.DISTINCT.loop.deep.sample.all %>% head(3)
df.DISTINCT.loop.deep.sample.all %>% count(resolution) # 31773/58992
# resolution     n
# 1         5K  6680
# 2        10K 12162
# 3        25K 12931

########################
# Exploratory Data analysis (EDA)  for df.DISTINCT.loop.deep.sample.all.lt.2mb (< 2Mb)
########################

df.DISTINCT.loop.deep.sample.all.lt.2mb.stats <- df.DISTINCT.loop.deep.sample.all.lt.2mb %>%
    summarise(
        Q1 = quantile(distance, 0.25), # 150,000
        Median = median(distance), # 195,000
        Q3 = quantile(distance, 0.75) # 375,000 (<-> 2,000,000)
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
        plot.title = element_text(hjust = 0.5) # title centering
    )

########################
# 2. CTCF
# 2-1. Exploratory Data analysis (EDA) : CTCF
########################
# Linux
# df.init.ctcf<-read.table(file="/home/pkim/Desktop/temp/enhancer/dropbox_enhancer_doosan/data/ctcf/fimo_2nd_trial_meme/fimo_2nd_trial.txt", header=TRUE, sep="\t")
df.init.ctcf <- read.table(file = "../data/ctcf/fimo_2nd_trial_meme/fimo_2nd_trial.txt", header = TRUE, sep = "\t")
df.init.ctcf %>% count() # 5767921

# Mac
# ctcf<-read.table(file="/home/pkim/Desktop/temp/enhancer/dropbox_enhancer_doosan/data/ctcf/fimo_2nd_trial_meme/fimo_2nd_trial.txt", sep="\t", col.names = c("chr", "start", "end", "strand", "length"), header = FALSE) %>%
# df.init.ctcf<-read.table(file="~/dropbox/K P/Gateway_to_Hao/enhancer/data/ctcf/fimo_2nd_trial_meme/fimo_2nd_trial.txt", header=TRUE, sep="\t")
df.init.ctcf <- read.table(file = "~/dropbox/Gateway_to_Hao/enhancer/data/ctcf/submission/E4/fimo_E4_submission_trial.txt", header = TRUE, sep = "\t") %>%
    dplyr::rename(chr = sequence_name, end = stop) %>%
    mutate(length = end - start)

df.init.ctcf %>% dim() # .4:6551641
df.init.ctcf %>% head(3) # chr    start      end strand length
df.init.ctcf %>% count(length) # length distribution: 9 ~ 39

# strand checking
df.init.ctcf %>% # +: 2891072, -: 2876849 = 5767921//+:3267352, -:3284289
    count(strand)
# removing dups including strand
df.init.ctcf %>%
    distinct(chr, start, end, strand) # .4: 3331233/6551641
# removing dups excluding strand
df.init.ctcf %>%
    distinct(chr, start, end) # .4: 3191859/6551641
# removing dups with all including length
df.init.ctcf %>%
    distinct() # .4: 3331233/6551641

########################
# 2. CTCF
# 2-2. CTCF data preprocessing: id and dedup GRange Obj. (df.DISTINCT.fimo.2nd.trial.ctcf/ df.DISTINCT.ctcf.2nd.fimo.GR)
########################
# CACHING for df.DISTINCT.fimo.2nd.trial.ctcf
cache_file_distinct_fimo_2nd_ctcf <- "../data/df.DISTINCT.fimo.2nd.trial.ctcf.rds"

if (file.exists(cache_file_distinct_fimo_2nd_ctcf)) {
    message("Loading cached DISTINCT fimo 2nd trial ctcf data from: ", cache_file_distinct_fimo_2nd_ctcf)
    df.DISTINCT.fimo.2nd.trial.ctcf <- readRDS(cache_file_distinct_fimo_2nd_ctcf)
} else {
    message("Processing and caching DISTINCT fimo 2nd trial ctcf data...")
    df.DISTINCT.fimo.2nd.trial.ctcf <- df.init.ctcf %>%
        distinct(chr, start, end) %>% # .4:3191859 ************** NO STRAND INFO
        mutate(start = as.numeric(start)) %>%
        mutate(end = as.numeric(end)) %>%
        mutate(ctcf_pos = as.numeric(round((start + end) / 2))) %>%
        mutate(id = str_c(chr, "_", start, "_", end, "_", ctcf_pos))

    saveRDS(df.DISTINCT.fimo.2nd.trial.ctcf, cache_file_distinct_fimo_2nd_ctcf)
}

df.DISTINCT.fimo.2nd.trial.ctcf %>% dim() # .4: 3191859/6551641 :0.4871847
df.DISTINCT.fimo.2nd.trial.ctcf %>% head(3) # chr    start      end ctcf_pos  id

# GRanges Obj.: df.DISTINCT.ctcf.2nd.fimo.GR
df.DISTINCT.ctcf.2nd.fimo.GR <- GRanges(
    seqnames = as.character(df.DISTINCT.fimo.2nd.trial.ctcf$chr),
    ranges = IRanges(
        start = df.DISTINCT.fimo.2nd.trial.ctcf$start,
        end = df.DISTINCT.fimo.2nd.trial.ctcf$end
    )
)

# metadata: id
mcols(df.DISTINCT.ctcf.2nd.fimo.GR)$id <- df.DISTINCT.fimo.2nd.trial.ctcf$id
mcols(df.DISTINCT.ctcf.2nd.fimo.GR)$ctcf_pos <- df.DISTINCT.fimo.2nd.trial.ctcf$ctcf_pos

########################
# 2. CTCF
# 2-3. overall distribution of CTCF on loops
########################
OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping.lt.2mb %>% head(3)
OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping.lt.2mb.GR <- creating_granges(OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping.lt.2mb) # utils_functions.R

index.distinct.ctcf.w.OVERALL.whole.loop <- findOverlaps(
    df.DISTINCT.ctcf.2nd.fimo.GR,
    OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping.lt.2mb.GR,
    type = "any",
    select = "all"
)

OVERALL.loop.for.ctcf.hits <- subjectHits(index.distinct.ctcf.w.OVERALL.whole.loop)
OVERALL.ctcf.on.loop.hits <- queryHits(index.distinct.ctcf.w.OVERALL.whole.loop)

df.ctcf.dist.result <- tibble(
    loop.id = OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping.lt.2mb$loop.id[OVERALL.loop.for.ctcf.hits],
    loop.start = OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping.lt.2mb$x0[OVERALL.loop.for.ctcf.hits],
    loop.end = OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping.lt.2mb$y3[OVERALL.loop.for.ctcf.hits],
    loop.res = OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping.lt.2mb$resolution[OVERALL.loop.for.ctcf.hits],
    ctcf.id = df.DISTINCT.fimo.2nd.trial.ctcf$id[OVERALL.ctcf.on.loop.hits],
    ctcf.pos = df.DISTINCT.fimo.2nd.trial.ctcf$ctcf_pos[OVERALL.ctcf.on.loop.hits],
) %>%
    mutate(across(where(is.numeric), ~ format(., scientific = FALSE)))

df.ctcf.dist.result %>% dim() # 40250853(any, w/o capping + lt 2mb)

# converting data types
relative.pos.df.ctcf.dist.result <- df.ctcf.dist.result %>%
    mutate(
        loop.start = as.numeric(loop.start),
        loop.end = as.numeric(loop.end),
        pos_coord = as.numeric(ctcf.pos),
        loop_length = (loop.end - loop.start),
        relative_pos = pos_coord - loop.start,
        value = relative_pos / (loop_length / 3) - 1
    ) %>%
    dplyr::select(loop.id, ctcf.id, value, loop_length, loop.res)

checking_component_distribution(relative.pos.df.ctcf.dist.result, "ctcf") # utils_functions.R
relative.pos.df.ctcf.dist.result %>% head(3)

########################
# 2. CTCF
# 2-3. overall distribution of CTCF on loops: figures
# 2-3-1. by CHROMOSOME
########################
plot.ctcf.hist <- plot_histogram(relative.pos.df.ctcf.dist.result) # utils_functions.R
plot.ctcf.dens <- plot_density(relative.pos.df.ctcf.dist.result) # utils_functions.R

saving_combined_plot(
    plot.ctcf.hist, plot.ctcf.dens, # utils_functions.R
    "figures/submission/lt2mb/overall_distribution_of_CTCF_by_resolution_wo_capping_lt2mb.pdf"
)

########################
# 2. CTCF
# 2-4. Distribution of CTCF at ends in a loop: for the number of CTCF used in filtering valid loops: figures
# so, the object should be used one without padding: df.DISTINCT.loop.deep.sample.all
# CACHING for df.overlapping.CTCF.w.BOTH.result
cache_file_overlapping_CTCF_BOTH <- "../data/df.overlapping.CTCF.w.BOTH.result.rds"

if (file.exists(cache_file_overlapping_CTCF_BOTH)) {
    message("Loading cached overlapping CTCF BOTH result from: ", cache_file_overlapping_CTCF_BOTH)
    df.overlapping.CTCF.w.BOTH.result <- readRDS(cache_file_overlapping_CTCF_BOTH)

    # We still need df.DISTINCT.loop.deep.sample.all.lt.2mb for later steps, so we define it here if loading from cache
    # Note: This assumes df.DISTINCT.loop.deep.sample.all is available (it should be from earlier steps)
    df.DISTINCT.loop.deep.sample.all.lt.2mb <- df.DISTINCT.loop.deep.sample.all %>%
        filter(distance < 2000000)
} else {
    message("Processing and caching overlapping CTCF BOTH result...")

    # 2-4-3. data processing for getting information of CTCF at ends in loops
    # loops only less than 2mb: 31019 from DISTINCT loops
    df.DISTINCT.loop.deep.sample.all.lt.2mb <- df.DISTINCT.loop.deep.sample.all %>%
        filter(distance < 2000000) # 31019/31773, only use less than 2mb

    df.DISTINCT.loop.deep.sample.all.lt.2mb %>% dim() # 31019
    df.DISTINCT.loop.deep.sample.all.lt.2mb %>% head(3) # 31019
    df.DISTINCT.loop.deep.sample.all.lt.2mb.up.GR <- creating_granges(df.DISTINCT.loop.deep.sample.all.lt.2mb,
        direction = "up"
    ) # utils_functions.R
    df.DISTINCT.loop.deep.sample.all.lt.2mb.up.GR # 31019

    # 1. CTCF + UPSTREAM (df.DISTINCT.loop.deep.sample.all.lt.2mb, df.DISTINCT.fimo.2nd.trial.ctcf)
    index.distinct.ctcf.w.up.loop <- findOverlaps(df.DISTINCT.ctcf.2nd.fimo.GR,
        df.DISTINCT.loop.deep.sample.all.lt.2mb.up.GR,
        type = "any",
        select = "all"
    )

    end.loop.up.ctcf.hits <- subjectHits(index.distinct.ctcf.w.up.loop)
    end.ctcf.up.hits <- queryHits(index.distinct.ctcf.w.up.loop)

    df.overlapping.CTCF.w.UPSTREAM.result <- tibble(
        up.loop.id = mcols(df.DISTINCT.loop.deep.sample.all.lt.2mb.up.GR)$loop.id[end.loop.up.ctcf.hits],
        end.up.distance = mcols(df.DISTINCT.loop.deep.sample.all.lt.2mb.up.GR)$end.distance[end.loop.up.ctcf.hits],
        distance = mcols(df.DISTINCT.loop.deep.sample.all.lt.2mb.up.GR)$distance[end.loop.up.ctcf.hits],
        resolution = mcols(df.DISTINCT.loop.deep.sample.all.lt.2mb.up.GR)$resolution[end.loop.up.ctcf.hits],
        ctcf.id = mcols(df.DISTINCT.ctcf.2nd.fimo.GR)$id[end.ctcf.up.hits],
        WHERE = "UP"
    ) %>%
        unite(ctcf.loop.up.id, up.loop.id, distance, ctcf.id, WHERE, sep = "|", remove = FALSE)

    df.overlapping.CTCF.w.UPSTREAM.result %>% dim() # .4 any lt2mb 1381954
    df.overlapping.CTCF.w.UPSTREAM.result %>% head()

    df.DISTINCT.loop.deep.sample.all.lt.2mb.down.GR <- creating_granges(df.DISTINCT.loop.deep.sample.all.lt.2mb,
        direction = "down"
    ) # utils_functions.R

    # 2. CTCF + DOWNSTREAM (df.DISTINCT.loop.deep.sample.all.lt.2mb, df.DISTINCT.fimo.2nd.trial.ctcf)
    index.distinct.ctcf.w.down.loop <- findOverlaps(df.DISTINCT.ctcf.2nd.fimo.GR,
        df.DISTINCT.loop.deep.sample.all.lt.2mb.down.GR,
        type = "any",
        select = "all"
    )

    end.loop.down.ctcf.hits <- subjectHits(index.distinct.ctcf.w.down.loop)
    end.ctcf.down.hits <- queryHits(index.distinct.ctcf.w.down.loop)

    df.overlapping.CTCF.w.DOWNSTREAM.result <- tibble(
        down.loop.id = mcols(df.DISTINCT.loop.deep.sample.all.lt.2mb.down.GR)$loop.id[end.loop.down.ctcf.hits],
        end.down.distance = mcols(df.DISTINCT.loop.deep.sample.all.lt.2mb.down.GR)$end.distance[end.loop.down.ctcf.hits],
        distance = mcols(df.DISTINCT.loop.deep.sample.all.lt.2mb.down.GR)$distance[end.loop.down.ctcf.hits],
        resolution = mcols(df.DISTINCT.loop.deep.sample.all.lt.2mb.down.GR)$resolution[end.loop.down.ctcf.hits],
        ctcf.id = mcols(df.DISTINCT.ctcf.2nd.fimo.GR)$id[end.ctcf.down.hits],
        WHERE = "DOWN"
    ) %>%
        unite(ctcf.loop.down.id, down.loop.id, distance, ctcf.id, WHERE, sep = "|", remove = FALSE)

    df.overlapping.CTCF.w.DOWNSTREAM.result %>% dim() # sub.4 any lt2mb 1393938
    df.overlapping.CTCF.w.DOWNSTREAM.result %>% count(end.down.distance)

    ########## bind_rows(UPSTREAM & DOWNSTREAM) -> BOTH
    df.overlapping.CTCF.w.BOTH.result <- bind_rows(
        df.overlapping.CTCF.w.UPSTREAM.result %>%
            mutate(loop.id = up.loop.id, end.distance = end.up.distance) %>%
            mutate(case.id = ctcf.loop.up.id) %>%
            dplyr::select(-c(up.loop.id, end.up.distance, ctcf.loop.up.id)),
        df.overlapping.CTCF.w.DOWNSTREAM.result %>%
            mutate(loop.id = down.loop.id, end.distance = end.down.distance) %>%
            mutate(case.id = ctcf.loop.down.id) %>%
            dplyr::select(-c(down.loop.id, end.down.distance, ctcf.loop.down.id))
    ) %>%
        mutate(chr = ifelse(WHERE == "UP", str_split_n(loop.id, "_", 1), str_split_n(loop.id, "_", 4))) %>% # utils_functions.R
        mutate(chr = factor(chr, levels = c(paste0("chr", 1:20), "chrX", "chrY"))) %>%
        mutate(WHERE = fct_relevel(WHERE, "UP", "DOWN")) %>%
        mutate(resolution = convert_to_resolution(end.distance)) # Using utility function # utils_functions.R

    saveRDS(df.overlapping.CTCF.w.BOTH.result, cache_file_overlapping_CTCF_BOTH)
}

df.overlapping.CTCF.w.BOTH.result %>% dim() # sub.4 any lt2mb 2775892
df.overlapping.CTCF.w.BOTH.result %>% head(3) # distance, resolution, ctcf.id, WHERE, loop.id, end.distance, case.id, chr
df.overlapping.CTCF.w.BOTH.result %>% colnames() # "distance" "resolution" "ctcf.id" "WHERE" "loop.id" "end.distance" "case.id" "chr"
df.overlapping.CTCF.w.BOTH.result %>% count(resolution)

df.overlapping.CTCF.w.BOTH.result %>% distinct(loop.id) # .4 any lt2mb 30,908

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
file.tss.list <- fs::dir_ls("/home/pkim/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/", regexp = "\\.txt$")
# Mac
file.tss.list <- fs::dir_ls("~/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss", regexp = "\\.txt$")

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
tss <- read.table(file = "/home/pkim/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/ucsc_start_codon.txt", sep = "\t", head = F)
# Mac
tss <- read.table(file = "~/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/ucsc_start_codon.txt", sep = "\t", head = F)
tss # 17849
tss %>% head()
head(tss)[, c(1, 4, 5, 7, 9)]
tss.select <- tss[, c(1, 4, 5, 7, 9)] # onlty take the relevant columns
tss.select
tss.select.sparate <- tss.select %>%
    separate(V9, into = c("gene_id", "transcript_id", "exon_number", "exon_id", "gene_name"), sep = "; ") %>%
    mutate(across(everything(), ~ gsub(".+\\s", "", .))) %>%
    mutate(gene_name = str_replace(gene_name, ";", ""))
tss.select.sparate %>% head()
tss.select.sparate %>% dim() # 17849
# tss.select.sparate %>% filter(gene_id != gene_name) # 0 rows, so gene_id and gene_name are consistent

########################
# 3. TSS
########################
#############################
# tss resource 2: deprecated (EXON used)
#############################
df.refgene.gtf <- read_tsv("~/dropbox/Gateway_to_Hao/workshop/2023_NIH_meeting/loop_N_tss/ucsc_refGene.gtf", # download from ucsc: https://hgdownload.soe.ucsc.edu/goldenPath/rn7/bigZips/genes/
    comment = "#",
    col_names = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"),
    col_types = cols(.default = "c")
) # char

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
# retrieving exon data 1 from RefSeq GTF
####################################################
# CACHING: df.refgene.gtf.for.exon processing
cache_file_refseq <- "../data/df_refgene_gtf_for_exon.rds"

if (file.exists(cache_file_refseq)) {
    message("Loading cached RefSeq exon data from: ", cache_file_refseq)
    df.refgene.gtf.for.exon <- readRDS(cache_file_refseq)
} else {
    message("Processing RefSeq exon data...")
    df.refgene.gtf.for.exon.raw <- df.refgene.gtf %>%
        filter(feature == "exon") %>%
        filter(chr %in% c(paste0("chr", 1:20), "chrX", "chrY"))

    # checking attribute keys
    refgene.exon.attribute.keys <- get_attribute_keys(df.refgene.gtf.for.exon.raw$attribute) # utils_functions.R
    # refgene.exon.attribute.keys
    # [1] "exon_id"       "exon_number"   "gene_id"       "gene_name"     "transcript_id"

    # adding columns from attribute
    df.refgene.gtf.for.exon.attribute <- df.refgene.gtf.for.exon.raw %>%
        bind_cols(df.refgene.gtf.for.exon.raw$attribute %>% map_dfr(~ extracting_attributes(.x, keys = refgene.exon.attribute.keys))) %>% # utils_functions.R
        dplyr::select(-c(source, feature, attribute, score, frame))

    # df.refgene.gtf.for.exon.attribute # 174,418
    # df.refgene.gtf.for.exon.attribute %>% # 174,418
    #   filter(gene_id != gene_name) # 0 rows, so gene_id and gene_name are consistent

    # filtering: get the last exon (largest exon_number) per gene_id
    df.refgene.gtf.for.exon <- df.refgene.gtf.for.exon.attribute %>%
        mutate(exon_number = as.numeric(exon_number)) %>% # convert to numeric
        group_by(gene_id) %>%
        slice_max(order_by = exon_number, n = 1, with_ties = FALSE) %>% # keep only the largest exon_number
        ungroup() %>%
        mutate(refseq_exon_id = str_c(chr, ":", start, ":", end, ":", strand, ":", gene_name, ":", gene_id, ":", exon_number))

    saveRDS(df.refgene.gtf.for.exon, cache_file_refseq)
}

df.refgene.gtf.for.exon %>% dim() # 17,488
df.refgene.gtf.for.exon %>% head(3)
df.refgene.gtf.for.exon %>%
    add_count(gene_id) %>%
    filter(n > 1) # should be 0

############################
# tss resource 3
############################
# TSS load
cache_file_ensembl_tss <- "../data/df_ensembl_gtf_for_tss_DISTINCT_geneid.rds"

if (file.exists(cache_file_ensembl_tss)) {
    message("Loading cached Ensembl TSS data from: ", cache_file_ensembl_tss)
    df.ensembl.gtf.for.tss.DISTINCT.geneid <- readRDS(cache_file_ensembl_tss)
} else {
    message("Processing Ensembl exon data...")
    df.ensembl.gtf <- read_tsv("~/dropbox/Gateway_to_Hao/enhancer/data/Rattus_norvegicus.mRatBN7.2.113.gtf", # download (https://ftp.ensembl.org/pub/release-113/gtf/rattus_norvegicus/)
        comment = "#",
        col_names = FALSE
    ) # 1,284,446 × 9

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
    df.ensembl.gtf %>%
        dplyr::select(X9) %>%
        head(3) # attribute

    colnames(df.ensembl.gtf) <- c(
        "chr", "source", "feature", "start", "end",
        "score", "strand", "frame", "attribute"
    )

    # retrieving tss data
    df.ensembl.gtf.for.tss <- df.ensembl.gtf %>%
        filter(feature == "start_codon") %>%
        filter(chr %in% c(as.character(1:20), "X", "Y")) %>%
        mutate(chr = str_c("chr", chr))

    df.ensembl.gtf.for.tss %>% dim() # 42925
    df.ensembl.gtf.for.tss %>% head(3)
    df.ensembl.gtf.for.tss %>%
        dplyr::select(attribute) %>%
        head(3)

    df.ensembl.gtf.for.tss %>% filter(str_detect(attribute, "ENSRNOG00000042691"))

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
        bind_cols(df.ensembl.gtf.for.tss$attribute %>% map_dfr(~ extracting_attributes(.x, keys = tss.attribute.keys))) %>% # utils_functions.R
        dplyr::select(-c(source, feature, attribute, score, frame, exon_number))
    df.ensembl.gtf.for.tss.attribute # 42,925

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
    df.ensembl.gtf.for.tss.DISTINCT.geneid <- bind_rows(
        df.ensembl.gtf.for.tss.attribute.filtered.tag.gnbt.geneid.NOdup,
        df.ensembl.gtf.for.tss.attribute.filtered.tag.gnbt.geneid.dup.1pick
    ) %>% # 21776
        dplyr::select(-n) %>%
        mutate(tss.id = paste(chr, start, end, strand, gene_id, gene_name, sep = ":")) %>%
        dplyr::select(chr, start, end, strand, gene_id, gene_name, tss.id)

    df.ensembl.gtf.for.tss.DISTINCT.geneid %>% # 21,725
        # count(chr, start, end, gene_id) #%>% # 21,725
        count(chr, start, end, strand, gene_id) # %>% # 21,725

    df.ensembl.gtf.for.tss.DISTINCT.geneid %>% count(gene_id) # 21,725

    # tss data integrity check
    df.ensembl.gtf.for.tss.DISTINCT.geneid # %>%
    # count(chr) %>% print(n = Inf)
    # head(3) # chr     start       end strand gene_id        gene_name
    # dim() # 21,725
    # distinct(gene_id) # 21,725
    # filter(str_detect(gene_id, 'LOC|RGD')) # 0

    saveRDS(df.ensembl.gtf.for.tss.DISTINCT.geneid, cache_file_ensembl_tss)
}

df.ensembl.gtf.for.tss.DISTINCT.geneid %>% dim() # 21,725
df.ensembl.gtf.for.tss.DISTINCT.geneid %>% head(2)

####################################################
# retrieving exon data 2 from Ensembl GTF
####################################################
# exon load
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
# CACHING: This step is slow, so saving/loading the result
cache_file <- "../data/df_ensembl_gtf_for_exon_attribute.rds"

if (file.exists(cache_file)) {
    message("Loading cached exon attribute data from: ", cache_file)
    df.ensembl.gtf.for.exon.attribute <- readRDS(cache_file)
} else {
    message("Processing exon attribute data (this may take a while)...")
    df.ensembl.gtf.for.exon.attribute <- df.ensembl.gtf.for.exon.raw %>%
        bind_cols(df.ensembl.gtf.for.exon.raw$attribute %>% map_dfr(~ extracting_attributes(.x, keys = exon.attribute.keys))) %>% # utils_functions.R
        dplyr::select(-c(source, feature, attribute, score, frame))

    message("Saving result to cache: ", cache_file)
    saveRDS(df.ensembl.gtf.for.exon.attribute, cache_file)
}
# 526,204

# CACHING for df.tss.ensembl
cache_file_tss_ensembl <- "../data/df.tss.ensembl.rds"

if (file.exists(cache_file_tss_ensembl)) {
    message("Loading cached df.tss.ensembl from: ", cache_file_tss_ensembl)
    df.tss.ensembl <- readRDS(cache_file_tss_ensembl)
} else {
    message("Processing and caching df.tss.ensembl...")

    # filtering
    df.ensembl.gtf.for.exon <- df.ensembl.gtf.for.exon.attribute %>% # 526,204
        filter(tag == "Ensembl_canonical") %>% # 1st filter
        filter(gene_biotype == "protein_coding") %>% # 2nd filter
        mutate(exon_number = as.numeric(exon_number)) %>% # convert to numeric
        group_by(gene_id) %>%
        slice_max(order_by = exon_number, n = 1, with_ties = FALSE) %>% # keep only the largest exon_number
        ungroup() %>%
        mutate(ensembl_exon_id = str_c(chr, ":", start, ":", end, ":", strand, ":", gene_id, ":", gene_name, ":", exon_number)) %>%
        dplyr::select(-c(gene_biotype, tag, transcript_biotype, transcript_version))

    df.tss.ensembl <- df.ensembl.gtf.for.tss.DISTINCT.geneid %>%
        left_join(df.ensembl.gtf.for.exon %>% dplyr::select(gene_id, ensembl_exon_id), by = "gene_id") %>% # 21,725
        left_join(df.refgene.gtf.for.exon %>% dplyr::select(gene_name, refseq_exon_id), by = "gene_name")

    # integrity check
    print(df.tss.ensembl %>% filter(!is.na(ensembl_exon_id) & !is.na(refseq_exon_id))) # 14,071/21,725
    print(df.tss.ensembl %>% filter(is.na(ensembl_exon_id) & is.na(refseq_exon_id))) # 1,360/21,725
    print(df.tss.ensembl %>% filter(is.na(ensembl_exon_id) & !is.na(refseq_exon_id))) # 0/21,725
    print(df.tss.ensembl %>% filter(!is.na(ensembl_exon_id) & is.na(refseq_exon_id))) # 6,294/21,725

    saveRDS(df.tss.ensembl, cache_file_tss_ensembl)
}
df.tss.ensembl %>% dim() # 21,725
df.tss.ensembl %>% head(2)

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
OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping.lt.2mb.GR # 30928

index.distinct.tss.w.OVERALL.whole.loop <- findOverlaps(
    df.tss.ensembl.GR,
    OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping.lt.2mb.GR,
    type = "any",
    select = "all"
)

index.distinct.tss.w.OVERALL.whole.loop # any: 359130/ w/o capping lt2mb: 208296/ ENSEMBL: 259900

OVERALL.loop.for.tss.hits <- subjectHits(index.distinct.tss.w.OVERALL.whole.loop)
OVERALL.tss.on.loop.hits <- queryHits(index.distinct.tss.w.OVERALL.whole.loop)

df.tss.dist.result <- tibble(
    loop.id = OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping.lt.2mb$loop.id[OVERALL.loop.for.tss.hits],
    loop.start = OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping.lt.2mb$x0[OVERALL.loop.for.tss.hits],
    loop.end = OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping.lt.2mb$y3[OVERALL.loop.for.tss.hits],
    loop.res = OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping.lt.2mb$resolution[OVERALL.loop.for.tss.hits],
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
    mutate(
        tss_start = as.numeric(tss_start),
        tss_end = as.numeric(tss_end),
        loop.start = as.numeric(loop.start),
        loop.end = as.numeric(loop.end)
    ) %>%
    mutate(pos_coord = (tss_start + tss_end) / 2) %>% # coordinate for midpoint
    mutate(loop_length = (loop.end - loop.start)) %>% # overall length of the loop
    mutate(relative_pos = pos_coord - loop.start) %>% # relative pos
    # mutate(value = (relative_pos/(loop_length / 2)) - 1) %>% for .5 distance
    mutate(value = (relative_pos / (loop_length / 3)) - 1) %>% # 1 distance
    mutate(resolution = case_when(
        loop.res == 5000 ~ "5K",
        loop.res == 10000 ~ "10K",
        loop.res == 25000 ~ "25K",
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
saving_combined_plot(
    plot.tss.hist, plot.tss.dens, # utils_functions.R
    "figures/submission/lt2mb/overall_distribution_of_TSS_by_resolution_wo_capping_lt2mb_ENSEMBL.pdf"
)

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

# CACHING for df.promoter.rn7.epd
cache_file_promoter_epd <- "../data/df.promoter.rn7.epd.rds"

if (file.exists(cache_file_promoter_epd)) {
    message("Loading cached df.promoter.rn7.epd from: ", cache_file_promoter_epd)
    df.promoter.rn7.epd <- readRDS(cache_file_promoter_epd)
} else {
    message("Processing and caching df.promoter.rn7.epd...")

    ############################
    # promoter resource 2
    ############################
    # Download: https://epd.expasy.org/ftp/epdnew/R_norvegicus/
    # CACHING: df.Rn_EPDnew_001_rn7 processing
    cache_file_epd_rn7 <- "../data/df_Rn_EPDnew_001_rn7.rds"

    if (file.exists(cache_file_epd_rn7)) {
        message("Loading cached EPD rn7 data from: ", cache_file_epd_rn7)
        df.Rn_EPDnew_001_rn7 <- readRDS(cache_file_epd_rn7)
    } else {
        message("Processing EPD rn7 data...")
        Rn_EPDnew_001_rn6.bed.raw <- read.table("~/dropbox/Gateway_to_Hao/enhancer/data/epdnew/001/Rn_EPDnew_001_rn6.bed", sep = " ") %>%
            dplyr::rename(
                seqnames = V1,
                start = V7,
                name = V4,
                score = V5,
                strand = V6
            ) %>%
            dplyr::select(-starts_with("V")) %>%
            mutate(start = as.numeric(start), end = start + 1, score = 1)

        # Rn_EPDnew_001_rn6.bed.raw # 12,601
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
        # Rn_EPDnew_001_rn6.bed.gr
        BiocIO::export(Rn_EPDnew_001_rn6.bed.gr, "/Users/pete/dropbox/Gateway_to_Hao/enhancer/data/epdnew/001/Rn_EPDnew_001_rn6.bed.gr", format = "BED")
        Rn_EPDnew_001_rn6.bed <- import("/Users/pete/dropbox/Gateway_to_Hao/enhancer/data/epdnew/001/Rn_EPDnew_001_rn6.bed.gr", format = "BED")
        # Rn_EPDnew_001_rn6.bed

        # chain file
        chain.rn6.to.rn7 <- import.chain("~/dropbox/Gateway_to_Hao/enhancer/data/epdnew/rn6ToRn7.over.chain")

        # 3. liftOver: rn6 → rn7
        Rn_EPDnew_001_rn7.list <- liftOver(Rn_EPDnew_001_rn6.bed, chain.rn6.to.rn7)
        # Rn_EPDnew_001_rn7.list # 12601

        # 4. GRangesList → GRanges
        Rn_EPDnew_001_rn7 <- unlist(Rn_EPDnew_001_rn7.list)

        # 5. BED export
        export(Rn_EPDnew_001_rn7, "~/dropbox/Gateway_to_Hao/enhancer/data/epdnew/001/Rn_EPDnew_001_rn7.bed", format = "BED")
        # export(Rn_EPDnew_001_rn7, "~/dropbox/Gateway_to_Hao/enhancer/data/epdnew/R_norvegicus_epdnew_rn7.bed", format = "BED")

        df.Rn_EPDnew_001_rn7 <- as_tibble(Rn_EPDnew_001_rn7)

        saveRDS(df.Rn_EPDnew_001_rn7, cache_file_epd_rn7)
    }

    # ENSEMBL ID with gene symbol
    df.gene.mapping.for.promoter <- read_tsv("~/dropbox/Gateway_to_Hao/enhancer/data/epdnew/001/db/promoter_ensembl.txt",
        col_names = c("promoter_id", "gene_id"),
        col_types = cols(
            promoter_id = col_character(),
            gene_id = col_character()
        )
    ) %>% # 12,793
        mutate(gene_id = if_else(str_detect(promoter_id, "Cfb_1"), "ENSRNOG00000051158.3", gene_id)) %>% # https://useast.ensembl.org/Rattus_norvegicus/Gene/Idhistory?g=ENSRNOG00000051158
        distinct() # 12,600

    # adding gene_id (ENSEMBL provieded by EPD)
    df.promoter.rn7.raw <- left_join(df.Rn_EPDnew_001_rn7, df.gene.mapping.for.promoter, by = c("name" = "promoter_id"))

    # adding exon information from ENSEMBL and RefSeq
    df.promoter.rn7.exon_id <- df.promoter.rn7.raw %>% # seqnames   start     end width strand name     score gene_id
        mutate(gene_id = if_else(str_detect(name, "AABR07053687"), "ENSRNOG00000015756", gene_id)) %>%
        dplyr::rename(promoter_id = name) %>%
        mutate(gene_name = str_split_n(promoter_id, "_", 1)) %>% # utils_functions.R
        left_join(df.ensembl.gtf.for.exon %>% dplyr::select(gene_id, ensembl_exon_id), by = "gene_id") %>%
        left_join(df.refgene.gtf.for.exon %>% dplyr::select(gene_name, refseq_exon_id), by = "gene_name")

    # Deduplication: keep only _1 promoters when seqnames, start, end, gene_id are identical
    df.promoter.rn7 <- df.promoter.rn7.exon_id %>%
        add_count(seqnames, start, end, gene_id, name = "dup_count") %>%
        filter(dup_count == 1 | str_detect(promoter_id, "_1$")) %>%
        dplyr::select(-dup_count)

    # df.promoter.rn7 (dups)
    df.promoter.rn7.epd <- df.promoter.rn7 %>%
        dplyr::rename(tss_start = start, tss_end = end) %>% # tss_start, tss_end : original coordinates of promoter
        mutate(
            start = tss_start - 40,
            end   = tss_start + 40
        ) %>%
        mutate(promoter.id = str_c(seqnames, ":", start, ":", end, ":", strand, ":", gene_id, ":", gene_name, ":", seqnames, ":", tss_start, ":", tss_end))

    saveRDS(df.promoter.rn7.epd, cache_file_promoter_epd)
}
df.promoter.rn7.epd
df.promoter.rn7.epd.GR <- GRanges(
    seqnames = df.promoter.rn7.epd$seqnames,
    ranges = IRanges(
        start = df.promoter.rn7.epd$start,
        end = df.promoter.rn7.epd$end
    )
)
# metadata
mcols(df.promoter.rn7.epd.GR) <- df.promoter.rn7.epd[, c("promoter.id", "gene_id", "gene_name", "ensembl_exon_id", "refseq_exon_id")]

########################
# 4. promoter
# 4-3. overall distribution of promoter on loops
########################
index.promoter.w.OVERALL.whole.loop <- findOverlaps(
    df.promoter.rn7.epd.GR,
    OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping.lt.2mb.GR,
    type = "any",
    select = "all"
)
index.promoter.w.OVERALL.whole.loop # w/o capping lt2mb ENSMBL 159,944

OVERALL.loop.for.promoter.hits <- subjectHits(index.promoter.w.OVERALL.whole.loop)
OVERALL.promoter.on.loop.hits <- queryHits(index.promoter.w.OVERALL.whole.loop)

df.promoter.dist.result <- data.frame(
    loop.id = OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping.lt.2mb$loop.id[OVERALL.loop.for.promoter.hits],
    loop.start = OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping.lt.2mb$x0[OVERALL.loop.for.promoter.hits],
    loop.end = OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping.lt.2mb$y3[OVERALL.loop.for.promoter.hits],
    loop.res = OVERALL.df.DISTINCT.loop.deep.sample.all.1.distance.wo.capping.lt.2mb$resolution[OVERALL.loop.for.promoter.hits],
    promoter_id = df.promoter.rn7.epd$promoter.id[OVERALL.promoter.on.loop.hits],
    promoter_start = df.promoter.rn7.epd$start[OVERALL.promoter.on.loop.hits],
    promoter_end = df.promoter.rn7.epd$end[OVERALL.promoter.on.loop.hits],
    promoter_gene_id = df.promoter.rn7.epd$gene_id[OVERALL.promoter.on.loop.hits]
) %>%
    mutate(across(where(is.numeric), ~ format(., scientific = FALSE)))

df.promoter.dist.result %>% head()
df.promoter.dist.result %>% dim() # w/o capping lt2mb ENSMBL 159,944

relative.pos.df.promoter.dist.result <- df.promoter.dist.result %>%
    mutate(
        promoter_start = as.numeric(promoter_start),
        promoter_end = as.numeric(promoter_end),
        loop.start = as.numeric(loop.start),
        loop.end = as.numeric(loop.end)
    ) %>%
    mutate(pos_coord = round((promoter_start + promoter_end) / 2)) %>%
    mutate(loop_length = (loop.end - loop.start)) %>% # x0, y3
    mutate(relative_pos = pos_coord - loop.start) %>% # relative position from loop start
    # mutate(value = (relative_pos / (loop_length / 2)) - 1) %>% # for padding .5x distance
    mutate(value = relative_pos / (loop_length / 3) - 1) %>% # for padding 1x distance
    mutate(resolution = case_when(
        loop.res == 5000 ~ "5K",
        loop.res == 10000 ~ "10K",
        loop.res == 25000 ~ "25K",
        TRUE ~ NA_character_
    )) %>%
    dplyr::select(loop.id, promoter_id, value, loop.res, promoter_gene_id)

relative.pos.df.promoter.dist.result %>% head()
relative.pos.df.promoter.dist.result %>% dim() # w/o capping lt2mb ENSMBL 159,944
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

saving_combined_plot(
    plot.promoter.hist, plot.promoter.dens, # utils_functions.R
    "figures/submission/lt2mb/overall_distribution_of_promoter_by_resolution_wo_capping_lt2mb_ENSEMBL.pdf"
)

#########################################
#########################################
# Figure for new ones
#########################################
#########################################

combining_and_save_plots(plot.ctcf.hist, plot.tss.hist, plot.promoter.hist, "histogram_combined_all_ENSEMBL.png") # utils_functions.R
combining_and_save_plots(plot.ctcf.dens, plot.tss.dens, plot.promoter.dens, "density_combined_all_ENSEMBL.png") # utils_functions.R
