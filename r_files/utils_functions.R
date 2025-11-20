################################################################################
# Utility Functions for enhancer_promoter_interaction.R
# Purpose: Extract repeated code patterns into reusable functions
################################################################################

library("tidyverse")
library("GenomicRanges")
library("fs")

################################################################################
# 0. Functions from funcs.R (originally in project_common_code)
################################################################################

#' Safely subset a vector by index
#' @param x Vector to subset
#' @param index Index to extract
#' @return Element at index or NA if index is out of bounds
subset_safely <- function(x, index) {
    if (length(x) < index) {
        return(NA_character_)
    }
    x[[index]]
}

#' Split string and extract nth element
#' @param string Character vector to split
#' @param pattern Pattern to split on
#' @param n Index of element to extract (1-indexed)
#' @return Character vector with nth element from each split
str_split_n <- function(string, pattern, n) {
    out <- str_split(string, pattern)
    vapply(out, subset_safely, character(1L), index = n)
}

#' Read BEDPE file
#' @param file Path to BEDPE file
#' @return Data frame with BEDPE data
bedpe_data_reader <- function(file) {
    read.csv(file, header = TRUE, sep = "\t")
}

#' Initialize BEDPE data frame from file list
#' @param file.list Character vector of file paths
#' @param deli_dir Directory delimiter for sample name extraction
#' @return Tibble with combined BEDPE data and calculated distances
init.bedpe.df <- function(file.list, deli_dir) {
    file.list %>%
        map_dfr(bedpe_data_reader, .id = "sample") %>%
        mutate(sample = str_split_n(sample, str_c(deli_dir, "\\/"), 2)) %>%
        mutate(sample = str_split_n(sample, "_", 1)) %>%
        mutate(distance = ((y2 + y1) - (x2 + x1)) / 2)
}

################################################################################
# 1. Resolution Conversion Functions
################################################################################

#' Convert end.distance to resolution label
#' @param end_distance Numeric vector of end distances
#' @return Factor vector with levels "5K", "10K", "25K"
convert_to_resolution <- function(end_distance) {
    resolution <- case_when(
        end_distance == 5000 ~ "5K",
        end_distance == 10000 ~ "10K",
        end_distance == 25000 ~ "25K",
        TRUE ~ NA_character_
    )
    factor(resolution, levels = c("5K", "10K", "25K"))
}

################################################################################
# 2. Loop ID Parsing Functions
################################################################################

#' Parse loop.id into components
#' @param loop_id Character vector of loop IDs (format: chr_x1_x2_chr_y1_y2_distance)
#' @return Tibble with columns: chr1, x1, x2, chr2, y1, y2, end.distance
parse_loop_id <- function(loop_id) {
    tibble(loop.id = loop_id) %>%
        separate(loop.id,
            into = c("chr1", "x1", "x2", "chr2", "y1", "y2", "end.distance"),
            sep = "_",
            convert = TRUE,
            remove = FALSE
        )
}

################################################################################
# 3. FindOverlaps Result Processing Functions
################################################################################

#' Extract overlap results from findOverlaps into a tibble
#' @param query_gr GRanges object (query)
#' @param subject_gr GRanges object (subject)
#' @param query_metadata Character vector of metadata column names to extract from query
#' @param subject_metadata Character vector of metadata column names to extract from subject
#' @param where_label Character label for WHERE column (e.g., "UP", "DOWN")
#' @return Tibble with overlap information
extract_overlap_results <- function(query_gr, subject_gr,
                                    query_metadata = NULL,
                                    subject_metadata = NULL,
                                    where_label = NULL) {
    # Find overlaps
    overlaps <- findOverlaps(query_gr, subject_gr, type = "any", select = "all")

    query_hits <- queryHits(overlaps)
    subject_hits <- subjectHits(overlaps)

    # Initialize result tibble
    result <- tibble()

    # Extract query metadata
    if (!is.null(query_metadata)) {
        for (col in query_metadata) {
            result[[paste0("query_", col)]] <- mcols(query_gr)[[col]][query_hits]
        }
    }

    # Extract subject metadata
    if (!is.null(subject_metadata)) {
        for (col in subject_metadata) {
            result[[paste0("subject_", col)]] <- mcols(subject_gr)[[col]][subject_hits]
        }
    }

    # Add WHERE label if provided
    if (!is.null(where_label)) {
        result$WHERE <- where_label
    }

    return(result)
}

################################################################################
# 4. Statistical Summary Functions
################################################################################

#' Calculate comprehensive statistics for a numeric vector
#' @param x Numeric vector
#' @return Tibble with min, Q1, median, mean, Q3, max, SD
calculate_stats <- function(x) {
    tibble(
        Min = min(x, na.rm = TRUE),
        Q1 = quantile(x, 0.25, na.rm = TRUE),
        Median = median(x, na.rm = TRUE),
        Mean = mean(x, na.rm = TRUE),
        Q3 = quantile(x, 0.75, na.rm = TRUE),
        Max = max(x, na.rm = TRUE),
        SD = sd(x, na.rm = TRUE)
    )
}

#' Calculate statistics grouped by one or more variables
#' @param df Data frame
#' @param value_col Name of the column to calculate statistics on
#' @param group_cols Character vector of column names to group by
#' @return Tibble with grouped statistics
calculate_grouped_stats <- function(df, value_col, group_cols = NULL) {
    if (is.null(group_cols)) {
        df %>%
            summarise(
                Min = min(.data[[value_col]], na.rm = TRUE),
                Q1 = quantile(.data[[value_col]], 0.25, na.rm = TRUE),
                Median = median(.data[[value_col]], na.rm = TRUE),
                Mean = mean(.data[[value_col]], na.rm = TRUE),
                Q3 = quantile(.data[[value_col]], 0.75, na.rm = TRUE),
                Max = max(.data[[value_col]], na.rm = TRUE),
                SD = sd(.data[[value_col]], na.rm = TRUE)
            )
    } else {
        df %>%
            group_by(across(all_of(group_cols))) %>%
            summarise(
                Min = min(.data[[value_col]], na.rm = TRUE),
                Q1 = quantile(.data[[value_col]], 0.25, na.rm = TRUE),
                Median = median(.data[[value_col]], na.rm = TRUE),
                Mean = mean(.data[[value_col]], na.rm = TRUE),
                Q3 = quantile(.data[[value_col]], 0.75, na.rm = TRUE),
                Max = max(.data[[value_col]], na.rm = TRUE),
                SD = sd(.data[[value_col]], na.rm = TRUE),
                .groups = "drop"
            )
    }
}

################################################################################
# 5. Plotting Helper Functions
################################################################################

#' Create a boxplot with statistics annotations
#' @param df Data frame
#' @param y_col Name of the y-axis column
#' @param fill_color Fill color for the boxplot
#' @param title Plot title
#' @param y_label Y-axis label
#' @return ggplot object
create_annotated_boxplot <- function(df, y_col, fill_color = "#4C9F70",
                                     title = NULL, y_label = NULL) {
    stats <- calculate_stats(df[[y_col]])

    p <- ggplot(df, aes(y = .data[[y_col]])) +
        geom_boxplot(outlier.shape = NA, fill = fill_color, color = "black") +
        geom_hline(yintercept = stats$Min, linetype = "dashed", color = "blue", linewidth = 0.6) +
        geom_hline(yintercept = stats$Median, linetype = "dotted", color = "red", linewidth = 0.8) +
        annotate("text",
            x = 0.2, y = stats$Min,
            label = paste0("Min: ", round(stats$Min)), color = "blue", hjust = 0
        ) +
        annotate("text",
            x = 0.2, y = stats$Median,
            label = paste0("Median: ", round(stats$Median)), color = "red", hjust = 0
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 0, hjust = 1),
            plot.title = element_text(hjust = 0.5),
            legend.position = "bottom"
        )

    if (!is.null(title)) p <- p + labs(title = title)
    if (!is.null(y_label)) p <- p + labs(y = y_label)

    return(p)
}

################################################################################
# 6. Relative Position Calculation Functions
################################################################################

#' Calculate relative position within loop
#' @param df Data frame with loop.start, loop.end, and position columns
#' @param pos_col Name of the position column
#' @param loop_start_col Name of the loop start column
#' @param loop_end_col Name of the loop end column
#' @param padding_factor Padding factor (default: 3 for 1 distance padding)
#' @return Data frame with added value column (relative position)
calculate_relative_position <- function(df, pos_col,
                                        loop_start_col = "loop.start",
                                        loop_end_col = "loop.end",
                                        padding_factor = 3) {
    df %>%
        mutate(
            loop_length = .data[[loop_end_col]] - .data[[loop_start_col]],
            relative_pos = .data[[pos_col]] - .data[[loop_start_col]],
            value = (relative_pos / (loop_length / padding_factor)) - 1
        )
}

################################################################################
# 7. Chromosome Filtering Functions
################################################################################

#' Filter data to standard chromosomes only
#' @param df Data frame with a chromosome column
#' @param chr_col Name of the chromosome column
#' @return Filtered data frame
filter_standard_chromosomes <- function(df, chr_col = "chr") {
    valid_chrs <- c(paste0("chr", 1:20), "chrX", "chrY")
    df %>% filter(.data[[chr_col]] %in% valid_chrs)
}

################################################################################
# 8. Count and Case Functions
################################################################################

#' Determine case (UP, DOWN, BOTH, NONE) for each loop based on WHERE column
#' @param df Data frame with loop.id and WHERE columns
#' @return Data frame with loop.id and case columns
determine_loop_case <- function(df) {
    df %>%
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
}

################################################################################
# 9. Data Binding Functions
################################################################################

#' Bind UP and DOWN results into a single data frame
#' @param df_up Data frame with upstream results
#' @param df_down Data frame with downstream results
#' @param up_id_col Name of the loop.id column in df_up
#' @param down_id_col Name of the loop.id column in df_down
#' @param up_distance_col Name of the end.distance column in df_up
#' @param down_distance_col Name of the end.distance column in df_down
#' @return Combined data frame
bind_up_down_results <- function(df_up, df_down,
                                 up_id_col = "up.loop.id",
                                 down_id_col = "down.loop.id",
                                 up_distance_col = "end.up.distance",
                                 down_distance_col = "end.down.distance") {
    bind_rows(
        df_up %>%
            mutate(
                loop.id = .data[[up_id_col]],
                end.distance = .data[[up_distance_col]]
            ) %>%
            dplyr::select(-all_of(c(up_id_col, up_distance_col))),
        df_down %>%
            mutate(
                loop.id = .data[[down_id_col]],
                end.distance = .data[[down_distance_col]]
            ) %>%
            dplyr::select(-all_of(c(down_id_col, down_distance_col)))
    ) %>%
        mutate(chr = ifelse(WHERE == "UP",
            str_split_n(loop.id, "_", 1),
            str_split_n(loop.id, "_", 4)
        )) %>%
        mutate(chr = factor(chr, levels = c(paste0("chr", 1:20), "chrX", "chrY"))) %>%
        mutate(WHERE = fct_relevel(WHERE, "UP", "DOWN")) %>%
        mutate(resolution = convert_to_resolution(end.distance))
}

message("Utility functions loaded successfully!")
