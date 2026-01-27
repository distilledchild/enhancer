################################################################################
# Utility Functions for enhancer_promoter_interaction.R
# Purpose: Extract repeated code patterns into reusable functions
################################################################################

library("tidyverse")
library("GenomicRanges")
library("fs")
library("ggplot2")
library("patchwork")
library("ggvenn")
library("pdftools")
library("magick")
library("tools")
library("cowplot")
library("circlize")

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

#' Compute inner distance statistics
#' @param df Data frame with x3 and y0 columns
#' @return List with ylim_vals (quantiles) and stats (min, median)
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

#' Save plot in both PDF and PNG formats
#' @param plot_obj ggplot object to save
#' @param filename_base Base filename (without extension)
#' @param output_dir Output directory path
#' @param width_in Width in inches
#' @param height_in Height in inches
#' @param scale_x X-axis scale factor
#' @param scale_y Y-axis scale factor
#' @param dpi DPI for PNG output
#' @param bg Background color
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
        dpi = dpi,
        bg = bg
    )
}

#' Save combined histogram and density plot
#' @param plot_left Left plot object
#' @param plot_right Right plot object
#' @param filename Output filename
#' @param width Plot width
#' @param height Plot height
saving_combined_plot <- function(plot_left, plot_right, filename, width = 11 * 0.8, height = 8.5 * 0.4) {
    pdf(file = filename, width = width, height = height)
    combined_plot <- plot_left | plot_right
    print(combined_plot)
    dev.off()
}

#' Create histogram plot
#' @param df Data frame
#' @param xvar X-axis variable (unquoted)
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @return ggplot object
plot_histogram <- function(df, xvar = value, xlab = "Relative Position to Loop", ylab = "Count") {
    ggplot(df, aes(x = {{ xvar }})) +
        geom_histogram(fill = "skyblue", color = "grey70", alpha = 0.7, bins = 200, linewidth = 0.1) +
        labs(x = xlab, y = ylab) +
        theme(plot.title = element_text(hjust = 0.5))
}

#' Create density plot
#' @param df Data frame
#' @param xvar X-axis variable (unquoted)
#' @param groupvar Grouping variable (unquoted)
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param colors Named vector of colors for groups
#' @return ggplot object
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

#' Plot inner distance boxplot with annotations
#' @param df Data frame with x3 and y0 columns
#' @param ylim_vals Y-axis limits
#' @param stats List with min_val and median_val
#' @param title_text Plot title
#' @return ggplot object
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

#' Combine three plots and save
#' @param p1 First plot
#' @param p2 Second plot
#' @param p3 Third plot
#' @param filename Output filename
#' @param output_dir Output directory
#' @param width Plot width
#' @param height Plot height
#' @param dpi DPI for output
combining_and_save_plots <- function(p1, p2, p3,
                                     filename = "histogram_combined_all.png",
                                     output_dir = "./figures/submission/lt2mb",
                                     width = 11, height = 8.5, dpi = 300) {
    # base theme
    base_theme <- theme_bw(base_size = 12) +
        theme(
            plot.title.position = "plot",
            panel.grid.minor = element_blank()
        )

    # theme
    p1 <- p1 + base_theme
    p2 <- p2 + base_theme
    p3 <- p3 + base_theme

    # combining plots
    fig_combined <- (p1 | p2 | p3) +
        plot_layout(ncol = 3, guides = "collect", widths = c(1, 1, 1)) &
        theme(
            legend.position = "bottom",
            legend.margin = margin(2, 6, 2, 6),
            plot.tag = element_text(face = "bold", size = 12)
        )

    # adding tag
    fig_combined <- fig_combined +
        plot_annotation(tag_levels = "a")

    # saving
    ggsave(file.path(output_dir, filename), fig_combined, width = width, height = height, dpi = dpi, bg = "white")
}

#' Create Venn diagram plot
#' @param ctcf_data Data frame with CTCF loops
#' @param promoter_data Data frame with promoter loops
#' @param tss_data Data frame with TSS loops
#' @param ctcf_label Label for CTCF (not used in current implementation)
#' @return ggplot object
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

################################################################################
# 10. GRanges Creation Functions
################################################################################

#' Create GRanges object from data frame
#' @param df Data frame with genomic coordinates
#' @param direction Direction: NULL (overall), "up" (upstream), or "down" (downstream)
#' @param chr.col Name of chromosome column
#' @param metadata.cols Metadata columns to include
#' @param use_anchor Whether to use anchor coordinates
#' @param point Whether to create GRanges at midpoint
#' @return GRanges object
creating_granges <- function(df, direction = NULL, chr.col = "chr1",
                             metadata.cols = NULL, use_anchor = FALSE,
                             point = FALSE) {
    # ---- 1. Determine start/end columns ----
    if (is.null(direction)) {
        start.col <- "x0"
        end.col <- "y3"
    } else if (direction == "up") {
        if (use_anchor) {
            start.col <- "x1"
            end.col <- "x2"
        } else {
            start.col <- "x0"
            end.col <- "x3"
        }
    } else if (direction == "down") {
        if (use_anchor) {
            start.col <- "y1"
            end.col <- "y2"
        } else {
            start.col <- "y0"
            end.col <- "y3"
        }
    } else {
        stop("Invalid direction. Use 'up', 'down', or NULL.")
    }

    # ---- 2. POINT MODE: create GRanges at midpoint ----
    if (point) {
        midpoint <- as.integer((df[[start.col]] + df[[end.col]]) / 2)

        gr <- GRanges(
            seqnames = df[[chr.col]],
            ranges = IRanges(start = midpoint, end = midpoint)
        )

        if (is.null(metadata.cols)) {
            metadata.cols <- setdiff(names(df), c(chr.col, start.col, end.col))
        }

        mcols(gr) <- df[, metadata.cols, drop = FALSE]
        return(gr)
    }

    # ---- 3. RANGE MODE ----
    gr <- GRanges(
        seqnames = df[[chr.col]],
        ranges = IRanges(
            start = as.integer(df[[start.col]]),
            end = as.integer(df[[end.col]])
        )
    )

    if (is.null(metadata.cols)) {
        metadata.cols <- setdiff(names(df), c(chr.col, start.col, end.col))
    }

    mcols(gr) <- df[, metadata.cols, drop = FALSE]
    return(gr)
}

################################################################################
# 11. PDF/Image Processing Functions
################################################################################

#' Merge PDF pages into a single image
#' @param pdf_path Path to input PDF
#' @param output_pdf_path Path for output PDF
#' @param output_png_path Path for output PNG
#' @param dpi DPI for rendering
#' @param stack Whether to stack images vertically
#' @return List with paths to output PDF and PNG
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

################################################################################
# 12. Component Distribution Analysis
################################################################################

#' Check component distribution by chromosome
#' @param data Data frame with loop and component information
#' @param component_name Name of component (e.g., "CTCF", "TSS")
#' @param output_dir Output directory for plots
checking_component_distribution <- function(data, component_name, output_dir = "./figures/submission/lt2mb") {
    chromosomes <- c(1:20, "X", "Y")
    plot_list <- list()

    for (chr in chromosomes) {
        tryCatch(
            {
                message("START: Processing chromosome: ", chr)

                data_chr <- data %>%
                    filter(str_detect(loop.id, paste0("_chr", chr, "_")))

                # Density plot
                plot_dens <- data_chr %>%
                    ggplot(aes(x = value)) +
                    geom_density(fill = "skyblue", color = "black", alpha = 0.5) +
                    ylim(c(0, 1)) +
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
            },
            error = function(e) {
                message("Error processing chromosome: ", chr)
                message("Error message: ", e$message)
            }
        )
    }

    pdf_path <- file.path(output_dir, paste0("overall_distribution_of_", component_name, "_by_chromosome_wo_capping_lt2mb.pdf"))
    png_path <- paste0(file_path_sans_ext(pdf_path), ".png")

    # saving PDF
    pdf(pdf_path, width = 11 * 0.8, height = 8.5 * 0.8)
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

################################################################################
# 13. Loop Validity and Functionality Functions
################################################################################

#' Add loop validity column (clear vs vague)
#' @param df Data frame with loop.id column
#' @return Data frame with loop_validity column added
adding_loop_validity <- function(df) {
    df %>%
        add_count(loop.id, name = "n_per_loop") %>%
        mutate(loop_validity = ifelse(n_per_loop == 1, "clear", "vague")) %>%
        dplyr::select(-n_per_loop)
}

#' Add gene location pattern column
#' @param df Data frame with loop.id, gene_chr, gene_start, gene_end columns
#' @return Data frame with loop coordinates and gene_LOC pattern column added
#' @details
#' Parses loop.id to extract coordinates and determines if gene is inside or outside loop
#' gene_LOC values:
#' - "2_INS": Both anchors have genes inside loop
#' - "INS_OUT": One anchor has gene inside, one outside
#' - "2_OUT": Both anchors have genes outside loop
adding_gene_location_pattern <- function(df) {
    df %>%
        mutate(
            # Parse loop.id: chr_x1_x2_chr_y1_y2_resolution
            loop_parts = str_split(loop.id, "_"),
            chr1 = map_chr(loop_parts, 1),
            x1 = as.numeric(map_chr(loop_parts, 2)),
            x2 = as.numeric(map_chr(loop_parts, 3)),
            chr2 = map_chr(loop_parts, 4),
            y1 = as.numeric(map_chr(loop_parts, 5)),
            y2 = as.numeric(map_chr(loop_parts, 6)),

            # Calculate midpoints
            mid_x = (x1 + x2) / 2,
            mid_y = (y1 + y2) / 2,

            # Check if gene is inside loop (between x1 and y2)
            gene_inside = case_when(
                is.na(gene_chr) | is.na(gene_start) | is.na(gene_end) ~ NA,
                gene_chr == chr1 & gene_start >= x1 & gene_end <= y2 ~ TRUE,
                TRUE ~ FALSE
            )
        ) %>%
        dplyr::select(-loop_parts) %>%
        # Add gene_LOC pattern per loop
        group_by(loop.id) %>%
        mutate(
            n_inside = sum(gene_inside == TRUE, na.rm = TRUE),
            n_outside = sum(gene_inside == FALSE, na.rm = TRUE),
            gene_LOC = case_when(
                n_inside == 2 ~ "2_INS", # Both anchors have genes inside
                n_inside == 1 & n_outside == 1 ~ "INS_OUT", # One inside, one outside
                n_outside == 2 ~ "2_OUT", # Both anchors have genes outside
                TRUE ~ NA_character_
            )
        ) %>%
        ungroup() %>%
        dplyr::select(-c(n_inside, n_outside))
}

#' Add final decision column for vague cases
#' @param df Data frame with loop_validity and functionality columns
#' @return Data frame with final decision column added
adding_final_decision <- function(df) {
    df %>%
        group_by(loop.id) %>%
        mutate(final = case_when(
            # For clear cases (only 1 row per loop.id)
            loop_validity == "clear" ~ "decided",

            # For vague cases (2 rows per loop.id), check functionality
            loop_validity == "vague" ~ {
                functional_count <- sum(functionality == "functional", na.rm = TRUE)

                if (functional_count == 2) {
                    "undecidable" # Both UP and DOWN are functional
                } else if (functional_count == 1) {
                    "decidable" # Only one (UP or DOWN) is functional
                } else {
                    "FALSE" # Neither is functional
                }
            },
            TRUE ~ NA_character_
        )) %>%
        ungroup()
}

#' Add TSS/Promoter location pattern column
#' @param df Data frame with loop.id and distance columns
#' @return Data frame with TSS_PRO_LOC column added
#' @details
#' TSS_PRO_LOC values:
#' - "0_0": Both anchors have distance = 0 (TSS/Promoter overlaps both anchors)
#' - "0_S": One anchor has distance = 0, the other has distance > 0
#' - "S_S": Both anchors have distance > 0 (TSS/Promoter separated from both anchors)
adding_tss_pro_location_pattern <- function(df) {
    df %>%
        group_by(loop.id) %>%
        mutate(
            n_zero_distance = sum(distance == 0),
            TSS_PRO_LOC = case_when(
                n_zero_distance == 2 ~ "0_0", # Both anchors overlap (distance = 0)
                n_zero_distance == 1 ~ "0_S", # One anchor overlaps, one separated
                n_zero_distance == 0 ~ "S_S", # Both anchors separated (distance > 0)
                TRUE ~ NA_character_
            )
        ) %>%
        ungroup() %>%
        dplyr::select(-n_zero_distance)
}

#' Extract gene coordinates from exon IDs
#' @param df Data frame with ensembl_exon_id and refseq_exon_id columns
#' @return Data frame with gene_chr, gene_start, gene_end, and gene_coord_source columns added
#' @details
#' Extracts chr, start, end from exon IDs (format: chr:start:end:gene_name:gene_id:exon_number)
#' Priority: ENSEMBL first, RefSeq as fallback
#' gene_coord_source tracks which source was used: "both", "ensembl", "refseq", or "none"
adding_gene_coord_from_exon_ids <- function(df) {
    df %>%
        mutate(
            # Extract from ENSEMBL exon ID (format: chr:start:end:gene_name:gene_id:exon_number)
            ensembl_chr = if_else(!is.na(ensembl_exon_id), str_split_n(ensembl_exon_id, ":", 1), NA_character_),
            ensembl_start = if_else(!is.na(ensembl_exon_id), as.numeric(str_split_n(ensembl_exon_id, ":", 2)), NA_real_),
            ensembl_end = if_else(!is.na(ensembl_exon_id), as.numeric(str_split_n(ensembl_exon_id, ":", 3)), NA_real_),

            # Extract from RefSeq exon ID (format: chr:start:end:gene_name:gene_name:exon_number)
            refseq_chr = if_else(!is.na(refseq_exon_id), str_split_n(refseq_exon_id, ":", 1), NA_character_),
            refseq_start = if_else(!is.na(refseq_exon_id), as.numeric(str_split_n(refseq_exon_id, ":", 2)), NA_real_),
            refseq_end = if_else(!is.na(refseq_exon_id), as.numeric(str_split_n(refseq_exon_id, ":", 3)), NA_real_),

            # Priority: Use ENSEMBL if available, otherwise use RefSeq
            gene_chr = if_else(!is.na(ensembl_chr), ensembl_chr, refseq_chr),
            gene_start = if_else(!is.na(ensembl_start), ensembl_start, refseq_start),
            gene_end = if_else(!is.na(ensembl_end), ensembl_end, refseq_end),

            # Track which source was used
            gene_coord_source = case_when(
                !is.na(ensembl_chr) & !is.na(refseq_chr) ~ "both",
                !is.na(ensembl_chr) ~ "ensembl",
                !is.na(refseq_chr) ~ "refseq",
                TRUE ~ "none"
            )
        ) %>%
        # Remove only intermediate columns, keep gene_chr, gene_start, gene_end
        dplyr::select(-c(
            ensembl_chr, ensembl_start, ensembl_end,
            refseq_chr, refseq_start, refseq_end
        ))
}

#' Add gene NA pattern column
#' @param df Data frame with loop.id and gene_chr columns
#' @return Data frame with gene_NA_pattern column added
#' @details
#' gene_NA_pattern values:
#' - "2_NA": Both anchors have NA gene coordinates
#' - "1_NA": One anchor has NA gene coordinates
#' - "0_NA": Both anchors have valid gene coordinates
adding_gene_na_pattern <- function(df) {
    df %>%
        group_by(loop.id) %>%
        mutate(
            n_na_gene = sum(is.na(gene_chr)),
            gene_NA_pattern = case_when(
                n_na_gene == 2 ~ "2_NA", # Both anchors have NA
                n_na_gene == 1 ~ "1_NA", # One anchor has NA
                n_na_gene == 0 ~ "0_NA", # Both anchors have valid coordinates
                TRUE ~ NA_character_
            )
        ) %>%
        ungroup() %>%
        dplyr::select(-n_na_gene)
}

################################################################################
# 14. Loop Overlap and Extraction Functions
################################################################################

#' Extract overlapping loops between datasets
#' @param ctcf_data CTCF data frame
#' @param promoter_data Promoter data frame
#' @param tss_data TSS data frame
#' @return List with ctcf_promoter_only, ctcf_tss_only, and ctcf_promoter_tss_overlap
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

################################################################################
# 15. GTF Attribute Extraction Functions
################################################################################

#' Get attribute keys from GTF attribute column
#' @param attribute_vector Character vector of GTF attributes
#' @return Sorted unique vector of attribute keys
get_attribute_keys <- function(attribute_vector) {
    attribute_vector %>%
        strsplit(";") %>%
        unlist() %>%
        trimws() %>%
        str_extract("^[^ ]+") %>%
        unique() %>%
        sort()
}

#' Extract attribute keys and values dynamically from GTF
#' @param attr_string GTF attribute string
#' @param keys Optional vector of keys to extract (if NULL, extracts all)
#' @return Tibble with extracted attribute values
extracting_attributes <- function(attr_string, keys = NULL) {
    # If keys are not provided, extract them from the attribute string
    if (is.null(keys)) {
        keys <- get_attribute_keys(attr_string)
    }

    # Helper function to safely extract a value for a given key
    safe_extract <- function(key, string) {
        match <- str_match(string, paste0(key, ' "([^"]*)?\"'))[, 2]
        if (is.na(match) || trimws(match) == "") {
            return(NA_character_)
        } else {
            return(match)
        }
    }

    # Create a named list of extracted values
    result <- setNames(
        lapply(keys, function(k) safe_extract(k, attr_string)),
        keys
    )

    # Convert to tibble
    as_tibble(result)
}

################################################################################
# 16. Loop Analysis Functions
################################################################################

#' Analyze loops by distance threshold
#' @param df Data frame with gene_id, loop.id, and distance columns
#' @param threshold_distance Maximum distance threshold
#' @param top_n_genes Number of top genes to plot
#' @param print_top_n Number of top genes to print
approach_2nd_analyze_loops_by_threshold <- function(df, threshold_distance = 2e5, top_n_genes = 80, print_top_n = 50) {
    # Task 1: Filtering and top N genes by loop count
    df.filtered <- df %>% filter(distance <= threshold_distance)
    df.gene_loop_count <- df.filtered %>%
        # filter(UP == "OK" | DOWN == "OK") %>%
        # filter(!is.na(gene_start)) %>% # in order to exclude NA values
        count(gene_id, sort = TRUE) # descend gene_id n 

    df.top_genes <- df.gene_loop_count %>%
        slice_max(n, n = top_n_genes) %>%
        arrange(n)

    plot <- ggplot(df.top_genes, aes(x = reorder(gene_id, n), y = n)) +
        geom_bar(stat = "identity", fill = "#1F78B4") +
        geom_text(aes(label = n), vjust = -0.3, color = "red", size = 3) +
        coord_flip() +
        labs(
            title = paste("Top", top_n_genes, "Genes by Number of Loops (distance <=", threshold_distance, ")"),
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
    slice_max(n, n = print_top_n)

    cat("Top", print_top_n, "genes:\n")
    print(top_genes, n = Inf)    

}

#' Plot and save filtering summary histogram
#' @param df_counts Data frame with count column
#' @param count_col Name of count column
#' @param output_prefix Prefix for output filenames
plotting_and_filtering_summary <- function(df_counts, count_col, output_prefix) {
    # 1. histogram
    p.hist <- ggplot(df_counts, aes(x = log2(.data[[count_col]]))) +
        geom_histogram() +
        facet_wrap(~ WHERE + resolution, scales = "free_y") +
        labs(
            y = "Count"
        )

    # 2. PDF
    pdf_path <- paste0("./figures/submission/lt2mb/histogram_number_of_", output_prefix, "_within_ends_of_loops_by_resolution_end.pdf")
    pdf(file = pdf_path, width = 11 * 0.8, height = 8.5 * 0.8)
    print(p.hist)
    dev.off()

    # 3. png
    png_path <- paste0("./figures/submission/lt2mb/histogram_number_of_", output_prefix, "_within_ends_of_loops_by_resolution_end.png")
    ggsave(filename = png_path, plot = p.hist, width = 11 * 0.8, height = 8.5 * 0.8, dpi = 300)
}

################################################################################
# 17. Circos Plotting Functions
################################################################################

#' Plot Circos diagram for a specific chromosome
#' @param chr Chromosome to plot
#' @note Requires df.circos.input.log.final.loop and resolution_colors in environment
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
            )
        }
        circos.clear()
    }
}

################################################################################
# 18. Distance Calculation Functions
################################################################################

#' Convert log2 value to actual distance
#' @param x Log2 value
#' @return Actual distance (2^x - 1)
get_distance <- function(x) {
    distance <- 2^x - 1
    return(distance)
}

################################################################################
# 19. Feature Binning Functions
################################################################################

#' Process feature bins for chromosome analysis
#' @param feature_df Data frame with feature information (chr, start columns)
#' @param chromosome_ends Data frame with chromosome end positions
#' @param bin_size Bin size in base pairs (default: 1e6)
#' @param label Label for the feature type
#' @return Data frame with binned feature counts
process_feature_bins <- function(feature_df, chromosome_ends, bin_size = 1e6, label = "Feature") {
    # col type
    chromosome_ends <- chromosome_ends %>% mutate(chr = str_remove(as.character(chr), "chr"))
    feature_df$chr <- as.character(feature_df$chr)

    result <- list()

    valid_chromosomes <- as.character(c(1:20, "X", "Y"))
    feature_df <- feature_df %>%
        mutate(chr = str_remove(chr, "chr")) %>%
        filter(chr %in% valid_chromosomes)

    for (current_chr in unique(feature_df$chr)) {
        chr_data <- feature_df %>% filter(chr == current_chr)
        chr_end <- chromosome_ends %>%
            filter(chr == current_chr) %>%
            pull(end)

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
            end = bin_edges[-1] - 1
        )

        if (tail(bins$end, 1) < chr_end) {
            bins <- rbind(bins, data.frame(start = bins$end[nrow(bins)] + 1, end = chr_end))
        }

        bins <- bins %>%
            mutate(
                chr = current_chr,
                Value = 0,
                Feature = label
            )

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

################################################################################
# 20. Dataframe comparison Functions
################################################################################

compare_df <- function(df1, df2) {
  cat("df1 rows:", nrow(df1), "\n")
  cat("df2 rows:", nrow(df2), "\n")
  
  if (nrow(df1) != nrow(df2)) {
    cat("Different number of rows!\n")
    return(FALSE)
  }
  
  # 정렬 후 비교
  df1_sorted <- df1 %>% arrange(across(everything()))
  df2_sorted <- df2 %>% arrange(across(everything()))
  
  result <- identical(df1_sorted, df2_sorted)
  cat("Identical:", result, "\n")
  
  return(result)
}

message("All utility functions loaded successfully!")


