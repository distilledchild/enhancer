library(fs)

################################################################################
# Shared Google Drive project paths
################################################################################

# Return the path of the R script currently being executed or sourced.
current_script_path <- function() {
  file.args <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )
  if (length(file.args) > 0L) {
    script.arg <- sub("^--file=", "", file.args[[1]])
    script.arg <- gsub("~+~", " ", script.arg, fixed = TRUE)
    return(normalizePath(
      script.arg,
      winslash = "/",
      mustWork = FALSE
    ))
  }

  frame.files <- vapply(
    sys.frames(),
    function(frame) {
      if (is.null(frame$ofile)) NA_character_ else as.character(frame$ofile)
    },
    character(1)
  )
  frame.files <- frame.files[!is.na(frame.files) & nzchar(frame.files)]
  if (length(frame.files) > 0L) {
    return(normalizePath(
      tail(frame.files, 1L),
      winslash = "/",
      mustWork = FALSE
    ))
  }
  NA_character_
}

# Resolve the shared Google Drive project and all revision-main directories.
resolve_enhancer_analysis_paths <- function(
  funcs.file,
  analysis.relative.path = file.path("revision", "revision_main")
) {
  funcs.file <- normalizePath(
    path.expand(funcs.file),
    winslash = "/",
    mustWork = TRUE
  )

  explicit.project.dir <- Sys.getenv("ENHANCER_PROJECT_DIR", unset = "")
  explicit.r.files.dir <- Sys.getenv("ENHANCER_R_FILES_DIR", unset = "")
  project.candidates <- unique(c(
    explicit.project.dir,
    if (nzchar(explicit.r.files.dir)) dirname(explicit.r.files.dir) else "",
    dirname(dirname(funcs.file)),
    Sys.glob(path.expand(
      "~/Library/CloudStorage/GoogleDrive-*/My Drive/research/enhancer"
    )),
    path.expand("~/Google Drive/My Drive/research/enhancer"),
    path.expand("~/Desktop/playground/enhancer")
  ))
  project.candidates <- project.candidates[
    nzchar(project.candidates) &
      file.exists(file.path(project.candidates, "r_files", "funcs.R"))
  ]
  if (length(project.candidates) == 0L) {
    stop(
      paste0(
        "Cannot locate the shared enhancer project. Make the Google Drive ",
        "folder available offline or set ENHANCER_PROJECT_DIR."
      ),
      call. = FALSE
    )
  }

  enhancer.project.dir <- normalizePath(
    project.candidates[[1]],
    winslash = "/",
    mustWork = TRUE
  )
  r.files.dir <- normalizePath(
    file.path(enhancer.project.dir, "r_files"),
    winslash = "/",
    mustWork = TRUE
  )
  analysis.dir <- path.expand(Sys.getenv(
    "RESUBMIT_ANALYSIS_DIR",
    unset = file.path(r.files.dir, analysis.relative.path)
  ))
  input.root <- path.expand(Sys.getenv(
    "RESUBMIT_INPUT_BUNDLE_DIR",
    unset = file.path(analysis.dir, "inputs")
  ))

  list(
    script.path = current_script_path(),
    enhancer.project.dir = enhancer.project.dir,
    r.files.dir = r.files.dir,
    revision.dir = file.path(r.files.dir, "revision"),
    analysis.dir = analysis.dir,
    input.root = input.root,
    data.dir = file.path(input.root, "data"),
    hiccups.loop.root = path.expand(Sys.getenv(
      "HICCUPS_LOOP_ROOT",
      unset = file.path(
        input.root,
        "hic",
        "2023A",
        "hic30_w_sb_options"
      )
    )),
    coord.cache.dir = file.path(analysis.dir, "cache_data"),
    output.dir = path.expand(Sys.getenv(
      "RESUBMIT_OUTPUT_DIR",
      unset = file.path(analysis.dir, "results")
    ))
  )
}

# func for str_split_n()
subset_safely <- function(x, index) {
  if (length(x) < index) {
    return(NA_character_)
  }
  x[[index]]
}

# func for str_split_n()
str_split_n <- function(string, pattern, n) {
  out <- str_split(string, pattern)
  vapply(out, subset_safely, character(1L), index = n)
}

str_split_first <- function(string, pattern) {
  str_split(string, pattern, simplify = TRUE)[, 1]
}

# Save one ggplot/patchwork object as publication-ready PDF and PNG files.
saving_plot_dual <- function(
  plot_obj,
  filename_base,
  output_dir = "./figures",
  width_in = 11,
  height_in = 8.5,
  scale_x = 1,
  scale_y = 1,
  dpi = 300,
  bg = "white"
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  final.width <- width_in * scale_x
  final.height <- height_in * scale_y

  ggplot2::ggsave(
    file.path(output_dir, paste0(filename_base, ".pdf")),
    plot = plot_obj,
    device = "pdf",
    width = final.width,
    height = final.height,
    units = "in",
    bg = bg
  )
  ggplot2::ggsave(
    file.path(output_dir, paste0(filename_base, ".png")),
    plot = plot_obj,
    device = "png",
    width = final.width,
    height = final.height,
    units = "in",
    dpi = dpi,
    bg = bg
  )

  invisible(file.path(
    output_dir,
    paste0(filename_base, c(".pdf", ".png"))
  ))
}

# file_list <- function(step, file_ext) {
#   # file.dir.5 <- path(step)
#   fs::dir_ls(path(step), regexp = str_c("\\.", file_ext, "$")
# }

# MiniMUGA output txt file data reader
txt_minimuga_data_reader <- function(file, colnames = FALSE) {
  # read_table2(file, skip = 1, col_names = colnames, guess_max = 2222)
  read_table2(file, skip = 1, col_names = colnames)
}

# MiniMUGA output csv file data reader (deprecated)
csv_minimuga_data_reader <- function(file) {
  # readr::read_csv(file, col_types = cols(.default = "c"))
  readr::read_csv(file, col_types = cols(chromosome = readr::col_character()))
}

# skip header part in a vcf file
vcf_header_skip <- function(file, pattern) {
  max(grep(pattern, read_lines(file)))
}

# vcf file data reader
vcf_data_reader <- function(file, pattern) {
  read_tsv(file, skip = vcf_header_skip(file, "^##"), col_name = TRUE, guess_max = 2222) %>%
    rename_at(ncol(.), ~"default")
}

# vcf file data reader for initial data after GLNexus
vcf_raw_data_reader <- function(file, pattern) {
  read_tsv(file, skip = vcf_header_skip(file, "^##"), col_types = cols("#CHROM" = readr::col_character()))
  # read_tsv(file, skip = vcf_header_skip(file, '^##'), col_name = TRUE, col_types = "cdcccdccccccccc")
}

# vcf file data reader for initial data splitted by chrs
vcf_by_chr_raw_data_reader <- function(file, pattern) {
  read_tsv(file, skip = vcf_header_skip(file, "^##"), col_name = TRUE, col_types = cols(.default = "d", pos = "d"))
}

# inital data of postition 151 provider
# params: file_list, deli_dir, meta_data
# file_list: file list
# deli_dir: dir name for delimiter
# meta_data: meta_data
initial_data_provider_151 <- function(file_list, deli_dir, meta_data) {
  file_list %>%
    map_dfr(vcf_data_reader, .id = "sample") %>%
    rename_at(vars(starts_with("#")), funs(str_replace(., "#", ""))) %>%
    filter(POS == 151 & FILTER != "RefCall") %>% # only 151 important
    mutate(CHROM = str_to_upper(CHROM)) %>%
    filter(str_detect(CHROM, regex("Q30", ignore_case = T))) %>% # only Q30 important
    filter(str_length(ALT) == 1 & str_length(REF) == 1) %>% # redundant condition
    filter(str_detect(CHROM, regex("snp", ignore_case = T))) %>% # only SNP markers are important
    # # mutate(sample = str_split_n(str_split_n(sample, str_c('\\/', deli_dir, '\\/'), 2), '\\.', 1)) %>%
    mutate(sample = str_split_n(sample, str_c(deli_dir, "\\/"), 2)) %>%
    mutate(sample = if_else(str_detect(sample, "_"), str_split_n(sample, "\\_", 1), str_split_n(sample, "\\.", 1))) %>%
    mutate(GT = if_else(str_detect(default, ":"), str_split_n(default, ":", 1), default)) %>%
    mutate(DP = if_else(str_detect(INFO, "="), as.integer(str_split_n(str_split_n(INFO, "DP=", 2), ";", 1)), as.integer(str_split_n(default, ":", 3)))) %>%
    select(-c(FORMAT, INFO)) %>%
    filter(!is.na(DP)) %>%
    mutate(number = str_split_n(sample, "Rat", 2)) %>%
    left_join(meta_data, by = c("number" = "number")) %>%
    mutate(column_name = paste0(.$sample, ".", .$strain, ".", .$sex)) %>%
    separate(CHROM, into = c("mq", "mstrain", "mvt", "mcoord"), sep = "[|_]", remove = FALSE) # prefix 'm' means 'marker' in order not to be confused with q, strain, variant type of data
}

# getting sequence region
# params: chrom, start, end, species
# TODO using ... for optional parameter
marker_seq_fetcher <- function(chrom, start, end, species, mask = NULL) {
  server <- "https://rest.ensembl.org"
  ext <- "/sequence/region/"

  r <- str_c(server, ext, species, "/", chrom, ":", start, "..", end, ":?", if_else(is.null(mask), "", str_c("mask=", mask))) %>%
    print() %>%
    GET(content_type("text/plain")) %>%
    content()

  return(r)
}

# setting markers across chromosome with coordinates
# using if-else for more readable
# should be used with purrr::map
# very strict input dataset (column names)
# TODO change column names to column index
marker_coord_setter <- function(...) {
  current.row.df <- tibble(...)

  if ((current.row.df$max_val - current.row.df$min_val) == 300) {
    # just one marker in a chrom

    current.row.df %>%
      mutate(Start = Init) %>%
      mutate(End = current.row.df$Start - 1) %>%
      bind_rows(current.row.df) %>%
      bind_rows(current.row.df %>% mutate(Start = current.row.df$End + 1) %>% mutate(End = current.row.df$Terminal))
  } else if (current.row.df$Start == current.row.df$min_val) {
    # 1st marker in a chrom

    current.row.df %>%
      mutate(Start = Init) %>%
      mutate(End = current.row.df$Start - 1) %>%
      bind_rows(current.row.df)
  } else if (current.row.df$End == current.row.df$max_val) {
    # the last marker in a chrom

    current.row.df %>%
      mutate(Start = NA) %>%
      mutate(End = current.row.df$End - 1) %>%
      bind_rows(current.row.df) %>%
      bind_rows(current.row.df %>% mutate(Start = current.row.df$End + 1) %>% mutate(End = current.row.df$Terminal))
  } else {
    # a marker
    current.row.df %>%
      mutate(Start = NA) %>%
      mutate(End = current.row.df$End - 1) %>%
      bind_rows(current.row.df)
  }
}

# convertor to XStringSet
converting_df_2_XStringSet <- function(dfx, type) {
  column_names <- colnames(dfx)

  if (!("name" %in% column_names) || !("sequence" %in% column_names)) {
    print("name or sequence column absent")
  } else {
    for (cname in column_names) {
      assign(cname, dfx[[cname]])
    }
  }

  names(sequence) <- name
  return(ifelse(toupper(type) == "D", DNAStringSet(sequence),
    ifelse(toupper(type) == "R", RNAStringSet(sequence), AAStringSet(sequence))
  ))
}

my_quantile <- function(x, probs) {
  tibble(qt = quantile(x, probs), probs = probs)
}


bedpe_data_reader <- function(file) {
  read.csv(file, header = T, sep = "\t")
}

init.bedpe.df <- function(file.list, deli_dir) {
  file.list %>%
    map_dfr(bedpe_data_reader, .id = "sample") %>%
    mutate(sample = str_split_n(sample, str_c(deli_dir, "\\/"), 2)) %>%
    mutate(sample = str_split_n(sample, "_", 1)) %>%
    mutate(distance = ((y2 + y1) - (x2 + x1)) / 2)
  # %>%
  #   view()
}

################################################################################
# 1. resubmission
################################################################################

# Confirm that an input table contains every column required by a downstream
# analysis step. Stops with an object-specific error listing missing columns.
check_required_columns <- function(df, required.columns, object.name) {
  missing.columns <- setdiff(required.columns, colnames(df))
  if (length(missing.columns) > 0) {
    stop(
      object.name,
      " is missing required column(s): ",
      paste(missing.columns, collapse = ", "),
      call. = FALSE
    )
  }
}

# Convert coordinate-like input to integer coordinates after checking for
# missing, non-finite, fractional, or out-of-range values. Factors are first
# converted to character so their level indices are never mistaken for bases.
as_integer_coordinate <- function(x, object.name) {
  if (is.factor(x)) {
    x <- as.character(x)
  }

  x.numeric <- suppressWarnings(as.numeric(x))

  invalid.coordinate <- (
    is.na(x.numeric) |
      !is.finite(x.numeric) |
      x.numeric != floor(x.numeric) |
      x.numeric > .Machine$integer.max
  )

  if (any(invalid.coordinate)) {
    stop(
      object.name,
      " contains a missing, non-integer, or out-of-range coordinate.",
      call. = FALSE
    )
  }

  as.integer(x.numeric)
}

# Validate intervals represented in BED convention: zero-based start and
# half-open end. This checks start >= 0 and end > start without transforming
# the coordinates. Returns invisibly when all intervals are valid.
validate_bed_half_open_intervals <- function(
  start,
  end,
  object.name
) {
  start <- as_integer_coordinate(start, paste0(object.name, " start"))
  end <- as_integer_coordinate(end, paste0(object.name, " end"))

  invalid.interval <- start < 0L | end <= start
  if (any(invalid.interval)) {
    stop(
      object.name,
      " must use valid 0-based, half-open intervals with ",
      "start >= 0 and end > start.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

# Validate intervals represented in the analysis convention: one-based start
# and inclusive end. This checks start >= 1 and end >= start and returns
# invisibly when all intervals are valid.
validate_one_based_inclusive_intervals <- function(
  start,
  end,
  object.name
) {
  start <- as_integer_coordinate(start, paste0(object.name, " start"))
  end <- as_integer_coordinate(end, paste0(object.name, " end"))

  invalid.interval <- start < 1L | end < start
  if (any(invalid.interval)) {
    stop(
      object.name,
      " must use valid 1-based, inclusive intervals with ",
      "start >= 1 and end >= start.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

# Convert a zero-based BED start to a one-based inclusive start. BED end
# coordinates remain unchanged and are therefore not accepted by this helper.
bed_start_to_one_based <- function(x, object.name) {
  x <- as_integer_coordinate(x, object.name)
  if (any(x < 0L)) {
    stop(object.name, " contains a negative BED start.", call. = FALSE)
  }
  x + 1L
}

# Standardise numeric and character resolution labels to 5K, 10K, or 25K.
# Unknown values are preserved so callers can detect or handle them explicitly.
normalise_resolution <- function(x) {
  case_when(
    as.character(x) %in% c("5000", "5K", "5k") ~ "5K",
    as.character(x) %in% c("10000", "10K", "10k") ~ "10K",
    as.character(x) %in% c("25000", "25K", "25k") ~ "25K",
    TRUE ~ as.character(x)
  )
}

# Return the resolution-adjusted minimum CTCF-motif count used in sensitivity
# analysis: 6 for 5K, 12 for 10K, and 30 for 25K loop anchors.
ctcf_resolution_adjusted_threshold <- function(x) {
  case_when(
    normalise_resolution(x) == "5K" ~ 6L,
    normalise_resolution(x) == "10K" ~ 12L,
    normalise_resolution(x) == "25K" ~ 30L,
    TRUE ~ NA_integer_
  )
}

# Count values that are explicitly TRUE while treating FALSE and missing values
# as not passing the criterion.
count_true <- function(x) {
  sum(x %in% TRUE, na.rm = TRUE)
}

# Calculate the percentage of values explicitly equal to TRUE. The result is
# rounded to the requested number of digits; an empty input returns NA.
percent_true <- function(x, digits = 1) {
  if (length(x) == 0) {
    return(NA_real_)
  }
  round(100 * mean(x %in% TRUE, na.rm = TRUE), digits)
}

# Divide a numerator by a valid positive denominator and otherwise return NA.
# Vectorised with dplyr::if_else for use inside mutate operations.
safe_ratio <- function(numerator, denominator) {
  if_else(
    !is.na(denominator) & denominator > 0,
    numerator / denominator,
    NA_real_
  )
}

################################################################################
# Revision: coordinate-preparation cache functions
################################################################################

# Return the complete set of normalized coordinate and coordinate-QC objects stored in the revision cache.
coordinate_cache_object_names <- function() {
  c(
    "df.analysis.input.files",
    "df.hic.source.provenance",
    "loop.file.metadata",
    "df.sample.loop.1based",
    "df.loop.pooled.support.summary",
    "df.loop.distinct",
    "df.loop.universe",
    "df.loop.source.count.check",
    "gr.ctcf.motif",
    "df.ctcf.fimo.summary",
    "df.transcript.ensembl.rn7.1based",
    "df.promoter.epd.rn7.1based",
    "gr.atac.rn7.1based",
    "df.coordinate.system.audit"
  )
}

# Build a named file-path vector for the RDS files that make up the coordinate-preparation cache.
coordinate_cache_paths <- function(
  cache.dir,
  object.names = coordinate_cache_object_names()
) {
  setNames(
    file.path(cache.dir, paste0(object.names, ".rds")),
    object.names
  )
}

# Report whether every expected coordinate-cache RDS file is present.
coordinate_cache_complete <- function(
  cache.dir,
  object.names = coordinate_cache_object_names()
) {
  all(file.exists(coordinate_cache_paths(cache.dir, object.names)))
}

# Save normalized coordinate objects as separate RDS files and write a cache manifest.
save_coordinate_cache_objects <- function(
  cache.dir,
  object.names = coordinate_cache_object_names(),
  envir = parent.frame()
) {
  dir.create(cache.dir, recursive = TRUE, showWarnings = FALSE)

  missing.objects <- object.names[
    !vapply(
      object.names,
      exists,
      logical(1L),
      envir = envir,
      inherits = FALSE
    )
  ]
  if (length(missing.objects) > 0L) {
    stop(
      "Cannot create coordinate cache; missing object(s): ",
      paste(missing.objects, collapse = ", "),
      call. = FALSE
    )
  }

  cache.paths <- coordinate_cache_paths(cache.dir, object.names)
  cache.objects <- mget(object.names, envir = envir, inherits = FALSE)

  for (object.name in object.names) {
    cache.path <- cache.paths[[object.name]]
    temporary.path <- paste0(cache.path, ".tmp")
    saveRDS(cache.objects[[object.name]], temporary.path)
    if (!file.rename(temporary.path, cache.path)) {
      unlink(temporary.path)
      stop("Cannot finalize coordinate cache file: ", cache.path, call. = FALSE)
    }
  }

  df.cache.manifest <- tibble(
    object_name = object.names,
    cache_file = unname(cache.paths),
    object_class = vapply(
      cache.objects,
      function(x) paste(class(x), collapse = ";"),
      character(1L)
    ),
    n_records = vapply(
      cache.objects,
      function(x) {
        if (is.data.frame(x)) {
          nrow(x)
        } else {
          length(x)
        }
      },
      integer(1L)
    ),
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
  readr::write_tsv(
    df.cache.manifest,
    file.path(cache.dir, "coordinate_cache_manifest.tsv")
  )

  invisible(df.cache.manifest)
}

# Load every normalized coordinate-cache object into the requested analysis environment.
load_coordinate_cache_objects <- function(
  cache.dir,
  object.names = coordinate_cache_object_names(),
  envir = parent.frame()
) {
  cache.paths <- coordinate_cache_paths(cache.dir, object.names)
  missing.files <- cache.paths[!file.exists(cache.paths)]
  if (length(missing.files) > 0L) {
    stop(
      "Coordinate cache is incomplete; missing file(s): ",
      paste(unname(missing.files), collapse = ", "),
      call. = FALSE
    )
  }

  for (object.name in object.names) {
    assign(
      object.name,
      readRDS(cache.paths[[object.name]]),
      envir = envir
    )
  }

  invisible(cache.paths)
}

# Summarise a correlation using complete finite pairs. Returns sample size,
# estimate, and p-value; constant vectors or fewer than three pairs return NA
# statistics rather than causing the analysis to fail.
safe_correlation_summary <- function(x, y, method = "pearson") {
  valid.rows <- complete.cases(x, y) & is.finite(x) & is.finite(y)
  x <- x[valid.rows]
  y <- y[valid.rows]

  if (
    length(x) < 3 ||
      isTRUE(all.equal(stats::sd(x), 0)) ||
      isTRUE(all.equal(stats::sd(y), 0))
  ) {
    return(tibble(
      n = length(x),
      estimate = NA_real_,
      p_value = NA_real_
    ))
  }

  test.result <- suppressWarnings(
    stats::cor.test(x, y, method = method, exact = FALSE)
  )

  tibble(
    n = length(x),
    estimate = unname(test.result$estimate),
    p_value = test.result$p.value
  )
}

# Fit loop count against log10 Hi-C contacts and append fitted loop counts and
# residuals. These residuals support exploratory depth sensitivity checks, not
# biological strain comparisons. Insufficient or invariant data yield NA.
add_depth_residuals <- function(df) {
  df <- df %>%
    mutate(log10_hic_contacts = log10(`Hi-C_Contacts`))

  valid.rows <- complete.cases(df$n_loops, df$log10_hic_contacts) &
    is.finite(df$n_loops) &
    is.finite(df$log10_hic_contacts)

  if (
    sum(valid.rows) < 3 ||
      isTRUE(all.equal(stats::sd(df$n_loops[valid.rows]), 0)) ||
      isTRUE(all.equal(stats::sd(df$log10_hic_contacts[valid.rows]), 0))
  ) {
    return(
      df %>%
        mutate(
          fitted_loop_count_by_hic_contacts = NA_real_,
          depth_residual_loop_count = NA_real_
        )
    )
  }

  depth.model <- stats::lm(
    n_loops ~ log10_hic_contacts,
    data = df[valid.rows, , drop = FALSE]
  )

  fitted.values <- rep(NA_real_, nrow(df))
  fitted.values[valid.rows] <- as.numeric(stats::predict(
    depth.model,
    newdata = df[valid.rows, , drop = FALSE]
  ))

  df %>%
    mutate(
      fitted_loop_count_by_hic_contacts = fitted.values,
      depth_residual_loop_count = n_loops - fitted_loop_count_by_hic_contacts
    )
}

# Read a ten-column narrowPeak file, validate its BED coordinates, convert only
# its starts from zero-based to one-based, and return a one-based inclusive
# GRanges object carrying the available peak statistics.
read_atac_narrowpeak <- function(file, sample.label = "Duttke2022_snATAC_PFC") {
  df.atac <- readr::read_tsv(
    file,
    col_names = c(
      "chr", "start", "end", "name", "score", "strand",
      "fc", "neglog10p", "neglog10q", "summit"
    ),
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  ) %>%
    mutate(
      start = as.integer(start),
      end = as.integer(end),
      score = suppressWarnings(as.numeric(score)),
      fc = suppressWarnings(as.numeric(fc)),
      neglog10q = suppressWarnings(as.numeric(neglog10q))
    )

  validate_bed_half_open_intervals(
    df.atac$start,
    df.atac$end,
    "ATAC narrowPeak intervals"
  )

  atac.start.one.based <- bed_start_to_one_based(
    df.atac$start,
    "ATAC narrowPeak start"
  )

  validate_one_based_inclusive_intervals(
    atac.start.one.based,
    df.atac$end,
    "Normalised ATAC intervals"
  )

  GRanges(
    seqnames = df.atac$chr,
    ranges = IRanges(
      start = atac.start.one.based,
      end = df.atac$end
    ),
    score = df.atac$score,
    fc = df.atac$fc,
    neglog10q = df.atac$neglog10q,
    sample = sample.label
  )
}

# Read full genomic CTCF FIMO output, retain the most significant prediction and
# multiplicity metadata for each exact genomic interval, validate one-based
# inclusive coordinates, and return a coordinate-distinct GRanges object.
read_ctcf_fimo <- function(file) {
  df.ctcf <- readr::read_tsv(
    file,
    comment = "#",
    col_types = cols(
      motif_id = col_character(),
      motif_alt_id = col_character(),
      sequence_name = col_character(),
      start = col_integer(),
      stop = col_integer(),
      strand = col_character(),
      score = col_double(),
      `p-value` = col_double(),
      `q-value` = col_double(),
      matched_sequence = col_character()
    ),
    show_col_types = FALSE
  )

  check_required_columns(
    df.ctcf,
    c(
      "motif_id", "motif_alt_id", "sequence_name", "start", "stop",
      "strand", "score", "p-value", "q-value", "matched_sequence"
    ),
    "CTCF FIMO table"
  )

  df.ctcf <- df.ctcf %>%
    transmute(
      motif_id = as.character(motif_id),
      motif_alt_id = as.character(motif_alt_id),
      chr = as.character(sequence_name),
      start = as_integer_coordinate(start, "CTCF FIMO start"),
      end = as_integer_coordinate(stop, "CTCF FIMO stop"),
      strand = as.character(strand),
      score = as.numeric(score),
      p_value = as.numeric(`p-value`),
      q_value = as.numeric(`q-value`),
      matched_sequence = as.character(matched_sequence)
    )

  if (
    any(is.na(df.ctcf$motif_id)) ||
      any(is.na(df.ctcf$score)) ||
      any(is.na(df.ctcf$p_value)) ||
      any(is.na(df.ctcf$q_value))
  ) {
    stop(
      "Full CTCF FIMO output contains missing motif IDs or score/p/q values.",
      call. = FALSE
    )
  }

  validate_one_based_inclusive_intervals(
    df.ctcf$start,
    df.ctcf$end,
    "CTCF FIMO intervals"
  )

  n.raw.predictions <- nrow(df.ctcf)
  n.raw.motif.ids <- n_distinct(df.ctcf$motif_id)

  # Sort each exact coordinate by p-value, q-value, descending score, motif ID,
  # and strand so selection of the cached representative row is deterministic.
  df.ctcf <- df.ctcf %>%
    arrange(chr, start, end, p_value, q_value, desc(score), motif_id, strand) %>%
    group_by(chr, start, end) %>%
    summarise(
      best_motif_id = dplyr::first(motif_id),
      best_motif_alt_id = dplyr::first(motif_alt_id),
      best_strand = dplyr::first(strand),
      best_score = dplyr::first(score),
      best_p_value = dplyr::first(p_value),
      best_q_value = dplyr::first(q_value),
      best_matched_sequence = dplyr::first(matched_sequence),
      n_fimo_predictions = n(),
      n_distinct_motif_ids = n_distinct(motif_id),
      n_distinct_strands = n_distinct(strand),
      predicted_on_both_strands = n_distinct(strand) > 1L,
      .groups = "drop"
    ) %>%
    mutate(
      ctcf_position = as.integer(round((start + end) / 2)),
      ctcf_id = str_c(chr, "_", start, "_", end, "_", ctcf_position)
    )

  df.ctcf.fimo.summary <- tibble(
    source_file = normalizePath(file, mustWork = TRUE),
    raw_fimo_predictions = n.raw.predictions,
    distinct_motif_ids = n.raw.motif.ids,
    exact_coordinate_distinct_intervals = nrow(df.ctcf),
    intervals_with_multiple_predictions = sum(
      df.ctcf$n_fimo_predictions > 1L
    ),
    intervals_with_multiple_motif_ids = sum(
      df.ctcf$n_distinct_motif_ids > 1L
    ),
    intervals_predicted_on_both_strands = sum(
      df.ctcf$predicted_on_both_strands
    ),
    coordinate_system = "FIMO_1_based_inclusive",
    cached_interval_definition = paste0(
      "Exact-coordinate interval; representative prediction selected by ",
      "minimum p-value, minimum q-value, maximum score, motif ID, and strand."
    )
  )

  gr.ctcf <- GRanges(
    seqnames = df.ctcf$chr,
    ranges = IRanges(
      start = df.ctcf$start,
      end = df.ctcf$end
    ),
    ctcf_id = df.ctcf$ctcf_id,
    ctcf_position = df.ctcf$ctcf_position,
    best_motif_id = df.ctcf$best_motif_id,
    best_motif_alt_id = df.ctcf$best_motif_alt_id,
    best_strand = df.ctcf$best_strand,
    best_score = df.ctcf$best_score,
    best_p_value = df.ctcf$best_p_value,
    best_q_value = df.ctcf$best_q_value,
    best_matched_sequence = df.ctcf$best_matched_sequence,
    n_fimo_predictions = df.ctcf$n_fimo_predictions,
    n_distinct_motif_ids = df.ctcf$n_distinct_motif_ids,
    n_distinct_strands = df.ctcf$n_distinct_strands,
    predicted_on_both_strands = df.ctcf$predicted_on_both_strands,
    source_coordinate_system = "FIMO_1_based_inclusive"
  )

  S4Vectors::metadata(gr.ctcf)$fimo_summary <- df.ctcf.fimo.summary
  gr.ctcf
}

################################################################################
# Revision: Hi-C loop-processing functions
################################################################################

# Parse one required numeric HiCCUPS field and fail on non-numeric source values.
parse_hiccups_numeric_field <- function(x, field.name) {
  x <- as.character(x)
  parsed <- suppressWarnings(as.numeric(x))
  invalid <- !is.na(x) & nzchar(x) & is.na(parsed)

  if (any(invalid)) {
    stop(
      "HiCCUPS field ", field.name,
      " contains non-numeric value(s): ",
      paste(head(unique(x[invalid]), 5L), collapse = ", "),
      call. = FALSE
    )
  }

  parsed
}

# Read one Juicer HiCCUPS BEDPE file and retain coordinates, provenance, and loop-quality fields.
read_hiccups_bedpe <- function(file, sample, strain) {
  header.line <- tryCatch(
    readLines(file, n = 1L, warn = FALSE),
    error = function(e) {
      stop(
        "Cannot read original HiCCUPS file: ", file, "\n",
        "If it is stored in Google Drive, select Available offline.\n",
        "Original error: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )
  column.names <- str_split(
    str_remove(header.line, "^#"),
    "\\t",
    simplify = TRUE
  ) %>%
    as.character()

  df.loop <- readr::read_tsv(
    file,
    comment = "#",
    col_names = column.names,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  )

  check_required_columns(
    df.loop,
    c(
      "chr1", "x1", "x2", "chr2", "y1", "y2",
      "name", "score", "strand1", "strand2", "color",
      "observed", "expectedBL", "expectedDonut", "expectedH", "expectedV",
      "fdrBL", "fdrDonut", "fdrH", "fdrV", "numCollapsed",
      "centroid1", "centroid2", "radius"
    ),
    paste0("HiCCUPS BEDPE file ", basename(file))
  )

  df.loop <- df.loop %>%
    transmute(
      sample = sample,
      strain = strain,
      source_file = file,
      chr1 = as.character(chr1),
      x1 = as_integer_coordinate(x1, "HiCCUPS x1"),
      x2 = as_integer_coordinate(x2, "HiCCUPS x2"),
      chr2 = as.character(chr2),
      y1 = as_integer_coordinate(y1, "HiCCUPS y1"),
      y2 = as_integer_coordinate(y2, "HiCCUPS y2"),
      hiccups_name = as.character(name),
      hiccups_score = as.character(score),
      hiccups_strand1 = as.character(strand1),
      hiccups_strand2 = as.character(strand2),
      hiccups_color = as.character(color),
      hiccups_observed = parse_hiccups_numeric_field(
        observed,
        "observed"
      ),
      hiccups_expected_bl = parse_hiccups_numeric_field(
        expectedBL,
        "expectedBL"
      ),
      hiccups_expected_donut = parse_hiccups_numeric_field(
        expectedDonut,
        "expectedDonut"
      ),
      hiccups_expected_h = parse_hiccups_numeric_field(
        expectedH,
        "expectedH"
      ),
      hiccups_expected_v = parse_hiccups_numeric_field(
        expectedV,
        "expectedV"
      ),
      hiccups_fdr_bl = parse_hiccups_numeric_field(fdrBL, "fdrBL"),
      hiccups_fdr_donut = parse_hiccups_numeric_field(
        fdrDonut,
        "fdrDonut"
      ),
      hiccups_fdr_h = parse_hiccups_numeric_field(fdrH, "fdrH"),
      hiccups_fdr_v = parse_hiccups_numeric_field(fdrV, "fdrV"),
      hiccups_num_collapsed = as_integer_coordinate(
        numCollapsed,
        "HiCCUPS numCollapsed"
      ),
      hiccups_centroid1_source = as_integer_coordinate(
        centroid1,
        "HiCCUPS centroid1"
      ),
      hiccups_centroid2_source = as_integer_coordinate(
        centroid2,
        "HiCCUPS centroid2"
      ),
      hiccups_radius_bp = as_integer_coordinate(
        radius,
        "HiCCUPS radius"
      )
    )

  validate_bed_half_open_intervals(
    df.loop$x1,
    df.loop$x2,
    paste0("HiCCUPS anchor1 intervals in ", basename(file))
  )
  validate_bed_half_open_intervals(
    df.loop$y1,
    df.loop$y2,
    paste0("HiCCUPS anchor2 intervals in ", basename(file))
  )

  df.loop %>%
    mutate(
      resolution_bp = x2 - x1,
      resolution = normalise_resolution(resolution_bp),
      loop_distance = as.integer(((y1 + y2) - (x1 + x2)) / 2),
      loop_id = str_c(
        chr1, x1, x2, chr2, y1, y2, resolution_bp,
        sep = "_"
      ),
      sample_loop_id = str_c(strain, loop_id, sep = "_")
    )
}

# Read all revision HiCCUPS files listed in metadata and validate the combined sample-level loop table.
read_hiccups_loop_files <- function(loop.file.metadata) {
  check_required_columns(
    loop.file.metadata,
    c("sample", "strain", "file"),
    "loop.file.metadata"
  )

  df.sample.loop.raw <- loop.file.metadata %>%
    dplyr::select(sample, strain, file) %>%
    pmap_dfr(
      function(sample, strain, file) {
        read_hiccups_bedpe(
          file = file,
          sample = sample,
          strain = strain
        )
      }
    )

  check_required_columns(
    df.sample.loop.raw,
    c(
      "sample", "strain", "source_file", "sample_loop_id", "loop_id",
      "chr1", "x1", "x2", "chr2", "y1", "y2",
      "resolution", "resolution_bp", "loop_distance"
    ),
    "df.sample.loop.raw"
  )

  validate_bed_half_open_intervals(
    df.sample.loop.raw$x1,
    df.sample.loop.raw$x2,
    "Raw HiCCUPS anchor1 intervals"
  )
  validate_bed_half_open_intervals(
    df.sample.loop.raw$y1,
    df.sample.loop.raw$y2,
    "Raw HiCCUPS anchor2 intervals"
  )

  if (
    nrow(df.sample.loop.raw) !=
      n_distinct(df.sample.loop.raw$sample_loop_id)
  ) {
    stop(
      "Original HiCCUPS files contain duplicate loop IDs within a sample.",
      call. = FALSE
    )
  }

  df.sample.loop.raw
}

# Convert HiCCUPS anchors and centroids to 1-based coordinates while preserving source quality fields.
normalize_hiccups_loop_coordinates <- function(
  df.sample.loop.raw,
  max.loop.distance = 2000000L
) {
  max.loop.distance <- as_integer_coordinate(
    max.loop.distance,
    "Maximum loop distance"
  )
  if (length(max.loop.distance) != 1L || max.loop.distance <= 0L) {
    stop("Maximum loop distance must be one positive integer.", call. = FALSE)
  }

  df.sample.loop.1based <- df.sample.loop.raw %>%
    transmute(
      sample,
      strain,
      source_file,
      sample_loop_id,
      loop_id,
      chr1 = as.character(chr1),
      start1 = bed_start_to_one_based(x1, "HiCCUPS anchor1 start"),
      end1 = as_integer_coordinate(x2, "HiCCUPS anchor1 end"),
      chr2 = as.character(chr2),
      start2 = bed_start_to_one_based(y1, "HiCCUPS anchor2 start"),
      end2 = as_integer_coordinate(y2, "HiCCUPS anchor2 end"),
      resolution = normalise_resolution(resolution),
      resolution_bp = as_integer_coordinate(
        resolution_bp,
        "HiCCUPS resolution"
      ),
      loop_distance = as_integer_coordinate(
        loop_distance,
        "HiCCUPS loop distance"
      ),
      passes_lt2mb = loop_distance < max.loop.distance,
      hiccups_name,
      hiccups_score,
      hiccups_strand1,
      hiccups_strand2,
      hiccups_color,
      hiccups_observed,
      hiccups_expected_bl,
      hiccups_expected_donut,
      hiccups_expected_h,
      hiccups_expected_v,
      hiccups_fdr_bl,
      hiccups_fdr_donut,
      hiccups_fdr_h,
      hiccups_fdr_v,
      hiccups_num_collapsed,
      hiccups_centroid1_source,
      hiccups_centroid2_source,
      hiccups_centroid1_1based = hiccups_centroid1_source + 1L,
      hiccups_centroid2_1based = hiccups_centroid2_source + 1L,
      hiccups_radius_bp
    )

  validate_one_based_inclusive_intervals(
    df.sample.loop.1based$start1,
    df.sample.loop.1based$end1,
    "Normalized HiCCUPS anchor1 intervals"
  )
  validate_one_based_inclusive_intervals(
    df.sample.loop.1based$start2,
    df.sample.loop.1based$end2,
    "Normalized HiCCUPS anchor2 intervals"
  )

  expected.anchor.width <- case_when(
    df.sample.loop.1based$resolution == "5K" ~ 5000L,
    df.sample.loop.1based$resolution == "10K" ~ 10000L,
    df.sample.loop.1based$resolution == "25K" ~ 25000L,
    TRUE ~ NA_integer_
  )
  observed.anchor1.width <- df.sample.loop.1based$end1 -
    df.sample.loop.1based$start1 + 1L
  observed.anchor2.width <- df.sample.loop.1based$end2 -
    df.sample.loop.1based$start2 + 1L

  if (
    any(is.na(expected.anchor.width)) ||
      any(observed.anchor1.width != expected.anchor.width) ||
      any(observed.anchor2.width != expected.anchor.width)
  ) {
    stop(
      "Normalized loop-anchor widths do not match HiCCUPS resolution.",
      call. = FALSE
    )
  }

  df.sample.loop.1based
}

# Summarize source-library support and descriptive HiCCUPS fields for each exact pooled loop call.
summarise_pooled_hiccups_support <- function(
  df.sample.loop.1based
) {
  df.sample.loop.1based %>%
    group_by(loop_id) %>%
    summarise(
      n_supporting_libraries = n_distinct(sample),
      n_supporting_strains = n_distinct(strain),
      supporting_samples = str_c(sort(unique(sample)), collapse = ";"),
      supporting_strains = str_c(sort(unique(strain)), collapse = ";"),
      hiccups_observed_min = min(hiccups_observed),
      hiccups_observed_median = median(hiccups_observed),
      hiccups_observed_max = max(hiccups_observed),
      hiccups_expected_bl_median = median(hiccups_expected_bl),
      hiccups_expected_donut_median = median(hiccups_expected_donut),
      hiccups_expected_h_median = median(hiccups_expected_h),
      hiccups_expected_v_median = median(hiccups_expected_v),
      hiccups_fdr_bl_min = min(hiccups_fdr_bl),
      hiccups_fdr_donut_min = min(hiccups_fdr_donut),
      hiccups_fdr_h_min = min(hiccups_fdr_h),
      hiccups_fdr_v_min = min(hiccups_fdr_v),
      hiccups_num_collapsed_median = median(hiccups_num_collapsed),
      hiccups_num_collapsed_max = max(hiccups_num_collapsed),
      hiccups_centroid1_1based_median = median(hiccups_centroid1_1based),
      hiccups_centroid2_1based_median = median(hiccups_centroid2_1based),
      hiccups_radius_bp_median = median(hiccups_radius_bp),
      hiccups_radius_bp_max = max(hiccups_radius_bp),
      .groups = "drop"
    )
}

# Collapse sample-level rows into exact pooled loop calls without treating nearby calls as identical loci.
build_pooled_hiccups_loop_resource <- function(
  df.sample.loop.1based
) {
  df.loop.pooled.support.summary <- summarise_pooled_hiccups_support(
    df.sample.loop.1based
  )

  df.loop.distinct <-
    df.sample.loop.1based %>%
    dplyr::select(
      loop_id,
      chr1, start1, end1,
      chr2, start2, end2,
      resolution, resolution_bp,
      loop_distance, passes_lt2mb
    ) %>%
    distinct(loop_id, .keep_all = TRUE) %>%
    left_join(df.loop.pooled.support.summary, by = "loop_id") %>%
    arrange(chr1, start1, end1, chr2, start2, end2)

  df.loop.universe <- df.loop.distinct %>%
    filter(passes_lt2mb)

  if (nrow(df.loop.universe) != n_distinct(df.loop.universe$loop_id)) {
    stop("The pooled loop resource contains duplicate loop IDs.", call. = FALSE)
  }

  list(
    support = df.loop.pooled.support.summary,
    distinct = df.loop.distinct,
    universe = df.loop.universe
  )
}

# Assign stable connected-component labels to an edge list without an external graph dependency.
assign_connected_components_from_edges <- function(
  n.nodes,
  from,
  to
) {
  n.nodes <- as.integer(n.nodes)
  parent <- seq_len(n.nodes)

  find.root <- function(node) {
    while (parent[[node]] != node) {
      node <- parent[[node]]
    }
    node
  }

  if (length(from) > 0L) {
    for (edge.index in seq_along(from)) {
      root.from <- find.root(from[[edge.index]])
      root.to <- find.root(to[[edge.index]])
      if (root.from != root.to) {
        parent[[max(root.from, root.to)]] <- min(root.from, root.to)
      }
    }
  }

  roots <- vapply(seq_len(n.nodes), find.root, integer(1L))
  component.order <- tibble(node = seq_len(n.nodes), root = roots) %>%
    group_by(root) %>%
    summarise(first_node = min(node), .groups = "drop") %>%
    arrange(first_node) %>%
    mutate(component_id = row_number())

  component.order$component_id[match(roots, component.order$root)]
}

# Group nearby exact HiCCUPS calls into all-pairs-constrained canonical loci.
build_approximate_loop_loci <- function(
  df.loop,
  merge.distance.bp = c("5K" = 20000L, "10K" = 20000L, "25K" = 50000L)
) {
  check_required_columns(
    df.loop,
    c(
      "loop_id", "chr1", "start1", "end1", "chr2", "start2", "end2",
      "resolution", "resolution_bp", "hiccups_centroid1_1based_median",
      "hiccups_centroid2_1based_median", "hiccups_radius_bp_max",
      "n_supporting_libraries", "supporting_samples"
    ),
    "Exact pooled HiCCUPS loop resource"
  )

  merge.distance.bp <- as.integer(merge.distance.bp)
  names(merge.distance.bp) <- c("5K", "10K", "25K")
  if (any(is.na(merge.distance.bp)) || any(merge.distance.bp <= 0L)) {
    stop("HiCCUPS merge-distance values must be positive integers.", call. = FALSE)
  }

  df.loop.ordered <- df.loop %>%
    arrange(chr1, start1, end1, chr2, start2, end2, resolution_bp, loop_id) %>%
    mutate(
      loop_index = row_number(),
      source_support_count = .data$n_supporting_libraries,
      locus_centroid1 = coalesce(
        hiccups_centroid1_1based_median,
        (start1 + end1) / 2
      ),
      locus_centroid2 = coalesce(
        hiccups_centroid2_1based_median,
        (start2 + end2) / 2
      ),
      locus_merge_distance_bp = unname(merge.distance.bp[resolution])
    )

  if (any(is.na(df.loop.ordered$locus_merge_distance_bp))) {
    stop("A loop resolution lacks a canonical-locus merge distance.", call. = FALSE)
  }

  # Use a broad anchor-1 window to enumerate nearby candidates efficiently;
  # the exact two-dimensional threshold is applied below.
  maximum.merge.distance.bp <- max(merge.distance.bp)
  gr.anchor1 <- GRanges(
    seqnames = df.loop.ordered$chr1,
    ranges = IRanges(
      start = pmax(
        1L,
        floor(df.loop.ordered$locus_centroid1 - maximum.merge.distance.bp)
      ),
      end = ceiling(
        df.loop.ordered$locus_centroid1 + maximum.merge.distance.bp
      )
    )
  )
  anchor1.hits <- findOverlaps(gr.anchor1, gr.anchor1, type = "any")

  df.candidate.edge <- tibble(
    loop_index_1 = queryHits(anchor1.hits),
    loop_index_2 = subjectHits(anchor1.hits)
  ) %>%
    filter(loop_index_1 < loop_index_2) %>%
    transmute(
      loop_index_1,
      loop_index_2,
      loop_id_1 = df.loop.ordered$loop_id[loop_index_1],
      loop_id_2 = df.loop.ordered$loop_id[loop_index_2],
      resolution_1 = df.loop.ordered$resolution[loop_index_1],
      resolution_2 = df.loop.ordered$resolution[loop_index_2],
      resolution_bp_1 = df.loop.ordered$resolution_bp[loop_index_1],
      resolution_bp_2 = df.loop.ordered$resolution_bp[loop_index_2],
      supporting_samples_1 =
        df.loop.ordered$supporting_samples[loop_index_1],
      supporting_samples_2 =
        df.loop.ordered$supporting_samples[loop_index_2],
      chr2_1 = df.loop.ordered$chr2[loop_index_1],
      chr2_2 = df.loop.ordered$chr2[loop_index_2],
      start1_1 = df.loop.ordered$start1[loop_index_1],
      end1_1 = df.loop.ordered$end1[loop_index_1],
      start1_2 = df.loop.ordered$start1[loop_index_2],
      end1_2 = df.loop.ordered$end1[loop_index_2],
      start2_1 = df.loop.ordered$start2[loop_index_1],
      end2_1 = df.loop.ordered$end2[loop_index_1],
      start2_2 = df.loop.ordered$start2[loop_index_2],
      end2_2 = df.loop.ordered$end2[loop_index_2],
      centroid1_1 = df.loop.ordered$locus_centroid1[loop_index_1],
      centroid1_2 = df.loop.ordered$locus_centroid1[loop_index_2],
      centroid2_1 = df.loop.ordered$locus_centroid2[loop_index_1],
      centroid2_2 = df.loop.ordered$locus_centroid2[loop_index_2],
      merge_distance_1 =
        df.loop.ordered$locus_merge_distance_bp[loop_index_1],
      merge_distance_2 =
        df.loop.ordered$locus_merge_distance_bp[loop_index_2]
    ) %>%
    filter(
      chr2_1 == chr2_2
    ) %>%
    mutate(
      anchor1_overlap_bp =
        pmin(end1_1, end1_2) - pmax(start1_1, start1_2) + 1,
      anchor2_overlap_bp =
        pmin(end2_1, end2_2) - pmax(start2_1, start2_2) + 1,
      centroid1_difference_bp = abs(centroid1_1 - centroid1_2),
      centroid2_difference_bp = abs(centroid2_1 - centroid2_2),
      centroid_distance_2d_bp = sqrt(
        centroid1_difference_bp^2 + centroid2_difference_bp^2
      ),
      centroid_threshold_bp = pmax(merge_distance_1, merge_distance_2),
      shares_source_library = map2_lgl(
        supporting_samples_1,
        supporting_samples_2,
        function(samples.1, samples.2) {
          samples.1 <- str_split(samples.1, fixed(";"), simplify = FALSE)[[1]]
          samples.2 <- str_split(samples.2, fixed(";"), simplify = FALSE)[[1]]
          length(intersect(samples.1, samples.2)) > 0L
        }
      ),
      normalized_centroid_distance =
        centroid_distance_2d_bp / centroid_threshold_bp
    ) %>%
    filter(
      !shares_source_library,
      centroid_distance_2d_bp <= centroid_threshold_bp
    ) %>%
    arrange(loop_index_1, loop_index_2)

  candidate.component.id <- assign_connected_components_from_edges(
    n.nodes = nrow(df.loop.ordered),
    from = df.candidate.edge$loop_index_1,
    to = df.candidate.edge$loop_index_2
  )

  # Partition each candidate component so every pair within a final locus has
  # a qualifying edge. This prevents transitive A-B-C chaining from merging A
  # and C when they are too far apart or were resolved separately in one library.
  qualifying.edge.distance <- setNames(
    df.candidate.edge$normalized_centroid_distance,
    str_c(
      pmin(df.candidate.edge$loop_index_1, df.candidate.edge$loop_index_2),
      pmax(df.candidate.edge$loop_index_1, df.candidate.edge$loop_index_2),
      sep = ":"
    )
  )
  final.cluster.id <- integer(nrow(df.loop.ordered))
  next.cluster.id <- 0L

  for (component in sort(unique(candidate.component.id))) {
    component.members <- which(candidate.component.id == component)
    member.order <- df.loop.ordered %>%
      filter(loop_index %in% component.members) %>%
      arrange(
        dplyr::desc(source_support_count),
        resolution_bp,
        loop_id
      ) %>%
      pull(loop_index)

    component.clusters <- list()
    for (member in member.order) {
      eligible.cluster <- which(vapply(
        component.clusters,
        function(cluster.members) {
          pair.keys <- str_c(
            pmin(member, cluster.members),
            pmax(member, cluster.members),
            sep = ":"
          )
          all(pair.keys %in% names(qualifying.edge.distance))
        },
        logical(1)
      ))

      if (length(eligible.cluster) == 0L) {
        component.clusters[[length(component.clusters) + 1L]] <- member
      } else {
        mean.distance <- vapply(
          eligible.cluster,
          function(cluster.index) {
            cluster.members <- component.clusters[[cluster.index]]
            pair.keys <- str_c(
              pmin(member, cluster.members),
              pmax(member, cluster.members),
              sep = ":"
            )
            mean(qualifying.edge.distance[pair.keys])
          },
          numeric(1)
        )
        chosen.cluster <- eligible.cluster[[which.min(mean.distance)]]
        component.clusters[[chosen.cluster]] <- c(
          component.clusters[[chosen.cluster]],
          member
        )
      }
    }

    for (cluster.members in component.clusters) {
      next.cluster.id <- next.cluster.id + 1L
      final.cluster.id[cluster.members] <- next.cluster.id
    }
  }

  df.loop.locus.map <- df.loop.ordered %>%
    transmute(
      loop_id,
      resolution,
      resolution_bp,
      approximate_loop_locus_id = str_c(
        "approx_locus_",
        str_pad(final.cluster.id, width = 6L, pad = "0")
      )
    ) %>%
    add_count(
      approximate_loop_locus_id,
      name = "n_exact_loop_calls_in_locus"
    ) %>%
    group_by(approximate_loop_locus_id) %>%
    mutate(n_resolutions_in_locus = n_distinct(resolution)) %>%
    ungroup() %>%
    mutate(is_multi_call_locus = n_exact_loop_calls_in_locus > 1L) %>%
    arrange(loop_id)

  df.loop.locus.edge <- df.candidate.edge %>%
    dplyr::select(
      loop_id_1,
      loop_id_2,
      resolution_1,
      resolution_2,
      anchor1_overlap_bp,
      anchor2_overlap_bp,
      centroid1_difference_bp,
      centroid2_difference_bp,
      centroid_distance_2d_bp,
      centroid_threshold_bp
    ) %>%
    left_join(
      df.loop.locus.map %>%
        dplyr::select(
          loop_id_1 = loop_id,
          approximate_loop_locus_id_1 = approximate_loop_locus_id
        ),
      by = "loop_id_1"
    ) %>%
    left_join(
      df.loop.locus.map %>%
        dplyr::select(
          loop_id_2 = loop_id,
          approximate_loop_locus_id_2 = approximate_loop_locus_id
        ),
      by = "loop_id_2"
    ) %>%
    filter(
      approximate_loop_locus_id_1 == approximate_loop_locus_id_2
    ) %>%
    transmute(
      approximate_loop_locus_id = approximate_loop_locus_id_1,
      dplyr::across(-c(
        approximate_loop_locus_id_1,
        approximate_loop_locus_id_2
      ))
    )

  df.loop.with.locus <- df.loop.ordered %>%
    left_join(
      df.loop.locus.map %>%
        dplyr::select(loop_id, approximate_loop_locus_id),
      by = "loop_id"
    )

  # Select a real source call as the representative medoid; supporting-library
  # count and finer resolution are used only to resolve exact distance ties.
  df.loop.locus.representative <- df.loop.with.locus %>%
    group_by(approximate_loop_locus_id) %>%
    mutate(
      locus_centroid1_median = median(locus_centroid1),
      locus_centroid2_median = median(locus_centroid2),
      distance_to_locus_median = sqrt(
        (locus_centroid1 - locus_centroid1_median)^2 +
          (locus_centroid2 - locus_centroid2_median)^2
      )
    ) %>%
    arrange(
      distance_to_locus_median,
      dplyr::desc(source_support_count),
      resolution_bp,
      loop_id,
      .by_group = TRUE
    ) %>%
    dplyr::slice(1L) %>%
    ungroup() %>%
    transmute(
      approximate_loop_locus_id,
      representative_loop_id = loop_id,
      representative_resolution = resolution,
      representative_chr1 = chr1,
      representative_start1 = start1,
      representative_end1 = end1,
      representative_chr2 = chr2,
      representative_start2 = start2,
      representative_end2 = end2
    )

  df.loop.locus.summary <- df.loop.with.locus %>%
    group_by(approximate_loop_locus_id) %>%
    summarise(
      n_exact_loop_calls = n(),
      n_resolutions = n_distinct(resolution),
      resolutions = str_c(sort(unique(resolution)), collapse = ";"),
      n_supporting_libraries = n_distinct(
        unlist(str_split(supporting_samples, fixed(";")))
      ),
      supporting_samples = str_c(
        sort(unique(unlist(str_split(supporting_samples, fixed(";"))))),
        collapse = ";"
      ),
      chr1 = dplyr::first(chr1),
      locus_anchor1_start = min(start1),
      locus_anchor1_end = max(end1),
      chr2 = dplyr::first(chr2),
      locus_anchor2_start = min(start2),
      locus_anchor2_end = max(end2),
      centroid1_1based_median = median(locus_centroid1),
      centroid2_1based_median = median(locus_centroid2),
      maximum_source_hiccups_radius_bp = max(hiccups_radius_bp_max),
      .groups = "drop"
    ) %>%
    left_join(
      df.loop.locus.representative,
      by = "approximate_loop_locus_id"
    ) %>%
    arrange(chr1, locus_anchor1_start, chr2, locus_anchor2_start)

  df.loop.locus.method <- tibble(
    analysis_role = "canonical_locus_sensitivity_main_exact_calls_preserved",
    pair_requirement = paste0(
      "same chromosome pair; no shared source library; two-dimensional ",
      "centroid distance <= max of resolution-specific HiCCUPS-like ",
      "merge distances (5K=", merge.distance.bp[["5K"]],
      ",10K=", merge.distance.bp[["10K"]],
      ",25K=", merge.distance.bp[["25K"]], " bp)"
    ),
    aggregation = paste0(
      "deterministic all-pairs-constrained partition; every call pair within ",
      "a locus must qualify, preventing transitive chaining"
    ),
    representative_rule = paste0(
      "source call nearest the locus median two-anchor centroid; ties resolved ",
      "by more supporting libraries, then finer resolution; exact calls and ",
      "provenance remain unchanged"
    ),
    n_exact_loop_calls = nrow(df.loop.ordered),
    n_qualifying_pair_edges = nrow(df.loop.locus.edge),
    n_approximate_loop_loci = nrow(df.loop.locus.summary),
    n_multi_call_loci = sum(df.loop.locus.summary$n_exact_loop_calls > 1L),
    maximum_calls_per_locus = max(df.loop.locus.summary$n_exact_loop_calls)
  )

  list(
    edge = df.loop.locus.edge,
    map = df.loop.locus.map,
    summary = df.loop.locus.summary,
    method = df.loop.locus.method
  )
}

# Count exact loop calls and approximate loop loci per gene after loop-gene deduplication.
summarise_gene_loop_locus_counts <- function(
  df.gene.loop.membership,
  df.loop.locus.map,
  evidence.definition
) {
  check_required_columns(
    df.gene.loop.membership,
    c("loop_id", "ensembl_gene_id", "gene_symbol"),
    "Gene-loop membership for locus sensitivity"
  )
  check_required_columns(
    df.loop.locus.map,
    c("loop_id", "approximate_loop_locus_id"),
    "Approximate loop-locus map"
  )

  df.membership.with.locus <- df.gene.loop.membership %>%
    filter(!is.na(ensembl_gene_id)) %>%
    distinct(loop_id, ensembl_gene_id, .keep_all = TRUE) %>%
    left_join(
      df.loop.locus.map %>%
        dplyr::select(loop_id, approximate_loop_locus_id),
      by = "loop_id"
    )

  if (any(is.na(df.membership.with.locus$approximate_loop_locus_id))) {
    stop(
      "One or more gene-loop records lack an approximate loop-locus ID.",
      call. = FALSE
    )
  }

  df.membership.with.locus %>%
    group_by(ensembl_gene_id) %>%
    summarise(
      gene_symbol = {
        symbol.values <- sort(unique(na.omit(gene_symbol)))
        if (length(symbol.values) == 0L) NA_character_ else symbol.values[[1]]
      },
      n_exact_loop_calls = n_distinct(loop_id),
      n_approximate_loop_loci = n_distinct(approximate_loop_locus_id),
      .groups = "drop"
    ) %>%
    mutate(evidence_definition = evidence.definition, .before = 1) %>%
    arrange(
      dplyr::desc(n_exact_loop_calls),
      dplyr::desc(n_approximate_loop_loci),
      gene_symbol
    )
}

# Summarize one Spearman comparison between two named gene-count columns.
summarise_gene_count_spearman <- function(
  df.gene.count,
  x.column,
  y.column,
  comparison.name
) {
  x <- df.gene.count[[x.column]]
  y <- df.gene.count[[y.column]]
  complete <- is.finite(x) & is.finite(y)

  if (
    sum(complete) < 3L ||
      n_distinct(x[complete]) < 2L ||
      n_distinct(y[complete]) < 2L
  ) {
    return(
      tibble(
        comparison = comparison.name,
        x_metric = x.column,
        y_metric = y.column,
        n_genes = sum(complete),
        spearman_rho = NA_real_,
        p_value = NA_real_
      )
    )
  }

  test.result <- suppressWarnings(
    cor.test(x[complete], y[complete], method = "spearman", exact = FALSE)
  )
  tibble(
    comparison = comparison.name,
    x_metric = x.column,
    y_metric = y.column,
    n_genes = sum(complete),
    spearman_rho = unname(test.result$estimate),
    p_value = test.result$p.value
  )
}

# Compare reconstructed HiCCUPS row counts with documented source counts and stop on any mismatch.
check_hiccups_loop_counts <- function(
  df.sample.loop.1based,
  df.loop.distinct,
  df.loop.universe,
  expected.n
) {
  expected.n <- as_integer_coordinate(expected.n, "Expected HiCCUPS counts")
  if (length(expected.n) != 3L) {
    stop("Expected HiCCUPS counts must contain exactly three values.", call. = FALSE)
  }

  df.loop.source.count.check <- tibble(
    object = c(
      "sample-level source rows",
      "pooled distinct loops",
      "pooled distinct loops shorter than 2 Mb"
    ),
    observed_n = c(
      nrow(df.sample.loop.1based),
      nrow(df.loop.distinct),
      nrow(df.loop.universe)
    ),
    expected_n = expected.n
  ) %>%
    mutate(matches_expected = observed_n == expected_n)

  if (any(!df.loop.source.count.check$matches_expected)) {
    stop(
      "Loop counts reconstructed from the original BEDPE files do not match ",
      "the documented source counts. Inspect df.loop.source.count.check.",
      call. = FALSE
    )
  }

  df.loop.source.count.check
}

################################################################################
# Revision: non-Hi-C coordinate-source readers
################################################################################

# Read transcript records from an Ensembl GTF file. GTF coordinates are already
# one-based and inclusive, so only chromosome-name normalisation is applied.
# Attribute parsing and strand-aware TSS construction are intentionally left to
# the next analysis step.
read_ensembl_gtf_transcripts <- function(file) {
  df.transcript <- readr::read_tsv(
    file,
    comment = "#",
    col_names = c(
      "source_chr", "source", "feature", "start", "end",
      "score", "strand", "frame", "attribute"
    ),
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  ) %>%
    filter(feature == "transcript") %>%
    mutate(
      transcript_start = as_integer_coordinate(
        start,
        "Ensembl GTF transcript start"
      ),
      transcript_end = as_integer_coordinate(
        end,
        "Ensembl GTF transcript end"
      ),
      chr = case_when(
        source_chr == "MT" ~ "chrM",
        str_starts(source_chr, "chr") ~ source_chr,
        TRUE ~ str_c("chr", source_chr)
      )
    )

  validate_one_based_inclusive_intervals(
    df.transcript$transcript_start,
    df.transcript$transcript_end,
    "Ensembl GTF transcript intervals"
  )

  if (any(!df.transcript$strand %in% c("+", "-"))) {
    stop(
      "Ensembl GTF transcript records contain an unsupported strand.",
      call. = FALSE
    )
  }

  df.transcript %>%
    transmute(
      chr,
      transcript_start,
      transcript_end,
      strand,
      source,
      feature,
      attribute,
      source_chr,
      source_coordinate_system = "GTF_1_based_inclusive",
      analysis_coordinate_system = "1_based_inclusive"
    )
}

################################################################################
# Revision: true TSS generation functions
################################################################################

# Extract one quoted value from a semicolon-delimited GTF attribute column.
extract_gtf_attribute <- function(attribute, attribute.name) {
  if (length(attribute.name) != 1L || is.na(attribute.name)) {
    stop("GTF attribute name must be one non-missing value.", call. = FALSE)
  }

  attribute.pattern <- str_c(
    "(?:^|;[[:space:]]*)",
    attribute.name,
    "[[:space:]]+\"([^\"]+)\""
  )
  str_match(attribute, attribute.pattern)[, 2]
}

# Generate one strand-aware one-base TSS record for every Ensembl transcript
# while retaining gene and transcript boundaries. GTF start/end remain in
# ascending genomic order on both strands; strand determines which endpoint is
# the TSS and the direction of transcription.
build_true_tss_annotation <- function(df.transcript) {
  check_required_columns(
    df.transcript,
    c(
      "chr", "transcript_start", "transcript_end", "strand",
      "source", "attribute"
    ),
    "Coordinate-normalized Ensembl transcript annotation"
  )

  df.true.tss <- df.transcript %>%
    transmute(
      chr = as.character(chr),
      transcript_start = as_integer_coordinate(
        transcript_start,
        "Ensembl transcript start"
      ),
      transcript_end = as_integer_coordinate(
        transcript_end,
        "Ensembl transcript end"
      ),
      strand = as.character(strand),
      gene_id = extract_gtf_attribute(attribute, "gene_id"),
      gene_version = extract_gtf_attribute(attribute, "gene_version"),
      gene_id_versioned = if_else(
        is.na(gene_version),
        gene_id,
        str_c(gene_id, ".", gene_version)
      ),
      gene_name = coalesce(
        extract_gtf_attribute(attribute, "gene_name"),
        gene_id
      ),
      gene_biotype = extract_gtf_attribute(attribute, "gene_biotype"),
      transcript_id = extract_gtf_attribute(attribute, "transcript_id"),
      transcript_version = extract_gtf_attribute(
        attribute,
        "transcript_version"
      ),
      transcript_id_versioned = if_else(
        is.na(transcript_version),
        transcript_id,
        str_c(transcript_id, ".", transcript_version)
      ),
      transcript_name = coalesce(
        extract_gtf_attribute(attribute, "transcript_name"),
        transcript_id
      ),
      transcript_biotype = extract_gtf_attribute(
        attribute,
        "transcript_biotype"
      ),
      is_ensembl_canonical = str_detect(
        attribute,
        fixed('tag "Ensembl_canonical"')
      ),
      true_tss_start = case_when(
        strand == "+" ~ transcript_start,
        strand == "-" ~ transcript_end,
        TRUE ~ NA_integer_
      ),
      true_tss_end = true_tss_start,
      transcript_length = transcript_end - transcript_start + 1L,
      true_tss_id = str_c(
        chr,
        true_tss_start,
        strand,
        transcript_id_versioned,
        sep = ":"
      ),
      source,
      source_coordinate_system = "GTF_1_based_inclusive",
      analysis_coordinate_system = "rn7_1_based_inclusive"
    )

  if (any(is.na(df.true.tss$gene_id))) {
    stop("One or more Ensembl transcript records lack gene_id.", call. = FALSE)
  }
  if (any(is.na(df.true.tss$transcript_id))) {
    stop(
      "One or more Ensembl transcript records lack transcript_id.",
      call. = FALSE
    )
  }
  if (any(is.na(df.true.tss$true_tss_start))) {
    stop("True TSS construction encountered an unsupported strand.", call. = FALSE)
  }
  if (
    any(
      df.true.tss$strand == "+" &
        df.true.tss$true_tss_start != df.true.tss$transcript_start
    ) ||
      any(
        df.true.tss$strand == "-" &
          df.true.tss$true_tss_start != df.true.tss$transcript_end
      )
  ) {
    stop(
      "Strand-aware TSS coordinates do not match transcript boundaries.",
      call. = FALSE
    )
  }
  if (
    nrow(df.true.tss) !=
      n_distinct(df.true.tss$transcript_id_versioned)
  ) {
    stop(
      "Ensembl transcript identifiers are not unique after versioning.",
      call. = FALSE
    )
  }

  validate_one_based_inclusive_intervals(
    df.true.tss$transcript_start,
    df.true.tss$transcript_end,
    "Ensembl transcript intervals retained with true TSS"
  )
  validate_one_based_inclusive_intervals(
    df.true.tss$true_tss_start,
    df.true.tss$true_tss_end,
    "Strand-aware true TSS intervals"
  )

  df.true.tss %>%
    arrange(chr, true_tss_start, strand, transcript_id_versioned)
}

# Convert transcript-level true TSS records into one-base stranded GRanges for overlap analysis.
create_true_tss_granges <- function(df.true.tss) {
  check_required_columns(
    df.true.tss,
    c(
      "chr", "true_tss_start", "true_tss_end", "strand",
      "true_tss_id", "gene_id", "gene_name",
      "transcript_id", "transcript_id_versioned",
      "transcript_start", "transcript_end", "is_ensembl_canonical"
    ),
    "Transcript-level true TSS annotation"
  )

  GRanges(
    seqnames = df.true.tss$chr,
    ranges = IRanges(
      start = df.true.tss$true_tss_start,
      end = df.true.tss$true_tss_end
    ),
    strand = df.true.tss$strand,
    true_tss_id = df.true.tss$true_tss_id,
    gene_id = df.true.tss$gene_id,
    gene_name = df.true.tss$gene_name,
    transcript_id = df.true.tss$transcript_id,
    transcript_id_versioned = df.true.tss$transcript_id_versioned,
    transcript_start = df.true.tss$transcript_start,
    transcript_end = df.true.tss$transcript_end,
    is_ensembl_canonical = df.true.tss$is_ensembl_canonical,
    source_coordinate_system = "GTF_1_based_inclusive",
    analysis_coordinate_system = "rn7_1_based_inclusive"
  )
}

# Convert coordinate-normalized EPD TSS positions into one-base stranded
# GRanges while preserving source-row order for lookup after overlap queries.
create_epd_tss_granges <- function(df.promoter) {
  check_required_columns(
    df.promoter,
    c("chr", "epd_tss_start", "epd_tss_end", "strand"),
    "Coordinate-normalized EPD TSS annotation"
  )

  validate_one_based_inclusive_intervals(
    df.promoter$epd_tss_start,
    df.promoter$epd_tss_end,
    "Coordinate-normalized EPD TSS intervals"
  )
  if (any(df.promoter$epd_tss_start != df.promoter$epd_tss_end)) {
    stop("EPD TSS records must be one-base intervals.", call. = FALSE)
  }

  GRanges(
    seqnames = df.promoter$chr,
    ranges = IRanges(
      start = df.promoter$epd_tss_start,
      end = df.promoter$epd_tss_end
    ),
    strand = df.promoter$strand
  )
}

# Expand one-base TSS positions symmetrically into promoter windows. Source
# metadata and row order are preserved so annotation lookup indices stay valid.
expand_tss_to_promoter_windows <- function(gr.tss, flank.bp = 1000L) {
  if (!inherits(gr.tss, "GRanges")) {
    stop("TSS input must be a GRanges object.", call. = FALSE)
  }
  if (any(width(gr.tss) != 1L)) {
    stop("Promoter windows must be constructed from one-base TSSs.", call. = FALSE)
  }

  flank.bp <- as_integer_coordinate(flank.bp, "Promoter-window flank")
  if (length(flank.bp) != 1L || flank.bp < 0L) {
    stop("Promoter-window flank must be one non-negative integer.", call. = FALSE)
  }

  gr.window <- gr.tss
  ranges(gr.window) <- IRanges(
    start = pmax(1L, start(gr.tss) - flank.bp),
    end = end(gr.tss) + flank.bp
  )
  mcols(gr.window)$promoter_window_flank_bp <- flank.bp
  gr.window
}

# Read the original EPD promoter-to-Ensembl mapping and apply the two mapping
# corrections retained from the previous analysis.
read_epd_gene_mapping <- function(mapping.file) {
  tryCatch(
    readr::read_tsv(
      mapping.file,
      col_names = c("epd_promoter_name", "gene_id"),
      col_types = cols(.default = col_character()),
      show_col_types = FALSE
    ),
    error = function(e) {
      stop(
        "Cannot read original EPD mapping file: ", mapping.file, "\n",
        "If it is stored in Google Drive, select Available offline.\n",
        "Original error: ", conditionMessage(e),
        call. = FALSE
      )
    }
  ) %>%
    mutate(
      gene_id = if_else(
        epd_promoter_name == "Cfb_1",
        "ENSRNOG00000051158.3",
        gene_id
      )
    ) %>%
    distinct()
}

# Add gene identifiers to an rn7 EPD TSS table, resolve duplicated promoter
# coordinates using the established _1 rule, and build an 81-bp promoter
# interval centered on each one-base TSS.
finalize_epd_promoters <- function(
  df.tss,
  mapping.file,
  source.coordinate.system,
  promoter.flank = 40L
) {
  check_required_columns(
    df.tss,
    c(
      "chr", "epd_tss_start", "epd_tss_end", "strand",
      "score", "epd_promoter_name"
    ),
    "EPD rn7 TSS table"
  )

  df.mapping <- read_epd_gene_mapping(mapping.file)

  df.promoter <- df.tss %>%
    left_join(df.mapping, by = "epd_promoter_name") %>%
    mutate(
      gene_id = if_else(
        str_detect(epd_promoter_name, "AABR07053687"),
        "ENSRNOG00000015756",
        gene_id
      ),
      gene_id = str_remove(gene_id, "\\.[0-9]+$"),
      gene_name = str_split_n(epd_promoter_name, "_", 1)
    )

  if (any(is.na(df.promoter$gene_id))) {
    stop("One or more EPD promoters lack an Ensembl gene ID.", call. = FALSE)
  }

  df.promoter <- df.promoter %>%
    add_count(
      chr,
      epd_tss_start,
      epd_tss_end,
      gene_id,
      name = "duplicate_count"
    ) %>%
    filter(
      duplicate_count == 1L |
        str_detect(epd_promoter_name, "_1$")
    ) %>%
    transmute(
      chr,
      promoter_start = pmax(1L, epd_tss_start - promoter.flank),
      promoter_end = epd_tss_start + promoter.flank,
      epd_tss_start,
      epd_tss_end,
      strand,
      score,
      epd_promoter_name,
      promoter_annotation_id = str_c(
        chr,
        promoter_start,
        promoter_end,
        strand,
        gene_id,
        gene_name,
        chr,
        epd_tss_start,
        epd_tss_end,
        sep = ":"
      ),
      gene_id,
      gene_name,
      source_coordinate_system = source.coordinate.system,
      analysis_coordinate_system = "rn7_1_based_inclusive"
    )

  validate_one_based_inclusive_intervals(
    df.promoter$promoter_start,
    df.promoter$promoter_end,
    "EPD promoter intervals"
  )
  validate_one_based_inclusive_intervals(
    df.promoter$epd_tss_start,
    df.promoter$epd_tss_end,
    "EPD TSS intervals"
  )

  if (any(!df.promoter$strand %in% c("+", "-"))) {
    stop(
      "EPD promoter records contain an unsupported strand.",
      call. = FALSE
    )
  }

  df.promoter
}

# Construct strand-aware rn6 EPD TSS coordinates directly from the original
# BED fields, validate them against EPD's promoter_coordinate.txt, and lift the
# one-base TSS intervals to rn7. In BED terms the TSS is thickStart + 1 on the
# plus strand and thickEnd on the minus strand.
read_epd_rn6_liftover_promoters <- function(
  bed.file,
  chain.file,
  coordinate.file,
  mapping.file,
  promoter.flank = 40L
) {
  df.bed <- tryCatch(
    readr::read_table(
      bed.file,
      col_names = FALSE,
      col_types = cols(.default = col_character()),
      show_col_types = FALSE
    ),
    error = function(e) {
      stop(
        "Cannot read original EPD rn6 BED file: ", bed.file, "\n",
        "If it is stored in Google Drive, select Available offline.\n",
        "Original error: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (ncol(df.bed) < 8L) {
    stop("The EPD rn6 BED must contain at least eight columns.", call. = FALSE)
  }

  df.rn6.tss <- df.bed %>%
    transmute(
      source_record_id = row_number(),
      source_chr = as.character(.data[["X1"]]),
      promoter_bed_start = as_integer_coordinate(
        .data[["X2"]],
        "EPD rn6 promoter BED start"
      ),
      promoter_bed_end = as_integer_coordinate(
        .data[["X3"]],
        "EPD rn6 promoter BED end"
      ),
      epd_promoter_name = as.character(.data[["X4"]]),
      score = suppressWarnings(as.numeric(.data[["X5"]])),
      strand = as.character(.data[["X6"]]),
      thick_start = as_integer_coordinate(
        .data[["X7"]],
        "EPD rn6 thickStart"
      ),
      thick_end = as_integer_coordinate(
        .data[["X8"]],
        "EPD rn6 thickEnd"
      ),
      rn6_tss = case_when(
        strand == "+" ~ thick_start + 1L,
        strand == "-" ~ thick_end,
        TRUE ~ NA_integer_
      )
    )

  validate_bed_half_open_intervals(
    df.rn6.tss$promoter_bed_start,
    df.rn6.tss$promoter_bed_end,
    "EPD rn6 promoter BED intervals"
  )

  if (any(is.na(df.rn6.tss$rn6_tss))) {
    stop("EPD rn6 records contain an unsupported strand.", call. = FALSE)
  }

  df.reported.coordinate <- readr::read_tsv(
    coordinate.file,
    col_names = c(
      "epd_promoter_name", "accession", "reported_tss", "reported_strand",
      "species", "promoter_type"
    ),
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  ) %>%
    transmute(
      epd_promoter_name,
      reported_tss = as_integer_coordinate(
        reported_tss,
        "EPD reported TSS"
      ),
      reported_strand
    )

  df.coordinate.check <- df.rn6.tss %>%
    left_join(df.reported.coordinate, by = "epd_promoter_name")

  if (
    any(is.na(df.coordinate.check$reported_tss)) ||
      any(df.coordinate.check$rn6_tss != df.coordinate.check$reported_tss) ||
      any(df.coordinate.check$strand != df.coordinate.check$reported_strand)
  ) {
    stop(
      "Strand-aware TSS coordinates derived from the EPD rn6 BED do not ",
      "match promoter_coordinate.txt.",
      call. = FALSE
    )
  }

  gr.rn6.tss <- GRanges(
    seqnames = df.rn6.tss$source_chr,
    ranges = IRanges(start = df.rn6.tss$rn6_tss, end = df.rn6.tss$rn6_tss),
    strand = df.rn6.tss$strand,
    source_record_id = df.rn6.tss$source_record_id,
    epd_promoter_name = df.rn6.tss$epd_promoter_name,
    score = df.rn6.tss$score,
    rn6_tss = df.rn6.tss$rn6_tss
  )

  chain.rn6.to.rn7 <- rtracklayer::import.chain(chain.file)
  gr.rn7.list <- rtracklayer::liftOver(gr.rn6.tss, chain.rn6.to.rn7)
  n.mappings <- S4Vectors::elementNROWS(gr.rn7.list)

  if (any(n.mappings > 1L)) {
    stop("One or more EPD TSS records map to multiple rn7 loci.", call. = FALSE)
  }

  gr.rn7.tss <- unlist(gr.rn7.list, use.names = FALSE)
  df.rn7.tss <- tibble(
    chr = as.character(seqnames(gr.rn7.tss)),
    epd_tss_start = start(gr.rn7.tss),
    epd_tss_end = end(gr.rn7.tss),
    strand = as.character(strand(gr.rn7.tss)),
    score = as.numeric(mcols(gr.rn7.tss)$score),
    epd_promoter_name = as.character(
      mcols(gr.rn7.tss)$epd_promoter_name
    )
  )

  df.promoter <- finalize_epd_promoters(
    df.tss = df.rn7.tss,
    mapping.file = mapping.file,
    source.coordinate.system = paste0(
      "EPD_rn6_BED_strand_aware_TSS_lifted_to_rn7; ",
      "rn7_1_based_inclusive"
    ),
    promoter.flank = promoter.flank
  )

  attr(df.promoter, "n_rn6_input") <- nrow(df.rn6.tss)
  attr(df.promoter, "n_rn7_lifted") <- length(gr.rn7.tss)
  attr(df.promoter, "n_failed_liftover") <- sum(n.mappings == 0L)
  attr(df.promoter, "failed_liftover_promoters") <-
    df.rn6.tss$epd_promoter_name[n.mappings == 0L]
  df.promoter
}

# Reconstruct the legacy promoter annotation from the previously exported rn7
# BED. This reader is retained only to quantify the coordinate correction; the
# revised primary promoter annotation comes from the rn6 liftOver reader above.
read_epd_rn7_bed_promoters <- function(
  bed.file,
  mapping.file,
  promoter.flank = 40L
) {
  df.bed <- tryCatch(
    readr::read_tsv(
      bed.file,
      comment = "#",
      col_names = FALSE,
      col_types = cols(.default = col_character()),
      show_col_types = FALSE
    ),
    error = function(e) {
      stop(
        "Cannot read derived EPD rn7 BED file: ", bed.file, "\n",
        "If it is stored in Google Drive, select Available offline.\n",
        "Original error: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (ncol(df.bed) < 6L) {
    stop("The EPD rn7 BED must contain at least six columns.", call. = FALSE)
  }

  df.rn7.tss <- df.bed %>%
    transmute(
      chr = as.character(.data[["X1"]]),
      bed_start = as_integer_coordinate(.data[["X2"]], "EPD BED start"),
      bed_end = as_integer_coordinate(.data[["X3"]], "EPD BED end"),
      epd_promoter_name = as.character(.data[["X4"]]),
      score = suppressWarnings(as.numeric(.data[["X5"]])),
      strand = as.character(.data[["X6"]])
    )

  validate_bed_half_open_intervals(
    df.rn7.tss$bed_start,
    df.rn7.tss$bed_end,
    "Derived EPD rn7 BED intervals"
  )

  df.rn7.tss <- df.rn7.tss %>%
    transmute(
      chr,
      epd_tss_start = bed_start_to_one_based(
        bed_start,
        "Derived EPD rn7 TSS start"
      ),
      epd_tss_end = as_integer_coordinate(
        bed_end,
        "Derived EPD rn7 TSS end"
      ),
      strand,
      score,
      epd_promoter_name
    )

  finalize_epd_promoters(
    df.tss = df.rn7.tss,
    mapping.file = mapping.file,
    source.coordinate.system = "derived_rn7_BED_0_based_half_open",
    promoter.flank = promoter.flank
  )
}

# Read the former cached EPD promoter object. This compatibility reader remains
# only until all downstream legacy-comparison code has moved to the raw BED
# reader above; it is not intended for the revised primary analysis.
read_epd_promoter_annotation <- function(file) {
  df.promoter <- readRDS(file)
  check_required_columns(
    df.promoter,
    c(
      "seqnames", "start", "end", "tss_start", "tss_end", "strand",
      "promoter_id", "gene_id", "gene_name", "promoter.id",
      "ensembl_exon_id", "refseq_exon_id"
    ),
    "EPD promoter annotation"
  )

  df.promoter <- df.promoter %>%
    transmute(
      chr = as.character(seqnames),
      promoter_start = as_integer_coordinate(
        start,
        "EPD promoter start"
      ),
      promoter_end = as_integer_coordinate(
        end,
        "EPD promoter end"
      ),
      epd_tss_start = as_integer_coordinate(
        tss_start,
        "EPD TSS start"
      ),
      epd_tss_end = as_integer_coordinate(
        tss_end,
        "EPD TSS end"
      ),
      strand = as.character(strand),
      epd_promoter_name = as.character(promoter_id),
      promoter_annotation_id = as.character(promoter.id),
      gene_id = as.character(gene_id),
      gene_name = as.character(gene_name),
      ensembl_exon_id = as.character(ensembl_exon_id),
      refseq_exon_id = as.character(refseq_exon_id),
      source_coordinate_system = "GRanges_1_based_inclusive",
      analysis_coordinate_system = "1_based_inclusive"
    )

  validate_one_based_inclusive_intervals(
    df.promoter$promoter_start,
    df.promoter$promoter_end,
    "EPD promoter intervals"
  )
  validate_one_based_inclusive_intervals(
    df.promoter$epd_tss_start,
    df.promoter$epd_tss_end,
    "EPD TSS intervals"
  )

  if (any(!df.promoter$strand %in% c("+", "-"))) {
    stop(
      "EPD promoter records contain an unsupported strand.",
      call. = FALSE
    )
  }

  df.promoter
}

# Convert either side of a normalised loop table into GRanges. The function
# validates one-based inclusive anchor coordinates and carries loop identity,
# resolution, and coordinate provenance as metadata.
create_loop_anchor_granges <- function(df, anchor.side) {
  if (anchor.side == "anchor1") {
    chr <- df$chr1
    start <- df$start1
    end <- df$end1
  } else if (anchor.side == "anchor2") {
    chr <- df$chr2
    start <- df$start2
    end <- df$end2
  } else {
    stop("Unknown anchor side: ", anchor.side, call. = FALSE)
  }

  validate_one_based_inclusive_intervals(
    start,
    end,
    paste0("Normalised loop ", anchor.side)
  )

  GRanges(
    seqnames = chr,
    ranges = IRanges(start = start, end = end),
    loop_id = df$loop_id,
    resolution = df$resolution,
    anchor_side = anchor.side,
    source_coordinate_system = "HiCCUPS_BEDPE_0_based_half_open",
    analysis_coordinate_system = "1_based_inclusive"
  )
}

################################################################################
# Revision: direct promoter/TSS-anchor overlap functions
################################################################################

# Convert normalized EPD promoter intervals into stranded GRanges while retaining promoter and gene identifiers.
create_epd_promoter_granges <- function(df.promoter) {
  check_required_columns(
    df.promoter,
    c(
      "chr", "promoter_start", "promoter_end",
      "epd_tss_start", "epd_tss_end", "strand",
      "epd_promoter_name", "promoter_annotation_id",
      "gene_id", "gene_name"
    ),
    "Normalized EPD promoter annotation"
  )

  validate_one_based_inclusive_intervals(
    df.promoter$promoter_start,
    df.promoter$promoter_end,
    "Normalized EPD promoter intervals"
  )

  GRanges(
    seqnames = df.promoter$chr,
    ranges = IRanges(
      start = df.promoter$promoter_start,
      end = df.promoter$promoter_end
    ),
    strand = df.promoter$strand,
    promoter_annotation_id = df.promoter$promoter_annotation_id,
    epd_promoter_name = df.promoter$epd_promoter_name,
    epd_tss_start = df.promoter$epd_tss_start,
    epd_tss_end = df.promoter$epd_tss_end,
    gene_id = df.promoter$gene_id,
    gene_name = df.promoter$gene_name,
    source_coordinate_system = df.promoter$source_coordinate_system,
    analysis_coordinate_system = df.promoter$analysis_coordinate_system
  )
}

# Return every interval overlap between one loop-anchor side and an annotation set without nearest-feature selection or deduplication.
find_direct_anchor_annotation_overlaps <- function(
  gr.anchor,
  gr.annotation
) {
  required.anchor.metadata <- c("loop_id", "resolution", "anchor_side")
  missing.anchor.metadata <- setdiff(
    required.anchor.metadata,
    colnames(mcols(gr.anchor))
  )
  if (length(missing.anchor.metadata) > 0L) {
    stop(
      "Loop-anchor GRanges lacks required metadata: ",
      paste(missing.anchor.metadata, collapse = ", "),
      call. = FALSE
    )
  }

  hits <- findOverlaps(
    gr.anchor,
    gr.annotation,
    ignore.strand = TRUE
  )

  if (length(hits) == 0L) {
    return(
      tibble(
        loop_id = character(),
        resolution = character(),
        anchor_side = character(),
        opposite_anchor_side = character(),
        anchor_chr = character(),
        anchor_start = integer(),
        anchor_end = integer(),
        anchor_width_bp = integer(),
        annotation_index = integer(),
        annotation_chr = character(),
        annotation_start = integer(),
        annotation_end = integer(),
        annotation_width_bp = integer(),
        annotation_strand = character(),
        direct_overlap_bp = integer()
      )
    )
  }

  query.index <- queryHits(hits)
  subject.index <- subjectHits(hits)
  gr.anchor.hit <- gr.anchor[query.index]
  gr.annotation.hit <- gr.annotation[subject.index]
  overlap.start <- pmax(start(gr.anchor.hit), start(gr.annotation.hit))
  overlap.end <- pmin(end(gr.anchor.hit), end(gr.annotation.hit))

  df.overlap <- tibble(
    loop_id = as.character(mcols(gr.anchor.hit)$loop_id),
    resolution = as.character(mcols(gr.anchor.hit)$resolution),
    anchor_side = as.character(mcols(gr.anchor.hit)$anchor_side),
    opposite_anchor_side = if_else(
      anchor_side == "anchor1",
      "anchor2",
      "anchor1"
    ),
    anchor_chr = as.character(seqnames(gr.anchor.hit)),
    anchor_start = start(gr.anchor.hit),
    anchor_end = end(gr.anchor.hit),
    anchor_width_bp = width(gr.anchor.hit),
    annotation_index = as.integer(subject.index),
    annotation_chr = as.character(seqnames(gr.annotation.hit)),
    annotation_start = start(gr.annotation.hit),
    annotation_end = end(gr.annotation.hit),
    annotation_width_bp = width(gr.annotation.hit),
    annotation_strand = as.character(strand(gr.annotation.hit)),
    direct_overlap_bp = overlap.end - overlap.start + 1L
  )

  if (any(df.overlap$direct_overlap_bp < 1L)) {
    stop(
      "A direct anchor-annotation hit has a non-positive overlap width.",
      call. = FALSE
    )
  }

  df.overlap
}

# Build one strict or +/-1-kb direct promoter/TSS evidence tier from normalized
# Ensembl and EPD annotations, then return the lossless records and summaries.
build_direct_promoter_tss_tier <- function(
  gr.loop.anchor.by.side,
  gr.true.tss.annotation,
  gr.epd.annotation,
  df.true.tss,
  df.epd.promoter,
  df.loop.universe,
  evidence.definition = c("strict", "primary_1kb"),
  promoter.window.flank.bp = NA_integer_
) {
  evidence.definition <- match.arg(evidence.definition)
  is.primary <- evidence.definition == "primary_1kb"

  true.source.lookup <- df.true.tss %>%
    mutate(annotation_index = row_number()) %>%
    select(
      annotation_index,
      true_tss_id,
      true_tss_start,
      true_tss_end,
      gene_id,
      gene_id_versioned,
      gene_name,
      gene_biotype,
      transcript_id,
      transcript_id_versioned,
      transcript_name,
      transcript_biotype,
      is_ensembl_canonical,
      transcript_start,
      transcript_end,
      source,
      source_coordinate_system,
      analysis_coordinate_system
    )

  epd.source.lookup <- df.epd.promoter %>%
    mutate(annotation_index = row_number()) %>%
    select(
      annotation_index,
      promoter_annotation_id,
      epd_promoter_name,
      epd_tss_start,
      epd_tss_end,
      gene_id,
      gene_name,
      source_coordinate_system,
      analysis_coordinate_system
    )

  true.lookup <- true.source.lookup %>%
    transmute(
      annotation_index,
      annotation_class = "true_TSS",
      base_annotation_id = true_tss_id,
      gene_id,
      gene_id_versioned,
      gene_name,
      gene_biotype,
      transcript_id,
      transcript_id_versioned,
      transcript_name,
      transcript_biotype,
      is_ensembl_canonical,
      transcript_start,
      transcript_end,
      promoter_annotation_id = NA_character_,
      epd_promoter_name = NA_character_,
      epd_tss_start = NA_integer_,
      epd_tss_end = NA_integer_,
      tss_start = true_tss_start,
      tss_end = true_tss_end,
      lookup_source = source,
      source_coordinate_system,
      analysis_coordinate_system
    )

  epd.lookup <- epd.source.lookup %>%
    transmute(
      annotation_index,
      annotation_class = "EPD_promoter",
      base_annotation_id = promoter_annotation_id,
      gene_id,
      gene_id_versioned = gene_id,
      gene_name,
      gene_biotype = NA_character_,
      transcript_id = NA_character_,
      transcript_id_versioned = NA_character_,
      transcript_name = NA_character_,
      transcript_biotype = NA_character_,
      is_ensembl_canonical = NA,
      transcript_start = NA_integer_,
      transcript_end = NA_integer_,
      promoter_annotation_id,
      epd_promoter_name,
      epd_tss_start,
      epd_tss_end,
      tss_start = epd_tss_start,
      tss_end = epd_tss_end,
      lookup_source = "EPDnew",
      source_coordinate_system,
      analysis_coordinate_system
    )

  build_source_evidence <- function(gr.annotation, df.lookup) {
    map_loop_anchors_dfr(
      gr.loop.anchor.by.side,
      find_direct_anchor_annotation_overlaps,
      gr.annotation = gr.annotation
    ) %>%
      left_join(df.lookup, by = "annotation_index")
  }

  df.evidence <- bind_rows(
    build_source_evidence(gr.true.tss.annotation, true.lookup),
    build_source_evidence(gr.epd.annotation, epd.lookup)
  ) %>%
    mutate(
      annotation_source = if (is.primary) {
        if_else(
          annotation_class == "true_TSS",
          str_c(lookup_source, "_GTF_TSS_plus_minus_1kb"),
          "EPDnew_TSS_plus_minus_1kb"
        )
      } else {
        if_else(
          annotation_class == "true_TSS",
          str_c(lookup_source, "_GTF_transcript"),
          "EPDnew_promoter"
        )
      },
      annotation_id = if (is.primary) {
        str_c(base_annotation_id, ":TSS_pm1kb")
      } else {
        base_annotation_id
      }
    ) %>%
    transmute(
      direct_evidence_id = str_c(
        loop_id, anchor_side, annotation_source, annotation_id, sep = "|"
      ),
      loop_id,
      resolution,
      anchor_side,
      opposite_anchor_side,
      anchor_chr,
      anchor_start,
      anchor_end,
      anchor_width_bp,
      annotation_class,
      annotation_source,
      annotation_id,
      annotation_chr,
      annotation_start,
      annotation_end,
      annotation_width_bp,
      annotation_strand,
      direct_overlap_bp,
      tss_start,
      tss_end,
      promoter_window_flank_bp = as.integer(promoter.window.flank.bp),
      gene_id,
      gene_id_versioned,
      gene_name,
      gene_biotype,
      transcript_id,
      transcript_id_versioned,
      transcript_name,
      transcript_biotype,
      is_ensembl_canonical,
      transcript_start,
      transcript_end,
      promoter_annotation_id,
      epd_promoter_name,
      epd_tss_start,
      epd_tss_end,
      annotation_source_coordinate_system = source_coordinate_system,
      analysis_coordinate_system
    ) %>%
    arrange(
      loop_id, anchor_side, gene_id, annotation_class,
      annotation_start, annotation_id
    )

  if (!is.primary) {
    df.evidence <- df.evidence %>%
      select(-tss_start, -tss_end, -promoter_window_flank_bp)
  }

  list(
    true_tss_lookup = true.source.lookup,
    epd_promoter_lookup = epd.source.lookup,
    true_tss_overlap = filter(df.evidence, annotation_class == "true_TSS"),
    epd_promoter_overlap = filter(
      df.evidence,
      annotation_class == "EPD_promoter"
    ),
    combined_overlap = df.evidence,
    summary = summarise_direct_promoter_tss_evidence(
      df.evidence,
      df.loop.universe
    )
  )
}

# Collapse a lossless direct promoter/TSS evidence table into consistent
# loop-anchor-gene, anchor-level, and loop-level summaries. Multiple genes are
# preserved; transcript and annotation duplicates are collapsed only within the
# same loop-anchor-gene assignment.
summarise_direct_promoter_tss_evidence <- function(
  df.evidence,
  df.loop.universe
) {
  check_required_columns(
    df.evidence,
    c(
      "direct_evidence_id", "loop_id", "resolution", "anchor_side",
      "opposite_anchor_side", "annotation_class", "gene_id", "gene_name",
      "transcript_id_versioned", "promoter_annotation_id"
    ),
    "Direct promoter/TSS evidence table"
  )
  check_required_columns(
    df.loop.universe,
    c(
      "loop_id", "resolution", "chr1", "start1", "end1",
      "chr2", "start2", "end2"
    ),
    "Pooled loop universe"
  )
  assert_analysis_unique_key(
    df.evidence,
    "direct_evidence_id",
    "Direct promoter/TSS evidence IDs are not unique."
  )

  df.gene.assignment <- df.evidence %>%
    group_by(
      loop_id,
      resolution,
      anchor_side,
      opposite_anchor_side,
      gene_id
    ) %>%
    summarise(
      gene_name = {
        gene.names <- sort(unique(na.omit(gene_name)))
        if (length(gene.names) == 0L) NA_character_ else gene.names[1]
      },
      n_direct_evidence_records = n(),
      has_direct_true_tss = any(annotation_class == "true_TSS"),
      has_direct_epd_promoter = any(annotation_class == "EPD_promoter"),
      n_true_tss_transcripts = n_distinct(
        transcript_id_versioned,
        na.rm = TRUE
      ),
      n_epd_promoters = n_distinct(
        promoter_annotation_id,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    arrange(loop_id, anchor_side, gene_id)

  df.anchor.index <- build_loop_anchor_index(df.loop.universe)

  df.anchor.count <- df.evidence %>%
    group_by(loop_id, resolution, anchor_side, opposite_anchor_side) %>%
    summarise(
      n_direct_evidence_records = n(),
      n_direct_genes = n_distinct(gene_id),
      n_true_tss_records = sum(annotation_class == "true_TSS"),
      n_true_tss_genes = n_distinct(
        gene_id[annotation_class == "true_TSS"]
      ),
      n_epd_promoter_records = sum(annotation_class == "EPD_promoter"),
      n_epd_promoter_genes = n_distinct(
        gene_id[annotation_class == "EPD_promoter"]
      ),
      .groups = "drop"
    )

  df.anchor.summary <- df.anchor.index %>%
    left_join(
      df.anchor.count,
      by = c(
        "loop_id", "resolution", "anchor_side", "opposite_anchor_side"
      )
    ) %>%
    mutate(
      across(
        c(
          n_direct_evidence_records,
          n_direct_genes,
          n_true_tss_records,
          n_true_tss_genes,
          n_epd_promoter_records,
          n_epd_promoter_genes
        ),
        ~ coalesce(.x, 0L)
      ),
      has_direct_true_tss = n_true_tss_records > 0L,
      has_direct_epd_promoter = n_epd_promoter_records > 0L,
      has_any_direct_promoter_tss = n_direct_evidence_records > 0L,
      has_both_direct_annotation_classes = (
        has_direct_true_tss & has_direct_epd_promoter
      )
    ) %>%
    arrange(loop_id, anchor_side)

  assert_analysis_row_count(
    df.anchor.summary,
    2L * nrow(df.loop.universe),
    "Direct-overlap summary did not retain both anchors of every loop."
  )

  df.anchor.wide <- df.anchor.summary %>%
    dplyr::select(
      loop_id,
      anchor_side,
      n_direct_evidence_records,
      n_direct_genes,
      n_true_tss_records,
      n_epd_promoter_records,
      has_direct_true_tss,
      has_direct_epd_promoter,
      has_any_direct_promoter_tss,
      has_both_direct_annotation_classes
    ) %>%
    pivot_wider(
      names_from = anchor_side,
      values_from = -loop_id,
      names_glue = "{.value}_{anchor_side}"
    )

  df.loop.gene.count <- df.gene.assignment %>%
    group_by(loop_id) %>%
    summarise(
      n_direct_genes_across_anchors = n_distinct(gene_id),
      .groups = "drop"
    )

  df.loop.summary <- df.loop.universe %>%
    left_join(df.anchor.wide, by = "loop_id") %>%
    left_join(df.loop.gene.count, by = "loop_id") %>%
    mutate(
      n_direct_genes_across_anchors = coalesce(
        n_direct_genes_across_anchors,
        0L
      ),
      n_direct_anchor_sides =
        as.integer(has_any_direct_promoter_tss_anchor1) +
          as.integer(has_any_direct_promoter_tss_anchor2),
      has_any_direct_promoter_tss = n_direct_anchor_sides > 0L,
      has_direct_promoter_tss_both_anchors = n_direct_anchor_sides == 2L,
      has_direct_true_tss_any_anchor = (
        has_direct_true_tss_anchor1 | has_direct_true_tss_anchor2
      ),
      has_direct_epd_promoter_any_anchor = (
        has_direct_epd_promoter_anchor1 |
          has_direct_epd_promoter_anchor2
      ),
      direct_anchor_assignment_class = case_when(
        n_direct_anchor_sides == 2L ~ "direct_both_anchors",
        has_any_direct_promoter_tss_anchor1 ~ "direct_anchor1_only",
        has_any_direct_promoter_tss_anchor2 ~ "direct_anchor2_only",
        TRUE ~ "no_direct_promoter_or_TSS"
      )
    )

  assert_analysis_condition(
    nrow(df.loop.summary) == nrow(df.loop.universe) &&
      !any(is.na(df.loop.summary$n_direct_anchor_sides)),
    "Direct-overlap loop summary did not preserve the pooled universe."
  )

  df.summary <- tibble(
    metric = c(
      "pooled_loops",
      "direct_true_TSS_overlap_records",
      "direct_EPD_promoter_overlap_records",
      "unique_loop_anchor_gene_assignments",
      "loops_with_any_direct_promoter_or_TSS",
      "loops_with_direct_promoter_or_TSS_at_one_anchor",
      "loops_with_direct_promoter_or_TSS_at_both_anchors",
      "loops_without_direct_promoter_or_TSS"
    ),
    n = c(
      nrow(df.loop.universe),
      sum(df.evidence$annotation_class == "true_TSS"),
      sum(df.evidence$annotation_class == "EPD_promoter"),
      nrow(df.gene.assignment),
      sum(df.loop.summary$has_any_direct_promoter_tss),
      sum(df.loop.summary$n_direct_anchor_sides == 1L),
      sum(df.loop.summary$has_direct_promoter_tss_both_anchors),
      sum(!df.loop.summary$has_any_direct_promoter_tss)
    )
  )

  list(
    gene_assignment = df.gene.assignment,
    anchor_index = df.anchor.index,
    anchor_count = df.anchor.count,
    anchor_summary = df.anchor.summary,
    anchor_wide = df.anchor.wide,
    loop_gene_count = df.loop.gene.count,
    loop_summary = df.loop.summary,
    summary = df.summary
  )
}

# Return every non-overlapping annotation within a fixed edge-to-edge distance of one loop-anchor side.
find_proximal_anchor_annotation_pairs <- function(
  gr.anchor,
  gr.annotation,
  max.distance = 200000L
) {
  required.anchor.metadata <- c("loop_id", "resolution", "anchor_side")
  missing.anchor.metadata <- setdiff(
    required.anchor.metadata,
    colnames(mcols(gr.anchor))
  )
  if (length(missing.anchor.metadata) > 0L) {
    stop(
      "Loop-anchor GRanges lacks required metadata: ",
      paste(missing.anchor.metadata, collapse = ", "),
      call. = FALSE
    )
  }

  max.distance <- as_integer_coordinate(
    max.distance,
    "Maximum proximal anchor-annotation distance"
  )
  if (length(max.distance) != 1L || max.distance < 1L) {
    stop(
      "Maximum proximal distance must be one positive integer.",
      call. = FALSE
    )
  }

  gr.anchor.search <- gr.anchor
  ranges(gr.anchor.search) <- IRanges(
    start = pmax(1L, start(gr.anchor) - max.distance),
    end = end(gr.anchor) + max.distance
  )

  hits <- findOverlaps(
    gr.anchor.search,
    gr.annotation,
    ignore.strand = TRUE
  )

  if (length(hits) == 0L) {
    return(
      tibble(
        loop_id = character(),
        resolution = character(),
        anchor_side = character(),
        opposite_anchor_side = character(),
        anchor_chr = character(),
        anchor_start = integer(),
        anchor_end = integer(),
        anchor_width_bp = integer(),
        annotation_index = integer(),
        annotation_chr = character(),
        annotation_start = integer(),
        annotation_end = integer(),
        annotation_width_bp = integer(),
        annotation_strand = character(),
        annotation_relative_to_anchor = character(),
        anchor_annotation_edge_distance_bp = integer(),
        proximal_max_distance_bp = integer()
      )
    )
  }

  query.index <- queryHits(hits)
  subject.index <- subjectHits(hits)
  gr.anchor.hit <- gr.anchor[query.index]
  gr.annotation.hit <- gr.annotation[subject.index]

  annotation.relative.to.anchor <- case_when(
    end(gr.annotation.hit) < start(gr.anchor.hit) ~ "left_of_anchor",
    start(gr.annotation.hit) > end(gr.anchor.hit) ~ "right_of_anchor",
    TRUE ~ "direct_overlap"
  )
  edge.distance <- case_when(
    annotation.relative.to.anchor == "left_of_anchor" ~
      start(gr.anchor.hit) - end(gr.annotation.hit),
    annotation.relative.to.anchor == "right_of_anchor" ~
      start(gr.annotation.hit) - end(gr.anchor.hit),
    TRUE ~ 0L
  )

  df.proximal <- tibble(
    loop_id = as.character(mcols(gr.anchor.hit)$loop_id),
    resolution = as.character(mcols(gr.anchor.hit)$resolution),
    anchor_side = as.character(mcols(gr.anchor.hit)$anchor_side),
    opposite_anchor_side = if_else(
      anchor_side == "anchor1",
      "anchor2",
      "anchor1"
    ),
    anchor_chr = as.character(seqnames(gr.anchor.hit)),
    anchor_start = start(gr.anchor.hit),
    anchor_end = end(gr.anchor.hit),
    anchor_width_bp = width(gr.anchor.hit),
    annotation_index = as.integer(subject.index),
    annotation_chr = as.character(seqnames(gr.annotation.hit)),
    annotation_start = start(gr.annotation.hit),
    annotation_end = end(gr.annotation.hit),
    annotation_width_bp = width(gr.annotation.hit),
    annotation_strand = as.character(strand(gr.annotation.hit)),
    annotation_relative_to_anchor = annotation.relative.to.anchor,
    anchor_annotation_edge_distance_bp = as.integer(edge.distance),
    proximal_max_distance_bp = max.distance
  ) %>%
    filter(
      anchor_annotation_edge_distance_bp >= 1L,
      anchor_annotation_edge_distance_bp <= proximal_max_distance_bp
    )

  if (
    any(df.proximal$annotation_relative_to_anchor == "direct_overlap") ||
      any(df.proximal$anchor_annotation_edge_distance_bp < 1L) ||
      any(
        df.proximal$anchor_annotation_edge_distance_bp >
          df.proximal$proximal_max_distance_bp
      )
  ) {
    stop(
      "The proximal table contains a direct or out-of-range annotation.",
      call. = FALSE
    )
  }

  df.proximal
}

# Classify positive anchor-edge distances into descriptive 10-kb, 50-kb, and 200-kb proximity bands.
classify_proximal_distance_tier <- function(distance) {
  distance <- as.integer(distance)
  case_when(
    distance >= 1L & distance <= 10000L ~ "01_1bp_to_10kb",
    distance > 10000L & distance <= 50000L ~ "02_gt10kb_to_50kb",
    distance > 50000L & distance <= 200000L ~ "03_gt50kb_to_200kb",
    TRUE ~ NA_character_
  )
}

################################################################################
# Revision: true-TSS-excluded ATAC functions
################################################################################

# Summarise overlap between every loop anchor and a reduced genomic feature set,
# retaining interval counts, covered bases, and both any-bp and >=minimum flags.
summarise_anchor_interval_overlap <- function(
  gr.anchor,
  gr.feature,
  minimum.overlap.bp = 50L
) {
  required.anchor.metadata <- c("loop_id", "resolution", "anchor_side")
  missing.anchor.metadata <- setdiff(
    required.anchor.metadata,
    colnames(mcols(gr.anchor))
  )
  if (length(missing.anchor.metadata) > 0L) {
    stop(
      "Loop-anchor GRanges lacks required metadata: ",
      paste(missing.anchor.metadata, collapse = ", "),
      call. = FALSE
    )
  }

  minimum.overlap.bp <- as.integer(minimum.overlap.bp)
  if (
    length(minimum.overlap.bp) != 1L ||
      is.na(minimum.overlap.bp) ||
      minimum.overlap.bp < 1L
  ) {
    stop("Minimum overlap must be one positive integer.", call. = FALSE)
  }

  gr.feature <- GenomicRanges::reduce(
    gr.feature,
    ignore.strand = TRUE
  )
  n.anchor <- length(gr.anchor)
  n.feature.intervals <- integer(n.anchor)
  n.feature.intervals.ge.minimum <- integer(n.anchor)
  feature.overlap.bp <- integer(n.anchor)

  if (n.anchor > 0L && length(gr.feature) > 0L) {
    hit.any <- findOverlaps(
      gr.anchor,
      gr.feature,
      ignore.strand = TRUE
    )

    if (length(hit.any) > 0L) {
      query.index <- queryHits(hit.any)
      n.feature.intervals <- tabulate(query.index, nbins = n.anchor)

      gr.intersection <- pintersect(
        gr.anchor[query.index],
        gr.feature[subjectHits(hit.any)],
        ignore.strand = TRUE
      )
      overlap.sum <- rowsum(
        width(gr.intersection),
        group = query.index,
        reorder = FALSE
      )
      feature.overlap.bp[as.integer(rownames(overlap.sum))] <-
        as.integer(overlap.sum[, 1])
    }

    hit.minimum <- findOverlaps(
      gr.anchor,
      gr.feature,
      minoverlap = minimum.overlap.bp,
      ignore.strand = TRUE
    )
    if (length(hit.minimum) > 0L) {
      n.feature.intervals.ge.minimum <- tabulate(
        queryHits(hit.minimum),
        nbins = n.anchor
      )
    }
  }

  tibble(
    anchor_index = seq_len(n.anchor),
    loop_id = as.character(mcols(gr.anchor)$loop_id),
    resolution = as.character(mcols(gr.anchor)$resolution),
    anchor_side = as.character(mcols(gr.anchor)$anchor_side),
    anchor_chr = as.character(seqnames(gr.anchor)),
    anchor_start = start(gr.anchor),
    anchor_end = end(gr.anchor),
    anchor_width_bp = width(gr.anchor),
    n_overlapping_feature_intervals = n.feature.intervals,
    n_overlapping_feature_intervals_ge_minimum =
      n.feature.intervals.ge.minimum,
    feature_overlap_bp = feature.overlap.bp,
    feature_overlap_fraction_of_anchor = safe_ratio(
      feature.overlap.bp,
      width(gr.anchor)
    ),
    has_feature_overlap_any = n.feature.intervals > 0L,
    has_feature_overlap_ge_minimum =
      n.feature.intervals.ge.minimum > 0L,
    has_feature_total_overlap_ge_minimum =
      feature.overlap.bp >= minimum.overlap.bp,
    minimum_overlap_bp = minimum.overlap.bp
  )
}

# Subtract reduced exclusion intervals from each anchor independently so that
# overlapping anchors never merge and every surviving non-excluded fragment
# retains its original loop and anchor identity.
subtract_exclusion_from_anchor_ranges <- function(
  gr.anchor,
  gr.exclusion
) {
  required.anchor.metadata <- c("loop_id", "resolution", "anchor_side")
  missing.anchor.metadata <- setdiff(
    required.anchor.metadata,
    colnames(mcols(gr.anchor))
  )
  if (length(missing.anchor.metadata) > 0L) {
    stop(
      "Loop-anchor GRanges lacks required metadata: ",
      paste(missing.anchor.metadata, collapse = ", "),
      call. = FALSE
    )
  }

  empty.result <- tibble(
    source_anchor_index = integer(),
    loop_id = character(),
    resolution = character(),
    anchor_side = character(),
    fragment_index = integer(),
    fragment_chr = character(),
    fragment_start = integer(),
    fragment_end = integer(),
    fragment_width_bp = integer()
  )
  if (length(gr.anchor) == 0L) {
    return(empty.result)
  }

  gr.exclusion <- GenomicRanges::reduce(
    gr.exclusion,
    ignore.strand = TRUE
  )
  exclusion.hit <- findOverlaps(
    gr.anchor,
    gr.exclusion,
    ignore.strand = TRUE
  )
  exclusion.index.by.anchor <- split(
    subjectHits(exclusion.hit),
    queryHits(exclusion.hit)
  )

  fragment.rows <- vector("list", length(gr.anchor))
  for (anchor.index in seq_along(gr.anchor)) {
    exclusion.index <- exclusion.index.by.anchor[[as.character(anchor.index)]]
    residual.range <- if (is.null(exclusion.index)) {
      ranges(gr.anchor[anchor.index])
    } else {
      IRanges::setdiff(
        ranges(gr.anchor[anchor.index]),
        ranges(gr.exclusion[exclusion.index])
      )
    }

    if (length(residual.range) == 0L) {
      next
    }

    fragment.rows[[anchor.index]] <- tibble(
      source_anchor_index = anchor.index,
      loop_id = as.character(mcols(gr.anchor)$loop_id[anchor.index]),
      resolution = as.character(
        mcols(gr.anchor)$resolution[anchor.index]
      ),
      anchor_side = as.character(
        mcols(gr.anchor)$anchor_side[anchor.index]
      ),
      fragment_index = seq_along(residual.range),
      fragment_chr = as.character(seqnames(gr.anchor)[anchor.index]),
      fragment_start = start(residual.range),
      fragment_end = end(residual.range),
      fragment_width_bp = width(residual.range)
    )
  }

  bind_rows(fragment.rows)
}

# Compare paired binary overlap states at the promoter and opposite anchors of
# the same loop and return rates, discordant cells, and a corrected McNemar test.
summarise_paired_binary_overlap <- function(
  df,
  promoter.column,
  candidate.column,
  resolution.label,
  comparison.label
) {
  check_required_columns(
    df,
    c(promoter.column, candidate.column),
    "Paired anchor-overlap table"
  )

  promoter.state <- as.logical(df[[promoter.column]])
  candidate.state <- as.logical(df[[candidate.column]])
  complete.pair <- !is.na(promoter.state) & !is.na(candidate.state)
  promoter.state <- promoter.state[complete.pair]
  candidate.state <- candidate.state[complete.pair]

  paired.table <- table(
    promoter = factor(promoter.state, levels = c(FALSE, TRUE)),
    candidate = factor(candidate.state, levels = c(FALSE, TRUE))
  )
  promoter.only <- unname(paired.table["TRUE", "FALSE"])
  candidate.only <- unname(paired.table["FALSE", "TRUE"])
  mcnemar.p <- if ((promoter.only + candidate.only) == 0L) {
    NA_real_
  } else {
    suppressWarnings(
      stats::mcnemar.test(paired.table, correct = TRUE)$p.value
    )
  }

  tibble(
    resolution = resolution.label,
    comparison = comparison.label,
    promoter_column = promoter.column,
    candidate_column = candidate.column,
    n_pairs_input = nrow(df),
    n_pairs_complete = length(promoter.state),
    n_both_positive = unname(paired.table["TRUE", "TRUE"]),
    n_promoter_only_positive = promoter.only,
    n_candidate_only_positive = candidate.only,
    n_neither_positive = unname(paired.table["FALSE", "FALSE"]),
    pct_promoter_positive = if (length(promoter.state) > 0L) {
      round(100 * mean(promoter.state), 1)
    } else {
      NA_real_
    },
    pct_candidate_positive = if (length(candidate.state) > 0L) {
      round(100 * mean(candidate.state), 1)
    } else {
      NA_real_
    },
    candidate_minus_promoter_percentage_points =
      pct_candidate_positive - pct_promoter_positive,
    mcnemar_p = mcnemar.p,
    test_note = paste0(
      "Paired anchors from the same single-direct-promoter/TSS loop; ",
      "continuity-corrected McNemar test."
    )
  )
}

################################################################################
# Revision: transcript-position functions
################################################################################

# Expand each loop-anchor-gene assignment to every annotated Ensembl transcript
# and flag how each transcript body, TSS, and TES relate to the loop boundaries.
annotate_assignment_transcript_positions <- function(
  df.assignment,
  df.transcript,
  df.loop,
  assignment.tier,
  promoter.window.flank.bp = 1000L
) {
  check_required_columns(
    df.assignment,
    c("loop_id", "resolution", "anchor_side", "gene_id"),
    "Promoter/TSS gene-assignment table"
  )
  check_required_columns(
    df.transcript,
    c(
      "chr", "gene_id", "gene_name", "transcript_id_versioned",
      "transcript_start", "transcript_end", "strand", "true_tss_start"
    ),
    "Strand-aware transcript annotation"
  )
  check_required_columns(
    df.loop,
    c(
      "loop_id", "chr1", "start1", "end1", "chr2", "start2", "end2"
    ),
    "Pooled loop table"
  )

  assignment.tier <- as.character(assignment.tier)
  if (length(assignment.tier) != 1L || is.na(assignment.tier)) {
    stop("Assignment tier must be one non-missing label.", call. = FALSE)
  }
  promoter.window.flank.bp <- as_integer_coordinate(
    promoter.window.flank.bp,
    "Transcript promoter-window flank"
  )
  if (
    length(promoter.window.flank.bp) != 1L ||
      promoter.window.flank.bp < 0L
  ) {
    stop(
      "Transcript promoter-window flank must be one non-negative integer.",
      call. = FALSE
    )
  }

  df.assignment.base <- df.assignment %>%
    mutate(
      assignment_tier = assignment.tier,
      assignment_id = str_c(
        assignment_tier,
        loop_id,
        anchor_side,
        gene_id,
        sep = "|"
      ),
      .before = 1
    )

  if (anyDuplicated(df.assignment.base$assignment_id)) {
    stop(
      "Loop-anchor-gene assignment identifiers are not unique.",
      call. = FALSE
    )
  }

  df.loop.lookup <- df.loop %>%
    transmute(
      loop_id,
      loop_chr = chr1,
      loop_is_intrachromosomal = chr1 == chr2,
      start1,
      end1,
      start2,
      end2,
      loop_span_start = pmin(start1, start2),
      loop_span_end = pmax(end1, end2),
      left_anchor_end = if_else(start1 <= start2, end1, end2),
      right_anchor_start = if_else(start1 <= start2, start2, start1),
      inter_anchor_start = left_anchor_end + 1L,
      inter_anchor_end = right_anchor_start - 1L,
      inter_anchor_interval_exists = inter_anchor_start <= inter_anchor_end
    )

  df.transcript.lookup <- df.transcript %>%
    transmute(
      gene_id,
      transcript_gene_name = gene_name,
      transcript_chr = chr,
      transcript_id_versioned,
      transcript_start,
      transcript_end,
      transcript_strand = strand,
      true_tss = true_tss_start,
      true_tes = if_else(
        transcript_strand == "+",
        transcript_end,
        transcript_start
      )
    ) %>%
    distinct()

  df.assignment.base %>%
    left_join(df.loop.lookup, by = "loop_id") %>%
    left_join(
      df.transcript.lookup,
      by = "gene_id",
      relationship = "many-to-many"
    ) %>%
    mutate(
      has_ensembl_transcript_annotation = !is.na(transcript_id_versioned),
      transcript_matches_loop_chromosome = (
        has_ensembl_transcript_annotation &
          loop_is_intrachromosomal &
          transcript_chr == loop_chr
      ),
      assigned_anchor_start = case_when(
        anchor_side == "anchor1" ~ start1,
        anchor_side == "anchor2" ~ start2,
        TRUE ~ NA_integer_
      ),
      assigned_anchor_end = case_when(
        anchor_side == "anchor1" ~ end1,
        anchor_side == "anchor2" ~ end2,
        TRUE ~ NA_integer_
      ),
      transcript_tss_overlaps_assigned_anchor = (
        transcript_matches_loop_chromosome &
          true_tss >= assigned_anchor_start &
          true_tss <= assigned_anchor_end
      ),
      transcript_promoter_window_start = if_else(
        has_ensembl_transcript_annotation,
        pmax(1L, true_tss - promoter.window.flank.bp),
        NA_integer_
      ),
      transcript_promoter_window_end = if_else(
        has_ensembl_transcript_annotation,
        true_tss + promoter.window.flank.bp,
        NA_integer_
      ),
      transcript_promoter_window_overlaps_assigned_anchor = (
        transcript_matches_loop_chromosome &
          transcript_promoter_window_end >= assigned_anchor_start &
          transcript_promoter_window_start <= assigned_anchor_end
      ),
      transcript_tes_overlaps_assigned_anchor = (
        transcript_matches_loop_chromosome &
          true_tes >= assigned_anchor_start &
          true_tes <= assigned_anchor_end
      ),
      transcript_overlaps_loop_span = (
        transcript_matches_loop_chromosome &
          transcript_end >= loop_span_start &
          transcript_start <= loop_span_end
      ),
      transcript_fully_within_loop_span = (
        transcript_matches_loop_chromosome &
          transcript_start >= loop_span_start &
          transcript_end <= loop_span_end
      ),
      transcript_overlaps_inter_anchor_interval = (
        transcript_matches_loop_chromosome &
          inter_anchor_interval_exists &
          transcript_end >= inter_anchor_start &
          transcript_start <= inter_anchor_end
      ),
      transcript_fully_within_inter_anchor_interval = (
        transcript_matches_loop_chromosome &
          inter_anchor_interval_exists &
          transcript_start >= inter_anchor_start &
          transcript_end <= inter_anchor_end
      ),
      transcript_crosses_left_loop_boundary = (
        transcript_matches_loop_chromosome &
          transcript_start < loop_span_start &
          transcript_end >= loop_span_start
      ),
      transcript_crosses_right_loop_boundary = (
        transcript_matches_loop_chromosome &
          transcript_start <= loop_span_end &
          transcript_end > loop_span_end
      ),
      transcript_spans_both_loop_boundaries = (
        transcript_matches_loop_chromosome &
          transcript_start < loop_span_start &
          transcript_end > loop_span_end
      ),
      transcript_tss_within_loop_span = (
        transcript_matches_loop_chromosome &
          true_tss >= loop_span_start &
          true_tss <= loop_span_end
      ),
      transcript_tes_within_loop_span = (
        transcript_matches_loop_chromosome &
          true_tes >= loop_span_start &
          true_tes <= loop_span_end
      ),
      transcript_position_class = case_when(
        !has_ensembl_transcript_annotation ~
          "no_Ensembl_transcript_annotation",
        !transcript_matches_loop_chromosome ~
          "transcript_on_different_chromosome",
        transcript_fully_within_inter_anchor_interval ~
          "fully_within_inter_anchor_interval",
        transcript_fully_within_loop_span ~
          "fully_within_loop_span_including_anchors",
        transcript_spans_both_loop_boundaries ~
          "spans_both_loop_boundaries",
        transcript_crosses_left_loop_boundary ~
          "crosses_left_loop_boundary",
        transcript_crosses_right_loop_boundary ~
          "crosses_right_loop_boundary",
        transcript_overlaps_loop_span ~ "partially_overlaps_loop_span",
        TRUE ~ "outside_loop_span"
      )
    )
}

# Collapse transcript-level position flags without selecting one isoform or
# using containment as a requirement for retaining a gene assignment.
summarise_assignment_transcript_positions <- function(
  df.position.detail,
  df.assignment
) {
  check_required_columns(
    df.position.detail,
    c(
      "assignment_id", "transcript_id_versioned",
      "has_ensembl_transcript_annotation",
      "transcript_fully_within_loop_span",
      "transcript_fully_within_inter_anchor_interval",
      "transcript_overlaps_loop_span",
      "transcript_tss_overlaps_assigned_anchor",
      "transcript_promoter_window_overlaps_assigned_anchor",
      "transcript_tes_within_loop_span",
      "transcript_crosses_left_loop_boundary",
      "transcript_crosses_right_loop_boundary"
    ),
    "Transcript-position detail table"
  )

  df.assignment.base <- df.assignment %>%
    mutate(
      assignment_id = str_c(
        assignment_tier,
        loop_id,
        anchor_side,
        gene_id,
        sep = "|"
      ),
      .before = 1
    )

  df.position.summary <- df.position.detail %>%
    group_by(assignment_id) %>%
    summarise(
      n_ensembl_transcripts = n_distinct(
        transcript_id_versioned[has_ensembl_transcript_annotation],
        na.rm = TRUE
      ),
      n_transcripts_tss_overlapping_assigned_anchor = sum(
        transcript_tss_overlaps_assigned_anchor,
        na.rm = TRUE
      ),
      n_transcript_promoter_windows_overlapping_assigned_anchor = sum(
        transcript_promoter_window_overlaps_assigned_anchor,
        na.rm = TRUE
      ),
      n_transcripts_fully_within_loop_span = sum(
        transcript_fully_within_loop_span,
        na.rm = TRUE
      ),
      n_transcripts_fully_within_inter_anchor_interval = sum(
        transcript_fully_within_inter_anchor_interval,
        na.rm = TRUE
      ),
      n_transcripts_overlapping_loop_span = sum(
        transcript_overlaps_loop_span,
        na.rm = TRUE
      ),
      n_transcripts_with_tes_in_loop_span = sum(
        transcript_tes_within_loop_span,
        na.rm = TRUE
      ),
      n_transcripts_crossing_left_loop_boundary = sum(
        transcript_crosses_left_loop_boundary,
        na.rm = TRUE
      ),
      n_transcripts_crossing_right_loop_boundary = sum(
        transcript_crosses_right_loop_boundary,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    mutate(
      has_ensembl_transcript_annotation = n_ensembl_transcripts > 0L,
      any_transcript_fully_within_loop_span =
        n_transcripts_fully_within_loop_span > 0L,
      all_transcripts_fully_within_loop_span = (
        n_ensembl_transcripts > 0L &
          n_transcripts_fully_within_loop_span == n_ensembl_transcripts
      ),
      any_transcript_fully_within_inter_anchor_interval =
        n_transcripts_fully_within_inter_anchor_interval > 0L,
      any_transcript_overlaps_loop_span =
        n_transcripts_overlapping_loop_span > 0L,
      any_transcript_tss_overlaps_assigned_anchor =
        n_transcripts_tss_overlapping_assigned_anchor > 0L,
      any_transcript_promoter_window_overlaps_assigned_anchor =
        n_transcript_promoter_windows_overlapping_assigned_anchor > 0L,
      any_transcript_tes_within_loop_span =
        n_transcripts_with_tes_in_loop_span > 0L,
      transcript_containment_summary = case_when(
        n_ensembl_transcripts == 0L ~
          "no_Ensembl_transcript_annotation",
        all_transcripts_fully_within_loop_span ~
          "all_transcripts_fully_within_loop_span",
        any_transcript_fully_within_loop_span ~
          "some_transcripts_fully_within_loop_span",
        any_transcript_overlaps_loop_span ~
          "overlap_without_full_transcript_containment",
        TRUE ~ "all_annotated_transcripts_outside_loop_span"
      )
    )

  df.assignment.base %>%
    left_join(df.position.summary, by = "assignment_id") %>%
    arrange(loop_id, anchor_side, gene_id)
}

# Select the promoter-side or opposite candidate-regulatory-side anchor from
# each directional promoter/TSS assignment and return the selected anchors as
# GRanges. WHERE determines whether anchor1 or anchor2 has the requested role.
create_anchor_granges <- function(df, anchor.role) {
  if (anchor.role == "promoter") {
    df.anchor <- bind_rows(
      df %>%
        filter(WHERE == "UP") %>%
        transmute(
          loop_id,
          chr = chr1,
          start = start1,
          end = end1,
          resolution,
          component,
          gene_id,
          gene_name,
          WHERE,
          anchor_side = "anchor1"
        ),
      df %>%
        filter(WHERE == "DOWN") %>%
        transmute(
          loop_id,
          chr = chr2,
          start = start2,
          end = end2,
          resolution,
          component,
          gene_id,
          gene_name,
          WHERE,
          anchor_side = "anchor2"
        )
    )
  } else if (anchor.role == "candidate_regulatory") {
    df.anchor <- bind_rows(
      df %>%
        filter(WHERE == "UP") %>%
        transmute(
          loop_id,
          chr = chr2,
          start = start2,
          end = end2,
          resolution,
          component,
          gene_id,
          gene_name,
          WHERE,
          anchor_side = "anchor2"
        ),
      df %>%
        filter(WHERE == "DOWN") %>%
        transmute(
          loop_id,
          chr = chr1,
          start = start1,
          end = end1,
          resolution,
          component,
          gene_id,
          gene_name,
          WHERE,
          anchor_side = "anchor1"
        )
    )
  } else {
    stop("Unknown anchor role: ", anchor.role, call. = FALSE)
  }

  GRanges(
    seqnames = df.anchor$chr,
    ranges = IRanges(
      start = df.anchor$start,
      end = df.anchor$end
    ),
    loop_id = df.anchor$loop_id,
    resolution = df.anchor$resolution,
    component = df.anchor$component,
    gene_id = df.anchor$gene_id,
    gene_name = df.anchor$gene_name,
    WHERE = df.anchor$WHERE,
    anchor_role = anchor.role,
    anchor_side = df.anchor$anchor_side
  )
}

# Calculate ordered pairwise loop sharing for all requested strains. For each
# pair, returns loop counts, intersection, union, Jaccard similarity, and
# Jaccard distance for the supplied resolution stratum.
create_pairwise_loop_overlap <- function(
  df.loop.presence,
  strain.levels,
  resolution.label = "ALL"
) {
  loop.sets <- split(df.loop.presence$loop_id, df.loop.presence$strain)
  pair.rows <- vector("list", length(strain.levels)^2)
  row.index <- 1L

  for (strain.1 in strain.levels) {
    loops.1 <- unique(loop.sets[[strain.1]])
    if (is.null(loops.1)) {
      loops.1 <- character(0)
    }

    for (strain.2 in strain.levels) {
      loops.2 <- unique(loop.sets[[strain.2]])
      if (is.null(loops.2)) {
        loops.2 <- character(0)
      }

      n.shared <- length(intersect(loops.1, loops.2))
      n.union <- length(union(loops.1, loops.2))

      pair.rows[[row.index]] <- tibble(
        resolution = resolution.label,
        strain1 = strain.1,
        strain2 = strain.2,
        n_loops_strain1 = length(loops.1),
        n_loops_strain2 = length(loops.2),
        n_shared_loops = n.shared,
        n_union_loops = n.union,
        jaccard_similarity = if (n.union > 0) n.shared / n.union else NA_real_,
        jaccard_distance = if (n.union > 0) 1 - (n.shared / n.union) else NA_real_
      )

      row.index <- row.index + 1L
    }
  }

  bind_rows(pair.rows)
}

# Summarise promoter-anchor and candidate-regulatory-anchor ATAC support for a
# loop set. Reports both all-loop percentages and the percentage among anchors
# that retain a non-TSS fragment eligible for regulatory ATAC testing.
summarise_atac_support <- function(df, resolution.label) {
  tibble(
    resolution = resolution.label,
    n_promoter_tss_candidates = nrow(df),
    n_promoter_anchor_atac = count_true(df$atac_promoter_anchor),
    pct_promoter_anchor_atac = percent_true(df$atac_promoter_anchor),
    n_candidate_regulatory_anchor_atac = count_true(
      df$atac_candidate_regulatory_anchor
    ),
    pct_candidate_regulatory_anchor_atac = percent_true(
      df$atac_candidate_regulatory_anchor
    ),
    n_candidate_regulatory_anchor_with_non_tss_fragment = count_true(
      df$candidate_regulatory_anchor_has_non_tss_fragment
    ),
    n_candidate_regulatory_anchor_non_tss_atac = count_true(
      df$atac_candidate_regulatory_anchor_non_tss
    ),
    pct_candidate_regulatory_anchor_non_tss_atac_all = percent_true(
      df$atac_candidate_regulatory_anchor_non_tss
    ),
    pct_candidate_regulatory_anchor_non_tss_atac_eligible = {
      eligible <- df$candidate_regulatory_anchor_has_non_tss_fragment %in% TRUE
      if (sum(eligible) == 0) {
        NA_real_
      } else {
        round(
          100 * mean(
            df$atac_candidate_regulatory_anchor_non_tss[eligible] %in% TRUE
          ),
          1
        )
      }
    }
  )
}

# Compare paired ATAC states at the promoter and candidate-regulatory anchors
# of the same loop. Returns the four paired cells and a continuity-corrected
# McNemar p-value; no discordant pairs produce an undefined p-value.
summarise_mcnemar <- function(df, resolution.label) {
  paired.table <- table(
    promoter_ATAC = factor(
      df$atac_promoter_anchor,
      levels = c(FALSE, TRUE)
    ),
    candidate_regulatory_ATAC = factor(
      df$atac_candidate_regulatory_anchor,
      levels = c(FALSE, TRUE)
    )
  )

  promoter.only <- unname(paired.table["TRUE", "FALSE"])
  candidate.only <- unname(paired.table["FALSE", "TRUE"])

  mcnemar.p <- if ((promoter.only + candidate.only) == 0) {
    NA_real_
  } else {
    suppressWarnings(stats::mcnemar.test(
      paired.table,
      correct = TRUE
    )$p.value)
  }

  tibble(
    resolution = resolution.label,
    n_pairs = nrow(df),
    n_both_atac = unname(paired.table["TRUE", "TRUE"]),
    n_promoter_only_atac = promoter.only,
    n_candidate_regulatory_only_atac = candidate.only,
    n_neither_atac = unname(paired.table["FALSE", "FALSE"]),
    mcnemar_p = mcnemar.p,
    test_note = paste0(
      "Paired promoter and candidate-regulatory anchors from the same loop; ",
      "continuity-corrected McNemar test."
    )
  )
}

# Read a long strain-pair genetic-distance table, accept a supported
# distance-column name, canonicalise pair order, and retain one finite value
# per unordered strain pair.
read_genetic_distance_long <- function(file) {
  df.distance <- readr::read_tsv(file, show_col_types = FALSE)
  required.pair.columns <- c("strain1", "strain2")
  check_required_columns(
    df.distance,
    required.pair.columns,
    "Genetic-distance table"
  )

  distance.column <- intersect(
    c(
      "genetic_distance",
      "plink2_ibs_distance",
      "ibs_distance",
      "distance"
    ),
    colnames(df.distance)
  )

  if (length(distance.column) == 0) {
    stop(
      "The genetic-distance table must contain genetic_distance, ",
      "plink2_ibs_distance, ibs_distance, or distance.",
      call. = FALSE
    )
  }

  df.distance %>%
    transmute(
      strain1 = as.character(strain1),
      strain2 = as.character(strain2),
      genetic_distance = as.numeric(.data[[distance.column[1]]]),
      pair_strain_a = pmin(strain1, strain2),
      pair_strain_b = pmax(strain1, strain2)
    ) %>%
    filter(
      strain1 != strain2,
      !is.na(genetic_distance),
      is.finite(genetic_distance)
    ) %>%
    distinct(pair_strain_a, pair_strain_b, .keep_all = TRUE)
}

################################################################################
# Revision: layered CTCF sensitivity and multi-gene resource functions
################################################################################

# Summarise one fixed per-anchor CTCF motif threshold while treating revised
# regulatory and promoter-promoter classes as independent annotation layers.
summarise_layered_ctcf_threshold <- function(
  df,
  threshold,
  resolution.label,
  regulatory.column = "revised_putative_regulatory_support",
  promoter.promoter.column = "revised_promoter_promoter_compatible"
) {
  check_required_columns(
    df,
    c(
      "ctcf_count_anchor1", "ctcf_count_anchor2",
      regulatory.column, promoter.promoter.column
    ),
    "Revised loop-evidence table"
  )

  threshold <- as.integer(threshold)
  passes.threshold <- (
    df$ctcf_count_anchor1 >= threshold &
      df$ctcf_count_anchor2 >= threshold
  )
  regulatory.state <- df[[regulatory.column]] %in% TRUE
  promoter.promoter.state <- df[[promoter.promoter.column]] %in% TRUE
  n.regulatory <- sum(regulatory.state)
  n.promoter.promoter <- sum(promoter.promoter.state)

  tibble(
    resolution = resolution.label,
    ctcf_min_threshold = threshold,
    n_loops_universe = nrow(df),
    n_passes_ctcf_threshold = sum(passes.threshold),
    pct_pooled_loops_passing_ctcf_threshold = round(
      100 * mean(passes.threshold),
      1
    ),
    n_revised_putative_regulatory = n.regulatory,
    n_revised_putative_regulatory_passing_ctcf = sum(
      regulatory.state & passes.threshold
    ),
    pct_revised_putative_regulatory_passing_ctcf = round(
      100 * safe_ratio(
        sum(regulatory.state & passes.threshold),
        n.regulatory
      ),
      1
    ),
    n_promoter_promoter_compatible = n.promoter.promoter,
    n_promoter_promoter_passing_ctcf = sum(
      promoter.promoter.state & passes.threshold
    ),
    pct_promoter_promoter_passing_ctcf = round(
      100 * safe_ratio(
        sum(promoter.promoter.state & passes.threshold),
        n.promoter.promoter
      ),
      1
    )
  )
}

# Compare the fixed >=6 motif rule with the anchor-width-scaled sensitivity
# rule within the revised putative-regulatory set and complete pooled universe.
summarise_layered_ctcf_resolution_rule <- function(
  df,
  resolution.label,
  regulatory.column = "revised_putative_regulatory_support"
) {
  check_required_columns(
    df,
    c(
      "passes_ctcf_ge6_both_anchors",
      "passes_ctcf_resolution_adjusted_both_anchors",
      "ctcf_resolution_adjusted_min_threshold",
      regulatory.column
    ),
    "Revised loop-evidence table"
  )

  fixed.pass <- df$passes_ctcf_ge6_both_anchors %in% TRUE
  adjusted.pass <-
    df$passes_ctcf_resolution_adjusted_both_anchors %in% TRUE
  regulatory.state <- df[[regulatory.column]] %in% TRUE
  n.regulatory <- sum(regulatory.state)
  adjusted.rule <- if (resolution.label == "ALL") {
    "5K>=6; 10K>=12; 25K>=30"
  } else {
    str_c(
      resolution.label,
      ">=",
      sort(unique(df$ctcf_resolution_adjusted_min_threshold))[1]
    )
  }

  tibble(
    resolution = resolution.label,
    resolution_adjusted_rule = adjusted.rule,
    n_loops_universe = nrow(df),
    n_fixed_ge6_both = sum(fixed.pass),
    pct_fixed_ge6_both = round(100 * mean(fixed.pass), 1),
    n_resolution_adjusted_both = sum(adjusted.pass),
    pct_resolution_adjusted_both = round(
      100 * mean(adjusted.pass),
      1
    ),
    n_pass_both_rules = sum(fixed.pass & adjusted.pass),
    n_fixed_ge6_only = sum(fixed.pass & !adjusted.pass),
    n_resolution_adjusted_only = sum(!fixed.pass & adjusted.pass),
    n_fail_both_rules = sum(!fixed.pass & !adjusted.pass),
    n_revised_putative_regulatory = n.regulatory,
    n_revised_putative_regulatory_fixed_ge6 = sum(
      regulatory.state & fixed.pass
    ),
    pct_revised_putative_regulatory_fixed_ge6 = round(
      100 * safe_ratio(
        sum(regulatory.state & fixed.pass),
        n.regulatory
      ),
      1
    ),
    n_revised_putative_regulatory_resolution_adjusted = sum(
      regulatory.state & adjusted.pass
    ),
    pct_revised_putative_regulatory_resolution_adjusted = round(
      100 * safe_ratio(
        sum(regulatory.state & adjusted.pass),
        n.regulatory
      ),
      1
    )
  )
}

# Count distinct loops per stable Ensembl gene within each revised gene set;
# multiple annotation sources or anchor records never inflate loop counts.
summarise_revised_gene_loop_counts <- function(df.gene.loop.membership) {
  check_required_columns(
    df.gene.loop.membership,
    c("gene_set", "loop_id", "ensembl_gene_id", "gene_symbol"),
    "Revised gene-loop membership table"
  )

  df.gene.loop.membership %>%
    filter(!is.na(ensembl_gene_id)) %>%
    distinct(gene_set, loop_id, ensembl_gene_id, .keep_all = TRUE) %>%
    group_by(gene_set, ensembl_gene_id) %>%
    summarise(
      gene_symbol = {
        symbol.values <- sort(unique(na.omit(gene_symbol)))
        if (length(symbol.values) == 0L) {
          NA_character_
        } else {
          symbol.values[1]
        }
      },
      n = n_distinct(loop_id),
      .groups = "drop"
    ) %>%
    arrange(gene_set, desc(n), gene_symbol, ensembl_gene_id)
}

# Select unique Ensembl genes from a named loop-derived gene set that meet a
# minimum loop-count threshold. Returns genes ordered by loop count and symbol.
prepare_go_gene_set <- function(df.gene.count, set.name, minimum.loops) {
  df.gene.count %>%
    filter(
      gene_set == set.name,
      n >= minimum.loops,
      !is.na(ensembl_gene_id)
    ) %>%
    distinct(ensembl_gene_id, gene_symbol, n) %>%
    arrange(desc(n), gene_symbol)
}

# Run rat Gene Ontology Biological Process enrichment for one selected gene
# set against an explicit Ensembl-gene universe. Missing packages or too few
# mapped genes cause a documented skip; successful results are written to TSV.
run_go_enrichment <- function(df.gene, universe.ensembl, set.name, out.dir) {
  required.packages <- c(
    "clusterProfiler",
    "org.Rn.eg.db",
    "AnnotationDbi"
  )
  missing.packages <- required.packages[
    !vapply(required.packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing.packages) > 0) {
    message(
      "Skipping GO for ",
      set.name,
      ": missing ",
      paste(missing.packages, collapse = ", "),
      "."
    )
    return(tibble())
  }

  ensembl.genes <- unique(na.omit(df.gene$ensembl_gene_id))
  if (length(ensembl.genes) < 5) {
    message(
      "Skipping GO for ",
      set.name,
      ": fewer than 5 Ensembl genes."
    )
    return(tibble())
  }

  df.gene.map <- AnnotationDbi::select(
    org.Rn.eg.db::org.Rn.eg.db,
    keys = ensembl.genes,
    keytype = "ENSEMBL",
    columns = c("ENTREZID", "SYMBOL")
  ) %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(ENTREZID, .keep_all = TRUE)

  df.universe.map <- AnnotationDbi::select(
    org.Rn.eg.db::org.Rn.eg.db,
    keys = unique(na.omit(universe.ensembl)),
    keytype = "ENSEMBL",
    columns = c("ENTREZID", "SYMBOL")
  ) %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(ENTREZID, .keep_all = TRUE)

  if (nrow(df.gene.map) < 5) {
    message(
      "Skipping GO for ",
      set.name,
      ": fewer than 5 mapped Entrez genes."
    )
    return(tibble())
  }

  go.result <- clusterProfiler::enrichGO(
    gene = df.gene.map$ENTREZID,
    universe = df.universe.map$ENTREZID,
    OrgDb = org.Rn.eg.db::org.Rn.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  )

  df.go.result <- as_tibble(go.result@result) %>%
    mutate(
      gene_set = set.name,
      input_ensembl_genes = length(ensembl.genes),
      mapped_entrez_genes = nrow(df.gene.map),
      universe_mapped_entrez_genes = nrow(df.universe.map)
    ) %>%
    arrange(p.adjust, pvalue)

  readr::write_tsv(
    df.go.result,
    file.path(out.dir, paste0("go_enrichment_BP_", set.name, ".tsv"))
  )

  df.go.result
}

# Stop with a descriptive message unless one scalar analysis invariant is true.
# Optionally prints a success message if condition evaluates to TRUE.
assert_analysis_condition <- function(condition, message, success.message = NULL) {
  if (length(condition) != 1L || is.na(condition) || !isTRUE(condition)) {
    stop(message, call. = FALSE)
  }

  if (!is.null(success.message) && nzchar(success.message)) {
    message(success.message)
  }

  invisible(TRUE)
}

# Require an analysis table to contain exactly the expected number of rows.
assert_analysis_row_count <- function(df, expected, message, success.message = NULL) {
  assert_analysis_condition(
    nrow(df) == as.integer(expected),
    message = message,
    success.message = success.message
  )
}

# Require one unique row for each key combination in an analysis table.
assert_analysis_unique_key <- function(df, columns, message) {
  check_required_columns(df, columns, "Analysis key validation")
  assert_analysis_condition(
    !anyDuplicated(df[columns]),
    message
  )
}

# Convert both sides of a normalized loop table into a named GRanges list.
create_loop_anchor_granges_by_side <- function(df.loop) {
  anchor.sides <- c("anchor1", "anchor2")
  set_names(
    map(anchor.sides, ~ create_loop_anchor_granges(df.loop, .x)),
    anchor.sides
  )
}

# Require a named two-element anchor list before applying shared side-wise logic.
validate_loop_anchor_granges_by_side <- function(gr.loop.anchor.by.side) {
  expected.sides <- c("anchor1", "anchor2")
  assert_analysis_condition(
    identical(names(gr.loop.anchor.by.side), expected.sides) &&
      all(map_lgl(gr.loop.anchor.by.side, ~ methods::is(.x, "GRanges"))),
    "Loop-anchor GRanges must be named anchor1 and anchor2."
  )

  invisible(TRUE)
}

# Combine the two named anchor GRanges in deterministic anchor1-then-anchor2 order.
combine_loop_anchor_granges_by_side <- function(gr.loop.anchor.by.side) {
  validate_loop_anchor_granges_by_side(gr.loop.anchor.by.side)
  do.call(c, unname(gr.loop.anchor.by.side))
}

# Convert a normalized loop table into one long row for each anchor side.
build_loop_anchor_index <- function(df.loop) {
  check_required_columns(
    df.loop,
    c(
      "loop_id", "resolution", "chr1", "start1", "end1",
      "chr2", "start2", "end2"
    ),
    "Loop anchor index input"
  )

  map_dfr(c("anchor1", "anchor2"), function(anchor.side) {
    side.number <- if (anchor.side == "anchor1") 1L else 2L
    opposite.side <- if (anchor.side == "anchor1") "anchor2" else "anchor1"

    df.loop %>%
      transmute(
        loop_id,
        resolution,
        anchor_side = anchor.side,
        opposite_anchor_side = opposite.side,
        anchor_chr = .data[[paste0("chr", side.number)]],
        anchor_start = .data[[paste0("start", side.number)]],
        anchor_end = .data[[paste0("end", side.number)]]
      )
  })
}

# Apply one interval operation to both loop-anchor GRanges and bind the rows.
map_loop_anchors_dfr <- function(gr.loop.anchor.by.side, .f, ...) {
  validate_loop_anchor_granges_by_side(gr.loop.anchor.by.side)

  map_dfr(gr.loop.anchor.by.side, .f, ...)
}

# Count feature overlaps on both anchor sides and return one wide loop table.
count_loop_anchor_feature_overlaps <- function(
  gr.loop.anchor.by.side,
  gr.feature,
  count.prefix,
  ...
) {
  validate_loop_anchor_granges_by_side(gr.loop.anchor.by.side)
  expected.sides <- names(gr.loop.anchor.by.side)

  map_dfr(expected.sides, function(anchor.side) {
    gr.anchor <- gr.loop.anchor.by.side[[anchor.side]]
    tibble(
      loop_id = as.character(mcols(gr.anchor)$loop_id),
      anchor_side = anchor.side,
      feature_count = countOverlaps(gr.anchor, gr.feature, ...)
    )
  }) %>%
    pivot_wider(
      names_from = anchor_side,
      values_from = feature_count,
      names_glue = paste0(count.prefix, "_{anchor_side}")
    )
}

# Write a named list of data frames using each list name as the output filename.
write_named_tsv_tables <- function(tables, output.dir) {
  output.names <- names(tables)
  assert_analysis_condition(
    !is.null(output.names) &&
      all(nzchar(output.names)) &&
      !anyDuplicated(output.names),
    "TSV output tables must have unique, non-empty filenames."
  )

  walk2(
    tables,
    output.names,
    ~ readr::write_tsv(.x, file.path(output.dir, .y))
  )

  invisible(file.path(output.dir, output.names))
}

# Resolve a filename/object-name registry into a named list of output tables.
resolve_output_table_registry <- function(registry, envir = parent.frame()) {
  check_required_columns(
    registry,
    c("output_file", "object_name"),
    "Output table registry"
  )
  assert_analysis_condition(
    all(nzchar(registry$output_file)) &&
      all(nzchar(registry$object_name)) &&
      !anyDuplicated(registry$output_file),
    "Output registry filenames and object names must be non-empty, and filenames must be unique."
  )

  missing.objects <- registry$object_name[
    !vapply(
      registry$object_name,
      exists,
      logical(1),
      envir = envir,
      inherits = TRUE
    )
  ]
  assert_analysis_condition(
    length(missing.objects) == 0L,
    paste0(
      "Output registry references missing object(s): ",
      paste(unique(missing.objects), collapse = ", ")
    )
  )

  set_names(
    mget(registry$object_name, envir = envir, inherits = TRUE),
    registry$output_file
  )
}

# Write a named list of character vectors using each name as the output filename.
write_named_line_vectors <- function(vectors, output.dir) {
  output.names <- names(vectors)
  assert_analysis_condition(
    !is.null(output.names) &&
      all(nzchar(output.names)) &&
      !anyDuplicated(output.names),
    "Text output vectors must have unique, non-empty filenames."
  )

  walk2(
    vectors,
    output.names,
    ~ readr::write_lines(as.character(.x), file.path(output.dir, .y))
  )

  invisible(file.path(output.dir, output.names))
}

# Revision-compatible loop-count analysis: deduplicate loop/gene pairs, count
# distinct interactions per gene, and return reusable multiple-interaction and
# plotting objects while preserving the legacy function name and core inputs.
approach_2nd_analyze_loops_by_threshold <- function(
  df,
  threshold_distance = 2e5,
  top_n_genes = 80L,
  print_top_n = 50L,
  minimum_interactions = 2L,
  create_plot = TRUE,
  print_results = TRUE
) {
  required.columns <- c("gene_id", "loop.id", "distance")
  missing.columns <- setdiff(required.columns, names(df))
  if (length(missing.columns) > 0L) {
    stop(
      "Loop-count input lacks required column(s): ",
      paste(missing.columns, collapse = ", "),
      call. = FALSE
    )
  }

  minimum_interactions <- as.integer(minimum_interactions)
  if (
    length(minimum_interactions) != 1L ||
      is.na(minimum_interactions) ||
      minimum_interactions < 1L
  ) {
    stop("minimum_interactions must be one positive integer.", call. = FALSE)
  }

  df.filtered <- df %>%
    filter(
      !is.na(gene_id),
      !is.na(loop.id),
      !is.na(distance),
      distance <= threshold_distance
    ) %>%
    distinct(gene_id, loop.id, .keep_all = TRUE)

  df.gene.loop.count <- df.filtered %>%
    count(gene_id, name = "n", sort = TRUE)

  df.multiple.interaction.genes <- df.gene.loop.count %>%
    filter(n >= minimum_interactions) %>%
    arrange(desc(n), gene_id)

  df.top.genes <- df.gene.loop.count %>%
    slice_max(n, n = top_n_genes, with_ties = TRUE) %>%
    arrange(n, gene_id)

  plot.object <- NULL
  if (create_plot && nrow(df.top.genes) > 0L) {
    plot.object <- ggplot(
      df.top.genes,
      aes(x = reorder(gene_id, n), y = n)
    ) +
      geom_col(fill = "#1F78B4") +
      geom_text(aes(label = n), vjust = -0.3, color = "red", size = 3) +
      coord_flip() +
      labs(
        title = paste(
          "Top",
          top_n_genes,
          "Genes by Number of Distinct Loops (distance <=",
          threshold_distance,
          ")"
        ),
        x = "Gene",
        y = "Number of distinct loops"
      ) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))

    print(plot.object)
  }

  total.loops <- df %>%
    filter(!is.na(loop.id)) %>%
    distinct(loop.id) %>%
    nrow()
  filtered.loops <- df.filtered %>%
    distinct(loop.id) %>%
    nrow()
  percent.retained <- if (total.loops == 0L) {
    NA_real_
  } else {
    round(100 * filtered.loops / total.loops, 2)
  }

  df.summary <- tibble(
    threshold_distance = threshold_distance,
    minimum_interactions = minimum_interactions,
    n_input_loops = total.loops,
    n_filtered_loops = filtered.loops,
    pct_loops_retained = percent.retained,
    n_genes = nrow(df.gene.loop.count),
    n_multiple_interaction_genes = nrow(df.multiple.interaction.genes),
    maximum_interactions_per_gene = if (nrow(df.gene.loop.count) == 0L) {
      0L
    } else {
      max(df.gene.loop.count$n)
    }
  )

  if (print_results) {
    message("Total loops: ", total.loops)
    message("Filtered loops: ", filtered.loops)
    message("Filtered loops retain ", percent.retained, "% of total loops")
    cat("Top", print_top_n, "genes:\n")
    print(
      df.gene.loop.count %>%
        slice_max(n, n = print_top_n, with_ties = TRUE),
      n = Inf
    )
  }

  invisible(
    list(
      filtered_loop_gene_pairs = df.filtered,
      gene_loop_count = df.gene.loop.count,
      multiple_interaction_genes = df.multiple.interaction.genes,
      top_genes = df.top.genes,
      summary = df.summary,
      plot = plot.object
    )
  )
}

################################################################################
# Revision: legacy midpoint analysis with coordinate-normalized inputs
################################################################################

# Rebuild the legacy TSS/promoter feature catalog from strand-aware Ensembl
# transcript TSSs and coordinate-normalized EPD promoters without legacy RDSs.
build_legacy_feature_catalog <- function(df.true.tss, df.promoter) {
  check_required_columns(
    df.true.tss,
    c(
      "chr", "true_tss_start", "true_tss_end", "strand", "true_tss_id",
      "gene_id", "gene_name", "gene_biotype", "transcript_id",
      "transcript_start", "transcript_end"
    ),
    "Strand-aware true TSS annotation"
  )
  check_required_columns(
    df.promoter,
    c(
      "chr", "promoter_start", "promoter_end", "strand",
      "promoter_annotation_id", "epd_promoter_name", "gene_id", "gene_name"
    ),
    "Coordinate-normalized EPD promoter annotation"
  )

  first_non_missing <- function(x) {
    x <- x[!is.na(x) & nzchar(as.character(x))]
    if (length(x) == 0L) NA_character_ else as.character(x[[1L]])
  }

  # EPD records are gene-level, so attach the full Ensembl gene span and retain
  # the number of transcripts used to construct that span.
  df.gene.span <- df.true.tss %>%
    group_by(chr, gene_id) %>%
    summarise(
      ensembl_gene_name = first_non_missing(gene_name),
      gene_biotype = first_non_missing(gene_biotype),
      gene_start = min(transcript_start, na.rm = TRUE),
      gene_end = max(transcript_end, na.rm = TRUE),
      n_gene_transcripts = n_distinct(transcript_id),
      .groups = "drop"
    )

  # Each Ensembl transcript contributes its true one-base TSS and its own
  # transcript interval for the legacy loop-containment check.
  df.tss.feature <- df.true.tss %>%
    transmute(
      chr,
      feature_start = as.integer(true_tss_start),
      feature_end = as.integer(true_tss_end),
      component_strand = strand,
      gene_id,
      gene_name,
      gene_biotype,
      transcript_id,
      gene_start = as.integer(transcript_start),
      gene_end = as.integer(transcript_end),
      n_gene_transcripts = 1L,
      component = "tss",
      component_source = "Ensembl_true_transcript_TSS",
      component_id = str_c(true_tss_id, "tss", sep = "|"),
      containment_unit = "transcript_interval",
      feature_priority = 2L
    )

  # EPD promoters lack transcript identifiers; use the union Ensembl gene span
  # for the same legacy positional screen and keep unmatched mappings explicit.
  df.promoter.feature <- df.promoter %>%
    left_join(df.gene.span, by = c("chr", "gene_id")) %>%
    transmute(
      chr,
      feature_start = as.integer(promoter_start),
      feature_end = as.integer(promoter_end),
      component_strand = strand,
      gene_id,
      gene_name = coalesce(gene_name, ensembl_gene_name, gene_id),
      gene_biotype,
      transcript_id = NA_character_,
      gene_start = as.integer(gene_start),
      gene_end = as.integer(gene_end),
      n_gene_transcripts = as.integer(n_gene_transcripts),
      component = "pro",
      component_source = "EPDnew_promoter",
      component_id = str_c(promoter_annotation_id, "pro", sep = "|"),
      containment_unit = "Ensembl_gene_span",
      feature_priority = 1L
    )

  df.feature <- bind_rows(df.tss.feature, df.promoter.feature) %>%
    filter(
      !is.na(chr),
      !is.na(feature_start),
      !is.na(feature_end),
      !is.na(gene_id),
      !is.na(component_id)
    ) %>%
    distinct(component_id, .keep_all = TRUE) %>%
    arrange(chr, feature_start, feature_end, feature_priority, component_id) %>%
    mutate(component_index = row_number(), .before = 1)

  validate_one_based_inclusive_intervals(
    df.feature$feature_start,
    df.feature$feature_end,
    "Legacy-compatible TSS/promoter features"
  )
  df.feature
}

# Convert the rebuilt legacy feature catalog into GRanges while retaining all
# metadata needed for directional nearest-feature assignment.
create_legacy_feature_granges <- function(df.feature) {
  check_required_columns(
    df.feature,
    c("chr", "feature_start", "feature_end", "component_id"),
    "Legacy-compatible feature catalog"
  )

  gr.feature <- GRanges(
    seqnames = df.feature$chr,
    ranges = IRanges(
      start = df.feature$feature_start,
      end = df.feature$feature_end
    ),
    strand = "*"
  )
  mcols(gr.feature) <- S4Vectors::DataFrame(
    df.feature %>% dplyr::select(-chr, -feature_start, -feature_end)
  )
  gr.feature
}

# Represent one normalized loop-anchor side by its midpoint and retain the
# midpoint between anchors used by the original directional eligibility rule.
create_legacy_anchor_midpoint_granges <- function(df.loop, anchor.side) {
  check_required_columns(
    df.loop,
    c(
      "loop_id", "resolution", "chr1", "start1", "end1",
      "chr2", "start2", "end2"
    ),
    "Coordinate-normalized loop table"
  )
  if (!anchor.side %in% c("anchor1", "anchor2")) {
    stop("Unknown legacy anchor side: ", anchor.side, call. = FALSE)
  }

  side.number <- if (anchor.side == "anchor1") 1L else 2L
  anchor.start <- df.loop[[paste0("start", side.number)]]
  anchor.end <- df.loop[[paste0("end", side.number)]]
  anchor.chr <- df.loop[[paste0("chr", side.number)]]
  midpoint1 <- floor((as.double(df.loop$start1) + df.loop$end1) / 2)
  midpoint2 <- floor((as.double(df.loop$start2) + df.loop$end2) / 2)
  anchor.midpoint <- floor((as.double(anchor.start) + anchor.end) / 2)
  loop.midpoint <- floor((midpoint1 + midpoint2) / 2)

  GRanges(
    seqnames = anchor.chr,
    ranges = IRanges(
      start = as.integer(anchor.midpoint),
      end = as.integer(anchor.midpoint)
    ),
    strand = "*",
    loop.id = df.loop$loop_id,
    resolution = df.loop$resolution,
    mid_loop = as.integer(loop.midpoint),
    anchor_interval_start = as.integer(anchor.start),
    anchor_interval_end = as.integer(anchor.end),
    WHERE = if (anchor.side == "anchor1") "UP" else "DOWN"
  )
}

# Find all equally nearest eligible features for one midpoint anchor side. The
# UP side admits features starting at or left of the loop midpoint; the DOWN
# side admits features ending at or right of that midpoint, matching the legacy
# directional screen without scanning every feature for every loop.
find_legacy_directional_nearest_assignments <- function(
  gr.anchor,
  gr.feature,
  where
) {
  if (!where %in% c("UP", "DOWN")) {
    stop("Legacy WHERE must be UP or DOWN.", call. = FALSE)
  }

  hit_pairs <- function(hits) {
    if (length(hits) == 0L) {
      return(tibble(query_index = integer(), subject_index = integer()))
    }
    tibble(
      query_index = queryHits(hits),
      subject_index = subjectHits(hits)
    )
  }

  # Only overlap plus the nearest feature on either side can be the nearest
  # eligible feature after applying the loop-midpoint directional restriction.
  df.pair <- bind_rows(
    hit_pairs(findOverlaps(gr.anchor, gr.feature, ignore.strand = TRUE)),
    hit_pairs(follow(
      gr.anchor,
      gr.feature,
      select = "all",
      ignore.strand = TRUE
    )),
    hit_pairs(precede(
      gr.anchor,
      gr.feature,
      select = "all",
      ignore.strand = TRUE
    ))
  ) %>%
    distinct(query_index, subject_index)

  if (nrow(df.pair) == 0L) {
    return(tibble())
  }

  feature.start <- start(gr.feature)[df.pair$subject_index]
  feature.end <- end(gr.feature)[df.pair$subject_index]
  loop.mid <- mcols(gr.anchor)$mid_loop[df.pair$query_index]
  eligible <- if (where == "UP") {
    feature.start <= loop.mid
  } else {
    feature.end >= loop.mid
  }
  df.pair <- df.pair[eligible, , drop = FALSE]

  if (nrow(df.pair) == 0L) {
    return(tibble())
  }

  df.pair$distance <- as.integer(distance(
    gr.anchor[df.pair$query_index],
    gr.feature[df.pair$subject_index],
    ignore.strand = TRUE
  ))

  df.anchor.metadata <- as_tibble(as.data.frame(mcols(gr.anchor)))[
    df.pair$query_index,
    ,
    drop = FALSE
  ]
  df.feature.metadata <- as_tibble(as.data.frame(mcols(gr.feature)))[
    df.pair$subject_index,
    ,
    drop = FALSE
  ]

  bind_cols(
    tibble(
      distance = df.pair$distance,
      loop_chr = as.character(seqnames(gr.anchor))[df.pair$query_index],
      loop_start = start(gr.anchor)[df.pair$query_index],
      loop_end = end(gr.anchor)[df.pair$query_index]
    ),
    df.anchor.metadata,
    df.feature.metadata
  ) %>%
    group_by(loop.id, WHERE) %>%
    filter(distance == min(distance)) %>%
    ungroup() %>%
    arrange(loop.id, WHERE, distance, feature_priority, gene_id, component_id)
}

# Replace the legacy manual LOC/exon-string tie edits with one deterministic
# nearest assignment per loop side while recording the original tie burden.
resolve_legacy_nearest_assignment_ties <- function(df.assignment, df.loop) {
  check_required_columns(
    df.assignment,
    c(
      "loop.id", "WHERE", "distance", "gene_id", "component_id",
      "component", "feature_priority", "gene_start", "gene_end"
    ),
    "Legacy nearest-feature assignments"
  )

  df.loop.lookup <- df.loop %>%
    transmute(
      loop.id = loop_id,
      loop_chr1 = chr1,
      x1 = start1,
      x2 = end1,
      loop_chr2 = chr2,
      y1 = start2,
      y2 = end2,
      loop_distance,
      resolution_bp
    )

  df.assignment %>%
    left_join(df.loop.lookup, by = "loop.id") %>%
    group_by(loop.id, WHERE) %>%
    mutate(
      tie.count.loop.where = n(),
      tie.gene.count.loop.where = n_distinct(gene_id),
      candidate_body_contained =
        !is.na(gene_start) &
        !is.na(gene_end) &
        loop_chr1 == loop_chr2 &
        loop_chr == loop_chr1 &
        gene_start >= x1 &
        gene_end <= y2
    ) %>%
    ungroup() %>%
    # Collapse multiple transcript/source records for the same tied gene first.
    group_by(loop.id, WHERE, gene_id) %>%
    arrange(
      desc(candidate_body_contained),
      feature_priority,
      component_id,
      .by_group = TRUE
    ) %>%
    dplyr::slice(1L) %>%
    ungroup() %>%
    # Resolve a remaining between-gene distance tie without manual gene removal.
    group_by(loop.id, WHERE) %>%
    arrange(
      desc(candidate_body_contained),
      feature_priority,
      gene_id,
      component_id,
      .by_group = TRUE
    ) %>%
    dplyr::slice(1L) %>%
    mutate(
      tie_resolution_rule = if_else(
        tie.count.loop.where == 1L,
        "unique_nearest",
        paste0(
          "contained_body_then_promoter_then_gene_id;nearest_rows=",
          tie.count.loop.where,
          ";nearest_genes=",
          tie.gene.count.loop.where
        )
      )
    ) %>%
    ungroup()
}
