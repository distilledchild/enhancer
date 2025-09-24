library(fs)

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
  read_tsv(file, skip = vcf_header_skip(file, '^##'), col_name = TRUE, guess_max = 2222) %>% 
  rename_at(ncol(.), ~"default")
}

# vcf file data reader for initial data after GLNexus
vcf_raw_data_reader <- function(file, pattern) {
  read_tsv(file, skip = vcf_header_skip(file, '^##'), col_types = cols('#CHROM' = readr::col_character()))
  # read_tsv(file, skip = vcf_header_skip(file, '^##'), col_name = TRUE, col_types = "cdcccdccccccccc")
}

# vcf file data reader for initial data splitted by chrs
vcf_by_chr_raw_data_reader <- function(file, pattern) {
  read_tsv(file, skip = vcf_header_skip(file, '^##'), col_name = TRUE, col_types = cols(.default = "d", pos = "d"))
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
    filter(POS == 151 & FILTER != 'RefCall') %>%                        # only 151 important
    mutate(CHROM = str_to_upper(CHROM)) %>% 
    filter(str_detect(CHROM, regex('Q30', ignore_case = T))) %>%        # only Q30 important
    filter(str_length(ALT) == 1 & str_length(REF) == 1) %>%             # redundant condition                                        
    filter(str_detect(CHROM, regex('snp', ignore_case = T))) %>%        # only SNP markers are important
    # # mutate(sample = str_split_n(str_split_n(sample, str_c('\\/', deli_dir, '\\/'), 2), '\\.', 1)) %>%
    mutate(sample = str_split_n(sample, str_c(deli_dir, '\\/'), 2)) %>%
    mutate(sample = if_else(str_detect(sample, "_"), str_split_n(sample, "\\_", 1), str_split_n(sample, '\\.', 1))) %>%
    mutate(GT = if_else(str_detect(default, ":"), str_split_n(default, ':', 1), default)) %>%
    mutate(DP = if_else(str_detect(INFO, "="), as.integer(str_split_n(str_split_n(INFO, 'DP=', 2), ';', 1)), as.integer(str_split_n(default, ':', 3)))) %>%
    select(-c(FORMAT, INFO)) %>%
    filter(!is.na(DP)) %>% 
    mutate(number = str_split_n(sample, 'Rat', 2)) %>% 
    left_join(meta_data, by = c("number" = "number")) %>%
    mutate(column_name = paste0(.$sample,".",.$strain, ".", .$sex)) %>% 
    separate(CHROM, into=c("mq", "mstrain", "mvt", "mcoord"), sep = '[|_]', remove=FALSE) #prefix 'm' means 'marker' in order not to be confused with q, strain, variant type of data
}

# getting sequence region
# params: chrom, start, end, species
# TODO using ... for optional parameter
marker_seq_fetcher <- function(chrom, start, end, species, mask=NULL) {
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
marker_coord_setter <- function(...){
                                
    current.row.df <- tibble(...)

    if((current.row.df$max_val - current.row.df$min_val) == 300){
        # just one marker in a chrom

        current.row.df %>% mutate(Start = Init) %>% mutate(End = current.row.df$Start - 1) %>% 
        bind_rows(current.row.df) %>% 
        bind_rows(current.row.df %>% mutate(Start = current.row.df$End + 1) %>% mutate(End = current.row.df$Terminal))

    }else if (current.row.df$Start == current.row.df$min_val){
        # 1st marker in a chrom

        current.row.df %>% mutate(Start = Init) %>% mutate(End = current.row.df$Start - 1) %>% 
        bind_rows(current.row.df)

    }else if (current.row.df$End == current.row.df$max_val){
        # the last marker in a chrom

        current.row.df %>% mutate(Start = NA) %>% mutate(End = current.row.df$End - 1) %>% 
        bind_rows(current.row.df) %>% 
        bind_rows(current.row.df %>% mutate(Start = current.row.df$End + 1) %>% mutate(End = current.row.df$Terminal))

    }else{
        # a marker
        current.row.df %>% mutate(Start = NA) %>% mutate(End = current.row.df$End - 1) %>% 
        bind_rows(current.row.df)
    }
}

# convertor to XStringSet
converting_df_2_XStringSet <- function(dfx, type){
  
  column_names <- colnames(dfx)
  
  if (!("name" %in% column_names) || !("sequence" %in% column_names)) {
    print("name or sequence column absent")
  } else {
    for(cname in column_names) {
      assign(cname, dfx[[cname]])
    }
  }
  
  names(sequence) <- name
  return (ifelse(toupper(type) == "D", DNAStringSet(sequence), 
         ifelse(toupper(type) == "R", RNAStringSet(sequence), AAStringSet(sequence))))
}

my_quantile <- function(x, probs) {
  tibble(qt = quantile(x, probs), probs = probs)
}


bedpe_data_reader <- function(file) {
  read.csv(file, header = T, sep = '\t')
}

init.bedpe.df <- function(file.list, deli_dir){
  file.list %>% 
    map_dfr(bedpe_data_reader, .id="sample") %>% 
    mutate(sample = str_split_n(sample, str_c(deli_dir, '\\/'), 2)) %>%
    mutate(sample = str_split_n(sample, '_', 1)) %>%
    mutate(distance = ((y2 + y1) - (x2 + x1))/2) 
  # %>%
  #   view()
}