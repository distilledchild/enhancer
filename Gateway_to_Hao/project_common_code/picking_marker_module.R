library(httr)
library(jsonlite)
library(xml2)
library(tidyverse) 
library(fs)
library(matrixStats)
library(RcppAlgos)
library(tidytidbits)
library(tidylog)
library(tidyselect)
library(pryr)
library(RevoUtilsMath)
library(parallel)

options(scipen = 100)

# retrieving #core 
avail_thread_n <- detectCores(TRUE) %>% print()
print(avail_thread_n)

SetMKL <-function(nCPU) {
	if(require("RevoUtilsMath")){
        setMKLthreads(nCPU)
    } else {
        print("Please install Microsoft R Open")
    }
}

### Check correctness with getMKLthreads()
GetMKL <-function() {

	if(require("RevoUtilsMath")){
        return(getMKLthreads())
    } else {
        print("Please install Microsoft R Open")
    }
}

getwd()

source(file.path('C:/Users/panju/Dropbox (UTHSC GGI)/WLIWMIB6Research', '/project_common_code/variables.R'))
source(file.path(root.dir, '/project_common_code/funcs.R'))

wdpath <- file.path(root.dir, oxy.project.ref.dir)
setwd(wdpath)

# listing up files (two ANN file of WLI and WMI)
file.dir <- path(file.path(root.dir, oxy.project.n.dir, "pipeline_output\\output\\bcftools_mpileup_output_vcf"))
file.dir
file.list = fs::dir_ls(file.dir, regexp = "\\.vcf$")
file.list

# needed for chromosomal end to end coordinate
rn6.ref <- read_tsv(file.path(root.dir, oxy.project.ref.dir, "rn6.txt")) %>% 
    rename_at(1, ~"Chr") %>%
    mutate(Start = 0)  %>% 
    rename_at(2, ~"End") %>% 
    select(Chr, Start, End) %>% 
    mutate_at(c(2:3), as.numeric) %>% 
    view()

########################################################################################################################################
distance_param <- 20000000
########################################################################################################################################
df.marker.picking.initial.data <- file.list %>% 
    map_dfr(vcf_raw_data_reader, .id = "sample") %>% 
    rename_at(vars(starts_with("#")), ~str_replace(., "#", "")) %>%
    mutate(sample = str_to_upper(str_split_n(str_split_n(sample, 'vcf\\/', 2), '\\.', 1))) %>%
    separate(sample, into=c("sample_file", "sstrain", "calls", "sq", "ann"), sep = '[|_]') %>%          
    filter(str_length(ALT) == 1 & str_length(REF) == 1) %>% 
    mutate(ANN = str_split_n(str_split_n(INFO, ";", 3), "\\|", 2)) %>% 
    mutate(WLI_IonProton_GT = str_split_n(WLI_IonProton, ":", 1)) %>% 
    mutate(WLI_chromium_GT = str_split_n(WLI_chromium, ":", 1)) %>% 
    mutate(WLI_xTen_GT = str_split_n(WLI_xTen, ":", 1)) %>% 
    mutate(WMI_IonProton_GT = str_split_n(WMI_IonProton, ":", 1)) %>% 
    mutate(WMI_chromium_GT = str_split_n(WMI_chromium, ":", 1)) %>% 
    mutate(WMI_xTen_GT = str_split_n(WMI_xTen, ":", 1)) %>% 
    mutate(WLI_IonProton_DP = as.numeric(str_split_n(WLI_IonProton, ":", 3))) %>% 
    mutate(WLI_chromium_DP = as.numeric(str_split_n(WLI_chromium, ":", 3))) %>% 
    mutate(WLI_xTen_DP = as.numeric(str_split_n(WLI_xTen, ":", 3))) %>% 
    mutate(WMI_IonProton_DP = as.numeric(str_split_n(WMI_IonProton, ":", 3))) %>% 
    mutate(WMI_chromium_DP = as.numeric(str_split_n(WMI_chromium, ":", 3))) %>% 
    mutate(WMI_xTen_DP = as.numeric(str_split_n(WMI_xTen, ":", 3))) %>% 
    select(-c(calls, ann, FILTER, INFO:WMI_xTen, WLI_IonProton_GT:WMI_xTen_DP)) %>% 
    arrange(as.integer(CHROM), as.integer(POS)) %>% 
    unite("ID", c(sstrain, ID, QUAL), remove = FALSE) %>% # new ID column
    view()

# non-intergenic 636 from WLI and WMI
ongoing.df <- df.marker.picking.initial.data %>% 
    # 모든 filter 선행에 두기
    filter(ANN != "intergenic_region") %>% 
    # 모든 filter 선행에 두기
    filter(sstrain == "WLI") %>% 
    group_by(CHROM) %>% 
    mutate(MSTART = min(POS), MEND = max(POS)) %>% # prefix: M : marker
    ungroup() %>% 
    left_join(rn6.ref %>% rename(REND = End) %>% rename(RSTART = Start) %>%  mutate(RSTART = RSTART + 1), by=c("CHROM"="Chr")) %>%  # 149 markers
    # pmap_dfr(marker_coord_setter) %>% 
    select(-c(sample_file, REF, ALT)) %>% 
    mutate(distance = POS - lag(POS)) %>%
    mutate(filter_distance = if_else(lag(CHROM) != CHROM | is.na(lag(CHROM)), 0, distance - distance_param)) %>% 
    mutate(cnddt_mrkr_grp = cumsum(filter_distance >= 0 )) %>% 
    group_by(cnddt_mrkr_grp) %>% 
    add_tally() %>% 
    mutate(max_pssbl_pick = if_else(n == 1, 1, ceiling((max(POS) - min(POS))/distance_param))) %>% 
    ungroup() %>% 
    mutate(cnddt_marker_set = if_else(n == 1 & max_pssbl_pick == 1, 'Y', 'N')) %>% 
    # filter(cnddt_mrkr_grp == 2 | cnddt_mrkr_grp == 11) %>% 
    filter(n != 1) %>% 
    # filter(max_pssbl_pick < 8) %>% 
    # filter(between(cnddt_mrkr_grp, 1, 3)) %>% 
    view()
    # dplyr::count(cnddt_mrkr_grp) %>% # 125 groups WLI
    # filter(n==1) %>% # 20 groups

###########################################################################
###########################################################################
distance_param = 10000000
distance_param = 900000
distance_param = 400000
distance_param = 7000000
distance_param = 15000000
distance_param
###########################################################################
priority_checker <- function(x) {
    
    print("priority_checker function start ======")

    # TODO 3 connecting this picking marker job to previous picked marker work
    # maximum range priority
    onlyone <- x %>% 
        filter(decision == 1) %>% 
        mutate(maxrange = rowSums(.[1:(match('decision', names(.)))])) %>% 
        view()
        # dplyr::slice_max(maxrange) %>% print(.)

    print("priority_checker function end ==========")

    # same number of strains
    # variant type
    # qual
}
distance_param
spacing_checker <- function(df) {
    # collect garbage
    gc()
    # seedlings
    set.seed (1)

    # if(df$n == 1) {
    #     df %>% mutate(picked_yn = 'Y') 
    #     return
    # }
    
    mxpick <- as.numeric(max(df$max_pssbl_pick))
    show <- NULL

    
    # TODO 6 data structure 바꾸기
    for (trial in mxpick:1) {
        print(is.atomic(df$ID))
        print(class(df$ID))
        print(dim(df)) # 22, 16 // 16, 16
        print(df) # 22, 16 // 16, 16
        print((trial)) 

        rm(show)
        
        # row ID value and qual_mean
        sho <- df$ID %>% 
            comboGeneral(m = trial, Parallel = TRUE, nThreads = avail_thread_n) %>% 
            as_tibble() %>% # for column names
            bind_cols(apply(., MARGIN=2, FUN=function(x){as.integer(str_split_n(x, '\\_', 6))})) %>%
            mutate(qual_mean = as.character(rowMeans(across(where(is.integer))))) %>% # mutate: new variable 'qual_mean' (character) with 178 unique values and 0% NA
            select_if(is.character) 
        
        # TODO 4 large n and r speed enhancement
        show <- sho %>% 
            select(-qual_mean) %>% # select: dropped one variable (qual_mean)
            apply(., 2, FUN=function(x){as.numeric(str_split_n(x, '\\_', 3))}) %>%
            {if(ncol(.) > 1) rowDiffs(.) else .} %>% 
            as_tibble() %>% 
            mutate(distance_range = rowSums(.[1:ncol(.)])) %>% # memory problem
            rowwise() %>% mutate(decision = if_else(trial == 1, 1, if_else(any(c_across(starts_with('V')) < distance_param), 0, 1))) %>% # memory problem
            bind_cols(sho %>% mutate(qual_mean = as.numeric(qual_mean))) 

        # iris %>%
        # slice(1:6) %>%
        # select(starts_with('Petal')) %>% 
        # rowwise() %>%
        # mutate(sum = sum(cur_data())) %>%
        # ungroup

        rm(sho)
        ls()

        # ifelse(show$decision == 0, rm(show), break)
        ifelse(show$decision == 0, print("keep looping"), break)
    }

    picked_set <- show %>% 
        filter(decision == 1) %>% 
        execute_if_else(mxpick == 1, arrange_at('qual_mean', desc), arrange(-distance_range, -qual_mean)) %>% 
        head(n = 1) %>% 
        select(-c(1:decision, qual_mean)) %>% 
        t() %>% 
        as_tibble() %>% 
        mutate(picked_yn = 'Y') %>% 
        rename(picked_markers = V1) %>% 
        right_join(df, by = c('picked_markers' = 'ID')) %>% 
        mutate(picked_yn = replace_na(picked_yn, 'N')) 

    rm(show)
    # TODO 1 additional priority
    # if(sum(show %>% pull(decision) == 1) > 1) show <- priority_checker(show) # only one case should be return
    return(picked_set)
}
###########################################################################

# marker.picked.df <- function() {
    first.output <- ongoing.df %>% 
        group_by(cnddt_mrkr_grp) %>%
        group_modify(~ spacing_checker(.x)) %>% 
        filter(picked_yn == 'Y') %>% 
        write_csv(file.path(root.dir, oxy.project.r.dir, "test1.csv"), append = FALSE, na = "NA", col_names = TRUE)
# } 

# TODO 2 FASTA file output
first.output %>% 
    filter(picked_yn == 'Y') %>% 
    ungroup() %>% 
    mutate(seq = mapply(marker_seq_fetcher, .$CHROM, .$POS - 150, .$POS + 150, species.rat)) %>% 
    # bind_cols(marker.picked.df %>% picked_yn == 'N' %>% mutate(seq = '')) %>% 
    # arrange(cnddt_mrkr_grp, POS) %>% 
    write_csv(file.path(root.dir, oxy.project.r.dir, "test2.csv"), append = FALSE, na = "NA", col_names = TRUE) %>% 
    view()

nth(marker.picked.df, 55) %>% view()

# TODO 5 ideogram work

file.path(root.dir, oxy.project.r.dir, "test1.csv")

SetMKL(avail_thread_n); cat(GetMKL(), "MKLthreads\n");
marker.picked.df()
