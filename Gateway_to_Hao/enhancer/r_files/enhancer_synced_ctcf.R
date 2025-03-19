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
setwd('/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files')
getwd()
source(file.path('/Users/PanjunKim/dropbox/Gateway_to_Hao/project_common_code/', 'variables.R'))
source(file.path('/Users/PanjunKim/dropbox/Gateway_to_Hao/project_common_code/', 'funcs.R'))

# source(file.path('/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/enhancer_synced_data_preparation.R'))
# source(file.path('/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/data_analysis.R'))
# source(file.path('/Users/PanjunKim/dropbox/Gateway_to_Hao/enhancer/r_files/data_analysis.R'))
