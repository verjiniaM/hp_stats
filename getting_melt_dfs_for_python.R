rm(list = ls())
# Load libraries 
library(lme4)
library(data.table)
library(emmeans)
library(tidyverse)
library(openxlsx)
library(r2glmm)

setwd('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/code/R/human_slice/hp_stats/')
source("funcs_human_stats.R")

data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/'
area <- 'temporal'
df_repatch <- fread(paste(data_dir, 'repatch_data_', area,'.csv', sep = ''), sep = ",")
df_slice <- fread(paste(data_dir, 'slice_data_', area,'.csv', sep = ''), sep = ",")

df_slice_inc_only <- fread(paste(data_dir, '10.07.25_slice_incubation_only_no_move.csv', sep = ''), sep = ",")
area <- 'temporal_inc_only' # or 'complete'
min_hrs_incubation = 15.9
max_hrs_after_op = 35
data_type <- paste('inc time ', min_hrs_incubation, max_hrs_after_op, 'hrs')
df_inc <- df_slice_inc_only[df_slice_inc_only$hrs_incubation > min_hrs_incubation, ]
df_inc <- df_inc[df_inc$hrs_after_OP < max_hrs_after_op, ]
# removing slices which were recorded on D1. Keeping really only incubation
df_inc <- df_inc[nchar(df_inc$slice) <= 2, ]
df_inc$slice_inc <- df_inc$slice

var_of_interest <- c('num_aps', 'IFF')
# define ID variables
ids_basic <- c('OP', 'day', 'treatment_r', 'hrs_after_OP')

df_list <- list(cell_ID_new = df_repatch, slice = df_slice, slice_inc = df_inc)
for (data_type in names(df_list)) {
  
  ids <- append(ids_basic, data_type, after = length(ids_basic))
  print(ids)

  df <- df_list[[data_type]]
  # df[[data_type]] <- factor(df[[data_type]])
  # df <- fix_df_slice_repatch(df)
  
  for (var_firing in var_of_interest) {
    melt_df <- get_melt_df_firing_all_inj(var_firing, df, ids)
  }}