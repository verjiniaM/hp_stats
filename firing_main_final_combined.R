rm(list = ls())
# Load libraries 
library(lme4)
library(data.table)
library(emmeans)
library(tidyverse)
library(openxlsx)
library(r2glmm)

setwd('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/code/R/human_slice/')
source("funcs_human_stats.R")

data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/'
area <- 'temporal'
df_repatch <- fread(paste(data_dir, 'repatch_data_', area,'.csv', sep = ''), sep = ",")
df_slice <- fread(paste(data_dir, 'slice_data_', area,'.csv', sep = ''), sep = ",")

var_of_interest <- c('num_aps', 'IFF')
# define ID variables
ids_basic <- c('OP', 'day', 'treatment_r', 'hrs_after_OP', 'patient_age')
ids_ext_b <- c('OP', 'day', 'treatment', 'hrs_group', 'age_group')

short_model_results <- data.frame()
ext_model_results <- data.frame()
R_sq_short <- data.frame()
R_sq_ext <- data.frame()
emm_df <- data.frame()
emm_CI_df <- data.frame()
ext_emm_df <- data.frame()
ext_emm_CI_df <- data.frame()

df_list <- list(cell_ID_new = df_repatch, slice = df_slice)
for (data_type in names(df_list)) {

  ids <- append(ids_basic, data_type, after = length(ids_basic))
  ids_ext <- append(ids_ext_b, data_type, after = length(ids_ext_b))
  
  df <- df_list[[data_type]]
  df[[data_type]] <- factor(df[[data_type]])
  df_ext <- df_ext_comparison_groups(df)
  df <- fix_df_slice_repatch(df)
  
  for (var_firing in var_of_interest) {
    
    melt_df <- get_melt_df_firing(var_firing, df, ids)
    
    # scaled df for model stability
    df_scaled <- melt_df 
    df_scaled$VAL <- as.vector(scale(df_scaled$VAL))
    df_scaled$inj_current <- as.vector(scale(df_scaled$inj_current))
    
    short_model_results <- get_results_anova_model_comparison(lmers_firing(df_scaled, var_firing, 'day', data_type),
                                                              short_model_results, data_type, var_firing, 'day')
    ext_model_results <- get_results_anova_model_comparison(lmers_firing_extended_no_age(df_scaled, var_firing, data_type),
                                                            ext_model_results, data_type, var_firing, 'hrs_after_OP')
    
    EMMS_short <- marg_effects_firing_short(df_scaled, data_type, var_firing, melt_df$VAL, 
                                            attributes(scale(melt_df$VAL)), emm_df, emm_CI_df)
    
    emm_CI_df <- EMMS_short$df_CIs
    
    EMMS_short <- marg_effects_firing_short(melt_df, data_type, var_firing, melt_df$VAL, 
                                            attributes(scale(melt_df$VAL)), emm_df, emm_CI_df)
    emm_df <- EMMS_short$df_emmeans
    
    # EMMS_short_change <-
    
    melt_df_ext <- get_melt_df_firing(var_firing, df_ext, ids_ext)
    EMMS_ext <- marg_effects_firing_ext(melt_df_ext, var_firing, data_type, ext_emm_df, ext_emm_CI_df)
    ext_emm_df <- EMMS_ext$df_emmeans
    ext_emm_CI_df <- EMMS_ext$df_CIs
    
    R_sq_short <- get_results_df_r2(df_scaled, R_sq_short, data_type, var_firing, 'firing short')
    R_sq_ext <- get_results_df_r2(df_scaled, R_sq_ext, data_type, var_firing, 'firing ext')

  }
  rm(ids, ids_ext)
}

save_dir <- '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/stats/'

wb <- createWorkbook()
addWorksheet(wb, "short_model_")
writeData(wb, "short_model_", short_model_results)
addWorksheet(wb, "ext_model_")
writeData(wb, "ext_model_", ext_model_results)
addWorksheet(wb, "short_emmeans_")
writeData(wb, "short_emmeans_", emm_df)
addWorksheet(wb, "short_emmeans_CIs")
writeData(wb, "short_emmeans_CIs", emm_CI_df)
addWorksheet(wb, "ext_emmeans")
writeData(wb, "ext_emmeans", ext_emm_df)
addWorksheet(wb, "ext_emmeans_CIs")
writeData(wb, "ext_emmeans_CIs", ext_emm_CI_df)
addWorksheet(wb, "short_expl_var")
writeData(wb, "short_expl_var", R_sq_short)
addWorksheet(wb, "ext_expl_var")
writeData(wb, "ext_expl_var", R_sq_ext)

saveWorkbook(wb, file = paste(save_dir, "firing_complete_",area,".xlsx", sep = ''))
