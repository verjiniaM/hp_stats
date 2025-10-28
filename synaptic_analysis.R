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
save_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/stats/events/'
df_single_events <- fread(paste(data_dir, '19.10.2025_single_events_QCed.csv', sep = ''), sep = ",")

event_type = 'spontan'

distribution_dict <- list(
  'IEI_Hz' = 'g',
  'abs_amp_barb_pA' = 'g')

short_anova_model_comp_df <- data.frame()
emm_df <- data.frame()
emm_CI_df <- data.frame()
R_sq_short <- data.frame()

df_single_events$source <- factor(df_single_events$source)
df_types = list('incubation' = '',
               'slice' = 'day')

for (data_type in names(df_types)) {
  df <- df_single_events[df_single_events$source ==  data_type, ]
  df <-  fix_df(df)
  
  # df = df %>%
  #   mutate(
  #      = scale(IEI_Hz),
  #      = scale(abs_amp_barb_pA)
  #   )
  
  for (var in names(distribution_dict)){
    id2 <- df_types[data_type]
    var_org <- df[[var]]
    attrs <- attributes(scale(df[[var]]))
    df[[var]] <- as.vector(scale(df[[var]]))
    
    short_anova_model_comp_df <- get_results_anova_model_comparison(lmers_synaptic(df, var, id2),
                                                                    short_anova_model_comp_df,
                                                                    paste(event_type, data_type, sep = '_'),
                                                                    var, id2)
    
    EMMS_short <- marg_effects_spontan(df, paste(event_type, data_type, sep = '_'), id2,
                                          var, var_org, attrs, emm_df, emm_CI_df)
    emm_df <- EMMS_short$df_emmeans
    emm_CI_df <- EMMS_short$df_CIs
    
    R_sq_short <- get_results_df_r2(df, R_sq_short, paste(event_type, data_type, sep = '_'), var, 'events')

  }
}

# # INTRINSIC Add data frames to different sheets
wb <- createWorkbook()
addWorksheet(wb, "short_model")
writeData(wb, "short_model", short_anova_model_comp_df)
addWorksheet(wb, "intrinsic_expl_var_short")
writeData(wb, "intrinsic_expl_var_short", R_sq_short)
addWorksheet(wb, "emmeans_short")
writeData(wb, "emmeans_short", emm_df)
addWorksheet(wb, "CI_short")
writeData(wb, "CI_short", emm_CI_df)

saveWorkbook(wb, file = paste(save_dir, event_type, "model_comparison.xlsx", sep = ''))


