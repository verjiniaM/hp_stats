rm(list = ls())
# Load libraries 
library(lme4)
library(ggplot2)
library(ggbeeswarm)
library(data.table)
library(emmeans)
library(dplyr)
library(tidyverse)
library(openxlsx)
library(r2glmm)

setwd('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/code/R/human_slice/')
source("funcs_human_stats.R")

data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/'
area <- 'temporal' # or 'complete'
df_repatch <- fread(paste(data_dir, 'repatch_data_', area,'.csv', sep = ''), sep = ",")
df_slice <- fread(paste(data_dir, 'slice_data_', area,'.csv', sep = ''), sep = ",")

distribution_dict <- list(
  'TH' = 'gamma',
  'rheo_ramp_c' = 'gaussian',
  'sag' = 'gaussian',
  'membra_time_constant_tau' = 'gaussian',
  'resting_potential' = 'gamma',
  'Rin' = 'gamma')

short_anova_model_comp_df <- data.frame()
ext_anova_model_comp_df <- data.frame()
emm_df <- data.frame()
emm_CI_df <- data.frame()
R_sq_short <- data.frame()
R_sq_ext <- data.frame()

df_list <- list(cell_ID_new = df_repatch, slice = df_slice)

for (data_type in names(df_list)) {
  
  df <- df_list[[data_type]]
  df[[data_type]] <- factor(df[[data_type]])
  df <- fix_df_slice_repatch(df)
  
  for (var in names(distribution_dict)) {
    print(paste('starting analysis of ', data_type, var))
    # saving org var; scaling var for better convergence
    var_org <- df[[var]]
    attrs <- attributes(scale(df[[var]]))
    df[[var]] <- as.vector(scale(df[[var]]))
    
    short_anova_model_comp_df <- get_results_anova_model_comparison(lmers(df, var, 'day', data_type),
                                                             short_anova_model_comp_df, data_type, var, 'day')
    
    ext_anova_model_comp_df <- get_results_anova_model_comparison(lmers_extended_2(df, var, 'hrs_after_OP', data_type),
                                                                     ext_anova_model_comp_df, data_type, var, 'hrs_after_OP')
    
    EMMS_short <- marg_effects_intr_short(df, data_type, var, var_org, attrs, emm_df, emm_CI_df)
    emm_df <- EMMS_short$df_emmeans
    emm_CI_df <- EMMS_short$df_CIs
    
    R_sq_short <- get_results_df_r2(df, R_sq_short, data_type, var, 'intr short')
    R_sq_ext <- get_results_df_r2(df, R_sq_ext, data_type, var, 'intr ext')
    
  }
}

# getting dfs for CIs
fig4_complete <- data.frame()
for (data_type in names(df_list)) {
  df <- df_list[[data_type]]
  df[[data_type]] <- factor(df[[data_type]])
  df <- df_ext_comparison_groups(df)
  for (var in names(distribution_dict)) {
    df$val <- df[[var]]
    df <- df[!is.na(df[[var]]), ]
    hrs_df <- as.data.frame(df %>%
                              group_by(treatment_word, day, hrs_group) %>%
                              summarize(mean = mean(val),
                                        median = median(val),
                                        SE = std.error(val),
                                        CI.L = confint(lm(val ~ 1), level=0.95)[1,1],
                                        CI.U = confint(lm(val ~ 1), level=0.95)[1,2]))
    hrs_df$group <- hrs_df$hrs_group
    hrs_df$hrs_group <- NULL
    hrs_df$data_type <- rep(data_type, nrow(hrs_df))
    hrs_df$param <- rep('hrs_after_OP', nrow(hrs_df))
    hrs_df$var <- rep(var, nrow(hrs_df))

    age_df <- as.data.frame(df %>%
                              group_by(treatment_word, day, age_group) %>%
                              summarize(mean = mean(val),
                                        median = median(val),
                                        SE = std.error(val),
                                        CI.L  = confint(lm(val ~ 1), level=0.95)[1,1],
                                        CI.U = confint(lm(val ~ 1), level=0.95)[1,2]))
    age_df$group <- age_df$age_group
    age_df$age_group <- NULL
    age_df$data_type <- rep(data_type, nrow(age_df))
    age_df$param <- rep('patient_age', nrow(age_df))
    age_df$var <- rep(var, nrow(age_df))

    fig4_complete <- rbind(fig4_complete, hrs_df, age_df)
  }
}

save_dir <- '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/stats/'
write.xlsx(fig4_complete, paste(save_dir, 'sum_data_age_hrs_', area,'.xlsx',sep = ''))

wb <- createWorkbook()
addWorksheet(wb, "short_model_")
writeData(wb, "short_model_", short_anova_model_comp_df)
addWorksheet(wb, "ext_model_")
writeData(wb, "ext_model_", ext_anova_model_comp_df)
addWorksheet(wb, "short_emmeans_")
writeData(wb, "short_emmeans_", emm_df)
addWorksheet(wb, "short_emmeans_CIs")
writeData(wb, "short_emmeans_CIs", emm_CI_df)
addWorksheet(wb, "short_expl_var")
writeData(wb, "short_expl_var", R_sq_short)
addWorksheet(wb, "ext_expl_var")
writeData(wb, "ext_expl_var", R_sq_ext)

saveWorkbook(wb, file = paste(save_dir, "intrinsic_complete_",area,".xlsx", sep = ''))

# OC plots ####
# Count unique slices per OP and treatment
df_s <- fix_df_slice_repatch(df_slice)

df_counts_slice <- df_s %>%
  group_by(OP, treatment) %>%
  summarise(unique_slice_count = n_distinct(slice), .groups = 'drop') %>%
  arrange(OP)

slice_per_op <- ggplot(df_counts_slice, aes(x = OP, y = unique_slice_count, fill = treatment)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Unique Slices per OP (stacked by treatment)", y = "Number of unique Slices") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/plots/QC_plots/slice/'
ggsave(filename = paste(save_dir, 'slice_repeat_slice_per_OP.png', sep = ''),
       plot = slice_per_op, # or replace with your plot object, e.g., p
       width = 8, height = 6, dpi = 300, bg = "white")

# VERY CRAZY PLOT ####
# # Create a dataframe with one row per OP, treatment, and slice
# df_expanded <- df_s %>%
#   distinct(OP, treatment, slice)
# 
# # Add a count column for stacking (all 1s, since each row is a unique slice)
# df_expanded$count <- 1
# 
# # Combine treatment and slice for unique fill color
# df_expanded$fill_group <- interaction(df_expanded$treatment, df_expanded$slice, drop = TRUE)
# 
# 
# # Suppose your treatments are "Ctrl" and "high K"
# treatments <- unique(df_expanded$treatment)
# 
# # Get unique slices for each treatment
# slices_ctrl <- unique(df_expanded$slice[df_expanded$treatment == "Ctrl"])
# slices_hik  <- unique(df_expanded$slice[df_expanded$treatment == "high K"])
# 
# # Create color palettes
# greens  <- colorRampPalette(brewer.pal(9, "Greens"))(length(slices_ctrl))
# purples <- colorRampPalette(brewer.pal(9, "Purples"))(length(slices_hik))
# 
# # Create a named vector for fill_group
# fill_colors <- c(
#   setNames(greens,  paste("Ctrl", slices_ctrl, sep = ".")),
#   setNames(purples, paste("high K", slices_hik, sep = "."))
# )
# 
# # Make sure fill_group is constructed as interaction(treatment, slice, sep = ".")
# df_expanded$fill_group <- interaction(df_expanded$treatment, df_expanded$slice, sep = ".")
# 
# colorslice_per_op <- ggplot(df_expanded, aes(x = OP, y = count, fill = fill_group)) +
#   geom_bar(stat = "identity", position = "stack") +
#   scale_fill_manual(values = fill_colors) +
#   labs(title = "Unique Slices per OP (stacked by treatment and slice)", y = "Number of unique Slices") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# print(colorslice_per_op)
# 
