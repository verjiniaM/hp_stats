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

df_slice_inc_only <- fread(paste(data_dir, '10.07.25_slice_incubation_only_no_move.csv', sep = ''), sep = ",")
area <- 'temporal_inc_only' # or 'complete'
min_hrs_incubation = 15.9
max_hrs_after_op = 35
data_type <- paste('inc time ', min_hrs_incubation, max_hrs_after_op, 'hrs')
df_slice_inc_only_f <- df_slice_inc_only[df_slice_inc_only$hrs_incubation > min_hrs_incubation, ]
df_slice_inc_only_f <- df_slice_inc_only_f[df_slice_inc_only_f$hrs_after_OP < max_hrs_after_op, ]
# removing slices which were recorded on D1. Keeping really only incubation
df_slice_inc_only_f <- df_slice_inc_only_f[nchar(df_slice_inc_only_f$slice) <= 2, ]

### INTRINSIC PROPERTIES ####
# distribution_dict <- list(
#   'TH' = 'gamma',
#   # 'rheo_ramp_c' = 'gaussian',
#   'sag' = 'gaussian',
#   'membra_time_constant_tau' = 'gaussian',
#   'resting_potential' = 'gamma',
#   'Rin' = 'gamma')
# 
# distribution_dict <- list(
#   'AP_heigth' = 'gamma',
#   'max_repol' = 'gaussian',
#   'max_depol' = 'gaussian',
#   'AP_halfwidth' = 'gamma')

distribution_dict <- list(
  'membra_time_constant_tau' = 'gamma',
  'resting_potential' = 'gaussian',
  'cap_adj' = 'gaussian')

short_anova_model_comp_df <- data.frame()
ext_anova_model_comp_df <- data.frame()
emm_df <- data.frame()
emm_CI_df <- data.frame()
emm_df_ext <- data.frame()
emm_CI_df_ext <- data.frame()
R_sq_short <- data.frame()
R_sq_ext <- data.frame()

df_slice_inc_only_f[['slice']] <- factor(df_slice_inc_only_f[['slice']])
# scaling of params; keeping original values
df <- fix_df(df_slice_inc_only_f)

for (var in names(distribution_dict)) {
  print(paste('starting analysis of ', var))
  
  # saving org var; scaling var for better convergence
  var_org <- df[[var]]
  attrs <- attributes(scale(df[[var]]))
  df[[var]] <- as.vector(scale(df[[var]]))
  
  short_anova_model_comp_df <- get_results_anova_model_comparison(lmers_slice_inc_only(df, var),
                                                                  short_anova_model_comp_df, data_type, var, 'day')
  
  ext_anova_model_comp_df <- get_results_anova_model_comparison(lmers_ext_inc_only(df, var, 'hrs_after_OP'),
                                                                ext_anova_model_comp_df, data_type, var, 'hrs_after_OP')
  
  EMMS_short <- marg_effects_intr_short_inc_only(df, '', var, var_org, attrs, emm_df, emm_CI_df)
  emm_df <- EMMS_short$df_emmeans
  emm_CI_df <- EMMS_short$df_CIs
  
  EMMS_ext <- marg_effects_intr_ext_inc_only(df, 'slice', var, var_org, attrs, emm_df_ext, emm_CI_df_ext)
  emm_df_ext <- EMMS_ext$df_emmeans
  emm_CI_df_ext <- EMMS_ext$df_CIs
  
  R_sq_short <- get_results_df_r2(df, R_sq_short, 'slice', var, 'intr_inc_only')
  R_sq_ext <- get_results_df_r2(df, R_sq_ext, 'slice', var, 'intr_ext_inc_only')
  
}

### FIRING ####

var_of_interest <- c('num_aps', 'IFF')
# define ID variables
ids <- c('OP', 'treatment_r', 'hrs_after_OP', 'slice')

short_model_results <- data.frame()
ext_model_results <- data.frame()
R_sq_short <- data.frame()
R_sq_ext <- data.frame()
emm_df <- data.frame()
emm_CI_df <- data.frame()

df_slice_inc_only_f[['slice']] <- factor(df_slice_inc_only_f[['slice']])
df <- fix_df(df_slice_inc_only_f)

for (var_firing in var_of_interest) {
  print(paste('starting analysis of ', var_firing))
  melt_df <- get_melt_df_firing(var_firing, df, ids)
  
  # scaled df for model stability
  df_scaled <- melt_df 
  df_scaled$VAL <- as.vector(scale(df_scaled$VAL))
  df_scaled$inj_current <- as.vector(scale(df_scaled$inj_current))
  
  short_model_results <- get_results_anova_model_comparison(lmers_firing_inc_only(df_scaled, var_firing, 'day', 'slice'),
                                                            short_model_results, data_type, var_firing, 'day')
  ext_model_results <- get_results_anova_model_comparison(lmers_firing_extended_inc_only(df_scaled, var_firing, 'slice'),
                                                          ext_model_results, data_type, var_firing, 'hrs_after_OP')
  
  EMMS_short <- marg_effects_firing_short_inc_only(df_scaled, 'slice', var_firing, melt_df$VAL, 
                                          attributes(scale(melt_df$VAL)), emm_df, emm_CI_df)
  emm_CI_df <- EMMS_short$df_CIs
  
  EMMS_short <- marg_effects_firing_short_inc_only(melt_df, 'slice', var_firing, melt_df$VAL, 
                                          attributes(scale(melt_df$VAL)), emm_df, emm_CI_df)
  emm_df <- EMMS_short$df_emmeans
  
  R_sq_short <- get_results_df_r2(df_scaled, R_sq_short, 'slice', var_firing, 'firing_short_inc_only')
  R_sq_ext <- get_results_df_r2(df_scaled, R_sq_ext, 'slice', var_firing, 'firing_ext_inc_only')
  
}

### SAVING STUFF ####
# save_dir <- '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/stats/inc_only/'
save_dir <- '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/figures/all_params_reference/add/'

# save the comparison table
# write.xlsx(anova_comp_df, paste(save_dir, 'model_com_slice_all_', area,'.xlsx',sep = ''))

# # INTRINSIC Add data frames to different sheets
wb <- createWorkbook()
# addWorksheet(wb, "cutoff_comparison")
# writeData(wb, "cutoff_comparison", anova_comp_df)
addWorksheet(wb, "intrinsic_short_model")
writeData(wb, "intrinsic_short_model", short_anova_model_comp_df)
addWorksheet(wb, "intrinsic_ext_only_hrs")
writeData(wb, "intrinsic_ext_only_hrs", ext_anova_model_comp_df)
addWorksheet(wb, "intrinsic_emmeans_short")
writeData(wb, "intrinsic_emmeans_short", emm_df)
addWorksheet(wb, "intrinsic_CI_short")
writeData(wb, "intrinsic_CI_short", emm_CI_df)
addWorksheet(wb, "intrinsic_emmeans_ext")
writeData(wb, "intrinsic_emmeans_ext", emm_df_ext)
addWorksheet(wb, "intrinsic_CI_ext")
writeData(wb, "intrinsic_CI_ext", emm_CI_df_ext)
addWorksheet(wb, "intrinsic_expl_var_short")
writeData(wb, "intrinsic_expl_var_short", R_sq_short)
addWorksheet(wb, "intrinsic_expl_var_ext")
writeData(wb, "intrinsic_expl_var_ext", R_sq_ext)

saveWorkbook(wb, file = paste(save_dir, "inc_only_intr_",max_hrs_after_op,".xlsx", sep = ''))

# FIRING
wb_fire <- createWorkbook()
addWorksheet(wb_fire, "firing_short_model")
writeData(wb_fire, "firing_short_model", short_model_results)
addWorksheet(wb_fire, "firing_ext_model")
writeData(wb_fire, "firing_ext_model", ext_model_results)
addWorksheet(wb_fire, "firing_short_emmeans")
writeData(wb_fire, "firing_short_emmeans", emm_df)
addWorksheet(wb_fire, "firing_short_CI")
writeData(wb_fire, "firing_short_CI", emm_CI_df)
addWorksheet(wb_fire, "exlp_var_short")
writeData(wb_fire, "exlp_var_short", R_sq_short)
addWorksheet(wb_fire, "exlp_var_ext")
writeData(wb_fire, "exlp_var_ext", R_sq_ext)
# saveWorkbook(wb_fire, file = paste(save_dir, "inc_only_firing_to_",max_hrs_after_op,".xlsx", sep = ''))


# ### comparisons to find last hrs incubation no change ####
# area <- 'temporal' # or 'complete'
# df_slice_all <- fread(paste(data_dir, 'slice_all.csv', sep = ''), sep = ",")
# 
# hist(df_slice_all$hrs_after_OP, 
#      main = "Histogram of Hours After Operation", 
#      xlab = "Hours After Operation", 
#      ylab = "Frequency", 
#      col = "blue", 
#      border = "black")
# 
# var <- 'TH'
# id2 <- 'hrs_after_OP'
# data_type <- 'slice'
# 
# df_slice_all[[data_type]] <- factor(df_slice_all[[data_type]])
# # scaling of params; keeping original values
# df_slice_all <- fix_df(df_slice_all)
# 
# anova_comp_df <- data.frame()
# for (hr in c(20, 25, 30, 35, 40, 45, 50)){
#   
#   df <- df_slice_all[df_slice_all$org_hrs_after_OP < hr, ]
#   
#   var_org <- df[[var]]
#   attrs <- attributes(scale(df[[var]]))
#   df[[var]] <- as.vector(scale(df[[var]]))
#   n_row_all <- nrow(df)
#   
#   #remove rows with nan values
#   df <- df[!is.na(df[[var]]), ]
#   del_nans <- n_row_all - nrow(df)
#   if (del_nans != 0){
#     paste('deleting', del_nans, 'number of rows with NANs')
#   }
#   
#   formula_full <- as.formula(paste(var, "~ treatment_r  +", id2, " + (1|OP) + (1|OP:", data_type, ")"))
#   formula_hrs_0 <- as.formula(paste(var, "~ treatment_r   + (1|OP) + (1|OP:", data_type, ")"))
#   
#   FullModel <- lmer(formula_full, data = df, REML = FALSE)
#   IncNullModel <- lmer(formula_hrs_0, data = df, REML = FALSE)
#   anova_full_inc_0 <- anova(FullModel, IncNullModel)
#   list(full_vs_inc_null = anova_full_inc_0)
#   
#   dt <- paste('slice under', hr, 'hrs incubation')
#   anova_comp_df <- get_results_anova_model_comparison(list(full_vs_inc_null = anova_full_inc_0), anova_comp_df, dt, var, id2)
# }
# 
# 
# ### QC plots ####
# # define df

df_counts <- df %>%
  group_by(OP, treatment) %>%
  summarise(unique_slice_count = n_distinct(slice), .groups = 'drop')

slice_per_op <- ggplot(df_counts, aes(x = OP, y = unique_slice_count, fill = treatment)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Unique Slices per OP (stacked by treatment)", y = "Number of unique Slices") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/plots/QC_plots/slice/'
ggsave(filename = paste(save_dir, 'inc_only_slices_per_op.png', sep = ''),
  plot = slice_per_op, # or replace with your plot object, e.g., p
  width = 8, height = 6, dpi = 300)

# ggplot(df, aes(x = OP, y = 'TH', fill = slice)) +
#   geom_bar(stat = "identity") +
#   labs(title = "Barplot of OP by Slice", y = "Value") +
#   theme_minimal()
# 
# # plotting IV curves from melt_df
# 
# # Calculate averages for inj_current and VAL grouped by treatment_r
# df_summary <- melt_df %>%
#   group_by(treatment_r, inj_current) %>%
#   summarise(mean_VAL = mean(VAL, na.rm = TRUE))
# 
# # Create the plot
# ggplot(melt_df, aes(x = inj_current, y = VAL, color = treatment_r)) +
#   geom_point(alpha = 0.5)+  # Scatter plot for individual points
#   geom_point(data = df_summary, aes(x = inj_current, y = mean_VAL), size = 3, shape = 21, fill = "white") +  # Dots for averages
#   geom_line(data = df_summary, aes(x = inj_current, y = mean_VAL, group = treatment_r)) +  # Line through averages
#   labs(title = "VAL vs inj_current by Treatment",
#        x = "Injection Current",
#        y = "VAL",
#        color = "Treatment") +
#   theme_minimal()
# 
