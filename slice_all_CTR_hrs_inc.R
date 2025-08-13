
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

df_slice_all <- fread(paste(data_dir, 'slice_all.csv', sep = ''), sep = ",")
# taking only 'Ctrl' treatment
df_slice_ctr <- df_slice_all %>%
  filter(day == "D1" | (day == "D2" & treatment == "Ctrl"))

# check for D1 and slice name >2
any(df_slice_ctr$day == "D1" & nchar(df_slice_ctr$slice) > 2)

# Plotting ####
df <- df_slice_ctr
df <- df %>%
  mutate(
    hrs_bin = cut(
      hrs_after_OP,
      breaks = c(5, 15, 25, 35, 45, 55),
      labels = c("5-15", "16-25", "26-35", "36-45", "46-55"),
      include.lowest = TRUE, right = TRUE
    )
  )

# Calculate mean and SE for each bin
summary_df <- df %>%
  group_by(hrs_bin) %>%
  summarise(
    mean_TH = mean(TH, na.rm = TRUE),
    se_TH = sd(TH, na.rm = TRUE) / sqrt(n()),
    n = n(),
    ci_lower = mean_TH - qt(0.975, n-1) * se_TH,
    ci_upper = mean_TH + qt(0.975, n-1) * se_TH
  )

# Plot with SE
p_se <- ggplot(df, aes(x = hrs_bin, y = TH)) +
  geom_jitter(width = 0.2, alpha = 0.6, color = "grey40") +
  geom_line(data = summary_df, aes(x = hrs_bin, y = mean_TH, group = 1), color = "orange") +
  geom_point(data = summary_df, aes(x = hrs_bin, y = mean_TH), color = "orange") +
  geom_errorbar(data = summary_df, aes(x = hrs_bin, ymin = mean_TH - se_TH, ymax = mean_TH + se_TH), 
                inherit.aes = FALSE, width = 0.2, color = "orange") +
  geom_ribbon(data = summary_df, aes(x = hrs_bin, ymin = mean_TH - se_TH, ymax = mean_TH + se_TH, group = 1),
              fill = "orange", alpha = 0.2, inherit.aes = FALSE)
labs(x = "Hours after OP (binned)", y = "TH", title = "TH vs hrs_after_OP (SE)") +
  theme_minimal()

p_ci <- ggplot(df, aes(x = hrs_bin, y = TH)) +
  geom_jitter(width = 0.2, alpha = 0.6, color = "grey40") +
  geom_line(data = summary_df, aes(x = hrs_bin, y = mean_TH, group = 1), color = "green") +
  geom_point(data = summary_df, aes(x = hrs_bin, y = mean_TH), color = "green") +
  geom_errorbar(data = summary_df, aes(x = hrs_bin, ymin = ci_lower, ymax = ci_upper), 
                inherit.aes = FALSE, width = 0.2, color = "green") +
  geom_ribbon(data = summary_df, aes(x = hrs_bin, ymin = ci_lower, ymax = ci_upper, group = 1),
              fill = "green", alpha = 0.2, inherit.aes = FALSE) +
  labs(x = "Hours after OP (binned)", y = "TH", title = "TH vs hrs_after_OP (95% CI)") +
  theme_minimal()

# analysis ####

distribution_dict <- list(
  'TH' = 'gamma',
  'rheo_ramp_c' = 'gaussian',
  'sag' = 'gaussian',
  'membra_time_constant_tau' = 'gaussian',
  'resting_potential' = 'gamma',
  'Rin' = 'gamma')

ext_anova_model_comp_df <- data.frame()
emm_df <- data.frame()
emm_CI_df <- data.frame()
R_sq_hrs <- data.frame()

df <- df_slice_ctr
df$slice <- factor(df$slice)
df <- fix_df_slice_repatch(df)

for (var in names(distribution_dict)) {
  print(paste('starting analysis of ', var))
  # saving org var; scaling var for better convergence
  var_org <- df[[var]]
  attrs <- attributes(scale(df[[var]]))
  df[[var]] <- as.vector(scale(df[[var]]))
  
  ext_anova_model_comp_df <- get_results_anova_model_comparison(lmers_extended_hrs_only(df, var, 'slice'),
                                                                ext_anova_model_comp_df, 'slice_all_CTR', var, 'hrs_after_OP')
  
  EMMS_short <- marg_effects_intr_hrs_only(df, var, var_org, attrs, emm_df, emm_CI_df)
  emm_df <- EMMS_short$df_emmeans
  emm_CI_df <- EMMS_short$df_CIs
  
}

save_dir <- '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/stats/'

wb <- createWorkbook()
addWorksheet(wb, "hrs_model_slice_CTR")
writeData(wb, "hrs_model_slice_CTR", ext_anova_model_comp_df)
addWorksheet(wb, "hrs_emmeans_")
writeData(wb, "hrs_emmeans_", emm_df)
addWorksheet(wb, "hrs_emmeans_CIs")
writeData(wb, "hrs_emmeans_CIs", emm_CI_df)

saveWorkbook(wb, file = paste(save_dir, "hrs_after_OP_slice_CTR.xlsx", sep = ''))

plot_save <- '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/figures/draft4_figs/from_r_for_comparison/'
ggsave(
  filename = paste(plot_save,  "SE_slice_all_CTR.png", sep = ''),   # Path and filename for the PNG
  plot = p_se,           # Or use your plot object, e.g., plot = p
  width = 8, height = 6, dpi = 300, # Adjust size and resolution as needed
  bg = "white"                  # Ensures a white background
)
ggsave(
  filename = paste(plot_save,  "CI_slice_all_CTR.png", sep = ''),   # Path and filename for the PNG
  plot = p_ci,           # Or use your plot object, e.g., plot = p
  width = 8, height = 6, dpi = 300, # Adjust size and resolution as needed
  bg = "white"                  # Ensures a white background
)



