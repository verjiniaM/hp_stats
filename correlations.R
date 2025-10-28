rm(list = ls())
# Load libraries 
library(lme4)
library(data.table)
library(emmeans)
library(tidyverse)
library(openxlsx)
library(r2glmm)
library('psych')
library(ggplot2)
library(gridExtra)


setwd('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/code/R/human_slice/hp_stats/')
source("funcs_human_stats.R")

analysis_type <- 'events' # 'events' or 'intrinsic' or 'intrinsic_time'
columns_of_i = c('hrs_after_OP',
                'Rin',
                'resting_potential',
                'rheo_ramp_c',
                'TH', 
                'AP_halfwidth',
                'max_depol',
                'max_repol',
                'source',
                'day',
                'treatment',
                'OP')

# separate by source
data_dir <- '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/'
save_dir <- paste(str_sub(data_dir, 1,89), 'stats/correlations/', sep = '')

# for events and intrinsic complete ####
if (analysis_type == 'events'){
  columns_of_i = c(columns_of_i, 'avg_amp_adj', 'avg_IEI_adj')
  
  df_single_events <- fread(paste(data_dir, '19.10.2025_single_events_QCed.csv', sep = ''), sep = ",")
  df_sum <- fread(paste(data_dir, '19.10.2025_events_summary_QCed.csv', sep = ''), sep = ",")
  df_sum  <- df_sum[, ..columns_of_i]
  df_list <- c('incubation', 'slice')
  
} else if (analysis_type == 'intrinsic'){
  df_repatch <- fread(paste(data_dir, 'repatch_data_temporal.csv', sep = ''), sep = ",")
  df_repatch$source <- 'intr_repatch'
  df_slice <- fread(paste(data_dir, 'slice_data_temporal.csv', sep = ''), sep = ",")
  df_slice$source <- 'intr_slice'
  df_slice_inc_only <- fread(paste(data_dir, '10.07.25_slice_incubation_only_no_move.csv', sep = ''), sep = ",")
  area <- 'temporal_inc_only' # or 'complete'
  min_hrs_incubation = 15.9
  max_hrs_after_op = 35
  data_type <- paste('inc time ', min_hrs_incubation, max_hrs_after_op, 'hrs')
  df_slice_inc_only_f <- df_slice_inc_only[df_slice_inc_only$hrs_incubation > min_hrs_incubation, ]
  df_slice_inc_only_f <- df_slice_inc_only_f[df_slice_inc_only_f$hrs_after_OP < max_hrs_after_op, ]
  # removing slices which were recorded on D1. Keeping really only incubation
  df_inc_only <- df_slice_inc_only_f[nchar(df_slice_inc_only_f$slice) <= 2, ]
  df_inc_only$source <- 'intr_incubation'
  
  df_list <- c('intr_incubation','intr_slice', 'intr_repatch')
  
  # create 1 summary df
  df_sum <- rbind(df_slice[, ..columns_of_i], 
                  df_repatch[, ..columns_of_i], 
                  df_inc_only[, ..columns_of_i])
}

# removing meta info columns
params_corr =  setdiff(columns_of_i, c('day', 'treatment', 'OP', 'source'))

conditions <- list('pre' = c('D1', ''),
                   'CTR' = c('D2','high K'), # using != for less loops
                   'HiK' = list('D2', 'Ctrl'))


for (data_type in df_list) {
  for (con in names(conditions)){
    
    if (str_detect(data_type, 'incubation') & con =='pre'){
      next
      }
    
      print(c(data_type, con))
      df <-  df_sum[df_sum$source == data_type &
                      df_sum$day ==  conditions[[con]][1] &
                      df_sum$treatment != conditions[[con]][2] ,]
      
      # calculate correlations
      # partial_corr_mixed(df, params_corr, '(1|OP)')
      per_cor <- pearson_corr(df, params_corr) # $ corrs, pvals, l_conf, u_conf
      
      plot_list <- list()
      plot_count <- 1
      
      for (i in 1:(length(params_corr)-1)) {
        for (j in (i+1):length(params_corr)) {
          param_x <- params_corr[i]
          param_y <- params_corr[j]
          
          # Get correlation values
          cor_val <- per_cor$corrs[param_x, param_y]
          p_val <- per_cor$pvals[param_x, param_y]
          lower_ci <- per_cor$l_conf[param_x, param_y]
          upper_ci <- per_cor$u_conf[param_x, param_y]
          
          # Create scatter plot with regression line and CI
          p <- ggplot(df, aes_string(x = param_x, y = param_y)) +
            geom_point(alpha = 0.6, size = 2) +
            geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
            labs(
              title = paste(param_x, "vs", param_y),
              subtitle = paste("r =", round(cor_val, 3), 
                               ", p =", round(p_val, 3),
                               "\nCI: [", round(lower_ci, 3), ",", round(upper_ci, 3), "]"),
              x = param_x,
              y = param_y
            ) +
            theme_minimal() +
            theme(plot.title = element_text(size = 10),
                  plot.subtitle = element_text(size = 8))
          
          plot_list[[plot_count]] <- p
          plot_count <- plot_count + 1
        }
      }
      
      n <-  length(params_corr)
      layout_matr <- get_layout_matr(n)
      layout_matr[layout_matr == 0] <- NA
      
      # Save combined plot
      combined_plot <- grid.arrange(grobs = plot_list, layout_matrix = layout_matr)
      
      # Save the plot
      plot_filename <- paste(save_dir,analysis_type, '/', data_type, "_", con, "_correlation.png", sep = '')
      ggsave(plot_filename, combined_plot, 
             width = n * 4, height = n * 3, dpi = 300)
      
      # saving the data
      wb <- createWorkbook()
      addWorksheet(wb, "cor_estimate")
      writeData(wb, "cor_estimate", per_cor$corrs)
      addWorksheet(wb, "p_vals")
      writeData(wb, "p_vals", per_cor$pvals)
      addWorksheet(wb, "lower_CI")
      writeData(wb, "lower_CI", per_cor$l_conf)
      addWorksheet(wb, "upper_CI")
      writeData(wb, "upper_CI", per_cor$u_conf)
      
      saveWorkbook(wb, file = paste(save_dir, analysis_type, '/', data_type, '_',con, ".xlsx", sep = ''))
  }
}

print(paste('all plots and corrlation tables for ', analysis_type, 'saved in', save_dir))






### Calculate correlations only hrs_incubation ####

params_corr =  setdiff(columns_of_i, c('day', 'treatment', 'OP', 'source'))

inc_0 <- fread(paste(data_dir, 'df_ctrl_0_inc.csv', sep = ''), sep = ",")
short_inc <- fread(paste(data_dir, 'df_ctrl_short_inc.csv', sep = ''), sep = ",")
long_inc <- fread(paste(data_dir, 'df_ctrl_long_inc.csv', sep = ''), sep = ",")

df_list <- list('no_inc' = inc_0,
                'short_inc' = short_inc,
                'long_inc' = long_inc)

for (data_type in names(df_list)) {
    
  df <- df_list[[data_type]][, ..params_corr]
  
  # calculate correlations
  # partial_corr_mixed(df, params_corr, '(1|OP)')
  pearson_cor <- pearson_corr(df, params_corr) # $ corrs, pvals, l_conf, u_conf
  # partial_cor <- partial_corr_mixed(df, params_corr, random_effect, '1|OP')
  
  plot_list <- list()
  plot_count <- 1
  
  for (i in 1:(length(params_corr)-1)) {
    for (j in (i+1):length(params_corr)) {
      param_x <- params_corr[i]
      param_y <- params_corr[j]
      
      # Get correlation values
      cor_val <- pearson_cor$corrs[param_x, param_y]
      p_val <- pearson_cor$pvals[param_x, param_y]
      lower_ci <- pearson_cor$l_conf[param_x, param_y]
      upper_ci <- pearson_cor$u_conf[param_x, param_y]
      
      # Create scatter plot with regression line and CI
      p <- ggplot(df, aes_string(x = param_x, y = param_y)) +
        geom_point(alpha = 0.6, size = 2) +
        geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
        labs(
          title = paste(param_x, "vs", param_y),
          subtitle = paste("r =", round(cor_val, 3), 
                           ", p =", round(p_val, 3),
                           "\nCI: [", round(lower_ci, 3), ",", round(upper_ci, 3), "]"),
          x = param_x,
          y = param_y
        ) +
        theme_minimal() +
        theme(plot.title = element_text(size = 10),
              plot.subtitle = element_text(size = 8))
      
      plot_list[[plot_count]] <- p
      plot_count <- plot_count + 1
    }
  }
  
  n <-  length(params_corr)
  layout_matr <- get_layout_matr(n)
  layout_matr[layout_matr == 0] <- NA
  layout_matr <- t(layout_matr)
  
  # Save combined plot
  combined_plot <- grid.arrange(grobs = plot_list, layout_matrix = layout_matr)
  
  # Save the plot
  plot_filename <- paste0(save_dir, analysis_type, '/', data_type, "_correlation.png")
  ggsave(plot_filename, combined_plot, 
         width = n * 4, height = n * 3, dpi = 300)
  
  wb <- createWorkbook()
  addWorksheet(wb, "cor_estimate")
  writeData(wb, "cor_estimate", pearson_cor$corrs)
  addWorksheet(wb, "p_vals")
  writeData(wb, "p_vals", pearson_cor$pvals)
  addWorksheet(wb, "lower_CI")
  writeData(wb, "lower_CI", pearson_cor$l_conf)
  addWorksheet(wb, "upper_CI")
  writeData(wb, "upper_CI", pearson_cor$u_conf)
  
  saveWorkbook(wb, file = paste(save_dir, analysis_type, '/', data_type, '.xlsx', sep = ''))

}