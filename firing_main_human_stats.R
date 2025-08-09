rm(list = ls())
# Load libraries 
library(lme4)
library(data.table)
library(emmeans)
library(dplyr)
library(tidyverse)
library(openxlsx)

setwd('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/code/R/human_slice/')
source("human_stats_functions.R")

data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/'
area <- 'temporal' # or 'complete'
df_repatch <- fread(paste(data_dir, 'repatch_data_', area,'.csv', sep = ''), sep = ",")
df_slice <- fread(paste(data_dir, 'slice_data_', area,'.csv', sep = ''), sep = ",")

id2 <- 'day'
var_of_interest <- c('num_aps', 'IFF')
# define ID variables
ids_basic <- c('OP', 'day', 'treatment_r', 'hrs_incubation', 'age_group')

anova_results <- data.frame()
marg_effects_contrasts <- data.frame()
marg_effects_conf_ints <- data.frame()
extended_results <- data.frame()
R_sq <- data.frame()

df_list <- list(cell_ID_new = df_repatch, slice = df_slice)
for (data_type in names(df_list)) {
  
  ids <- append(ids_basic, data_type, after = length(ids_basic))
  #get rid of NANs in amplitude
  # remove TTX
  df <- df_list[[data_type]]
  df[[data_type]] <- factor(df[[data_type]])
  df <- subset(df, treatment != 'TTX')
  df$day <- factor(df$day)
  df$day_cat <- df$day # keep the categorical data in a column
  df$day <- ifelse(df$day_cat == "D1", 0, 1)
  df$treatment_r <- as.numeric(df$treatment_r)
  df$inj_current <- as.numeric(df$inj_current)
  df$hrs_after_OP <- as.numeric(df$hrs_after_OP)
  df$patient_age <- as.numeric(df$patient_age)
  #grouping patient ages for more power
  df <- df %>% 
    mutate(age_group = case_when(
      patient_age >= 11 & patient_age <= 20 ~ 1,
      patient_age >= 21 & patient_age <= 30 ~ 2,
      patient_age >= 31 & patient_age <= 40 ~ 3,
      patient_age >= 41 & patient_age <= 50 ~ 4,
      patient_age >= 51 & patient_age <= 61 ~ 5,
      patient_age >= 61 & patient_age <= 80 ~ 6,
      TRUE ~ NA_real_  # Assign NA for ages outside the specified ranges
    ))
  setnames(df, colnames(df), make.names(colnames(df))) # make col names without spaces

  if ("frontal" %in% df$area) {
    print("frontal is present in the area column")
  }
  
  for (var_firing in var_of_interest) {
    
    col_names_inj <- grep(var_firing, colnames(df), value = TRUE)
    col_names_inj <- col_names_inj[6:length(col_names_inj)]
  
    # create new col_names
    new_col_names <- c()
    for (name in col_names_inj){
      
      stop_at <- nchar(name) - 15
      new_name <- substr(name, 4, stop_at)
      
      if (var_firing == 'IFF'){
        stop_at <- nchar(name) - 11
        new_name <- substr(name, 4, stop_at)
      }
      
      df[[new_name]] <- df[[name]]
      new_col_names <- append(new_col_names, new_name, after = length(new_col_names))
    }
    
    col_to_keep <- c(ids, new_col_names)
    melt_df <- df[, ..col_to_keep]
    melt_df <- melt(melt_df, id = ids, variable = 'inj_current', value = var_firing)
    melt_df$inj_current <- as.numeric(as.character(melt_df$inj_current))
    # melt_df[[var_firing]] <- as.numeric(melt_df[[var_firing]])
    
    # # calculate pearson for each value, see when drops
    # con_table <- data.frame()
    # for (val in unique(melt_df$inj_current)[-c(1,2)]){
    #   df_cor <- melt_df[inj_current <= val]
    #   x <- df_cor$inj_current
    #   y <-  df_cor[[var_firing]]
    #   correlation <- cor(x, y, method = "pearson")
    #   
    #   con_table<-  rbind(con_table, data.frame(
    #     up_to = val,
    #     cor = correlation
    #   ))
    # }
    # print(con_table)
    
    # take the inj_values where the pearson correlation is of highest value
    # inj_to_keep <- con_table$up_to[con_table$cor == max(con_table$cor)]
    # print(paste(var_firing, 'keeping up to inj current of', inj_to_keep, 'pA', data_type))
    
    # decidede to keep up to 600 for all
    inj_to_keep <- 600
    melt_df <- melt_df[inj_current <= inj_to_keep]
    
    # remove 0 values, aka where no APs fired
    melt_df <- melt_df[melt_df[[var_firing]] != 0]
    melt_df$inj_current <- as.numeric(melt_df$inj_current)
    print(paste(data_type, var_firing ,(prod(dim(melt_df)) / 7)))
    
    # Q-Q plot distributions to find distribution
    # for (col in col_names_inj){
    #   df[[col]] <- as.numeric(df[[col]])
    #   plot_distributions(df, col, 'day', 'repatch')
    # }
    
    results_firing <- lmers_firing(melt_df, var_firing, id2, data_type)
    
    # get results table from model comparison
    for (name in names(results_firing)) {
      anova_results <-  rbind(anova_results, data.frame(
        data_type = rep(data_type, 2),
        DV = rep(var_firing,2),
        ID2 = rep(id2,2),
        model_comparison = rownames(results_firing[[name]]),
        logLik = results_firing[[name]][['logLik']],
        deviance = results_firing[[name]][['deviance']],
        p_vals = results_firing[[name]][['Pr(>Chisq)']]))
    }
  
    marg_effects <- marginal_effects_day_firing (melt_df, var_firing, data_type)
    
    marg_effects_contrasts <- rbind(marg_effects_contrasts, data.frame(
      data_type = rep(data_type, length(marg_effects$test$contrast)),
      IV = rep(id2, length(marg_effects$test$contrast)), #repeat as long as the contrasts are
      DV = rep(var_firing, length(marg_effects$test$contrast)),
      contrast = marg_effects$test$contrast,
      se = marg_effects$test$SE,
      p_vales = marg_effects$test$p.value))
  
    marg_effects_conf_ints <- rbind(marg_effects_conf_ints , data.frame(
      data_type = rep(data_type, length(marg_effects$effects[[id2]])),
      IV = marg_effects$effects[[id2]], 
      DV = rep(var_firing, length(marg_effects$effects[[id2]])),
      treatment = marg_effects$effects$treatment_r,
      lower_CI = marg_effects$effects$lower.CL,
      upper_CI = marg_effects$effects$upper.CL))
  
    
    formula_full <- as.formula(paste(var_firing, "~ inj_current + treatment_r * age_group + hrs_incubation + (1|OP) + (1|OP:", data_type, ")"))
    FullModel <-lmer(formula_full, data = melt_df, REML = FALSE)
    
    sum_model <- summary(FullModel)
    extended_results <- rbind(extended_results, data.frame(
      data_type = rep(data_type, 6),
      DV = rep(var_firing,6),
      param = rownames(sum_model$coefficients),
      estimate = sum_model$coefficients[,1],
      SE = sum_model$coefficients[,2],
      dfs = sum_model$coefficients[,3],
      t_val = sum_model$coefficients[,4],
      p_val = sum_model$coefficients[,5]))
    
    # extended_results <- lmers_firing_extended(melt_df, var_firing, data_type) # get results table from model comparison
    
    r2 <- r2beta(FullModel)
    R_sq <- rbind(R_sq, data.frame(
      data_type = rep(data_type, 6),
      DV = rep(var_firing,6),
      param = r2$Effect,
      R_squared = r2$Rsq,
      upper_CL = r2$upper.CL,
      lower_CL = r2$lower.CL))
    
    # for (name in names(extended_results)) {
    #   ext_anova_results <-  rbind(ext_anova_results, data.frame(
    #     data_type = rep(data_type, 2),
    #     DV = rep(var_firing,2),
    #     ID2 = rep('hrs',2),
    #     model_comparison = rownames(extended_results [[name]]),
    #     logLik = extended_results [[name]][['logLik']],
    #     deviance = extended_results [[name]][['deviance']],
    #     p_vals = extended_results [[name]][['Pr(>Chisq)']]))
    #   }
    
  }
  rm(ids)
}

save_dir <- '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/stats/'
# write.xlsx(anova_results, paste(save_dir, 'model_comaprisons_firing_',area ,'.xlsx',sep = ''))
# write.xlsx(marg_effects_contrasts, paste(save_dir, 'contrasts_firing_',area ,'.xlsx',sep = ''))
# write.xlsx(marg_effects_conf_ints, paste(save_dir, 'conf_intervals_firing_',area ,'.xlsx',sep = ''))
#write.xlsx(ext_anova_results, paste(save_dir, 'extended_firing_',area ,'.xlsx',sep = ''))
# write.xlsx(extended_results, paste(save_dir, '2_extended_firing_',area ,'.xlsx',sep = ''))
# write.xlsx(R_sq, paste(save_dir, 'explained_var_firing_', area,'.xlsx',sep = ''))

