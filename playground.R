#### ff ####

## Checking if for rheobase the hrs_incubation better explain the results
id2 <- 'day'
var_of_interest <- c('num_aps', 'IFF')
df_list <- list(cell_ID_new = df_repatch, slice = df_slice)

for (data_type in names(df_list)) {
  df <- df_list[[data_type]]
  df <- subset(df, treatment != 'TTX')
  df$day <- as.factor(df$day)
  setnames(df, colnames(df), make.names(colnames(df))) # make col names without spaces
  var <- 'rheo_ramp_c'
  results_day <- lmers(df, var, 'day', data_type)
  results_hrs <- lmers(df, var, 'hrs_incubation', data_type)
  contr_day <- marginal_effects_day(df, var, 'gaussian', data_type)
  contr_hrs <- marginal_effects_hrs(df, var, 'gaussian', data_type)
  
  # get results table from model comparison
  for (name in names(results_day)) {
    anova_results <-  rbind(anova_results, data.frame(
      data_type = rep(data_type, 2),
      DV = rep(var,2),
      ID2 = rep('day',2),
      model_comparison = rownames(results_day[[name]]),
      logLik = results_day[[name]][['logLik']],
      deviance = results_day[[name]][['deviance']],
      p_vals = results_day[[name]][['Pr(>Chisq)']]))
  }
  for (name in names(results_hrs)) {
    anova_results <-  rbind(anova_results, data.frame(
      data_type = rep(data_type, 2),
      DV = rep(var,2),
      ID2 = rep('hrs',2),
      model_comparison = rownames(results_hrs[[name]]),
      logLik = results_hrs[[name]][['logLik']],
      deviance = results_hrs[[name]][['deviance']],
      p_vals = results_hrs[[name]][['Pr(>Chisq)']]))
  }
  
  contr_list <- list(day = contr_day, hrs = contr_hrs)
  for (contr in names(contr_list)){
    marg_effects <- contr_list[[contr]]
    marg_effects_contrasts <- rbind(marg_effects_contrasts, data.frame(
      data_type = rep(data_type, length(marg_effects$test$contrast)),
      IV = rep(contr, length(marg_effects$test$contrast)),#repeat as long as the contrasts are
      DV = rep(var, length(marg_effects$test$contrast)),
      contrast = marg_effects$test$contrast,
      se = marg_effects$test$SE,
      p_vales = marg_effects$test$p.value))
  }

}





data_d1 <- melt_df[melt_df$day == 'D1']
data_d2 <- melt_df[melt_df$day == 'D2']

# basic plot for checking data
ggplot(data_d2, aes(x = inj_current, y = num_aps, color = treatment)) +
  # geom_point(size = 3) +  # Add points
  stat_summary(fun =median, geom = "point")+
  labs(
    title = "Number of Action Potentials vs Injection Current",
    x = "Injection Current (pA)",
    y = "Number of Action Potentials"
  ) +
  theme_minimal()




# variance 

distribution_dict <- list(
  'TH' = 'gamma',
  'rheo_ramp_c' = 'gaussian')

id2 <- 'day'
var_of_interest <- c('num_aps', 'IFF')

anova_results <- data.frame()
marg_effects_contrasts <- data.frame()
marg_effects_conf_ints <- data.frame()

df_list <- list(cell_ID_new = df_repatch, slice = df_slice)

formula_full <- as.formula(paste(var, "~ treatment *",id2,"+ (1|OP) + (1|OP:", data_type, ")"))

# random interceptt for age
# random intercept for sex
# hrs_incubation - continuous
formula_full <- as.formula(paste(var, "~ treatment *",id2,"+ (1|patient_age) +(1|OP) + (1|OP:", data_type, ")"))





# FUNCS FOR REPATCH ONLY
# funcs analysis intrinsic excitability
# keeping only treatment and the interaction with day
# how the models were initially designed

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

setwd('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/code/R/human_slice/')
source("human_stats_functions.R")

data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/'
area <- 'temporal' # or 'complete'
df_repatch <- fread(paste(data_dir, 'repatch_data_', area,'.csv', sep = ''), sep = ",")
df_slice <- fread(paste(data_dir, 'slice_data_', area,'.csv', sep = ''), sep = ",")


lmers_repatch <- function(df, var, id2, data_type) {
  # make positive
  if (mean(df[[var]], na.rm = TRUE) < 0) {
    df[[var]] = -df[[var]] 
  }
  #remove rows with nan values
  df <- df[!is.na(df[[var]]), ]
  
  formula_full <- as.formula(paste(var, "~ treatment *", id2, "+ (1|OP) + (1|OP:", data_type, ")"))
  formula_tr_0 <- as.formula(paste(var, "~ ", id2 ," + (1|OP) + (1|OP:", data_type, ")"))
  formula_day_0 <- as.formula(paste(var, "~ treatment + (1|OP) + (1|OP:", data_type, ")"))
  formula_interaction_0 <- as.formula(paste(var, "~ treatment +",id2, "+ (1|OP) + (1|OP:", data_type, ")"))
  
  
  FullModel <- lmer(formula_full, data = df, REML = FALSE)
  TreatmentNullModel <-lmer(formula_tr_0, data = df, REML = FALSE)
  dayNullModel <- lmer(formula_day_0, data = df, REML = FALSE)
  interactionNullModel <- lmer(formula_interaction_0, data = df, REML = FALSE)
  
  anova_full_tr_0 <- anova(FullModel, TreatmentNullModel)
  anova_full_day_0 <- anova(FullModel, dayNullModel)
  anova_full_interaction_0 <- anova(FullModel, interactionNullModel)
  
  return(list(
    full_vs_treat_null = anova_full_tr_0, 
    full_vs_day_null = anova_full_day_0,
    full_vs_interact_null = anova_full_interaction_0))
}

glmers_repatch <- function(df, var, id2, data_type) {
  # make positive
  if (mean(df[[var]], na.rm = TRUE) < 0) {
    df[[var]] = -df[[var]] 
  }
  # remove rows with nans
  df <- df[!is.na(df[[var]]), ]
  
  formula_full <- as.formula(paste(var, "~ treatment * ",id2,"+ (1|OP) + (1|OP:", data_type, ")"))
  formula_tr_0 <- as.formula(paste(var, "~", id2, " + (1|OP) + (1|OP:", data_type, ")"))
  formula_day_0 <- as.formula(paste(var, "~ treatment + (1|OP) + (1|OP:", data_type, ")"))
  formula_interaction_0 <- as.formula(paste(var, "~ treatment +",id2, "+ (1|OP) + (1|OP:", data_type, ")"))
  
  FullModel <- glmer(formula_full, data = df,
                     family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  TreatmentNullModel <- glmer(formula_tr_0, data = df,
                              family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  dayNullModel <- glmer(formula_day_0, data = df,
                        family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  interactionNullModel <- glmer(formula_interaction_0, data = df,
                        family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  
  anova_full_tr_0 <- anova(FullModel, TreatmentNullModel)
  anova_full_day_0 <- anova(FullModel, dayNullModel)
  anova_full_interaction_0 <- anova(FullModel, interactionNullModel)
  
  return(list(
    full_vs_treat_null = anova_full_tr_0, 
    full_vs_day_null = anova_full_day_0,
    full_vs_interact_null = anova_full_interaction_0))
}


distribution_dict <- list(
  'TH' = 'gamma',
  'rheo_ramp_c' = 'gaussian',
  'sag' = 'gaussian')

id2 <- 'day'
df_list <- list(cell_ID_new = df_repatch, slice = df_slice)


data_type <- 'cell_ID_new'

df <- df_list[[data_type]]
df <- subset(df, treatment != 'TTX')
df$day <- factor(df$day, levels = c('D1', 'D2'), labels = c('D1', 'D2'))
df$treatment <- factor(df$treatment, levels = c('Ctrl', 'high K'), labels = c('Ctrl', 'high K'))
setnames(df, colnames(df), make.names(colnames(df))) # make col names without spaces

anova_results <- data.frame()

for (var in names(distribution_dict)) {
  
  # choose the correct model
  if (distribution_dict[[var]] == 'gaussian'){
    results <- lmers_repatch(df, var, id2, data_type)
  } else if (distribution_dict[[var]] == 'gamma'){
    results <- glmers_repatch(df, var, id2, data_type)
  }
  
  # get results table from model comparison
  for (name in names(results)) {
    anova_results <-  rbind(anova_results, data.frame(
      data_type = rep(data_type, 2),
      DV = rep(var,2),
      ID2 = rep(id2,2),
      model_comparison = rownames(results[[name]]),
      logLik = results[[name]][['logLik']],
      deviance = results[[name]][['deviance']],
      p_vals = results[[name]][['Pr(>Chisq)']]))
  }
}

save_dir <- '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/stats/'
write.xlsx(anova_results, paste(save_dir, 'REPATCH_SIGN_model_comaprisons_intrinsic_', area,'.xlsx',sep = ''))




# EXPLAINED VARIANCE IN IFF and NUM_APS


df <- df_repatch
df$day_cat <- df$day # keep the categorical data in a column
df$day <- ifelse(df$day_cat == "D1", 0, 1)

# Create a scatterplot
ggplot(df, aes(x = day, y = TH, color = factor(day))) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("blue", "red"), labels = c("D1 (0)", "D2 (1)")) +
  labs(title = "Scatterplot with Day Coded as Numeric",
       x = "X-axis Label",
       y = "Y-axis Label",
       color = "Day")

ggplot(df, aes(x = day_cat, y = TH, color = factor(day))) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("blue", "red"), labels = c("D1 (0)", "D2 (1)")) +
  labs(title = "Scatterplot with Day Coded as Numeric",
       x = "X-axis Label",
       y = "Y-axis Label",
       color = "Day")



##### 29.05  - plotting effects understanding stuff

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

setwd('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/code/R/human_slice/')
source("human_stats_functions.R")

data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/'
area <- 'temporal' # or 'complete'
# df_repatch <- fread(paste(data_dir, 'repatch_data_', area,'.csv', sep = ''), sep = ",")
df_slice <- fread(paste(data_dir, 'slice_data_', area,'.csv', sep = ''), sep = ",")

### choosing params
var <- 'TH'
data_type <- 'slice'
id2 <- 'day'
df <- df_slice

### df preparation 
df$day <- factor(df$day)
df$treatment_r <- as.numeric(df$treatment_r)
df$OP <- factor(df$OP)
df$slice <- factor(df$slice)
df <- within(df, sample <- factor(OP:slice))
df$day_cat <- df$day # keep the categorical data in a column
df$day <- ifelse(df$day_cat == "D1", 0, 1)
setnames(df, colnames(df), make.names(colnames(df))) 

df[[var]] = -df[[var]]

### from glmers function

formula_full <- as.formula(paste(var, "~ treatment_r + ",id2,"+ (1|OP) + (1|OP:", data_type, ")"))
formula_tr_0 <- as.formula(paste(var, "~", id2, " + (1|OP) + (1|OP:", data_type, ")"))
formula_day_0 <- as.formula(paste(var, "~ treatment_r + (1|OP) + (1|OP:", data_type, ")"))

FullModel <- glmer(formula_full,data = df, 
                   family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
TreatmentNullModel <- glmer(formula_tr_0,data = df, 
                            family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
DayNullModel <- glmer(formula_day_0,  data = df,
                      family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))


anova_full_tr_0 <- anova(FullModel, TreatmentNullModel)
anova_full_day_0 <- anova(FullModel, DayNullModel)

### plot
(plot <- ggplot(df, aes(x = treatment_r, y = TH , colour = OP)) +
    # facet_wrap(~OP, nrow=2) +   # a panel for each mountain range
    # geom_point(alpha = 0.5) +
    theme_classic() +
    geom_line(data = cbind(df, pred = predict(FullModel)), aes(y = pred), size = 1) +  # adding predicted line from mixed model 
    theme(legend.position = "none",
          panel.spacing = unit(2, "lines"))  # adding space between panels
)





### save those plots

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

setwd('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/code/R/human_slice/')
source("human_stats_functions.R")

data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/'
area <- 'temporal' # or 'complete'
df_repatch <- fread(paste(data_dir, 'repatch_data_', area,'.csv', sep = ''), sep = ",")
df_slice <- fread(paste(data_dir, 'slice_data_', area,'.csv', sep = ''), sep = ",")

distribution_dict <- list(
  'TH' = 'gamma',
  'resting_potential' = 'gamma',
  'membra_time_constant_tau' = 'gaussian',
  'Rin' = 'gamma',
  'rheo_ramp_c' = 'gaussian',
  'sag' = 'gaussian')

id2 <- 'day'


df_list <- list(cell_ID_new = df_repatch, slice = df_slice)
for (data_type in names(df_list)) {
  df <- df_list[[data_type]]
  df[[data_type]] <- factor(df[[data_type]])
  df <- fix_df(df)
  for (var in names(distribution_dict)) {
    plot_distributions(df, var, 'day', data_type)
  }
}

### exteded contrast ####
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
area <- 'temporal' # or 'complete'
df_slice <- fread(paste(data_dir, 'slice_data_', area,'.csv', sep = ''), sep = ",")
df_repatch <- fread(paste(data_dir, 'repatch_data_', area,'.csv', sep = ''), sep = ",")

df <- df_ext_comparison_groups(df_repatch)
df$treatment <- as.factor(df$treatment)
# df$day_cat <- df$day # keep the categorical data in a column
# df$day <- ifelse(df$day_cat == "D1", 0, 1)

data_type <- 'cell_ID_new'
var <- 'TH'

formula_full <- as.formula(paste(var, "~ treatment  * hrs_group + (1|OP) + (1|OP:", data_type, ")"))
FullModel <- lmer(formula_full, data = df, control =  lmerControl(optimizer = "Nelder_Mead"))

summary(FullModel)

# here for TH maybe an interesting contrast start end and treatment
# no interaction effect of hrs_after_OP and treatment




### play ####

df_slice <- fread(paste(data_dir, 'slice_data_', area,'.csv', sep = ''), sep = ",")

ids_basic <- c('OP', 'day', 'treatment_r', 'hrs_after_OP', 'patient_age')
ids_ext_b <- c('OP', 'day', 'treatment', 'hrs_group', 'age_group')

ids <- append(ids_basic, data_type, after = length(ids_basic))
ids_ext <- append(ids_ext_b, data_type, after = length(ids_ext_b))

data_type <- 'slice'
var_firing <- 'num_aps'

df <- df_list[[data_type]]
df[[data_type]] <- factor(df[[data_type]])
df_ext <- df_ext_comparison_groups(df)
df <- fix_df(df)

melt_df <- get_melt_df_firing(var_firing, df, ids)

# scaled df for model stability
df_scaled <- melt_df 
df_scaled$VAL <- as.vector(scale(df_scaled$VAL))
df_scaled$inj_current <- as.vector(scale(df_scaled$inj_current))

df_scaled <- melt_df
ar_org <- melt_df$VAL
attrs <- attributes(scale(melt_df$VAL))

formula_full <- as.formula(paste("VAL ~ inj_current + treatment_r + day  + (1|OP) + (1|OP:", data_type, ")"))
FullModel <- lmer(formula_full, data = df_scaled)

# Post-hoc analysis with Estimated Marginal Means ####
MarginalEffects <- emmeans(object = FullModel, 
                           specs =  ~ inj_current + treatment_r + day,
                           at = list(inj_current = c(200, 400, 600)),
                           tran = "link") # no transformation 
MarginalEffects



#### f ####
# EMMS_short_change <- 

short_emms_difference <- function(df, var_firing, data_type){
  formula_full <- as.formula(paste("VAL ~ inj_current + treatment_r + day  + (1|OP) + (1|OP:", data_type, ")"))
  FullModel <- lmer(formula_full, data = df)
  
  # Post-hoc analysis with Estimated Marginal Means ####
  MarginalEffects <- emmeans(object = FullModel, 
                             specs =  ~ inj_current + treatment_r + day,
                             at = list(inj_current = c(200, 400, 600)))
  
  ContrastMatrix <- diag(8)
  gr_day <- MarginalEffects@grid$day
  gr_tr <- MarginalEffects@grid$treatment_r
  gr_inj <- MarginalEffects@grid$inj_current
  
  MarginalEffectTest <- contrast(object = MarginalEffects, type = "response",
                                 method = list("400 to 600 CTR" = (ContrastMatrix[gr_inj==400&gr_day==0&gr_tr==0,] - ContrastMatrix[gr_inj==600&gr_day==0&gr_tr==0,]) -
                                                 (ContrastMatrix[gr_inj==400&gr_day==1&gr_tr==0,] - ContrastMatrix[gr_inj==600&gr_day==1&gr_tr==0,]),
                                               "400 to 600 HiK" = (ContrastMatrix[gr_inj==400&gr_day==0&gr_tr==0,] - ContrastMatrix[gr_inj==600&gr_day==0&gr_tr==0,]) -
                                                 (ContrastMatrix[gr_inj==400&gr_day==1&gr_tr==1,] - ContrastMatrix[gr_inj==600&gr_day==1&gr_tr==1,])),
                                 adjust="BH")
  
}




#### for plotting in python ####
library("plotrix")
setwd('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/code/R/human_slice/')
source("funcs_human_stats.R")

data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/'
area <- 'temporal' # or 'complete'
df_repatch <- fread(paste(data_dir, 'repatch_data_', area,'.csv', sep = ''), sep = ",")
df_slice <- fread(paste(data_dir, 'slice_data_', area,'.csv', sep = ''), sep = ",")

df_repatch$cell_ID_new <- factor(df_repatch$cell_ID_new)
# df_repatch <- fix_df(df_repatch)
df_repatch <- df_ext_comparison_groups(df_repatch)

df_slice$slice <- factor(df_slice$slice)
# df_slice <- fix_df(df_slice)
df_slice <- df_ext_comparison_groups(df_slice)

df_list <- list(cell_ID_new = df_repatch, slice = df_slice)
hrs_complete <- data.frame()
age_complete <- data.frame()
for (data_type in names(df_list)) {
  df <- df_list[[data_type]]
  df[[data_type]] <- factor(df[[data_type]])
  df <- df_ext_comparison_groups(df)
  
  hrs_df <- as.data.frame(df %>%
    group_by(treatment_word, day, hrs_group) %>%
    summarize(mean_hrs = mean(hrs_after_OP),
              median_hrs = median(hrs_after_OP),
              SE_hrs = std.error(hrs_after_OP),
              CI_hrs = round(confint(lm(hrs_after_OP ~ 1), level=0.95),3)))
  hrs_df$data_type <- rep(data_type, nrow(hrs_df))
  hrs_complete <- rbind(hrs_complete, hrs_df)
  
  age_df <- as.data.frame(df %>%
     group_by(treatment_word, day, age_group) %>%
     summarize(means_age = mean(patient_age),
               median_age = median(patient_age),
               SE_age = std.error(patient_age),
               CI_age = round(confint(lm(patient_age ~ 1), level=0.95),3)))
  
  age_df$data_type <- rep(data_type, nrow(age_df))
  age_complete <- rbind(age_complete, age_df)
}



#### CALCULATING SUMMARY DF FOR PLOTTING FIRING ####
library(plotrix )
setwd('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/code/R/human_slice/')
source("funcs_human_stats.R")

data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/'
# area <- 'temporal' # or 'complete'
# df_repatch <- fread(paste(data_dir, 'repatch_data_', area,'.csv', sep = ''), sep = ",")
# df_slice <- fread(paste(data_dir, 'slice_data_', area,'.csv', sep = ''), sep = ",")

df_slice_inc_only <- fread(paste(data_dir, 'slice_all.csv', sep = ''), sep = ",")
filter_hrs = 31
data_type <- paste('slice under ', filter_hrs, 'hrs')
df_slice_inc_only_f <- df_slice_inc_only[df_slice_inc_only$hrs_after_OP < filter_hrs, ]

var_of_interest <- c('num_aps', 'IFF')
# define ID variables
ids <- c('OP', 'day', 'treatment', 'hrs_after_OP', 'slice')

df_slice_inc_only_f[['slice']] <- factor(df_slice_inc_only_f[['slice']])
df <- fix_df(df_slice_inc_only_f)
CIs_complete <- data.frame()

for (var_firing in var_of_interest){
  melt_df <- get_melt_df_firing(var_firing, df, ids)
  
  CIs_df <- as.data.frame(melt_df %>%
                            group_by(treatment, day, inj_current) %>%
                            summarize(means = mean(VAL),
                                      median = median(VAL),
                                      SE = std.error(VAL),
                                      CI_L = confint(lm(VAL ~ 1), level=0.95)[1],
                                      CI_U = confint(lm(VAL ~ 1), level=0.95)[2]))
  
  CIs_df$data_type <- rep(data_type, nrow(CIs_df))
  CIs_df$var_firing <- rep(var_firing, nrow(CIs_df))
  CIs_complete <- rbind(CIs_complete, CIs_df)
}

save_dir <- '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/stats/sum_data_CIs_SEs_for_plotting/'
# write.xlsx(CIs_complete, paste(save_dir, 'firing_plot_CIs_', area,'.xlsx',sep = ''))
# write.xlsx(CIs_complete, paste(save_dir, 'firing_plot_CIs_inc_only.xlsx',sep = ''))


### plotting firing properties incubation only



df_slice_inc_only_f[['slice']] <- factor(df_slice_inc_only_f[['slice']])
df <- fix_df(df_slice_inc_only_f)
melt_df <- get_melt_df_firing(var_firing, df, ids)
var_firing <- 'IFF'
ids <- c('OP', 'day', 'treatment_r', 'hrs_after_OP', 'slice')

melt_df <- get_melt_df_firing(var_firing, df, ids)
# Calculate means for each treatment/day/injection combination
df_means <- melt_df %>%
  group_by(treatment_r, day, inj_current) %>%
  summarise(mean_VAL = mean(VAL, na.rm = TRUE), .groups = "drop")


df_means_2 <- melt_df %>% group_by(treatment_r, day, inj_current) %>%
            summarize(means = mean(VAL),
                      median = median(VAL),
                      SE = std.error(VAL),
                      CI_L = confint(lm(VAL ~ 1), level=0.95)[1],
                      CI_U = confint(lm(VAL ~ 1), level=0.95)[2])

summary_df <- fread('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/stats/inc_only/summary_df_inc_only_plot.csv')
# Plot function

plot_1 <- ggplot(melt_df, aes(x = inj_current, y = VAL, color = interaction(treatment_r, day))) +
  # geom_point(alpha = 0.3, position = position_jitter(width = 0.2)) +
  geom_line(data = df_means, aes(x = inj_current, y = mean_VAL, group = interaction(treatment_r, day)), size = 1) +  # Lines through means
  geom_line(data = df_means_2, aes(x = inj_current, y = means, group = interaction(treatment_r, day)), size = 1) +
  geom_errorbar(data = df_means_2, aes(x = inj_current, ymin = means - SE, ymax = means + SE), width = 0.2, alpha = 0.5)+
  # geom_ribbon(data = df_summary, aes(x = inj, ymin = CI_L, ymax = CI_U, fill = interaction(treatment, day)), alpha = 0.3)
  # facet_wrap(~ treatment_r + day, scales = "free_y") +
  labs(
    x = "Current Injection (pA)",
    y = "Firing Frequency (Hz)",
    title = paste("Firing Properties - Error Type:")
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

plot_1




#### CALCULATING SUMMARY DF FOR PLOTTING FIRING ####
rm(list = ls())
library(plotrix )
setwd('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/code/R/human_slice/')
source("funcs_human_stats.R")

data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/inc_only_for_python_plot_delete_soon.csv'
area <- 'temporal' # or 'complete'
# df_repatch <- fread(paste(data_dir, 'repatch_data_', area,'.csv', sep = ''), sep = ",")
# df_slice <- fread(paste(data_dir, 'slice_data_', area,'.csv', sep = ''), sep = ",")

#df_slice_inc_only <- fread(paste(data_dir, 'slice_all.csv', sep = ''), sep = ",")
df_slice_inc_only <- fread(data_dir, sep = ",")
filter_hrs = 35
data_type <- paste('slice incubation only')
df_slice_inc_only_f <- df_slice_inc_only[df_slice_inc_only$hrs_after_OP < filter_hrs, ]

var_of_interest <- c('num_aps', 'IFF')
# define ID variables
ids <- c('OP', 'treatment_r', 'hrs_after_OP', 'slice')

df_slice_inc_only_f[['slice']] <- factor(df_slice_inc_only_f[['slice']])
df <- fix_df(df_slice_inc_only_f)
CIs_complete <- data.frame()

for (var_firing in var_of_interest){
  
  melt_df <- get_melt_df_firing(var_firing, df, ids)
  
  CIs_df <- as.data.frame(melt_df %>%
                            group_by(treatment_r, inj_current) %>%
                            summarize(means = mean(VAL),
                                      median = median(VAL),
                                      SE = std.error(VAL),
                                      CI_L = confint(lm(VAL ~ 1), level=0.95)[1],
                                      CI_U = confint(lm(VAL ~ 1), level=0.95)[2],
                                      .groups = "drop"))

  CIs_df$data_type <- rep(data_type, nrow(CIs_df))
  CIs_df$var_firing <- rep(var_firing, nrow(CIs_df))
  CIs_complete <- rbind(CIs_complete, CIs_df)
}

save_dir <- '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/stats/sum_data_CIs_SEs_for_plotting/'
write.xlsx(CIs_complete, paste(save_dir, 'firing_plot_temporary_inc_only_', area,'.xlsx',sep = ''))



# plot here

# plotting IV curves
# Ensure treatment_r is treated as a factor (dichotomous variable)
ids <- c('OP', 'day', 'treatment_r', 'hrs_after_OP', 'slice')
melt_df <- get_melt_df_firing('num_aps', df_slice_inc_only_f, ids)

melt_df$treatment_r <- as.factor(melt_df$treatment_r)

# Calculate averages for inj_current and VAL grouped by treatment_r
df_summary <- melt_df %>%
  group_by(treatment_r, inj_current) %>%
  summarise(mean_VAL = mean(VAL, na.rm = TRUE), .groups = "drop")

# Create the plot
ggplot(melt_df, aes(x = inj_current, y = VAL, color = treatment_r)) +
  geom_point(alpha = 0.5) +  # Scatter plot for individual points
  geom_point(data = df_summary, aes(x = inj_current, y = mean_VAL), size = 3, shape = 21, fill = "white") +  # Dots for averages
  geom_line(data = df_summary, aes(x = inj_current, y = mean_VAL, group = treatment_r)) +  # Line through averages
  labs(title = "VAL vs inj_current by Treatment",
       x = "Injection Current",
       y = "VAL",
       color = "Treatment") +
  theme_minimal()















