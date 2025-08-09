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

# AIS ANALYSIS ####
df_ais <- fread(paste(data_dir, 'AIS_all_data','.csv', sep = ''), sep = ",")

ais_anova_df <- data.frame()
R_sq <- data.frame()
df_ais <- fix_ais_df(df_ais)

ais_anova_df <- get_results_anova_model_comparison(lmers_ais(df_ais),
                                                   ais_anova_df, 'slice', 'AIS_len', 'hrs_incubation')

FullModel <- lmer(Length ~  treatment_r  + day + (1|OP), data = df_ais, REML = FALSE)

hrs_inc_attr <- attributes(scale(df_ais$org_hrs_incubation))
hrs_0 <- (0  - hrs_inc_attr$`scaled:center`) / hrs_inc_attr$`scaled:scale`
hrs_18 <- (18 - hrs_inc_attr$`scaled:center`) / hrs_inc_attr$`scaled:scale`

# EMMeans comparison
MarginalEffects <- emmeans(object = FullModel, 
                           specs =  ~ treatment_r + day,
                           type="response", tran = "log") #log transforms

#plot(MarginalEffects)
ContrastMatrix <- diag(4)
gr_tr <- MarginalEffects@grid$treatment_r
gr_day <- MarginalEffects@grid$day

### Contrast and adjustment for multiple testing (Benjamini-Hochberg) 
MarginalEffect_tr_test <- contrast(object = MarginalEffects, type = "response",
                                   method = list("(Ctrl 0) / (Ctrl 18)" = ContrastMatrix[gr_tr==0&gr_hrs==hrs_0,] - ContrastMatrix[gr_tr==0&gr_hrs==hrs_18,],
                                                 "(high K 0) / (high K 18)" = ContrastMatrix[gr_tr==0&gr_hrs==hrs_0,] - ContrastMatrix[gr_tr==1&gr_hrs==hrs_18,],
                                                 "(Ctrl 18) / (high K 18)" = ContrastMatrix[gr_tr==0&gr_hrs==hrs_18,] - ContrastMatrix[gr_tr==1&gr_hrs==hrs_18,]),
                                   adjust="BH")

MarginalEffect_tr_test <- contrast(object = MarginalEffects, type = "response",
                                   method = list("(Ctrl 0) / (Ctrl 18)" = ContrastMatrix[gr_tr==0&gr_day==0,] - ContrastMatrix[gr_tr==0&gr_day==1,],
                                                 "(high K 0) / (high K 18)" = ContrastMatrix[gr_tr==0&gr_day==0,] - ContrastMatrix[gr_tr==1&gr_day==1,],
                                                 "(Ctrl 18) / (high K 18)" = ContrastMatrix[gr_tr==0&gr_day==1,] - ContrastMatrix[gr_tr==1&gr_day==1,]),
                                   adjust="BH")

### Plot marginal effects
MarginalPlot <- plot(MarginalEffect_tr_test, colors = "grey10") +
  labs(title = paste(var, "Effects of Treatment and Day"))+
  ylab(label = "")+
  geom_vline(xintercept = 1, linetype="dashed")+
  xlab(label = "Ratio")+
  theme_classic()
MarginalPlot

m_effects <- summary(MarginalEffects)
emm_df <- rbind(emm_df , data.frame(
  data_type = rep(data_type, length(m_effects$day)),
  IV =  m_effects$day, 
  DV = rep(var, length(m_effects$day)),
  treatment =  m_effects$treatment_r,
  response = m_effects$response * (attrs$'scaled:scale') + attrs$'scaled:center',
  SE = m_effects$SE * (attrs$'scaled:scale') + attrs$'scaled:center',
  lower_CI =  m_effects$lower.CL * (attrs$'scaled:scale') + attrs$'scaled:center',
  upper_CI =  m_effects$upper.CL * (attrs$'scaled:scale') + attrs$'scaled:center'))


m_contrasts <- summary(MarginalEffect_tr_test)
emm_CI_df <- rbind(emm_CI_df, data.frame(
  data_type = rep(data_type, length(m_contrasts$contrast)),
  IV = rep('day', length(m_contrasts$contrast)), #repeat as long as the contrasts are
  DV = rep(var, length(m_contrasts$contrast)),
  ratio = m_contrasts$ratio,
  contrast = m_contrasts$contrast,
  se = m_contrasts$SE * (attrs$'scaled:scale') + attrs$'scaled:center',
  p_vales = m_contrasts$p.value))

r2 <- r2beta(FullModel)
R_sq <- data.frame(
  DV = rep('AIS_length',5),
  param = r2$Effect,
  R_squared = r2$Rsq,
  upper_CL = r2$upper.CL,
  lower_CL = r2$lower.CL)

save_dir <- '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/stats/AIS_MEA/'
# write.xlsx(ais_anova_df, paste(save_dir, 'AIS_model_comaprisons_intrinsic_.xlsx',sep = ''))
# write.xlsx(R_sq, paste(save_dir, 'AIS_explained_variance.xlsx',sep = ''))


## MEA SIMPLE ANALYSIS ####
df_mea <- fread(paste(data_dir, 'summary_pMEA_r.csv', sep = ''), sep = ",")

# fixing the dsataframe for analysis
df_mea$treatment_r <- df_mea$treatment # incubation solution
df_mea$treatment_r <- ifelse(df_mea$treatment == "Ctrl", 0, 1)
df_mea <- df_mea[!is.na(df_mea[["value"]]), ]
df_mea$value_org <- df_mea$value
df_mea = df_mea %>%
  mutate(value = scale(value))

mea_anova_results <- data.frame()
mea_R_sq <- data.frame()

mea_anova_results <- get_results_anova_model_comparison(lmers_mea_simple(df_mea),
                                                        mea_anova_results, 'slice', 'spikes/electrode(log)', 'trVScondition')

MEAFullModel <- lmer(value ~ treatment_r * patient_age + (1|OP),
                      data = df_mea, REML = FALSE)

MEA_r2 <- r2beta(MEAFullModel)
MEA_R_sq <- data.frame(
  DV = rep('spikes/electrode(log)',4 ),
  param = MEA_r2$Effect,
  R_squared = MEA_r2$Rsq,
  upper_CL = MEA_r2$upper.CL,
  lower_CL = MEA_r2$lower.CL)

# write.xlsx(mea_anova_results, paste(save_dir, 'MEA_model_comaprisons_intrinsic_.xlsx',sep = ''))
# write.xlsx(MEA_R_sq, paste(save_dir, 'MEA_explained_variance.xlsx',sep = ''))

### MEA - unnecessary=ily complicated
# df_mea <- fread(paste(data_dir, 'summary_pMEA_R.csv', sep = ''), sep = ",")
# 
# mea_anova_results <- data.frame()
# mea_R_sq <- data.frame()
# 
# df_mea <- df_mea[df_mea$Condition != 'washout', ]
# 
# df_mea$treatment_r <- df_mea$treatment # incubation solution
# df_mea$treatment_r <- ifelse(df_mea$treatment == "Ctrl", 0, 1) #
# 
# df_mea$condition_r <- df_mea$Condition # condition baseline - just after incubation, effects of time. Washout -> effects of treatment  solution 
# df_mea$condition_r <- ifelse(df_mea$Condition == "BL", 0, 1)

# mea_results <- glmers_mea(df_mea)
# 
# for (name in names(mea_results)) {
#   mea_anova_results <-  rbind(mea_anova_results, data.frame(
#     DV = rep('spikes/electrode(log)', 2),
#     ID2 = rep('trVScondition', 2),
#     model_comparison = rownames(mea_results[[name]]),
#     logLik = mea_results[[name]][['logLik']],
#     deviance = mea_results[[name]][['deviance']],
#     p_vals = mea_results[[name]][['Pr(>Chisq)']]))
# }
# 
# MEAFullModel <- glmer(value ~ condition_r + treatment_r + patient_age + (1|OP) + (1|OP:slice),
#                    data = df_mea, family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
# 
# MEA_r2 <- r2beta(MEAFullModel)
# MEA_R_sq <- data.frame(
#   DV = rep('spikes/electrode(log)',4 ),
#   param = MEA_r2$Effect,
#   R_squared = MEA_r2$Rsq,
#   upper_CL = MEA_r2$upper.CL,
#   lower_CL = MEA_r2$lower.CL)



# QC plots ####

# Count unique slices per OP and treatment
df_counts_ais <- df_ais %>%
  group_by(OP, treatment) %>%
  summarise(unique_slice_count = n_distinct(Slice), .groups = 'drop')

AIS_slice_per_op <- ggplot(df_counts_ais, aes(x = OP, y = unique_slice_count, fill = treatment)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Unique Slices per OP (stacked by treatment)", y = "Number of unique Slices") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/plots/QC_plots/slice/'
ggsave(filename = paste(save_dir, 'AIS_slice_per_OP.png', sep = ''),
       plot = AIS_slice_per_op, # or replace with your plot object, e.g., p
       width = 8, height = 6, dpi = 300, bg = "white")

# same for MEA data
# Count unique slices per OP and treatment
df_counts_mea <-  df_mea %>% 
  group_by(OP, treatment) %>%
  summarise(unique_slice_count = n_distinct(slice), .groups = 'drop')

MEA_slice_per_op <- ggplot(df_counts_mea, aes(x = OP, y = unique_slice_count, fill = treatment)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Unique Slices per OP (stacked by treatment)", y = "Number of unique Slices") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = paste(save_dir, 'MEA_slice_per_OP.png', sep = ''),
       plot = MEA_slice_per_op, # or replace with your plot object, e.g., p
       width = 8, height = 6, dpi = 300, bg = "white")
  
# ggplot(data, aes(x = Treatment, y = Length, color = Slice)) +
#   geom_jitter(position = position_dodge(width = 0.7), size = 2, alpha = 0.7) +
#   labs(title = "Axon Lengths by Condition and Slice",
#        x = "Condition",
#        y = "Axon Length (µm)",
#        color = "Slice (File)") +
#   theme_minimal() +
#   theme(legend.position = "right")