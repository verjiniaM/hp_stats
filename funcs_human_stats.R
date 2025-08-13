# Load libraries 

plot_distributions <- function(df, var, id2, data_type) {
  #remove rows with nan values
  df <- df[!is.na(df[[var]]), ]
  # if 0 in columns, change to a very small value
  if (any(df[[var]] == 0, na.rm = TRUE)) {
    df[[var]][!is.na(df[[var]]) & df[[var]] == 0] <- 1e-6
  }
  
  if (mean(df[[var]], na.rm = TRUE) < 0) {
    df[[var]] = -df[[var]] 
  }
  formula <- as.formula(paste(var, "~ treatment_r +", id2,"+ (1|OP) + (1|OP:",data_type,")"))
  
  # Fit the model
  FullModelLmer <- lmer(formula, data = df)
  FullModel_gamma <- glmer(formula, data = df, family = Gamma(link = "log"), 
                           control = glmerControl(nAGQ0initStep = 0))
  FullModel_poisson <- glmer(formula, data = df, family = poisson(link = "log"), 
                             control = glmerControl(nAGQ0initStep = 0))
  
  df[,`:=`(`Scaled Gaussian Residuals`=scale(residuals(FullModelLmer)),
           `Scaled Gamma Residuals`=scale(residuals(FullModel_gamma)),
           `Scaled Poisson Residuals`=scale(residuals(FullModel_poisson))),]
  
  # reshaping of the data
  ResidualDT <- melt(df, measure.vars = c("Scaled Gaussian Residuals", "Scaled Gamma Residuals",
                                          "Scaled Poisson Residuals"),
                     id.vars = c("OP", "treatment_r"))
  
  # plot
  plt1 <- ggplot(data = ResidualDT, aes(sample = value))+
    stat_qq()+
    stat_qq_line()+
    facet_wrap(~variable)+
    scale_x_continuous(name = "Theoretical quantiles")+
    scale_y_continuous(name = "Scaled Residuals")+
    ggtitle(paste("Q-Q Plot for", data_type, var)) +
    theme_classic()
  
  print(plt1)
}

# funcs analysis intrinsic excitability
lmers <- function(df, var, id2, data_type) {
  # make positive
  if (mean(df[[var]], na.rm = TRUE) < 0) {
    df[[var]] = -df[[var]] 
  }
  #remove rows with nan values
  df <- df[!is.na(df[[var]]), ]
  
  formula_full <- as.formula(paste(var, "~ treatment_r +", id2, "+ (1|OP) + (1|OP:", data_type, ")"))
  formula_tr_0 <- as.formula(paste(var, "~ ", id2 ," + (1|OP) + (1|OP:", data_type, ")"))
  formula_day_0 <- as.formula(paste(var, "~ treatment_r + (1|OP) + (1|OP:", data_type, ")"))
  #formula_interact_0 <- as.formula(paste(var, "~ treatment_r +", id2, "+ (1|OP) + (1|OP:", data_type, ")"))

  FullModel <- lmer(formula_full,data = df, REML = FALSE)
  TreatmentNullModel <-lmer(formula_tr_0, data = df, REML = FALSE)
  DayNullModel <- lmer(formula_day_0, data = df, REML = FALSE)
  #InteractNullModel <- lmer(formula_interact_0, data = df, REML = FALSE)
  
  anova_full_tr_0 <- anova(FullModel, TreatmentNullModel)
  anova_full_day_0 <- anova(FullModel, DayNullModel)
  #anova_full_interact_0 <- anova(FullModel, InteractNullModel)

  return(list(
    full_vs_treat_null = anova_full_tr_0, 
    full_vs_day_null = anova_full_day_0))
    # full_vs_interact_null = anova_full_interact_0))
}

glmers <- function(df, var, id2, data_type) {
  # make positive
  if (mean(df[[var]], na.rm = TRUE) < 0) {
    df[[var]] = -df[[var]] 
  }
  # remove rows with nans
  df <- df[!is.na(df[[var]]), ]
  
  formula_full <- as.formula(paste(var, "~ treatment_r + ",id2,"+ (1|OP) + (1|OP:", data_type, ")"))
  formula_tr_0 <- as.formula(paste(var, "~", id2, " + (1|OP) + (1|OP:", data_type, ")"))
  formula_day_0 <- as.formula(paste(var, "~ treatment_r + (1|OP) + (1|OP:", data_type, ")"))
  # formula_interact_0 <- as.formula(paste(var, "~ treatment_r +", id2, "+ (1|OP) + (1|OP:", data_type, ")"))

  FullModel <- glmer(formula_full,data = df, 
                     family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  TreatmentNullModel <- glmer(formula_tr_0,data = df, 
                              family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  DayNullModel <- glmer(formula_day_0,  data = df,
                        family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  #InteractNullModel <- glmer(formula_interact_0,  data = df,
                             # family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  
  anova_full_tr_0 <- anova(FullModel, TreatmentNullModel)
  anova_full_day_0 <- anova(FullModel, DayNullModel)
  #anova_full_interact_0 <- anova(FullModel, InteractNullModel)
  
  return(list(
    full_vs_treat_null = anova_full_tr_0, 
    full_vs_day_null = anova_full_day_0))
    #full_vs_interact_null = anova_full_interact_0))
}

lmers_simple <- function(df, var, data_type) {
  # make positive
  if (mean(df[[var]], na.rm = TRUE) < 0) {
    df[[var]] = -df[[var]] 
  }
  #remove rows with nan values
  df <- df[!is.na(df[[var]]), ]
  
  formula_full <- as.formula(paste(var, "~ treatment_r * hrs_incubation + (1|OP) + (1|OP:", data_type, ")"))
  formula_tr_0 <- as.formula(paste(var, "~ (1|OP) + (1|OP:", data_type, ")"))
  formula_hrs_0 <- as.formula(paste(var, "~ treatment_r + (1|OP) + (1|OP:", data_type, ")"))
  formula_interact_0 <- as.formula(paste(var, "~ treatment_r + hrs_incubation+ (1|OP) + (1|OP:", data_type, ")"))

  FullModel <- lmer(formula_full, data = df, REML = FALSE)
  TreatmentNullModel <-lmer(formula_tr_0, data = df, REML = FALSE)
  HrsNullModel <-lmer(formula_hrs_0, data = df, REML = FALSE)
  InteractNullModel <-lmer(formula_interact_0, data = df, REML = FALSE)
  
  anova_full_tr_0 <- anova(FullModel, TreatmentNullModel)
  anova_full_hrs_0 <- anova(FullModel, HrsNullModel)
  anova_full_interact_0 <- anova(FullModel, InteractNullModel)

  return(list(full_vs_treat_null = anova_full_tr_0,
              full_vs_hrs_null = anova_full_hrs_0,
              full_vs_interact_null = anova_full_interact_0))
  }

glmers_simple <- function(df, var, data_type) {
  # make positive
  if (mean(df[[var]], na.rm = TRUE) < 0) {
    df[[var]] = -df[[var]] 
  }
  # remove rows with nans
  df <- df[!is.na(df[[var]]), ]
  
  formula_full <- as.formula(paste(var, "~ treatment_r * hrs_incubation + (1|OP) + (1|OP:", data_type, ")"))
  formula_tr_0 <- as.formula(paste(var, "~ (1|OP) + (1|OP:", data_type, ")"))
  formula_hrs_0 <- as.formula(paste(var, "~ treatment_r + (1|OP) + (1|OP:", data_type, ")"))
  formula_interact_0 <- as.formula(paste(var, "~ treatment_r + hrs_incubation+ (1|OP) + (1|OP:", data_type, ")"))
  
  FullModel <- glmer(formula_full,
                     data = df,
                     family = Gamma(link = "log"), 
                     control = glmerControl(nAGQ0initStep = 0))
  TreatmentNullModel <- glmer(formula_tr_0,
                              data = df,
                              family = Gamma(link = "log"),
                              control = glmerControl(nAGQ0initStep = 0))
  HrsNullModel <- glmer(formula_hrs_0,
                        data = df,
                        family = Gamma(link = "log"),
                        control = glmerControl(nAGQ0initStep = 0))
  InteractNullModel <- glmer(formula_interact_0,
                             data = df,
                             family = Gamma(link = "log"),
                             control = glmerControl(nAGQ0initStep = 0))
  
  
  anova_full_tr_0 <- anova(FullModel, TreatmentNullModel)
  anova_full_hrs_0 <- anova(FullModel, HrsNullModel)
  anova_full_interact_0 <- anova(FullModel, InteractNullModel)
  
  return(list(full_vs_treat_null = anova_full_tr_0,
              full_vs_hrs_null = anova_full_hrs_0,
              full_vs_interact_null = anova_full_interact_0))
}

marginal_effects_simple <- function(df, var, distr, data_type) {
  
  if (mean(df[[var]], na.rm = TRUE) < 0) {
    df[[var]] = -df[[var]] 
  }
  
  # remove rows with nans
  df <- df[!is.na(df[[var]]), ]
  
  formula_full <- as.formula(paste(var, "~ treatment_r + (1|OP) + (1|OP:", data_type, ")"))
  if (distr == 'gaussian'){
    FullModel <- lmer(formula_full,
                      data = df,
                      control =  lmerControl(optimizer = "Nelder_Mead"))
  } else if (distr == 'gamma'){
    FullModel <- glmer(formula_full,
                       data = df, 
                       family = Gamma(link = "log"), 
                       control = glmerControl(nAGQ0initStep = 0))
  }
  
  # Post-hoc analysis with Estimated Marginal Means ####
  MarginalEffects <- emmeans(object = FullModel, 
                             specs =  ~ treatment_r,
                             type="response", tran = "log")
  
  #plot(MarginalEffects)
  ContrastMatrix <- diag(3)
  gr_tr <- MarginalEffects@grid$treatment_r
  
  ### Contrast and adjustment for multiple testing (Benjamini-Hochberg) ####
  MarginalEffectTest <- contrast(object = MarginalEffects, type = "response",
                                 method = list("(BL) / (Ctrl)" = ContrastMatrix[gr_tr=="0",] - ContrastMatrix[gr_tr=="0",],
                                               # "(Ctrl BL) / (Ctrl mid)" = ContrastMatrix[gr_tr=="0",] - ContrastMatrix[gr_tr=="mid Ctrl",],
                                               # "(Ctrl BL) / (Ctrl long)" = ContrastMatrix[gr_tr=="0",] - ContrastMatrix[gr_tr=="long Ctrl",],
                                               "(BL) / (high K)" = ContrastMatrix[gr_tr=="0",] - ContrastMatrix[gr_tr=="1",],
                                               # "(high K BL) / (high K mid)" = ContrastMatrix[gr_tr=="0",] - ContrastMatrix[gr_tr=="mid high K",],
                                               "(Ctrl D2) / (high K D2)" = ContrastMatrix[gr_tr=="0",] - ContrastMatrix[gr_tr=="1",]),
                                 adjust="BH")
  # print(MarginalEffectTest)
  
  ### Plot marginal effects ####
  MarginalPlot <- plot(MarginalEffectTest, colors = "grey10") +
    labs(title = paste(var, "Effects of Treatment and Day"))+
    ylab(label = "")+
    geom_vline(xintercept = 1, linetype="dashed")+
    xlab(label = "Ratio")+
    theme_classic()
  
  return(list(
    effects = summary(MarginalEffects),
    test = summary(MarginalEffectTest)))
}

marginal_effects_day <- function(df, var, distr, data_type) {
  
  if (mean(df[[var]], na.rm = TRUE) < 0) {
    df[[var]] = -df[[var]] 
  }
  
  # remove rows with nans
  df <- df[!is.na(df[[var]]), ]
  
  formula_full <- as.formula(paste(var, "~ treatment_r + day + (1|OP) + (1|OP:", data_type, ")"))
  if (distr == 'gaussian'){
    FullModel <- lmer(formula_full,
                      data = df,
                      control =  lmerControl(optimizer = "Nelder_Mead"))
  } else if (distr == 'gamma'){
    FullModel <- glmer(formula_full,
                       data = df, 
                       family = Gamma(link = "log"), 
                       control = glmerControl(nAGQ0initStep = 0))
  }
  
  # Post-hoc analysis with Estimated Marginal Means ####
  MarginalEffects <- emmeans(object = FullModel, 
                             specs =  ~ treatment_r + day,
                             type="response", tran = "log")
  
  #plot(MarginalEffects)
  
  ContrastMatrix <- diag(4)
  gr_day <- MarginalEffects@grid$day
  gr_tr <- MarginalEffects@grid$treatment_r
  
  ### Contrast and adjustment for multiple testing (Benjamini-Hochberg) ####
  MarginalEffectTest <- contrast(object = MarginalEffects, type = "response",
                                 method = list("(Ctrl D1) / (Ctrl D2)" = ContrastMatrix[gr_tr==0&gr_day==0,] - ContrastMatrix[gr_tr==0&gr_day==1,],
                                               "(high K D1) / (high K D2)" = ContrastMatrix[gr_tr==0&gr_day==0,] - ContrastMatrix[gr_tr==1&gr_day==1,],
                                               "(Ctrl D2) / (high K D2)" = ContrastMatrix[gr_tr==0&gr_day==1,] - ContrastMatrix[gr_tr==1&gr_day==1,]),
                                 adjust="BH")
  # print(MarginalEffectTest)
  
  ### Plot marginal effects ####
  MarginalPlot <- plot(MarginalEffectTest, colors = "grey10") +
    labs(title = paste(var, "Effects of Treatment and Day"))+
    ylab(label = "")+
    geom_vline(xintercept = 1, linetype="dashed")+
    xlab(label = "Ratio")+
    theme_classic()
  MarginalPlot
  
  return(list(
    emmeans_effects = summary(MarginalEffects),
    contrasts_test = summary(MarginalEffectTest)))
}

# functions for the analysis of IFF and num_aps
lmers_firing <- function(df, var, id2, data_type) {
  # id2  = day or hrs_incubation
  # var = num_aps or IFF
  # make positive
  # if (mean(df[[var]], na.rm = TRUE) < 0) {
  #   df[[var]] = -df[[var]] 
  # }
  #remove rows with nan values
  # df <- df[!is.na(df$VAL), ]
  
  formula_full <- as.formula(paste("VAL~ treatment_r + inj_current + ", id2," + (1|OP) + (1|OP:", data_type, ")"))
  formula_tr_day_0 <- as.formula(paste("VAL ~ inj_current + (1|OP) + (1|OP:", data_type, ")"))
  formula_tr_0 <- as.formula(paste("VAL ~ inj_current +" ,id2, "+ (1|OP) + (1|OP:", data_type, ")"))
  formula_day_0 <- as.formula(paste(" VAL ~ inj_current + treatment_r + (1|OP) + (1|OP:", data_type, ")"))

  FullModel <- lmer(formula_full, data = df, REML = FALSE)
  TreatmentNullModel <- lmer(formula_tr_0, data = df, REML = FALSE)
  DayNullModel <- lmer(formula_day_0, data = df, REML = FALSE)
  NullModel_tr_day <- lmer(formula_tr_day_0, data = df, REML = FALSE)
  
  anova_full_0_tr_day <- anova(FullModel, NullModel_tr_day)
  anova_full_tr_0 <- anova(FullModel, TreatmentNullModel)
  anova_full_day_0 <- anova(FullModel, DayNullModel)
  
  return(list(
    full_vs_tr_day_null = anova_full_0_tr_day, # only injection explains the data as well as when other params included
    full_vs_treat_null = anova_full_tr_0,
    full_vs_day_null = anova_full_day_0))
}

marginal_effects_day_firing <- function(df, var, data_type){
  
  formula_full <- as.formula(paste(var, "~ inj_current + day + treatment_r + (1|OP) + (1|OP:", data_type, ")"))
  
  FullModel <- lmer(formula_full,
                    data = df,
                    control =  lmerControl(optimizer = "Nelder_Mead"))
  
  # Post-hoc analysis with Estimated Marginal Means ####
  MarginalEffects <- emmeans(object = FullModel, 
                             specs =  ~ day + treatment_r + inj_current,
                             # uncomment in case comparison at specifc inj_values
                             at = list(inj_current = c(200, 400, 600)),
                             type="response", tran = "log",
                             pbkrtest.limit = 5000)
  
  
  ContrastMatrix <- diag(12)
  gr_day <- MarginalEffects@grid$day
  gr_tr <- MarginalEffects@grid$treatment_r
  gr_inj <- MarginalEffects@grid$inj_current
  
  ### Contrast and adjustment for multiple testing (Benjamini-Hochberg) ####
  MarginalEffectTest <- contrast(object = MarginalEffects, type = "response",
                                 method = list("200 pA (Ctrl D1) / (Ctrl D2)" = ContrastMatrix[gr_tr==0&gr_day==0&gr_inj=='200',] - ContrastMatrix[gr_tr==0&gr_day==1&gr_inj=='200',],
                                               "400 pA (Ctrl D1) / (Ctrl D2)" = ContrastMatrix[gr_tr==0&gr_day==0&gr_inj=='400',] - ContrastMatrix[gr_tr==0&gr_day==1&gr_inj=='400',],
                                               "600 pA (Ctrl D1) / (Ctrl D2)" = ContrastMatrix[gr_tr==0&gr_day==0&gr_inj=='600',] - ContrastMatrix[gr_tr==0&gr_day==1&gr_inj=='600',],
                                               "200 pA (high K D1) / (high K D2)" = ContrastMatrix[gr_tr==0&gr_day==0&gr_inj=='200',] - ContrastMatrix[gr_tr==1&gr_day==1&gr_inj=='200',],
                                               "400 pA (high K D1) / (high K D2)" = ContrastMatrix[gr_tr==0&gr_day==0&gr_inj=='400',] - ContrastMatrix[gr_tr==1&gr_day==1&gr_inj=='400',],
                                               "600 pA (high K D1) / (high K D2)" = ContrastMatrix[gr_tr==0&gr_day==0&gr_inj=='600',] - ContrastMatrix[gr_tr==1&gr_day==1&gr_inj=='600',]),
                                               # "600 pA (Ctrl D2) / (high K D2)" = ContrastMatrix[gr_tr=="Ctrl"&gr_day=='D2'&gr_inj=='600',] - ContrastMatrix[gr_tr=="high K"&gr_day=='D2'&gr_inj=='600',]),
                                 adjust="BH")
  
  return(list(
    effects = summary(MarginalEffects),
    test = summary(MarginalEffectTest)))
}

marginal_effects_day_firing <- function(df, var, data_type){
  
  formula_full <- as.formula(paste(var, "~ inj_current + day + treatment_r + (1|OP) + (1|OP:", data_type, ")"))
  
  FullModel <- lmer(formula_full,
                    data = df,
                    control =  lmerControl(optimizer = "Nelder_Mead"))
  
  # Post-hoc analysis with Estimated Marginal Means ####
  MarginalEffects <- emmeans(object = FullModel, 
                             specs =  ~ day + treatment_r + inj_current,
                             # uncomment in case comparison at specifc inj_values
                             at = list(inj_current = c(200, 400, 600)),
                             type="response", tran = "log",
                             pbkrtest.limit = 5000)
  
  
  ContrastMatrix <- diag(12)
  gr_day <- MarginalEffects@grid$day
  gr_tr <- MarginalEffects@grid$treatment_r
  gr_inj <- MarginalEffects@grid$inj_current
  
  ### Contrast and adjustment for multiple testing (Benjamini-Hochberg) ####
  MarginalEffectTest <- contrast(object = MarginalEffects, type = "response",
                                 method = list("200 pA (Ctrl D1) / (Ctrl D2)" = ContrastMatrix[gr_tr==0&gr_day==0&gr_inj=='200',] - ContrastMatrix[gr_tr==0&gr_day==1&gr_inj=='200',],
                                               "400 pA (Ctrl D1) / (Ctrl D2)" = ContrastMatrix[gr_tr==0&gr_day==0&gr_inj=='400',] - ContrastMatrix[gr_tr==0&gr_day==1&gr_inj=='400',],
                                               "600 pA (Ctrl D1) / (Ctrl D2)" = ContrastMatrix[gr_tr==0&gr_day==0&gr_inj=='600',] - ContrastMatrix[gr_tr==0&gr_day==1&gr_inj=='600',],
                                               "200 pA (high K D1) / (high K D2)" = ContrastMatrix[gr_tr==0&gr_day==0&gr_inj=='200',] - ContrastMatrix[gr_tr==1&gr_day==1&gr_inj=='200',],
                                               "400 pA (high K D1) / (high K D2)" = ContrastMatrix[gr_tr==0&gr_day==0&gr_inj=='400',] - ContrastMatrix[gr_tr==1&gr_day==1&gr_inj=='400',],
                                               "600 pA (high K D1) / (high K D2)" = ContrastMatrix[gr_tr==0&gr_day==0&gr_inj=='600',] - ContrastMatrix[gr_tr==1&gr_day==1&gr_inj=='600',]),
                                 # "600 pA (Ctrl D2) / (high K D2)" = ContrastMatrix[gr_tr=="Ctrl"&gr_day=='D2'&gr_inj=='600',] - ContrastMatrix[gr_tr=="high K"&gr_day=='D2'&gr_inj=='600',]),
                                 adjust="BH")
  
  return(list(
    effects = summary(MarginalEffects),
    test = summary(MarginalEffectTest)))
}

lmers_extended <- function(df, var, id2, data_type) {
  # make positive
  if (mean(df[[var]], na.rm = TRUE) < 0) {
    df[[var]] = -df[[var]] 
  }
  #remove rows with nan values
  df <- df[!is.na(df[[var]]), ]
  
  formula_full <- as.formula(paste(var, "~ treatment_r +", id2, "+ patient_age + (1|OP) + (1|OP:", data_type, ")"))
  formula_tr_0 <- as.formula(paste(var, "~ ", id2 ,"+ patient_age + (1|OP) + (1|OP:", data_type, ")"))
  formula_hrs_0 <- as.formula(paste(var, "~ treatment_r + patient_age + (1|OP) + (1|OP:", data_type, ")"))
  formula_age_0 <- as.formula(paste(var, "~ treatment_r +", id2, "  + (1|OP) + (1|OP:", data_type, ")"))
  #formula_sex_0 <- as.formula(paste(var, "~ treatment_r +", id2, "+ patient_age + (1|OP) + (1|OP:", data_type, ")"))
  
  FullModel <-    lmer(formula_full, data = df, REML = FALSE)
  TrtNullModel <- lmer(formula_tr_0, data = df, REML = FALSE)
  IncNullModel <- lmer(formula_hrs_0, data = df, REML = FALSE)
  AgeNullModel <- lmer(formula_age_0, data = df, REML = FALSE)
  #SexNullModel <- lmer(formula_sex_0, data = df, REML = FALSE)
  
  anova_full_tr_0 <- anova(FullModel, TrtNullModel)
  anova_full_inc_0 <- anova(FullModel, IncNullModel)
  anova_full_age_0 <- anova(FullModel, AgeNullModel)
  #anova_full_sex_0 <- anova(FullModel, SexNullModel)
  
  return(list(
    full_vs_treat_null = anova_full_tr_0, 
    full_vs_inc_null = anova_full_inc_0,
    full_vs_age_null = anova_full_age_0))
    #full_vs_sex_null = anova_full_sex_0))
}

glmers_extended <- function(df, var, id2, data_type) {
  # make positive
  if (mean(df[[var]], na.rm = TRUE) < 0) {
    df[[var]] = -df[[var]] 
  }
  # remove rows with nans
  df <- df[!is.na(df[[var]]), ]
  
  formula_full <- as.formula(paste(var, "~ treatment_r +", id2, "+ patient_age + (1|OP) + (1|OP:", data_type, ")"))
  formula_tr_0 <- as.formula(paste(var, "~ ", id2 ,"+ patient_age + (1|OP) + (1|OP:", data_type, ")"))
  formula_hrs_0 <- as.formula(paste(var, "~ treatment_r + patient_age  + (1|OP) + (1|OP:", data_type, ")"))
  formula_age_0 <- as.formula(paste(var, "~ treatment_r +", id2, " + (1|OP) + (1|OP:", data_type, ")"))
  #formula_sex_0 <- as.formula(paste(var, "~ treatment_r +", id2, "+ patient_age + (1|OP) + (1|OP:", data_type, ")"))
  
  FullModel <- glmer(formula_full, data = df,
                     family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  TrtNullModel <- glmer(formula_tr_0, data = df,
                        family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  IncNullModel <- glmer(formula_hrs_0, data = df,
                        family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  AgeNullModel <- glmer(formula_age_0, data = df,
                        family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  #SexNullModel <- glmer(formula_sex_0, data = df,
   #                    family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  
  anova_full_tr_0 <- anova(FullModel, TrtNullModel)
  anova_full_inc_0 <- anova(FullModel, IncNullModel)
  anova_full_age_0 <- anova(FullModel, AgeNullModel)
  #anova_full_sex_0 <- anova(FullModel, SexNullModel)
  
  return(list(
    full_vs_treat_null = anova_full_tr_0, 
    full_vs_inc_null = anova_full_inc_0,
    full_vs_age_null = anova_full_age_0))
    #full_vs_sex_null = anova_full_sex_0))

}

# extended_2 - have the interaction term patient_age * treatment_r
lmers_extended_2 <- function(df, var, id2, data_type) {
  # make positive
  if (mean(df[[var]], na.rm = TRUE) < 0) {
    df[[var]] = -df[[var]] 
  }
  #remove rows with nan values
  df <- df[!is.na(df[[var]]), ]
  
  formula_full <- as.formula(paste(var, "~ treatment_r * patient_age +", id2, " + (1|OP) + (1|OP:", data_type, ")"))
  formula_tr_0 <- as.formula(paste(var, "~  ", id2 ,"+ patient_age + (1|OP) + (1|OP:", data_type, ")"))
  formula_hrs_0 <- as.formula(paste(var, "~  + treatment_r * patient_age  + (1|OP) + (1|OP:", data_type, ")"))
  formula_age_0 <- as.formula(paste(var, "~  + treatment_r +", id2, " + (1|OP) + (1|OP:", data_type, ")"))
  formula_interaction_0 <- as.formula(paste(var, "~ treatment_r + patient_age +", id2, " + (1|OP) + (1|OP:", data_type, ")"))

  FullModel <-    lmer(formula_full, data = df, REML = FALSE)
  TrtNullModel <- lmer(formula_tr_0, data = df, REML = FALSE)
  IncNullModel <- lmer(formula_hrs_0, data = df, REML = FALSE)
  AgeNullModel <- lmer(formula_age_0, data = df, REML = FALSE)
  InteractNullModel <- lmer(formula_interaction_0, data = df, REML = FALSE)
  
  anova_full_tr_0 <- anova(FullModel, TrtNullModel)
  anova_full_inc_0 <- anova(FullModel, IncNullModel)
  anova_full_age_0 <- anova(FullModel, AgeNullModel)
  anova_full_interact_0 <- anova(FullModel, InteractNullModel)
  
  return(list(
    full_vs_treat_null = anova_full_tr_0, 
    full_vs_inc_null = anova_full_inc_0,
    full_vs_age_null = anova_full_age_0,
    full_vs_interact_null = anova_full_interact_0))
}

lmers_extended_hrs_only <- function(df, var, data_type) {
  # make positive
  if (mean(df[[var]], na.rm = TRUE) < 0) {
    df[[var]] = -df[[var]] 
  }
  #remove rows with nan values
  df <- df[!is.na(df[[var]]), ]
  
  formula_full <- as.formula(paste(var, "~  hrs_after_OP + (1|OP) + (1|OP:", data_type, ")"))
  # formula_tr_0 <- as.formula(paste(var, "~ hrs_after_OP + (1|OP) + (1|OP:", data_type, ")"))
  formula_hrs_0 <- as.formula(paste(var, "~ (1|OP) + (1|OP:", data_type, ")"))

  FullModel <-    lmer(formula_full, data = df, REML = FALSE)
  # TrtNullModel <- lmer(formula_tr_0, data = df, REML = FALSE)
  IncNullModel <- lmer(formula_hrs_0, data = df, REML = FALSE)
  
  #anova_full_tr_0 <- anova(FullModel, TrtNullModel)
  anova_full_inc_0 <- anova(FullModel, IncNullModel)
  
  return(list(full_vs_inc_null = anova_full_inc_0)) 
    # full_vs_treat_null = anova_full_tr_0))
}

glmers_extended_2 <- function(df, var, id2, data_type) {
  # make positive
  if (mean(df[[var]], na.rm = TRUE) < 0) {
    df[[var]] = -df[[var]] 
  }
  # remove rows with nans
  df <- df[!is.na(df[[var]]), ]
  
  formula_full <- as.formula(paste(var, "~ treatment_r * patient_age +", id2, " + (1|OP) + (1|OP:", data_type, ")"))
  formula_tr_0 <- as.formula(paste(var, "~ ", id2 ,"+ patient_age + (1|OP) + (1|OP:", data_type, ")"))
  formula_hrs_0 <- as.formula(paste(var, "~ treatment_r * patient_age  + (1|OP) + (1|OP:", data_type, ")"))
  formula_age_0 <- as.formula(paste(var, "~ treatment_r +", id2, " + (1|OP) + (1|OP:", data_type, ")"))
  formula_interaction_0 <- as.formula(paste(var, "~ treatment_r + patient_age +", id2, " + (1|OP) + (1|OP:", data_type, ")"))
  
  FullModel <- glmer(formula_full, data = df,
                     family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  TrtNullModel <- glmer(formula_tr_0, data = df,
                        family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  IncNullModel <- glmer(formula_hrs_0, data = df,
                        family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  AgeNullModel <- glmer(formula_age_0, data = df,
                        family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  InteractNullModel <- glmer(formula_interaction_0, data = df,
                             family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  
  anova_full_tr_0 <- anova(FullModel, TrtNullModel)
  anova_full_inc_0 <- anova(FullModel, IncNullModel)
  anova_full_age_0 <- anova(FullModel, AgeNullModel)
  anova_full_interact_0 <- anova(FullModel, InteractNullModel)
  
  return(list(
    full_vs_treat_null = anova_full_tr_0, 
    full_vs_inc_null = anova_full_inc_0,
    full_vs_age_null = anova_full_age_0,
    full_vs_interact_null = anova_full_interact_0))
  
}

marginal_effects_hrs <- function(df, var, distr, data_type) {

  if (mean(df[[var]], na.rm = TRUE) < 0) {
    df[[var]] = -df[[var]]
  }

  # remove rows with nans
  df <- df[!is.na(df[[var]]), ]

  formula_full <- as.formula(paste(var, "~ treatment_r + hrs_incubation + (1|OP) + (1|OP:", data_type, ")"))
  if (distr == 'gaussian'){
    FullModel <- lmer(formula_full,
                       data = df,
                       control =  lmerControl(optimizer = "Nelder_Mead"))
  } else if (distr == 'gamma'){
    FullModel <- glmer(formula_full,
                       data = df,
                       family = Gamma(link = "log"),
                       control = glmerControl(nAGQ0initStep = 0))
  }


  # Post-hoc analysis with Estimated Marginal Means ####
  ## dependent var per treatment for time (0, 18, 21, 24)
  MarginalEffects <- emmeans(object = FullModel,
                             specs =  ~ treatment_r + hrs_incubation,
                             at = list(hrs_incubation = c(0, 18, 21, 24)),
                             type="response", tran = "log")

  #plot(MarginalEffects)

  ContrastMatrix <- diag(12)
  gr_inc <- MarginalEffects@grid$hrs_incubation
  gr_tr <- MarginalEffects@grid$treatment_r

  ### Contrast and adjustment for multiple testing (Benjamini-Hochberg) ####
  MarginalEffectTest <- contrast(object = MarginalEffects, type = "response",
                                 method = list("(Ctrl 0) / (Ctrl 18)" = ContrastMatrix[gr_tr=="BL"&gr_inc==0,] - ContrastMatrix[gr_tr=="CTR"&gr_inc==18,],
                                               "(Ctrl 0) / (Ctrl 21)" = ContrastMatrix[gr_tr=="BL"&gr_inc==0,] - ContrastMatrix[gr_tr=="CTR"&gr_inc==21,],
                                               "(Ctrl 0) / (Ctrl 24)" = ContrastMatrix[gr_tr=="BL"&gr_inc==0,] - ContrastMatrix[gr_tr=="CTR"&gr_inc==24,],
                                               "(high K 0) / (high K 18)" = ContrastMatrix[gr_tr=="BL"&gr_inc==0,] - ContrastMatrix[gr_tr=="HiK"&gr_inc==18,],
                                               "(high K 0) / (high K 21)" = ContrastMatrix[gr_tr=="BL"&gr_inc==0,] - ContrastMatrix[gr_tr=="HiK"&gr_inc==21,],
                                               "(high K 0) / (high K 24)" = ContrastMatrix[gr_tr=="BL"&gr_inc==0,] - ContrastMatrix[gr_tr=="HiK"&gr_inc==24,]),
                                 adjust="BH")
  # print(MarginalEffectTest)

  # ### Plot marginal effects ####
  # MarginalPlot <- plot(MarginalEffectTest, colors = "grey10") +
  #    labs(title = paste(var, " Treatment and Hours of incubation Effects"))+
  #    ylab(label = "")+
  #    geom_vline(xintercept = 1, linetype="dashed")+
  #    xlab(label = "Ratio")+
  #    theme_classic()

  return(list(
    effects = summary(MarginalEffects),
    test = summary(MarginalEffectTest)))
}

lmers_extended_no_patient_age <- function(df, var, id2, data_type) {
  # make positive
  if (mean(df[[var]], na.rm = TRUE) < 0) {
    df[[var]] = -df[[var]] 
  }
  #remove rows with nan values
  df <- df[!is.na(df[[var]]), ]

  formula_full <- as.formula(paste(var, "~ treatment_r +", id2, " + (1|OP) + (1|OP:", data_type, ")"))
  formula_tr_0 <- as.formula(paste(var, "~  ", id2 ," + (1|OP) + (1|OP:", data_type, ")"))
  formula_hrs_0 <- as.formula(paste(var, "~  + treatment_r  + (1|OP) + (1|OP:", data_type, ")"))
  #formula_interaction_0 <- as.formula(paste(var, "~ treatment_r +", id2, " + (1|OP) + (1|OP:", data_type, ")"))

  FullModel <-    lmer(formula_full, data = df, REML = FALSE)
  TrtNullModel <- lmer(formula_tr_0, data = df, REML = FALSE)
  IncNullModel <- lmer(formula_hrs_0, data = df, REML = FALSE)
  #InteractNullModel <- lmer(formula_interaction_0, data = df, REML = FALSE)
  
  anova_full_tr_0 <- anova(FullModel, TrtNullModel)
  anova_full_inc_0 <- anova(FullModel, IncNullModel)
  #anova_full_interact_0 <- anova(FullModel, InteractNullModel)
  
  return(list(
    full_vs_treat_null = anova_full_tr_0, 
    full_vs_inc_null = anova_full_inc_0))
    # full_vs_interact_null = anova_full_interact_0))
}

glmers_extended_no_patient_age <- function(df, var, id2, data_type) {
  # make positive
  if (mean(df[[var]], na.rm = TRUE) < 0) {
    df[[var]] = -df[[var]] 
  }
  # remove rows with nans
  df <- df[!is.na(df[[var]]), ]
  
  formula_full <- as.formula(paste(var, "~ treatment_r * ", id2, " + (1|OP) + (1|OP:", data_type, ")"))
  formula_tr_0 <- as.formula(paste(var, "~ ", id2 ," + (1|OP) + (1|OP:", data_type, ")"))
  formula_hrs_0 <- as.formula(paste(var, "~ treatment_r + (1|OP) + (1|OP:", data_type, ")"))
  # formula_interaction_0 <- as.formula(paste(var, "~ treatment_r  + ", id2, " + (1|OP) + (1|OP:", data_type, ")"))
  
  FullModel <- glmer(formula_full, data = df,
                     family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  TrtNullModel <- glmer(formula_tr_0, data = df,
                        family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  IncNullModel <- glmer(formula_hrs_0, data = df,
                        family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  # InteractNullModel <- glmer(formula_interaction_0, data = df,
  #                            family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  # 
  anova_full_tr_0 <- anova(FullModel, TrtNullModel)
  anova_full_inc_0 <- anova(FullModel, IncNullModel)
  # anova_full_interact_0 <- anova(FullModel, InteractNullModel)
  
  return(list(
    full_vs_treat_null = anova_full_tr_0, 
    full_vs_inc_null = anova_full_inc_0))
   #  full_vs_interact_null = anova_full_interact_0))
}

# functions for the analysis of IFF and num_aps
lmers_firing_extended <- function(df, var, data_type) {

  # if (mean(df[[var]], na.rm = TRUE) < 0) {
  #   df[[var]] = -df[[var]] 
  # }
  #remove rows with nan values
  # df <- df[!is.na(df[[var]]), ]
  
  formula_full <- as.formula(paste("VAL ~ inj_current + treatment_r * patient_age + hrs_after_OP + (1|OP) + (1|OP:", data_type, ")"))
  formula_tr_0 <- as.formula(paste("VAL ~ inj_current + hrs_after_OP + patient_age + (1|OP) + (1|OP:", data_type, ")"))
  formula_hrs_0 <- as.formula(paste("VAL ~ inj_current + treatment_r * patient_age  + (1|OP) + (1|OP:", data_type, ")"))
  formula_age_0 <- as.formula(paste("VAL ~ inj_current + treatment_r + hrs_after_OP + (1|OP) + (1|OP:", data_type, ")"))
  formula_interaction_0 <- as.formula(paste("VAL ~ inj_current + treatment_r + patient_age + hrs_after_OP + (1|OP) + (1|OP:", data_type, ")"))
  
  FullModel <- lmer(formula_full, data = df, REML = FALSE)
  TreatmentNullModel <- lmer(formula_tr_0, data = df, REML = FALSE)
  HrsNullModel <- lmer(formula_hrs_0, data = df, REML = FALSE)
  AgeNullModel <- lmer(formula_age_0, data = df, REML = FALSE)
  InteractNullModel <- lmer(formula_interaction_0, data = df, REML = FALSE)
  
  anova_full_tr_0 <- anova(FullModel, TreatmentNullModel)
  anova_full_hrs_0 <- anova(FullModel, HrsNullModel)
  anova_full_age_0 <- anova(FullModel, AgeNullModel)
  anova_full_interact_0 <- anova(FullModel, InteractNullModel)
  
  return(list(
    full_vs_treat_null = anova_full_tr_0,
    full_vs_day_null = anova_full_hrs_0,
    full_vs_age_0 = anova_full_age_0, # only injection explains the data as well as when other params included
    full_vs_interact_null = anova_full_interact_0))
}

lmers_ais <- function(df) {
  
  FullModel <- lmer(Length ~  treatment_r * patient_age + hrs_incubation + (1|OP) + (1|OP:Slice),
                     data = df, REML = FALSE)
  TreatmentNullModel <- lmer(Length ~  patient_age + hrs_incubation + (1|OP) + (1|OP:Slice),
                              data = df, REML = FALSE)
  HrsNullModel <- lmer(Length ~ treatment_r * patient_age + (1|OP) + (1|OP:Slice),  
                        data = df, REML = FALSE)
  AgeNullModel <- lmer(Length ~ treatment_r + hrs_incubation  + (1|OP) + (1|OP:Slice),
                        data = df, REML = FALSE)
  InteractNullModel <- lmer(Length ~ hrs_incubation + treatment_r + patient_age  + (1|OP) + (1|OP:Slice),
                             data = df, REML = FALSE)
  
  anova_full_tr_0 <- anova(FullModel, TreatmentNullModel)
  anova_full_hrs_0 <- anova(FullModel, HrsNullModel)
  #anova_full_soma_0 <- anova(FullModel, SomaNullModel)
  anova_full_age_0 <- anova(FullModel, AgeNullModel)
  anova_full_interact_0 <- anova(FullModel, InteractNullModel)
  
  return(list(
    full_vs_treat_null = anova_full_tr_0, 
    full_vs_hrs_null = anova_full_hrs_0,
    #full_vs_soma_null = anova_full_soma_0,
    full_vs_age_null = anova_full_age_0,
    full_vs_interact_null = anova_full_interact_0))
}

lmers_firing_extended_no_age <- function(df, var, data_type) {
  
  formula_full <- as.formula(paste("VAL ~ inj_current + treatment_r  + hrs_after_OP + (1|OP) + (1|OP:", data_type, ")"))
  formula_tr_0 <- as.formula(paste("VAL ~ inj_current + hrs_after_OP  + (1|OP) + (1|OP:", data_type, ")"))
  formula_hrs_0 <- as.formula(paste("VAL ~ inj_current + treatment_r   + (1|OP) + (1|OP:", data_type, ")"))

  FullModel <- lmer(formula_full, data = df, REML = FALSE)
  TreatmentNullModel <- lmer(formula_tr_0, data = df, REML = FALSE)
  HrsNullModel <- lmer(formula_hrs_0, data = df, REML = FALSE)

  anova_full_tr_0 <- anova(FullModel, TreatmentNullModel)
  anova_full_hrs_0 <- anova(FullModel, HrsNullModel)

  return(list(
    full_vs_treat_null = anova_full_tr_0,
    full_vs_day_null = anova_full_hrs_0))
}

lmers_ais <- function(df) {
  
  FullModel <- lmer(Length ~  treatment_r * patient_age + hrs_incubation + (1|OP) + (1|OP:Slice),
                    data = df, REML = FALSE)
  TreatmentNullModel <- lmer(Length ~  patient_age + hrs_incubation + (1|OP) + (1|OP:Slice),
                             data = df, REML = FALSE)
  HrsNullModel <- lmer(Length ~ treatment_r * patient_age + (1|OP) + (1|OP:Slice),  
                       data = df, REML = FALSE)
  AgeNullModel <- lmer(Length ~ treatment_r + hrs_incubation  + (1|OP) + (1|OP:Slice),
                       data = df, REML = FALSE)
  InteractNullModel <- lmer(Length ~ hrs_incubation + treatment_r + patient_age  + (1|OP) + (1|OP:Slice),
                            data = df, REML = FALSE)
  
  anova_full_tr_0 <- anova(FullModel, TreatmentNullModel)
  anova_full_hrs_0 <- anova(FullModel, HrsNullModel)
  #anova_full_soma_0 <- anova(FullModel, SomaNullModel)
  anova_full_age_0 <- anova(FullModel, AgeNullModel)
  anova_full_interact_0 <- anova(FullModel, InteractNullModel)
  
  return(list(
    full_vs_treat_null = anova_full_tr_0, 
    full_vs_hrs_null = anova_full_hrs_0,
    #full_vs_soma_null = anova_full_soma_0,
    full_vs_age_null = anova_full_age_0,
    full_vs_interact_null = anova_full_interact_0))
}

glmers_ais <- function(df) {
  # remove rows with nans
  df <- df[!is.na(df[["Length"]]), ]

  FullModel <- glmer(Length ~  treatment_r * patient_age + hrs_incubation + (1|OP) + (1|OP:Slice),
                     data = df, family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  TreatmentNullModel <- glmer(Length ~  patient_age + hrs_incubation + (1|OP) + (1|OP:Slice),
                              data = df, family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  HrsNullModel <- glmer(Length ~ treatment_r * patient_age + (1|OP) + (1|OP:Slice),  
                        data = df, family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  AgeNullModel <- glmer(Length ~ treatment_r + hrs_incubation  + (1|OP) + (1|OP:Slice),
                        data = df, family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  InteractNullModel <- glmer(Length ~ hrs_incubation + treatment_r + patient_age  + (1|OP) + (1|OP:Slice),
                            data = df, family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  
  anova_full_tr_0 <- anova(FullModel, TreatmentNullModel)
  anova_full_hrs_0 <- anova(FullModel, HrsNullModel)
  #anova_full_soma_0 <- anova(FullModel, SomaNullModel)
  anova_full_age_0 <- anova(FullModel, AgeNullModel)
  anova_full_interact_0 <- anova(FullModel, InteractNullModel)
  
  return(list(
    full_vs_treat_null = anova_full_tr_0, 
    full_vs_hrs_null = anova_full_hrs_0,
    #full_vs_soma_null = anova_full_soma_0,
    full_vs_age_null = anova_full_age_0,
   full_vs_interact_null = anova_full_interact_0))
}

glmers_mea <- function(df) {
  
  # remove rows with nans
  df <- df[!is.na(df[["value"]]), ]
  
  FullModel <- glmer(value ~ condition_r + treatment_r + patient_age  + (1|OP) + (1|OP:slice),
                     data = df, family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  TreatmentNullModel <- glmer(value ~ condition_r + patient_age  + (1|OP) + (1|OP:slice),
                              data = df, family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  CondNullModel <- glmer(value ~ treatment_r + patient_age + (1|OP) + (1|OP:slice),
                        data = df, family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))
  AgeNullModel <- glmer(value ~ treatment_r + condition_r  + (1|OP) + (1|OP:slice),
                        data = df, family = Gamma(link = "log"), control = glmerControl(nAGQ0initStep = 0))

  anova_full_tr_0 <- anova(FullModel, TreatmentNullModel)
  anova_full_cond_0 <- anova(FullModel, CondNullModel)
  anova_full_age_0 <- anova(FullModel, AgeNullModel)
  #anova_full_interact_0 <- anova(FullModel, InteractNullModel)
  
  return(list(
    full_vs_treat_null = anova_full_tr_0, 
    full_vs_cond_null = anova_full_cond_0,
    full_vs_age_null = anova_full_age_0))
    #full_vs_interact_null = anova_full_interact_0))
}

lmers_mea_simple <- function(df) {
  FullModel <- lmer(value ~  treatment_r + (1|OP),
                    data = df, REML = FALSE)
  TreatmentNullModel <- lmer(value ~  (1|OP),
                             data = df, REML = FALSE)
  
  anova_full_tr_0 <- anova(FullModel, TreatmentNullModel)

  return(list(full_vs_treat_null = anova_full_tr_0))
}

lmers_mea_extended <- function(df) {
  FullModel <- lmer(value ~  treatment_r * patient_age  + (1|OP),
                     data = df, REML = FALSE)
  TreatmentNullModel <- lmer(value ~ patient_age  + (1|OP),
                              data = df, REML = FALSE)
  AgeNullModel <- lmer(value ~ treatment_r  + (1|OP),
                        data = df, REML = FALSE)
  InteractNullModel <- lmer(value ~ treatment_r + patient_age + (1|OP),
                       data = df, REML = FALSE)
  
  anova_full_tr_0 <- anova(FullModel, TreatmentNullModel)
  anova_full_age_0 <- anova(FullModel, AgeNullModel)
  anova_full_interact_0 <- anova(FullModel, InteractNullModel)
  
  return(list(
    full_vs_treat_null = anova_full_tr_0, 
    full_vs_age_null = anova_full_age_0,
    full_vs_interact_null = anova_full_interact_0))
}

################################### FINAL #########################

fix_df <- function(df){
  
  df$OP <- factor(df$OP)
  df$day <- factor(df$day)
  df$day_cat <- df$day # keep the categorical data in a column
  df$day <- ifelse(df$day_cat == "D1", 0, 1)
  df$treatment_r <- ifelse(df$treatment == "Ctrl", 0, 1)
  df$treatment_r <- as.numeric(df$treatment_r)
  df$hrs_after_OP <- as.numeric(df$hrs_after_OP)
  # df$patient_age <- as.numeric(df$patient_age)
  setnames(df, colnames(df), make.names(colnames(df))) # make col names without spaces
  
  df$org_patient_age <- df$patient_age
  df$org_hrs_after_OP <- df$hrs_after_OP
  df = df %>%
    mutate(
      hrs_after_OP = scale(hrs_after_OP),
      # patient_age = scale(patient_age)
    )
  
  return(df)
}

fix_df_slice_repatch <- function(df){
  
  df$OP <- factor(df$OP)
  df$day <- factor(df$day)
  df$day_cat <- df$day # keep the categorical data in a column
  df$day <- ifelse(df$day_cat == "D1", 0, 1)
  df$treatment_r <- as.numeric(df$treatment_r)
  df$hrs_after_OP <- as.numeric(df$hrs_after_OP)
  # df$patient_age <- as.numeric(df$patient_age)
  setnames(df, colnames(df), make.names(colnames(df))) # make col names without spaces
  
  df$org_patient_age <- df$patient_age
  df$org_hrs_after_OP <- df$hrs_after_OP
  df = df %>%
    mutate(
      hrs_after_OP = scale(hrs_after_OP),
      # patient_age = scale(patient_age)
    )
  
  return(df)
}

fix_ais_df <- function(df_ais){
  
  # remove rows with nans
  df_ais <- df_ais[!is.na(df_ais[["Length"]]), ]
  
  setnames(df_ais, colnames(df_ais), make.names(colnames(df_ais)))
  df_ais <- df_ais[df_ais$patient_age >= 12, ]
  df_ais <- df_ais[df_ais$area != 'frontal', ]
  df_ais$day <- ifelse(df_ais$hrs_incubation == 0, 0, 1)
  
  df_ais$Length <- as.numeric(df_ais$Length)
  df_ais$treatment_r <- as.numeric(df_ais$treatment_r)
  df_ais$hrs_incubation <- as.numeric(df_ais$hrs_incubation)
  df_ais$patient_age <- as.numeric(df_ais$patient_age)
  
  df_ais$OP <- factor(df_ais$OP)
  df_ais$Soma <- factor(df_ais$Soma)
  df_ais$Slice <- as.factor(df_ais$Slice)
  df_ais$treatment <- as.factor(df_ais$treatment)# Grouping factor for slices
  
  df_ais$org_hrs_incubation <- df_ais$hrs_incubation
  df_ais$org_patient_age <- df_ais$patient_age
  df_ais$org_length <- df_ais$Length

  df_ais = df_ais %>%
    mutate(
      hrs_incubation = scale(hrs_incubation),
      patient_age = scale(patient_age),
      Length = scale(Length)
    )
  return(df_ais)
}

get_melt_df_firing <- function(var_firing, df, ids){
  setnames(df, colnames(df), make.names(colnames(df)))
  # get inj column names
  col_names_inj <- grep(var_firing, colnames(df), value = TRUE) # all relevant column names
  col_names_inj <- col_names_inj[6:length(col_names_inj)] # only positive injs
  
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
  melt_df <- melt(melt_df, id = ids, variable = 'inj_current', value = 'VAL')
  melt_df$inj_current <- as.numeric(as.character(melt_df$inj_current))
  
  # decidede to keep up to 600 pA
  inj_to_keep <- 600
  melt_df <- melt_df[inj_current <= inj_to_keep]
  
  # remove 0 values, aka where no APs fired
  melt_df <- melt_df[melt_df$VAL != 0]
 
  return(melt_df) 
}

df_ext_comparison_groups <- function(df){
  # have treatment as a word
  df$treatment_word <- ifelse(df$treatment_r == 0, 'CTR', 'HiK')
  
  # grouping patient ages for more power
  df <- df %>%
    mutate(age_group = case_when(
      patient_age >= 11 & patient_age <= 20 ~ 15,
      patient_age >= 21 & patient_age <= 30 ~ 25,
      patient_age >= 31 & patient_age <= 40 ~ 35,
      patient_age >= 41 & patient_age <= 50 ~ 45,
      patient_age >= 51 & patient_age <= 80 ~ 60))
  
  df <- df %>%
    mutate(hrs_group = case_when(
      hrs_after_OP >= 5 & hrs_after_OP <= 15.990 ~ 10,
      hrs_after_OP >= 16 & hrs_after_OP <= 25.999 ~ 20, # short
      hrs_after_OP >= 26 & hrs_after_OP <= 35.999 ~ 30, # middle
      hrs_after_OP >= 36 & hrs_after_OP <= 45.999 ~ 40,
      hrs_after_OP >= 46 & hrs_after_OP <= 70 ~ 50)) # long
    
  df$age_group <- as.numeric(df$age_group)
  df$hrs_group <- as.numeric(df$hrs_group)
  return(df)
}

df_ext_comparison_groups_inc_only <- function(df){
  # have treatment as a word
  df$treatment_word <- ifelse(df$treatment_r == 0, 'CTR', 'HiK')
  
  df <- df %>%
    mutate(hrs_group = case_when(
      org_hrs_after_OP >= 5 & org_hrs_after_OP <= 9.999 ~ 10,
      org_hrs_after_OP >= 10 & org_hrs_after_OP <= 14.999 ~ 15, # short
      org_hrs_after_OP >= 15 & org_hrs_after_OP <= 19.999 ~ 20, # middle
      org_hrs_after_OP >= 20 & org_hrs_after_OP <= 24.999 ~ 25,
      org_hrs_after_OP >= 25 & org_hrs_after_OP <= 35 ~ 30)) # long
  
  df$hrs_group <- as.numeric(df$hrs_group)
  return(df)
}

df_hrs_groups_all_slice <- function(df){
  # have treatment as a word
  df$treatment_word <- ifelse(df$treatment_r == 0, 'CTR', 'HiK')
  
  df <- df %>%
    mutate(hrs_group = case_when(
      org_hrs_after_OP >= 5 & org_hrs_after_OP <= 14.999 ~ 10,
      org_hrs_after_OP >= 15 & org_hrs_after_OP <= 24.999 ~ 20, # short
      org_hrs_after_OP >= 25 & org_hrs_after_OP <= 34.999 ~ 30, # middle
      org_hrs_after_OP >= 35 & org_hrs_after_OP <= 44.999 ~ 40,
      org_hrs_after_OP >= 45 & org_hrs_after_OP <= 51 ~ 47)) # long
  
  df$hrs_group <- as.numeric(df$hrs_group)
  return(df)
}

get_results_df_glm <- function(result_df, sum_model, data_type, var, distr){
    
  rep_num = length(rownames(sum_model$coefficients))
  if (distr == 'gaussian'){
    result_df <- rbind(result_df, data.frame(
      data_type = rep(data_type, rep_num),
      DV = rep(var, rep_num),
      param = rownames(sum_model$coefficients),
      estimate = sum_model$coefficients[,1],
      SE = sum_model$coefficients[,2],
      dfs = sum_model$coefficients[,3],
      t_val = sum_model$coefficients[,4],
      p_val = sum_model$coefficients[,5]))
    
  } else if (distr == 'gamma'){
    
    result_df <- rbind(result_df, data.frame(
      data_type = rep(data_type, rep_num),
      DV = rep(var, rep_num),
      param = rownames(sum_model$coefficients),
      estimate = sum_model$coefficients[,1],
      SE = sum_model$coefficients[,2],
      dfs = rep('no',rep_num),
      t_val = sum_model$coefficients[,3],
      p_val = sum_model$coefficients[,4]))
  }
  return(result_df)
}

get_results_anova_model_comparison <- function(results, anova_df, data_type, var, id2){
  # get results table from model comparison
  num_rep <- length(rownames(results[[names(results)[1]]]))
  
  for (name in names(results)) {
    anova_df <-  rbind(anova_df, data.frame(
      data_type = rep(data_type, num_rep),
      DV = rep(var, num_rep),
      ID2 = rep(id2, num_rep),
      model_comparison = rownames(results[[name]]),
      logLik = results[[name]][['logLik']],
      deviance = results[[name]][['deviance']],
      ch_sq = results[[name]][['Chisq']],
      p_vals = results[[name]][['Pr(>Chisq)']]))
  }
  return(anova_df)
}

get_results_df_r2 <- function(df, r2_df, data_type, var, param){

  if (param == 'intr short') {
    formula_full <- as.formula(paste(var, "~ treatment_r + day + (1|OP) + (1|OP:", data_type, ")"))
  } else if (param == 'intr_inc_only'){
    formula_full <- as.formula(paste(var, "~ treatment_r + (1|OP)"))
  } else if (param == 'intr ext') {
    formula_full <- as.formula(paste(var, "~ treatment_r + hrs_after_OP + patient_age + (1|OP) + (1|OP:", data_type, ")"))
  } else if (param == 'firing short') {
    formula_full <- as.formula(paste("VAL ~ inj_current * treatment_r + day  + (1|OP) + (1|OP:", data_type, ")"))
  } else if (param == 'firing ext'){
    formula_full <- as.formula(paste("VAL ~ inj_current + treatment_r * patient_age + hrs_after_OP + (1|OP) + (1|OP:", data_type, ")"))
  } else if (param == 'firing_short_inc_only') {
    formula_full <- as.formula("VAL ~ inj_current + treatment_r + (1|OP)")
  } else if (param == 'intr_ext_inc_only'){
    formula_full <- as.formula(paste(var, "~ treatment_r + hrs_after_OP + (1|OP)"))
  } else if (param == 'firing_ext_inc_only'){
    formula_full <- as.formula(paste("VAL ~ inj_current + treatment_r + hrs_after_OP + (1|OP)"))
  } else if (param == 'slice_all_CTR'){
    formula_full <- as.formula(paste(var, "~ hrs_after_OP + (1|OP) + (1|OP:", data_type, ")"))
  }
  
  
  FullModel <- lmer(formula_full, data = df, control =  lmerControl(optimizer = "Nelder_Mead"))
  r2 <- r2beta(FullModel)
  
  num_reps <- length(r2$Effect)
  r2_df <- rbind(r2_df, data.frame(
    data_type = rep(data_type, num_reps),
    DV = rep(var, num_reps),
    param = r2$Effect,
    R_squared = r2$Rsq,
    upper_CL = r2$upper.CL,
    lower_CL = r2$lower.CL))
  
  return(r2_df)
}

marg_effects_intr_short <- function(df, data_type, var, var_org, attrs,emm_df, emm_CI_df){
  
  formula_full <- as.formula(paste(var, "~ treatment_r + day + (1|OP) + (1|OP:", data_type, ")"))

  FullModel <- lmer(formula_full, data = df, control =  lmerControl(optimizer = "Nelder_Mead"))
  
  # Post-hoc analysis with Estimated Marginal Means ####
  MarginalEffects <- emmeans(object = FullModel, 
                             specs =  ~ treatment_r + day,
                             type="response", tran = "log") #log transforms
  
  #plot(MarginalEffects)
  ContrastMatrix <- diag(4)
  gr_day <- MarginalEffects@grid$day
  gr_tr <- MarginalEffects@grid$treatment_r
  
  ### Contrast and adjustment for multiple testing (Benjamini-Hochberg) ####
  MarginalEffect_tr_test <- contrast(object = MarginalEffects, type = "response",
                                 method = list("(Ctrl D1) / (Ctrl D2)" = ContrastMatrix[gr_tr==0&gr_day==0,] - ContrastMatrix[gr_tr==0&gr_day==1,],
                                               "(high K D1) / (high K D2)" = ContrastMatrix[gr_tr==0&gr_day==0,] - ContrastMatrix[gr_tr==1&gr_day==1,],
                                               "(Ctrl D2) / (high K D2)" = ContrastMatrix[gr_tr==0&gr_day==1,] - ContrastMatrix[gr_tr==1&gr_day==1,]),
                                 adjust="BH")

  # print(MarginalEffectTest)
  
  ### Plot marginal effects ####
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
  

  return(list(
    df_emmeans = emm_df, 
    df_CIs = emm_CI_df))
}

marg_effects_firing_short <- function(df, data_type, var, var_org, attrs, emm_df, emm_CI_df){

  formula_full <- as.formula(paste("VAL ~ inj_current + treatment_r + day  + (1|OP) + + (1|OP:", data_type, ")"))
  FullModel <- lmer(formula_full, data = df)
  
  # Post-hoc analysis with Estimated Marginal Means ####
  MarginalEffects <- emmeans(object = FullModel, 
                             specs =  ~ inj_current + treatment_r + day,
                             at = list(inj_current = c(200, 400, 600)),
                             type="response", tran = "log") #log transforms
  #plot(MarginalEffects)
  ContrastMatrix <- diag(12)
  gr_day <- MarginalEffects@grid$day
  gr_tr <- MarginalEffects@grid$treatment_r
  gr_inj <- MarginalEffects@grid$inj_current
  
  ### Contrast and adjustment for multiple testing (Benjamini-Hochberg) ####
  MarginalEffectTest <- contrast(object = MarginalEffects, type = "response",
                                 method = list("200 pA (Ctrl D1) / (Ctrl D2)" = ContrastMatrix[gr_tr==0&gr_day==0&gr_inj==200,] - ContrastMatrix[gr_tr==0&gr_day==1&gr_inj==200,],
                                               "400 pA (Ctrl D1) / (Ctrl D2)" = ContrastMatrix[gr_tr==0&gr_day==0&gr_inj==400,] - ContrastMatrix[gr_tr==0&gr_day==1&gr_inj==400,],
                                               "600 pA (Ctrl D1) / (Ctrl D2)" = ContrastMatrix[gr_tr==0&gr_day==0&gr_inj==600,] - ContrastMatrix[gr_tr==0&gr_day==1&gr_inj==600,],
                                               "200 pA (high K D1) / (high K D2)" = ContrastMatrix[gr_tr==0&gr_day==0&gr_inj==200,] - ContrastMatrix[gr_tr==1&gr_day==1&gr_inj==200,],
                                               "400 pA (high K D1) / (high K D2)" = ContrastMatrix[gr_tr==0&gr_day==0&gr_inj==400,] - ContrastMatrix[gr_tr==1&gr_day==1&gr_inj==400,],
                                               "600 pA (high K D1) / (high K D2)" = ContrastMatrix[gr_tr==0&gr_day==0&gr_inj==600,] - ContrastMatrix[gr_tr==1&gr_day==1&gr_inj==600,]),
                                 # "600 pA (Ctrl D2) / (high K D2)" = ContrastMatrix[gr_tr=="Ctrl"&gr_day=='D2'&gr_inj=='600',] - ContrastMatrix[gr_tr=="high K"&gr_day=='D2'&gr_inj=='600',]),
                                 adjust="BH")
  
  m_contrasts <- summary(MarginalEffectTest)
  emm_CI_df <- rbind(emm_CI_df, data.frame(
    data_type = rep(data_type, length(m_contrasts$contrast)),
    IV = rep('day', length(m_contrasts$contrast)), #repeat as long as the contrasts are
    DV = rep(var, length(m_contrasts$contrast)),
    ratio = m_contrasts$ratio,
    contrast = m_contrasts$contrast,
    se = m_contrasts$SE * (attrs$'scaled:scale') + attrs$'scaled:center',
    p_vales = m_contrasts$p.value))
  
  # introducing this so that the emmean estimated and CIs are correct
  # some mistake due to the 3 predicotrs occurs when log scaled
  MarginalEffects_report <- emmeans(object = FullModel, 
                             specs =  ~ inj_current + treatment_r + day,
                             at = list(inj_current = c(200, 400, 600)),
                             tran = "link")
  
  m_effects <- summary(MarginalEffects_report)
  emm_df <- rbind(emm_df , data.frame(
    data_type = rep(data_type, length(m_effects$day)),
    IV =  m_effects$day,
    DV = rep(var, length(m_effects$day)),
    treatment =  m_effects$treatment_r,
    inj = m_effects$inj_current,
    response = m_effects$emmean,
    SE = m_effects$SE,
    lower_CI =  m_effects$lower.CL,
    upper_CI =  m_effects$upper.CL))
  
  # m_effects <- summary(MarginalEffects_report)
  # emm_df <- rbind(emm_df , data.frame(
  #   data_type = rep(data_type, length(m_effects$day)),
  #   IV =  m_effects$day,
  #   DV = rep(var, length(m_effects$day)),
  #   treatment =  m_effects$treatment_r,
  #   response = m_effects$response * (attrs$'scaled:scale') + attrs$'scaled:center',
  #   SE = m_effects$SE * (attrs$'scaled:scale') + attrs$'scaled:center',
  #   lower_CI =  m_effects$lower.CL * (attrs$'scaled:scale') + attrs$'scaled:center',
  #   upper_CI =  m_effects$upper.CL * (attrs$'scaled:scale') + attrs$'scaled:center'))
  
  return(list(
    df_emmeans = emm_df, 
    df_CIs = emm_CI_df))
}

marg_effects_firing_ext <- function(melt_df, var_firing, data_type, emm_df_ext, emm_CI_df_ext){
  
  fire_formula_full <- as.formula(paste("VAL ~ inj_current + treatment * hrs_group + (1|OP) + (1|OP:", data_type, ")"))
  fire_FullModel <- lmer(fire_formula_full, data = melt_df, control =  lmerControl(optimizer = "Nelder_Mead"))
  
  # Post-hoc analysis with Estimated Marginal Means ####
  MarginalEffects <- emmeans(object = fire_FullModel, 
                             specs =  ~ inj_current +  treatment * hrs_group,
                             at = list(inj_current = c(400, 600), hrs_group = c(29,45)))
  
  ContrastMatrix <- diag(8)
  gr_inj <- MarginalEffects@grid$inj_current
  gr_hrs <- MarginalEffects@grid$hrs_group
  gr_tr <- MarginalEffects@grid$treatment
  
  MarginalEffectTest <- contrast(object = MarginalEffects, type = "response",
                                 method = list("400 pA after 29 CTR vs HiK" = ContrastMatrix[gr_inj==400&gr_hrs==29&gr_tr=='Ctrl',] - ContrastMatrix[gr_inj==400&gr_hrs==29&gr_tr=='high K',],
                                               "600 pA after 29 h CTR vs HiK" = ContrastMatrix[gr_inj==600&gr_hrs==29&gr_tr=='Ctrl',] - ContrastMatrix[gr_inj==600&gr_hrs==29&gr_tr=='high K',],
                                               "400 pA after 29 to 45 CTR vs HiK" = (ContrastMatrix[gr_inj==400&gr_hrs==29&gr_tr=='Ctrl',] - ContrastMatrix[gr_inj==400&gr_hrs==45&gr_tr=='Ctrl',]) -
                                                 (ContrastMatrix[gr_inj==400&gr_hrs==29&gr_tr=='high K',] - ContrastMatrix[gr_inj==400&gr_hrs==45&gr_tr=='high K',]),
                                               "600 pA after 29 to 45 CTR vs HiK" = (ContrastMatrix[gr_inj==600&gr_hrs==29&gr_tr=='Ctrl',] - ContrastMatrix[gr_inj==600&gr_hrs==45&gr_tr=='Ctrl',]) -
                                                 (ContrastMatrix[gr_inj==600&gr_hrs==29&gr_tr=='high K',] - ContrastMatrix[gr_inj==600&gr_hrs==45&gr_tr=='high K',])),
                                 adjust="BH")
  # append all to a table
  m_effects <- summary(MarginalEffects)
  emm_df_ext <- rbind(emm_df_ext, data.frame(
    data_type = rep(data_type, length(m_effects$treatment)),
    DV = rep(var_firing, length(m_effects$treatment)),
    inj =  m_effects$inj_current, 
    hrs_group = (m_effects$hrs_group),
    treatment =  m_effects$treatment,
    response = m_effects$emmean,
    SE = m_effects$SE,
    lower_CI =  m_effects$lower.CL,
    upper_CI =  m_effects$upper.CL))
  
  # append all toa table
  m_contrasts <- summary(MarginalEffectTest)
  emm_CI_df_ext <- rbind(emm_CI_df_ext, data.frame(
    data_type = rep(data_type, length(m_contrasts$contrast)),
    IV = rep('tr*hrs_after_OP_group', length(m_contrasts$contrast)), #repeat as long as the contrasts are
    DV = rep(var_firing, length(m_contrasts$contrast)),
    ratio = m_contrasts$t.ratio,
    contrast = m_contrasts$contrast,
    estimate = m_contrasts$estimate,
    se = m_contrasts$SE,
    p_vales = m_contrasts$p.value))
  
  return(list(
    df_emmeans = emm_df_ext, 
    df_CIs = emm_CI_df_ext))
}

marg_effects_intr_ext_no_pat_age <- function(df, data_type, var, var_org, attrs, emm_df, emm_CI_df){
  
  # create groups of hours of incubation
  df <- df_ext_comparison_groups_inc_only(df)
  
  if (data_type ==''){
    formula_full <- as.formula(paste(var, "~ treatment_r + hrs_group + (1|OP)"))
  }
  else {
    formula_full <- as.formula(paste(var, "~ treatment_r + hrs_group + (1|OP) + (1|OP:", data_type, ")"))
  }
  
  FullModel <- lmer(formula_full, data = df, control =  lmerControl(optimizer = "Nelder_Mead"))
  
  # Post-hoc analysis with Estimated Marginal Means ####
  MarginalEffects <- emmeans(object = FullModel, 
                             specs =  ~ treatment_r + hrs_group,
                             type="response", tran = "log",
                             at = list(hrs_group = c(10, 25))) #log transforms
  
  #plot(MarginalEffects)
  ContrastMatrix <- diag(4)
  gr_hrs <- MarginalEffects@grid$hrs_group
  gr_tr <- MarginalEffects@grid$treatment_r
  
  ### Contrast and adjustment for multiple testing (Benjamini-Hochberg) ####
  MarginalEffect_tr_test <- contrast(object = MarginalEffects, type = "response",
                                     method = list("(Ctrl 10 hrs) / (Ctrl 25hrs)" = ContrastMatrix[gr_tr==0&gr_hrs==10,] - ContrastMatrix[gr_tr==0&gr_hrs==25,],
                                                   "(high K 10hrs) / (high K 25 hrs)" = ContrastMatrix[gr_tr==0&gr_hrs==10,] - ContrastMatrix[gr_tr==1&gr_hrs==25,],
                                                   "(Ctrl 25 hrs) / (high K 25 hrs)" = ContrastMatrix[gr_tr==0&gr_hrs==25,] - ContrastMatrix[gr_tr==1&gr_hrs==25,]),
                                     adjust="BH")
  
  # print(MarginalEffectTest)
  
  ### Plot marginal effects ####
  MarginalPlot <- plot(MarginalEffect_tr_test, colors = "grey10") +
    labs(title = paste(var, "Effects of Treatment and Day"))+
    ylab(label = "")+
    geom_vline(xintercept = 1, linetype="dashed")+
    xlab(label = "Ratio")+
    theme_classic()
  MarginalPlot
  
  m_effects <- summary(MarginalEffects)
  emm_df <- rbind(emm_df , data.frame(
    data_type = rep(data_type, length(m_effects$hrs_group)),
    IV =  m_effects$hrs_group, 
    DV = rep(var, length(m_effects$hrs_group)),
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
  
  
  return(list(
    df_emmeans = emm_df, 
    df_CIs = emm_CI_df))
}

marg_effects_intr_hrs_only <- function(df, var, var_org, attrs, emm_df, emm_CI_df){
  
  # create groups of hours of incubation
  df <- df_hrs_groups_all_slice(df)

  formula_full <- as.formula(paste(var, "~ hrs_group + (1|OP) + (1|OP:slice)"))
  FullModel <- lmer(formula_full, data = df, control =  lmerControl(optimizer = "Nelder_Mead"))
  
  # Post-hoc analysis with Estimated Marginal Means ####
  MarginalEffects <- emmeans(object = FullModel, 
                             specs =  ~ hrs_group,
                             type="response", tran = "log",
                             at = list(hrs_group = c(10, 30, 47))) #log transforms
  
  #plot(MarginalEffects)
  ContrastMatrix <- diag(3)
  gr_hrs <- MarginalEffects@grid$hrs_group
  
  ### Contrast and adjustment for multiple testing (Benjamini-Hochberg) ####
  MarginalEffect_tr_test <- contrast(object = MarginalEffects, type = "response",
                                     method = list("(Ctrl 10 hrs) / (Ctrl 30hrs)" = ContrastMatrix[gr_hrs==10,] - ContrastMatrix[gr_hrs==30,],
                                                   "(Ctr 10 hrs) / (Ctrl 47 hrs)" = ContrastMatrix[gr_hrs==10,] - ContrastMatrix[gr_hrs==47,]),
                                     adjust="BH")
  
  m_effects <- summary(MarginalEffects)
  emm_df <- rbind(emm_df , data.frame(
    data_type = rep("slice all CTR", length(m_effects$hrs_group)),
    IV =  m_effects$hrs_group, 
    DV = rep(var, length(m_effects$hrs_group)),
    response = m_effects$response * (attrs$'scaled:scale') + attrs$'scaled:center',
    SE = m_effects$SE * (attrs$'scaled:scale') + attrs$'scaled:center',
    lower_CI =  m_effects$lower.CL * (attrs$'scaled:scale') + attrs$'scaled:center',
    upper_CI =  m_effects$upper.CL * (attrs$'scaled:scale') + attrs$'scaled:center'))
  
  m_contrasts <- summary(MarginalEffect_tr_test)
  emm_CI_df <- rbind(emm_CI_df, data.frame(
    data_type = rep("slice all CTR", length(m_contrasts$contrast)),
    IV = rep('day', length(m_contrasts$contrast)), #repeat as long as the contrasts are
    DV = rep(var, length(m_contrasts$contrast)),
    ratio = m_contrasts$ratio,
    contrast = m_contrasts$contrast,
    se = m_contrasts$SE * (attrs$'scaled:scale') + attrs$'scaled:center',
    p_vales = m_contrasts$p.value))
  
  return(list(
    df_emmeans = emm_df, 
    df_CIs = emm_CI_df))
}


### FUNCS FOR INC ONLY ####

lmers_slice_inc_only <- function(df, var) {
  # make positive
  # df$treatment_r <- ifelse(df$treatment == "Ctrl", 0, 1)
  if (mean(df[[var]], na.rm = TRUE) < 0) {
    df[[var]] = -df[[var]] 
  }
  #remove rows with nan values
  df <- df[!is.na(df[[var]]), ]
  
  formula_full <- as.formula(paste(var, "~ treatment_r + (1|OP)"))
  formula_tr_0 <- as.formula(paste(var, "~ (1|OP)"))
  
  FullModel <- lmer(formula_full,data = df, REML = FALSE)
  TreatmentNullModel <-lmer(formula_tr_0, data = df, REML = FALSE)
  
  anova_full_tr_0 <- anova(FullModel, TreatmentNullModel)
  
  return(list(full_vs_treat_null = anova_full_tr_0))
}

lmers_ext_inc_only <- function(df, var, id2) {
  # make positive
  if (mean(df[[var]], na.rm = TRUE) < 0) {
    df[[var]] = -df[[var]] 
  }
  #remove rows with nan values
  df <- df[!is.na(df[[var]]), ]

  formula_full <- as.formula(paste(var, "~ treatment_r +", id2, " + (1|OP)"))
  formula_tr_0 <- as.formula(paste(var, "~  ", id2 ," + (1|OP)"))
  formula_hrs_0 <- as.formula(paste(var, "~  + treatment_r  + (1|OP)"))
  
  FullModel <-    lmer(formula_full, data = df, REML = FALSE)
  TrtNullModel <- lmer(formula_tr_0, data = df, REML = FALSE)
  IncNullModel <- lmer(formula_hrs_0, data = df, REML = FALSE)
  #InteractNullModel <- lmer(formula_interaction_0, data = df, REML = FALSE)
  
  anova_full_tr_0 <- anova(FullModel, TrtNullModel)
  anova_full_inc_0 <- anova(FullModel, IncNullModel)
  #anova_full_interact_0 <- anova(FullModel, InteractNullModel)
  
  return(list(
    full_vs_treat_null = anova_full_tr_0, 
    full_vs_inc_null = anova_full_inc_0))
  # full_vs_interact_null = anova_full_interact_0))
}

marg_effects_intr_short_inc_only <- function(df, data_type, var, var_org, attrs, emm_df, emm_CI_df){

  formula_full <- as.formula(paste(var, "~ treatment_r + (1|OP)"))
  FullModel <- lmer(formula_full, data = df)
  
  # Post-hoc analysis with Estimated Marginal Means ####
  MarginalEffects <- emmeans(object = FullModel, 
                             specs =  ~ treatment_r,
                             type="response", tran = "log") #log transforms
  
  #plot(MarginalEffects)
  ContrastMatrix <- diag(2)
  gr_tr <- MarginalEffects@grid$treatment_r
  
  ### Contrast and adjustment for multiple testing (Benjamini-Hochberg) ####
  MarginalEffect_tr_test <- contrast(object = MarginalEffects, type = "response",
                                     method = list("(Ctrl) / (high K)" = ContrastMatrix[gr_tr==0,] - ContrastMatrix[gr_tr==1,]),
                                     adjust="BH")
  
  # print(MarginalEffectTest)
  
  ### Plot marginal effects ####
  MarginalPlot <- plot(MarginalEffect_tr_test, colors = "grey10") +
    labs(title = paste(var, "Effects of Treatment and Day"))+
    ylab(label = "")+
    geom_vline(xintercept = 1, linetype="dashed")+
    xlab(label = "Ratio")+
    theme_classic()
  MarginalPlot
  
  m_effects <- summary(MarginalEffects)
  emm_df <- rbind(emm_df , data.frame(
    data_type = rep(data_type, length(m_effects$treatment_r)),
    DV = rep(var, length(m_effects$treatment_r)),
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
  
  
  return(list(
    df_emmeans = emm_df, 
    df_CIs = emm_CI_df))
}

marg_effects_mea <- function(df, data_type, attrs, emm_df, emm_CI_df){
  
  formula_full <- as.formula(paste("value ~ treatment_r + (1|OP)"))
  FullModel <- lmer(formula_full, data = df)
  
  # Post-hoc analysis with Estimated Marginal Means ####
  MarginalEffects <- emmeans(object = FullModel, 
                             specs =  ~ treatment_r,
                             type="response", tran = "log") #log transforms
  
  #plot(MarginalEffects)
  ContrastMatrix <- diag(2)
  gr_tr <- MarginalEffects@grid$treatment_r
  
  ### Contrast and adjustment for multiple testing (Benjamini-Hochberg) ####
  MarginalEffect_tr_test <- contrast(object = MarginalEffects, type = "response",
                                     method = list("(Ctrl) / (high K)" = ContrastMatrix[gr_tr==0,] - ContrastMatrix[gr_tr==1,]),
                                     adjust="BH")
  
  m_effects <- summary(MarginalEffects)
  emm_df <- rbind(emm_df , data.frame(
    data_type = rep(data_type, length(m_effects$treatment_r)),
    treatment =  m_effects$treatment_r,
    response = m_effects$response * (attrs$'scaled:scale') + attrs$'scaled:center',
    SE = m_effects$SE * (attrs$'scaled:scale') + attrs$'scaled:center',
    lower_CI =  m_effects$lower.CL * (attrs$'scaled:scale') + attrs$'scaled:center',
    upper_CI =  m_effects$upper.CL * (attrs$'scaled:scale') + attrs$'scaled:center'))
  
  m_contrasts <- summary(MarginalEffect_tr_test)
  emm_CI_df <- rbind(emm_CI_df, data.frame(
    data_type = rep(data_type, length(m_contrasts$contrast)),
    ratio = m_contrasts$ratio,
    contrast = m_contrasts$contrast,
    se = m_contrasts$SE * (attrs$'scaled:scale') + attrs$'scaled:center',
    p_vales = m_contrasts$p.value))
  
  
  return(list(
    df_emmeans = emm_df, 
    df_CIs = emm_CI_df))
}


marg_effects_intr_ext_inc_only <- function(df, data_type, var, var_org, attrs, emm_df, emm_CI_df){
  df$treatment_r <- ifelse(df$treatment == "Ctrl", 0, 1)
  
  # create groups of hours of incubation
  df <- df_ext_comparison_groups_inc_only(df)
  
  formula_full <- as.formula(paste(var, "~ treatment_r + hrs_group + (1|OP)"))
  FullModel <- lmer(formula_full, data = df, control =  lmerControl(optimizer = "Nelder_Mead"))
  
  # Post-hoc analysis with Estimated Marginal Means ####
  MarginalEffects <- emmeans(object = FullModel, 
                             specs =  ~ treatment_r + hrs_group,
                             type="response", tran = "log",
                             at = list(hrs_group = c(10, 25))) #log transforms
  
  #plot(MarginalEffects)
  ContrastMatrix <- diag(4)
  gr_hrs <- MarginalEffects@grid$hrs_group
  gr_tr <- MarginalEffects@grid$treatment_r
  
  ### Contrast and adjustment for multiple testing (Benjamini-Hochberg) ####
  MarginalEffect_tr_test <- contrast(object = MarginalEffects, type = "response",
                                     method = list("(Ctrl 10 hrs) / (Ctrl 25hrs)" = ContrastMatrix[gr_tr==0&gr_hrs==10,] - ContrastMatrix[gr_tr==0&gr_hrs==25,],
                                                   "(high K 10hrs) / (high K 25 hrs)" = ContrastMatrix[gr_tr==1&gr_hrs==10,] - ContrastMatrix[gr_tr==1&gr_hrs==25,],
                                                   "(Ctrl 25 hrs) / (high K 25 hrs)" = ContrastMatrix[gr_tr==0&gr_hrs==25,] - ContrastMatrix[gr_tr==1&gr_hrs==25,]),
                                     adjust="BH")
  
  # print(MarginalEffectTest)
  
  ### Plot marginal effects ####
  MarginalPlot <- plot(MarginalEffect_tr_test, colors = "grey10") +
    labs(title = paste(var, "Effects of Treatment and Day"))+
    ylab(label = "")+
    geom_vline(xintercept = 1, linetype="dashed")+
    xlab(label = "Ratio")+
    theme_classic()
  MarginalPlot
  
  m_effects <- summary(MarginalEffects)
  emm_df <- rbind(emm_df, data.frame(
    data_type = rep(data_type, length(m_effects$hrs_group)),
    IV =  m_effects$hrs_group, 
    DV = rep(var, length(m_effects$hrs_group)),
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
  
  
  return(list(
    df_emmeans = emm_df, 
    df_CIs = emm_CI_df))
}

lmers_firing_inc_only <- function(df, var, id2, data_type) {
  # df$treatment_r <- ifelse(df$treatment == "Ctrl", 0, 1)
  
  formula_full <- as.formula(paste("VAL~ treatment_r + inj_current + (1|OP)"))
  formula_tr_0 <- as.formula(paste("VAL ~ inj_current + (1|OP)"))

  FullModel <- lmer(formula_full, data = df, REML = FALSE)
  TreatmentNullModel <- lmer(formula_tr_0, data = df, REML = FALSE)

  anova_full_tr_0 <- anova(FullModel, TreatmentNullModel)

  return(list(full_vs_treat_null = anova_full_tr_0))
}

lmers_firing_extended_inc_only <- function(df, var, data_type) {
  # df$treatment_r <- ifelse(df$treatment == "Ctrl", 0, 1)
  
  formula_full <- as.formula(paste("VAL ~ inj_current + treatment_r + hrs_after_OP + (1|OP)"))
  formula_tr_0 <- as.formula(paste("VAL ~ inj_current + hrs_after_OP  + (1|OP)"))
  formula_hrs_0 <- as.formula(paste("VAL ~ inj_current + treatment_r  + (1|OP)"))
  
  FullModel <- lmer(formula_full, data = df, REML = FALSE)
  TreatmentNullModel <- lmer(formula_tr_0, data = df, REML = FALSE)
  HrsNullModel <- lmer(formula_hrs_0, data = df, REML = FALSE)
  
  anova_full_tr_0 <- anova(FullModel, TreatmentNullModel)
  anova_full_hrs_0 <- anova(FullModel, HrsNullModel)
  
  return(list(
    full_vs_treat_null = anova_full_tr_0,
    full_vs_day_null = anova_full_hrs_0))
}

marg_effects_firing_short_inc_only <- function(df, data_type, var, var_org, attrs, emm_df, emm_CI_df){
  # df$treatment_r <- ifelse(df$treatment == "Ctrl", 0, 1)
  
  formula_full <- as.formula(paste("VAL ~ inj_current + treatment_r + (1|OP)"))
  FullModel <- lmer(formula_full, data = df)
  
  # Post-hoc analysis with Estimated Marginal Means ####
  MarginalEffects <- emmeans(object = FullModel, 
                             specs =  ~ inj_current + treatment_r,
                             at = list(inj_current = c(200, 400, 600)),
                             type="response", tran = "log") #log transforms
  #plot(MarginalEffects)
  ContrastMatrix <- diag(6)
  gr_tr <- MarginalEffects@grid$treatment_r
  gr_inj <- MarginalEffects@grid$inj_current
  
  ### Contrast and adjustment for multiple testing (Benjamini-Hochberg) ####
  MarginalEffectTest <- contrast(object = MarginalEffects, type = "response",
                                 method = list("200 pA (Ctrl) / (high K)" = ContrastMatrix[gr_tr==0&gr_inj==200,] - ContrastMatrix[gr_tr==1&gr_inj==200,],
                                               "400 pA (Ctrl) / (high K" = ContrastMatrix[gr_tr==0&gr_inj==400,] - ContrastMatrix[gr_tr==1&gr_inj==400,],
                                               "600 pA (Ctrl) / (high K)" = ContrastMatrix[gr_tr==0&gr_inj==600,] - ContrastMatrix[gr_tr==1&gr_inj==600,]),
                                 adjust="BH")
  
  m_contrasts <- summary(MarginalEffectTest)
  emm_CI_df <- rbind(emm_CI_df, data.frame(
    DV = rep(var, length(m_contrasts$contrast)),
    ratio = m_contrasts$ratio,
    contrast = m_contrasts$contrast,
    se = m_contrasts$SE * (attrs$'scaled:scale') + attrs$'scaled:center',
    p_vales = m_contrasts$p.value))
  
  # introducing this so that the emmean estimated and CIs are correct
  # some mistake due to the 3 predicotrs occurs when log scaled
  MarginalEffects_report <- emmeans(object = FullModel, 
                                    specs =  ~ inj_current + treatment_r,
                                    at = list(inj_current = c(200, 400, 600)),
                                    tran = "link")
  
  m_effects <- summary(MarginalEffects_report)
  emm_df <- rbind(emm_df , data.frame(
    DV = rep(var, length(m_contrasts$contrast)),
    treatment =  m_effects$treatment_r,
    inj = m_effects$inj_current,
    response = m_effects$emmean,
    SE = m_effects$SE,
    lower_CI =  m_effects$lower.CL,
    upper_CI =  m_effects$upper.CL))
  
  return(list(
    df_emmeans = emm_df, 
    df_CIs = emm_CI_df))
}


#### CALCULATING SUMMARY DF FOR PLOTTING FIRING ####
# library(plotrix )
# 
# data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/'
# # area <- 'temporal' # or 'complete'
# # df_repatch <- fread(paste(data_dir, 'repatch_data_', area,'.csv', sep = ''), sep = ",")
# # df_slice <- fread(paste(data_dir, 'slice_data_', area,'.csv', sep = ''), sep = ",")
# 
# df_slice_inc_only <- fread(paste(data_dir, 'slice_all.csv', sep = ''), sep = ",")
# filter_hrs = 31
# data_type <- paste('slice under ', filter_hrs, 'hrs')
# df_slice_inc_only_f <- df_slice_inc_only[df_slice_inc_only$hrs_after_OP < filter_hrs, ]
# 
# var_of_interest <- c('num_aps', 'IFF')
# # define ID variables
# ids <- c('OP', 'day', 'treatment', 'hrs_after_OP', 'patient_age', 'slice')
# 
# df_slice_inc_only_f[['slice']] <- factor(df_slice_inc_only_f[['slice']])
# df <- fix_df(df_slice_inc_only_f)
# CIs_complete <- data.frame()
# 
# for (var_firing in var_of_interest){
#   melt_df <- get_melt_df_firing(var_firing, df, ids)
#   
#   CIs_df <- as.data.frame(melt_df %>%
#                             group_by(treatment, day, inj_current) %>%
#                             summarize(means = mean(VAL),
#                                       median = median(VAL),
#                                       SE = std.error(VAL),
#                                       CI_L = confint(lm(VAL ~ 1), level=0.95)[1],
#                                       CI_U = confint(lm(VAL ~ 1), level=0.95)[2]))
#   
#   CIs_df$data_type <- rep(data_type, nrow(CIs_df))
#   CIs_df$var_firing <- rep(var_firing, nrow(CIs_df))
#   CIs_complete <- rbind(CIs_complete, CIs_df)
# }
# 
# # save_dir <- '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/stats/sum_data_CIs_SEs_for_plotting/'
# # write.xlsx(CIs_complete, paste(save_dir, 'firing_plot_CIs_', area,'.xlsx',sep = ''))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # # maybe helpful to get quick overview
# # # t %>% takes the output of the expression on its left and passes it as the first argument to the function on its right
# # df_repatch %>% 
# #   group_by(treatment, day) %>%
# #   summarise(
# #     count = n(),
# #     mean_TH = mean(TH),
# #     mean_RMP = mean(resting_potential)
# #   )
# # 
# # # equivalent to aggregate(df_repatch$TH, by=list(df_repatch$day, df_repatch$treatment), FUN=mean)
# 
# 
# 