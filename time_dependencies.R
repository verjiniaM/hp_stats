setwd("/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/figures/analysis_try_01Feb2022/")

rm(list = ls())

scaleFUN <- function(x) sprintf("%.2f", x) # adds the desired number of elements after the decimal point
#function that allows to add axis text on 2 lines
addline_format <- function(x,...){
  gsub('\\s\\s','\n',x)
}

theme_set(theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(color = "black"),
                axis.text=element_text(size=14), axis.title=element_text(size=16),
                legend.title = element_text(size = 14), legend.text = element_text(size = 12),
                strip.text.x = element_text(size = 16),
                plot.title = element_text(size = rel(1.8), hjust = 0.5, margin=margin(0,0,15,0))))

#desce_stat_tab <- "/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/summary_data_tables/2022-03-12_complete+times.csv"
desce_stat_tab <- "/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/OP220228/data_tables/OP220228_Intrinsic_and_synaptic_properties.csv"
data_tab <- read.csv(desce_stat_tab, header = T, na.strings=c("",".","NA")) #use ";" as separator when saved from excel


#factorize string data
data_tab$tissue_source <- as.factor(data_tab$tissue_source)
data_tab$OP <- as.factor(data_tab$OP)
data_tab$patcher <- as.factor(data_tab$patcher)
data_tab$filename <- as.factor(data_tab$filename)
data_tab$slice <- as.factor(data_tab$slice)
data_tab$day <- as.factor(data_tab$day)
data_tab$treatment <- as.factor(data_tab$treatment)
data_tab$cell_ID<- as.factor(data_tab$cell_ID)
data_tab$repatch <- as.factor(data_tab$repatch)
data_tab$age_cat <- with(data_tab, ifelse(patient_age<15,yes = "J",no = "A"))

#selecting adult/J data
data_adult_all <- data_tab %>% filter(age_cat == "A")

#removing missing columns
data_adult <- data_adult_all %>% drop_na(treatment)
data_adult <- data_adult %>% drop_na(day)

#Repatch data
data_adult_BA <- data_adult %>% filter(repatch == "n",
                                            max_spikes < 100,
                                            day != "",
                                            Rs < 30,
                                            max_repol > -120,
                                            capacitance < 20000)
data_adult_BA <- data_adult %>% filter(
                                       max_spikes < 100,
                                       day != "",
                                       Rs < 30,
                                       max_repol > -120,
                                       capacitance < 20000)

y_names <- list("mV", "pA", "mV", "mV", "F",
                "mV/ms", "mV/ms", " ","tau", "mV","Ω","Ω")

titles <- list("AP amplitude", "Rheobase", "AP threshold", "Vm", "Capacitance",
               "Upstroke", "Downstroke", " Max number of spikes", "Membrane time constant",
               "Resting potential", "Series resistance", "Input resistance") 

#filtering
data_ctrl <- data_adult_BA %>% filter(treatment == 'Ctrl') #select treatment

for (j in 16:27){
  data_plot <- cbind(data_ctrl[3:15], data_ctrl[j])
  names(data_plot)[14] <- "y_var"
  
  plot1 <- ggplot(data = data_plot, aes(x = hrs_after_op, y = y_var))+
    geom_point(aes(color=OP))+
    labs(x="Times after OP (hrs)", y= y_names[j-15])+
    #geom_smooth(method = "lm", se = TRUE, aes(color = "Linear"), size=.7)+
    #scale_color_manual(name = "Model fit", breaks = "Linear", values = c("Linear" = "blue"))+
    #theme(legend.justification = "top")+ #theme(legend.position = c(0.5, 0.6))
    #theme(legend.position = "none")+
    ggtitle(titles[j-15])
  
  #save_folder <- '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/analysis/plots_lab_seminar_16.March.22/time_dep/'
  save_folder <- '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/Human_tissue_meetings/OP220228_'
  plot_name_png <- paste(save_folder, titles[j-15], "OP.png", sep="")
  ggsave(plot_name_png, plot1 , device = "png", width = 15, height = 10, units = "cm",  dpi = 400) 
}

data_resting_repatch <-  data_ctrl %>% filter(resting_potential < -50, resting_potential > -110)

plot2 <- ggplot(data = data_resting_repatch, aes(x = hrs_after_op, y = resting_potential))+
  geom_point()+
  labs(x="Times after OP (hrs)", y= "mV")+
  geom_smooth(method = "lm", se = TRUE, aes(color = "Linear"), size=.7)+
  #scale_color_manual(name = "Model fit", breaks = "Linear", values = c("Linear" = "blue"))+
  #theme(legend.justification = "top")+ #theme(legend.position = c(0.5, 0.6))
  theme(legend.position = "none")+
  ggtitle("Resting potential")

save_folder <- '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/analysis/plots_lab_seminar_16.March.22/time_dep/'
#save_folder <- '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/Human_tissue_meetings/OP220228_'
plot_name_png <- paste(save_folder, "Resting potential", "OP.png", sep="")
ggsave(plot_name_png, plot2 , device = "png", width = 15, height = 10, units = "cm",  dpi = 400) 

#For performing statistical tests if neccessary

time <- data_ctrl$hrs_after_op
TH <- data_ctrl$TH
summary(lm(TH ~ time))

ampl <- data_ctrl$AP_heigth
summary(lm(ampl ~ time))

Rin <- data_ctrl$Rin
summary(lm(Rin ~ time))

down <- data_ctrl$max_repol
summary(lm(down ~ time))

up <- data_ctrl$max_depol
summary(lm(up ~ time))

spikes <- data_ctrl$max_spikes
summary(lm(spikes ~ time))
  
  
  
  

  
  
  