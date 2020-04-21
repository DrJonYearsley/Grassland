###################################################################################
## This script fits a linear model of the soil mositure resilience measure with  ##
## the fixed effect Day (days since the end of the simulated drought; covariate) ##
## and treatment (Control vs. Drought) and with random effects Day (see above)   ##
## and Farm (1-5) nested in Region (Border, Cork, Dublin, Limerick)              ##
## The response is the ratio of soil moisture in each plot to the average of the ##
## control plots at the same site (farm) and time step.                          ##
##                                                                               ##
## Author of the script:                                                         ##
## Maja Ilic M.Ilic@qub.ac.uk                                                    ##
## first modified: 03 Mar 2020                                                   ##
## last modified: 05 Mar 2020                                                    ##
###################################################################################

#================================================
# Clear objects from the workspace ----

rm(list = ls())

#================================================
# Set working directory ----

# Maja's desktop PC

setwd("C:/Users/3054311/Documents/My Documents/Grassland project/05_Resilience")

data.dir <- paste0(getwd(),"/Data/")
figures.dir <- paste0(getwd(),"/Figures/SOIL MOISTURE/")
script.dir <- paste0(getwd(),"/Script/")

#================================================
# Packages ----

library(ggplot2)
library(lubridate)
library(dplyr)
library(tidyr)
library(stringr)
library(lme4)
library(rcompanion)
library(car)
library(performance)
library(insight)
library(see)
library(grid)
library(cowplot)
library(patchwork)
library(gridExtra)
library(effects)
library(afex)
library(RColorBrewer)

#================================================
# Set directory for data and plots ----

dir.create(paste0(figures.dir,"LMER model soil moisture day trt"))

mydir.data <- data.dir
mydir <- paste0(figures.dir,"LMER model soil moisture day trt")

#================================================
# Get raw data for soil moisture ----

raw.data.dir <- "C:/Users/3054311/Documents/My Documents/Grassland project/01_Probes_soil_moisture_and_temp/Data/"

load(paste0(raw.data.dir,"Cleaned_Soil_Moisture_Data_from_Loggers.RData"))

# Events 

load(file = paste0(raw.data.dir,"All events.RData"))

#================================================
# Calculate soil moisture ratio ----
# Idea from Paul: take soil moisture values within a certain time window
# For days 8, 16, 32 and 84, take the soil moisture data since the last sampling event
# E.g. for day 8, take soil moisture between day 0 and day 8, for day 16 take soil moisture between day 8 and day 16 and so on
# For day 0, take the soil moisture from the entire period of the simulated drought
# Daytime chosen: 12:00 on each sampling day
# Most of the following code is copy + pasted from the function "moisture.mov.win"

Sys.setenv(tz='UTC')  ## change timezone to UTC

sites <- unique(df.moisture$Site)  ## extract all sites

daytime.hr <- 12  ## 12:00 noon

for (i in 1:length(sites)) {
  
  mysite <- df.moisture %>% 
    filter(Site == sites[i])
  
  ## Get all dates of sampling days starting at specific day time
  
  day0 <- as.POSIXct(mdy(events$EndDrought[events$Site == sites[i]])) + hours(daytime.hr)
  day8 <- as.POSIXct(mdy(events$Day8[events$Site == sites[i]])) + hours(daytime.hr)
  day16 <- as.POSIXct(mdy(events$Day16[events$Site == sites[i]])) + hours(daytime.hr)
  day32 <- as.POSIXct(mdy(events$Day32[events$Site == sites[i]])) + hours(daytime.hr)
  day64 <- as.POSIXct(mdy(events$Day64[events$Site == sites[i]])) + hours(daytime.hr)
  
  start.drought <- as.POSIXct(mdy(events$StartDrought[events$Site == sites[i]])) + hours(daytime.hr)
  drought.duration <- as.numeric(day0 - start.drought)
  
  ## Create intervals for each sampling day of specific duration
  
  day0.interval <- day0 - hours(0:(24*drought.duration-1))
  day8.interval <- day8 - hours(0:(24*8-1))
  day16.interval <- day16 - hours(0:(24*8-1))
  day32.interval <- day32 - hours(0:(24*16-1))
  day64.interval <- day64 - hours(0:(24*32-1))
    
  date.intervals <- c(day0.interval,
                      day8.interval,
                      day16.interval,
                      day32.interval,
                      day64.interval)
  
  moisture.site <- mysite %>% 
    filter(date.time %in% date.intervals) %>% 
    mutate(Day = case_when(date.time %in% day0.interval ~ 0,
                           date.time %in% day8.interval ~ 8,
                           date.time %in% day16.interval ~ 16,
                           date.time %in% day32.interval ~ 32,
                           date.time %in% day64.interval ~ 64))
  
  if (i == 1) {
    df.window <- moisture.site
  }
  
  if (i > 1) {
    df.window <- rbind(df.window,moisture.site)
  }
  
  rm(day0,day8,day16,day32,day64,
     day0.interval,day8.interval,day16.interval,day32.interval,day64.interval,
     date.intervals,mysite,moisture.site)
}

# Calculate the average soil moisture over the period

avg.window <- df.window %>% 
  group_by(Day, Site, Region, Farm, Plot, site_ID, treatment) %>% 
  summarize(Avg.Moisture = mean(Moisture, na.rm = T))

# Calculate Cmean (baseline) for each site

baseline.window <- avg.window %>% 
  filter(treatment == "C") %>% 
  group_by(Day, Site, Region, Farm, site_ID) %>% 
  summarize(control_mean = mean(Avg.Moisture, na.rm = T))

# Join datasets 

moisture.window <- left_join(avg.window, baseline.window)

# Calculate ratio of soil moisture in each plot to the average soil moisture in the control plots in the respective site

moisture.window$Ratio <- moisture.window$Avg.Moisture/moisture.window$control_mean

# Save data

write.table(moisture.window, paste0(mydir.data,"Soil moisture resilience.csv"), sep = ",", row.names = F)
save(moisture.window, file = paste0(mydir.data,"Soil_moisture_resilience.RData"))

#================================================
# Plot ratios per sampling day ----

moisture.window$Day.Title <- paste("Day",moisture.window$Day)
moisture.window$Day.Title <- as.factor(moisture.window$Day.Title)
moisture.window$Day.Title <- factor(moisture.window$Day.Title, 
                                    levels = c("Day 0", "Day 8", "Day 16", "Day 32", "Day 64"))

#================================================
# Drought treatment ----

drought.window <- moisture.window %>% 
  filter(treatment == "D")

ylab.D <- expression(paste(bold("Ratio"~D[bolditalic("i")]~":"~"Control"[bolditalic("mean")])))

title.D <- paste0("Soil moisture ratio in drought plots, starting at 12-00")
figure.title1 <- paste0(mydir,"/Soil moisture ratio in drought plots, starting at 12-00 V1.png")
figure.title2 <- paste0(mydir,"/Soil moisture ratio in drought plots, starting at 12-00 V2.png")

# Panels: Day

g1 <- ggplot(drought.window, aes(x = Site, y = Ratio)) +
  geom_boxplot(aes(color = Region, fill = Region),
               alpha = 0.5) +
  facet_wrap(~ Day.Title, nrow = 5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(0.4, "cm")) +
  labs(x = "Day", y = ylab.D, title = title.D,
       subtitle = "Day 0: Drought period\nAll other days: soil moisture between two sampling events") +
  scale_y_continuous(limits = c(0,2)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50")

ggsave(figure.title1,
       g1, width = 18, height = 24, units = "cm")

# Panels: Region ~ Farm

g2 <- ggplot(drought.window, aes(x = as.factor(Day), y = Ratio)) +
  geom_boxplot(aes(color = Region, fill = Region),
               alpha = 0.5) +
  facet_grid(Region ~ Farm) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        strip.text = element_text(size = 12),
        legend.position = "none",
        panel.spacing = unit(0.4, "cm")) +
  labs(x = "Day", y = ylab.D, title = title.D,
       subtitle = "Day 0: Drought period\nAll other days: soil moisture between two sampling events") +
  scale_y_continuous(limits = c(0,2)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50")

ggsave(figure.title2,
       g2, width = 24, height = 22, units = "cm")

#================================================
# Control treatment ----

control.window <- moisture.window %>% 
  filter(treatment == "C")

ylab.C <- expression(paste(bold("Ratio"~C[bolditalic("i")]~":"~"Control"[bolditalic("mean")])))

title.C <- paste0("Soil moisture ratio in control plots, starting at 12-00")
figure.title3 <- paste0(mydir,"/Soil moisture ratio in control plots, starting at 12-00 V1.png")
figure.title4 <- paste0(mydir,"/Soil moisture ratio in control plots, starting at 12-00 V2.png")

# Panels: Day

g3 <- ggplot(control.window, aes(x = Site, y = Ratio)) +
  geom_boxplot(aes(color = Region, fill = Region),
               alpha = 0.5) +
  facet_wrap(~ Day.Title, nrow = 5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(0.4, "cm")) +
  labs(x = "Day", y = ylab.C, title = title.C,
       subtitle = "Day 0: Drought period\nAll other days: soil moisture between two sampling events") +
  scale_y_continuous(limits = c(0,2)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50")

ggsave(figure.title3,
       g3, width = 18, height = 24, units = "cm")

# Panels: Region ~ Farm

g4 <- ggplot(control.window, aes(x = as.factor(Day), y = Ratio)) +
  geom_boxplot(aes(color = Region, fill = Region),
               alpha = 0.5) +
  facet_grid(Region ~ Farm) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        strip.text = element_text(size = 12),
        legend.position = "none",
        panel.spacing = unit(0.4, "cm")) +
  labs(x = "Day", y = ylab.C, title = title.C,
       subtitle = "Day 0: Drought period\nAll other days: soil moisture between two sampling events") +
  scale_y_continuous(limits = c(0,2)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50")

ggsave(figure.title4,
       g4, width = 24, height = 22, units = "cm")

#================================================
# Functions for plots ----

# Histogram with Gaussian curve, Shapiro's test p-value, median, mean, variance and SD of x

source("C:/Users/3054311/Documents/My Documents/Grassland project/Functions/Function ggHistNorm.R")

## qqPlot made with ggplot

source("C:/Users/3054311/Documents/My Documents/Grassland project/Functions/Function ggQQplot.R")

#================================================
# Prepare data for the model

resilience_soil <- moisture.window[,c(1,3,4,5,7,10)]

resilience_soil$treatment <- gsub("C", "Control", 
                                as.character(resilience_soil$treatment))

resilience_soil$treatment <- gsub("D", "Drought", 
                                as.character(resilience_soil$treatment))

resilience_soil$treatment <- as.factor(resilience_soil$treatment)
resilience_soil$Region <- as.factor(resilience_soil$Region)
resilience_soil$Farm <- as.factor(resilience_soil$Farm)

resilience_soil <- resilience_soil %>% rename(plot_meanC_ratio = Ratio)

resilience_soil <- resilience_soil[!is.na(resilience_soil$plot_meanC_ratio),]

# Define response variable 

resp.variable <- "Soil moisture ratio"

#================================================
# Model assumptions ----

# Check normality of the data (soil moisture ratio) ----

check.norm <- ggHistNorm(resilience_soil, resilience_soil$plot_meanC_ratio, resp.variable, 0.05)
plot.norm <- check.norm[1][[1]]

plot.qq <- ggQQplot(resilience_soil, resilience_soil$plot_meanC_ratio, resp.variable)

# Treatment Control

C_resilience <- resilience_soil %>% filter(treatment == "Control")

check.norm.C <- ggHistNorm(C_resilience, C_resilience$plot_meanC_ratio, paste(resp.variable,"- Control"), 0.05)
plot.norm.C <- check.norm.C[1][[1]]

plot.qq.C <- ggQQplot(C_resilience, C_resilience$plot_meanC_ratio, paste(resp.variable,"- Control"))

# Treatment Drought

D_resilience <- resilience_soil %>% filter(treatment == "Drought")

check.norm.D <- ggHistNorm(D_resilience, D_resilience$plot_meanC_ratio, paste(resp.variable,"- Drought"), 0.05)
plot.norm.D <- check.norm.D[1][[1]]

plot.qq.D <- ggQQplot(D_resilience, D_resilience$plot_meanC_ratio, paste(resp.variable,"- Drought"))

## Check variance homogeneity ----

ylab <- expression(paste(bold("Ratio"~"Plot"[bolditalic("i")]~":"~"Control"[bolditalic("mean")])))

# ~ Day

OBS_var_day <- leveneTest(plot_meanC_ratio ~ as.factor(Day), data = resilience_soil)

if(OBS_var_day$`Pr(>F)`[1] < 0.001){
  main.day.0 <- "p < 0.001"
}

if(OBS_var_day$`Pr(>F)`[1] >= 0.001 & OBS_var_day$`Pr(>F)`[1] < 0.01){
  main.day.0 <- "p < 0.01"
}

if(OBS_var_day$`Pr(>F)`[1] >= 0.01 & OBS_var_day$`Pr(>F)`[1] < 0.05){
  main.day.0 <- "p < 0.05"
}

if(OBS_var_day$`Pr(>F)`[1] >= 0.05 ){
  main.day.0 <- parse(text = paste0('p == ', round(OBS_var_day$`Pr(>F)`[1], digits = 3)))
}

box.obs.day <- ggplot(resilience_soil, aes(x = as.factor(Day), 
                                           y = plot_meanC_ratio, 
                                           fill = as.factor(Day))) +
  geom_boxplot(outlier.shape = 21) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12, face = "bold")) +
  scale_fill_brewer(palette = "BuPu") +
  labs(x = "\nDays since the end of the simulated drought",
       y = ylab,
       title = "Soil moisture ratio ~ Day",
       subtitle = main.day.0)

# ~ treatment

OBS_var_trt <- leveneTest(plot_meanC_ratio ~ treatment, data = resilience_soil)

if(OBS_var_trt$`Pr(>F)`[1] < 0.001){
  main.trt.0 <- "p < 0.001"
}

if(OBS_var_trt$`Pr(>F)`[1] >= 0.001 & OBS_var_trt$`Pr(>F)`[1] < 0.01){
  main.trt.0 <- "p < 0.01"
}

if(OBS_var_trt$`Pr(>F)`[1] >= 0.01 & OBS_var_trt$`Pr(>F)`[1] < 0.05){
  main.trt.0 <- "p < 0.05"
}

if(OBS_var_trt$`Pr(>F)`[1] >= 0.05 ){
  main.trt.0 <- parse(text = paste0('p == ', round(OBS_var_trt$`Pr(>F)`[1], digits = 3)))
}

box.obs.trt <- ggplot(resilience_soil, aes(x = treatment, y = plot_meanC_ratio, fill = treatment)) +
  geom_boxplot(outlier.shape = 21) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12, face = "bold")) +
  scale_fill_manual(values = c("blue","red")) +
  labs(x = "\nTreatment",
       y = ylab,
       title = "Soil moisture ratio ~ treatment",
       subtitle = main.trt.0)

# ~ Region

OBS_var_region <- leveneTest(plot_meanC_ratio ~ Region, data = resilience_soil)

if(OBS_var_region$`Pr(>F)`[1] < 0.001){
  main.region.0 <- "p < 0.001"
}

if(OBS_var_region$`Pr(>F)`[1] >= 0.001 & OBS_var_region$`Pr(>F)`[1] < 0.01){
  main.region.0 <- "p < 0.01"
}

if(OBS_var_region$`Pr(>F)`[1] >= 0.01 & OBS_var_region$`Pr(>F)`[1] < 0.05){
  main.region.0 <- "p < 0.05"
}

if(OBS_var_region$`Pr(>F)`[1] >= 0.05 ){
  main.region.0 <- parse(text = paste0('p == ', round(OBS_var_region$`Pr(>F)`[1], digits = 3)))
}

box.obs.region <- ggplot(resilience_soil, aes(x = Region, y = plot_meanC_ratio, fill = Region)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.7) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12, face = "bold")) +
  labs(x = "\nRegion",
       y = ylab,
       title = "Soil moisture ratio ~ Region",
       subtitle = main.region.0)

# ~ Farm

OBS_var_farm <- leveneTest(plot_meanC_ratio ~ Farm, data = resilience_soil)

if(OBS_var_farm$`Pr(>F)`[1] < 0.001){
  main.farm.0 <- "p < 0.001"
}

if(OBS_var_farm$`Pr(>F)`[1] >= 0.001 & OBS_var_farm$`Pr(>F)`[1] < 0.01){
  main.farm.0 <- "p < 0.01"
}

if(OBS_var_farm$`Pr(>F)`[1] >= 0.01 & OBS_var_farm$`Pr(>F)`[1] < 0.05){
  main.farm.0 <- "p < 0.05"
}

if(OBS_var_farm$`Pr(>F)`[1] >= 0.05 ){
  main.farm.0 <- parse(text = paste0('p == ', round(OBS_var_farm$`Pr(>F)`[1], digits = 3)))
}

box.obs.farm <- ggplot(resilience_soil, aes(x = Farm, y = plot_meanC_ratio, fill = Farm)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.7) +
  scale_fill_brewer(palette = "GnBu") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12, face = "bold")) +
  labs(x = "\nFarm",
       y = ylab,
       title = "Soil moisture ratio ~ Farm",
       subtitle = main.farm.0)

# ~ Region * Farm

OBS_var_region_farm <- leveneTest(plot_meanC_ratio ~ Region*Farm, data = resilience_soil)

if(OBS_var_region_farm$`Pr(>F)`[1] < 0.001){
  main.region.farm.0 <- "p < 0.001"
}

if(OBS_var_region_farm$`Pr(>F)`[1] >= 0.001 & OBS_var_region_farm$`Pr(>F)`[1] < 0.01){
  main.region.farm.0 <- "p < 0.01"
}

if(OBS_var_region_farm$`Pr(>F)`[1] >= 0.01 & OBS_var_region_farm$`Pr(>F)`[1] < 0.05){
  main.region.farm.0 <- "p < 0.05"
}

if(OBS_var_region_farm$`Pr(>F)`[1] >= 0.05 ){
  main.region.farm.0 <- parse(text = paste0('p == ', round(OBS_var_region_farm$`Pr(>F)`[1], digits = 3)))
}

box.obs.region.farm <- ggplot(resilience_soil, aes(x = Region, y = plot_meanC_ratio, fill = Farm)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.7) +
  scale_fill_brewer(palette = "GnBu") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12, face = "bold")) +
  labs(x = "\nRegion",
       y = ylab,
       title = "Soil moisture ratio ~ Region * Farm",
       subtitle = main.region.farm.0)

#================================================
# Plot model assumptions ----
# Normality

title.model <- "lmer(plot_meanC_ratio ~ Day*treatment + (1|Day) + (1|Region/Farm))"

figure.title1 <- paste0(mydir,"/Model assumptions - Soil moisture ratio vs Day and Treatment V1.png")

final.plot1 <- grid.arrange(plot.norm,
                            plot.qq,
                            plot.norm.C,
                            plot.qq.C,
                            plot.norm.D,
                            plot.qq.D) 

final.plot1 <- grid.arrange(final.plot1,
                            top = tableGrob(t(title.model),
                                            theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                   base_colour = "black",
                                                                   base_size = 16)))

ggsave(figure.title1,
       grid.arrange(arrangeGrob(final.plot1,
                                top = textGrob(t(resp.variable),
                                               gp = gpar(fontsize = 18, fontface = "bold")))),
       
       width = 25, height = 32, units = "cm")

# Homogeneity

figure.title2 <- paste0(mydir,"/Model assumptions - Soil moisture ratio vs Day and Treatment V2.png")

final.plot2 <- grid.arrange(arrangeGrob(box.obs.day,
                                        box.obs.trt,
                                        box.obs.region,
                                        box.obs.farm, ncol = 2),
                            box.obs.region.farm, nrow = 2,
                            heights = c(2.2,1)) 

final.plot2 <- grid.arrange(final.plot2,
                            top = tableGrob(t(title.model),
                                            theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                   base_colour = "black",
                                                                   base_size = 16)))

ggsave(figure.title2,
       grid.arrange(arrangeGrob(final.plot2,
                                top = textGrob(t(resp.variable),
                                               gp = gpar(fontsize = 18, fontface = "bold")))),
       
       width = 25, height = 32, units = "cm")

#================================================
# Fit a model ----

# Fit a model of the ratio of the plot to the mean of the control plots as 
# a function Day since end of drought (numeric) and treatment (Control vs Drought)
# Random factors included: Day (1|day) and Farm nested in Region (1|Region/Farm)

lmer_full_mod <- lmer(plot_meanC_ratio ~ Day*treatment + (1|Day) + (1|Region/Farm),
                      data = resilience_soil)

#================================================
# Add fitted values and residuals to the raw data ----

resilience_mod <- data.frame(resilience_soil, "Fitted" = fitted(lmer_full_mod))  
resilience_mod$Residuals <- residuals(lmer_full_mod)

#================================================
# Run anova() and extract the results ----

aov_mod <- as.data.frame(anova(lmer_full_mod))
aov_mod$Term <- rownames(aov_mod)
aov_mod <- aov_mod[,c(7,1:6)]
names(aov_mod)[2:7] <- c("Sum.Sq","Mean.Sq","NumDF","DenDF","F.value","p.value")
rownames(aov_mod) <- 1:nrow(aov_mod)

#================================================
# Extract model coefficients ----

# Fixed

coeff_fixef <- summary(lmer_full_mod)$coefficients
rownames_fixef <- data.frame("Term" = rownames(coeff_fixef))
dimnames(coeff_fixef)[[2]] <- c("Estimate","Std.Error","DF","t.value","p.value")
fixef_out <- cbind(rownames_fixef, coeff_fixef)
rownames(fixef_out) <- rownames(rownames_fixef)

# Random: Farm:Region (Farm nested in Region) ----

coeff_random.1 <- ranef(lmer_full_mod)$`Farm:Region`
rownames_random.1 <- data.frame(rownames(coeff_random.1))
random_out.1 <- rownames_random.1 %>% separate(rownames.coeff_random.1., c("Farm", "Region"), ":")
random_out.1 <- data.frame(random_out.1, "Intercept" = coeff_random.1$`(Intercept)`)
random_out.1$Intercept_C <- random_out.1$Intercept + fixef_out$Estimate[fixef_out$Term == "(Intercept)"] 
random_out.1$Intercept_D <- random_out.1$Intercept + fixef_out$Estimate[fixef_out$Term == "(Intercept)"] + fixef_out$Estimate[fixef_out$Term == "treatmentDrought"]

# Random: Region ----

coeff_random.2 <- ranef(lmer_full_mod)$Region
rownames_random.2 <- rownames(coeff_random.2)
random_out.2 <- data.frame("Region" = rownames_random.2, "Intercept" = coeff_random.2$`(Intercept)`)
random_out.2$Intercept_1 <- random_out.2$Intercept + fixef_out$Estimate[fixef_out$Term == "(Intercept)"] 

# Random: Day ----

coeff_random.3 <- ranef(lmer_full_mod)$Day
rownames_random.3 <- rownames(coeff_random.3)
random_out.3 <- data.frame("day" = rownames_random.3, "Intercept" = coeff_random.3$`(Intercept)`)
random_out.3$Intercept_1 <- random_out.3$Intercept + fixef_out$Estimate[fixef_out$Term == "(Intercept)"] 

#================================================
# Model validation ----

# Observed vs fitted values ----

# ~ Day

obs.fit.day <- ggplot(resilience_mod, aes(x = Fitted, y = plot_meanC_ratio, fill = as.factor(Day))) +
  geom_point(shape = 21, size = 2, alpha = 0.5, color = "grey50") +
  facet_wrap(~ Day, ncol = 5) +
  geom_abline(slope = 1, intercept = 0) +
  theme_minimal() +
  theme(strip.background = element_rect(color = "grey50"),
        legend.position = "none",
        axis.title = element_text(size = 12, face = "bold")) +
  scale_fill_brewer(palette = "BuPu") +
  labs(x = "\nFitted",
       y = ylab,
       title = "Observed vs. Fitted",
       subtitle = "~ Day")

# ~ treatment

obs.fit.trt <- ggplot(resilience_mod, aes(x = Fitted, y = plot_meanC_ratio, 
                                          color = treatment, fill = treatment)) +
  geom_point(shape = 21, size = 2, alpha = 0.5) +
  facet_wrap(~ treatment, ncol = 2) +
  geom_abline(slope = 1, intercept = 0) +
  theme_minimal() +
  theme(strip.background = element_rect(color = "grey50"),
        legend.position = "none",
        axis.title = element_text(size = 12, face = "bold")) +
  scale_fill_manual(values = c("blue","red")) +
  scale_color_manual(values = c("blue","red")) +
  labs(x = "\nFitted",
       y = ylab,
       title = "Observed vs. Fitted",
       subtitle = "~ treatment")

# ~ Region

obs.fit.region <- ggplot(resilience_mod, aes(x = Fitted, y = plot_meanC_ratio, 
                                             color = Region, fill = Region)) +
  geom_point(shape = 21, size = 2, alpha = 0.5) +
  facet_wrap(~ Region, ncol = 2) +
  geom_abline(slope = 1, intercept = 0) +
  theme_minimal() +
  theme(strip.background = element_rect(color = "grey50"),
        legend.position = "none",
        axis.title = element_text(size = 12, face = "bold")) +
  labs(x = "\nFitted",
       y = ylab,
       title = "Observed vs. Fitted",
       subtitle = "~ Region")

# ~ Farm

obs.fit.farm <- ggplot(resilience_mod, aes(x = Fitted, y = plot_meanC_ratio, fill = Farm)) +
  geom_point(shape = 21, size = 2, alpha = 0.5, color = "grey50") +
  scale_fill_brewer(palette = "GnBu") +
  facet_wrap(~ Farm, ncol = 3) +
  geom_abline(slope = 1, intercept = 0) +
  theme_minimal() +
  theme(strip.background = element_rect(color = "grey50"),
        legend.position = "none",
        axis.title = element_text(size = 12, face = "bold")) +
  labs(x = "\nFitted",
       y = ylab,
       title = "Observed vs. Fitted",
       subtitle = "~ Farm")

# ~ Region * Farm

obs.fit.region.farm <- ggplot(resilience_mod, aes(x = Fitted, y = plot_meanC_ratio, 
                                          color = Region, fill = Region)) +
  geom_point(shape = 21, size = 2, alpha = 0.5) +
  facet_grid(Region ~ Farm) +
  geom_abline(slope = 1, intercept = 0) +
  theme_minimal() +
  theme(strip.background = element_rect(color = "grey50"),
        legend.position = "none",
        axis.title = element_text(size = 12, face = "bold")) +
  labs(x = "\nFitted",
       y = ylab,
       title = "Observed vs. Fitted",
       subtitle = "~ Region * Farm")

#================================================
# Collinearity ----

result.coll <- check_collinearity(lmer_full_mod)

result.coll[which(result.coll$VIF < 5),"Correlation"] <- "low"
result.coll[which(result.coll$VIF >= 5 & result.coll$VIF < 10),"Correlation"] <- "moderate"
result.coll[which(result.coll$VIF >= 10),"Correlation"] <- "high"

result.coll$Correlation <- as.factor(result.coll$Correlation)
result.coll$Correlation <- factor(result.coll$Correlation, levels = c("low", "moderate", "high"))

mycol.coll <- c(rgb(39, 174, 96, max = 255),
                rgb(230, 126, 34, max = 255),
                rgb(228, 26, 28, max = 255))

corr.levels <- unique(result.coll$Correlation)

mycol <- c()

if ("low" %in% corr.levels){
  mycol <- mycol.coll[1]
}

if ("moderate" %in% corr.levels){
  mycol <- c(mycol,mycol.coll[2])
}

if ("high" %in% corr.levels){
  mycol <- c(mycol,mycol.coll[3])
}

plot.coll <- ggplot(result.coll, aes(x = Parameter, y = VIF)) +
  geom_bar(stat = "identity", width = 0.7, aes(fill = Correlation)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_fill_manual(values = mycol) +
  labs(x = "\nParameter", y = "VIF", title = "Check for Multicollinearity", subtitle = "")

#================================================
# Residuals ----
# Binned residuals
# Function binned_residuals produces a plot if not saved to an object
# However, that plot can't be further modified
# Therefore, I save the output to an object and recreated the plot 

result.binned <- binned_residuals(lmer_full_mod)
resid.inside.err <- sum(result.binned$group == "yes")/nrow(result.binned)*100
resid.inside.err <- round(resid.inside.err, digits = 2)

plot.binned <- ggplot(result.binned, aes(x = xbar*100, y = ybar)) +
  geom_ribbon(aes(ymin = -Inf, ymax = -se),
              color = "grey80", fill = "grey95", alpha = 0.5) +
  geom_ribbon(aes(ymin = se, ymax = +Inf),
              color = "grey80", fill = "grey95", alpha = 0.5) +
  geom_hline(yintercept = 0, color = "grey80") +
  geom_point(aes(color = group), size = 3) +
  theme_bw() +
  scale_color_brewer(palette = "Set1") + 
  labs(x = paste0("\nEstimated probability of ", resp.variable), 
       y = "Average residual", 
       title = "Binned residuals", 
       subtitle = paste0(resid.inside.err, "% of the residuals are inside the error bounds."))

# ~ Fitted

result.heteroscedasticity <- check_heteroscedasticity(lmer_full_mod)

if(result.heteroscedasticity[1] < 0.001){
  p.res.fit <- "p < 0.001"
}

if(result.heteroscedasticity[1] >= 0.001 & result.heteroscedasticity[1] < 0.01){
  p.res.fit <- "p < 0.01"
}

if(result.heteroscedasticity[1] >= 0.01 & result.heteroscedasticity[1] < 0.05){
  p.res.fit <- "p < 0.05"
}

if(result.heteroscedasticity[1] >= 0.05 ){
  p.res.fit <- parse(text = paste0('p == ', round(result.heteroscedasticity[1], digits = 3)))
}

plot.res.fit <- ggplot(resilience_mod, aes(x = Fitted, y = Residuals)) +
  geom_point(size = 2, color = rgb(44, 62, 80, max = 255)) +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1, color = rgb(228, 26, 28, max = 255), se = F) +
  labs(x = "\nFitted",
       y = "Residuals",
       title = "Residuals vs. Fitted",
       subtitle = p.res.fit)

# ~ Day

RESID_var_day <- leveneTest(Residuals ~ as.factor(Day), data = resilience_mod)

if(RESID_var_day$`Pr(>F)`[1] < 0.001){
  main.day <- "p < 0.001"
}

if(RESID_var_day$`Pr(>F)`[1] >= 0.001 & RESID_var_day$`Pr(>F)`[1] < 0.01){
  main.day <- "p < 0.01"
}

if(RESID_var_day$`Pr(>F)`[1] >= 0.01 & RESID_var_day$`Pr(>F)`[1] < 0.05){
  main.day <- "p < 0.05"
}

if(RESID_var_day$`Pr(>F)`[1] >= 0.05 ){
  main.day <- parse(text = paste0('p == ', round(RESID_var_day$`Pr(>F)`[1], digits = 3)))
}

box.day <- ggplot(resilience_mod, aes(x = as.factor(Day), 
                                      y = Residuals, 
                                      fill = as.factor(Day))) +
  geom_boxplot(outlier.shape = 21) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12, face = "bold")) +
  scale_fill_brewer(palette = "BuPu") +
  labs(x = "\nDays since the end of the simulated drought",
       y = "Residuals",
       title = "Residuals ~ Day",
       subtitle = main.day)

# ~ treatment

RESID_var_trt <- leveneTest(Residuals ~ treatment, data = resilience_mod)

if(RESID_var_trt$`Pr(>F)`[1] < 0.001){
  main.trt <- "p < 0.001"
}

if(RESID_var_trt$`Pr(>F)`[1] >= 0.001 & RESID_var_trt$`Pr(>F)`[1] < 0.01){
  main.trt <- "p < 0.01"
}

if(RESID_var_trt$`Pr(>F)`[1] >= 0.01 & RESID_var_trt$`Pr(>F)`[1] < 0.05){
  main.trt <- "p < 0.05"
}

if(RESID_var_trt$`Pr(>F)`[1] >= 0.05 ){
  main.trt <- parse(text = paste0('p == ', round(RESID_var_trt$`Pr(>F)`[1], digits = 3)))
}

box.trt <- ggplot(resilience_mod, aes(x = treatment, y = Residuals, fill = treatment)) +
  geom_boxplot(outlier.shape = 21) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12, face = "bold")) +
  scale_fill_manual(values = c("blue","red")) +
  labs(x = "\ntreatment",
       y = "Residuals",
       title = "Residuals ~ treatment",
       subtitle = main.trt)

# ~ Region

RESID_var_region <- leveneTest(Residuals ~ Region, data = resilience_mod)

if(RESID_var_region$`Pr(>F)`[1] < 0.001){
  main.region <- "p < 0.001"
}

if(RESID_var_region$`Pr(>F)`[1] >= 0.001 & RESID_var_region$`Pr(>F)`[1] < 0.01){
  main.region <- "p < 0.01"
}

if(RESID_var_region$`Pr(>F)`[1] >= 0.01 & RESID_var_region$`Pr(>F)`[1] < 0.05){
  main.region <- "p < 0.05"
}

if(RESID_var_region$`Pr(>F)`[1] >= 0.05 ){
  main.region <- parse(text = paste0('p == ', round(RESID_var_region$`Pr(>F)`[1], digits = 3)))
}

box.region <- ggplot(resilience_mod, aes(x = Region, y = Residuals, fill = Region)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.7) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12, face = "bold")) +
  labs(x = "\nRegion",
       y = "Residuals",
       title = "Residuals ~ Region",
       subtitle = main.region)

# ~ Farm

RESID_var_farm <- leveneTest(Residuals ~ Farm, data = resilience_mod)

if(RESID_var_farm$`Pr(>F)`[1] < 0.001){
  main.farm <- "p < 0.001"
}

if(RESID_var_farm$`Pr(>F)`[1] >= 0.001 & RESID_var_farm$`Pr(>F)`[1] < 0.01){
  main.farm <- "p < 0.01"
}

if(RESID_var_farm$`Pr(>F)`[1] >= 0.01 & RESID_var_farm$`Pr(>F)`[1] < 0.05){
  main.farm <- "p < 0.05"
}

if(RESID_var_farm$`Pr(>F)`[1] >= 0.05 ){
  main.farm <- parse(text = paste0('p == ', round(RESID_var_farm$`Pr(>F)`[1], digits = 3)))
}

box.farm <- ggplot(resilience_mod, aes(x = Farm, y = Residuals, fill = Farm)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.7) +
  scale_fill_brewer(palette = "GnBu") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12, face = "bold")) +
  labs(x = "\nFarm",
       y = "Residuals",
       title = "Residuals ~ Farm",
       subtitle = main.farm)

# ~ Region * Farm

RESID_var_region_farm <- leveneTest(Residuals ~ Region*Farm, data = resilience_mod)

if(RESID_var_region_farm$`Pr(>F)`[1] < 0.001){
  main.region.farm <- "p < 0.001"
}

if(RESID_var_region_farm$`Pr(>F)`[1] >= 0.001 & RESID_var_region_farm$`Pr(>F)`[1] < 0.01){
  main.region.farm <- "p < 0.01"
}

if(RESID_var_region_farm$`Pr(>F)`[1] >= 0.01 & RESID_var_region_farm$`Pr(>F)`[1] < 0.05){
  main.region.farm <- "p < 0.05"
}

if(RESID_var_region_farm$`Pr(>F)`[1] >= 0.05 ){
  main.region.farm <- parse(text = paste0('p == ', round(RESID_var_region_farm$`Pr(>F)`[1], digits = 3)))
}

box.region.farm <- ggplot(resilience_mod, aes(x = Region, y = Residuals, fill = Farm)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.7) +
  scale_fill_brewer(palette = "GnBu") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12, face = "bold")) +
  labs(x = "\nRegion",
       y = "Residuals",
       title = "Residuals ~ Region * Farm",
       subtitle = main.region.farm)

#================================================
# Check for normal distribution of residuals ----

check.resid.norm <- ggHistNorm(resilience_mod, resilience_mod$Residuals, "Residuals", 0.05)
plot.norm.resid <- check.resid.norm[1][[1]]

plot.qq.resid <- ggQQplot(resilience_mod, resilience_mod$Residuals, "Residuals")

# Treatment Control

C_mod <- resilience_mod %>% filter(treatment == "Control")

check.resid.norm.C <- ggHistNorm(C_mod, C_mod$Residuals, "Residuals - Control", 0.05)
plot.norm.resid.C <- check.resid.norm.C[1][[1]]

plot.qq.resid.C <- ggQQplot(C_mod, C_mod$Residuals, "Residuals - Control")

# Treatment Control

D_mod <- resilience_mod %>% filter(treatment == "Drought")

check.resid.norm.D <- ggHistNorm(D_mod, D_mod$Residuals, "Residuals - Drought", 0.05)
plot.norm.resid.D <- check.resid.norm.D[1][[1]]

plot.qq.resid.D <- ggQQplot(D_mod, D_mod$Residuals, "Residuals - Drought")

#================================================
# Normality of Random Effects ----

result.mod <- check_model(lmer_full_mod)

# ~ day

REQQ.day <- result.mod$REQQ$Day

reqq.plot.day <- ggplot(REQQ.day, aes(x = x, y = y)) +
  geom_point(size = 2, color = rgb(44, 62, 80, max = 255)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) + 
  geom_smooth(method = "lm", color = rgb(22, 160, 133, max = 255), size = 1, se = F) +
  theme_minimal() +
  labs(x = "Theoretical Quantiles",
       y = "RE Quantiles",
       title = "Normality of Random Effects",
       subtitle = "Day")

# ~ Region

REQQ.region <- result.mod$REQQ$Region

reqq.plot.region <- ggplot(REQQ.region, aes(x = x, y = y)) +
  geom_point(size = 2, color = rgb(44, 62, 80, max = 255)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) + 
  geom_smooth(method = "lm", color = rgb(22, 160, 133, max = 255), size = 1, se = F) +
  theme_minimal() +
  labs(x = "Theoretical Quantiles",
       y = "RE Quantiles",
       title = "Normality of Random Effects",
       subtitle = "Region")

# ~ Farm:Region

REQQ.farm.region <- result.mod$REQQ$"Farm:Region"

reqq.plot.farm.region <- ggplot(REQQ.farm.region, aes(x = x, y = y)) +
  geom_point(size = 2, color = rgb(44, 62, 80, max = 255)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) + 
  geom_smooth(method = "lm", color = rgb(22, 160, 133, max = 255), size = 1, se = F) +
  theme_minimal() +
  labs(x = "Theoretical Quantiles",
       y = "RE Quantiles",
       title = "Normality of Random Effects",
       subtitle = "Farm:Region")

#================================================
# Plot observed vs. fitted ----

figure.title3 <- paste0(mydir,"/Observed vs Fitted V1.png")

final.plot3 <- grid.arrange(obs.fit.day,
                            obs.fit.trt,
                            obs.fit.region,
                            obs.fit.farm, 
                            ncol = 2,
                            top = tableGrob(t(title.model),
                                            theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                   base_colour = "black",
                                                                   base_size = 16)))

ggsave(figure.title3,
       grid.arrange(arrangeGrob(final.plot3,
                                top = textGrob(t(resp.variable),
                                               gp = gpar(fontsize = 18, fontface = "bold")))),
       width = 25, height = 27, units = "cm")

# Farm:Region

obs.fit.region.farm <- grid.arrange(obs.fit.region.farm,
                                    top = tableGrob(t(title.model),
                                                    theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                           base_colour = "black",
                                                                           base_size = 16)))

ggsave(paste0(mydir,"/Observed vs Fitted V2.png"),
       grid.arrange(arrangeGrob(obs.fit.region.farm,
                                top = textGrob(t(resp.variable),
                                               gp = gpar(fontsize = 18, fontface = "bold")))),
       width = 20, height = 20, units = "cm")

#================================================
# Plot residuals vs. fitted ----

figure.title4 <- paste0(mydir,"/Residuals V1.png")

final.plot4 <- grid.arrange(plot.norm.resid,
                            plot.qq.resid,
                            plot.norm.resid.C,
                            plot.qq.resid.C,
                            plot.norm.resid.D,
                            plot.qq.resid.D,
                            ncol = 2,
                            top = tableGrob(t(title.model),
                                            theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                   base_colour = "black",
                                                                   base_size = 16)))

ggsave(figure.title4,
       grid.arrange(arrangeGrob(final.plot4,
                                top = textGrob(t(resp.variable),
                                               gp = gpar(fontsize = 18, fontface = "bold")))),
       width = 25, height = 32, units = "cm")

#================================================
# Plot residuals vs. random effects ----

figure.title5 <- paste0(mydir,"/Residuals V2.png")

final.plot5 <- grid.arrange(arrangeGrob(box.day,
                                        box.trt,
                                        box.region,
                                        box.farm, ncol = 2),
                            box.region.farm, nrow = 2,
                            heights = c(2.2,1),
                            top = tableGrob(t(title.model),
                                            theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                   base_colour = "black",
                                                                   base_size = 16)))

ggsave(figure.title5,
       grid.arrange(arrangeGrob(final.plot5,
                                top = textGrob(t(resp.variable),
                                               gp = gpar(fontsize = 18, fontface = "bold")))),
       width = 25, height = 32, units = "cm")

#================================================
# Plot residuals and collinearity ----

figure.title6 <- paste0(mydir,"/Residuals V3.png")

final.plot6 <- grid.arrange(plot.res.fit,
                            plot.binned,
                            reqq.plot.day,
                            reqq.plot.region,
                            reqq.plot.farm.region,
                            plot.coll,
                            ncol = 2,
                            top = tableGrob(t(title.model),
                                            theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                   base_colour = "black",
                                                                   base_size = 16)))

ggsave(figure.title6,
       grid.arrange(arrangeGrob(final.plot6,
                                top = textGrob(t(resp.variable),
                                               gp = gpar(fontsize = 18, fontface = "bold")))),
       width = 25, height = 32, units = "cm")

#================================================
# Plot data - Model output ----
# Panels: Region ~ Farm

# Raw data

g1 <- ggplot(resilience_mod, 
             aes(x = Day, y = plot_meanC_ratio)) +
  geom_point(aes(fill = treatment, color = treatment),
             size = 2, alpha = 0.2) +
  facet_grid(Region ~ Farm) +
  geom_smooth(aes(y = Fitted, color = treatment),
              method = "lm", se = F) +
  theme_minimal() +
  theme(strip.background = element_rect(),
        panel.spacing = unit(1, "lines"),
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 12, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 14, face = "bold")) +
  scale_fill_manual(values = c("blue","red")) +
  scale_color_manual(values = c("blue","red")) +
  scale_x_continuous(limits = c(0,64), breaks = c(0,8,16,32,64)) +
  labs(x = "Days since the end of the simulated drought",
       y = ylab,
       title = "Raw data",
       subtitle = "lm(fitted data ~ day)")

# Line plot

g2 <- ggplot(resilience_mod, 
             aes(x = Day, y = Fitted)) +
  facet_grid(Region ~ Farm) +
  geom_smooth(aes(color = treatment),
              method = "lm", se = F) +
  theme_minimal() +
  theme(strip.background = element_rect(),
        panel.spacing = unit(1, "lines"),
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 12, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 14, face = "bold")) +
  scale_fill_manual(values = c("blue","red")) +
  scale_color_manual(values = c("blue","red")) +
  scale_x_continuous(limits = c(0,64), breaks = c(0,8,16,32,64)) +
  labs(x = "Days since the end of the simulated drought",
       title = "",
       subtitle = "lm(fitted data ~ day)")

## Panels: Region

# Raw data

g3 <- ggplot(resilience_mod, 
             aes(x = Day, y = plot_meanC_ratio)) +
  geom_point(aes(fill = treatment, color = treatment),
             size = 3, alpha = 0.2) +
  facet_wrap(~ Region) +
  geom_smooth(aes(y = Fitted, color = treatment, linetype = Farm),
              method = "lm", se = F, size = 0.6) +
  geom_smooth(aes(y = Fitted, color = treatment),
              method = "lm", se = F, size = 1) +
  theme_minimal() +
  theme(strip.background = element_rect(),
        panel.spacing = unit(1, "lines"),
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 12, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 14, face = "bold")) +
  scale_fill_manual(values = c("blue","red")) +
  scale_color_manual(values = c("blue","red")) +
  scale_linetype_manual(values = c(2,3,4,5,6)) +
  scale_x_continuous(limits = c(0,64), breaks = c(0,8,16,32,64)) +
  labs(x = "Days since the end of the simulated drought",
       y = ylab,
       title = "Raw data",
       subtitle = "lm(fitted data ~ day)")

# Line plot

g4 <-  ggplot(resilience_mod, 
              aes(x = Day, y = Fitted)) +
  facet_wrap(~ Region) +
  geom_smooth(aes(color = treatment, linetype = Farm),
              method = "lm", se = F, size = 0.7) +
  geom_smooth(aes(color = treatment),
              method = "lm", se = F, size = 1) +
  theme_minimal() +
  theme(strip.background = element_rect(),
        panel.spacing = unit(1, "lines"),
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 12, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 14, face = "bold")) +
  scale_fill_manual(values = c("blue","red")) +
  scale_color_manual(values = c("blue","red")) +
  scale_linetype_manual(values = c(2,3,4,5,6)) +
  scale_x_continuous(limits = c(0,64), breaks = c(0,8,16,32,64)) +
  labs(x = "Days since the end of the drought",
       title = "",
       subtitle = "lm(fitted data ~ day)")

# All four plots together

figure.title7 <- paste0(mydir,"/Model output.png")

final.plot7 <- grid.arrange(g1,g3,g2,g4,
                            ncol = 2,
                            top = tableGrob(t(title.model),
                                            theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                   base_colour = "black",
                                                                   base_size = 16)))


ggsave(figure.title7,
       grid.arrange(arrangeGrob(final.plot7,
                                top = textGrob(t(resp.variable),
                                               gp = gpar(fontsize = 18, fontface = "bold")))),
       width = 40, height = 32, units = "cm")

# Only the last plot

figure.title8 <- paste0(mydir,"/Model output 2.png")

g4 <- g4 + labs(title = title.model)

ggsave(figure.title8, g4, width = 9, height = 7)

#================================================
# Statistics
 
AIC(lmer_full_mod)                                  ## -116.9684
model_performance(lmer_full_mod)$R2_marginal        ## 0.2085482
model_performance(lmer_full_mod)$R2_conditional     ## 0.3771294

# write.table(aov_mod, paste0(mydir.data, "/LMER Anova output - Soil moisture ratio vs Day and Trt.csv"), 
#             sep = ",", row.names = F)
# 
# write.table(fixef_out, paste0(mydir.data, "/LMER Model output fixef - Soil moisture ratio vs Day and Trt.csv"), 
#             sep = ",", row.names = F)
# 
# write.table(random_out.1, paste0(mydir.data, "/LMER Model output random Farm Region - Soil moisture ratio vs Day and Trt.csv"), 
#             sep = ",", row.names = F)
# 
# write.table(random_out.2, paste0(mydir.data, "/LMER Model output random Region - Soil moisture ratio vs Day and Trt.csv"), 
#             sep = ",", row.names = F)
# 
# write.table(random_out.3, paste0(mydir.data, "/LMER Model output random Day - Soil moisture ratio vs Day and Trt.csv"), 
#             sep = ",", row.names = F)
