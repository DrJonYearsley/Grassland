# Author: Maja Ilic
# Last update: November 2019

# This is a script to calculate soil moisture ratios Drought : Control from the Grassland Resilience 
# Year 1 All Ireland Drought Experiment. In the previous script (05), I calculated running ratio
# (ratio of drought to mean control at each time point, resolution: 1 hour). In this script, 
# I am extracting the soil mositure data from drought and control treatments on/around sampling days
# 0, 8, 16, 32 and 64. In the next step, I calculate again the ratio of drought to mean control.
# To get the ratio for each of the sampling days, three options were suggested:
# I. 12:00 noon on the day of sampling
# II. Average over 24 h, starting at 12:00 noon
# III. Average over 48 h, starting at 12:00 noon
# Starting at 00:00 (midnight) does not make sense for Day 0 (lids on till noon? Check with Lupe!) 

# See also email from Willson Gaul (UCD) on 18/11/2019

#===========================================
# Clear objects from the workspace

rm(list = ls())

#===========================================
# Set working directory

# Maja's desktop PC

setwd("C:/Users/3054311/Documents/My Documents/Grassland project/01_Probes_soil_moisture_and_temp")

data.dir <- paste0(getwd(),"/Data/")
figures.dir <- paste0(getwd(),"/Figures/")
script.dir <- paste0(getwd(),"/Script/")

#===========================================
# Packages

library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(lubridate)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)

#================================================
# Load data

# Soil moisture

load(paste0(data.dir,"Cleaned_Soil_Moisture_Data_from_Loggers.RData"))
df <- df.moisture   ## shorter
rm(df.moisture)
names(df)

# Events

load(file = paste0(data.dir,"All events.RData"))

#================================================
# Create new folder for data and plots

dir.create(paste0(data.dir,"Soil moisture ratios DtoC per site and sampling day"))
mydir.data <- paste0(data.dir,"Soil moisture ratios DtoC per site and sampling day/")

dir.create(paste0(figures.dir,"Soil moisture ratios DtoC per site and sampling day"))
mydir <- paste0(figures.dir,"Soil moisture ratios DtoC per site and sampling day/")

#================================================

##--------------------##
##                    ##
##     12:00 noon     ##
##                    ##
##--------------------##

# Extract soil moisture at 12:00 noon on all sampling days for each site

sites <- unique(df$Site)

for (i in 1:length(sites)) {
  
  mysite <- df %>% 
    filter(Site == sites[i])
  
  day0 <- mdy(events$EndDrought[events$Site == sites[i]])
  day8 <- mdy(events$Day8[events$Site == sites[i]])
  day16 <- mdy(events$Day16[events$Site == sites[i]])
  day32 <- mdy(events$Day32[events$Site == sites[i]])
  day64 <- mdy(events$Day64[events$Site == sites[i]])
  
  dates <- c(day0,day8,day16,day32,day64)
  
  moisture.site <- mysite %>% 
    filter(date(date.time) %in% dates) %>% 
    filter(hour(date.time) == 12) %>% 
    mutate(Day = case_when(date(date.time) == day0 ~ 0,
                           date(date.time) == day8 ~ 8,
                           date(date.time) == day16 ~ 16,
                           date(date.time) == day32 ~ 32,
                           date(date.time) == day64 ~ 64))
  
  if (i == 1) {
    df.12.noon <- moisture.site
  }
  
  if (i > 1) {
    df.12.noon <- rbind(df.12.noon,moisture.site)
  }
  
  rm(day0,day8,day16,day32,day64,
     dates,mysite,moisture.site)
  
}

## Save data

write.table(df.12.noon, paste0(mydir.data,"Soil moisture on sampling days - 12 noon.csv"), sep = ",", row.names = F)
save(df.12.noon, file = paste0(mydir.data,"Soil_moisture_on_sampling_days_12_noon.RData"))

#================================================
## Calculate ratio for each Drought replicate as Ratio = Di / Cmean for every sampling day and every site

# Calculate Cmean (baseline) for each site

baseline.12.noon <- df.12.noon %>% 
  filter(treatment == "C") %>% 
  group_by(Day, Site, Region, Farm, site_ID) %>% 
  summarize(control_mean = mean(Moisture, na.rm = T))

# Join datasets 

moisture.12.noon <- left_join(df.12.noon, baseline.12.noon)

# Calculate ratio

moisture.12.noon$Ratio <- moisture.12.noon$Moisture/moisture.12.noon$control_mean

#================================================
# Plot ratios per sampling day

moisture.12.noon$Day.Title <- paste("Day",moisture.12.noon$Day)
moisture.12.noon$Day.Title <- as.factor(moisture.12.noon$Day.Title)
moisture.12.noon$Day.Title <- factor(moisture.12.noon$Day.Title, 
                                     levels = c("Day 0", "Day 8", "Day 16", "Day 32", "Day 64"))

# Only drought treatment

drought.12.noon <- moisture.12.noon %>% 
  filter(treatment == "D")

ylab <- expression(paste(bold("Ratio"~D["i"]~":"~C["mean"])))

g1 <- ggplot(drought.12.noon, aes(x = Site, y = Ratio)) +
  geom_boxplot(aes(color = Region, fill = Region),
               alpha = 0.5) +
  facet_wrap(~ Day.Title, nrow = 5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(0.4, "cm")) +
  labs(x = "Date", y = ylab, title = "Soil moisture ratio 12:00 noon") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  scale_y_continuous(limits = c(0,3.5))

ggsave(paste0(mydir,"Soil moisture ratio 12-00 noon V1.png"),
       g1, width = 18, height = 24, units = "cm")

####

g2 <- ggplot(drought.12.noon, aes(x = as.factor(Day), y = Ratio)) +
  geom_boxplot(aes(color = Region, fill = Region),
               alpha = 0.5) +
  facet_grid(Region ~ Farm) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        strip.text = element_text(size = 12),
        legend.position = "none",
        panel.spacing = unit(0.4, "cm")) +
  labs(x = "Date", y = ylab, title = "Soil moisture ratio 12:00 noon") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  scale_y_continuous(limits = c(0,3.5))

ggsave(paste0(mydir,"Soil moisture ratio 12-00 noon V2.png"),
       g2, width = 24, height = 18, units = "cm")

#================================================
# Plot ratios per sampling day regardless of site (see Willson's script)

ylab2 <- expression(paste(bold("Ratio"~"Plot"["i"]~":"~"Plot"["control mean"])))

g3 <- ggplot(data = moisture.12.noon, 
             aes(x = Day, y = Ratio, 
                 color = treatment, fill = treatment)) + 
  geom_point() + 
  geom_smooth(method = "loess", alpha = 0.2) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0))) +
  labs(y = ylab2)

ggsave(paste0(mydir,"Soil moisture ratio 12-00 noon V3.png"), 
       g3, width = 8, height = 5)

#================================================
# Plot ratios per sampling day for each region and farm

g4 <- ggplot(data = moisture.12.noon, 
             aes(x = Day, y = Ratio, 
                 color = treatment, fill = treatment)) + 
  geom_point() + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        strip.text = element_text(size = 12),
        panel.spacing = unit(0.4, "cm")) +
  geom_smooth(method = "loess", alpha = 0.2) +
  facet_grid(Region ~ Farm) +
  labs(y = ylab2) +
  scale_y_continuous(limits = c(0,4))

ggsave(paste0(mydir,"Soil moisture ratio 12-00 noon V4.png"), 
       g4, width = 24, height = 18, units = "cm")

#================================================

##---------------------##
##                     ##
##     24 h period     ##
##                     ##
##---------------------##

# Extract soil moisture from 12:00 noon to 11:00 am next day
# on all sampling days for each site

Sys.setenv(tz='UTC')  ## change timezone to UTC

sites <- unique(df$Site)

for (i in 1:length(sites)) {
  
  mysite <- df %>% 
    filter(Site == sites[i])
  
  ## Get all dates of sampling days exactly at 12:00 noon
  
  day0 <- as.POSIXct(mdy(events$EndDrought[events$Site == sites[i]])) + hours(12)
  day8 <- as.POSIXct(mdy(events$Day8[events$Site == sites[i]])) + hours(12)
  day16 <- as.POSIXct(mdy(events$Day16[events$Site == sites[i]])) + hours(12)
  day32 <- as.POSIXct(mdy(events$Day32[events$Site == sites[i]])) + hours(12)
  day64 <- as.POSIXct(mdy(events$Day64[events$Site == sites[i]])) + hours(12)
  
  ## Create intervals of 24 hours for each sampling day
  # This also works with interval(day0, day0 + hours(23)) and using %within%, see lubridate cheatsheet 
  
  day0.24h <- day0 + hours(0:23)
  day8.24h <- day8 + hours(0:23)
  day16.24h <- day16 + hours(0:23)
  day32.24h <- day32 + hours(0:23)
  day64.24h <- day64 + hours(0:23)
  
  dates <- c(day0.24h,day8.24h,day16.24h,day32.24h,day64.24h)
  
  moisture.site <- mysite %>% 
    filter(date.time %in% dates) %>% 
    mutate(Day = case_when(date.time %in% day0.24h ~ 0,
                           date.time %in% day8.24h ~ 8,
                           date.time %in% day16.24h ~ 16,
                           date.time %in% day32.24h ~ 32,
                           date.time %in% day64.24h ~ 64))
  
  if (i == 1) {
    df.24h <- moisture.site
  }
  
  if (i > 1) {
    df.24h <- rbind(df.24h,moisture.site)
  }
  
  rm(day0,day8,day16,day32,day64,
     day0.24h,day8.24h,day16.24h,day32.24h,day64.24h,
     dates,mysite,moisture.site)
}

## Save data

write.table(df.24h, paste0(mydir.data,"Soil moisture on sampling days - 24 h.csv"), sep = ",", row.names = F)
save(df.24h, file = paste0(mydir.data,"Soil_moisture_on_sampling_days_24_h.RData"))

#================================================
# Calculate the average soil moisture over the 24 h period

avg.24h <- df.24h %>% 
  group_by(Day, Site, Region, Farm, Plot, site_ID, treatment) %>% 
  summarize(Avg.Moisture = mean(Moisture, na.rm = T))

## Calculate ratio for each Drought replicate as Ratio = Di / Cmean for every sampling day and every site

# Calculate Cmean (baseline) for each site

baseline.24h <- avg.24h %>% 
  filter(treatment == "C") %>% 
  group_by(Day, Site, Region, Farm, site_ID) %>% 
  summarize(control_mean = mean(Avg.Moisture, na.rm = T))

# Join datasets 

moisture.24h <- left_join(avg.24h, baseline.24h)

# Calculate ratio

moisture.24h$Ratio <- moisture.24h$Avg.Moisture/moisture.24h$control_mean

#================================================
# Plot ratios per sampling day

moisture.24h$Day.Title <- paste("Day",moisture.24h$Day)
moisture.24h$Day.Title <- as.factor(moisture.24h$Day.Title)
moisture.24h$Day.Title <- factor(moisture.24h$Day.Title, 
                                 levels = c("Day 0", "Day 8", "Day 16", "Day 32", "Day 64"))

# Only drought treatment

drought.24h <- moisture.24h %>% 
  filter(treatment == "D")

ylab <- expression(paste(bold("Ratio"~D["i"]~":"~C["mean"])))

g5 <- ggplot(drought.24h, aes(x = Site, y = Ratio)) +
  geom_boxplot(aes(color = Region, fill = Region),
               alpha = 0.5) +
  facet_wrap(~ Day.Title, nrow = 5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(0.4, "cm")) +
  labs(x = "Date", y = ylab, title = "Soil moisture ratio 24 h period, starting at 12:00 noon") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  scale_y_continuous(limits = c(0,3.5))

g5

# ggsave(paste0(mydir,"Soil moisture ratio 24 h V1.png"),
#        g5, width = 18, height = 24, units = "cm")

####

g6 <- ggplot(drought.24h, aes(x = as.factor(Day), y = Ratio)) +
  geom_boxplot(aes(color = Region, fill = Region),
               alpha = 0.5) +
  facet_grid(Region ~ Farm) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        strip.text = element_text(size = 12),
        legend.position = "none",
        panel.spacing = unit(0.4, "cm")) +
  labs(x = "Date", y = ylab, title = "Soil moisture ratio 24 h period, starting at 12:00 noon") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  scale_y_continuous(limits = c(0,3.5))

g6

# ggsave(paste0(mydir,"Soil moisture ratio 24 h V2.png"),
#        g6, width = 24, height = 18, units = "cm")

#================================================
# Plot ratios per sampling day regardless of site (see Willson's script)

ylab2 <- expression(paste(bold("Ratio"~"Plot"["i"]~":"~"Plot"["control mean"])))

g7 <- ggplot(data = moisture.24h, 
             aes(x = Day, y = Ratio, 
                 color = treatment, fill = treatment)) + 
  geom_point() + 
  geom_smooth(method = "loess", alpha = 0.2) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0))) +
  labs(y = ylab2)

g7

# ggsave(paste0(mydir,"Soil moisture ratio 24 h V3.png"), 
#        g3, width = 8, height = 5)

#================================================
# Plot ratios per sampling day for each region and farm

g8 <- ggplot(data = moisture.24h, 
             aes(x = Day, y = Ratio, 
                 color = treatment, fill = treatment)) + 
  geom_point() + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        strip.text = element_text(size = 12),
        panel.spacing = unit(0.4, "cm")) +
  geom_smooth(method = "loess", alpha = 0.2) +
  facet_grid(Region ~ Farm) +
  labs(y = ylab2) +
  scale_y_continuous(limits = c(0,4))

g8

# ggsave(paste0(mydir,"Soil moisture ratio 24 h V4.png"), 
#        g8, width = 24, height = 18, units = "cm")

#================================================

##---------------------##
##                     ##
##     48 h period     ##
##                     ##
##---------------------##

# Extract soil moisture from 12:00 noon to 11:00 am two days later
# on all sampling days for each site

#Sys.setenv(tz='UTC')  ## change timezone to UTC

sites <- unique(df$Site)

for (i in 1:length(sites)) {
  
  mysite <- df %>% 
    filter(Site == sites[i])
  
  ## Get all dates of sampling days exactly at 12:00 noon
  
  day0 <- as.POSIXct(mdy(events$EndDrought[events$Site == sites[i]])) + hours(12)
  day8 <- as.POSIXct(mdy(events$Day8[events$Site == sites[i]])) + hours(12)
  day16 <- as.POSIXct(mdy(events$Day16[events$Site == sites[i]])) + hours(12)
  day32 <- as.POSIXct(mdy(events$Day32[events$Site == sites[i]])) + hours(12)
  day64 <- as.POSIXct(mdy(events$Day64[events$Site == sites[i]])) + hours(12)
  
  ## Create intervals of 48 hours for each sampling day
  # This also works with interval(day0, day0 + hours(23)) and using %within%, see lubridate cheatsheet 
  
  day0.48h <- day0 + hours(0:47)
  day8.48h <- day8 + hours(0:47)
  day16.48h <- day16 + hours(0:47)
  day32.48h <- day32 + hours(0:47)
  day64.48h <- day64 + hours(0:47)
  
  dates <- c(day0.48h,day8.48h,day16.48h,day32.48h,day64.48h)
  
  moisture.site <- mysite %>% 
    filter(date.time %in% dates) %>% 
    mutate(Day = case_when(date.time %in% day0.48h ~ 0,
                           date.time %in% day8.48h ~ 8,
                           date.time %in% day16.48h ~ 16,
                           date.time %in% day32.48h ~ 32,
                           date.time %in% day64.48h ~ 64))
  
  if (i == 1) {
    df.48h <- moisture.site
  }
  
  if (i > 1) {
    df.48h <- rbind(df.48h,moisture.site)
  }
  
  rm(day0,day8,day16,day32,day64,
     day0.48h,day8.48h,day16.48h,day32.48h,day64.48h,
     dates,mysite,moisture.site)
}

## Save data

write.table(df.48h, paste0(mydir.data,"Soil moisture on sampling days - 48 h.csv"), sep = ",", row.names = F)
save(df.48h, file = paste0(mydir.data,"Soil_moisture_on_sampling_days_48_h.RData"))

#================================================
# Calculate the average soil moisture over the 48 h period

avg.48h <- df.48h %>% 
  group_by(Day, Site, Region, Farm, Plot, site_ID, treatment) %>% 
  summarize(Avg.Moisture = mean(Moisture, na.rm = T))

## Calculate ratio for each Drought replicate as Ratio = Di / Cmean for every sampling day and every site

# Calculate Cmean (baseline) for each site

baseline.48h <- avg.48h %>% 
  filter(treatment == "C") %>% 
  group_by(Day, Site, Region, Farm, site_ID) %>% 
  summarize(control_mean = mean(Avg.Moisture, na.rm = T))

# Join datasets 

moisture.48h <- left_join(avg.48h, baseline.48h)

# Calculate ratio

moisture.48h$Ratio <- moisture.48h$Avg.Moisture/moisture.48h$control_mean

#================================================
# Plot ratios per sampling day

moisture.48h$Day.Title <- paste("Day",moisture.48h$Day)
moisture.48h$Day.Title <- as.factor(moisture.48h$Day.Title)
moisture.48h$Day.Title <- factor(moisture.48h$Day.Title, 
                                 levels = c("Day 0", "Day 8", "Day 16", "Day 32", "Day 64"))

# Only drought treatment

drought.48h <- moisture.48h %>% 
  filter(treatment == "D")

ylab <- expression(paste(bold("Ratio"~D["i"]~":"~C["mean"])))

g9 <- ggplot(drought.48h, aes(x = Site, y = Ratio)) +
  geom_boxplot(aes(color = Region, fill = Region),
               alpha = 0.5) +
  facet_wrap(~ Day.Title, nrow = 5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(0.4, "cm")) +
  labs(x = "Date", y = ylab, title = "Soil moisture ratio 48 h period, starting at 12:00 noon") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  scale_y_continuous(limits = c(0,3.5))

g9

# ggsave(paste0(mydir,"Soil moisture ratio 48 h V1.png"),
#        g9, width = 18, height = 24, units = "cm")

####

g10 <- ggplot(drought.48h, aes(x = as.factor(Day), y = Ratio)) +
  geom_boxplot(aes(color = Region, fill = Region),
               alpha = 0.5) +
  facet_grid(Region ~ Farm) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        strip.text = element_text(size = 12),
        legend.position = "none",
        panel.spacing = unit(0.4, "cm")) +
  labs(x = "Date", y = ylab, title = "Soil moisture ratio 48 h period, starting at 12:00 noon") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  scale_y_continuous(limits = c(0,3.5))

g10

# ggsave(paste0(mydir,"Soil moisture ratio 48 h V2.png"),
#        g10, width = 24, height = 18, units = "cm")

#================================================
# Plot ratios per sampling day regardless of site (see Willson's script)

ylab2 <- expression(paste(bold("Ratio"~"Plot"["i"]~":"~"Plot"["control mean"])))

g11 <- ggplot(data = moisture.48h, 
             aes(x = Day, y = Ratio, 
                 color = treatment, fill = treatment)) + 
  geom_point() + 
  geom_smooth(method = "loess", alpha = 0.2) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0))) +
  labs(y = ylab2)

g11

# ggsave(paste0(mydir,"Soil moisture ratio 48 h V3.png"), 
#        g11, width = 8, height = 5)

#================================================
# Plot ratios per sampling day for each region and farm

g12 <- ggplot(data = moisture.48h, 
             aes(x = Day, y = Ratio, 
                 color = treatment, fill = treatment)) + 
  geom_point() + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        strip.text = element_text(size = 12),
        panel.spacing = unit(0.4, "cm")) +
  geom_smooth(method = "loess", alpha = 0.2) +
  facet_grid(Region ~ Farm) +
  labs(y = ylab2) +
  scale_y_continuous(limits = c(0,4))

g12

# ggsave(paste0(mydir,"Soil moisture ratio 48 h V4.png"), 
#        g12, width = 24, height = 18, units = "cm")
