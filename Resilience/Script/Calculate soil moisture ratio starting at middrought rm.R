###################################################################################
## This script calculates the soil moisture ratio for each sampling event with   ##
## a given window width. The idea is to start at the middrought (grass trimming) ##
## and calculate the soil moisture ratio based on the average soil moisture      ##
## between two grass trimming (sampling) events. This would then correspond to   ##
## the same time window as the plant biomass (repeated measures). This means,    ##
## soil moisture ratio for day 0 is calculated based on average soil moisture    ##
## between middrought (12:00 noon) and day 0 (12:00 noon), soil moisture ratio   ##
## for day 8 is calculated based on the average soil moisture between day 0      ##
## (12:00 noon) and day 8 (12:00) and so on.                                     ##
##                                                                               ##
## Author of the script:                                                         ##
## Maja Ilic M.Ilic@qub.ac.uk                                                    ##
## first modified: 24 Mar 2020                                                   ##
## last modified: 21 Apr 2020                                                    ##
###################################################################################

#================================================
# Clear objects from the workspace ----

rm(list = ls())

#================================================
# Set working directory ----

# Maja's desktop PC

setwd("C:/Users/3054311/Documents/My Documents/Grassland project/05_Resilience")

data.dir <- paste0(getwd(),"/Data/")

#================================================
# Packages ----

library(ggplot2)
library(lubridate)
library(dplyr)
library(tidyr)
library(stringr)

#================================================
# Get raw data for soil moisture ----

raw.data.dir <- "C:/Users/3054311/Documents/My Documents/Grassland project/01_Probes_soil_moisture_and_temp/Data/"

load(paste0(raw.data.dir,"Cleaned_Soil_Moisture_Data_from_Loggers.RData"))

# Events 

load(file = paste0(raw.data.dir,"All events.RData"))

# Middrought dates

base_mid_dates <- read.csv("C:/Users/3054311/Documents/My Documents/Grassland project/BaselineAndMiddroughtDates.csv", 
                           sep = ",", header = T)

mid_dates <- base_mid_dates %>% filter(Event == "Middrought") %>% droplevels()

#================================================
# Calculate soil moisture ratio ----
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
  
  mid.drought <- as.POSIXct(dmy(mid_dates$Date[mid_dates$Site == sites[i]])) + hours(daytime.hr)
  mid.day0.duration <- as.numeric(day0 - mid.drought)
  
  ## Create intervals for each sampling day of specific duration
  
  day0.interval <- day0 - hours(0:(24*mid.day0.duration-1))
  day8.interval <- day8 - hours(0:(24*8-1))
  day16.interval <- day16 - hours(0:(24*16-1))
  day32.interval <- day32 - hours(0:(24*32-1))
  day64.interval <- day64 - hours(0:(24*64-1))
  
  ## Extract soil moisture for each interval
  
  date.intervals <- c(day0.interval,
                      day8.interval,
                      day16.interval,
                      day32.interval,
                      day64.interval)
  
  moisture.site <- mysite %>% 
    filter(date.time %in% date.intervals) %>% 
    mutate(day = case_when(date.time %in% day0.interval ~ 0,
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
     date.intervals,mid.drought,mid.day0.duration,mysite,moisture.site)
}

# Calculate the average soil moisture over the period

avg.window <- df.window %>% 
  group_by(day, Site, Region, Farm, Plot, site_ID, treatment) %>% 
  summarize(Avg.Moisture = mean(Moisture, na.rm = T))

# Calculate Cmean (baseline) for each site

baseline.window <- avg.window %>% 
  filter(treatment == "C") %>% 
  group_by(day, Site, Region, Farm, site_ID) %>% 
  summarize(control_mean = mean(Avg.Moisture, na.rm = T))

# Join datasets 

moisture.window <- left_join(avg.window, baseline.window)

# Calculate ratio of soil moisture in each plot to the average soil moisture in the control plots in the respective site

moisture.window$plot_meanC_ratio <- moisture.window$Avg.Moisture/moisture.window$control_mean

# Rename

resilience_soil_rm <- moisture.window

# Rename C to Control and D to drought

resilience_soil_rm$treatment <- gsub("C", "Control", 
                                     as.character(resilience_soil_rm$treatment))

resilience_soil_rm$treatment <- gsub("D", "Drought", 
                                     as.character(resilience_soil_rm$treatment))

# Change Region, Farm and treatment to factor

resilience_soil_rm$Region <- as.factor(resilience_soil_rm$Region)
resilience_soil_rm$Farm <- as.factor(resilience_soil_rm$Farm)
resilience_soil_rm$treatment <- as.factor(resilience_soil_rm$treatment)

#================================================
# Load rainfall data (MERA) ----  

load(file = "C:/Users/3054311/Documents/My Documents/Grassland project/04_MERA_rainfall_and_temp/Data/MERA_Accumulated_percipitation_during_drought.RData")
rain.acc$Farm <- as.factor(rain.acc$Farm)

# Combine rainfall with soil moisture data

resilience_soil_rm <- left_join(resilience_soil_rm,
                                rain.acc[,c("Region","Farm","Acc.Rain.D")],
                                by = c("Region","Farm"))

# Rename Acc.Rain.D to rain_out

names(resilience_soil_rm)[names(resilience_soil_rm) == "Acc.Rain.D"] <- "acc_rain"

# Replace accumulated rainfall with zero for all drought plots

resilience_soil_rm$acc_rain[resilience_soil_rm$treatment == "Drought"] <- 0

resilience_soil_rm$Region <- as.factor(resilience_soil_rm$Region)

#================================================
# Save data ----

write.table(resilience_soil_rm, paste0(data.dir,"Resilience_Soil_with_Rainfall_rm.csv"), sep = ",", row.names = F)
save(resilience_soil_rm, file = paste0(data.dir,"Resilience_Soil_with_Rainfall_rm.RData"))

#================================================
# Plot ratios per region and farm ----

ylab <- expression(paste(bold("Ratio"~"Plot"[bolditalic("i")]~":"~"Control"[bolditalic("mean")])))

g1 <- ggplot(resilience_soil_rm,aes(x = day, y = plot_meanC_ratio)) + 
  geom_point(aes(color = treatment, fill = treatment),
             size = 3, shape = 21, alpha = 0.5) +
  facet_grid(Region ~ Farm) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(),
        panel.spacing = unit(0.6, "lines"),
        axis.text = element_text(size = 10),
        axis.title.y = element_text(size = 14, margin = margin(0,15,0,0)),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        plot.title = element_text(size = 16, face = "bold")) +
  scale_color_manual(values = c("blue","red")) +
  scale_fill_manual(values = c("blue","red")) +
  scale_y_continuous(limits = c(0,3), exp = c(0,0)) +
  scale_x_continuous(limits = c(-1,65), breaks = c(0,8,16,32,64)) +
  labs(x = "Days since the end of the simulated drought",
       y = ylab,
       title = "Soil moisture ratio",
       subtitle = "Starting at mid-drought (rm)") 

g2 <- ggplot(resilience_soil_rm,aes(x = day, y = plot_meanC_ratio)) + 
  geom_smooth(aes(color = treatment),
              method = "lm", se = F) +
  facet_grid(Region ~ Farm) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(),
        panel.spacing = unit(0.6, "lines"),
        axis.text = element_text(size = 10),
        axis.title.y = element_text(size = 14, margin = margin(0,15,0,0)),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        plot.title = element_text(size = 16, face = "bold")) +
  scale_color_manual(values = c("blue","red")) +
  scale_y_continuous(limits = c(0,3), exp = c(0,0)) +
  scale_x_continuous(limits = c(-1,65), breaks = c(0,8,16,32,64)) +
  labs(x = "Days since the end of the simulated drought",
       y = ylab,
       title = "Soil moisture ratio",
       subtitle = "Starting at mid-drought (rm)") 

ggsave("Figures/Soil moisture ratio - Middrought (rm) V1.png", g1, width = 10, height = 8)
ggsave("Figures/Soil moisture ratio - Middrought (rm) V2.png", g2, width = 10, height = 8)
