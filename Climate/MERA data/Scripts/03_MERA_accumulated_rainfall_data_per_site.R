# Author: Maja Ilic
# Last update: 29 Feb 2020

# This script is used to calculate accumulated rainfall for each farm and region during the simulated drought period
# Data used is MERA rainfall data (Jun - Aug 2017) for each site from the Grassland Resilience 
# Year 1 All Ireland Drought Experiment.

#===========================================
# Clear objects from the workspace ----

rm(list = ls())

#===========================================
# Set working directory ----

# Maja's desktop PC

setwd("C:/Users/3054311/Documents/My Documents/Grassland project/04_MERA_rainfall_and_temp")

data.dir <- paste0(getwd(),"/Data/")
figures.dir <- paste0(getwd(),"/Figures/")
script.dir <- paste0(getwd(),"/Script/")

#===========================================
# Packages ----

library(ggplot2)
library(maps)
library(mapdata)
library(RColorBrewer)
library(scales)
library(dplyr)
library(lubridate)
library(purrr)

#===========================================
# Import extracted MERA rainfall data for all sites ---- 

load(file = paste0(data.dir,"MERA_rain_all_sites_Apr_Oct_2017.RData"))
head(rain.2017)
str(rain.2017)

#===========================================
# Import events data ----

load(file = "C:/Users/3054311/Documents/My Documents/Grassland project/01_Probes_soil_moisture_and_temp/Data/All events.RData")

#===========================================
# Accumulated rainfall during the simulated drought ----
# For each site, calculate how much rain was "kept out" during the simulated drought 
# Additionally calculate stepwise accumulation of the rainfall
# All sites using for-loop

sites <- unique(rain.2017$Site)

Sys.setenv(tz='UTC')  ## change timezone to UTC

for (site in sites) {
  
  # Start and end date of the simulated drought
  start.drought <- as.POSIXct(mdy(events$StartDrought[events$Site == site]))
  end.drought <- as.POSIXct(mdy(events$EndDrought[events$Site == site]))
  
  # Subset rainfall data for the current site
  mysite.rain.drought <- rain.2017 %>% 
    filter(Site == site & dataDate >= start.drought & dataDate <= end.drought) %>% 
    droplevels()
  
  # Calculate stepwise accumulated rainfall (growing daily precipitation, GDP)
  mysite.rain.drought$GDP <- accumulate(mysite.rain.drought$Rainfall,`+`)
  
  if (site == sites[1]) {
    sites.rain.drought <- mysite.rain.drought
  }
  
  if (site != sites[1]) {
    sites.rain.drought <- rbind(sites.rain.drought,mysite.rain.drought)
  }
  
}

# Accumulated rain per site
rain.acc <- sites.rain.drought %>% 
  group_by(Site,Region,Farm) %>% 
  summarize(Acc.Rain.D = sum(Rainfall))

# Average accumulated rain per region
rain.acc <- rain.acc %>% 
  group_by(Region) %>% 
  mutate(Avg.Acc.Rain.D = mean(Acc.Rain.D))

#===========================================
# Export data ----

write.table(sites.rain.drought,paste0(data.dir,"MERA Percipitation during drought.csv"), sep = ",", row.names = F)
save(sites.rain.drought,file = paste0(data.dir,"MERA_Percipitation_during_drought.RData"))
write.table(rain.acc,paste0(data.dir,"MERA Accumulated percipitation during drought.csv"), sep = ",", row.names = F)
save(rain.acc,file = paste0(data.dir,"MERA_Accumulated_percipitation_during_drought.RData"))

#===========================================
# Plot accumulated rainfall ----

ylab1 <- expression(paste(bold("Precipitation (kg"~m^"-2"*")")))
title1 <- "Accumulated precipitation\nduring the simulated drought period\n"

rain.acc$Region <- factor(rain.acc$Region, levels = c("Border","Cork","Dublin","Limerick"))

g1 <- ggplot(rain.acc, aes(x = Site, y = Acc.Rain.D)) +
  geom_point(aes(fill = Region, color = Region), size = 5, alpha = 0.5, shape = 25) +
  geom_hline(aes(yintercept = Avg.Acc.Rain.D, color = Region), size = 1, linetype = "dashed") +
  facet_wrap(~ Region, ncol = 4, scales = "free_x") +
  theme_minimal() +
  theme(strip.background = element_rect(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.position = "none") +
  scale_y_continuous(limits = c(140,260), breaks = c(140, 180, 220, 260)) +
  labs(x = "Region", y = ylab1, title = title1)

ggsave(paste0(figures.dir, "Accumulated precipitation during the simulated drought period per site - MERA2.png"),
       g1,
       width = 18, height = 13, units = "cm")

#===========================================
# Plot growing daily precipitation

ylab2 <- expression(paste(bold("Growing Daily Precipitation (kg"~m^"-2"*")")))
title2 <- "Growing daily precipitation\nduring the simulated drought period\n"

sites.rain.drought$Region <- factor(sites.rain.drought$Region, levels = c("Border","Cork","Dublin","Limerick"))

g2 <- ggplot(sites.rain.drought, aes(x = dataDate)) +
  geom_line(aes(y = Rainfall), size = 0.8) +
  geom_line(aes(y = GDP), color = "blue", size = 1) +
  facet_grid(Region ~ Farm, scales = "free_x") +
  theme_test() +
  theme(strip.background = element_rect(),
        axis.text.x = element_text(size = 12, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.position = "none") +
  labs(x = "Date", y = ylab2, title = title2)

ggsave(paste0(figures.dir, "Growing Daily Precipitation during the simulated drought period per site - MERA.png"),
       g2,
       width = 24, height = 20, units = "cm")
