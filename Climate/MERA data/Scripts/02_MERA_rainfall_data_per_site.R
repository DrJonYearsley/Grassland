# Author: Maja Ilic
# Last update: 20 Feb 2020

# This script is for plotting MERA rainfall data (Jun - Aug 2017) for each site from the Grassland Resilience 
# Year 1 All Ireland Drought Experiment.

#===========================================
# Clear objects from the workspace

rm(list = ls())

#===========================================
# Set working directory

# Maja's desktop PC

setwd("C:/Users/3054311/Documents/My Documents/Grassland project/04_MERA_rainfall_and_temp")

data.dir <- paste0(getwd(),"/Data/")
figures.dir <- paste0(getwd(),"/Figures/")
script.dir <- paste0(getwd(),"/Script/")

#===========================================
# Packages

library(ggplot2)
library(maps)
library(mapdata)
library(RColorBrewer)
library(scales)
library(dplyr)
library(lubridate)

#===========================================
# Import extracted MERA rainfall data for all sites  

load(file = paste0(data.dir,"MERA_rain_all_sites_Apr_Oct_2017.RData"))
str(rain.2017)

#===========================================
# Plot rainfall data for each site
# Data used: MERA rainfall data, 0-33 hours step (according to Jon)
# All sites using facet_grid

Sys.setenv(tz='UTC')  ## change timezone to UTC

min.date <- as.POSIXct(min(rain.2017$dataDate))
max.date <- as.POSIXct(max(rain.2017$dataDate))

main.title <- c("Percipitation per site in year 2017 - MERA, step 0-33 hrs")
ylab <- expression(paste(bold("Precipitation (kg"~m^"-2"*")")))

g <- ggplot() +
  geom_line(data = rain.2017, 
            aes(x = as.POSIXct(dataDate), y = Rainfall),
            color = "black") +
  facet_grid(Region ~ Farm) +
  scale_x_datetime(limits = c(min.date,max.date),
                   date_breaks = "1 months",
                   date_labels = "%b") +
  scale_y_continuous(limits = c(-5,65), breaks = c(0,20,40,60), exp = c(0,0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 14, face = "bold")) +
  labs(x = "Date", y = ylab)

ggsave(paste0(figures.dir,"Precipitation per site in year 2017 - MERA step 0-33 hrs.png"),
       g,
       width = 24, height = 16, units = "cm")

#===========================================
# All sites using for-loop

# Events

load(file = "C:/Users/3054311/Documents/My Documents/Grassland project/01_Probes_soil_moisture_and_temp/Data/All events.RData")

sites <- unique(rain.2017$Site)

for (site in sites){
  
  mysite.MERA <- rain.2017 %>% 
    filter(Site == site) %>% 
    droplevels()
  
  start.drought <- as.POSIXct(mdy(events$StartDrought[events$Site == site]), tz = "GMT") + hours(12)
  end.drought <- as.POSIXct(mdy(events$EndDrought[events$Site == site]), tz = "GMT") + hours(12)
  day8 <- as.POSIXct(mdy(events$Day8[events$Site == site]), tz = "GMT") + hours(12)
  day16 <- as.POSIXct(mdy(events$Day16[events$Site == site]), tz = "GMT") + hours(12)
  day32 <- as.POSIXct(mdy(events$Day32[events$Site == site]), tz = "GMT") + hours(12)
  day64 <- as.POSIXct(mdy(events$Day64[events$Site == site]), tz = "GMT") + hours(12)
  
  ## 
  
  title <- paste(mysite.MERA$Region[1],mysite.MERA$Farm[1])
  figure.title <- paste("Rainfall for",title,"in year 2017 - MERA step 0-33 hrs.png")
  
  ##
  
  plot.site <- ggplot() + 
    geom_rect(data = mysite.MERA,
              xmin = start.drought,
              xmax = end.drought,
              ymin = -Inf,
              ymax = +Inf,
              fill = "gray80") +
    geom_line(data = mysite.MERA, 
              aes(x = as.POSIXct(dataDate), y = Rainfall)) +
    geom_vline(xintercept = day8, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = day16, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = day32, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = day64, linetype = "dashed", color = "gray50") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
          axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
    labs(x = "Date", y = ylab, title = title) + 
    scale_x_datetime(limits = c(min.date,max.date),
                     date_breaks = "1 months",
                     date_labels = "%b") +
    scale_y_continuous(limits = c(-5,65), breaks = c(0,20,40,60), exp = c(0,0))
  
  ##
  
  ggsave(paste0(figures.dir,figure.title),
         plot.site,
         width = 15, height = 12, units = "cm")
  
}  
