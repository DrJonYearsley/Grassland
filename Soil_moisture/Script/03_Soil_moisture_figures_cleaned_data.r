# Authors: Mark Emmerson and Maja Ilic
# Last update: November 2019

# This is a script to explore and plot soil moisture data from the Grassland Resilience 
# Year 1 All Ireland Drought Experiment.

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
library(lubridate)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)
library(cowplot)

#===========================================
# Load data and add column treatment

load(paste0(data.dir,"Cleaned_Soil_Moisture_Data_from_Loggers.RData"))
df <- df.moisture   ## shorter
rm(df.moisture)
names(df)

dir.create(paste0(figures.dir,"Soil moisture per site"))
mydir <- paste0(figures.dir,"Soil moisture per site/")

#================================================
# Theme for ggplots

min.date <- min(df$date.time)
max.date <- max(df$date.time)

mytheme <- theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        panel.spacing = unit(1, "lines"),
        legend.position = "none")

#================================================
# Add events (start drought, end drought, days 8. 16. 32 and 64)

load(file = paste0(data.dir,"All events.RData"))

##------------------##
##                  ##
##     for-loop     ##
##                  ##
##------------------##

regions <- unique(df$Region)
sites <- unique(df$Site)

for (site in sites){
  
  mysite <- df %>% 
    filter(Site == site) %>% 
    droplevels()
  
  title <- paste(mysite$Region[1],mysite$Farm[1])

  ## Get all important dates (only day, month and year) and add 11 hours ( = 12:00:00)
  
  start.drought <- as.POSIXct(mdy(events$StartDrought[events$Site == site]), tz = "GMT") + 11*60*60
  end.drought <- as.POSIXct(mdy(events$EndDrought[events$Site == site]), tz = "GMT") + 11*60*60
  day8 <- as.POSIXct(mdy(events$Day8[events$Site == site]), tz = "GMT") + 11*60*60
  day16 <- as.POSIXct(mdy(events$Day16[events$Site == site]), tz = "GMT") + 11*60*60
  day32 <- as.POSIXct(mdy(events$Day32[events$Site == site]), tz = "GMT") + 11*60*60
  day64 <- as.POSIXct(mdy(events$Day64[events$Site == site]), tz = "GMT") + 11*60*60
  
  g <- ggplot(mysite,
              aes(x = date.time, y = Moisture)) +
    geom_rect(aes(xmin = start.drought,
                  xmax = end.drought,
                  ymin = -Inf,
                  ymax = +Inf),
              fill = "gray80") +
    geom_line(aes(color = treatment)) +
    facet_wrap(~ Plot) +
    labs(x = "Date", y = "Soil Moisture (%)", title = title) +
    mytheme +
    scale_color_manual(values = c("red","blue")) +
    scale_y_continuous(limits = c(-2,102), breaks = c(0,25,50,75,100)) +
    scale_x_datetime(limits = c(min.date,max.date), 
                     breaks = c(start.drought,
                                end.drought,
                                day64),
                     date_labels = "%b %d") +
    geom_vline(aes(xintercept = day8),
               linetype = "dashed", color = "gray50") +
    geom_vline(aes(xintercept = day16),
               linetype = "dashed", color = "gray50") +
    geom_vline(aes(xintercept = day32),
               linetype = "dashed", color = "gray50") +
    geom_vline(aes(xintercept = day64),
               linetype = "dashed", color = "gray50")
  
  ggsave(paste0(mydir,site," - Soil moisture.png"), g,
         width = 8, height = 6)
  
}
