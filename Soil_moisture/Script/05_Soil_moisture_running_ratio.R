# Author: Maja Ilic
# Last update: November 2019

# This is a script to calculate soil moisture ratios Drought : Control from the Grassland Resilience 
# Year 1 All Ireland Drought Experiment.

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
# Create new folder for plots

dir.create(paste0(figures.dir,"Soil moisture ratios DtoC per site"))
mydir <- paste0(figures.dir,"Soil moisture ratios DtoC per site/")

#================================================

##-----------------------##
##                       ##
##     Running ratio     ##
##                       ##
##-----------------------##

## Ratio calculated at all time points that are in the original soil probe data

## Idea: use the mean of Control treatments as "baseline"

## Calculate ratio for each Drought replicate as Ratio = Di / Cmean

# Calculate the baseline (mean of control treatments) for each time point and site

baseline <- df %>% 
  filter(treatment == "C") %>% 
  group_by(date.time, Site, Region, Farm, site_ID) %>% 
  summarize(control_mean = mean(Moisture, na.rm = T))

# Extract drought treatment from the raw data

drought <- df %>% 
  filter(treatment == "D")

# Join datasets 

moisture <- left_join(drought, baseline)

# Calculate running ratio

moisture$Ratio <- moisture$Moisture/moisture$control_mean

# Plot ratios per site

sites <- unique(moisture$Site)

ylab <- expression(paste(bold("Ratio"~D["i"]~":"~C["mean"])))

min.date <- min(moisture$date.time)
max.date <- max(moisture$date.time)

# Inspecting the ratios showed that there are some outliers
# I extracted these and exported them in an extra file
# In total, 1750 outliers were found (threshold = 2)
# of which 1632 were from the site L4 (Limerick 4)

outliers <- moisture %>% 
  filter(Ratio > 2)

write.table(outliers, paste0(data.dir,"Soil moisture ratio outliers.csv"), 
            sep = ",", row.names = F)

# for-loop

for (site in sites){
  
  mysite <- moisture %>% 
    filter(Site == site) %>% 
    droplevels()
  
  title <- paste(mysite$Region[1],mysite$Farm[1])

  start.drought <- as.POSIXct(mdy(events$StartDrought[events$Site == site]), tz = "GMT") + 11*60*60
  end.drought <- as.POSIXct(mdy(events$EndDrought[events$Site == site]), tz = "GMT") + 11*60*60
  day8 <- as.POSIXct(mdy(events$Day8[events$Site == site]), tz = "GMT") + 11*60*60
  day16 <- as.POSIXct(mdy(events$Day16[events$Site == site]), tz = "GMT") + 11*60*60
  day32 <- as.POSIXct(mdy(events$Day32[events$Site == site]), tz = "GMT") + 11*60*60
  day64 <- as.POSIXct(mdy(events$Day64[events$Site == site]), tz = "GMT") + 11*60*60
  
  single.panel <- ggplot(mysite, aes(x = date.time, y = Ratio)) +
    geom_rect(xmin = start.drought,
              xmax = end.drought,
              ymin = -Inf,
              ymax = +Inf,
              fill = "gray80") +
    geom_line(aes(linetype = Plot)) + 
    geom_vline(xintercept = day8, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = day16, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = day32, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = day64, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = 1, color = "red") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 12),
          axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
          axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          legend.title = element_text(size = 14, face = "bold")) +
    labs(x = "Date", y = ylab, title = title) +
    scale_x_datetime(limits = c(min.date,max.date),
                     breaks = c(start.drought, 
                                end.drought),
                     date_labels = "%b %d") +
    scale_y_continuous(limits = c(0,2), breaks = c(0,0.5,1,1.5,2))
  
  figure.title <- paste0(mydir,site," - Soil moisture running ratio.png")
  
  ggsave(figure.title,
         single.panel,
         width = 17, height = 12, units = "cm")
  
  # Remove axis titles and legend and assign each plot to a named object
  
  single.panel <- single.panel +
    theme(axis.title = element_blank(),
          legend.position = "none")
  
  assign(site,single.panel)
  
  rm(single.panel)

}
  