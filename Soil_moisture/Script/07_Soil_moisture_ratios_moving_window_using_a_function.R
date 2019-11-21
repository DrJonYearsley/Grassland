# Author: Maja Ilic
# Last update: November 2019

# This is a script to calculate soil moisture ratios Drought : Control 
# from the Grassland Resilience Year 1 All Ireland Drought Experiment. 
# In one of the previous scripts (05), I calculated running ratio
# (ratio of drought to mean control at each time point, resolution: 1 hours). 
# In this script, I am extracting the soil mositure data from drought and control treatments 
# on/around sampling days 0, 8, 16, 32 and 64. 
# In the next step, I calculate again the ration of drought to mean control.
# To get the ratio for each of the sampling days, several options were suggested (see script 06).

# In this script, I define a function that takes five inputs:
# 1. Data (e.g. soil moisture)
# 2. Dates (e.g. dates of the sampling days in each site)
# 3. Time of the day (e.g. 12:00 noon) set as the starting point of the moving window (in 00:00 - 23:00)
# 4. Duration of the period (width of the moving window), e.g. 24 h
# 5. Limits (range) for the y-axis

# Outputs of the function are the following:
# 1. Extracted soil moisture data for all days (or intervals/windows) set 
# 2. Averaged soil moisture over the specified period/window
# 3. ggplot object g1 (panels: sampling days)
# 4. ggplot object g2 (Pnales: Region and Farm)

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
# Set directory for pdata and lots

mydir.data <- paste0(data.dir,"Soil moisture ratios DtoC per site and sampling day/")
mydir <- paste0(figures.dir,"Soil moisture ratios DtoC per site and sampling day/")

#================================================
# Use function moisture.mov.win

source(paste0(script.dir,"Function moisture.mov.win.R"))

## Use function

# Start: 12:00 noon

window.0h <- moisture.mov.win(data = df,
                              dates = events,
                              daytime.hr = 12,
                              duration.hr = 0,
                              ylim = c(0,3.5))

window.6h <- moisture.mov.win(data = df,
                              dates = events,
                              daytime.hr = 12,
                              duration.hr = 6,
                              ylim = c(0,3.5))

window.12h <- moisture.mov.win(data = df,
                              dates = events,
                              daytime.hr = 12,
                              duration.hr = 12,
                              ylim = c(0,3.5))

window.18h <- moisture.mov.win(data = df,
                              dates = events,
                              daytime.hr = 12,
                              duration.hr = 18,
                              ylim = c(0,3.5))

window.24h <- moisture.mov.win(data = df,
                              dates = events,
                              daytime.hr = 12,
                              duration.hr = 24,
                              ylim = c(0,3.5))

window.36h <- moisture.mov.win(data = df,
                              dates = events,
                              daytime.hr = 12,
                              duration.hr = 36,
                              ylim = c(0,3.5))

window.48h <- moisture.mov.win(data = df,
                              dates = events,
                              daytime.hr = 12,
                              duration.hr = 48,
                              ylim = c(0,3.5))

# Start: 16:00

window.0h.16 <- moisture.mov.win(data = df,
                              dates = events,
                              daytime.hr = 16,
                              duration.hr = 0,
                              ylim = c(0,3.5))

window.6h.16 <- moisture.mov.win(data = df,
                              dates = events,
                              daytime.hr = 16,
                              duration.hr = 6,
                              ylim = c(0,3.5))

window.12h.16 <- moisture.mov.win(data = df,
                               dates = events,
                               daytime.hr = 16,
                               duration.hr = 12,
                               ylim = c(0,3.5))

window.18h.16 <- moisture.mov.win(data = df,
                               dates = events,
                               daytime.hr = 16,
                               duration.hr = 18,
                               ylim = c(0,3.5))

window.24h.16 <- moisture.mov.win(data = df,
                               dates = events,
                               daytime.hr = 16,
                               duration.hr = 24,
                               ylim = c(0,3.5))

window.36h.16 <- moisture.mov.win(data = df,
                               dates = events,
                               daytime.hr = 16,
                               duration.hr = 36,
                               ylim = c(0,3.5))

window.48h.16 <- moisture.mov.win(data = df,
                               dates = events,
                               daytime.hr = 16,
                               duration.hr = 48,
                               ylim = c(0,3.5))
