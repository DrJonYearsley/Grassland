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
# Clear objects from the workspace ----

rm(list = ls())

#===========================================
# Set working directory ----

# Maja's desktop PC

setwd("C:/Users/3054311/Documents/My Documents/Grassland project/01_Probes_soil_moisture_and_temp")

data.dir <- paste0(getwd(),"/Data/")
figures.dir <- paste0(getwd(),"/Figures/")
script.dir <- paste0(getwd(),"/Script/")

#===========================================
# Packages ----

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
# Load data ----

# Soil moisture

load(paste0(data.dir,"Cleaned_Soil_Moisture_Data_from_Loggers.RData"))

# Events

load(file = paste0(data.dir,"All events.RData"))

#================================================
# Create new folder for data and plots ----

dir.create(paste0(figures.dir,"Soil moisture effect size boxplot"))
dir.create(paste0(figures.dir,"Soil moisture effect size V1"))
dir.create(paste0(figures.dir,"Soil moisture effect size V2"))

mydir <- paste0(figures.dir)

# Create a list of figure titles to be used within the mov.win.diff() function

mydir.plot1 <- paste0(mydir,"/Soil moisture effect size boxplot/")
mydir.plot2 <- paste0(mydir,"/Soil moisture effect size V1/")
mydir.plot3 <- paste0(mydir,"/Soil moisture effect size V2/")

myplot <- c("plot1", "plot2", "plot3")
plot.dir <- c(mydir.plot1, mydir.plot2, mydir.plot3)

list.figures <- setNames(as.list(plot.dir), myplot)

#================================================
# Load function mov.win.diff ----

source("C:/Users/3054311/Documents/My Documents/Grassland project/Functions/Function mov.win.diff.2.R")

#================================================
# Plot soil moisture ES

daytime.hr <- 12
width <- seq(-192,0,1)

for (duration.hr in width) {
  
  # Use the mov.win.diff function to extract relevant soil moisture data for the given time inverval
  # and calculate the effect size (for the days 0, 32 and 64)
  
  df_soil <- mov.win.diff.2(data = df.moisture,
                            dates = events,
                            daytime.hr = daytime.hr,
                            duration.hr = duration.hr,
                            doplot = TRUE,
                            list.figures = list.figures)
}
