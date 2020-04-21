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
# 4. ggplot object g2 (panels: Region and Farm)

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
library(patchwork)
library(lubridate)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(RColorBrewer)

#================================================
# Load data ----

# Soil moisture

load(paste0(data.dir,"Cleaned_Soil_Moisture_Data_from_Loggers.RData"))

# Events

load(file = paste0(data.dir,"All events.RData"))

#================================================
# Set directory for plots ----

# dir.create(paste0(figures.dir,"Soil moisture moving window plots"))
mydir <- paste0(figures.dir,"Soil moisture moving window plots/")

#================================================
# For each site, plot the moving window on the top of the raw data

sites <- unique(df.moisture$Site)
exp.plots <- unique(df.moisture$Plot)

for (i in 1:length(sites)){
  
  # Change timezone to UTC
  
  Sys.setenv(tz='UTC')  
  
  # Get sampling days 0, 8, 16, 32 and 64 for the specific site
  
  day0 <- as.POSIXct(mdy(events$EndDrought[events$Site == sites[i]])) + hours(12)
  day8 <- as.POSIXct(mdy(events$Day8[events$Site == sites[i]])) + hours(12)
  day16 <- as.POSIXct(mdy(events$Day16[events$Site == sites[i]])) + hours(12)
  day32 <- as.POSIXct(mdy(events$Day32[events$Site == sites[i]])) + hours(12)
  day64 <- as.POSIXct(mdy(events$Day64[events$Site == sites[i]])) + hours(12)
  
  exp.days <- data.frame("Day" = c(0,8,16,32,64), 
                         "date.time" = c(day0,day8,day16,day32,day64))
  
  # Define the start of the plot as sampling day 0 (12:00) minus 8 days
  
  start.plot <- day0 - days(8) 
  
  # Extract soil moisture data for the specific site within the time window day 0 - 8 days and day 64 post-drought
  
  mysite <- df.moisture %>% 
    filter(Site == sites[i] & date.time >= start.plot & date.time <= day64) %>% 
    droplevels()
  
  for (j in 1:length(exp.plots)) {
    
    # Previous data exploration showed that at certain window widths, less data is included due to missing data
    # These window widths are:
    # - 1 h           Limerick 2, Day 32
    # - 16 h          Limerick 2, Day 0
    # - 39 h          Border 1, Day 64
    # - 80 h          Dublin 5, Day 64
    # Total window width: -192 h (= 8 days)
    
    # For each experimental plot (C1 - D3), find the lowest and highest soil moisture within the given window range
    # for each sampling day and set it as lower and upper border of the moving window rectangle
    
    # First, extract soil moisture data for the specific experimental plot
    
    myexp.plot <- mysite %>% 
      filter(Plot == exp.plots[j]) %>% 
      droplevels()
    
    # Plot data
    
    title.plot <- paste("Plot",exp.plots[j])
    
    g <- ggplot(myexp.plot, aes(x = date.time, y = Moisture)) +
      theme_minimal() +
      theme(axis.title = element_blank(),
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
      labs(title = title.plot) +
      scale_y_continuous(limits = c(0,105), breaks = c(0,20,40,60,80,100), exp = c(0,0))
    
    # Set colors for moving window rectangles
    
    mycol <- colorRampPalette(brewer.pal(8, "Set1"))(8)
    
    # Plot moving window rectangles for each sampling day and each time step as provided above
    
    for (d in c(0,8,16,32,64)) {
      
      m <- 3
      
      for (h in c(192,81,40,17,2,0)) {
        
        xmin <- exp.days$date.time[exp.days$Day == d] - hours(h)
        xmax <- exp.days$date.time[exp.days$Day == d]
        
        win.range <- myexp.plot %>% filter(date.time >= xmin & date.time <= xmax)
        
        ymin <- min(win.range$Moisture, na.rm = T)
        ymax <- max(win.range$Moisture, na.rm = T)
        
        if(is.na(ymin)){
          ymin <- -Inf
        }
        
        if(is.na(ymax)){
          ymax <- Inf
        }
        
        g <- g + annotate("rect",
                          xmin = xmin,
                          xmax = xmax,
                          ymin = ymin,
                          ymax = ymax,
                          color = mycol[m],
                          fill = mycol[m],
                          alpha = 0.5)
        
        # Increase m
        
        m <- m + 1
        
      }
    }
    
    g <- g + geom_line() +
      geom_vline(xintercept = day0, color = "grey50", linetype = "dashed") +
      geom_vline(xintercept = day8, color = "grey50", linetype = "dashed") +
      geom_vline(xintercept = day16, color = "grey50", linetype = "dashed") +
      geom_vline(xintercept = day32, color = "grey50", linetype = "dashed") +
      geom_vline(xintercept = day64, color = "grey50", linetype = "dashed")
    
    assign(paste0("plot.",exp.plots[j]), g)
    rm(g)
  }
  
  # Final plot
  
  plot.final <- plot_grid(plot.C1,plot.C2,plot.C3,
                          NULL,NULL,NULL,
                          plot.D1,plot.D2,plot.D3,
                          ncol = 3,
                          labels = NULL,
                          rel_heights = c(10,1,10))
  
  y.grob <- textGrob("Soil Moisture (%)", 
                     gp = gpar(fontsize = 14, fontface = "bold"), rot = 90)
  
  x.grob <- textGrob("Date", 
                     gp = gpar(fontsize = 14, fontface = "bold"))
  
  title.site <- paste(mysite$Region[1],mysite$Farm[1])
  
  title.grob <- textGrob(title.site, 
                         gp = gpar(fontsize = 18, fontface = "bold"))
  
  figure.title <- paste0(mydir,title.site,".png")

  ggsave(figure.title,
         grid.arrange(arrangeGrob(plot.final, left = y.grob, bottom = x.grob, top = title.grob)),
         width = 14, height = 10)
}
