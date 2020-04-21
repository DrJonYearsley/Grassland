# Authors: Jon Yearsley and Maja Ilic
# Last update: 20 Feb 2020

# This script imports MERA data (previously extracted from the GRIB1 files using grib_ls and grib_get_data functions). 
# Next, rainfall data is extracted for each site from the Grassland Resilience 
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
library(stringr)
library(lubridate)
library(sp)
library(dplyr)
library(maps)
library(mapdata)

#===========================================
# Define all files to be imported ----

all_files <- dir(data.dir)

myfiles <- all_files[grep(pattern = "TotalPrecip", x = all_files)]

#===========================================
# Use a for-loop to import and read all data at once ----

setwd(paste0(getwd(),"/Data"))  
Sys.setenv(tz='UTC')  ## change timezone to UTC
ac <- as.character
an <- as.numeric

# Coordinates of each site ----

sites <- read.csv("C:/Users/3054311/Documents/My Documents/Grassland project/GPS data sites Grassland.csv", sep = ";", header = T)

# Open a graphical window for the png-file ----
# (used later for the map of the Island or Ireland with all sites and MERA coordinates used to extract rainfall data for 2017)

nrow.png <- ncol.png <- ceiling(sqrt(length(myfiles)))
width.png <- ncol.png * 6
height.png <- nrow.png * 6

png(paste0(figures.dir,"Island of Ireland map with all sites and MERA coordinates for rainfall 2017.png"),
    width = width.png, height = height.png, units = "cm", res = 600)
par(mfrow = c(nrow.png,ncol.png), mar = c(1,1,1,1))

# Run for-loop to extract data ----

for (m in 1:length(myfiles)) {

  # Import MERA rainfall data for 2017  

  rain <- read.csv(myfiles[m], sep = "", header = T)
  names(rain) <- c("Latitude", "Longitude", "Rainfall", "dataDate", "dataTime", "validityDate", "validityTime")
  
  exclude <- which(rain$Latitude == "Latitude,")
  
  if (length(exclude) > 0) {
    rain <- rain[-exclude,]
  }
  
  # Get proper date format
  
  rain$dataDate <- ymd(rain$dataDate)
  rain$validityDate <- ymd(rain$validityDate)
  
  # Extract month and year for each dataset (useful later)
  
  month_year <- paste0(lubridate::month(rain$dataDate[1], label = T),"_",year(rain$dataDate[1]))
  month_year2 <- paste(lubridate::month(rain$dataDate[1], label = T),year(rain$dataDate[1]))
  
  # Change Latitude, Longitude, Rainfall, dataTime and validityTime to numerical values
  
  rain$Latitude <- an(ac(rain$Latitude))
  rain$Longitude <- an(ac(rain$Longitude))
  rain$Rainfall <- an(ac(rain$Rainfall))
  rain$dataTime <- an(ac(rain$dataTime))
  rain$validityTime <- an(ac(rain$validityTime))
  
  # Substract 360° from the Longitude
  
  rain$Longitude[which(rain$Longitude > 180)] <- rain$Longitude[which(rain$Longitude > 180)] - 360
  
  # Subset coordinates for Ireland
  
  rain <- rain %>% 
    filter(Longitude >= -11 & Longitude <= -5)
  
  # All unique combinations of latitude and longitude (length(combinations) is 79047)
  
  combinations <- rain %>% select(Latitude, Longitude) %>% distinct() %>% arrange(Latitude, Longitude)
  
  # For each site, find the closest lat and lon coordinates in the MERA dataset
  # This is a bit tricky, as not all combinations of each lat and lon are present
  # Also, the distribution of lat and log is not quite "linear"
  # Therefore, I came up with the following approach:
  # First, calculate deviation from each single lon and lat combination from actual lon and lat per site
  # Second, find the smallest absolute deviation (given as sum(diff_lat,diff_lon))
  # Thrid, extract lon and lat combination with the smalles deviation from actual lon and lat per site
  
  sites$Lat_MERA <- rep(NA,nrow(sites))
  sites$Lon_MERA <- rep(NA,nrow(sites))
  
  for (k in 1:nrow(sites)) {
  
    diff_lat <- combinations$Latitude - sites$Latitude[k]
    diff_lon <- combinations$Longitude - sites$Longitude[k]
    
    diff_sum <- abs(diff_lat) + abs(diff_lon)
    min.index <- which.min(diff_sum)
    
    sites$Lat_MERA[k] <- combinations$Latitude[min.index]
    sites$Lon_MERA[k] <- combinations$Longitude[min.index]
    
  }
  
  if (sum(is.na(sites$Lat_MERA)) > 0) {
    cat("Problem with latitude coordinates for", month_year2)
  }
  
  if (sum(is.na(sites$Lon_MERA)) > 0) {
    cat("Problem with longitude coordinates for", month_year2)
  }
  
  # Plot IOI map 
  
  IOI.map <- map("worldHires", c("Ireland","UK:Northern Ireland"), 
                 fill = T, col = "grey80")
  title(main = month_year2)
  map.axes(las = 1)
  map.cities(label = T, minpop = 90000, cex = 1.5)
  map.cities(capitals = 1)
  
  # Plot actual lon and lat per site
  
  points(sites$Latitude ~ sites$Longitude, pch = 21, col = "red", bg = "yellow", cex = 1.2)
  
  # Plot MERA lon and lat per site
  
  points(sites$Lat_MERA ~ sites$Lon_MERA, pch = 24, col = "green", bg = "cyan", cex = 1.2)

  # Extract MERA data for each site
  
  for (l in 1:nrow(sites)) {
    
    rain.site <- rain %>% 
      filter(rain$Latitude == sites$Lat_MERA[l] & 
             rain$Longitude == sites$Lon_MERA[l])
    
    rain.site$Site <- sites$Site[l]
    rain.site$Region <- sites$Region[l]
    rain.site$Farm <- sites$Farm[l]
    
    if (l == 1) {
      rain.sites <- rain.site
    }
    
    if (l > 1) {
      rain.sites <- rbind(rain.sites,rain.site)
    }
  }

  # Save as .csv and .Rdata
  
  file.name.csv <- paste0(data.dir,"MERA_rain_all_sites_",month_year,".csv")
  file.name.RData <- paste0(data.dir,"MERA_rain_all_sites_",month_year,".RData")

  write.table(rain.sites, file.name.csv, sep = ",", row.names = F)
  save(rain.sites, file = file.name.RData)
  
  if (m == 1) {
    rain.2017 <- rain.sites
  }
  
  if (m > 1) {
    rain.2017 <- rbind(rain.2017,rain.sites)
  }

}

dev.off()

write.table(rain.2017, paste0(data.dir,"MERA_rain_all_sites_Apr_Oct_2017.csv"), sep = ",", row.names = F)
save(rain.2017, file = paste0(data.dir,"MERA_rain_all_sites_Apr_Oct_2017.RData"))

#===========================================
# Plot data

ggplot(rain.2017, aes(x = dataDate, y = Rainfall)) +
  geom_line() +
  facet_grid(Region ~ Farm)
