# Authors: Mark Emmerson and Maja Ilic
# Last update: November 2019

# This is a script to extract and clean the soil moisture and tempererature data 
# from the Grassland Resilience Year 1 All Ireland Drought Experiment.

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

# Raw data folder

mydir <- paste0(data.dir,"Raw data/")

#===========================================
# Packages

library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)

#===========================================
# Load each of the files in turn
# Add new columns Site, Region and Farm

#================================
# Border - Armagh
# Note: there is no data for the sites B2 and B4

B1 <- read.csv(paste0(mydir,"B1.csv"), header = T)
B3 <- read.csv(paste0(mydir,"B3.csv"), header = T)
B5 <- read.csv(paste0(mydir,"B5.csv"), header = T)

B1$Site <- "B1"
B3$Site <- "B3"
B5$Site <- "B5"

B1$Region <- "Border"
B3$Region <- "Border"
B5$Region <- "Border"

B1$Farm <- 1
B3$Farm <- 3
B5$Farm <- 5

#----------------
# Cork

C1 <- read.csv(paste0(mydir,"C1.csv"), header = T)
C2 <- read.csv(paste0(mydir,"C2.csv"), header = T)
C3 <- read.csv(paste0(mydir,"C3.csv"), header = T)
C4 <- read.csv(paste0(mydir,"C4.csv"), header = T)
C5 <- read.csv(paste0(mydir,"C5.csv"), header = T)

C1$Site <- "C1"
C2$Site <- "C2"
C3$Site <- "C3"
C4$Site <- "C4"
C5$Site <- "C5"

C1$Region <- "Cork"
C2$Region <- "Cork"
C3$Region <- "Cork"
C4$Region <- "Cork"
C5$Region <- "Cork"

C1$Farm <- 1
C2$Farm <- 2
C3$Farm <- 3
C4$Farm <- 4
C5$Farm <- 5

#----------------
# Dublin

D1 <- read.csv(paste0(mydir,"D1.csv"), header = T)
D2 <- read.csv(paste0(mydir,"D2.csv"), header = T)
D3 <- read.csv(paste0(mydir,"D3.csv"), header = T)
D4 <- read.csv(paste0(mydir,"D4.csv"), header = T)
D5 <- read.csv(paste0(mydir,"D5.csv"), header = T)

D1$Site <- "D1"
D2$Site <- "D2"
D3$Site <- "D3"
D4$Site <- "D4"
D5$Site <- "D5"

D1$Region <- "Dublin"
D2$Region <- "Dublin"
D3$Region <- "Dublin"
D4$Region <- "Dublin"
D5$Region <- "Dublin"

D1$Farm <- 1
D2$Farm <- 2
D3$Farm <- 3
D4$Farm <- 4
D5$Farm <- 5

#----------------
# Limerick

L1 <- read.csv(paste0(mydir,"L1.csv"), header = T)
L2 <- read.csv(paste0(mydir,"L2.csv"), header = T)
L3 <- read.csv(paste0(mydir,"L3.csv"), header = T)
L4 <- read.csv(paste0(mydir,"L4.csv"), header = T)
L5 <- read.csv(paste0(mydir,"L5.csv"), header = T)

L1$Site <- "L1"
L2$Site <- "L2"
L3$Site <- "L3"
L4$Site <- "L4"
L5$Site <- "L5"

L1$Region <- "Limerick"
L2$Region <- "Limerick"
L3$Region <- "Limerick"
L4$Region <- "Limerick"
L5$Region <- "Limerick"

L1$Farm <- 1
L2$Farm <- 2
L3$Farm <- 3
L4$Farm <- 4
L5$Farm <- 5

#================================================

##--------------##
##              ##
##   Cleaning   ##
##              ##
##--------------##

# Create a list with all datasets

df.list <- list("B1" = B1, "B3" = B3, "B5" = B5,
                "C1" = C1, "C2" = C2, "C3" = C3, "C4" = C4, "C5" = C5,
                "D1" = D1, "D2" = D2, "D3" = D3, "D4" = D4, "D5" = D5,
                "L1" = L1, "L2" = L2, "L3" = L3, "L4" = L4, "L5" = L5)

# Remove all individual datasets

rm(B1,B3,B5,
   C1,C2,C3,C4,C5,
   D1,D2,D3,D4,D5,
   L1,L2,L3,L4,L5)

# Remove the row containing the units from each dataset (where necessary)
# Additionally, remove all rows with empty cells in "Label"

for (i in 1:length(df.list)){
  if (df.list[[i]][1,1] == "Units") {
    df.list[[i]] <- df.list[[i]][-1,]
    cat("Removed first row from",names(df.list)[[i]],"\n")
  }
  
  empty <- which(df.list[[i]]$Label == "")
  
  if(length(empty) > 0){
    df.list[[i]] <- df.list[[i]][-empty,]
    cat("Removed",length(empty),"empty rows from",names(df.list)[[i]],"\n")
  }
}


#================================================
# Find all the NANs, -INF and +INF error values in each dataset and replace with NA's

for (i in 1:length(df.list)){
  ## replace #NAN
  
  ind1 <- df.list[[i]] == "#NAN"
  df.list[[i]][ind1] <- NA
  
  if(!is.na(sum(ind1,na.rm = T))){
    cat("Replaced",sum(ind1,na.rm = T),"#NAN with NAs in",names(df.list)[[i]],"\n")
  }
  
  # replace #-INF
  
  ind2 <- df.list[[i]] == "#-INF"
  df.list[[i]][ind2] <- NA
  
  if(!is.na(sum(ind2,na.rm = T))){
    cat("Replaced",sum(ind2,na.rm = T),"#-INF with NAs in",names(df.list)[[i]],"\n")
  }
  
  # replace #+INF
  
  ind3 <- df.list[[i]] == "#+INF"
  df.list[[i]][ind3] <- NA
  
  if(!is.na(sum(ind3,na.rm = T))){
    cat("Replaced",sum(ind3,na.rm = T),"#+INF with NAs in",names(df.list)[[i]],"\n")
  }
}

#================================================
# Change the names of the column headings
# Change the column date.time to POSIXct
# Change all factors to numeric variables

an <- as.numeric   ## shorter
ac <- as.character ## shorter

# See file Logger_connections.xlsx
# Set 1: connections as planned

set1 <- c("B1", "B5", 
          "C2", "C3", "C4", "C5", 
          "D2", "D3", "D4", "D5", 
          "L1", "L2", "L3", "L5")

col.headings.1 <- c("date.time", "Power", 
                    "MoistureD1", "TemperatureD1", "MoistureC1", "TemperatureC1", "MoistureD2", "TemperatureD2",
                    "MoistureC2", "TemperatureC2", "MoistureC3", "TemperatureC3", "MoistureD3", "TemperatureD3",
                    "Site", "Region", "Farm")

## Set 2: Exceptions for B3 and L4

set2 <- c("B3","L4")

col.headings.2 <- c("date.time", "Power", 
                    "MoistureD1", "TemperatureD1", "MoistureC1", "TemperatureC1", "MoistureD2", "TemperatureD2",
                    "MoistureC3", "TemperatureC3", "MoistureC2", "TemperatureC2", "MoistureD3", "TemperatureD3",
                    "Site", "Region", "Farm")

## Set 3: Exceptions for C1 and D1

set3 <- c("C1","D1")

col.headings.3 <- c("date.time", "Power", 
                    "MoistureD1", "TemperatureD1", "MoistureC1", "TemperatureC1", "MoistureD2", "TemperatureD2",
                    "MoistureC2", "TemperatureC2", "MoistureD3", "TemperatureD3", "MoistureC3", "TemperatureC3",
                    "Site", "Region", "Farm")

dir.create(paste0(data.dir,"Cleaned Soil Moisture Data from Loggers"))
dir.create(paste0(data.dir,"Extracted Temperature Data from Loggers"))

#================================================
## for-loop

for (i in 1:length(df.list)){
  
  ## Change headings depending on the cable connections
  
  if(names(df.list)[i] %in% set1){
    names(df.list[[i]]) <- col.headings.1
  }

  if(names(df.list)[i] %in% set2){
    names(df.list[[i]]) <- col.headings.2
  }
  
  if(names(df.list)[i] %in% set3){
    names(df.list[[i]]) <- col.headings.3
  }
  
  cat("\nDataset",names(df.list)[[i]],"\n")
  
  ## Change date.time to class POSIXct
  
  if (i %in% c(1,4:18)){  ## excluding B3 and B5
    df.list[[i]]$date.time <- dmy_hms(df.list[[i]]$date.time)
  }
  
  if (i %in% c(2,3)){    ## B3 and B5
    df.list[[i]]$date.time <- dmy_hm(df.list[[i]]$date.time)
  }
  
  cat("Number of NAs in date.time:",sum(is.na(df.list[[i]]$date.time)),"\n")
  
  ## Change all variables from factor to numeric
  
  df.list[[i]]$Power <- an(ac(df.list[[i]]$Power))
  df.list[[i]]$MoistureD1 <- an(ac(df.list[[i]]$MoistureD1))
  df.list[[i]]$MoistureD2 <- an(ac(df.list[[i]]$MoistureD2))
  df.list[[i]]$MoistureD3 <- an(ac(df.list[[i]]$MoistureD3))
  df.list[[i]]$TemperatureD1 <- an(ac(df.list[[i]]$TemperatureD1))
  df.list[[i]]$TemperatureD2 <- an(ac(df.list[[i]]$TemperatureD2))
  df.list[[i]]$TemperatureD3 <- an(ac(df.list[[i]]$TemperatureD3))
  df.list[[i]]$MoistureC1 <- an(ac(df.list[[i]]$MoistureC1))
  df.list[[i]]$MoistureC2 <- an(ac(df.list[[i]]$MoistureC2))
  df.list[[i]]$MoistureC3 <- an(ac(df.list[[i]]$MoistureC3))
  df.list[[i]]$TemperatureC1 <- an(ac(df.list[[i]]$TemperatureC1))
  df.list[[i]]$TemperatureC2 <- an(ac(df.list[[i]]$TemperatureC2))
  df.list[[i]]$TemperatureC3 <- an(ac(df.list[[i]]$TemperatureC3))
  df.list[[i]]$site_ID <- paste0(df.list[[i]]$Region,df.list[[i]]$Farm)
  
  cat("Number of NAs: ",sum(is.na(df.list[[i]])),"\n")
  str(df.list[[i]])
  
  ## Rearrange the data (change format from wide to long)
  
  ##-------------------##
  ##   Soil moisture   ##
  ##-------------------##
  
  # Extract soil moisture and rearrange from wide to long format
  # Use the index of columns you want to keep
  
  site.moisture <- gather(df.list[[i]][,c(1,3,5,7,9,11,13,15:18)], 
                          key = "Plot",
                          value = "Moisture",
                          -date.time,-Site,-Region,-Farm,-site_ID)
  site.moisture$Plot <- str_remove(site.moisture$Plot,"Moisture")
  
  ##-----------------##
  ##   Temperature   ##
  ##-----------------##
  
  # Extract temperature and rearrange from wide to long format
  # Use the index of columns you want to keep
  
  site.temperature <- gather(df.list[[i]][,c(1,4,6,8,10,12,14,15:18)], 
                             key = "Plot",
                             value = "Temperature",
                             -date.time,-Site,-Region,-Farm,-site_ID)
  site.temperature$Plot <- str_remove(site.temperature$Plot,"Temperature")
  
  ##----------------##
  ##    Save data   ##
  ##----------------##
  
  ## Save each dataset in a new folder
  
  write.table(site.moisture, 
              file = paste0(data.dir,"Cleaned Soil Moisture Data from Loggers/",
                            names(df.list)[[i]]," - Soil Moisture Data from Loggers - cleaned.csv"),
              sep = ",",
              row.names = F)
  
  write.table(site.temperature, 
              file = paste0(data.dir,"Extracted Temperature Data from Loggers/",
                            names(df.list)[[i]]," - Temperature Data from Loggers.csv"),
              sep = ",",
              row.names = F)
  
  ##----------------##
  ##    Join data   ##
  ##----------------##
  
  ## Join all sites into one dataframe
  
  if (i == 1) {
    df.moisture <- site.moisture
    df.temperature <- site.temperature
  }
  
  if (i > 1) {
    df.moisture <- rbind(df.moisture,site.moisture)
    df.temperature <- rbind(df.temperature,site.temperature)
  }

}

#================================================
# Add column treatment

df.moisture$treatment <- str_remove(df.moisture$Plot,"[1,2,3]")
df.temperature$treatment <- str_remove(df.temperature$Plot,"[1,2,3]")

#================================================
# Save each (cleaned) dataset in a new folder

# Soil Moisture

write.table(df.moisture, paste0(data.dir,"Cleaned_Soil_Moisture_Data_from_Loggers.csv"),
            sep = ",", row.names = F)
save(df.moisture, file = paste0(data.dir,"Cleaned_Soil_Moisture_Data_from_Loggers.RData"))

# Temperature

write.table(df.temperature, paste0(data.dir,"Extracted_Temperature_Data_from_Loggers.csv"),
            sep = ",", row.names = F)
save(df.temperature, file = paste0(data.dir,"Extracted_Temperature_Data_from_Loggers.RData"))