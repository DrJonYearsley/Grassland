# Author: Maja Ilic
# Last update: November 2019

# This is a script to extract dates from all events from the Grassland Resilience 
# Year 1 All Ireland Drought Experiment

#===========================================
# Clear objects from the workspace

rm(list = ls())

#===========================================
# Set working directory

# Maja's desktop PC

setwd("C:/Users/3054311/Documents/My Documents/Grassland project/01_Probes_soil_moisture_and_temp")

data.dir <- paste0(getwd(),"/Data/")

#===========================================
# Packages

library(ggplot2)
library(lubridate)
library(dplyr)
library(tidyr)

#===========================================
# Load probe data (soil moisture) 

load(paste0(data.dir,"Cleaned_Soil_Moisture_Data_from_Loggers.RData"))
names(df.moisture)

#================================================
# Datalogger start and end dates per site and plot

datalogger.dates <- df.moisture %>% 
  group_by(Site,Plot) %>% 
  summarize(Start.Date = lubridate::date(min(date.time)),
            End.Date = lubridate::date(max(date.time)))

write.table(datalogger.dates, "DataLogger start and end dates per site and plot.csv",
            sep = ",", row.names = F)

# Datalogger start and end dates per site

datalogger.dates2 <- df.moisture %>% 
  group_by(Site) %>% 
  summarize(Start.Date = lubridate::date(min(date.time)),
            End.Date = lubridate::date(max(date.time)))

write.table(datalogger.dates2, "DataLogger start and end dates per site.csv",
            sep = ",", row.names = F)

#================================================
# Every event plotted for each site

events <- read.csv(paste0(data.dir,"ExperimentDatesMI.csv"), sep = ",", header = T)
head(events)
str(events)

events.long <- gather(events, key = "Event", value = "Date",
                      -Region, -Site, -Latitude, -Longitude)

events.long$Date <- mdy(events.long$Date)

##

DataLogger <- c("DataloggerStart","DataloggerStart_actual","DataloggerEnd_actual")
Drought <- c("StartDrought","EndDrought")
Sampling <- c("Day8","Day16","Day32","Day64")
Flight <- c("Flight1","Flight2","Flight3","Flight4","Flight5","Flight6","Flight7")

events.long[which(events.long$Event %in% DataLogger),"Group"] <- "DataLogger"
events.long[which(events.long$Event %in% Drought),"Group"] <- "Drought"
events.long[which(events.long$Event %in% Sampling),"Group"] <- "Sampling"
events.long[which(events.long$Event %in% Flight),"Group"] <- "Flight"

str(events.long)

events.long$Event <- as.factor(events.long$Event)
events.long$Event <- factor(events.long$Event, levels = c(DataLogger,Drought,Sampling,Flight))
events.long$Group <- as.factor(events.long$Group)
events.long$Group <- factor(events.long$Group, levels = c("DataLogger","Drought","Sampling","Flight"))

str(events.long)

save(events,events.long,file = paste0(data.dir,"All events.RData"))

## Plot

events.plot <- ggplot(events.long, 
                      aes(x = Date, 
                          y = factor(Site, rev(levels(factor(Site)))))) + 
  geom_point(aes(color = Event, fill = Event, shape = Group), 
             size = 4, alpha = 0.5, stroke = 1.2) +
  scale_shape_manual(values = c(21,22,24,4)) + 
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold")) +
  geom_hline(aes(yintercept = 5.5), color = "gray30") +
  geom_hline(aes(yintercept = 10.5), color = "gray30") +
  geom_hline(aes(yintercept = 15.5), color = "gray30") +
  guides(colour = guide_legend(order = 1),
         fill = guide_legend(order = 1),
         shape = guide_legend(order = 2)) +
  ylab("Site")

events.plot

ggsave(paste0(figures.dir,"All events per site.png"), events.plot,
       width = 25, height = 18, units = "cm")

##

events.plot2 <- ggplot(events.long, 
                       aes(x = Date, 
                           y = factor(Site, rev(levels(factor(Site)))))) +
  geom_point(aes(color = Event, fill = Event, shape = Event), 
             size = 4, alpha = 0.5, stroke = 1.2) +
  scale_shape_manual(values = c(21,21,21,
                                22,22,
                                24,24,24,24,
                                4,4,4,4,4,4,4)) + 
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold")) +
  geom_hline(aes(yintercept = 5.5), color = "gray30") +
  geom_hline(aes(yintercept = 10.5), color = "gray30") +
  geom_hline(aes(yintercept = 15.5), color = "gray30") +
  ylab("Site")

events.plot2

ggsave(paste0(figures.dir,"All events per site 2.png"), events.plot2,
       width = 25, height = 18, units = "cm")

