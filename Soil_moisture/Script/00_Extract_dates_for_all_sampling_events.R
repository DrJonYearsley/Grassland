#############################################################################################
## This is a script to extract dates for all sampling events from the Grassland Resilience ##
## Year 1 All Ireland Drought Experiment                                                   ##
## The dates are derived from Lupe's Excel sheet (FIELD SCHEDULE WITH SITES - Internship)  ##
## Note: all dates are from the sheet "UPDATED SCHEDULE 19.7.17"                           ##          
##                                                                                         ##
## Author: Maja Ilic                                                                       ##
## Last update: 16 Mar 2020                                                                ##
#############################################################################################

#===========================================
# Clear objects from the workspace ----

rm(list = ls())

#===========================================
# Set working directory ----

# Maja's desktop PC

setwd("C:/Users/3054311/Documents/My Documents/Grassland project/")

#===========================================
# Packages ----

library(ggplot2)
library(lubridate)
library(dplyr)
library(tidyr)
library(formattable)
library(htmltools)
library(webshot)

#===========================================
# Load sampling dates ----

dates <- read.csv("SamplingDates.csv", sep = ",", header = T)
head(dates)
str(dates)

dates$Date <- dmy(dates$Date)

#================================================
# Every sampling event plotted for each site ----

sampling.plot <- ggplot(dates, 
                      aes(x = Date, 
                          y = factor(Site, rev(levels(factor(Site)))))) + 
  geom_point(aes(color = as.factor(Day), fill = as.factor(Day)), 
             size = 4, shape = 21, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 14, face = "bold")) +
  geom_hline(aes(yintercept = 5.5), color = "gray30") +
  geom_hline(aes(yintercept = 10.5), color = "gray30") +
  geom_hline(aes(yintercept = 15.5), color = "gray30") +
  guides(colour = guide_legend(title = "Day", order = 1),
         fill = guide_legend(title = "Day",order = 1)) +
  scale_x_date(date_breaks = "1 weeks",
               date_labels = "%d-%b") +
  labs(y = "Site", title = "Grassland Resilience - Sampling events")

ggsave("All sampling events per site.png", 
       sampling.plot,
       width = 25, height = 18, units = "cm")

#================================================
# Calculate the difference between each sampling event to Day 0 ----

EndDrought <- dates %>% 
  group_by(Site) %>% 
  filter(Day == 0)

names(EndDrought)[names(EndDrought) == "Date"] <- "Day0"

dates <- left_join(dates,EndDrought[,c("Region","Farm","Site","site_ID","Day0")])

dates$DayDiff <- dates$Date - dates$Day0

mytable <- formattable(dates[which(dates$Day != dates$DayDiff),], 
                       align = c("l","c","c","l","c","r","c","r"))
mytable

# webshot::install_phantomjs()

export_formattable <- function(f, file, width = "100%", height = NULL, 
                               background = "white", delay = 0.2)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}

export_formattable(mytable, "MismatchingDates.png") 
