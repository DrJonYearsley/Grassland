# Author: Maja Ilic
# Last update: November 2019

# This is a script to explore and plot average soil moisture data from the Grassland Resilience 
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

dir.create(paste0(figures.dir,"Average soil moisture per site and treatment"))
mydir <- paste0(figures.dir,"Average soil moisture per site and treatment/")

#================================================
# Calculate the average drought and control soil moistures per time point and site 

##########################
##                      ##
##   Average moisture   ##
##                      ##
##########################

avg.moisture <- df %>% 
  group_by(date.time, Site, Region, Farm, treatment) %>% 
  summarize(Mean = mean(Moisture, na.rm = T),
            SD = sd(Moisture, na.rm = T),
            n = sum(!is.na(Moisture))) 

head(avg.moisture)
str(avg.moisture)

save(avg.moisture, file = paste0(data.dir,"Average_soil_moisture_per_site_and_treatment.RData"))

#================================================
# x limits

min.date <- min(df$date.time)
max.date <- max(df$date.time)

#================================================
# Plot average moisture +- SD

p1 <- ggplot(avg.moisture, aes(x = date.time, y = Mean)) +
  geom_line(aes(color = treatment)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD,fill = treatment), 
              alpha = 0.3) +
  facet_grid(Region ~ Farm) +
  labs(x = "Date", y = "Soil Moisture (%)", title = "Average moisture per site and treatment using raw data") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 14, face = "bold")) +
  scale_color_manual(values = c("red","blue"), name = "Treatment") +
  scale_fill_manual(values = c("red","blue"), name = "Treatment") +
  scale_y_continuous(limits = c(-2,102), breaks = c(0,25,50,75,100)) +
  scale_x_datetime(limits = c(min.date,max.date))

p1

ggsave(paste0(mydir,"All plots - Average soil moisture using raw data.png"),
       p1,
       width = 24, height = 18, units = "cm")

#================================================
# Visualizing mean soil moisture with missing data (less than 3 datapoints per time point)

##----------------------------------##
##                                  ##
##     Plot with missing values     ##
##   for all sites using for-loop   ##
##                                  ##
##----------------------------------##

sites <- unique(avg.moisture$Site)

Site <- c()
All.Rep.C <- c()
Missing.Rep.C <- c()
All.Rep.D <- c()
Missing.Rep.D <- c()

m <- 1

for (site in sites){
  
  mysite <- avg.moisture %>% 
    filter(Site == site) %>% 
    droplevels()
  
  control.avg <- mysite %>% 
    filter(treatment == "C") %>% 
    droplevels()
  
  drought.avg <- mysite %>% 
    filter(treatment == "D") %>% 
    droplevels()
  
  if((nrow(control.avg) + nrow(drought.avg)) != nrow(mysite)){
    cat("Problem with site",site,"\nControl and Drought subsets don't match")
  }
  
  Site[m] <- site
  
  start.drought <- as.POSIXct(mdy(events$StartDrought[events$Site == site]), tz = "GMT") + 11*60*60
  end.drought <- as.POSIXct(mdy(events$EndDrought[events$Site == site]), tz = "GMT") + 11*60*60
  day8 <- as.POSIXct(mdy(events$Day8[events$Site == site]), tz = "GMT") + 11*60*60
  day16 <- as.POSIXct(mdy(events$Day16[events$Site == site]), tz = "GMT") + 11*60*60
  day32 <- as.POSIXct(mdy(events$Day32[events$Site == site]), tz = "GMT") + 11*60*60
  day64 <- as.POSIXct(mdy(events$Day64[events$Site == site]), tz = "GMT") + 11*60*60
  
  
  ##-------------------##
  ##   Control plots   ##
  ##-------------------##
  
  # Extract all rows with n = 3 (all replicates)
  all.rep.C <- control.avg %>% 
    filter(n == 3)
  
  # Extract all rows with n < 3 (missing values, i.e. replicates)
  missing.rep.C <- control.avg %>% 
    filter(n != 3)
  
  if((nrow(all.rep.C) + nrow(missing.rep.C)) != nrow(control.avg)){
    cat("Problem with site",site,"\nControl subsets don't match")
  }
  
  All.Rep.C[m] <- nrow(all.rep.C)
  Missing.Rep.C[m] <- nrow(missing.rep.C)
  
  # Left join with all dates from the original subset
  # This will create NAs for all dates that are not present in all.rep and missing.rep
  
  all.rep.C <- left_join(control.avg[,1], all.rep.C)
  missing.rep.C <- left_join(control.avg[,1], missing.rep.C)
  
  ## Plot data
  
  title.C <- paste(mysite$Region[1],mysite$Farm[1],"- Control")
  
  plot.c <- ggplot() + 
    geom_rect(data = mysite,
              xmin = start.drought,
              xmax = end.drought,
              ymin = -Inf,
              ymax = +Inf,
              fill = "gray80") +
    geom_line(data = all.rep.C, 
              aes(x = date.time, y = Mean),
              size = 1, color = "red") +
    geom_line(data = missing.rep.C, 
              aes(x = date.time, y = Mean),
              size = 0.8, color = "indianred1") +
    geom_ribbon(data = all.rep.C,
                aes(x = date.time, 
                    ymin = Mean - SD,
                    ymax = Mean + SD),
                fill = "red", alpha = 0.6) +
    geom_ribbon(data = missing.rep.C,
                aes(x = date.time, 
                    ymin = Mean - SD,
                    ymax = Mean + SD),
                fill = "indianred1", alpha = 0.4) + 
    geom_vline(xintercept = day8, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = day16, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = day32, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = day64, linetype = "dashed", color = "gray50") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 12),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
    labs(x = "", y = "", title = title.C) +
    scale_x_datetime(limits = c(min.date,max.date)) +
    scale_y_continuous(limits = c(-2,102), breaks = c(0,25,50,75,100))
  
  
  ##-------------------##
  ##   Drought plots   ##
  ##-------------------##
  
  # Extract all rows with n = 3 (all replicates)
  all.rep.D <- drought.avg %>% 
    filter(n == 3)
  
  # Extract all rows with n < 3 (missing values, i.e. replicates)
  missing.rep.D <- drought.avg %>% 
    filter(n != 3)
  
  if((nrow(all.rep.D) + nrow(missing.rep.D)) != nrow(drought.avg)){
    cat("Problem with site",site,"\nDrought subsets don't match")
  }
  
  All.Rep.D[m] <- nrow(all.rep.D)
  Missing.Rep.D[m] <- nrow(missing.rep.D)
  
  # Left join with all dates from the original subset
  # This will create NAs for all dates that are not present in all.rep and missing.rep
  
  all.rep.D <- left_join(drought.avg[,1], all.rep.D)
  missing.rep.D <- left_join(drought.avg[,1], missing.rep.D)
  
  # Plot data
  
  title.D <- paste(mysite$Region[1],mysite$Farm[1],"- Drought")
    
  plot.d <- ggplot() +
    geom_rect(data = mysite,
              xmin = start.drought,
              xmax = end.drought,
              ymin = -Inf,
              ymax = +Inf,
              fill = "gray80") +
    geom_line(data = all.rep.D, 
              aes(x = date.time, y = Mean),
              size = 1, color = "blue") +
    geom_line(data = missing.rep.D, 
              aes(x = date.time, y = Mean),
              size = 0.8, color = "dodgerblue1") +
    geom_ribbon(data = all.rep.D,
                aes(x = date.time, 
                    ymin = Mean - SD,
                    ymax = Mean + SD),
                fill = "blue", alpha = 0.5) +
    geom_ribbon(data = missing.rep.D,
                aes(x = date.time, 
                    ymin = Mean - SD,
                    ymax = Mean + SD),
                fill = "dodgerblue1", alpha = 0.5) + 
    geom_vline(xintercept = day8, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = day16, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = day32, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = day64, linetype = "dashed", color = "gray50") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 12),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
    labs(x = "", y = "", title = title.D) +
    scale_x_datetime(limits = c(min.date,max.date)) +
    scale_y_continuous(limits = c(-2,102), breaks = c(0,25,50,75,100))
  
  
  ##------------------##
  ##   Joined plots   ##
  ##------------------##
  
  p <- plot_grid(plot.c, plot.d,
                 ncol = 2,
                 labels = NULL)
  
  y.grob <- textGrob("Soil Moisture (%), mean \u00B1 s.d.", 
                     gp = gpar(fontsize = 14, fontface = "bold"), rot = 90)
  
  x.grob <- textGrob("Date", 
                     gp = gpar(fontsize = 14, fontface = "bold"))
  
  # add to plot
  
  figure.title.1 <- paste0(mydir,site," - Average soil moisture joined.png")
  
  ggsave(figure.title.1,
         grid.arrange(arrangeGrob(p, left = y.grob, bottom = x.grob)),
         width = 22, height = 10, units = "cm")
  
  
  ##------------------##
  ##   Single panel   ##
  ##------------------##
  
  main.title <- c(mysite$Region[1]," ",mysite$Farm[1]," - ","Control"," vs ","Drought")
  title.colors <- c("black","black","black","black","red","black","blue")
  
  single.panel <- ggplot() + 
    geom_rect(data = mysite,
              xmin = start.drought,
              xmax = end.drought,
              ymin = -Inf,
              ymax = +Inf,
              fill = "gray80") +
    geom_line(data = all.rep.C, 
              aes(x = date.time, y = Mean),
              size = 1, color = "red") +
    geom_line(data = missing.rep.C, 
              aes(x = date.time, y = Mean),
              size = 0.8, color = "indianred1") +
    geom_ribbon(data = all.rep.C,
                aes(x = date.time, 
                    ymin = Mean - SD,
                    ymax = Mean + SD),
                fill = "red", alpha = 0.6) +
    geom_ribbon(data = missing.rep.C,
                aes(x = date.time, 
                    ymin = Mean - SD,
                    ymax = Mean + SD),
                fill = "indianred1", alpha = 0.4) +
    geom_line(data = all.rep.D, 
              aes(x = date.time, y = Mean),
              size = 1, color = "blue") +
    geom_line(data = missing.rep.D, 
              aes(x = date.time, y = Mean),
              size = 0.8, color = "dodgerblue1") +
    geom_ribbon(data = all.rep.D,
                aes(x = date.time, 
                    ymin = Mean - SD,
                    ymax = Mean + SD),
                fill = "blue", alpha = 0.5) +
    geom_ribbon(data = missing.rep.D,
                aes(x = date.time, 
                    ymin = Mean - SD,
                    ymax = Mean + SD),
                fill = "dodgerblue1", alpha = 0.5) + 
    geom_vline(xintercept = day8, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = day16, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = day32, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = day64, linetype = "dashed", color = "gray50") +
    mytheme +
    labs(x = "Date", 
         y = "Soil Moisture (%), mean \u00B1 s.d.") +
    scale_x_datetime(limits = c(min.date,max.date), 
                     breaks = c(start.drought,
                                end.drought),
                     date_labels = "%b %d") +
    scale_y_continuous(limits = c(-2,102), breaks = c(0,25,50,75,100))
  
  figure.title.2 <- paste0(mydir,site," - Average soil moisture.png")
  
  ggsave(figure.title.2,
         grid.arrange(single.panel,
                      top = tableGrob(t(main.title),
                                      theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                             base_colour = title.colors,
                                                             base_size = 16))),
         width = 15, height = 12, units = "cm") 
  
  # Remove axis titles and assign each plot to a named object
  
  title <- paste(mysite$Region[1],mysite$Farm[1])
  
  single.panel <- single.panel +
    theme(axis.title = element_blank(),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
    labs(title = title)
  
  assign(site,single.panel)
  
  rm(single.panel)
  
  # Adjust running index m
  
  m <- m + 1
}

#================================================
# Export all missing values

mydata <- data.frame(Site,All.Rep.C,Missing.Rep.C,All.Rep.D,Missing.Rep.D)
write.table(mydata, "Number of full replicated treatments per site moisture.csv", sep = ",", row.names = F)

#================================================
# Final plot

plot.final <- plot_grid(B1,NULL,B3,NULL,B5,
                        C1,C2,C3,C4,C5,
                        D1,D2,D3,D4,D5,
                        L1,L2,L3,L4,L5,
                        ncol = 5,
                        labels = NULL)

y.grob <- textGrob("Soil Moisture (%), mean \u00B1 s.d.", 
                   gp = gpar(fontsize = 14, fontface = "bold"), rot = 90)

x.grob <- textGrob("Date", 
                   gp = gpar(fontsize = 14, fontface = "bold"))

figure.title.3 <- paste0(mydir,"Average soil moisture - all sites.png")

ggsave(figure.title.3,
       grid.arrange(arrangeGrob(plot.final, left = y.grob, bottom = x.grob)),
       width = 18, height = 14)