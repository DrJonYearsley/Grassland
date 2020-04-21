# Author: Maja Ilic
# Last update: December 2019

# This is a script to explore datalogger temperature data from the Grassland Resilience 
# Year 1 All Ireland Drought Experiment. 
# Outliers are removed and the cleaned data is exported. 
# Plots are generated for each site and plot.

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

library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)

#===========================================

load(file = paste0(data.dir,"Extracted_Temperature_Data_from_Loggers.RData"))

#================================================
# Explore data to identify outliers for each site

# Start and end date

min.date <- min(df.temperature$date.time)
max.date <- max(df.temperature$date.time)

# Create directory to save plots

dir.create(paste0(figures.dir,"Temperature per site and plot"))
mydir <- paste0(figures.dir,"Temperature per site and plot/")

#================================================

##---------------##
##               ##
##   Plot data   ##
##               ##
##---------------##

## Lines

p1 <- ggplot(df.temperature, aes(x = date.time, y = Temperature)) +
  geom_line(aes(color = treatment, linetype = Plot)) +
  facet_grid(Region ~ Farm, scales = "free_y") +
  labs(x = "Date", y = "Temperature (°C)", title = "Temperature per site and plot") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold")) +
  scale_x_datetime(limits = c(min.date,max.date)) +
  geom_hline(aes(yintercept = 30), linetype = "dashed")

p1

ggsave(paste0(mydir,"Temperature per site and plot - extracted data from DataLoggers.png"),
       p1,
       width = 13, height = 10)

## Boxplots

p2 <- ggplot(df.temperature, aes(x = treatment, y = Temperature)) +
  geom_boxplot(aes(fill = treatment),
               outlier.shape = 21) +
  facet_grid(Region ~ Farm, scales = "free_y") +
  labs(x = "Date", y = "Temperature (°C)", title = "Temperature per site and plot") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold")) +
  geom_hline(aes(yintercept = 30), linetype = "dashed")

p2

ggsave(paste0(mydir,"Temperature per site and plot - extracted data from DataLoggers - boxplots.png"),
       p2,
       width = 13, height = 10)

## Boxplots with raw data

p3 <- ggplot(df.temperature, aes(x = treatment, y = Temperature, fill = treatment, color = treatment)) +
  geom_point(size = 2, shape = 21, alpha = 0.5,
             position = position_jitterdodge()) +
  facet_grid(Region ~ Farm, scales = "free_y") +
  labs(x = "Date", y = "Temperature (°C)", title = "Temperature per site and plot") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold")) +
  geom_hline(aes(yintercept = 30), linetype = "dashed")

p3

ggsave(paste0(mydir,"Temperature per site and plot - extracted data from DataLoggers - jitter.png"),
       p3,
       width = 13, height = 10)

#================================================
# Cleaning the data - step 1
# This has been done per site individually, by visual inspection of the raw data
# All outliers are saved in the vector "out" and replaced with "NA" in the original dataset

##--------##
##   B1   ##
##--------##

# remove temp over 30 °C

out <- which(df.temperature$Site == "B1" & df.temperature$Temperature > 30)
temp <- df.temperature$Temperature[out]
site <- rep("B1",length(out))

df.outlier <- data.frame(site,out,temp)
df.outlier.all <- df.outlier

df.temperature$Temperature[out] <- NA

##--------##
##   B3   ##
##--------##

# remove temp over 25 °C

out <- which(df.temperature$Site == "B3" & df.temperature$Temperature > 25)
temp <- df.temperature$Temperature[out]
site <- rep("B3",length(out))

df.outlier <- data.frame(site,out,temp)
df.outlier.all <- rbind(df.outlier.all, df.outlier)

df.temperature$Temperature[out] <- NA

##--------##
##   B5   ##
##--------##

# remove temp over 30 °C

out <- which(df.temperature$Site == "B5" & df.temperature$Temperature > 30)
temp <- df.temperature$Temperature[out]
site <- rep("B5",length(out))

df.outlier <- data.frame(site,out,temp)
df.outlier.all <- rbind(df.outlier.all, df.outlier)

df.temperature$Temperature[out] <- NA

##--------##
##   C2   ##
##--------##

# remove temp over 50 °C

out <- which(df.temperature$Site == "C2" & df.temperature$Temperature > 50)
temp <- df.temperature$Temperature[out]
site <- rep("C2",length(out))

df.outlier <- data.frame(site,out,temp)
df.outlier.all <- rbind(df.outlier.all, df.outlier)

df.temperature$Temperature[out] <- NA

## Here, data will bel cleaned later on as there are still some possible outliers included in the raw data

##--------##
##   C3   ##
##--------##

# remove temp over 30 °C and below 0 °C

out <- which(df.temperature$Site == "C3" & df.temperature$Temperature > 30 | df.temperature$Site == "C3" &  df.temperature$Temperature < 0)
temp <- df.temperature$Temperature[out]
site <- rep("C3",length(out))

df.outlier <- data.frame(site,out,temp)
df.outlier.all <- rbind(df.outlier.all, df.outlier)

df.temperature$Temperature[out] <- NA

##--------##
##   C4   ##
##--------##

# remove temp below 0 °C

out <- which(df.temperature$Site == "C4" &  df.temperature$Temperature < 0)
temp <- df.temperature$Temperature[out]
site <- rep("C4",length(out))

df.outlier <- data.frame(site,out,temp)
df.outlier.all <- rbind(df.outlier.all, df.outlier)

df.temperature$Temperature[out] <- NA

##--------##
##   D2   ##
##--------##

# remove temp below 0 °C

out <- which(df.temperature$Site == "D2" &  df.temperature$Temperature < 0)
temp <- df.temperature$Temperature[out]
site <- rep("D2",length(out))

df.outlier <- data.frame(site,out,temp)
df.outlier.all <- rbind(df.outlier.all, df.outlier)

df.temperature$Temperature[out] <- NA                        

##--------##
##   L2   ##
##--------##

# remove temp over 30 °C and below 0°

out <- which(df.temperature$Site == "L2" & df.temperature$Temperature > 30 | df.temperature$Site == "L2" &  df.temperature$Temperature < 0)
temp <- df.temperature$Temperature[out]
site <- rep("L2",length(out))

df.outlier <- data.frame(site,out,temp)
df.outlier.all <- rbind(df.outlier.all, df.outlier)

df.temperature$Temperature[out] <- NA

## Export outliers

df.outlier.all$date.time <- df.temperature$date.time[df.outlier.all$out]

write.table(df.outlier.all, "All temperature outliers - visual inspection.csv", sep = ",", row.names = F)

#================================================
# Plot cleaned data

## Lines

p4 <- ggplot(df.temperature, aes(x = date.time, y = Temperature)) +
  geom_line(aes(color = treatment, linetype = Plot)) +
  facet_grid(Region ~ Farm, scales = "free_y") +
  labs(x = "Date", y = "Temperature (°C)", title = "Temperature per site and plot") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold")) +
  scale_x_datetime(limits = c(min.date,max.date)) +
  geom_hline(aes(yintercept = 30), linetype = "dashed")

p4

ggsave(paste0(mydir,"Temperature per site and plot - cleaned.png"),
       p4,
       width = 13, height = 10)

## Boxplots

p5 <- ggplot(df.temperature, aes(x = treatment, y = Temperature)) +
  geom_boxplot(aes(fill = treatment),
               outlier.shape = 21) +
  facet_grid(Region ~ Farm, scales = "free_y") +
  labs(x = "Date", y = "Temperature (°C)", title = "Temperature per site and plot") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold")) +
  geom_hline(aes(yintercept = 30), linetype = "dashed")

p5

ggsave(paste0(mydir,"Temperature per site and plot - boxplots - cleaned.png"),
       p5,
       width = 13, height = 10)

## Boxplots with raw data

p6 <- ggplot(df.temperature, aes(x = treatment, y = Temperature, fill = treatment, color = treatment)) +
  geom_point(size = 2, shape = 21, alpha = 0.5,
             position = position_jitterdodge()) +
  facet_grid(Region ~ Farm, scales = "free_y") +
  labs(x = "Date", y = "Temperature (°C)", title = "Temperature per site and plot") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold"))

p6

ggsave(paste0(mydir,"Temperature per site and plot - jitter - cleaned.png"),
       p6,
       width = 13, height = 10)

#================================================
# Cleaning the data - step 2
# Plot temperature for each plot in different color, to identify plots with additional potential outliers

ggplot(df.temperature, aes(x = date.time, y = Temperature)) +
  geom_line(aes(color = Plot)) +
  facet_grid(Region ~ Farm, scales = "free_y") +
  labs(x = "Date", y = "Temperature (°C)", title = "Temperature per site and plot") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold")) +
  scale_x_datetime(limits = c(min.date,max.date)) +
  geom_hline(aes(yintercept = 30), linetype = "dashed")

# It seems like there might be some additional outliers in Cork 2, Dublin 5 and Limerick 2
# Let's plot these separatelly

Cork2 <- df.temperature %>% 
  filter(Site == "C2") %>% 
  droplevels()

Dublin5 <- df.temperature %>% 
  filter(Site == "D5") %>% 
  droplevels()

Limerick2 <- df.temperature %>% 
  filter(Site == "L2") %>% 
  droplevels()

################
##            ##
##   Cork 2   ##
##            ##
################

cork2.plot1 <- ggplot(Cork2, aes(x = date.time, y = Temperature)) +
  geom_line(aes(color = Plot)) +
  facet_wrap(~Plot) +
  labs(x = "Date", y = "Temperature (°C)", title = "Cork 2 - Temperature per plot") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold")) +
  scale_x_datetime(limits = c(min.date,max.date)) +
  geom_hline(aes(yintercept = 30), linetype = "dashed")

cork2.plot1

ggsave(paste0(mydir,"Cork 2 - Temperature per plot - cleaned 1.png"),
       cork2.plot1,
       width = 8, height = 6)

## C3 and all drought plots look good, no additional potential outliers are present
## C1 has additional outliers over 30 °C:

out <- which(df.temperature$Site == "C2" & df.temperature$Plot == "C1" & df.temperature$Temperature > 30)
temp <- df.temperature$Temperature[out]
site <- rep("C2",length(out))
plot <- rep("C1",length(out))

df.outlier <- data.frame(site,out,temp,plot)
df.outlier.all <- df.outlier

df.temperature$Temperature[out] <- NA

# Plot again

df.temperature %>% 
  filter(Site == "C2" & Plot == "C1") %>% 
  ggplot(aes(x = date.time, y = Temperature)) + 
  geom_line() +
  scale_y_continuous(limits = c(0,30)) +
  theme_bw() ## looks much more reasonable!

## C2 has also some outliers, that seem to follow the pattern of other control plots, but are much lower
## Trying to identify the whole period by looking for largest differences from timepoint t to timepoint t+1
## I will also exclude all NAs from the dataset for this

Cork2.C2 <- Cork2 %>% 
  filter(Plot == "C2" & !is.na(Temperature)) %>% 
  droplevels()

diff <- c()

for (i in 1:(nrow(Cork2.C2)-1)) {
  
  diff[i] <- Cork2.C2$Temperature[i+1] - Cork2.C2$Temperature[i]
  
}

# Find the highest difference and check the temperature and date for t and t+1

diff[which.max(abs(diff))]
which.max(abs(diff))
Cork2.C2$Temperature[which.max(abs(diff))]
Cork2.C2$Temperature[which.max(abs(diff))+1]
Cork2.C2$date.time[which.max(abs(diff))]
Cork2.C2$date.time[which.max(abs(diff))+1]

df.temperature %>% 
  filter(Site == "C2" & Plot == "C2") %>% 
  ggplot(aes(x = date.time, y = Temperature)) + 
  geom_line() +
  scale_y_continuous(limits = c(0,30)) +
  theme_bw() +
  geom_vline(xintercept = Cork2.C2$date.time[which.max(abs(diff))], col = "red") 

# Right border of the outlier period identified!

right.border <- Cork2.C2$date.time[which.max(abs(diff))]

# Exclude the highest difference and search for the next highest

diff <- diff[-which.max(abs(diff))]

diff[which.max(abs(diff))]
which.max(abs(diff))
Cork2.C2$Temperature[which.max(abs(diff))]
Cork2.C2$Temperature[which.max(abs(diff))+1]
Cork2.C2$date.time[which.max(abs(diff))]
Cork2.C2$date.time[which.max(abs(diff))+1]

df.temperature %>% 
  filter(Site == "C2" & Plot == "C2") %>% 
  ggplot(aes(x = date.time, y = Temperature)) + 
  geom_line() +
  scale_y_continuous(limits = c(0,30)) +
  theme_bw() +
  geom_vline(xintercept = Cork2.C2$date.time[which.max(abs(diff))], col = "red") 

# This seems to be a bit after the actual left border of the outlier period
# Therefore, I repeated the procedure:

diff <- diff[-which.max(abs(diff))]

diff[which.max(abs(diff))]
which.max(abs(diff))
Cork2.C2$Temperature[which.max(abs(diff))]
Cork2.C2$Temperature[which.max(abs(diff))+1]
Cork2.C2$date.time[which.max(abs(diff))]
Cork2.C2$date.time[which.max(abs(diff))+1]

df.temperature %>% 
  filter(Site == "C2" & Plot == "C2") %>% 
  ggplot(aes(x = date.time, y = Temperature)) + 
  geom_line() +
  scale_y_continuous(limits = c(0,30)) +
  theme_bw() +
  geom_vline(xintercept = Cork2.C2$date.time[which.max(abs(diff))], col = "red") 

# Left border of the outlier period identified!

left.border <- Cork2.C2$date.time[which.max(abs(diff))]

# Next, I replaced all outliers withing this period with NA and added the outliers to the outliers dataframe

out <- which(df.temperature$Site == "C2" & df.temperature$Plot == "C2" & df.temperature$date.time > left.border & df.temperature$date.time <= right.border)
temp <- df.temperature$Temperature[out]
site <- rep("C2",length(out))
plot <- rep("C2",length(out))

df.outlier <- data.frame(site,out,temp,plot)
df.outlier.all <- rbind(df.outlier.all, df.outlier)

df.temperature$Temperature[out] <- NA

# Plot again

df.temperature %>% 
  filter(Site == "C2" & Plot == "C2") %>% 
  ggplot(aes(x = date.time, y = Temperature)) + 
  geom_line() +
  scale_y_continuous(limits = c(0,30)) +
  theme_bw() ## all potential outliers were removed!

# Plot all data for Cork 2 again

Cork2.clean <- df.temperature %>% 
  filter(Site == "C2") %>% 
  droplevels()

cork2.plot2 <- ggplot(Cork2.clean, aes(x = date.time, y = Temperature)) +
  geom_line(aes(color = Plot)) +
  facet_wrap(~Plot) +
  labs(x = "Date", y = "Temperature (°C)", title = "Cork 2 - Temperature per plot") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold")) +
  scale_x_datetime(limits = c(min.date,max.date)) +
  geom_hline(aes(yintercept = 30), linetype = "dashed") + 
  scale_y_continuous(limits = c(0,30))

cork2.plot2

ggsave(paste0(mydir,"Cork 2 - Temperature per plot - cleaned 2.png"),
       cork2.plot2,
       width = 8, height = 6)


##################
##              ##
##   Dublin 5   ##
##              ##
##################

dublin5.plot1 <- ggplot(Dublin5, aes(x = date.time, y = Temperature)) +
  geom_line(aes(color = Plot)) +
  facet_wrap(~Plot) +
  labs(x = "Date", y = "Temperature (°C)", title = "Cork 2 - Temperature per plot") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold")) +
  scale_x_datetime(limits = c(min.date,max.date)) +
  geom_hline(aes(yintercept = 25), linetype = "dashed")

dublin5.plot1

ggsave(paste0(mydir,"Dublin 5 - Temperature per plot - cleaned 1.png"),
       dublin5.plot1,
       width = 8, height = 6)

## Here, it is difficult to say, if the temperature in control 1 was correctly recorded
## Therefore, no data will be removed except for obvious outliers at the beginning and the end of the measurement
## Find the date and time for the first recorded temperature for each plot (should be the same across plots)

start.date.time <- Dublin5 %>% 
  filter(!is.na(Temperature)) %>% 
  group_by(Plot) %>% 
  summarize(Start = min(date.time))

start.date.time

# Check temperature on these days in each plot
# For start, this is easy, as it is the same date & time for each plot

start.temp <- Dublin5 %>% 
  filter(date.time == start.date.time$Start[1])

start.temp  

# For the end of the measuring period, only the plot drought 3 has potential outliers towards the end of the recording period

Dublin5 %>% 
  filter(Plot == "D3") %>% 
  ggplot(aes(x = date.time, y = Temperature)) + 
  geom_line() +
  scale_y_continuous(limits = c(0,30)) +
  theme_bw() +
  geom_hline(yintercept = 15, color = "red")

# This is not the most perfect solution and is rather a "lucky guess", but I decided here 
# to exclude all values over 15°C in the late September

Dublin5$date.time[max(which(Dublin5$Temperature < 15))]
Dublin5$Temperature[max(which(Dublin5$Temperature < 15))]

# Save this as outlier.border

outlier.border <- Dublin5$date.time[max(which(Dublin5$Temperature < 15))]

# Exclude all first recordings from each plot
# Then, exclude outliers from the plot drought 3

## C1

out <- which(df.temperature$Site == "D5" & df.temperature$Plot == "C1" & df.temperature$date.time == start.date.time$Start[1])
temp <- df.temperature$Temperature[out]
site <- rep("D5",length(out))
plot <- rep("C1",length(out))  

df.outlier <- data.frame(site,out,temp,plot)
df.outlier.all <- rbind(df.outlier.all, df.outlier)

df.temperature$Temperature[out] <- NA

## C2

out <- which(df.temperature$Site == "D5" & df.temperature$Plot == "C2" & df.temperature$date.time == start.date.time$Start[2])
temp <- df.temperature$Temperature[out]
site <- rep("D5",length(out))
plot <- rep("C2",length(out))  

df.outlier <- data.frame(site,out,temp,plot)
df.outlier.all <- rbind(df.outlier.all, df.outlier)

df.temperature$Temperature[out] <- NA

## C3

out <- which(df.temperature$Site == "D5" & df.temperature$Plot == "C3" & df.temperature$date.time == start.date.time$Start[3])
temp <- df.temperature$Temperature[out]
site <- rep("D5",length(out))
plot <- rep("C3",length(out))  

df.outlier <- data.frame(site,out,temp,plot)
df.outlier.all <- rbind(df.outlier.all, df.outlier)

df.temperature$Temperature[out] <- NA

## D1

out <- which(df.temperature$Site == "D5" & df.temperature$Plot == "D1" & df.temperature$date.time == start.date.time$Start[4])
temp <- df.temperature$Temperature[out]
site <- rep("D5",length(out))
plot <- rep("D1",length(out))  

df.outlier <- data.frame(site,out,temp,plot)
df.outlier.all <- rbind(df.outlier.all, df.outlier)

df.temperature$Temperature[out] <- NA

## D2

out <- which(df.temperature$Site == "D5" & df.temperature$Plot == "D2" & df.temperature$date.time == start.date.time$Start[5])
temp <- df.temperature$Temperature[out]
site <- rep("D5",length(out))
plot <- rep("D2",length(out))  

df.outlier <- data.frame(site,out,temp,plot)
df.outlier.all <- rbind(df.outlier.all, df.outlier)

df.temperature$Temperature[out] <- NA

## D3
out <- which(df.temperature$Site == "D5" & df.temperature$Plot == "D3" & df.temperature$date.time == start.date.time$Start[6] 
             | df.temperature$Site == "D5" & df.temperature$Plot == "D3" & df.temperature$date.time > outlier.border)
temp <- df.temperature$Temperature[out]
site <- rep("D5",length(out))
plot <- rep("D3",length(out))   

df.outlier <- data.frame(site,out,temp,plot)
df.outlier.all <- rbind(df.outlier.all, df.outlier)

df.temperature$Temperature[out] <- NA

# Plot all data for Dublin 5 again

Dublin5.clean <- df.temperature %>% 
  filter(Site == "D5") %>% 
  droplevels()

dublin5.plot2 <- ggplot(Dublin5.clean, aes(x = date.time, y = Temperature)) +
  geom_line(aes(color = Plot)) +
  facet_wrap(~Plot) +
  labs(x = "Date", y = "Temperature (°C)", title = "Dublin 5 - Temperature per plot") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold")) +
  scale_x_datetime(limits = c(min.date,max.date)) +
  geom_hline(aes(yintercept = 30), linetype = "dashed") + 
  scale_y_continuous(limits = c(0,30))

dublin5.plot2

ggsave(paste0(mydir,"Dublin 5 - Temperature per plot - cleaned 2.png"),
       dublin5.plot2,
       width = 8, height = 6)


####################
##                ##
##   Limerick 2   ##
##                ##
####################

limerick2.plot1 <- ggplot(Limerick2, aes(x = date.time, y = Temperature)) +
  geom_line(aes(color = Plot)) +
  facet_wrap(~Plot) +
  labs(x = "Date", y = "Temperature (°C)", title = "Limerick 2 - Temperature per plot") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold")) +
  scale_x_datetime(limits = c(min.date,max.date)) +
  geom_hline(aes(yintercept = 30), linetype = "dashed")

limerick2.plot1

ggsave(paste0(mydir,"Limerick 2 - Temperature per plot - cleaned 1.png"),
       limerick2.plot1,
       width = 8, height = 6)

# Control 1 does not have any additional outliers
# All other plots seem to have some outliers starting in August, when the recording obviously failed
# All these outliers are approximatelly below 10 °C
# Add the column month to Limerick2 first

Limerick2$month <- as.numeric(substr(Limerick2$date.time, start = 6, stop = 7))

limerick2.below10 <- Limerick2 %>% 
  filter(Temperature < 10 & month >= 8) 

limerick2.below10

# Find the earliest date and time for each plot where the values below 10 occur

earliest.date <- limerick2.below10 %>% 
  group_by(Plot) %>% 
  summarize(Date = min(date.time))

earliest.date

# Exclude all outliers

## C1

out <- which(df.temperature$Site == "L2" & df.temperature$Plot == "C1" & df.temperature$date.time >= earliest.date$Date[1])
temp <- df.temperature$Temperature[out]
site <- rep("L2",length(out))
plot <- rep("C1",length(out))

df.outlier <- data.frame(site,out,temp,plot)
df.outlier.all <- rbind(df.outlier.all, df.outlier)

df.temperature$Temperature[out] <- NA

## C2

out <- which(df.temperature$Site == "L2" & df.temperature$Plot == "C2" & df.temperature$date.time >= earliest.date$Date[2])
temp <- df.temperature$Temperature[out]
site <- rep("L2",length(out))
plot <- rep("C2",length(out))

df.outlier <- data.frame(site,out,temp,plot)
df.outlier.all <- rbind(df.outlier.all, df.outlier)

df.temperature$Temperature[out] <- NA

## C3

out <- which(df.temperature$Site == "L2" & df.temperature$Plot == "C3" & df.temperature$date.time >= earliest.date$Date[3])
temp <- df.temperature$Temperature[out]
site <- rep("L2",length(out))
plot <- rep("C3",length(out))

df.outlier <- data.frame(site,out,temp,plot)
df.outlier.all <- rbind(df.outlier.all, df.outlier)

df.temperature$Temperature[out] <- NA

## D1

out <- which(df.temperature$Site == "L2" & df.temperature$Plot == "D1" & df.temperature$date.time >= earliest.date$Date[4])
temp <- df.temperature$Temperature[out]
site <- rep("L2",length(out))
plot <- rep("D1",length(out))

df.outlier <- data.frame(site,out,temp,plot)
df.outlier.all <- rbind(df.outlier.all, df.outlier)

df.temperature$Temperature[out] <- NA

## D2

out <- which(df.temperature$Site == "L2" & df.temperature$Plot == "D2" & df.temperature$date.time >= earliest.date$Date[5])
temp <- df.temperature$Temperature[out]
site <- rep("L2",length(out))
plot <- rep("D2",length(out))

df.outlier <- data.frame(site,out,temp,plot)
df.outlier.all <- rbind(df.outlier.all, df.outlier)

df.temperature$Temperature[out] <- NA

## D3

out <- which(df.temperature$Site == "L2" & df.temperature$Plot == "D3" & df.temperature$date.time >= earliest.date$Date[6])
temp <- df.temperature$Temperature[out]
site <- rep("L2",length(out))
plot <- rep("D3",length(out))

df.outlier <- data.frame(site,out,temp,plot)
df.outlier.all <- rbind(df.outlier.all, df.outlier)

df.temperature$Temperature[out] <- NA

# Plot all data for Limerick 2 again

Limerick2.clean <- df.temperature %>% 
  filter(Site == "L2") %>% 
  droplevels()

limerick2.plot2 <- ggplot(Limerick2.clean, aes(x = date.time, y = Temperature)) +
  geom_line(aes(color = Plot)) +
  facet_wrap(~Plot) +
  labs(x = "Date", y = "Temperature (°C)", title = "Limerick 2 - Temperature per plot") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold")) +
  scale_x_datetime(limits = c(min.date,max.date)) +
  geom_hline(aes(yintercept = 30), linetype = "dashed") + 
  scale_y_continuous(limits = c(0,30))

limerick2.plot2  ## looks much more reasonable!

ggsave(paste0(mydir,"Limerick 2 - Temperature per plot - cleaned 2.png"),
       limerick2.plot2,
       width = 8, height = 6)

#================================================
# Export outliers

df.outlier.all$date.time <- df.temperature$date.time[df.outlier.all$out]

write.table(df.outlier.all, "All temperature outliers - visual inspection of each plot.csv", sep = ",", row.names = F)

#================================================
# Save new data

save(df.temperature, file = paste0(data.dir,"Temperature data - outliers removed.RData"))

#================================================
# Plot cleaned data

## Lines

p7 <- ggplot(df.temperature, aes(x = date.time, y = Temperature)) +
  geom_line(aes(color = treatment, linetype = Plot)) +
  facet_grid(Region ~ Farm, scales = "free_y") +
  labs(x = "Date", y = "Temperature (°C)", title = "Temperature per site and plot") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold")) +
  scale_x_datetime(limits = c(min.date,max.date)) + 
  scale_y_continuous(limits = c(0,30))
  
p7

ggsave(paste0(mydir,"Temperature per site and plot - cleaned final.png"),
       p7,
       width = 13, height = 10)

## Boxplots

p8 <- ggplot(df.temperature, aes(x = treatment, y = Temperature)) +
  geom_boxplot(aes(fill = treatment),
               outlier.shape = 21) +
  facet_grid(Region ~ Farm, scales = "free_y") +
  labs(x = "Date", y = "Temperature (°C)", title = "Temperature per site and plot") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold")) + 
  scale_y_continuous(limits = c(0,30))

p8

ggsave(paste0(mydir,"Temperature per site and plot - boxplots - cleaned final.png"),
       p8,
       width = 13, height = 10)

## Scatterplot of raw data with jitter

p9 <- ggplot(df.temperature, aes(x = treatment, y = Temperature, fill = treatment, color = treatment)) +
  geom_point(size = 2, shape = 21, alpha = 0.5,
             position = position_jitterdodge()) +
  facet_grid(Region ~ Farm, scales = "free_y") +
  labs(x = "Date", y = "Temperature (°C)", title = "Temperature per site and plot") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold")) + 
  scale_y_continuous(limits = c(0,30))

p9

ggsave(paste0(mydir,"Temperature per site and plot - jitter - cleaned final.png"),
       p9,
       width = 13, height = 10)


#================================================
# Add events (start drought, end drought, day 64)

## THIS PART IS NOT FINISHED YET!!!

##------------------##
##                  ##
##     for-loop     ##
##                  ##
##------------------##

events <- read.csv("ExperimentDatesMI.csv", sep = ",", header = T)

regions <- unique(df.temperature$Region)
sites <- unique(df.temperature$Site)

for (site in sites){
  mysite <- df.temperature %>% 
    filter(Site == site) %>% 
    droplevels()
  
  title <- paste(mysite$Region[1],mysite$No.Site[1])
  
  start.drought <- as.POSIXct(mdy(events$StartDrought[events$Site == site]), tz = "GMT") + 11*60*60
  end.drought <- as.POSIXct(mdy(events$EndDrought[events$Site == site]), tz = "GMT") + 11*60*60
  day64 <- as.POSIXct(mdy(events$Day64[events$Site == site]), tz = "GMT") + 11*60*60
  
  g <- ggplot(mysite,
              aes(x = date.time, y = Temperature)) +
    geom_line(aes(color = Plot)) +
    facet_wrap(~ Plot) +
    labs(x = "Date", y = "Temperature (°C)", title = title) +
    mytheme +
    scale_y_continuous(limits = c(-1,41), breaks = c(0,10,20,30,40)) +
    scale_x_datetime(limits = c(min.date,max.date), 
                     breaks = c(start.drought,
                                end.drought,
                                day64),
                     date_labels = "%b %d") +
    geom_vline(aes(xintercept = start.drought,
                   linetype = Plot), color = "gray50") +
    geom_vline(aes(xintercept = end.drought,
                   linetype = Plot), color = "gray50") +
    geom_vline(aes(xintercept = day64),
               linetype = "dashed", color = "gray50") +
    scale_linetype_manual(values = c(2,2,2,2,2,2))
    #scale_linetype_manual(values = c(0,0,0,2,2,2))
  
  ggsave(paste0(mydir,site," - Soil moisture with events.png"), g,
         width = 8, height = 6)
  
}
