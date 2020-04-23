####################################################################
## Calculating Bray-Curtis Dissimilarity using Invertebrates data ##
## from the Grassland Resilience Experiment in 2017               ##
## Data provided by Matthew Magilton                              ##
##                                                                ##
## Idea: for each Region and Farm (i.e. each site), calculate an  ##
## average community composition for all control plots. Use this  ##
## as the reference community and compare it to each plot (both   ##
## control and drought) using the Bray-Curtis Dissimilarity       ##
## If BC = 0, the compared communities are equal                  ##
## If BC = 1, the compared communities do not share any species   ##
##                                                                ##
## Author of the script:                                          ##
## Maja Ilic (QUB)                                                ##
##                                                                ##
## First modified: 10 Feb 2020                                    ##
## Last modified: 23 Apr 2020                                     ##
####################################################################

#===========================================
# Clear objects from the workspace ####

rm(list = ls())

#===========================================
# Set working directory ####

# Maja's desktop PC

setwd("C:/Users/3054311/Documents/My Documents/Grassland project/07_Invertebrates")

data.dir <- paste0(getwd(),"/Data/")
figures.dir <- paste0(getwd(),"/Figures/")
script.dir <- paste0(getwd(),"/Script/")

#======================================
# Packages ####

library(ggplot2)      ## for plots
library(dplyr)        ## for easier data management
library(vegan)        ## for Bray-Curtis

#======================================
# Import data ####

df <- read.csv(paste0(data.dir,"Invertebrates.csv"), sep = ",", header = T)  ## read in the dataset

# Add column Region

df[which(df$loc == "bo"),"Region"] <- "Border"
df[which(df$loc == "co"),"Region"] <- "Cork"
df[which(df$loc == "du"),"Region"] <- "Dublin"
df[which(df$loc == "li"),"Region"] <- "Limerick"

# Change "site" to "Farm"

names(df)[which(names(df) == "site")] <- "Farm"

# Change "treat" to "treatment" 

names(df)[which(names(df) == "treat")] <- "treatment"

# Change "c" to "Control" and "d" to "Drought"

df$treatment <- gsub("c", "Control", 
                     as.character(df$treatment))

df$treatment <- gsub("d", "Drought", 
                     as.character(df$treatment))

# Add column "Plot"

df$Plot <- paste0(substr(df$treatment,1,1),df$rep)

# Add column site_ID

df$site_ID <- paste0(df$Region, df$Farm)

# Export the raw data with correct names of the columns (exclude the column "moisture")

write.table(df[,c(40,2,42,3,41,5:9,11:39)], 
            paste0(data.dir,"Grassland Resilience 2017 - Invertebrates raw data.csv"), 
            sep = ",", row.names = F)

# Extract only the relevant columns and put them in a logical order 

df.new <- df[,c(40,2,42,3,41,5,13:39)]

#======================================
#### Prior to calculations ####

# Prior to the Bray-Curtis calculation, find out which plots at which sites have missing data
# (no sample available) or for which no invertebrates were found (all counts are zero)
# For the later, it has to be differentiated between "real soil species only" and "all speices"

df.new$Total.Soil <- rowSums(df.new[,c(7:17)], na.rm = F)  
df.new$Total.All <- rowSums(df.new[,c(7:33)], na.rm = F) 

soil.zero <- df.new %>% filter(Total.Soil == 0)  
all.zero <- df.new %>% filter(Total.All == 0)

# soil.zero and all.zero are exactly the same

# Finding lines with all NAs is easy, as given the structure of the dataset
# as for each unavailable sample, NAs were inserted for every order
# Therefore, it is enough to check one single order for NAs:

all.NA <- df.new %>% filter(is.na(earth))

# Combine all.zero and all.NA datasets together and export them 

zero.or.NA <- rbind(all.zero,all.NA) 
zero.or.NA$CODE <- paste(zero.or.NA$site_ID, zero.or.NA$Plot, zero.or.NA$day, sep = ".")
unique(zero.or.NA$CODE)

# There are in total 9 samples with all zero data and 43 samples that were not available

#======================================
# Calculate Bray-Curtis Dissimilarity using a for-loop ####
# The for-loop contains 2 nested for-loops:
# Outer for-loop loops through all sites
# Inner for-loop loops trough all days
# Because two running indeces are used (i for sites and j for days), 
# I included the running index "m" that gets increased by 1 after every inner loop process is done

# Save all sites and all days as single vector to be used later in the for-loop

sites <- unique(df.new$site_ID)
days <- unique(df.new$day)

# Include running index "m", set it to 1 

m <- 1

for (i in 1:length(sites)){   ## outer loop for sites
  
  for (j in 1:length(days)){  ## inner loop for days
  
    # Extract data for one single site and one single day
    
    mysite <- df.new %>% 
      filter(site_ID == sites[i] & day == days[j])
    
    # Add row names to your data frame, useful later for the Bray-Curtis calculation
    
    rownames(mysite) <- paste(mysite$site_ID, mysite$Plot, mysite$day, sep = ".")
    
    # Extract control plots
    
    control_plots <- mysite %>% 
      filter(treatment == "Control")
  
    # Calculate average community composition (average abundance of each group/column)
    # Important: only columns with species are included 
    # ref.comm.all => includes all species
    # ref.comm.soil => includes only the "real soil species"
    
    ref.comm.all <- control_plots[,c(7:33)] %>% summarise_each(funs( mean( .,na.rm = TRUE)))
    rownames(ref.comm.all) <- "REF.COMM.ALL"
    
    ref.comm.soil <- control_plots[,c(7:17)] %>% summarise_each(funs( mean( .,na.rm = TRUE)))
    rownames(ref.comm.soil) <- "REF.COMM.SOIL"
    
    # Add the reference community to your site
    
    z.all <- rbind(mysite[,c(7:33)], ref.comm.all)
    z.soil <- rbind(mysite[,c(7:17)], ref.comm.soil)
    
    # Calculate Bray-Curtis Dissimilary using the function vegdist
    # Here, there are different options for the actual calculation
    # binary = "TRUE": performs presence/absence standardization before analysis using decostand.
    # I only did the calculations without
    
    # Matrix of all Bray-Curtis Dissimilarities
    
    calc.BC.all <- as.matrix(vegdist(z.all, method = "bray", binary = F, na.rm = T))
    calc.BC.soil <- as.matrix(vegdist(z.soil, method = "bray", binary = F, na.rm = T))
    
    # Extract just the last column containing the pairwise comparison of each plot 
    # to the average control community composition
    
    BC.mean.C.all <- calc.BC.all[1:8,"REF.COMM.ALL"]
    BC.mean.C.soil <- calc.BC.soil[1:8,"REF.COMM.SOIL"]
    
    # Add the BCs to the mysite dataset
    
    mysite$BC_meanC_all <- BC.mean.C.all
    mysite$BC_meanC_soil <- BC.mean.C.soil
    
    # This is optional: add the names from the BC vectors to mysite 
    # to be sure you extracted the correct comparions of single communities to REF.COMM.
    
    mysite$community_all <- names(BC.mean.C.all)
    mysite$community_soil <- names(BC.mean.C.soil)
    
    # Add all chunks of the dataset together (combine all sites again)
    
    if (m == 1){
      df.sites <- mysite
    }
    
    if (m > 1){
      df.sites <- rbind(df.sites,mysite)
    }
    
    # Increase m (the running index) by 1
    
    m <- m + 1
  }
}

#======================================
# Export data ####

write.table(df.sites, paste0(data.dir,"Invertebrates with Bray Curtis.csv"), sep = ",", row.names = F)

#======================================
# Plot the results ####

# All species included
# Panels: Region ~ Farm

g1 <- ggplot(df.sites, aes(x = day, y = BC_meanC_all)) +
  geom_point(aes(color = treatment, fill = treatment),
             size = 3, shape = 21, alpha = 0.5) +
  geom_smooth(method = "lm", se = F, aes(color = treatment)) +
  facet_grid(Region ~ Farm) +
  theme_minimal() +
  theme(panel.spacing.y = unit(0.5, "lines"),
        strip.background = element_rect(),
        axis.title.x = element_text(size = 12, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 12, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(face = "bold")) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,65), breaks = c(0,32,64)) +
  labs(x = "Days since the end of the drought",
       y = "Bray-Curtis Dissimilariy Index",
       title = "All species included",
       subtitle = "0: communities are equal;   1: communities do not share any species")

ggsave(paste0(figures.dir,"Invertebrates Bray Curtis all species V1.png"), g1, width = 10, height = 7)

# Panels: ~ Region

g2 <- ggplot(df.sites, aes(x = day, y = BC_meanC_all)) +
  geom_point(aes(color = treatment, fill = treatment),
             size = 3, shape = 21, alpha = 0.5) +
  geom_smooth(method = "lm", se = F, aes(color = treatment)) +
  facet_wrap(~ Region, ncol = 4) +
  theme_minimal() +
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_rect(),
        axis.title.x = element_text(size = 12, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 12, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(face = "bold")) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,65), breaks = c(0,32,64)) +
  labs(x = "Days since the end of the drought",
       y = "Bray-Curtis Dissimilariy Index",
       title = "All species included",
       subtitle = "0: communities are equal;   1: communities do not share any species")

ggsave(paste0(figures.dir,"Invertebrates Bray Curtis all species V2.png"), g2, width = 10, height = 6)

# "Real soil species"
# Panels: Region ~ Farm

g3 <- ggplot(df.sites, aes(x = day, y = BC_meanC_soil)) +
  geom_point(aes(color = treatment, fill = treatment),
             size = 3, shape = 21, alpha = 0.5) +
  geom_smooth(method = "lm", se = F, aes(color = treatment)) +
  facet_grid(Region ~ Farm) +
  theme_minimal() +
  theme(panel.spacing.y = unit(0.5, "lines"),
        strip.background = element_rect(),
        axis.title.x = element_text(size = 12, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 12, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(face = "bold")) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,65), breaks = c(0,32,64)) +
  labs(x = "Days since the end of the drought",
       y = "Bray-Curtis Dissimilariy Index",
       title = "Only \"real soil species\" included",
       subtitle = "0: communities are equal;   1: communities do not share any species")

ggsave(paste0(figures.dir,"Invertebrates Bray Curtis soil species V1.png"), g3, width = 10, height = 7)

# Panels: ~ Region

g4 <- ggplot(df.sites, aes(x = day, y = BC_meanC_soil)) +
  geom_point(aes(color = treatment, fill = treatment),
             size = 3, shape = 21, alpha = 0.5) +
  geom_smooth(method = "lm", se = F, aes(color = treatment)) +
  facet_wrap(~ Region, ncol = 4) +
  theme_minimal() +
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_rect(),
        axis.title.x = element_text(size = 12, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 12, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(face = "bold")) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,65), breaks = c(0,32,64)) +
  labs(x = "Days since the end of the drought",
       y = "Bray-Curtis Dissimilariy Index",
       title = "Only \"real soil species\" included",
       subtitle = "0: communities are equal;   1: communities do not share any species")

ggsave(paste0(figures.dir,"Invertebrates Bray Curtis soil species V2.png"), g4, width = 10, height = 6)

# Both BC indices combined (only lines)

# Panels: Region ~ Farm

g5 <- ggplot(df.sites) +
  geom_smooth(aes(x = day, y = BC_meanC_all, color = treatment),
              method = "lm", se = F, size = 1) +
  geom_smooth(aes(x = day, y = BC_meanC_soil, color = treatment),
              method = "lm", se = F, size = 0.8, linetype = "dashed") +
  facet_grid(Region ~ Farm) +
  theme_minimal() +
  theme(panel.spacing.y = unit(0.5, "lines"),
        strip.background = element_rect(),
        axis.title.x = element_text(size = 12, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 12, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(face = "bold")) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,65), breaks = c(0,32,64)) +
  labs(x = "Days since the end of the drought",
       y = "Bray-Curtis Dissimilariy Index",
       title = "Solid: all species;   dashed: only real soil species",
       subtitle = "0: communities are equal;   1: communities do not share any species")

ggsave(paste0(figures.dir,"Invertebrates Bray Curtis V1.png"), g5, width = 10, height = 7)

# Panels: ~ Region

g6 <- ggplot(df.sites) +
  geom_smooth(aes(x = day, y = BC_meanC_all, color = treatment),
              method = "lm", se = F, size = 1) +
  geom_smooth(aes(x = day, y = BC_meanC_soil, color = treatment),
              method = "lm", se = F, size = 0.8, linetype = "dashed") +
  facet_wrap(~ Region, ncol = 4) +
  theme_minimal() +
  theme(panel.spacing.y = unit(0.5, "lines"),
        strip.background = element_rect(),
        axis.title.x = element_text(size = 12, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 12, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(face = "bold")) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,65), breaks = c(0,32,64)) +
  labs(x = "Days since the end of the drought",
       y = "Bray-Curtis Dissimilariy Index",
       title = "Solid: all species;   dashed: only real soil species",
       subtitle = "0: communities are equal;   1: communities do not share any species")

ggsave(paste0(figures.dir,"Invertebrates Bray Curtis V2.png"), g6, width = 10, height = 6)

# Heatmap
# All species included

g9 <- ggplot(df.sites, aes(x = as.factor(day), y = Plot, fill = BC_meanC_all)) +
  geom_tile() +
  facet_wrap(~ site_ID, scales = "free_x") +
  theme_minimal() +
  theme(panel.spacing.y = unit(0.5, "lines"),
        strip.background = element_rect(),
        axis.title.x = element_text(size = 12, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 12, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(face = "bold")) +
  scale_fill_viridis_c(name = "Bray-Curtis") +
  geom_hline(yintercept = 4.5, color = "white", size = 1) +
  labs(x = "Days since the end of the drought",
       y = "Bray-Curtis Dissimilariy Index",
       title = "All species included",
       subtitle = "0: communities are equal;   1: communities do not share any species")

ggsave(paste0(figures.dir,"Invertebrates Bray Curtis Heatmap.png"), g9, width = 10, height = 8)

#======================================
# Plot soil properties

df.control <- df %>% filter(treatment == "Control")
df.drought <- df %>% filter(treatment == "Drought")

# Soil organic matter

# Control plots 

om.C <- ggplot(df.control, aes(x = as.factor(day), y = soil_om, fill = Region, color = Region)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_point(size = 2, shape = 21, alpha = 0.7) +
  facet_grid(Region ~ Farm) + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 12, color = "black", face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 12, color = "black", face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_y_continuous(limits = c(0,80)) +
  labs(x = "Day since the end of the simulated drought", y = "Soil organic matter (%)",
       title = "Soil organic matter - Control plots")

ggsave("Soil organic matter - Control plots.png", om.C, width = 10, height = 7.5)

om.C2 <- om.C + scale_y_continuous(limits = c(0,45)) +
  labs(title = "Soil organic matter without outlier in Border 2, Day 64 - Control plots")
om.C2

ggsave("Soil organic matter wo outlier - Control plots.png", om.C2, width = 10, height = 7.5)

# Drought plots 

om.D <- ggplot(df.drought, aes(x = as.factor(day), y = soil_om, fill = Region, color = Region)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_point(size = 2, shape = 21, alpha = 0.7) +
  facet_grid(Region ~ Farm) + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 12, color = "black", face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 12, color = "black", face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_y_continuous(limits = c(0,45)) +
  labs(x = "Day since the end of the simulated drought", y = "Soil organic matter (%)",
       title = "Soil organic matter - Drought plots")

ggsave("Soil organic matter - Drought plots.png", om.D, width = 10, height = 7.5)

# Control and Drought plots in one figure

om.plot <- ggplot(df, aes(x = as.factor(day), y = soil_om, fill = Region, color = Region)) +
  geom_boxplot(outlier.shape = NA, aes(alpha = treatment)) +
  geom_point(size = 2, shape = 21, aes(alpha = treatment), position = position_dodge(0.75)) +
  facet_grid(Region ~ Farm) + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 12, color = "black", face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 12, color = "black", face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_alpha_manual(values = c(0.2,0.6), name = "Treatment") +
  scale_y_continuous(limits = c(0,45)) +
  labs(x = "Day since the end of the simulated drought", y = "Soil organic matter (%)",
       title = "Soil organic matter - Control & Drought plots",
       subtitle = "Outlier in Control plot in Border 2, Day 64, was excluded from plotting")

ggsave("Soil organic matter - All plots.png", om.plot, width = 10, height = 7.75)

# Soil pH

# Control plots

pH.C <- ggplot(df.control, aes(x = as.factor(day), y = soil_ph, fill = Region, color = Region)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_point(size = 2, shape = 21, alpha = 0.7) +
  facet_grid(Region ~ Farm) + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 12, color = "black", face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 12, color = "black", face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_y_continuous(limits = c(4.9,9.1)) +
  labs(x = "Day since the end of the simulated drought", y = "Soil pH",
       title = "Soil pH - Control plots") +
  geom_hline(aes(yintercept = 7), linetype = "dashed")

ggsave("Soil pH - Control plots.png", pH.C, width = 10, height = 7.5)

# Drought plots

pH.D <- ggplot(df.drought, aes(x = as.factor(day), y = soil_ph, fill = Region, color = Region)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_point(size = 2, shape = 21, alpha = 0.7) +
  facet_grid(Region ~ Farm) + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 12, color = "black", face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 12, color = "black", face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_y_continuous(limits = c(4.9,9.1)) +
  labs(x = "Day since the end of the simulated drought", y = "Soil pH",
       title = "Soil pH - Drought plots") +
  geom_hline(aes(yintercept = 7), linetype = "dashed")

ggsave("Soil pH - Drought plots.png", pH.D, width = 10, height = 7.5)

# Control and Drought plots in one figure

pH.plot <- ggplot(df, aes(x = as.factor(day), y = soil_ph, fill = Region, color = Region)) +
  geom_boxplot(outlier.shape = NA, aes(alpha = treatment)) +
  geom_point(size = 2, shape = 21, aes(alpha = treatment), position = position_dodge(0.75)) +
  facet_grid(Region ~ Farm) + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 12, color = "black", face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 12, color = "black", face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_alpha_manual(values = c(0.2,0.6), name = "Treatment") +
  scale_y_continuous(limits = c(4.9,9.1)) +
  labs(x = "Day since the end of the simulated drought", y = "Soil pH",
       title = "Soil pH - Control & Drought plots") +
  geom_hline(aes(yintercept = 7), linetype = "dashed")

ggsave("Soil pH - All plots.png", pH.plot, width = 10, height = 7.5)

#======================================
# Correlate soil organic matter with soil moisture

# Load soil moisture data

load(file = "C:/Users/3054311/Documents/My Documents/Grassland project/05_Resilience/Data/Resilience_Soil_with_Rainfall_non_rm.RData")
resilience_soil_non_rm <- resilience_soil_non_rm[!is.na(resilience_soil_non_rm$plot_meanC_ratio),]
resilience_soil_non_rm <- resilience_soil_non_rm %>% filter(day == 0 | day == 32 | day == 64)

df$Farm <- as.factor(df$Farm)

df.soil <- left_join(resilience_soil_non_rm[,c("Region","Farm","Plot","treatment","day","Avg.Moisture")],
                     df[,c("Region","Farm","Plot","treatment","day","soil_om")],
                     by = c("Region","Farm","Plot","treatment","day"))

g.soil <- ggplot(df.soil,aes(x = Avg.Moisture, y = soil_om)) +
  geom_point(stat = "identity", size = 3, shape = 21, alpha = 0.5,
             aes(fill = treatment, color = treatment)) +
  facet_grid(Region ~ Farm) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 12, color = "black", face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 12, color = "black", face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_fill_manual(values = c("red","blue"), name = "Treatment") +
  scale_color_manual(values = c("red","blue"), name = "Treatment") +
  labs(x = "Average soil moisture (%)", y = "Soil organic matter (%)") +
  geom_smooth(method = "lm", se = F, color = "black", size = 0.8) +
  geom_smooth(method = "lm", se = F, aes(color = treatment), size = 0.6)

ggsave("Soil moisture vs soil organic matter.png", g.soil, width = 10, height = 7.5) 
