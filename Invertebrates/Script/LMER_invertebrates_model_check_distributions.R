###################################################################################
## This script fits linear models to the invertebrate Bray-Curtis Dissimilarity  ##
## in the Grassland Resilience experiment with fixed factors "day" (days after   ##
## the end of the drought), effect size of soil moisture (calculated as soil     ##
## moisture in each plot C1-D3 minus the average soil moisture in control plots),##
## treatment (control vs. drought) and Region (Border, Cork, Dublin, Limerick)   ##
## Random effects: Farm (1-5); and Farm nested in Region                         ##
## The response is the invertebrate BC Index (each plot compared to an average   ##
## invertebrate community composition in control plots at each site and each     ##
## sampling day (9, 32, 64), i.e. reference community)                           ##
##                                                                               ##
## The script has be modified in order to test for the best moving window width  ##
## used to calculate the effect size of the soil moisture                        ##
## It uses previously written function that calculates the soil moisture ES      ##
## for any given timepoint of the sampling day (e.g. 12:00 noon) and any given   ##
## window width (e.g. 24 h, 48 h etc., negative if in the past)                  ##
## The aim is to find the moving window width that explains most of the variance ##
## and to find the best model by comparing all models with each other            ##
##                                                                               ##
## Additionally, the distribution of the residuals and response variable is      ##
## checked using check_distribution from the performace and see package          ##
##                                                                               ##
## Author of the modified script:                                                ##
## Maja Ilic M.Ilic@qub.ac.uk                                                    ##
## first modified: 25 Feb 2020                                                   ##
## last modified: 25 Feb 2020                                                    ##
###################################################################################

#===========================================
# Clear objects from the workspace ----

rm(list = ls())

#===========================================
# Set working directory ----

# Maja's desktop PC

setwd("C:/Users/3054311/Documents/My Documents/Grassland project/07_Invertebrates")

data.dir <- paste0(getwd(),"/Data/")
figures.dir <- paste0(getwd(),"/Figures/")
script.dir <- paste0(getwd(),"/Script/")

#===========================================
#Packages ----

library(ggplot2)
library(lubridate)
library(dplyr)
library(tidyr)
library(stringr)
library(lme4)
library(rcompanion)
library(car)
library(performance)
library(insight)
library(see)
library(cowplot)
library(patchwork)
library(gridExtra)
library(effects)
library(afex)
library(RColorBrewer)

#===========================================
# Get raw data for soil moisture ----

raw.data.dir <- "C:/Users/3054311/Documents/My Documents/Grassland project/01_Probes_soil_moisture_and_temp/Data/"

load(paste0(raw.data.dir,"Cleaned_Soil_Moisture_Data_from_Loggers.RData"))

# Events ----

load(file = paste0(raw.data.dir,"All events.RData"))

#================================================
# Set directory for data and plots ----

# dir.create(paste0(data.dir,"LMER model invertebrates check distributions"))
# dir.create(paste0(figures.dir,"LMER model invertebrates check distributions"))

mydir.data <- paste0(data.dir,"LMER model invertebrates check distributions/")
mydir <- paste0(figures.dir,"LMER model invertebrates check distributions/")

#================================================
# Load function mov.win.diff ----

source(paste0(script.dir,"Function mov.win.diff.R"))

#================================================
# Import invertebrate data: Bray-Curtis Dissimilarity ----

inverts_data <- read.csv(paste0(data.dir,"Invertebrates with Bray Curtis.csv"), sep = ",", header = T)

# Extract only relevant columns

df_inverts <- inverts_data[,c(1:6,36:37)]

# Change the column "BC_meanC_all" to "inverts_BC_all"

names(df_inverts)[which(names(df_inverts) == "BC_meanC_all")] <- "inverts_BC_all"

# Change the column "BC_meanC_soil" to "inverts_BC_soil"

names(df_inverts)[which(names(df_inverts) == "BC_meanC_soil")] <- "inverts_BC_soil"

# Remove plots C4 and D4

df_inverts <- df_inverts %>% filter(Plot != "C4" & Plot != "D4")

# Remove Border 2 and 4

df_inverts <- df_inverts %>% filter(site_ID != "Border2" & site_ID != "Border4")

#================================================
# Run a for-loop for soil moisture data ----

timepoint <- 12
width <- seq(-192,0,1)

m <- 1

resp.variable <- "inverts_BC_all"

for (daytime.hr in timepoint) {
  
  for (duration.hr in width) {
    
    # Extract soil moisture data ----
    # Use the mov.win.diff function to extract relevant soil moisture data for the given time inverval
    # and calculate the effect size (for the days 0, 32 and 64)
    
    df_soil <- mov.win.diff(data = df.moisture,
                            dates = events,
                            daytime.hr = daytime.hr,
                            duration.hr = duration.hr,
                            doplot = FALSE,
                            list.figures = NULL)
    
    #================================================
    # Combine effect size (ES) for soil moisture and BC for invertebrates ----
    
    names(df_soil)[names(df_soil) == "plot_meanC_ES"] <- "soil_ES"
    
    df_joined <- full_join(df_inverts[, c("Region", "Farm", "Plot", 
                                          "site_ID", "day", 
                                          "treatment", "inverts_BC_all")], 
                           df_soil[, c("Region", "Farm", "Plot", 
                                       "site_ID", "day", 
                                       "treatment", "soil_ES")], 
                           by = c("Region", "Farm", "Plot", 
                                  "site_ID", "day", 
                                  "treatment"))
    
    # Exclude all rows with NAs (columns BC_meanC_all and soil_ES) 
    
    df_joined <- df_joined[which(!is.na(df_joined$inverts_BC_all)),]
    df_joined <- df_joined[which(!is.na(df_joined$soil_ES)),]
    
    # Change region, farm and treatment factors
    
    df_joined$Region <- as.factor(df_joined$Region)
    df_joined$Farm <- as.factor(df_joined$Farm)
    df_joined$treatment <- as.factor(df_joined$treatment)
    
    #================================================
    
    ##################
    ##              ##
    ##  Fit models  ## ----
    ##              ##
    ##################
    
    lmer_day <- lmer(inverts_BC_all ~ day + (1|Region/Farm),
                     data = df_joined)
    
    lmer_soil <- lmer(inverts_BC_all ~ soil_ES + (1|Region/Farm),
                      data = df_joined)
    
    lmer_trt <- lmer(inverts_BC_all ~ treatment + (1|Region/Farm),
                     data = df_joined)
    
    lmer_region <- lmer(inverts_BC_all ~ Region + (1|Farm),
                        data = df_joined)
    
    lmer_day_soil <- lmer(inverts_BC_all ~ day*soil_ES + (1|Region/Farm),
                          data = df_joined)
    
    lmer_day_trt <- lmer(inverts_BC_all ~ day*treatment + (1|Region/Farm),
                         data = df_joined)
    
    lmer_day_region <- lmer(inverts_BC_all ~ day*Region + (1|Farm),
                            data = df_joined)
    
    lmer_soil_trt <- lmer(inverts_BC_all ~ soil_ES*treatment + (1|Region/Farm),
                          data = df_joined)
    
    lmer_soil_region <- lmer(inverts_BC_all ~ soil_ES*Region + (1|Farm),
                             data = df_joined)
    
    lmer_trt_region <- lmer(inverts_BC_all ~ treatment*Region + (1|Farm),
                            data = df_joined)
    
    lmer_day_soil_trt <- lmer(inverts_BC_all ~ day*soil_ES*treatment + (1|Region/Farm),
                              data = df_joined)
    
    lmer_day_soil_region <- lmer(inverts_BC_all ~ day*soil_ES*Region + (1|Farm),
                                 data = df_joined)
    
    lmer_day_trt_region <- lmer(inverts_BC_all ~ day*treatment*Region + (1|Farm),
                                data = df_joined)
    
    lmer_soil_trt_region <- lmer(inverts_BC_all ~ soil_ES*treatment*Region + (1|Farm),
                                 data = df_joined)
    
    lmer_day_soil_trt_region <- lmer(inverts_BC_all ~ day*soil_ES*treatment*Region + (1|Farm),
                                     data = df_joined)
    
    #================================================
    # Check the distribution of residuals and response in all models ----
    
    result_day <- as.data.frame(check_distribution(lmer_day))
    result_soil <- as.data.frame(check_distribution(lmer_soil))
    result_trt <- as.data.frame(check_distribution(lmer_trt))
    result_region <- as.data.frame(check_distribution(lmer_region))
    result_day_soil <- as.data.frame(check_distribution(lmer_day_soil))
    result_day_trt <- as.data.frame(check_distribution(lmer_day_trt))
    result_day_region <- as.data.frame(check_distribution(lmer_day_region))
    result_soil_trt <- as.data.frame(check_distribution(lmer_soil_trt))
    result_soil_region <- as.data.frame(check_distribution(lmer_soil_region))
    result_trt_region <- as.data.frame(check_distribution(lmer_trt_region))
    result_day_soil_trt <- as.data.frame(check_distribution(lmer_day_soil_trt))
    result_day_soil_region <- as.data.frame(check_distribution(lmer_day_soil_region))
    result_day_trt_region <- as.data.frame(check_distribution(lmer_day_trt_region))
    result_soil_trt_region <- as.data.frame(check_distribution(lmer_soil_trt_region))
    result_day_soil_trt_region <- as.data.frame(check_distribution(lmer_day_soil_trt_region))
    
    result_day$Model <- "lmer_day"
    result_soil$Model <- "lmer_soil"
    result_trt$Model <- "lmer_trt"
    result_region$Model <- "lmer_region"
    result_day_soil$Model <- "lmer_day_soil"
    result_day_trt$Model <- "lmer_day_trt"
    result_day_region$Model <- "lmer_day_region"
    result_soil_trt$Model <- "lmer_soil_trt"
    result_soil_region$Model <- "lmer_soil_region"
    result_trt_region$Model <- "lmer_trt_region"
    result_day_soil_trt$Model <- "lmer_day_soil_trt"
    result_day_soil_region$Model <- "lmer_day_soil_region"
    result_day_trt_region$Model <- "lmer_day_trt_region"
    result_soil_trt_region$Model <- "lmer_soil_trt_region"
    result_day_soil_trt_region$Model <- "lmer_day_soil_trt_region"
    
    result_all_mod <- rbind(result_day,
                            result_soil,
                            result_trt,
                            result_region,
                            result_day_soil,
                            result_day_trt,
                            result_day_region,
                            result_soil_trt,
                            result_soil_region,
                            result_trt_region,
                            result_day_soil_trt,
                            result_day_soil_region,
                            result_day_trt_region,
                            result_soil_trt_region,
                            result_day_soil_trt_region)
    
    result_all_mod$Daytime <- daytime.hr
    result_all_mod$Duration <- duration.hr
    result_all_mod$Trial <- m
    
    if (m == 1){
      result_final <- result_all_mod
    }
    
    if (m > 1){
      result_final <- rbind(result_final,result_all_mod)
    }
    
    #================================================
    # Plot results
    
    result_all_mod_long <- gather(result_all_mod, key = "Group", value = "Probability",
                                  -Distribution, -Model, -Daytime, -Duration, -Trial)
    
    result_all_mod_long$Group[which(result_all_mod_long$Group == "p_Residuals")] <- "Residuals"
    result_all_mod_long$Group[which(result_all_mod_long$Group == "p_Response")] <- "Response"
    
    result_all_mod_long$Model <- as.factor(result_all_mod_long$Model)
    result_all_mod_long$Model <- factor(result_all_mod_long$Model, levels = unique(result_all_mod_long$Model))
    
    g <- ggplot(result_all_mod_long, aes(x = Distribution, y = Probability*100)) +
      geom_point(size = 2.5, aes(color = Group), position = position_dodge(0.9)) +
      geom_linerange(aes(x = Distribution,
                         ymin = 0, ymax = Probability*100,
                         color = Group),
                   size = 1, position = position_dodge(0.9)) +
      facet_wrap(~ Model) +
      theme_minimal() +
      theme(strip.background = element_rect(),
            plot.title = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 12, face = "bold")) +
      scale_color_brewer(palette = "Set1") +
      coord_flip() +
      labs(x = "\nDistribution",
           y = "Probability (%)\n",
           title = "LMER - Invertebrates Bray-Curtis Dissimilarity",
           subtitle = paste("Daytime: 12:00 h\nWindow width:",duration.hr,"h"))
    
    ggsave(paste0(mydir,"LMER Inverts BC check distributions - Window width ",duration.hr," h.png"),
           g, width = 35, height = 30, units = "cm")
    
    #================================================
    # Increase m ----
    
    m <- m + 1
    
    #================================================
    # Remove all models ----
    
    rm(lmer_day,
       lmer_soil,
       lmer_trt,
       lmer_region,
       lmer_day_soil,
       lmer_day_trt,
       lmer_day_region,
       lmer_soil_trt,
       lmer_soil_region,
       lmer_trt_region,
       lmer_day_soil_trt,
       lmer_day_soil_region,
       lmer_day_trt_region,
       lmer_soil_trt_region,
       lmer_day_soil_trt_region)
  }
}

#================================================
# Export results ----

write.table(result_final, paste0(mydir.data, "LMER Inverts BC check distributions.csv"), 
            sep = ",", row.names = F)
    
#================================================
# Plot results ----
# Remove all distributions that are zero for all models and all window widths

result_final_long <- gather(result_final, key = "Group", value = "Probability",
                            -Distribution, -Model, -Daytime, -Duration, -Trial)

result_final_long$Group[which(result_final_long$Group == "p_Residuals")] <- "Residuals"
result_final_long$Group[which(result_final_long$Group == "p_Response")] <- "Response"

sums <- result_final_long %>% 
  group_by(Distribution) %>% 
  summarize(Sum = sum(Probability))

non_zero <- sums$Distribution[which(sums$Sum > 0)]

result_final_relevant <- result_final_long %>% 
  filter(Distribution %in% non_zero)

result_mean <- result_final_relevant %>% 
  group_by(Model,Distribution,Group) %>% 
  summarize(Mean = mean(Probability,na.rm = T))

summary.g <- ggplot(result_mean, aes(x = Distribution, y = Mean*100)) +
  geom_point(size = 3, aes(color = Group, fill = Group), 
             position = position_dodge(0.9)) +
  geom_linerange(aes(x = Distribution,
                     ymin = 0, ymax = Mean*100,
                     color = Group),
                 size = 1, position = position_dodge(0.9)) +
  facet_wrap(~ Model) +
  theme_minimal() +
  theme(strip.background = element_rect(),
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12, face = "bold")) +
  scale_color_brewer(palette = "Set1") +
  coord_flip() +
  labs(x = "\nDistribution",
       y = "Probability (%)\n",
       title = "LMER - Invertebrates Bray-Curtis Dissimilarity",
       subtitle = paste("Daytime: 12:00 h\nAverage probability across all window widths from -192 to 0 h"))

ggsave(paste0(figures.dir, "LMER Inverts BC check distributions - average probability.png"),
       summary.g, width = 35, height = 30, units = "cm")
