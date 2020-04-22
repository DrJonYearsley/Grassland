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
## Author of the modified script:                                                ##
## Maja Ilic M.Ilic@qub.ac.uk                                                    ##
## first modified: 20 Feb 2020                                                   ##
## last modified: 20 Feb 2020                                                    ##
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

dir.create(paste0(data.dir,"LMER model invertebrates comparison"))
dir.create(paste0(figures.dir,"LMER model invertebrates comparison"))

mydir.data <- paste0(data.dir,"LMER model invertebrates comparison/")
mydir <- paste0(figures.dir,"LMER model invertebrates comparison/")

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
    # Compare all models ----
    
    result <- compare_performance(lmer_day,
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
    
    df.results <- as.data.frame(result)
    df.results$Daytime <- 12
    df.results$Duration <- duration.hr
    df.results$Trial <- m
    
    if (m == 1) {
      df.results.final <- df.results
    }
    
    if (m > 1) {
      df.results.final <- rbind(df.results.final, df.results)
    }
    
    myplot <- plot(result)
    
    #================================================
    # Compare models with only one fixed factor ----
    
    result.1.fixef <- compare_performance(lmer_day,
                                          lmer_soil,
                                          lmer_trt,
                                          lmer_region)
    
    plot.1.fixef <- plot(result.1.fixef)
    
    #================================================
    # Compare models with 2 fixed factors ----
    
    result.2.fixef <- compare_performance(lmer_day_soil,
                                          lmer_day_trt,
                                          lmer_day_region,
                                          lmer_soil_trt,
                                          lmer_soil_region,
                                          lmer_trt_region)
    
    plot.2.fixef <- plot(result.2.fixef)
    
    #================================================
    # Compare models with 3 fixed factor ----
    
    result.3.fixef <- compare_performance(lmer_day_soil_trt,
                                          lmer_day_soil_region,
                                          lmer_day_trt_region,
                                          lmer_soil_trt_region)
    
    plot.3.fixef <- plot(result.3.fixef)
    
    #================================================
    # Compare models with Farm nested in Region as random effect ----
    
    result.farm.region <- compare_performance(lmer_day,
                                              lmer_soil,
                                              lmer_trt,
                                              lmer_day_soil,
                                              lmer_day_trt,
                                              lmer_soil_trt,
                                              lmer_day_soil_trt)
    
    plot.farm.region <- plot(result.farm.region)
    
    #================================================
    # Compare models with Farm as random effect ----
    
    result.farm <- compare_performance(lmer_region,
                                       lmer_day_region,
                                       lmer_soil_region,
                                       lmer_trt_region,
                                       lmer_day_soil_region,
                                       lmer_day_trt_region,
                                       lmer_soil_trt_region,
                                       lmer_day_soil_trt_region)
    
    plot.farm <- plot(result.farm)
    
    #================================================
    # Plot the comparison ----
    
    title <- paste("Window width:",duration.hr,"h")
    figure.title <- paste0(mydir, "Model comparison/LMER model comparison for Invertebrates, window width ",duration.hr," h.png")
    
    final.plot <- grid.arrange(myplot,
                               plot.1.fixef,
                               plot.2.fixef,
                               plot.3.fixef,
                               plot.farm.region,
                               plot.farm,
                               ncol = 2,
                               top = tableGrob(t(title),
                                               theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                      base_colour = "black",
                                                                      base_size = 16)))
    
    ggsave(figure.title,
           final.plot <- grid.arrange(final.plot,
                                      top = tableGrob(t("lmer(inverts_BC_all ~ fixed_effects + random_effects)"),
                                                      theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                             base_colour = "black",
                                                                             base_size = 16))),
           width = 30, height = 30, units = "cm")
    
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

write.table(df.results.final, paste0(mydir.data, "LMER model invertebrates comparison.csv"), 
            sep = ",", row.names = F)
    
#================================================
# Plot results ----

mycolors <- colorRampPalette(brewer.pal(10, "Spectral"))(length(unique(df.results.final$Model)))
show_col(mycolors)

df.results.final$Model <- as.factor(df.results.final$Model)
df.results.final$Model <- factor(df.results.final$Model, levels = unique(df.results.final$Model))

# AIC vs. window width 

g1 <- ggplot(df.results.final, aes(x = Duration, y = AIC, fill = Model, color = Model)) +
  geom_point(size = 3, shape = 21, alpha = 0.3) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_fill_manual(values = mycolors) +
  scale_color_manual(values = mycolors) +
  labs(x = "Window width (h)", 
       title = "LMER - Invertebrates Bray-Curtis Dissimilarity",
       subtitle = "Daytime: 12:00 h") +
  geom_vline(xintercept = -0.5, color = "grey80", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = -15.5, color = "grey80", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = -38.5, color = "grey80", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = -79.5, color = "grey80", linetype = "dashed", size = 0.5)

ggsave(paste0(mydir,"LMER model invertebrates comparison - AIC vs. window width.png"),
       g1, width = 20, height = 15, units = "cm")

# R2 marginal vs. window width

g2 <- ggplot(df.results.final, aes(x = Duration, y = R2_marginal, fill = Model, color = Model)) +
  geom_point(size = 3, shape = 21, alpha = 0.3) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_fill_manual(values = mycolors) +
  scale_color_manual(values = mycolors) +
  labs(x = "Window width (h)", y = expression(paste(bold("Marginal")~bolditalic("R"^"2"))),
       title = "LMER - Invertebrates Bray-Curtis Dissimilarity",
       subtitle = "Daytime: 12:00 h") +
  geom_vline(xintercept = -0.5, color = "grey80", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = -15.5, color = "grey80", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = -38.5, color = "grey80", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = -79.5, color = "grey80", linetype = "dashed", size = 0.5)

ggsave(paste0(mydir,"LMER model invertebrates comparison - R2 marginal vs. window width.png"),
       g2, width = 20, height = 15, units = "cm")

# R2 conditional vs. window width

g3 <- ggplot(df.results.final, aes(x = Duration, y = R2_conditional, fill = Model, color = Model)) +
  geom_point(size = 3, shape = 21, alpha = 0.3) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_fill_manual(values = mycolors) +
  scale_color_manual(values = mycolors) +
  labs(x = "Window width (h)", y = expression(paste(bold("Conditional")~bolditalic("R"^"2"))),
       title = "LMER - Invertebrates Bray-Curtis Dissimilarity",
       subtitle = "Daytime: 12:00 h") +
  geom_vline(xintercept = -0.5, color = "grey80", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = -15.5, color = "grey80", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = -38.5, color = "grey80", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = -79.5, color = "grey80", linetype = "dashed", size = 0.5)

ggsave(paste0(mydir,"LMER model invertebrates comparison - R2 conditional vs. window width.png"),
       g3, width = 20, height = 15, units = "cm")

# AIC vs. window width and R2 marginal

g4 <- ggplot(df.results.final, aes(x = Duration, y = AIC, fill = R2_marginal, color = R2_marginal)) +
  geom_point(size = 3, shape = 21, alpha = 0.3) +
  facet_wrap(~ Model) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = c(0.875,0.1),
        strip.background = element_rect()) +
  scale_fill_viridis_c(name = expression(paste(bold("Marginal")~bolditalic("R"^"2")))) +
  scale_color_viridis_c(name = expression(paste(bold("Marginal")~bolditalic("R"^"2")))) +
  labs(x = "Window width (h)", y = "AIC",
       title = "LMER - Invertebrates Bray-Curtis Dissimilarity",
       subtitle = "Daytime: 12:00 h") +
  geom_vline(xintercept = -0.5, color = "grey80", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = -15.5, color = "grey80", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = -38.5, color = "grey80", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = -79.5, color = "grey80", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, color = "grey80", size = 0.5) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5,
                               direction = "horizontal", barwidth = unit(5, "cm")),
         color = guide_colorbar(title.position = "top", title.hjust = 0.5,
                                direction = "horizontal", barwidth = unit(5, "cm")))

ggsave(paste0(mydir,"LMER model invertebrates comparison - AIC vs. window width and mar R2.png"),
       g4, width = 25, height = 20, units = "cm")

# AIC vs. window width and R2 conditional

g5 <- ggplot(df.results.final, aes(x = Duration, y = AIC, fill = R2_conditional, color = R2_conditional)) +
  geom_point(size = 3, shape = 21, alpha = 0.3) +
  facet_wrap(~ Model) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = c(0.875,0.1),
        strip.background = element_rect()) +
  scale_fill_viridis_c(name = expression(paste(bold("Conditional")~bolditalic("R"^"2")))) +
  scale_color_viridis_c(name = expression(paste(bold("Conditional")~bolditalic("R"^"2")))) +
  labs(x = "Window width (h)", y = "AIC",
       title = "LMER - Invertebrates Bray-Curtis Dissimilarity",
       subtitle = "Daytime: 12:00 h") +
  geom_vline(xintercept = -0.5, color = "grey80", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = -15.5, color = "grey80", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = -38.5, color = "grey80", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = -79.5, color = "grey80", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, color = "grey80", size = 0.5) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5,
                               direction = "horizontal", barwidth = unit(5, "cm")),
         color = guide_colorbar(title.position = "top", title.hjust = 0.5,
                                direction = "horizontal", barwidth = unit(5, "cm")))

ggsave(paste0(mydir,"LMER model invertebrates comparison - AIC vs. window width and con R2.png"),
       g5, width = 25, height = 20, units = "cm")

# AIC vs. R2 marginal

g6 <- ggplot(df.results.final, aes(x = R2_marginal, y = AIC, fill = Model, color = Model)) +
  geom_point(size = 3, shape = 21, alpha = 0.3) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_fill_manual(values = mycolors) +
  scale_color_manual(values = mycolors) +
  labs(x = expression(paste(bold("Marginal")~bolditalic("R"^"2"))),
       title = "LMER - Invertebrates Bray-Curtis Dissimilarity",
       subtitle = "Daytime: 12:00 h")

ggsave(paste0(mydir,"LMER model invertebrates comparison - AIC vs. R2 marginal.png"),
       g6, width = 20, height = 15, units = "cm")

# AIC vs. R2 conditional

g7 <- ggplot(df.results.final, aes(x = R2_conditional, y = AIC, fill = Model, color = Model)) +
  geom_point(size = 3, shape = 21, alpha = 0.3) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_fill_manual(values = mycolors) +
  scale_color_manual(values = mycolors) +
  labs(x = expression(paste(bold("Conditional")~bolditalic("R"^"2"))),
       title = "LMER - Invertebrates Bray-Curtis Dissimilarity",
       subtitle = "Daytime: 12:00 h")

ggsave(paste0(mydir,"LMER model invertebrates comparison - AIC vs. R2 conditional.png"),
       g7, width = 20, height = 15, units = "cm")
