################################
## This script fits a linear model of the resilience measures for the experiment
## The response is the ratio of each measure in each plot to the control plots
## at the same site (farm) and time step.
## 
## Willson Gaul willson.gaul@ucdconnect.ie
## created: 21 Nov 2019
## last modified: 21 Jan 2020
###############################

# Get resilience metrics for different levels of the ecosystem.  Eventually we
# will do this by reading the files from the Dropbox final data folder, but
# for now I will use a mix of sourcing scripts and loading preliminary data files.

source("./plant_biomass_resilience.R")
resilience_soil <- read_csv("./Soil_moisture/Data/Soil moisture ratios DtoC per site and sampling day/Avg soil RH 24 h, 12-00.csv")

# make treatment values for match for all datasets
resilience_plant$treatment <- gsub("0", "Control", 
                                   as.character(resilience_plant$treatment))
resilience_plant$treatment <- gsub("1", "Drought", 
                                   as.character(resilience_plant$treatment))
resilience_soil$treatment <- gsub("C", "Control", 
                                  as.character(resilience_soil$treatment))
resilience_soil$treatment <- gsub("D", "Drought", 
                                  as.character(resilience_soil$treatment))

# make timestep values match for all datasets
resilience_plant$timestep <- gsub("Day.", "", resilience_plant$timestep)

# get variable classes to match
resilience_soil$Day <- as.character(resilience_soil$Day)

# add metric column to soil moisture
resilience_soil$metric <- "soil moisture"

## join measures from different ecosystem levels into a single data frame
resilience_joined <- full_join(resilience_plant[, c("Region", "Farm", "CODE", 
                                                    "site_ID", "timestep", 
                                                    "treatment", "metric", 
                                                    "plot_meanC_ratio")], 
                               resilience_soil[, c("Region", "Farm", "Plot", 
                                                   "site_ID", "Day", 
                                                   "treatment", "metric",
                                                   "Ratio")], 
                               by = c("Region", "Farm", "CODE" = "Plot", 
                                      "site_ID", "timestep" = "Day", 
                                      "treatment", "metric",
                                      "plot_meanC_ratio" = "Ratio"))
# make a column with numeric days since end of drought
resilience_joined$day <- as.numeric(as.character(resilience_joined$timestep))

# take a quick look to make sure the join worked correctly
ggplot(data = resilience_joined[!is.na(resilience_joined$day), ], 
       aes(x = as.numeric(day), 
           y = plot_meanC_ratio, 
           color = factor(as.character(treatment)))) + 
  geom_smooth() + 
  geom_point() +
  facet_wrap(~ metric)


### Fit a model --------------------------------------------------------------
## fit a model of the ratio of the plot to the mean of the control plots as 
## a function of treatment (binary) and metric (plant biomass, soil moisture, 
## etc.), and Day since end of drought (numeric)
## 
## This is not meant to be a real final model.  Rather this is just to 
## demonstrate how we could use a single outcome (the ratio of each plot to the
## mean of the control plots at its site) for all our measures. 
## 
## Eventually, we might want to add random effects (e.g. for site or region) and
## possibly other covariates.

res_mod <- glm(plot_meanC_ratio ~ as.numeric(day)*factor(treatment)*factor(metric), 
               data = resilience_joined)





