#########################
## Grassland plant biomass resilience measures
## 
## TODO:  - add control_i to control mean ratios
## 
## Willson Gaul, Lupe Leon-Sanchez, Hannah White
## created: 14 Nov 2019
## last modified: 14 Nov 2019
##########################
rm(list = ls())
# RM = repeated measurement
# Can drop the final duplicate time and dry weight rm columns

# for TREATMENT, 1 = Drought, 0 = control

# 8 plots per site, 5 sites per region, 4 regions, 7 time points

setwd("~/Documents/UCD/Grassland_project/")

library(Hmisc)
library(tidyverse)

BIOMASA <- read_csv("~/Dropbox/Grassland/Data/EMPIRICAL DATA FIELD EXPERIMENT/DATASET BIOMASS SOIL MOISTURE EXPERIMENT.csv")

BIOMASA <- BIOMASA[, colnames(BIOMASA) %nin% c("Time_1", "DRY WEIGHT RM_1")]

BIOMASA$TREATMENT <- as.character(BIOMASA$TREATMENT)
BIOMASA$Farm <- as.character(BIOMASA$Farm)
BIOMASA$CODE <- as.character(BIOMASA$CODE)

# Make a unique site identifier
BIOMASA$site_ID <- paste0(as.character(BIOMASA$Region), as.character(BIOMASA$Farm))

# make column with days elapsed since end of drought (for plotting)
BIOMASA$elapsed_days <- gsub("Day.", "", BIOMASA$Time)
BIOMASA$elapsed_days[BIOMASA$elapsed_days == "Baseline"] <- NA
BIOMASA$elapsed_days[BIOMASA$elapsed_days == "mid.drought"] <- NA
BIOMASA$elapsed_days <- as.numeric(as.character(BIOMASA$elapsed_days))

# make a dataframe to hold resilience metric results
resilience <- select(BIOMASA, Region, Farm, TREATMENT, CODE, Time, elapsed_days, 
                     site_ID, `DRY WEIGHT NON RM`)

#####################################################################
### calculate control means at each site to use -------------------------------
calc_control_mean <- function(timestep, data, variable) {
  ## calculate the mean of controls plots for variable at a single timestep
  df <- BIOMASA %>%
    filter(TREATMENT == "0" & Time == timestep) %>%
    select(site_ID, eval(variable), Time) %>%
    group_by(site_ID) %>%
    mutate(control_mean_dry_weight_non_rm = mean(`DRY WEIGHT NON RM`, na.rm = T)) %>%
    select(-`DRY WEIGHT NON RM`) %>%
    unique()
  df
}

# calculate mean of control dry mass weights for all sites and time steps
control_means <- lapply(unique(BIOMASA$Time), FUN = calc_control_mean, 
                        data = BIOMASA, variable = "DRY WEIGHT NON RM")
names(control_means) <- unique(BIOMASA$Time)

control_means <- bind_rows(control_means)
## end calculate control means ----------------------------------------------
###########################################################################

resilience <- left_join(resilience, control_means)

#### Resistance -------------------------------------------------------------
# group by Farm (equivalent to site) and Time, then calculate mean of controls 
# in the Farm & Time groups as the demoninator in the resistance equation.
# control_mean <- BIOMASA %>%
#   filter(TREATMENT == "0" & Time == "Day.0") %>%
#   select(site_ID, `DRY WEIGHT NON RM`) %>%
#   group_by(site_ID) %>%
#   mutate(control_mean_dry_weight_non_rm = mean(`DRY WEIGHT NON RM`, na.rm = T)) %>%
#   select(-`DRY WEIGHT NON RM`) %>%
#   unique()

# calculate resistance at each plot
# resist <- left_join(BIOMASA, control_means) %>%
#   filter(TREATMENT == "1" & Time == "Day.0") %>%
#   select(-`SOIL MOISTURE`, -`DRY WEIGHT RM`, -Farm, -Region) %>% 
#   # group_by(site_ID) %>% 
#   mutate(resistance = `DRY WEIGHT NON RM` / control_mean_dry_weight_non_rm) %>%
#   select(-TREATMENT, -`DRY WEIGHT NON RM`)
# 
# resilience <- left_join(resilience, resist)
### end resistance calculation -----------------------------------------------

### Recovery & Resistance ------------------------------------------------------
## Use ratio at day 64 instead of rate?
# For recovery, we are using non-repeated measures of biomass

# calculate the ratio of the biomass value for each plot to the mean value for
# control plots at that site and time step
resilience$plot_meanC_ratio <- resilience$`DRY WEIGHT NON RM` / 
  resilience$control_mean_dry_weight_non_rm



ggplot(data = resilience, aes(x = elapsed_days, y = plot_meanC_ratio, 
                              color = factor(as.character(TREATMENT)))) + 
  geom_point() + 
  geom_smooth(method = "loess")

### end recovery -------------------------------------------------------------
