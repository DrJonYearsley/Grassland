#########################
## Grassland plant biomass resilience measures
## 
## inputs:  - 'DATASET BIOMASS SOIL MOISTURE EXPERIMENT.csv'
## outpus:  - resilience: a data frame with the ratio of plant biomass in each
##              plot to the mean of plant biomass in control plots at the same 
##              site and time step.
## 
## TODO:  - add ratio for repeated measures samples
##
## NOTES: For calculating reslilience measures for other datasets, I suggest
##    getting those datasets into a dataframe with colum names that match those
##    used here.  
## 
## Willson Gaul, Lupe Leon-Sanchez, Hannah White
## created: 14 Nov 2019
## last modified: 15 Nov 2019
##########################
rm(list = ls())
# RM means repeated measurement
# Can drop the final duplicate time and dry weight rm columns that are in the 
# csv file.

# for TREATMENT, 1 = Drought, 0 = control

# 8 plots per site, 5 sites per region, 4 regions, 7 time points = 1120

setwd("~/Documents/Data_Analysis/UCD/Grassland/")

library(Hmisc)
library(tidyverse)

BIOMASA <- read_csv("~/Dropbox/Grassland/Data/FinalExptData/DATASET BIOMASS SOIL MOISTURE EXPERIMENT.csv")

## define function to calculate mean of controls
calc_control_mean <- function(ts, data, control_val, metric_name) {
  ## calculate the mean of controls plots for variable at a single timestep
  ##
  ## This function requires data to have the following column names:
  ##  timestep - indicating the sampling time (e.g. day 0)
  ##  variable - this is the variable for which the mean of control plots is 
  ##              desired (e.g. plant_biomass)
  ##  treatment - giving the control/treatment label for each observation
  ##  site_ID - unique site (farm) identifier (site not plot)
  ##  
  ##  ARGS: ts - character string with a single value from the timestep column 
  ##            at which the means should be calculated
  ##        data - data frame containing the columns listed above
  ##        control_val - character string giving the value used to indicate
  ##            controls in the treatment column
  ##        metric_name - character string with the name of the variable being
  ##            measured (the original name of variable, e.g. plant_biomass)
  df <- BIOMASA %>%
    filter(treatment == eval(control_val) & timestep == ts) %>%
    select(site_ID, variable, timestep) %>%
    group_by(site_ID) %>%
    summarise(control_mean = mean(variable, na.rm = T), 
              timestep = unique(!!ts), 
              metric = !!metric_name) 
  df
} ## end function definition --------------------------------------------------

# Plant biomass data format prep ---------------------------------------------
# Make a unique site identifier
BIOMASA$site_ID <- paste0(as.character(BIOMASA$Region), as.character(BIOMASA$Farm))
# standardize column names
names(BIOMASA)[names(BIOMASA) == "Time"] <- "timestep" # standardize column names
names(BIOMASA)[names(BIOMASA) == "TREATMENT"] <- "treatment" # standardize column names


# make column with days elapsed since end of drought (for plotting)
BIOMASA$elapsed_days <- gsub("Day.", "", BIOMASA$timestep)
BIOMASA$elapsed_days[BIOMASA$elapsed_days == "Baseline"] <- NA
BIOMASA$elapsed_days[BIOMASA$elapsed_days == "mid.drought"] <- NA
BIOMASA$elapsed_days <- as.numeric(as.character(BIOMASA$elapsed_days))
# copy the variable of interest to a column named 'variable' as required by function
BIOMASA$variable <- BIOMASA$`DRY WEIGHT NON RM` 
## end plant biomass data prep -----------------------------------------------


### calculate control means at each site to use -------------------------------
# make a dataframe to hold resilience metric results
# this only needs to be done for the first variable (e.g. plant biomass). 
# Results for subsequent variables (e.g. EVI, soil moister) can be joined to
# this data frame
resilience <- select(BIOMASA, Region, Farm, treatment, CODE, timestep, elapsed_days, 
                     site_ID, `DRY WEIGHT NON RM`)

# calculate mean of control dry mass weights for all sites and time steps
control_means <- lapply(unique(BIOMASA$timestep), FUN = calc_control_mean, 
                        data = BIOMASA, control_val = "0", 
                        metric_name = "DRY WEIGHT NON RM")

control_means <- bind_rows(control_means) # should be 140 rows
## end calculate control means ----------------------------------------------
###########################################################################

resilience <- left_join(resilience, control_means)

### Recovery & Resistance ------------------------------------------------------
# group by Farm (equivalent to site) and Time, then calculate mean of controls 
# in the Farm & Time groups as the demoninator in the resistance equation.
# 
# We are using non-repeated measures of biomass

# calculate the ratio of the biomass value for each plot to the mean value for
# control plots at that site and time step
resilience$plot_meanC_ratio <- resilience$`DRY WEIGHT NON RM` / 
  resilience$control_mean


ggplot(data = resilience, aes(x = elapsed_days, y = plot_meanC_ratio, 
                              color = factor(as.character(treatment)))) + 
  geom_point() + 
  geom_smooth(method = "loess") #+ 
  # facet_wrap(~site_ID)

### end recovery & resistance calculation -------------------------------------
########## end plant biomass resilience calculation --------------------------

# rename data frame so it is clear this is the plant biomass data
resilience_plant <- resilience
