#########################
## Grassland soil moisture resilience measures
## 
## inputs:  - 'DATASET SOIL MOISTURE EXPERIMENT.csv'
## outpus:  - resilience: a data frame with the ratio of soil in each
##            plot to the mean of plant biomass in control plots at the same 
##            site and time step.
## 
## TODO:  
##
## NOTES: For calculating resilience measures for other datasets, I suggest
##    getting those datasets into a dataframe with column names that match those
##    used here.  
## 
## Willson Gaul, Lupe Leon-Sanchez, Hannah White, Maja Ilic
## created: 14 Nov 2019
## last modified: 20 Nov 2019 (Maja)
##########################

#===========================================
# Clear objects from the workspace

rm(list = ls())

# for TREATMENT, D = Drought, C = control

# 6 plots per site, 5 sites per region, 4 regions, 5 time points = 600
# Two sites are missing, B2 and B4, therefore, expected number of observations (ratios) is 540

#===========================================
# Set working directory and create folders for output data and figures

# Maja's desktop PC

setwd("C:/Users/3054311/Documents/My Documents/Grassland project/01_Probes_soil_moisture_and_temp")

data.dir <- paste0(getwd(),"/Data/")
figures.dir <- paste0(getwd(),"/Figures/")

dir.create(paste0(data.dir,"Soil moisture resilience Willson"))
mydir.data <- paste0(data.dir,"Soil moisture resilience Willson/")

dir.create(paste0(figures.dir,"Soil moisture resilience Willson"))
mydir <- paste0(figures.dir,"Soil moisture resilience Willson/")

#===========================================
# Packages

library(Hmisc)
library(tidyverse)

#===========================================
# Load soil moisture data
# This data is generated in the scripts 06 and 07 for soil moisture
# Extracted is soil moisture data on all sampling days,
# within a period of 24 h, starting at 12:00 noon 
# Note that previous on using the function, the average soil moisture
# within the period of 24 h has tu be calculated

df <- read.csv(paste0(data.dir,
                            "Soil moisture ratios DtoC per site and sampling day/Soil RH 24 h, 12-00.csv"),
                     sep = ",", header = T)

MOISTURE <- df %>% 
  group_by(Site,Region,Farm,site_ID,Plot,treatment,Day) %>% 
  summarize(variable = mean(Moisture, na.rm = T)) %>% 
  ungroup()

# start soil moisture data prep -------------
names(MOISTURE)
names(MOISTURE)[which(names(MOISTURE) == "Day")] <- "elapsed_days"

MOISTURE$timestep <- paste("day",MOISTURE$elapsed_days)
MOISTURE$metric_name <- "Moisture"

head(MOISTURE)
# end soil moisture data prep -------------

#===========================================
# Define function to calculate mean of controls

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
  ##             at which the means should be calculated
  ##        data - data frame containing the columns listed above
  ##        control_val - character string giving the value used to indicate
  ##                      controls in the treatment column
  ##        metric_name - character string with the name of the variable being
  ##                      measured (the original name of variable, e.g. plant_biomass)
  
  df <- MOISTURE %>%
    filter(treatment == eval("C") & timestep == ts) %>%
    select(site_ID, variable, timestep) %>%
    group_by(site_ID) %>%
    summarise(control_mean = mean(variable, na.rm = T), 
              timestep = unique(!!ts), 
              metric = !!metric_name) 
  df
  
} ## end function definition --------------------------------------------------

### calculate control means at each site to use -------------------------------
# make a dataframe to hold resilience metric results
# this only needs to be done for the first variable (e.g. plant biomass). 
# Results for subsequent variables (e.g. EVI, soil moister) can be joined to
# this data frame

# Replaced CODE by Plot

resilience <- MOISTURE %>% 
  select(Region, Farm, treatment, Plot, timestep, elapsed_days, site_ID, variable)

# calculate mean of control dry mass weights for all sites and time steps

control_means <- lapply(unique(MOISTURE$timestep), FUN = calc_control_mean, 
                        data = MOISTURE, control_val = "C", 
                        metric_name = "Moisture")

control_means <- bind_rows(control_means) # should be 18 * 5 = 90 rows

## only 88 rows: no data for day 64, sites Border1 (B1) and Dublin5 (D5)

## end calculate control means ----------------------------------------------
###########################################################################

resilience <- left_join(resilience, control_means)   ## Maja: why do we keep the control replicates?

### Recovery & Resistance ------------------------------------------------------
# group by Farm (equivalent to site) and Time, then calculate mean of controls 
# in the Farm & Time groups as the denominator in the resistance equation.

# calculate the ratio of the soil moisture value for each plot to the mean value for
# control plots at that site and time step

resilience$plot_meanC_ratio <- resilience$variable / 
  resilience$control_mean


g1 <- ggplot(data = resilience, aes(x = elapsed_days, y = plot_meanC_ratio, 
                              color = factor(as.character(treatment)))) + 
  geom_point() + 
  geom_smooth(method = "loess") #+ 
  # facet_wrap(~site_ID)

ggsave(paste0(mydir,"Soil resilience Willson.png"), g1, width = 8, height = 5)

## Produce similar graph like in Maja's other script
# For this, extract only drought treatment

drought <- resilience %>% 
  filter(treatment == "D")

ylab <- expression(paste(bold("Ratio"~D["i"]~":"~C["mean"])))

g2 <- ggplot(drought, aes(x = as.factor(elapsed_days), y = plot_meanC_ratio)) +
  geom_boxplot(aes(color = Region, fill = Region),
               alpha = 0.5) +
  facet_grid(Region ~ Farm) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        strip.text = element_text(size = 12),
        legend.position = "none",
        panel.spacing = unit(0.4, "cm")) +
  labs(x = "Date", y = ylab, title = "Soil resilience") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  scale_y_continuous(limits = c(0,3.5))

ggsave(paste0(mydir,"Soil resilience Willson 2.png"), g2, 
       width = 24, height = 18, units = "cm")

### end recovery & resistance calculation -------------------------------------
########## end plant biomass resilience calculation --------------------------

## Save data

write.table(resilience, paste0(mydir.data,"Soil resilience Willson.csv"), sep = ",", row.names = F)
save(resilience, file = paste0(mydir.data,"Soil resilience Willson.RData"))
