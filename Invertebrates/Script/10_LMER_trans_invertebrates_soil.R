###################################################################################
## This script fits a linear model of the invertebrate Bray-Curtis Dissimilarity ##
## in the Grassland Resilience experiment with fixed factors "day" (days after   ##
## the end of the drought) and effect size of soil moisture (calculated as soil  ##
## moisture in each plot C1-D3 minus the average soil moisture in control plots) ## 
## and random effects Region (Border, Cork, Dublin, Limerick) and Farm (1-5)     ##
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
##                                                                               ##
## Addition to this script: raw data (invertebrates BC) is transformed in order  ##
## to achieve the normal distribution of the raw data (and residuals)            ##
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
# Packages ----

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

#===========================================
# Get raw data for soil moisture ----

raw.data.dir <- "C:/Users/3054311/Documents/My Documents/Grassland project/01_Probes_soil_moisture_and_temp/Data/"

load(paste0(raw.data.dir,"Cleaned_Soil_Moisture_Data_from_Loggers.RData"))

# Events

load(file = paste0(raw.data.dir,"All events.RData"))

#================================================
# Set directory for data and plots ----

# dir.create(paste0(data.dir,"LMER model trans invertebrates soil"))
# dir.create(paste0(figures.dir,"LMER model trans invertebrates soil"))
# dir.create(paste0(figures.dir,"LMER model trans invertebrates soil/Model validation"))
# dir.create(paste0(figures.dir,"LMER model trans invertebrates soil/Model output"))
# dir.create(paste0(figures.dir,"LMER model trans invertebrates soil/Residuals"))

mydir.data <- paste0(data.dir,"LMER model trans invertebrates soil")
mydir <- paste0(figures.dir,"LMER model trans invertebrates soil")

#================================================
# Load function mov.win.diff ----

source(paste0(script.dir,"Function mov.win.diff.R"))

#================================================
# Functions for plots ----

# Histogram with Gaussian curve, Shapiro's test p-value, median, mean, variance and SD of x

source(paste0(script.dir,"Function ggHistNorm.R"))

## qqPlot made with ggplot

source(paste0(script.dir,"Function ggQQplot.R"))

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
# Transform invert BC to get a normal distribution ----

asinTransform <- function(p) { asin(sqrt(p)) }

df_inverts$inverts_BC_all_arcsin <- asinTransform(df_inverts$inverts_BC_all)
df_inverts$inverts_BC_all_log <- log(df_inverts$inverts_BC_all + 1)
df_inverts$inverts_BC_all_sqrt <- sqrt(df_inverts$inverts_BC_all)
df_inverts$inverts_BC_all_cbrt <- (df_inverts$inverts_BC_all)^(1/3)
df_inverts$inverts_BC_all_4rt <- (df_inverts$inverts_BC_all)^(1/4)

ggHistNorm(df_inverts, df_inverts$inverts_BC_all, "Bray-Curtis Dissimilarity for Invertebrates - raw data")
ggHistNorm(df_inverts, df_inverts$inverts_BC_all_arcsin, "Bray-Curtis Dissimilarity for Invertebrates - arcsin transformation")
ggHistNorm(df_inverts, df_inverts$inverts_BC_all_log, "Bray-Curtis Dissimilarity for Invertebrates - log(x+1) transformation")
ggHistNorm(df_inverts, df_inverts$inverts_BC_all_sqrt, "Bray-Curtis Dissimilarity for Invertebrates - sqrt transformation")
ggHistNorm(df_inverts, df_inverts$inverts_BC_all_cbrt, "Bray-Curtis Dissimilarity for Invertebrates - cube root transformation")
ggHistNorm(df_inverts, df_inverts$inverts_BC_all_4rt, "Bray-Curtis Dissimilarity for Invertebrates - 4th root transformation")

ggQQplot(df_inverts, df_inverts$inverts_BC_all, "Bray-Curtis Dissimilarity for Invertebrates - raw data")
ggQQplot(df_inverts, df_inverts$inverts_BC_all_arcsin, "Bray-Curtis Dissimilarity for Invertebrates - arcsin transformation")
ggQQplot(df_inverts, df_inverts$inverts_BC_all_log, "Bray-Curtis Dissimilarity for Invertebrates - log(x+1) transformation")
ggQQplot(df_inverts, df_inverts$inverts_BC_all_sqrt, "Bray-Curtis Dissimilarity for Invertebrates - sqrt transformation")
ggQQplot(df_inverts, df_inverts$inverts_BC_all_cbrt, "Bray-Curtis Dissimilarity for Invertebrates - cube root transformation")
ggQQplot(df_inverts, df_inverts$inverts_BC_all_4rt, "Bray-Curtis Dissimilarity for Invertebrates - 4th root transformation")

#================================================
# Run a for-loop for soil moisture data ----

timepoint.hr <- c()
width.hr <- c()
obs_norm <- c()
obs_var_region <- c()
obs_var_farm <- c()
obs_var_region_farm <- c()
resid_norm <- c()
resid_var_region <- c()
resid_var_farm <- c()
resid_var_region_farm <- c()
mod_AIC <- c()
mod_Rsq_marginal <- c()
mod_Rsq_conditional <- c()

timepoint <- 12
width <- seq(-192,0,1)

m <- 1

resp.variable <- "inverts_BC_all_cbrt"

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
    # Combine effect size (ES) for soil moisture and BC for invertebrates (cube root transformed) ----
    # Keep untransformed data as well
    
    names(df_soil)[names(df_soil) == "plot_meanC_ES"] <- "soil_ES"
    
    df_joined <- full_join(df_inverts[, c("Region", "Farm", "Plot", 
                                          "site_ID", "day", 
                                          "treatment", "inverts_BC_all", "inverts_BC_all_cbrt")], 
                           df_soil[, c("Region", "Farm", "Plot", 
                                       "site_ID", "day", 
                                       "treatment", "soil_ES")], 
                           by = c("Region", "Farm", "Plot", 
                                  "site_ID", "day", 
                                  "treatment"))
    
    # Exclude all rows with NAs (columns BC_meanC_all and soil_ES) ----
    
    df_joined <- df_joined[which(!is.na(df_joined$inverts_BC_all_cbrt)),]
    df_joined <- df_joined[which(!is.na(df_joined$soil_ES)),]
    
    # Change region, farm and treatment factors ----
    
    df_joined$Region <- as.factor(df_joined$Region)
    df_joined$Farm <- as.factor(df_joined$Farm)
    df_joined$treatment <- as.factor(df_joined$treatment)
    
    #================================================
    
    ###################
    ##               ##
    ##  Fit a model  ## ----
    ##               ##
    ###################
    
    # Fit a model of the ratio of the plot to the mean of the control plots as 
    # a function Day since end of drought (numeric) and the effect size of soil moisture (numeric)
    # Random factors included: Farm nested in Region (1|Region/Farm)
    
    lmer_full_mod <- lmer(inverts_BC_all_cbrt ~ day*soil_ES + (1|Region/Farm),
                          data = df_joined)
    
    #================================================
    # Add fitted values and residuals to the raw data ----
    
    df_mod <- data.frame(df_joined, "Fitted" = fitted(lmer_full_mod)^3)  # back-transform
    df_mod$Residuals <- residuals(lmer_full_mod)
    
    #================================================
    # Run anova() and extract the results ----
    
    aov_mod <- as.data.frame(anova(lmer_full_mod))
    aov_mod$Term <- rownames(aov_mod)
    aov_mod <- aov_mod[,c(7,1:6)]
    names(aov_mod)[2:7] <- c("Sum.Sq","Mean.Sq","NumDF","DenDF","F.value","p.value")
    rownames(aov_mod) <- 1:nrow(aov_mod)
    aov_mod$Daytime <- daytime.hr
    aov_mod$Duration <- duration.hr
    aov_mod$Trial <- m
    
    if(m == 1){
      aov_mod_final <- aov_mod
    }
    
    if(m > 1){
      aov_mod_final <- rbind(aov_mod_final,aov_mod)
    }
    
    #================================================
    # Extract model coefficients ----
    
    # Fixed
    
    coeff_fixef <- summary(lmer_full_mod)$coefficients
    rownames_fixef <- data.frame("Term" = rownames(coeff_fixef))
    dimnames(coeff_fixef)[[2]] <- c("Estimate","Std.Error","DF","t.value","p.value")
    fixef_out <- cbind(rownames_fixef, coeff_fixef)
    rownames(fixef_out) <- rownames(rownames_fixef)
    fixef_out$Back_trans_coeff <- fixef_out$Estimate^3
    fixef_out$Daytime <- daytime.hr
    fixef_out$Duration <- duration.hr
    fixef_out$Trial <- m
    
    if(m == 1){
      fixef_out_final <- fixef_out
    }
    
    if(m > 1){
      fixef_out_final <- rbind(fixef_out_final,fixef_out)
    }
    
    # Random: Farm:Region (Farm nested in Region) ----
    
    coeff_random.1 <- ranef(lmer_full_mod)$`Farm:Region`
    rownames_random.1 <- data.frame(rownames(coeff_random.1))
    random_out.1 <- rownames_random.1 %>% separate(rownames.coeff_random.1., c("Farm", "Region"), ":")
    random_out.1 <- data.frame(random_out.1, "Intercept" = coeff_random.1$`(Intercept)`)
    random_out.1$Intercept_1 <- random_out.1$Intercept + fixef_out$Estimate[fixef_out$Term == "(Intercept)"]  
    random_out.1$Back_trans_intercept <- random_out.1$Intercept_1^3
    random_out.1$Daytime <- daytime.hr
    random_out.1$Duration <- duration.hr
    random_out.1$Trail <- m
    
    if(m == 1){
      random_out_final.1 <- random_out.1
    }
    
    if(m > 1){
      random_out_final.1 <- rbind( random_out_final.1, random_out.1)
    }
    
    # Random: Farm ----
    
    coeff_random.2 <- ranef(lmer_full_mod)$Region
    rownames_random.2 <- rownames(coeff_random.2)
    random_out.2 <- data.frame("Region" = rownames_random.2, "Intercept" = coeff_random.2$`(Intercept)`)
    random_out.2$Intercept_1 <- random_out.2$Intercept + fixef_out$Estimate[fixef_out$Term == "(Intercept)"] 
    random_out.2$Back_trans_intercept <- random_out.2$Intercept_1^3
    random_out.2$Daytime <- daytime.hr
    random_out.2$Duration <- duration.hr
    random_out.2$Trail <- m
    
    if(m == 1){
      random_out_final.2 <- random_out.2
    }
    
    if(m > 1){
      random_out_final.2 <- rbind( random_out_final.2, random_out.2)
    }
    
    #================================================
    # Define figure titles for the plots ----
    
    if (daytime.hr < 10) {
      title <- paste0("LMER Model ",duration.hr, " h period, starting at 0",daytime.hr,":00")
      figure.title0 <- paste0(mydir,"/Model validation/Model assumptions ",duration.hr, " h, 0",daytime.hr,"-00.png")
      figure.title1 <- paste0(mydir,"/Model validation/Model validation ",duration.hr, " h, 0",daytime.hr,"-00.png")
      figure.title2.1 <- paste0(mydir,"/Residuals/Residuals ",duration.hr, " h, 0",daytime.hr,"-00.png")
      figure.title2.2 <- paste0(mydir,"/Residuals/Residuals fixed effect ",duration.hr, " h, 0",daytime.hr,"-00.png")
      figure.title2.3 <- paste0(mydir,"/Residuals/Residuals random effect ",duration.hr, " h, 0",daytime.hr,"-00.png")
      figure.title3.1 <- paste0(mydir,"/Model output/Model output A",duration.hr, " h, 0",daytime.hr,"-00.png")
      figure.title3.2 <- paste0(mydir,"/Model output/Model output B",duration.hr, " h, 0",daytime.hr,"-00.png")
    }
    
    if (daytime.hr >= 10) {
      title <- paste0("LMER Model ",duration.hr, " h period, starting at ",daytime.hr,":00")
      figure.title0 <- paste0(mydir,"/Model validation/Model assumptions ",duration.hr, " h, ",daytime.hr,"-00.png")
      figure.title1 <- paste0(mydir,"/Model validation/Model validation ",duration.hr, " h, ",daytime.hr,"-00.png")
      figure.title2.1 <- paste0(mydir,"/Residuals/Residuals ",duration.hr, " h, ",daytime.hr,"-00.png")
      figure.title2.2 <- paste0(mydir,"/Residuals/Residuals fixed effect ",duration.hr, " h, ",daytime.hr,"-00.png")
      figure.title2.3 <- paste0(mydir,"/Residuals/Residuals random effect ",duration.hr, " h, ",daytime.hr,"-00.png")
      figure.title3.1 <- paste0(mydir,"/Model output/Model output A",duration.hr, " h, ",daytime.hr,"-00.png")
      figure.title3.2 <- paste0(mydir,"/Model output/Model output B",duration.hr, " h, ",daytime.hr,"-00.png")
    }
    
    #================================================
    ## Model assumptions ----
    
    # Check normality of the data (invert_BC_all) ----
    
    check.norm <- ggHistNorm(df_mod, df_mod$inverts_BC_all_cbrt, "Invertebrates Bray-Curtis Dissimilarity Index")
    plot.norm <- check.norm[1][[1]]
    
    plot.qq <- ggQQplot(df_mod, df_mod$inverts_BC_all_cbrt, "Invertebrates Bray-Curtis Dissimilarity Index")
    
    ## Check variance homogeneity ----
    
    # ~ Region
    
    OBS_var_region <- leveneTest(inverts_BC_all_cbrt ~ Region, data = df_mod)
    
    if(OBS_var_region$`Pr(>F)`[1] < 0.001){
      main.region.0 <- "p < 0.001"
    }
    
    if(OBS_var_region$`Pr(>F)`[1] >= 0.001 & OBS_var_region$`Pr(>F)`[1] < 0.01){
      main.region.0 <- "p < 0.01"
    }
    
    if(OBS_var_region$`Pr(>F)`[1] >= 0.01 & OBS_var_region$`Pr(>F)`[1] < 0.05){
      main.region.0 <- "p < 0.05"
    }
    
    if(OBS_var_region$`Pr(>F)`[1] >= 0.05 ){
      main.region.0 <- parse(text = paste0('p == ', round(OBS_var_region$`Pr(>F)`[1], digits = 3)))
    }
    
    box.obs.region <- ggplot(df_mod, aes(x = Region, y = inverts_BC_all_cbrt, fill = Region)) +
      geom_boxplot(outlier.shape = 21, alpha = 0.7) +
      theme_minimal() +
      theme(legend.position = "none") +
      labs(x = "Region",
           y = "Bray-Curtis Dissimilarity Index\nInvertebrates",
           title = "Invertebrates Bray-Crutis Dissimilarity Index ~ Region",
           subtitle = main.region.0)
    
    # ~ Farm
    
    OBS_var_farm <- leveneTest(inverts_BC_all_cbrt ~ Farm, data = df_mod)
    
    if(OBS_var_farm$`Pr(>F)`[1] < 0.001){
      main.farm.0 <- "p < 0.001"
    }
    
    if(OBS_var_farm$`Pr(>F)`[1] >= 0.001 & OBS_var_farm$`Pr(>F)`[1] < 0.01){
      main.farm.0 <- "p < 0.01"
    }
    
    if(OBS_var_farm$`Pr(>F)`[1] >= 0.01 & OBS_var_farm$`Pr(>F)`[1] < 0.05){
      main.farm.0 <- "p < 0.05"
    }
    
    if(OBS_var_farm$`Pr(>F)`[1] >= 0.05 ){
      main.farm.0 <- parse(text = paste0('p == ', round(OBS_var_farm$`Pr(>F)`[1], digits = 3)))
    }
    
    box.obs.farm <- ggplot(df_mod, aes(x = Farm, y = inverts_BC_all_cbrt, fill = Farm)) +
      geom_boxplot(outlier.shape = 21, alpha = 0.7) +
      scale_fill_brewer(palette = "GnBu") +
      theme_minimal() +
      theme(legend.position = "none") +
      labs(x = "Farm",
           y = "Bray-Curtis Dissimilarity Index\nInvertebrates",
           title = "Invertebrates Bray-Crutis Dissimilarity Index ~ Farm",
           subtitle = main.farm.0)
    
    # ~ Region * Farm
    
    OBS_var_region_farm <- leveneTest(inverts_BC_all_cbrt ~ Region*Farm, data = df_mod)
    
    if(OBS_var_region_farm$`Pr(>F)`[1] < 0.001){
      main.region.farm.0 <- "p < 0.001"
    }
    
    if(OBS_var_region_farm$`Pr(>F)`[1] >= 0.001 & OBS_var_region_farm$`Pr(>F)`[1] < 0.01){
      main.region.farm.0 <- "p < 0.01"
    }
    
    if(OBS_var_region_farm$`Pr(>F)`[1] >= 0.01 & OBS_var_region_farm$`Pr(>F)`[1] < 0.05){
      main.region.farm.0 <- "p < 0.05"
    }
    
    if(OBS_var_region_farm$`Pr(>F)`[1] >= 0.05 ){
      main.region.farm.0 <- parse(text = paste0('p == ', round(OBS_var_region_farm$`Pr(>F)`[1], digits = 3)))
    }
    
    box.obs.region.farm <- ggplot(df_mod, aes(x = Region, y = inverts_BC_all_cbrt, fill = Farm)) +
      geom_boxplot(outlier.shape = 21, alpha = 0.7) +
      scale_fill_brewer(palette = "GnBu") +
      theme_minimal() +
      theme(legend.position = "none") +
      labs(x = "Region",
           y = "Bray-Curtis Dissimilarity Index\nInvertebrates",
           title = "Invertebrates Bray-Crutis Dissimilarity Index ~ Region * Farm",
           subtitle = main.region.farm.0)
    
    #================================================
    ## Model validation ----
    
    # Observed vs fitted values (fitted values were back transformed!)
    
    # ~ Region
    
    obs.fit.region <- ggplot(df_mod, aes(x = Fitted, y = inverts_BC_all, color = Region, fill = Region)) +
      geom_point(shape = 21, size = 2, alpha = 0.5) +
      facet_wrap(~ Region, ncol = 2) +
      geom_abline(slope = 1, intercept = 0) +
      theme_minimal() +
      theme(strip.background = element_rect(color = "grey50"),
            legend.position = "none") +
      labs(x = "Fitted",
           y = "Bray-Curtis Dissimilarity Index\nInvertebrates",
           title = "Observed vs. Fitted",
           subtitle = "~ Region")
    
    # ~ Farm
    
    obs.fit.farm <- ggplot(df_mod, aes(x = Fitted, y = inverts_BC_all, fill = Farm)) +
      geom_point(shape = 21, size = 2, alpha = 0.5, color = "grey50") +
      scale_fill_brewer(palette = "GnBu") +
      facet_wrap(~ Farm, ncol = 3) +
      geom_abline(slope = 1, intercept = 0) +
      theme_minimal() +
      theme(strip.background = element_rect(color = "grey50"),
            legend.position = "none") +
      labs(x = "Fitted",
           y = "Bray-Curtis Dissimilarity Index\nInvertebrates",
           title = "Observed vs. Fitted",
           subtitle = "~ Farm")
    
    # ~ Region * Farm
    
    obs.fit.region.farm <- ggplot(df_mod, aes(x = Fitted, y = inverts_BC_all, 
                                              color = Region, fill = Region)) +
      geom_point(shape = 21, size = 2, alpha = 0.5) +
      facet_grid(Region ~ Farm) +
      geom_abline(slope = 1, intercept = 0) +
      theme_minimal() +
      theme(strip.background = element_rect(color = "grey50"),
            legend.position = "none") +
      labs(x = "Fitted",
           y = "Bray-Curtis Dissimilarity Index\nInvertebrates",
           title = "Observed vs. Fitted",
           subtitle = "~ Region * Farm")
    
    #================================================
    # Collinearity ----
    
    result.coll <- check_collinearity(lmer_full_mod)
    
    result.coll[which(result.coll$VIF < 5),"Correlation"] <- "low"
    result.coll[which(result.coll$VIF >= 5 & result.coll$VIF < 10),"Correlation"] <- "moderate"
    result.coll[which(result.coll$VIF >= 10),"Correlation"] <- "high"
    
    result.coll$Correlation <- as.factor(result.coll$Correlation)
    result.coll$Correlation <- factor(result.coll$Correlation, levels = c("low", "moderate", "high"))
    
    mycol.coll <- c(rgb(39, 174, 96, max = 255),
                    rgb(230, 126, 34, max = 255),
                    rgb(228, 26, 28, max = 255))
    
    corr.levels <- unique(result.coll$Correlation)
    
    mycol <- c()
    
    if ("low" %in% corr.levels){
      mycol <- mycol.coll[1]
    }
    
    if ("moderate" %in% corr.levels){
      mycol <- c(mycol,mycol.coll[2])
    }
    
    if ("high" %in% corr.levels){
      mycol <- c(mycol,mycol.coll[3])
    }
    
    plot.coll <- ggplot(result.coll, aes(x = Parameter, y = VIF)) +
      geom_bar(stat = "identity", width = 0.7, aes(fill = Correlation)) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
      scale_fill_manual(values = mycol) +
      labs(x = "Parameter", y = "VIF", title = "Check for Multicollinearity", subtitle = "")
    
    #================================================
    
    #################
    ##             ##
    ##  Residuals  ## ----
    ##             ##
    #################
    
    # Binned residuals ----
    # Function binned_residuals produces a plot if not saved to an object
    # However, that plot can't be further modified
    # Therefore, I save the output to an object and recreated the plot 
    
    result.binned <- binned_residuals(lmer_full_mod)
    resid.inside.err <- sum(result.binned$group == "yes")/nrow(result.binned)*100
    resid.inside.err <- round(resid.inside.err, digits = 2)
    
    plot.binned <- ggplot(result.binned, aes(x = xbar*100, y = ybar)) +
      geom_ribbon(aes(ymin = -Inf, ymax = -se),
                  color = "grey80", fill = "grey95", alpha = 0.5) +
      geom_ribbon(aes(ymin = se, ymax = +Inf),
                  color = "grey80", fill = "grey95", alpha = 0.5) +
      geom_hline(yintercept = 0, color = "grey80") +
      geom_point(aes(color = group), size = 3) +
      theme_bw() +
      scale_color_brewer(palette = "Set1") + 
      labs(x = paste0("Estimated probability of ", resp.variable), 
           y = "Average residual", 
           title = "Binned residuals", 
           subtitle = paste0(resid.inside.err, "% of the residuals are inside the error bounds."))
    
    # ~ Fitted
    
    result.heteroscedasticity <- check_heteroscedasticity(lmer_full_mod)
    
    if(result.heteroscedasticity[1] < 0.001){
      p.res.fit <- "p < 0.001"
    }
    
    if(result.heteroscedasticity[1] >= 0.001 & result.heteroscedasticity[1] < 0.01){
      p.res.fit <- "p < 0.01"
    }
    
    if(result.heteroscedasticity[1] >= 0.01 & result.heteroscedasticity[1] < 0.05){
      p.res.fit <- "p < 0.05"
    }
    
    if(result.heteroscedasticity[1] >= 0.05 ){
      p.res.fit <- parse(text = paste0('p == ', round(result.heteroscedasticity[1], digits = 3)))
    }
    
    plot.res.fit <- ggplot(df_mod, aes(x = Fitted, y = Residuals)) +
      geom_point(size = 2, color = rgb(44, 62, 80, max = 255)) +
      theme_minimal() +
      geom_smooth(method = "loess", size = 1, color = rgb(228, 26, 28, max = 255), se = F) +
      labs(x = "Fitted",
           y = "Residuals",
           title = "Residuals vs. Fitted",
           subtitle = p.res.fit)

    # ~ Region
    
    RESID_var_region <- leveneTest(Residuals ~ Region, data = df_mod)
    
    if(RESID_var_region$`Pr(>F)`[1] < 0.001){
      main.region <- "p < 0.001"
    }
    
    if(RESID_var_region$`Pr(>F)`[1] >= 0.001 & RESID_var_region$`Pr(>F)`[1] < 0.01){
      main.region <- "p < 0.01"
    }
    
    if(RESID_var_region$`Pr(>F)`[1] >= 0.01 & RESID_var_region$`Pr(>F)`[1] < 0.05){
      main.region <- "p < 0.05"
    }
    
    if(RESID_var_region$`Pr(>F)`[1] >= 0.05 ){
      main.region <- parse(text = paste0('p == ', round(RESID_var_region$`Pr(>F)`[1], digits = 3)))
    }
    
    box.region <- ggplot(df_mod, aes(x = Region, y = Residuals, fill = Region)) +
      geom_boxplot(outlier.shape = 21, alpha = 0.7) +
      theme_minimal() +
      theme(legend.position = "none") +
      labs(x = "Region",
           y = "Residuals",
           title = "Residuals ~ Region",
           subtitle = main.region)
    
    # ~ Farm
    
    RESID_var_farm <- leveneTest(Residuals ~ Farm, data = df_mod)
    
    if(RESID_var_farm$`Pr(>F)`[1] < 0.001){
      main.farm <- "p < 0.001"
    }
    
    if(RESID_var_farm$`Pr(>F)`[1] >= 0.001 & RESID_var_farm$`Pr(>F)`[1] < 0.01){
      main.farm <- "p < 0.01"
    }
    
    if(RESID_var_farm$`Pr(>F)`[1] >= 0.01 & RESID_var_farm$`Pr(>F)`[1] < 0.05){
      main.farm <- "p < 0.05"
    }
    
    if(RESID_var_farm$`Pr(>F)`[1] >= 0.05 ){
      main.farm <- parse(text = paste0('p == ', round(RESID_var_farm$`Pr(>F)`[1], digits = 3)))
    }
    
    box.farm <- ggplot(df_mod, aes(x = Farm, y = Residuals, fill = Farm)) +
      geom_boxplot(outlier.shape = 21, alpha = 0.7) +
      scale_fill_brewer(palette = "GnBu") +
      theme_minimal() +
      theme(legend.position = "none") +
      labs(x = "Farm",
           y = "Residuals",
           title = "Residuals ~ Farm",
           subtitle = main.farm)
    
    # ~ Region * Farm
    
    RESID_var_region_farm <- leveneTest(Residuals ~ Region*Farm, data = df_mod)
    
    if(RESID_var_region_farm$`Pr(>F)`[1] < 0.001){
      main.region.farm <- "p < 0.001"
    }
    
    if(RESID_var_region_farm$`Pr(>F)`[1] >= 0.001 & RESID_var_region_farm$`Pr(>F)`[1] < 0.01){
      main.region.farm <- "p < 0.01"
    }
    
    if(RESID_var_region_farm$`Pr(>F)`[1] >= 0.01 & RESID_var_region_farm$`Pr(>F)`[1] < 0.05){
      main.region.farm <- "p < 0.05"
    }
    
    if(RESID_var_region_farm$`Pr(>F)`[1] >= 0.05 ){
      main.region.farm <- parse(text = paste0('p == ', round(RESID_var_region_farm$`Pr(>F)`[1], digits = 3)))
    }
    
    box.region.farm <- ggplot(df_mod, aes(x = Region, y = Residuals, fill = Farm)) +
      geom_boxplot(outlier.shape = 21, alpha = 0.7) +
      scale_fill_brewer(palette = "GnBu") +
      theme_minimal() +
      theme(legend.position = "none") +
      labs(x = "Region",
           y = "Residuals",
           title = "Residuals ~ Region * Farm",
           subtitle = main.region.farm)
    
    #================================================
    # Check for normal distribution of residuals ----
    
    check.resid.norm <- ggHistNorm(df_mod, df_mod$Residuals, "Residuals", 0.05)
    plot.norm.resid <- check.resid.norm[1][[1]]
    
    plot.qq.resid <- ggQQplot(df_mod, df_mod$Residuals, "Residuals")
    
    #================================================
    # Normality of Random Effects ----
    
    result.mod <- check_model(lmer_full_mod)
    
    # ~ Region
    
    REQQ.region <- result.mod$REQQ$Region
    
    reqq.plot.region <- ggplot(REQQ.region, aes(x = x, y = y)) +
      geom_point(size = 2, color = rgb(44, 62, 80, max = 255)) +
      geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) + 
      geom_smooth(method = "lm", color = rgb(22, 160, 133, max = 255), size = 1, se = F) +
      theme_minimal() +
      labs(x = "Theoretical Quantiles",
           y = "RE Quantiles",
           title = "Normality of Random Effects",
           subtitle = "Region")
    
    # ~ Farm:Region
    
    REQQ.farm.region <- result.mod$REQQ$"Farm:Region"
    
    reqq.plot.farm.region <- ggplot(REQQ.farm.region, aes(x = x, y = y)) +
      geom_point(size = 2, color = rgb(44, 62, 80, max = 255)) +
      geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) + 
      geom_smooth(method = "lm", color = rgb(22, 160, 133, max = 255), size = 1, se = F) +
      theme_minimal() +
      labs(x = "Theoretical Quantiles",
           y = "RE Quantiles",
           title = "Normality of Random Effects",
           subtitle = "Farm:Region")
    
    #================================================
    
    #################################################
    ##                                             ##
    ##  Plot raw data vs fixed and random effects  ## ----
    ##                                             ##
    #################################################
    
    final.plot0 <- grid.arrange(plot.norm,
                                plot.qq,
                                box.obs.region,
                                box.obs.farm,
                                box.obs.region.farm,
                                ncol = 2,
                                top = tableGrob(t("lmer(inverts_BC_all_cbrt ~ day*soil_ES + (1|Region/Farm)"),
                                                theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                       base_colour = "black",
                                                                       base_size = 12)))
    
    ggsave(figure.title0,
           final.plot0 <- grid.arrange(final.plot0,
                                       top = tableGrob(t(title),
                                                       theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                              base_colour = "black",
                                                                              base_size = 16))),
           width = 25, height = 32, units = "cm")
    
    #================================================
    
    ################################################
    ##                                            ##
    ##  Plot observed vs fitted and collinearity  ## ----
    ##                                            ##
    ################################################
    
    final.plot1 <- grid.arrange(obs.fit.region,
                                obs.fit.farm,
                                obs.fit.region.farm,
                                plot.coll,
                                ncol = 2,
                                top = tableGrob(t("lmer(inverts_BC_all_cbrt ~ day*soil_ES + (1|Region/Farm)"),
                                                theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                       base_colour = "black",
                                                                       base_size = 12)))
    
    ggsave(figure.title1,
           final.plot1 <- grid.arrange(final.plot1,
                                       top = tableGrob(t(title),
                                                       theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                              base_colour = "black",
                                                                              base_size = 16))),
           width = 30, height = 26, units = "cm")
    
    #================================================
    
    ################################
    ##                            ##
    ##  Plot residuals vs fitted  ## ----
    ##                            ##
    ################################

    final.plot2 <- grid.arrange(plot.norm.resid,
                                plot.qq.resid,
                                plot.res.fit,
                                plot.binned,
                                ncol = 2,
                                top = tableGrob(t("lmer(inverts_BC_all_cbrt ~ day*soil_ES + (1|Region/Farm)"),
                                                theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                       base_colour = "black",
                                                                       base_size = 12)))
    
    ggsave(figure.title2.1,
           final.plot2 <- grid.arrange(final.plot2,
                                     top = tableGrob(t(title),
                                                     theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                            base_colour = "black",
                                                                            base_size = 16))),
           width = 25, height = 25, units = "cm")
    
    #================================================
    
    ########################################
    ##                                    ##
    ##  Plot residuals vs random effects  ## ----
    ##                                    ##
    ########################################
    
    final.plot3 <- grid.arrange(arrangeGrob(box.region, box.farm, ncol = 2),
                                box.region.farm, nrow = 2) 
    
    final.plot3 <- grid.arrange(final.plot3,
                                arrangeGrob(reqq.plot.region, reqq.plot.farm.region, ncol = 2),
                                nrow = 2,
                                heights = c(2.2,1),
                                top = tableGrob(t("lmer(inverts_BC_all_cbrt ~ day*soil_ES + (1|Region/Farm)"),
                                                theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                       base_colour = "black",
                                                                       base_size = 12)))
    
    ggsave(figure.title2.3,
           final.plot3 <- grid.arrange(final.plot3,
                                       top = tableGrob(t(title),
                                                       theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                              base_colour = "black",
                                                                              base_size = 16))),
           width = 25, height = 32, units = "cm")
    
    #================================================
    
    ################################
    ##                            ##
    ##  Plot data - Model output  ## ----
    ##                            ##
    ################################
    
    ylab <- expression(paste(bold("Invertebrates Bray-Curtis Dissimilarity")))
    
    ## Panels: Region ~ Farm
    # Fitted data
    
    g1 <- ggplot(df_mod, 
                 aes(x = day, y = Fitted)) +
      geom_point(aes(fill = Region, color = Region),
                 size = 2, alpha = 0.2) +
      facet_grid(Region ~ Farm) +
      geom_smooth(aes(color = Region),
                  method = "lm", se = F) +
      theme_minimal() +
      theme(strip.background = element_rect(),
            panel.spacing = unit(1, "lines")) +
      scale_shape_manual(values = c(21,22,23,24,25)) +
      labs(x = "Days since the end of the drought",
           title = "Fitted data back transformed",
           subtitle = "lm(fitted data ~ day)")
    
    # Raw data
    
    g2 <- ggplot(df_mod, 
                 aes(x = day, y = inverts_BC_all)) +
      geom_point(aes(fill = Region, color = Region),
                 size = 2, alpha = 0.2) +
      facet_grid(Region ~ Farm) +
      geom_smooth(aes(y = Fitted, color = Region),
                  method = "lm", se = F) +
      theme_minimal() +
      theme(strip.background = element_rect(),
            panel.spacing = unit(1, "lines")) +
      scale_shape_manual(values = c(21,22,23,24,25)) +
      labs(x = "Days since the end of the drought",
           title = "Raw (untransformed) data",
           subtitle = "lm(fitted data ~ day)")
    
    ## Panels: Region
    # Fitted data

    g3 <- ggplot(df_mod, 
                 aes(x = day, y = Fitted)) +
      geom_point(aes(fill = Region, color = Region),
                 size = 3, alpha = 0.2) +
      facet_wrap(~ Region) +
      geom_smooth(aes(color = Region, linetype = Farm),
                  method = "lm", se = F, size = 0.7) +
      geom_smooth(method = "lm", se = F,
                  color = "black", size = 1) +
      theme_minimal() +
      theme(strip.background = element_rect(),
            panel.spacing = unit(1, "lines")) +
      scale_shape_manual(values = c(21,22,23,24,25)) +
      labs(x = "Days since the end of the drought",
           title = "Fitted data back transformed",
           subtitle = "lm(fitted data ~ day)")
    
    # Raw data
    
    g4 <- ggplot(df_mod, 
                 aes(x = day, y = inverts_BC_all)) +
      geom_point(aes(fill = Region, color = Region),
                 size = 3, alpha = 0.2) +
      facet_wrap(~ Region) +
      geom_smooth(aes(y = Fitted, color = Region, linetype = Farm),
                  method = "lm", se = F, size = 0.7) +
      geom_smooth(aes(y = Fitted),
                  method = "lm", se = F,
                  color = "black", size = 1) +
      theme_minimal() +
      theme(strip.background = element_rect(),
            panel.spacing = unit(1, "lines")) +
      scale_shape_manual(values = c(21,22,23,24,25)) +
      labs(x = "Days since the end of the drought",
           title = "Raw (untransformed) data",
           subtitle = "lm(fitted data ~ day)")
    
    # All four plots together
    
    final.plot4 <- grid.arrange(g2,g4,g1,g3,
                                ncol = 2,
                                top = tableGrob(t("lmer(inverts_BC_all_cbrt ~ day*soil_ES + (1|Region/Farm)"),
                                                theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                       base_colour = "black",
                                                                       base_size = 12)))
    
    
    ggsave(figure.title3.1,
           final.plot4 <- grid.arrange(final.plot4,
                                       top = tableGrob(t(title),
                                                       theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                              base_colour = "black",
                                                                              base_size = 16))),
           width = 40, height = 30, units = "cm")
    
    #================================================
    ## Panels: Day
    # Fitted data ~ soil moisture ----
    
    g5 <- ggplot(df_mod, 
                 aes(x = soil_ES, y = Fitted)) +
      geom_point(aes(fill = Region, color = Region),
                 size = 3, alpha = 0.2) +
      facet_grid(Region ~ day) +
      geom_smooth(aes(color = Region, linetype = Farm),
                  method = "lm", se = F, size = 0.7) +
      geom_smooth(method = "lm", se = F,
                  color = "black", size = 1) +
      theme_minimal() +
      theme(strip.background = element_rect(),
            panel.spacing = unit(1, "lines")) +
      scale_shape_manual(values = c(21,22,23,24,25)) +
      labs(x = "Soil moisture - Effect size",
           title = "Fitted data back transformed",
           subtitle = "lm(fitted data ~ soil moisture effect size)")
    
    # Raw data
    
    g6 <- ggplot(df_mod, 
                 aes(x = soil_ES, y = inverts_BC_all_cbrt)) +
      geom_point(aes(fill = Region, color = Region),
                 size = 3, alpha = 0.2) +
      facet_grid(Region ~ day) +
      geom_smooth(aes(y = Fitted, color = Region, linetype = Farm),
                  method = "lm", se = F, size = 0.7) +
      geom_smooth(aes(y = Fitted),
                  method = "lm", se = F,
                  color = "black", size = 1) +
      theme_minimal() +
      theme(strip.background = element_rect(),
            panel.spacing = unit(1, "lines")) +
      scale_shape_manual(values = c(21,22,23,24,25)) +
      labs(x = "Soil moisture - Effect size",
           title = "Raw (untransformed) data",
           subtitle = "lm(fitted data ~ soil moisture effect size)")
    
    # Combine plots
    
    final.plot5 <- grid.arrange(g6,g5,
                                ncol = 1,
                                top = tableGrob(t("lmer(inverts_BC_all_cbrt ~ day*soil_ES + (1|Region/Farm)"),
                                                theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                       base_colour = "black",
                                                                       base_size = 12)))
    
    
    ggsave(figure.title3.2,
           final.plot5 <- grid.arrange(final.plot5,
                                       top = tableGrob(t(title),
                                                       theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                              base_colour = "black",
                                                                              base_size = 16))),
           width = 24, height = 30, units = "cm")
    
    #================================================
    
    #######################
    ##                   ##
    ##  Save statistics  ## ----
    ##                   ##
    #######################
    
    timepoint.hr[m] <- daytime.hr
    width.hr[m] <- duration.hr
    obs_norm[m] <- check.norm[2][[1]]
    obs_var_region[m] <- OBS_var_region$`Pr(>F)`[1]
    obs_var_farm[m] <- OBS_var_farm$`Pr(>F)`[1]
    obs_var_region_farm[m] <- OBS_var_region_farm$`Pr(>F)`[1]
    resid_norm[m] <- check.resid.norm[2][[1]]
    resid_var_region[m] <- RESID_var_region$`Pr(>F)`[1]
    resid_var_farm[m] <- RESID_var_farm$`Pr(>F)`[1]
    resid_var_region_farm[m] <- RESID_var_region_farm$`Pr(>F)`[1]
    mod_AIC[m] <- AIC(lmer_full_mod)
    
    ##########################
    ##                      ##
    ##  Calculation of Rsq  ## ----
    ##                      ##
    ##########################
    
    # For details, see Nakagawa and Schielzeth, 2012
    # https://doi.org/10.1111/j.2041-210x.2012.00261.x
    # Approach later modified by Johnson, 2014
    # https://doi.org/10.1111/2041-210X.12225
    # var(f):     variance of the fixed effects
    # var(r):     variance of the radnom effects
    # var(e):     variance of the model residuals
    
    # Marginal
    # var(f) / [var(f) + var(r) + var(e)]
    
    mod_Rsq_marginal[m] <- model_performance(lmer_full_mod)$R2_marginal
      
    # Conditional
    # [var(f) + var(r)] / [var(f) + var(r) + var(e)]
      
    mod_Rsq_conditional[m] <- model_performance(lmer_full_mod)$R2_conditional
    
    # Increase m
    
    m <- m + 1
    
  }
}

#================================================
# Cobine all summary vectors into a data frame

mod_summary <- data.frame(timepoint.hr, width.hr, 
                          obs_norm, obs_var_region, obs_var_farm, obs_var_region_farm,
                          resid_norm, resid_var_region, resid_var_farm, resid_var_region_farm,
                          mod_AIC, mod_Rsq_marginal, mod_Rsq_conditional)

#================================================
# Export summary and output of the models

write.table(mod_summary, paste0(mydir.data, "/LMER Model summary - cbrt inverts BC vs soil moisture.csv"), 
            sep = ",", row.names = F)

write.table(aov_mod_final, paste0(mydir.data, "/LMER Anova output - cbrt inverts BC vs soil moisture.csv"), 
            sep = ",", row.names = F)

write.table(fixef_out_final, paste0(mydir.data, "/LMER Model output fixef - cbrt inverts BC vs soil moisture.csv"), 
            sep = ",", row.names = F)

write.table(random_out_final.1, paste0(mydir.data, "/LMER Model output random Region - cbrt inverts BC vs soil moisture.csv"), 
            sep = ",", row.names = F)

write.table(random_out_final.2, paste0(mydir.data, "/LMER Model output random Farm Region - cbrt inverts BC vs soil moisture.csv"), 
            sep = ",", row.names = F)

#================================================
# Explore the summary

par(mfrow = c(1,1), oma = c(0,0,0,0))

mar.Rsq <- expression(paste(bold("Marginal")~bolditalic("R"^"2")))
con.Rsq <- expression(paste(bold("Conditional")~bolditalic("R"^"2")))

# AIC vs. Window width and marginal R2

gp1 <- ggplot(mod_summary, aes(x = width.hr, y = mod_AIC)) +
  geom_point(aes(color = mod_Rsq_marginal, fill = mod_Rsq_marginal), shape = 21, size = 4, alpha = 0.5) +
  scale_color_viridis_c(name = mar.Rsq) +
  scale_fill_viridis_c(name = mar.Rsq) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold")) + 
  labs(x = expression(paste(bold("Window width (h)"))),
       y = expression(paste(bold("AIC"))),
       title = "LMER - Invertebrates Bray-Curtis Dissimilarity (cube root)",
       subtitle = "Fixed: day * soil moisture effect size\nRandom: Region/Farm\nDaytime: 12:00 h") +
  geom_vline(xintercept = -0.5, color = "grey80", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = -15.5, color = "grey80", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = -38.5, color = "grey80", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = -79.5, color = "grey80", linetype = "dashed", size = 0.5)

ggsave(paste0(mydir, "/LMER Model validation cbrt Inverts vs Soil Moisture - AIC vs Window width and mar R2.png"),
       gp1, width = 20, height = 16, units = "cm")

# AIC vs. Window width and conditional R2

gp2 <- ggplot(mod_summary, aes(x = width.hr, y = mod_AIC)) +
  geom_point(aes(color = mod_Rsq_conditional, fill = mod_Rsq_conditional), shape = 21, size = 4, alpha = 0.5) +
  scale_color_viridis_c(name = con.Rsq) +
  scale_fill_viridis_c(name = con.Rsq) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold")) + 
  labs(x = expression(paste(bold("Window width (h)"))),
       y = expression(paste(bold("AIC"))),
       title = "LMER - Invertebrates Bray-Curtis Dissimilarity (cube root)",
       subtitle = "Fixed: day * soil moisture effect size\nRandom: Region/Farm\nDaytime: 12:00 h") +   
  geom_vline(xintercept = -0.5, color = "grey80", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = -15.5, color = "grey80", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = -38.5, color = "grey80", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = -79.5, color = "grey80", linetype = "dashed", size = 0.5)

ggsave(paste0(mydir, "/LMER Model validation cbrt Inverts vs Soil Moisture - AIC vs Window width and con R2.png"),
       gp2, width = 20, height = 16, units = "cm")

# Both combined

gp1 <- gp1 + theme(legend.position = "bottom") +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5,
                               direction = "horizontal", barwidth = unit(5, "cm")),
         color = guide_colorbar(title.position = "top", title.hjust = 0.5,
                                direction = "horizontal", barwidth = unit(5, "cm")))

gp2 <- gp2 + theme(legend.position = "bottom") +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5,
                               direction = "horizontal", barwidth = unit(5, "cm")),
         color = guide_colorbar(title.position = "top", title.hjust = 0.5,
                                direction = "horizontal", barwidth = unit(5, "cm")))

ggsave(paste0(mydir, "/LMER Model validation cbrt Inverts vs Soil Moisture - AIC vs Window width R2.png"),
       gp1 | gp2, width = 28, height = 16, units = "cm")

# Marginal R2 vs. AIC

gp3 <- ggplot(mod_summary, aes(x = mod_AIC, y = mod_Rsq_marginal, 
                               fill = width.hr, color = width.hr)) +
  geom_point(size = 4, shape = 21, alpha = 0.5) +
  scale_fill_viridis_c(name = "Window width (h)") +
  scale_color_viridis_c(name = "Window width (h)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) + 
  labs(x = expression(paste(bold("AIC"))),
       y = mar.Rsq,
       title = "LMER - Invertebrates Bray-Curtis Dissimilarity (cube root)",
       subtitle = "Fixed: day * soil moisture effect size\nRandom: Region/Farm\nDaytime: 12:00 h")

ggsave(paste0(mydir, "/LMER Model validation cbrt Inverts vs Soil Moisture - Marginal R2 vs AIC.png"),
       gp3, width = 20, height = 16, units = "cm")

# Conditional R2 vs. AIC

gp4 <- ggplot(mod_summary, aes(x = mod_AIC, y = mod_Rsq_conditional,                                
                               fill = width.hr, color = width.hr)) +
  geom_point(size = 4, shape = 21, alpha = 0.5) +
  scale_fill_viridis_c(name = "Window width (h)") +
  scale_color_viridis_c(name = "Window width (h)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) + 
  labs(x = expression(paste(bold("AIC"))),
       y = con.Rsq,
       title = "LMER - Invertebrates Bray-Curtis Dissimilarity (cube root)",
       subtitle = "Fixed: day * soil moisture effect size\nRandom: Region/Farm\nDaytime: 12:00 h")

ggsave(paste0(mydir, "/LMER Model validation cbrt Inverts vs Soil Moisture - Conditional R2 vs AIC.png"),
       gp4, width = 20, height = 16, units = "cm")

#################
##             ##
##  Residuals  ##
##             ##
#################

ylab <- expression(paste(bolditalic("P"),bold("-value")))

# Shapiro.test

p1 <- ggplot(mod_summary, aes(x = width.hr, y = resid_norm)) +
  geom_point(aes(color = mod_Rsq_conditional, fill = mod_Rsq_conditional, size = mod_AIC), 
             shape = 21, alpha = 0.5) +
  scale_color_viridis_c(name = con.Rsq) +
  scale_fill_viridis_c(name = con.Rsq) + 
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  labs(x = "Window width (h)", y = ylab, 
       title = "LMER - Invertebrates Bray-Curtis Dissimilarity (cube root)\nShapiro's Test Residuals",
       subtitle = "Fixed: day * soil moisture effect size\nRandom: Region/Farm\nDaytime: 12:00 h")

ggsave(paste0(mydir, "/LMER Model validation cbrt Inverts vs Soil Moisture - Shapiro test.png"),
       p1, width = 20, height = 16, units = "cm")

# leveneTest ~ Region

p2 <- ggplot(mod_summary, aes(x = width.hr, y = resid_var_region)) +
  geom_point(aes(color = mod_Rsq_conditional, fill = mod_Rsq_conditional, size = mod_AIC), 
             shape = 21, alpha = 0.5) +
  scale_color_viridis_c(name = con.Rsq) +
  scale_fill_viridis_c(name = con.Rsq) + 
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  labs(x = "Window width (h)", y = ylab, 
       title = "LMER - Invertebrates Bray-Curtis Dissimilarity (cube root)\nResiduals ~ Region",
       subtitle = "Fixed: day * soil moisture effect size\nRandom: Region/Farm\nDaytime: 12:00 h")

ggsave(paste0(mydir, "/LMER Model validation cbrt Inverts vs Soil Moisture - Levene Test Region.png"),
       p2, width = 20, height = 16, units = "cm")

# leveneTest ~ Farm

p3 <- ggplot(mod_summary, aes(x = width.hr, y = resid_var_farm)) +
  geom_point(aes(color = mod_Rsq_conditional, fill = mod_Rsq_conditional, size = mod_AIC), 
             shape = 21, alpha = 0.5) +
  scale_color_viridis_c(name = con.Rsq) +
  scale_fill_viridis_c(name = con.Rsq) + 
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  labs(x = "Window width (h)", y = ylab, 
       title = "LMER - Invertebrates Bray-Curtis Dissimilarity (cube root)\nResiduals ~ Farm",
       subtitle = "Fixed: day * soil moisture effect size\nRandom: Region/Farm\nDaytime: 12:00 h")

ggsave(paste0(mydir, "/LMER Model validation Inverts vs Soil Moisture - Levene Test Farm.png"),
       p3, width = 20, height = 16, units = "cm")

# leveneTest ~ Region * Farm

p4 <- ggplot(mod_summary, aes(x = width.hr, y = resid_var_region_farm)) +
  geom_point(aes(color = mod_Rsq_conditional, fill = mod_Rsq_conditional, size = mod_AIC), 
             shape = 21, alpha = 0.5) +
  scale_color_viridis_c(name = con.Rsq) +
  scale_fill_viridis_c(name = con.Rsq) + 
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  labs(x = "Window width (h)", y = ylab, 
       title = "LMER - Invertebrates Bray-Curtis Dissimilarity (cube root)\nResiduals ~ Region * Farm",
       subtitle = "Fixed: day * soil moisture effect size\nRandom: Region/Farm\nDaytime: 12:00 h")

ggsave(paste0(mydir, "/LMER Model validation cbrt Inverts vs Soil Moisture - Levene Test Region Farm.png"),
       p4, width = 20, height = 16, units = "cm")

################
##            ##
##  p-values  ##
##            ##
################

names(mod_summary)[1:2] <- c("Daytime","Duration")

mod_output <- left_join(aov_mod_final,
                        mod_summary[,c("Daytime","Duration","mod_AIC","mod_Rsq_marginal","mod_Rsq_conditional")], 
                        by = c("Daytime","Duration"))

mod_output$Term <- as.factor(mod_output$Term)
mod_output$Term <- factor(mod_output$Term, levels = unique(mod_output$Term))

# p-values vs. window width

p5 <- ggplot(mod_output, aes(x = Duration, y = p.value)) +
  geom_point(aes(color = mod_Rsq_marginal, fill = mod_Rsq_marginal), 
             shape = 21, alpha = 0.3, size = 3) +
  facet_wrap(~ Term, ncol = 3) +
  scale_color_viridis_c(name = mar.Rsq) +
  scale_fill_viridis_c(name = mar.Rsq) + 
  theme_minimal() +
  theme(strip.background = element_rect(),
        plot.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  labs(x = "Window width (h)", y = ylab, 
       title = "LMER - Invertebrates Bray-Curtis Dissimilarity (cube root)",
       subtitle = "Fixed: day * soil moisture effect size\nRandom: Region/Farm\nDaytime: 12:00 h") +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_y_continuous(limits = c(0,0.8), 
                     breaks = c(0,0.05,0.2,0.4,0.6,0.8),
                     labels = c("0",expression(paste(bold("0.05"))),"0.2","0.4","0.6","0.8"))

ggsave(paste0(mydir, "/LMER Model validation cbrt Inverts vs Soil Moisture - P values vs Window width.png"),
       p5, width = 24, height = 14, units = "cm")

# p-values as barchart (significant vs non significant)

mod_output[which(mod_output$p.value < 0.05), "Significant"] <- "yes"
mod_output[which(mod_output$p.value >= 0.05), "Significant"] <- "no"

mod_sig <- mod_output %>% 
  group_by(Term,Significant) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n) * 100)

p6 <- ggplot(mod_sig, aes(x = Term, y = freq)) +
  geom_bar(stat = "identity", width = 0.8, aes(fill = Significant)) +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(strip.background = element_rect(),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold")) +
  labs(x = "Term", y = "Frequency (%)", 
       title = "LMER - Invertebrates Bray-Curtis Dissimilarity (cube root)",
       subtitle = "Fixed: day * soil moisture effect size\nRandom: Region/Farm\nDaytime: 12:00 h")

ggsave(paste0(mydir, "/LMER Model validation cbrt Inverts vs Soil Moisture - Significance of the terms.png"),
       p6, width = 18, height = 15, units = "cm")
