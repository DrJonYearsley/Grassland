###################################################################################
## This script fits a linear model of the invertebrates resilience measure with  ##
## the fixed effect day (days since the end of the simulated drought; covariate) ##
## and accumulated rain, which was kept out from the drought plots during the    ##
## simulated drought (covariate) and with random effect Farm (1-5) nested in     ##
## Region (Border, Cork, Dublin, Limerick)                                       ##
## The response is the invertebrate Bray-Curtis Dissimilarity Index, calculated  ##
## by comparing the invertebrate community in each plot with the average         ##
## invertebrate community of the control plots within the same site and time     ##
## step (reference community). BC-Index of 0 indicates that the two compared     ##
## communities are equal and share all species, while BC-Index of 1 indicates    ##
## that the two compared communities do not share any species                    ##
##                                                                               ##
## Author of the script:                                                         ##
## Maja Ilic M.Ilic@qub.ac.uk                                                    ##
## first modified: 06 Mar 2020                                                   ##
## last modified: 06 Mar 2020                                                    ##
###################################################################################

#================================================
# Clear objects from the workspace ----

rm(list = ls())

#================================================
# Set working directory ----

# Maja's desktop PC

setwd("C:/Users/3054311/Documents/My Documents/Grassland project/05_Resilience")

data.dir <- paste0(getwd(),"/Data/")
figures.dir <- paste0(getwd(),"/Figures/")
script.dir <- paste0(getwd(),"/Script/")

#================================================
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
library(grid)
library(cowplot)
library(patchwork)
library(gridExtra)
library(effects)
library(afex)
library(RColorBrewer)

#================================================
# Set directory for data and plots ----

dir.create(paste0(figures.dir,"LMER model inverts BC day rain out 2"))

mydir.data <- data.dir
mydir <- paste0(figures.dir,"LMER model inverts BC day rain out 2")

#================================================
# Import invertebrate data: Bray-Curtis Dissimilarity

inverts_data <- read.csv("C:/Users/3054311/Documents/My Documents/Grassland project/07_Invertebrates/Data/Invertebrates with Bray Curtis.csv", sep = ",", header = T)

# Extract only relevant columns

resilience_inverts <- inverts_data[,c(1:6,36)]

# Change the column "BC_meanC_all" to "inverts_BC_all"

names(resilience_inverts)[which(names(resilience_inverts) == "BC_meanC_all")] <- "inverts_BC_all"

# Change the column Farm to a factor

resilience_inverts$Farm <- as.factor(resilience_inverts$Farm)

# Remove any NAs

resilience_inverts <- resilience_inverts %>% filter(!is.na(inverts_BC_all))

# Define response variable 

resp.variable <- "Invertebrates Bray-Curtis Dissimilarity Index"

#================================================
# Load rainfall data (MERA) ----

load(file = "C:/Users/3054311/Documents/My Documents/Grassland project/04_MERA_rainfall_and_temp/Data/MERA_Accumulated_percipitation_during_drought.RData")
rain.acc$Farm <- as.factor(rain.acc$Farm)

# Combine rainfall with inverts data

resilience_inverts <- left_join(resilience_inverts,
                                rain.acc[,c("Region","Farm","Acc.Rain.D")],
                                by = c("Region","Farm"))

# Rename Acc.Rain.D to rain_out

names(resilience_inverts)[names(resilience_inverts) == "Acc.Rain.D"] <- "rain_out"

# Replace accumulated rainfall with zero for all drought plots

resilience_inverts$rain_out[resilience_inverts$treatment == "Drought"] <- 0

resilience_inverts$Region <- as.factor(resilience_inverts$Region)

#================================================
# Functions for plots ----

# Histogram with Gaussian curve, Shapiro's test p-value, median, mean, variance and SD of x

source("C:/Users/3054311/Documents/My Documents/Grassland project/Functions/Function ggHistNorm.R")

## qqPlot made with ggplot

source("C:/Users/3054311/Documents/My Documents/Grassland project/Functions/Function ggQQplot.R")

#================================================
# Fit a model ----

# Fit a model of the invertebrates Bray-Curtis Dissimilarity Index as 
# a function Day since end of drought (numeric) and accumulated rain kept out 
# of the drought plots during the simulated drought 
# Random factors included: Farm nested in Region (1|Region/Farm)

lmer_full_mod <- lmer(inverts_BC_all ~ day*rain_out + (1|Region/Farm),
                      data = resilience_inverts)

#================================================
# Add fitted values and residuals to the raw data ----

resilience_mod <- data.frame(resilience_inverts, "Fitted" = fitted(lmer_full_mod))  
resilience_mod$Residuals <- residuals(lmer_full_mod)

#================================================
# Run anova() and extract the results ----

aov_mod <- as.data.frame(anova(lmer_full_mod))
aov_mod$Term <- rownames(aov_mod)
aov_mod <- aov_mod[,c(7,1:6)]
names(aov_mod)[2:7] <- c("Sum.Sq","Mean.Sq","NumDF","DenDF","F.value","p.value")
rownames(aov_mod) <- 1:nrow(aov_mod)

#================================================
# Extract model coefficients ----

# Fixed

coeff_fixef <- summary(lmer_full_mod)$coefficients
rownames_fixef <- data.frame("Term" = rownames(coeff_fixef))
dimnames(coeff_fixef)[[2]] <- c("Estimate","Std.Error","DF","t.value","p.value")
fixef_out <- cbind(rownames_fixef, coeff_fixef)
rownames(fixef_out) <- rownames(rownames_fixef)

# Random: Farm:Region (Farm nested in Region) ----

coeff_random.1 <- ranef(lmer_full_mod)$`Farm:Region`
rownames_random.1 <- data.frame(rownames(coeff_random.1))
random_out.1 <- rownames_random.1 %>% separate(rownames.coeff_random.1., c("Farm", "Region"), ":")
random_out.1 <- data.frame(random_out.1, "Intercept" = coeff_random.1$`(Intercept)`)
random_out.1$Intercept_1 <- random_out.1$Intercept + fixef_out$Estimate[fixef_out$Term == "(Intercept)"] 

# Random: Region ----

coeff_random.2 <- ranef(lmer_full_mod)$Region
rownames_random.2 <- rownames(coeff_random.2)
random_out.2 <- data.frame("Region" = rownames_random.2, "Intercept" = coeff_random.2$`(Intercept)`)
random_out.2$Intercept_1 <- random_out.2$Intercept + fixef_out$Estimate[fixef_out$Term == "(Intercept)"] 

#================================================
# Model validation ----

ylab <- "Invertebrates Bray-Curtis\nDissimilarity Index"

# Observed vs fitted values ----

# ~ day

obs.fit.day <- ggplot(resilience_mod, aes(x = Fitted, y = inverts_BC_all, fill = as.factor(day))) +
  geom_point(shape = 21, size = 2, alpha = 0.5, color = "grey50") +
  facet_wrap(~ day, ncol = 5) +
  geom_abline(slope = 1, intercept = 0) +
  theme_minimal() +
  theme(strip.background = element_rect(color = "grey50"),
        legend.position = "none",
        axis.title = element_text(size = 12, face = "bold")) +
  scale_fill_brewer(palette = "BuPu") +
  labs(x = "\nFitted",
       y = ylab,
       title = "Observed vs. Fitted",
       subtitle = "~ day")

# ~ rain_out

rain.lab <- expression(paste(bold("Accumulated Precipitation (kg "*m^"-2"*")")))

midpoint <- min(resilience_inverts$rain_out[resilience_inverts$rain_out > 0], na.rm = T)

obs.fit.rain <- ggplot(resilience_mod, aes(x = Fitted, 
                                           y = inverts_BC_all, 
                                           fill = rain_out,
                                           color = rain_out)) +
  geom_point(shape = 21, size = 2, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0) +
  theme_minimal() +
  theme(strip.background = element_rect(color = "grey50"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "top") +
  scale_fill_gradient2(name = rain.lab, midpoint = midpoint, 
                       low = "red", mid = "lightblue", high = "blue") +
  scale_color_gradient2(name = rain.lab, midpoint = midpoint, 
                        low = "red", mid = "lightblue", high = "blue") +
  labs(x = "\nFitted",
       y = ylab,
       title = "Observed vs. Fitted",
       subtitle = "~ rain_out") +
  guides(fill = guide_colorbar(title.position = "top", barwidth = unit(6,"cm")),
         color = guide_colorbar(title.position = "top", barwidth = unit(6,"cm")))

# ~ Region

obs.fit.region <- ggplot(resilience_mod, aes(x = Fitted, y = inverts_BC_all, 
                                             color = Region, fill = Region)) +
  geom_point(shape = 21, size = 2, alpha = 0.5) +
  facet_wrap(~ Region, ncol = 2) +
  geom_abline(slope = 1, intercept = 0) +
  theme_minimal() +
  theme(strip.background = element_rect(color = "grey50"),
        legend.position = "none",
        axis.title = element_text(size = 12, face = "bold")) +
  labs(x = "\nFitted",
       y = ylab,
       title = "Observed vs. Fitted",
       subtitle = "~ Region")

# ~ Farm

obs.fit.farm <- ggplot(resilience_mod, aes(x = Fitted, y = inverts_BC_all, fill = Farm)) +
  geom_point(shape = 21, size = 2, alpha = 0.5, color = "grey50") +
  scale_fill_brewer(palette = "GnBu") +
  facet_wrap(~ Farm, ncol = 3) +
  geom_abline(slope = 1, intercept = 0) +
  theme_minimal() +
  theme(strip.background = element_rect(color = "grey50"),
        legend.position = "none",
        axis.title = element_text(size = 12, face = "bold")) +
  labs(x = "\nFitted",
       y = ylab,
       title = "Observed vs. Fitted",
       subtitle = "~ Farm")

# ~ Region * Farm

obs.fit.region.farm <- ggplot(resilience_mod, aes(x = Fitted, y = inverts_BC_all, 
                                          color = Region, fill = Region)) +
  geom_point(shape = 21, size = 2, alpha = 0.5) +
  facet_grid(Region ~ Farm) +
  geom_abline(slope = 1, intercept = 0) +
  theme_minimal() +
  theme(strip.background = element_rect(color = "grey50"),
        legend.position = "none",
        axis.title = element_text(size = 12, face = "bold")) +
  labs(x = "\nFitted",
       y = ylab,
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
  labs(x = "\nParameter", y = "VIF", title = "Check for Multicollinearity", subtitle = "")

#================================================
# Residuals ----
# Binned residuals
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
  labs(x = paste0("\nEstimated probability of ", resp.variable), 
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

plot.res.fit <- ggplot(resilience_mod, aes(x = Fitted, y = Residuals)) +
  geom_point(size = 2, color = rgb(44, 62, 80, max = 255)) +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1, color = rgb(228, 26, 28, max = 255), se = F) +
  labs(x = "\nFitted",
       y = "Residuals",
       title = "Residuals vs. Fitted",
       subtitle = p.res.fit)

# ~ day

RESID_var_day <- leveneTest(Residuals ~ as.factor(day), data = resilience_mod)

if(RESID_var_day$`Pr(>F)`[1] < 0.001){
  main.day <- "p < 0.001"
}

if(RESID_var_day$`Pr(>F)`[1] >= 0.001 & RESID_var_day$`Pr(>F)`[1] < 0.01){
  main.day <- "p < 0.01"
}

if(RESID_var_day$`Pr(>F)`[1] >= 0.01 & RESID_var_day$`Pr(>F)`[1] < 0.05){
  main.day <- "p < 0.05"
}

if(RESID_var_day$`Pr(>F)`[1] >= 0.05 ){
  main.day <- parse(text = paste0('p == ', round(RESID_var_day$`Pr(>F)`[1], digits = 3)))
}

box.day <- ggplot(resilience_mod, aes(x = as.factor(day), 
                                      y = Residuals, 
                                      fill = as.factor(day))) +
  geom_boxplot(outlier.shape = 21) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12, face = "bold")) +
  scale_fill_brewer(palette = "BuPu") +
  labs(x = "\nDays since the end of the simulated drought",
       y = "Residuals",
       title = "Residuals ~ day",
       subtitle = main.day)

# ~ Region

RESID_var_region <- leveneTest(Residuals ~ Region, data = resilience_mod)

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

box.region <- ggplot(resilience_mod, aes(x = Region, y = Residuals, fill = Region)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.7) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12, face = "bold")) +
  labs(x = "\nRegion",
       y = "Residuals",
       title = "Residuals ~ Region",
       subtitle = main.region)

# ~ Farm

RESID_var_farm <- leveneTest(Residuals ~ Farm, data = resilience_mod)

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

box.farm <- ggplot(resilience_mod, aes(x = Farm, y = Residuals, fill = Farm)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.7) +
  scale_fill_brewer(palette = "GnBu") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12, face = "bold")) +
  labs(x = "\nFarm",
       y = "Residuals",
       title = "Residuals ~ Farm",
       subtitle = main.farm)

# ~ Region * Farm

RESID_var_region_farm <- leveneTest(Residuals ~ Region*Farm, data = resilience_mod)

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

box.region.farm <- ggplot(resilience_mod, aes(x = Region, y = Residuals, fill = Farm)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.7) +
  scale_fill_brewer(palette = "GnBu") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12, face = "bold")) +
  labs(x = "\nRegion",
       y = "Residuals",
       title = "Residuals ~ Region * Farm",
       subtitle = main.region.farm)

#================================================
# Check for normal distribution of residuals ----

check.resid.norm <- ggHistNorm(resilience_mod, resilience_mod$Residuals, "Residuals", 0.05)
plot.norm.resid <- check.resid.norm[1][[1]]

plot.qq.resid <- ggQQplot(resilience_mod, resilience_mod$Residuals, "Residuals")

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
# Plot observed vs. fitted ----

title.model <- "lmer(inverts_BC_all~ day*rain_out + (1|Region/Farm))"

figure.title3 <- paste0(mydir,"/Observed vs Fitted V1.png")

final.plot3 <- grid.arrange(obs.fit.day,
                            obs.fit.rain,
                            obs.fit.region,
                            obs.fit.farm, 
                            ncol = 2,
                            top = tableGrob(t(title.model),
                                            theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                   base_colour = "black",
                                                                   base_size = 16)))

ggsave(figure.title3,
       grid.arrange(arrangeGrob(final.plot3,
                                top = textGrob(t(resp.variable),
                                               gp = gpar(fontsize = 18, fontface = "bold")))),
       width = 25, height = 27, units = "cm")

# Farm:Region

obs.fit.region.farm <- grid.arrange(obs.fit.region.farm,
                                    top = tableGrob(t(title.model),
                                                    theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                           base_colour = "black",
                                                                           base_size = 16)))

ggsave(paste0(mydir,"/Observed vs Fitted V2.png"),
       grid.arrange(arrangeGrob(obs.fit.region.farm,
                                top = textGrob(t(resp.variable),
                                               gp = gpar(fontsize = 18, fontface = "bold")))),
       width = 20, height = 20, units = "cm")

#================================================
# Plot residuals vs. random effects ----

figure.title4 <- paste0(mydir,"/Residuals V1.png")

final.plot4.A <- arrangeGrob(plot.norm.resid,
                             plot.qq.resid, 
                             ncol = 2)

final.plot4.B <- arrangeGrob(box.day,
                             box.region,
                             box.farm, 
                             ncol = 3)

final.plot4 <- grid.arrange(final.plot4.A,
                            final.plot4.B,
                            box.region.farm, 
                            nrow = 3,
                            top = tableGrob(t(title.model),
                                            theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                   base_colour = "black",
                                                                   base_size = 16)))

ggsave(figure.title4,
       grid.arrange(arrangeGrob(final.plot4,
                                top = textGrob(t(resp.variable),
                                               gp = gpar(fontsize = 18, fontface = "bold")))),
       width = 25, height = 32, units = "cm")

#================================================
# Plot residuals and collinearity ----

figure.title5 <- paste0(mydir,"/Residuals V2.png")

final.plot5 <- grid.arrange(plot.res.fit,
                            plot.binned,
                            reqq.plot.region,
                            reqq.plot.farm.region,
                            plot.coll,
                            ncol = 2,
                            top = tableGrob(t(title.model),
                                            theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                   base_colour = "black",
                                                                   base_size = 16)))

ggsave(figure.title5,
       grid.arrange(arrangeGrob(final.plot5,
                                top = textGrob(t(resp.variable),
                                               gp = gpar(fontsize = 18, fontface = "bold")))),
       width = 25, height = 32, units = "cm")

#================================================
# Plot data - Model output ----
# Panels: Region ~ Farm

# Raw data vs. day

g1 <- ggplot(resilience_mod, 
             aes(x = day, y = inverts_BC_all)) +
  geom_point(aes(fill = rain_out, color = rain_out),
             size = 2, shape = 21, alpha = 0.2) +
  facet_grid(Region ~ Farm) +
  geom_smooth(aes(y = Fitted),
              method = "lm", se = F, color = "black") +
  theme_minimal() +
  theme(strip.background = element_rect(),
        panel.spacing = unit(1, "lines"),
        legend.position = "bottom",
        axis.title.x = element_text(size = 12, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 12, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 14, face = "bold")) +
  scale_fill_gradient2(name = rain.lab, midpoint = midpoint, 
                       low = "red", mid = "lightblue", high = "blue") +
  scale_color_gradient2(name = rain.lab, midpoint = midpoint, 
                        low = "red", mid = "lightblue", high = "blue") +
  scale_x_continuous(limits = c(0,64), breaks = c(0,8,16,32,64)) +
  labs(x = "Days since the end of the simulated drought",
       y = ylab, 
       title = "Raw data",
       subtitle = "lm(fitted data ~ day)") +
  guides(fill = guide_colorbar(title.position = "top", barwidth = unit(6,"cm")),
         color = guide_colorbar(title.position = "top", barwidth = unit(6,"cm")))

## Panels: Region

# Raw data

g2 <- ggplot(resilience_mod, 
             aes(x = day, y = inverts_BC_all)) +
  geom_point(aes(fill = rain_out, color = rain_out),
             size = 3, shape = 21, alpha = 0.2) +
  facet_wrap(~ Region) +
  geom_smooth(aes(y = Fitted, linetype = Farm),
              method = "lm", se = F, size = 0.6, color = "black") +
  theme_minimal() +
  theme(strip.background = element_rect(),
        panel.spacing = unit(1, "lines"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", margin = margin(15,0,0,0)),
        axis.title.y = element_text(size = 12, face = "bold", margin = margin(0,15,0,0)),
        plot.title = element_text(size = 14, face = "bold")) +
  scale_linetype_manual(values = c(2,3,4,5,6)) +
  scale_fill_gradient2(name = rain.lab, midpoint = midpoint, 
                       low = "red", mid = "lightblue", high = "blue") +
  scale_color_gradient2(name = rain.lab, midpoint = midpoint, 
                        low = "red", mid = "lightblue", high = "blue") +
  scale_x_continuous(limits = c(0,64), breaks = c(0,8,16,32,64)) +
  labs(x = "Days since the end of the simulated drought",
       y = ylab, 
       title = "Raw data",
       subtitle = "lm(fitted data ~ day)") +
  guides(fill = guide_colorbar(title.position = "top", barwidth = unit(6,"cm")),
         color = guide_colorbar(title.position = "top", barwidth = unit(6,"cm")),
         linetype = guide_legend(title.position = "top"))

# Both plots together

figure.title6 <- paste0(mydir,"/Model output.png")

final.plot6 <- grid.arrange(g1,g2,
                            ncol = 2,
                            top = tableGrob(t(title.model),
                                            theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                   base_colour = "black",
                                                                   base_size = 16)))


ggsave(figure.title6,
       grid.arrange(arrangeGrob(final.plot6,
                                top = textGrob(t(resp.variable),
                                               gp = gpar(fontsize = 18, fontface = "bold")))),
       width = 30, height = 19, units = "cm")

#================================================
# Statistics ----
 
AIC(lmer_full_mod)                                  ## -272.4079
model_performance(lmer_full_mod)$R2_marginal        ## 0.1945697
model_performance(lmer_full_mod)$R2_conditional     ## 0.2225412

write.table(aov_mod, paste0(mydir.data, "/LMER Anova output - ",resp.variable," vs Day and Rain Out 2.csv"), 
            sep = ",", row.names = F)

write.table(fixef_out, paste0(mydir.data, "/LMER Model output fixef - ",resp.variable," vs Day and Rain Out 2.csv"), 
            sep = ",", row.names = F)

write.table(random_out.1, paste0(mydir.data, "/LMER Model output random Farm Region - ",resp.variable," vs Day and Rain Out 2.csv"), 
            sep = ",", row.names = F)

write.table(random_out.2, paste0(mydir.data, "/LMER Model output random Region - ",resp.variable," vs Day and Rain Out 2.csv"), 
            sep = ",", row.names = F)
