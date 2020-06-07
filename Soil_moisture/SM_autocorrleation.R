# Soil Moisture time series
#
# Look at correlation between changes in soil moisture and 
# lagged rainfall
# The fit some time series models. AR1 looks appropriate.
# Make sure that rainfall is appropriately lagged so that 
# difference in soil moisture from day t-1 and day t corresponds 
# to rainfall on day t-1 (not day t)
#
# Jon Yearsley  Jon.Yearsley@cd.ie
# 5th June 2020
# **************************************************

rm(list=ls())
library(forecast)
library(ggplot2)


# Load soil moisture data
sm = read.table('~/Dropbox/Grassland/Data/Climate/04_MERA_rainfall_and_temp/Rainfall and soil moisture per site 2017.csv', header=T, sep=',')
# Create day of year
sm$doy = as.numeric(format(as.Date(sm$Date),"%j"))

# Compute lagged rainfall and soil moisture by 1 day, 
# so that rainfall corresponds 
create_df = TRUE
for (s in unique(sm$Site)) {
  for (p in unique(sm$Plot)) {
    tmp = subset(sm, Site==s & Plot==p)
    if (nrow(tmp)!=0) {
      tmp$Rainfall_lagged = dplyr::lead(tmp$Rainfall, 1)
      tmp$Moisture_1day = dplyr::lag(tmp$Avg.Daily.Moisture, 1)
      tmp$diff_Moisture = c(NA,diff(tmp$Avg.Daily.Moisture))
      
      if (create_df) {
        sm_new = tmp
        create_df = FALSE
      } else {
        sm_new = rbind(sm_new, tmp)
      }
    }
  }
}


# Correlation between change in soil moisture and rainfall
cors = data.frame(Region=NULL, 
                  Site=NULL, 
                  treatment=NULL, 
                  lag = NULL, 
                  cor=NULL)
for (s in unique(sm$Site)) {
  for (t in unique(sm$treatment)) {
    
    tmp = subset(sm, Site==s & treatment==t)
    
    cor.df = data.frame(Region=unique(tmp$Region),
                        Site=s,
                        treatment=t,
                        lag = seq(-5,5),
                        cor = NA)
    
    for (l in 1:nrow(cor.df)) {
      if (cor.df$lag[l]<0) {
        cor.df$cor[l] = cor(dplyr::lead(tmp$Rainfall,abs(cor.df$lag[l])),
                             c(diff(tmp$Avg.Daily.Moisture),NA),
                             use = 'complete.obs')
      } else {
        cor.df$cor[l] = cor(dplyr::lag(tmp$Rainfall,cor.df$lag[l]),
                            c(diff(tmp$Avg.Daily.Moisture),NA),
                             use = 'complete.obs')
      }
    }
    
    cors = rbind(cors,cor.df)
    
  }
}

# Plot correlation between lagged rainfall and difference in soil moisture
ggplot(data=cors, aes(x=lag, y=cor, colour=treatment)) +
  geom_point() + 
  theme_bw() + 
  ggtitle('Correlaton of change in soil moisture versus lagged rainfall')




# **************************************************
# Do some time series modelling (AR1 seems appropriate)

library(nlme)

# Look at one plot to start with
tmp = subset(sm_new, Site=='C2' & Plot=='C3')
#tmp$Rainfall_lag = dplyr::lag(tmp$Rainfall, n=1)

# Fit an AR1 model using arima function
m.ar1 = arima(tmp$Avg.Daily.Moisture, 
              order = c(1,0,0),
              xreg=tmp$Rainfall_lagged)

tmp$fit.ar = fitted(m.ar1)

# Fit model without rainfall being lagged
m.ar1b = arima(tmp$Avg.Daily.Moisture, 
               order = c(1,0,0),
               xreg = tmp$Rainfall)
tmp$fit.arb = fitted(m.ar1b)

m.ar1c = arima(tmp$Avg.Daily.Moisture, 
              order = c(1,0,0))
tmp$fit.arc = fitted(m.ar1c)

summary(m.ar1)

# Plot model fit versus data for three models. The strong autocorrelation dominates
ggplot(data=tmp, aes(x=doy, y=Avg.Daily.Moisture)) + 
  geom_point() +
  geom_path(aes(y=fit.ar), colour='red') + 
  geom_path(aes(y=fit.arb), colour='orange') + 
  geom_path(aes(y=fit.arc), colour='blue') + 
  ggtitle('arima (red), arima (rainfall not lagged, orange), no covariate (blue) and data (black)') +
  theme_bw()


# ******************************
# Recreate this model using gls from nlme
m.gls1 = gls(Avg.Daily.Moisture~Rainfall_lagged, 
            correlation = corAR1(form=~doy, fixed=FALSE),
            data=tmp,
            na.action = 'na.exclude')

summary(m.gls1)

# Manually created fitted values with AR1 process added
b0 = coef(m.gls1)[1]
b1 = coef(m.gls1)[2]
p1 = intervals(m.gls1)$corStruct[2]  # AR1 parameter

fit = b0 + b1*tmp$Rainfall_lagged
phi = diag(1,length(fit)-1)
# diag(phi[-1,]) <- -p1
# fit2 = solve(phi, fit[-1]-b0) + b0

# Add on AR1 part (lagged effect of soil moisture)
tmp$fit.gls = fit + c(0, (tmp$Avg.Daily.Moisture[-length(fit)]-b0)*p1)


# Plot nlme model predictions next to data and arima predictions
ggplot(data=tmp, aes(x=doy, y=Avg.Daily.Moisture)) + 
  geom_point() +
  geom_path(aes(y=fit.gls), colour='blue')+
  geom_path(aes(y=fit.ar), colour='red') + 
  ggtitle('arima (red), gls (blue) and data (black)') +
  theme_bw()

#######################################################
# Fit model for all control plots
m.gls2 = gls(Avg.Daily.Moisture~Rainfall+factor(Site)+factor(Plot), 
            correlation = corAR1(form=~doy|Site/Plot, fixed=FALSE),
            data=subset(sm_new,treatment=='Control'),
            subset=treatment=='Control',
            na.action = 'na.exclude')

summary(m.gls2)
# Strong effect of rainfall


# ***********************************
# Fit model looking at differences in soil moisture ----
# Just look at one site
library(mgcv)

# Fit AR1 model
m.ar = gls(diff_Moisture~1+I(log(1+Rainfall)), 
               correlation = corAR1(form=~doy, fixed=FALSE),
               data=tmp,
               na.action = 'na.exclude')

# Fit GAM model
m.gam = gam(diff_Moisture~s(I(log(Rainfall+1))), 
             data=tmp,
             na.action = 'na.exclude')


summary(m.gam)

# Fit model without rainfall
m0 = update(m.gam,.~1+I(log(1+Rainfall)))


tmp$fit.gam = fitted(m.gam)
tmp$fit.ar = fitted(m.ar)
tmp$fit.lm = fitted(m0)

# Plot gam model predictions next to data and arima predictions
ggplot(data=tmp, aes(x=doy, y=diff_Moisture)) + 
  geom_point() +
  geom_path(aes(y=fit.gam), colour='blue') +
  geom_path(aes(y=fit.lm), colour='red') + 
  geom_path(aes(y=fit.ar), colour='orange') + 
  geom_linerange(aes(ymin=0,ymax=Rainfall/max(Rainfall)*5), colour='blue', size=2) +
  ggtitle('Linear model (red), gam (blue), AR1 (orange) and data (black)') +
  theme_bw() + 
  ylim(c(-4,10))

#plot(m.gam, residuals=TRUE, pch=20)


# Looks like a linear model performs well

# Model all sites with a linear model and log(Rainfall+1) as covariate

library(lme4)
m.lmfull = lmer(diff_Moisture~1+I(log(Rainfall+1))+(1|Region), 
               data=subset(sm_new,treatment=='Control'),
               subset=treatment=='Control',
               na.action = 'na.exclude')

summary(m.lmfull)


check_model(m.lmfull)
model_performance(m.lmfull)
