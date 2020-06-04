# Look at temporal autocorrleation in soil moisture versus rainfall
#
# Jon Yearsley  Jon.Yearsley@cd.ie
# 4th June 2020
# **************************************************

rm(list=ls())
library(forecast)

# Load soil moisture data
sm = read.table('~/Dropbox/Grassland/Data/Climate/04_MERA_rainfall_and_temp/Rainfall and soil moisture per site 2017.csv', header=T, sep=',')

# Read in rainfall data

sm$doy = as.numeric(format(as.Date(sm$Date),"%j"))

tmp = subset(sm,  Plot=='C2' & Site=='D1')

pacf(ts(tmp$Avg.Daily.Moisture))

m = auto.arima(ts(tmp$Avg.Daily.Moisture), xreg=tmp$Rainfall, 
               d=0,D=0, max.p=10,max.q=10, max.P=10, max.Q=10, 
               start.p = 9, start.q=10,start.P=10, start.Q=10,
               trace=TRUE, max.order=10, stepwise=TRUE)

plot(forecast(m, xreg=tmp$Rainfall))
points(tmp$Avg.Daily.Moisture)

par(mfrow=c(2,1))
plot(tmp$doy,fitted(m),t='l')
points(tmp$doy,tmp$Avg.Daily.Moisture,t='p')

plot(tmp$doy,tmp$Rainfall,t='l',col='red')
par(mfrow=c(1,1))

library(nlme)

m.lme = gls(Avg.Daily.Moisture~Rainfall, 
            correlation = corAR1(0.4,form=~doy, fixed=FALSE),
            data=tmp,
            na.action = 'na.omit')

summary(m.lme)

tmp$fit = predict(m.lme) + residuals(m.lme)

phi = diag(1,nrow(tmp))
diag(phi[-1,]) <- 0.8538


tmp$fit2 = solve(phi,residuals(m.lme)) + tmp$fit


library(ggplot2)

ggplot(data=tmp, aes(x=doy, y=Avg.Daily.Moisture)) + 
  geom_point() +
  geom_path(aes(y=fit))


qqnorm(residuals(m.lme))



