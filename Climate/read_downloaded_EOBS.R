# This script reads in the raw EOBS gridded climate data 
# for 1995-2016 and pulls out the spatial region for Ireland 
# (longitude: -5.5 -- -11, latitude: 51 -- 55.4)

# Most recent EOBS data is 0.1 deg resolution (roughly 10km)
# Earlier data is 0.25 deg resolution.

# The data can be downloaded from 
# http://www.ecad.eu/download/ensembles/downloadchunks.php
# 
# Variables in this data set are:
# tg = daily mean temperature (units = deg Celsius)
# tn = daily minimum temperature (units = deg Celsius)
# tx = daily maximum temperature (units = deg Celsius)
# rr = daily precipitation sum (units = mm)
# pp = daily averaged sea level pressure  (units = hecto Pascals)
# elevation = mean elevation (m)
#
#  Jon Yearsley   6th June 2017 (jon.yearsley@ucd.ie)

# Aug 2019: Updated for new EOBS data
###########################################


library(ncdf4)

rm(list=ls())

# Names of climate variables to import (if data exists)
vars = data.frame(name=c('tg','tn','tx','rr','pp'),
                 desc=c('Mean daily temp',
                        'Min daily temp',
                        'Max daily temp',
                        'Daily rainfall (mm)',
                        'Average daily pressure'),
                 available=TRUE,
                 stringsAsFactors = FALSE)

eobsDir = '~/WorkFiles/Data/Climate/eobs_0.1deg' # Directory containing data
outfile = 'eobs_0.1deg.RData'
outDir = '~/WorkFiles/Data/Climate'

files = list.files(path = eobsDir,
                   pattern = paste(vars$name,collapse='|'),
                   full.names = TRUE)

for (f in 1:length(vars$name)) {
  if (!any(grepl(pattern = vars$name[f],files))) {  
  vars$available[f] = FALSE
  }
}

# Loop around each variable. Import the variable and the accompanying standard error
inds = which(vars$available)

for (f in 1:length(inds)) {
    fileInd2 = grep(pattern = vars$name[inds[f]],files)  
    fileInd1 = grep(pattern = paste0(vars$name[inds[f]],'[[:print:]]*(spread|stderr)'),files)  
    nc.se = nc_open(files[fileInd1])
    nc = nc_open(files[fileInd2[fileInd2!=fileInd1]])
    
    # Cut out long-lat of Ireland
    indLong = which(nc$dim$longitude$val > -11 & nc$dim$longitude$val < -5.5)
    indLat = which(nc$dim$latitude$val > 51 & nc$dim$latitude$val < 55.4) # Use max latitude of 55.3 to remove Scotland
    d = ncvar_get(nc,varid=vars$name[inds[f]],start=c(min(indLong),min(indLat),1), count=c(length(indLong), length(indLat), -1), verbose = T, collapse_degen=FALSE)
    
    # Cut out long-lat of Ireland
    indLong2 = which(nc.se$dim$longitude$val > -11 & nc.se$dim$longitude$val < -5.5)
    indLat2 = which(nc.se$dim$latitude$val > 51 & nc.se$dim$latitude$val < 55.4) # Use max latitude of 55.3 to remove Scotland
    d.se = ncvar_get(nc.se,varid=vars$name[inds[f]],start=c(min(indLong2),min(indLat2),1), count=c(length(indLong2), length(indLat2), -1), verbose = T, collapse_degen=FALSE)
    
    
    ############
    if (f==1) {
      # Use first file to define the dates, longitudes and latitudes
      longitude = nc$dim$longitude$val[indLong]
      latitude = nc$dim$latitude$val[indLat]
      date = as.Date(nc$dim$time$val,  origin=as.Date('1950-01-01'))  # Time units are number of days since 01/01/1950
      
      n1 = length(longitude)
      n2 = length(latitude)
      n3 = length(date)
      var_all = array(NA, dim=c(n1*n2*n3, 2*length(inds)))
    }
    
    nc_close(nc)
    nc_close(nc.se)
    
    # Put data into one large array
    var_all[,f] = d
    var_all[,length(inds)+f] = d.se
}
# Finally import elevation data
elev_file = list.files(path = eobsDir,
                       pattern = 'elev',
                       full.names = TRUE)

nc = nc_open(elev_file)
indLong = which(nc$dim$longitude$val > -11 & nc$dim$longitude$val < -5.5)
indLat = which(nc$dim$latitude$val > 51 & nc$dim$latitude$val < 55.4) # Use max latitude of 55.3 to remove Scotland
d.elev = ncvar_get(nc,varid='elevation',start=c(min(indLong),min(indLat)), count=c(length(indLong), length(indLat)), verbose = T, collapse_degen=FALSE)

# Create a big data frame with all the data
df = data.frame(longitude = rep(longitude, times=n2*n3), 
                latitude = rep(latitude, times=n3, each=n1),
                date = rep(date, each=n1*n2), 
                elevation=rep(d.elev,times=n3),
                var=var_all)
names(df) <- c('longitude','latitude','date','elevation',
               as.character(vars$name[vars$available]),
               paste(as.character(vars$name[vars$available]),'.se',sep=''))

# Remove grid squares with no climate data
eobs_data = subset(df, apply(df[1:200,-c(1:4)],MARGIN=1,FUN=function(x){any(is.finite(x))}))

# Save the combined data
save(eobs_data, file=file.path(outDir,outfile))


########### Some quick validation plots
# Plot one variable as a check
image(longitude, latitude, d.elev)  #Plot elevation
image(longitude, latitude, apply(d, MARGIN=c(1,2), FUN=mean, na.rm=T))  # Plot a climate variable

# Plot a time series
d.sub = subset(eobs_data, longitude < -9.62 & longitude > -9.63 & latitude>51.62 & latitude<51.63)
library(ggplot2)
ggplot(data=d.sub, aes(x=date,y=tg)) +  geom_line()
