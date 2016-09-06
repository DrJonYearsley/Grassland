# Import the processed MODIS data for pastures in Ireland (aggregated at the hectad scale)
# Look at the relationship between spatial mean and spatial standard deviation within a hectad
#
# MODIS data are from the Terra platform,(version 5)
# Only MODIS data with a 'useable' quality score have been used. 
#
# The agregated data were created by a MODIS raster with a spatial resolution of 250m. 
# Each hectad therefore has a maximum of 1600 (40x40) cells that could be aggregated over. 
# The actual number of 'pasture' cells 
#
# Jon Yearsley   Aug 2016

rm(list=ls())
library(raster)
library(sp)
library(rgdal)

#setwd('~/Dropbox/SFIDEL_Finalise/MODIS_prelim/')
setwd('~/WorkFiles/PeopleStuff/Resilience_MarkJack/')

# Import the MODIS data
load('modis_1km_par.RData')

# Look at standard deviation as a function of number of pasture cells
plot(x=evi.ncell, y=evi.sd, xlab='Number of cells in a hectad', ylab='EVI standard deviation', pch='.')

# Mean vs SD relationship for each year with ncell>500

yearList=unique(year)
for(yr in yearList) {
  x=evi.mean[,,year==yr]
  y=evi.sd[,,year==yr]
  ind = evi.ncell[,,year==yr]>100
  
  plot(x[ind], y[ind], xlab='EVI Mean', ylab='EVI s.d.', pch='.',main=paste('Year', yr))
#  plot(x[ind], y[ind]/x[ind], xlab='EVI Mean', ylab='EVI Coef. Var.', pch=19, main=paste('Year', yr))
}


season = list(c(2:4),c(5:7),c(8:10),c(11,12,1))
seasonName = c('Feb-Apr','May-Jul','Aug-Oct','Nov-Jan')
yearList=c(2005:2008)
for(s in 1:4) {
  ind = as.numeric(month)%in%season[[s]] & as.numeric(year)%in%yearList
  x=apply(evi.mean[,,ind],c(1,2),mean, na.rm=T)
  y=apply(evi.mean[,,ind],c(1,2),sd, na.rm=T)
  
  indPlot = apply(evi.ncell[,,ind],c(1,2),min,na.rm=T)>100
  plot(x[indPlot], y[indPlot], xlab='EVI temporal mean of hectad mean', ylab='EVI temporal s.d. of hectad mean', pch=19, main=paste('Season', seasonName[s],', Years',min(yearList),'-',max(yearList)))
#  plot(x[ind], y[ind]/x[ind], xlab='EVI Mean', ylab='EVI Coef. Var.', pch=19, main=paste('Season', seasonName[s]))
}

yearList=c(2009:2014)
corResults = array(NA,dim=c(1,length(yearList)))
seasonID=4
for(yr in 1:length(yearList)) {
  ind = as.numeric(month)%in%season[[seasonID]] & as.numeric(year)%in%yearList[yr]
  
  # Calculate temporal variation within the season and months specified by ind
  x=apply(evi.mean[,,ind],c(1,2),mean, na.rm=T)
  y=apply(evi.mean[,,ind],c(1,2),sd, na.rm=T)

  # Only plot data from hectads with more than 10 valid MODIS pixels
  indPlot = apply(evi.ncell[,,ind],c(1,2),min,na.rm=T)>10
  
  
  corResults[yr] = cor(c(x[indPlot]),c(y[indPlot]),
                       use='complete.obs',method='spearman')
  
  plot(x[indPlot], y[indPlot],
       xlab='EVI temporal mean of hectad mean',
       ylab='EVI temporal s.d. of hectad mean', pch='.',
       main=paste('Season', seasonName[seasonID],', Years',yearList[yr]))
  #  plot(x[ind], y[ind]/x[ind], xlab='EVI Mean', ylab='EVI Coef. Var.', pch=19, main=paste('Season', seasonName[s]))
}



# Create a spatial plot of the evi variability
ind=as.numeric(year)%in%c(2009:2014) & as.numeric(month)%in%season[[1]]
x = raster(apply(evi.mean[,,ind],c(1,2),sd,na.rm=T), crs=projection)
extent(x) <- rasterExtent
plot(x,xlab='Eastings',ylab='Northings',main='Temporal SD of hectad mean EVI')
# Add on the coeatline of Ireland
ie = readOGR(dsn='./Data', layer='country')
ie.grid = spTransform(ie, CRS("+init=epsg:29903"))   # Transform to Irish Grid TM75
plot(ie.grid, add=T)

# Create a temporal plot of evi in a few cells
ind = as.numeric(year)%in%c(2009:2014)



x = evi.mean[230,140,]
plot(modis.dates[ind],x[ind],type='l',xlab='Date (years)',ylab='Mean EVI within the hectad')
lines(modis.dates[ind],evi.mean[200,140,ind],col='blue')
lines(modis.dates[ind],evi.mean[26,32,ind],col='darkred')
lines(modis.dates[ind],evi.mean[30,29,ind],col='darkgreen')
lines(modis.dates[ind],evi.mean[10,19,ind],col='darkorange')
