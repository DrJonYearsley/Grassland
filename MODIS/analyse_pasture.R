# Script to read in processed MODIS data for pasture in Ireland
# and analyse the variance and mean within hectads
#
# Jon Yearsley Aug 2016

library(raster)


rm(list=ls())
#setwd("/media/jon/3TB/jon/PeopleStuff/Resilience_MarkJack")

# Scaling for aggregation. Raw data is at 250mx250m
#scaleFactor=4  # Scale to 1km
scaleFactor=40 # Scale to 10km (hectads)
filename_out = 'modis_10km.RData'

# Load MODIS data
evi.files = list.files(path='./Data/MODIS/',pattern='EVI_pasture_([0-9]{4})_([0-9]{2})_([0-9]{2}).grd',full.names=T)
ndvi.files = list.files(path='./Data/MODIS/',pattern='NDVI_pasture_([0-9]{4})_([0-9]{2})_([0-9]{2}).grd',full.names=T)

evi.dates = strptime(substr(evi.files,start=27, stop=36), format="%Y_%m_%d")
date.order=order(evi.dates)

modis.dates = evi.dates[date.order]
julian.date = julian(modis.dates)
day = format(modis.dates, '%d')
month = format(modis.dates, '%m')
year = format(modis.dates, '%Y')



for (f in 1:length(evi.files)) {
#  print(modis.dates[f])
  # Read in processed MODIS data (subsetted to pasture and CRS set to Irish TM75, rounded to nearest hectad)
  evi = raster(evi.files[date.order[f]])

  evi.hectad.mean =  aggregate(evi, fact=scaleFactor, fun=mean, na.rm=T)
  evi.hectad.sd =  aggregate(evi, fact=scaleFactor, fun=sd, na.rm=T)
  evi.hectad.ncell =  aggregate(evi, fact=scaleFactor, fun=function(x, na.rm=T){sum(!is.na(x), na.rm=na.rm)})

  ndvi = raster(ndvi.files[date.order[f]])

  ndvi.hectad.mean =  aggregate(ndvi, fact=scaleFactor, fun=mean, na.rm=T)
  ndvi.hectad.sd =  aggregate(ndvi, fact=scaleFactor, fun=sd, na.rm=T)
  ndvi.hectad.ncell =  aggregate(ndvi, fact=scaleFactor, fun=function(x, na.rm=T){sum(!is.na(x), na.rm=na.rm)})
  
    
  if (f==1) {
    rasterExtent = extent(evi.hectad.mean)

    evi.mean = array(NA,dim=c(nrow(evi.hectad.mean),ncol(evi.hectad.mean),length(modis.dates)))
    evi.sd = array(NA,dim=c(nrow(evi.hectad.mean),ncol(evi.hectad.mean),length(modis.dates)))
    evi.ncell = array(NA,dim=c(nrow(evi.hectad.mean),ncol(evi.hectad.mean),length(modis.dates)))
    
    ndvi.mean = array(NA,dim=c(nrow(evi.hectad.mean),ncol(evi.hectad.mean),length(modis.dates)))
    ndvi.sd = array(NA,dim=c(nrow(evi.hectad.mean),ncol(evi.hectad.mean),length(modis.dates)))
    ndvi.ncell = array(NA,dim=c(nrow(evi.hectad.mean),ncol(evi.hectad.mean),length(modis.dates)))

    coord = coordinates(evi.hectad.mean)
    eastings = t(array(coord[,1],dim=rev(dim(evi.mean[,,1]))))
    northings = t(array(coord[,2],dim=rev(dim(evi.mean[,,1]))))
    projection = proj4string(evi.hectad.mean)
  }
  
  evi.mean[,,f] = as.matrix(evi.hectad.mean)
  evi.sd[,,f] = as.matrix(evi.hectad.sd)
  evi.ncell[,,f] = as.matrix(evi.hectad.ncell)

  ndvi.mean[,,f] = as.matrix(ndvi.hectad.mean)
  ndvi.sd[,,f] = as.matrix(ndvi.hectad.sd)
  ndvi.ncell[,,f] = as.matrix(ndvi.hectad.ncell)
}

save(evi.mean, evi.sd, evi.ncell, ndvi.mean, ndvi.sd, ndvi.ncell, eastings, northings, projection, modis.dates, julian.date, day,month,year, rasterExtent,file=filename_out )

