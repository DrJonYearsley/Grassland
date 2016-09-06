# Script to read in processed MODIS data for pasture in Ireland
# and analyse the variance and mean within hectads
#
# Jon Yearsley Aug 2016

library(raster)
library(foreach)
library(doParallel)
#library(doMC)


# Register cluster with 2 nodes
cl<-makeCluster(2)
registerDoParallel(cl)


rm(list=ls())
setwd("/media/jon/3TB/jon/PeopleStuff/Resilience_MarkJack")

# Scaling for aggregation. Raw data is at 250mx250m
scaleFactor=40  # Scale to 10km
filename_out = 'test.RData'

# Load MODIS data
evi.files = list.files(path='./Data/MODIS',pattern='EVI_pasture_([0-9]{4})_([0-9]{2})_([0-9]{2}).grd',full.names=T)
ndvi.files = list.files(path='./Data/MODIS/',pattern='NDVI_pasture_([0-9]{4})_([0-9]{2})_([0-9]{2}).grd',full.names=T)

dateInd = c(regexpr(evi.files[1],pattern='([0-9]{4})_([0-9]{2})_([0-9]{2})',fixed=F))
evi.dates = strptime(substr(evi.files,start=dateInd, stop=dateInd+9), format="%Y_%m_%d")

date.order=order(evi.dates)

modis.dates = evi.dates[date.order]
julian.date = julian(modis.dates)
day = format(modis.dates, '%d')
month = format(modis.dates, '%m')
year = format(modis.dates, '%Y')


# # Read in Ireland coastline
# ie = readOGR(dsn='Data', layer='country')
# ie.grid = spTransform(ie, CRS=CRS("+init=epsg:29903"))   # Transform to Irish Grid TM75

# Get the coords of the new aggregated data
tmp = aggregate(raster(evi.files[1]), fact=scaleFactor, fun=mean, na.rm=T)
coord = coordinates(tmp)
eastings = t(array(coord[,1],dim=rev(dim(tmp)[1:2])))
northings = t(array(coord[,2],dim=rev(dim(tmp)[1:2])))
projection = proj4string(tmp)
rasterExtent = extent(tmp)


# evi.mean.list <- foreach (f = icount(3), .packages='raster', .inorder=T) %dopar% {
#   # Read in processed MODIS data (subsetted to pasture and CRS set to Irish TM75, rounded to nearest hectad)
#   as.matrix(aggregate(raster(evi.files[f]), fact=scaleFactor, fun=mean, na.rm=T))
# }

# Aggregate EVI data
evi.mean.list <- foreach (f = icount(length(evi.files)), .packages='raster', .inorder=T) %dopar% {
  # Read in processed MODIS data (subsetted to pasture and CRS set to Irish TM75, rounded to nearest hectad)
  as.matrix(aggregate(raster(evi.files[f]), fact=scaleFactor, fun=mean, na.rm=T))
}
evi.sd.list <- foreach (f = icount(length(evi.files)), .packages='raster', .inorder=T) %dopar% {
  # Read in processed MODIS data (subsetted to pasture and CRS set to Irish TM75, rounded to nearest hectad)
  as.matrix(aggregate(raster(evi.files[f]), fact=scaleFactor, fun=sd, na.rm=T))
}
evi.ncell.list <- foreach (f = icount(length(evi.files)), .packages='raster', .inorder=T) %dopar% {
  # Read in processed MODIS data (subsetted to pasture and CRS set to Irish TM75, rounded to nearest hectad)
  as.matrix(aggregate(raster(evi.files[f]), fact=scaleFactor, fun=function(x, na.rm=T){sum(!is.na(x), na.rm=na.rm)}))
}

# Aggregate NDVI data
ndvi.mean.list <- foreach (f = icount(length(ndvi.files)), .packages='raster', .inorder=T) %dopar% {
  # Read in processed MODIS data (subsetted to pasture and CRS set to Irish TM75, rounded to nearest hectad)
  as.matrix(aggregate(raster(ndvi.files[f]), fact=scaleFactor, fun=mean, na.rm=T))
}
ndvi.sd.list <- foreach (f = icount(length(ndvi.files)), .packages='raster', .inorder=T) %dopar% {
  # Read in processed MODIS data (subsetted to pasture and CRS set to Irish TM75, rounded to nearest hectad)
  as.matrix(aggregate(raster(ndvi.files[f]), fact=scaleFactor, fun=sd, na.rm=T))
}
ndvi.ncell.list <- foreach (f = icount(length(ndvi.files)), .packages='raster', .inorder=T) %dopar% {
  # Read in processed MODIS data (subsetted to pasture and CRS set to Irish TM75, rounded to nearest hectad)
  as.matrix(aggregate(raster(ndvi.files[f]), fact=scaleFactor, fun=function(x, na.rm=T){sum(!is.na(x), na.rm=na.rm)}))
}

# Put the results in a big matrix
dimMatrix=c(dim(evi.mean.list[[1]]), length(modis.dates))
evi.mean = array(NA, dim=dimMatrix)
evi.sd = array(NA, dim=dimMatrix)
evi.ncell = array(NA, dim=dimMatrix)
for (f in 1:length(evi.files)) {
  evi.mean[,,f] = evi.mean.list[[date.order[f]]]
  evi.sd[,,f] = evi.sd.list[[date.order[f]]]
  evi.ncell[,,f] = evi.ncell.list[[date.order[f]]]
}

ndvi.mean = array(NA, dim=dimMatrix)
ndvi.sd = array(NA, dim=dimMatrix)
ndvi.ncell = array(NA, dim=dimMatrix)
for (f in 1:length(evi.files)) {
  ndvi.mean[,,f] = ndvi.mean.list[[date.order[f]]]
  ndvi.sd[,,f] = ndvi.sd.list[[date.order[f]]]
  ndvi.ncell[,,f] = ndvi.ncell.list[[date.order[f]]]
}


save(evi.mean, evi.sd, evi.ncell, ndvi.mean, ndvi.sd, ndvi.ncell, eastings, northings, projection, modis.dates, julian.date, day,month,year, rasterExtent,file=filename_out )


stopCluster(cl)
