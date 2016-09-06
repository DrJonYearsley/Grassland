# Script to read in MODIS data, crop to Ireland, rescale NDVI and EVI 
# and remove poor quality pixels
#
# Jon Yearsley Aug 2016

library(rgdal)
library(raster)
library(gdalUtils)

rm(list=ls())
setwd("/media/jon/3TB/jon/PeopleStuff/Resilience_MarkJack")

# Load CORINE 2012 raster data and crop to Ireland
corine = raster('Data/CORINE_2012_raster_Europe/g100_clc12_V18_5.tif')
corine2 = crop(corine, extent(2891919, 3366575 , 3230104 , 3808712))
corine.crs = crs(corine)

# Read in Ireland coastline
ie = readOGR(dsn='Data', layer='country')
ie.lc = spTransform(ie, CRS=crs(corine2))

corine.ie = mask(corine2, ie.lc)

writeRaster(corine.ie,file='./Data/CORINE_IE',format='raster',overwrite=TRUE)
