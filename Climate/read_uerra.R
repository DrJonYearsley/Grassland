# Import uerra data from a GRIB file
#
# This script imports GRIB data
#
# To do: 
# select the pixels that correspond to the experimental sites
#
# ****************************

setwd('~/WorkFiles/MEGA/Projects/GrasslandResilience/MERA_Data/')

rm(list=ls())
library(raster)
library(rgdal)

# Import UERRA temperature data for August 2017. 
# There are four data points per day (i.e. 124 data points per pixel)
temp = stack('uerra_5km_temp_2m_Aug_2017.grib')

# Import outline of Ireland and convert to UERRA CRS 
# (Lambert Conformal Conic)
ir = readOGR(dsn='../country.shp', layer='country')
ir2 = spTransform(ir, CRSobj = crs(temp[[1]]))

# Crop raster stack to Ireland
temp_ir = crop(temp, extent(ir2))

# Project UERRA data onto WGS84 (i.e. Lat Long) and Irish Grid (TM75)
temp_ir_wgs = projectRaster(temp_ir[[1]], crs=CRS("+init=epsg:4326"))
temp_ir_TM75 = projectRaster(temp_ir[[1]], crs=CRS("+init=epsg:29903"))

# Plot data in WGS84 CRS
plot(temp_ir_wgs)

