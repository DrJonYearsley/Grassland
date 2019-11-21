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
ir2 = spTransform(ir, CRSobj = crs(d[[1]]))

# Crop raster stack to Ireland
temp_ir = crop(temp, extent(ir2))

# Project UERRA data onto WGS84 (i.e. Lat Long) and Irish Grid (TM75)
temp_ir_wgs = projectRaster(temp_ir[[1]], crs=CRS("+init=epsg:4326"))
temp_ir_TM75 = projectRaster(temp_ir[[1]], crs=CRS("+init=epsg:29903"))

# Plot data in WGS84 CRS
plot(temp_ir_wgs)

########################################################
# Transform coords to Lambert Conformal Conic  (MESCAN data)
# Central meridian: 8
#Standard parallel1: 50
# Standard parallel2: 50
# Latitude of origin: 50
# Earth assumed spherical with radius: 6371229m
# Latitude of first grid point in degrees: 20.292
# Longitude of first grid point in degrees: 342.514
GR_lcc = spTransform(GR, CRS("+proj=lcc +lon_0=8 +lat_0=48 +lat_1=48 +lat_2=48 +R=6371229 +datum=WGS84"))

# Transform coords to Lambert Conformal Conic  (MESCAN data)
# Central meridian: 8
#Standard parallel1: 50
# Standard parallel2: 50
# Latitude of origin: 50
# Earth assumed spherical with radius: 6371229m
# Latitude of first grid point in degrees: 20.292
# Longitude of first grid point in degrees: 342.514
GR_lcc = spTransform(GR, CRS("+proj=lcc +lon_0=8 +lat_0=50 +lat_1=50 +lat_2=50 +R=6371229 +datum=WGS84"))

crs_old = crs(dsub)
crs(dsub) = CRS("+proj=lcc +lon_0=8 +lat_0=48 +lat_1=48 +lat_2=48 +R=6371229 +datum=WGS84")

ir = readOGR(dsn='../country.shp', layer='country')
europe = readOGR(dsn='../Europe/Europe_coastline_poly.shp', layer='Europe_coastline_poly')
spain = readOGR(dsn='../Europe/ESP_adm0.shp', layer='ESP_adm0')

ir2 = spTransform(ir, CRSobj = crs(tmp))
spain2 = spTransform(spain, CRSobj = crs(dsub))

plot(dsub)
plot(spain2, add=T)

plot(gb)
plot(spain, add=T)
plot(spain2, axes=T)







