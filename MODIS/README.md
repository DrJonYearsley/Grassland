Management and analysis of the NDVI and EVI products from MODIS

**process_CORINE.R** reads in the CORINE ratser file for Europe and crops the data to be just Ireland

**process_MODIS.R** this is the main script for processing MODIS data. It performs the following steps:

1. read in downloaded MODIS data (the MOD13Q1 product in geoTIFF format), 
2. crop MODIS data to Ireland, 
3. rescale NDVI and EVI to be on a scale 0-1
4. remove poor quality pixels
5. reproject data onto CORINE raster coordinate system
6. select pixels correpsonding to pasture (corine code 18)
7. reproject onto Irish grid (TM75) coordinate system
8. rescale extent so that grid aligns with hectads (to allow aggregation)
9. saves the processed raster data for NDVI, EVI and the quality data

The final data that is saved to file will have a spatial resolution of approximately 250m

**analyse_pasture_parallel.R** reads in the processed MODIS data (EVI and NDVI) for Ireland and calculates summary statistics (e.g. mean, variance, number of pixels) at a specified spatial scale. The script is set up to run on a cluster if necessary.
