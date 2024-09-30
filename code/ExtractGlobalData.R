### Extract global data for ice-free land area for representative analysis ###
### Sophie von Fromm ###
### 2023-10-27 ###

library(tidyverse)
library(terra)
library(sf)
library(ncdf4)

## Load global soil map to define area of interest (aoi)
GSM_dir <- "C:/Users/f0076db/Documents/PhD/AfSIS_GlobalData/USDA_GlobalSoilMap/so2015v2.tif"
GSM_raster <- terra::rast(GSM_dir)

plot(GSM_raster)
summary(GSM_raster)
print(GSM_raster)

#Remove values 0 to 4: water bodies, shifting sands, rocky land, ice/glacier, salt
GSM_raster_NA <- terra::classify(GSM_raster, 
                                 cbind(-Inf, 4, NA))
print(GSM_raster_NA)
plot(GSM_raster_NA)

crs(GSM_raster_NA) <- "+proj=longlat +datum=WGS84 +no_defs"

# crop out Antarctica
aoi <- extent(c(xmin = -180, xmax = 180, ymin = -58, ymax = 90))

GSM_raster_crop <- terra::crop(GSM_raster_NA, aoi)

plot(GSM_raster_crop)

## Load global data of interest
# MAT
MAT_dir <- "C:/Users/f0076db/Documents/PhD/AfSIS_GlobalData/wc2.0_30s_bio/wc2.0_bio_30s_01.tif"
MAT_global <- terra::rast(MAT_dir)
MAT_global_crop <- terra::crop(MAT_global, aoi)

plot(MAT_global_crop)

# MAP
MAP_dir <- "C:/Users/f0076db/Documents/PhD/AfSIS_GlobalData/wc2.0_30s_bio/wc2.0_bio_30s_12.tif"
MAP_global <- terra::rast(MAP_dir)
MAP_global_crop <- terra::crop(MAP_global, aoi)

# PET
PET_dir <- "D:/Sophie/PhD/AfSIS_GlobalData/global-et0_annual.tif/et0_yr/et0_yr.tif"
PET_global <- terra::rast(PET_dir)
PET_global_crop <- terra::crop(PET_global, aoi)

PET_resample <- terra::resample(PET_global_crop, MAT_global_crop)

climate_raster <- rast(list(MAT = MAT_global_crop, 
                            MAP = MAP_global_crop,
                            PET = PET_resample))

print(climate_raster)

# GPP
read_nc_fun <- function(path, pattern = "*.nc"){
  list.files(path, pattern, full.names = TRUE) %>% 
    map(~rast(.))
}

GPP_global <- read_nc_fun(path = "C:/Users/f0076db/Documents/PhD/AfSIS_GlobalData/Fluxcom_GPP")
print(GPP_global)

#calculate yearly means for each year
Global_calc <- map(GPP_global, ~calc(., fun = mean))

#Calculate long-term mean
GPP_mean <- (Global_calc[[1]]+Global_calc[[2]]+Global_calc[[3]]+
               Global_calc[[4]]+Global_calc[[5]]+Global_calc[[6]]+
               Global_calc[[7]]+Global_calc[[8]]+Global_calc[[9]]+
               Global_calc[[10]]+Global_calc[[11]]+Global_calc[[12]])/12

GPP_resample <- raster::resample(GPP_mean, climate_raster)

climate_GPP_raster <- stack(list(climate_raster, GPP = GPP_resample))


## Crop and mask global data to aoi
GSM_resample <- raster::resample(GSM_raster_crop, climate_GPP_raster)

climate_GPP_mask <- raster::mask(climate_GPP_raster, GSM_resample)

### Extract global values
# set.seed(42)
set.seed(123)
climate_GPP_subset <- sampleRandom(climate_GPP_mask, size = 1000000, xy = TRUE) %>% 
  as.data.frame()

skimr::skim_without_charts(climate_GPP_subset)

summary(climate_GPP_subset$MAT)

write_csv(climate_GPP_subset, file = paste0("./data/GlobalClimateGPPSubset_1000000_", 
                                            Sys.Date(), ".csv"))

