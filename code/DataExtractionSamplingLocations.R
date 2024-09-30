### Extract global data for soil sampling locations ###
### Sophie von Fromm ###
### 2024-09-29 ###

library(tidyverse)
library(raster)
library(sf)

## israd
library(ISRaD)

ISRaD_dir <- "H:/PostDoc/ISRaD/ISRaD_data_files"
ISRaD_extra <- ISRaD.getdata(directory = ISRaD_dir,
                             dataset = "full", extra = TRUE,
                             force_download = FALSE)

inc_data <- ISRaD.flatten(ISRaD_extra, 'incubation') %>% 
  drop_na(inc_14c) %>% 
  dplyr::select(pro_long, pro_lat, inc_14c) %>% 
  tibble() %>% 
  mutate(lyr_14c = NA_integer_,
         frc_14c = NA_integer_,
         ist_14c = NA_integer_)

skimr::skim_without_charts(inc_data)
    
lyr_data <- ISRaD.flatten(ISRaD_extra, 'layer') %>% 
  drop_na(lyr_14c) %>% 
  dplyr::select(pro_long, pro_lat, lyr_14c) %>% 
  tibble() %>% 
  mutate(inc_14c = NA_integer_,
         frc_14c = NA_integer_,
         ist_14c = NA_integer_)

skimr::skim_without_charts(lyr_data)

frc_data <- ISRaD.flatten(ISRaD_extra, 'fraction') %>% 
  drop_na(frc_14c) %>% 
  dplyr::select(pro_long, pro_lat, frc_14c) %>% 
  tibble() %>% 
  mutate(inc_14c = NA_integer_,
         lyr_14c = NA_integer_,
         ist_14c = NA_integer_)

skimr::skim_without_charts(frc_data)

ist_data <- ISRaD.flatten(ISRaD_extra, 'interstitial') %>% 
  drop_na(ist_14c) %>% 
  dplyr::select(pro_long, pro_lat, ist_14c) %>% 
  tibble() %>% 
  mutate(inc_14c = NA_integer_,
         lyr_14c = NA_integer_,
         frc_14c = NA_integer_)

skimr::skim_without_charts(ist_data)

israd <- rbind(lyr_data, frc_data, inc_data, ist_data) %>% 
  rename(lat_dec_deg = pro_lat,
         long_dec_deg = pro_long) %>% 
  mutate(source = "israd",
         row_id = row_number(),
         label = "radiocarbon") %>% 
  unite(ID, c(row_id, source), remove = FALSE) %>% 
  dplyr::select(ID, lat_dec_deg, long_dec_deg, source, label)


head(israd)
skimr::skim_without_charts(israd)

israd %>% 
  ggplot(aes(y = lat_dec_deg, x = long_dec_deg)) +
  geom_point() +
  theme_classic()

israd_distinct <- israd %>% 
  distinct(lat_dec_deg, long_dec_deg, source, label)

## wosis
wosis <- read_csv("./data/wosis_profiles.csv") %>% 
  mutate(source = "wosis",
         label = "C stock",
         row_id = row_number()) %>% 
  unite(ID, c(row_id, source), remove = FALSE) %>% 
  dplyr::select(ID, lat_dec_deg, long_dec_deg, source, label)

head(wosis)
skimr::skim_without_charts(wosis)

wosis %>% 
  ggplot(aes(y = lat_dec_deg, x = long_dec_deg)) +
  geom_point() +
  theme_classic()

wosis_distinct <- wosis %>% 
  distinct(lat_dec_deg, long_dec_deg, source, label)

## iscn
iscn <- read_csv("./data/iscn_profiles.csv") %>% 
  mutate(source = "iscn",
         label = "soil properties",
         row_id = row_number()) %>% 
  unite(ID, c(row_id, source), remove = FALSE) %>% 
  dplyr::select(ID, lat_dec_deg, long_dec_deg, source, label)

head(iscn)
skimr::skim_without_charts(iscn)

iscn %>% 
  ggplot(aes(y = lat_dec_deg, x = long_dec_deg)) +
  geom_point() +
  theme_classic()

iscn_distinct <- iscn %>% 
  distinct(lat_dec_deg, long_dec_deg, source, label)

## maom
maom <- read_csv("./data/MOC_synthesis.csv") %>% 
  mutate(source = "maom",
         label = "MAOM",
         row_id = row_number()) %>% 
  rename(lat_dec_deg = Lat,
         long_dec_deg = Lon) %>% 
  drop_na(lat_dec_deg) %>% 
  unite(ID, c(row_id, source), remove = FALSE) %>% 
  dplyr::select(ID, lat_dec_deg, long_dec_deg, source, label)
  
head(maom)
skimr::skim_without_charts(maom)

maom %>% 
  ggplot(aes(y = lat_dec_deg, x = long_dec_deg)) +
  geom_point() +
  theme_classic()

maom_distinct <- maom %>% 
  distinct(lat_dec_deg, long_dec_deg, source, label)

## srdb
srdb <- read_csv("./data/srdb-data-V5.csv") %>% 
  mutate(source = "srdb",
         label = "respiration",
         row_id = row_number()) %>% 
  rename(lat_dec_deg = Latitude,
         long_dec_deg = Longitude) %>% 
  drop_na(lat_dec_deg, long_dec_deg) %>% 
  unite(ID, c(row_id, source), remove = FALSE) %>% 
  dplyr::select(ID, lat_dec_deg, long_dec_deg, source, label)

head(srdb)
skimr::skim_without_charts(srdb)

srdb %>% 
  ggplot(aes(y = lat_dec_deg, x = long_dec_deg)) +
  geom_point() +
  theme_classic()

srdb_distinct <- srdb %>% 
  distinct(lat_dec_deg, long_dec_deg, source, label)

## ossl
ossl <- read_csv("./data/ossl_soilsite_L0_v1.2.csv") %>% 
  mutate(source = "ossl",
         label = "spectral",
         row_id = row_number()) %>% 
  rename(lat_dec_deg = latitude.point_wgs84_dd,
         long_dec_deg = longitude.point_wgs84_dd) %>% 
  drop_na(lat_dec_deg, long_dec_deg) %>% 
  unite(ID, c(row_id, source), remove = FALSE) %>% 
  dplyr::select(ID, lat_dec_deg, long_dec_deg, source, label)

head(ossl)
skimr::skim_without_charts(ossl)

ossl %>% 
  ggplot(aes(y = lat_dec_deg, x = long_dec_deg)) +
  geom_point() +
  theme_classic()

ossl_distinct <- ossl %>% 
  distinct(lat_dec_deg, long_dec_deg, source, label)

## timeseries
timeseries <- read_csv("./data/SOC_repeated_latlong_2024May10.csv") %>% 
  mutate(source = "timeseries",
         label = "time series",
         row_id = row_number()) %>% 
  rename(lat_dec_deg = lat,
         long_dec_deg = long) %>% 
  drop_na(lat_dec_deg, long_dec_deg) %>% 
  unite(ID, c(row_id, source), remove = FALSE) %>% 
  dplyr::select(ID, lat_dec_deg, long_dec_deg, source, label)

head(timeseries)
skimr::skim_without_charts(timeseries)

timeseries %>% 
  ggplot(aes(y = lat_dec_deg, x = long_dec_deg)) +
  geom_point() +
  theme_classic()

timeseries_distinct <- timeseries %>% 
  distinct(lat_dec_deg, long_dec_deg, source, label)

## Soils1000
Soils1000 <- read_csv("./data/1000Soils_Metadata_Site_Mastersheet_v1.csv") %>% 
  mutate(source = "soils1000",
         label = "microbial biomass",
         row_id = row_number()) %>% 
  rename(lat_dec_deg = Lat,
         long_dec_deg = Long) %>% 
  drop_na(lat_dec_deg, long_dec_deg) %>% 
  unite(ID, c(row_id, source), remove = FALSE) %>% 
  dplyr::select(ID, lat_dec_deg, long_dec_deg, source, label)

head(Soils1000)
skimr::skim_without_charts(Soils1000)

Soils1000 %>% 
  ggplot(aes(y = lat_dec_deg, x = long_dec_deg)) +
  geom_point() +
  theme_classic()

Soils1000_distinct <- Soils1000 %>% 
  distinct(lat_dec_deg, long_dec_deg, source, label)

## Guillaume
Guillaume <- read_csv("./data/Guillaume_microbial_data.csv") %>% 
  mutate(source = "Guillaume",
         label = "microbial biomass") %>% 
  rename(lat_dec_deg = latitude,
         long_dec_deg = longitude, 
         ID = gID) %>% 
  drop_na(lat_dec_deg, long_dec_deg) %>% 
  dplyr::select(ID, lat_dec_deg, long_dec_deg, source, label)

head(Guillaume)
skimr::skim_without_charts(Guillaume)

Guillaume %>% 
  ggplot(aes(y = lat_dec_deg, x = long_dec_deg)) +
  geom_point() +
  theme_classic()

Guillaume_distinct <- Guillaume %>% 
  distinct(lat_dec_deg, long_dec_deg, source, label)

# Merge both microbial datasets
microbial_distinct <- rbind(Soils1000_distinct, Guillaume_distinct) %>% 
  tibble()

## Necromass
necromass <- read_csv("./data/patel_necromass.csv") %>% 
  mutate(source = "necromass",
         label = "microbial necromass",
         row_id = row_number()) %>% 
  rename(lat_dec_deg = lat,
         long_dec_deg = lon) %>% 
  drop_na(lat_dec_deg, long_dec_deg) %>% 
  unite(ID, c(row_id, source), remove = FALSE) %>% 
  dplyr::select(ID, lat_dec_deg, long_dec_deg, source, label)

head(necromass)
skimr::skim_without_charts(necromass)

necromass %>% 
  ggplot(aes(y = lat_dec_deg, x = long_dec_deg)) +
  geom_point() +
  theme_classic()

necromass_distinct <- necromass %>% 
  distinct(lat_dec_deg, long_dec_deg, source, label)

## Merge all data
all_data <- wosis_distinct %>% 
  full_join(iscn_distinct) %>% 
  full_join(israd_distinct) %>% 
  full_join(maom_distinct) %>% 
  full_join(srdb_distinct) %>% 
  full_join(ossl_distinct) %>%
  full_join(timeseries_distinct) %>%
  full_join(microbial_distinct) %>%
  full_join(necromass_distinct) %>%
  # Remove Antarctica: 
  filter(lat_dec_deg > -58)

names(all_data)

## Extract global data
# MAT
MAT_dir <- "C:/Users/f0076db/Documents/PhD/AfSIS_GlobalData/wc2.0_30s_bio/wc2.0_bio_30s_01.tif"
MAT_global <- raster(MAT_dir)
MAT_all <- raster::extract(MAT_global, cbind(all_data$long_dec_deg,
                                             all_data$lat_dec_deg))

# MAP
MAP_dir <- "C:/Users/f0076db/Documents/PhD/AfSIS_GlobalData/wc2.0_30s_bio/wc2.0_bio_30s_12.tif"
MAP_global <- raster(MAP_dir)
MAP_all <- raster::extract(MAP_global, cbind(all_data$long_dec_deg,
                                             all_data$lat_dec_deg))

# PET
PET_dir <- "C:/Users/f0076db/Documents/PhD/AfSIS_GlobalData/global-et0_annual.tif/et0_yr/et0_yr.tif"
PET_global <- raster(PET_dir) 
PET_resample <- raster::resample(PET_global, MAT_global)

PET_all <- raster::extract(PET_resample, cbind(all_data$long_dec_deg,
                                             all_data$lat_dec_deg))

## GPP (from Fluxcom)
# GPP from FLUXCOM
read_nc_fun <- function(path, pattern = "*.nc"){
  list.files(path, pattern, full.names = TRUE) %>% 
    map(~brick(.))
}

GPP_global <- read_nc_fun(path = "C:/Users/f0076db/Documents/PhD/AfSIS_GlobalData/Fluxcom_GPP")

#calculate yearly means for each year
Global_calc <- map(GPP_global, ~calc(., fun = mean))

#Calculate long-term mean
GPP_mean <- (Global_calc[[1]]+Global_calc[[2]]+Global_calc[[3]]+
               Global_calc[[4]]+Global_calc[[5]]+Global_calc[[6]]+
               Global_calc[[7]]+Global_calc[[8]]+Global_calc[[9]]+
               Global_calc[[10]]+Global_calc[[11]]+Global_calc[[12]])/12

# Match resolution
GPP_resample <- raster::resample(GPP_mean, MAT_global)

GPP_all <- raster::extract(GPP_resample, cbind(all_data$long_dec_deg,
                                           all_data$lat_dec_deg))


data_all_climate <- cbind(all_data, MAT = MAT_all, MAP = MAP_all, PET = PET_all, 
                          GPP = GPP_all)

skimr::skim_without_charts(data_all_climate)

write_csv(data_all_climate, file = paste0("./data/all_data_combined_", 
                                            Sys.Date(), ".csv"))
