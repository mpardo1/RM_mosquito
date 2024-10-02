rm(list = ls())
library(eurostat)
library(terra)
library(tidyverse)
library(data.table)
library(sf)
library(giscoR)
source("funcR0.R")

# Load population density from: https://sedac.ciesin.columbia.edu/data/collection/gpw-v4
Path <- "data/gpw_v4_population_density_rev11_2020_2pt5_min.tif"
pop <- rast(Path)
plot(log(pop))

# Crop the extent
exact_extent <- c(xmin = -25, xmax = 40, ymin = 25, ymax = 75)
pop_eu <- crop(pop,exact_extent)
plot(log(pop_eu))

# Remove non European countries
sf_eu <- gisco_get_countries(year = "2020", region = "Europe")
sf_eu <- sf_eu[sf_eu$CNTR_ID != "BY" &
                 sf_eu$CNTR_ID != "RU"  & sf_eu$CNTR_ID != "IS", ]
plot(sf_eu[,"CNTR_ID"])
pop_eu <- terra::crop(pop_eu, sf_eu) %>% terra::mask(., sf_eu)
writeRaster(pop_eu,paste0("output/pop_eu.tif"),overwrite=TRUE)

# climate data europe ----------------------------------------
# https://cds.climate.copernicus.eu/datasets/reanalysis-cerra-land?tab=overview
rain_eu <- rast("data/monthly_mean_rain_2020_europe.grib") # Total precipitation
rain_eu <- tapp(rain_eu, "month", mean) # Compute monthly mean
plot(rain_eu[[3]])

temp_eu <- rast("data/monthly_mean_temp_2020_europe.grib") # Skin temperature
temp_eu <- tapp(temp_eu, "month", mean) # Compute monthly mean
plot(temp_eu[[7]])

# Load map CMIP6 to have the same coord sys for present and future maps
Path <- paste0("output/245_tmin2041-2060.tif") # File from CMIP6
tmin_w = rast(Path)
exact_extent <- c(xmin = -25, xmax = 40, ymin = 25, ymax = 75)
tmin_w <- crop(tmin_w, exact_extent)

# change coordinate system to crop ---------------------------
temp_eu <- terra::project(temp_eu,tmin_w, method = "average")
pop_eu <- terra::project(pop_eu,tmin_w, method = "near")
rain_eu <- terra::project(rain_eu,tmin_w, method = "average")

plot(pop_eu[[1]])
plot(temp_eu[[1]])
plot(rain_eu[[5]])

# Create a grid of longitude and latitude values
lon <- seq(from = xmin(temp_eu), to = xmax(temp_eu),
           by = res(temp_eu)[1])
lat <- seq(from = ymin(temp_eu), to = ymax(temp_eu),
           by = res(temp_eu)[2])
grid_points <- expand.grid(lon = lon, lat = lat)
# extract values as df --------------------------------------
temp <- terra::extract(temp_eu,
                       grid_points, xy =TRUE)
colnames(temp)[2:13] <- c(1:12)
temp <- reshape::melt(temp[,c(1:13)],id.vars = "ID")
colnames(temp) <- c("id", "month", "tmean")
rain <- terra::extract(rain_eu,
                       grid_points, xy =TRUE)
colnames(rain)[2:13] <- c(1:12)
rain <- reshape::melt(rain[,c(1:13)],id.vars = "ID")
colnames(rain) <- c("id", "month", "prec")

# transform rain into mm per squared km
rain$prec <- rain$prec*1000
pop <- terra::extract(pop_eu,
                      grid_points, xy =TRUE)[,c(1:2)]
colnames(pop) <- c("id", "pop")

# join df -----------------------------------------------------
clim <- temp %>% left_join(rain)
clim_pop <- setDT(clim %>% left_join(pop))

# create df europe climate monthly data -----------------------
grid_points$id <- c(1:nrow(grid_points))
clim_pop <- clim_pop %>% left_join(grid_points)

# compute R0 --------------------------------------------------
clim_pop[, R0_alb := mapply(R0_func_alb, tmean, prec, pop)]
clim_pop[, R0_aeg := mapply(R0_func_aeg, tmean, prec, pop)]
# clim_pop[, R0_jap := mapply(R0_func_jap, tmean, prec, pop)]
saveRDS(clim_pop,
        paste0("output/eu_clim_same_coords_",2020,".Rds"))
# clim_pop <- readRDS(paste0("~/INVASIBILITY_THRESHOLD/data/ERA5/Europe/eu_clim_",2020,".Rds"))

clim_pop$bool_alb <- ifelse(clim_pop$R0_alb>1,1,0)
clim_pop$bool_aeg <- ifelse(clim_pop$R0_aeg>1,1,0)
clim_pop$bool_jap <- ifelse(clim_pop$R0_jap>1,1,0)

# Group by location
clim_pop <- clim_pop[,.(sum_alb = sum(bool_alb),
                        sum_aeg = sum(bool_aeg),
                        sum_jap = sum(bool_jap)), by = list(id)]

# add the lon lat -------------------------------------
grid_points$id <- c(1:nrow(grid_points))
saveRDS(grid_points, "~/output/grid_points.Rds")

clim_pop <- clim_pop %>% left_join(grid_points)

ggplot(clim_pop,
       aes(x = lon, y = lat,
           fill = sum_alb)) +
  geom_raster()

# save the df
saveRDS(clim_pop,
        paste0("output/eu_R0_fitfuture_clim_",2020,".Rds"))

clim_pop <- readRDS("output/eu_R0_fitfuture_clim_2020.Rds")
