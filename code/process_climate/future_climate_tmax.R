# Code that extracts the future weather from CMIP6
# Tutorial: https://geofabio.com/2022/12/13/modelamiento-con-ecocrop-para-identificar-impacto-del-cambio-climatico-sobre-el-cultivo-de-cafe/
# Load libraries ----------------------------------------------------------
library(terra)
library(geodata)

g <- gc(reset = T)
rm(list = ls())

# Download 30s ------------------------------------------------------------
# A valid Shared Socio-economic Pathway code: "126", "245", "370" or "585".
# path = 'tmpr_245'  path = 'tmpr_370'  path = 'tmpr_585'
# (optimistic: SSP245; middle of the road: SSP370; and pessimistic: SSP585)
dataset = 'ACCESS-CM2'
path_dir <-'tmpr_2060_2445'
ssp = '245'

dataset <- c("ACCESS-CM2",  "AWI-CM-1-1-MR", 
             "BCC-CSM2-MR", "CanESM5", "CanESM5-CanOE", 
             "CMCC-ESM2", "CNRM-CM6-1", "CNRM-CM6-1-HR",
             "CNRM-ESM2-1", "EC-Earth3-Veg", "EC-Earth3-Veg-LR",
             "GISS-E2-1-G", "GISS-E2-1-H",
             "INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR",
             "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", 
             "MRI-ESM2-0", "UKESM1-0-LL")


# tmax ---------------------------------------------------------------------
# Future predictions tmax 2041-2060
time = '2061-2080'
prec_w1 <- geodata::cmip6_world(model = dataset[1],
                                ssp = ssp, time = time,
                                var = 'tmax', path = path_dir, res = 2.5)

for(i in c(2:length(dataset))){
  print(paste0("dataset:",dataset[i]))
  prec_w2 <- geodata::cmip6_world(model = dataset[i],
                                  ssp = ssp, time = time,
                                  var = 'tmax', path = path_dir, res = 2.5)
  prec_w1 <- mean(prec_w1, prec_w2)
  print("mean done")
}

writeRaster(prec_w1,paste0("output/ssp245_tmax",
                           time,".tif"))

# Future predictions tmax 2041-2060
time = '2041-2060'
prec_w1 <- geodata::cmip6_world(model = dataset[1],
                                ssp = ssp, time = time,
                                var = 'tmax', path = path_dir, res = 2.5)

for(i in c(2:length(dataset))){
  print(paste0("dataset:",dataset[i]))
  prec_w2 <- geodata::cmip6_world(model = dataset[i],
                                  ssp = ssp, time = time,
                                  var = 'tmax', path = path_dir, res = 2.5)
  prec_w1 <- mean(prec_w1, prec_w2)
  print("mean done")
}

writeRaster(prec_w1,paste0("output/ssp245_tmax",
                           time,".tif"))
