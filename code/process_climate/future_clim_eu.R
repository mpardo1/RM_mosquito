# Code that extracts the future weather from CMIP6
# Tutorial: https://geofabio.com/2022/12/13/modelamiento-con-ecocrop-para-identificar-impacto-del-cambio-climatico-sobre-el-cultivo-de-cafe/
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, eurostat,tidyverse, sf, RColorBrewer,
               geodata, data.table,giscoR)

g <- gc(reset = T)
rm(list = ls())

# Functions R_M -----------------------------------------------------------
source("funcR0.R")

# Download data -----------------------------------------------------------
SHP_0 <- get_eurostat_geospatial(resolution = 10, 
                                 nuts_level = 0, 
                                 year = 2016)

# Load future climate data average all data sets------------------------------------------------
# code in ~/INVASIBILITY_THRESHOLD/code/future-climate/avg_all_dataset_future_clim.R 
time = "2061-2080"
var = "prec"
Path <- paste0("output/ssp245_",
                var,time,".tif")
prec_w = rast(Path)
plot(prec_w[[2]])

var = "tmin"
Path <- paste0("output/ssp245_",
               var,time,".tif")
tmin_w = rast(Path)
plot(tmin_w[[2]])

var = "tmax"
Path <- paste0("output/ssp245_",
               var,"_mean",time,".tif")
tmax_w = rast(Path)
plot(tmax_w[[2]])

# Crop the raster to the exact extent ------------------------
exact_extent <- c(xmin = -25, xmax = 40, ymin = 25, ymax = 75)
tmin_w <- crop(tmin_w, exact_extent)
tmax_w <- crop(tmax_w, exact_extent)
prec_w <- crop(prec_w, exact_extent)

# read population density europe -----------------------------
Path <- "data/GHS_POP_E2030_GLOBE_R2023A_reprojected_bilinear_025.tif"
pop_eu <- rast(Path)
# plot(log(pop_eu[[1]]))
pop_eu <- terra::project(pop_eu,tmin_w, method = "near")

# Create a grid of longitude and latitude values
lon <- seq(from = xmin(tmin_w), to = xmax(tmin_w),
           by = res(tmin_w)[1])
lat <- seq(from = ymin(tmin_w), to = ymax(tmin_w),
           by = res(tmin_w)[2])
grid_points <- expand.grid(lon = lon, lat = lat)

# extract values as df --------------------------------------
# plot(tmin_w[[1]])
tmin_w <- terra::extract(tmin_w,
                       grid_points, xy =TRUE)
colnames(tmin_w)[2:13] <- as.character(c(1:12))
tmin_w <- reshape::melt(tmin_w[,c(1:13)],id.vars = "ID")
colnames(tmin_w) <- c("id", "month", "tmin")

# plot(tmax_w[[1]])
tmax_w <- terra::extract(tmax_w,
                         grid_points, xy =TRUE)
colnames(tmax_w)[2:13] <- as.character(c(1:12))
tmax_w <- reshape::melt(tmax_w[,c(1:13)],id.vars = "ID")
colnames(tmax_w) <- c("id", "month", "tmax")

# plot(prec_w[[1]])
prec_w <- terra::extract(prec_w,
                         grid_points, xy =TRUE)
colnames(prec_w)[2:13] <- as.character(c(1:12))
prec_w <- reshape::melt(prec_w[,c(1:13)],id.vars = "ID")
colnames(prec_w) <- c("id", "month", "prec")
prec_w$prec <- prec_w$prec/(as.numeric(days_in_month(as.Date(paste0("2020-",prec_w$month,"-01"),
                                                  format ="%Y-%m-%d"))))

# Extract popilation density
# Remove non European countries
sf_eu <- gisco_get_countries(year = "2020", region = "Europe")
sf_eu <- sf_eu[sf_eu$CNTR_ID != "BY" &
                 sf_eu$CNTR_ID != "RU"  & sf_eu$CNTR_ID != "IS", ]
plot(sf_eu[,"CNTR_ID"])
pop_eu <- terra::crop(pop_eu, sf_eu) %>% terra::mask(., sf_eu)
pop <- terra::extract(pop_eu,
                         grid_points, xy =TRUE)
pop <- pop[,c(1:2)]
colnames(pop) <- c("id", "pop")

# Join all df -------------------------------------------
clim_df <- setDT(tmin_w) %>% left_join(setDT(tmax_w))
rm(tmin_w,tmax_w)
clim_df$tmean <- (clim_df$tmin + clim_df$tmax)/2
clim_df$tmin <- NULL
clim_df$tmax <- NULL

clim_df <- setDT(clim_df) %>% left_join(setDT(prec_w))
clim_df <- clim_df %>% left_join(setDT(pop))
rm(prec_w,pop)

# compute R0 --------------------------------------------------
clim_df[, R0_alb := mapply(R0_func_alb, tmean, prec, pop)]
clim_df[, R0_aeg := mapply(R0_func_aeg, tmean, prec, pop)]
clim_df[, R0_jap := mapply(R0_func_jap, tmean, prec, pop)]

# plots seasonal ----------------------------------------------
grid_points$id <- c(1:nrow(grid_points))
clim_df <- clim_df %>% left_join(grid_points)
 saveRDS(clim_df,
        paste0("output/eu_alb_aeg_",time,"_.Rds"))
 
clim_df <- readRDS(paste0("output/eu_alb_aeg_",time,"_.Rds"))
library(latex2exp)
month_n =5
plot_5_aeg <- ggplot(clim_df[which(clim_df$month == month_n),],
       aes(x = lon, y = lat, fill = R0_aeg)) +
  geom_raster() +
  scale_fill_distiller(na.value = "white",
                       palette = "Spectral",
                       name = TeX("$R_M$"),
                       limits = c(0,6))+
  ylim(c(25,75)) + xlim(c(-30,40)) +
  ggtitle(month_n) +
  theme_minimal() + coord_fixed() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "italic"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.margin = unit(c(0, 0, 0, 0), "null"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks.length = unit(0, "null"),
        axis.ticks.margin = unit(0, "null")) 
plot_11_alb

library(ggpubr)
ggarrange(plot_5_alb, plot_8_alb,plot_11_alb,
          plot_5_aeg, plot_8_aeg,plot_11_aeg,
          ncol = 3,nrow=2, common.legend= TRUE)

# aggregate by year -------------------------------------------
clim_df$bool_alb <- ifelse(clim_df$R0_alb>1,1,0)
clim_df$bool_aeg <- ifelse(clim_df$R0_aeg>1,1,0)
clim_df$bool_jap <- ifelse(clim_df$R0_jap>1,1,0)

clim_df <- clim_df[,.(sum_alb = sum(bool_alb),
                        sum_aeg = sum(bool_aeg),
                        sum_jap = sum(bool_jap)), by = list(id)]

# plot  sum months ------------------------------------------------
grid_points$id <- c(1:nrow(grid_points))
clim_df <- clim_df %>% left_join(grid_points)
saveRDS(clim_df,paste0("output/summon_eu_",time,".Rds"))
clim_df <- readRDS(paste0("output/summon_eu_",time,".Rds"))

library(RColorBrewer)
name_pal = "RdYlBu"
display.brewer.pal(11, name_pal)
pal <- rev(brewer.pal(11, name_pal))
pal[11]
pal[12] = "#74011C"
pal[13] = "#4B0011"
letsize = 16
alb <- ggplot(clim_df) +
  geom_raster(aes(x = lon, y = lat, 
                fill = as.factor(sum_alb)),alpha = 1) +
  scale_fill_manual(values = pal,
                    name = "Nº suitable \n months",
                    limits = factor(seq(0,12,1)),
                    na.value = "white") +
  ylim(c(25,75)) + xlim(c(-30,40)) +
  xlab("") + ylab("") +
  theme_minimal() + coord_fixed() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "italic"),
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "null"),
    panel.margin = unit(c(0, 0, 0, 0), "null"),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks.length = unit(0, "null"),
    axis.ticks.margin = unit(0, "null"), 
    legend.position = "rigth")


library(ggpubr)
legend_only <- get_legend(ggplot(clim_df) +
  geom_tile(aes(x = lon, y = lat, 
                fill = as.factor(sum_alb)),alpha = 1) +
  scale_fill_manual(values = pal,
                    name = "Nº suitable \n months",
                    limits = factor(seq(0,12,1)),
                    na.value = "white") +
  theme_minimal() +
  theme(legend.position = "right"))


plot(legend_only)
aeg <- ggplot(clim_df, aes(x = lon, y = lat,
                    fill = as.factor(sum_aeg))) +
  geom_tile() +
  scale_fill_manual(values = pal,
                    name = "Nº suitable \n months",
                    limits = factor(seq(0,12,1)),
                    na.value = "white")+
  xlab("") + ylab("") +
  ylim(c(25,75)) + xlim(c(-30,40)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "italic"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.margin = unit(c(0, 0, 0, 0), "null"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks.length = unit(0, "null"),
        axis.ticks.margin = unit(0, "null"),
        legend.text = element_text(14), 
        legend.position = "none" )

# jap <- ggplot(clim_df, aes(x = lon, y = lat,
#                     fill = as.factor(sum_jap))) +
#   geom_tile() +
#   xlab("") + ylab("") +
#   scale_fill_manual(values = pal,
#                     name = "Nº months\n suitable",
#                     limits = factor(seq(0,12,1)),
#                     na.value = "white")+
#   ylim(c(25,75)) + xlim(c(-30,40)) +
#   ggtitle("Aedes japonicus 2061-2080") +
#   theme_minimal() +
#   theme(plot.title = element_text(hjust = 0.5,
#                                   face = "italic"),
#         panel.background = element_rect(fill = "transparent", colour = NA),
#         plot.background = element_rect(fill = "transparent", colour = NA),
#         panel.grid = element_blank(),
#         panel.border = element_blank(),
#         plot.margin = unit(c(0, 0, 0, 0), "null"),
#         panel.margin = unit(c(0, 0, 0, 0), "null"),
#         axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank(),
#         axis.line = element_blank(),
#         axis.ticks.length = unit(0, "null"),
#         axis.ticks.margin = unit(0, "null"),
#         legend.text = element_text(14) ) +
#   guides(fill = guide_legend(nrow = 1),
#          label.position = "top")

library(ggpubr)
# ggarrange(alb, aeg, jap,
#           ncol = 3, common.legend = TRUE)
ggarrange(alb, aeg,
          ncol = 2, common.legend = TRUE)
