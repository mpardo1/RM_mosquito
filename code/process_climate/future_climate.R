# Code that extracts the future weather from CMIP6
# Tutorial: https://geofabio.com/2022/12/13/modelamiento-con-ecocrop-para-identificar-impacto-del-cambio-climatico-sobre-el-cultivo-de-cafe/
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, tidyverse, glue, sf, RColorBrewer,
               ggspatial, hrbrthemes, showtext, rnaturalearthdata,
               rnaturalearth, extrafont, geodata, data.table, raster,mapSpain)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Font --------------------------------------------------------------------
font_add_google(family = 'Fira Sans', name = 'Fira Sans Condensed')
showtext_auto()

# Functions R_M -----------------------------------------------------------
source("~/INVASIBILITY_THRESHOLD/code/funcR0.R")

# Download data -----------------------------------------------------------
esp0 <- geodata::gadm(country = 'ESP', level = 0,
                      path = 'tmpr')
plot(esp0)
vars <- c('prec', 'tmax', 'tmin')

# Download 30s ------------------------------------------------------------
# A valid Shared Socio-economic Pathway code: "126", "245", "370" or "585".
# path = 'tmpr_245'  path = 'tmpr_370'  path = 'tmpr_585'
# (optimistic: SSP245; middle of the road: SSP370; and pessimistic: SSP585)

 # time = '2061-2080'
 time = "2041-2060"
 var = "prec"
 Path <- paste0("~/INVASIBILITY_THRESHOLD/data/future-climate/",
                var,"_mean",time,".tif")
 prec_w = rast(Path)
 plot(prec_w[[2]])
 
 var = "tmin"
 Path <- paste0("~/INVASIBILITY_THRESHOLD/data/future-climate/",
                var,"_mean",time,".tif")
 tmin_w = rast(Path)
 plot(tmin_w[[2]])
 
 var = "tmax"
 Path <- paste0("~/INVASIBILITY_THRESHOLD/data/future-climate/",
                var,"_mean",time,".tif")
 tmax_w = rast(Path)

# Change coordinate system to WGS84 ---------------------------------------
coord_sys <- crs("+proj=longlat +datum=WGS84")
crs(esp0) <- coord_sys
crs(prec_w) <- coord_sys
crs(tmin_w) <- coord_sys
crs(tmax_w) <- coord_sys

# Extract by mask ---------------------------------------------------------
prec_esp <- terra::crop(prec_w, esp0) %>% terra::mask(., esp0)
# plot(prec_esp[[9]])
tmax_esp <- terra::crop(tmax_w, esp0) %>% terra::mask(., esp0)
tmin_esp <- terra::crop(tmin_w, esp0) %>% terra::mask(., esp0)
plot(tmax_esp[[10]])
plot(tmin_esp[[10]])

# resolution
res_clim <- res(prec_esp)

prec_df <- terra::as.data.frame(prec_esp, xy = TRUE)
tmax_df <- terra::as.data.frame(tmax_esp, xy = TRUE)
tmin_df <- terra::as.data.frame(tmin_esp, xy = TRUE)

# Visualized raster with temp ---------------------------------------------
colnames(prec_df) <- c("x", "y","Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",
                       "Aug", "Sep", "Oct", "Nov", "Dec")
ggplot(prec_df) +
  geom_tile(aes(x = x,y=y,fill=Jan))

# Transform into a sf object ----------------------------------------------
prec_sf <- st_as_sf(prec_df, coords = c('x', 'y'), remove = F)
prec_sf <- prec_sf[,c(3:ncol(prec_sf))]
colnames(prec_sf) <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",
                       "Aug", "Sep", "Oct", "Nov", "Dec", "geometry")
tmax_sf <- st_as_sf(tmax_df, coords = c('x', 'y'), remove = F)
tmax_sf <- tmax_sf[,c(3:ncol(tmax_sf))]
colnames(tmax_sf) <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",
                       "Aug", "Sep", "Oct", "Nov", "Dec", "geometry")
tmin_sf <- st_as_sf(tmin_df, coords = c('x', 'y'), remove = F)
tmin_sf <- tmin_sf[,c(3:ncol(tmin_sf))]
colnames(tmin_sf) <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",
                       "Aug", "Sep", "Oct", "Nov", "Dec", "geometry")

# Spain map with NATCODE --------------------------------------------------
library(mapSpain)
esp_can <- esp_get_munic_siane(moveCAN = FALSE)
can_box <- esp_get_can_box()
esp_can$NATCODE <- as.numeric(paste0("34",esp_can$codauto,
                                     esp_can$cpro,
                                     esp_can$LAU_CODE))
# Transform in to the same coord. system
esp_can <- st_transform(esp_can, crs = coord_sys)

# Save data for future use --------------------------------------------------
saveRDS(esp_can, "~/INVASIBILITY_THRESHOLD/data/future-climate/esp_can.Rds")
saveRDS(prec_sf, "~/INVASIBILITY_THRESHOLD/data/future-climate/prec_sf.Rds")

# Intersect data from climate to mapSpain -----------------------------------
st_crs(prec_sf) <- st_crs(esp_can)
st_crs(tmax_sf) <- st_crs(esp_can)
st_crs(tmin_sf) <- st_crs(esp_can)
geometry <-  esp_can[,"geometry"]

# Read intersection compute in intersect_raster_sf.R ------------------------
intersect_df <- readRDS("~/INVASIBILITY_THRESHOLD/data/future-climate/intersect_df_parall.Rds")
intersect_df <- do.call(rbind, intersect_df)
colnames(intersect_df) <- c("points", "geometry")

# Join columns NATCODE y lon-lat --------------------------------------------
intersect_df$geom <- prec_sf[intersect_df$point,"geometry"]
esp_can$geometry <- NULL
esp_can <- as.data.frame(esp_can)
intersect_df$NATCODE <- as.numeric(esp_can[intersect_df$geometry,"NATCODE"])

# Free memory ---------------------------------------------------------------
rm(prec_df, prec_esp, prec_w,
   tmin_df, tmin_esp, tmin_w,
   tmax_df, tmax_esp, tmax_w, geometry)

# Add to prec_sf NATCODE ----------------------------------------------------
prec_sf$ind <- seq(1,nrow(prec_sf),1)
prec_sf$geometry <- NULL
prec_sf <- intersect_df %>%
  left_join(prec_sf, by = join_by(points == ind))
prec_sf$geom <- NULL

# Add to prec_sf NATCODE ----------------------------------------------------
tmin_sf$ind <- seq(1,nrow(tmin_sf),1)
tmin_sf$geometry <- NULL
tmin_sf <- intersect_df %>%
  left_join(tmin_sf, by = join_by(points == ind))
tmin_sf$geom <- NULL

# Add to prec_sf NATCODE ----------------------------------------------------
tmax_sf$ind <- seq(1,nrow(tmax_sf),1)
tmax_sf$geometry <- NULL
tmax_sf <- intersect_df %>%
  left_join(tmax_sf, by = join_by(points == ind))
tmax_sf$geom <- NULL

# Check with map visually ---------------------------------------------------
esp_can <- esp_get_munic_siane(moveCAN = TRUE)
can_box <- esp_get_can_box()
esp_can$NATCODE <- as.numeric(paste0("34",esp_can$codauto,
                                     esp_can$cpro,
                                     esp_can$LAU_CODE))
df_plot <- esp_can %>% left_join(tmax_sf , by = join_by(NATCODE == NATCODE))
ggplot(df_plot) +
  geom_sf(aes(fill = Oct), colour = NA) +
  geom_sf(data = can_box) +
  scale_fill_viridis_c(option = "magma")

list_NA_NATCODE <- df_plot[which(is.na(df_plot$Oct)), "NATCODE"]
list_NA_NATCODE$geometry.x <- NULL
saveRDS(list_NA_NATCODE, "~/INVASIBILITY_THRESHOLD/data/list_NA_futureclim.Rds")

# Population density 2022------------------------------------------------------
Path <- "~/INVASIBILITY_THRESHOLD/data/pop/pobmun22.csv"
pop22 <- read.csv(Path, sep = ",")
pop22$cmun <- ifelse(pop22$CMUN<10, paste0("00",pop22$CMUN),
                     ifelse(pop22$CMUN<100, paste0("0",pop22$CMUN),
                            as.character(pop22$CMUN)))
pop22$cpro <- ifelse(pop22$CPRO<10,
                     paste0("0",pop22$CPRO),as.character(pop22$CPRO))
pop <- esp_can %>% left_join(pop22)
pop[which(is.na(pop$POB22)),"POB22"] <- 0
pop$area <- as.numeric(st_area(pop))/1000000
pop$dens <- pop$POB22/pop$area
pop <- pop[, c("dens","NATCODE")]
pop$geometry <- NULL

# Load population density
Path <- "/home/marta/INVASIBILITY_THRESHOLD/data/pop/dens_fut.Rds"
pop <- readRDS(Path)
colnames(pop) <- c("NATCODE", "dens")

# Set to province data muni with no data -------------------------------------
prov_muni <- setDT(esp_can[,c("NATCODE", "cpro")])
prov_muni$geometry <- NULL

rm_NA <- function(tmax_sf){
  tmax_sf <- prov_muni %>% left_join(tmax_sf)
  tmax_sf_prov <- tmax_sf[which(is.na(tmax_sf$Jan)==FALSE ),
                          .(Jan1 = mean(Jan), Feb1 = mean(Feb),
                             Mar1 = mean(Mar), Apr1 = mean(Apr),
                             May1 = mean(May), Jun1 = mean(Jun),
                             Jul1 = mean(Jul), Aug1 = mean(Aug),
                             Sep1 = mean(Sep), Oct1 = mean(Oct),
                             Nov1 = mean(Nov), Dec1 = mean(Dec)),
                     by = list(cpro)]

  tmax_sf <- tmax_sf %>% left_join(tmax_sf_prov)
  tmax_sf$Jan <- ifelse(is.na(tmax_sf$Jan), tmax_sf$Jan1,
                        tmax_sf$Jan )
  tmax_sf$Feb <- ifelse(is.na(tmax_sf$Feb), tmax_sf$Feb1,
                        tmax_sf$Feb )
  tmax_sf$Mar <- ifelse(is.na(tmax_sf$Mar), tmax_sf$Mar1,
                        tmax_sf$Mar )
  tmax_sf$Apr <- ifelse(is.na(tmax_sf$Apr), tmax_sf$Apr1,
                        tmax_sf$Apr )
  tmax_sf$May <- ifelse(is.na(tmax_sf$May), tmax_sf$May1,
                        tmax_sf$May )
  tmax_sf$Jun <- ifelse(is.na(tmax_sf$Jun), tmax_sf$Jun1,
                        tmax_sf$Jun )
  tmax_sf$Jul <- ifelse(is.na(tmax_sf$Jul), tmax_sf$Jul1,
                        tmax_sf$Jul )
  tmax_sf$Aug <- ifelse(is.na(tmax_sf$Aug), tmax_sf$Aug1,
                        tmax_sf$Aug )
  tmax_sf$Sep <- ifelse(is.na(tmax_sf$Sep), tmax_sf$Sep1,
                        tmax_sf$Sep )
  tmax_sf$Oct <- ifelse(is.na(tmax_sf$Oct), tmax_sf$Oct1,
                        tmax_sf$Oct )
  tmax_sf$Nov <- ifelse(is.na(tmax_sf$Nov), tmax_sf$Nov1,
                        tmax_sf$Nov )
  tmax_sf$Dec <- ifelse(is.na(tmax_sf$Dec), tmax_sf$Dec1,
                        tmax_sf$Dec )

  tmax_sf <- as.data.frame(tmax_sf)
  if(nrow(tmax_sf[which(is.na(tmax_sf$Aug)),]) == 0){
    print("No NA")
  }
  return(tmax_sf)
}

tmax_sf[which(is.na(tmax_sf$Oct)),"NATCODE"]
tmax_sf <- rm_NA(tmax_sf)
tmin_sf <- rm_NA(tmin_sf)
prec_sf <- rm_NA(prec_sf)


# Function to compute plot and df monthly --------------------------------------
plot_month <- function(month, esp){
  # Join df by months ----------------------------------------------------------
  tmin_Jan <- tmin_sf[,c("NATCODE", month)]
  colnames(tmin_Jan) <- c("NATCODE", "tmin")
  tmax_Jan <- tmax_sf[,c("NATCODE", month)]
  colnames(tmax_Jan) <- c("NATCODE", "tmax")
  prec_Jan <- prec_sf[,c("NATCODE", month)]
  colnames(prec_Jan) <- c("NATCODE", "prec")
  
  # Group by NATCODE there are more than 1 raster point at each municipality
  tmin_Jan <- tmin_Jan %>% group_by(NATCODE) %>%
    summarise(tmin = mean(tmin))
  tmax_Jan <- tmax_Jan %>% group_by(NATCODE) %>%
    summarise(tmax = mean(tmax))
  prec_Jan <- prec_Jan %>% group_by(NATCODE) %>%
    summarise(prec = mean(prec))
  
  df_Jan <- setDT(prec_Jan %>% left_join(tmax_Jan) %>% 
                    left_join(tmin_Jan) %>% left_join(pop))
  df_Jan$tmean <- (df_Jan$tmin + df_Jan$tmax)/2
  
  # Extract month number
  month_num <- ifelse(month == "Jan", 1, ifelse(
                        month == "Feb", 2, ifelse(
                          month == "Mar", 3, ifelse(
                            month == "Apr", 4, ifelse(
                              month == "May", 5, ifelse(
                                month == "Jun", 6, ifelse(
                                  month == "Jul", 7, ifelse(
                                    month == "Aug", 8, ifelse(
                                      month == "Sep", 9, ifelse(
                                        month == "Oct", 10, ifelse(
                                          month == "Nov", 11, 12)))))))))))
  
  # The rainfall it is accumulated months, but we need an average per day.
  df_Jan$prec <- df_Jan$prec/as.numeric(days_in_month(month_num))
  df_Jan[, R0_alb := mapply(R0_func_alb, tmean, prec, dens)]
  df_Jan[, R0_aeg := mapply(R0_func_aeg, tmean, prec, dens)]
  df_Jan[, R0_jap := mapply(R0_func_jap, tmean, prec, dens)]
  
  NA_df <- df_Jan[which(is.na(R0_alb)),]
  
  df_Jan_p <- esp_can %>% left_join(df_Jan)
  if(esp == "alb"){
    plot <- ggplot(df_Jan_p) +
      geom_sf(data = can_box) +
      geom_sf(aes(fill = R0_alb), colour = NA) +
      coord_sf(datum = NA) +
      scale_fill_distiller(palette = "Spectral",
                           limits = c(min(df_Jan_p$R0_alb),max(df_Jan_p$R0_alb)),
                           name = TeX("$R_M$")) +
      ggtitle(as.character(month_num)) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = c(0.1,0.6))
  }else{
    plot <- ggplot(df_Jan_p) +
      geom_sf(data = can_box) +
      geom_sf(aes(fill = R0_aeg), colour = NA) +
      coord_sf(datum = NA) +
      scale_fill_distiller(palette = "Spectral",
                           limits = c(min(df_Jan_p$R0_aeg),max(df_Jan_p$R0_aeg)),
                           name = TeX("$R_M$")) +
      ggtitle(as.character(month_num)) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = c(0.1,0.6))
  }
  
  
  return(list(plot,df_Jan))
}

library(latex2exp)

# Plot whole year monthly ----------------------------------------------------
list_month = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",
               "Aug", "Sep", "Oct", "Nov", "Dec")
ind = 3
# plot_3 <- plot_month(list_month[ind], "alb")[[1]]
plot_3 <- plot_month(list_month[ind], "aeg")[[1]]
plot_11
library(ggpubr)
ggarr <- ggarrange(plot_3,plot_4,plot_5,
                   plot_6,plot_7,plot_8,
                   plot_9,plot_10,plot_11,
                   nrow=3,ncol = 3, common.legend = TRUE)

ggarr

# Produce a monthly df for all year ------------------------------------------
lmon <- list("Feb", "Mar", "Apr", "May", "Jun", "Jul", 
             "Aug", "Sep", "Oct", "Nov", "Dec")
library(latex2exp)
df_y <- plot_month("Jan","alb")[[2]]
df_y$month <- "Jan"
for(i in c(1:length(lmon))){
  df_aux <- plot_month(lmon[[i]],"alb")[[2]]
  df_aux$month <- lmon[[i]]
  df_y <- rbind(df_y,df_aux)
}

Path <- paste0("~/INVASIBILITY_THRESHOLD/output/ERA5/temp/2020/monthly_clim_2060.Rds")
saveRDS(df_y, Path)

# Create df for num months suitable -----------------------------------------
df_y$bool_alb <- ifelse(df_y$R0_alb <1,0,1)
df_y$bool_aeg <- ifelse(df_y$R0_aeg <1,0,1)
df_y$bool_jap <- ifelse(df_y$R0_jap <1,0,1)

df_g <- df_y %>% group_by(NATCODE) %>%
  summarise(alb = sum(bool_alb),
            aeg = sum(bool_aeg),
            jap = sum(bool_jap),
            avg_alb = mean(R0_alb),
            avg_aeg = mean(R0_aeg),
            avg_jap = mean(R0_jap)
            )

library(RColorBrewer)
name_pal = "RdYlBu"
display.brewer.pal(11, name_pal)
pal <- rev(brewer.pal(11, name_pal))
pal[11]
pal[12] = "#74011C"
pal[13] = "#4B0011"
ggplot(df_g) +
  geom_sf(data = can_box) +
  geom_sf(aes(fill = as.factor(alb)), colour = NA) +
  coord_sf(datum = NA) +
  scale_fill_manual(values = pal,
                       limits = factor(seq(0,12,1)),
                       name = TeX("$R_M$")) +
  theme_bw()

# Save file
Path <- paste0("~/INVASIBILITY_THRESHOLD/output/ERA5/temp/2020/clim_2060.Rds")
saveRDS(df_g, Path)

df_g <- readRDS(Path)

# Create plot months suitable -----------------------------------------
esp_can <- esp_get_munic_siane(moveCAN = TRUE)
can_box <- esp_get_can_box()
esp_can$NATCODE <- as.numeric(paste0("34",esp_can$codauto,
                                     esp_can$cpro,
                                     esp_can$LAU_CODE))

df_g <- esp_can %>% left_join(df_g)
name_pal = "RdYlBu"
display.brewer.pal(11, name_pal)
pal <- rev(brewer.pal(11, name_pal))
pal[11]
pal[12] = "#74011C"
pal[13] = "#4B0011"
letsize = 16

library(latex2exp)
# Albopictus
alb <- ggplot(df_g) + 
  geom_sf(aes(fill = as.factor(alb), alpha = 0.4*as.factor(alb)), colour = NA) +
  geom_sf(data = can_box) + coord_sf(datum = NA) +
  scale_fill_manual(values = pal, limits = c(1:12),
                    name = TeX("$R_M>1$ \n months")) +
  theme_bw()

# Aegipty
aeg <- ggplot(df_g) + 
  geom_sf(aes(fill = as.factor(aeg)), colour = NA) +
  geom_sf(data = can_box) + coord_sf(datum = NA) +
  scale_fill_manual(values = pal, limits = c(1:12),
                    name = TeX("$R_M>1$ \n months")) +
  theme_bw()

# Japonicus
jap <- ggplot(df_g) + 
  geom_sf(aes(fill = as.factor(jap)), colour = NA) +
  geom_sf(data = can_box) + coord_sf(datum = NA) +
  scale_fill_manual(values = pal, limits = c(1:12),
                    name = TeX("$R_M>1$ \n months")) +
  theme_bw()

# Join all species ------------------------------------------------
# Extract the legend
library(ggpubr)
legend_only <- get_legend(alb +
                            theme(legend.position = "top"))
ggarrange(alb + ggtitle(expression(italic("Ae. albopictus")))+
            theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5)),
          aeg + ggtitle(expression(italic("Ae. aegypti")))+
            theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5)),
          jap + ggtitle(expression(italic("Ae. japonicus")))+
            theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5)),
          legend_only,
          ncol=2, nrow = 2)

# Compare a month 3 species --------------------------------------------------
