rm(list=ls())
library(mcera5)
library(mapSpain)
library(sf)
library(ggplot2)
library(tidyverse)
library(parallel)
library(data.table)
library("viridis")
library("gganimate")

source("~/INVASIBILITY_THRESHOLD/code/funcR0.R")
#----------------------------------------------------------------------#
# Check hatching rate japonicus
vec <- seq(0,20,0.1)
out <- sapply(vec, h_f_jap, land = 0.4)
df_h <- data.frame(vec, out)
ggplot(df_h) + geom_line(aes(vec,out)) +
  xlab("Monthly average rainfall, mm") +
  ylab("Hatching rate") + theme_bw()


## Read the data for the R0 computed daily:
# year = 2022
# Path <- paste0("/home/marta/INVASIBILITY_THRESHOLD/output/mcera5/process_Daily_ERA5_daily_mcera_",year,".Rds")
# year = 2004
# Path <- paste0("/home/marta/INVASIBILITY_THRESHOLD/output/mcera5/process_hourly_daily_ERA5_daily_mcera_",
#                year,".Rds")
year = 2020
Path <- paste0("~/INVASIBILITY_THRESHOLD/output/ERA5/temp/2020/clim_",year,".Rds")
# saveRDS(dt_weather,Path)
df_group <- setDT(readRDS(Path))
df_group$id <- 1
test <- df_group[,.(n = sum(id)), by = list(NATCODE)]
min(test$n)
max(test$n)
nrow(test)

# Load Spanish municipality maps
esp_can <- esp_get_munic_siane(moveCAN = TRUE)
can_box <- esp_get_can_box()
esp_can$NATCODE <- as.numeric(paste0("34",esp_can$codauto,
                                     esp_can$cpro,
                                     esp_can$LAU_CODE))
# esp_can$cpro <- as.numeric(esp_can$cpro)
# write_sf(esp_can, "~/Dengue/data/esp_muni.shp")
df_group$month <- lubridate::month(df_group$date)

# Check rainfall ----------------------------------------------------
df_group_y  <- df_group[,.(prec = sum(prec1)), by = list(NATCODE)]
df_group_y  <- esp_can %>% left_join(df_group_y )
ggplot(df_group_y) +  geom_sf(aes(fill = prec), color = NA) + 
  scale_fill_viridis_c()

# Add population density ---------------------------------------------
Path <- "/home/marta/INVASIBILITY_THRESHOLD/data/pop/pobmun20.csv"
pop22 <- read.csv(Path, sep = ";")
# Path <- "/home/marta/INVASIBILITY_THRESHOLD/data/pop/pobmun04.csv"
# pop22 <- read.csv(Path, sep = ";")
pop22$cmun <- ifelse(pop22$CMUN<10, paste0("00",pop22$CMUN),
                     ifelse(pop22$CMUN<100, paste0("0",pop22$CMUN),
                            as.character(pop22$CMUN)))
pop22$cpro <- ifelse(pop22$CPRO<10,
                     paste0("0",pop22$CPRO),as.character(pop22$CPRO))
pop22$POB22 <- gsub('\\.','',pop22$TOTAL)
pop22$POB22 <- as.numeric(pop22$POB22)
pop22$POB22 <- as.numeric(pop22$POB20)
esp_can <- esp_can %>% left_join(pop22, by = c("cpro","cmun") )
test_pop <- esp_can[which(is.na(esp_can$POB22)),
                    c("name", "NATCODE", "POB22")]

test_pop$geometry <- NULL
nrow(test_pop)
esp_can[which(is.na(esp_can$POB22)),"POB22"] <- 0

# Transform squared meters (output st_area) to squared km
esp_can$area <- as.numeric(st_area(esp_can))/1000000
esp_can$dens <- esp_can$POB22/esp_can$area

# Path <- "~/INVASIBILITY_THRESHOLD/output/pop/pop22.Rds"
# saveRDS(esp_can_pop, Path)

# Join landcover ---------------------------------------------------
path <- "~/INVASIBILITY_THRESHOLD/data/landcov_fracs.Rds"
landcover <- readRDS(path)
landcover <- landcover[landcover$value == 23 ,c("NATCODE", "freq")]
landcover[is.na(landcover$freq) , "freq"] <- 0
esp_can <- esp_can %>% left_join(landcover)
esp_can[is.na(esp_can$freq) , "freq"] <- 0
ggplot(esp_can) + geom_sf(aes(fill = freq), color = NA) +
  scale_fill_viridis_c()

# Set to data table and join with climateic variables -------------
esp_can <- setDT(esp_can[,c("NATCODE","name", "dens", "freq")])
esp_can$geometry <- NULL
# Join pop and landcover data with climate conditions
df_group <- df_group %>% left_join(esp_can)

# Test weather data -------------------------------------------------
esp_can <- esp_get_munic_siane(moveCAN = TRUE)
can_box <- esp_get_can_box()
esp_can$NATCODE <- as.numeric(paste0("34",esp_can$codauto,
                                     esp_can$cpro,
                                     esp_can$LAU_CODE))
df_day <- df_group[which(df_group$date == as.Date("2020-07-05")),]
df_day <- esp_can %>% left_join(df_day)
ggplot(df_day) + 
  geom_sf(aes(fill = prec1 ), color = NA) +
  scale_fill_viridis_c(option = "magma") + 
  theme_bw()

ggplot(df_day) + 
  geom_sf(aes(fill = freq ), color = NA) +
  scale_fill_viridis_c(option = "magma") + 
  theme_bw()

hist(df_group$prec1)
rm(df_day)

# Compute R_M for each specie and each day ------------------------------
## Es el prec1 no el prec !! IMPORTANTEE!!!
df_group$prec <- as.numeric(df_group$prec1)
df_group[, R0_dai_alb := mapply(R0_func_alb, tmean, prec, dens)]
df_group[, R0_dai_aeg := mapply(R0_func_aeg, tmean, prec, dens)]
df_group[, R0_dai_jap := mapply(R0_func_jap, tmean, prec, freq)]

# Aggregate monthly ------------------------------------------------------
df_group_mon <- df_group[, .(tmean = mean(tmean),
                             tmin = min(tmean),
                             tmax = max(tmean),
                             prec = sum(prec), 
                             precmean = mean(prec), 
                             dens = min(dens),
                             freq = min(freq),
                             dens1 = max(dens),
                             R0_dai_alb = mean(R0_dai_alb),
                             R0_dai_aeg = mean(R0_dai_aeg),
                             R0_dai_jap = mean(R0_dai_jap)), 
                         by=list(NATCODE,month)]

df_group_mon[, R0_mon_alb := mapply(R0_func_alb, tmean, precmean, dens)]
df_group_mon[, R0_mon_aeg := mapply(R0_func_aeg, tmean, precmean, dens)]
df_group_mon[, R0_mon_jap := mapply(R0_func_jap, tmean, precmean, freq)]
df_group_mon[, R0_mon_alb_min := mapply(R0_func_alb, tmin, precmean, dens)]
df_group_mon[, R0_mon_aeg_min := mapply(R0_func_aeg, tmin, precmean, dens)]
df_group_mon[, R0_mon_jap_min := mapply(R0_func_jap, tmin, precmean, freq)]
df_group_mon[, R0_mon_alb_max := mapply(R0_func_alb, tmax, precmean, dens)]
df_group_mon[, R0_mon_aeg_max := mapply(R0_func_aeg, tmax, precmean, dens)]
df_group_mon[, R0_mon_jap_max := mapply(R0_func_jap, tmax, precmean, freq)]


# Maps number of months that > 1
df_group_mon$bool_R0_alb <- ifelse(df_group_mon$R0_mon_alb < 1,0,1)
df_group_mon$bool_R0_aeg <- ifelse(df_group_mon$R0_mon_aeg < 1,0,1)
df_group_mon$bool_R0_jap <- ifelse(df_group_mon$R0_mon_jap < 1,0,1)
df_group_mon$bool_R0_alb_min <- ifelse(df_group_mon$R0_mon_alb_min < 1,0,1)
df_group_mon$bool_R0_aeg_min <- ifelse(df_group_mon$R0_mon_aeg_min < 1,0,1)
df_group_mon$bool_R0_jap_min <- ifelse(df_group_mon$R0_mon_jap_min < 1,0,1)
df_group_mon$bool_R0_alb_max <- ifelse(df_group_mon$R0_mon_alb_max < 1,0,1)
df_group_mon$bool_R0_aeg_max <- ifelse(df_group_mon$R0_mon_aeg_max < 1,0,1)
df_group_mon$bool_R0_jap_max <- ifelse(df_group_mon$R0_mon_jap_max < 1,0,1)
df_group_mon$bool_R0_alb_dai <- ifelse(df_group_mon$R0_dai_alb < 1,0,1)
df_group_mon$bool_R0_aeg_dai <- ifelse(df_group_mon$R0_dai_aeg < 1,0,1)
df_group_mon$bool_R0_jap_dai <- ifelse(df_group_mon$R0_dai_jap < 1,0,1)

# Save DF -------------------------------------------------------------
Path <- paste0("~/INVASIBILITY_THRESHOLD/output/ERA5/temp/R0_clim_monthly",
               year,".Rds")
saveRDS(df_group_mon, Path)

## Test:
ggplot(df_group_mon) + geom_point(aes(tmean, bool_R0_jap))
df_group_mon[which(df_group_mon$bool_R0_jap == 0 &
                     df_group_mon$tmean < 15),]

# R0_func_jap(13.7,56,0.6)
# Plots months---------------------------------------------------------
library(RColorBrewer)
library(ggpubr)
library("latex2exp")
plot_months <- function(df, month){
  df1 <- df[which(df$month == month),]
  # Create a palette function using colorRampPalette
  plot <- ggplot(df1) +
    geom_sf(aes(fill = R0), colour = NA) +
    geom_sf(data = can_box) + coord_sf(datum = NA) +
    scale_fill_distiller(palette = "Spectral",
                         limits = c(min(df$R0),max(df$R0)),
                         name = TeX("$R_M$")) +
    ggtitle(as.character(month)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = c(0.1,0.6))
  return(plot)
}

# Esp can shapefile --------------------------------------------------
esp_can <- esp_get_munic_siane(moveCAN = TRUE)
can_box <- esp_get_can_box()
esp_can$NATCODE <- as.numeric(paste0("34",esp_can$codauto,
                                     esp_can$cpro,
                                     esp_can$LAU_CODE))

## Para hacer un cuadrado con seis plots cambio el numero de 
# month y el nombre del plot y lo hago para cada especie.
df_group_mon <- esp_can %>% left_join(df_group_mon, by = "NATCODE")
# df_group_mon$R0 <- df_group_mon$R0_mon_aeg
# month = 7
# plot_7 <- plot_months(df_group_mon,month)
# plot_7
# ggarr <- ggarrange(plot_3,plot_4,plot_5,
#           plot_6,plot_7,plot_8,
#           plot_9,plot_10,plot_11,
#           nrow=3,ncol = 3, common.legend = TRUE)
# 
# ggarr
# 
# Path <- paste0("~/Documentos/PHD/2023/INVASIBILITY/Plots/MS/Maps/2020/WholeyearJAP",year,".pdf")
# dev.copy2pdf(file=Path, width = 7, height = 5)
# ggsave(Path, plot = ggarr)

# Check the Guadalquivir region ------------------------------------------------
# unique(df_group_mon$ine.prov.name)
# df_jaen <- df_group_mon[which(df_group_mon$ine.prov.name == "Jaén"), ]
# df_jaen <- df_group_mon[which(df_jaen$name %in% unique(df_jaen$name)[1:10]), ]
# df_val <- df_group_mon[which(df_group_mon$ine.prov.name == "Valencia/València"), ]
# df_val <- df_val[which(df_val$name %in% unique(df_val$name)[1:10]), ]
# 
# ggarrange(ggplot(df_val) +
#             geom_line(aes(month, tmean), color = "red") +
#             geom_line(aes(month, precmean), color = "blue") +
#             geom_point(aes(month, R0_mon_jap), color = "brown", size = 0.5) +
#             # geom_point(aes(month, R0_mon_aeg), color = "black", size = 0.5) +
#             theme(legend.position = "none") +
#             geom_hline(yintercept = 14.3, color = "blue") +
#             theme_bw() +
#             geom_hline(yintercept = 15.2, color = "red")  +
#             theme(legend.position = "none") +ggtitle("Valencia"),
#           ggplot(df_jaen) +
#             geom_line(aes(month, tmean), color = "red") +
#             geom_line(aes(month, precmean), color = "blue") +
#             geom_point(aes(month, R0_mon_jap), color = "brown", size = 0.5) +
#             # geom_point(aes(month, R0_mon_aeg), color = "black", size = 0.5) +
#             geom_hline(yintercept = 14.3, color = "blue") +
#             geom_hline(yintercept = 15.2, color = "red") +
#             theme_bw() +
#             theme(legend.position = "none") + ggtitle("Jaen"))

# check rainfall
ggplot(df_group_mon[which(df_group_mon$month == 8),]) +
  geom_sf(aes(fill = precmean), color = NA) +
  scale_fill_viridis_c()

## ----------PLOT ANNUAL AVERAGE SEASON----------#
## Group by year:
df_group_mon$geometry <- NULL
df_group_mon <- setDT(df_group_mon)
df_group_y <- df_group_mon[,.(tmean = mean(tmean),
                              prec = mean(precmean),
                              dens =min(dens),
                              freq = min(freq),
                              R0_an_alb = mean(R0_mon_alb),
                              R0_an_aeg = mean(R0_mon_aeg),
                              R0_an_jap = mean(R0_mon_jap),
                              R0_sum_alb = sum(bool_R0_alb),
                              R0_sum_aeg = sum(bool_R0_aeg),
                              R0_sum_jap = sum(bool_R0_jap),
                              R0_sum_alb_min = sum(bool_R0_alb_min),
                              R0_sum_aeg_min = sum(bool_R0_aeg_min),
                              R0_sum_jap_min = sum(bool_R0_jap_min),
                              R0_sum_alb_max = sum(bool_R0_alb_max),
                              R0_sum_aeg_max = sum(bool_R0_aeg_max),
                              R0_sum_jap_max = sum(bool_R0_jap_max),
                              R0_sum_alb_dai = sum(bool_R0_alb_dai),
                              R0_sum_aeg_dai = sum(bool_R0_aeg_dai),
                              R0_sum_jap_dai = sum(bool_R0_jap_dai)),
                           by = list(NATCODE)]

df_group_y[, R0_anual_alb := mapply(R0_func_alb, tmean, prec, dens)]
df_group_y[, R0_anual_aeg := mapply(R0_func_aeg, tmean, prec, dens)]
df_group_y[, R0_anual_jap := mapply(R0_func_jap, tmean, prec, freq)]

# Save DF -------------------------------------------------------------
Path <- paste0("~/INVASIBILITY_THRESHOLD/output/ERA5/temp/2020/R0_clim_",
               year,".Rds")
saveRDS(df_group_y, Path)

# Plot R0 annual avg temp ---------------------------------------------
df_group_y <- esp_can %>% left_join(df_group_y)
ggplot(df_group_y) +
  geom_sf(aes(fill = R0_an_jap), colour = NA) +
  geom_sf(data = can_box) + coord_sf(datum = NA) +
  scale_fill_distiller(palette = "Spectral",
                       name = TeX("$R_M$")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

# Save RDS for validation ---------------------------------------------
Path <- paste0("~/INVASIBILITY_THRESHOLD/output/R0/datasets/R0_",year,".Rds")
saveRDS(df_group_y,Path)

# Whole map group by number of months suitable
library(ggpubr)
library(RColorBrewer)
name_pal = "RdYlBu"
display.brewer.pal(11, name_pal)
pal <- rev(brewer.pal(11, name_pal))
pal[11]
pal[12] = "#74011C"
pal[13] = "#4B0011"
letsize = 16
plot_summonths <- function(df){
  num_colors <- 13
  # Create a palette function using colorRampPalette
  my_palette <- colorRampPalette(c("#faf0ca","#f95738", "#732c2c"))
  
  colors <- c("#43A2CA", "#7BCCC4", "#BAE4BC", "#F0F9E8",
              "#FFF7EC","#FEE8C8","#FDD49E","#FDBB84",
              "#FC8D59","#EF6548","#D7301F", "#B30000",
              "#7F0000") 
  ggplot(df) +
    geom_sf(aes(fill = as.factor(R0)), colour = NA) +
    geom_sf(data = can_box) + coord_sf(datum = NA) +
    scale_fill_manual(values = pal,
                      name = "Nº months\n suitable",
                      limits = factor(seq(0,12,1))) +
    theme_bw() 
}

## Computed with sum months maps 
df_group_y$R0 <- df_group_y$R0_sum_alb
plot_sum_alb <- plot_summonths(df_group_y)
plot_sum_alb

Path <- paste0("~/Documentos/PHD/2023/INVASIBILITY/Plots/MS/Maps/",year,"/New_AlboSum",year,".pdf")
dev.copy2pdf(file=Path, width = 7, height = 5)

df_group_y$R0 <- df_group_y$R0_sum_aeg
plot_sum_aeg <- plot_summonths(df_group_y)
plot_sum_aeg
Path <- paste0("~/Documentos/PHD/2023/INVASIBILITY/Plots/MS/Maps/",year,"/New_Sum",year,".pdf")
dev.copy2pdf(file=Path, width = 7, height = 5)

df_group_y$R0 <- df_group_y$R0_sum_jap
plot_sum_jap <- plot_summonths(df_group_y)
plot_sum_jap
Path <- paste0("~/Documentos/PHD/2023/INVASIBILITY/Plots/MS/Maps/",year,"/New_JapSum",year,".pdf")
dev.copy2pdf(file=Path, width = 7, height = 5)

ggplot(df_group_y) +
  geom_sf(aes(fill = prec), color = NA) +
  scale_fill_viridis_c()

# Extract the legend
legend_only <- get_legend(plot_sum_alb +
                            theme(legend.position = "top"))
ggarrange(plot_sum_alb + ggtitle(expression(italic("Ae. albopictus")))+
            theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5)),
          plot_sum_aeg + ggtitle(expression(italic("Ae. aegypti")))+
            theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5)),
          plot_sum_jap + ggtitle(expression(italic("Ae. japonicus")))+
            theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5)),
          legend_only,
          ncol=2, nrow = 2)

## Test if variables make sense:
# df_group_y <- esp_can %>% left_join(df_group_y)
ggplot(df_group_y) + 
  geom_sf(aes(fill = tmean), color = NA) +
  scale_fill_viridis_c()

ggplot(df_group_y) + 
  geom_sf(aes(fill = prec), color = NA) +
  scale_fill_viridis_c()

ggplot(df_group_y) + 
  geom_sf(aes(fill = dens), color = NA) +
  scale_fill_viridis_c()
