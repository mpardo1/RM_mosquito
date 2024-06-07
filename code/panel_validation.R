### Code to do a validation between female traps counts and R0
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
library(RColorBrewer)
library("latex2exp")
source("~/RM_mosquito/code/funcR0.R")

# Load BG count data
Path <- "~/RM_mosquito/data/BG_count_data_albopictus.RData"
load(Path)
city_count <- gi_min_model_pred %>% group_by(city) %>% 
  summarise(n=n())
# Process data from Catu with trap data ----------------------------
trap_data <- setDT(gi_min_model_pred[,c("trap_name", "province",
                                        "city", "start_date","end_date",
                                        "females","precipitation", "mean_temperature",
                                        "population", "pred", "l7precipitation",
                                        "l14precipitation", "l21precipitation", 
                                        "trapping_effort", "pred")])
rm(gi_min_model_pred)
trap_data$prec7 <- trap_data$l7precipitation/7
trap_data$prec14 <- trap_data$l14precipitation/14
trap_data$trapping_effort <- as.numeric(trap_data$trapping_effort)
trap_data$females_daily <- trap_data$females/trap_data$trapping_effort
trap_data[which(trap_data$city == "La Bisbal de l'Empordà"), "population"] <- 10859
trap_data[, R0_alb := mapply(R0_func_alb, mean_temperature, prec14, population)]
trap_data$female_norm <- trap_data$females/max(trap_data$females)
trap_data$pred_norm <- trap_data$pred/max(trap_data$pred)
trap_data$female_norm <- trap_data$females
trap_data$pred_norm <- trap_data$pred

# Normalized R0 between 0 and 1 to compare with counts --------------
trap_data$R0_alb_norm <- trap_data$R0_alb/max(trap_data$R0_alb)

# DF of specific cities ---------------------------------------------------
list_cit <- unique(trap_data$city)
# list_cit <- list("Blanes","Tordera", "Palafolls", "Lloret de Mar")
trap_data_filt <- trap_data[which(trap_data$city %in% list_cit)] %>%
  group_by(city, start_date) %>% 
  summarise(female = sum(females),
            R0_alb = mean(R0_alb))

# DF compute maximum -------------------------------------------------------
trap_citi_max <- trap_data_filt %>%
  group_by(city) %>% 
  summarise(female_max= max(female),
            R0_alb_max = max(R0_alb))

# DF compute norm female and R0  --------------------------------------------
trap_data_filt <- trap_data_filt %>% left_join(trap_citi_max)
trap_data_filt$female_norm <- trap_data_filt$female/trap_data_filt$female_max
trap_data_filt$R0_alb_norm <- trap_data_filt$R0_alb/trap_data_filt$R0_alb_max

ggplot(trap_data_filt[trap_data_filt$city == list_cit[[1]],]) +
  geom_point(aes(R0_alb,female, color = city))

# Plot with geom_smooth 
library("latex2exp")
list_cit_1 <- list_cit[c(1,2,3,4)]
name_pal = "Set1"
display.brewer.pal(length(list_cit), name_pal)
pal <- brewer.pal(length(list_cit), name_pal)
plot_corr <- ggplot(trap_data_filt[trap_data_filt$city %in% list_cit_1,],
                    aes(R0_alb, log10(female), color = city))+
  geom_point(alpha = 0.7) + geom_smooth(method = "lm", alpha = 0.3) +
  scale_color_manual(values = pal) + theme_bw() + xlab(TeX("$R_M$")) +
  theme(legend.position = c(0.2,0.8), text = element_text(size = 14)) +
  ylab("Count female mosquito trap (logarithmic scale)")
plot_corr

# Map spain 2020 Rm albopictus:
Path <- "~/RM_mosquito/data/R0_avg_2003-2020.Rds"
df_2020 <- setDT(readRDS(Path))
df_2020 <- df_2020[,c("NATCODE", "sum_alb")]
colnames(df_2020) <-c ("NATCODE", "Alb_2020")

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
  library(mapSpain)
  esp_can <- esp_get_munic_siane(moveCAN = TRUE)
  can_box <- esp_get_can_box()
  esp_can$NATCODE <- as.numeric(paste0("34",esp_can$codauto,
                                       esp_can$cpro,
                                       esp_can$LAU_CODE))
  df <- esp_can %>% left_join(df)
  num_colors <- 13
  # Create a palette function using colorRampPalette
  my_palette <- colorRampPalette(c("#faf0ca","#f95738", "#732c2c"))
  
  colors <- c("#43A2CA", "#7BCCC4", "#BAE4BC", "#F0F9E8",
              "#FFF7EC","#FEE8C8","#FDD49E","#FDBB84",
              "#FC8D59","#EF6548","#D7301F", "#B30000",
              "#7F0000") 
  ggplot(df) +
    geom_sf(aes(fill = as.factor(R0)),
            colour = NA) +
    geom_sf(data = can_box) +
    coord_sf(datum = NA) +
    scale_fill_manual(values = pal,
                      name = "Nº suitable \n months",
                      limits = factor(seq(0,12,1))) +
    theme_minimal()  +
    theme(legend.position = "right",
          legend.text = element_text(14)) 
}

# Albopictus ---------------------------------------------------------
# 2004
df_2020$R0 <- df_2020$Alb_2020
plot_2020 <- plot_summonths(df_2020)
plot_2020
df_2020$R0 <- NULL

# Presence absence data Albopictus Spain:
Path <- "~/RM_mosquito/data/PA_Oct_2023_albopictus.csv"
df_pa <- read.csv(Path)

# Map Spain --------------------------------------------------------------------
esp_can <- esp_get_munic_siane(moveCAN = TRUE)
esp_can$NATCODE <- as.numeric(paste0("34",
                                     esp_can$codauto,
                                     esp_can$cpro,
                                     esp_can$LAU_CODE))
can_box <- esp_get_can_box()

# Plot PA Albopictus -----------------------------------------------------------
prov_esp <- esp_get_prov_siane()
df_pa <- esp_can %>% left_join(df_pa)
gir_star <- as.data.frame(st_centroid(esp_can[which(esp_can$name == "Blanes"),"geometry"]))
lonlat_gir <- gir_star %>%
  mutate(long = unlist(map(gir_star$geometry,1)),
         lat = unlist(map(gir_star$geometry,2)))

# map ccaa spain -----------------------------------------------------------
ccaa <- esp_get_ccaa(ccaa = c(
  "Catalunya",
  "Comunidad Valenciana",
  "Aragón",
  "País Vasco","Andalucía"
))

# create palette -------------------------------------------------------------
name_pal = "Set1"
display.brewer.pal(7, name_pal)
pal <- brewer.pal(7, name_pal)
pal <- pal[c(1,3:5,7)]
pal[length(pal)] <- "#08B2E3"
pal[4] <- "#ECA400"
pal <- c("#301A4B", "#ECA400", "#B3001B", "#08B2E3", "#79B473")
non_pa_col <- "#EBEBEB"
pa_col <-  "#8CAADF"
df_pa[is.na(df_pa$PA),"PA"] <- 0
PA_alb <- ggplot() +
  geom_sf(data =df_pa, aes(fill = as.factor(PA)),
          # lwd = 0.03, alpha = 0.7, color = "#6E6E6E") +
          color = NA, alpha = 0.7) +
  geom_sf(data = ccaa,
          aes(fill = codauto, color = codauto),
          alpha = 0, lwd = 0.5) +
  geom_sf(data = can_box) + coord_sf(datum = NA) +
  scale_color_manual(values = pal, name = " ",guide = "none") +
  theme_minimal() +
  scale_fill_manual(values = c(non_pa_col,rep(pa_col,6)),
                    name = " ", limits = c("0","1")) +
  geom_point(data = lonlat_gir, aes(long,lat),
             color = "red", shape = 8) +
  rremove("xlab") + rremove("ylab") +
  theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
        legend.position = c(0.1,0.8)) 
PA_alb

### Comparison between presence absence and number of months R0>1
df_pa <- df_pa[, c("NATCODE", "PA")]
df_pa$geometry <- NULL
df_group_m <- df_2020[,c("NATCODE",
                              "Alb_2020")]
colnames(df_group_m) <- c("NATCODE", "R0_sum_alb")
df_pa <- df_pa %>% left_join(df_group_m)

# Extract data for CAT
NATCODE_CAT <- esp_can[which(esp_can$ine.ccaa.name == "Cataluña"),"NATCODE"]
NATCODE_CAT$geometry <- NULL
df_pa_CAT <- df_pa[which(as.numeric(df_pa$NATCODE) %in% 
                           as.numeric(NATCODE_CAT$NATCODE)),]

# Prop presence df ------------------------------------------
df_num_months <- function(df_pa_CAT){
  df_r_1 <- df_pa_CAT[which(df_pa_CAT$PA == 1),] %>% 
    group_by(R0_sum_alb) %>% summarize(num_1 = n())
  df_r_0 <- df_pa_CAT[which(df_pa_CAT$PA == 0),] %>% 
    group_by(R0_sum_alb) %>% summarize(num_0 = n())
  
  df_r <- merge(df_r_1 ,df_r_0,
                all.x = TRUE, all.y = TRUE)
  df_r$prop_1 <- ifelse(is.na(df_r$num_0),1,
                        ifelse(is.na(df_r$num_1),
                               0,df_r$num_1/(df_r$num_1+df_r$num_0)))
  df_r$num_1 <- ifelse(is.na(df_r$num_1),0,df_r$num_1)
  df_r$num_0 <- ifelse(is.na(df_r$num_0),0,df_r$num_0)
  df_r$sum_muni <- df_r$num_1 + df_r$num_0
  return(df_r)
  
}

# Filter df -----------------------------------------------------
df_pa <- df_pa[, c("NATCODE", "PA")]
df_pa$geometry <- NULL

# Join map with PA data -----------------------------------------
df_group_m <- df_group_tot
df_pa <- df_pa %>% left_join(df_group_m)

# compute a boxplot For a specific region or set of regions
list_ccaa = c("Cataluña" ,  
              "Comunitat Valenciana", 
              "País Vasco","Andalucía","Aragón")
# list_ccaa = c("Cataluña")
# df_pa_filt <- df_pa[which(esp_can$ine.ccaa.name %in% list_ccaa ),]
# ggplot(df_pa_filt) +
#   geom_violin(aes(as.factor(PA), R0_sum_alb, fill = as.factor(PA)), 
#                                              alpha = 0.3) +
#   geom_boxplot(aes(as.factor(PA), R0_sum_alb, fill = as.factor(PA)), alpha = 0.2,
#                width = 0.2) +
#   ylab("Sum of suitable months") +
#   xlab("") + 
#   scale_fill_viridis_d(name = "") + 
#   scale_color_viridis_d(name = "") +
#   theme_bw()
# 
# wc <- wilcox.test(df_pa_filt[which(df_pa_filt$PA == 0),"R0_sum_alb"],
#                   df_pa_filt[which(df_pa_filt$PA == 1),"R0_sum_alb"])
# 
# # For whole Spain
# ggplot(df_pa) +
#   # geom_violin(aes(as.factor(PA), R0_sum_alb, fill = as.factor(PA)), 
#   #             alpha = 0.3) +
#   geom_boxplot(aes(as.factor(PA), R0_sum_alb, fill = as.factor(PA)), alpha = 0.2,
#                width = 0.2) +
#   geom_jitter(aes(as.factor(PA), R0_sum_alb, fill = as.factor(PA)), alpha = 0.2,
#                width = 0.2) +
#   ylab("Sum of suitable months") +
#   xlab("") + 
#   scale_fill_viridis_d(name = "") + 
#   scale_color_viridis_d(name = "") +
#   theme_bw()
# 
# wc <- wilcox.test(df_pa[which(df_pa$PA == 0),"R0_sum_alb"],
#                   df_pa[which(df_pa$PA == 1),"R0_sum_alb"])
# 
# unique(esp_can$ine.ccaa.name)

# Check particular ccaa
# colors <- c("#43A2CA", "#7BCCC4", "#BAE4BC", "#F0F9E8",
#             "#FFF7EC","#FEE8C8","#FDD49E","#FDBB84",
#             "#FC8D59","#EF6548","#D7301F", "#B30000",
#             "#7F0000") 
# df_pa_aux <- df_pa_CAT
# df_pa_aux <- esp_can %>% left_join(df_pa_CAT)
# unique(df_pa_CAT$ine.ccaa.name)
# ggplot(df_pa_aux[df_pa_aux$ine.ccaa.name == "País Vasco",])+
#   geom_sf(aes(fill = as.factor(R0_sum_alb)), color = NA) +
#   geom_sf(aes(color = as.factor(PA)), alpha= 0.3) +
#   scale_fill_manual(values=colors)

# Func to compute plot PA prop vs summonths --------------------
plot_sum_p <- function(ccaa){
  NATCODE_CAT <- esp_can[which(esp_can$ine.ccaa.name == ccaa),"NATCODE"]
  NATCODE_CAT$geometry <- NULL
  df_pa_CAT <- df_pa[which(as.numeric(df_pa$NATCODE) %in% 
                             as.numeric(NATCODE_CAT$NATCODE)),]
  
  df_sum_CAT <- df_num_months(df_pa_CAT)
  df_sum_CAT$ccaa <- ccaa
  
  ggplot(df_sum_CAT) +
    geom_point(aes(R0_sum_alb,prop_1)) +
    xlab("Nº months suitable") + 
    ylab("Proportion of municipalities with presence") +
    ylim(c(0,1)) +
    ggtitle(ccaa) +    
    scale_x_continuous(breaks = seq(1,12,1)) +
    theme_bw()
}

# DF with ccaa names ----------------------------------------------
esp_can_ccaa <- esp_can[,c("NATCODE", "ine.ccaa.name")]
df_pa_ccaa <- esp_can %>% left_join(df_pa)

# Join more than one ccaa ---------------------------------------------
list_ccaa = c("Cataluña" ,  
              "Comunitat Valenciana", 
              "País Vasco","Andalucía","Aragón")

NATCODE_CAT <- esp_can[which(esp_can$ine.ccaa.name %in% list_ccaa ),c("NATCODE", "ine.ccaa.name")]

# Add PA data --------------------------------------------------------
NATCODE_CAT$geometry <- NULL
df_pa_CAT <- df_pa[which(as.numeric(df_pa$NATCODE) %in% 
                           as.numeric(NATCODE_CAT$NATCODE)),] %>% left_join(NATCODE_CAT)

# Compute df with prop PA --------------------------------------------
df_sum_CAT <- data.frame()
for(i in c(1:length(list_ccaa))){
  df_aux <- df_num_months(df_pa_CAT[which(df_pa_CAT$ine.ccaa.name == list_ccaa[i] ),])
  df_aux$ccaa_n <- list_ccaa[i]
  df_sum_CAT <- rbind(df_aux,df_sum_CAT)
}

# Plot color related to ccaa ----------------------------------------
name_pal = "Set1"
display.brewer.pal(length(list_ccaa), name_pal)
# pal <- rev(brewer.pal(length(list_ccaa), name_pal))
plot_ccaa <- ggplot(df_sum_CAT) +
  geom_line(aes(R0_sum_alb,prop_1, color =ccaa_n)) +
  geom_point(aes(R0_sum_alb,prop_1, color =ccaa_n,
                 size = sum_muni), alpha = 0.6) +
  xlab("Nº months suitable") + 
  ylab("Proportion of municipalities with presence") +
  ylim(c(0,1)) + 
  scale_color_manual(name = "", values = pal,
                     labels = c("Andalusia","Aragon","Catalonia" ,  
                                "Valencian C.", 
                                "Basque C.")) +
  scale_size_continuous(name = "Nº municipalities",
                        breaks = c(5,100,200),
                        labels =c("<5","[5,100]",
                                  ">100")) +
  scale_x_continuous(breaks = seq(1,12,1)) +
  theme_bw() +
  theme(legend.position = c(0.12,0.72),
        text = element_text(size = 14)) 
plot_ccaa

# Plot all ccaa selected together ----------------------------------------
df_group <- df_sum_CAT %>% group_by(R0_sum_alb) %>%
  summarise(num_1 = sum(num_1),
            num_0 = sum(num_0),
            prop_1 = sum(num_1)/(sum(num_0)+sum(num_1)),
            sum_muni = sum(sum_muni))
plot_tot <- ggplot(df_group) +
  geom_line(aes(R0_sum_alb,prop_1)) +
  geom_point(aes(R0_sum_alb,prop_1,
                 size = sum_muni), alpha = 0.6) +
  xlab("Nº months suitable") + 
  ylab("Proportion of municipalities with presence") +
  ylim(c(0,1)) + 
  scale_size_continuous(name = "Number municipalities",
                        breaks = c(5,100,500,1000,1500),
                        labels =c("[0,5]","[5,100]",
                                  "[100,500]","[500,1000]",
                                  "[1000, 1500]")) +
  scale_x_continuous(breaks = seq(1,12,1)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = c(0.2,0.6)) 
plot_tot
