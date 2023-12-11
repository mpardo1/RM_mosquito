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

source("~/INVASIBILITY_THRESHOLD/code/funcR0.R")
## Data Catu ------------------------------------------------------
Path <- "~/INVASIBILITY_THRESHOLD/data/Datos_Catu/gi_min_model_pred.RData"
load(Path)
unique(gi_min_model_pred$trap_name)

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
trap_data[, R0_alb := mapply(R0_func_alb, mean_temperature, prec7, population)]
trap_data$female_norm <- trap_data$females/max(trap_data$females)
trap_data$pred_norm <- trap_data$pred/max(trap_data$pred)
trap_data$female_norm <- trap_data$females
trap_data$pred_norm <- trap_data$pred

# Compare model sampling effort versus count females ----------------
ggplot(trap_data) +
  geom_point(aes(start_date, female_norm), color = "red")+
  geom_point(aes(start_date, pred_norm), color = "blue")

# Normalized R0 between 0 and 1 to compare with counts --------------
trap_data$R0_alb_norm <- trap_data$R0_alb/max(trap_data$R0_alb)

# DF of specific cities ---------------------------------------------------
list_cit <- list("Blanes","Tordera", "Palafolls", "Lloret de Mar")
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

ggplot(trap_data_filt) +
  geom_point(aes(R0_alb_norm,female_norm, color = city))

# Fit exp city with nls -------------------------------------------------------
df_out_fem <- data.frame()
for(i in c(1:length(list_cit))){
  df_aux1 <- trap_data_filt[which(trap_data_filt$city == list_cit[[i]]),]
  Fitting_fem <- nls(female ~ exp(cont1*R0_alb)*cont,
                     data = df_aux1,
                     start = list(cont = 0.2, cont1 = 0.12))
  
  summary(Fitting_fem)
  
  mod <- function(te){
    t0 <- as.numeric(Fitting_fem$m$getPars()[1])
    tm <- as.numeric(Fitting_fem$m$getPars()[2])
    t0*exp(tm*te)
  }
  
  vec <- seq(0,max(df_aux1$R0_alb)+1,0.01)
  df_aux <- data.frame(temp_ae <- vec,
                       fem <- sapply(vec, mod))
  df_aux$cit <- list_cit[[i]]
  colnames(df_aux) <- c("vec", "out", "cit")
  df_out_fem <- rbind(df_aux, df_out_fem)
  colnames(df_out_fem) <- c("vec", "out", "cit")
}


name_pal = "Set1"
display.brewer.pal(length(list_cit), name_pal)
pal <- rev(brewer.pal(length(list_cit), name_pal))
sizelet = 14
plot_counts <- ggplot(data = trap_data_filt) +
  geom_point(aes(x = R0_alb, 
                 y = female, color = city), size = 0.8, alpha = 0.6) +
  geom_line(data = df_out_fem,
            aes(vec , out, color = cit),
            lwd = 0.8) +
  scale_color_manual(values = pal,
                     labels = c(TeX("Blanes"),
                                TeX("Lloret de Mar"),
                                TeX("Palafolls"),
                                TeX("Tordera")), name = "") +
  xlab(TeX("$R_M$")) +
  ylab("Number of females") +
  xlim(c(0,max(trap_data_filt$R0_alb)+0.3)) +
  ylim(c(0,max(trap_data_filt$female))) +
  theme_bw() + 
  theme(text = element_text(size = sizelet),
        legend.position = c(0.2,0.65))
plot_counts

# exponential function with lm---------------------------------------------------
df_out_fem <- data.frame()
for(i in c(1:length(list_cit))){
  df_aux1 <- trap_data_filt[which(trap_data_filt$city == list_cit[[i]]),]
  model <- lm(log(female)~ R0_alb, data =df_aux1)
  vec <- seq(0,max(df_aux1$R0_alb),0.01)
  df_out <- data.frame(vec = vec, 
                       # out =  exp(coef(model)[1])*(exp(df_aux1$R0_alb*coef(model)[2])))
                       out =  coef(model)[1] + vec*coef(model)[2])
  df_out$cit <- list_cit[[i]]
  df_out$rsq <- summary(model)$adj.r.squared
  colnames(df_out) <- c("vec", "out", "cit","rsq")
  df_out_fem <- rbind(df_out, df_out_fem)
  colnames(df_out_fem) <- c("vec", "out", "cit","rsq")
  ggplot(df_out)+
    geom_line(aes(vec,out)) +
    geom_point(data = df_aux1 ,
              aes(R0_alb,female)) +
    ylim(c(0,max(df_aux1$female)))
}

# Plot with dots and fit exp ---------------------------------------------
name_pal = "Set1"
display.brewer.pal(length(list_cit), name_pal)
pal <- rev(brewer.pal(length(list_cit), name_pal))
unique(df_out_fem[,c("cit", "rsq")])
sizelet = 14
plot_counts <- ggplot(data = trap_data_filt) +
  geom_point(aes(x = R0_alb, 
                 y = log(female), color = city), size = 0.8, alpha = 0.6) +
  geom_line(data = df_out_fem,
            aes(vec , out, color = cit),
            lwd = 0.8) +
  scale_color_manual(values = pal,
                     labels = c(TeX("Blanes, $R^2$: 0.37 "),
                                TeX("Lloret de Mar, $R^2$: 0.59"),
                                TeX("Palafolls, $R^2$: 0.32"),
                                TeX("Tordera, $R^2$: 0.31")),
                     name = "") +
  xlab(TeX("$R_M$")) +
  ylab("Number of females (ln)") +
  # xlim(c(0,max(trap_data_filt$R0_alb)+0.1)) +
  # ylim(c(0,max(trap_data_filt$female))) +
  theme_bw() + 
  theme(text = element_text(size = sizelet),
        legend.position = c(0.25,0.8))
plot_counts

# map spain 2020 Rm albopictus:
year = 2020
Path <- paste0("~/INVASIBILITY_THRESHOLD/output/ERA5/temp/2020/R0_clim_",
               year,".Rds")
df_2020 <- setDT(readRDS(Path))
df_2020 <- df_2020[,c("NATCODE", "R0_sum_alb", "R0_sum_aeg", "R0_sum_jap")]
colnames(df_2020) <-c ("NATCODE", "Alb_2020", "Aeg_2020", "Jap_2020")

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

# Load Data --------------------------------------------------------------------
## year = 2022, we will use this year to validate the data
year = 2020
Path <- paste0("~/INVASIBILITY_THRESHOLD/output/R0/datasets/R0_",year,".Rds")
# check with no h (no rainfall and human density) validation
# Path <- paste0("INVASIBILITY_THRESHOLD/output/R0/datasets/R0_H2020.Rds")
df_group_tot <- readRDS(Path)
head(df_group_tot)

## Data Catu BG traps
Path <- "~/INVASIBILITY_THRESHOLD/data/Datos_Catu/gi_min_model_pred.RData"
load(Path)
unique(gi_min_model_pred$trap_name)

# Presence absence data Albopictus Spain:
Path = "~/INVASIBILITY_THRESHOLD/data/PA/Albopictus_Spain_Pa.csv"
df_pa <- read.csv(Path)
# add cadiz
df_pa[which(df_pa$NATCODE == "34011111012"), "PA"] <- 1
df_pa[which(df_pa$name == "Moaña"), "PA"] <- 1
df_pa[which(df_pa$name == "Vigo"), "PA"] <- 1
df_pa[which(df_pa$name == "Redondela"), "PA"] <- 1

# Read csv OCT PA
Path = "~/INVASIBILITY_THRESHOLD/data/PA/ALBO_PA_2023_OCT.csv"
df_pa <- read.csv(Path)
df_pa$PA <- ifelse(df_pa$A_PRIM_DET_OFICIAL != 0 |
                     df_pa$A_PRIM_DET_CITSCI != 0,1,0)
write.csv(df_pa, "~/INVASIBILITY_THRESHOLD/data/PA/ALBO_PA_102023.csv")

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

df_group_m <- df_group_tot[,c("NATCODE",
                              "R0_sum_alb",
                              "R0_sum_alb_min",
                              "R0_sum_alb_max",
                              "R0_sum_alb_dai")]
df_pa <- df_pa %>% left_join(df_group_m)

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

unique(esp_can$ine.ccaa.name)

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

NATCODE_CAT <- esp_can[which(esp_can$ine.ccaa.name %in% list_ccaa ),"NATCODE"]

# Add PA data --------------------------------------------------------
NATCODE_CAT$geometry <- NULL
df_pa_CAT <- df_pa[which(as.numeric(df_pa$NATCODE) %in% 
                           as.numeric(NATCODE_CAT$NATCODE)),]

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
  theme(legend.position = c(0.17,0.6),
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

# Plot all ccaa ----------------------------------------
df_aux <- df_num_months(df_pa)
df_group <- df_aux %>% group_by(R0_sum_alb) %>%
  summarise(num_1 = sum(num_1),
            num_0 = sum(num_0),
            prop_1 = sum(num_1)/(sum(num_0)+sum(num_1)),
            sum_muni = sum(sum_muni))
sum(df_group$sum_muni)
plot_esp <- ggplot(df_group) +
  geom_line(aes(R0_sum_alb,prop_1)) +
  geom_point(aes(R0_sum_alb,prop_1,
                 size = sum_muni), alpha = 0.6) +
  xlab("Nº months suitable") + 
  ylab("Proportion of municipalities with presence") +
  ylim(c(0,1)) + 
  scale_size_continuous(name = "Number municipalities",
                        breaks = c(100,500,1000,2000,3000, 4000),
                        labels =c("[0,100]","[100,500]",
                                  "[500,1000]","[1000,2000]",
                                  "[2000, 3000]",
                                  "[3000, 4000]")) +
  scale_x_continuous(breaks = seq(1,12,1)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = c(0.15,0.6)) 
plot_esp

# Join min, max, mean toguether ---------------------------------------
NATCODE_CAT$geometry <- NULL
df_pa_CAT <- df_pa[which(df_pa$ine.ccaa.name == "Cataluña"),]

library(ggpubr)
leg <- get_legend(plot_2020)
as_ggplot(leg)
ggarr <- ggarrange( PA_alb+ ggtitle("A"),
                    plot_2020 + ggtitle("B") +
                      theme(legend.position = "none"),
                    leg,
                    ncol = 3,
                   nrow = 1, widths = c(1,1,0.2))

ggarr1 <- ggarrange(plot_ccaa+ ggtitle("C"),
                   plot_counts + ggtitle("D"),
                   widths = c(1.2,1),
                   ncol = 2, nrow = 1)
ggarrange(ggarr, ggarr1, ncol = 1)

# Other way to plot -----------------------------------------
ggarrange(ggarr ,
          plot_ccaa+ ggtitle("C"),
          plot_counts + ggtitle("D"),
          ncol = 1, nrow = 3, heights = c(1,1,1))

ggarrange(ggarr,
          plot_ccaa + ylab("Proportion presence")
          + ggtitle("C"),
          plot_min_max + ylab("Proportion presence")
          + ggtitle("D"),
          ncol = 1, nrow = 3)

