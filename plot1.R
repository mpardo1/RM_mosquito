## Code that Compare the PA data for albopictus in comparison with
# the number o months in which R0>1 and the avg R0
rm(list=ls())
library(mapSpain)
library(ggplot2)
library("ggpubr")
library(data.table)
source("~/INVASIBILITY_THRESHOLD/code/funcR0.R")

# R_M ---------------------------------------------------------------------
vec <- seq(5,40,0.001)
aegypti <- sapply(vec,R0_func_aeg, hum = 500,rain = 8)
albopictus <- sapply(vec,R0_func_alb, hum = 500,rain = 8) 
japonicus <- sapply(vec,R0_func_jap, hum = 500,rain =8) 


df_out <- data.frame(vec, aegypti = aegypti,
                     albopictus = albopictus)
                     # ,
                     # japonicus = japonicus)
df_out <- reshape2::melt( df_out, id.vars = "vec")

# checks for the text :
esp = "aegypti"
min(df_out[which(df_out$variable == esp &
                   df_out$value >1), "vec"])
max(df_out[which(df_out$variable == esp &
                   df_out$value >1), "vec"])
max_r <- max(df_out[which(df_out$variable == esp &
                            df_out$value >1), "value"])
df_out[which(df_out$value == max_r), "vec"]

library(RColorBrewer)
name_pal = "Set1"
display.brewer.pal(3, name_pal)
pal <- brewer.pal(3, name_pal)[2:3]
letsize = 16
library("latex2exp")
plot_temp <- ggplot(df_out) + 
  geom_line(aes(vec,value, color=variable), size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  ylab(TeX("$R_M$")) + 
  scale_color_manual(name = "", values =pal,
                     labels = c(expression(italic("Ae. aegypti")),
                                expression(italic("Ae. albopictus")))) +
                                # ,
                                # expression(italic("Ae. japonicus")))) +
  xlab("Temperature (Cº)") +
  scale_x_continuous(breaks = seq(5,41,4)) +
  theme_bw() + theme(legend.position = c(0.18,0.75),
                     text = element_text(size = letsize),
                     legend.text.align = 0)
# 
ggplot(df_out[which(df_out$variable == "albopictus"),]) +
  geom_line(aes(vec,value), size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  ylab(TeX("$R_M$")) + theme_bw() +
  theme(legend.position = c(0.18,0.75),
        text = element_text(size = letsize),
        legend.text.align = 0) +
  xlab("Temperature (Cº)") + 
  ggtitle(expression(italic("Ae. albopictus")))


# rm as a function of rainfall ------------------------------------------
vec <- seq(0,16,0.001)
rain <- sapply(vec,R0_func_aeg, hum = 0, Te = 17.5)
df_rain <- data.frame(vec, rain)
plot_rain <- ggplot(df_rain) + 
  geom_line(aes(vec,rain)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  xlab("Rainfall (mm)") + ylab(TeX("$R_M$")) + 
  theme_bw() + theme(text = element_text(size = letsize))

# rm as a function of human density ------------------------------------------
vec <- seq(0,1000,0.1)
hum <- sapply(vec,R0_func_aeg, rain = 0, Te = 17.5)
df_hum <- data.frame(vec, hum)
plot_hum <- ggplot(df_hum) + 
  geom_line(aes(vec,hum)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  xlab("Human density") + ylab(TeX("$R_M$")) + 
  theme_bw() + theme(text = element_text(size = letsize))

# join all the plots ------------------------------------------------------
ggarrange(plot_temp + ggtitle("A"),
          plot_rain + rremove("ylab")+ ggtitle("B"),
          plot_hum + rremove("ylab")+ ggtitle("C"),
          ncol = 3,
          widths = c(1,0.7,0.7))


# check one parameter dependecy ------------------------------------------
# vec <- seq(0,40,0.1)
# hum <- sapply(vec,a_f_alb)
# df_hum <- data.frame(vec, hum)
# ggplot(df_hum) + 
#   geom_line(aes(vec,hum)) +
#   xlab("Human density") + ylab(TeX("$R_M$")) + 
#   theme_bw() + theme(text = element_text(size = letsize))

# compute phase space for rainfall
library(data.table)
temp <- seq(0,40,length.out = 1000)
rain <- seq(0,16,length.out = 1000)
df_clim <- setDT(expand.grid(temp =temp, rain = rain, hum = 0))

df_clim[, R0_alb := mapply(R0_func_alb, temp, rain, hum)]
df_clim[, R0_aeg := mapply(R0_func_aeg, temp, rain, hum)]

library(latex2exp)
max(df_clim$R0_alb)
alb <- ggplot(df_clim,aes(temp, rain, fill = R0_alb)) +
  geom_raster() + 
  geom_contour(aes(z = R0_alb),breaks = 1,
               color = "black", linetype = "dashed") +
  scale_fill_distiller(palette = "Spectral",
                       name = TeX("$R_M$"), 
                       limits = c(0, max(df_clim$R0_alb))) +
  xlab("Temperature") + ylab("Rainfall") +
  xlim(c(10,38)) +
  theme_bw() + theme(text = element_text(size = letsize))


aeg <- ggplot(df_clim,aes(temp, rain, fill = R0_aeg)) +
  geom_raster() + 
  geom_contour(aes(z = R0_aeg),breaks = 1,
               color = "black", linetype = "dashed") +
  scale_fill_distiller(palette = "Spectral",
                       name = TeX("$R_M$"), 
                       limits = c(0, max(df_clim$R0_alb))) +
  xlab("Temperature") + ylab("Rainfall") +
  xlim(c(10,38)) +
  theme_bw() +  theme(text = element_text(size = letsize))

# compute phase space for human density
temp <- seq(0,40,length.out = 1000)
hum <- seq(0,800,length.out = 1000)
df_clim <- setDT(expand.grid(temp =temp, rain = 0, hum = hum))

df_clim[, R0_alb := mapply(R0_func_alb, temp, rain, hum)]
df_clim[, R0_aeg := mapply(R0_func_aeg, temp, rain, hum)]

library(latex2exp)
max(df_clim$R0_alb)
alb2 <- ggplot(df_clim,aes(temp, hum, fill = R0_alb)) +
  geom_raster() + 
  geom_contour(aes(z = R0_alb),breaks = 1,
               color = "black", linetype = "dashed") +
  scale_fill_distiller(palette = "Spectral",
                       name = TeX("$R_M$"), 
                       limits = c(0, max(df_clim$R0_alb))) +
  xlab("Temperature") + ylab("Human density") +
  xlim(c(10,38)) +
  theme_bw() + theme(text = element_text(size = letsize))


aeg2 <- ggplot(df_clim,aes(temp, hum, fill = R0_aeg)) +
  geom_raster() + 
  geom_contour(aes(z = R0_aeg),breaks = 1,
               color = "black", linetype = "dashed") +
  scale_fill_distiller(palette = "Spectral",
                       name = TeX("$R_M$"), 
                       limits = c(0, max(df_clim$R0_alb))) +
  xlab("Temperature") + ylab("Human density") +
  xlim(c(10,38)) +
  theme_bw() +  theme(text = element_text(size = letsize))
ggarrange(alb,aeg,alb2, aeg2, nrow = 2, ncol = 2,
          common.legend = TRUE)

# check if the maximun does not move
temp <- seq(0,40,length.out = 100)
hum <- seq(0,800,length.out = 100)
rain <- seq(0,16,length.out = 100)
df_clim <- setDT(expand.grid(temp =temp, rain = rain, hum = hum))

df_clim[, R0_alb := mapply(R0_func_alb, temp, rain, hum)]
df_clim[, R0_aeg := mapply(R0_func_aeg, temp, rain, hum)]

df_clim[which(df_clim$R0_alb == max(df_clim$R0_alb)), "temp"]
df_clim[which(df_clim$R0_aeg == max(df_clim$R0_aeg)), "temp"]

# test the influence of the constant e0, the one that weight the ingluence of 
# rainfall and human density
# R0 function by temperature:
R0_func_alb <- function(Te, rain, hum, erat){
  # Constants:
  e0 = 1.5
  evar = 0.05
  eopt = 8
  efac = 0.01
  edens = 0.01
  
  h <- (1-erat)*(((1+e0)*exp(-evar*(rain-eopt)^2))/(exp(-evar*(rain - eopt)^2) + e0)) +
    erat*(edens/(edens + exp(-efac*hum)))
  if(is.na(Te) | is.na(rain) | is.na(hum)){
    R0 <- NA
  }else{
    a <- a_f_alb(Te)
    f <- (1/2)*TFD_f_alb(Te)
    deltaa <- lf_f_alb(Te)
    dE <- dE_f_alb(Te)
    probla <- pLA_f_alb(Te)
    deltaE = 0.1
    R0 <- ((f*a*deltaa)*probla*(h*dE/(h*dE+deltaE)))^(1/3)
  }
  return(R0)
}

erat_test <- seq(0,1,0.01)
tempcte =  20
dt_erat <- data.table(rain = 8,hum = 0,Te = tempcte, erat =erat_test)
dt_erat1 <- data.table(rain = 8,hum = 100,Te = tempcte, erat =erat_test)
dt_erat2 <- data.table(rain = 8,hum = 250,Te = tempcte, erat =erat_test)
dt_erat3 <- data.table(rain = 8,hum = 500,Te = tempcte, erat =erat_test)
dt_erat <- rbind(dt_erat,dt_erat1,dt_erat2,dt_erat3)
dt_erat[, R0_alb := mapply(R0_func_alb, Te, rain, hum,erat)]
ggplot(dt_erat) + 
  geom_point(aes(erat,R0_alb , color = as.factor(hum))) +
  theme_bw()

dt_erat <- data.table(rain = 0,hum = 500,Te = tempcte, erat =erat_test)
dt_erat1 <- data.table(rain = 4,hum = 500,Te = tempcte, erat =erat_test)
dt_erat2 <- data.table(rain = 8,hum = 500,Te = tempcte, erat =erat_test)
dt_erat3 <- data.table(rain = 12,hum = 500,Te = tempcte, erat =erat_test)
dt_erat <- rbind(dt_erat,dt_erat1,dt_erat2,dt_erat3)
dt_erat[, R0_alb := mapply(R0_func_alb, Te, rain, hum,erat)]
ggplot(dt_erat) + 
  geom_point(aes(erat,R0_alb , color = as.factor(rain))) +
  theme_bw()
