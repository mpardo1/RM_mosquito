## Code to produce Figure 1 
rm(list=ls())
library(mapSpain)
library(ggplot2)
library("ggpubr")
library(data.table)
source("~/RM_mosquito/code/funcR0.R")

# R_M ---------------------------------------------------------------------
vec <- seq(5,40,0.001)
aegypti <- sapply(vec,R0_func_aeg, hum = 500,rain = 8)
albopictus <- sapply(vec,R0_func_alb, hum = 500,rain = 8) 

df_out <- data.frame(vec,
                     aegypti = aegypti,
                     albopictus = albopictus)
df_out <- reshape2::melt( df_out, id.vars = "vec")

# checks for the text :
esp = "albopictus"
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
  theme_bw() + theme(legend.position = c(0.2,0.75),
                     text = element_text(size = letsize),
                     legend.text.align = 0)

plot_temp

# rm as a function of rainfall ------------------------------------------
vec <- seq(0,16,0.001)
temp_opt <- 15
aegypti <- sapply(vec,R0_func_aeg, hum = 0, Te = temp_opt)
albopictus <- sapply(vec,R0_func_alb, hum = 0, Te = temp_opt)
df_rain <- data.frame(vec, albopictus, aegypti)
df_rain <- reshape2::melt(df_rain, id.vars = "vec")
plot_rain <- ggplot(df_rain) + 
  geom_line(aes(vec,value, color = variable), size = 1) +
  # geom_hline(yintercept = 1, linetype = "dashed", color = "red") + 
  scale_color_manual(name = "", values =pal,
                     labels = c(expression(italic("Ae. aegypti")),
                                expression(italic("Ae. albopictus")))) +
  # ,
  # expression(italic("Ae. japonicus")))) +
  xlab("Rainfall (mm)") + ylab(TeX("$R_M$")) + 
  theme_bw() + theme(text = element_text(size = letsize),
                     legend.position = "none")
plot_rain

# rm as a function of human density ------------------------------------------
vec <- seq(0,1000,0.1)
aegypti <- sapply(vec,R0_func_aeg, rain = 0, Te = temp_opt)
albopictus <- sapply(vec,R0_func_alb, rain = 0, Te = temp_opt)
df_hum <- data.frame(vec,albopictus, aegypti)
df_hum <- reshape2::melt(df_hum, id.vars = "vec")
plot_hum <- ggplot(df_hum) + 
  geom_line(aes(vec,value, color = variable), size = 1) +
  # geom_hline(yintercept = 1, linetype = "dashed", color = "red") + 
  scale_color_manual(name = "", values =pal,
                     labels = c(expression(italic("Ae. aegypti")),
                                expression(italic("Ae. albopictus")))) +
  # ,
  # expression(italic("Ae. japonicus")))) +
  xlab("Human density") + ylab(TeX("$R_M$")) + 
  theme_bw() + theme(text = element_text(size = letsize),
                     legend.position = "none")
plot_hum

# join all the plots ------------------------------------------------------
ggarrange(plot_temp + ggtitle("a)"),
          plot_rain + rremove("ylab")+ ggtitle("b)"),
          plot_hum + rremove("ylab")+ ggtitle("c)"),
          ncol = 3,
          widths = c(1,0.7,0.7))

