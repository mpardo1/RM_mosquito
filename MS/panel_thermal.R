## Compute the thermal responses for Aedes Japonicus from literature data
# compare also different functions and compute its AIC value
rm(list= ls())
library(thermPerf)
library(ggplot2)
library(tidyverse)
library(nls2)

## Data taken from https://parasitesandvectors.biomedcentral.com/articles/10.1186/s13071-018-2659-1
## https://www.researchgate.net/publication/235430511_The_ecology_of_the_exotic_mosquito_Ochlerotatus_Finlaya_japonicus_japonicus_Theobald_1901_Diptera_Culicidae_and_an_examination_of_its_role_in_the_West_Nile_virus_cycle_in_New_Jersey
Path <- "~/INVASIBILITY_THRESHOLD/data/japonicus/adult_larva_lifespan.csv"
Japonicus <- read.csv(Path)

head(Japonicus)
Japonicus$lifespan <- Japonicus$Age_adult_death_mean_female  
plot_deltaA <- ggplot(Japonicus) + 
  geom_point(aes(Temp,lifespan)) + theme_bw()
plot_deltaA

## Linear Fit
Fitting_deltaA <- nls(lifespan ~ cont*Temp + cont1,
                      data = Japonicus,
                      start = list(cont = 0.001, cont1 = 0.0))
summary(Fitting_deltaA)

mod <- function(te){
  c <- as.numeric(Fitting_deltaA$m$getPars()[1])
  c1 <- as.numeric(Fitting_deltaA$m$getPars()[2])
  c*te+c1
}
# Quadratic Fit
Fitting_deltaA_quad <- nls(lifespan ~ cont*Temp^2 + cont1*Temp +cont2,
                           data = Japonicus,
                           start = list(cont = 0.001, cont1 = 0.0, cont2 = 0.0))

summary(Fitting_deltaA_quad)

# Compute the AIC value
AIC(Fitting_deltaA,Fitting_deltaA_quad)

# mod <- function(te){
#   c <- as.numeric(Fitting_deltaA$m$getPars()[1])
#   c1 <- as.numeric(Fitting_deltaA$m$getPars()[2])
#   c2 <- as.numeric(Fitting_deltaA$m$getPars()[3])
#   c*te^2+c1*te+c2
# }

## Exponential Fit
# Fitting_deltaA <- nls(lifespan ~ cont*exp(cont1+cont2*Temp),
#                       data = Japonicus,
#                       start = list(cont = 15, cont1 = 2.1, cont2 = -0.01))
# 
# summary(Fitting_deltaA)

# mod <- function(te){
#   c <- as.numeric(Fitting_deltaA$m$getPars()[1])
#   c1 <- as.numeric(Fitting_deltaA$m$getPars()[2])
#   c2 <- as.numeric(Fitting_deltaA$m$getPars()[3])
#   c*te^2+c1*te+c2
# }
vec <- seq(0,45,0.001)
df_out_deltaA <- data.frame(temp_ae = vec,
                            deltaA_jap <- sapply(vec, mod))
colnames(df_out_deltaA) <- c("temp_ae","deltaA_jap")
df_out_deltaA[which(df_out_deltaA$deltaA_jap < 0 ),2] <- 0
# df_out_deltaA[which(df_out_deltaA$temp_ae < 5 ),2] <- 0
# df_out_deltaA[which(df_out_deltaA$temp_ae > 34 ),2] <- 0
plotdeltaA <- ggplot(df_out_deltaA) +
  geom_line(aes(temp_ae,deltaA_jap), size = 0.7) +
  geom_point(data = Japonicus,
             aes(x = Temp,y = lifespan),
             size = 0.9, color = "red") +
  xlim(c(0,45))  +
  ylab("Adult Life Span") + xlab("Temperature (Cº)") +
  theme_bw()
plotdeltaA

####----------- +/- SD-------------####
# Mean - sd
mod_min <- function(te){
  c <- as.numeric(Fitting_deltaA$m$getPars()[1]) - 
    summary(Fitting_deltaA)$coefficients[1,2]
  c1 <- as.numeric(Fitting_deltaA$m$getPars()[2])-
    summary(Fitting_deltaA)$coefficients[2,2]
  c*te+c1
}

vec <- seq(0,45,0.001)
df_out_deltaA_min <- data.frame(temp_ae = vec,
                                deltaA_jap <- sapply(vec, mod_min))
colnames(df_out_deltaA_min) <- c("temp_ae","deltaA_jap")
# df_out_deltaA_min[which(df_out_deltaA_min$deltaA_jap < 0),2] <- 0
df_out_deltaA_min$group <- "min"

ggplot(df_out_deltaA_min) + 
  geom_line(aes(temp_ae,deltaA_jap))

### Mean + sd
mod_max <- function(te){
  c <- as.numeric(Fitting_deltaA$m$getPars()[1]) +
    summary(Fitting_deltaA)$coefficients[1,2]
  c1 <- as.numeric(Fitting_deltaA$m$getPars()[2]) + 
    summary(Fitting_deltaA)$coefficients[2,2]
  c*te+c1
}

vec <- seq(0,45,0.001)
df_out_deltaA_max <- data.frame(temp_ae = vec,
                                deltaA_jap <- sapply(vec, mod_max))
colnames(df_out_deltaA_max) <- c("temp_ae","deltaA_jap")
# df_out_deltaA_max[which(df_out_deltaA_max$deltaA_jap < 0),2] <- 0
df_out_deltaA_max$group <- "max"
df_out_deltaA$group <- "mean"

# Plot all three curves together
df_out_deltaA <- rbind(df_out_deltaA_max,
                       df_out_deltaA_min,
                       df_out_deltaA)

plotdeltaA <- ggplot(df_out_deltaA) +
  geom_line(aes(temp_ae,deltaA_jap,
                color = group,
                group = group, 
                alpha = group), size = 0.7) +
  geom_point(data = Japonicus,aes(x = Temp,y = lifespan),
             size = 0.9, color = "black") +
  scale_color_manual(values=c("red", "blue", "red")) + 
  scale_alpha_manual(values = c(0.5,1,0.5)) +
  xlim(c(5,35)) + ylim(c(0,80)) +
  guides( color =FALSE, alpha = FALSE) +
  ylab("Adult life span") + xlab("Temperature (Cº)") +
  theme_bw() 

plotdeltaA 

# add grey ribbon --------------------------------------------
df_out_deltaA_w <- reshape(df_out_deltaA, idvar ="temp_ae" ,
        timevar = "group", direction = "wide")

plotdeltaA_w <- ggplot(df_out_deltaA_w, aes(x=temp_ae,y=deltaA_jap.min)) +
  geom_line(aes(y=deltaA_jap.min), color = "grey", size = 0.7, alpha=0.5) +
  geom_line(aes(y=deltaA_jap.max), color = "grey", size = 0.7, alpha=0.5) +
  geom_line(aes(y=deltaA_jap.mean), color = "blue", size = 0.7) +
  geom_ribbon(data = df_out_deltaA_w,aes(ymin=deltaA_jap.min,
                  ymax=deltaA_jap.max), fill="grey", alpha=0.5) +
  geom_point(data = Japonicus,aes(x = Temp,y = lifespan),
             size = 0.9, color = "black") +
  xlim(c(5,33)) +  ylim(c(-20,80)) + 
  guides( color =FALSE, alpha = FALSE) +
  ylab("Adult life span") + xlab("Temperature (Cº)") +
  theme_bw() 
plotdeltaA_w  

###----------------------------------------------
# Path <- "~/INVASIBILITY_THRESHOLD/data/japonicus/japonicus_temp_developmenttime.csv"
# developL <- read.csv(Path)
# head(developL)
# 
# n_points <- 8
# r1 <- rnorm(n_points,as.numeric(gsub(",", ".",developL$First_instar_mean[1])),
#                    as.numeric(gsub(",", ".",developL$First_instar_sd[1])) )
# r2 <- rnorm(n_points,as.numeric(gsub(",", ".",developL$First_instar_mean[2])),
#             as.numeric(gsub(",", ".",developL$First_instar_sd[2])) )
# r3 <- rnorm(n_points,as.numeric(gsub(",", ".",developL$First_instar_mean[3])),
#             as.numeric(gsub(",", ".",developL$First_instar_sd[3])) )
# r4 <- rnorm(n_points,as.numeric(gsub(",", ".",developL$First_instar_mean[4])),
#             as.numeric(gsub(",", ".",developL$First_instar_sd[4])) )
# r5 <- rnorm(n_points,as.numeric(gsub(",", ".",developL$First_instar_mean[5])),
#             as.numeric(gsub(",", ".",developL$First_instar_sd[5])) )
# 
# developL <- data.frame(Temp = sort(rep(developL[,1],n_points)),
#                        First_instar_mean = c(r1,r2,r3,r4,r5))
#
# saveRDS(developL, "~/INVASIBILITY_THRESHOLD/data/japonicus/developL.Rds")
# developL$First_instar_mean <- 1/developL$First_instar_mean

###### Random sample from a gaussian distribution, since if we take the mean
# there are only 4 points in the data, and the fit it is overfitting
# We do it once since each iteration the curves are different. So we save the random sample.
# tiene sentido como lo estoy congiendo por que si ves el texto asociado pone que
# first instar es desde huevo a first instar, es decir cuanto tarda en hatch.
developL <- readRDS( "~/INVASIBILITY_THRESHOLD/data/japonicus/developL.Rds")

plot_dE <- ggplot(developL) + 
  geom_point(aes(Temp,First_instar_mean)) + theme_bw()
plot_dE

# Compute the fit for the Briere function:
Fitting_dE <- nls(First_instar_mean ~ cont*Temp*(Temp-cont1)*(cont2-Temp)^(1/2) ,
                  data = developL,
                  start = list(cont = 0.00035, cont1 = 9.5, cont2 = 36))

summary(Fitting_dE)

# Compute the fit for the Quadratic function:
Fitting_dE_quad <- nls(First_instar_mean ~ cont*(Temp-cont1)*(Temp - cont2) ,
                       data = developL,
                       start = list(cont = 0.00035, cont1 = 9.5, cont2 = 36))

summary(Fitting_dE_quad)

# Compute the fit for the Linear function:
Fitting_dE_lin <- nls(First_instar_mean ~ cont*Temp + cont1 ,
                      data = developL,
                      start = list(cont = 0.00035, cont1 = 0))

summary(Fitting_dE_lin)

# Compute the AIC value for the three models:
AIC(Fitting_dE_lin,Fitting_dE,Fitting_dE_quad)

mod <- function(te){
  c <- as.numeric(Fitting_dE$m$getPars()[1])
  c1 <- as.numeric(Fitting_dE$m$getPars()[2])
  c2 <- as.numeric(Fitting_dE$m$getPars()[3])
  c*te*(te-c1)*(c2-te)^(1/2)
}


vec <- seq(0,45,0.001)
df_out_dE  <- data.frame(temp_ae = vec,
                         dE_jap <- sapply(vec, mod))
colnames(df_out_dE) <- c("temp_ae","dE_jap")
df_out_dE[which(df_out_dE$dE_jap < 0),2] <- 0
df_out_dE$group <- "mean"
plotdE <- ggplot(df_out_dE) +
  geom_line(aes(temp_ae,dE_jap), size = 0.7) +
  geom_point(data = developL,aes(Temp,First_instar_mean), size = 0.9, color = "red") +
  xlim(c(0,40)) + ylim(c(0,0.7)) +  
  ylab("Develop rate from Egg to Larva") + xlab("Temperature (Cº)") +
  theme_bw()
plotdE

####----------- +/- SD-------------####
# Mean - sd
mod_min <- function(te){
  c <- as.numeric(Fitting_dE$m$getPars()[1]) - summary(Fitting_dE)$coefficients[1,2]
  c1 <- as.numeric(Fitting_dE$m$getPars()[2])- summary(Fitting_dE)$coefficients[2,2]
  c2 <- as.numeric(Fitting_dE$m$getPars()[3])- summary(Fitting_dE)$coefficients[3,2]
  c*te*(te-c1)*(c2-te)^(1/2)
}

vec <- seq(0,45,0.001)
df_out_dE_min  <- data.frame(temp_ae = vec,
                             dE_jap <- sapply(vec, mod_min))
colnames(df_out_dE_min) <- c("temp_ae","dE_jap")
df_out_dE_min[which(df_out_dE_min$dE_jap < 0),2] <- 0
df_out_dE_min$group <- "min"

### Mean + sd
mod_max <- function(te){
  c <- as.numeric(Fitting_dE$m$getPars()[1]) + summary(Fitting_dE)$coefficients[1,2]
  c1 <- as.numeric(Fitting_dE$m$getPars()[2]) + summary(Fitting_dE)$coefficients[2,2]
  c2 <- as.numeric(Fitting_dE$m$getPars()[3]) + summary(Fitting_dE)$coefficients[3,2]
  c*te*(te-c1)*(c2-te)^(1/2)
}

vec <- seq(0,45,0.001)
df_out_dE_max  <- data.frame(temp_ae = vec,
                             dE_jap <- sapply(vec, mod_max))
colnames(df_out_dE_max) <- c("temp_ae","dE_jap")
df_out_dE_max[which(df_out_dE_max$dE_jap < 0),2] <- 0
df_out_dE_max$group <- "max"

# Plot all three curves together
df_out_dE <- rbind(df_out_dE_max,df_out_dE_min, df_out_dE)
plotdE <- ggplot(df_out_dE) +
  geom_line(aes(temp_ae,dE_jap,
                color = group,
                group = group, 
                alpha = group), size = 0.7) +
  geom_point(data = developL,aes(Temp,First_instar_mean),
             size = 0.9, color = "black") +
  scale_color_manual(values=c("red", "blue", "red")) + 
  scale_alpha_manual(values = c(0.5,1,0.5)) +
  xlim(c(5,36)) + ylim(c(0,0.7)) + 
  guides( color =FALSE, alpha = FALSE) +
  ylab("Develop rate from Egg to Larva") + xlab("Temperature (Cº)") +
  theme_bw()
plotdE

# ####-----------CI-------------####
# new.data <- data.frame(Temp=seq(5, 35.9, by = 0.1))
# interval <- as_tibble(predFit(Fitting_dE, newdata = new.data,
#                               interval = "confidence", level= 0.9)) %>%
#   mutate(Temp = new.data$Temp)
# 
# p1 <-  ggplot(data = developL,
#               aes(x = Temp,y = First_instar_mean)) +
#   geom_point(size = 0.7)
# 
# plotdeltaA <- p1 +
#   geom_line(data = df_out_dE, aes(temp_ae,dE_jap), size = 0.7) +
#   geom_ribbon(data=interval, aes(x=Temp, ymin=lwr, ymax=upr),
#               alpha=0.5, inherit.aes=F, fill="blue") +
#   xlim(c(0,40))  + ylim(c(0,0.5)) +
#   ylab("Adult mortality rate") + xlab("Temperature (Cº)") +
#   theme_bw()
# 
# plotdeltaA

# add grey ribbon --------------------------------------------
df_out_dE_w <- reshape(df_out_dE, idvar ="temp_ae" ,
                           timevar = "group", direction = "wide")

plotdEjap_w <- ggplot(df_out_dE_w, aes(x=temp_ae,y=dE_jap.min)) +
  geom_line(aes(y=dE_jap.min), color = "grey", size = 0.7, alpha=0.5) +
  geom_line(aes(y=dE_jap.max), color = "grey", size = 0.7, alpha=0.5) +
  geom_line(aes(y=dE_jap.mean), color = "blue", size = 0.7) +
  geom_ribbon(data = df_out_dE_w,aes(ymin=dE_jap.min,
                                         ymax=dE_jap.max), fill="grey", alpha=0.5) +
  geom_point(data = developL,aes(Temp,First_instar_mean),
             size = 0.9, color = "black") +
  xlim(c(5,36)) + ylim(c(0,0.7)) + 
  guides( color =FALSE, alpha = FALSE) +
  ylab("Develop rate from Egg to Larva") + xlab("Temperature (Cº)") +
  theme_bw()
plotdEjap_w  

#--------------------------------------------------------
# Paper Germany:
Path <- "~/INVASIBILITY_THRESHOLD/data/japonicus/adult_larva_lifespan.csv"
Japonicus <- read.csv(Path)
head(Japonicus)
## Aunque ponga male es female
# Coger esta variable tiene sentido por que en el experimento cuentan como dia cero
# cuando la larva tiene como maximo 24h. Es decir este tiempo es de larva a adulto
Japonicus$FemaledL <- 1/(Japonicus$Age_emergence_female_mean) 
Japonicus <- Japonicus[,c("Temp","FemaledL")]
# Thesis Jamesina
# Japonicus <- data.frame(Temp <- numeric(), FemaledL <- numeric())
Japonicus <- rbind(Japonicus, c(10,1/140.8))
Japonicus <- rbind(Japonicus, c(16,1/84))
Japonicus <- rbind(Japonicus, c(22,1/31.3))
Japonicus <- rbind(Japonicus, c(28,1/17))
Japonicus <- rbind(Japonicus, c(34,0))
colnames(Japonicus) <-  c("Temp", "FemaledL")
plot_dL <- ggplot(Japonicus) + 
  geom_point(aes(Temp,FemaledL)) + theme_bw()
plot_dL

# Briere function fit    
Fitting_dL <- nls(FemaledL ~ cont*Temp*(Temp-cont1)*(cont2-Temp)^(1/2) ,
                  data = Japonicus, algorithm = "port",
                  start = list(cont = 0.0035, cont1 = 9.5, cont2 = 36), 
                  lower=c(7e-05,min(developL$Temp)-5,max(developL$Temp)+0.1), upper=c(1,min(developL$Temp)-0.1,max(developL$Temp)+4))

summary(Fitting_dL)

# Linear function fit
Fitting_dL_lin <- nls(FemaledL ~ cont*Temp + cont1,
                      data = Japonicus,
                      start = list(cont = 0.00035, cont1 = 0))
summary(Fitting_dL_lin)

# Compute the AIC value
AIC(Fitting_dL,Fitting_dL_lin)
mod <- function(te){
  c <- as.numeric(Fitting_dL$m$getPars()[1])
  c1 <- as.numeric(Fitting_dL$m$getPars()[2])
  c2 <- as.numeric(Fitting_dL$m$getPars()[3])
  c*te*(te-c1)*(c2-te)^(1/2)
}

vec <- seq(0,45,0.001)
df_out_dL  <- data.frame(temp_ae = vec,
                         dL_jap <- sapply(vec, mod))

colnames(df_out_dL) <- c("temp_ae","dL_jap")
df_out_dL[which(df_out_dL$dL_jap < 0),2] <- 0

plotdL <- ggplot(df_out_dL) +
  geom_line(aes(temp_ae,dL_jap), size = 0.7) +
  geom_point(data = Japonicus,aes(Temp,FemaledL), size = 0.9, color = "red") +
  xlim(c(0,45)) + 
  ylab("Develop rate from Larva to Adult") + xlab("Temperature (Cº)") +
  theme_bw()
plotdL

####----------- +/- SD-------------####
# Mean - sd
mod_min <- function(te){
  c <- as.numeric(Fitting_dL$m$getPars()[1]) - 
    summary(Fitting_dL)$coefficients[1,2]
  c1 <- as.numeric(Fitting_dL$m$getPars()[2])- 
    summary(Fitting_dL)$coefficients[2,2]
  c2 <- as.numeric(Fitting_dL$m$getPars()[3])- 
    summary(Fitting_dL)$coefficients[3,2]
  c*te*(te-c1)*(c2-te)^(1/2)
}

vec <- seq(0,45,0.001)
df_out_dL_min  <- data.frame(temp_ae = vec,
                             dL_jap <- sapply(vec, mod_min))

colnames(df_out_dL_min) <- c("temp_ae","dL_jap")
df_out_dL_min[which(df_out_dL_min$dL_jap < 0),2] <- 0
df_out_dL_min$group <- "min"

### Mean + sd
mod_max <- function(te){
  c <- as.numeric(Fitting_dL$m$getPars()[1]) + 
    summary(Fitting_dL)$coefficients[1,2]
  c1 <- as.numeric(Fitting_dL$m$getPars()[2]) +
    summary(Fitting_dL)$coefficients[2,2]
  c2 <- as.numeric(Fitting_dL$m$getPars()[3]) + 
    summary(Fitting_dL)$coefficients[3,2]
  c*te*(te-c1)*(c2-te)^(1/2)
}

vec <- seq(0,45,0.001)
df_out_dL_max  <- data.frame(temp_ae = vec,
                             dL_jap <- sapply(vec, 
                                              mod_max))

colnames(df_out_dL_max) <- c("temp_ae","dL_jap")
df_out_dL_max[which(df_out_dL_max$dL_jap < 0),2] <- 0
df_out_dL_max$group <- "max"
df_out_dL$group <- "mean"

# Plot all three curves together
df_out_dL <- rbind(df_out_dL_min,
                   df_out_dL_max,
                   df_out_dL)
plotdL <- ggplot(df_out_dL) +
  geom_line(aes(temp_ae,dL_jap,
                color = group,
                group = group, 
                alpha = group), size = 0.7) +
  geom_point(data = Japonicus,aes(Temp,FemaledL),
             size = 0.9, color = "black") +
  scale_color_manual(values=c("red", "blue", "red")) + 
  scale_alpha_manual(values = c(0.5,1,0.5)) +
  xlim(c(0,45)) + 
  ylab("Develop rate from Larva to Adult") +
  guides(color = FALSE, alpha = FALSE) +
  xlab("Temperature (Cº)") + ylab("Develop rate from Egg to Larva") + xlab("Temperature (Cº)") +
  theme_bw()
plotdL

# add grey ribbon --------------------------------------------
df_out_dL_w <- reshape(df_out_dL, idvar ="temp_ae" ,
                       timevar = "group", direction = "wide")

plotdL_w <- ggplot(df_out_dL_w, aes(x=temp_ae,y=dL_jap.min)) +
  geom_line(aes(y=dL_jap.min), color = "grey", size = 0.7, alpha=0.5) +
  geom_line(aes(y=dL_jap.max), color = "grey", size = 0.7, alpha=0.5) +
  geom_line(aes(y=dL_jap.mean), color = "blue", size = 0.7) +
  geom_ribbon(data = df_out_dL_w,aes(ymin=dL_jap.min,
                                     ymax=dL_jap.max), fill="grey", alpha=0.5) +
  geom_point(data = Japonicus,aes(Temp,FemaledL),
             size = 0.9, color = "black") +
  xlim(c(0,45)) + 
  ylab("Develop rate from Larva to Adult") +
  guides(color = FALSE, alpha = FALSE) +
  xlab("Temperature (Cº)") + ylab("Develop rate from Egg to Larva") + xlab("Temperature (Cº)") +
  theme_bw()
plotdL_w 

###----------------------------------------------
## Paper Germany:
#https://parasitesandvectors.biomedcentral.com/articles/10.1186/s13071-018-2659-1
Lmortality <- data.frame(Temp = c(0,5,10,12,14,15,17,19,20,23,25,26,27,28,29,31),
                         mean_mort_perc = c(100,99.5,16,38.5,18,15,19,29.5,11.3,48.5,13.8,6,41.5,12.5,70.5,87.5),
                         sd_mort_perc = c(0,1.1,5.5,14.2,4.8,7.9,9.6,12.4,6.5,27.6,8.4,5.2,31.1,7.7,22.2,6.4))
Lmortality$mean_mort_perc = Lmortality$mean_mort_perc/100
Lmortality[nrow(Lmortality) +1,] <- c(10,0.5,0)
Lmortality[nrow(Lmortality) +1,] <- c(16,1-0.6,0)
Lmortality[nrow(Lmortality) +1,] <- c(22,1-0.27,0)
Lmortality[nrow(Lmortality) +1,] <- c(28,1-0.33,0)
Lmortality[nrow(Lmortality) +1,] <- c(34,1,0)
Lmortality[nrow(Lmortality) +1,] <- c(40,1,0)

head(Lmortality)

plot_deltaL <- ggplot(Lmortality) + 
  geom_point(aes(Temp,mean_mort_perc)) + theme_bw()
plot_deltaL

## Quadratic normal fit
Fitting_deltaL <- nls(mean_mort_perc ~ cont*Temp^2 + cont1*Temp + cont2,
                      data = Lmortality,
                      start = list(cont = 15, cont1 = -20,
                                   cont2 = 30))

summary(Fitting_deltaL)

# Linear function fit
Fitting_deltaL_lin <- nls(mean_mort_perc ~ cont*Temp + cont1,
                          data = Lmortality,
                          start = list(cont = 0, cont1 = 0))

summary(Fitting_deltaL_lin)

# AIC for the two models:
AIC(Fitting_deltaL,Fitting_deltaL_lin)

mod <- function(te){
  c <- as.numeric(Fitting_deltaL$m$getPars()[1])
  c1 <- as.numeric(Fitting_deltaL$m$getPars()[2])
  c2 <- as.numeric(Fitting_deltaL$m$getPars()[3])
  c*te^2+c1*te+c2
}

vec <- seq(0,45,0.001)
df_out_deltaL <- data.frame(temp_ae = vec,
                            deltaL_jap <- sapply(vec, mod))

colnames(df_out_deltaL) <- c("temp_ae","deltaL_jap")
df_out_deltaL[which(df_out_deltaL$deltaL_jap < 0),2] <- 0

plotdeltaL <- ggplot(df_out_deltaL) +
  geom_line(aes(temp_ae,deltaL_jap), size = 0.8) +
  geom_point(data = Lmortality,aes(Temp,mean_mort_perc),
             size = 0.9, color = "red") +
  ylab("Larva mortality rate") + xlab("Temperature (Cº)") +
  theme_bw()
plotdeltaL

####----------- +/- SD-------------####
# Mean - sd
mod_min <- function(te){
  c <- as.numeric(Fitting_deltaL$m$getPars()[1]) - 
    summary(Fitting_deltaL)$coefficients[1,2]
  c1 <- as.numeric(Fitting_deltaL$m$getPars()[2])- 
    summary(Fitting_deltaL)$coefficients[2,2]
  c2 <- as.numeric(Fitting_deltaL$m$getPars()[3])- 
    summary(Fitting_deltaL)$coefficients[3,2]
  c*te^2+c1*te+c2
}

vec <- seq(0,45,0.001)
df_out_deltaL_min <- data.frame(temp_ae = vec,
                                deltaL_jap <- sapply(vec,
                                                     mod_min))

colnames(df_out_deltaL_min) <- c("temp_ae","deltaL_jap")
df_out_deltaL_min$group <- "min"

### Mean + sd
mod_max <- function(te){
  c <- as.numeric(Fitting_deltaL$m$getPars()[1]) + 
    summary(Fitting_deltaL)$coefficients[1,2]
  c1 <- as.numeric(Fitting_deltaL$m$getPars()[2])+ 
    summary(Fitting_deltaL)$coefficients[2,2]
  c2 <- as.numeric(Fitting_deltaL$m$getPars()[3])+ 
    summary(Fitting_deltaL)$coefficients[3,2]
  c*te^2+c1*te+c2
}

vec <- seq(0,45,0.001)
df_out_deltaL_max <- data.frame(temp_ae = vec,
                                deltaL_jap <- sapply(vec,
                                                     mod_max))

colnames(df_out_deltaL_max) <- c("temp_ae","deltaL_jap")
df_out_deltaL_max[which(df_out_deltaL_max$deltaL_jap < 0),2] <- 0
df_out_deltaL_max$group <- "max"
df_out_deltaL$group <- "mean"

# Plot all three curves together
df_out_deltaL <- rbind(df_out_deltaL_min,
                       df_out_deltaL_max,
                       df_out_deltaL)

plotdeltaL <- ggplot(df_out_deltaL) +
  geom_line(aes(temp_ae,deltaL_jap,
                color = group,
                group = group, 
                alpha = group), size = 0.8) +
  geom_point(data = Lmortality,aes(Temp,mean_mort_perc),
             size = 0.9, color = "black") +
  ylab("Larva mortality rate") + xlab("Temperature (Cº)") +
  scale_color_manual(values=c("red", "blue", "red")) + 
  scale_alpha_manual(values = c(0.5,1,0.5)) +
  xlim(c(5,35)) + ylim(c(-1,2)) + 
  guides(color = FALSE, alpha = FALSE) +
  theme_bw()
plotdeltaL

# add grey ribbon --------------------------------------------
df_out_deltaL_w <- reshape(df_out_deltaL, idvar ="temp_ae" ,
                       timevar = "group", direction = "wide")

plotdeltaL_w <- ggplot(df_out_deltaL_w, aes(x=temp_ae,y=deltaL_jap.min)) +
  geom_line(aes(y=deltaL_jap.min), color = "grey", size = 0.7, alpha=0.5) +
  geom_line(aes(y=deltaL_jap.max), color = "grey", size = 0.7, alpha=0.5) +
  geom_line(aes(y=deltaL_jap.mean), color = "blue", size = 0.7) +
  geom_ribbon(data = df_out_deltaL_w,aes(ymin=deltaL_jap.min,
                                     ymax=deltaL_jap.max), fill="grey", alpha=0.5) +
  geom_point(data = Lmortality,aes(Temp,mean_mort_perc),
             size = 0.9, color = "black") +
  ylab("Larva mortality rate") + xlab("Temperature (Cº)") +
  xlim(c(5,35)) + ylim(c(-0.6,2.3)) + 
  guides(color = FALSE, alpha = FALSE) +
  theme_bw()
plotdeltaL_w 

###----------------------------------------------
library(ggpubr)
sizelet = 14
ggarrange(plotdE  +
            theme(text = element_text(size = sizelet)) +
            xlab("") ,
          plotdL +
            theme(text = element_text(size = sizelet)) +
            ylab("Development rate from Egg to Larva") +
            xlab(""),
          plotdeltaL + ylim(c(0,1.3)) +
            theme(text = element_text(size = sizelet)),
          plotdeltaA  +
            theme(text = element_text(size = sizelet)))

# Albopictus and Aegypti --------------------------------
# Data frame data taken from Tun-Lin et al 2001 
# https://resjournals.onlinelibrary.wiley.com/doi/full/10.1046/j.1365-2915.2000.00207.x
# https://academic.oup.com/jme/article/56/6/1661/5505308

df_aeg <- data.frame(temp = c(10,15,20,25,27,30,35,40,12,14,16),
                     proportion_surv = c(0,0.235,0.9,0.88,0.93,0.88,0.67,0,0.23,0.87,0.93))

plot_aeg <- ggplot(df_aeg) + 
  geom_point(aes(temp,proportion_surv)) + theme_bw()

Fitting_aeg <- nls(proportion_surv ~ (-cont*(temp-Tmin)*(temp - Tmax)),
                   data = df_aeg,
                   start = list(cont = 0.001, Tmin = 5, Tmax = 20))

summary(Fitting_aeg)

Fitting_aeg_lin <- nls(proportion_surv ~ cont*temp+cont1,
                       data = df_aeg,
                       start = list(cont = 0.001, cont1 = 0))

summary(Fitting_aeg_lin)

AIC(Fitting_aeg,Fitting_aeg_lin)

mod <- function(te){
  t0 <- as.numeric(Fitting_aeg$m$getPars()[2])
  tm <- as.numeric(Fitting_aeg$m$getPars()[3])
  c <- as.numeric(Fitting_aeg$m$getPars()[1])
  (-c*(te-t0)*(te - tm))
}

vec <- seq(0,45,0.01)
df_out_aeg <- data.frame(temp_ae = vec, life_span_ae <- sapply(vec, mod))
colnames(df_out_aeg) <- c("temp_ae", "life_span_ae")
df_out_aeg[which(df_out_aeg$life_span_ae < 0),2] <- 0
df_out_aeg$group <- "mean"

###---------------+/- SD---------------------######
## Mean - SD
mod_min <- function(te){
  t0 <- as.numeric(Fitting_aeg$m$getPars()[2]) - 
    summary(Fitting_aeg)$coefficients[2,2]
  tm <- as.numeric(Fitting_aeg$m$getPars()[3]) - 
    summary(Fitting_aeg)$coefficients[3,2]
  c <- as.numeric(Fitting_aeg$m$getPars()[1]) - 
    summary(Fitting_aeg)$coefficients[1,2]
  (-c*(te-t0)*(te - tm))
}

vec <- seq(0,45,0.01)
df_out_aeg_min <- data.frame(temp_ae = vec,
                             life_span_ae <- sapply(vec, mod_min))
colnames(df_out_aeg_min) <- c("temp_ae", "life_span_ae")
df_out_aeg_min[which(df_out_aeg_min$life_span_ae < 0),2] <- 0
df_out_aeg_min$group <- "min"

plot(df_out_aeg_min$temp_ae,df_out_aeg_min$life_span_ae)

## Mean - SD
mod_max <- function(te){
  t0 <- as.numeric(Fitting_aeg$m$getPars()[2]) + 
    summary(Fitting_aeg)$coefficients[2,2]
  tm <- as.numeric(Fitting_aeg$m$getPars()[3]) + 
    summary(Fitting_aeg)$coefficients[3,2]
  c <- as.numeric(Fitting_aeg$m$getPars()[1]) + 
    summary(Fitting_aeg)$coefficients[1,2]
  (-c*(te-t0)*(te - tm))
}

vec <- seq(0,45,0.01)
df_out_aeg_max <- data.frame(temp_ae = vec,
                             life_span_ae <- sapply(vec, mod_max))
colnames(df_out_aeg_max) <- c("temp_ae", "life_span_ae")
df_out_aeg_max[which(df_out_aeg_max$life_span_ae < 0),2] <- 0
df_out_aeg_max$group <- "max"
plot(df_out_aeg_max$temp_ae,df_out_aeg_max$life_span_ae)

# Plot all three curves together
df_out_aeg <- rbind(df_out_aeg_max,
                    df_out_aeg_min,
                    df_out_aeg)

plotaeg <- ggplot(df_out_aeg) +
  geom_line(aes(temp_ae,life_span_ae,
                color = group,
                group = group, 
                alpha = group), size = 0.7) +
  geom_point(data = df_aeg,aes(temp,proportion_surv),
             size = 0.9, color = "black") +
  scale_color_manual(values=c("red", "blue", "red")) + 
  scale_alpha_manual(values = c(0.5,1,0.5)) +
  xlim(c(5,42)) + ylim(c(0,1.3)) +
  guides( color =FALSE, alpha = FALSE) +
  ylab("Development rate from Larvae to Adult") + xlab("Temperature (Cº)") +
  theme_bw() 
plotaeg

# add ribbon -----------------------------------------------------
df_out_aeg_w <- reshape(df_out_aeg, idvar ="temp_ae" ,
                           timevar = "group", direction = "wide")

plotaeg_w <- ggplot(df_out_aeg_w, aes(x=temp_ae,y=life_span_ae.min)) +
  geom_line(aes(y=life_span_ae.min), color = "grey", size = 0.7, alpha=0.5) +
  geom_line(aes(y=life_span_ae.max), color = "grey", size = 0.7, alpha=0.5) +
  geom_line(aes(y=life_span_ae.mean), color = "blue", size = 0.7) +
  geom_ribbon(data = df_out_aeg_w,aes(ymin=life_span_ae.min,
                                         ymax=life_span_ae.max), fill="grey",
              alpha=0.5) +
  geom_point(data = df_aeg,aes(temp,proportion_surv),
             size = 0.9, color = "black") +
  xlim(c(5,42)) + ylim(c(0,1.3)) +
  guides( color =FALSE, alpha = FALSE) +
  ylab("Development rate from Larvae to Adult") + xlab("Temperature (Cº)") +
  theme_bw() 
plotaeg_w 

# Data frame data taken from Delatte et al 2009. 
# https://academic.oup.com/jme/article/46/1/33/902827?login=false
df_albo <- data.frame(temp = c(5,10,15,20,25,30,35,40),
                      proportion_surv = c(0,0,4/8,62/80,61/80,54/80,3.5/120,0))

ggplot(df_albo) + 
  geom_point(aes(temp,proportion_surv)) + theme_bw()

### Survival curve Laura:
# Path <- "~/Documentos/PHD/2023/Laura proj/DATA_LARVA_TEM.csv"
# df <- read.csv(file = Path)
# plot_alb <- ggplot(df) + geom_point(aes(temp_chamber,total_lived))
# Path <- "~/INVASIBILITY_THRESHOLD/data/pLA_Laura.Rds"
# albo_lau <- readRDS(Path)[,c(1,5)]
# colnames(albo_lau) <- colnames(df_albo)
# df_albo <- rbind(df_albo, albo_lau)

ggplot(df_albo) + geom_point(aes(temp,proportion_surv))

# My data nls:
Fitting <- nls(proportion_surv ~ (-cont*(temp-Tmin)*(temp - Tmax)),
               data = df_albo,
               start = list(cont = 0.001, Tmin = 5, Tmax = 20))

summary(Fitting)

Fitting_Lin <- nls(proportion_surv ~ cont*temp+cont1,
               data = df_albo,
               start = list(cont = 0.001, cont1 = 0))

summary(Fitting_Lin)
AIC(Fitting_Lin,Fitting)

mod <- function(te){
  t0 <- as.numeric(Fitting$m$getPars()[2])
  tm <- as.numeric(Fitting$m$getPars()[3])
  c <- as.numeric(Fitting$m$getPars()[1])
  (-c*(te-t0)*(te - tm))
}

vec <- seq(0,40,0.001)
df_alb <- data.frame(temp = vec, life_span <- sapply(vec, mod))
colnames(df_alb) <- c("temp", "life_span")
df_alb[which(df_alb$life_span < 0),2] <- 0

###------------ Mean +/- SD -------------###
### Mean+SD
mod_max <- function(te){
  t0 <- as.numeric(Fitting$m$getPars()[2]) +
    summary(Fitting)$coefficients[2,2]
  tm <- as.numeric(Fitting$m$getPars()[3]) +
    summary(Fitting)$coefficients[3,2]
  c <- as.numeric(Fitting$m$getPars()[1]) +
    summary(Fitting)$coefficients[1,2]
  (-c*(te-t0)*(te - tm))
}

vec <- seq(0,40,0.001)
df_alb_max <- data.frame(temp = vec, life_span <- sapply(vec, mod_max))
colnames(df_alb_max) <- c("temp", "life_span")
df_alb_max[which(df_alb_max$life_span < 0),2] <- 0
df_alb_max$group <- "max"

### Mean-SD
mod_min <- function(te){
  t0 <- as.numeric(Fitting$m$getPars()[2]) -
    summary(Fitting)$coefficients[2,2]
  tm <- as.numeric(Fitting$m$getPars()[3]) -
    summary(Fitting)$coefficients[3,2]
  c <- as.numeric(Fitting$m$getPars()[1]) -
    summary(Fitting)$coefficients[1,2]
  (-c*(te-t0)*(te - tm))
}

vec <- seq(0,40,0.001)
df_alb_min <- data.frame(temp = vec, life_span <- sapply(vec, mod_min))
colnames(df_alb_min) <- c("temp", "life_span")
df_alb_min[which(df_alb_min$life_span < 0),2] <- 0
df_alb_min$group <- "min"
df_alb$group <- "mean"

# Plot all three curves together
df_alb <- rbind(df_alb_max,
                df_alb_min,
                df_alb)

plotalb <- ggplot(df_alb) +
  geom_line(aes(temp,life_span,
                color = group,
                group = group, 
                alpha = group), size = 0.7) +
  geom_point(data = df_albo,aes(temp,proportion_surv),
             size = 0.9, color = "black") +
  scale_color_manual(values=c("red", "blue", "red")) + 
  scale_alpha_manual(values = c(0.5,1,0.5)) +
  xlim(c(5,40)) +
  guides( color =FALSE, alpha = FALSE) +
  ylab("Prob from Larva to Adult") + xlab("Temperature (Cº)") +
  theme_bw() 
plotalb

# add ribbon -----------------------------------------------------
df_alb_w <- reshape(df_alb, idvar ="temp" ,
                        timevar = "group", direction = "wide")

plotalb_w <- ggplot(df_alb_w, aes(x=temp,y=life_span.min)) +
  geom_line(aes(y=life_span.min), color = "grey", size = 0.7, alpha=0.5) +
  geom_line(aes(y=life_span.max), color = "grey", size = 0.7, alpha=0.5) +
  geom_line(aes(y=life_span.mean), color = "blue", size = 0.7) +
  geom_ribbon(data = df_alb_w,aes(ymin=life_span.min,
                                      ymax=life_span.max), fill="grey",
              alpha=0.5) +
  geom_point(data = df_albo,aes(temp,proportion_surv),
             size = 0.9, color = "black") +
  xlim(c(5,40)) +
  guides( color =FALSE, alpha = FALSE) +
  ylab("Prob from Larva to Adult") + xlab("Temperature (Cº)") +
  theme_bw() 
plotalb_w 

# Join all the plots in one figure
library(ggpubr)
ggarrange(plotalb  + ggtitle("Aedes Albopictus") +
            theme(text = element_text(size = 15)) +
            ylab("Probability from Larvae to Adult"),
          plotaeg + ggtitle("Aedes Aegypti") +
            theme(text = element_text(size = 15)) + 
            ylab(""))

#--------------------Egg Development Aedes Albopictus--------------------#
# ## Info taken Table1: https://www.scielo.br/j/rsp/a/dvPQ8QMr7Y687hPJxsTxjDg/abstract/?lang=en
df_dE <- data.frame(temp = c(5,15,20,25,30,35),
                     develop_rate_mu = c(11,7.4,2.9,4.5,6.7,7.1),
                     develop_rate_sd = c(1.3,1.8,0.4,0.7,0.7,0.8)
                     )

n = 8
r1 <- rnorm(n, df_dE$develop_rate_mu[1],
            df_dE$develop_rate_sd[1] )
r2 <- rnorm(n, df_dE$develop_rate_mu[2],
            df_dE$develop_rate_sd[2] )
r3 <- rnorm(n, df_dE$develop_rate_mu[3],
            df_dE$develop_rate_sd[3] )
r4 <- rnorm(n, df_dE$develop_rate_mu[4],
            df_dE$develop_rate_sd[4] )
r5 <- rnorm(n, df_dE$develop_rate_mu[5],
            df_dE$develop_rate_sd[5] )
r6 <- rnorm(n, df_dE$develop_rate_mu[6],
            df_dE$develop_rate_sd[6] )

df_dE <- data.frame(temp = sort(rep(df_dE[,1],n)),
                        develop_rate = c(r1,r2,r3,r4,r5,r6))

df_dE$develop_rate <- 1/df_dE$develop_rate
saveRDS(df_dE, "~/INVASIBILITY_THRESHOLD/data/df_dE.Rds")

###### Random sample from a Gaussian distribution, since if we take the mean
# thre are only 4 points in the data, and the fit it is overfitting
# We do it once since each iteration the curves are different. So we save the random sample.
df_dE <- readRDS( "~/INVASIBILITY_THRESHOLD/data/df_dE.Rds")
df_dE$develop_rate
plot_dE <- ggplot(df_dE) + 
  geom_point(aes(temp,develop_rate)) + theme_bw()
plot_dE

# Fitting_dE <- nls(develop_rate ~ c*temp*(temp-c1)*(c2-temp)^(1/2),
#                   data = df_dE,
#                   start = list(c = 0.001, c1 = 5 , c2 = 30))


Fitting_dE <- nls(develop_rate ~ c*(temp-c1)*(temp - c2),
                  data = df_dE,
                  start = list(c = 0.001, c1 = 5 , c2 = 30))

summary(Fitting_dE)

Fitting_dE_lin <- nls(develop_rate ~ c*temp+c1,
                      data = df_dE,
                      start = list(c = 0.001, c1 = 0))

summary(Fitting_dE_lin)

AIC(Fitting_dE,Fitting_dE_lin)

mod <- function(te){
  c <- as.numeric(Fitting_dE$m$getPars()[1])
  c1 <- as.numeric(Fitting_dE$m$getPars()[2])
  c2 <- as.numeric(Fitting_dE$m$getPars()[3])
  # c*te*(te-c1)*(c2-te)^(1/2)
  c*(te-c1)*(te-c2)
}

vec <- seq(5,40,0.01)
df_out <- data.frame(temp = vec,
                     devep_rate <- sapply(vec, mod))
df_out[which(df_out$devep_rate < 0),2] <- 0
df_out$group <- "mean"

##--------- Mean - SD ----------###
mod_min <- function(te){
  c <- as.numeric(Fitting_dE$m$getPars()[1]) - 
    summary(Fitting_dE)$coefficients[1,2]
  c1 <- as.numeric(Fitting_dE$m$getPars()[2]) - 
    summary(Fitting_dE)$coefficients[2,2]
  c2 <- as.numeric(Fitting_dE$m$getPars()[3]) - 
    summary(Fitting_dE)$coefficients[3,2]
  # c*te*(te-c1)*(c2-te)^(1/2)
  c*(te-c1)*(te-c2)
}

df_out_min <- data.frame(temp = vec,
                         devep_rate = sapply(vec, mod_min))
df_out_min[which(df_out_min$devep_rate < 0),2] <- 0
df_out_min$group <- "min"

##--------- Mean + SD ----------###
mod_max <- function(te){
  c <- as.numeric(Fitting_dE$m$getPars()[1]) + 
    summary(Fitting_dE)$coefficients[1,2]
  c1 <- as.numeric(Fitting_dE$m$getPars()[2]) + 
    summary(Fitting_dE)$coefficients[2,2]
  c2 <- as.numeric(Fitting_dE$m$getPars()[3]) + 
    summary(Fitting_dE)$coefficients[3,2]
  # c*te*(te-c1)*(c2-te)^(1/2)
  c*(te-c1)*(te-c2)
}

df_out_max <- data.frame(temp = vec, 
                         devep_rate = sapply(vec, mod_max))
df_out_max[which(df_out_max$devep_rate < 0),2] <- 0
df_out_max$group <- "max"
colnames(df_out) <- colnames(df_out_min)
# Plot all three curves together
df_out <- rbind(df_out_max,
                df_out_min,
                df_out)

plotdE <- ggplot(df_out) +
  geom_line(aes(temp,devep_rate,
                color = group,
                group = group, 
                alpha = group), size = 0.7) +
  geom_point(data = df_dE, aes(temp,develop_rate),
             size = 0.9, color = "black") +
  scale_color_manual(values=c("red", "blue", "red")) + 
  scale_alpha_manual(values = c(0.5,1,0.5)) +
  xlim(c(5,40)) +
  guides( color =FALSE, alpha = FALSE) +
  ylab("Egg development time") + xlab("Temperature (Cº)") +
  theme_bw() 

plotdE 

# add ribbon -----------------------------------------------------
df_out_w <- reshape(df_out, idvar ="temp" ,
                    timevar = "group", direction = "wide")
df_out_w$devep_rate.min <- ifelse(is.na(df_out_w$devep_rate.min),0, df_out_w$devep_rate.min)
plotdE_w <- ggplot(df_out_w, aes(x=temp,y=devep_rate.max)) +
  geom_line(aes(y=devep_rate.min), color = "grey", size = 0.7, alpha=0.5) +
  geom_line(aes(y=devep_rate.max), color = "grey", size = 0.7, alpha=0.5) +
  geom_line(aes(y=devep_rate.mean), color = "blue", size = 0.7) +
  geom_ribbon(data = df_out_w,aes(ymin=devep_rate.max,
                                  ymax=devep_rate.min), fill="grey",
              alpha=0.5) +
  geom_point(data = df_dE, aes(temp,develop_rate),
             size = 0.9, color = "black") +
  xlim(c(5,40)) +  
  guides( color =FALSE, alpha = FALSE) +
  ylab("Egg development time") + xlab("Temperature (Cº)") +
  theme_bw() 
plotdE_w 

### --------------Development Egg Aegypti -----------------------#
## https://pubmed.ncbi.nlm.nih.gov/19274388/
df_dE_aeg <- data.frame(temp=c(16,22,25,28,31,35),
                        hours_mean = c(489.3,98.3,77.4,61.6,48.4,50.3),
                        days_sd = c(0.6,0.7,0.8,1.2,0.5,0.3))

n = 8
r1 <- rnorm(n, df_dE_aeg$hours_mean[1],
            df_dE_aeg$days_sd[1] )
r2 <- rnorm(n, df_dE_aeg$hours_mean[2],
            df_dE_aeg$days_sd[2] )
r3 <- rnorm(n, df_dE_aeg$hours_mean[3],
            df_dE_aeg$days_sd[3] )
r4 <- rnorm(n, df_dE_aeg$hours_mean[4],
            df_dE_aeg$days_sd[4] )
r5 <- rnorm(n, df_dE_aeg$hours_mean[5],
            df_dE_aeg$days_sd[5] )
r6 <- rnorm(n, df_dE_aeg$hours_mean[6],
            df_dE_aeg$days_sd[6] )

df_dE_aeg <- data.frame(temp = sort(rep(df_dE_aeg[,1],n)),
                        hour_devep = c(r1,r2,r3,r4,r5,r6))

df_dE_aeg$develop_rate <- 1/(df_dE_aeg$hour_devep/24)
df_dE_aeg <- rbind(df_dE_aeg, c(36,0,0,0.001))

# saveRDS(df_dE_aeg, "~/INVASIBILITY_THRESHOLD/data/df_dE_aeg.Rds")

## I use the random sample fixed, other wise the param would change.
# df_dE_aeg <- readRDS("~/INVASIBILITY_THRESHOLD/data/df_dE_aeg.Rds")


ggplot(df_dE_aeg) + 
  geom_point(aes(temp,develop_rate))

Fitting_dE_aeg <- nls(develop_rate ~ c*temp*(temp-c1)*(c2-temp)^(1/2),
                      data = df_dE_aeg, algorithm = "port",
                      start = list(c = 0.0003, c1 = 14, c2 = 40), 
                      lower=c(0.00001,7,30), upper = c(0.04,15,50))

summary(Fitting_dE_aeg)

Fitting_dE_aeg_lin <- nls(develop_rate ~ c*temp + c1,
                          data = df_dE_aeg,
                          start = list(c = 0.0003, c1 = 0))

summary(Fitting_dE_aeg_lin)

AIC(Fitting_dE_aeg, Fitting_dE_aeg_lin)

mod <- function(te){
  c <- as.numeric(Fitting_dE_aeg$m$getPars()[1])
  c1 <- as.numeric(Fitting_dE_aeg$m$getPars()[2])
  c2 <- as.numeric(Fitting_dE_aeg$m$getPars()[3])
  c*te*(te-c1)*(c2-te)^(1/2)
}

vec <- seq(0,45,0.001)
df_out_aeg <- data.frame(temp = vec,
                         devep_rate = sapply(vec, mod))
df_out_aeg[which(df_out_aeg$devep_rate < 0),2] <- 0

####------Mean +/- SD -----------###
# Mean + SD
mod_max <- function(te){
  c <- as.numeric(Fitting_dE_aeg$m$getPars()[1])  + 
    summary(Fitting_dE_aeg)$coefficients[1,2]
  c1 <- as.numeric(Fitting_dE_aeg$m$getPars()[2])  + 
    summary(Fitting_dE_aeg)$coefficients[2,2]
  c2 <- as.numeric(Fitting_dE_aeg$m$getPars()[3])  + 
    summary(Fitting_dE_aeg)$coefficients[3,2]
  c*te*(te-c1)*(c2-te)^(1/2)
}

df_out_aeg_max <- data.frame(temp = vec,
                             devep_rate = sapply(vec, mod_max))
df_out_aeg_max[which(df_out_aeg_max$devep_rate < 0),2] <- 0
df_out_aeg_max$group = "max"

## MEan - SD
mod_min <- function(te){
  c <- as.numeric(Fitting_dE_aeg$m$getPars()[1])  - 
    summary(Fitting_dE_aeg)$coefficients[1,2]
  c1 <- as.numeric(Fitting_dE_aeg$m$getPars()[2])  - 
    summary(Fitting_dE_aeg)$coefficients[2,2]
  c2 <- as.numeric(Fitting_dE_aeg$m$getPars()[3])  - 
    summary(Fitting_dE_aeg)$coefficients[3,2]
  c*te*(te-c1)*(c2-te)^(1/2)
}

df_out_aeg_min <- data.frame(temp = vec,
                             devep_rate = sapply(vec, mod_min))
df_out_aeg_min[which(df_out_aeg_min$devep_rate < 0),2] <- 0
df_out_aeg_min$group = "min"
df_out_aeg$group = "mean"


df_out_aeg <- rbind(df_out_aeg_max,
                    df_out_aeg_min,
                    df_out_aeg)

plotdE_aeg <- ggplot(df_out_aeg) +
  geom_line(aes(temp,devep_rate,
                color = group,
                group = group, 
                alpha = group), size = 0.7) +
  geom_point(data =  df_dE_aeg, aes(temp,develop_rate),
             size = 0.9, color = "black") +
  scale_color_manual(values=c("red", "blue", "red")) + 
  scale_alpha_manual(values = c(0.5,1,0.5)) +
  xlim(c(5,40)) +
  guides( color =FALSE, alpha = FALSE) +
  ylab("Egg development time") + xlab("Temperature (Cº)") +
  theme_bw()
plotdE_aeg

# add ribbon -----------------------------------------------------
df_out_aeg_w <- reshape(df_out_aeg, idvar ="temp" ,
                    timevar = "group", direction = "wide")
df_out_aeg_w$devep_rate.min <- ifelse(is.na(df_out_aeg_w$devep_rate.min),0, df_out_aeg_w$devep_rate.min)
plotdE_aeg_w <- ggplot(df_out_aeg_w, aes(x=temp,y=devep_rate.max)) +
  geom_line(aes(y=devep_rate.min), color = "grey", size = 0.7, alpha=0.5) +
  geom_line(aes(y=devep_rate.max), color = "grey", size = 0.7, alpha=0.5) +
  geom_line(aes(y=devep_rate.mean), color = "blue", size = 0.7) +
  geom_ribbon(data = df_out_aeg_w,aes(ymin=devep_rate.max,
                                  ymax=devep_rate.min), fill="grey",
              alpha=0.5) +
  geom_point(data =  df_dE_aeg, aes(temp,develop_rate),
             size = 0.8, color = "black") +
  xlim(c(5,40)) +
  guides( color =FALSE, alpha = FALSE) +
  ylab("Egg development time") + xlab("Temperature (Cº)") +
  theme_bw()
plotdE_aeg_w 

# Join thermal responses all species --------------------
library(latex2exp)
sizelet = 12
ggarrange( plotalb_w  +
             theme(text = element_text(size = sizelet)) +
             ylab(TeX("Prob. from Larva to Adult, $p_{LA}$")) +
             xlab("") + ylim(c(0,1.3)) +
             ggtitle(expression(italic("Ae. albopictus"))) +
             ylab(TeX("Prob. from Larva to Adult, $p_{LA}$"))+ 
             theme(text = element_text(size = sizelet)),
           plotaeg_w  + 
             theme(text = element_text(size = sizelet)) +
             ylab(TeX("Prob. from Larva to Adult, $p_{LA}$")) +
             xlab("") + ylim(c(0,1.3))  +
             ggtitle(expression(italic("Ae. aegypti"))) +
             theme(text = element_text(size = sizelet)),
           plotdE_w + 
             ylab(TeX("Egg development rate, $d_E$")) +
             theme(text = element_text(size = sizelet)),
           plotdE_aeg_w + ylim(c(0,0.55)) + ylab("") +
             ylab(TeX("Egg development rate, $d_E$")) +
             theme(text = element_text(size = sizelet)),
           ncol = 2, nrow = 2)


# Larva development rate --------------------------------------

### --------------Development Egg Aegypti -----------------------#
## Delatte
  df_dL_alb <- data.frame(temp=c(15,20,25,30,35),
                          days_mean = c(35,14.4,10.4,8.8,12.3),
                          days_sd = c(0.9,0.4,0.7,0.6,0.7))
  ggplot(df_dL_alb) +
    geom_point(aes(temp,days_mean))
  
  n = 5
  r1 <- rnorm(n, df_dL_alb$days_mean[1],
              df_dL_alb$days_sd[1] )
  r2 <- rnorm(n, df_dL_alb$days_mean[2],
              df_dL_alb$days_sd[2] )
  r3 <- rnorm(n, df_dL_alb$days_mean[3],
              df_dL_alb$days_sd[3] )
  r4 <- rnorm(n, df_dL_alb$days_mean[4],
              df_dL_alb$days_sd[4] )
  r5 <- rnorm(n, df_dL_alb$days_mean[5],
              df_dL_alb$days_sd[5] )
 
  df_dL_alb <- data.frame(temp = sort(rep(df_dL_alb[,1],n)),
                          days_devep = c(r1,r2,r3,r4,r5))
 
  ggplot(df_dL_alb) +
    geom_point(aes(temp,develop_rate))
  
  df_dL_alb$develop_rate <- 1/(df_dL_alb$days_devep)
  df_dL_alb <- rbind(df_dL_alb, c(5,0,0,0))
  df_dL_alb <- rbind(df_dL_alb, c(10,0,0,0))
  ggplot(df_dL_alb) +
    geom_point(aes(temp,develop_rate))
  saveRDS(df_dL_alb, "~/INVASIBILITY_THRESHOLD/data/df_dL_alb.Rds")

## I use the random sample fixed, other wise the param would change.
  df_dL_alb <- readRDS("~/INVASIBILITY_THRESHOLD/data/df_dL_alb.Rds")

Fitting_dL_alb <- nls(develop_rate ~ c*temp*(temp-c1)*(c2-temp)^(1/2),
                      data = df_dL_alb, algorithm = "port",
                      start = list(c = 0.0003, c1 = 14, c2 = 40))

summary(Fitting_dL_alb)

Fitting_dL_alb_lin <- nls(develop_rate ~ c*temp + c1,
                          data = df_dL_alb,
                          start = list(c = 0.0003, c1 = 0))

summary(Fitting_dL_alb_lin)

AIC(Fitting_dL_alb, Fitting_dL_alb_lin)

mod <- function(te){
  c <- as.numeric(Fitting_dL_alb$m$getPars()[1])
  c1 <- as.numeric(Fitting_dL_alb$m$getPars()[2])
  c2 <- as.numeric(Fitting_dL_alb$m$getPars()[3])
  c*te*(te-c1)*(c2-te)^(1/2)
}

vec <- seq(0,45,0.001)
df_out_alb <- data.frame(temp = vec,
                         devep_rate = sapply(vec, mod))


# df_out_aeg[which(df_out_aeg$devep_rate < 0),2] <- 0

plotdL_alb <- ggplot(df_out_alb) +
  geom_line(aes(temp,devep_rate), size = 0.7) +
  geom_point(data =  df_dL_alb, aes(temp,develop_rate),
             size = 0.9, color = "black") +
  xlim(c(5,40)) +
  guides( color =FALSE, alpha = FALSE) +
  ylab("Larva development time") + xlab("Temperature (Cº)") +
  theme_bw()
plotdL_alb
