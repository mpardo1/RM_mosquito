## Compute the thermal responses for Aedes aegypti and albopictus from literature data
# compare also different functions and compute its AIC value
rm(list= ls())
library(ggplot2)
library(tidyverse)
# library(nls2)

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
df_dE <- readRDS( "~/RM_mosquito/data/df_dE.Rds")
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
df_dE_aeg <- readRDS("~/RM_mosquito/data/df_dE_aeg.Rds")
ggplot(df_dE_aeg) + 
  geom_point(aes(temp,develop_rate))

# Do a briere fitting
Fitting_dE_aeg <- nls(develop_rate ~ c*temp*(temp-c1)*(c2-temp)^(1/2),
                      data = df_dE_aeg, algorithm = "port",
                      start = list(c = 0.0003, c1 = 14, c2 = 40), 
                      lower=c(0.00001,7,30), upper = c(0.04,15,50))

summary(Fitting_dE_aeg)

# Do a linear fitting
Fitting_dE_aeg_lin <- nls(develop_rate ~ c*temp + c1,
                          data = df_dE_aeg,
                          start = list(c = 0.0003, c1 = 0))

summary(Fitting_dE_aeg_lin)

# Compute the AIC value to compare models
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

# Store in a data frame and add column with measure
df_out_aeg_min <- data.frame(temp = vec,
                             devep_rate = sapply(vec, mod_min))
df_out_aeg_min[which(df_out_aeg_min$devep_rate < 0),2] <- 0
df_out_aeg_min$group = "min"
df_out_aeg$group = "mean"

# Join all data frames
df_out_aeg <- rbind(df_out_aeg_max,
                    df_out_aeg_min,
                    df_out_aeg)

# Plot results
# plotdE_aeg <- ggplot(df_out_aeg) +
#   geom_line(aes(temp,devep_rate,
#                 color = group,
#                 group = group, 
#                 alpha = group), size = 0.7) +
#   geom_point(data =  df_dE_aeg, aes(temp,develop_rate),
#              size = 0.9, color = "black") +
#   scale_color_manual(values=c("red", "blue", "red")) + 
#   scale_alpha_manual(values = c(0.5,1,0.5)) +
#   xlim(c(5,40)) +
#   guides( color =FALSE, alpha = FALSE) +
#   ylab("Egg development time") + xlab("Temperature (Cº)") +
#   theme_bw()
# plotdE_aeg

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

# Development Larvae albopictus -----------------------------------------------
## Delatte.https://pubmed.ncbi.nlm.nih.gov/19274388/
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


# Egg mortality rate ------------------------------------------
# Data from paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2944657/
# Aedes Albopictus deltaE -------------------------------------
Egg_mortality_alb <- data.frame(
  month = numeric(0), temp = numeric(0), rel_hum = numeric(0), prop_mort =numeric(0)
)

# Plot month 1
# Egg_mortality_alb[nrow(Egg_mortality_alb)+1,] <- c(1, 22,95,0.35)
# Egg_mortality_alb[nrow(Egg_mortality_alb)+1,] <- c(1, 24,95,0.05)
# Egg_mortality_alb[nrow(Egg_mortality_alb)+1,] <- c(1, 26,95,0.1)
Egg_mortality_alb[nrow(Egg_mortality_alb)+1,] <- c(1, 22,75,0.4)
Egg_mortality_alb[nrow(Egg_mortality_alb)+1,] <- c(1, 24,75,0.25)
Egg_mortality_alb[nrow(Egg_mortality_alb)+1,] <- c(1, 26,75,0.2)
# Egg_mortality_alb[nrow(Egg_mortality_alb)+1,] <- c(1, 22,55,0.1)
# Egg_mortality_alb[nrow(Egg_mortality_alb)+1,] <- c(1, 24,55,0.2)
# Egg_mortality_alb[nrow(Egg_mortality_alb)+1,] <- c(1, 26,55,0.39)
# Egg_mortality_alb[nrow(Egg_mortality_alb)+1,] <- c(1, 22,25,0.61)
# Egg_mortality_alb[nrow(Egg_mortality_alb)+1,] <- c(1, 24,25,0.38)
# Egg_mortality_alb[nrow(Egg_mortality_alb)+1,] <- c(1, 26,25,0.62)

ggplot(Egg_mortality_alb) +
  geom_point(aes(temp,prop_mort, 
                 color=as.factor(rel_hum),
                 shape = as.factor(month))) + theme_bw()

# Egg mortality albo 
# Paper https://www.scielo.br/j/rsp/a/dvPQ8QMr7Y687hPJxsTxjDg/abstract/?lang=en#
Egg_mortality_alb_2 <- data.frame(
  temp = numeric(0), rel_hum = character(0), prop_mort =numeric(0)
)

# Plot month 1
Egg_mortality_alb_2[nrow(Egg_mortality_alb_2)+1,] <- c(15, "70-85",1-0.66)
Egg_mortality_alb_2[nrow(Egg_mortality_alb_2)+1,] <- c(20,"70-85",1-0.82)
Egg_mortality_alb_2[nrow(Egg_mortality_alb_2)+1,] <- c(25, "70-85",1-0.735)
Egg_mortality_alb_2[nrow(Egg_mortality_alb_2)+1,] <- c(30, "70-85",1-0.80)

ggplot(Egg_mortality_alb_2) +
  geom_point(aes(temp,prop_mort)) + theme_bw()


# Egg mortality albo paper 2
# Paper Table 2 https://oaktrust.library.tamu.edu/bitstream/handle/1969.1/ETD-TAMU-2508/DICKERSON-DISSERTATION.pdf?sequence=1&isAllowed=y
Egg_mortality_alb_3 <- data.frame(
  temp = numeric(0), rel_hum = character(0), prop_mort =numeric(0)
)

# Plot month 1
Egg_mortality_alb_3[nrow(Egg_mortality_alb_3)+1,] <- c(15, NA ,1-0.56)
Egg_mortality_alb_3[nrow(Egg_mortality_alb_3)+1,] <- c(21,NA ,1-0.75)
Egg_mortality_alb_3[nrow(Egg_mortality_alb_3)+1,] <- c(27,NA ,1-0.74)
Egg_mortality_alb_3[nrow(Egg_mortality_alb_3)+1,] <- c(32,NA ,1-0.53)
Egg_mortality_alb_3[nrow(Egg_mortality_alb_3)+1,] <- c(35,NA ,1-0.52)

# Join data sets to fit thermal response
Egg_mortality_alb$paper <- "Juliano et al 2002, Oecologia."
Egg_mortality_alb_2$paper <- "Calado et al 2002, Scielo."
Egg_mortality_alb_3$paper <- "Dickerson et al 2007, Thesis"
Egg_mortality_alb <- Egg_mortality_alb[,c(2:5)] 
Egg_mortality_alb$rel_hum <- as.character(Egg_mortality_alb$rel_hum)
Egg_mort_alb <- rbind(Egg_mortality_alb, Egg_mortality_alb_2, Egg_mortality_alb_3)
Egg_mort_alb$prop_mort <- as.numeric(Egg_mort_alb$prop_mort)
Egg_mort_alb$temp <- as.numeric(Egg_mort_alb$temp)
ggplot(Egg_mort_alb) +
  geom_point(aes(temp, prop_mort, color = paper)) + 
  theme_bw()

# Do the fitting with nls
Fitting_alb_quad <- nls(prop_mort ~ (cont*temp*temp+Tmin*temp+Tmax),
                   data = Egg_mort_alb,
                   start = list(cont = 0.001, Tmin = 0.001, Tmax = 0.001))

summary(Fitting_alb_quad)

# Fit for a linear model
Fitting_alb <- nls(prop_mort ~ cont*temp+cont1,
                       data = Egg_mort_alb,
                       start = list(cont = 0.001, cont1 = 0))

summary(Fitting_alb)

# Compute the AIC value for the two models
AIC(Fitting_alb,Fitting_alb_quad)

# Plot the model
mod <- function(te){
  c <- as.numeric(Fitting_alb$m$getPars()[1])
  t0 <- as.numeric(Fitting_alb$m$getPars()[2])
  c*te+t0
}

# Plot the model Quad
mod <- function(te){
  c <- as.numeric(Fitting_alb_quad$m$getPars()[1])
  t0 <- as.numeric(Fitting_alb_quad$m$getPars()[2])
  t1 <- as.numeric(Fitting_alb_quad$m$getPars()[3])
  c*te^2+t0*te+t1
}

vec <- seq(0,45,0.01)
df_out_alb <- data.frame(temp_ae = vec, life_span_ae <- sapply(vec, mod))
colnames(df_out_alb) <- c("temp_ae", "deltaE")
df_out_alb$group <- "mean"

ggplot(df_out_alb) +
  geom_line(aes(temp_ae, deltaE)) +
  geom_point(data = Egg_mort_alb, aes(temp, prop_mort), color = "red") + theme_bw()

# Check minimum
df_out_alb[(df_out_alb$deltaE == min(df_out_alb$deltaE)), "temp_ae"]

# +/- SD---------------------
## Mean - SD
mod_min <- function(te){
  t0 <- as.numeric(Fitting_alb$m$getPars()[2]) - 
    summary(Fitting_alb)$coefficients[2,2]
  c <- as.numeric(Fitting_alb$m$getPars()[1]) - 
    summary(Fitting_alb)$coefficients[1,2]
  c*te+t0
}

mod_min <- function(te){
  t0 <- as.numeric(Fitting_alb_quad$m$getPars()[2]) - 
    summary(Fitting_alb_quad)$coefficients[2,2]
  c <- as.numeric(Fitting_alb_quad$m$getPars()[1]) - 
    summary(Fitting_alb_quad)$coefficients[1,2]
  t1 <- as.numeric(Fitting_alb_quad$m$getPars()[3]) - 
    summary(Fitting_alb_quad)$coefficients[3,2]
  c*te^2+t0*te+t1
}

vec <- seq(0,45,0.01)
df_out_alb_min <- data.frame(temp_ae = vec,
                             deltaE = sapply(vec, mod_min))
df_out_alb_min$group <- "min"

plot(df_out_alb_min$temp_ae,df_out_alb_min$deltaE)

## Mean - SD
mod_max <- function(te){
  t0 <- as.numeric(Fitting_alb$m$getPars()[2]) + 
    summary(Fitting_alb)$coefficients[2,2]
  c <- as.numeric(Fitting_alb$m$getPars()[1]) + 
    summary(Fitting_alb)$coefficients[1,2]
  c*te+t0
}

mod_max <- function(te){
  t0 <- as.numeric(Fitting_alb_quad$m$getPars()[2]) + 
    summary(Fitting_alb_quad)$coefficients[2,2]
  c <- as.numeric(Fitting_alb_quad$m$getPars()[1]) +
    summary(Fitting_alb_quad)$coefficients[1,2]
  t1 <- as.numeric(Fitting_alb_quad$m$getPars()[3]) + 
    summary(Fitting_alb_quad)$coefficients[3,2]
  c*te^2+t0*te+t1
}
vec <- seq(0,45,0.01)
df_out_alb_max <- data.frame(temp_ae = vec,
                             deltaE= sapply(vec, mod_max))
df_out_alb_max$group <- "max"
plot(df_out_alb_max$temp_ae,df_out_alb_max$deltaE)

# Plot all three curves together
df_out_alb <- rbind(df_out_alb_max,
                    df_out_alb_min,
                    df_out_alb)

df_out_w <- reshape(df_out_alb, idvar ="temp_ae" ,
                    timevar = "group", direction = "wide")
plotdeltaE_w <- ggplot(df_out_w, aes(x=temp_ae,y=deltaE.max)) +
  geom_line(aes(y=deltaE.min), color = "grey", size = 0.7, alpha=0.5) +
  geom_line(aes(y=deltaE.max), color = "grey", size = 0.7, alpha=0.5) +
  geom_line(aes(y=deltaE.mean), color = "blue", size = 0.7) +
  geom_ribbon(data = df_out_w,aes(ymin=deltaE.max,
                                  ymax=deltaE.min), fill="grey",
              alpha=0.5) +
  geom_point(data = Egg_mort_alb,aes(temp,prop_mort),
             size = 0.9, color = "black") +
  xlim(c(5,37)) + ylim(c(-2,3)) +
  guides( color =FALSE, alpha = FALSE) +
  ylab("Egg mortality rate") + xlab("Temperature (Cº)") +
  theme_bw() 
plotdeltaE_w 

# Aedes Aegypti deltaE ------------------------------------------------
Egg_mortality_aeg <- data.frame(
  month = numeric(0), temp = numeric(0), rel_hum = numeric(0), prop_mort =numeric(0)
)

# Plot month 1
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(1, 22,95,0.02)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(1, 24,95,0.3)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(1, 26,95,0.01)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(1, 22,75,0.05)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(1, 24,75,0.03)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(1, 26,75,0.05)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(1, 22,55,0.02)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(1, 24,55,0.01)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(1, 26,55,0.1)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(1, 22,25,0.08)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(1, 24,25,0.02)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(1, 26,25,0.05)

# Plot month 2
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(2, 22,95,0)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(2, 24,95,0.01)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(2, 26,95,0.01)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(2, 22,75,0)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(2, 24,75,0)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(2, 26,75,0.05)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(2, 22,55,0)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(2, 24,55,0.01)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(2, 26,55,0.05)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(2, 22,25,0.1)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(2, 24,25,0.01)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(2, 26,25,0.01)

# Month 3
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(3, 22,95,0.1)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(3, 24,95,0.35)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(3, 26,95,0.3)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(3, 22,75,0.1)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(3, 24,75,0.1)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(3, 26,75,0.35)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(3, 22,55,0.15)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(3, 24,55,0.2)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(3, 26,55,0.7)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(3, 22,25,0.55)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(3, 24,25,0.75)
Egg_mortality_aeg[nrow(Egg_mortality_aeg)+1,] <- c(3, 26,25,0.95)

ggplot(Egg_mortality_aeg) +
  geom_point(aes(temp,prop_mort, 
                 color=as.factor(rel_hum),
                 shape = as.factor(month))) + theme_bw()

# Join all in one df
Egg_mortality_alb$esp <- "Albopictus"
Egg_mortality_aeg$esp <- "Aegypti"
df_egg_mort <- rbind(Egg_mortality_alb,Egg_mortality_aeg)

# Egg mortality aeg
# Paper: https://academic.oup.com/jme/article/51/1/97/863390?login=false
Egg_mortality_aeg_2 <- data.frame(
  temp = numeric(0), rel_hum = character(0), prop_mort =numeric(0)
)

# Add manually values from table 2
Egg_mortality_aeg_2[nrow(Egg_mortality_aeg_2)+1,] <- c(12, "70-80",1-0.43)
Egg_mortality_aeg_2[nrow(Egg_mortality_aeg_2)+1,] <- c(14,"70-80",1-0.53)
Egg_mortality_aeg_2[nrow(Egg_mortality_aeg_2)+1,] <- c(16, "70-80",1-0.60)
Egg_mortality_aeg_2[nrow(Egg_mortality_aeg_2)+1,] <- c(18, "70-80",1-0.55)
Egg_mortality_aeg_2[nrow(Egg_mortality_aeg_2)+1,] <- c(20, "70-80",1-0.72)

ggplot(Egg_mortality_aeg_2) +
  geom_point(aes(temp,prop_mort)) + theme_bw()

# Egg mortality aeg
# Paper: https://www.scielo.br/j/mioc/a/6QjKLKZYL5Yr8KFfwGGqSPc/?format=pdf&lang=en#:~:text=Temperatures%20tested%20ranged%20between%2012,tropical%20and%20subtropical%20world%20regions
Egg_mortality_aeg_3 <- data.frame(
  temp = numeric(0), rel_hum = character(0), prop_mort =numeric(0)
)

# Add manually values from unique table 
Egg_mortality_aeg_3[nrow(Egg_mortality_aeg_3)+1,] <- c(16, "NA",1-0.81)
Egg_mortality_aeg_3[nrow(Egg_mortality_aeg_3)+1,] <- c(22,"NA",1-0.94)
Egg_mortality_aeg_3[nrow(Egg_mortality_aeg_3)+1,] <- c(25, "NA",1-0.96)
Egg_mortality_aeg_3[nrow(Egg_mortality_aeg_3)+1,] <- c(28, "NA",1-0.93)
Egg_mortality_aeg_3[nrow(Egg_mortality_aeg_3)+1,] <- c(31, "NA",1-0.83)
Egg_mortality_aeg_3[nrow(Egg_mortality_aeg_3)+1,] <- c(35, "NA",1-0.49)
Egg_mortality_aeg_3[nrow(Egg_mortality_aeg_3)+1,] <- c(36, "NA",1)

ggplot(Egg_mortality_aeg_3) +
  geom_point(aes(temp,prop_mort)) + theme_bw()

# Join data sets to fit thermal response
Egg_mortality_aeg$paper <- "Juliano et al 2002, Oecologia"
Egg_mortality_aeg_2$paper <- "Byttebier et al 2014, J.Med.Entomol."
Egg_mortality_aeg_3$paper <- "Farnesi et al 2009 Scielo."
Egg_mortality_aeg_1 <- Egg_mortality_aeg[Egg_mortality_aeg$rel_hum == 75,
                                         c("temp", "rel_hum", "prop_mort","paper")]
Egg_mortality_aeg_1$rel_hum <- as.character(Egg_mortality_aeg_1$rel_hum)
Egg_mort_aeg <- rbind(Egg_mortality_aeg_1, Egg_mortality_aeg_2, Egg_mortality_aeg_3)
Egg_mort_aeg$prop_mort <- as.numeric(Egg_mort_aeg$prop_mort)
Egg_mort_aeg$temp <- as.numeric(Egg_mort_aeg$temp)
ggplot(Egg_mort_aeg) +
  geom_point(aes(temp, prop_mort)) + 
  theme_bw()

# Do the fitting with nls
Fitting_aeg <- nls(prop_mort ~ (cont*temp*temp+Tmin*temp+Tmax),
                   data = Egg_mort_aeg,
                   start = list(cont = 0.001, Tmin = 0.001, Tmax = 0.001))

summary(Fitting_aeg)

# Fit for a linear model
Fitting_aeg_lin <- nls(prop_mort ~ cont*temp+cont1,
                       data = Egg_mort_aeg,
                       start = list(cont = 0.001, cont1 = 0))

summary(Fitting_aeg_lin)

# Compute the AIC value for the two models
AIC(Fitting_aeg,Fitting_aeg_lin)

# Plot the model
mod <- function(te){
  c <- as.numeric(Fitting_aeg$m$getPars()[1])
  t0 <- as.numeric(Fitting_aeg$m$getPars()[2])
  tm <- as.numeric(Fitting_aeg$m$getPars()[3])
  c*te*te+t0*te+tm
}

vec <- seq(0,45,0.01)
df_out_aeg <- data.frame(temp_ae = vec, life_span_ae <- sapply(vec, mod))
colnames(df_out_aeg) <- c("temp_ae", "deltaE")
df_out_aeg$group <- "mean"

# Check minimum
df_out_aeg[(df_out_aeg$deltaE == min(df_out_aeg$deltaE)), "temp_ae"]

###---------------+/- SD---------------------######
## Mean - SD
mod_min <- function(te){
  t0 <- as.numeric(Fitting_aeg$m$getPars()[2]) - 
    summary(Fitting_aeg)$coefficients[2,2]
  tm <- as.numeric(Fitting_aeg$m$getPars()[3]) - 
    summary(Fitting_aeg)$coefficients[3,2]
  c <- as.numeric(Fitting_aeg$m$getPars()[1]) - 
    summary(Fitting_aeg)$coefficients[1,2]
  c*te*te+t0*te+tm
}

vec <- seq(0,45,0.01)
df_out_aeg_min <- data.frame(temp_ae = vec,
                             deltaE = sapply(vec, mod_min))
df_out_aeg_min$group <- "min"

plot(df_out_aeg_min$temp_ae,df_out_aeg_min$deltaE)

## Mean - SD
mod_max <- function(te){
  t0 <- as.numeric(Fitting_aeg$m$getPars()[2]) + 
    summary(Fitting_aeg)$coefficients[2,2]
  tm <- as.numeric(Fitting_aeg$m$getPars()[3]) + 
    summary(Fitting_aeg)$coefficients[3,2]
  c <- as.numeric(Fitting_aeg$m$getPars()[1]) + 
    summary(Fitting_aeg)$coefficients[1,2]
  c*te*te+t0*te+tm
}

vec <- seq(0,45,0.01)
df_out_aeg_max <- data.frame(temp_ae = vec,
                            deltaE= sapply(vec, mod_max))
df_out_aeg_max$group <- "max"
plot(df_out_aeg_max$temp_ae,df_out_aeg_max$deltaE)

# Plot all three curves together
df_out_aeg <- rbind(df_out_aeg_max,
                    df_out_aeg_min,
                    df_out_aeg)

df_out_w <- reshape(df_out_aeg, idvar ="temp_ae" ,
                    timevar = "group", direction = "wide")
plotdeltaE_w_aeg <- ggplot(df_out_w, aes(x=temp_ae,y=deltaE.max)) +
  geom_line(aes(y=deltaE.min), color = "grey", size = 0.7, alpha=0.5) +
  geom_line(aes(y=deltaE.max), color = "grey", size = 0.7, alpha=0.5) +
  geom_line(aes(y=deltaE.mean), color = "blue", size = 0.7) +
  geom_ribbon(data = df_out_w,aes(ymin=deltaE.max,
                                  ymax=deltaE.min), fill="grey",
              alpha=0.5) +
  geom_point(data = Egg_mort_aeg,aes(temp,prop_mort),
             size = 0.9, color = "black") +
  xlim(c(5,37)) + ylim(c(-1.7,3.1)) + 
  guides( color =FALSE, alpha = FALSE) +
  ylab("Egg mortality rate") + xlab("Temperature (Cº)") +
  theme_bw()
plotdeltaE_w_aeg 

# Join thermal responses all species --------------------
library(latex2exp)
sizelet = 12
ggarrange( plotalb_w  +
             theme(text = element_text(size = sizelet)) +
             ylab(TeX("Prob. from Larva to Adult, $p_{LA}$")) +
             ylim(c(0,1.3)) +
             ggtitle(TeX("a) \\textit{Ae. albopictus}")) +
             ylab(TeX("Prob. from Larva to Adult, $p_{LA}$"))+ 
             theme(text = element_text(size = sizelet)),
           plotaeg_w  + 
             theme(text = element_text(size = sizelet)) +
             ylab(TeX("Prob. from Larva to Adult, $p_{LA}$")) +
              ylim(c(0,1.3))  +
             ggtitle(TeX("b) \\textit{Ae. aegypti}")) +
             theme(text = element_text(size = sizelet)),
           plotdeltaE_w + 
             ylab(TeX("Egg mortality rate, $\\delta_{E}$")) +
             ggtitle("c)") +
             theme(text = element_text(size = sizelet)),
           plotdeltaE_w_aeg  + 
             theme(text = element_text(size = sizelet)) +
             ylab(TeX("Egg mortality rate, $\\delta_{E}$")) +
             ggtitle("d)") +
             theme(text = element_text(size = sizelet)),
           plotdE_w + 
             ylab(TeX("Egg development rate, $d_E$")) +
             ggtitle("e)") +
             theme(text = element_text(size = sizelet)),
           plotdE_aeg_w + ylim(c(0,0.55)) + ylab("") +
             ggtitle("f)") +
             ylab(TeX("Egg development rate, $d_E$")) +
             theme(text = element_text(size = sizelet)),
           ncol = 2, nrow = 3)


