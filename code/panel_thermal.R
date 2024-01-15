## Compute the thermal responses for Aedes albopictus a aegypti from literature data
# compare also different functions and compute its AIC value
rm(list= ls())
library(thermPerf)
library(ggplot2)
library(tidyverse)
library(nls2)

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
# ## Info taken Table2:https://academic.oup.com/jme/article/46/1/33/902827?login=true
df_dE <- data.frame(temp = c(5,10,15,20,25,30,35),
                     develop_rate_mu = c(11,2,7.4,2.9,4.5,6.7,7.1),
                     develop_rate_sd = c(1.3,0,1.8,0.4,0.7,0.7,0.8)
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
plot_dE <- ggplot(df_dE) + 
  geom_point(aes(temp,develop_rate)) + theme_bw()
plot_dE

# Fitting_dE <- nls(develop_rate ~ c*temp*(temp-c1)*(c2-temp)^(1/2),
#                   data = df_dE,
#                   start = list(c = 0.001, c1 = 5 , c2 = 30))

# Fit the data with a specific functional form 
Fitting_dE <- nls(develop_rate ~ c*(temp-c1)*(temp - c2),
                  data = df_dE,
                  start = list(c = 0.001, c1 = 5 , c2 = 30))

summary(Fitting_dE)

Fitting_dE_lin <- nls(develop_rate ~ c*temp+c1,
                      data = df_dE,
                      start = list(c = 0.001, c1 = 0))

summary(Fitting_dE_lin)

# Compute the AIC function 
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

# Plot the data
ggplot(df_dE_aeg) + 
  geom_point(aes(temp,develop_rate))

Fitting_dE_aeg <- nls(develop_rate ~ c*temp*(temp-c1)*(c2-temp)^(1/2),
                      data = df_dE_aeg, algorithm = "port",
                      start = list(c = 0.0003, c1 = 14, c2 = 40), 
                      lower=c(0.00001,7,30), upper = c(0.04,15,50))

summary(Fitting_dE_aeg)

# Fit the data 
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
