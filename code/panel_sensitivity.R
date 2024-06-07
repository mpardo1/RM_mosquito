### Code to do the sensitivity analysis for the R0 panel manuscipt
rm(list=ls())
library(tidyverse)
library(ggplot2)
library(dplyr)
library("ggpubr")
library(viridis)
library(gdata)
library("data.table")
library('numDeriv')
library("ggbreak")
library(latex2exp)
source("~/RM_mosquito/code/funcR0.R")

# dRM/dX*dX/dT; X each param -------------------------------------------
# Main functions
Briere_df <- function(cte, tmin, tmax, temp){
  outp <- cte*(2*temp - tmin)*(tmax - temp)^(1/2) - ((1/2)*cte*(temp^2-tmin*temp)*(tmax-temp)^(-(1/2)))
  
  return(outp)
}

Quad_df <- function(cte, tmin, tmax, temp){
  outp <- -cte*(2*temp - (tmax + tmin))
  
  return(outp)
}


Lin_df <- function(cte, tmin, tmax, temp){
  outp <- cte
  
  return(outp)
}

QuadN_df <- function(cte, tmin, tmax, temp){
  outp <- 2*cte*temp + tmin
  
  return(outp)
}
# For albopictus ----------------------------------------------
a_f_alb <- function(temp){Briere_func(0.000193,10.25,38.32,temp)} # Biting rate
TFD_f_alb <- function(temp){Briere_func(0.0488,8.02,35.65,temp)} # Fecundity
pEA_f_alb <- function(temp){Quad_func(0.002663,6.668,38.92,temp)} # Survival probability Egg-Adult
lf_f_alb <- function(temp){Quad_func(1.43,13.41,31.51,temp)} # Adult life span
dE_f_alb <- function(temp){Quad_func(0.00071,1.73,40.51,temp)} # Adult life span
deltaE_f_alb <- function(temp){QuadN_func(0.0019328,-0.091868,1.3338874,temp)} # Egg mortality rate
a_df_alb <- function(temp){Briere_df(0.000193,10.25,38.32,temp)} # Biting rate
TFD_df_alb <- function(temp){Briere_df(0.0488,8.02,35.65,temp)} # Fecundity
pEA_df_alb <- function(temp){Quad_df(0.002663,6.668,38.92,temp)} # Survival probability Egg-Adult
lf_df_alb <- function(temp){Quad_df(1.43,13.41,31.51,temp)} # Adult life span
dE_df_alb <- function(temp){Quad_df(0.00071,1.73,40.51,temp)} # Egg Development Rate
deltaE_df_alb <- function(temp){QuadN_df(0.0019328,-0.091868,1.3338874,temp)} # Egg mortality Rate

R0_dfunc_alb <- function(rain,hum,Te,var){
  a <- a_f_alb(Te)
  f <- TFD_f_alb(Te)
  deltaa <- lf_f_alb(Te)
  probla <- pEA_f_alb(Te)
  dE <- dE_f_alb(Te)
  h <- h_f(hum,rain)
  deltE <- deltaE_f_alb(Te)
  R0 <- f*deltaa*a*probla*((h*dE)/(h*dE+deltE))
  dffT <- TFD_df_alb(Te)
  dfaT <- a_df_alb(Te)
  dfdeltaAT <- lf_df_alb(Te)
  dfplaT <- pEA_df_alb(Te)
  dfdET <- dE_df_alb(Te)
  dfdeltaET <- deltaE_df_alb(Te)
  dffaR0 <- (1/3)*((R0)^(-2/3))*((deltaa*h*dE*probla)/(h*dE+deltE))*(dffT*a+f*dfaT) # derivative fa 
  dffR0 <- (1/3)*((R0)^(-2/3))*((deltaa*a*h*dE*probla)/(h*dE+deltE))*dffT
  dfaR0 <- (1/3)*((R0)^(-2/3))*((deltaa*f*h*dE*probla)/(h*dE+deltE))*dfaT
  dfdeltAR0 <- (1/3)*((R0)^(-2/3))*((f*a*h*dE*probla)/(h*dE+deltE))*dfdeltaAT
  dfpLAR0 <- (1/3)*((R0)^(-2/3))*((deltaa*a*h*dE*f)/(h*dE+deltE))*dfplaT
  dfdeltaER0 <- (1/3)*((R0)^(-2/3))*(-(deltaa*a*h*dE*f*probla)/(h*dE+deltE)^2)*dfdeltaET
  dfdER0 <- (1/3)*((R0)^(-2/3))*((deltaa*a*f*
                                    probla)*((h*(h*dE+deltE)- (h*dE*h))/(h*dE+deltE)^2))*dfdET
  dfR0 <- dffR0 + dfaR0 + dfdeltAR0 + dfpLAR0 + dfdER0
  dfR0 <- ifelse(var == "RM",dfR0,
                 ifelse(var == "a",dfaR0,
                        ifelse(var == "f",dffR0,
                               ifelse(var == "deltaA",dfdeltAR0,
                                      ifelse(var == "pLA",dfpLAR0,
                                             ifelse(var == "deltaE",dfpLAR0,
                                                    ifelse(var == "fa",dffaR0,dfdER0)))))))
  if(var == "deltaA" | var == "RM"){
    
  }else{
    dfR0 <-ifelse(is.na(dfR0),0,dfR0)
  }
  
  return(dfR0)
}

# set to cte the human density an rainfall -------------
rain_cte <- 8
hum_cte = 500
vec = seq(5,35,0.001)
var_list = c("RM","a","f","deltaA","pLA", "dE", "deltaE")
var_list = c("RM","fa","deltaA","pLA", "dE", "deltaE")

# compute the derivative for each param ----------------
df_dT <- data.frame()
for(i in c(1:length(var_list))){
  R0df_ana <- sapply(vec, function(x){R0_dfunc_alb(rain_cte,hum_cte,
                                                   x,var_list[i])} )
  df_RM <- data.frame(vec= vec, out =R0df_ana)
  df_RM$var <- var_list[i]
  df_dT <- rbind(df_dT,df_RM)
}

# Plot all curves together -----------------------------------------
library(RColorBrewer)
name_pal = "Dark2"
# display.brewer.pal(8, name_pal)
# pal <- brewer.pal(8, name_pal)
# pal <- c("#00798c","#d1495b", "#edae49", "#66a182",
#          "#1C0da2", "#8d96a3", "#000000")
pal <- c("#00798c","#d1495b", "#edae49", "#66a182",
          "#8d96a3", "#000000")
# col_a = pal[1]
# col_f = pal[2]
col_fa = pal[1]
col_deltaA = pal[3]
col_dE = pal[4]
col_deltaE = pal[5]
col_pLA = pal[2]
col_R = "#000000"


# Plot
df_alb <- ggplot(df_dT) +
  geom_line(aes(vec,out, color =var), size = 1) +
  ylim(c(-2,2)) + theme_bw() +
  xlab("Temperature") + ylab(TeX("Derivative, $dR_M/dT$")) +
  # scale_color_manual(name = "",
  #                    values = c(col_a,col_dE,col_deltaA,
  #                               col_deltaE,
  #                               col_f,col_pLA,col_R),
  #                    labels = c("a", expression(d[E]),
  #                               expression(delta[A]),expression(delta[E]),
  #                               "f", expression(paste(p[LA])),
  #                               expression(R[M]))) +
  scale_color_manual(name = "",
                     values = c(col_dE,col_deltaA,
                                col_deltaE,col_fa,
                                col_pLA,col_R),
                     labels = c(expression(d[E]),
                                expression(delta[A]),expression(delta[E]),
                                "fa", 
                                expression(paste(p[LA])),
                                expression(R[M]))) +
  theme(legend.key.size = unit(0.5, 'cm'),
        legend.key.width = unit(1.5, 'cm'),
        legend.key.height = unit(0.5, 'cm'),
        legend.position = c(0.1,0.8)) 

df_alb

# Aegypti ----------------------------------------------------------
EFD_f_aeg <- function(temp){Briere_func(0.00856,14.58,34.61,temp)} # Fecundity
pLA_f_aeg <- function(temp){Quad_func(0.004186,9.373,40.26,temp)} # Survival probability Egg-Adult
MDR_f_aeg <- function(temp){Briere_func(0.0000786,11.36,39.17,temp)} # Mosquito Development Rate
lf_f_aeg <- function(temp){Quad_func(0.148,9.16,37.73,temp)} # Adult life span
dE_f_aeg <- function(temp){Briere_func(0.0003775 ,14.88,37.42,temp)} # Egg development rate
deltaE_f_aeg <- function(temp){QuadN_func(0.004475,-0.210787,2.552370,temp)} # Egg mortality rate
a_df_aeg <- function(temp){Briere_df(0.000202,13.35,40.08,temp)} # Biting rate
pLA_df_aeg <- function(temp){Quad_df(0.004186,9.373,40.26,temp)} # Survival probability Egg-Adult
lf_df_aeg <- function(temp){Quad_df(0.148,9.16,37.73,temp)} # Adult life span
dE_df_aeg <- function(temp){Briere_df(0.0003775 ,14.88,37.42,temp)} # Mosquito Development Rate
deltaE_df_aeg <- function(temp){QuadN_df(0.004475,-0.210787,2.552370,temp)} # Mosquito Development Rate
EFD_df_aeg <- function(temp){Briere_df(0.00856,14.58,34.61,temp)} # Derivative fecundity

R0_dfunc_aeg <- function(rain,hum,Te,var){
  a <- 1#a_f_aeg(Te)
  f <- EFD_f_aeg(Te)
  deltaa <- lf_f_aeg(Te)
  probla <- pLA_f_aeg(Te)
  dE <- dE_f_aeg(Te)
  h <- h_f(hum,rain)
  deltE <- deltaE_f_aeg(Te)
  R0 <- f*deltaa*a*probla*((h*dE)/(h*dE+deltE))
  dfaT <- a_df_aeg(Te)
  dffaT <- EFD_df_aeg(Te)
  dfdeltaAT <- lf_df_aeg(Te)
  dfplaT <- pLA_df_aeg(Te)
  dfdET <- dE_df_aeg(Te)
  dfdeltaET <- deltaE_df_aeg(Te)
  dfaR0 <- (1/3)*((R0)^(-2/3))*((deltaa*f*h*dE*probla)/(h*dE+deltE))*dfaT
  dffaR0 <- (1/3)*((R0)^(-2/3))*((deltaa*h*dE*probla)/(h*dE+deltE))*dffaT
  dfdeltAR0 <- (1/3)*((R0)^(-2/3))*((f*a*h*dE*probla)/(h*dE+deltE))*dfdeltaAT
  dfpLAR0 <- (1/3)*((R0)^(-2/3))*((deltaa*a*h*dE*f)/(h*dE+deltE))*dfplaT
  dfdeltaER0 <- (1/3)*((R0)^(-2/3))*(-(deltaa*a*h*dE*f*probla)/(h*dE+deltE)^2)*dfdeltaET
  dfdER0 <- (1/3)*((R0)^(-2/3))*((deltaa*a*f*
                                    probla)*((h*(h*dE+deltE)- (h*dE*h))/(h*dE+deltE)^2))*dfdET
  dfR0 <- dfaR0 + dfdeltAR0 + dfpLAR0 + dfdER0
  dfR0 <- ifelse(var == "RM",dfR0,
                 ifelse(var == "a",dfaR0,
                        ifelse(var == "deltaA",dfdeltAR0,
                               ifelse(var == "pLA",dfpLAR0,
                                      ifelse(var == "deltaE",dfdeltaER0,
                                             ifelse(var == "fa",dffaR0,dfdER0))))))
  if(var == "deltaA"|var == "a"|var == "RM"|var == "dE"|var == "deltaE"){
    
  }else{
    dfR0 <-ifelse(is.na(dfR0),0,dfR0)
  }
  
  return(dfR0)
}

vec = seq(5,40,0.0001)
# var_list = c("RM","a","deltaA","pLA", "dE", "deltaE")
var_list = c("RM","fa","deltaA","pLA", "dE", "deltaE")

df_dT <- data.frame()
for(i in c(1:length(var_list))){
  R0df_ana <- sapply(vec, function(x){R0_dfunc_aeg(rain_cte,hum_cte,
                                                   x,var_list[i])} )
  df_RM <- data.frame(vec= vec, out =R0df_ana)
  df_RM$var <- var_list[i]
  df_dT <- rbind(df_dT,df_RM)
}

# Plot all curves together -----------------------------------------
# ggplot(df_dT[df_dT$var == "dE",]) +
#   geom_line(aes(vec,out, color =var), size =1)
df_aeg <- ggplot(df_dT) +
  geom_line(aes(vec,out, color =var), size =0.8) +
  ylim(c(-6,6)) + theme_bw() +
  xlab("Temperature") + ylab(TeX("Derivative, $dR_M/dT$")) +
  # scale_color_manual(name = "",
  #                    values = c(col_a,col_dE,col_deltaA,col_deltaE,col_pLA,col_R),
  #                    labels = c("a",TeX("$ d_E$"),TeX(" $ \\delta_A$"),TeX(" $ \\delta_E$"),
  #                               TeX( " $ p_{LA}$"), TeX( " $ R_M$") )) +
  scale_color_manual(name = "",
                     values = c(col_dE,col_deltaA,col_deltaE,col_fa,col_pLA,col_R),
                     labels = c(TeX("$ d_E$"),TeX(" $ \\delta_A$"),TeX(" $ \\delta_E$"),
                                "fa",
                                TeX( " $ p_{LA}$"), TeX( " $ R_M$") )) +
  theme(legend.key.size = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm')) 
df_aeg

# Sensitivy Fixed one param --------------------------------------------
source("~/INVASIBILITY_THRESHOLD/code/funcR0.R")

# Albopictus -----------------------------------------------------
vec <- seq(5,40,0.001)
albopictus <- sapply(vec,R0_func_alb, hum = 500,rain = 8) 


# Fixed bitting rate
cte <- max(sapply(seq(1,40,0.01),a_f_alb))

# R0 function by temperature:
R0_func_alb <- function(Te, rain, hum){
  if(is.na(Te) | is.na(rain) | is.na(hum)){
    R0 <- NA
  }else{
    a <- cte #a_f_alb(Te)
    f <- (1/2)*TFD_f_alb(Te)
    deltaa <- lf_f_alb(Te)
    dE <- dE_f_alb(Te)
    probla <- pLA_f_alb(Te)
    h <- h_f(hum,rain)
    deltaE = deltaE_f_alb(Te)#0.1
    
    R0 <- ((f*a*deltaa)*probla*((h*dE)/(h*dE+deltaE)))^(1/3)
  }
  return(R0)
}

# Run model with fixed param
vec <- seq(5,40,0.001)
albopictus_a_cte <- sapply(vec,R0_func_alb, hum = 500,rain = 8) 

# Fixed fa
cte <- max(sapply(seq(1,40,0.01),function(x){a_f_alb(x)*TFD_f_alb(x)}))

# R0 function by temperature:
R0_func_alb <- function(Te, rain, hum){
  if(is.na(Te) | is.na(rain) | is.na(hum)){
    R0 <- NA
  }else{
    a <- a_f_alb(Te)
    f <- (1/2)*TFD_f_alb(Te)
    deltaa <- lf_f_alb(Te)
    dE <- dE_f_alb(Te)
    probla <- pLA_f_alb(Te)
    h <- h_f(hum,rain)
    deltaE = deltaE_f_alb(Te)#0.1
    
    R0 <- ((cte*deltaa)*probla*((h*dE)/(h*dE+deltaE)))^(1/3)
  }
  return(R0)
}

# Run model with fixed param
vec <- seq(5,40,0.001)
albopictus_fa_cte <- sapply(vec,R0_func_alb, hum = 500,rain = 8) 

# Fixed bitting rate
cte <- (max(sapply(seq(1,40,0.01),TFD_f_alb)))

# R0 function by temperature:
R0_func_alb <- function(Te, rain, hum){
  if(is.na(Te) | is.na(rain) | is.na(hum)){
    R0 <- NA
  }else{
    a <- a_f_alb(Te)
    f <- (1/2)*cte#TFD_f_alb(Te)
    deltaa <- lf_f_alb(Te)
    dE <- dE_f_alb(Te)
    probla <- pLA_f_alb(Te)
    h <- h_f(hum,rain)
    deltaE = deltaE_f_alb(Te)#0.1
    
    R0 <- ((f*a*deltaa)*probla*((h*dE)/(h*dE+deltaE)))^(1/3)
  }
  return(R0)
}

# Run model with fixed param
vec <- seq(5,40,0.001)
albopictus_f_cte <- sapply(vec,R0_func_alb, hum = 500,rain = 8) 

# Fixed adult mortality rate
cte <- (max(sapply(seq(1,40,0.01),lf_f_alb)))

# R0 function by temperature:
R0_func_alb <- function(Te, rain, hum){
  if(is.na(Te) | is.na(rain) | is.na(hum)){
    R0 <- NA
  }else{
    a <- a_f_alb(Te)
    f <- (1/2)*TFD_f_alb(Te)
    deltaa <- cte#lf_f_alb(Te)
    dE <- dE_f_alb(Te)
    probla <- pLA_f_alb(Te)
    h <- h_f(hum,rain)
    deltaE = deltaE_f_alb(Te)#0.1
    
    R0 <- ((f*a*deltaa)*probla*((h*dE)/(h*dE+deltaE)))^(1/3)
  }
  return(R0)
}

# Run model with fixed param
vec <- seq(5,40,0.001)
albopictus_lf_cte <- sapply(vec,R0_func_alb, hum = 500,rain = 8) 

# Fixed Egg development rate
cte <- (max(sapply(seq(1,40,0.01),dE_f_alb)))

# R0 function by temperature:
R0_func_alb <- function(Te, rain, hum){
  if(is.na(Te) | is.na(rain) | is.na(hum)){
    R0 <- NA
  }else{
    a <- a_f_alb(Te)
    f <- (1/2)*TFD_f_alb(Te)
    deltaa <- lf_f_alb(Te)
    dE <- cte#dE_f_alb(Te)
    probla <- pLA_f_alb(Te)
    h <- h_f(hum,rain)
    deltaE = deltaE_f_alb(Te)#0.1
    
    R0 <- ((f*a*deltaa)*probla*((h*dE)/(h*dE+deltaE)))^(1/3)
  }
  return(R0)
}

# Run model with fixed param
vec <- seq(5,40,0.001)
albopictus_dE_cte <- sapply(vec,R0_func_alb, hum = 500,rain = 8) 

# Fixed adult probability from Larvae to adult
cte <- (max(sapply(seq(1,40,0.01),pLA_f_alb)))

# R0 function by temperature:
R0_func_alb <- function(Te, rain, hum){
  if(is.na(Te) | is.na(rain) | is.na(hum)){
    R0 <- NA
  }else{
    a <- a_f_alb(Te)
    f <- (1/2)*TFD_f_alb(Te)
    deltaa <- lf_f_alb(Te)
    dE <- dE_f_alb(Te)
    probla <- cte#pLA_f_alb(Te)
    h <- h_f(hum,rain)
    deltaE = deltaE_f_alb(Te)#0.1
    
    R0 <- ((f*a*deltaa)*probla*((h*dE)/(h*dE+deltaE)))^(1/3)
  }
  return(R0)
}

# Run model with fixed param
vec <- seq(5,40,0.001)
albopictus_pLA_cte <- sapply(vec,R0_func_alb, hum = 500,rain = 8) 

# Fixed adult Egg mortality rate
cte <- (min(sapply(seq(1,40,0.01),deltaE_f_alb)))

# R0 function by temperature:
R0_func_alb <- function(Te, rain, hum){
  if(is.na(Te) | is.na(rain) | is.na(hum)){
    R0 <- NA
  }else{
    a <- a_f_alb(Te)
    f <- (1/2)*TFD_f_alb(Te)
    deltaa <- lf_f_alb(Te)
    dE <- dE_f_alb(Te)
    probla <- pLA_f_alb(Te)
    h <- h_f(hum,rain)
    deltaE = cte#deltaE_f_alb(Te)#0.1
    
    R0 <- ((f*a*deltaa)*probla*((h*dE)/(h*dE+deltaE)))^(1/3)
  }
  return(R0)
}

# Run model with fixed param
vec <- seq(5,40,0.001)
albopictus_deltaE_cte <- sapply(vec,R0_func_alb, hum = 500,rain = 8) 

# Create data frame with all vecs
# df_cte <- data.frame(vec,  
#                      albopictus_a_cte, albopictus_dE_cte,
#                      albopictus_f_cte, albopictus_lf_cte,
#                      albopictus_pLA_cte, albopictus_deltaE_cte,
#                      albopictus)
df_cte <- data.frame(vec,  
                     albopictus_dE_cte,
                     albopictus_fa_cte, albopictus_lf_cte,
                     albopictus_pLA_cte, albopictus_deltaE_cte,
                     albopictus)
# colnames(df_cte) <- c("Temperature",  "a", "dE",
#                       "f", "deltaA", "pLA", "deltaE", "No cte")
colnames(df_cte) <- c("Temperature",  "dE",
                      "fa", "deltaA", "pLA", "deltaE", "No cte")
df_cte <- reshape2::melt( df_cte, id.vars = "Temperature")
library(RColorBrewer)
name_pal = "Set1"
# display.brewer.pal(7, name_pal)
# pal <- brewer.pal(7, name_pal)

letsize = 16
library("latex2exp")
plot_temp_alb <- ggplot(df_cte) + 
  geom_line(aes(Temperature,value, color=variable), size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  ylab(TeX("$R_M$")) + 
  # scale_color_manual(name = "",
  #                    values = c(col_a,col_dE,col_f,col_deltaA,col_deltaE,col_pLA,col_R),
  #                    labels = c("a",TeX("$ d_E$"),"f",TeX(" $ \\delta_A$"),TeX( " $ p_{LA}$"),
  #                               TeX(" $ \\delta_E$"), TeX( "Original") )) +
  scale_color_manual(name = "",
                     values = c(col_dE,col_fa,col_deltaA,col_pLA,col_deltaE,col_R),
                     labels = c(TeX("$ d_E$"),"fa",TeX(" $ \\delta_A$"),TeX( " $ p_{LA}$"),
                                TeX(" $ \\delta_E$"), TeX( "Original") )) +
  xlab("Temperature (Cº)") +
  theme_bw() + theme(legend.position = c(0.18,0.75),
                     text = element_text(size = letsize),
                     legend.text.align = 0)

plot_temp_alb

# Aegypti ---------------------------------------------------------------
source("~/RM_mosquito/code/funcR0.R") # Reload
vec <- seq(5,40,0.001)
aegypti <- sapply(vec,R0_func_aeg, hum = 500,rain = 8) 

# Fixed bitting rate
cte <- max(sapply(seq(1,40,0.01),a_f_aeg))

# R0 function by temperature:
R0_func_aeg <- function(Te, rain,hum){
  if(is.na(Te) | is.na(rain) | is.na(hum)){
    R0 <- NA
  }else{
    a <- cte#a_f_aeg(Te)
    f <- 40#EFD_f_aeg(Te) #40
    deltaa <- lf_f_aeg(Te)
    dE <- dE_f_aeg(Te)
    probla <- pLA_f_aeg(Te)
    h <- h_f(hum,rain)
    deltaE = deltaE_f_aeg(Te)
    R0 <- ((f*a*deltaa)*probla*((h*dE)/(h*dE+deltaE)))^(1/3)
  }
  return(R0)
}

# Run model with fixed param
vec <- seq(5,40,0.001)
aegypti_a_cte <- sapply(vec,R0_func_aeg, hum = 500,rain = 8) 

# Fixed fa
cte <- max(sapply(seq(1,40,0.01),EFD_f_aeg))

# R0 function by temperature:
R0_func_aeg <- function(Te, rain,hum){
  if(is.na(Te) | is.na(rain) | is.na(hum)){
    R0 <- NA
  }else{
    a <- 1#a_f_aeg(Te)
    f <- EFD_f_aeg(Te) #40
    deltaa <- lf_f_aeg(Te)
    dE <- dE_f_aeg(Te)
    probla <- pLA_f_aeg(Te)
    h <- h_f(hum,rain)
    deltaE = deltaE_f_aeg(Te)
    R0 <- ((cte*deltaa)*probla*((h*dE)/(h*dE+deltaE)))^(1/3)
  }
  return(R0)
}

# Run model with fixed param
vec <- seq(5,40,0.001)
aegypti_fa_cte <- sapply(vec,R0_func_aeg, hum = 500,rain = 8) 

# Fixed adult mortality rate
cte <- (max(sapply(seq(1,40,0.01),lf_f_aeg)))

# R0 function by temperature:
R0_func_aeg <- function(Te, rain,hum){
  if(is.na(Te) | is.na(rain) | is.na(hum)){
    R0 <- NA
  }else{
    a <- 1#a_f_aeg(Te)
    f <- EFD_f_aeg(Te) #40
    deltaa <- cte#lf_f_aeg(Te)
    dE <- dE_f_aeg(Te)
    probla <- pLA_f_aeg(Te)
    h <- h_f(hum,rain)
    deltaE = deltaE_f_aeg(Te)
    R0 <- ((f*a*deltaa)*probla*((h*dE)/(h*dE+deltaE)))^(1/3)
  }
  return(R0)
}

# Run model with fixed param
vec <- seq(5,40,0.001)
aegypti_lf_cte <- sapply(vec,R0_func_aeg, hum = 500,rain = 8) 

# Fixed Egg development rate
cte <- (max(sapply(seq(1,40,0.01),dE_f_aeg)))

# R0 function by temperature:
R0_func_aeg <- function(Te, rain,hum){
  if(is.na(Te) | is.na(rain) | is.na(hum)){
    R0 <- NA
  }else{
    a <- 1#a_f_aeg(Te)
    f <- EFD_f_aeg(Te) #40
    deltaa <- lf_f_aeg(Te)
    dE <- cte#dE_f_aeg(Te)
    probla <- pLA_f_aeg(Te)
    h <- h_f(hum,rain)
    deltaE = deltaE_f_aeg(Te)
    R0 <- ((f*a*deltaa)*probla*((h*dE)/(h*dE+deltaE)))^(1/3)
  }
  return(R0)
}

# Run model with fixed param
vec <- seq(5,40,0.001)
aegypti_dE_cte <- sapply(vec,R0_func_aeg, hum = 500,rain = 8) 

# Fixed adult probability from Larvae to adult
cte <- (max(sapply(seq(1,40,0.01),pLA_f_aeg)))

# R0 function by temperature:
R0_func_aeg <- function(Te, rain,hum){
  if(is.na(Te) | is.na(rain) | is.na(hum)){
    R0 <- NA
  }else{
    a <- 1#a_f_aeg(Te)
    f <- EFD_f_aeg(Te) #40
    deltaa <- lf_f_aeg(Te)
    dE <- dE_f_aeg(Te)
    probla <- cte#pLA_f_aeg(Te)
    h <- h_f(hum,rain)
    deltaE = deltaE_f_aeg(Te)
    R0 <- ((f*a*deltaa)*probla*((h*dE)/(h*dE+deltaE)))^(1/3)
  }
  return(R0)
}

# Run model with fixed param
vec <- seq(5,40,0.001)
aegypti_pLA_cte <- sapply(vec,R0_func_aeg, hum = 500,rain = 8) 

# Fixed adult Egg mortality rate
cte <- (min(sapply(seq(1,40,0.01),deltaE_f_aeg)))

# R0 function by temperature:
R0_func_aeg <- function(Te, rain,hum){
  if(is.na(Te) | is.na(rain) | is.na(hum)){
    R0 <- NA
  }else{
    a <- 1#a_f_aeg(Te)
    f <- EFD_f_aeg(Te) #40
    deltaa <- lf_f_aeg(Te)
    dE <- dE_f_aeg(Te)
    probla <- pLA_f_aeg(Te)
    h <- h_f(hum,rain)
    deltaE = cte#deltaE_f_aeg(Te)
    R0 <- ((f*a*deltaa)*probla*((h*dE)/(h*dE+deltaE)))^(1/3)
  }
  return(R0)
}

# Run model with fixed param
vec <- seq(5,40,0.001)
aegypti_deltaE_cte <- sapply(vec,R0_func_aeg, hum = 500,rain = 8) 

# Create data frame with all vecs
# df_cte <- data.frame(vec,  
#                      aegypti_a_cte, aegypti_dE_cte,
#                      aegypti_lf_cte,
#                      aegypti_pLA_cte, aegypti_deltaE_cte,
#                      aegypti)
# colnames(df_cte) <- c("Temperature",  "a", "dE",
#                       "deltaA", "pLA", "deltaE", "No cte")
df_cte <- data.frame(vec,  
                     aegypti_fa_cte, aegypti_dE_cte,
                     aegypti_lf_cte,
                     aegypti_pLA_cte, aegypti_deltaE_cte,
                     aegypti)
colnames(df_cte) <- c("Temperature",  "fa", "dE",
                      "deltaA", "pLA", "deltaE", "No cte")
df_cte <- reshape2::melt( df_cte, id.vars = "Temperature")

# Plot the results
letsize = 16
library("latex2exp")
plot_temp_aeg <- ggplot(df_cte) + 
  geom_line(aes(Temperature,value, color=variable), size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  ylab(TeX("$R_M$")) + 
  # scale_color_manual(name = "",
  #                    values = c(col_a,col_dE,col_deltaA,col_pLA,col_deltaE,col_R),
  #                    labels = c("a",TeX("$ d_E$"),TeX(" $ \\delta_A$"),
  #                               TeX(" $ \\delta_E$"),TeX( " $ p_{LA}$"), TeX( "Original") )) +
  scale_color_manual(name = "",
                     values = c(col_fa,col_dE,col_deltaA,col_pLA,col_deltaE,col_R),
                     labels = c("fa",TeX("$ d_E$"),TeX(" $ \\delta_A$"),
                                TeX( " $ p_{LA}$"),TeX(" $ \\delta_E$"), TeX( "Original") )) +
  xlab("Temperature (Cº)") +
  theme_bw() + theme(legend.position = c(0.18,0.75),
                     text = element_text(size = letsize),
                     legend.text.align = 0)

plot_temp_aeg

# Arrange a panel for SM ------------------------------------------------
sizelet = 13
legend_only <- get_legend(df_alb +
                            theme(legend.position = "right",
                                  legend.text = element_text(size = sizelet),
                                  text = element_text(size = sizelet),
                                  legend.key.size = unit(0.5, 'cm'),
                                  legend.key.width = unit(1.5, 'cm'),
                                  legend.key.height = unit(0.8, 'cm')))

ggarrange(df_alb +  xlab("") +
            ylim(c(-1.5,1.2)) +
            ggtitle(expression(paste("a) ",italic("Ae. Albopictus"))))+
            theme(text = element_text(size = sizelet)),
          df_aeg +  xlab("") +
            ylim(c(-1,1)) +
            ggtitle(expression(paste("b) ",
                                     italic("Ae. Aegypti")))) +
            theme(legend.position = "none",
                  text = element_text(size = sizelet)),
          plot_temp_alb + ggtitle("c)") +
            theme(text = element_text(size = sizelet)),
          plot_temp_aeg + ggtitle("d)") +
            theme(text = element_text(size = sizelet)),
          ncol=2,nrow = 2,  heights = c(1,0.8),
          common.legend = TRUE)

# Check difference when fecundity is EFD, not a*cte ---------------------
# R0 function by temperature:
R0_func_aeg <- function(Te, rain,hum){
  if(is.na(Te) | is.na(rain) | is.na(hum)){
    R0 <- NA
  }else{
    a <- 1#a_f_aeg(Te)
    f <- EFD_f_aeg(Te) #40
    deltaa <- lf_f_aeg(Te)
    dE <- dE_f_aeg(Te)
    probla <- pLA_f_aeg(Te)
    h <- h_f(hum,rain)
    deltaE = cte#deltaE_f_aeg(Te)
    R0 <- ((f*a*deltaa)*probla*((h*dE)/(h*dE+deltaE)))^(1/3)
  }
  return(R0)
}

# Run model with fixed param 
vec <- seq(5,40,0.001)
aegypti_EFD <- sapply(vec,R0_func_aeg, hum = 500,rain = 8) 

# Create data frame with all vecs
df_cte1 <- data.frame(vec, aegypti, aegypti_EFD)
df_cte1 <- reshape2::melt( df_cte1, id.vars = "vec")
ggplot(df_cte1) + 
  geom_line(aes(vec,value, color=variable), size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  ylab(TeX("$R_M$")) + theme_bw() +  xlab("Temperature") + 
  scale_color_manual(name= "", values = c("#edae49", "#66a182"),
                     labels = c(expression(italic("Cte fecundity")),
                                expression(italic("Fecundity depending on temp")))) +
  theme(legend.position = c(0.18,0.75),
        text = element_text(size = letsize),
        legend.text.align = 0)
