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
source("~/INVASIBILITY_THRESHOLD/code/funcR0.R")
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
a_df_alb <- function(temp){Briere_df(0.000193,10.25,38.32,temp)} # Biting rate
TFD_df_alb <- function(temp){Briere_df(0.0488,8.02,35.65,temp)} # Fecundity
pEA_df_alb <- function(temp){Quad_df(0.002663,6.668,38.92,temp)} # Survival probability Egg-Adult
lf_df_alb <- function(temp){Quad_df(1.43,13.41,31.51,temp)} # Adult life span
dE_df_alb <- function(temp){Quad_df(0.00071,1.73,40.51,temp)} # Mosquito Development Rate

R0_dfunc_alb <- function(rain,hum,Te,var){
  a <- a_f_alb(Te)
  f <- TFD_f_alb(Te)
  deltaa <- lf_f_alb(Te)
  probla <- pEA_f_alb(Te)
  dE <- dE_f_alb(Te)
  h <- h_f(hum,rain)
  deltE = 0.1
  R0 <- f*deltaa*a*probla*((h*dE)/(h*dE+deltE))
  dffT <- TFD_df_alb(Te)
  dfaT <- a_df_alb(Te)
  dfdeltaAT <- lf_df_alb(Te)
  dfplaT <- pEA_df_alb(Te)
  dfdET <- dE_df_alb(Te)
  dffR0 <- (1/3)*((R0)^(-2/3))*((deltaa*a*h*dE*probla)/(h*dE+deltE))*dffT
  dfaR0 <- (1/3)*((R0)^(-2/3))*((deltaa*f*h*dE*probla)/(h*dE+deltE))*dfaT
  dfdeltAR0 <- (1/3)*((R0)^(-2/3))*((f*a*h*dE*probla)/(h*dE+deltE))*dfdeltaAT
  dfpLAR0 <- (1/3)*((R0)^(-2/3))*((deltaa*a*h*dE*f)/(h*dE+deltE))*dfplaT
  dfdER0 <- (1/3)*((R0)^(-2/3))*((deltaa*a*f*
                                    probla)*((h*(h*dE+deltE)- (h*dE*h))/(h*dE+deltE)^2))*dfdET
  dfR0 <- dffR0 + dfaR0 + dfdeltAR0 + dfpLAR0 + dfdER0
  dfR0 <- ifelse(var == "RM",dfR0,
                 ifelse(var == "a",dfaR0,
                        ifelse(var == "f",dffR0,
                               ifelse(var == "deltaA",dfdeltAR0,
                                      ifelse(var == "pLA",dfpLAR0,dfdER0)))))
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
var_list = c("RM","a","f","deltaA","pLA", "dE")

# compute the derivative for each param ----------------
df_dT <- data.frame()
for(i in c(1:length(var_list))){
  R0df_ana <- sapply(vec, function(x){R0_dfunc_alb(rain_cte,hum_cte,
                                                   x,var_list[i])} )
  df_RM <- data.frame(vec= vec, out =R0df_ana)
  df_RM$var <- var_list[i]
  df_dT <- rbind(df_dT,df_RM)
}

# Add artificially to get color in the legend ---------------------
# df_out_deltaL <- data.frame(vec = 5,
#                             out = 0,
#                             var = "deltaL")
# df_out_dL <- data.frame(vec = 5,
#                         out = 0,
#                         var = "dL")
# df_dT <- rbind(df_dT,df_out_dL,df_out_deltaL)

# Plot all curves together -----------------------------------------
library(RColorBrewer)
name_pal = "Dark2"
display.brewer.pal(8, name_pal)
pal <- brewer.pal(8, name_pal)

col_a = pal[1]
col_f = pal[2]
col_lf = pal[3]
# col_deltaL = pal[4]
# col_dL = pal[5]
col_dE = pal[4]
col_pLA = pal[6]
col_deltaA = pal[8]
col_R = "#000000"

# Plot
df_alb <- ggplot(df_dT) +
  geom_line(aes(vec,out, color =var), size = 1) +
  ylim(c(-10,10)) + theme_bw() +
  xlab("Temperature") + ylab(TeX("Derivative, $dR_M/dT$")) +
  scale_color_manual(name = "",
                     values = c(col_a,col_dE,col_deltaA,
                                col_f,col_pLA,col_R),
                     labels = c("a", expression(d[E]),
                                expression(delta[A]),
                                "f", expression(paste(p[LA])),
                                expression(R[M]))) +
  theme(legend.key.size = unit(0.5, 'cm'),
        legend.key.width = unit(1.5, 'cm'),
        legend.key.height = unit(0.5, 'cm'),
        legend.position = c(0.1,0.8)) 

# Aegypti ----------------------------------------------------------
a_f_aeg <- function(temp){Briere_func(0.000202,13.35,40.08,temp)} # Biting rate
pEA_f_aeg <- function(temp){Quad_func(0.004186,9.373,40.26,temp)} # Survival probability Egg-Adult
lf_f_aeg <- function(temp){Quad_func(0.148,9.16,37.73,temp)} # Adult life span
dE_f_aeg <- function(temp){Briere_func(0.0003775 ,14.88,37.42,temp)} # Adult life span
a_df_aeg <- function(temp){Briere_df(0.000202,13.35,40.08,temp)} # Biting rate
pEA_df_aeg <- function(temp){Quad_df(0.004186,9.373,40.26,temp)} # Survival probability Egg-Adult
lf_df_aeg <- function(temp){Quad_df(0.148,9.16,37.73,temp)} # Adult life span
dE_df_aeg <- function(temp){Briere_df(0.0003775 ,14.88,37.42,temp)} # Mosquito Development Rate

R0_dfunc_aeg <- function(rain,hum,Te,var){
  a <- a_f_aeg(Te)
  f <- 40
  deltaa <- lf_f_aeg(Te)
  probla <- pEA_f_aeg(Te)
  dE <- dE_f_aeg(Te)
  h <- h_f(hum,rain)
  deltE = 0.1
  R0 <- f*deltaa*a*probla*((h*dE)/(h*dE+deltE))
  dfaT <- a_df_aeg(Te)
  dfdeltaAT <- lf_df_aeg(Te)
  dfplaT <- pEA_df_aeg(Te)
  dfdET <- dE_df_aeg(Te)
  dfaR0 <- (1/3)*((R0)^(-2/3))*((deltaa*f*h*dE*probla)/(h*dE+deltE))*dfaT
  dfdeltAR0 <- (1/3)*((R0)^(-2/3))*((f*a*h*dE*probla)/(h*dE+deltE))*dfdeltaAT
  dfpLAR0 <- (1/3)*((R0)^(-2/3))*((deltaa*a*h*dE*f)/(h*dE+deltE))*dfplaT
  dfdER0 <- (1/3)*((R0)^(-2/3))*((deltaa*a*f*
                                    probla)*((h*(h*dE+deltE)- (h*dE*h))/(h*dE+deltE)^2))*dfdET
  dfR0 <- dfaR0 + dfdeltAR0 + dfpLAR0 + dfdER0
  dfR0 <- ifelse(var == "RM",dfR0,
                 ifelse(var == "a",dfaR0,
                        ifelse(var == "deltaA",dfdeltAR0,
                               ifelse(var == "pLA",dfpLAR0,dfdER0))))
  if(var == "deltaA"|var == "a"|var == "RM"|var == "dE"){
    
  }else{
    dfR0 <-ifelse(is.na(dfR0),0,dfR0)
  }
  
  return(dfR0)
}

vec = seq(5,40,0.0001)
var_list = c("RM","a","deltaA","pLA", "dE")
df_dT <- data.frame()
for(i in c(1:length(var_list))){
  R0df_ana <- sapply(vec, function(x){R0_dfunc_aeg(rain_cte,hum_cte,
                                                   x,var_list[i])} )
  df_RM <- data.frame(vec= vec, out =R0df_ana)
  df_RM$var <- var_list[i]
  df_dT <- rbind(df_dT,df_RM)
}

# Plot all curves together -----------------------------------------
ggplot(df_dT[df_dT$var == "dE",]) +
  geom_line(aes(vec,out, color =var), size =1)
df_aeg <- ggplot(df_dT) +
  geom_line(aes(vec,out, color =var), size =0.8) +
  ylim(c(-6,6)) + theme_bw() +
  xlab("Temperature") + ylab(TeX("Derivative, $dR_M/dT$")) +
  scale_color_manual(name = "",
                     values = c(col_a,col_dE,col_deltaA,col_pLA,col_R),
                     labels = c("a",TeX("$ d_E$"),TeX(" $ \\delta_A$"),
                                TeX( " $ p_{LA}$"), TeX( " $ R_M$") )) +
  theme(legend.key.size = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm')) 

# # Japonicus ----------------------------------------------------------
# dE_f_jap <- function(temp){Briere_func(0.0002859,6.360,35.53 ,temp)} # Mosquito Development Rate
# dL_f_jap <- function(temp){Briere_func(7.000e-05,9.705e+00,3.410e+01,temp)} # Survival probability Egg-Adult
# lf_f_jap <- function(temp){Lin_func(-2.5045,82.6525,temp)} # Adult life span
# deltaL_f_jap <- function(temp){QuadN_func(0.0021476,-0.0806067 ,1.0332455,temp)} # Adult life span
# dE_df_jap <- function(temp){Briere_df(0.0002859,6.360,35.53,temp)} # Biting rate
# dL_df_jap <- function(temp){Briere_df(7.000e-05,9.705e+00,3.410e+01,temp)} # Survival probability Egg-Adult
# lf_df_jap <- function(temp){Lin_df(-2.5045,82.6525,temp)} # Adult life span
# deltaL_df_jap <- function(temp){QuadN_df(0.0021476,-0.0806067 ,1.0332455,temp)} # Mosquito Development Rate
#
# R0_dfunc_jap <- function(rain,hum,Te,var){
#   a <- 0.35
#   f <- 40 #183/2
#   deltaa <- lf_f_jap(Te)
#   deltaL <- deltaL_f_jap(Te)
#   deltE = 0.1
#   dE <- dE_f_jap(Te)
#   dL <- dL_f_jap(Te)
#   h <- h_f(hum,rain)
#   R0 <- ((f*a*deltaa)*(dL/(dL+deltaL))*((h*dE)/(h*dE+deltE)))
#   dfdeltaAT <- lf_df_jap(Te)
#   dfdeltaLT <- deltaL_df_jap(Te)
#   dfdET <- dE_df_jap(Te)
#   dfdLT <- dL_df_jap(Te)
#   
#   dfdeltaAR0 <- (1/3)*((R0)^(-2/3))*((f*a)*(dL/(dL+deltaL))*((h*dE)/(h*dE+deltE)))*dfdeltaAT
#   dfdLR0 <- (1/3)*((R0)^(-2/3))*((f*a*deltaa)*(((dL+deltaL)-dL)/(dL+deltaL)^2)*((h*dE)/(h*dE+deltE)))*dfdLT
#   dfdeltaLR0 <- (1/3)*((R0)^(-2/3))*((f*a*deltaa)*(-dL/(dL+deltaL)^2)*((h*dE)/(h*dE+deltE)))*dfdeltaLT
#   dfdER0 <- (1/3)*((R0)^(-2/3))*(deltaa*a*f*
#                                    (dL/(dL+deltaL))*((h*(h*dE+deltE)- (h*dE*h))/(h*dE+deltE)^2))*dfdET
#   dfR0 <- dfdeltaAR0 + dfdeltaLR0 + dfdLR0 + dfdER0
#   dfR0 <- ifelse(var == "RM",dfR0,
#                  ifelse(var == "deltaA",dfdeltaAR0,
#                         ifelse(var == "deltaL",dfdeltaLR0,
#                                ifelse(var == "dL",dfdLR0,dfdER0))))
#   if(var == "RM" | var == "dL"){
#     
#   }else{
#     dfR0 <-ifelse(is.na(dfR0),0,dfR0)
#   }
#   return(dfR0)
# }
#
# vec = seq(5,40,0.01)
# var_list = c("RM","deltaA","dL", "dE","deltaL")
# df_dT <- data.frame()
# for(i in c(1:length(var_list))){
#   R0df_ana <- sapply(vec, function(x){R0_dfunc_jap(rain_cte,hum_cte,
#                                                    x,var_list[i])} )
#   df_RM <- data.frame(vec= vec, out =R0df_ana)
#   df_RM$var <- var_list[i]
#   df_dT <- rbind(df_dT,df_RM)
# }
#
# # Plot all curves together -----------------------------------------
# df_jap <- ggplot(df_dT) +
#   geom_line(aes(vec,out, color =var), size = 1) +
#   ylim(c(-1,1)) + theme_bw() +
#   xlab("Temperature") + ylab(TeX("Derivative, $dR_M/dT$")) +
#   scale_color_manual(name = "",
#                      values = c(col_dE,col_deltaA,col_deltaL,col_dL,col_R)) +
#   theme(legend.key.size = unit(1, 'cm'),
#         legend.key.width = unit(1, 'cm')) 

# Join all plots ---------------------------------------------------
sizelet = 13
legend_only <- get_legend(df_alb +
                            theme(legend.position = "right",
                                  legend.text = element_text(size = sizelet),
                                  text = element_text(size = sizelet),
                                  legend.key.size = unit(0.5, 'cm'),
                                  legend.key.width = unit(1.5, 'cm'),
                                  legend.key.height = unit(0.8, 'cm')))
# ggarrange(df_alb + ylab("dx/dT") + xlab("") +
#             ylim(c(-1,1.2)) +
#             ggtitle(expression(paste("A         ",italic("Ae. Albopitus"))))+
#             theme(legend.position = "none",
#                   text = element_text(size = sizelet)),
#           df_aeg + ylab("") + xlab("") +
#             ylim(c(-1.5,1)) +
#             ggtitle(expression(paste("B           ",
#                                      italic("Ae. Aegypti")))) +
#             theme(legend.position = "none",
#                   text = element_text(size = sizelet)),
#           df_jap  + ylab("dx/dT") +
#             ggtitle(expression(paste("C        ",
#                                      italic("Ae. Japonicus"))))+
#             ylim(c(-1,1)) +
#             theme(legend.position = "none",
#                   text = element_text(size = sizelet)),
#           legend_only,
#           ncol=2, nrow = 2)

ggarrange(df_alb + ylab("dx/dT") + xlab("") +
            ylim(c(-1.5,1.2)) +
            ggtitle(expression(paste("A         ",italic("Ae. Albopitus"))))+
            theme(text = element_text(size = sizelet)),
          df_aeg + ylab("") + xlab("") +
            ylim(c(-1.5,1.2)) +
            ggtitle(expression(paste("B           ",
                                     italic("Ae. Aegypti")))) +
            theme(legend.position = "none",
                  text = element_text(size = sizelet)),
          ncol=2, widths = c(1,1), common.legend = TRUE)