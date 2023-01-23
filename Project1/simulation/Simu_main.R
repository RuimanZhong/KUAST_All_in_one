library(raster)
library(RandomFields)
library(spatstat)
library(stars)
library(sf)
library(ggplot2)
library(INLA)
library(maptools)
library(ggpubr)
setwd("~/Documents/All_in_one/All_in_one/Project1/simulation")
source('samplegenerator.r')
source("fnCreateMesh.R")
source('fnCheckInputsMelding.R')
source("fnCheckInputsDown.R")
source('fnPredictMelding.R')
source('fnPredictDown.R')
source('simulate_fun.R')
source('~/Documents/All_in_one/All_in_one/visualization.R')
xlim = c (0,1)
ylim = c(0,1)

xlim0 = c(0,1)
ylim0 = c(0,1)
by = 0.01


# simulation_domain -------------------------------------------------------


win = owin(xrange = xlim, yrange = ylim)
boundaryregion <- sf::st_as_sf(win)
x <- seq(from = xlim[1] + (by / 2), to = xlim[2] - (by / 2), by = by)
y <- seq(from = ylim[1] + (by / 2), to = ylim[2] - (by / 2), by = by)
coord <- expand.grid(x = x, y = y)
coop_sf <- sf::st_as_sf(coord, coords = c('x','y'))
dppoint <- coop_sf %>% st_join(boundaryregion, left = FALSE)

# simulation_scenarios ----------------------------------------------------


sig.err = 0

beta0 = 0
beta1 = 0

#####################################
# INI SIMULATIONS
#####################################

# Number iterations
N <- 100

# Common parameters
nu <- 1

intercept <- 10


# no_ocv
#----
#
# Scenario


# Scenario
sce <- "sce13"
var = 1; scale = 0.01; unif = F;beta = 2;
for(pnum in c(10,50,100)){
  for(anum in c(2,10)){
    scenario <- paste0(sce, "pnum", pnum, "anum", anum)
    fnSimulateCov(scenario, N, var, scale, nu, intercept, beta, pnum, anum, unif)
    fnCalculateMSE(scenario, N)
  }}

# Scenario
sce <- "sce14"
var = 1; scale = 0.1; unif = F;beta = 2;
for(pnum in c(10,50,100)){
  for(anum in c(2,10)){
    scenario <- paste0(sce, "pnum", pnum, "anum", anum)
    fnSimulateCov(scenario, N, var, scale, nu, intercept, beta, pnum, anum, unif)
    fnCalculateMSE(scenario, N)
  }}

# Scenario
sce <- "sce15"
var = 4; scale = 0.01; unif = F;beta = 2;
for(pnum in c(10,50,100)){
  for(anum in c(2,10)){
    scenario <- paste0(sce, "pnum", pnum, "anum", anum)
    fnSimulateCov(scenario, N, var, scale, nu, intercept, beta, pnum, anum, unif)
    fnCalculateMSE(scenario, N)
  }}

# Scenario
sce <- "sce16"
var = 4; scale = 0.1; unif = F;beta = 2;
for(pnum in c(10,50,100)){
  for(anum in c(2,10)){
    scenario <- paste0(sce, "pnum", pnum, "anum", anum)
    fnSimulateCov(scenario, N, var, scale, nu, intercept, beta, pnum, anum, unif)
    fnCalculateMSE(scenario, N)
  }}



# Results
# Figure 5, Figure 6, Figure 7, Figure 8
#----
a = 2
a= 10
N = 100

varlist = rep(c(1,1,4,4,1,1,4,4),2)
scalelist = rep(c(0.01, 0.1), 8)
betalist = rep(c(0,2),each = 8)
preferlist = rep(c(0,0,0,0,1,1,1,1),2)
index = c(1:4, 1:4, 5:8, 5:8)
final <- NULL
for(i in c(1:4, 9:12)){
  sce <- paste0('sce',i)
  var = varlist[i]
  scale = scalelist[i]
  beta = betalist[i]
  for(pnum in c(10, 50, 100)){
    for(anum in c(2,10)){
      scenario <- paste0(sce, "pnum", pnum, "anum", anum)
      results <- read.csv(paste0("results/mse/", scenario, "MSE", N, ".csv"))
      print(scenario)
      final <- rbind(final,
                       data.frame(MSE = results[,1], Point = as.factor(pnum), Method = 'Melding',Var = var,Prefer = 0, Scenario = paste(anum^2, 'area','scale =', scale,'betac =',beta, paste0('(S',index[i],')'))),
                       data.frame(MSE = results[,2], Point = as.factor(pnum), Method = 'Downscaler', Var = var,Prefer = 0, Scenario = paste(anum^2, 'area','scale =', scale,'betac =',beta, paste0('(S',index[i],')'))))
      
    }}
}

final_1 <- dplyr::filter(final, Var == 1,Prefer== 0)
fnMSEplot(final_1)
final_2 <- dplyr::filter(final, Var == 4,Prefer== 0)
fnMSEplot(final_2)
final <- NULL
for(i in c(5:8, 13:16)){
  sce <- paste0('sce',i)
  var = varlist[i]
  scale = scalelist[i]
  beta = betalist[i]
  for(pnum in c(10, 50, 100)){
    for(anum in c(2,10)){
      scenario <- paste0(sce, "pnum", pnum, "anum", anum)
      results <- read.csv(paste0("results/mse/", scenario, "MSE", N, ".csv"))
      print(scenario)
      final <- rbind(final,
                     data.frame(MSE = results[,1],Method = results[,2], Point = as.factor(pnum), Prefer = 1,Var = var, Scenario = paste(anum^2, 'area','scale =', scale,'betac =',beta, paste0('(S',index[i],')') )))
      
    }}
}
final_m <- final %>% dplyr:: filter(Method == 'Meld')%>% mutate(Method = 'Melding')
final <- final%>% dplyr:: filter(Method == 'Downscaler')%>% rbind(final_m)
final_3 <- dplyr::filter(final, Var == 1)
g = fnMSEplot(final_3)
g 
final_4 <- dplyr::filter(final, Var == 4)
fnMSEplot(final_4)

final = rbind(final_1,final_2, final_3,final_4)
final_11 <- dplyr::filter(final, Var == 1, Method == 'Melding') %>% mutate(Prefer = as.factor(Prefer))
final_12 <- dplyr::filter(final, Var == 1, Method == 'Downscaler') %>% mutate(Prefer = as.factor(Prefer))
fnPreplot(final_11)
fnPreplot(final_12)

# Figure 1

# generating surface from simulate_fun.R
theme_set(theme_minimal())
theme_replace(text = element_text(size = 15))
fnfigure1(r,p1,a1,meld[[1]], 'Prediction')

# Figure 10

