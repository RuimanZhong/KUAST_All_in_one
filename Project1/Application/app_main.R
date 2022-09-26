library(ggplot2)
library(stringr)
source('Application/fnCreateCV.R')
source('Application/fnPredictUK_wCov.R')
source('spatialM/fnPredictMelding.R')
source('spatialM/fnPredictDown.R')


k = 5
rep = 3
mse1 = NULL
mse2 = NULL
mse3 = NULL
mse4 = NULL
mse11 = NULL
mse21 = NULL
mse31 = NULL
mse41 = NULL
for(j in 1: rep){
  flag_SCV = fnCreateCV(depoint_cov, k, type = 'standard')
  mse1 = c(mse1,fnCrossEVa(flag_SCV)[[1]])
  mse2 = c(mse2, fnCrossEVa(flag_SCV)[[2]])
  mse3 = c(mse3, fnCrossEVa(flag_SCV)[[3]])
  mse4 = c(mse4, fnCrossEVa(flag_SCV)[[4]])
  mse11 = c(mse11,fnCrossEVa(flag_SCV)[[5]])
  mse21 = c(mse21, fnCrossEVa(flag_SCV)[[6]])
  mse31 = c(mse31, fnCrossEVa(flag_SCV)[[7]])
  mse41 = c(mse41, fnCrossEVa(flag_SCV)[[8]])
}

fnCrossEVa = function(flag_SCV){
  mse_meld = NULL
  mse_geo = NULL
  mse_area = NULL
  mse_down = NULL
  mse_meldcov = NULL
  mse_geocov = NULL
  mse_areacov = NULL
  mse_downcov = NULL
  for(i in 1:length(flag_SCV)){
   label = flag_SCV[[i]]
   train = depoint_cov[-label,]
   test = depoint_cov[label,]
    
  meld = fnPredictMelding(train, dearea_cov, dppoint = test, boundaryregion = boundaryregion)[[1]]
  geo = fnPredictMelding(train,dearea = NULL, dppoint = test, boundaryregion = boundaryregion)[[1]]
  area = fnPredictMelding(NULL,dearea_cov, dppoint = test, 
                              boundaryregion = boundaryregion)[[1]]
  down = fnPredictDown(train,dearea_cov, dppoint = test, boundaryregion = boundaryregion)[[1]]
  mse_meld = c(mse_meld, Metrics::mse(test$value, meld$pred_mean))
  mse_geo = c(mse_geo, Metrics::mse(test$value, geo$pred_mean))
  mse_area = c(mse_area, Metrics::mse(test$value, area$pred_mean))
  mse_down = c(mse_down, Metrics::mse(test$value, down$pred_mean))
  
  meldcov = fnPredictMeldingUK_1(train, dearea_cov, dppoint = test, boundaryregion = boundaryregion)[[1]]
  geocov = fnPredictMeldingUK_1(train,dearea = NULL, dppoint = test, boundaryregion = boundaryregion)[[1]]
  areacov = fnPredictMeldingUK_1(NULL,dearea_cov, dppoint = test, 
                          boundaryregion = boundaryregion)[[1]]
  downcov = fnPredictDownUK_3(train,dearea_cov, dppoint = test, boundaryregion = boundaryregion)[[1]]
  mse_meldcov = c(mse_meldcov, Metrics::mse(test$value, meldcov$pred_mean))
  mse_geocov = c(mse_geocov, Metrics::mse(test$value, geocov$pred_mean))
  mse_areacov = c(mse_areacov, Metrics::mse(test$value, areacov$pred_mean))
  mse_downcov = c(mse_downcov, Metrics::mse(test$value, downcov$pred_mean))
  }
  mse1 = mean(mse_meld)
  mse2 = mean(mse_geo)
  mse3 = mean(mse_area)
  mse4 = mean(mse_down)
  mse11 = mean(mse_meldcov)
  mse21 = mean(mse_geocov)
  mse31 = mean(mse_areacov)
  mse41 = mean(mse_downcov)
  return(list(mse1,mse2,mse3,mse4,mse11,mse21,mse31,mse41))
}


# No CV
#----
meld1 = fnPredictMelding(depoint, dearea, dppoint = depoint, boundaryregion = boundaryregion)
down1 = fnPredictDown(depoint,dearea, dppoint = depoint, boundaryregion = boundaryregion)
meld2 = fnPredictMeldingUK_1(depoint_cov, dearea_cov, dppoint = depoint_cov, boundaryregion = boundaryregion)
down2 = fnPredictDownUK_3 (depoint_cov,dearea_cov, dppoint = depoint_cov, boundaryregion = boundaryregion)
Metrics::mse(depoint$value, meld1[[1]]$pred_mean)
Metrics::mse(depoint_cov$value, meld2[[1]]$pred_mean)
Metrics::mse(depoint$value, down1[[1]]$pred_mean)
Metrics::mse(depoint_cov$value, down2[[1]]$pred_mean)

# Covariates Application 
#----

depoint <-  st_read("depoint_cov.shp")%>%
  st_transform(crsproj)

dearea <- st_read("dearea_cov.shp")%>%
  st_transform(crsproj)



# model selection
meld_1 = fnPredictMeldingUK_1(depoint_cov, dearea_cov, dppoint = dpconsurface, boundaryregion = boundaryregion)
meld_2 = fnPredictMeldingUK_2(de1, de2, dppoint = de1, boundaryregion = boundaryregion)[[3]]
meld_3 = fnPredictMeldingUK_3(de1, de2, dppoint = de1, boundaryregion = boundaryregion)[[3]]
meld_4 = fnPredictMeldingUK_4(de1, de2, dppoint = de1, boundaryregion = boundaryregion)[[3]]


for(j in 1: rep){
  flag_SCV = fnCreateCV(de1,n = 0.8, type = 'block')
  mse1 = c(mse1,fnCrossEVaCov(flag_SCV)[[1]])
  mse2 = c(mse2, fnCrossEVaCov(flag_SCV)[[2]])
  mse3 = c(mse3, fnCrossEVaCov(flag_SCV)[[3]])
}

fnCrossEVaCov = function(flag){
  mse_meld = NULL
  mse_geo = NULL
  mse_area = NULL
  mse_down = NULL
  for(i in 1:length(flag_SCV)){
    label = flag_SCV[[i]]
    train = depoint_cov[-label,]
    test = depoint_cov[label,]
  
    
    meld = fnPredictMeldingUK_1(train, dearea_cov, dppoint = test, boundaryregion = boundaryregion)[[1]]
    geo = fnPredictMeldingUK_1(train,dearea = NULL, dppoint = test, boundaryregion = boundaryregion)[[1]]
    area = fnPredictMeldingUK_1(NULL,dearea_cov, dppoint = test, 
                            boundaryregion = boundaryregion)[[1]]
    down = fnPredictDownUK_3(train, dearea_cov, dppoint = test, boundaryregion = boundaryregion)[[1]]
    mse_meld = c(mse_meld, Metrics::mse(test$value, meld$pred_mean) )
    mse_geo = c(mse_geo, Metrics::mse(test$value, geo$pred_mean))
    mse_area = c(mse_area, Metrics::mse(test$value, area$pred_mean))
    mse_down = c(mse_down, Metrics::mse(test$value, down$pred_mean))
  }
  mse1 = mean(mse_meld)
  mse2 = mean(mse_geo)
  mse3 = mean(mse_area)
  return(list(mse1,mse2,mse3))
}

mse1
mse2
mse3
mse4