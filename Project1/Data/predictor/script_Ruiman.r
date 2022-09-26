## Analysing Ruiman's data
library(maptools) # To deal with shp boundary
library(spatstat)
library(sf)

## Reading edge data and transforming it
uk_bd <- readShapeSpatial("area_bound/area_boud.shp")
proj4string(uk_bd) <- CRS("+init=epsg:4326")
SANC <- spTransform(uk_bd, "+proj=utm +zone=29 +nord +units=km")

# Building a field of view
W <- as.owin(SANC) # First a suitable window
plot(W)

#############
#Streets #LOAD FIRST the folder
setwd("road_coord")
temp <- list.files(pattern="*.csv")
myfiles <- lapply(temp, read.csv)

LinePP <- function(DF, win = W){
  coordinates(DF) <- c("X", "Y")
  proj4string(DF) <- CRS("+init=epsg:4326")
  projF <- spTransform(DF, "+proj=utm +zone=29 +nord +units=km")
  DG <- as.data.frame(projF)
  unique.ppp(ppp(x = DG$X, y = DG$Y, window = win))
}
ALLpp <- as.solist(lapply(myfiles, LinePP))
ALLpp <- superimpose(ALLpp)
ALLpp <- unique(ALLpp)
Allppt <- rthin(ALLpp, 0.04)
plot(Allppt, size = 0.02, cols = "black", pch = 20)

#################
# Dummy point pattern
Dummy <- rpoispp(4e8, win = W)
NNmarks <- nncross(Dummy, Allppt)$dist
#thresh <- quantile(NNmarks, 0.6) #This is to cut off the extreme values
#NNmarks[NNmarks > thresh] <- thresh
marks(Dummy) <- NNmarks
Cov <- Smooth.ppp(Dummy, diggle = T, sigma = 0.12e-4 , dimyx = 512)
plot(Cov, main = "Covariate")
#################
save(Cov, file = "CovariateRuiman.RData")
