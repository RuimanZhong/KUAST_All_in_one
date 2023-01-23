#-----------------------------
# tiff data project 2
#-----------------------------
library(raster)
library(sf)
library(sp)

str <- 'gbr_ppp_2020_UNadj_constrained.tiff'
r <- raster(str)
sr <- "CRS('+init=EPSG:4326')" 
raster::crs(r) <- "EPSG:4326"
plot(r)
spol <- rasterToPolygons(r, dissolve = F)
population <- st_as_sf(spol)
colnames(dearea) <- c('value','geometry')