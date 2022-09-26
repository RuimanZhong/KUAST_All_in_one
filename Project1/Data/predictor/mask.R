library(raster)
library(sf)
library(rgeoboundaries)
library(ggplot2)
library(viridis)
library(stars)
library(tidyverse)

theme_set(theme_minimal())
options(ggplot2.continuous.colour="viridis")
options(ggplot2.continuous.fill = "viridis")

crsproj <- 2158

boundaryregion <- geoboundaries('United Kingdom')
boundaryregion <- st_transform(boundaryregion, 4326)
setwd("~/Documents/project /TBC/project 1/spatialM /areal_data")
depoint <- read.csv("pointdata.csv")
depoint <- depoint[, c(2, 3, 4)]
colnames(depoint) <- c('value', 'y', 'x')

depoint <- depoint %>% st_as_sf(coords = c("x", "y"), dim = "XY") %>%
  st_set_crs(4326) %>%
  st_transform(2158)

setwd("~/Documents/project /TBC/project 1/raw_data/sdei-global-annual-gwr-pm2-5-modis-misr-seawifs-aod-v4-gl-03-2016-geotiff")
str <- 'p25_2016.tif'
glo_pm <- raster(str)
sr <- "+proj=utm +zone=29 +ellps=GRS80 +datum=NAD83 +units=m +no_defs" 

rr <- mask(crop(glo_pm, boundaryregion, snap = 'out'), boundaryregion)
rr <- projectRaster(rr, crs = sr)
# plot(rr)
fa <- 60
ra <- raster::aggregate(rr, fact = fa, fun = mean, na.rm = T)
plot(ra)

boundaryregion <- st_transform(boundaryregion, crsproj)
raa = mask(ra,boundaryregion)
plot(raa)

spol <- rasterToPolygons(ra, dissolve = F)
dearea <- st_as_sf(spol)
colnames(dearea) <- c('value','geometry')

dearea  <- st_transform(dearea, 2158)
depoint <- st_transform(depoint, 2158)
#r_df <- as.data.frame(ra, xy = TRUE, na.rm = TRUE)

# Covariates 
#----
setwd("~/Documents/project /TBC/project 1/raw_data/predictor")
dearea  <- st_transform(dearea, 4326)
depoint <- st_transform(depoint, 4326)

library(ncdf4)
library(mapview)

fnCreateRaster <- function(filename,savename,var){
  nc_data <- nc_open(filename)
  # Save the print(nc) dump to a text file
  {
    sink(savename)
    print(nc_data)
    sink()
  }
  lon <- ncvar_get(nc_data, "longitude")
  lat <- ncvar_get(nc_data, "latitude", verbose = F)
  t <- ncvar_get(nc_data, "time")
  value <- ncvar_get(nc_data, var) # store the data in a 3-dimensional array
  
  fillvalue <- ncatt_get(nc_data, var, "_FillValue")
  nc_close(nc_data)
  value[value == fillvalue$value] <- NA
  
  sum =  value[ , ,1]
  for(i in 2:12){
    sum <- sum + value[ , ,i]
  }
  annual_mean <- sum/12
  r <- raster(t(annual_mean), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), 
                        crs=CRS("EPSG:4326"))
  return(r)
}

# temp
dearea_cov <- st_read('dearea_cov.shp')
df <- read_ncdf("temp.nc" )
r_temp_2016 <- fnCreateRaster('temp.nc','tem_meata.txt','t2m')
area_bound <- st_as_sf(st_union(dearea_cov))

temp_uk <- mask(r_temp_2016,area_bound)
plot(temp_uk)
spol <- rasterToPolygons(temp_uk, dissolve = F)
temp_sf <- st_as_sf(spol)
colnames(temp_sf) <- c('temp','geometry')


df <- read_ncdf("precipation.nc" )
r_precipation_2016 <- fnCreateRaster('precipation.nc','pre_meata.txt','tp')
precipation_uk <- mask(r_precipation_2016,area_bound)
plot(precipation_uk)
spol <- rasterToPolygons(precipation_uk, dissolve = F)
precipation_sf <- st_as_sf(spol)
colnames(precipation_sf) <- c('precipation','geometry')

de_new <- dearea_cov
de_new$group <- 1:169
de_new <- de_new %>%
  st_join(precipation_sf)

de_new$area <- as.numeric(st_area(de_new))

de_new <- de_new%>%
  group_by(group,value)%>%
  summarise(area_sum = sum(area),
            precipation = weighted.mean(precipation,area/sum(area)))

dearea_cov<- de_new%>%
  st_join(temp_sf)%>%
  group_by(group,value)%>%
  summarise(area_sum = sum(area),
            temp_mean = weighted.mean(temp,area/sum(area)))

dearea_combine$temp_mean <- dearea_cov$temp_mean-273.15

dearea_cov <- dearea_combine[,c(2,4,5,6)]

# point

depoint_cov <- depoint_cov[,c(1,6)] %>%
  st_join(temp_sf) %>%
  st_join(precipation_sf)

depoint_cov$temp_mn <- depoint_cov$temp - 273.15
depoint_cov$road_length_1000 <- depoint_cov$r__1500
scaled_road<- scale(depoint_cov$road_length_1000)
depoint_cov$road_scale <- scaled_road
dearea_cov$road_scale <- scale(dearea_cov$road_length_1000)
# road length 
# area
library(rlang)
library(osmdata)

fnHighwayL <- function(box){
  highways <- opq(bbox = box) %>%
    add_osm_feature(key = "highway", value = c("motorway", "trunk", "primary", "secondary")) %>% osmdata_sf()
  if(!is.null(highways$osm_lines)){
    road <- sum(st_length(highways$osm_lines))
  } else{
    road <- 0
  }
  
  return(road)
}

road_length <- NULL
for (i in length(road_length)+1 : nrow(dearea_combine)){
  print(i)
  box <- unname(st_bbox(dearea_combine[i,]))
  road <- fnHighwayL(box)
  road_length <- c(road_length, road)
}

dearea_cov$road_length_1000 <- road_length

# buffered 
# Here use UTM and meters
buffer_margin_1000 = units::set_units(1000, m)
buffered_point_utm <- st_buffer(depoint_cov, buffer_margin_1000)
buffered_point <- st_transform(buffered_point_utm,4326)
road_length <- NULL
for (i in length(road_length)+1 : nrow(buffered_point)){
  print(i)
  box <- unname(st_bbox(buffered_point[i,]))
  road <- fnHighwayL(box)
  road_length <- c(road_length, road)
}

depoint_cov$road_length_1000 <- road_length

# Areal road proportion
#----
fnCreateAreaProp <- function(radius, route, area){
  
  utm_route = st_transform(route, 2158)
  buffer_margin = units::set_units(radius, m)
  road_prop <- utm_route%>%
    st_buffer(buffer_margin) %>% 
    st_transform(4326)%>%
    st_union() %>%
    st_area()/st_area(area)
  
  return(road_prop)
}
road_prop <- NULL
dearea_combine <-dearea_cov

for (i in length(road_prop)+1 : nrow(dearea_combine)){
  print(i)
  box <- unname(st_bbox(dearea_combine[i,]))
  
  highways <- opq(bbox = box) %>%
    add_osm_feature(key = "highway", value = c("motorway", "trunk", "primary", "secondary")) %>% 
    osmdata_sf()
  
  route <- highways$osm_lines
 if(!is.null(route)){prop <- fnCreateAreaProp(1000,route,dearea_combine[i,])} else {prop <- 0}

  road_prop<- c(road_prop, prop)
}

dearea_cov$road_prop <- road_prop
# point within buffered route
#----
fnFindRoadIn <-function(radius,box){
  highways <- opq(bbox = box) %>%
    add_osm_feature(key = "highway", value = c("motorway", "trunk", "primary", "secondary")) %>% 
    osmdata_sf()
  buffer_margin = units::set_units(radius, m)
  if(!is.null(highways$osm_lines)){
    road_buffered <- highways$osm_lines%>%
      st_transform(4326) %>%
      st_buffer(buffer_margin) %>% 
      st_transform(4326)
    
    road_in_point <- as.numeric(nrow(st_intersection(depoint_cov, road_buffered)))/as.numeric(nrow(road_buffered))
    
  }else {road_in_point <- 0}
  
  return(road_in_point)
}
radius <- 1000
road_in <- NULL
# find each point within which area, 
# since we can not download the road information of whole UK, 
# so I do it area by area
point_in_area <- st_within(depoint_cov,dearea_cov)
for(i in length(road_in)+1 : nrow(point_in_area)){
  print(i)
  box <- unname(st_bbox(dearea_cov[point_in_area[[i]],]))
  road_in_point <- fnFindRoadIn(radius,box)
  road_in <- c(road_in,road_in_point)
}

depoint_cov$road_in_prop <- road_in
#----
st_write(de1, "depoint_cov.shp")
st_write(dearea_cov, "dearea_cov.shp")

de1 <-  st_read("depoint_cov.shp")%>%
  st_transform(crsproj)

de2 <- st_read("dearea_cov.shp")%>%
  st_transform(crsproj)

de1<- de1[,-3]
#----

highways <- opq(bbox = box) %>%
  add_osm_feature(key = "highway", value = c("motorway", "trunk", "primary", "secondary")) %>% 
  osmdata_sf()
buffer_margin = units::set_units(radius, m)
if(!is.null(highways$osm_lines)){
  road_buffered <- highways$osm_lines%>%
    st_transform(4326) %>%
    st_buffer(buffer_margin) %>% 
    st_transform(4326)}
  


# data description 
#----
gplot_data <- function(x, maxpixels = 50000)  {
  x <- raster::sampleRegular(x, maxpixels, asRaster = TRUE)
  coords <- raster::xyFromCell(x, seq_len(raster::ncell(x)))
  ## Extract values
  dat <- utils::stack(as.data.frame(raster::getValues(x))) 
  names(dat) <- c('value', 'variable')
  
  dat <- tibble::as_tibble(data.frame(coords, dat))
  
  if (!is.null(levels(x))) {
    dat <- dplyr::left_join(dat, levels(x)[[1]], 
                            by = c("value" = "ID"))
  }
  dat
}

rr_t <- gplot_data(rr)
raa_t <- gplot_data(raa)

 p1 = ggplot()+
  geom_tile(data = dplyr::filter(rr_t, !is.na(value)), 
            aes(x = x, y = y,fill = value)) +
  geom_tile(data = dplyr::filter(raa_t, !is.na(value)), 
            aes(x = x, y = y, fill = value))  +
   geom_sf(data = boundaryregion, inherit.aes = FALSE, fill = NA)+
   labs(
     title = 'Areal Data',
     y = 'Longitude',
     x = 'Latitude',
     fill = "PM 2.5"
   )+theme_bw()+theme_minimal()
 
p2 = ggplot(data = boundaryregion) + geom_sf() +
  geom_sf(data = depoint, aes(col = value))+
  labs(
    title = 'Point Data',
    y = 'Longitude',
    x = 'Latitude',
    color = "PM 2.5"
  )+theme_bw() + theme_minimal()

require(gridExtra)
plot1 <- p1
plot2 <- p2
grid.arrange(plot1, plot2, ncol=2)

setwd("~/Documents/project /TBC/project 1/results")
png('raw_data.png',type = 'cairo', res = 300, width = 1600, height = 960)
grid.arrange(plot1, plot2, ncol=2)
dev.off()

# Prediction Mask
result <- results[,c(6,7,8,9)]
example_points <- as(result, "Spatial")

# Generate empty raster layer and rasterize points
raster <- raster(resolution = 15000, crs = sr, 
                    xmn = 500000, xmx = 1500000, ymn = 5000000, ymx = 7000000) %>% 
  rasterize(example_points, ., fun = "last", background = 0)

plot(raster)

# Covariates raster

# Plot
#----
crsproj <- 2158
boundaryregion <- geoboundaries('United Kingdom')
boundaryregion <- st_transform(boundaryregion, crsproj)
sr <- '+proj=utm +zone=29 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m
+no_defs'

plot(r_new)

class(Cov)
road_raster <- raster(Cov)
crs(road_raster) <- "+proj=utm +zone=29 +nord +units=km"
road_raster <- projectRaster(road_raster, crs = sr)
plot(road_raster)

r_new <- projectRaster(precipation_uk, crs = sr)
r_df <- as.data.frame(r_new, xy = TRUE, na.rm = TRUE)
r_df$layer <- as.vector(scale(r_df$layer,center = T))
p1 <- ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = layer)) +
  geom_sf(data = boundaryregion, inherit.aes = FALSE,colour = 'white',fill = NA) +
  scale_fill_viridis() +
  labs(
    title = "Total Precipitation",
    fill = "Precipitation",
    subtitle = "units:m"
  )+ theme(axis.title.x=element_blank(),
           axis.ticks.x=element_blank(),
           axis.title.y=element_blank(),
           axis.ticks.y=element_blank())

r_new <- projectRaster(temp_uk, crs = sr)
r_df <- as.data.frame(r_new, xy = TRUE, na.rm = TRUE)
r_df$layer <- as.vector(scale(r_df$layer,center = T))
p2 <- ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = layer-273.15)) +
  geom_sf(data = boundaryregion, inherit.aes = FALSE,colour = 'white',fill = NA) +
  scale_fill_viridis() +
  labs(
    title = "Temperature",
    fill = "Temperature",
    subtitle = "units:Celsius "
  ) + theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.y=element_blank())

r_df <- as.data.frame(road_raster, xy = TRUE, na.rm = TRUE)
r_df$layer <- as.vector(log(r_df$layer))

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(r_df$layer, probs=c(.25, .95), na.rm = na.rm)
  y <- r_df[r_df$layer < 2.2, ]
  y
}
r_df <- remove_outliers(r_df)

p3 <- ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = layer)) +
  geom_sf(data = boundaryregion, inherit.aes = FALSE,colour = 'white',fill = NA) +
  scale_fill_viridis() +
  labs(
    title = "Road Distance",
    fill = "Distance",
    subtitle = "units:log(Kilometers) "
  ) + theme(axis.title.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.ticks.y=element_blank())

p3


library(ggpubr)
ggarrange(p1, p2,p3, ncol = 3, nrow = 1)
