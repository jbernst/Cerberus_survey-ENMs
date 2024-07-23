### - Cerberus ENM - ###
### - script by: Justin Matthew Bernstein - ###
### - date: 3 February 2024 - ###
# Note: Due to size, bioclim layers are not hosted on Github.
### --- clear up your Global Environment & Devises --- ###
rm(list = ls())
dev.off()

# load libraries
library(raster) # for raster management
library(rasterVis) # for visualizing raster data
library(rgdal) # for raster management
library(sf) # for raster management
library(sp) # for shapefile management
library(spThin) # for spatial data processing
library(ecospat) # for niche statistics and decomposition
library(ENMTools) # for niche statistic and decomposition
library(humboldt) # for niche statistics and decomposition
library(tcltk) # to use the tcltk progress bar in quartz
library(maptools) # for mapping
library(tidyverse) # for data wrangling and mapping functions
library(dismo) # for niche modeling
library(rJava)
options(java.parameters = "-Xmx32g" ) # turn this on to give the maximum memory# for niche modeling
library(ENMeval) # for niche modeling parameter evaluation
library(maxnet) # for niche modeling parameter evaluation
library(rmaxent) # to calculate MESS plots
library(gridExtra) # for combining plots
library(rworldmap)

getwd()
setwd("D:/Documents/Publications/Cerberus-Autecology/Cerberus_ENM/")

## spatial data processing
## part 1. species observations wrangling
# load species data
cerb <- read.csv(file = "../Correlation-analysis/Cerberus_records.csv", header = T)

# make data backup
data -> data.master

# remove duplicate records
#data.dup = duplicated(data[, c("x", "y")])
#unique(data.dup) # verify if there are any dups
# data.clean <- data[!data.dup, ] # data is already clean for dups

# verify the species data 
# get map data
# base R map 
library(maptools)
# attach map daata
data("wrld_simpl")

# plot the Philippines
shp <- readOGR("D:/Documents/Publications/Cerberus-Autecology/PHL_adm/PHL_adm0.shp")
plot(shp)


# add species observation points
points(x = cerb$lon, y = cerb$lat, col = "#46804E", pch = 20, cex = 1.5)

###############################################
## Import climate data
###############################################

# import clim data + elevation
raster_preswc <- list.files("/Documents/R/wc2.1_2.5m_bio/", full.names = T, pattern = ".tif") ## wc2
# print/confirm list
raster_preswc

# if need to remove, for example, elevation (#20 on the list)
raster_preswc <- raster_preswc[-20] ## wc2
raster_preswc <- raster_preswc[c(1,12,14,17,4,7,8,10,11)]

# create a raster stack using the list you created
predictors <- stack(raster_preswc) ## wc2
# verify rasters
# plot(predictors[[1]], main = names(predictors)[1])
plot(predictors$bio_2, main = "Bio 2")
# from predictors, choose the first one. then take the name of predictor 1


# crop environmental layers
# 1. create a projection extent around the species data
# oder is (xmn, xmx, ymn, ymx), use the data point summary() to get the max/min for lon/lat
# 10 degrees = ~1110`km
geo.ext.sqbuff <- extent(shp)
# crop the predictors to a narrower geo.extent
env <- crop(x = predictors, y = geo.ext.sqbuff)
## or use predictors.pr object name to match
# verify
plot(env$bio_2, col = viridis::turbo(99, alpha = 0.85, direction = -1), main = "Bio2 Projection Extent", xlab = "longitude", ylab = "latitude")

# add species observation points
points(x = cerb$lon, y = cerb$lat, col = "blue", pch = 20, cex = 0.75)
# 
# ## LGM projection data
# raster_pastwc <- list.files("/Volumes/ANGELO4/ASC_GIS/Layers/Climate/Climate/Past/wc_1_CCSM_LGM_2-5m/", full.names = T, pattern = ".tif") ## wc1
# # confirm list
# raster_pastwc
# # if need to remove, for example, elevation (#20 on the list)
# raster_pastwc <- raster_pastwc[c(1,12,14,2,4,7,8,10)]
# 
# # create a raster stack using the list you created
# p.predictors <- stack(raster_pastwc) ## wc1
# # verify rasters
# # plot(predictors[[1]], main = names(predictors)[1])
# plot(p.predictors$bio_2, main = "Bio 2")
# 
# # crop environmental layers of LGM
# # 1. create a projection extent around the species data
# # oder is (xmn, xmx, ymn, ymx), use the data point summary() to get the max/min for lon/lat
# # 10 degrees = ~1110`km
# geo.ext.sqbuff <- extent(90, 142, -13, 29)
# # crop the predictors to a narrower geo.extent
# env1 <- crop(x = p.predictors, y = geo.ext.sqbuff)
# ## or use predictors.pr object name to match
# # verify
# plot(env1$bio_2, col = viridis::turbo(99, alpha = 0.85, direction = -1), main = "Bio2 Projection Extent LGM", xlab = "longitude", ylab = "latitude")
# 

# 2. crate calibration extents
# WITH A POINT BUFFER -- use package sf

# Make our occs into a sf object -- as the coordinate reference system (crs) for these 
# points is WGS84, a geographic crs (lat/lon) and the same as our envs rasters, we specify it 
# as the RasterStack's crs.
occs.cerb.sf <- sf::st_as_sf(cerb, coords = c("lon", "lat"), crs = raster::crs(env)) ## point here to the original data.frame for species

# Now, we project our point data to an equal-area projection, which converts our 
# degrees to meters, which is ideal for buffering (the next step). 
# We use the typical Eckert IV projection.
eckertIV <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
occs.cerb.sf <- sf::st_transform(occs.cerb.sf, crs = eckertIV)

# Buffer all occurrences by 500 km, union the polygons together 
# (for visualization), and convert back to a form that the raster package 
# can use. Finally, we reproject the buffers back to WGS84 (lat/lon).
# We choose 500 km here to avoid sampling the Caribbean islands.
occs.buf <- sf::st_buffer(occs.cerb.sf, dist = 500000) %>% 
  sf::st_union() %>% 
  sf::st_sf() %>%
  sf::st_transform(crs = raster::crs(env))
plot(env$bio_2, main = names(env)[1])
points(cerb)

# To add sf objects to a plot, use add = TRUE
plot(occs.buf, border = "blue", lwd = 3, add = TRUE)

# Crop environmental rasters to match the study extent
envs.bg.cerb <- raster::crop(env, occs.buf)
# Next, mask the rasters to the shape of the buffers
envs.bg.cerb <- raster::mask(envs.bg.cerb, occs.buf)
plot(envs.bg.cerb[[1]], main = names(env)[1])
points(cerb)
plot(occs.buf, border = "blue", lwd = 3, add = TRUE)

# Randomly sample 10,000 background points from one background extent raster 
# (only one per cell without replacement). Note: Since the raster has <10,000 pixels, 
# you'll get a warning and all pixels will be used for background. We will be sampling 
# from the biome variable because it is missing some grid cells, and we are trying to 
# avoid getting background points with NA. If one raster in the stack has NAs where the
# other rasters have data, ENMeval internally converts these cells to NA.
bg.cerb <- dismo::randomPoints(envs.bg.cerb[[7]], n = 10000) %>% as.data.frame()
colnames(bg.cerb) <- colnames(cerb)

# Notice how we have pretty good coverage (every cell).
plot(envs.bg.cerb[[1]])
points(bg.cerb, pch = 20, cex = 0.2)


#################################################################
##  ENMeval model testing                                      ##
#################################################################

####### compare a range of models using ENMeval


# set the parameter tunings
tune.args <- list(fc = c("L", "LQ", "H", "LQH", "LQHP"), rm = 1:5)

eeval.cerb <- ENMevaluate(occs = cerb, envs = env, bg = bg.cerb, 
                          algorithm = 'maxnet', partitions = 'block', 
                          tune.args = tune.args)
# save the results
write.csv(eeval.cerb@results, file = "cerb.eval.res1.csv", row.names = F)
# best tuning
# fc = LQ, rm = 1



############################################################
## Run Maxent model using custom parameters
############################################################


list.files()
#dir.create("cerb_mxnt")
setwd("cerb_mxnt/")

# set maxent parameters for MATA from ENMeval fc = LQHP, rm = 2
args <- c("linear=true", "quadratic=true", 
          "product=false", "hinge=false", "threshold=false", 
          "betamultiplier=1")

# make single run model using ALL species points
cerb.cur.dist <- maxent(envs.bg.cerb, cerb, args = args, path = "./")

# NOTES: if response curves & jacknife plots are needed, add args = c(-J", "-P", ...)
# all args must be added without spaces
# threshold features are ommited by default in maxent to produce smoother/simpler models that are likely more realistic *use threshold with caution

# check variable contribution
plot(cerb.cur.dist, main = "C schneiderii present")
# check response curves
response(cerb.cur.dist, main = "C schneiderii present")


# make present model prediction
cerb.cur.mod <- predict(cerb.cur.dist, env, args = c("outputformat=cloglog"), progress = "text")
# NOTES: may change outputformat options "=raw" or "=logistic" or "=cloglog"

# calculate Boyce index
ecospat.boyce(cerb.cur.mod, cerb, window.w = "default", res = 100, PEplot = T)

# performance statistics
# AUC = 0.832; BI = 0.931

# plot original model
plot(raster, col = viridis::turbo(n = 99, alpha = 0.8, direction = 1), main = "C schneiderii (fc = LQ, rm = 1 :: AUC = 0.832; BI = 0.931)", xlab = "longitude", ylab = "latitude")

View(cerb.cur.mod)

# save raster for plotting
list.files()
#dir.create("cerb_plots")
setwd("cerb_plots/")
writeRaster(cerb.cur.mod, filename = "cerbENM_present.tif")

# extract suitability at occurence points
setwd("D:/Documents/Publications/Cerberus-Autecology/Cerberus_ENM/cerb_mxnt/cerb_plots")
raster <- raster("cerbENM_present.tif")
plot(raster)
occurences <- cerb
occurrence_points <- SpatialPoints(cerb[, c("lon", "lat")])
extracted_values <- extract(raster, occurrence_points)
View(occurrence_points)

extracted_values2 <- cbind(cerb[, c("lon", "lat")], new_column = extracted_values)
View(extracted_values2)
write.csv(extracted_values2, file = "cerb_coordinates_suitability.csv")

#NEEDS WORK
##### linear regressions

### PPT correlations
occurences.water <- read.csv("water-averages-coordinates.csv")
occurrence_points.water <- SpatialPoints(occurences.water[, c("lon", "lat")])
extracted_values <- extract(raster, occurrence_points.water)
water <- read.csv("water-averages.csv")

# check for normality
shapiro.test(water$PPT) # p = 0.0001 = not normally distributed

# transform them
ppt.log <- log(water$PPT)

# no significant different found
lm.ppt <- lm(extracted_values~water$PPT)
summary(lm.ppt)

plot(extracted_values~water$PPT)
plot(extracted_values~ppt.log)

### mangroves correlations

# import clim data + elevation
raster_mangroves <- list.files("/Documents/Publications/Cerberus-Autecology/Cerberus_ENM/GMW-12_1996-2020_v3.0/", full.names = T, pattern = ".tif") ## wc2


# import shape file of mangroves
setwd("~/Publications/Cerberus-Autecology/Cerberus_ENM/")
mangrove.all <- readOGR("Philippines_GMW_v3_2020.shp")
plot(mangrove.all)
points(x = cerb$lon, y = cerb$lat, col = "#46804E", pch = 20, cex = 1.5)
cerb

# convert cerb points to spatial points
sp_points <- st_as_sf(cerb, coords = c('lon',"lat"))#make points spatial
st_crs(sp_points) = 4326 # Give the points a coordinate reference system (CRS); can find this in properties of the layer in QGIS

sp_points <- st_transform(sp_points, crs = st_crs(countries)) # Match the point and polygon CRS if they aren't anymore (st_transform matches the projections)

polyinter <- st_intersection(points_shp, mangrove.all)


#----TEST-----


rst_phil_mgrv_high_res <- rasterize(mangrove.all, rst_template_high_res)
plot(rst_phil_mgrv_high_res, col = "red", legend = FALSE, xlab = "Longitude", ylab = "Latitude")
plot(rasterToPolygons(rst_phil_mgrv_high_res), add = TRUE, col = "red")
plot(mangrove.all, add = TRUE)

## rasterize

rst_template <- raster(ncols = 4500, nrows = 4500, 
                        crs = projection(mangrove.all), 
                        ext = extent(shp))

rst_phil_mgrv <- rasterize(mangrove.all, rst_template)
plot(rst_phil_mgrv, col = "red", legend = FALSE, xlab = "lon", ylab = "lat")
plot(shp)
plot(rst_phil_mgrv, col = "red", legend = FALSE, xlab = "lon", ylab = "lat", add = TRUE)

# change back to a shape file
plot(rasterToPolygons(rst_phil_mgrv), add = TRUE, col = "blue")

## raster template2
rst_template2 <- raster(ncols = 300, nrows = 400, 
                       crs = projection(mangrove.all), 
                       ext = extent(shp))

## rasterize
rst_phil_mgrv2 <- rasterize(mangrove.all, rst_template2)
plot(rst_phil_mgrv2, col = "red", legend = FALSE, xlab = "lon", ylab = "lat")
writeRaster(rst_phil_mgrv2, filename = "mangroves_2020.tif")


plot(shp)
plot(rst_phil_mgrv2, col = "red", legend = FALSE, xlab = "lon", ylab = "lat", add = TRUE)
plot(rasterToPolygons(rst_phil_mgrv2), add = TRUE, col = "red")



## raster template3
rst_template <- raster(ncols = 100, nrows = 400, 
                       crs = projection(mangrove.all), 
                       ext = extent(shp))
extent_shp <- extent(mangrove.all)

ncols <- 100
nrows <- 400

cellsize_x <- (extent_shp@xmax - extent_shp@xmin)/ncols
cellsize_y <- (extent_shp@ymax - extent_shp@ymin)/nrows

rst_template_high_res <- raster(ncols = ncols, nrows = nrows, 
                                xmn = extent_shp@xmin, xmx = extent_shp@xmax,
                                ymn = extent_shp@ymin, ymx = extent_shp@ymax,
                                crs = projection(mangrove.all), 
                                res = c(cellsize_x, cellsize_y))
#----TEST-----






# if mangrove data is a raster, can crop it to the shape file and then mask
# obj1 <- crop(raster, shape)
# obj2 <- mask(ohb1, shape)

# if need to remove, for example, elevation (#20 on the list)
raster_preswc <- raster_preswc[-20] ## wc2
raster_preswc <- raster_preswc[c(1,12,14,17,4,7,8,10,11)]

# create a raster stack using the list you created
predictors <- stack(raster_preswc) ## wc2
# verify rasters
# plot(predictors[[1]], main = names(predic





# create a raster stack using the list you created
occurences <- cerb
occurrence_points <- SpatialPoints(cerb[, c("lon", "lat")])

result <- data.frame(lon = occurrence_points$lon, lat = occurrence_points$lat)

# Loop through each raster in 'mangroves' and extract values
for (raster_file in raster_mangroves) {
  # Load the raster
  current_raster <- raster(raster_file)
  
  # Resample if needed (optional)
  # current_raster <- resample(current_raster, target_raster)
  
  # Extract values
  values <- extract(current_raster, occurrence_points)
  
  # Create a new column in 'result' for the current raster
  result[paste0("Value_", basename(raster_file))] <- values
}

View(result)

extracted_values <- extract(mangroves, occurrence_points)




