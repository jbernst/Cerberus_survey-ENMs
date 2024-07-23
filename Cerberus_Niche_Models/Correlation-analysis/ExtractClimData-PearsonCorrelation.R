# load libraries
library(raster) # for raster management
library(rgdal)
library(ecospat) # for niche statistics and analysis
library(maptools) # for mapping
library(tidyverse) # for data wranging aand mapping functions
library(sf)
library(heatmaply)

rm(list=ls())

setwd

setwd("D:/Documents/Publications/Cerberus-Autecology/Correlation-analysis/")

##Load your spatial data

#Read in coordinate data for all species
sp<- read.csv(file = "Cerberus_records.csv", header = T)
View(sp)
colnames(sp) <- c("lon", "lat")
sp.xy<- sp[c("lon", "lat")] # this is the xy dataset

## create species data.frame with extracted climate info 
# import clim data
raster_files <- list.files("D:/Documents/R/wc2.1_2.5m_bio/", full.names = T, pattern = ".tif") #wc2
#raster_files <- list.files("D:/Documents/R/wc2-5/", full.names = T, pattern = ".bil") #wc1

# create a raster stack using the list you created
predictors <- stack(raster_files)
# verify full clim data
#plot(predictors$bio_1)

## create geographic extent around the species data
# oder is (xmn, xmx, ymn, ymx), use the data point summary() to get the max/min for lon/lat
  # 10 degrees = ~1110`km

geo.ext.sqbuff <- extent(min(sp$lon)-5, max(sp$lon)+5, min(sp$lat)-5, max(sp$lat)+3)



# crop the predictors to a narrower geo.extent
predictors <- crop(x = predictors, y = geo.ext.sqbuff)
# verify
plot(predictors$bio_1)
points(sp.xy[,c(1,2)], col = "blue", pch = 21)

# extract predictor (raster) data into species points
ext.sp <- raster::extract(predictors, sp.xy)


# cbind() the species point data with the raster values
clim.sp <- cbind(sp.xy, ext.sp)
str(clim.sp)

# add the occurrence variable to the species file
clim.sp$species_occ <- 1
str(clim.sp)

# METHOD 1: BUFFERS - create spatial buffers around the points for each species
## based on a point buffer of 5  degrees #### This is preferred ####

# convert the species points into sf object

#sp.sf <- st_as_sf(sp.xy[, c(1:2)], coords = c("lon","lat"), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

sp.sf <- st_as_sf(sp.xy[, c(1:2)], coords = c("lon","lat"))

summary(sp.sf) # examine 
class(sp.sf)  # examine 

#crs(predictors) <- raster::crs(sp.sf) # match the predictor & species CRS.

# create a buffer on the species points with a 5 degree distance, then unite all buffer circles and convert to sf
sp.buf <- sf::st_buffer(sp.sf, dist = 5) %>% sf::st_union() %>% sf::st_sf()
# NOTE: it will show a message that it is not projected, but that is OK
#  verify the full predictors
plot(predictors[[1]], main = names(predictors)[1])
# add the species points 
points(sp.xy[, c(1:2)])
# add the buffer 
# NOTE: use add = TRUE to include in the current plot
plot(sp.buf, border = "blue", lwd = 3, add = TRUE)
# crop the predictors based on the buffer
predictors.an <- crop(predictors, sp.buf)
# create a mask to remove the area outside the predictors
predictors.an <- raster::mask(predictors.an, sp.buf)
# verify predictors were properly cropped
plot(predictors.an$bio_1)

# create XX thousand random background points
  # Randomly sample 10,000 background points from one background extent raster 
  # (only one per cell without replacement). NOTE: Since the raster has <10,000 pixels, 
  # you'll get a warning and all pixels will be used for background.
bg.an <- dismo::randomPoints(predictors.an[[1]], n = 10000) %>% as.data.frame()
colnames(bg.an) <- colnames(sp.xy)

# Notice how we have pretty good coverage (every cell).
plot(predictors.an[[1]])
points(bg.an, pch = 20, cex = 0.2)

# extract raster stack values into species points
  # point file must have only two variables x & y
  # extract raster data to background points
  # if you  have tidyverse on, then use the function raster::extract() because extract() will conflict

ext.bg <- raster::extract(predictors.an, bg.an)
  # this produces a large matrix of only the raster values for each point
  # cbind() the species point data with the raster values
ext.bg.an <- cbind(bg.an, ext.bg)
str(ext.bg.an)
# add the occurrence variable to the species file and verify
ext.bg.an$species_occ <- 0
## join the presence vs background data.frames
Cerberus.f <- rbind(clim.sp, ext.bg.an)

# verify all data is clean
summary(Cerberus.f)
  # remove NAs
# create a new data.frame with complete records
# tow ways to do it: base R
Cerberus.f <- subset(Cerberus.f, !is.na(bio_1) & !is.na(bio_10))
summary(Cerberus.f)

# save data.frame
write.csv(Cerberus.f, file = "./Cerberus_final.csv", row.names = F)

#### Pearson Correlation ## by Angelo Soto-Centeno

# read file
clim <- read.csv(file = "Cerberus_final.csv", header = T) #species width data file (swd)

# inspect data
head(clim)
# set up correlation in data.rame 
# the characters in brackets [rows,column] indicate portion of data.frame to do correlation on
corclim <- cor(clim[,3:21])  
# run Pearson Correlation
corclim
# shows statistics for variables
summary(lm(bio_18 ~ bio_11, data =clim))
# show the correlation plots 
plot(clim[, 3:21]) 
# save correlation matrix
write.csv(corclim, file = "Cerberus_corClim1.csv", row.names = F)

spclim <- read.csv("Cerberus_final.csv", header = T)
head(spclim)

#Instead of the above, can also use heatmaply to make an interactive heatmap of correlated variables
heatmaply(x = cor(spclim[,3:21]),
          xlab = "Evironment", ylab = "Environment", k_col = 2, k_row = 2)
View(spclim)
#Check which are less than or equal to -0.8 or greater than or equal to 0.8
#run the heatmaply function again only on the uncorrelated ones and you can get 
#another heatmap one for ones that are uncorrelated
View(spclim)
heatmaply(x = cor(spclim[,c(3,14,16,19,6,9,10,12,13)]),
          xlab = "Evironment", ylab = "Environment", k_col = 2, k_row = 2)

#Another way to visualize the heatmap is below
cor.test.p <- function(x){
  FUN <- function(x, y) cor.test(x, y)[["p.value"]]
  z <- outer(
    colnames(x), 
    colnames(x), 
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}

r <- cor(spclim[,3:21])

p <- cor.test.p(spclim[,3:21])

heatmaply_cor(
  r,
  node_type = "scatter",
  point_size_mat = -log10(p), 
  point_size_name = "-log10(p-value)",
  label_names = c("x", "y", "Correlation")
)
