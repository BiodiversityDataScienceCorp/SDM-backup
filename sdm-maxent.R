
#Mila's code, for reference

# dependencies

install.packages("dismo")
install.packages("maptools")
install.packages("tidyverse")
install.packages("maxnet")
install.packages("rJava")

library(dismo)
library(maptools)
library(tidyverse)
library(maxnet)
library(rJava)

# world map data
data(wrld_simpl)
#n.america <- wrld_simpl %>% filter(ISO3 == "USA" | ISO3 == "MEX"| ISO3 == "CAN")

# read occurence data
ranaData <- read_csv("ranaData.csv")
ranaData <- ranaData %>% dplyr::select(longitude, latitude)
ranaDataNotCoords <- read_csv("ranaData.csv") %>% dplyr::select(longitude, latitude)

# convert one data set (ranaData) to spatial
ranaDataSpatialPts <- SpatialPoints(ranaData, proj4string = CRS("+proj=longlat"))

# climate data: use get data only once
# getData("worldclim", var="bio", res=2.5) #not the correct name I think for the variable

clim_list <- list.files(path = "wc2-5/", pattern = ".bil$", 
                        full.names = T)  # '..' leads to the path above the folder where the .rmd file is located

# stacking the bioclim variables to process them at one go
clim <- raster::stack(clim_list)

plot(clim[[1]]) # show that it is the first env layer ( = ?)
plot(ranaDataSpatialPts, add = TRUE) 


# determine geographic extent of our data
geographic.extent <- extent(x = ranaDataSpatialPts)


# Create pseudo-absence points (making them up, using 'background' approach)
# first we need a raster layer to make the points up on, just picking 1
bil.files <- list.files(path = "wc2-5", 
                        pattern = "*.bil$", 
                        full.names = TRUE)
mask <- raster(bil.files[1])

# Random points for background (same number as our observed points)
set.seed(19470909) # seed set so we get the same background points each time we run this code 

# raster package will complain about not having coordinate reference system,
# so we suppress that warning
# extf makes the etent of the random points slightly larger than the given ext
background.points <- randomPoints(mask = mask, 
                                  n = .5 * nrow(ranaDataNotCoords),
                                  ext = geographic.extent, 
                                  extf = 1.25,
                                  warn = 0)
# add col names (can click and see right now they are x and y)
colnames(background.points) <- c("longitude", "latitude")

# separate presence data into training model and testing model
# randomly select 50% of data for training
set.seed(98)
selected <- sample(1:nrow(ranaDataNotCoords), (0.5 * nrow(ranaDataNotCoords)))
occ_train <- ranaDataNotCoords[selected, ]  # this is the selection to be used for model training
occ_test <- ranaDataNotCoords[-selected, ]  # this is the opposite of the selection which will be used for model testing

# Data for observation sites (presence and background), with climate data
# run each separately, wait until finished at 2.5 gb RAM
occ.train.values <- raster::extract(x = clim, y = occ_train) # why are there NA values ? all of the ranaD. has values
occ.test.values <- raster::extract(x = clim, y = occ_test) # why are there NA values ? all of the ranaD. has values
absence.values <- raster::extract(x = clim, y = background.points)

occ.train.values2 <- na.omit(occ.train.values)
occ.test.values2 <- na.omit(occ.test.values)
absence.values2 <- na.omit(absence.values)

testing.Nas <- cbind(occ_train, occ.train.values)
testing.Nas2 <- cbind(occ_test, occ.test.values)



# Create data frame with presence training data and backround points (0 = abs, 1 = pres)
presence.absence.vector <- c(rep(1, nrow(occ.train.values2)), rep(0, nrow(absence.values)))
presence.absence.train.env.data <- as.data.frame(rbind(occ.train.values2, absence.values)) 
# Mila: dimensions of these objects good bc = 2n(.5 * ranaData?)

## with dismo (java)
ranaModelDismo <- dismo::maxent(x = presence.absence.train.env.data, ## env conditions
                                p = presence.absence.vector,   ## 1:presence or 0:absence
                                path=paste0("maxent_outputs"), ## folder for maxent output; 
                                # if we do not specify a folder R will put the results in a temp file, 
                                # and it gets messy to read those. . .
                                args=c("responsecurves") ## parameter specification
)

# view the maxent model 
ranaModelDismo
plot(ranaModelDismo)
response(ranaModelDismo)

predictExtent <- 1.25 * geographic.extent # choose here what is reasonable for your pts


# graph it
geographicArea <- crop(clim, predictExtent)

ranaPredictPlot <- raster::predict(ranaModelDismo, geographicArea) 

plot(ranaPredictPlot)



## with maxnet
# train Maxent with tabular data
ranaModel <- maxnet::maxnet(data = presence.absence.train.env.data, ## env conditions
                            p = presence.absence.vector)
## parameter specification
#)

ranaModel
plot(ranaModel)
#response(ranaModel)
#maxnet::response.plot(ranaModel)

# the maxent functions runs a model in the default settings. To change these parameters,
# you have to tell it what you want...i.e. response curves or the type of features
















#light evaluation
# Evaluation indices include AUC, TSS, Sensitivity, Specificity, etc
mod_eval_train <- dismo::evaluate(p = occ.train.values2, a = absence.values2, model = ranaModelDismo) 
print(mod_eval_train)

mod_eval_test <- dismo::evaluate(p = occ.test.values2, a = absence.values2, model = ranaModelDismo)
print(mod_eval_test)  # training AUC may be higher than testing AUC


# calculate thresholds of models
thd1 <- threshold(mod_eval_train, "no_omission")  # 0% omission rate 
thd2 <- threshold(mod_eval_train, "spec_sens")  # highest TSS

plot(ranaPredictPlot, main='Maxent, raw values')
plot(wrld_simpl, add=TRUE, border='dark grey')
points(occ_train, pch='+')
points(occ_test, col = "pink")

plot(ranaPredictPlot > thd2, main='presence/absence')
plot(wrld_simpl, add=TRUE, border='dark grey')
points(occ_train, pch='+')
points(occ_test, col = "pink")



# model evaluation (cross validation, very fancy but good)

# plotting points that are above the previously calculated
# thresholded value
plot(ranaPredictPlot >= thd1)


# load the function that prepares parameters for maxent
source("../code/Appendix2_prepPara.R")

mod1_autofeature <- maxent(x=presence.absence.train.env.data[c("bio1","bio4","bio11")], 
                           ## env conditions, here we selected only 3 predictors
                           p=presence.absence.vector,
                           ## 1:presence or 0:absence
                           path="../output/maxent_outputs1_auto",
                           ## this is the folder you will find manxent output
                           args=prepPara(userfeatures=NULL) ) 
## default is autofeature

# or select Linear& Quadratic features
mod1_lq <- maxent(x=pder[c("bio1","bio4","bio11")],
                  p=pa,
                  path=paste0("../output/maxent_outputs1_lq"),
                  args=prepPara(userfeatures="LQ") ) 
## default is autofeature, here LQ represents Linear& Quadratic
## (L-linear, Q-Quadratic, H-Hinge, P-Product, T-Threshold)







#
# extracting env conditions for training occ from the raster
# stack; a data frame is returned (i.e multiple columns)
p <- raster::extract(clim, occ_train)
# env conditions for testing occ
p_test <- raster::extract(clim, occ_test)
# extracting env conditions for background
a <- raster::extract(clim, bg)


pa <- c(rep(1, nrow(p)), rep(0, nrow(a)))

# (rep(1,nrow(p)) creating the number of rows as the p data
# set to have the number '1' as the indicator for presence;
# rep(0,nrow(a)) creating the number of rows as the a data
# set to have the number '0' as the indicator for absence;
# the c combines these ones and zeros into a new vector that
# can be added to the Maxent table data frame with the
# environmental attributes of the presence and absence
# locations
pder <- as.data.frame(rbind(p, a))


# train Maxent with tabular data
mod <- maxent(x=pder, ## env conditions
              p=pa,   ## 1:presence or 0:absence
              
              path=paste0("maxent_outputs"), ## folder for maxent output; 
              # if we do not specify a folder R will put the results in a temp file, 
              # and it gets messy to read those. . .
              args=c("responsecurves") ## parameter specification
)
# the maxent functions runs a model in the default settings. To change these parameters,
# you have to tell it what you want...i.e. response curves or the type of features

# view the maxent model 
mod
plot(mod)
response(mod)

# predict to dataset

# first crop environment to the local species range +/- 10 degrees
model.extent<-extent(min(ranaData$longitude)-30,max(ranaData$longitude)+30,min(ranaData$latitude)-30,max(ranaData$latitude)+30)
modelEnv=crop(currentEnv,model.extent)


land <- crop(clim, (25+extent(ranaData)))
r <- raster::predict(mod, land)



plot(r)








############ mila cleaning to here

# Determine minimum and maximum values of latidue and longitudegitude

#' Finds minimum and maximum latidue and longitudegitude
MinMaxCoordinates <- function(x, padding = 0.1) {
  # If passed a single data.frame, wrap in a list
  if(class(x) == "data.frame") {
    x <- list(x)
  }
  
  # Establish starting min/max values
  max.latitude <- -90
  min.latitude <- 90
  max.longitude <- -180
  min.longitude <- 180
  
  # Iterate over all elements of list x and find min/max values
  for (i in 1:length(x)) {
    max.latitude <- ceiling(max(x[[i]]$latitude, max.latitude))
    min.latitude <- floor(min(x[[i]]$latitude, min.latitude))
    max.longitude <- ceiling(max(x[[i]]$longitude, max.longitude))
    min.longitude <- floor(min(x[[i]]$longitude, min.longitude))
  }
  
  if (padding > 0) {
    # Pad the values a bit so we don't end up with straight line distribution edges
    longitude.pad <- padding * (max.longitude - min.longitude)
    latitude.pad <- padding * (max.latitude - min.latitude)
    max.latitude <- max.latitude + latitude.pad
    min.latitude <- min.latitude - latitude.pad
    max.longitude <- max.longitude + longitude.pad
    min.longitude <- min.longitude - longitude.pad
  }
  
  # Format results and return
  min.max.coords <- c(min.longitude, max.longitude, min.latitude, max.latitude)
  names(min.max.coords) <- c("min.longitude", "max.longitude", "min.latitude", "max.latitude")
  return(min.max.coords)
}

# Determine geographic extent of our data
min.max <- MinMaxCoordinates(x = ranaDataNotCoords, padding = 1)
geographic.extent <- extent(x = min.max)


# Create pseudo-absence points (making them up, using 'background' approach)
# raster.files <- list.files(path = paste0(system.file(package = "dismo"), "/ex"),
#                            pattern = "grd", full.names = TRUE)
# mask <- raster(raster.files[1])
bil.files <- list.files(path = "wc2-5", 
                        pattern = "*.bil$", 
                        full.names = TRUE)
mask <- raster(bil.files[1])

# Random points for background (same number as our observed points)
set.seed(19470909)
background.extent <- extent(x = MinMaxCoordinates(x = ranaDataNotCoords, padding = 0.0))
# raster package will complain about not having coordinate reference system,
# so we suppress that warning
background.points <- suppressWarnings(randomPoints(mask = mask, 
                                                   n = nrow(ranaDataNotCoords), 
                                                   ext = background.extent, 
                                                   extf = 1.25))
colnames(background.points) <- c("longitude", "latitude")

# Data for observation sites (presence and background)
presence.values <- raster::extract(x = bioclim_data, y = ranaData)
absence.values <- raster::extract(x = bioclim_data, y = background.points)

# get the same random sample for training and testing
set.seed(1)

# randomly select 50% for training
selected <- sample(1:nrow(ranaDataNotCoords), nrow(ranaDataNotCoords) * 0.5)

occ_train <- ranaDataNotCoords[selected, ]  # this is the selection to be used for model training
occ_test <- ranaDataNotCoords[-selected, ]  # this is the opposite of the selection which will be used for model testing



# extracting env conditions for training occ from the raster
# stack; a data frame is returned (i.e multiple columns)
p <- raster::extract(clim, occ_train)
# env conditions for testing occ
p_test <- raster::extract(clim, occ_test)
# extracting env conditions for background
a <- raster::extract(clim, bg)


pa <- c(rep(1, nrow(p)), rep(0, nrow(a)))

# (rep(1,nrow(p)) creating the number of rows as the p data
# set to have the number '1' as the indicator for presence;
# rep(0,nrow(a)) creating the number of rows as the a data
# set to have the number '0' as the indicator for absence;
# the c combines these ones and zeros into a new vector that
# can be added to the Maxent table data frame with the
# environmental attributes of the presence and absence
# locations
pder <- as.data.frame(rbind(p, a))


# train Maxent with tabular data
mod <- maxent(x=pder, ## env conditions
              p=pa,   ## 1:presence or 0:absence
              
              path=paste0("maxent_outputs"), ## folder for maxent output; 
              # if we do not specify a folder R will put the results in a temp file, 
              # and it gets messy to read those. . .
              args=c("responsecurves") ## parameter specification
)
# the maxent functions runs a model in the default settings. To change these parameters,
# you have to tell it what you want...i.e. response curves or the type of features

# view the maxent model 
mod
plot(mod)
response(mod)

# predict to dataset

# first crop environment to the local species range +/- 10 degrees
model.extent<-extent(min(ranaData$longitude)-30,max(ranaData$longitude)+30,min(ranaData$latitude)-30,max(ranaData$latitude)+30)
modelEnv=crop(currentEnv,model.extent)


land <- crop(clim, (25+extent(ranaData)))
r <- raster::predict(mod, land)



plot(r)


# example 1, project to study area [raster]
ped1 <- raster::predict(mod,studyArea )  # studyArea is the clipped rasters 
plot(ped1)  # plot the continuous prediction


# this creates a 4-decimal-degree buffer around the
# occurrence data
occ_buff <- raster::buffer(foo, 10, dissolve = T)

# plot the first element ([[1]]) in the raster stack
plot(clim[[1]],  
     xlim = c(min_longitude, max_longitude),
     ylim = c(min_latitude, max_latitude))

plot(ranaData, add = T, col = "red")  # adds occurrence data to the plot
plot(occ_buff, col = "blue")  # adds buffer polygon to the plot

# crop study area to a manageable extent (rectangle shaped)
studyArea <- crop(clim,extent(occ_buff))  

# the 'study area' created by extracting the buffer area from the raster stack
studyArea <- mask(studyArea,occ_buff)
# output will still be a raster stack, just of the study area

# save the new study area rasters as ascii
writeRaster(studyArea,
            # a series of names for output files
            filename=paste0("studyarea/",names(studyArea),".asc"), 
            format="ascii", ## the output format
            bylayer=TRUE, ## this will save a series of layers
            overwrite=T)

set.seed(1) 
bg <- sampleRandom(x=studyArea,
                   size=10000,
                   na.rm=T, #removes the 'Not Applicable' points  
                   sp=T) # return spatial points 

plot(studyArea[[1]])
# add the background points to the plotted raster
plot(bg,add=T) 
# add the occurrence data to the plotted raster
plot(occ_final,add=T,col="red")

# Determine geographic extent of our data
max_latitude <- ceiling(max(ranaData$latidue))
min_latitude <- floor(min(ranaData$latidue))
max_longitude <- ceiling(max(ranaData$longitudegitude))
min_longitude <- floor(min(ranaData$longitudegitude))
geographic_extent <- extent(x = c(min_longitude, max_longitude, min_latitude, max_latitude))

# Load the data to use for our base map
data(wrld_simpl)

## observations 

# Plot the base map
plot(wrld_simpl, 
     xlim = c(min_longitude, max_longitude),
     ylim = c(min_latitude, max_latitude),
     axes = TRUE, 
     col = "grey95")

# Add the points for individual observation
points(x = ranaData$longitudegitude, 
       y = ranaData$latidue, 
       col = "olivedrab", 
       pch = 20, 
       cex = 0.75)
# And draw a little box around the graph
box()

# Crop bioclim data to geographic extent 
bioclim_data <- crop(x = bioclim_data, y = geographic_extent)

# view one layer
plot(bioclim_data, )













# with bioclim
# Build species distribution model (I need absence points)
bc_model <- bioclim(x = bioclim_data, p = ranaData)

predict_presence <- dismo::predict(object = bc_model, 
                                   x = bioclim_data)

# Plot base map
plot(wrld_simpl, 
     xlim = c(min_longitude, max_longitude),
     ylim = c(min_latitude, max_latitude),
     axes = TRUE, 
     col = "grey95")

# Add model probabilities
plot(predict_presence, add = TRUE)

# Redraw those country borders
plot(wrld_simpl, add = TRUE, border = "grey5")




