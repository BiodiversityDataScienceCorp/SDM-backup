# one file test


# contents of setup.R:

#########################################################

# Setup script for required data and package installation
# Jeffrey C. Oliver
# jcoliver@email.arizona.edu
# 2017-11-02

################################################################################
# SUMMARY
# Install dependencies
# Download climate data

################################################################################
# Install dependencies
required <- c("raster", "sp", "dismo", "maptools", "spocc", "tidyverse")
install.packages(required)
library("spocc")
library(tidyverse)

# Make sure packages all installed
successful <- required %in% rownames(installed.packages())
unsuccessful <- required[!successful]

if (length(unsuccessful) > 0) {
  unsuccessful.string <- paste0(unsuccessful, collapse = ", ")
  stop(paste0("One or more required packages could not be installed: ", 
              unsuccessful.string))
}

################################################################################
# Download climate data

# Make sure data directory is writable
if (file.access(names = "data") != 0) {
  stop(paste0("You do not have sufficient write access to data directory.\n"))
}

# Make sure raster package was installed and load it
if (!require(package = "raster")) {
  stop("Setup requires the raster package, which does not appear to be available.\n")
}

# Download bioclim data
message("Downloading climate data from WorldClim")
bioclim.data <- getData(name = "worldclim",
                        var = "bio",
                        res = 2.5, # Could try for better resolution, 0.5, but would then need to provide lat & long...
                        path = "data/")

# Unzip forecast data
message("Extracting forecast climate data (this may take a moment)")
forecast.archives <- list.files(path = "data/cmip5/2_5m", 
                                pattern = "*.zip$",
                                full.names = TRUE)
forecast.data <- lapply(X = forecast.archives, FUN = unzip)
# unzip(zipfile = "data/cmip5/2_5m/forecast-data.zip")

# NOPE archive is too big (> 100 MB) for GitHub. But there might be a solution
# GitHub large file storage https://git-lfs.github.com/
# Better yet, just make a few (4?) smaller archives

# Downloading of forecast data deprecated to avoid dependency on troublesome 
# rgdal. Instead use 
#   `forecast.data <- raster::stack(x = "data/cmip5/2_5m/forecast-raster.gri")`
# when forecast data are needed
# Download forecast data
# See https://link.springer.com/article/10.1007/s00382-014-2418-8
# for recommendations of the model to use
# forecast.data <- getData(name = "CMIP5", # forecast
#                          var = "bio", # bioclim
#                          res = 2.5,
#                          path = "data/",
#                          model = "GD", # GFDL-ESM2G
#                          rcp = "45", # CO2 increase 4.5
#                          year = 70) # 2070
# For those interested, the workaround was:
#  1. With rgdal installed, use the getData code as above
#  2. With the `bioclim.data` object in memory, run 
#      `names(forecast.data) <- names(bioclim.data)`
#  3. Save the file as a raster using `raster::writeRaster` with default 
#      format (.gri)
#  4. Compress the resultant .gri and .grd files via
#     `zip(zipfile = "data/cmip5/2_5m/forecast-data.zip", 
#          files = c("data/cmip5/2_5m/forcast-raster.grd", 
#                    "data/cmip5/2_5m/forcast-raster.gri"))`
#  6. Include that zip archive under version control, but not the .gri and .grd
#      files
#  5. Update this file (scripts/setup.R) to unzip the archive, inflating the 
#      .grd and .gri files (directory structure was preserved by `zip` command)

# Clean up workspace
rm(required, successful, unsuccessful, bioclim.data, forecast.archives, forecast.data)

message("Setup complete.")

####################################################################

# end setup.R


# query gbif data, and save as CSV

########

ranaQuery<-occ(query='Rana draytonii', from="gbif")

#navigate object
ranaData<-ranaQuery$gbif$data$Rana_draytonii



################ 5. Save to CSV #####################

# first, ensure all data is character data
#df <- apply(df,2,as.character)

ranaData <- apply(ranaData,2,as.character)

# use write.csv to write the data frame to 'data' directory
# make sure the file name matches what you indicated in step 3 on line 14
write.csv(ranaData, "data/rana.csv")




#Now run though sdm-single code:
########################################



# Script to run contemporary species distribution model for a single species
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2017-09-07

rm(list = ls())

################################################################################
# SETUP
# Gather path information
# Load dependancies

# Things to set:
infile <- "data/rana.csv"
outprefix <- "rana"
outpath <- "output/"

# Make sure the input file exists
if (!file.exists(infile)) {
  stop(paste0("Cannot find ", infile, ", file does not exist.\n"))
}

# Make sure the input file is readable
if (file.access(names = infile, mode = 4) != 0) {
  stop(paste0("You do not have sufficient access to read ", infile, "\n"))
}

# Make sure the output path ends with "/" (and append one if it doesn't)
if (substring(text = outpath, first = nchar(outpath), last = nchar(outpath)) != "/") {
  outpath <- paste0(outpath, "/")
}

# Make sure directories are writable
required.writables <- c("data", outpath)
write.access <- file.access(names = required.writables)
if (any(write.access != 0)) {
  stop(paste0("You do not have sufficient write access to one or more directories. ",
              "The following directories do not appear writable: \n",
              paste(required.writables[write.access != 0], collapse = "\n")))
}

# Load dependancies, keeping track of any that fail
required.packages <- c("raster", "sp", "dismo", "maptools")
missing.packages <- character(0)
for (one.package in required.packages) {
  if (!suppressMessages(require(package = one.package, character.only = TRUE))) {
    missing.packages <- cbind(missing.packages, one.package)
  }
}

if (length(missing.packages) > 0) {
  stop(paste0("Missing one or more required packages. The following packages are required for run-sdm: ", paste(missing.packages, sep = "", collapse = ", ")), ".\n")
}

source(file = "src/sdm-functions.R")

################################################################################
# ANALYSES
# Prepare data
# Run species distribution modeling
# Combine results from butterflies and plants

# Prepare data
prepared.data <- PrepareData(file = infile)

# Run species distribution modeling
sdm.raster <- SDMRaster(data = prepared.data)

################################################################################
# PLOT
# Determine size of plot
# Plot to pdf file

# Add small value to all raster pixels so plot is colored correctly
sdm.raster <- sdm.raster + 0.00001

# Determine the geographic extent of our plot
xmin <- extent(sdm.raster)[1]
xmax <- extent(sdm.raster)[2]
ymin <- extent(sdm.raster)[3]
ymax <- extent(sdm.raster)[4]

sdf <- rasterToPoints(sdm.raster, spatial = TRUE)
#sdffuture <- rasterToPoints(sdm.raster.future, spatial = TRUE)
# Then to a 'conventional' dataframe
rasterDF  <- data.frame(sdf)

# removes absence data
sdmRasterDF<-rasterDF %>% subset(layer>1)

wrld<-ggplot2::map_data("world", c("mexico", "canada"))
#az<-map_data("county", "arizona")


ggplot(prepared.data) +
  geom_tile(data = sdmRasterDF , aes(x = x, y = y, color="Current", fill="green"), show.legend=TRUE, alpha=0.1, col="green")+ 
#  geom_tile(data= sdmRasterDFfuture, aes(x=x, y=y, color="Future", fill="blue"), show.legend=TRUE, alpha=0.1, col="blue") +
  geom_point(aes(x=lon, y=lat, color='red'),  size=1.5) +
  borders("state", xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
  #geom_polygon(data=az, mapping=aes(x=long, y=lat,group = group), fill = NA, colour = "grey60") +
 scale_size_area() +
 coord_quickmap() +
  coord_fixed(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) 
#  labs(title="Current and Future SDM Predictions in Arizona", x="longitude", y="latitude")+ 
#  scale_fill_identity(name = 'Species Distribution Model', guide = 'legend',labels = c('Future', "Current"))+
 # scale_colour_manual(name = 'Occurrences', 
#                      values =c('red'='#FF5733'), labels = c('Pineneedle milkweed'))

ggsave(plot.file.sdm, presentFuture)



# Plot the model; save to pdf
plot.file <- paste0(outpath, outprefix, "-single-prediction.pdf")
pdf(file = plot.file, useDingbats = FALSE)

# Load in data for map borders
data(wrld_simpl)

# Draw the base map
plot(wrld_simpl, xlim = c(xmin, xmax), ylim = c(ymin, ymax), axes = TRUE, col = "gray95", 
     main = paste0(gsub(pattern = "_", replacement = " ", x = outprefix), " - current"))

# Add the model rasters
plot(sdm.raster, legend = FALSE, add = TRUE)

# Redraw the borders of the base map
plot(wrld_simpl, xlim = c(xmin, xmax), ylim = c(ymin, ymax), add = TRUE, border = "gray10", col = NA)

# Add bounding box around map
box()

# Stop re-direction to PDF graphics device
dev.off()

# Let user know analysis is done.
message(paste0("\nAnalysis complete. Map image written to ", plot.file, "."))

rm(list = ls())




