###############################################
#### Caffea Arabica: Modelin and Prediction ###
###############################################

#' Species distribution models (SDM) using caffea arabica data from GBIF
#' ============================================
#' May 2019 (24.06.2019)
#' Magdalena Halbgewachs
#' R version 3.5.1
#' Occurrence data source: [gbif](http://www.gbif.org)  
#' Environmental data source: [worldclim](http://www.worldclim.org)  
#' Algorithms: GAM, RF
#' 
#' Overview
#' ===============================
#' - Set working directory and load packages
#' - Import data
#' - Data preprocessing
#'  * Presence records for which environmental data is available
#'  * Collinearity
#' - Modelling
#' - Model evaluation
#'  * Model performance on test data
#'  * Variable importance
#'  * Response functions
#' - Model predictions
#' 

############################################################################################################
#' Set working directory and load packages
setwd("C:\\02_Studium\\02_Master\\02_Semester_2\\MET1_Spatial_modelling_and_prediction\\caffea_arabica")

library(rms)
#install.packages("raster")
library(raster)
library(mgcv) # gam
#' We assume that the file "varImpBiomod.R" is in the working directory
source("varImpBiomod.R") 
library(randomForest)
library(dismo)
library(rgdal)
library(ellipse)
library(randomForest)
library(rJava)
library(XML)
library(sp)

#################################################################################################################
#' Import data: Environment and species occurrences
#' Import shape of the study area (Ethiopia)

study_area <- readOGR("./ethiopia_shp", "ETH_outline")
plot(study_area)

#' Download and read bioclim variables
#' -----------------------------------
#' The function getData will load the data from the working directory, if available. If they are not available,
#' the function will attempt to download the data from the internet.
#' For a definition of the variables, see [http://www.worldclim.org/bioclim](http://www.worldclim.org/bioclim)
#' 
#' Variable | Description
#' -------- | -----------
#' BIO1 | Annual Mean Temperature
#' BIO2 | Mean Diurnal Range (Mean of monthly (max temp - min temp))
#' BIO3 | Isothermality (BIO2/BIO7) (* 100)
#' BIO4 | Temperature Seasonality (standard deviation *100)
#' BIO5 | Max Temperature of Warmest Month
#' BIO6 | Min Temperature of Coldest Month
#' BIO7 | Temperature Annual Range (BIO5-BIO6)
#' BIO8 | Mean Temperature of Wettest Quarter
#' BIO9 | Mean Temperature of Driest Quarter
#' BIO10 | Mean Temperature of Warmest Quarter
#' BIO11 | Mean Temperature of Coldest Quarter
#' BIO12 | Annual Precipitation
#' BIO13 | Precipitation of Wettest Month
#' BIO14 | Precipitation of Driest Month
#' BIO15 | Precipitation Seasonality (Coefficient of Variation)
#' BIO16 | Precipitation of Wettest Quarter
#' BIO17 | Precipitation of Driest Quarter
#' BIO18 | Precipitation of Warmest Quarter
#' BIO19 | Precipitation of Coldest Quarter
#' 
#' getdata(climatevariable)

bio <- raster::getData("worldclim", var = "bio", res = 2.5)

#' Plot the first raster layer, i.e. annual mean temperature
plot(raster(bio, 1))
# Add the outline of the study area
plot(study_area, add=TRUE)

#' Crop to study area extent (with a 3 degree buffer in each direction)
biocrop <- crop(bio, extent(study_area) +3)

#' Plot the first raster layer of the cropped climate data
plot(raster(biocrop, 7))
plot(study_area, add=TRUE)

###########################################################################################################
#' ==========================================
#' Read occurrence points
#' ==========================================

# Download species location data from gbif
# species <- gbif("Syncerus", "caffer", ext = extent(bio), sp = TRUE, removeZeros = TRUE)
species0 <- gbif('Coffea arabica L.')
species <- subset(species0,select=c("lat","lon"))
species <- na.omit(species)
coordinates(species) <- c("lon", "lat")  # set spatial coordinates
  
  
# Add projection information
proj4string(species) <- CRS("+proj=longlat +datum=WGS84")
# Save species records in mif-format (preserves full column names)
#writeOGR(species, "./GIS/Synceruscaffer", "Synceruscaffer", driver="MapInfo File", dataset_options="FORMAT=MIF")

plot(raster(biocrop, 1))
plot(study_area, add=TRUE)
plot(species, add = TRUE)

############################################################################################################
#' ==============================
#' Data preprocessing
#' ==============================

#' Select species records for which environmental information is available
#' -------------------------------
species <- species[complete.cases(extract(biocrop, species)), ]

#' Collinearity
#' -----------------------------
#' ### Visual inspection of collinearity ###
cm <- cor(getValues(bio), use = "complete.obs")
plotcorr(cm, col=ifelse(abs(cm) > 0.7, "red", "grey"))

#' ### Select an uncorrelated subset of environmental variables ###
env <- subset(biocrop, c("bio3", "bio5", "bio6", "bio12", "bio18"))
#env <- subset(biocrop, c("bio1", "bio4", "bio5", "bio12", "bio14", "bio15")) #selfdefined

############################################################################################################
#' ==========================================
#' Sampling of (pseudo-)absence points
#' ==========================================
#' The function randomPoints in package dismo allows to 
#' randomly select a certain number of random points,
#' and to adjust the probability of selecting a cell
#' according to its size, which is relevant in lat-lon-grids,
#' where cells are of differing size

#' Selecting 2000 random background points, excluding cells where
#' the species is present
set.seed(2)
background <- randomPoints(env, 2000, species)
#' Select only one presence record in each cell of the environmental layer
presence <- gridSample(species, env, n = 1)

#' combining the presence and background points, adding a 
#' column "species" that contains the information about presence (1)
#' and background (0)
fulldata <- SpatialPointsDataFrame(rbind(presence, background),
                                   data = data.frame("species" = rep(c(1,0), 
                                                                     c(nrow(presence), nrow(background)))),
                                   match.ID = FALSE,
                                   proj4string = CRS(projection(env)))
#' Add information of environmental conditions at point locations
fulldata@data <- cbind(fulldata@data, extract(env, fulldata))

# Split data set into a training and test data set
set.seed(2)
fold <- kfold(fulldata, k = 5)
traindata <- fulldata[fold != 1, ]
testdata <- fulldata[fold == 1, ]

############################################################################################################
#' ==========================================
#' Evaluate model on test data
#' ==========================================
#' We can now use a range of statistical methods to estimate the
#' probability of species occurrence.
#' Unfortunately, there are often subtle differences in how the models
#' are specified and in which data formats are useable

varnames <- c("bio3", "bio5", "bio6", "bio12", "bio18")  #paper
#varnames <- c("bio1", "bio4", "bio5", "bio12", "bio14", "bio15") #selfdenied


#' ==========================================
#' GAM algorithm (Generalized additive models)
#' ==========================================
gammodel <- gam(species ~ s(bio3) + s(bio5) + s(bio6) + s(bio12) + s(bio18),
                family="binomial", data=traindata)
#gammodel <- gam(species ~ s(bio1) + s(bio4) + s(bio5) + s(bio12) + s(bio14) + s(bio15),
                #family="binomial", data=traindata)  #selfdefined
summary(gammodel)

plot(gammodel)

# a) Predict to test data
gamtest <- predict(gammodel, newdata = testdata, type = "response")
# b) Calculate performance indices
val.prob(gamtest, testdata[["species"]])

# Variable importance
gamimp <- varImpBiomod(gammodel, varnames,
                       traindata)
barplot(100 * gamimp/sum(gamimp), ylab = "Variable importance (%)")

# Response functions --> how good are the variables selected?
plot(gammodel, pages = 1)

# png("gammodel_resp.png", 800, 800)
# plot(gammodel, pages = 1)
# dev.off()

# Prediction map
gammap <- predict(env, gammodel, type = "response")

plot(gammap)
plot(study_area, add=TRUE)

#' ==========================================
## Random forest
#' ==========================================
# randomForest requires the dependent variable to be a factor
# if we want to do classification
rftraindata <- as(traindata, "data.frame")
rftraindata$species <- factor(rftraindata$species)

# TODO: check proper settings of random forest algorithm
rfmodel <- randomForest(species ~ bio1 + bio4 + bio5 + bio12 + bio14 + bio15, data = rftraindata)

# Evaluate model on test data
# a) Predict to test data
rftest <- predict(rfmodel, newdata = testdata, type = "prob")[,2]
# b) Calculate performance indices
val.prob(rftest, testdata[["species"]])

# Variable importance
rfImp <- importance(rfmodel)
varImpPlot(rfmodel)

# Response functions
par(mfrow=c(3,2))
for (i in seq_along(varnames)) {
  partialPlot(rfmodel, rftraindata, varnames[i], xlab = varnames[i], main="")  
}

# Prediction map
rfmap <- predict(env, rfmodel, type = "prob", index = 2)
par(mfrow=c(1, 1))
plot(rfmap)
plot(study_area, add=TRUE)


############################################################################
# prediction
############################################################################

#load datasets
#stack all rasters of one scenario and year
#crop layer to study area extent with a buffer of 3
#rename layer --> must be equal to actual climate data

rcp2_2050_list <- list.files(path="D:\\01_Uni\\02_Master\\MET1_Modeling_Prediction\\RCP\\mpi_esm_lr_rcp2_6_2050s_bio_30s_r1i1p1_b4_asc\\bio_b4", pattern = ".asc", full.names = TRUE)
rcp2_2050 <- raster::stack(rcp2_2050_list)
rcp2_2050_crop <- crop(rcp2_2050, extent(study_area) + 3)
names(rcp2_2050_crop) <- gsub("_", "", names(rcp2_2050_crop))
plot(rcp2_2050_crop[[1]])
plot(study_area, add=T)


rcp2_2080_list <- list.files(path="D:\\01_Uni\\02_Master\\MET1_Modeling_Prediction\\RCP\\mpi_esm_lr_rcp2_6_2080s_bio_30s_r1i1p1_b4_asc\\bio_b4", pattern = ".asc", full.names = TRUE)
rcp2_2080 <- raster::stack(rcp2_2080_list)
rcp2_2080_crop <- crop(rcp2_2080, extent(study_area) + 3)
names(rcp2_2080_crop) <- gsub("_", "", names(rcp2_2080_crop))
plot(rcp2_2080_crop[[1]])
plot(study_area, add=T)

rcp4_2050_list <- list.files(path="D:\\01_Uni\\02_Master\\MET1_Modeling_Prediction\\RCP\\mpi_esm_lr_rcp4_5_2050s_bio_30s_r1i1p1_b4_asc\\bio_b4", pattern = ".asc", full.names = TRUE)
rcp4_2050 <- raster::stack(rcp4_2050_list)
rcp4_2050_crop <- crop(rcp4_2050, extent(study_area) + 3)
names(rcp4_2050_crop) <- gsub("_", "", names(rcp4_2050_crop))
plot(rcp4_2050_crop[[1]])
plot(study_area, add=T)

rcp4_2080_list <- list.files(path="D:\\01_Uni\\02_Master\\MET1_Modeling_Prediction\\RCP\\mpi_esm_lr_rcp4_5_2080s_bio_30s_r1i1p1_b4_asc\\bio_b4", pattern = ".asc", full.names = TRUE)
rcp4_2080 <- raster::stack(rcp4_2080_list)
rcp4_2080_crop <- crop(rcp4_2080, extent(study_area) + 3)
names(rcp4_2080_crop) <- gsub("_", "", names(rcp4_2080_crop))
plot(rcp4_2080_crop[[1]])
plot(study_area, add=T)

rcp8_2050_list <- list.files(path="D:\\01_Uni\\02_Master\\MET1_Modeling_Prediction\\RCP\\mpi_esm_lr_rcp8_5_2050s_bio_30s_r1i1p1_b4_asc\\bio_b4", pattern = ".asc", full.names = TRUE)
rcp8_2050 <- raster::stack(rcp8_2050_list)
rcp8_2050_crop <- crop(rcp8_2050, extent(study_area) + 3)
names(rcp8_2050_crop) <- gsub("_", "", names(rcp8_2050_crop))
plot(rcp8_2050_crop[[1]])
plot(study_area, add=T)

rcp8_2080_list <- list.files(path="D:\\01_Uni\\02_Master\\MET1_Modeling_Prediction\\RCP\\mpi_esm_lr_rcp8_5_2080s_bio_30s_r1i1p1_b4_asc\\bio_b4", pattern = ".asc", full.names = TRUE)
rcp8_2080 <- raster::stack(rcp8_2080_list)
rcp8_2080_crop <- crop(rcp8_2080, extent(study_area) + 3)
names(rcp8_2080_crop) <- gsub("_", "", names(rcp8_2080_crop))
plot(rcp8_2080_crop[[1]])
plot(study_area, add=T)

#-----------------------------------------
# GAM model for future climate scenarios
#-----------------------------------------
#select subset of uncorrelated variables
#use the previous calculates GAM model for future prediction

env_rcp2_2050 <- subset(rcp2_2050_crop, c("bio3", "bio5", "bio6", "bio12", "bio18"))#c("bio1", "bio4", "bio5", "bio12", "bio14", "bio15"))
gammap_rcp2_2050 <- predict(env_rcp2_2050, gammodel, type = "response")
plot(gammap_rcp2_2050)
plot(study_area, add=T)
#writeRaster(class, filename = paste0(binary_out_folder, "/class_", date, ".tif"), format = "GTiff", overwrite = T)

env_rcp2_2080 <- subset(rcp2_2080_crop, c("bio3", "bio5", "bio6", "bio12", "bio18"))
gammap_rcp2_2080 <- predict(env_rcp2_2080, gammodel, type = "response")
plot(gammap_rcp2_2080)
plot(study_area, add=T)

env_rcp4_2050 <- subset(rcp4_2050_crop, c("bio3", "bio5", "bio6", "bio12", "bio18"))
gammap_rcp4_2050 <- predict(env_rcp4_2050, gammodel, type = "response")
plot(gammap_rcp4_2050)
plot(study_area, add=T)

env_rcp4_2080 <- subset(rcp4_2080_crop, c("bio3", "bio5", "bio6", "bio12", "bio18"))
gammap_rcp4_2080 <- predict(env_rcp4_2080, gammodel, type = "response")
plot(gammap_rcp4_2080)
plot(study_area, add=T)

env_rcp8_2050 <- subset(rcp8_2050_crop, c("bio3", "bio5", "bio6", "bio12", "bio18"))
gammap_rcp8_2050 <- predict(env_rcp8_2050, gammodel, type = "response")
plot(gammap_rcp8_2050)
plot(study_area, add=T)

env_rcp8_2080 <- subset(rcp8_2080_crop, c("bio3", "bio5", "bio6", "bio12", "bio18"))
gammap_rcp8_2080 <- predict(env_rcp8_2080, gammodel, type = "response")
plot(gammap_rcp8_2080)
plot(study_area, add=T)

#-------------------------------------------------------------------------------------------
#creating binary map with treshold 0.7 -> 70% probability of growing caffea arabica
class_out_folder <- "C:\\02_Studium\\02_Master\\02_Semester_2\\MET1_Spatial_modelling_and_prediction\\Caffea_arabica\\data\\ggRGB"
binaryMap <- function(raster, threshold) {
  bin <- raster
  bin[raster <= threshold] <- NA
  bin[raster > threshold] <- bin[raster]
  return(bin)
}

threshold <- 0.7

#binary map
ras <- brick(gammap)
class_current <- binaryMap(ras, threshold)
proj4string(class_current) <- CRS("+proj=longlat +datum=WGS84")
class_current <- raster::resample(class_current, class_rcp2_2050)
class_current <- raster::mask(class_current, study_area)
plot(class_current)
plot(study_area, add=T)
writeRaster(class_current, filename = "class_current", format = "GTiff", overwrite = T)
#raster to polygon
pol_class_current <- rasterToPolygons(class_current, dissolve = T)
plot(pol_class_current)
writeOGR(pol_class_current, ".", "pol_class_current", driver="ESRI Shapefile")

ras <- brick(gammap_rcp2_2050)
class_rcp2_2050 <- binaryMap(ras, threshold)
proj4string(class_rcp2_2050) <- CRS("+proj=longlat +datum=WGS84")
class_rcp2_2050 <- raster::mask(class_rcp2_2050, study_area)
plot(class_rcp2_2050)
plot(study_area, add=T)
writeRaster(class_rcp2_2050, filename = "class_rcp2_2050", format = "GTiff", overwrite = T)
pol_class_rcp2_2050 <- rasterToPolygons(class_rcp2_2050, dissolve = T)
plot(pol_class_rcp2_2050)
writeOGR(pol_class_rcp2_2050, ".", "pol_class_rcp2_2050", driver="ESRI Shapefile")


ras <- brick(gammap_rcp2_2080)
class_rcp2_2080 <- binaryMap(ras, threshold)
proj4string(class_rcp2_2080) <- CRS("+proj=longlat +datum=WGS84")
class_rcp2_2080 <- raster::mask(class_rcp2_2080, study_area)
plot(class_rcp2_2080)
plot(study_area, add=T)
writeRaster(class_rcp2_2080, filename = "class_rcp2_2080", format = "GTiff", overwrite = T)
pol_class_rcp2_2080 <- rasterToPolygons(class_rcp2_2080, dissolve = T)
plot(pol_class_rcp2_2080)
writeOGR(pol_class_rcp2_2080, ".", "pol_class_rcp2_2080", driver="ESRI Shapefile")


ras <- brick(gammap_rcp4_2050)
class_rcp4_2050 <- binaryMap(ras, threshold)
proj4string(class_rcp4_2050) <- CRS("+proj=longlat +datum=WGS84")
class_rcp4_2050 <- raster::mask(class_rcp4_2050, study_area)
plot(class_rcp4_2050)
plot(study_area, add=T)
writeRaster(class_rcp4_2050, filename = "class_rcp4_2050", format = "GTiff", overwrite = T)
pol_class_rcp4_2050 <- rasterToPolygons(class_rcp4_2050, dissolve = T)
plot(pol_class_rcp4_2050)
writeOGR(pol_class_rcp4_2050, ".", "pol_class_rcp4_2050", driver="ESRI Shapefile")


ras <- brick(gammap_rcp4_2080)
class_rcp4_2080 <- binaryMap(ras, threshold)
proj4string(class_rcp4_2080) <- CRS("+proj=longlat +datum=WGS84")
class_rcp4_2080 <- raster::mask(class_rcp4_2080, study_area)
plot(class_rcp4_2080)
plot(study_area, add=T)
writeRaster(class_rcp4_2080, filename = "class_rcp4_2080", format = "GTiff", overwrite = T)
pol_class_rcp4_2080 <- rasterToPolygons(class_rcp4_2080, dissolve = T)
plot(pol_class_rcp4_2080)
writeOGR(pol_class_rcp4_2080, ".", "pol_class_rcp4_2080", driver="ESRI Shapefile")


ras <- brick(gammap_rcp8_2050)
class_rcp8_2050 <- binaryMap(ras, threshold)
proj4string(class_rcp8_2050) <- CRS("+proj=longlat +datum=WGS84")
class_rcp8_2050 <- raster::mask(class_rcp8_2050, study_area)
plot(class_rcp8_2050)
plot(study_area, add=T)
writeRaster(class_rcp8_2050, filename = "class_rcp8_2050", format = "GTiff", overwrite = T)
pol_class_rcp8_2050 <- rasterToPolygons(class_rcp8_2050, dissolve = T)
plot(pol_class_rcp8_2050)
writeOGR(pol_class_rcp8_2050, ".", "pol_class_rcp8_2050", driver="ESRI Shapefile")


ras <- brick(gammap_rcp8_2080)
class_rcp8_2080 <- binaryMap(ras, threshold)
proj4string(class_rcp8_2080) <- CRS("+proj=longlat +datum=WGS84")
class_rcp8_2080 <- raster::mask(class_rcp8_2080, study_area)
plot(class_rcp8_2080)
plot(study_area, add=T)
writeRaster(class_rcp8_2080, filename = "class_rcp8_2080", format = "GTiff", overwrite = T)
pol_class_rcp8_2080 <- rasterToPolygons(class_rcp8_2080, dissolve = T)
plot(pol_class_rcp8_2080)
writeOGR(pol_class_rcp8_2080, ".", "pol_class_rcp8_2080", driver="ESRI Shapefile")

#----------------------------------------------------
#calculate area of binary results

calc_area <- function(raster, res_x, res_y){
  #define default value for res_x and res_y
  if (missing(res_x) & missing(res_y)){
    res_x <- 0.92
    res_y <- 0.92
  }
  vals <- getValues(raster)
  cells <- sum(vals ==1, na.rm = T)
  area <-  sum(vals ==1, na.rm = T) * res_x * res_y
  return(area)
}

area_current <- calc_area(class_current)  #pixel 4.6km
area_rcp2_2050 <- calc_area(class_rcp2_2050)  #pixel 0.92km
area_rcp2_2080 <- calc_area(class_rcp2_2080)
area_rcp4_2050 <- calc_area(class_rcp4_2050)
area_rcp4_2080 <- calc_area(class_rcp4_2080)
area_rcp8_2050 <- calc_area(class_rcp8_2050)
area_rcp8_2080 <- calc_area(class_rcp8_2080)
  
  
#------------------------------------------------
# Plot per each scenario
library(RStoolbox)

rcp2_list <- list.files(path="C:\\02_Studium\\02_Master\\02_Semester_2\\MET1_Spatial_modelling_and_prediction\\Caffea_arabica\\data\\ggRGB\\class_rcp2", full.names = T)
stack_rcp2 <- raster::stack(rcp2_list)
ggRGB(stack_rcp2, r = 3, g = 2, b = 1, alpha=0.5)+
  #xlab("test")+
  ggtitle("Prediction in Caffea Arabica for RCP2")+
  theme_void()  # Empty theme without axis lines and texts
 
rcp4_list <- list.files(path="C:\\02_Studium\\02_Master\\02_Semester_2\\MET1_Spatial_modelling_and_prediction\\Caffea_arabica\\data\\ggRGB\\class_rcp4", full.names = T)
stack_rcp4 <- raster::stack(rcp4_list)
plot_rcp4 <- ggRGB(stack_rcp4, r = 3, g = 2, b = 1)
plot_rcp4

writeRaster(stack_rcp4, filename="RGB_RCP4.tif", options="INTERLEAVE=BAND")

rcp8_list <- list.files(path="C:\\02_Studium\\02_Master\\02_Semester_2\\MET1_Spatial_modelling_and_prediction\\Caffea_arabica\\data\\ggRGB\\class_rcp8", full.names = T)
stack_rcp8 <- raster::stack(rcp8_list)
plot_rcp8 <- ggRGB(stack_rcp8, r = 3, g = 2, b = 1)
plot_rcp8

gammap_res <- raster::resample(gammap, class_rcp2_2050)
gam_rcp2_stack <- stack(gammap_res, gammap_rcp2_2050, gammap_rcp2_2080)
writeRaster(gam_rcp2_stack, filename="RGB_RCP2_full.tif", options="INTERLEAVE=BAND")
ggRGB(gam_rcp2_stack, r = 3, g = 2, b = 1)

#-------------------------------------------------
#elevation data
setwd("D:\\01_Uni\\02_Master\\MET1_Modeling_Prediction")
#install.packages("elevatr")
library(elevatr)
data(study_area)
elevation <- get_elev_raster(gammap, z = 9)
plot(elevation)


srtm <- raster::getData('SRTM', lon=8, lat=39)
plot(srtm)
plot(study_area, add=T)


require(devtools)
install.packages("rjson")
install_github("ropensci/geonames")
library(rjson)
library(geonames)
topo30 <- GNgtopo30(lat=39,lng=8)







# TODO: 
#       calculate area of binary results
#       elevation miteinbeziehen (corrplot höhe-wahrscheinlichkeit occurence) > Höhenlinine auch plotten https://plot.ly/r/contour-plots/
#       https://cran.r-project.org/web/packages/elevatr/vignettes/introduction_to_elevatr.html
#       https://gis.stackexchange.com/questions/293993/plotting-and-analyzing-extracted-elevation-data-in-r
#       https://earthdata.nasa.gov/: srtm, topo30 gtopo 1km, aspect slope, resample
#       https://www.rdocumentation.org/packages/geonames/versions/0.999/topics/GNgtopo30
#       gdal merge

#       doku where and how to download bioclim data step by step
#       Korrelation mit Kaffeekirschenkäfer? Osttimor

#       andere bios nehmen? --> paper
#       crop all to extent!


#https://search.earthdata.nasa.gov/data/retrieve/0375455176




