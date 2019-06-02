###############################################
#### Caffea Arabica: Modeling and Prediction ###
###############################################

#' Species distribution model (SDM) using caffea arabica data from GBIF
#' ============================================
#' May 2019 (24.06.2019)
#' Magdalena Halbgewachs
#' R version 3.5.1
#' Occurrence data source: [gbif](http://www.gbif.org)  
#' Environmental data source: [worldclim](http://www.worldclim.org)  
#' Algorithms: GAM
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
setwd("C:\\02_Studium\\02_Master\\02_Semester_2\\MET1_Spatial_modelling_and_prediction\\Caffea_arabica")

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
#' ==========================================
#' Import data: Environment and species occurrences
#' Import shape of the study area (Ethiopia)
#' #' ==========================================

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
#' BIO20 | Digital elevation model
#' BIO21 | Slope
#' getdata(climatevariable)

bio <- raster::getData("worldclim", var = "bio", res = 2.5)                           #climate data from Bioclim
gtopo30 <- raster("D:\\01_Uni\\02_Master\\MET1_Modeling_Prediction\\gt30e020n40.tif") #elevation data from USGS Earth Explorer
names(gtopo30) <- gsub("gt30e020n40","bio20", names(gtopo30))

#' Plot the first raster layer, i.e. annual mean temperature
plot(raster(bio, 1))
# Add the outline of the study area
plot(study_area, add=TRUE)

#' Crop to study area extent (with a 3 degree buffer in each direction)
biocrop <- crop(bio, extent(study_area) +3)
gtopo30_crop <- crop(gtopo30, extent(study_area) +3)

#' resample gtopo dataset to pixelsize of bioclilm data
gtopo30_res <- raster::resample(gtopo30_crop, biocrop)

#' calculate slope out of gtopo30 data
slope_out_folder <- "D:\\01_Uni\\02_Master\\MET1_Modeling_Prediction"
slope <- terrain(gtopo30_res, opt='slope', unit='degrees', neighbors=8, filename = paste0(slope_out_folder, "/slope.tif"))

#' assign same coordinate system
proj4string(biocrop) <- CRS("+proj=longlat +datum=WGS84")
proj4string(gtopo30_res) <- CRS("+proj=longlat +datum=WGS84")

#' stack datasets together
biocrop <- stack(biocrop, gtopo30_res)

#' Plot the first raster layer of the cropped and stacked climate data
plot(raster(biocrop, 20))
plot(study_area, add=TRUE)

###########################################################################################################
#' ==========================================
#' Read occurrence points
#' ==========================================

# Download species location data from gbif
species0 <- gbif('Coffea arabica L.')
species <- subset(species0,select=c("lat","lon"))
species <- na.omit(species)
coordinates(species) <- c("lon", "lat")  # set spatial coordinates


# Add projection information
proj4string(species) <- CRS("+proj=longlat +datum=WGS84")
# Save species records in mif-format (preserves full column names)
#writeOGR(species, "./GIS/Synceruscaffer", "Synceruscaffer", driver="MapInfo File", dataset_options="FORMAT=MIF")

plot(raster(biocrop, 20))
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
env <- subset(biocrop, c("bio3", "bio5", "bio6", "bio12", "bio18", "bio20"))

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

varnames <- c("bio3", "bio5", "bio6", "bio12", "bio18", "bio20")

#' ==========================================
#' GAM algorithm (Generalized additive models)
#' ==========================================
gammodel <- gam(species ~ s(bio3) + s(bio5) + s(bio6) + s(bio12) + s(bio18) + s(bio20),
                family="binomial", data=traindata)
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


############################################################################################################
############################################################################################################
#' ==========================================
#' Prediction
#' ==========================================
#load datasets
rcp2_2050_list <- list.files(path="D:\\01_Uni\\02_Master\\MET1_Modeling_Prediction\\RCP\\mpi_esm_lr_rcp2_6_2050s_bio_30s_r1i1p1_b4_asc\\bio_b4", pattern = ".asc", full.names = TRUE)
rcp2_2080_list <- list.files(path="D:\\01_Uni\\02_Master\\MET1_Modeling_Prediction\\RCP\\mpi_esm_lr_rcp2_6_2080s_bio_30s_r1i1p1_b4_asc\\bio_b4", pattern = ".asc", full.names = TRUE)
rcp4_2050_list <- list.files(path="D:\\01_Uni\\02_Master\\MET1_Modeling_Prediction\\RCP\\mpi_esm_lr_rcp4_5_2050s_bio_30s_r1i1p1_b4_asc\\bio_b4", pattern = ".asc", full.names = TRUE)
rcp4_2080_list <- list.files(path="D:\\01_Uni\\02_Master\\MET1_Modeling_Prediction\\RCP\\mpi_esm_lr_rcp4_5_2080s_bio_30s_r1i1p1_b4_asc\\bio_b4", pattern = ".asc", full.names = TRUE)
rcp8_2050_list <- list.files(path="D:\\01_Uni\\02_Master\\MET1_Modeling_Prediction\\RCP\\mpi_esm_lr_rcp8_5_2050s_bio_30s_r1i1p1_b4_asc\\bio_b4", pattern = ".asc", full.names = TRUE)
rcp8_2080_list <- list.files(path="D:\\01_Uni\\02_Master\\MET1_Modeling_Prediction\\RCP\\mpi_esm_lr_rcp8_5_2080s_bio_30s_r1i1p1_b4_asc\\bio_b4", pattern = ".asc", full.names = TRUE)

#stack all rasters of one scenario and year
rcp2_2050 <- raster::stack(rcp2_2050_list)
rcp2_2080 <- raster::stack(rcp2_2080_list)
rcp4_2050 <- raster::stack(rcp4_2050_list)
rcp4_2080 <- raster::stack(rcp4_2080_list)
rcp8_2050 <- raster::stack(rcp8_2050_list)
rcp8_2080 <- raster::stack(rcp8_2080_list)

#crop layer to study area extent with a buffer of 3
rcp2_2050_crop <- crop(rcp2_2050, extent(study_area) + 3)
rcp2_2080_crop <- crop(rcp2_2080, extent(study_area) + 3)
rcp4_2050_crop <- crop(rcp4_2050, extent(study_area) + 3)
rcp4_2080_crop <- crop(rcp4_2080, extent(study_area) + 3)
rcp8_2050_crop <- crop(rcp8_2050, extent(study_area) + 3)
rcp8_2080_crop <- crop(rcp8_2080, extent(study_area) + 3)

#rename layer --> must be equal to actual climate data
names(rcp2_2050_crop) <- gsub("_", "", names(rcp2_2050_crop))
names(rcp2_2080_crop) <- gsub("_", "", names(rcp2_2080_crop))
names(rcp4_2050_crop) <- gsub("_", "", names(rcp4_2050_crop))
names(rcp4_2080_crop) <- gsub("_", "", names(rcp4_2080_crop))
names(rcp8_2050_crop) <- gsub("_", "", names(rcp8_2050_crop))
names(rcp8_2080_crop) <- gsub("_", "", names(rcp8_2080_crop))

#resample gtopo dataset to pixelsize of RCP data
gtopo30_res_rcp <- raster::resample(gtopo30_crop, rcp2_2050_crop)

#stack scenarios togehter with elevatin data (gtopo30)
rcp2_2050_crop <- stack(rcp2_2050_crop, gtopo30_res_rcp)
rcp2_2080_crop <- stack(rcp2_2080_crop, gtopo30_res_rcp)
rcp4_2050_crop <- stack(rcp4_2050_crop, gtopo30_res_rcp)
rcp4_2080_crop <- stack(rcp4_2080_crop, gtopo30_res_rcp)
rcp8_2050_crop <- stack(rcp8_2050_crop, gtopo30_res_rcp)
rcp8_2080_crop <- stack(rcp8_2080_crop, gtopo30_res_rcp)

#plot RCP scenarios (BIO20)
plot(rcp2_2050_crop[[20]])
plot(study_area, add=T)


#list_rcp <- as.list(rcp2_2050_list, rcp2_2080_list, rcp4_2050_list, rcp4_2080_list, rcp8_2050_list, rcp8_2080_list)

############################################################################################################
#' ==========================================
#' GAM algorithm for future climate scenarios
#' ==========================================
#' 
gammap_out_folder <- "D:\\01_Uni\\02_Master\\MET1_Modeling_Prediction\\data\\gammap"

#select subset of uncorrelated variables
env_rcp2_2050 <- subset(rcp2_2050_crop, c("bio3", "bio5", "bio6", "bio12", "bio18", "bio20"))
env_rcp2_2080 <- subset(rcp2_2080_crop, c("bio3", "bio5", "bio6", "bio12", "bio18", "bio20"))
env_rcp4_2050 <- subset(rcp4_2050_crop, c("bio3", "bio5", "bio6", "bio12", "bio18", "bio20"))
env_rcp4_2080 <- subset(rcp4_2080_crop, c("bio3", "bio5", "bio6", "bio12", "bio18", "bio20"))
env_rcp8_2050 <- subset(rcp8_2050_crop, c("bio3", "bio5", "bio6", "bio12", "bio18", "bio20"))
env_rcp8_2080 <- subset(rcp8_2080_crop, c("bio3", "bio5", "bio6", "bio12", "bio18", "bio20"))

#use the previous calculates GAM model for future prediction
gammap_rcp2_2050 <- predict(env_rcp2_2050, gammodel, type = "response")
gammap_rcp2_2080 <- predict(env_rcp2_2080, gammodel, type = "response")
gammap_rcp4_2050 <- predict(env_rcp4_2050, gammodel, type = "response")
gammap_rcp4_2080 <- predict(env_rcp4_2080, gammodel, type = "response")
gammap_rcp8_2050 <- predict(env_rcp8_2050, gammodel, type = "response")
gammap_rcp8_2080 <- predict(env_rcp8_2080, gammodel, type = "response")

#plot model results
plot(gammap_rcp4_2050)
plot(study_area, add=T)

#resample/mask current climate to pixel size of future scenarios (!!)
gammap_res <- raster::resample(gammap, gammap_rcp2_2050)
gammap_mask <- raster::mask(gammap_res, study_area)
plot(gammap_mask)
plot(study_area, add=T)

#stack RCP scenario with current climate
gam_rcp2_stack <- stack(gammap_res, gammap_rcp2_2050, gammap_rcp2_2080)
gam_rcp4_stack <- stack(gammap_res, gammap_rcp4_2050, gammap_rcp4_2080)
gam_rcp8_stack <- stack(gammap_res, gammap_rcp8_2050, gammap_rcp8_2080)

#mask
gam_rcp2_stack <- raster::mask(gam_rcp2_stack, study_area)
gam_rcp4_stack <- raster::mask(gam_rcp4_stack, study_area)
gam_rcp8_stack <- raster::mask(gam_rcp8_stack, study_area)

#write gammaps scenarios (stacked)
writeRaster(gammap_mask, filename = paste0(gammap_out_folder, "/model_current.tif"), options="INTERLEAVE=BAND", overwrite=T)
writeRaster(gam_rcp2_stack, filename = paste0(gammap_out_folder, "/RGB_RCP2_full.tif"), options="INTERLEAVE=BAND", overwrite=T)
writeRaster(gam_rcp4_stack, filename = paste0(gammap_out_folder, "/RGB_RCP4_full.tif"), options="INTERLEAVE=BAND", overwrite=T)
writeRaster(gam_rcp8_stack, filename = paste0(gammap_out_folder, "/RGB_RCP8_full.tif"), options="INTERLEAVE=BAND", overwrite=T)


############################################################################################################
#' ==========================================
#' defining threshold for data vizualising
#' ==========================================
#creating binary map with treshold 0.7 -> 70% probability of growing caffea arabica
class_out_folder <- "D:\\01_Uni\\02_Master\\MET1_Modeling_Prediction\\data\\class"
binaryMap <- function(raster, threshold) {
  bin <- raster
  bin[raster <= threshold] <- 0
  bin[raster > threshold] <- 1
  return(bin)
}

threshold <- 0.7

#brick GAM models
ras_current <- brick(gammap)
ras_rcp2_2050 <- brick(gammap_rcp2_2050)
ras_rcp2_2080 <- brick(gammap_rcp2_2080)
ras_rcp4_2050 <- brick(gammap_rcp4_2050)
ras_rcp4_2080 <- brick(gammap_rcp4_2080)
ras_rcp8_2050 <- brick(gammap_rcp8_2050)
ras_rcp8_2080 <- brick(gammap_rcp8_2080)

#create scenarios with threshold
class_current <- binaryMap(ras_current, threshold)
class_rcp2_2050 <- binaryMap(ras_rcp2_2050, threshold)
class_rcp2_2080 <- binaryMap(ras_rcp2_2080, threshold)
class_rcp4_2050 <- binaryMap(ras_rcp4_2050, threshold)
class_rcp4_2080 <- binaryMap(ras_rcp4_2080, threshold)
class_rcp8_2050 <- binaryMap(ras_rcp8_2050, threshold)
class_rcp8_2080 <- binaryMap(ras_rcp8_2080, threshold)

#projection (?)
proj4string(class_current) <- CRS("+proj=longlat +datum=WGS84")
proj4string(class_rcp2_2050) <- CRS("+proj=longlat +datum=WGS84")
proj4string(class_rcp2_2080) <- CRS("+proj=longlat +datum=WGS84")
proj4string(class_rcp4_2050) <- CRS("+proj=longlat +datum=WGS84")
proj4string(class_rcp4_2080) <- CRS("+proj=longlat +datum=WGS84")
proj4string(class_rcp8_2050) <- CRS("+proj=longlat +datum=WGS84")
proj4string(class_rcp8_2080) <- CRS("+proj=longlat +datum=WGS84")

#resample current climate to pixel size of future scenarios (!!)
class_current <- raster::resample(class_current, class_rcp2_2050)
plot(class_current)

#stack RCP scenario with current climate
class_rcp2_stack <- stack(class_current, class_rcp2_2050, class_rcp2_2080)
class_rcp4_stack <- stack(class_current, class_rcp4_2050, class_rcp4_2080)
class_rcp8_stack <- stack(class_current, class_rcp8_2050, class_rcp8_2080)

#mask
class_rcp2_stack <- raster::mask(class_rcp2_stack, study_area)
class_rcp4_stack <- raster::mask(class_rcp4_stack, study_area)
class_rcp8_stack <- raster::mask(class_rcp8_stack, study_area)

#write raster
writeRaster(class_rcp2_stack, filename = paste0(class_out_folder, "/RGB_RCP2_class.tif"), options="INTERLEAVE=BAND", overwrite=T)
writeRaster(class_rcp4_stack, filename = paste0(class_out_folder, "/RGB_RCP4_class.tif"), options="INTERLEAVE=BAND", overwrite=T)
writeRaster(class_rcp8_stack, filename = paste0(class_out_folder, "/RGB_RCP8_class.tif"), options="INTERLEAVE=BAND", overwrite=T)


plot(class_rcp2_2050)
plot(study_area, add=T)


############################################################################################################
#' ==========================================
#' plotting results for each scenario as ggRGB
#' ==========================================
library(RStoolbox)

#plot of each RCP scenario with full data range
ggRGB(gam_rcp2_stack, r = 1, g = 2, b = 3)+
  #xlab("test")+
  ggtitle("Prediction in Caffea Arabica for RCP2")+
  theme_void()  # Empty theme without axis lines and texts

ggRGB(gam_rcp4_stack, r = 1, g = 2, b = 3)+
  #xlab("test")+
  ggtitle("Prediction in Caffea Arabica for RCP4")+
  theme_void()  # Empty theme without axis lines and texts

ggRGB(gam_rcp8_stack, r = 1, g = 2, b = 3)+
  #xlab("test")+
  ggtitle("Prediction in Caffea Arabica for RCP8")+
  theme_void()  # Empty theme without axis lines and texts

#-------------------
#plot of each scenario with threshold 0.7
ggRGB(class_rcp2_stack, r = 3, g = 2, b = 1)+
  #xlab("test")+
  ggtitle("Prediction in Caffea Arabica for RCP2, threshold 0.7")+
  theme_void()  # Empty theme without axis lines and texts

ggRGB(class_rcp4_stack, r = 3, g = 2, b = 1)+
  #xlab("test")+
  ggtitle("Prediction in Caffea Arabica for RCP4, threshold 0.7")+
  theme_void()  # Empty theme without axis lines and texts

ggRGB(class_rcp8_stack, r = 3, g = 2, b = 1)+
  #xlab("test")+
  ggtitle("Prediction in Caffea Arabica for RCP8, threshold 0.7")+
  theme_void()  # Empty theme without axis lines and texts

############################################################################################################
#' ==========================================
#' calculate area of binary results
#' ==========================================

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

xres(class_rcp2_2050)
area_current <- calc_area(class_current)
area_rcp2_2050 <- calc_area(class_rcp2_2050)  #pixel 0.92km
area_rcp2_2080 <- calc_area(class_rcp2_2080)
area_rcp4_2050 <- calc_area(class_rcp4_2050)
area_rcp4_2080 <- calc_area(class_rcp4_2080)
area_rcp8_2050 <- calc_area(class_rcp8_2050)
area_rcp8_2080 <- calc_area(class_rcp8_2080)



#------------------------------------------------------------------------------------------
#install.packages("SSDM")
library(SSDM)
#install.packages("Rmpi")
library(Rmpi)

stack_env <- raster::stack(env)
esdm_traindata <- as(traindata, "data.frame")
esdm_traindata$species <- factor(esdm_traindata$species)


ESDM <- ensemble_modelling(algorithms=c('RF', 'GAM', 'ANN', 'MAXENT'), Occurrences = esdm_traindata, Pcol='species', Env = stack_env,
                           rep = 1, Xcol = 'lon', Ycol = 'lat', method = "pSSDM", verbose = FALSE,
                           endemism = "WEI")

ESDM <- ensemble_modelling(algorithms=c('RF', 'GAM', 'ANN', 'MAXENT'), Occurrences = esdm_traindata, Pcol='species', Env = stack_env,
                           rep = 1, Xcol = 'lon', Ycol = 'lat')
plot(ESDM)
plot(study_area, add=T)

