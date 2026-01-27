.libPaths( c( "~/projects/def-jsunday/nikkim/sdms/packages" , .libPaths() ))

## increase amount of memory
options(java.parameters = "-Xmx32g" )

## for cluster 
library(tidyverse)
library(data.table)
library(dismo)
library(predicts)
library(terra)
library(stringr)
set.seed(123)
select = dplyr::select

## read in background samples
bkgd <- read.csv("background_bioclim.csv")
nrow(bkgd)

## read in weather variables 
files = list.files("weather-data", full = T)
files = files[str_detect(files, "date")]
stacks = list(length = length(files))
for(i in 1:length(files)) {
  stacks[[i]] = rast(files[i])
}
names(stacks) = str_split_fixed(files, "weather-data/", 2)[,2]
names(stacks) = str_split_fixed(names(stacks), ".tif", 2)[,1]

##################
#    FIT SDMs    #
##################
## get occurrence files 
occ_files <- list.files("occ-data", full.names = T)
length(occ_files)

## get species list 
species = str_split_fixed(occ_files, "occ-data/", 2)[,2]
species = str_split_fixed(species, "_distinct", 2)[,1]

## for each species 
sp = 226
while(sp <= length(occ_files)) {
  ## read in occurrences + background data
  occ <- read.csv(occ_files[sp])
  
  ## make a folder to save everything in
  path = paste0("outputs/sdms/", species[sp])
  if(!file.exists(path)){
    dir.create(path)
  }
  
  ## mark the presences in the bkgd data 
  occ = select(occ, -species)
  occ$presence = 1
  
  data = left_join(bkgd, occ)
  
  data$presence = ifelse(is.na(data$presence),0,data$presence)
  
  data <- select(data, presence, x, y, everything())
  
  ## split predictors / response
  weather <- data[, 6:ncol(data)]
  presences <- data$presence
  
  ## fit maxent model for 10-fold cross validation
  path_cv = paste0(path, "/cross-validation")
  maxent_mod_cv <- maxent(x = weather, p = presences, 
                          path = path_cv,
                          #factors = "month",
                          args = c("nowarnings",
                                   "noprefixes",
                                   "responsecurves",
                                   "noaskoverwrite",
                                   "noremoveduplicates",
                                   "nothreshold",
                                   "nohinge",
                                   "replicates=10",
                                   "replicatetype=crossvalidate"))
  
  ## get AUC
  AUC = mean(maxent_mod_cv@results[grep("Test.AUC", rownames(maxent_mod_cv@results)), ])
  
  ## if mean AUC > 0.5 
  if(AUC > 0.5) {
    ## fit maxent model with all data
    maxent_mod <- maxent(x = weather, p = presences, 
                         path = path,
                         #factors = "month",
                         args = c("nowarnings",
                                  "noprefixes",
                                  "responsecurves",
                                  "noaskoverwrite",
                                  "noremoveduplicates",
                                  "nothreshold",
                                  "nohinge"))
    
    
    
    ## make predictions for each month from 1966 - 2023
    ## make a folder to save everything in
    path = paste0("outputs/sdms/", species[sp], "/predictions/")
    if(!file.exists(path)){
      dir.create(path)
    }
    for(i in 1:length(stacks)) {
      if(i == 1) {
        pred = predict(stacks[[i]], maxent_mod, na.rm = T)
      }
      else {
        pred = c(pred, predict(stacks[[i]], maxent_mod, na.rm = T))
      }
      print(i)
    }
    names(pred) = names(stacks)
    writeRaster(pred, paste0(path, species[sp], "_predictions_1966-2024.tif"), overwrite = TRUE)
  }
  sp = sp + 1
}
  
