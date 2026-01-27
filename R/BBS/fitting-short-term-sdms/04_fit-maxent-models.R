## try recreating a MaxEnt models from Bateman et al.
library(data.table)
library(tidyverse)
library(terra)
select = dplyr::select

## increase amount of memory
options(java.parameters = "-Xmx32g" )
library(dismo)

## read in background samples
bkgd <- read.csv("/Volumes/NIKKI/background_bioclim.csv")
nrow(bkgd)

### make a prediction
################################################### 
## read in all environmental layers 
avgtmp_6m <- rast("outputs/data-processed/weather-vars/avgtmp_6m.tiff")
avgtmp_12m <- rast("outputs/data-processed/weather-vars/avgtmp_12m.tiff")
avgtmp_36m <- rast("outputs/data-processed/weather-vars/avgtmp_36m.tiff")
tmpseas_6m <- rast("outputs/data-processed/weather-vars/tmpseas_6m.tiff")
tmpseas_12m <- rast("outputs/data-processed/weather-vars/tmpseas_12m.tiff")
tmpseas_36m <- rast("outputs/data-processed/weather-vars/tmpseas_36m.tiff")
tmx_6m <- rast("outputs/data-processed/weather-vars/tmx_6m.tiff")
tmx_12m<-rast("outputs/data-processed/weather-vars/tmx_12m.tiff")
tmx_36m<-rast("outputs/data-processed/weather-vars/tmx_36m.tiff")
tmn_6m<-rast("outputs/data-processed/weather-vars/tmn_6m.tiff")
tmn_12m<-rast("outputs/data-processed/weather-vars/tmn_12m.tiff")
tmn_36m<-rast("outputs/data-processed/weather-vars/tmn_36m.tiff")
annualpre_6m <- rast("outputs/data-processed/weather-vars/annualpre_6m.tiff")
annualpre_12m <- rast("outputs/data-processed/weather-vars/annualpre_12m.tiff")
annualpre_36m <- rast("outputs/data-processed/weather-vars/annualpre_36m.tiff")
prewet_12m <- rast("outputs/data-processed/weather-vars/prewet_12m.tiff")
prewet_36m <- rast("outputs/data-processed/weather-vars/prewet_36m.tiff")
predry_12m <- rast("outputs/data-processed/weather-vars/predry_12m.tiff")
predry_36m<- rast("outputs/data-processed/weather-vars/predry_36m.tiff")
preseas_6m <- rast("outputs/data-processed/weather-vars/preseas_6m.tiff")
preseas_12m <- rast("outputs/data-processed/weather-vars/preseas_12m.tiff")
preseas_36m <- rast("outputs/data-processed/weather-vars/preseas_36m.tiff")

list = list(avgtmp_6m, avgtmp_12m, avgtmp_36m,
            tmpseas_6m, tmpseas_12m, tmpseas_36m,
            tmx_6m, tmx_12m, tmx_36m,
            tmn_6m, tmn_12m, tmn_36m,
            annualpre_6m, annualpre_12m, annualpre_36m,
            prewet_12m, prewet_36m,
            predry_12m, predry_36m,
            preseas_6m, preseas_12m, preseas_36m)
names(list) = list("avgtmp_6m", "avgtmp_12m", "avgtmp_36m",
                   "tmpseas_6m", "tmpseas_12m", "tmpseas_36m",
                   "tmx_6m", "tmx_12m", "tmx_36m",
                   "tmn_6m", "tmn_12m", "tmn_36m",
                   "annualpre_6m", "annualpre_12m", "annualpre_36m",
                   "prewet_12m", "prewet_36m",
                   "predry_12m", "predry_36m",
                   "preseas_6m", "preseas_12m", "preseas_36m")

## get index of 1966
names = str_split_fixed(names(avgtmp_6m), "_",3)[,3]
y1 = first(which(substr(names, 1, 4) == 1966))

## apply to each clim variable 
i = 1
while(i <= length(list)) {
  
  ## read in layer 
  lyr <- list[[i]]
  
  ## subset to data after 1966
  lyr = lyr[[y1:nlyr(lyr)]]
  
  list[[i]] <- lyr
  
  ## make stack to use for prediction for each month x year combination
  names(lyr)
  
  i = i + 1
}


## make a raster stack for each year x month 
## variables must be in order that they are in background data 
order = colnames(bkgd)[5:26]
list = list[order]
list[[1]]

## number of layers
n <- nlyr(list[[1]])

## create list to hold the new stacks
stacks <- vector("list", n)

for(i in 1:n) {
  ## extract layer i from each raster and stack them
  layers_i <- lapply(list, \(x) x[[i]])
  stacks[[i]] <- rast(layers_i)
}

## name them by date 
names(stacks) <- paste0("date_", str_split_fixed(names(list[[1]]), "_", 3)[,3])

first = stacks[[1]]
first
plot(first[[1]])

plot(annualpre_12m[[y1]])

## save:
for (i in 1:n) {
  writeRaster(stacks[[i]], paste0("outputs/data-processed/weather-vars/", names(stacks)[i], ".tif"), overwrite = TRUE)
}

saveRDS(stacks, "outputs/data-processed/weather-vars/stacks.rds")

##################
#    FIT SDMs    #
##################
## read in weather data 

## get occurrence files 
occ_files <- list.files("/Volumes/NIKKI/gbif-occurrences_filtered", full.names = T)

length(occ_files)

## get species list 
species = str_split_fixed(occ_files, "gbif-occurrences_filtered/", 2)[,2]
species = str_split_fixed(species, "_distinct", 2)[,1]

## for each species 
sp = 1
sp = which(occ_files == "/Volumes/NIKKI/gbif-occurrences_filtered/Pyrocephalus_rubinus_distinct.csv")
while(sp <= length(occ_files)) {
  
  ## read in occurrences + background data
  occ <- read.csv(occ_files[sp])
  
  ## make a folder to save everything in
  path = paste0("/Volumes/NIKKI/sdms/", species[sp])
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
  
  set.seed(123)
  
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
    path = paste0("/Volumes/NIKKI/sdms/", species[sp], "/predictions/")
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

## check which species don't have sdms:

sdms = list.files("/Volumes/NIKKI/sdms/")

missing = species[which(!species %in% sdms)]



## compare to Bateman et al. predictions 
## read in their data
files = list.files("data-raw/short-term-sdms/Short-Term_Climate_Variability_Models_10/Passerella_iliaca/outputs", 
                   full = T)
files = files[which(str_detect(files, "\\.asc"))]

## check that all can be read properly
no_errors = c()
for(i in 1:length(files)) {
  tryCatch({
    ## attempt to read file and access data
    raster_data <- rast(files[i])
    
    summary(values(raster_data))
    
    no_errors = append(no_errors, i)
  }, error = function(e) {
    print(paste("Error reading file:", i))
  })
}

## get rid of files with errors 
files = files[no_errors]

## read in rasters 
all = rast(files)

## read in our data
all_ours = rast("/Volumes/NIKKI/sdms/Passerella_iliaca/predictions/Passerella_iliaca_predictions_1966-2024.tif")

names(all)
names(all_ours)

## change extent
all_ours = crop(all_ours, all)
all = crop(all, all_ours)

all = resample(all, all_ours)

mask = ifel(is.na(all[[1]]), NA, 1)
plot(mask)

all_ours = mask(all_ours, mask)

which(names(all) == "198207")
which(names(all_ours) == "date_1982_7")

plot(all[[132]]) 
plot(all_ours[[68]] > 0.3)

ggplot() +
  geom_spatraster(data = all[[132]]) +
  scale_fill_continuous(limits = c(0,1))

ggplot() +
  geom_spatraster(data = all_ours[[168]]) +
  scale_fill_continuous(limits = c(0,1))



plot(all[[67]])
plot(all_ours[[4]])

plot(all[[67]] - all_ours[[4]])


## plot bbs data on top



plot(pred[[60]])



