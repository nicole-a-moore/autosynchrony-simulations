## combine all bioclim layers into a data frame & extract values at all coordinates of the background layer 
library(tidyverse)
library(terra)
extract = terra::extract

#########################################
##         PREP BIOCLIM LAYERS         ##
#########################################
files = list.files("outputs/data-processed/weather-vars", full.names = T)
files = files[-which(str_detect(files, "mask"))]
files = files[-which(str_detect(files, ".csv"))]
files = files[-which(str_detect(files, "date"))]

i = 1
while(i <= length(files)) {
  
  ## read in layer 
  lyr <- rast(files[[i]])
  
  ## covert to xyz dataframe 
  lyr = as.data.frame(lyr, xy = TRUE)
  
  ## get rid of NA columns
  lyr = lyr[, colSums(is.na(lyr)) != nrow(lyr)]
  
  ## gather and create columns for each month/year
  lyr_name =  paste(str_split_fixed(names(lyr)[3], "\\_", 4)[,1], str_split_fixed(names(lyr)[3], "\\_", 4)[,2], sep = "_")
  lyr = lyr %>%
    gather(key = "layer", value = lyr_name, colnames(.[,3:ncol(.)])) %>%
    mutate(year = as.integer(str_split_fixed(layer, "\\_", 4)[,3]),
           month = as.integer(str_split_fixed(layer, "\\_", 4)[,4])) %>%
    select(-layer)
  colnames(lyr)[which(colnames(lyr) == "lyr_name")] = lyr_name
  
  ## join to other bioclim data
  if(i == 1) {
    bioclim = lyr
  }
  else {
   bioclim = left_join(bioclim, lyr, by = join_by("x", "y", "year", "month"))
  }
  
  i = i + 1
  print(i)
}

## save 
write.csv(bioclim, "outputs/data-processed/weather-vars/all_vars.csv", row.names = F)


#########################################
##         MAKE BACKGROUND LAYER       ##
#########################################
## get names of occurrence files 
occ_files = list.files("/Volumes/NIKKI/gbif-occurrences_filtered", full.names = T)

## read in all bioclim data 
bioclim <- read.csv("outputs/data-processed/weather-vars/all_vars.csv")

i = 1
bkgd <- c()
while(i <= length(occ_files)) {
  
  ## read in occurrences 
  occ <- read.csv(occ_files[i])
  
  ## make background layer 
  bkgd <- rbind(bkgd, select(occ, -species)) %>%
    distinct()
  
  
  print(i)
  i = i + 1
}

## add climate data to background points:
bkgd <- left_join(bkgd, bioclim, by = c("year", "month", "x", "y"))

## write:
write.csv(bkgd, "/Volumes/NIKKI/background_bioclim.csv", row.names = F)




