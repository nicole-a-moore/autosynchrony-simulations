## extract values of seasonal bioclimatic variables at GBIF occurrences 
library(tidyverse)
library(terra)
extract = terra::extract

## match date and year of occurrence to proper month of bioclimatic variables:

## data format output:
## species, lon, lat, month, bioclim_6m, bioclim_xx, etc.

## prepare background data (combine all species-level data into one dataframe and keep distinct rows)

## prep occurrences 
################################################
occ_df <- read.csv("outputs/data-processed/gbif-occurrences/occ_df.csv")

## prep bioclim layers  
################################################
files = list.files("outputs/data-processed/weather-vars", full.names = T)
files = files[-which(str_detect(files, "mask"))]

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


## join to occurrence data   
################################################
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
  
  ## left join bioclim data to occurrences 
  occ <- left_join(occ, bioclim, by = c("year", "month", "x", "y"))
  
  ## save
  write.csv(occ, paste0("/Volumes/NIKKI/occurrences-x-bioclim/", 
                        str_replace_all(unique(occ$species), "\\ ", "\\_"),
                        "_bioclim.csv"), row.names = F)
  print(i)
  i = i + 1
}

## add climate data to background points:
bkgd <- left_join(bkgd, bioclim, by = c("year", "month", "x", "y"))

## write:
write.csv(bkgd, "/Volumes/NIKKI/occurrences-x-bioclim/background_bioclim.csv", row.names = F)




