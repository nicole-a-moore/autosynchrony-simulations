## download and filter GBIF occurrences for breeding bird survey birds 
options(warn = 1)
library(tidyverse)
library(rgbif)
library(data.table)
library(ggplot2)
library(tidyterra)
library(terra)
library(CoordinateCleaner)
library(countrycode)
theme_set(theme_bw())
extract = terra::extract

## start with: Passerella_iliaca

## wipe old occurrences 
file.remove(list.files("outputs/data-processed/gbif-occurrences/", full.names = TRUE))


####################################
#    DOWNLOAD GBIF OCCURRENCES     #
####################################
## download occurrences for Passerella_iliaca
## from 1951 - 2024 (to match temperature data I have downloaded)

## create list of landbird species from bbs to download occurrences for
## for now: just ones that passed the spatial filter
bbs <- read.csv("outputs/data-processed/BBS_spatial-filter.csv")
sp = unique(bbs$genus_sp)
rm("bbs")

## get backbone for species of interest
backbones = lapply(sp, FUN = name_backbone)
names(backbones) <- sp

## make occurrence request for each species
dkeys <- c()
i = 1
while(i <= length(backbones)) {
  backbone = backbones[[i]]
  
  request <- occ_download(
    pred_in("speciesKey", backbone$speciesKey),
    pred_in("basisOfRecord", c("HUMAN_OBSERVATION", "OBSERVATION", "OCCURRENCE")), # basis of records for downloading GBIF data
    ## filter to 1951-06-16 to 2024-12-16
    pred_gte("year", 1951),
    pred_lte("year", 2024),
    ## filter to April, May, June, July
    pred_gte("month", 4), 
    pred_lte("month", 7), 
    pred("hasCoordinate", TRUE), # must have coordinates
    pred("hasGeospatialIssue", FALSE), # and no geospatial issues
    pred("occurrenceStatus", "PRESENT"), 
    format = "SIMPLE_CSV",
    user = Sys.getenv('GBIF_USER'), 
    pwd = Sys.getenv('GBIF_PWD'),
    email = "nicole.moore@mail.mcgill.ca")
  
  occ_download_wait(request)
  
  # save information about request
  tmp = request
  tmp.save <- attributes(tmp)
  tmp.save <- data.frame(download_key = tmp[1],
                         created = tmp.save$created,
                         download_link = tmp.save$downloadLink,
                         citation = tmp.save$citation,
                         format = tmp.save$format,
                         user = tmp.save$user,
                         email = tmp.save$email)
  request <- tmp.save
  
  ## write out info about request
  write.csv(request, 
            paste0("outputs/data-processed/gbif-requests/gbif_request_", 
                   sub(" ", "-", backbone$canonicalName),
                   ".csv"), row.names = F)
  
  occ_download_get(request$download_key, 
                   path = "outputs/data-processed/gbif-occurrences",
                   overwrite = TRUE)
  
  ## make list of download keys
  dkeys <- append(dkeys, request$download_key)
  
  print(paste("finished species no. ", i))
  i = i + 1
}

## add species name to download keys and save 
dkeys = data.frame(download_keys = dkeys, canonicalName = sp)
write.csv(dkeys, "outputs/data-processed/gbif-requests/gbif-download-keys.csv", row.names = F)


####################################
#     FILTER GBIF OCCURRENCES      #
####################################
## read in download list
dkeys = read.csv("outputs/data-processed/gbif-requests/gbif-download-keys.csv")

## read in climate data mask
mask <- rast("outputs/data-processed/weather-vars/mask.tiff")

## loop through files
lessthan30 <- c()
f = 1
while(f <= nrow(dkeys)) {
  
  ## unzip and read in raw occurrences 
  unzip(paste0("/Volumes/NIKKI/gbif-occurrences/", dkeys$download_keys[f], ".zip"), exdir = "/Volumes/NIKKI/gbif-occurrences/") # unzip file
  occ <- fread(paste0("/Volumes/NIKKI/gbif-occurrences/", dkeys$download_keys[f], ".csv"))
  
  ## 1. use coordinate cleaner to flag issues with occurrences 
  #############################################################
  ## convert country code from ISO2c to ISO3c
  occ$countryCode <-  countrycode(occ$countryCode, 
                                  origin =  'iso2c',
                                  destination = 'iso3c')
  
  ## convert to dataframe
  occ = as.data.frame(occ)
  
  flags <- clean_coordinates(x = occ, 
                             lon = "decimalLongitude", 
                             lat = "decimalLatitude",
                             species = "species",
                             tests = c("capitals", "centroids"),
                             capitals_rad = 1000,
                             centroids_rad = 1000) 
  
  ## exclude flagged data from occurrences
  occ <- occ[flags$.summary,]
  
  ## 2. filter to locations with CRU TS data
  ###########################################
  occ <- select(occ, decimalLongitude, decimalLatitude, species, year, month) %>%
    mutate(gbif_lon = decimalLongitude, gbif_lat = decimalLatitude)
  occ_vect = vect(occ, geom = c("decimalLongitude", "decimalLatitude"))
  crs(occ_vect) = crs(mask)
  
  ## remove occurrences outside of mask 
  vals <- extract(mask, occ_vect)
  occ_vect <- occ_vect[!is.na(vals[,2]), ]
  
  ## 3. get cell coordinates within CRU TS data
  #############################################
  ## get coordinates of raster cells of occurrences
  cell_ids <- cellFromXY(mask, geom(occ_vect)[,c(3,4)])
  cell_xy <- xyFromCell(mask, cell_ids)
  occ_vect$x <- cell_xy[,1]
  occ_vect$y <- cell_xy[,2]
  
  occ = as.data.frame(occ_vect)
  
  ## save only occurrences in distinct months x years x cells 
  occ_distinct = select(occ, -gbif_lon, -gbif_lat) %>%
    distinct()
  
  ## 4. filter to species with more than 30 occurrences 
  ######################################################
  if(nrow(occ_distinct) > 30) {
    
    ## 5. filter to species that Nikki and Maxence love
    ######################################################
    
    ## 6. save
    #######################################
    ## write out 
    # write.csv(occ, paste0("/Volumes/NIKKI/gbif-occurrences_filtered/", 
    #                                str_replace_all(dkeys$canonicalName[f], "\\ ", "\\_"), ".csv"), row.names = F)
    write.csv(occ_distinct, paste0("/Volumes/NIKKI/gbif-occurrences_filtered/", 
                                  str_replace_all(dkeys$canonicalName[f], "\\ ", "\\_"), "_distinct.csv"), row.names = F)
  }
  else {
    lessthan30 <- append(lessthan30, dkeys$canonicalName[f])
  }
  
  ## remove unzipped file 
  file.remove(paste0("/Volumes/NIKKI/gbif-occurrences/", dkeys$download_keys[f], ".csv"))
  
  f = f + 1
}



## how many distinct occurrences? 
nrow(occ_distinct) 

## what dates are they from? 
occ_distinct %>%
  ggplot(aes(x = year)) +
  geom_bar()

## where are they?
# download data with geodata's world function to use for our base map
world_map <- geodata::world(resolution = 3,
                   path = "data/")

# crop to area of interest
# store boundaries in a single extent object
geographic_extent <- ext(x = c(min(occ$gbif_lon) - 3, 
                               max(occ$gbif_lon) + 3, 
                               min(occ$gbif_lat) - 3, 
                               max(occ$gbif_lat) + 3))

map <- crop(x = world_map, y = geographic_extent)

map %>%
  ggplot() +
  geom_spatvector() +
  geom_point(data = occ,
             aes(x = gbif_lon, y = gbif_lat),
             size = 0.1) +
  facet_wrap(~month)

map %>%
  ggplot() +
  geom_spatvector() +
  geom_point(data = occ_distinct,
             aes(x = x, y = y),
             size = 0.1) +
  facet_wrap(~month)





####################################
#       MAKE BACKGROUND LAYER      #
####################################
## target-group background
## make background layer of all unique lat x lon x year x month combinations where any of the species in the sample 




