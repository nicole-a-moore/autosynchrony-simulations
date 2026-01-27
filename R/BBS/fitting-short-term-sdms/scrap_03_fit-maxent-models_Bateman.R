## try recreating a MaxEnt models from Bateman et al.
library(data.table)
library(tidyverse)
select = dplyr::select

## increase amount of memory
options(java.parameters = "-Xmx32g" )
library(dismo)

##################
#    FIT SDM     #
##################
## read in occurrences + climate data at occurrences for Passerella iliaca
occ <- read.csv("data-raw/short-term-sdms/Short-Term_Climate_Variability_Models_10/Passerella_iliaca/weather.occur.csv")

## read in background samples
bkgd <- read.csv("data-raw/short-term-sdms/weather.bkgd.csv")

nrow(bkgd)
nrow(occ)

## check they are the same env data
occ[100,]
bkgd[which(bkgd$lat == occ$lat[100] &bkgd$lon == occ$lon[100]),]

## mark the presences in the bkgd data 
occ = select(occ, -species)
bkgd = select(bkgd, -species)

occ$presence = 1

data = left_join(bkgd, occ)

data$presence = ifelse(is.na(data$presence),0,data$presence)

## check to make sure it worked
length(which(data$presence == 1)) == nrow(occ)
nrow(data) == nrow(bkgd)

data <- select(data, presence, lon, lat, everything())

## idea: use all presences to train maxent model for training and verify with bbs data 
# Regularization values: linear/quadratic/product: 0.050, categorical: 0.250, threshold: 1.000, hinge: 0.500
# Feature types used: product linear quadratic
# responsecurves: true
# outputdirectory: outputs
# samplesfile: weather.occur.csv
# environmentallayers: weather.bkgd.csv
# warnings: false
# askoverwrite: false
# removeduplicates: false
# threshold: false
# hinge: false
# autorun: true
# visible: false

## fit the model
weather = data[,c(5:ncol(data))]
presences = data[,c(1)]
#weather$month = as.factor(weather$month)
  
maxent_mod <- maxent(x = weather, p = presences, 
                     path = "outputs/data-processed/maxent-results/",
                     #factors = "month",
                     args = c("nowarnings",
                              "noprefixes",
                              "responsecurves",
                              "noaskoverwrite",
                              "noremoveduplicates",
                              "nothreshold",
                              "nohinge"))

## predict 
pred <- predict(maxent_mod, bioclim_data)  

weather 




## plot predictions - suitability
suit = map %>%
  ggplot() +
  geom_spatvector() +
  geom_spatraster(data = pred) +
  scale_fill_viridis_c(na.value = "transparent") +
  geom_point(data = occ, aes(x = decimalLongitude, y = decimalLatitude),
             size = 0.1, alpha = 0.5) +
  labs(title = dkeys$canonicalName[index[i]], 
       fill = "Suitability", x = "Longitude", y = "Latitude")

## save 
ggsave(suit, path = "figures", filename = paste0("suitability_", dkeys$canonicalName[1], ".png"),
       width = 6, height = 6)

## plot predictions - p/a
pa_pred <- pred > threshold

pa_plot = map %>%
  ggplot() +
  geom_spatvector() +
  geom_spatraster(data = pa_pred) +
  scale_fill_viridis_c(na.value = "transparent") +
  geom_point(data = occ, aes(x = decimalLongitude, y = decimalLatitude),
             size = 0.1, alpha = 0.5) +
  labs(title = dkeys$canonicalName[index[i]], 
       fill = "Pr/ab", x = "Longitude", y = "Latitude")

## save 
ggsave(pa_plot, path = "figures", filename = paste0("presabs_", dkeys$canonicalName[1], ".png"),
       width = 6, height = 6)


