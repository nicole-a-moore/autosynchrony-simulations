## calculate 22 seasonal variables with time lags from CRU TS
## following Bateman et al. 
library(data.table)
library(tidyverse)
library(terra)
select = dplyr::select

## STEPS:
###################################
## 1. download monthly CRU TS data
#      tmp - mean monthly temperature (°C)
#      tmn - minimum monthly temperature (°C)
#      tmx - maximum monthly temperature (°C)
#      pre - monthly precipitation (mm)
######################################################################
## 2. calculate the following 8 bioclimatic variables:
#      mean monthly temperature (°C)
#      temperature seasonality (standard deviation × 100)
#      maximum temperature of the warmest month (°C)
#      minimum temperature of the coldest month (°C)
#      annual precipitation (mm)
#      precipitation of the wettest quarter (mm)
#      precipitation of the driest quarter (mm)
#      precipitation seasonality (coefficient of variation, CV).
## for each 6 month, 12 month, and 36 month period proceeding each month 
## ** precip of driest and wetest quarter not calculated for 6 month period
############################################################################
 

dates = paste(rep(1951:2024, each = 12), rep(1:12, n = 888/12), sep = "_")
breeding_months = which(str_split_fixed(dates, "\\_", 2)[,2] %in% c(4,5,6,7))

## make shapefile for countries
countries <- necountries::countries(c("Canada", "United States of America"), part = TRUE)
countries <- sf::st_union(countries)
countries <- vect(countries)


##############################
##           TMP            ##
##############################
## read in tmp layers
####################################################################
tmp = list.files("data-raw/CRU_TS/tmp", full.names = T)
tmp = tmp[order(tmp)]
tmp = rast(tmp, subds = "tmp")
nlyr(tmp) ## 888

## mean monthly temperature (°C)
###################################
avgtmp_6m <- roll(tmp, n = 6, fun = "mean", circular = FALSE, type = "to")
avgtmp_12m <- roll(tmp, n = 12, fun = "mean", circular = FALSE, type = "to")
avgtmp_36m <- roll(tmp, n = 36, fun = "mean", circular = FALSE, type = "to")

## temperature seasonality (standard deviation × 100)
#####################################################
tmpseas_6m <- roll(tmp, n = 6, fun = "sd", circular = FALSE, type = "to")*100
tmpseas_12m <- roll(tmp, n = 12, fun = "sd", circular = FALSE, type = "to")*100
tmpseas_36m <- roll(tmp, n = 36, fun = "sd", circular = FALSE, type = "to")*100

## change layer names to include new variable name, date
names(avgtmp_6m) = paste0("avgtmp_6m_", dates)
names(avgtmp_12m) = paste0("avgtmp_12m_", dates)
names(avgtmp_36m) = paste0("avgtmp_36m_", dates)
names(tmpseas_6m) = paste0("tmpseas_6m_", dates)
names(tmpseas_12m) = paste0("tmpseas_12m_", dates)
names(tmpseas_36m) = paste0("tmpseas_36m_", dates)

## get rid of layers for non-breeding season months (4,5,6,7)
avgtmp_6m = avgtmp_6m[[breeding_months]]
avgtmp_12m = avgtmp_12m[[breeding_months]]
avgtmp_36m = avgtmp_36m[[breeding_months]]
tmpseas_6m = tmpseas_6m[[breeding_months]]
tmpseas_12m = tmpseas_12m[[breeding_months]]
tmpseas_36m = tmpseas_36m[[breeding_months]]

## crop to north america
avgtmp_6m = mask(avgtmp_6m, countries)
avgtmp_12m = mask(avgtmp_12m, countries)
avgtmp_36m = mask(avgtmp_36m, countries)
tmpseas_6m = mask(tmpseas_6m, countries)
tmpseas_12m = mask(tmpseas_12m, countries)
tmpseas_36m = mask(tmpseas_36m, countries)

## save tmp layers
######################
writeRaster(avgtmp_6m, "outputs/data-processed/weather-vars/avgtmp_6m.tiff", overwrite = T)
writeRaster(avgtmp_12m, "outputs/data-processed/weather-vars/avgtmp_12m.tiff", overwrite = T)
writeRaster(avgtmp_36m, "outputs/data-processed/weather-vars/avgtmp_36m.tiff", overwrite = T)
writeRaster(tmpseas_6m, "outputs/data-processed/weather-vars/tmpseas_6m.tiff", overwrite = T)
writeRaster(tmpseas_12m, "outputs/data-processed/weather-vars/tmpseas_12m.tiff", overwrite = T)
writeRaster(tmpseas_36m, "outputs/data-processed/weather-vars/tmpseas_36m.tiff", overwrite = T)

## remove tmp objects from memory:
rm("tmp", "avgtmp_6m", "avgtmp_12m", "avgtmp_36m", "tmpseas_6m", "tmpseas_12m", "tmpseas_36m")

##############################
##       TMX  &  TMN        ##
##############################
## read in tmx layers
####################################################################
tmx = list.files("data-raw/CRU_TS/tmx", full.names = T)
tmx = tmx[order(tmx)]
tmx = rast(tmx, subds = "tmx")
nlyr(tmx) ## 888

## maximum temperature of the warmest month (°C)
#####################################################
tmx_6m <- roll(tmx, n = 6, fun = "max", circular = FALSE, type = "to")
tmx_12m <- roll(tmx, n = 12, fun = "max", circular = FALSE, type = "to")
tmx_36m <- roll(tmx, n = 36, fun = "max", circular = FALSE, type = "to")

## read in tmn layers
####################################################################
tmn = list.files("data-raw/CRU_TS/tmn", full.names = T)
tmn = tmn[order(tmn)]
tmn = rast(tmn, subds = "tmn")
nlyr(tmn) ## 888

## minimum temperature of the coldest month (°C)
#####################################################
tmn_6m <- roll(tmn, n = 6, fun = "min", circular = FALSE, type = "to")
tmn_12m <- roll(tmn, n = 12, fun = "min", circular = FALSE, type = "to")
tmn_36m <- roll(tmn, n = 36, fun = "min", circular = FALSE, type = "to")

## change layer names to include new variable name, date
names(tmx_6m) = paste0("tmx_6m_", dates)
names(tmx_12m) = paste0("tmx_12m_", dates)
names(tmx_36m) = paste0("tmx_36m_", dates)
names(tmn_6m) = paste0("tmn_6m_", dates)
names(tmn_12m) = paste0("tmn_12m_", dates)
names(tmn_36m) = paste0("tmn_36m_", dates)

## get rid of layers for non-breeding season months (4,5,6,7)
tmn_6m = tmn_6m[[breeding_months]]
tmn_12m = tmn_12m[[breeding_months]]
tmn_36m = tmn_36m[[breeding_months]]
tmx_6m = tmx_6m[[breeding_months]]
tmx_12m = tmx_12m[[breeding_months]]
tmx_36m = tmx_36m[[breeding_months]]

## crop to north america
tmn_6m = mask(tmn_6m, countries)
tmn_12m = mask(tmn_12m, countries)
tmn_36m = mask(tmn_36m, countries)
tmx_6m = mask(tmx_6m, countries)
tmx_12m = mask(tmx_12m, countries)
tmx_36m = mask(tmx_36m, countries)

## save tmx and tmn layers
######################
writeRaster(tmx_6m, "outputs/data-processed/weather-vars/tmx_6m.tiff", overwrite = T)
writeRaster(tmx_12m, "outputs/data-processed/weather-vars/tmx_12m.tiff", overwrite = T)
writeRaster(tmx_36m, "outputs/data-processed/weather-vars/tmx_36m.tiff", overwrite = T)
writeRaster(tmn_6m, "outputs/data-processed/weather-vars/tmn_6m.tiff", overwrite = T)
writeRaster(tmn_12m, "outputs/data-processed/weather-vars/tmn_12m.tiff", overwrite = T)
writeRaster(tmn_36m, "outputs/data-processed/weather-vars/tmn_36m.tiff", overwrite = T)

## remove tmp objects from memory:
rm("tmx", "tmn", "tmn_6m", "tmn_12m", "tmn_36m", "tmx_6m", "tmx_12m","tmx_36m")

##############################
##            PRE           ##
##############################
## read in pre layers
####################################################################
pre = list.files("data-raw/CRU_TS/pre", full.names = T)
pre = pre[order(pre)]
pre = rast(pre, subds = "pre")
nlyr(pre) ## 888

## annual precipitation (mm)
#####################################################
annualpre_6m <- roll(pre, n = 6, fun = "sum", circular = FALSE, type = "to")
annualpre_12m <- roll(pre, n = 12, fun = "sum", circular = FALSE, type = "to")
annualpre_36m <- roll(pre, n = 36, fun = "sum", circular = FALSE, type = "to")

## change layer names to include new variable name, date
names(annualpre_6m) = paste0("annualpre_6m_", dates)
names(annualpre_12m) = paste0("annualpre_12m_", dates)
names(annualpre_36m) = paste0("annualpre_36m_", dates)

## get rid of layers for non-breeding season months (4,5,6,7)
annualpre_6m = annualpre_6m[[breeding_months]]
annualpre_12m = annualpre_12m[[breeding_months]]
annualpre_36m = annualpre_36m[[breeding_months]]

## crop to north america
annualpre_6m = mask(annualpre_6m, countries)
annualpre_12m = mask(annualpre_12m, countries)
annualpre_36m = mask(annualpre_36m, countries)

## save annual pre layers
###########################
writeRaster(annualpre_6m, "outputs/data-processed/weather-vars/annualpre_6m.tiff", overwrite = T)
writeRaster(annualpre_12m, "outputs/data-processed/weather-vars/annualpre_12m.tiff", overwrite = T)
writeRaster(annualpre_36m, "outputs/data-processed/weather-vars/annualpre_36m.tiff", overwrite = T)

## remove tmp objects from memory:
rm("annualpre_6m", "annualpre_12m", "annualpre_36m")

## precipitation of the wettest quarter (mm)
#####################################################
## calculate precipitation over quarters 
n_months <- nlyr(pre)
n_years <- n_months / 12
quarters <- lapply(1:n_years, function(i){
 
 start <- (i - 1) * 12
 
 q1 <- sum(pre[[start + 1:3]])
 q2 <- sum(pre[[start + 4:6]])
 q3 <- sum(pre[[start + 7:9]])
 q4 <- sum(pre[[start + 10:12]])
 
 c(q1, q2, q3, q4)
})

## covert list to spat rast 
pre_quarters <- do.call(c, quarters)

## calculate wetest quarter
prewet_12m <- roll(pre_quarters, n = 4, fun = "max", circular = FALSE, type = "to")
prewet_36m <- roll(pre_quarters, n = 4, fun = "max", circular = FALSE, type = "to")

## repeat each by 3 
prewet_12m = rep(prewet_12m, each = 3)
prewet_36m = rep(prewet_36m, each = 3)

## precipitation of the driest quarter (mm)
#####################################################
predry_12m <- roll(pre_quarters, n = 4, fun = "min", circular = FALSE, type = "to")
predry_36m <- roll(pre_quarters, n = 4, fun = "min", circular = FALSE, type = "to")

## repeat each by 3 
predry_12m = rep(predry_12m, each = 3)
predry_36m = rep(predry_36m, each = 3)

## change layer names to include new variable name, date
names(prewet_12m) = paste0("prewet_12m_", dates)
names(prewet_36m) = paste0("prewet_36m_", dates)
names(predry_12m) = paste0("predry_12m_", dates)
names(predry_36m) = paste0("predry_36m_", dates)

## get rid of layers for non-breeding season months (4,5,6,7)
prewet_12m = prewet_12m[[breeding_months]]
prewet_36m = prewet_36m[[breeding_months]]
predry_12m = predry_12m[[breeding_months]]
predry_36m = predry_36m[[breeding_months]]

## crop to north america
prewet_12m = mask(prewet_12m, countries)
prewet_36m = mask(prewet_36m, countries)
predry_12m = mask(predry_12m, countries)
predry_36m = mask(predry_36m, countries)

## save quarterly pre layers
###########################
writeRaster(prewet_12m, "outputs/data-processed/weather-vars/prewet_12m.tiff", overwrite = T)
writeRaster(prewet_36m, "outputs/data-processed/weather-vars/prewet_36m.tiff", overwrite = T)
writeRaster(predry_12m, "outputs/data-processed/weather-vars/predry_12m.tiff", overwrite = T)
writeRaster(predry_36m, "outputs/data-processed/weather-vars/predry_36m.tiff", overwrite = T)

## remove tmp objects from memory:
rm("prewet_12m", "prewet_36m", "predry_12m", "predry_36m", "pre_quarters")


## precipitation seasonality (coefficient of variation, CV).
#####################################################
presd_6m <- roll(pre, n = 6, fun = "sd", circular = FALSE, type = "to")
presd_12m <- roll(pre, n = 12, fun = "sd", circular = FALSE, type = "to")
presd_36m <- roll(pre, n = 36, fun = "sd", circular = FALSE, type = "to")

premean_6m <- roll(pre, n = 6, fun = "mean", circular = FALSE, type = "to")
premean_12m <- roll(pre, n = 12, fun = "mean", circular = FALSE, type = "to")
premean_36m <- roll(pre, n = 36, fun = "mean", circular = FALSE, type = "to")

preseas_6m<- (presd_6m/premean_6m)*100
preseas_12m <- (presd_12m/premean_12m)*100
preseas_36m <- (presd_36m/premean_36m)*100

## change layer names to include new variable name, date
names(preseas_6m) = paste0("preseas_6m_", dates)
names(preseas_12m) = paste0("preseas_12m_", dates)
names(preseas_36m) = paste0("preseas_36m_", dates)

## get rid of layers for non-breeding season months (4,5,6,7)
preseas_6m = preseas_6m[[breeding_months]]
preseas_12m = preseas_12m[[breeding_months]]
preseas_36m = preseas_36m[[breeding_months]]

## crop to north america
preseas_6m = mask(preseas_6m, countries)
preseas_12m = mask(preseas_12m, countries)
preseas_36m = mask(preseas_36m, countries)

## save seasonal pre layers
###########################
writeRaster(preseas_6m, "outputs/data-processed/weather-vars/preseas_6m.tiff", overwrite = T)
writeRaster(preseas_12m, "outputs/data-processed/weather-vars/preseas_12m.tiff", overwrite = T)
writeRaster(preseas_36m, "outputs/data-processed/weather-vars/preseas_36m.tiff", overwrite = T)

## remove tmp objects from memory:
rm("pre", "preseas_6m", "preseas_12m", "preseas_36m", "presd_6m", "presd_12m", "presd_36m",
   "premean_6m", "premean_12m", "premean_36m")


## make a data mask:
avgtmp_6m <- rast("outputs/data-processed/weather-vars/avgtmp_6m.tiff")
plot(avgtmp_6m[[6]])

mask = ifel(is.na(avgtmp_6m[[6]]), NA, 1)
plot(mask)
writeRaster(mask, "outputs/data-processed/weather-vars/mask.tiff", overwrite = T)







### garbage:
#######################

## check
# plot(annualpre_36m[[1]])
# plot(sum(pre[[1:36]]))

## save all the rasters
r_list = list(avgtmp_6m, avgtmp_12m, avgtmp_36m,
              tmpseas_6m, tmpseas_12m, tmpseas_36m,
              tmx_6m, tmx_12m, tmx_36m,
              tmn_6m, tmn_12m, tmn_36m,
              annualpre_6m, annualpre_12m, annualpre_36m,
              prewet_12m, prewet_36m,
              predry_12m, predry_36m,
              preseas_6m, preseas_12m, preseas_36m)

saveRDS(r_list, "outputs/data-processed/weather-vars/r_list.rds")

## get rid of NA layers 
avgtmp_6m = avgtmp_6m[[6:nlyr(avgtmp_6m)]]
avgtmp_12m = avgtmp_12m[[12:nlyr(avgtmp_12m)]]
avgtmp_36m = avgtmp_36m[[36:nlyr(avgtmp_36m)]]
tmpseas_6m = avgtmp_6m[[6:nlyr(tmpseas_6m)]]
tmpseas_12m = avgtmp_12m[[12:nlyr(tmpseas_12m)]]
tmpseas_36m = avgtmp_36m[[36:nlyr(tmpseas_36m)]]
tmx_6m = avgtmp_6m[[6:nlyr(tmx_6m)]]
tmx_12m = avgtmp_12m[[12:nlyr(tmx_12m)]]
tmx_36m = avgtmp_36m[[36:nlyr(tmx_36m)]]
tmn_6m = avgtmp_6m[[6:nlyr(tmn_6m)]]
tmn_12m = avgtmp_12m[[12:nlyr(tmn_12m)]]
tmn_36m = avgtmp_36m[[36:nlyr(tmn_36m)]]
annualpre_6m = annualpre_6m[[6:nlyr(annualpre_6m)]]
annualpre_12m = annualpre_12m[[12:nlyr(annualpre_12m)]]
annualpre_36m = annualpre_36m[[36:nlyr(annualpre_36m)]]
prewet_12m = prewet_12m[[12:nlyr(prewet_12m)]]
prewet_36m = annualpre_36m[[36:nlyr(prewet_36m)]]
predry_12m = prewet_12m[[12:nlyr(predry_12m)]]
predry_36m = annualpre_36m[[36:nlyr(predry_36m)]]
preseas_6m = preseas_6m[[6:nlyr(preseas_6m)]]
preseas_12m = preseas_12m[[12:nlyr(preseas_12m)]]
preseas_36m = annualpre_36m[[36:nlyr(preseas_36m)]]

