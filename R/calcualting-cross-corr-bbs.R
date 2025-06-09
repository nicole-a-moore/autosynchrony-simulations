## measuring pair-wise cross-correlation in temperature between each bbs route 
library(tidyverse)
library(parallel)
library(foreach)
library(doParallel)
library(terra)
library(raster)
library(sf)
theme_set(theme_bw())

## set up parallel processing:
starts <- rep(100, 40)
fx <- function(nstart) kmeans(Boston, 4, nstart=nstart)
numCores <- detectCores() ## detect cores
numCores

registerDoParallel(numCores)  ## use multicore, set to the number of our cores

## get filenames
files <- readRDS("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/tas_filenames.rds")

## keep only files 1960 and onwards
files = files[which(str_detect(files, "1960")):length(files)]

## make date vector
dates <- paste(rep(seq(1960, 2021), each = 365), ".", rep(seq(1:365), 62), sep = "")

## get index of first day of 1966
index = which(as.character(dates) == "1966.1")

dates = dates[index:length(dates)]

temp = rast(files[1])
file = 2
while(file <= length(files)) {
  temp = c(temp, rast(files[file]))
  print(file)
  file = file + 1
}

## crop to dates past 1966
be = temp[[index:dim(temp)[3]]]
dim(be)[3] == length(dates) ## make sure we have the right number of days 
rm("temp")

names(be) = dates
saveRDS(be, "outputs/data-processed/be.rds")

be <- readRDS("outputs/data-processed/be.rds")

## crop to North America
ext(be) <- c(-180, 180, -90, 90)
na = mask(be, countries, inverse = FALSE, updatevalue = NA, touches = TRUE)
plot(na[[1]])
rm("be")

na = crop(na, c(-180,-35,20,84))

## save
writeRaster(na, "outputs/data-processed/BerkeleyEarth_1966_NorthAmerica.tif", overwrite = TRUE)

## now loop through chunks and calculate anomaly 
md = str_split_fixed(dates, "\\.", 2)[,2]
md <- unique(md)

na <- stack("outputs/data-processed/BerkeleyEarth_1966_NorthAmerica.tif")

list = list()
i = 1
while(i <= length(md)) {
  layers <- which(str_split_fixed(dates, "\\.", 2)[,2] == md[i])
  list[[i]] <- mean(na[[layers]], na.rm = TRUE)
  print(i)
  i = i + 1
}
climatology_stack <- stack(list)
names(climatology_stack) <- md

saveRDS(climatology_stack, "outputs/data-processed/climatology_stack.rds")

climatology_stack <- readRDS("outputs/data-processed/climatology_stack.rds")


## read in bbs route data
bbs <- read.csv("outputs/data-processed/BBS_all.csv")

## filter to unique routes
bbs <- bbs %>%
  mutate(route = paste(Route, CountryNum, StateNum, sep = "_")) %>%
  select(route, Longitude, Latitude) %>%
  unique()

## figure out which cells each route falls within 
bbs$cell = cellFromXY(na, xy = as.matrix(bbs[,2:3]))

cells = unique(bbs$cell)
length(cells) ## 1438 unique cells

## get pairwise combination of each
pw <- expand.grid(cell1_old = cells, cell2_old = cells) %>%
  filter(cell1_old != cell2_old) %>%
  mutate(cell1 = ifelse(cell1_old > cell2_old, cell1_old, cell2_old),
         cell2 = ifelse(cell2_old > cell1_old, cell1_old, cell2_old)) %>%
  select(cell1, cell2) %>%
  distinct()

## and calculate pearson correlation for pairwise combination of bbs routes 
log_file <- "outputs/data-processed/log.txt"
file.create(log_file)

result = c()
i = 1
result = foreach(i = 1:nrow(pw), .combine = 'rbind') %dopar% {
  write(paste("Processing i =", i), file = log_file, append = TRUE)
  
  ## get time series at each cell
  ext = terra::extract(na, c(pw$cell1[i], pw$cell2[i]))
  
  ## get seasonal climatologies and subtract
  climat = terra::extract(climatology_stack, c(pw$cell1[i], pw$cell2[i]))
  for(x in 1:55) {
    climat = cbind(climat, climat[,1:365])
  }
  ext = ext - climat
  
  pc = cor(as.numeric(ext[1,]), as.numeric(ext[2,]), method = 'pearson')
  
  ## return the pearson corr
  return(c(pc))
}
result = as.data.frame(result)
