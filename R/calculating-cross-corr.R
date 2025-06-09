## testing out ways to measure synchrony

## take an environmental grid as an example
filename = "outputs/data-processed/env-grids/shift_grid90_p0_beta0_r1.2_K100_d0.1_icp0.1_L1500_reps100.tif"
env = rast(filename)


## measure: avg. cross correlation between each time series and it's 8 neighbours
## (should be 1)

## first-difference the time series values
diffs = diff(be)
df_diffs = as.data.frame(diffs, xy = TRUE)

## and calculate pearson correlation for each neighbourhood of 8 cells
pearson_corr <- function(diffs, df_diffs, focal) {
  neighbours <- as.data.frame(adjacent(diffs, focal, directions = "8", pairs = TRUE, include = FALSE, symmetrical = FALSE))
  pc = c()
  for(i in 1:length(neighbours$to)) {
    cur = neighbours$to[i]
    pc = append(pc, cor(as.numeric(df_diffs[focal,c(3:ncol(df_diffs))]), as.numeric(df_diffs[cur,c(3:ncol(df_diffs))]), method = 'pearson'))
  }
  neighbours$pearson_cor = pc

  # Return the average value of the surrounding 8 cells
  return(c(focal, mean(neighbours$pearson_cor)))
}

result = c()
for(i in 1:nrow(df_diffs)) {
  result = rbind(result, pearson_corr(diffs, df_diffs, focal = i))
}
result = as.data.frame(result)

## rasterize result
result$x = df_diffs$x
result$y = df_diffs$y

result <- dplyr::select(result, c(x, y, V2))

cors = rast(result, type = "xyz")

plot(cors)


mean(values(cors))

## next: try on Berkeley Earth data across bbs sampling period (1966 - 2023)
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

## calculate synchrony of each cell with its adjacent cells
## first-difference the time series values
diffs = diff(be)
rm("be")
saveRDS(diffs, "outputs/data-processed/diffs.rds")
diffs = readRDS("outputs/data-processed/diffs.rds")

## keep only cells in North America 
countries <- necountries::countries(c("Canada", "United States of America"), part = TRUE)
countries <- st_union(countries)
countries <- vect(countries)

ext(diffs) <- c(-180, 180, -90, 90)
na = mask(diffs, countries, inverse = FALSE, updatevalue = NA, touches = TRUE)
plot(na[[1]])
diffs = na

## and calculate pearson correlation for each neighbourhood of 8 cells
log_file <- "outputs/data-processed/log.txt"
file.create(log_file)

non_na = which(!is.na(diffs[,,1]))
result = c()
i = 1
result = foreach(i = 1:length(non_na), .combine = 'rbind') %dopar% {
  write(paste("Processing i =", i), file = log_file, append = TRUE)
  
  ## get cell neighbours
  focal = non_na[i]
  neighbours <- as.data.frame(adjacent(diffs, focal, directions = "8", pairs = TRUE, include = FALSE, symmetrical = FALSE))
  
  ## extract ts 
  ext = terra::extract(diffs, c(focal, neighbours$to))
  
  pc = c()
  for(t in 1:length(neighbours$to)) {
    pc = append(pc, cor(as.numeric(ext[1,]), as.numeric(ext[t,]), method = 'pearson'))
  }
  neighbours$pearson_cor = pc

  # Return the average value of the surrounding 8 cells
 return(c(focal, mean(neighbours$pearson_cor, na.rm = TRUE)))
}
result = as.data.frame(result)

## rasterize result
cors = diffs[[1]]
cors[result$V1] <- result$V2

plot(cors)
  
names(cors) = "mean_cross_corr"
cors %>%
  ggplot(aes(fill = mean_cross_corr)) +
  geom_spatraster(data = cors) +
  scale_fill_viridis_c(na.value = "transparent") +
  scale_x_continuous(limits = c(-180,-50)) +
  scale_y_continuous(limits = c(20,90)) +
  labs(fill = "Avg. cross\ncorrelation")

## save 
writeRaster(cors, "outputs/data-processed/BerkeleyEarth_cross-correlation.tif", 
              overwrite = TRUE)


## now do the same for seasonally-detrended temperatures 
## first, create raster of seasonally detrended temperatures 
## read in Berkeley Earth data 
be <- readRDS("outputs/data-processed/be.rds")

## crop to North America
ext(be) <- c(-180, 180, -90, 90)
na = mask(be, countries, inverse = FALSE, updatevalue = NA, touches = TRUE)
plot(na[[1]])
rm("be")

na = crop(na, c(-180,-35,20,84))

## save
writeRaster(na, "outputs/data-processed/BerkeleyEarth_1966_NorthAmerica.tif", overwrite = TRUE)

# ## create output directory for block-wise climatology rasters
# dir.create("outputs/data-processed/temp_climatology", showWarnings = FALSE)
# 
# ## save in spatial chunks 
# row = nrow(be)/3 
# col = ncol(be)/3 
# for(i in 1:3) {
#   c = 1
#   while(c <= 3) {
#     e <- ext(col - ncol(be)/3, col,
#              row - 150, row - 90)
#     temp <- crop(be, e)
#     
#     writeRaster(temp, paste0("outputs/data-processed/temp_climatology/be-chunk_x", 
#                              e[1],"-",e[2],"_y", e[3],"_",e[4],".tif") , overwrite = TRUE)
#     row = row + nrow(be)/3
#     print(paste0("c: ",c))
#     c = c + 1
#   }
#   row = nrow(be)/3
#   col = col + ncol(be)/3
#   
#   print(paste0("i: ", i))
# }
# 
# rm("be")

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

## seasonally detrend by subtracting climatology from each day
# i = 1
# list = list()
# while(i <= dim(na)[3]/365) {
#   # ((i-1)*365+1):(i*365)
#   list <- append(list, na[[((i-1)*365+1):(i*365)]] - climatology_stack)
#   print(i)
#   i = i + 1
# }



# files = list.files("outputs/data-processed/temp_climatology", full.names = T)
# file = 1
# while(file <= length(files)) {
#   ## read rast
#   r = rast(files[file])
#   
#   ## aggregate raster to daily climatology
#   climatology <- lapply(md, function(d) {
#     layers <- which(str_split_fixed(dates, "\\.", 2)[,2] == d)
#     mean(r[[layers]], na.rm = TRUE)
#   })
#   climatology_stack <- rast(climatology)
#   names(climatology_stack) <- md
#   
#   ## seasonally detrend by subtracting climatology from each day
#   r = r - rep(climatology_stack, dim(r)[3]/365)
#   
#   ## save 
#   writeRaster(r, paste0(str_split_fixed(files[file], ".tif", 2)[1,1], "_s-detrended",".tif") , overwrite = TRUE)
#   
#   print(file)
#   file = file + 1
# }
# 
# ## recombine detrended time series 
# files = list.files("outputs/data-processed/temp_climatology", full.names = T)
# files = files[which(str_detect(files, "s-det"))]
# 
# file = 1
# while(file <= length(files)) {
#   ## read rast
#   r = stack(files[file])
#   
#   ## combine with others
#   if(file == 1) {
#     r_all = r
#   }
#   else {
#     r_all = mosaic(r_all, r, fun = "mean")
#   }
#   
#   print(file)
#   file = file + 1
# }
# ## save
# writeRaster(r_all, "outputs/data-processed/be_s-detrended.tif", overwrite = TRUE)
# 
# r_all<- rast("outputs/data-processed/be_s-detrended.tif")

## calculate synchrony of each cell with its adjacent cells

## and calculate pearson correlation for each neighbourhood of 8 cells
log_file <- "outputs/data-processed/log.txt"
file.create(log_file)

non_na = which(!is.na(na[,,1]))
result = c()
i = 1
result = foreach(i = 1:length(non_na), .combine = 'rbind') %dopar% {
  write(paste("Processing i =", i), file = log_file, append = TRUE)
  
  ## get cell neighbours
  focal = non_na[i]
  neighbours <- as.data.frame(adjacent(na, focal, directions = "8", pairs = TRUE, 
                                       include = FALSE, symmetrical = FALSE))
  
  ## extract ts 
  ext = terra::extract(na, c(focal, neighbours$to))
  
  ## get seasonal climatologies and subtract
  climat = terra::extract(climatology_stack, c(focal, neighbours$to))
  for(x in 1:55) {
    print(x)
    climat = cbind(climat, climat[,1:365])
  }
  ext = ext - climat
  
  pc = c()
  for(t in 1:length(neighbours$to)) {
    pc = append(pc, cor(as.numeric(ext[1,]), as.numeric(ext[t,]), method = 'pearson'))
  }
  neighbours$pearson_cor = pc
  
  # Return the average value of the surrounding 8 cells
  return(c(focal, mean(neighbours$pearson_cor, na.rm = TRUE)))
}
result = as.data.frame(result)

## rasterize result
cors = na[[1]]
cors[result$V1] <- result$V2

plot(cors)

names(cors) = "mean_cross_corr"
cors %>%
  ggplot(aes(fill = mean_cross_corr)) +
  geom_spatraster(data = cors) +
  scale_fill_viridis_c(na.value = "transparent") +
  scale_x_continuous(limits = c(-180,-50)) +
  scale_y_continuous(limits = c(20,90)) +
  labs(fill = "Avg. cross\ncorrelation")

## save 
writeRaster(cors, "outputs/data-processed/BerkeleyEarth_cross-correlation_seasonally-detrended.tif", 
            overwrite = TRUE)



