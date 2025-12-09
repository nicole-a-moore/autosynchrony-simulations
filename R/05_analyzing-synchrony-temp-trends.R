## analyzing how much variation in warming trends there is across space from 1966-2023
library(terra)
library(raster)

## read in BE data across North America from 1966-now
## as raster:
na <- stack("outputs/data-processed/BerkeleyEarth_1966_NorthAmerica.tif")
plot(na[[1]])

## as 3d array:
array <- readRDS("outputs/data-processed/BE_all_array.rds")

## loop through each cell, fit linear regression to time series, and calculate & save slope
slopes = array[,,1]
x = 1
while(x <= dim(array)[2]) {
  
  y = 1
  while(y <= dim(array)[1]) {
    
    ## get local time series 
    ts <- data.frame(time = 1:(length(array[1,1,])), 
                     temp = array[y,x,])
    if(!any(is.nan(ts$temp))) {
      ## run linear regression for grid cell
      lm_out <- lm(ts, formula = temp ~ time)
      
      ## get slope of regression and save 
      slopes[y,x] <- lm_out$coefficients[2]
    }
    else {
      slopes[y,x] = NA
    }
    
    y = y + 1
  }
  x = x + 1
}

hist(slopes)

max(slopes, na.rm = T)
## 0.0001711523 degrees per day
## 0.06247059 degrees per year 

## plot raster showing variation in slopes across space
r = raster(slopes)
plot(r*365*10)

 ## interesting - some variation but not a lot 
# 0.2 to 0.6deg C per 10 years, but lots of spatial autocorrelation in warming trends 


## maybe need to look at the spatial autocorrelation of anomalies from a baseline period?
## make date vector from 1966 - 1980
dates <- paste(rep(seq(1966, 1979), each = 365), ".", rep(seq(1:365), 14), sep = "")

## now loop through chunks and calculate anomaly 
md = str_split_fixed(dates, "\\.", 2)[,2]
md <- unique(md)

clim = array[,,1]
x = 1
while(x <= dim(clim)[2]) {
  y = 1
  while(y <= dim(clim)[1]) {
    
    i = 1
    while(i <= length(md)) {
      layers <- which(str_split_fixed(dates, "\\.", 2)[,2] == md[i])
      clim[y,x] <- mean(array[y,x,layers], na.rm = TRUE)
      i = i + 1
    }
    
    y = y + 1
  }
  print(x)
  x = x + 1
}

plot(rast(clim))

saveRDS(clim, "outputs/data-processed/clim_array_1966-1979.rds")

climatology_stack <- readRDS("outputs/data-processed/climatology_stack.rds")
clim_array = as.array(climatology_stack)

## repeat climatologies for each year 
for(x in 1:56) {
  if(x == 1) {
    clim_array_all = clim_array
  }
  else {
    clim_array_all = abind::abind(clim_array_all, clim_array)
  }
}
dim(clim_array_all)[3] == dim(array)[3]

## subtract climatology from each temp
anom = array - clim_array_all

plot(raster(clim_array_all[,,1]))
plot(raster(anom[,,1]))
plot(raster(array[,,1]))
plot(raster(clim_array_all[,,1] + anom[,,1]))

## look at anomalies over time (e.g., the first year)
library(terra)

ggplot() +
  geom_spatraster(data = rast(anom[,,190])) +
  scale_fill_gradient2(high = "red", low = "blue", mid = "white", midpoint = 0,
                       limits = c(-20,20))




