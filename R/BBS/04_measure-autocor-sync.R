## measure synchrony and autocorrelation for each bb in analysis 
library(tidyverse)
library(terra)
library(sf)
select = dplyr::select

## read in clean bb data subset 
bb = read.csv("outputs/data-processed/bbs_clean-subset.csv")

## get canada and us map for plotting
countries <- necountries::countries(c("Canada", "United States of America"), part = TRUE)
countries <- st_union(countries)

#######################################
##          autocorrelation          ##
#######################################
## read in analysis of Berkeley Earth data
se <- readRDS("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/BerkeleyEarth_noise-colour.rds")

## select data on low spectral exponent measured across entire time period
se =  se %>%
  select(lat, lon, s_spec_exp_PSD_low_all) %>%
  unique()

## create raster from data frame 
se$lon = se$lon - 180
r = se %>%
  select(lon, lat, s_spec_exp_PSD_low_all) %>%
  rast()

plot(r)

## subset to unique locations and make into sf 
routes = bb %>%
  select(route, Longitude, Latitude) %>%
  distinct() %>%
  st_as_sf(., coords = c("Longitude", "Latitude"))

## extract temporal autocorrelation at each unique route in the bbs data 
st_crs(routes) = st_crs(countries)
beta <- terra::extract(r, routes)
routes$beta = beta$s_spec_exp_PSD_low_all

routes %>% 
  ggplot(aes(colour = beta)) +
  geom_sf(inherit.aes = FALSE, data = countries) +
  geom_sf(size = 0.5) +
  scale_x_continuous(limits = c(-130,-60)) +
  scale_y_continuous(limits = c(25,60)) +
  scale_colour_viridis_c() +
  labs(colour = "")

min(routes$beta, na.rm = T) ## 0.3
max(routes$beta, na.rm = T) ## 1.4

## join to bb data 
routes = cbind(routes, st_coordinates(routes)) %>%
  st_drop_geometry() %>%
  rename("Latitude" = Y, "Longitude" = X)

bb <- left_join(bb, routes)

## turn into sf
bb <- st_as_sf(bb, coords = c("Longitude", "Latitude"))
st_crs(bb) = st_crs(countries)

## for each species, make a plot of the autocorrelation across its bbs routes 
## and calculate mean, max, min

sp = 1
while(sp <= length(unique(bb$AOU))) {
  
  ## filter to data on one species 
  sp_sub = filter(bb, AOU == unique(bb$AOU)[sp]) %>%
    filter(!is.na(TotalAbd)) %>%
    ## get rid of routes where the species was never seen
    group_by(route) %>%
    mutate(sum_abd = sum(TotalAbd)) %>%
    ungroup() %>%
    filter(sum_abd >= 1) 
  
  stats = sp_sub %>%
    summarize(mean_beta = round(mean(beta, na.rm = T), digits = 2), 
              min_beta =  round(min(beta, na.rm = T), digits = 2), 
              max_beta = round(max(beta, na.rm = T), digits = 2)) %>%
    st_drop_geometry()
  
  plot = sp_sub %>% 
    ggplot(aes(colour = beta)) +
    geom_sf(inherit.aes = FALSE, data = countries) +
    geom_sf(size = 0.1) +
    scale_x_continuous(limits = c(-165,-55)) +
    scale_y_continuous(limits = c(25,70)) +
    scale_colour_viridis_c(limits = c(0.2, 1.5)) +
    labs(colour = "", title = unique(sp_sub$genus_sp),
         subtitle = paste0("mean: ", stats$mean_beta, " range: ", stats$min_beta, " - ", stats$max_beta)) 
  
  ## save plot
  ggsave(plot, path = "outputs/figures/autocorrelation", 
         filename = paste0(unique(sp_sub$genus_sp), ".png"), 
         width = 6, height = 4)
  
  sp = sp + 1
}

bb_nogeom = st_drop_geometry(bb)

## calculate stats for all and plot
stats = bb_nogeom %>%
  filter(!is.na(TotalAbd), !TotalAbd == 0) %>% ## filter to only routes where each species was seen 
  select(AOU, genus_sp, beta, route) %>% 
  distinct() %>% ## get rid of duplicated routes within species 
  group_by(AOU, genus_sp) %>% ## group by species and calculate mean 
  summarize(mean_beta = round(mean(beta, na.rm = T), digits = 2), 
            min_beta =  round(mean(beta, na.rm = T), digits = 2), 
            max_beta = round(max(beta, na.rm = T), digits = 2)) 

stats %>%
  ggplot(aes(x = mean_beta)) +
  geom_histogram()

## save
write.csv(stats, "outputs/data-processed/bbs-autocorrelation-stats.csv", row.names = F)

#################################
##          synchrony          ##
#################################
## calculate pairwise cross-correlation in temperature across sites during the BBS period  
library(tidyverse)
library(parallel)
library(foreach)
library(doParallel)
library(terra)
library(raster)
library(sf)
theme_set(theme_bw())
select <- dplyr::select

## set up parallel processing:
starts <- rep(100, 40)
fx <- function(nstart) kmeans(Boston, 4, nstart=nstart)
numCores <- detectCores() ## detect cores
numCores

registerDoParallel(numCores)  ## use multicore, set to the number of our cores

## read in clean bb data subset 
bb = read.csv("outputs/data-processed/bbs_clean-subset.csv") 

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
countries <- vect(countries)
ext(be) <- c(-180, 180, -90, 90)
na = mask(be, countries, inverse = FALSE, updatevalue = NA, touches = TRUE)
plot(na[[1]])
rm("be")

na = crop(na, c(-180,-35,20,84))

na = stack(na)

## save
writeRaster(na, "outputs/data-processed/BerkeleyEarth_1966_NorthAmerica.tif", overwrite = TRUE)

na <- stack("outputs/data-processed/BerkeleyEarth_1966_NorthAmerica.tif")

## now loop through chunks and calculate anomaly 
md = str_split_fixed(dates, "\\.", 2)[,2]
md <- unique(md)

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

## get unique routes 
routes = bb %>%
  select(route, Longitude, Latitude) %>%
  distinct()

## figure out which cell each route falls within 
routes$cell = cellFromXY(climatology_stack[[1]], cbind(routes$Longitude, routes$Latitude))

cells = unique(routes$cell)
length(cells) ## 1423 unique cells

## make routes into sf points df
routes = routes %>%
  st_as_sf(., coords = c("Longitude", "Latitude")) 
st_crs(routes) = st_crs(countries)

## for each route, get list of cells within certain radius 
int_all <- c()
for(rt in 1:nrow(routes)) {

  ## make buffer around the cell
  buff = st_buffer(routes[rt,], dist = 5) ## of 5 degree radius
  #plot(buff)
  
  ## figure out which other routes are within that buffer 
  int = st_intersection(buff, routes) 
  
  ## get unique pairs of cells to compare 
  int <- int %>% 
    select(route, cell, route.1, cell.1) %>% 
    st_drop_geometry() %>%
    distinct()
  
  ## add to list 
  int_all <- rbind(int, int_all)
  
  print(rt)
}

int_all <- int_all %>%
  rename("cell1" = cell, "cell2" = cell.1, "route1" = route, "route2" = route.1)

## get each unique pairwise combination of cells to calculate pearson correlation between 
pw <- int_all %>%
  select(cell1, cell2) %>%
  distinct()

## create 3d array of temperatures across north america to make ts extractions faster 
file = 1
while(file <= length(files)) {
  
  ## read in raster
  temp = rast(files[file])
  
  ## fix extent
  ext(temp) <- c(-180, 180, -90, 90)
  #plot(temp[[1]])
  
  ## crop data to north america
  temp = mask(temp, countries, inverse = FALSE, updatevalue = NA, touches = TRUE)
  temp = crop(temp, c(-180,-35,20,84))
  #plot(temp[[1]])
  
  ## convert data to array and store 
  r_array = as.array(temp)
  
  if (file == 1) {
    all_array = r_array
  }
  else {
    all_array = abind::abind(all_array, r_array, along = 3)
  }
  
  print(file)
  file = file + 1
}

## crop to only dates past 1966
all_array = all_array[,,index:dim(all_array)[3]]

## save the array
saveRDS(all_array, "outputs/data-processed/BE_all_array.rds")

all_array <- readRDS("outputs/data-processed/BE_all_array.rds")

## make climatology array
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
dim(clim_array_all)[3] == dim(all_array)[3]

## and calculate pearson correlation for pairwise combination of bbs routes 
# log_file <- "outputs/data-processed/log.txt"
# file.create(log_file)

## get row and column of each cell 
temp = climatology_stack[[1]]
values(temp) = 1
spts <- rasterToPoints(temp)

rownames(all_array) = unique(spts[,2])
colnames(all_array) = unique(spts[,1])

result = c()
i = 1
while(i <= nrow(pw))  {
 # write(paste("Processing i =", i), file = log_file, append = TRUE)
  
  rowcol = xyFromCell(climatology_stack, c(pw$cell1[i])) %>%
    rbind(., xyFromCell(climatology_stack, c(pw$cell2[i])))
  
  ## get time series at each cell
  ts1 = all_array[which(rownames(all_array) == rowcol[1,2]), which(colnames(all_array) == rowcol[1,1]),]  # note: row = y, col = x
  ts2 = all_array[which(rownames(all_array) == rowcol[2,2]), which(colnames(all_array) == rowcol[2,1]),]
  
  ## get seasonal climatologies at each cell 
  clim1 =  clim_array_all[which(rownames(all_array) == rowcol[1,2]), which(colnames(all_array) == rowcol[1,1]),]
  clim2 =  clim_array_all[which(rownames(all_array) == rowcol[2,2]), which(colnames(all_array) == rowcol[2,1]),]
  
  ## and subtract climatology from temp
  ts1 = ts1 - clim1
  ts2 = ts2 - clim2
  
  ## calculate pairwise pearson correlation between detrended time series 
  pc = cor(ts1, ts2, method = 'pearson')
  
  ## return the pearson corr
  result = append(result, pc)
  
  print(round((i/102689)*100, 4))
  
  i = i + 1
}

## add to df
pw$pearson_corr = as.vector(result)

## visualize 
hist(pw$pearson_corr)

## add back route info
pear_cor <- left_join(int_all, pw)

## get rid of routes falling in cells with no temp data 
pear_cor <- pear_cor[which(!is.na(pear_cor$pearson_corr)),]

## save 
write.csv(pear_cor, "outputs/data-processed/synchrony_all.csv", row.names = T)

## for each species, calculate mean cross correlation between temperature at each route where species was seen and raster cells containing routes where the species was seen within a 5 degree radius
# (if species was seen multiple times in a nearby cell, only count once)
# (filter out routes in same cell)

## get list of each species and routes where it was seen

## filter to unique combinations of routes + species 
routes_bb = bb %>%
  filter(!is.na(TotalAbd)) %>%
  ## get rid of routes where the species was never seen
  group_by(AOU, route) %>%
  mutate(sum_abd = sum(TotalAbd)) %>%
  ungroup() %>%
  filter(sum_abd >= 1) %>%
  select(route, AOU, genus_sp) %>%
  distinct() %>%
  rename("route1" = route)

## for each species + route, get list of routes within 5 deg
routes_5deg = left_join(routes_bb, pear_cor, relationship = "many-to-many")

## filter to unique pairs of cells, get rid of focal cell comparison to self, group by focal cell, and calculate mean
routes_5deg = routes_5deg %>%
  select(-route2) %>%
  distinct() %>%
  filter(cell1 != cell2) %>%
  group_by(genus_sp, AOU, route1, cell1) %>%
  summarize(mean_pearson_corr = mean(pearson_corr, na.rm = T)) %>% 
  rename(route = route1)

hist(routes_5deg$mean_pearson_corr)  

## plot range of mean pearson corr per route for each species 
routes_5deg %>%
  ggplot(aes(x = mean_pearson_corr, fill = AOU)) +
  geom_histogram() +
  theme(legend.position = "none")

## plot stats per species
## calculate stats for each species 
stats <- routes_5deg %>%
  group_by(AOU, genus_sp) %>%
  summarize(mean_mean_pearson_corr = mean(mean_pearson_corr),
            p95 = quantile(mean_pearson_corr, c(0.95), na.rm = T),
            p5 = quantile(mean_pearson_corr, c(0.05),  na.rm = T)) 

stats %>%
  arrange(mean_mean_pearson_corr) %>% 
  mutate(AOU = factor(AOU, levels = c(unique(AOU)), ordered = T)) %>%
  ggplot(aes(x = AOU, y = mean_mean_pearson_corr)) +
  geom_pointrange(aes(ymax = p95, ymin = p5)) 

## save stats
write.csv(stats, "outputs/data-processed/bbs-synchrony-stats.csv", row.names = F)
  

## add geometry 
locs = select(bb, route, Latitude, Longitude) %>% 
  group_by(route) %>%
  distinct()

routes_5deg <- left_join(routes_5deg, locs, relationship = "many-to-many")

routes_5deg <- vect(routes_5deg, geom = c("Longitude", "Latitude")) %>%
  st_as_sf()
st_crs(routes_5deg) = st_crs(countries)

## make plots for each species 
sp = 1
while(sp <= length(unique(routes_5deg$AOU))) {
  
  sp_sub = filter(routes_5deg, AOU == unique(routes_5deg$AOU)[sp])
  
  stats = sp_sub %>%
    summarize(mean_p = round(mean(mean_pearson_corr, na.rm = T), digits = 2), 
              min_p =  round(min(mean_pearson_corr, na.rm = T), digits = 2), 
              max_p = round(max(mean_pearson_corr, na.rm = T), digits = 2)) %>%
    st_drop_geometry()
  
  plot = sp_sub %>% 
    ggplot(aes(colour = mean_pearson_corr)) +
    geom_sf(inherit.aes = FALSE, data = countries) +
    geom_sf(size = 0.1) +
    scale_x_continuous(limits = c(-165,-55)) +
    scale_y_continuous(limits = c(25,70)) +
    scale_colour_viridis_c(limits = c(0.6, 1)) +
    labs(colour = "", title = unique(sp_sub$genus_sp),
         subtitle = paste0("mean: ", stats$mean_p, " range: ", stats$min_p, " - ", stats$max_p)) 
  
  ## save plot
  ggsave(plot, path = "outputs/figures/synchrony", 
         filename = paste0(unique(sp_sub$genus_sp), ".png"), 
         width = 6, height = 4)
  
  sp = sp + 1
}



## what if we calculate it just at the range edge?
## try for northern cardinal 
## get mean pearson correlation across routes falling within the 75 percentile of latitude?



## filter to unique combinations of routes + species 
routes_bb = bb %>%
  select(route, AOU, Latitude, Longitude) %>%
  distinct() %>%
  rename("route1" = route)

## for each species + route, get list of routes within 5 deg
routes_5deg = left_join(routes_bb, pear_cor, relationship = "many-to-many")

## filter to routes at each species range edge (75th percentile lat)
bb_p75 <- bb %>%
  select(AOU, route, Latitude, Longitude) %>%
  distinct() %>%
  group_by(AOU) %>%
  summarise(p75 = quantile(Latitude, 0.75)) %>%
  left_join(routes_5deg, .)

bb_p75 <- bb_p75 %>%
  filter(Latitude >= p75) %>% 
  select(-route2) %>%
  distinct() %>%
  filter(cell1 != cell2) %>%
  group_by(AOU, cell1) %>%
  summarize(mean_pearson_corr = mean(pearson_corr, na.rm = T)) 

## plot stats per species
bb_p75 %>%
  group_by(AOU) %>%
  summarize(mean_mean_pearson_corr = mean(mean_pearson_corr),
            p95 = quantile(mean_pearson_corr, c(0.95), na.rm = T),
            p5 = quantile(mean_pearson_corr, c(0.05),  na.rm = T)) %>%
  arrange(mean_mean_pearson_corr) %>% 
  mutate(AOU = factor(AOU, levels = c(unique(AOU)), ordered = T)) %>%
  ggplot(aes(x = AOU, y = mean_mean_pearson_corr)) +
  geom_pointrange(aes(ymax = p95, ymin = p5)) 

#########################################
##          climate velocity           ##
#########################################
library(VoCC)

## read in Berkeley Earth daily air temperature from 1966 onwards
na <- stack("outputs/data-processed/BerkeleyEarth_1966_NorthAmerica.tif")

## read in as array also
na_arr = tiff::readTIFF("outputs/data-processed/BerkeleyEarth_1966_NorthAmerica.tif")

## crop to 1980 onwards 
## make date vector
dates <- paste(rep(seq(1966, 2021), each = 365), ".", rep(seq(1:365), 62), sep = "")
## get index of first day of 1980
index = which(as.character(dates) == "1980.1")

na_arr = na_arr[,,index:dim(na_arr)[3]]


## calculate temporal trend and mean temp in each cell 
n = dim(na_arr)[3]
temp_trends = means = na_arr[,,1]
x = 1
while(x <= dim(na_arr)[1]) {
  
  y = 1
  while(y <= dim(na_arr)[2]) {
    
    if(!is.nan(na_arr[x,y,1])) {
      ## get time series
      df = data.frame(temp = na_arr[x,y,], time = 1:n)
      
      ## fit regression
      lm = lm(temp ~ time, data = df)
      
      ## save slope 
      temp_trends[x,y] = lm$coefficients[2]
      
      ## save mean 
      means[x,y] = mean(df$temp, na.rm = T)
    }
    y = y + 1
  }
  
  print(x)
  x = x + 1
}

## translate to deg C per year
temp_trends = temp_trends*365

## rasterize
temp_trend_rast = raster(temp_trends, na)
plot(temp_trend_rast)
means_rast = raster(means, na)
plot(means_rast)

## calculate the spatial gradient from mean temperatures 
spat_grad = spatGrad(means_rast)
plot(spat_grad)

## calculate climate velocity 
cv = gVoCC(tempTrend = temp_trend_rast, spatGrad = spat_grad)
plot(cv)

## get latitudinal velocity component 
cos_angle = cos(cv$voccAng)
cv_lat = temp_trend_rast*cos_angle

## units are degrees latitude/year
## convert to km/y (approximately)
cv_lat = cv_lat*111

plot(cv_lat)
mean(values(cv_lat), na.rm = T) # 0.08 km/y
max(values(cv_lat), na.rm = T) # 7.99 km/y
min(values(cv_lat), na.rm = T) # -8.07 km/y 
sd(values(cv_lat), na.rm = T) # 3.57

hist(cv_lat)
hist(shifts$slope*111)

## save climate velocity raster
writeRaster(cv_lat, "outputs/data-processed/cv_north-america_1980-2021.tif")


## calculate across the cardinal's range
botw = vect("outputs/data-processed/BBS_BOTW.shp")

## filter to cardinal 
range = botw[which(botw$sci_name == "Cardinalis cardinalis"),]

cv_lat = rast(cv_lat)

crs(cv_lat) = crs(range) 

test = mask(cv_lat, range)
plot(test)

mean(values(test), na.rm = T)

shifts$slope[which(shifts$genus_sp == 'Cardinalis cardinalis' & shifts$range_edge == "leading")]*110

##########################################
##          standard deviation          ##
##########################################
## measure standard deviation of anomalies and see how much variation there is across space
## also just plot some and look at them 
## how long are these good and bad periods? how good and bad?
library(ncdf4)

## read in seasonally detrended temperature data
s_filenames = readRDS("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/BerkeleyEarth_sp_files.rds")
s_filenames <- str_replace_all(s_filenames, "spatial_temps", 's-detrended')

## make date vector
## first date: 1880/1/1
## last date: 2021/12/31
dates = seq(ymd("1880-01-01"), ymd("2021-12-31"), by = "day")
## remove leap year dates 
dates = dates[-which(str_split_fixed(dates, "-", 2)[,2] == "02-29")]

## make lat long vectors 
lat <- seq(from = 89.5, to = -89.5, length.out = 180) 
lon <-seq(from = 0.5, to = 359.5, length.out = 360) 

lon_index = 0
lat_index = 90
count = 1
while(count <= length(s_filenames)) {
  
  ## get lat and lon bounds to extract in between:
  lon_bound1 <- lon_index 
  lon_bound2 <- lon_index + 60
  lat_bound1 <- lat_index
  lat_bound2 <- lat_index - 60
  
  ## open file 
  s_open = nc_open(s_filenames[count])
  s_detrended_tas = ncvar_get(s_open, "var1_1")
  nc_close(s_open)
  
  ## plot map 
  #test = rast(s_detrended_tas)
  #plot(test[[1]])
  
  ## change to xyz dataframe
  t = reshape::melt(s_detrended_tas)
  
  ## rename
  t = rename(t, "y" = X1, "x" = X2, "time" = X3, "anomaly" = value)
  
  ## fix date
  key = data.frame(time = 1:length(dates), date = dates)
  t = left_join(t, key) %>%
    select(-time)
  
  ## fix latitude
  lat_key = data.frame(y = 1:60, lat = seq(lat_bound1-0.5, lat_bound2))
  t = left_join(t, lat_key)
  
  ## fix longitude
  lon_key = data.frame(x = 1:60, lon = seq(lon_bound1+0.5, lon_bound2)-180)
  t = left_join(t, lon_key) 
    
  ## rearrange 
  t = t %>% select(-x, -y) %>% select(lon, lat, date, anomaly)
  
  test_t = filter(t, date %in% unique(dates)[1])
  
  test = rast(test_t, type = "xyz")
  plot(test)

  
  ## look at one with high autocorrelation versus low
  ## read in autocorrelation file 
  se = readRDS("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/BerkeleyEarth_noise-colour.rds")
  
  se %>%
    select(lat, lon, s_spec_exp_PSD_low_all) %>%
    unique()
  
  ## find highest
  sub = se %>%
    filter(lat < max(t$lat), lon < max(t$lon),
           lat > min(t$lat), lon > min(t$lon)) 
  sub[which(sub$s_spec_exp_PSD_low_all == max(sub$s_spec_exp_PSD_low_all)),]
  
  ## get one time series and plot it 
  one = filter(t, t$lat == 56.5)
  one = filter(one, one$lon == -121.5)
  
  one %>%
    ggplot(aes(x = date, y = anomaly)) +
    geom_line(linewidth = 0.1) +
    geom_hline(yintercept = 0, colour = "red")
  
  ## zoom in on a year  
  one %>%
    filter(date >= "1980-01-01") %>%
    filter(date <= "1981-01-01") %>%
    ggplot(aes(x = date, y = anomaly)) +
    geom_line(linewidth = 0.1) +
    geom_hline(yintercept = 0, colour = "red")
  
  ## calculate sd 
  sd(one$anomaly, na.rm = TRUE) ## 2.71

  ## find lowest
  sub[which(sub$s_spec_exp_PSD_low_all == min(sub$s_spec_exp_PSD_low_all)),]
  
  ## get one time series and plot it 
  one = filter(t, t$lat == 52.5)
  one = filter(one, one$lon == -177.5)
  
  one %>%
    ggplot(aes(x = date, y = anomaly)) +
    geom_line(linewidth = 0.1) +
    geom_hline(yintercept = 0, colour = "red")
  
  ## zoom in on a year  
  one %>%
    filter(date >= "1980-01-01") %>%
    filter(date <= "1981-01-01") %>%
    ggplot(aes(x = date, y = anomaly)) +
    geom_line(linewidth = 0.1) +
    geom_hline(yintercept = 0, colour = "red")
  
  ## calculate sd 
  sd(one$anomaly, na.rm = TRUE) ## 2.02
  
  
  ## create raster from data frame 
  se$lon = se$lon - 180
  r = se %>%
    select(lon, lat, s_spec_exp_PSD_low_all) %>%
    rast()
  
  plot(r)
  
  
  
  
  
  ## advance lat and long indecies to move to next chunk
  if (count == 6) {
    lat_index <- lat_index + 60
    lon_index <- 0
  }
  else if (count == 12) {
    lat_index <- lat_index + 60
    lon_index <- 0
  }
  else {
    lon_index <- lon_index + 60
  }
  
  spec_exp_df <- c()
  
  count = count + 1
}



