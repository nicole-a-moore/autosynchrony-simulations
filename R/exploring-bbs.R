## playing with BBS data
library(tidyverse)
library(tidyterra)
library(sf)
library(terra)

#####################################################
##                  read in metadata               ##
#####################################################
## read in weather data
weather = read.csv("data-raw/BBS/Weather.csv")
## filter out routes with RunType == 0
weather = filter(weather, RunType != 0)
## make code for route
weather$code = paste(weather$Route, weather$CountryNum, weather$StateNum, weather$Year, sep = "_")

## read in species list
sp <- read.csv("data-raw/BBS/SpeciesList.csv")
nrow(sp) == length(unique(sp$English_Common_Name)) ## 1 row per sp

## read in route data 
routes <- read.csv("data-raw/BBS/Routes.csv") %>%
  select(CountryNum, StateNum, Route, RouteName, Latitude, Longitude)

## see how many routes were sampled in all years
min(unique(weather$Year)) ## 1966
max(unique(weather$Year)) ## 2023
length(unique(weather$Year)) ## 57 years overall!!

weather %>% 
  group_by(Route, CountryNum, StateNum) %>%
  tally(.) %>% 
  arrange(-n) %>% 
  ggplot(aes(x = n)) +
  geom_histogram()

########################################################
## read in species-level survey data from 1966 - 1977 ##
########################################################
folder = "data-raw/BBS/States"
files = list.files(folder, full.names = T)

file = 1
while(file <= length(files)) {
  
  cur = read.csv(files[[file]])
  
  ## filter to only routes with runtype == 1
  cur <- cur %>%
    mutate(code = paste(Route, CountryNum, StateNum, Year, sep = "_")) %>%
    filter(code %in% weather$code)
  
  ## rename "SpeciesTotal" to "TotalAbd"
  cur <- rename(cur, "TotalAbd" = SpeciesTotal)
  
  ## get rid of Count columns and StopTotal
  cols = c(which(str_detect(colnames(cur), "Count")), which(str_detect(colnames(cur), "StopTotal")))
  cur <- cur[,-c(cols)]
  cur <- left_join(cur, routes) ## add route info
  
  if(file == 1) {
    cur_all = cur
  } 
  else {
    cur_all <- rbind(cur_all, cur)
  }

  print(file)
  file = file + 1
}

#######################################################
## read in species-level survey data from 1997 - now ##
#######################################################
folder = "data-raw/BBS/50-StopData"
files = list.files(folder, full.names = T)

file = 1
while(file <= length(files)) {
  
  cur = read.csv(files[[file]])
  
  ## filter to only routes with runtype == 1
  cur <- cur %>%
    mutate(code = paste(Route, CountryNum, StateNum, Year, sep = "_")) %>%
    filter(code %in% weather$code)

  ## calculate total abundance across all 50 plots across transect 
  cols = which(str_detect(colnames(cur), "Stop"))
  cur$TotalAbd = apply(cur, 1, FUN = function(x) {
    return(sum(as.numeric(x[cols])))
  })
  
  cur <- cur[,-c(cols)]
  cur <- left_join(cur, routes) ## add route info
  
  ## plot some examples
  # ebb %>%
  #   filter(code == unique(.$code)[8]) %>%
  #   ggplot(aes(x = Year, y = TotalAbd)) +
  #   geom_point() + 
  #   geom_line(aes(group = sequence))

  cur_all <- rbind(cur_all, cur)
 
  print(file)
  file = file + 1
}

## save 
write.csv(cur_all, "outputs/data-processed/BBS_all.csv", row.names = F)

cur_all <- read.csv("outputs/data-processed/BBS_all.csv")

## create a 'study extent' polygon that encompasses all sites in the bbs
codes = paste(cur_all$Route, cur_all$CountryNum, cur_all$StateNum, sep = "_") %>% unique()
routes$code = paste(routes$Route, routes$CountryNum, routes$StateNum, sep = "_") 
routes = routes[which(routes$code %in% codes),]

routes_vec = vect(routes, geom = c("Longitude", "Latitude"))

poly <- convHull(routes_vec)

plot(poly)
points(routes_vec)

## crop by Canada and US
countries <- necountries::countries(c("Canada", "United States of America"), part = TRUE)
countries <- st_union(countries)
countries <- vect(countries)

poly <- crop(countries, poly)

plot(poly)

## save
writeVector(poly, "outputs/data-processed/bbs_study-extent.shp", overwrite = T)

## filter to data on a single species 
ebb_all <- filter(cur_all, AOU == 7660) # eastern blue bird 

## add absence data 
# group by route, add years that were sampled on that route but where species wasn't found
sampling_scheme <- weather %>% 
  select(Route, CountryNum, StateNum, Year, code, RouteDataID, RPID)

ebb_all = ebb_all %>%
  left_join(sampling_scheme, .) %>%
  mutate(TotalAbd = ifelse(is.na(TotalAbd), 0, TotalAbd)) %>%
  mutate(temp = paste(Route, CountryNum, StateNum)) %>%
  group_by(temp) %>%
  tidyr::fill(., RPID, AOU, Latitude, Longitude, RouteName, .direction = "updown") %>%
  ungroup() %>%
  filter(!is.na(AOU)) %>% ## get rid of rows representing 0 abd in plots where spp was never seen once
  select(-temp) %>%
  distinct()

## plot across space 
ebb_all = vect(ebb_all, geom = c("Longitude", "Latitude"))

## bring in North America outline 
# get map data
crs(ebb_all) <- as.character(crs(countries))

ebb_all = st_as_sf(ebb_all)

min(ebb_all$Year)
year = 1975
ebb_all %>%
  filter(!TotalAbd == 0 & Year == year) %>%
  ggplot(aes(colour = TotalAbd)) +
  geom_sf(inherit.aes = FALSE, data = countries) +
  geom_sf(size = 1, data = filter(ebb_all, TotalAbd == 0 & Year == year), inherit.aes = FALSE, 
          fill = "transparent", linewidth = 0.0001, shape = 1) +
  geom_sf(size = 1) +
  scale_colour_gradient(trans = "log") +
  scale_x_continuous(limits = c(-130,-60)) +
  scale_y_continuous(limits = c(25,60))
  
## actually, can we use gappy data for the time series analysis?
## establishment = appearance after 3+ consecutive years of sampling without finding spp. 
## extinction = disappearance for 3+ consecutive years of sampling after presence 

## for the estimation of the range edges, I think we need to use only cites that were sampled every year
## could we agreggate sites within an x km radius to fill in gappy data? 


###################################################
##          get consecutive series               ##
###################################################
## count the number of consecutive periods each route was sampled 
consec_count = weather %>% 
  arrange(Route, CountryNum, StateNum, Year) %>%  
  group_by(Route, CountryNum, StateNum) %>%
  mutate(deltaLag = Year - lag(Year, 1)) %>% 
  group_by(Route, CountryNum, StateNum, Year) %>%
  mutate(sequence1 = case_when(is.na(deltaLag) | deltaLag > 1 ~ 1,
                               TRUE ~ 2)) %>% 
  ungroup() %>% 
  mutate(sequence = cumsum(sequence1==1)) %>% 
  select(-deltaLag, -sequence1) %>% 
  group_by(Route, CountryNum, StateNum, sequence) %>%
  tally() 

consec_count %>%
  ggplot(aes(x = n)) +
  geom_histogram()

length(unique(consec_count$Route[which(consec_count$n >= 10)])) ## 349 routes
length(which(consec_count$n >= 10)) #3701 time series

## filter to bouts of at least 10 consecutive years of samples 
consec_count <- consec_count %>% 
  mutate(code_temp = paste(Route, CountryNum, StateNum, sequence, sep = "_")) 
## %>% filter(n >= 10)

gt_10 = weather %>% 
  arrange(Route, CountryNum, StateNum, Year) %>%  
  group_by(Route, CountryNum, StateNum) %>%
  mutate(deltaLag = Year - lag(Year, 1)) %>% 
  group_by(Route, CountryNum, StateNum, Year) %>%
  mutate(sequence1 = case_when(is.na(deltaLag) | deltaLag > 1 ~ 1,
                               TRUE ~ 2)) %>% 
  ungroup() %>% 
  mutate(sequence = cumsum(sequence1==1)) %>% 
  select(-deltaLag, -sequence1) %>% 
  mutate(code_temp = paste(Route, CountryNum, StateNum, sequence, sep = "_")) %>%
  filter(code_temp %in% consec_count$code_temp) %>%
  select(-code_temp) %>%
  mutate(code = paste(Route, CountryNum, StateNum, sep = "_")) %>%
  select(RouteDataID, CountryNum, StateNum, Route, Year, ObsN, TotalSpp, sequence, code)

## plot some examples
gt_10 %>%
  filter(code == unique(.$code)[4]) %>%
  ggplot(aes(x = Year, y = TotalSpp)) +
  geom_point() + 
  geom_line(aes(group = sequence))

## combine 
gt_10 <- left_join(gt_10, routes)

## filter species-level data to only consecutive series 
gt_10 <- cur_all %>%
  mutate(code = paste(Route, CountryNum, StateNum, sep = "_")) %>%
  left_join(gt_10, .)

gt_10 %>%
  ggplot(aes(x = Year, y = TotalAbd, colour = code)) +
  geom_point() + 
  geom_line(aes(group = sequence)) +
  theme(legend.position = "none")


## filter to bb 
ebb_all <- filter(gt_10, AOU == 7660) %>%
  mutate(code = paste(Route, CountryNum, StateNum, sep = "_"))

ebb_all %>%
  filter(code == unique(.$code)[166]) %>%
  ggplot(aes(x = Year, y = TotalAbd, colour = code)) +
  geom_point() + 
  geom_line(aes(group = sequence)) +
  theme(legend.position = "none")


## read in botw ranges 
botw = st_read("outputs/data-processed/BBS_BOTW.shp")

## crop by botw study area 
poly = vect("outputs/data-processed/bbs_study-extent.shp")

botw_int = st_intersection(botw, poly)
split <- c()
for(i in 1:nrow(botw_int)) {
  one = botw_int[i,]
  
  # Calculate the latitudinal midpoint
  midpoint_lat <- st_coordinates(st_centroid(one))[2]
  
  ## make two separate polygons by cutting with bbox 
  bbox = st_bbox(one)
  bbox_upper = st_bbox(c(xmin = as.numeric(bbox$xmin), xmax = as.numeric(bbox$xmax), 
                         ymax = as.numeric(bbox$ymax), ymin = midpoint_lat))
  bbox_lower = st_bbox(c(xmin = as.numeric(bbox$xmin), xmax = as.numeric(bbox$xmax), 
                         ymax = midpoint_lat, ymin = as.numeric(bbox$ymin)))

  ## turn into linestrings
  one = st_cast(one, "MULTILINESTRING")
  
  ## crop
  upper = st_crop(one, bbox_upper)
  lower = st_crop(one, bbox_lower)

  upper$half = "upper"
  lower$half = "lower"
  
  split = rbind(split, upper, lower)
  
  print(i)
}

## initialize a raster to hold counts
r <- rast(extent = ext(botw))
line_raster_lower <-line_raster_upper <- r
values(line_raster_lower) <- values(line_raster_upper) <- 0

## rasterize each line segment individually and count overlaps
for (i in seq_len(nrow(split))) {
  ## rasterize one segment
  line_raster_cur <- rasterize(vect(split[i,]), r, field=1, touches=TRUE)
  
  ## add to the count raster
  if(split$half[i] == "upper") {
    line_raster_upper <- mosaic(line_raster_upper, line_raster_cur, fun = "sum")
  }
  else {
    line_raster_lower <- mosaic(line_raster_lower, line_raster_cur, fun = "sum")
  }
  
  print(i)
}

## plot result
plot(line_raster_lower, main = "Number of trailing range edges per cell")
plot(line_raster_upper, main = "Number of leading range edges per cell")

## crop by study area
line_raster_lower = mask(line_raster_lower, poly)
line_raster_upper = mask(line_raster_upper, poly)

## turn 0 to NA
line_raster_lower[line_raster_lower == 0] = NA
line_raster_upper[line_raster_upper == 0] = NA

ggplot() +
  geom_spatraster(data = line_raster_lower) +
  geom_sf(inherit.aes = FALSE, data = countries, fill = "transparent") +
  scale_x_continuous(limits = c(-180,-50)) +
  scale_y_continuous(limits = c(0,90)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_viridis_c(na.value = "transparent")

ggplot() +
  geom_spatraster(data = line_raster_upper) +
  geom_sf(inherit.aes = FALSE, data = countries, fill = "transparent") +
  scale_x_continuous(limits = c(-180,-50)) +
  scale_y_continuous(limits = c(0,90)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_viridis_c(na.value = "transparent")

## plot the range edges on top of each other 
upper = split[which(split$half == "upper"),]
upper = mask(vect(upper), poly)
ggplot() +
  geom_spatvector(data = upper[1:10,], colour = "red", linewidth = 0.1) +
  geom_spatvector(data = poly, fill = "transparent") +
  geom_sf(inherit.aes = FALSE, data = countries, fill = "transparent") +
  scale_x_continuous(limits = c(-180,-50)) +
  scale_y_continuous(limits = c(0,90)) +
  theme_bw() +
  theme(panel.grid = element_blank()) 

lower = split[which(split$half == "lower"),]
lower = mask(vect(lower), poly)
ggplot() +
  geom_spatvector(data = lower[1:10,], colour = "red", linewidth = 0.1) +
  geom_spatvector(data = poly, fill = "transparent") +
  geom_sf(inherit.aes = FALSE, data = countries, fill = "transparent") +
  scale_x_continuous(limits = c(-180,-50)) +
  scale_y_continuous(limits = c(0,90)) +
  theme_bw() +
  theme(panel.grid = element_blank()) 


