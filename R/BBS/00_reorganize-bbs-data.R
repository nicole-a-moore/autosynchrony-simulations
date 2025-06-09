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

#####################################################################
##          get breeding and resident range polygons               ##
#####################################################################
### see where species in the BBS have their range edges in the study area
botw = st_read("/Users/nikkimoore/Documents/bioshifts-traits/data-raw/large-data/BirdsOfTheWorld/BOTW/BOTW.gdb", "All_Species")

species = unique(paste0(sp$Genus, " ", sp$Species))

## filter BOTW to species list
botw = botw[which(botw$sci_name %in% species),]

## filter to breeding and resident range
botw = botw[which(botw$seasonal %in% c(1, 2)),]

## combine for each species 
botw_cb <- botw %>%
  group_by(sci_name) %>% 
  filter(n() > 1) %>%
  summarize(Shape = st_union(st_make_valid(st_cast(Shape, to = "MULTIPOLYGON")))) %>%
  ungroup()

botw_singles = botw %>%
  group_by(sci_name) %>% 
  filter(n() == 1) %>%
  ungroup() %>%
  select(sci_name, Shape) %>%
  mutate(Shape = st_make_valid(st_cast(Shape, to = "MULTIPOLYGON")))

## combine
botw = rbind(botw_cb, botw_singles)

## save 
st_write(botw, "outputs/data-processed/BBS_BOTW.shp", append = FALSE)

