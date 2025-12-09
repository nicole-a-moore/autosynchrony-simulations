## filtering BBS data
library(tidyverse)
library(tidyterra)
library(sf)
library(terra)
theme_set(theme_bw())
select = dplyr::select

###################################################
##                add absences/NAs               ##
###################################################
## read in bbs data for all spp
cur_all <- read.csv("outputs/data-processed/BBS_all.csv")

## read in species list
sp <- read.csv("data-raw/BBS/SpeciesList.csv")
nrow(sp) == length(unique(sp$English_Common_Name)) ## 1 row per sp

## filter to landbirds
sp = sp[which(!sp$Order %in% c("Anseriformes", "Charadriiformes", "Galliformes", "Pelecaniformes", "Suliformes", "Gaviiformes", "Ciconiiformes", "Podicipediformes", "Gruiformes", "Phaethontiformes", "Procellariiformes")),]
length(unique(sp$AOU)) # 500

## get rid of hybrids 
sp = filter(sp, !str_detect(sp$Species, " x "))

## get rid of unidentified species 
sp = filter(sp, !str_detect(sp$English_Common_Name, "unid."))

cur_all = cur_all[which(cur_all$AOU %in% sp$AOU),]
length(unique(cur_all$AOU)) # 440

## read in study area extent
extent = vect("outputs/data-processed/bbs_study-extent.shp")

## read in info about each route
weather = read.csv("data-raw/BBS/Weather.csv")
## filter out routes with RunType == 0
weather = filter(weather, RunType != 0)
## make code for route
weather$code = paste(weather$Route, weather$CountryNum, weather$StateNum, weather$Year, sep = "_")

## get canada and us map for plotting
countries <- necountries::countries(c("Canada", "United States of America"), part = TRUE)
countries <- st_union(countries)
countries <- vect(countries)

## group species and subspecies together by recoding AOU
sp <- filter(sp, AOU %in% cur_all$AOU)

sp$genus_sp = paste0(sp$Genus, " ", sp$Species)

subsp = filter(sp, str_split_fixed(sp$Species, " ", 2)[,2] != "")
sp = filter(sp, str_split_fixed(sp$Species, " ", 2)[,2] == "")

subsp$genus_sp = paste0(subsp$Genus, " ", str_split_fixed(subsp$Species, " ", 2)[,1])
subsp$old_AOU = subsp$AOU
subsp <- select(subsp, -AOU)
## give subspecies a new AOU
subsp = subsp %>% group_by(genus_sp) %>%
  mutate(AOU = paste("subsp_gr_", cur_group_id()))
sp$old_AOU = NA

sp <- rbind(sp, subsp)

length(unique(sp$AOU)) ## 433 species after regrouping 

key = select(sp, old_AOU, AOU) %>%
  filter(!is.na(old_AOU))

## change AOU of subspecies to AOU of greater species 
old_data = filter(cur_all, AOU %in% key$old_AOU)
cur_all = filter(cur_all, !AOU %in% key$old_AOU)

old_data$old_AOU = old_data$AOU
old_data <- select(old_data, -AOU)
new_data <- left_join(old_data, key) %>%
  select(-old_AOU)

cur_all <- rbind(cur_all, new_data)

## make sure we have the right number of species 
length(unique(cur_all$AOU)) == length(unique(sp$AOU))

nrow(cur_all) ## 9306792

## group by species, route, and year and sum abundance to consolidate abundance counts within a species 
cur_all <- cur_all %>%
  mutate(route = paste(Route, CountryNum, StateNum, sep = "_")) %>%
  group_by(AOU, route, Year) %>%
  mutate(TotalAbd = sum(TotalAbd)) %>%
  ungroup() %>%
  distinct()

## save 
write.csv(cur_all, "outputs/data-processed/BBS_songbirds.csv", row.names = FALSE)

cur_all <- read.csv("outputs/data-processed/BBS_songbirds.csv")

## read in route data 
routes <- read.csv("data-raw/BBS/Routes.csv") %>%
  select(CountryNum, StateNum, Route, RouteName, Latitude, Longitude)

## add absence data 
## add rows of 0 abundance for species if route was sampled in that year but species wasn't seen
sampling_scheme <- weather %>% 
  select(Route, CountryNum, StateNum, Year, code, RouteDataID, RPID) %>%
  left_join(routes) %>%
  mutate(route = paste(Route, CountryNum, StateNum, sep = "_")) 

## make list of all routes x all species 
grid = expand.grid(route = unique(sampling_scheme$route), AOU = unique(cur_all$AOU))
sampling_scheme <- left_join(sampling_scheme, grid, relationship = "many-to-many")

bbs_clean = cur_all %>%
  left_join(sampling_scheme, .) %>%
  group_by(AOU) %>%
  mutate(TotalAbd = ifelse(is.na(TotalAbd), 0, TotalAbd))   

length(which(is.na(bbs_clean$TotalAbd))) ## 0 NA
length(which(bbs_clean$TotalAbd == 0)) ## 49212413 zeros

## add NA for missed years that each route wasn't sampled 
year_key = bbs_clean %>% 
  group_by(route) %>%
  summarise(min_year = min(Year),
         max_year = max(Year)) %>%
  select(route, min_year, max_year) %>%
  distinct()

na_key <- c()
for(route in 1:length(unique(year_key$route))) {
  sub = year_key[which(year_key$route == unique(year_key$route)[route]),]
  na_key <- rbind(na_key, data.frame(route = unique(sub$route), 
                                     Year = seq(sub$min_year[1], sub$max_year[1])))
}

bbs_clean = left_join(na_key, bbs_clean) %>%
  group_by(route) %>%
  tidyr::fill(., Route, CountryNum, StateNum, RouteDataID, RPID, RouteName, Latitude,
              Longitude, AOU, .direction = "updown") %>% select(-code)
  
length(which(bbs_clean$TotalAbd == 0)) ## 49212413 zeros
length(which(is.na(bbs_clean$TotalAbd))) ## 44913 NAs

## save
write.csv(bbs_clean, "outputs/data-processed/bbs_landbirds_with_absence.csv", row.names = FALSE)

####################################################
##      filter species by inclusion criteria      ##
####################################################
bbs_clean <- read.csv("outputs/data-processed/bbs_landbirds_with_absence.csv")

## get rid of species with:
## - fewer than 50 occurrences 
## - less than 10 years of occurrences 
## - for whom BBS covers less than 70% of breeding and/or resident distribution in north america
## then, for now, just look at birds with both range edges in study area (according to Rushing et al. Table 1)

## GET RID OF SPECIES 
## with fewer than 50 occurrences across any route 
bbs_key = bbs_clean %>%
  group_by(AOU) %>%
  mutate(n_ind = sum(TotalAbd, na.rm = T)) %>% 
  ungroup() %>%
  filter(n_ind >= 50) %>%
  select(-n_ind)

length(unique(bbs_key$AOU)) ## 393 species remain

sp_to_keep = unique(bbs_key$AOU)

## with fewer than 10 years of presence across any route
bbs_key = bbs_clean %>%
  filter(TotalAbd != 0, !is.na(TotalAbd)) %>%
  group_by(AOU) %>%
  mutate(n_years = length(unique(Year))) %>% 
  ungroup() %>% 
  filter(n_years >= 10) %>%
  select(-n_years)

sp_to_keep <- sp_to_keep[which(sp_to_keep %in% bbs_key$AOU)]
length(unique(sp_to_keep)) ## 391 species remain

## that are not present across at least 10 routes
bbs_key = bbs_clean %>%
  filter(TotalAbd != 0, !is.na(TotalAbd)) %>%
  group_by(AOU) %>%
  mutate(n_routes = length(unique(route))) %>% 
  ungroup() %>% 
  filter(n_routes >= 10) %>%
  select(-n_routes)

sp_to_keep <- sp_to_keep[which(sp_to_keep %in% bbs_key$AOU)]
length(unique(sp_to_keep)) ## 366 species remain

## GET RID OF SITES 
## remove sites with only 1 year sampled 
bbs_key <- bbs_clean %>%
  filter(!is.na(TotalAbd)) %>%
  select(route, Year) %>%
  distinct() %>%
  group_by(route) %>%
  filter(n() > 1)

routes_to_keep = unique(bbs_key$route) 
length(routes_to_keep) ## 5141

## filter dataset to species to keep and routes to keep
bbs_clean = filter(bbs_clean, route %in% routes_to_keep, AOU %in% sp_to_keep)

name_key = select(sp, AOU, genus_sp) %>%
  distinct()

bbs_clean = left_join(bbs_clean, name_key)

## stats:
length(unique(bbs_clean$genus_sp)) ## 366 species
length(unique(bbs_clean$route)) ## 5141 routes

## save 
write.csv(bbs_clean, "outputs/data-processed/bbs_clean-subset.csv", row.names = F)
  



##### garbage

## read in Table 1 list from rushing et all
rushing = read.csv("data-raw/Rushing_et_al_2020_Table1.csv")

## plot some random species on a map

## add coordinates to dataframe
bbs_vect = vect(bbs_clean, geom = c("Longitude", "Latitude"))
bbs_clean = st_as_sf(bbs_vect)
st_crs(bbs_clean) = st_crs(countries)

onesp = bbs_clean %>%
  filter(AOU == unique(bbs_clean$AOU)[15]) 

onesp %>%
  filter(!TotalAbd == 0) %>%
  ggplot(aes(colour = TotalAbd)) +
  geom_sf(inherit.aes = FALSE, data = countries) +
  geom_sf(size = 1, data = filter(onesp, TotalAbd == 0), inherit.aes = FALSE, 
          fill = "transparent", linewidth = 0.0001, shape = 1) +
  geom_sf(size = 1) +
  scale_colour_gradient(trans = "log", limits = c(1, 191)) +
  scale_x_continuous(limits = c(-165,-55)) +
  scale_y_continuous(limits = c(25,80)) +
  labs(colour = "Total abundance", title = unique(onesp$genus_sp))



## then, map each bird's sampled routes on top of their distributions 
#botw = st_read("outputs/data-processed/BBS_BOTW.shp")

## then, calculate 10 northernmost and southernmost latitudinal cells where > 1 individual was present 
## track over time 
## (see if you need to group into time periods rather than using single years)


## try 






