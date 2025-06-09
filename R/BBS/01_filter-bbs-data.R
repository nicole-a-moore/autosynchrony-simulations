## filtering BBS data
library(tidyverse)
library(tidyterra)
library(sf)
library(terra)
theme_set(theme_bw())

###################################################
##                add absences/NAs               ##
###################################################
## read in bbs data for all spp
cur_all <- read.csv("outputs/data-processed/BBS_all.csv")

## read in species list
sp <- read.csv("data-raw/BBS/SpeciesList.csv")
nrow(sp) == length(unique(sp$English_Common_Name)) ## 1 row per sp

## filter to songbirds
sp = sp[which(sp$Order == "Passeriformes"),]
length(unique(sp$AOU)) # 354

## get rid of hybrids 
sp = filter(sp, !str_detect(sp$Species, " x "))

## get rid of unidentified species 
sp = filter(sp, !str_detect(sp$English_Common_Name, "unid."))

cur_all = cur_all[which(cur_all$AOU %in% sp$AOU),]
length(unique(cur_all$AOU)) # 309

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

length(unique(sp$AOU)) ## 303 species after regrouping 

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

nrow(cur_all)

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

## add absence data 
# group by route, add years that were sampled on that route but where species wasn't found
sampling_scheme <- weather %>% 
  select(Route, CountryNum, StateNum, Year, code, RouteDataID, RPID)

## split by species
split = cur_all %>%
  group_by(AOU) %>% 
  group_split(.) 

## for each species, add absence data
bbs_clean = c()
for(bird in 1:length(split)) {
  cur = split[[bird]]
  
  cur = cur %>%
    left_join(sampling_scheme, .) %>%  
    mutate(TotalAbd = ifelse(is.na(TotalAbd), 0, TotalAbd)) %>% 
    mutate(temp = paste(Route, CountryNum, StateNum)) %>%
    group_by(temp) %>%
    tidyr::fill(., RPID, AOU, Latitude, Longitude, RouteName, .direction = "updown") %>%
    ungroup() %>%
    filter(!is.na(AOU)) %>% ## get rid of rows representing 0 abd in plots where spp was never seen once
    select(-temp) %>%
    distinct()
  
  ## for each route, add in rows with NA for years that were not sampled within the sampling period
  cur <- cur %>%
    select(-code) %>%
    mutate(route = paste(Route, CountryNum, StateNum, sep = "_")) 
  
  bbs_key <- cur %>% 
    group_by(route) %>%
    mutate(first_year = min(Year), last_year = max(Year)) %>%
    select(route, first_year, last_year) %>%
    st_drop_geometry() %>%
    distinct()
  
  for(i in 1:nrow(bbs_key)) {
    df <- data.frame(Year = seq(bbs_key$first_year[i], bbs_key$last_year[i], by = 1), route = bbs_key$route[i])
    
    if(i==1) {
      df_all = df
    }
    else {
      df_all = rbind(df_all, df)
    }
  }
  
  t = left_join(df_all, cur) %>%
    group_by(route) %>% ## group by route and species 
    tidyr::fill(colnames(.)[which(colnames(.) != "TotalAbd")], .direction = "updown") %>%
    as.data.frame() 
  
  ## combine to the rest
  bbs_clean = rbind(bbs_clean, t)
  
  print(bird)
}

## save
write.csv(bbs_clean, "outputs/data-processed/bbs_songbirds_with_absence.csv", row.names = FALSE)

bbs_clean <- read.csv("outputs/data-processed/bbs_songbirds_with_absence.csv")


## get rid of species with:
## - fewer than 50 occurrences 
## - less than 10 years of occurrences 
## - for whom BBS covers less than 70% of breeding and/or resident distribution in north america
## then, for now, just look at birds with both range edges in study area (according to Rushing et al. Table 1)

## fewer than 50 occurrences
bbs_key = bbs_clean %>%
  group_by(AOU) %>%
  mutate(n_ind = sum(TotalAbd, na.rm = TRUE)) %>% 
  ungroup() %>%
  filter(n_ind >= 50) %>%
  select(-n_ind)

length(unique(bbs_key$AOU)) ## 291 species remain

## more than 10 years
bbs_key = bbs_key %>%
  filter(!is.na(TotalAbd), TotalAbd != 0) %>%
  group_by(AOU) %>%
  mutate(n_years = length(unique(Year))) %>% 
  ungroup() %>% 
  filter(n_years >= 10) %>%
  select(-n_years)

length(unique(bbs_key$AOU)) ## 288 species remain

## across at least 10 routes
bbs_key = bbs_key %>%
  group_by(AOU) %>%
  mutate(n_routes = length(unique(route))) %>% 
  ungroup() %>% 
  filter(n_routes >= 10) %>%
  select(-n_routes)

length(unique(bbs_key$AOU)) ## 270 species remain

bbs_key <- bbs_key %>%
  select(AOU, route) %>%
  distinct()

bbs_clean = filter(bbs_clean, paste(bbs_clean$route, bbs_clean$AOU) %in% paste(bbs_key$route, bbs_key$AOU))
  
## read in Table 1 list from rushing et all
rushing = read.csv("data-raw/Rushing_et_al_2020_Table1.csv")

## plot some random species on a map
name_key = select(sp, AOU, genus_sp) %>%
  distinct()

bbs_clean = left_join(bbs_clean, name_key)

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









