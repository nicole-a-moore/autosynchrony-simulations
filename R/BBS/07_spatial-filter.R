## filtering to only data where a northern range edge can be accurately measured  
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(ggridges)
theme_set(theme_bw())

## read in clean bbs data
bbs = read.csv("outputs/data-processed/bbs_clean-subset.csv")

## filter out data for years where route wasn't sampled
bbs <- filter(bbs, !is.na(TotalAbd))

## filter to northern cardinal as example
nc <- filter(bbs, AOU == 5930) 

## northern cardinal 5930
## bc chickadee 7350 

## issues:
############################################
## 1. number of samples increases over time
nc %>%
  ggplot(aes(x = Year)) +
  geom_bar()
## see Shoo et al. 2006 for effect of uneven sampling over time on bias in shift detection

## 2. increase in the latitude of routes sampled over time 
nc %>%
  ggplot(aes(x = Year, y = Latitude)) +
  geom_point(alpha = 0.1)

## 3. fewer samples at higher latitudes
nc %>%
  ggplot(aes(y = Year, x = Latitude, group = Year)) +
  geom_density_ridges(alpha = 0.5, aes(x = Latitude), rel_min_height = 0.01, fill = "#F8766D") +
  coord_flip() 

## 4. lower detectability per unit effort (sampled route) at higher latitude 
## e.g., northern cardinal
nc %>%
  mutate(Latitude_bin = cut(Latitude, breaks=seq(0, 90, by = 5))) %>%
  group_by(Latitude_bin) %>%
  mutate(Detectability = length(which(TotalAbd != 0))/n(),
         n = n()) %>%
  ungroup() %>%
  ggplot(aes(x = Latitude_bin, y = Detectability, colour = n)) +
  geom_point() +
  labs(colour = "Number of samples\nin latitudinal bin")


## solution: 
############################################
## 1. for each species, classify data into 10 degree longitude bins
## 2. pool data for years between 1966-1980 together 
## 3. filter out years of observations where there are less than 50 observations/bin/year 
## 4. filter out routes in longitudinal bands with less than 50 presences, years with less than 10 presences per year 
## 5. measure northern edge of a) species presence and b) sampled routes in each longitudinal bin each year (e.g., 95th percentile of latitude)
## 6. get rid of observations where less than 10 sites were sampled north of the northern edge of presence within the bin that year

## effectively, this should:
##  - subset to species whose northern range edge is well within sampling domain
##  - filter to only years where sampling domain was far enough from range edge to detect a shift 

## apply the steps to each species
########################################
## for each species, take presences over all years and calculate latitudinal span 
bounds = bbs %>%
  group_by(AOU) %>%
  filter(TotalAbd != 0) %>%
  mutate(Num_breaks = floor((max(Longitude) - min(Longitude))/10)) %>%
  ungroup()

bounds_new <- c()
i = 1
while(i <= length(unique(bounds$AOU))) {
  temp = filter(bounds, AOU == unique(bounds$AOU)[i])
  temp_wthabs <- filter(bbs, AOU == unique(bounds$AOU)[i])
  
  ## get min and max longitude of presences
  x_min <- min(temp$Longitude)
  x_max <- max(temp$Longitude)
  
  if(unique(temp$Num_breaks) == 0) {
    x_min = floor((x_min + 0.25) * 2) / 2 - 0.25
    x_max = ceiling((x_max + 0.25) * 2) / 2 - 0.25
    
    temp_wthabs$Longitude_bin = rep(paste0("(",x_min, ",", x_max, "]"),  n = nrow(temp))
  }
  else if(unique(temp$Num_breaks) == 1) {
    x_min = floor((x_min + 0.25) * 2) / 2 - 0.25
    x_max = ceiling((x_max + 0.25) * 2) / 2 - 0.25
    
    temp_wthabs$Longitude_bin = cut(temp_wthabs$Longitude, breaks = c(x_min,
                                                                      floor(((x_min + (x_max - x_min)/2) + 0.25) * 2) / 2 - 0.25,
                                                               x_max))
  }
  else {
    ## get full longitude range
    total_range <- x_max - x_min
    
    # get num of full 10-wide bins
    n_full <- unique(temp$Num_breaks)
    
    ## calc remainder that does not fit into 10-wide bins
    remainder <- total_range - n_full * 10
    
    first_width = last_width = remainder/2
    
    first_break = x_min + first_width + 10
    last_break = last(seq(first_break + 10, by = 10, length.out = unique(temp$Num_breaks)-1)) + last_width 
    
    ## round to nearest 0.25 or 0.75
    first_break = floor((first_break + 0.25) * 2) / 2 - 0.25
    last_break = ceiling((last_break + 0.25) * 2) / 2 - 0.25
    
    ## if less than 5 degrees remain, add to previous break
    if(last_width <= 5) {
      
      first_val = floor(((first_break - first_width - 10) + 0.25) * 2) / 2 - 0.25
      last_val = last(seq(first_break + 10, by = 10, length.out = unique(temp$Num_breaks)-2)) + last_width
      last_val = ceiling((last_val + 0.25) * 2) / 2 - 0.25 + 10
    
      ## calculate breaks
      breaks <- c(first_val, first_break, 
                  seq(first_break + 10, by = 10, length.out = unique(temp$Num_breaks)-2), 
                  last_val)
    }
    else {
      breaks <- c(first_val, floor((x_min + 0.25) * 2) / 2 - 0.25,  first_break, 
                  seq(first_break + 10, by = 10, length.out = unique(temp$Num_breaks)-1), 
                  (ceiling(((last(seq(first_break + 10, by = 10, 
                                     length.out = unique(temp$Num_breaks)-1)) + 0.25) * 2) + last_width) / 2 - 0.25))
    }
   
    
    ## cut breaks
    temp_wthabs$Longitude_bin = cut(temp_wthabs$Longitude, breaks = breaks, dig.lab = 4)
  }
 
  bounds_new <- rbind(bounds_new, temp_wthabs)                                
  print(i)
  i = i + 1
}
  
## join back to all data 
bounds_new = bounds_new %>%
  select(AOU, Longitude_bin, route) %>%
  distinct()

bbs_sp <- left_join(bbs, bounds_new, relationship = "many-to-many")

## save 
write.csv(bbs_sp, "outputs/data-processed/bbs_with-longitude-breaks.csv", row.names = F)

bbs_sp <- read.csv("outputs/data-processed/bbs_with-longitude-breaks.csv")

## 2. pool data for years between 1966-1980 together 
bbs_sp = mutate(bbs_sp, Year = ifelse(Year <= 1980, 1980, Year)) %>%
  ## get rid of duplicate presence/absences per route per year 
  group_by(AOU, Year, route) %>%
  mutate(TotalAbd = max(TotalAbd)) %>%
  ungroup() %>%
  select(-RouteDataID) %>%
  distinct() 

## 3. filter out years of longitudinal bins where there are less than 50 observations/bin/year 
bbs_sp_filtered = bbs_sp %>%  
  group_by(AOU, Longitude_bin, Year) %>%
  mutate(routes_sampled_per_year = n()) %>%
  ungroup() %>%
  filter(routes_sampled_per_year >= 50) %>%
  ungroup()

## 4. filter out routes in longitudinal bins with less than 50 presences overall, count number of presences in latitudinal bin each year 
bbs_sp_filtered = group_by(bbs_sp_filtered, AOU, Longitude_bin) %>%
  mutate(n_presences = sum(TotalAbd != 0)) %>% 
  filter(n_presences >= 50) %>%
  select(-n_presences) %>%
  group_by(AOU, Longitude_bin, Year) %>%
  mutate(n_presences = sum(TotalAbd != 0)) 

## removed: remove years with less than x presences 
  # group_by(AOU, Longitude_bin, Year) %>%
  # mutate(n_presences = sum(TotalAbd != 0)) %>% 
  # filter(n_presences >= 5) %>%
  # select(-n_presences) 

## 5. measure northern edge of a) species presence and b) sampled routes in each longitudinal bin each year (e.g., mean latitude of 5 most northern cells)

range_edges = bbs_sp_filtered %>%
  filter(TotalAbd != 0) %>% 
  arrange(AOU, Longitude_bin, -Latitude) %>% 
  ungroup() %>%
  slice(1:5, .by = c("AOU", "Longitude_bin", "Year")) %>%
  group_by(AOU, Longitude_bin, Year) %>%
  summarize(presence_edge = mean(Latitude)) 
sample_edges = bbs_sp_filtered %>%
  arrange(AOU, Longitude_bin, -Latitude) %>% 
  ungroup() %>%
  slice(1:5, .by = c("AOU", "Longitude_bin", "Year")) %>%
  group_by(AOU, Longitude_bin, Year) %>%
  summarize(sample_edge = mean(Latitude)) 

bbs_sp_filtered <- left_join(bbs_sp_filtered, range_edges) %>%
  left_join(., sample_edges)

## old: use 95th percentile instead of mean 5
# bbs_sp = bbs_sp %>%
#   group_by(AOU, Year, Longitude_bin) %>%
#   mutate(presence_edge = quantile(Latitude[which(TotalAbd != 0)], 0.95),
#          sample_edge = quantile(Latitude, 0.95)) %>%
#   filter(!is.na(presence_edge)) %>%
#   ungroup()

## 6. get rid of observations where less than 10 sites were sampled north of the northern edge of presence within the bin that year
bbs_sp_filtered = bbs_sp_filtered %>%
  group_by(AOU, Year, Longitude_bin) %>%
  ## count number of routes sampled beyond presence edge
  mutate(n_north_of_edge = length(which(Latitude > presence_edge))) %>% 
  ## get rid of years where less than 20 sites were sampled north of the edge of presence 
  filter(n_north_of_edge >= 20) %>% 
  ungroup()

# 7. get rid of observations more than 1 presences are less than 1 degree latitude below or are higher in latitude than sampling edge within the bin that year
bbs_sp_filtered <- bbs_sp_filtered %>%
  group_by(AOU, Year, Longitude_bin) %>%
  mutate(num_above_sample_edge = length(which(Latitude > sample_edge & TotalAbd != 0))) %>% 
  filter(num_above_sample_edge <= 1) %>%
  ungroup() %>%
  select(-num_above_sample_edge) %>%
  gather(key = "type", value = "Latitude_of_edge", c(sample_edge, presence_edge))

## plot the data that's left (for northern cardinal)
bbs_sp_filtered %>%
  filter(AOU == 5930) %>%
  ggplot(aes(x = Year, y = Latitude_of_edge, colour = type)) +
  geom_point() +
  facet_wrap(~Longitude_bin) +
  labs(colour = "Edge type")

## how many species left?
length(unique(bbs_sp_filtered$AOU)) ## 322 species 

## add back filtered out data but flag
bbs_sp_filtered$keep = 1

bbs_sp = left_join(bbs_sp, bbs_sp_filtered)

##save data
write.csv(bbs_sp, "outputs/data-processed/BBS_spatial-filter.csv", row.names = F)


## map data to see if we are happy with the criteria
############################################
##          map range edge shift          ##
############################################
bbs_sp <- read.csv("outputs/data-processed/BBS_spatial-filter.csv") %>%
  filter(keep == 1)

## get canada and us map for plotting
countries <- necountries::countries(c("Canada", "United States of America"), part = TRUE)
countries <- st_union(countries)
# countries <- vect(countries)

## get lakes for plotting
lakes <- rnaturalearth::ne_download(scale = 110, 
                                    type = 'lakes', 
                                    category = 'physical') %>% 
  sf::st_as_sf(lakes110, crs = 4269)
## crop to north america
lakes = st_transform(lakes, st_crs(countries))
lakes = st_intersection(lakes, countries)

## get ecoregions for plotting
ecoregions = read_sf("data-raw/Ecoregions_LevelI/na_cec_eco_l1/NA_CEC_Eco_Level1.shp")

## crop to north america
ecoregions = st_transform(ecoregions, st_crs(countries))
ecoregions = st_intersection(ecoregions, countries)

## loop through each species 
sp = 1
while(sp <= length(unique(bbs$AOU))) {
  
  ## get sp name
  genus_sp = unique(bbs$genus_sp[which(bbs$AOU == unique(bbs$AOU)[sp])])
  
  ## subset filtered bbs routes to species 
  bbs_filtered <- filter(bbs_sp, AOU == unique(bbs$AOU)[sp]) 
  
  ## check if there is data for the species 
  data_left = nrow(bbs_filtered) != 0
  
  ## subset all bbs routes to species 
  bbs_all <- filter(bbs, AOU == unique(bbs$AOU)[sp]) 
  
  ## pool data for years between 1966-1980 together 
  bbs_all =  mutate(bbs_all, Year = ifelse(Year <= 1980, 1980, Year)) %>%
    ## get rid of duplicate presence/absences per route per year 
    group_by(Year, route) %>%
    mutate(TotalAbd = max(TotalAbd)) %>%
    ungroup() %>%
    select(-RouteDataID) %>%
    distinct() 
  
  ## add identifier for unique route + year
  bbs_all$route_year = paste(bbs_all$route, bbs_all$Year, sep = "_")
  bbs_filtered$route_year = paste(bbs_filtered$route, bbs_filtered$Year, sep = "_")
  
  ## turn data frames into sf data frame
  bbs_all <- vect(bbs_all, geom = c("Longitude", "Latitude")) %>%
    st_as_sf(.)
  bbs_filtered <- vect(bbs_filtered, geom = c("Longitude", "Latitude")) %>%
    st_as_sf(.)
  st_crs(bbs_all) = st_crs(countries)
  st_crs(bbs_filtered) = st_crs(countries)
  
  ## add coordinates back to df as columns
  bbs_all <- cbind(bbs_all, st_coordinates(bbs_all)) %>%
    rename("Latitude" = Y, "Longitude" = X)
  
  if(nrow(bbs_filtered) != 0) {
    bbs_filtered <- cbind(bbs_filtered, st_coordinates(bbs_filtered)) %>%
      rename("Latitude" = Y, "Longitude" = X)
    
    ## make columns for min and max longitude 
    bbs_filtered = bbs_filtered %>%
      mutate(xmin = str_split_fixed(Longitude_bin, ",", 2)[,1], 
             xmax = str_split_fixed(Longitude_bin, ",", 2)[,2]) %>%
      mutate(xmin = as.numeric(str_remove_all(xmin, "\\(")), 
             xmax = as.numeric(str_remove_all(xmax, "\\]")))
  }
 
  ## create directory for plots
  if(data_left) {
    if(!dir.exists(paste0("outputs/figures/range-edge-shifts_after-filtering/species-with-data/", 
                          str_replace_all(genus_sp, " ", "\\_")))) {
      dir.create(paste0("outputs/figures/range-edge-shifts_after-filtering/species-with-data/", 
                        str_replace_all(genus_sp, " ", "\\_")))
    }
    # if(!dir.exists(paste0("outputs/figures/range-edge-shifts_after-filtering/species-with-data/ecoregions/", 
    #                       str_replace_all(genus_sp, " ", "\\_")))) {
    #   dir.create(paste0("outputs/figures/range-edge-shifts_after-filtering/species-with-data/ecoregions/", 
    #                     str_replace_all(genus_sp, " ", "\\_")))
    # }
    path = "outputs/figures/range-edge-shifts_after-filtering/species-with-data/"
    path2 = "outputs/figures/range-edge-shifts_after-filtering/species-with-data/ecoregions/"
  }
  else {
    if(!dir.exists(paste0("outputs/figures/range-edge-shifts_after-filtering/species-no-data/", 
                          str_replace_all(genus_sp, " ", "\\_")))) {
      dir.create(paste0("outputs/figures/range-edge-shifts_after-filtering/species-no-data/", 
                        str_replace_all(genus_sp, " ", "\\_")))
    }
    # if(!dir.exists(paste0("outputs/figures/range-edge-shifts_after-filtering/species-no-data/ecoregions/", 
    #                       str_replace_all(genus_sp, " ", "\\_")))) {
    #   dir.create(paste0("outputs/figures/range-edge-shifts_after-filtering/species-no-data/ecoregions/", 
    #                     str_replace_all(genus_sp, " ", "\\_")))
    # }
    path = "outputs/figures/range-edge-shifts_after-filtering/species-no-data/"
    path2 = "outputs/figures/range-edge-shifts_after-filtering/species-no-data/ecoregions/"
  }
  
  ## plot data in each year and save 
  for(y in 1:length(unique(bbs_all$Year))) {
    
    year = unique(bbs_all$Year)[order(unique(bbs_all$Year))][y]
    
    if(nrow(bbs_filtered) != 0) {
      ## if data for that year
      if(length(which(bbs_all$Year == year & bbs_all$TotalAbd != 0)) != 0) {
        
        spread = spread(bbs_filtered, key = type, value = Latitude_of_edge)  %>%
          filter(Year == year) %>%
          select(Longitude_bin, xmax, xmin, sample_edge, presence_edge) %>% 
          st_drop_geometry() %>%
          distinct()
        
         p = bbs_all %>%
          ## filter to data from current year
          filter(Year == year) %>%
          ggplot(aes()) +
          geom_sf(inherit.aes = FALSE, data = countries) +
          geom_sf(data = lakes, mapping = aes(geometry = geometry), color = "black",
                   fill = "white") +
          ## all routes
          geom_sf(size = 0.5, colour = "darkgrey", alpha = 0.8, stroke = 0.1, shape = 4) +
          ## absences in areas being used for analysis
          geom_sf(size = 0.5, data = filter(bbs_filtered, TotalAbd == 0 & Year == year), inherit.aes = FALSE,
                  stroke = 0.1, shape = 4, colour = "black") +
          ## presences that aren't used in analysis
          geom_sf(data = filter(bbs_all, TotalAbd != 0 & Year == year),
                  colour = "red", alpha = 0.6, size = 0.1) +
          scale_colour_gradient(trans = "log", limits = c(1, max(bbs_filtered$TotalAbd, na.rm = T))) +
          scale_x_continuous(limits = c(-165,-48.5)) +
          scale_y_continuous(limits = c(25,70)) +
          labs(colour = "Total abundance", title = year) +
          ## presences in areas being used for analysis
          geom_sf(size = 0.1, data = filter(bbs_filtered, TotalAbd != 0 & Year == year), inherit.aes = FALSE,
                  aes(colour = TotalAbd)) +
          geom_segment(data = spread, inherit.aes = F,
                       aes(x = xmax, xend = xmin, y = presence_edge), colour = "red") +
          geom_segment(data = spread, inherit.aes = F,
                       aes(x = xmax, xend = xmin, y = sample_edge), colour = "darkgrey") +
          labs(x = "", y = "")
         
         # ## plot with ecoregions behind 
         # p2 = bbs_all %>%
         #   ## filter to data from current year
         #   filter(Year == year) %>%
         #   ggplot(aes()) +
         #   geom_sf(inherit.aes = FALSE, data = countries) +
         #   geom_sf(data = lakes, mapping = aes(geometry = geometry), color = "black",
         #           fill = "white") +
         #   geom_sf(data = ecoregions, mapping = aes(geometry = geometry, fill = NA_L1NAME, colour = NA_L1NAME), 
         #           alpha = 0.25) +
         #   ## all routes 
         #   geom_sf(size = 0.5, colour = "darkgrey", alpha = 0.8, stroke = 0.1, shape = 4) +
         #   geom_sf(data = filter(bbs_all, TotalAbd != 0 & Year == year),
         #           colour = "red", alpha = 0.6, size = 0.1) +
         # 
         #   scale_x_continuous(limits = c(-165,-48.5)) +
         #   scale_y_continuous(limits = c(25,70)) +
         #   labs(title = year, colour = 'Ecoregion', fill = 'Ecoregion') +
         #   ## presences in areas being used for analysis 
         #   geom_segment(data = spread, inherit.aes = F,
         #                aes(x = xmax, xend = xmin, y = presence_edge), colour = "red") +
         #   geom_segment(data = spread, inherit.aes = F,
         #                aes(x = xmax, xend = xmin, y = sample_edge), colour = "darkgrey") +
         #   labs(x = "", y = "") +
         #   theme(legend.position = "none")
      }
      ## else if no data for that year 
      else {
        p = bbs_all %>%
          ## filter to data from current year
          filter(Year == year) %>%
          ggplot() +
          geom_sf(inherit.aes = FALSE, data = countries) +
          geom_sf(size = 0.1, colour = "darkgrey", alpha = 0.6) +
          scale_x_continuous(limits = c(-165,-48.5)) +
          scale_y_continuous(limits = c(25,70)) +
          labs(title = year) +
          geom_sf(size = 0.1, data = filter(bbs_all, TotalAbd != 0 & Year == year), inherit.aes = FALSE,
                  colour = "red", alpha = 0.6, size = 0.1) +
          scale_colour_gradient(trans = "log", limits = c(1, max(bbs_all$TotalAbd, na.rm = T)))+
          labs(x = "", y = "")
        
        # p2 = bbs_all %>% ## filter to data from current year
        #   filter(Year == year) %>%
        #   ggplot(aes()) +
        #   geom_sf(inherit.aes = FALSE, data = countries) +
        #   geom_sf(data = lakes, mapping = aes(geometry = geometry), color = "black",
        #           fill = "white") +
        #   geom_sf(data = ecoregions, mapping = aes(geometry = geometry, fill = NA_L1NAME, colour = NA_L1NAME), 
        #           alpha = 0.25) +
        #   ## all routes 
        #   geom_sf(size = 0.5, colour = "darkgrey", alpha = 0.8, stroke = 0.1, shape = 4) +
        #   geom_sf(data = filter(bbs_all, TotalAbd != 0 & Year == year),
        #           colour = "red", alpha = 0.6, size = 0.1) +
        #   
        #   scale_x_continuous(limits = c(-165,-48.5)) +
        #   scale_y_continuous(limits = c(25,70)) +
        #   labs(title = year, colour = 'Ecoregion', fill = 'Ecoregion') +
        #   theme(legend.position = "none")
        
      }
    }
    ## if no range edges to be measured, plot presences anyways 
    else {
      p = bbs_all %>%
        ## filter to data from current year
        filter(Year == year) %>%
        ggplot() +
        geom_sf(inherit.aes = FALSE, data = countries) +
        geom_sf(size = 0.1, colour = "darkgrey", alpha = 0.6) +
        scale_x_continuous(limits = c(-165,-48.5)) +
        scale_y_continuous(limits = c(25,70)) +
        labs(title = year) +
        geom_sf(size = 0.1, data = filter(bbs_all, TotalAbd != 0 & Year == year), inherit.aes = FALSE,
                colour = "red", alpha = 0.6, size = 0.1) +
        scale_colour_gradient(trans = "log", limits = c(1, max(bbs_all$TotalAbd, na.rm = T))) +
        labs(x = "", y = "")
      
      # p2 = bbs_all %>% ## filter to data from current year
      #   filter(Year == year) %>%
      #   ggplot(aes()) +
      #   geom_sf(inherit.aes = FALSE, data = countries) +
      #   geom_sf(data = lakes, mapping = aes(geometry = geometry), color = "black",
      #           fill = "white") +
      #   geom_sf(data = ecoregions, mapping = aes(geometry = geometry, fill = NA_L1NAME, colour = NA_L1NAME), 
      #           alpha = 0.25) +
      #   ## all routes 
      #   geom_sf(size = 0.5, colour = "darkgrey", alpha = 0.8, stroke = 0.1, shape = 4) +
      #   geom_sf(data = filter(bbs_all, TotalAbd != 0 & Year == year),
      #           colour = "red", alpha = 0.6, size = 0.1) +
      #   scale_x_continuous(limits = c(-165,-48.5)) +
      #   scale_y_continuous(limits = c(25,70)) +
      #   labs(title = year, colour = 'Ecoregion', fill = 'Ecoregion') +
      #   theme(legend.position = "none")
    }
    ## save 
    ggsave(p, path = paste0(path, str_replace_all(genus_sp, " ", "\\_")),
           filename = paste0("BBS_", str_replace_all(genus_sp," ", "\\_"), "_",
                             year, ".png"), width = 6, height = 3)
    # ggsave(p2, path = paste0(path2, str_replace_all(genus_sp, " ", "\\_")), 
    #        filename = paste0("BBS_ecoregion_", str_replace_all(genus_sp," ", "\\_"), "_", 
    #                          year, ".png"), width = 6, height = 3)
  }
  sp = sp + 1
  
  print(paste0("On species number: ", sp))
}


## plot one for committee meeting document 

spread = spread(bbs_filtered, key = type, value = Latitude_of_edge)  %>%
  select(Longitude_bin, xmax, xmin, sample_edge, presence_edge) %>% 
  st_drop_geometry() %>%
  distinct()

temp =  bbs_filtered %>%
  group_by(Route, route, RouteName) %>%
  mutate(TotalAbd_max = max(TotalAbd)) %>%
  ungroup() %>%
  filter(TotalAbd != 0) %>%
  select(-TotalAbd) %>%
  distinct() 

example = bbs_filtered %>%
  ggplot(aes()) +
  geom_sf(inherit.aes = FALSE, data = countries) +
  geom_sf(data = lakes, mapping = aes(geometry = geometry), color = "black",
          fill = "white") +
  ## all routes
  #geom_sf(size = 0.5, colour = "darkgrey", alpha = 0.8, stroke = 0.1, shape = 4) +
  ## absences in areas being used for analysis
  geom_sf(size = 0.5, data = filter(bbs_all), inherit.aes = FALSE,
          stroke = 0.1, shape = 4, colour = "grey80") +
  geom_sf(size = 0.5, data = filter(bbs_filtered, TotalAbd == 0), inherit.aes = FALSE,
          stroke = 0.1, shape = 4, colour = "grey50") +
  scale_colour_viridis_c(trans = "log", breaks = c(1, 5, 50, 500)) +
  scale_x_continuous(limits = c(-165,-48.5)) +
  scale_y_continuous(limits = c(25,70)) +
  labs(colour = "Total\nabundance", x = "Longitude", y = "Latitude") +
  ## presences in areas being used for analysis
  geom_sf(size = 0.1, data = filter(temp), inherit.aes = FALSE,
          aes(colour = TotalAbd_max)) +
  geom_segment(data = spread, inherit.aes = F,
               aes(x = xmax, y = 60, yend = 25), colour = "black") +
  geom_segment(data = spread, inherit.aes = F,
               aes(x = xmin, y = 60, yend = 25), colour = "black") +
  theme(panel.grid = element_blank())

ggsave(example, path = "outputs/figures",
       filename = "example-for-committee-meeting.png", width = 6, height = 3)





## next improvements/steps:
## - find a way to split species ranges into longitudinal chunks in a way that makes sense for each distribution
## - for longitudinal chunks that have discontinuous bits of the range inside, get rid of more southernly part
## - add back non-songbirds 
## - rethink whether I need to remove years with less than 10 presences (benefits of keeping less certain data and having more versus removing)

## plot range edges and sampling edge over time for each species 
############################################
##         plot range edge shift          ##
############################################
bbs_sp <- read.csv("outputs/data-processed/BBS_spatial-filter.csv")

length(unique(bbs_sp$genus_sp))

## get some stats
## how many range edge observations remain per species per longitudinal bin?
stats = bbs_sp %>%
  select(AOU, Longitude_bin, Year) %>%
  distinct() %>%
  group_by(Longitude_bin, AOU) %>%
  tally() 

hist(stats$n) ## lots have 1, but lots also have 40-ish
max(stats$n) ## some have up to 43

## get rid of observations in bins with less than 10 years of sampling 
bbs_sp = bbs_sp %>%
  group_by(Longitude_bin, AOU) %>%
  mutate(n_obs_bin = length(unique(Year))) %>% 
  filter(n_obs_bin >= 10) %>% 
  select(-n_obs_bin)
  
## plot range edge shifts over time
r_edges = bbs_sp %>%
  select(AOU, genus_sp, Longitude_bin, Year, type, Latitude_of_edge) %>%
  distinct() 

r_edges %>%
  filter(genus_sp == "Catherpes mexicanus") %>%
  mutate(group = paste(AOU, type, Longitude_bin)) %>%
  ggplot(aes(x = Year, y = Latitude_of_edge, colour = type, group = group)) +
  geom_point(size = 0.5) +
  geom_smooth(method = 'lm', se = F, linewidth = 0.2) +
  facet_wrap(~Longitude_bin)

## measure shifts with lms
library(broom)
lms = r_edges %>%
  filter(type == "presence_edge") %>% # filter to range edges
  group_by(AOU, genus_sp, Longitude_bin) %>% # group by species and longitude bin
  do(tidy(lm(Latitude_of_edge ~ Year, data = .), conf.int = TRUE)) %>% 
  ungroup()

## plot slopes 
slopes = lms %>%
  filter(term == "Year") %>%
  mutate(estimate = estimate*111) ## convert to km/y

slopes %>%
  ggplot(aes(x = estimate)) +
  geom_histogram() +
  geom_vline(xintercept = 0) +
  geom_vline(aes(xintercept = mean(estimate)), colour = 'red')

slopes %>%
  group_by(Longitude_bin) %>%
  mutate(mean = mean(estimate)) %>%
  ungroup() %>%
  ggplot(aes(x = estimate)) +
  geom_histogram() +
  facet_wrap(~Longitude_bin) +
  geom_vline(xintercept = 0) +
  geom_vline(aes(xintercept = mean), colour = 'red')

## max
max(slopes$estimate)
## 26.1 km/y

## min
min(slopes$estimate)
## -25.5 km/y 

## look at the crazy ones
max = slopes[which(slopes$estimate == max(slopes$estimate)),]
r_edges %>%
  filter(AOU == max$AOU, Longitude_bin == max$Longitude_bin) %>%
  mutate(group = paste(AOU, type, Longitude_bin)) %>%
  ggplot(aes(x = Year, y = Latitude_of_edge, colour = type, group = group)) +
  geom_point(size = 0.5) +
  geom_smooth(method = 'lm', se = F, linewidth = 0.2) 
## house finch (introduced)
max$genus_sp
max$Longitude_bin

## NOTE: get rid of introduced species

min = slopes[which(slopes$estimate == min(slopes$estimate)),]
r_edges %>%
  filter(AOU == min$AOU, Longitude_bin == min$Longitude_bin) %>%
  mutate(group = paste(AOU, type, Longitude_bin)) %>%
  ggplot(aes(x = Year, y = Latitude_of_edge, colour = type, group = group)) +
  geom_point(size = 0.5) +
  geom_smooth(method = 'lm', se = F, linewidth = 0.2) 
min$genus_sp
min$Longitude_bin


## look at relationship between number of points used in regression and outlier-ness
r_edges %>%
  group_by(Longitude_bin, AOU) %>%
  mutate(n_obs_bin = length(unique(Year))) %>% 
  ungroup() %>%
  select(n_obs_bin, AOU, Longitude_bin) %>%
  distinct() %>%
  left_join(., slopes) %>%
  ggplot(aes(x = n_obs_bin, y = estimate)) +
  geom_point(size = 0.5) +
  geom_smooth(method = 'lm', se = F, linewidth = 0.2) 




## read in latitudinal climate velocity 
cv_lat <- rast("outputs/data-processed/cv_north-america_1980-2021.tif")

## crop to longitudinal bounds 
ext = ext(-180, -50, 18, 84)

xmins = seq(from = -180, to = -60, by = 10)
xmaxs = seq(from = -170, to = -50, by = 10)

polys = c()
for(i in 1:length(xmins)) {
  polys = st_as_sf(vect(ext(xmins[i], xmaxs[i], 18, 84))) %>%
    rbind(polys, .)
}

plot(polys)

st_crs(polys) = st_crs(cv_lat)

plot(cv_lat, add = T)

## calcualate mean climate velocity within each longitudinal bound 
cv_means = extract(cv_lat, polys, fun = mean, na.rm = TRUE)

polys$cv_mean = cv_means$`cv_north-america_1980-2021`

polys %>%
  ggplot() +
  geom_sf(aes(fill = cv_mean))

hist(cv_lat)
hist(slopes$estimate)

  