## script to analyze the monthly predictions from sdms of breeding bird survey
library(tidyverse)
library(terra)
library(tidyterra)

## things to measure:
## 1. shift of suitability edge 
## 2. autocorrelation of suitability
## 3. spatial synchrony of suitability
## 4. change in suitability over time 
## 5. standard deviation of suitability 

############################################
##          analyze suitability shift     ##
############################################
## make list of sdm output filenames/directories 
sdms = list.files("/Volumes/NIKKI/sdms", full.names = T)

files = c()
species = c()
for(i in 1:length(sdms)) {
  files = append(files, list.files(list.files(sdms[i], full.names = T)[7], full.names = T))
  species = append(species, str_split_fixed(sdms[i], "sdms/", 2)[,2])
}

## read in list of species with study divided into longitudinal blocks
bbs <- read.csv("outputs/data-processed/BBS_spatial-filter.csv")

## filter to species with occurrence data 
bbs <- filter(bbs, genus_sp %in% str_replace_all(species, "_", " "))
length(unique(bbs$genus_sp))

## get rid of absences outside of longitudes of presence
bbs <- filter(bbs, !is.na(Longitude_bin))

## for each species, measure aspects of suitability within each latitudinal block 
all_results = c()
sp = 1
while(sp <= length(files)) {
  ## get species name 
  cur_sp = species[sp]
  
  ## get longitudinal bins
  lon_bins <- bbs %>%
    filter(genus_sp == str_replace_all(cur_sp, "_", " ")) %>%
    select(Longitude_bin) %>%
    distinct() %>%
    mutate(min = str_split_fixed(Longitude_bin, ",", 2)[,1],
           max = str_split_fixed(Longitude_bin, ",", 2)[,2]) %>%
    mutate(min = as.numeric(as.character(str_remove(min, "\\("))),
           max = as.numeric(as.character(str_remove(max, "\\]"))), 
           bin = 1:nrow(.))
  
  ## get sdm suggested minimum training presence threshold
  maxent = read.csv(paste0("/Volumes/NIKKI/sdms/", cur_sp, "/maxentResults.csv"))
  threshold = maxent$Maximum.training.sensitivity.plus.specificity.area
  
  ## read in sdm predictions 
  pred = rast(files[sp])
  
  ## change extent
  pred = crop(pred, c(-170, -50, 24.5, 88))
  
  # ggplot() +
  #   geom_spatraster(data = pred[[104]]) +
  #   scale_fill_continuous(limits = c(0,1), na.value = "transparent")
  
  ## convert to xy data frame
  pred = as.data.frame(pred, xy = TRUE)
  
  ## assign each cell to a longitudinal bin
  bin = c()
  for(i in 1:nrow(pred)) {
    bin = append(bin, ifelse(length(which(pred$x[i] >= lon_bins$min & pred$x[i] <= lon_bins$max)) == 0, 
                             NA, 
                             which(pred$x[i] >= lon_bins$min & pred$x[i] <= lon_bins$max)))
  }
  pred$bin = bin
  pred = left_join(pred, lon_bins, by = "bin")

  ## get rid of data outside of longitudinal bands
  pred <- pred[!is.na(pred$Longitude_bin),]
  
  ## split by longitudinal bin and calculate suitability edge in each month 
  pred = pred %>%
    arrange(Longitude_bin) 
  
  groups <- pred %>%
    group_by(Longitude_bin) %>%
    group_split()
  
  p95 <- lapply(groups, FUN = function(x) {
    p95 = c()
    for(i in 3:(ncol(x)-4)) {
      cur = as.data.frame(x[,c(1:2, i)])
      cur = cur[which(!is.na(cur[,3])),]
      p95 = append(p95, modi::weighted.quantile(cur$y, w = cur[,3], 0.95))
    }
    return(p95)
  })
  
  p75 <- lapply(groups, FUN = function(x) {
    p75 = c()
    for(i in 3:(ncol(x)-4)) {
      cur = as.data.frame(x[,c(1:2, i)])
      cur = cur[which(!is.na(cur[,3])),]
      p75 = append(p75, modi::weighted.quantile(cur$y, w = cur[,3], 0.75))
    }
    return(p75)
  })
  
  ## use sdm suggested threshold:
  pthresh <- lapply(groups, FUN = function(x) {
    pthresh = c()
    for(i in 3:(ncol(x)-4)) {
      cur = as.data.frame(x[,c(1:2, i)])
      cur = cur[which(!is.na(cur[,3])),]
      ## get rid of cells with suitability less than threshold
      cur = cur[cur[,3] >= threshold,]
      ## calculate max latitude within bin that falls within threshold 
      pthresh = append(pthresh, max(cur$y))
    }
    return(pthresh)
  })
  
  ## calculate across entire distribution
  pre = pred %>%
    select(-bin, -Longitude_bin, -min, -max) 
    
  p95_sp = c()
  for(i in 3:(ncol(pred)-4)) {
    cur = as.data.frame(pred[,c(1:2, i)])
    cur = cur[which(!is.na(cur[,3])),]
    p95_sp = append(p95_sp, modi::weighted.quantile(cur$y, w = cur[,3], 0.95))
  }
  
  ## make dataframe of results
  results = data.frame(genus_sp = rep(cur_sp, length(unlist(p95))),
                       date = rep(colnames(pred)[3:(ncol(pred)-4)], length(groups)),
                       year = rep(as.numeric(as.character(paste0(rep(1966:2024, each = 4), ".", 
                                                                 rep(4:7, length.out = 236)))), length(groups)),
             Longitude_bin = rep(unique(pred$Longitude_bin), each = length(colnames(pred)[3:(ncol(pred)-4)])),
             p95 = unlist(p95),
             p95_sp = rep(p95_sp, n = length(groups)),
             p75 = unlist(p75),
             threshold = rep(threshold, length.out = length(unlist(p95))),
             pthresh = unlist(pthresh))
  
  ## plot
  # results %>%
  #   ggplot(aes(x = year, y = p95, colour = Longitude_bin)) +
  #   geom_line()
  # 
  # results %>%
  #   filter(!is.infinite(pthresh)) %>%
  #   ggplot(aes(x = year, y = pthresh, colour = Longitude_bin)) +
  #   geom_line()
  #
  # results %>%
  #   select(year, p95_sp) %>%
  #   distinct() %>% 
  #   ggplot(aes(x = year, y = p95_sp)) +
  #   geom_line() +
  #   geom_smooth(method = "lm")
  
  ## join to all results
  all_results = rbind(all_results, results)
  
  ## go to next species 
  sp = sp + 1
  print(sp)
}

all_results = all_results %>%
  mutate(min = str_split_fixed(Longitude_bin, ",", 2)[,1],
         max = str_split_fixed(Longitude_bin, ",", 2)[,2]) %>%
  mutate(min = as.numeric(as.character(str_remove(min, "\\("))),
         max = as.numeric(as.character(str_remove(max, "\\]")))) %>%
  mutate(Year =  as.numeric(as.character(substr(year, 1, 4))))

## save 
write.csv(all_results, "outputs/data-processed/suitability_edges.csv", row.names = F)


## plot against range shift rates in each longitudinal band
all_results <- read.csv("outputs/data-processed/suitability_edges.csv")

bbs_edges <- bbs %>%
  filter(keep == 1) %>%
  spread(key = type, value = Latitude_of_edge) 

bbs_edges = bbs_edges %>%
  group_by(genus_sp, AOU, Year, Longitude_bin) %>% ## group by year and species and longitude bin
  filter(TotalAbd != 0) %>% ## get rid of absences 
  filter(n() >= 30) %>% ## filter to species x longitude bins x years with more than 30 observations 
  ## calculate range edge metrics for each year for each species:
  mutate(p95_w = modi::weighted.quantile(Latitude, w = TotalAbd, 0.95),
         p95 = quantile(Latitude, 0.95)) %>%
  select(genus_sp, AOU, Year, Longitude_bin, sample_edge, presence_edge, p95, p95_w, TotalAbd, Latitude) %>%
  distinct() %>%
  group_by(genus_sp, AOU, Year) %>% ## group by year and species
  ## calculate range edge metrics across whole range edge:
  mutate(p95_sp_w = modi::weighted.quantile(Latitude, w = TotalAbd, 0.95),
         p95_sp = quantile(Latitude, 0.95)) %>%
  select(-TotalAbd, -Latitude) %>%
  distinct() %>%
  group_by(genus_sp, Longitude_bin) %>% ## group by species and latitudinal bin 
  filter(n() >= 10) %>% ## filter to longitudinal bins with more than 10 years of range edge estimates 
  gather(key = "range_edge", value = "latitude", c(p95_w, p95, p95_sp_w, p95_sp)) %>%
  ungroup() %>%
  mutate(min = str_split_fixed(Longitude_bin, ",", 2)[,1],
         max = str_split_fixed(Longitude_bin, ",", 2)[,2]) %>%
  mutate(min = as.numeric(as.character(str_remove(min, "\\("))),
         max = as.numeric(as.character(str_remove(max, "\\]"))))

length(unique(bbs_edges$genus_sp)) ## 162

## intraspecific
bbs_edges %>%
 # filter(genus_sp == "Cardinalis cardinalis") %>%
  mutate(east_west = ifelse(bbs_edges$max <= -101, "West", 
                            "East")) %>%
  filter(range_edge == "p95") %>%
  mutate(group = paste(genus_sp, Longitude_bin)) %>%
  ggplot(aes(x = Year, y = latitude, colour = range_edge, group = group)) +
  geom_line() +
  geom_smooth(method = "lm", se = F, colour = "black", size = 0.1) +
  labs(x = "Year", y = "Latitude") +
  facet_wrap(~east_west)

## overall trend
bbs_edges %>%
  # filter(genus_sp == "Cardinalis cardinalis") %>%
  mutate(east_west = ifelse(bbs_edges$max <= -101, "West", 
                            "East")) %>%
  filter(range_edge == "p95_sp") %>%
  mutate(group = paste(genus_sp, Longitude_bin)) %>%
  ggplot(aes(x = Year, y = latitude, colour = range_edge, group = group)) +
  geom_line() +
  geom_smooth(method = "lm", se = F, colour = "black", size = 0.1) +
  labs(x = "Year", y = "Latitude") +
  facet_wrap(~east_west)


## next measure shifts using lms

## fit slopes to suitability shifts 
library(broom)

## measure over same years that bbs estimates of range edge are available
key = bbs_edges %>%
  group_by(genus_sp, Longitude_bin) %>%
  summarize(min_year = min(Year),
            max_year = max(Year))

all_results = all_results %>% 
  mutate(genus_sp = str_replace_all(genus_sp, "_", " ")) %>%
  left_join(., key) %>%
  filter(year >= min_year & year <= max_year)

## intraspecific
suit_trend <- all_results %>%
  group_by(genus_sp, Longitude_bin, min, max) %>%
  do(tidy(lm(data = ., p95 ~ year))) %>%
  filter(term == "year")

edge_trend = bbs_edges %>%
  filter(range_edge == "p95") %>%
  group_by(genus_sp, Longitude_bin, min, max) %>% 
  filter(n() >= 10) %>% ## filter to only groups with 10 or more observations
  do(tidy(lm(data = ., latitude ~ Year))) %>%
  filter(term == "Year")

## overall trend
suit_trend_sp <- all_results %>%
 # filter(str_detect(year, ".6")) %>% ## ONLY JUNE
  select(genus_sp, p95_sp, year) %>%
  distinct() %>%
  group_by(genus_sp) %>%
  do(tidy(lm(data = ., p95_sp ~ year))) %>%
  filter(term == "year")

edge_trend_sp = bbs_edges %>%
  filter(range_edge == "p95_sp") %>%
  select(genus_sp, range_edge, latitude, Year) %>%
  distinct() %>%
  group_by(genus_sp) %>% 
  filter(n() >= 10) %>% ## filter to only groups with 10 or more observations
  do(tidy(lm(data = ., latitude ~ Year))) %>%
  filter(term == "Year")

## plot the range edge shift and suitability shift in each longitude band for a given species:
species = unique(bbs_edges$genus_sp)[32]
bbs_edges %>%
  filter(genus_sp == species) %>%
  filter(range_edge == "p95") %>%
  ggplot(aes(x = Year, y = latitude, colour = range_edge)) +
  geom_line() +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Year", y = "Latitude", title = species) +
  # geom_line(data = all_results %>% filter(genus_sp == "Cardinalis_cardinalis"), 
  #           aes(x = year, y = p75), inherit.aes = F) +
  # geom_smooth(data = all_results %>% filter(genus_sp == "Cardinalis_cardinalis"), 
  #             method = "lm", aes(x = year, y = p75), inherit.aes = F, colour = "black", se = F) +
  geom_line(data = all_results %>% filter(genus_sp == str_replace(species, " ", " "), year >= 1980),
            aes(x = year, y = p95), inherit.aes = F) +
  geom_smooth(data = all_results %>% filter(genus_sp == str_replace(species, " ", " "), year >= 1980),
              method = "lm", aes(x = year, y = p95), inherit.aes = F, colour = "black", se = F) +
  # geom_line(data = df, aes(x = year, y = y), inherit.aes = F) +
  # geom_line(data = df2, aes(x = year, y = y), inherit.aes = F, colour = "red") +
  facet_wrap(~min)  + ## this orders from west to east 
  scale_y_continuous(limits = c(27, 83))

## overall trend
bbs_edges %>%
  filter(genus_sp == species) %>%
  filter(range_edge == "p95_sp") %>%
  ggplot(aes(x = Year, y = latitude, colour = range_edge)) +
  geom_line() +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Year", y = "Latitude", title = species) +
  # geom_line(data = all_results %>% filter(genus_sp == "Cardinalis_cardinalis"), 
  #           aes(x = year, y = p75), inherit.aes = F) +
  # geom_smooth(data = all_results %>% filter(genus_sp == "Cardinalis_cardinalis"), 
  #             method = "lm", aes(x = year, y = p75), inherit.aes = F, colour = "black", se = F) +
  geom_line(data = all_results %>% filter(genus_sp == str_replace(species, " ", " "), year >= 1980),
            aes(x = year, y = p95_sp), inherit.aes = F) +
  geom_smooth(data = all_results %>% filter(genus_sp == str_replace(species, " ", " "), year >= 1980),
              method = "lm", aes(x = year, y = p95_sp), inherit.aes = F, colour = "black", se = F) +
  # geom_line(data = df, aes(x = year, y = y), inherit.aes = F) +
  # geom_line(data = df2, aes(x = year, y = y), inherit.aes = F, colour = "red") +
  scale_y_continuous(limits = c(27, 83))

bbs_edges %>%
  #filter(genus_sp == species) %>%
  filter(range_edge == "p95") %>%
  left_join(., all_results) %>% 
  ggplot(aes(x = p95, y = latitude)) +
  geom_point() +
  labs(y = "Lat of range edge", x = "Lat of suitability edge") +
  geom_abline(intercept = 0, slope = 1, colour = 'red')

suit = suit_trend %>%
  rename("suit_trend" = estimate, "suit_error" = std.error) %>%
  select(genus_sp, Longitude_bin, suit_trend, suit_error, min, max) 

edge = edge_trend %>%
  rename("edge_trend" = estimate, "edge_error" = std.error) %>%
  select(genus_sp, Longitude_bin, edge_trend, edge_error, min, max) %>%
  mutate(genus_sp = str_replace_all(genus_sp, "\\ ", " ")) 

suit_sp = suit_trend_sp %>%
  rename("suit_trend_sp" = estimate, "suit_error_sp" = std.error) %>%
  select(genus_sp, suit_trend_sp, suit_error_sp) 

edge_sp = edge_trend_sp %>%
  rename("edge_trend_sp" = estimate, "edge_error_sp" = std.error) %>%
  select(genus_sp,edge_trend_sp, edge_error_sp) %>%
  mutate(genus_sp = str_replace_all(genus_sp, "\\ ", " ")) 

trends = left_join(edge, suit)
trends_sp = left_join(edge_sp, suit_sp)

## convert to km per year 
trends$edge_trend = trends$edge_trend*111
trends$suit_trend = trends$suit_trend*111
trends_sp$edge_trend_sp = trends_sp$edge_trend_sp*111
trends_sp$suit_trend_sp = trends_sp$suit_trend_sp*111

## plot suitability shift versus range edge shift trend against each other
trends %>%
  ggplot(aes(x = suit_trend, y = edge_trend)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(method = "lm") 

trends_sp %>%
  ggplot(aes(x = suit_trend_sp, y = edge_trend_sp)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(method = "lm") 

trends %>%
  gather(key = "type", value = "trend", c(suit_trend, edge_trend)) %>%
  ggplot(aes(x = type, y = trend)) +
  geom_violin()

trends_sp %>%
  gather(key = "type", value = "trend", c(suit_trend_sp, edge_trend_sp)) %>%
  ggplot(aes(x = type, y = trend)) +
  geom_violin()

mean(trends$suit_trend, na.rm = T) ## mean suitability edge shift: 2.2 km/y N
mean(trends$edge_trend, na.rm = T) ## mean edge shift: 1.3 km/y N
mean(trends$suit_trend - trends$edge_trend, na.rm = T) ## mean difference: 0.92 km/y N

mean(trends_sp$suit_trend_sp, na.rm = T) ## mean suitability edge shift: 2.2 km/y N
mean(trends_sp$edge_trend_sp, na.rm = T) ## mean edge shift: 2.6 km/y N
mean(trends_sp$suit_trend_sp - trends_sp$edge_trend_sp, na.rm = T) ## mean difference: -0.43 km/y N

trends %>%
  gather(key = "type", value = "trend", c(suit_trend, edge_trend)) %>%
  ggplot(aes(x = trend, fill = type)) +
  geom_histogram(bins = 50, position = position_dodge())

trends_sp %>%
  gather(key = "type", value = "trend", c(suit_trend_sp, edge_trend_sp)) %>%
  ggplot(aes(x = trend, fill = type)) +
  geom_histogram(bins = 50, position = position_dodge())

hist(trends$suit_trend - trends$edge_trend)
  
trends %>%
  mutate(east_west = ifelse(max <= -101, "West", 
                            "East")) %>%
  gather(key = "type", value = "trend", c(suit_trend, edge_trend)) %>%
  ggplot(aes(x = trend, fill = east_west)) +
  geom_histogram(bins = 50, position = position_dodge())

trends %>%
  select(genus_sp, max) %>% 
  distinct() %>%
  left_join(trends_sp, .) %>%
  mutate(east_west = ifelse(max <= -101, "West", 
                            "East")) %>%
  gather(key = "type", value = "trend", c(suit_trend_sp, edge_trend_sp)) %>%
  ggplot(aes(x = trend, fill = east_west)) +
  geom_histogram(bins = 50, position = position_dodge())



## need to get rid of non-native specices 
# - Streptopelia decaocto

## also need to make stricter too close to sampling edge criteria 

## to measure 1 per species, need to come up with a criteria of how to define edge 




## test cross correlation

edge_cc = bbs_edges %>%
  filter(range_edge == "p95",
         genus_sp == "Pipilo erythrophthalmus", 
         min == "-105.2") %>%
  distinct()

suit_cc = all_results %>%
  filter(genus_sp == "Pipilo erythrophthalmus", 
         min == "-105.2") %>%
  filter(str_detect(year, ".6")) %>%
  mutate(Year = as.integer(as.character(str_split_fixed(year, "\\.", 2)[,1])))

cc = left_join(edge_cc, suit_cc) %>%
  filter(!is.na(.$latitude), !is.na(.$p95))

library(testcorr)
cc.test(cc$p95, cc$latitude, max.lag = 5)
## negative lag: Y leads X
## positive lag: X leads Y
## lag 0: instantanous coupling


test = seewave::coh(cc$p95, cc$latitude, f = 1)












