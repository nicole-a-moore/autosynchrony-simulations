## performing stratified permutation on BBS data 
library(tidyverse)
library(sf)

## idea:
# -keep the set of sampling events (years × route locations) fixed
# -randomly shuffle presence/absence among years within latitude strata
# -measure range edge shifts and generate a null distribution 
# -compare to observed range edge shift
## this preserves the latitudinal distribution of survey effort in each year and the overall prevalence within latitudinal strata, while removing any temporal signal in species occurrence
## e.g., if high-latitude routes were sampled more in recent years (as is true), the permutation preserves this 

## read in clean bbs data
bbs =  read.csv("outputs/data-processed/bbs_clean-subset.csv")

## filter to northern cardinal for now
nc <- filter(bbs, AOU == 5930) 

## northern cardinal 5930
## bc chickadee 7350 

## calculate range shift metrics
## leading edge: mean latitude of 5 northernmost presences, abd weighted 95th quantile
## trailing edge: mean latitude of 5 southernmost presences, qbs weighted 5th quantile
mean_max_lat = nc %>%
  filter(!TotalAbd == 0 & !is.na(TotalAbd)) %>% ## get rid of absences 
  group_by(genus_sp, AOU, Year) %>% ## group by year and species 
  ## calculate range edge metrics for each year for each species:
  arrange(-Latitude) %>%
  slice_head(n = 5) %>% 
  summarise(mean_max_lat = mean(Latitude), 
            max_lat = max(Latitude)) %>% 
  st_drop_geometry()

mean_min_lat = nc %>%
  filter(!TotalAbd == 0 & !is.na(TotalAbd)) %>% ## get rid of absences 
  group_by(genus_sp, AOU, Year) %>% ## group by year and species 
  ## calculate range edge metrics for each year for each species:
  arrange(Latitude) %>%
  slice_head(n = 5) %>% 
  summarise(mean_min_lat = mean(Latitude),
            min_lat = min(Latitude)) %>% 
  st_drop_geometry()

quantiles <-  nc %>%
  filter(!TotalAbd == 0 & !is.na(TotalAbd)) %>% ## get rid of absences 
  group_by(genus_sp, AOU, Year) %>% ## group by year and species 
  ## calculate range edge metrics for each year for each species:
  summarise(p95_w = modi::weighted.quantile(Latitude, w = TotalAbd, 0.95),
            p5_w = modi::weighted.quantile(Latitude, w = TotalAbd, 0.05),
            p95 = quantile(Latitude, 0.95), 
            p5 = quantile(Latitude, 0.05)) %>% ## calculate range edge metrics for that year       
  ungroup() %>% ## ungroup
  st_drop_geometry() ## drop geometry

range_shift_sum <- left_join(quantiles, mean_max_lat) %>%
  left_join(., mean_min_lat)

## plot range edges over time 
range_shift_sum %>%
  gather(key = "range_edge_metric", value = "latitude", c(p95_w, p95, p5_w, p5, mean_min_lat, mean_max_lat, max_lat, min_lat)) %>%
  mutate(range_edge = ifelse(range_edge_metric %in% c("max_lat", "p95", "p95_w", "mean_max_lat"), "Leading edge", "Trailing edge")) %>%
  ggplot(aes(x = Year, y = latitude, colour = range_edge_metric, group = range_edge_metric)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", se = F, size = 0.5) +
  facet_wrap(~range_edge) +
  scale_colour_discrete(labels = c("Maximum latitude", "Mean latitude of 5 most Northern occurrences", "Mean latitude of 5 most Southern occurrences",
                                 "Minimum latitude", "5th percentile of presence",  "Abundance-weighted 5th percentile of presence",
                                 "95th percentile of presence",  "Abundance-weighted 95th percentile of presence")) +
  labs(colour = "Range edge definition")
## :-)

## now do the same, but for sites sampled over time 
## plot latitude of sampled sites over time
nc %>%
  select(Latitude, Longitude, Year, TotalAbd) %>%
  distinct() %>%
  mutate(is_present = TotalAbd != 0) %>%
  ggplot(aes(x = Year, y = Latitude, colour = is_present)) +
  geom_point(alpha = 0.1) +
  theme_bw() +
  labs(colour = "") +
  scale_colour_discrete(labels = c("Absence", "Presence"))

## calculate shift in sampling of sites over time 
mean_max_lat_samples = nc %>%
  ## keep all sites (even where abundance was 0)
  group_by(genus_sp, AOU, Year) %>% ## group by year and species 
  ## calculate range edge metrics for each year for each species:
  arrange(-Latitude) %>%
  slice_head(n = 5) %>% 
  summarise(mean_max_lat = mean(Latitude), 
            max_lat = max(Latitude)) %>% 
  st_drop_geometry()

mean_min_lat_samples = nc %>%
  filter(!TotalAbd == 0 & !is.na(TotalAbd)) %>% ## get rid of absences 
  group_by(genus_sp, AOU, Year) %>% ## group by year and species 
  ## calculate range edge metrics for each year for each species:
  arrange(Latitude) %>%
  slice_head(n = 5) %>% 
  summarise(mean_min_lat = mean(Latitude),
            min_lat = min(Latitude)) %>% 
  st_drop_geometry()

quantiles_samples <- nc %>%
  group_by(genus_sp, AOU, Year) %>% ## group by year and species 
  ## calculate range edge metrics for each year for each species:
  summarise(p95 = quantile(Latitude, 0.95), 
            p5 = quantile(Latitude, 0.05)) %>% ## calculate range edge metrics for that year       
  ungroup() %>% ## ungroup
  st_drop_geometry() ## drop geometry

range_shift_sum_samples <- left_join(quantiles_samples, mean_max_lat_samples) %>%
  left_join(., mean_min_lat_samples)

## prepare to join dataframes 
range_shift_sum <- range_shift_sum %>%
  mutate(type = "range_edge") %>%
  gather(key = "range_edge_metric", value = "latitude", c(p95_w, p95, p5_w, p5, mean_min_lat, mean_max_lat, max_lat, min_lat)) %>%
  mutate(range_edge = ifelse(range_edge_metric %in% c("max_lat", "p95", "p95_w", "mean_max_lat"), "Leading edge", "Trailing edge")) 
  
range_shift_sum_samples <- range_shift_sum_samples %>%
  mutate(type = "sample_edge") %>%
  gather(key = "range_edge_metric", value = "latitude", c(p95, p5, mean_min_lat, mean_max_lat, max_lat, min_lat)) %>%
  mutate(range_edge = ifelse(range_edge_metric %in% c("max_lat", "p95", "p95_w", "mean_max_lat"), "Leading edge", "Trailing edge")) 

## plot sample shift together with range edge shift 
# New facet label names for supp variable
labs <- c("Maximum latitude","Mean latitude of 5 most\nNorthern occurrences","95th percentile of presence","Abundance-weighted 95th\npercentile of presence")
names(labs) <- c("max_lat", "mean_max_lat", "p95", "p95_w")

rbind(range_shift_sum, range_shift_sum_samples) %>%
  filter(range_edge == "Leading edge") %>%
  ggplot(aes(x = Year, y = latitude, colour = type, group = type)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", se = F, size = 0.5) +
  facet_wrap(~ range_edge_metric, labeller = labeller(range_edge_metric = labs)) +
  theme() +
  labs(fill = "") +
  scale_colour_discrete(labels = c("Range edge", "Sampling edge"))

## plot 5 max occurrences 
max5_occ <- nc %>%
  filter(!TotalAbd == 0 & !is.na(TotalAbd)) %>% ## get rid of absences 
  group_by(genus_sp, AOU, Year) %>% ## group by year and species 
  ## calculate range edge metrics for each year for each species:
  arrange(-Latitude) %>%
  slice_head(n = 5) %>%
  mutate(type = "Range edge")

max5_samp <- nc %>%
  group_by(genus_sp, AOU, Year) %>% ## group by year and species 
  ## calculate range edge metrics for each year for each species:
  arrange(-Latitude) %>%
  slice_head(n = 5) %>%
  mutate(type = "Sampling edge")

rbind(max5_occ, max5_samp) %>% 
  ggplot(aes(x = Year, y = Latitude, colour = type, group = type)) +
  geom_point(size = 0.5) +
  theme() +
  labs(fill = "") +
  scale_colour_discrete(labels = c("Range edge", "Sampling edge"))







## GARBAGE
##########################
## define latitude strata 
## start with every 5 degrees
nc$lat_bin <- cut(nc$Latitude, 
                   breaks = seq(from = 0, to = 90, by = 5))

nc %>%
  ggplot(aes(x = lat_bin)) + geom_bar()

## randomly shuffle presence/absence among years within latitude strata
nc = nc %>%
  group_by(lat_bin, Year) %>%
  mutate(TotalAbd_null = sample(TotalAbd, size = n(), replace = F))

## scrambled latitudinal presences:
nc %>%
  select(Latitude, Longitude, Year, TotalAbd_null) %>%
  distinct() %>%
  mutate(is_present = TotalAbd_null != 0) %>%
  filter(is_present) %>%
  ggplot(aes(x = Year, y = Latitude, colour = is_present)) +
  geom_point(alpha = 0.1) +
  theme_bw()

## versus real presences:
nc %>%
  select(Latitude, Longitude, Year, TotalAbd) %>%
  distinct() %>%
  mutate(is_present = TotalAbd != 0) %>%
  filter(is_present) %>%
  ggplot(aes(x = Year, y = Latitude, colour = is_present)) +
  geom_point(alpha = 0.1) +
  theme_bw()

## calculate null range edge shift over time 
mean_max_lat = nc %>%
  filter(!TotalAbd_null == 0 & !is.na(TotalAbd_null)) %>% ## get rid of absences 
  group_by(genus_sp, AOU, Year) %>% ## group by year and species 
  ## calculate range edge metrics for each year for each species:
  arrange(-Latitude) %>%
  slice_head(n = 5) %>% 
  summarise(mean_max_lat = mean(Latitude), 
            max_lat = max(Latitude)) %>% 
  st_drop_geometry()

mean_min_lat = nc %>%
  filter(!TotalAbd_null == 0 & !is.na(TotalAbd_null)) %>% ## get rid of absences 
  group_by(genus_sp, AOU, Year) %>% ## group by year and species 
  ## calculate range edge metrics for each year for each species:
  arrange(Latitude) %>%
  slice_head(n = 5) %>% 
  summarise(mean_min_lat = mean(Latitude),
            min_lat = min(Latitude)) %>% 
  st_drop_geometry()

quantiles <-  nc %>%
  filter(!TotalAbd_null == 0 & !is.na(TotalAbd_null)) %>% ## get rid of absences 
  group_by(genus_sp, AOU, Year) %>% ## group by year and species 
  ## calculate range edge metrics for each year for each species:
  summarise(p95_w = modi::weighted.quantile(Latitude, w = TotalAbd_null, 0.95),
            p5_w = modi::weighted.quantile(Latitude, w = TotalAbd_null, 0.05),
            p95 = quantile(Latitude, 0.95), 
            p5 = quantile(Latitude, 0.05)) %>% ## calculate range edge metrics for that year       
  ungroup() %>% ## ungroup
  st_drop_geometry() ## drop geometry

range_shift_sum_null <- left_join(quantiles, mean_max_lat) %>%
  left_join(., mean_min_lat)

range_shift_sum_null$type = "null"
range_shift_sum$type = "observed"

## plot range edges over time 
range_shift_sum_null %>%
  rbind(., range_shift_sum) %>%
  gather(key = "range_edge_metric", value = "latitude", c(p95_w, p95, p5_w, p5, mean_min_lat, mean_max_lat, max_lat, min_lat)) %>%
  mutate(range_edge = ifelse(range_edge_metric %in% c("max_lat", "p95", "p95_w", "mean_max_lat"), "Leading edge", "Trailing edge")) %>%
  ggplot(aes(x = Year, y = latitude, colour = range_edge_metric, group = range_edge_metric)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", se = F, size = 0.2) +
  facet_wrap(~type) 
## :-)


## try scrambling across all latitudes each year instead of binning
## want to see if the increase in sampling effort at higher latitudes over time is causing increased presence at high latitudes over time
## randomly shuffle presence/absence among years within latitude strata
nc = nc %>%
  group_by(Year) %>%
  mutate(TotalAbd_null = sample(TotalAbd, size = n(), replace = F))

## scrambled latitudinal presences:
nc %>%
  select(Latitude, Longitude, Year, TotalAbd_null) %>%
  distinct() %>%
  mutate(is_present = TotalAbd_null != 0) %>%
  filter(is_present) %>%
  ggplot(aes(x = Year, y = Latitude, colour = is_present)) +
  geom_point(alpha = 0.1) +
  theme_bw()

## calculate null range edge shift over time 
mean_max_lat = nc %>%
  filter(!TotalAbd_null == 0 & !is.na(TotalAbd_null)) %>% ## get rid of absences 
  group_by(genus_sp, AOU, Year) %>% ## group by year and species 
  ## calculate range edge metrics for each year for each species:
  arrange(-Latitude) %>%
  slice_head(n = 5) %>% 
  summarise(mean_max_lat = mean(Latitude), 
            max_lat = max(Latitude)) %>% 
  st_drop_geometry()

mean_min_lat = nc %>%
  filter(!TotalAbd_null == 0 & !is.na(TotalAbd_null)) %>% ## get rid of absences 
  group_by(genus_sp, AOU, Year) %>% ## group by year and species 
  ## calculate range edge metrics for each year for each species:
  arrange(Latitude) %>%
  slice_head(n = 5) %>% 
  summarise(mean_min_lat = mean(Latitude),
            min_lat = min(Latitude)) %>% 
  st_drop_geometry()

quantiles <-  nc %>%
  filter(!TotalAbd_null == 0 & !is.na(TotalAbd_null)) %>% ## get rid of absences 
  group_by(genus_sp, AOU, Year) %>% ## group by year and species 
  ## calculate range edge metrics for each year for each species:
  summarise(p95_w = modi::weighted.quantile(Latitude, w = TotalAbd_null, 0.95),
            p5_w = modi::weighted.quantile(Latitude, w = TotalAbd_null, 0.05),
            p95 = quantile(Latitude, 0.95), 
            p5 = quantile(Latitude, 0.05)) %>% ## calculate range edge metrics for that year       
  ungroup() %>% ## ungroup
  st_drop_geometry() ## drop geometry

range_shift_sum_null <- left_join(quantiles, mean_max_lat) %>%
  left_join(., mean_min_lat)

range_shift_sum_null$type = "null"
range_shift_sum$type = "observed"

## plot range edges over time 
range_shift_sum_null %>%
  rbind(., range_shift_sum) %>%
  gather(key = "range_edge_metric", value = "latitude", c(p95_w, p95, p5_w, p5, mean_min_lat, mean_max_lat, max_lat, min_lat)) %>%
  mutate(range_edge = ifelse(range_edge_metric %in% c("max_lat", "p95", "p95_w", "mean_max_lat"), "Leading edge", "Trailing edge")) %>%
  ggplot(aes(x = Year, y = latitude, colour = range_edge_metric, group = range_edge_metric)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", se = F, size = 0.2) +
  facet_wrap(~type) 
## :-)


nc %>%
  group_by(Year, Latitude) %>%
  mutate(total_abundance = sum(TotalAbd, na.rm = T)) %>%
  ungroup() %>%
  filter(total_abundance !=0) %>%
  ggplot(aes(x = Latitude, y = total_abundance)) +
  geom_point() +
  facet_wrap(~Year) 




