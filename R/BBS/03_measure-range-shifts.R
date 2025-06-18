## measure range edge shifts for each species
## note: need to ensure later that range edges are actually range edges 
library(tidyverse)
library(sf)
library(terra)
select = dplyr::select

## read in clean bbs data
bb =  read.csv("outputs/data-processed/bbs_clean-subset.csv")

## get canada and us map for plotting
countries <- necountries::countries(c("Canada", "United States of America"), part = TRUE)
countries <- st_union(countries)
# countries <- vect(countries)

## for each species:
## - plot range edges over time
## - measure range edge shifts over time

## plot distribution of latitude of sampled sites over time
bb %>%
  filter(!is.na(TotalAbd)) %>% 
  #mutate(year = ifelse(Year <= 1993, 1993, Year)) %>%
  select(Latitude, Longitude, Year) %>%
  distinct() %>%
  ggplot(aes(x = Year, y = Latitude)) +
  geom_point()
## 1980 is when highish-latitude sites (>55 deg N) begin being sampled consistently
## 1993 is when high-latitude sites (>65 deg N) begin being sampled consistently
## to avoid detecting range shift because of shift in sampling latitude, group observations before 1980/1993 together and begin measuring shift from there
## if species' leading range edge occurs at less than 55 deg N in 1980, group obs before 1980 and start ts there
## otherwise if species' leading range edge occurs at greater than 55 deg N in 1980, group obs before 19993 and start ts there

## even better: make sure there was enough sampling within the longitudinal domain of the range   

## for each year, calculate mean latitude of 5 most northern samples within the longitudinal bounds of each species for data before 1980
range_shift_sum <- bb %>%
  filter(!TotalAbd == 0 & !is.na(TotalAbd))## get rid of absence data and NAs (i.e., years where routes weren't sampled)

lon_edges <- range_shift_sum %>%
  filter(Year <= 1980) %>%
  group_by(genus_sp) %>%
  summarise(min_lon = quantile(Longitude, 0.05), 
            max_lon = quantile(Longitude, 0.95)) 

lat_edges = range_shift_sum %>%
  filter(Year <= 1980) %>%
  group_by(genus_sp) %>%
  arrange(-Latitude) %>%
  slice_head(n = 5) %>%
  summarize(n_edge = mean(Latitude))

## for each species, filter to samples within longitudinal bounds and calculate mean latitude of 5 highest latitude samples
sampling_sum <- bb %>%
  filter(Year <= 1980) %>%
  filter(!is.na(TotalAbd)) %>%
  select(Year, Latitude, Longitude) %>%
  distinct()

lats = c()
for(sp in 1:length(unique(lon_edges$genus_sp))) {
  ## filter to lat and lon mins 
  sampling_sum_sub <- sampling_sum %>%
    filter(Longitude >= lon_edges$min_lon[sp], Longitude <= lon_edges$max_lon[sp])
  
  ## filter to years before 1980, calculate mean latitude of 5 highest latitude samples
  mean_max_lat = sampling_sum_sub %>%
    filter(Year <= 1980) %>%
    arrange(-Latitude) %>%
    slice_head(n = 5) %>%
    summarize(mean_max_lat = mean(Latitude))
  
  ## calculate mean latitude of 5 highest latitude samples before 1980
  # mean_max_lat = sampling_sum_sub %>%
  #   summarise(n_edge_sample = quantile(Latitude, 0.95)) 
  
  ## save 
  lats = append(lats, as.numeric(mean_max_lat))
}
lon_edges$n_edge_sample = lats

## compare to latitudinal edge of sampling within longitudinal bounds 
n_edges <- left_join(lon_edges, lat_edges)

################################################
##          measure range edge shift          ##
################################################
range_shift_sum <- bb %>%
  filter(!TotalAbd == 0 & !is.na(TotalAbd))## get rid of absence data and NAs (i.e., years where routes weren't sampled)

range_shift_sum <- left_join(range_shift_sum, n_edges) 

range_shift_sum <- range_shift_sum %>%
  mutate(start_year = ifelse(n_edge + 1 <= n_edge_sample, 1980, 1993)) %>% ## decide what year to start at based on northern range edge position in 1980 relative to sampling 
  ## if latitudinal edge of sampling before 1980 is more 10 degree latitude higher than range edge, start measuring shift at 1980
  ## otherwise start at 1993
  mutate(year = ifelse(Year <= start_year, start_year, Year)) %>% ## combine data from start year and before together
  group_by(genus_sp, AOU, year) %>% ## group by year and species 
  ## calculate range edge metrics for each year for each species:
  summarise(p95_w = modi::weighted.quantile(Latitude, w = TotalAbd, 0.95),
            p5_w = modi::weighted.quantile(Latitude, w = TotalAbd, 0.05),
            p95 = quantile(Latitude, 0.95), 
            p5 = quantile(Latitude, 0.05)) %>% ## calculate range edge metrics for that year       
  ungroup() %>% ## ungroup
  st_drop_geometry() ## drop geometry

## calculate shift in 95 and 5 percentile latitude of sampled routes where each species was seen over time
## if there is a shift in 'range edges' of sampling, cannot be sure range shift observed is not caused by change in sampling latitude 
sampling_sum <- bb %>%
  filter(!is.na(TotalAbd)) 

sampling_sum <- left_join(sampling_sum, n_edges)

temp = sampling_sum[which(!sampling_sum$n_edge + 1 <= sampling_sum$n_edge_sample),] %>%
  select(genus_sp, n_edge, n_edge_sample, max_lon, min_lon) %>%
  distinct()

sampling_sum <- sampling_sum %>%
  mutate(start_year = ifelse(n_edge + 1 <= n_edge_sample, 1980, 1993)) %>% ## decide what year to start at based on northern range edge position in 1980 relative to sampling 
  ## if latitudinal edge of sampling before 1980 is more 1 degree latitude higher than range edge, start measuring shift at 1980
  ## otherwise start at 1993
  mutate(year = ifelse(Year <= start_year, start_year, Year)) %>% ## combine data from start year and before together %>%
  group_by(genus_sp, AOU, year) %>% ## group by year and species 
  ## calculate range edge metrics for each year for each species:
  summarise(p95 = quantile(Latitude, 0.95), 
            p5 = quantile(Latitude, 0.05),
            max_lat = max(Latitude), 
            min_lat = min(Latitude)) %>% ## calculate range edge metrics for that year       
  ungroup() %>% ## ungroup
  st_drop_geometry() ## drop geometry

## NOTE: is it better to calculate this at a study level? i.e., to include routes sampled that species was not ever seen in? technically this is absence data 

## plot over time
range_shift_sum %>%
  gather(key = "range_edge", value = "latitude", c(p95_w, p95, p5_w, p5)) %>%
  ggplot(aes(x = year, y = latitude, colour = range_edge, group = AOU)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", se = F, size = 0.2) +
  facet_wrap(~range_edge) 

sampling_sum %>%
  gather(key = "range_edge", value = "latitude", c(p95, p5)) %>%
  ggplot(aes(x = year, y = latitude, colour = range_edge, group = AOU)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", se = F, size = 0.2) +
  facet_wrap(~range_edge) 

## plot in species / sample pairs
sampling_sum$type = "edge of sampling"
range_shift_sum$type = "edge of range"

range_shift_sum %>%
  left_join(., select(sampling_sum, max_lat, min_lat, AOU, year)) %>%
  gather(key = "range_edge", value = "latitude", c(p95, p5)) %>%
  filter(genus_sp %in% unique(genus_sp)[26:50]) %>%
  mutate(group = paste(type, range_edge)) %>%
  ggplot(aes(x = year, y = latitude, colour = range_edge, group = group)) +
  geom_point(size = 0.4) +
  geom_smooth(method = "lm", se = F, size = 0.5) +
  facet_wrap(~genus_sp) +
  geom_point(aes(x = year, y = max_lat), colour = "black", size = 0.4) +
  geom_point(aes(x = year, y = min_lat), colour = "black", size = 0.4)

## calculate shift rates by fitting linear regressions 
lm_p95 <- range_shift_sum %>%
  filter(!is.na(p95)) %>%
  group_by(genus_sp) %>%
  filter(n() > 1) %>%
  do(broom::tidy(lm(p95 ~ year, data = .), conf.int = TRUE)) %>% 
  ungroup()

lm_p5 <- range_shift_sum %>%
  filter(!is.na(p5)) %>%
  group_by(genus_sp) %>%
  filter(n() > 1) %>%
  do(broom::tidy(lm(p5 ~ year, data = .), conf.int = TRUE)) %>% 
  ungroup()

lm_p95$range_edge = "leading"
lm_p5$range_edge = "trailing"

lms = rbind(lm_p95, lm_p5)

lms <- mutate(lms, term = ifelse(term != "(Intercept)", "slope", "intercept")) %>%
  select(term, estimate, genus_sp, range_edge) %>%
  spread(key = "term", value = "estimate") 

## plot histogram of shifts 
## units are degrees latitude per y 
lms %>%
  ggplot(aes(x = slope, fill = range_edge)) +
  geom_histogram() +
  facet_wrap(~range_edge)
## leading edge shifts more pronounced, likely because many of the trailing edge shifts are not real trailing edges 
## do higher latitude range edges have greater shifts? if there is a relationship, it is maybe a sign that patterns are maybe confounded by sampling 
## plot in relation to predicted latitude when year was 1980
lms$lat_1980 = lms$slope*1980 + lms$intercept

lms %>%
  ggplot(aes(x = lat_1980, y = slope, fill = range_edge)) +
  geom_point() +
  facet_wrap(~range_edge)

## look at really fast ones
lms$genus_sp[which(lms$slope >= 0.2)]

## Empidonax flaviventris - shift measurement starts at 1993 now but still fast
## Myadestes townsendi - same thing

## write
write.csv(lms, "outputs/data-processed/bbs-range-shift-rates.csv", row.names = F)


############################################
##          map range edge shift          ##
############################################

sp = 1
while(sp <= length(unique(bb$AOU))) {
  
  ## filter to data on just one species 
  sp_df = bb[which(bb$AOU == unique(bb$AOU)[sp]),] 
  
  ## make data frame an sf data frame
  sp_df <- vect(sp_df, geom = c("Longitude", "Latitude")) %>%
    st_as_sf(.)
  
  st_crs(sp_df) = st_crs(countries)
  
  ## add coordinates back to df as columns
  sp_df <- cbind(sp_df, st_coordinates(sp_df)) %>%
    rename("Latitude" = Y, "Longitude" = X)
  
  ## make df of just presence
  sp_pres <- sp_df %>%
    filter(!is.na(TotalAbd) & TotalAbd != 0) %>%
    ## get rid of routes where the species was never seen
    group_by(route) %>%
    mutate(sum_abd = sum(TotalAbd)) %>%
    ungroup() %>%
    filter(sum_abd >= 1) 
  
  ## create directory for plots
  if(!dir.exists(paste0("outputs/figures/range-edge-shifts/", 
                        str_replace_all(unique(sp_df$genus_sp), " ", "\\_")))) {
    dir.create(paste0("outputs/figures/range-edge-shifts/", 
                      str_replace_all(unique(sp_df$genus_sp), " ", "\\_")))
  }
  
  ## figure out which year to start measuring shift
  ## if species' mean max latitude occurs more than than 1 degrees from the mean max latitude if sampling within the longitudinal bounds of the species' presence in 1980, group obs before 1980 and start ts there
  ## otherwise, group obs before 19993 and start ts there
  n_edges_sub = n_edges[which(n_edges$genus_sp == unique(sp_pres$genus_sp)),]
  
  if(min(sp_pres$Year) <= 1980) {
    start_year = as.numeric(ifelse(n_edges_sub$n_edge + 1 <= n_edges_sub$n_edge_sample, 1980, 1993))
  }
  else {
    start_year = 1993
  }
  
  ## plot data in each year and save 
  for(y in 1:length(unique(sp_df$Year))) {
    
    year = unique(sp_df$Year)[order(unique(sp_df$Year))][y]
    
    ## if data for that year
    if(length(which(sp_df$Year == year & sp_df$TotalAbd != 0 & !is.na(sp_df$TotalAbd))) != 0) {
      
      if(year > start_year) {
        ## calculate 95th percentile of occupied cells
        cur <- sp_pres %>%
          filter(Year == year) %>%
          mutate(p95_w = modi::weighted.quantile(Latitude, w = TotalAbd, 0.95),
                 p5_w = modi::weighted.quantile(Latitude, w = TotalAbd, 0.05)) %>%
          mutate(p95 = quantile(Latitude, 0.95), 
                 p5 = quantile(Latitude, 0.05)) 
        
        p = cur %>%
          ggplot(aes(colour = TotalAbd)) +
          geom_sf(inherit.aes = FALSE, data = countries) +
          geom_sf(size = 0.5, data = filter(sp_df, TotalAbd == 0 & Year == year), inherit.aes = FALSE, 
                  fill = "transparent", stroke = 0.1, shape = 1, alpha = 0.5) +
          geom_sf(size = 0.1) +
          scale_colour_gradient(trans = "log", limits = c(1, max(sp_df$TotalAbd, na.rm = T))) +
          scale_x_continuous(limits = c(-165,-55)) +
          scale_y_continuous(limits = c(25,70)) +
          labs(colour = "Total abundance", title = year) +
          geom_hline(aes(yintercept = unique(p95)), colour = "red") +
          geom_hline(aes(yintercept = unique(p95_w)), colour = "blue") +
          geom_hline(aes(yintercept = unique(p5)), colour = "red") +
          geom_hline(aes(yintercept = unique(p5_w)), colour = "blue")
        
        ggsave(p, path = paste0("outputs/figures/range-edge-shifts/", 
                                str_replace_all(unique(sp_df$genus_sp), " ", "\\_")), 
               filename = paste0("BBS_", str_replace_all(unique(sp_df$genus_sp)," ", "\\_"), "_", 
                                 year, ".png"), width = 4, height = 3)
      }
      else if (year == start_year){
        min_year = min(sp_df$Year)
        
        ## calculate 95th percentile of occupied cells
        cur <- sp_pres %>%
          filter(Year <= start_year) %>%
          mutate(p95_w = modi::weighted.quantile(Latitude, w = TotalAbd, 0.95),
                 p5_w = modi::weighted.quantile(Latitude, w = TotalAbd, 0.05)) %>%
          mutate(p95 = quantile(Latitude, 0.95), 
                 p5 = quantile(Latitude, 0.05)) 
        
        absences = filter(sp_df, TotalAbd == 0 & Year <= year) %>%
          select(geometry) %>%
          unique()
        
        p = cur %>%
          ggplot(aes(colour = TotalAbd)) +
          geom_sf(inherit.aes = FALSE, data = countries) +
          geom_sf(size = 0.5, data = absences, inherit.aes = FALSE,
                  fill = "transparent", stroke = 0.1, shape = 1, alpha = 0.5) +
          geom_sf(size = 0.1) +
          scale_colour_gradient(trans = "log", limits = c(1, max(sp_df$TotalAbd, na.rm = T))) +
          scale_x_continuous(limits = c(-165,-55)) +
          scale_y_continuous(limits = c(25,70)) +
          labs(colour = "Total abundance", title = year) +
          geom_hline(aes(yintercept = unique(p95)), colour = "red") +
          geom_hline(aes(yintercept = unique(p95_w)), colour = "blue") +
          geom_hline(aes(yintercept = unique(p5)), colour = "red") +
          geom_hline(aes(yintercept = unique(p5_w)), colour = "blue")
        
        ggsave(p, path = paste0("outputs/figures/range-edge-shifts/", 
                                str_replace_all(unique(sp_df$genus_sp), " ", "\\_")), 
               filename = paste0("BBS_", 
                                 str_replace_all(unique(sp_df$genus_sp)," ", "\\_"),"_", 
                                 min_year, "to", start_year, ".png"), width = 4, height = 3)
      }
    }
  }
  
  print(paste0("Species: ", sp))
  sp = sp + 1
}
