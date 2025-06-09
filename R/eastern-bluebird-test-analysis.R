## test analysis of whether autocorrelation + synchrony explain patterns of extinction/establishment in the bbs data
## using eastern bluebird as an example
library(tidyverse)
library(tidyterra)
library(sf)
library(terra)
theme_set(theme_bw())

###########################################
##                bbs data               ##
###########################################
## read in bbs data for all spp
cur_all <- read.csv("outputs/data-processed/BBS_all.csv")

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

## filter to data on a only eastern bluebird
ebb_all <- filter(cur_all, AOU == 7660) 

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

## turn bluebird data into vect
ebb_all = vect(ebb_all, geom = c("Longitude", "Latitude"))

## bring in North America outline 
# get map data
crs(ebb_all) <- as.character(crs(countries))

ebb_all = st_as_sf(ebb_all)

## plot data in each year and save 
for(y in 1:length(unique(ebb_all$Year))) {
  year = unique(ebb_all$Year)[order(unique(ebb_all$Year))][y]
  p = ebb_all %>%
    filter(!TotalAbd == 0 & Year == year) %>%
    ggplot(aes(colour = TotalAbd)) +
    geom_sf(inherit.aes = FALSE, data = countries) +
    geom_sf(size = 1, data = filter(ebb_all, TotalAbd == 0 & Year == year), inherit.aes = FALSE, 
            fill = "transparent", linewidth = 0.0001, shape = 1) +
    geom_sf(size = 1) +
    scale_colour_gradient(trans = "log", limits = c(1, 191)) +
    scale_x_continuous(limits = c(-130,-60)) +
    scale_y_continuous(limits = c(25,60)) +
    labs(colour = "Total abundance", title = year)
  
  ggsave(p, path = "outputs/figures/eastern-bluebird", 
         filename = paste0("BBS_eastern-bluebird_", year, ".png"), width = 4, height = 3)
}

## for each route, add in rows with NA for years that were not sampled within the sampling period
ebb_all <- ebb_all %>%
  select(-code) %>%
  mutate(route = paste(Route, CountryNum, StateNum, sep = "_")) 

ebb_key <- ebb_all %>% 
  group_by(route) %>%
  mutate(first_year = min(Year), last_year = max(Year)) %>%
  select(route, first_year, last_year) %>%
  st_drop_geometry() %>%
  distinct()

for(i in 1:nrow(ebb_key)) {
  df <- data.frame(Year = seq(ebb_key$first_year[i], ebb_key$last_year[i], by = 1), route = ebb_key$route[i])
  
  if(i==1) {
    df_all = df
  }
  else {
    df_all = rbind(df_all, df)
  }
}

t = left_join(df_all, ebb_all) %>%
  group_by(route) %>%
  tidyr::fill(colnames(.)[which(colnames(.) != "TotalAbd")], .direction = "updown") %>%
  as.data.frame() %>%
  select(-geometry)

geom = ebb_all %>%
  select(route) %>%
  unique()

## add geometry 
t <- left_join(t, geom)

## now, count extinctions, establishments and whether cell is new establishment or not
## split by routes
split = t %>%
  group_by(route) %>%
  group_split(.)

list = lapply(split, FUN = function(df) {
  ts = df$TotalAbd
  zero_count = 0
  ext = 0
  est = 0
  max_zero_count <- 0
  for(i in 1:length(ts)) {
    if(i == 1) {
      if(ts[i] == 0) {
        zero_count = 1
      }
    }
    else if(is.na(ts[i]))  {
      if(zero_count > max_zero_count) {
        max_zero_count = zero_count
      }
      zero_count = 0
      if(i == length(ts)) {
        if(max_zero_count >= 2) {
          ext = ext + 1
        }
      }
    }
    else if(ts[i] == 0) {
      zero_count = zero_count + 1
      if(zero_count > max_zero_count) {
        max_zero_count = zero_count
      }
      if(i == length(ts)) {
        if(max_zero_count >= 2) {
          ext = ext + 1
        }
      }
    }
    else if(ts[i] >= 1 & max_zero_count >= 2) {
      if(!any(ts[which(!is.na(ts[1:(i-1)]))] != 0) & first(ts) == 0) {
        ext = 0
      }
      else {
        ext = ext + 1
      }
      est = est + 1
      zero_count = 0
      max_zero_count = 0
    }
    else if(ts[i] >= 1 & max_zero_count <= 2) {
      zero_count = 0
      max_zero_count = 0
    }
  }
  
  if(any(ts[which(!is.na(ts))[1:2]] >= 1)) {
    new_est = FALSE
  }
  else {
    new_est = TRUE
  }
  
  ext_est = df %>%
    mutate(n_est = est, n_ext = ext, new_est = new_est, n_ts = length(ts))
  
  return(ext_est)
})

## bind them all 
list <- bind_rows(list)

## plot some time series 
temp = list %>%
  arrange(Route, CountryNum, StateNum, Year) %>%  
  group_by(Route, CountryNum, StateNum) %>%
  mutate(deltaLag = Year - lag(Year, 1)) %>% 
  group_by(Route, CountryNum, StateNum, Year) %>%
  mutate(sequence1 = case_when(is.na(deltaLag) | deltaLag > 1 ~ 1,
                               TRUE ~ 2)) %>% 
  ungroup() %>% 
  mutate(sequence = cumsum(sequence1==1)) %>% 
  select(-deltaLag, -sequence1) %>% 
  ungroup() 
temp %>%
  filter(route == unique(.$route)[2]) %>%
  ggplot(aes(x = Year, y = TotalAbd)) +
  geom_point() +
  geom_line(aes(group = sequence)) +
  geom_vline(aes(xintercept = min(Year)), colour = "red") +
  geom_vline(aes(xintercept = max(Year)), colour = "red")

temp %>%
  filter(route == unique(.$route)[2]) %>% View

## plot histogram 
list %>%
  select(route, n_ext, n_est) %>%
  distinct() %>% 
  gather(key = "measure", value = "count", c(n_ext, n_est)) %>%
  mutate(measure = ifelse(measure == "n_ext", "Extinctions", "Establishments")) %>%
  ggplot(aes(x = count, fill = measure)) +
  geom_histogram() +
  facet_wrap(~measure) + 
  theme(legend.position = "none") +
  labs(x = "No. of local events", y = "Frequency")

list %>%
  ggplot(aes(x = new_est)) +
  geom_bar() +
  theme(legend.position = "none") +
  labs(x = "No. of routes", y = "New estabishment?")

list = vect(st_as_sf(list))

## plot no. extinctions and no. establishments against autocorrelation at site
###########################################
##          autocorrelation data         ##
###########################################
## read in analysis of Berkeley Earth data
se <- readRDS("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/BerkeleyEarth_noise-colour.rds")

## select data on low spectral exponent measured across entire time period
se =  se %>%
  select(lat, lon, s_spec_exp_PSD_low_all) %>%
  unique()

se %>%
  ggplot(aes(x = lon, y = lat, fill = s_spec_exp_PSD_low_all)) +
  geom_raster() +
  coord_fixed() +
  scale_fill_gradient(high = "#FF5C7A", low = "#FFFFFF", na.value = "grey") +
  labs(fill = "Spectral exponent") +
  theme(panel.background = element_rect(fill = "#ECECEC", colour = "#ECECEC"))

## extract temporal autocorrelation at each point the northern cardinal was sampled
se$lon = se$lon - 180
r = se %>%
  select(lon, lat, s_spec_exp_PSD_low_all) %>%
  rast()
crs(r) <- crs(ebb_all)

plot(r)

## extract spectral exponent values across each route where eastern bluebird was sampled
beta <- extract(r, list)
list$beta = beta$s_spec_exp_PSD_low_all

list %>% 
  ggplot(aes(colour = beta)) +
  geom_sf(inherit.aes = FALSE, data = countries) +
  geom_sf(size = 0.5) +
  scale_x_continuous(limits = c(-130,-60)) +
  scale_y_continuous(limits = c(25,60)) +
  scale_colour_viridis_c() +
  labs(colour = "")

df = as.data.frame(list) %>%
  group_by(route) %>%
  mutate(n_ts = length(unique(Year))) %>% ## calculate time series length 
  ungroup() %>% 
  select(route, beta, n_ext, n_est, new_est, n_ts) %>%
  distinct()

df %>%
  gather(key = "measure", value = "count", c(n_ext, n_est)) %>%
  mutate(measure = ifelse(measure == "n_ext", "Extinctions", "Establishments")) %>%
  ggplot(aes(x = beta, y = count)) +
  geom_point(position = position_jitter(width = 0.01, height = 0.25)) +
  geom_smooth(method = "lm") +
  facet_wrap(~measure) + 
  labs(x = "Spectral exponent", y = "Frequency")
## extinction and establishment increase with increase autocorrelation :-)
df %>%
  mutate(new_est = ifelse(new_est == FALSE, 0, 1)) %>%
  ggplot(aes(x = beta, y = new_est)) +
  geom_point(position = position_jitter(width = 0.01, height = 0.25)) +
  geom_smooth(method = "lm") + 
  labs(x = "Spectral exponent", y = "New establishment?")
  
## so does probability of new establishment in a cell ?

## note for later:
# longer time series = more chance of seeing extinction/establishment? control for this 
df %>%
  gather(key = "measure", value = "count", c(n_ext, n_est)) %>%
  mutate(measure = ifelse(measure == "n_ext", "Extinctions", "Establishments")) %>%
  ggplot(aes(x = n_ts, y = count)) +
  geom_point(position = position_jitter(width = 0.01, height = 0.25)) +
  geom_smooth(method = "lm") + 
  labs(y = "No. of local events", x = "Length of time series (yrs)") +
  facet_wrap(~measure)

## next: separate into range centre vs. range edge populations
## see if there is:
## - more est and ext at the range edges 
## - clearer effect of autocorrelation at the range edges  

## read in birds of the world range 
#botw = st_read("/Users/nikkimoore/Documents/bioshifts-traits/data-raw/large-data/BirdsOfTheWorld/BOTW/BOTW.gdb")

## filter to eastern bluebird
ebb_range = botw[which(botw$sci_name == "Sialia sialis"),]

## save 
st_write(ebb_range, "outputs/data-processed/ebb_range.gdb", overwrite = T)

ebb_range = st_read("outputs/data-processed/ebb_range.gdb")

# 1 - resident
# 2 - breeding
# 3 - non-breeding
# 4 - passage
# 5 - uncertain

ebb_range = ebb_range %>%
  mutate(range_type = ifelse(seasonal == 1, "Resident", 
                             ifelse(seasonal == 2, "Breeding",
                                    ifelse(seasonal == 3, "Non-breeding",
                                           ifelse(seasonal == 4, "Migration",
                                                  ifelse(seasonal == 5, "Rare/uncertain", NA))))))

## plot all together 
ebb_range %>%
  ggplot(aes(fill = range_type)) +
  geom_sf(inherit.aes = FALSE, data = countries, fill = "transparent") +
  geom_sf() +
  scale_x_continuous(limits = c(-130,-60)) +
  theme(panel.grid = element_blank())
  
## plot routes where ebb was seen on top of just breeding + resident range
list %>%
  filter(!is.na(TotalAbd), TotalAbd != 0) %>%
  ggplot(aes()) +
  geom_sf(inherit.aes = FALSE, data = countries) +
  geom_sf(data = ebb_range, inherit.aes = FALSE, aes(fill = range_type)) +
  # geom_sf(size = 1, data = filter(list, TotalAbd == 0 & Year == year), inherit.aes = FALSE, 
  #         fill = "transparent", linewidth = 0.0001, shape = 1) +
  geom_sf(size = 0.5) +
  scale_colour_gradient(trans = "log") +
  scale_x_continuous(limits = c(-130,-60)) +
  scale_y_continuous(limits = c(25,60)) 

## combine all range bits
union = st_union(ebb_range[which(ebb_range$range_type %in% c("Breeding", "Resident")),])
plot(union)

## buffer edge 
buff = st_buffer(union, 2)
plot(buff)

negbuff = st_buffer(union, -2)
plot(negbuff)

edge = mask(list, vect(negbuff), inverse = T)
plot(edge)

centre =  mask(list, vect(negbuff))
plot(centre)

centre$range_position = "centre"
edge$range_position = "edge"

range = rbind(centre, edge)

## plot range edge vs. centre points 
range %>%
  ggplot(aes(colour = range_position)) +
  geom_sf(inherit.aes = FALSE, data = countries) +
  geom_sf(size = 0.5) +
  scale_x_continuous(limits = c(-130,-60)) +
  scale_y_continuous(limits = c(25,60))

## plot ext. and est. on map - see if higher at range edge 
range %>%
  select(n_ext, route) %>%
  distinct() %>%
  ggplot(aes(colour = n_ext)) +
  geom_sf(inherit.aes = FALSE, data = countries) +
  geom_sf(size = 0.5) +
  scale_x_continuous(limits = c(-130,-60)) +
  scale_y_continuous(limits = c(25,60)) +
  labs(colour = "No. extinctions")
  
range %>%
  select(n_est, route) %>%
  distinct() %>%
  ggplot(aes(colour = n_est)) +
  geom_sf(inherit.aes = FALSE, data = countries) +
  geom_sf(size = 0.5) +
  scale_x_continuous(limits = c(-130,-60)) +
  scale_y_continuous(limits = c(25,60)) +
  labs(colour = "No. establishments")

range %>%
  select(new_est, route) %>%
  distinct() %>%
  ggplot(aes(colour = new_est)) +
  geom_sf(inherit.aes = FALSE, data = countries) +
  geom_sf(size = 0.5) +
  scale_x_continuous(limits = c(-130,-60)) +
  scale_y_continuous(limits = c(25,60)) +
  labs(colour = "New establishment?")


range = as.data.frame(range) %>%
  group_by(route) %>%
  mutate(n_ts = length(unique(Year))) %>% ## calculate time series length 
  ungroup() %>% 
  select(route, beta, n_ext, n_est, new_est, n_ts, range_position) %>%
  distinct()

## learn some thing about edge vs. centre routes 
range %>%
  gather(key = "measure", value = "count", c(n_ext, n_est)) %>%
  mutate(measure = ifelse(measure == "n_ext", "Extinctions", "Establishments")) %>%
  select(route, measure, count, range_position) %>%
  distinct() %>%
  ggplot(aes(x = count, fill = range_position)) +
  geom_histogram(position = "dodge") +
  facet_wrap(~measure) + 
  labs(x = "No. of local events", y = "Frequency")
## edge routes have greater frequency of local extinction and establishment
range %>%
  ggplot(aes(x = new_est, fill = range_position)) +
  geom_bar(position = "dodge") +
  theme(legend.position = "none") +
  labs(x = "No. of routes", y = "New estabishment?")

range %>%
  gather(key = "measure", value = "count", c(n_ext, n_est)) %>%
  mutate(measure = ifelse(measure == "n_ext", "Extinctions", "Establishments")) %>%
  ggplot(aes(x = beta, y = count, colour = range_position)) +
  geom_point(position = position_jitter(width = 0.01, height = 0.25),
             size = 0.1) +
  geom_smooth(method = "lm") +
  facet_grid(measure~range_position) +
  labs(x = "Spectral exponent", y = "No. of local events")

## extinction and establishment increase with increasing autocorrelation across edge and centre :-)
range %>%
  mutate(new_est = ifelse(new_est == TRUE, 1, 0)) %>%
  ggplot(aes(x = beta, y = new_est, colour = range_position)) +
  geom_point() +
  geom_point(position = position_jitter(width = 0.01, height = 0.25),
                           size = 0.1) +
  geom_smooth(method = "lm") +
  facet_wrap(~range_position) +
  labs(x = "Spectral exponent", y = "New establishment?")
## so does probability of new establishment in a cell 
range %>%
  ggplot(aes(x = new_est, fill = range_position)) +
  geom_histogram(position = "dodge")
## greater number of new establishments at the range edge
range %>%
  ggplot(aes(x = new_est, fill = range_position)) +
  geom_histogram(position = "dodge")

c = length(which(range$range_position == "centre"))
e = length(which(range$range_position == "edge"))

c_new_est = length(which(range$range_position == "centre" & range$new_est == 1))/length(which(range$range_position == "centre"))
e_new_est = length(which(range$range_position == "edge" & range$new_est == 1))/length(which(range$range_position == "edge"))
c_new_est ## 20% of cells experience new establishment in range centre 
e_new_est ## 53% of cells experience new establishment at range edge

c_no_ext = length(which(range$range_position == "centre" & range$n_ext == 0))/length(which(range$range_position == "centre"))
e_no_ext = length(which(range$range_position == "edge" & range$n_ext == 0))/length(which(range$range_position == "edge"))
c_no_ext ## 65% of cells experience no extinction in range centre 
e_no_ext ## 29% of cells experience no extinction at range edge

##### try measuring range shift 
## measure shift in abundance weighted mean across latitude

## add coordinates to dataframe
ebb_df <- cbind(ebb_all, st_coordinates(ebb_all))

## see if sampling routes have shifted latitude and longitude over time 
ebb_df %>%
  select(route, Year, X, Y) %>%
  unique() %>%
  select(-route) %>%
  group_by(Year) %>% 
  summarize(mean_lat = mean(Y), mean_lon = mean(X)) %>%
  ggplot(aes(x = Year, y = mean_lat)) +
  geom_point() +
  geom_line() +
  geom_smooth(method = 'lm')

ebb_df %>%
  select(route, Year, X, Y) %>%
  unique() %>%
  select(-route) %>%
  group_by(Year) %>% 
  summarize(mean_lat = mean(Y), mean_lon = mean(X)) %>%
  ggplot(aes(x = Year, y = mean_lon)) +
  geom_point() +
  geom_line() +
  geom_smooth(method = 'lm')

## for each year, calculate: (latitude x abundance) / total abundance
ebb_df %>%
  select(route, TotalAbd, Year, X, Y) %>%
  unique() %>%
  select(-route) %>%
  group_by(Year) %>% 
  mutate(TotalAbd_year = sum(TotalAbd)) %>%
  mutate(weighted_lat = sum(Y*TotalAbd)/TotalAbd_year)  %>%
  ungroup() %>%
  select(Year, weighted_lat) %>%
  unique() %>%
  ggplot(aes(x = Year, y = weighted_lat)) +
  geom_point() +
  geom_smooth(method = "lm")

ebb_df %>%
  select(route, TotalAbd, Year, X, Y) %>%
  unique() %>%
  group_by(Year) %>% 
  mutate(TotalAbd_year = sum(TotalAbd)) %>%
  select(-route) %>%
  mutate(weighted_lat = sum(Y*TotalAbd)/TotalAbd_year) %>%
  ungroup() %>%
  select(Year, weighted_lat) %>%
  unique() %>%
  ggplot(aes(x = Year, y = weighted_lat)) +
  geom_point() +
  geom_smooth(method = "lm") 


ebb_df %>%
  select(route, TotalAbd, Year, X, Y) %>%
  unique() %>%
  select(-route) %>%
  filter(TotalAbd != 0) %>%
  group_by(Year) %>% 
  mutate(p95 = quantile(Y, c(0.95)),
         p5 = quantile(Y, c(0.05))) %>%
  mutate(b4 = ifelse(Year <= 1980, "pre-warming", "post-warming")) %>%
  ungroup() %>%
  select(Year, p95, p5, b4) %>%
  unique() %>%
  ggplot(aes(x = Year, y = p95)) +
  geom_point() +
  geom_smooth(method = "lm", aes(group = b4)) + 
  geom_point(aes(y = p5), colour = "red") 

## figure out how to:
## 1. measure rate of climate change (e.g., expected shift)
## 2. define range edge 
## 3. test for effect of synchrony / autocorrelation
## -> within species? (e.g., multiple measures of edge shift across species range)
## -> across species? (e.g., one measure of edge shift, synchrony, autococorrelation per species per edge)
## 4. what to do when trailing edge isn't actual trailing edge (bound by continent)
## 5. do we want to measure synchrony in population dynamics? or in the environment? we are assuming that synchrony in the environment leads to synchrony in pop dynamic, but why not just directly measure?
## 6. what about using population decline/increase at range edges as metric?
## shift in optima instead of edges?  
## if species wasn't present before 1980, consider it a new establishment. if it was and disappears, consider it an extinction
## WHAT ABOUT MEASURING SHIFTS/LAGS IN NICHE SPACE? 


## try measuring occupancy rate across latitude 
ebb_all = cbind(ebb_all, st_coordinates(ebb_all))

## bin into latitudinal bins
ebb_all$lat_bin = round(ebb_all$Y, digits = 0)

lat_occ = ebb_all %>%
  arrange(lat_bin, Year) %>%
  group_by(lat_bin, Year) %>%
  summarize(Nt_total = sum(TotalAbd),
            occ =  sum(TotalAbd >= 1),
            n_obs = n()) %>%
  mutate(occ_rate = occ/n_obs)

lat_occ %>%
  filter(Year %in% c(1980, 1990, 2001, 2015, 2022)) %>%
  ggplot(aes(x = lat_bin, y = occ_rate, colour = Year, group = Year)) +
  geom_point(size = 0.5) +
  geom_smooth()

## plot occupancy at one latitude over time
lat_occ %>%
  filter(lat_bin == 52) %>%
  ggplot(aes(x = Year, y = occ_rate)) +
  geom_point(size = 0.5) +
  geom_smooth()

## try plotting metrics by latitudinal bin? 
