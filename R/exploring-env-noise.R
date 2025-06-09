## exploring temporal autocorrelation
library(sf)
library(tidyverse)
theme_set(theme_bw())

## read in analysis of Berkeley Earth data
data <- readRDS("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/BerkeleyEarth_noise-colour.rds")

## plot low spectral exponent measured across entire time period on a map
sub = data %>%
  select(lat, lon, s_spec_exp_PSD_low_all) %>%
  unique()

sub %>%
  ggplot(aes(x = lon, y = lat, fill = s_spec_exp_PSD_low_all)) +
  geom_raster() +
  coord_fixed() +
  scale_fill_gradient(high = "#FF5C7A", low = "#FFFFFF", na.value = "grey", limits = c(-0.25, 3)) +
  labs(fill = "Spectral exponent") +
  theme(panel.background = element_rect(fill = "#ECECEC", colour = "#ECECEC"))

## see how much variation in spectral exponent there is across the Northern Cardinal sampling sites in the Breeding Bird Survey
bbs <- read.csv("outputs/data-processed/BBS_all.csv")

## read in weather data
weather = read.csv("data-raw/BBS/Weather.csv")
## filter out routes with RunType == 0
weather = filter(weather, RunType != 0)
## make code for route
weather$code = paste(weather$Route, weather$CountryNum, weather$StateNum, weather$Year, sep = "_")

## filter to sites where eastern blue bird was seen
nc <- filter(bbs, AOU == 7660)

## add absence data 
# group by route, add years that were sampled on that route but where species wasn't found
sampling_scheme <- weather %>% 
  select(Route, CountryNum, StateNum, Year, code, RouteDataID, RPID)

nc = nc %>%
  left_join(sampling_scheme, .) %>%
  mutate(TotalAbd = ifelse(is.na(TotalAbd), 0, TotalAbd)) %>%
  mutate(temp = paste(Route, CountryNum, StateNum)) %>%
  group_by(temp) %>%
  tidyr::fill(., RPID, AOU, Latitude, Longitude, RouteName, .direction = "updown") %>%
  ungroup() %>%
  filter(!is.na(AOU)) %>% ## get rid of rows representing 0 abd in plots where spp was never seen once
  select(-temp)

## plot across space 
nc = vect(nc, geom = c("Longitude", "Latitude"))

## bring in North America outline 
# get map data
countries <- necountries::countries(c("Canada", "United States of America"), part = TRUE)
crs(nc) <- as.character(crs(countries))

nc = st_as_sf(nc)

min(nc$Year)
year = 1975
nc %>%
  filter(!TotalAbd == 0 & Year == year) %>%
  ggplot(aes(colour = TotalAbd)) +
  geom_raster(inherit.aes = FALSE, data = sub, 
              aes(x = lon - 180, y = lat, fill = s_spec_exp_PSD_low_all)) +
  geom_sf(inherit.aes = FALSE, data = countries, fill = "transparent") +
  geom_sf(size = 1, data = filter(nc, TotalAbd == 0 & Year == year), inherit.aes = FALSE, 
          fill = "transparent", linewidth = 0.0001, shape = 1) +
  geom_sf(size = 1) +
  scale_colour_gradient(trans = "log") +
  scale_x_continuous(limits = c(-130,-60)) +
  scale_y_continuous(limits = c(25,60)) +
  scale_fill_gradient(high = "#FF5C7A", low = "#FFFFFF", na.value = "grey", limits = c(-0.25, 1.5)) +
  labs(fill = "Spectral exponent") 


## extract temporal autocorrelation at each point the northern cardinal was sampled
sub$lon = sub$lon - 180
r = sub %>%
  select(lon, lat, s_spec_exp_PSD_low_all) %>%
  rast()
crs(r) <- crs(nc)

plot(r)

beta <- extract(r, nc)
nc$beta = beta$s_spec_exp_PSD_low_all

## histogram of beta 
nc %>%
  ggplot(aes(x = beta)) +
  geom_histogram()
min(nc$beta, na.rm = TRUE) ## ranges from 0.56
max(nc$beta, na.rm = TRUE) ## to 1.22
max(nc$beta, na.rm = TRUE) - min(nc$beta, na.rm = TRUE) ## 0.66 range 

## let's see within the study extent of the bbs 
poly <- vect("outputs/data-processed/bbs_study-extent.shp")

beta_bbs <- mask(r, poly)
plot(beta_bbs)

min(values(beta_bbs), na.rm = T) ## ranges from 0.3
max(values(beta_bbs), na.rm = T) ## to 1.4
max(values(beta_bbs), na.rm = TRUE) - min(values(beta_bbs), na.rm = TRUE) ## 1.1 range 


