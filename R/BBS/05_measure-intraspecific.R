## try measuring range edge shift, synchrony, and autocorrelation at multiple points along Cardinal's range edge
library(tidyverse)
library(terra)
library(raster)
library(sf)
theme_set(theme_bw())

## read in clean bbs data
bb =  read.csv("outputs/data-processed/bbs_clean-subset.csv")

## get canada and us map for plotting
countries <- necountries::countries(c("Canada", "United States of America"), part = TRUE)
countries <- st_union(countries)
# countries <- vect(countries)

## filter to Cardinalis cardinalis
cc = bb[which(bb$genus_sp == "Cardinalis cardinalis"),] 

## make data frame an sf data frame
cc <- vect(cc, geom = c("Longitude", "Latitude")) %>%
  st_as_sf(.)

st_crs(cc) = st_crs(countries)



## start with a simple method 
## divide extent of study into 10 deg longitude blocks 
## measure shift separately within each chunk
min(bb$Longitude)
max(bb$Longitude)

ext = ext(-180, -50, 18, 84)

xmins = seq(from = -180, to = -60, by = 10)
xmaxs = seq(from = -170, to = -50, by = 10)

polys = c()
for(i in 1:length(xmins)) {
  polys = st_as_sf(vect(ext(xmins[i], xmaxs[i], 18, 84))) %>%
    rbind(polys, .)
}

plot(polys)

st_crs(polys) = st_crs(countries)

#######################################
##          range shift          ##
#######################################
## plot 
absences = filter(cc, TotalAbd == 0) %>%
  select(geometry) %>%
  unique()

cc %>%
  filter(TotalAbd != 0) %>%
  ggplot(aes(colour = TotalAbd)) +
  geom_sf(inherit.aes = FALSE, data = countries) +
  geom_sf(size = 0.5, data = absences, inherit.aes = FALSE,
          fill = "transparent", stroke = 0.1, shape = 1, alpha = 0.5) +
  geom_sf(size = 0.1) +
  scale_colour_gradient(trans = "log", limits = c(1, max(cc$TotalAbd, na.rm = T))) +
  scale_x_continuous(limits = c(-165,-55)) +
  scale_y_continuous(limits = c(25,70)) +
  labs(colour = "Total abundance") +
  geom_sf(inherit.aes = F, aes(), polys, fill = "transparent") +
  theme(panel.grid = element_blank())
 
## measure range shift within each section
polys$ID = 1:nrow(polys)

## assign points to polygons 
points_in_poly <- terra::extract(vect(polys), vect(cc))

cc$polyID = points_in_poly$ID

cc %>%
  filter(TotalAbd != 0) %>%
  ggplot(aes(colour = polyID)) +
  geom_sf(inherit.aes = FALSE, data = countries) +
  geom_sf(size = 0.5, data = absences, inherit.aes = FALSE,
          fill = "transparent", stroke = 0.1, shape = 1, alpha = 0.5) +
  geom_sf(size = 0.1) +
  scale_x_continuous(limits = c(-165,-55)) +
  scale_y_continuous(limits = c(25,70)) +
  labs(colour = "Total abundance") +
  geom_sf(inherit.aes = F, aes(), polys, fill = "transparent") +
  theme(panel.grid = element_blank())

## get rid of polygons that have less than 100 non-zero abundances 
cc = cc %>%
  group_by(genus_sp, polyID) %>%
  mutate(count = length(which(TotalAbd != 0))) %>%
  filter(count >= 100) %>%
  select(-count) %>%
  ungroup()
  

## now measure range shift in each 
## group by polygon polyID, then:
## for each year, calculate mean latitude of 5 most northern samples within the longitudinal bounds of each species for data before 1980
range_shift_sum <- cc %>%
  filter(!TotalAbd == 0 & !is.na(TotalAbd)) ## get rid of absence data and NAs (i.e., years where routes weren't sampled)

## add coordinates back to df as columns
range_shift_sum <- cbind(range_shift_sum, st_coordinates(range_shift_sum)) %>%
  rename("Latitude" = Y, "Longitude" = X)

## and max/mins
mean_max_lat = range_shift_sum %>%
  mutate(start_year = 1980) %>% ## decide what year to start at based on northern range edge position in 1980 relative to sampling 
  ## if latitudinal edge of sampling before 1980 is more 10 degree latitude higher than range edge, start measuring shift at 1980
  ## otherwise start at 1993
  mutate(year = ifelse(Year <= start_year, start_year, Year)) %>% ## combine data from start year and before together
  group_by(genus_sp, AOU, year, polyID) %>%
  arrange(-Latitude, by_group = T) %>%
  slice_head(n = 5) %>% 
  summarise(mean_max_lat = mean(Latitude))

mean_min_lat = range_shift_sum %>%
  mutate(start_year = 1980) %>% ## decide what year to start at based on northern range edge position in 1980 relative to sampling 
  ## if latitudinal edge of sampling before 1980 is more 10 degree latitude higher than range edge, start measuring shift at 1980
  ## otherwise start at 1993
  mutate(year = ifelse(Year <= start_year, start_year, Year)) %>% ## combine data from start year and before together
  group_by(genus_sp, AOU, year, polyID) %>%
  arrange(Latitude, by_group = T) %>%
  slice_head(n = 5) %>% 
  summarise(mean_min_lat = mean(Latitude))

range_shift_quants = range_shift_sum %>%
  mutate(start_year = 1980) %>% ## decide what year to start at based on northern range edge position in 1980 relative to sampling 
  ## if latitudinal edge of sampling before 1980 is more 10 degree latitude higher than range edge, start measuring shift at 1980
  ## otherwise start at 1993
  mutate(year = ifelse(Year <= start_year, start_year, Year)) %>% ## combine data from start year and before together
  group_by(genus_sp, AOU, year, polyID) %>% ## group by year and species and polygon
  ## calculate range edge metrics for each year for each species:
  summarise(p95_w = modi::weighted.quantile(Latitude, w = TotalAbd, 0.95),
            p5_w = modi::weighted.quantile(Latitude, w = TotalAbd, 0.05),
            p95 = quantile(Latitude, 0.95), 
            p5 = quantile(Latitude, 0.05)) %>%  
  ungroup() %>% ## ungroup
  st_drop_geometry()  ## drop geometry

range_shift_sum = left_join(st_drop_geometry(mean_min_lat), st_drop_geometry(mean_max_lat)) %>%
  left_join(., st_drop_geometry(range_shift_quants))

## plot over time
range_shift_sum %>%
  gather(key = "range_edge", value = "latitude", c(p95_w, p95, p5_w, p5, mean_max_lat, mean_min_lat)) %>%
  ggplot(aes(x = year, y = latitude, colour = polyID, group = polyID)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", se = F, size = 0.2) +
  facet_wrap(~range_edge) 

## make example plot for cardinal
## plot over time
plot = range_shift_sum %>%
  gather(key = "range_edge", value = "latitude", c(p95_w, p95, p5_w, p5, mean_min_lat, mean_max_lat)) %>%
  mutate(group = paste(range_edge, polyID)) %>%
  ggplot(aes(x = year, y = latitude, colour = range_edge, group = group)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", se = F, size = 0.2) +
  scale_color_manual(values = c("green", "green", "red", "blue", "red", "blue")) +
  labs(x = "Year", y = "Latitude") 

ggsave(plot, path = "outputs/figures/range-edge-shifts/Cardinalis_cardinalis/", 
       filename = "plot-year-vs-latitude_intraspecific.png", width = 5, height = 3)

## make info about polys
polys$xmins = xmins
polys$xmaxs = xmaxs
polys_df = polys %>% 
  st_drop_geometry()

range_shift_sum <- left_join(range_shift_sum, polys_df, by = c("polyID" = 'ID'))

## plot 
for(y in 1980:2023) {
  
  ## plot 
  absences = filter(cc, TotalAbd == 0, Year == y) %>%
    select(geometry) %>%
    unique()
  
  p = cc %>%
    filter(TotalAbd != 0, Year == y) %>%
    ggplot(aes(colour = TotalAbd)) +
    geom_sf(inherit.aes = FALSE, data = countries) +
    geom_sf(size = 0.5, data = absences, inherit.aes = FALSE,
            fill = "transparent", stroke = 0.1, shape = 1, alpha = 0.5) +
    geom_sf(size = 0.1) +
    scale_x_continuous(limits = c(-165,-55)) +
    scale_y_continuous(limits = c(25,70)) +
    labs(colour = "Total abundance") +
    geom_sf(inherit.aes = F, aes(), polys, fill = "transparent") +
    theme(panel.grid = element_blank()) +
    scale_colour_gradient(trans = "log", limits = c(1, max(cc$TotalAbd, na.rm = T))) +
    geom_segment(data = filter(range_shift_sum, year == y), inherit.aes = F,
                 aes(x = xmins, xend = xmaxs, y = p95), colour = "red") +
    geom_segment(data = filter(range_shift_sum, year == y), inherit.aes = F,
                 aes(x = xmins, xend = xmaxs, y = p5), colour = "red") +
    geom_segment(data = filter(range_shift_sum, year == y), inherit.aes = F,
                 aes(x = xmins, xend = xmaxs, y = p95_w), colour = "blue") +
    geom_segment(data = filter(range_shift_sum, year == y), inherit.aes = F,
               aes(x = xmins, xend = xmaxs, y = p5_w), colour = "blue") +
    geom_segment(data = filter(range_shift_sum, year == y), inherit.aes = F,
                 aes(x = xmins, xend = xmaxs, y = mean_max_lat), colour = "green") +
    geom_segment(data = filter(range_shift_sum, year == y), inherit.aes = F,
                 aes(x = xmins, xend = xmaxs, y = mean_min_lat), colour = "green") +
    labs(x = "", y = "")
  
  ggsave(p, path = paste0("outputs/figures/range-edge-shifts/", 
                          str_replace_all(unique(cc$genus_sp), " ", "\\_")), 
         filename = paste0("intraspecific_BBS_", str_replace_all(unique(cc$genus_sp)," ", "\\_"), "_", 
                           y, ".png"), width = 6, height = 3)
  
}

## calculate shift rates by fitting linear regressions 
lm_max <- range_shift_sum %>%
  group_by(genus_sp, polyID) %>%
  filter(n() > 1) %>%
  do(broom::tidy(lm(mean_max_lat ~ year, data = .), conf.int = TRUE)) %>% 
  ungroup()

lm_min <- range_shift_sum %>%
  group_by(genus_sp, polyID) %>%
  filter(n() > 1) %>%
  do(broom::tidy(lm(mean_min_lat ~ year, data = .), conf.int = TRUE)) %>% 
  ungroup()

lm_max$range_edge = "leading"
lm_min$range_edge = "trailing"

lms = rbind(lm_max, lm_min)

lms <- mutate(lms, term = ifelse(term != "(Intercept)", "slope", "intercept")) %>%
  select(term, estimate, genus_sp, range_edge, polyID) %>%
  spread(key = "term", value = "estimate") 

## plot histogram of shifts 
## units are degrees latitude per y 
lms %>%
  ggplot(aes(x = slope, fill = range_edge)) +
  geom_histogram() +
  facet_wrap(~range_edge)


#######################################
##          autocorrelation          ##
#######################################
## read in analysis of Berkeley Earth data
se <- readRDS("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/BerkeleyEarth_noise-colour.rds")

## select data on low spectral exponent measured across entire time period
se =  se %>%
  select(lat, lon, s_spec_exp_PSD_low_all) %>%
  unique()

## create raster from data frame 
se$lon = se$lon - 180
r = se %>%
  select(lon, lat, s_spec_exp_PSD_low_all) %>%
  rast()

plot(r)

## draw buffer around cardinal presence across all years
## get presence points
hull_pts = cc %>%
  filter(!TotalAbd == 0 & !is.na(TotalAbd)) %>%
  select(geometry) %>%
  unique() 

## draw buffered points
buffer = buffer(vect(hull_pts), width = 100000)
plot(buffer)

## dissolve indiivudal polygons
buffer <- st_as_sf(buffer)
buffer = st_union(buffer)
plot(buffer)

## smooth over
buffer = buffer(vect(buffer), width = 100000)
plot(buffer)
plot(hull_pts, add = T)

## crop by North America's edge
buffer_na = crop(buffer, vect(countries))

## cut out middle
buffer_inv = buffer(buffer_na, width = -500000)
plot(buffer_na)
plot(buffer_inv, add = T)


## measure autocorrelation within buffer, separately for each polygon 
plot(r)

r_masked = mask(r, buffer_na) %>%
  mask(., buffer_inv, inv = T)

plot(r_masked)

r_masked_2  = mask(r, buffer_na)

plot(r_masked_2)

## crop buffer by polygons and split into leading and trailing edge
buffer_polys = c()
for(i in 1:nrow(polys)) {
  ## perform crop
  cropped = crop(buffer_na, vect(polys[i,]))
  
  ## if there is a range edge within that polygon  
  if(nrow(cropped) != 0) {
    ## split into leading and trailing edge at mid latitude using bbox
    mid_lat = as.vector(ext(cropped)[3] + ext(cropped)[4]) / 2
    bbox_leading = ext(cropped)
    bbox_leading[3] = mid_lat ## ymin
    bbox_trailing = ext(cropped)
    bbox_trailing[4] = mid_lat ## ymax
    
    leading = crop(cropped, bbox_leading)
    trailing = crop(cropped, bbox_trailing)
    # plot(leading)
    # plot(trailing)
    
    ## combine and turn into sf
    leading = st_as_sf(leading) %>%
      mutate(range_edge = "leading")
    trailing = st_as_sf(trailing) %>%
      mutate(range_edge = "trailing")
    
    cropped = rbind(leading, trailing)
    
    ## save leading and trailing halves, cropped by polygon 
    buffer_polys = cropped %>%
      mutate(polyID = i) %>%
      rbind(buffer_polys, .) 
  }

}

## plot result 
plot(buffer_polys)

## change crs
st_crs(buffer_polys) = st_crs(countries)
crs(r) = crs(vect(buffer_polys))

## calculate autocorrelation in each buffer polygon
mean_betas <- c()
for(i in 1:nrow(buffer_polys)) {
  mean_beta = mean(values(mask(r_masked, vect(buffer_polys[i,]))), na.rm = T)
  
  mean_betas = append(mean_betas, mean_beta)
}
mean_betas

buffer_polys$mean_beta = mean_betas


## join to range shift data
lms_join = left_join(lms, st_drop_geometry(buffer_polys))

## plot
lms_join %>%
  ggplot(aes(x = mean_beta, y = slope, colour = range_edge)) +
  geom_point() +
  geom_smooth(method = "lm")


## do the same thing with synchrony
## read in pairwise synchrony between a route and all routes within a 5 degree radius

sync = read.csv("outputs/data-processed/synchrony_all.csv")

sync = sync %>%
  select(cell1, route1, cell2, pearson_corr) %>%
  distinct() %>%
  group_by(route1) %>%
  summarise(mean_pearson_corr = mean(pearson_corr, na.rm = T)) 

## bind route info to synchrony measurements 
routes <- cc %>%
  filter(!TotalAbd == 0 & !is.na(TotalAbd)) %>% ## get rid of absence data and NAs (i.e., years where routes weren't sampled)
 left_join(., sync, by = c("route" = "route1")) 

plot(routes[13])

## mask by buffers 
routes = mask(vect(routes), buffer_inv, inverse = T)
plot(routes)

## calculate mean synchrony in each polygon at each range edge, just like for autocorrelation
mean_pearsons <- c()
for(i in 1:nrow(buffer_polys)) {
  masked = mask(routes, vect(buffer_polys[i,]))
                
  mean_pearson = mean(masked$mean_pearson_corr, na.rm = T)
  
  mean_pearsons = append(mean_pearsons, mean_pearson)
}
mean_pearsons

buffer_polys$mean_pearsons = mean_pearsons


## join to range shift data
lms_join = left_join(lms_join, st_drop_geometry(buffer_polys))

## plot
lms_join %>%
  ggplot(aes(x = mean_beta, y = slope, colour = mean_pearsons)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~range_edge)

lms_join %>%
  ggplot(aes(x = mean_pearsons, y = slope, colour = mean_pearsons)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~range_edge)


## other species with expansion
## Carolina wren: Thryothorus ludovicianus
## Tufted titmouse: Baeolophus bicolor
## Red-bellied woodpecker: Melanerpes carolinus

## what do now?
## does significance of shifts matter?
sd = range_shift_sum %>%
  group_by(polyID) %>%
  mutate(sd_min = sd(mean_min_lat), 
         sd_max = sd(mean_max_lat)) %>% 
  left_join(st_drop_geometry(buffer_polys))


sd %>%
  filter(range_edge == "leading") %>%
  ggplot(aes(x = mean_beta, y = sd_max, colour = mean_pearsons)) +
  geom_point()
## higher synchrony + autocorrelation = more variable range edge 





## garbage 

## get BOTW range
botw = vect("outputs/data-processed/BBS_BOTW.shp")

## filter to cardinal 
range = botw[which(botw$sci_name == "Cardinalis cardinalis"),]

## buffer BOTW range edge 
buffered_range = buffer(range, width = 1)

plot(range)
plot(buffered_range)

## crop to study area extent 
range = crop(range, extent(cc))

## get max and min longitudes 
ext = ext(range)
xmin = ext[1]
xmax = ext[2]

xmin
xmax

## divide into 10 deg longitude chunks 


