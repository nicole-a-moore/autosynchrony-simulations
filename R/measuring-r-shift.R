## try measuring shift in latitudinal suitability across a species range 
library(sf)
library(raster)


## 1. get a species' range polygon - Cardinalis cardinalis
## read in botw ranges 
botw = st_read("outputs/data-processed/BBS_BOTW.shp")

#botw = botw[which(botw$sci_name == "Cardinalis cardinalis"),]
botw = botw[which(botw$sci_name == "Calamospiza melanocorys"),]

plot(botw)

## 2. create a theoretical response curve that describes how suitability (r) changes with temperature
## read in Berkeley Earth air temperature data as 3d array
be = stack("outputs/data-processed/BerkeleyEarth_1966_NorthAmerica.tif")
be_array = tiff::readTIFF("outputs/data-processed/BerkeleyEarth_1966_NorthAmerica.tif")
#plot(rast(be_array[,,1]))

## see distribution temps across the species range in 1966
## mask by species range 

range_rast = rasterize(botw, be[[1]])
#plot(range_rast)

range_array = as.array(range_rast)[,,1]
#plot(rast(range_array[,,1]))
range_array = replicate(dim(be_array)[3], range_array)

## set cells not in species range to na
sp_array = be_array
sp_array[which(is.na(range_array))] <- NA
#plot(rast(sp_array[,,10]))

## plot temps across species range in 1966
temps = sp_array[,,1:365]
 
hist(temps)
min(temps, na.rm = T)
max(temps, na.rm = T)
mean(temps, na.rm = T)

## theoretical optimum: 5 degrees 
## max suitability: 1.2
## min suitability: -0.25
## skewed response 
curve = fGarch::dsnorm(-60:60, mean = 10, sd = 9, xi = 1/1.5, log = FALSE)
max = max(curve)
curve = curve/max*(1.2+0.25)-0.25
plot(x = -60:60, y = curve) 

## see if extremes match the species distribution
min((-60:60)[which(curve > 0)])
max((-60:60)[which(curve > 0)])
## good enough for now 

## 3. plot change in latitudinal position of mean, 5th and 95th percentile of suitability > 0 across latitude across time

## translate temperature to suitability using curve 
suit_array = be_array
max = max(fGarch::dsnorm(suit_array[,,], mean = 12, sd = 10, xi = 1/1.5, log = FALSE), na.rm = T)
suit_array[,,] <- fGarch::dsnorm(suit_array[,,], mean = 12, sd = 10, xi = 1/1.5, log = FALSE)/max*(1.2+0.25)-0.25
plot(rast(be_array[,,1]))
plot(rast(suit_array[,,1]))

## for each day, plot histogram of mean suitability > 0 across latitude 
r = be[[1]]
values(r) = 1
coords = as.data.frame(rasterToPoints(r))
rownames(suit_array) = as.numeric(unique(coords$y))
colnames(suit_array) = as.numeric(unique(coords$x))

for(t in 1:dim(suit_array)[3]) {
  if(t == 1) {
    df_suit <- data.frame(time = t, 
                          lat = as.numeric(rownames(suit_array)),
                          mean_suit = sapply(1:nrow(suit_array), FUN = function (x) {
                            mean(suit_array[x,1:ncol(suit_array),t], na.rm = T)}))
  }
  else {
    df_suit <- rbind(df_suit,  data.frame(time = t, 
                                          lat = as.numeric(rownames(suit_array)),
                                          mean_suit = sapply(1:nrow(suit_array), FUN = function (x) {
                                            mean(suit_array[x,1:ncol(suit_array),t], na.rm = T)})))
  }
  print(t)
}

df_suit %>%
  filter(time == 100) %>%
  ggplot(aes(x = lat, y = mean_suit)) +
  geom_point() +
  geom_hline(yintercept = 0, colour = "red")

## plot shift in latitude of mean/max suitability > 0 over time
df_suit %>%
  group_by(time) %>%
  filter(mean_suit > 0) %>%
  summarise(max_suit = max(mean_suit)) %>%
  ungroup() %>%
  left_join(., df_suit, by = c("max_suit" = "mean_suit", "time")) %>%
  mutate(group = ifelse(time <= 5110, "1980", "after")) %>%
  ggplot(aes(x = time, y = lat)) +
  geom_line() +
  geom_smooth(method = "lm", aes(group = group)) + 
  scale_y_continuous(limits = c(min(df_suit$lat), max(df_suit$lat))) 

## calculate in windows of 100 time points
for(t in 50:(dim(suit_array)[3]-50)) {
  if(t == 50) {
    df_suit_mw = data.frame(time = t, 
               lat = as.numeric(rownames(suit_array)),
               mean_suit = sapply(1:nrow(suit_array), FUN = function (x) {
                 mean(suit_array[x,1:ncol(suit_array),(t-49):(t+49)], na.rm = T)}))
  }
  else {
    df_suit_mw = rbind(df_suit_mw, 
                       data.frame(time = t, 
                                  lat = as.numeric(rownames(suit_array)),
                                  mean_suit = sapply(1:nrow(suit_array), FUN = function (x) {
                                    mean(suit_array[x,1:ncol(suit_array),(t-49):(t+49)], na.rm = T)})))
  }
  print(t)
}

df_suit_mw %>%
  filter(time == 100) %>%
  ggplot(aes(x = lat, y = mean_suit)) +
  geom_point() +
  geom_hline(yintercept = 0, colour = "red")

## plot shift in latitude of mean/max/min/peak suitability > 0 over time
df_suit_mw %>%
  group_by(time) %>%
  summarise(p95 = quantile(.$lat[which(mean_suit >= 0)], 0.95),
            p5 = quantile(.$lat[which(mean_suit >= 0)], 0.05),
            peak_r = mean(.$lat[which(mean_suit == max(mean_suit, na.rm = T))])) %>%
  gather(key = "measurement", value = "lat", c(p95, p5, peak_r)) %>%
  mutate(group = paste0(measurement, "_", ifelse(time <= 5110, "1980", "after"))) %>%
  ggplot(aes(x = time, y = lat, colour = measurement, group = measurement)) +
  geom_line(alpha = 0.5) +
  geom_smooth(method = "lm", aes(group = group)) +
  scale_y_continuous(limits = c(min(df_suit_mw$lat), max(df_suit_mw$lat))) 

## calculate rate of shift in parameters 
params = df_suit_mw %>%
  group_by(time) %>%
  summarise(p95 = quantile(.$lat[which(mean_suit >= 0)], 0.95),
            p5 = quantile(.$lat[which(mean_suit >= 0)], 0.05),
            peak_r = mean(.$lat[which(mean_suit == max(mean_suit, na.rm = T))])) 

dates = seq(as.Date("1966/1/1"), as.Date("2021/12/31"), "days")
dates = dates[!str_detect(dates,"-02-29")]

params$date = seq(as.Date("1966/1/1"), as.Date("2021/12/31"), "days")[50:(length(dates)-50)]

## fit lm
peak_lm = lm(peak_r ~ date,
             data = params)

p95_lm = lm(p95 ~ date,
            data = params)

p5_lm = lm(p5 ~ date,
          data = params)

summary(peak_lm)
## 0.00007755 degrees latitude per day 
summary(p95_lm)
## 0.00007983 degrees latitude per day 
summary(p5_lm)
## 0.00001008 degrees latitude per day 

## fit to baseline and shifting period separately
bl = params[1:which(params$date == "1980-01-01"),]
sh = params[(which(params$date == "1980-01-01")+1):nrow(params),]

bl_peak_lm = lm(peak_r ~ date,
             data = bl)

bl_p95_lm = lm(p95 ~ date,
            data = bl)

bl_p5_lm = lm(p5 ~ date,
           data = bl)

summary(bl_peak_lm)
summary(bl_p95_lm)
summary(bl_p5_lm)

sh_peak_lm = lm(peak_r ~ date,
                data = sh)

sh_p95_lm = lm(p95 ~ date,
               data = sh)

sh_p5_lm = lm(p5 ~ date,
              data = sh)

summary(sh_peak_lm)
summary(sh_p95_lm)
summary(sh_p5_lm)



## separate by east west and see if pattern differs
## divide at the 100th meridian 
west = which(colnames(suit_array) == "-100.5")
east = which(colnames(suit_array) == "-99.5")

## calculate in windows of 100 time points
for(t in 50:(max(df_suit$time)-50)) {
  if(t == 50) {
    df_suit_region = data.frame(time = t, 
                            lat = as.numeric(rownames(suit_array)),
                            mean_suit = sapply(1:nrow(suit_array), FUN = function (x) {
                              mean(suit_array[x,east:ncol(suit_array),(t-49):(t+49)], na.rm = T)}),
                            region = "east")
    
    df_suit_region = rbind(df_suit_region, 
                       data.frame(time = t, 
                                  lat = as.numeric(rownames(suit_array)),
                                  mean_suit = sapply(1:nrow(suit_array), FUN = function (x) {
                                    mean(suit_array[x,1:west,(t-49):(t+49)], na.rm = T)}),
                                  region = "west"))
    
  }
  else {
    df_suit_region = rbind(df_suit_region, 
                       data.frame(time = t, 
                                  lat = as.numeric(rownames(suit_array)),
                                  mean_suit = sapply(1:nrow(suit_array), FUN = function (x) {
                                    mean(suit_array[x,east:ncol(suit_array),(t-49):(t+49)], na.rm = T)}),
                                  region = "east"))
    
    df_suit_region = rbind(df_suit_region, 
                       data.frame(time = t, 
                                  lat = as.numeric(rownames(suit_array)),
                                  mean_suit = sapply(1:nrow(suit_array), FUN = function (x) {
                                    mean(suit_array[x,1:west,(t-49):(t+49)], na.rm = T)}),
                                  region = "west"))
  }
  print(t)
}



## plot shift in latitude of mean/max/min/peak suitability > 0 over time
df_suit_region %>%
  group_by(time, region) %>%
  summarise(p95 = quantile(.$lat[which(mean_suit >= 0)], 0.95),
            p5 = quantile(.$lat[which(mean_suit >= 0)], 0.05),
            peak_r = mean(.$lat[which(mean_suit == max(mean_suit, na.rm = T))])) %>%
  gather(key = "measurement", value = "lat", c(p95, p5, peak_r)) %>%
  mutate(group = paste0(measurement, "_", ifelse(time <= 5110, "1980", "after"))) %>%
  ggplot(aes(x = time, y = lat, colour = measurement, group = measurement)) +
  geom_line(alpha = 0.5) +
  geom_smooth(method = "lm", aes(group = group)) +
  scale_y_continuous(limits = c(min(df_suit_region$lat), max(df_suit_region$lat))) +
  facet_wrap(~region)


## what about just within the species range 
## translate temperature to suitability using curve 
suit_array = sp_array
max = max(fGarch::dsnorm(suit_array[,,], mean = 12, sd = 10, xi = 1/1.5, log = FALSE), na.rm = T)
suit_array[,,] <- fGarch::dsnorm(suit_array[,,], mean = 12, sd = 10, xi = 1/1.5, log = FALSE)/max*(1.2+0.25)-0.25
plot(rast(be_array[,,1]))
plot(rast(suit_array[,,1]))

## for each day, plot histogram of mean suitability > 0 across latitude 
r = be[[1]]
values(r) = 1
coords = as.data.frame(rasterToPoints(r))
rownames(suit_array) = as.numeric(unique(coords$y))
colnames(suit_array) = as.numeric(unique(coords$x))

## calculate in windows of 100 time points
for(t in 50:(max(df_suit$time)-50)) {
  if(t == 50) {
    df_suit_range = data.frame(time = t, 
                            lat = as.numeric(rownames(suit_array)),
                            mean_suit = sapply(1:nrow(suit_array), FUN = function (x) {
                              mean(suit_array[x,1:ncol(suit_array),(t-49):(t+49)], na.rm = T)}))
  }
  else {
    df_suit_range = rbind(df_suit_range, 
                       data.frame(time = t, 
                                  lat = as.numeric(rownames(suit_array)),
                                  mean_suit = sapply(1:nrow(suit_array), FUN = function (x) {
                                    mean(suit_array[x,1:ncol(suit_array),(t-49):(t+49)], na.rm = T)})))
  }
  print(t)
}

## plot shift in latitude of mean/max/min/peak suitability > 0 over time
df_suit_range %>%
  group_by(time) %>%
  summarise(p95 = quantile(.$lat[which(mean_suit >= 0)], 0.95),
            p5 = quantile(.$lat[which(mean_suit >= 0)], 0.05),
            peak_r = mean(.$lat[which(mean_suit == max(mean_suit, na.rm = T))])) %>%
  gather(key = "measurement", value = "lat", c(p95, p5, peak_r)) %>%
  mutate(group = paste0(measurement, "_", ifelse(time <= 5110, "1980", "after"))) %>%
  ggplot(aes(x = time, y = lat, colour = measurement, group = measurement)) +
  geom_line(alpha = 0.5) +
  geom_smooth(method = "lm", aes(group = group)) +
  scale_y_continuous(limits = c(min(df_suit_region$lat), max(df_suit_region$lat))) 


## maybe average just across longitudes where the species is present?


## compare to magnitude of range edge shift
lms <- read.csv("outputs/data-processed/bbs-range-shift-rates.csv")
card = lms[which(lms$genus_sp == "Calamospiza melanocorys"),]

## leading edge 
card

card$slope[which(card$range_edge == "leading")]
## 0.02965894 deg latitude per year
## compare:
p95_lm$coefficients
## convert to same units
0.00007982575*365
## 0.0291364 omg HAHA this has to be coincidence

card$slope[which(card$range_edge == "trailing")]
## -0.009154661 deg latitude per year
## compare:
p5_lm$coefficients
## convert to same units
0.00001008522*365
## 0.003681105


## plot 
c = bb %>%
  filter(genus_sp == "Calamospiza melanocorys") %>%
  filter(!TotalAbd == 0 & !is.na(TotalAbd)) %>% ## get rid of absence data and NAs (i.e., years where routes weren't sampled)
  group_by(Year) %>%
  summarise(p95_w = modi::weighted.quantile(Latitude, w = TotalAbd, 0.95),
            p5_w = modi::weighted.quantile(Latitude, w = TotalAbd, 0.05),
            p95 = quantile(Latitude, 0.95), 
            p5 = quantile(Latitude, 0.05)) %>% ## calculate range edge metrics for that year       
  ungroup() %>% ## ungroup
  gather(key = "range_edge", value = "latitude", c(p95_w, p95, p5_w, p5)) %>%
  mutate(date = as.Date(paste(Year, "-06-01", sep = "")))

c %>%
  ggplot(aes(x = date, y = latitude, colour = range_edge)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", se = F, size = 0.2)

df_time = data.frame(time = 1:20440, date = dates)

df_suit_year %>%
  left_join(., df_time) %>%
  group_by(time) %>%
  summarise(p95 = quantile(.$lat[which(mean_suit >= 0.65)], 0.95),
            p5 = quantile(.$lat[which(mean_suit >= 0.65)], 0.05),
            peak_r = mean(.$lat[which(mean_suit == max(mean_suit, na.rm = T))]), 
            date = date) %>%
  gather(key = "measurement", value = "lat", c(p95, p5, peak_r)) %>%
  mutate(group = paste0(measurement, "_", ifelse(time <= 5110, "1980", "after"))) %>%
  ggplot(aes(x = date, y = lat, colour = measurement, group = measurement)) +
  geom_line(alpha = 0.5) +
  geom_smooth(method = "lm") +
  scale_y_continuous(limits = c(min(df_suit_mw$lat), max(df_suit_mw$lat))) +
  geom_point(data = c, aes(x = date, y = latitude, colour = range_edge), inherit.aes = F, size = 0.5) 






## calculate in windows of 100 time points
for(t in 183:(max(df_suit$time)-182)) {
  if(t == 183) {
    df_suit_year = data.frame(time = t, 
                               lat = as.numeric(rownames(suit_array)),
                               mean_suit = sapply(1:nrow(suit_array), FUN = function (x) {
                                 mean(suit_array[x,1:ncol(suit_array),(t-182):(t+182)], na.rm = T)}))
  }
  else {
    df_suit_year = rbind(df_suit_year, 
                          data.frame(time = t, 
                                     lat = as.numeric(rownames(suit_array)),
                                     mean_suit = sapply(1:nrow(suit_array), FUN = function (x) {
                                       mean(suit_array[x,1:ncol(suit_array),(t-182):(t+182)], na.rm = T)})))
  }
  print(t)
}






#### NEXT STEPS: ######
## - find a way to solidly estimate a temperature response curve for each species 
## - figure out where across North America to measure the suitability shift for each species 
## - keep working on filtering sp list 
##    - change to Martins et al. sp list instead of limiting to songbirds?
##    - independently calculate proportion of range covered by bbs or use Martins et al. 
## - keep working on deciding start_dates





