library(tidyverse)
library(terra)
library(tidyterra)

## data notes:
## - only projected for breeding season months (april, may, june, july) 


## read in data
files = list.files("data-raw/short-term-sdms/Short-Term_Climate_Variability_Models_04/Cardinalis_cardinalis/outputs", 
                   full = T)

files = files[which(str_detect(files, "asc"))]


## try first
all = rast(files)

plot(all[[176]])

## see one time series 
pt = vect(data.frame(lon = -120, lat = 40))
crs(pt) = crs(all)

ggplot() +
  geom_spatraster(data = all[[1]]) +
  scale_fill_continuous(na.value = "transparent") +
  geom_spatvector(data = pt, colour = "red")
  
ts <- extract(all, pt)
date = as.numeric(colnames(ts)[-1])
year = paste(substr(date, 1,4))
month = paste(substr(date, 5,6))
date = paste0(substr(date, 1, 4), "-", paste(substr(date, 5,6)), "-01")
date = as_date(date)

ts = data.frame(date = date, suitability = unlist(ts[1,2:ncol(ts)]),
                year = year, month = month)

ts %>% 
  ggplot(aes(x = date, y = suitability)) +
  geom_point() +
  geom_line() +
  geom_smooth(method = "lm", inherit.aes = F, aes(x = date, y = suitability)) 


spectral_exponent_calculator_PSD(ts$suitability, l = nrow(ts))

## add in 0s for months in between breeding season to preseve temporal structure 
key = data.frame(expand.grid(year = as.character(1950:2011), as.character(c(paste0(0, c(1:9)), 10:12))))

new_ts = left_join(key, ts) %>%
  mutate(date = as_date(paste(year, month, "01", sep = "-"))) %>%
  mutate(suitability = ifelse(is.na(suitability), 0, suitability))

new_ts %>% 
  ggplot(aes(x = date, y = suitability)) +
  geom_point() +
  geom_line() +
  geom_smooth(method = "lm", inherit.aes = F, aes(x = date, y = suitability)) 

spec_exp(new_ts$suitability, l = nrow(ts))

## try adding a yearly signal and seeing if spectral exponent changes 
ts2 = new_ts
ts2 = arrange(ts2, date)
ts2$suitability = new_ts$suitability + rep(c(rep(5, 36), rep(0, 36)), 100)[1:length(new_ts$suitability)]

ts2 %>% 
  ggplot(aes(x = date, y = suitability)) +
  geom_point() +
  geom_line() +
  geom_smooth(method = "lm", inherit.aes = F, aes(x = date, y = suitability)) 

spec_exp(ts2$suitability)
# :-) it does


## how much variation in the spectral exponent of suitability is there for the cardinal?
## upscale data for now
all_resampled = aggregate(all, 5, mean)

plot(all[[1]])
plot(all_resampled[[1]])

## turn raster stack to df
cardinal = as.data.frame(all_resampled, xy = TRUE, na.rm = FALSE)

## get sample dates 
date = as.numeric(colnames(cardinal)[-c(1,2)])
year = paste(substr(date, 1,4))
month = paste(substr(date, 5,6))

## add zeroes for months with no data
key = data.frame(expand.grid(year = as.character(1950:2011), month = as.character(c(paste0(0, c(1:9)), 10:12))))

## change col names 
colnames(cardinal)[3:ncol(cardinal)] = paste0("c", colnames(cardinal)[3:ncol(cardinal)])

## get rid of cells with all na 
cardinal = cardinal %>%
  filter(if_any(c195004:c201107, ~ !is.na(.)))

## calculate monthly climatology in each cell
april = which(substr(colnames(cardinal), 6,7) == "04")
april_clims = sapply(1:nrow(cardinal), FUN = function(x) {
  mean(as.numeric(cardinal[x,april]), na.rm = T)
})
may = which(substr(colnames(cardinal), 6,7) == "05")
may_clims = sapply(1:nrow(cardinal), FUN = function(x) {
  mean(as.numeric(cardinal[x,may]), na.rm = T)
})
june = which(substr(colnames(cardinal), 6,7) == "06")
june_clims = sapply(1:nrow(cardinal), FUN = function(x) {
  mean(as.numeric(cardinal[x,june]), na.rm = T)
})
july = which(substr(colnames(cardinal), 6,7) == "07")
july_clims = sapply(1:nrow(cardinal), FUN = function(x) {
  mean(as.numeric(cardinal[x,july]), na.rm = T)
})

i = 1
se = se_detrended = se_zeros = suitability_trends = c()
while(i <= nrow(cardinal)) {
  
  ## get suitability
  ts = as.numeric(cardinal[i,3:ncol(cardinal)])

  ts = data.frame(date = date, suitability = ts,
                  year = year, month = month)
  
  ## plot
  # ts %>%
  #   mutate(date = as_date(paste(year, month, "01", sep = "-"))) %>%
  #   ggplot(aes(x = year, y = suitability)) +
  #   geom_point()

  ## if not all na 
  if(!all(is.na(ts$suitability))) {
    
    ## calculate spectral exponent
    se = append(se, spec_exp(ts$suitability))
    
    ## subtract climatology, linearly detrend, and calculate spectral exponent
    ## get climatology 
    clim = rep(c(april_clims[i],may_clims[i],june_clims[i], july_clims[i]), length.out = length(date))
    
    ## plot
    # data.frame(year = year, month = month, clim = clim) %>%
    #   mutate(date = as_date(paste(year, month, "01", sep = "-"))) %>%
    #   ggplot(aes(x = date, y = clim)) +
    #   geom_point()
    
    ## subtract from suitability
    ts$suitability_detrended = ts$suitability - clim
    
    ## plot
    # ts%>%
    #   mutate(date = as_date(paste(year, month, "01", sep = "-"))) %>%
    #   ggplot(aes(x = date, y = suitability_detrended)) +
    #   geom_point()
    
    ## linearly detrend by fitting linear model and extracting residuals
    ## save slope and intercept
    lm = lm(suitability_detrended ~ date,
            data = ts)
    suitability_trends = rbind(suitability_trends, tidy(lm) %>% mutate(i = i)) 
    
    ## plot
    # ts%>%
    #   ggplot(aes(x = date, y = suitability_detrended)) +
    #   geom_point() +
    #   geom_abline(slope = lm$coefficients[[2]], intercept = lm$coefficients[[1]], colour = "red")
    
    ## set residuals as suitability
    ts$suitability_detrended = lm$residuals
    
    ## calculate spectral exponent
    se_detrended = append(se_detrended, spec_exp(ts$suitability_detrended))
    
    # ## plot
    # ts%>%
    #   mutate(date = as_date(paste(year, month, "01", sep = "-"))) %>%
    #   ggplot(aes(x = date, y = suitability_detrended)) +
    #   geom_point()
    
    # add 0s for months not sampled
    ts = left_join(key, ts, by = c("year", "month")) %>%
      mutate(date = as_date(paste(year, month, "01", sep = "-"))) %>%
      mutate(suitability = ifelse(is.na(suitability), 0, suitability)) %>%
      arrange(date)
    
    ## calculate spectral exponent
    se_zeros = append(se_zeros, spec_exp(ts$suitability))
  } 
  else {
    se = append(se, NA)
    se_detrended = append(se_detrended, NA)
    se_zeros = append(se_zeros, NA)
    suitability_trends = rbind(suitability_trends, data.frame(i = i, 
                                                              term = c(NA,NA),
                                                              estimate = c(NA,NA),
                                                              std.error = c(NA,NA),
                                                              statistic = c(NA,NA),
                                                              p.value = c(NA,NA))) 
  }
  
  i = i + 1
  print(i)
}
hist(se)
hist(se_detrended)
hist(se_zeros)




## 1. plot se of suitability on a map - are there spatial patterns? is it just the places with highest suitability?
## 2. look at one with se = 0 - does it look like white noise? 
## 3. what about for another species? how different is this?


## plot for now
temp = cardinal
temp$se = se
temp$se_detrended = se_detrended
temp$se_zeros = se_zeros

temp = dplyr::select(temp, x, y, se, se_detrended, se_zeros)

## add suitability trend 
trends = suitability_trends %>%
  filter(term == "date")
temp$trends = trends$estimate

r_temp = rast(temp)
plot(-r_temp$se)
plot(-r_temp$se_detrended)
plot(-r_temp$se_zeros)
plot(r_temp$trends)

plot(all[[1]])

## try for another species 





## another species 
## read in data
files = list.files("data-raw/short-term-sdms/Short-Term_Climate_Variability_Models_10/Passerella_iliaca/outputs", 
                   full = T)
files = files[which(str_detect(files, "\\.asc"))]

## check that all can be read properly
no_errors = c()
for(i in 1:length(files)) {
  tryCatch({
    ## attempt to read file and access data
    raster_data <- rast(files[i])
    
    summary(values(raster_data))
    
    no_errors = append(no_errors, i)
  }, error = function(e) {
    print(paste("Error reading file:", i))
  })
}

## get rid of files with errors 
files = files[no_errors]

## read in rasters 
all = rast(files)

plot(all[[100]])

## upscale data for now
all_resampled = aggregate(all, 5, mean)

plot(all[[1]])
plot(all_resampled[[1]])

## see one time series 
pt = vect(data.frame(lon = -115.125, lat = 33.16667))
crs(pt) = crs(all)

ggplot() +
  geom_spatraster(data = all[[1]]) +
  scale_fill_continuous(na.value = "transparent") +
  geom_spatvector(data = pt, colour = "red")

ts <- extract(all_resampled, pt)
date = as.numeric(colnames(ts)[-1])
year = paste(substr(date, 1,4))
month = paste(substr(date, 5,6))
date = paste0(substr(date, 1, 4), "-", paste(substr(date, 5,6)), "-01")
date = as_date(date)

ts = data.frame(date = date, suitability = unlist(ts[1,2:ncol(ts)]),
                year = year, month = month)

ts %>% 
  ggplot(aes(x = date, y = suitability)) +
  geom_point() +
  geom_line() +
  geom_smooth(method = "lm", inherit.aes = F, aes(x = date, y = suitability)) 

## add in 0s for months in between breeding season to preseve temporal structure 
key = data.frame(expand.grid(year = as.character(1950:2011), as.character(c(paste0(0, c(1:9)), 10:12))))

new_ts = left_join(key, ts) %>%
  mutate(date = as_date(paste(year, month, "01", sep = "-"))) %>%
  mutate(suitability = ifelse(is.na(suitability), 0, suitability))

new_ts %>% 
  ggplot(aes(x = date, y = suitability)) +
  geom_point() +
  geom_line() +
  geom_smooth(method = "lm", inherit.aes = F, aes(x = date, y = suitability)) 

spec_exp(new_ts$suitability)

## turn raster stack to df
df = as.data.frame(all_resampled, xy = TRUE, na.rm = FALSE)

## get sample dates 
date = as.numeric(colnames(df)[-c(1,2)])
year = paste(substr(date, 1,4))
month = paste(substr(date, 5,6))

## add zeroes for months with no data
key = data.frame(expand.grid(year = as.character(1950:2011), month = as.character(c(paste0(0, c(1:9)), 10:12))))

## change col names 
colnames(df)[3:ncol(df)] = paste0("c", colnames(df)[3:ncol(df)])

## get rid of cells with all na 
df = df %>%
  filter(if_any(c195004:c201107, ~ !is.na(.)))

## calculate monthly climatology in each cell
april = which(substr(colnames(df), 6,7) == "04")
april_clims = sapply(1:nrow(df), FUN = function(x) {
  mean(as.numeric(df[x,april]), na.rm = T)
})
may = which(substr(colnames(df), 6,7) == "05")
may_clims = sapply(1:nrow(df), FUN = function(x) {
  mean(as.numeric(df[x,may]), na.rm = T)
})
june = which(substr(colnames(df), 6,7) == "06")
june_clims = sapply(1:nrow(df), FUN = function(x) {
  mean(as.numeric(df[x,june]), na.rm = T)
})
july = which(substr(colnames(df), 6,7) == "07")
july_clims = sapply(1:nrow(df), FUN = function(x) {
  mean(as.numeric(df[x,july]), na.rm = T)
})

i = 1
se = se_detrended = se_zeros = suitability_trends = c()
while(i <= nrow(df)) {
  
  ## get suitability
  ts = as.numeric(df[i,3:ncol(df)])
  
  ts = data.frame(date = date, suitability = ts,
                  year = year, month = month)
  
  ## plot
  # ts %>%
  #   mutate(date = as_date(paste(year, month, "01", sep = "-"))) %>%
  #   ggplot(aes(x = year, y = suitability)) +
  #   geom_point()

  ## if not all na 
  if(!all(is.na(ts$suitability))) {
    
    ## calculate spectral exponent
    se = append(se, spec_exp(ts$suitability))
    
    ## subtract climatology, linearly detrend, and calculate spectral exponent
    ## get climatology 
    clim = rep(c(april_clims[i],may_clims[i],june_clims[i], july_clims[i]), length.out = length(date))
    
    ## plot
    # data.frame(year = year, month = month, clim = clim) %>%
    #   mutate(date = as_date(paste(year, month, "01", sep = "-"))) %>%
    #   ggplot(aes(x = date, y = clim)) +
    #   geom_point()

    ## subtract from suitability
    ts$suitability_detrended = ts$suitability - clim
    
    ## plot
    # ts%>%
    #   mutate(date = as_date(paste(year, month, "01", sep = "-"))) %>%
    #   ggplot(aes(x = date, y = suitability_detrended)) +
    #   geom_point()
    
    ## linearly detrend by fitting linear model and extracting residuals
    ## save slope and intercept
    lm = lm(suitability_detrended ~ date,
            data = ts)
    suitability_trends = rbind(suitability_trends, tidy(lm) %>% mutate(i = i)) 
    
    ## plot
    # ts%>%
    #   ggplot(aes(x = date, y = suitability_detrended)) +
    #   geom_point() +
    #   geom_abline(slope = lm$coefficients[[2]], intercept = lm$coefficients[[1]], colour = "red")
    
    ## set residuals as suitability
    ts$suitability_detrended = lm$residuals
    
    ## calculate spectral exponent
    se_detrended = append(se_detrended, spec_exp(ts$suitability_detrended))
    
    # ## plot
    # ts%>%
    #   mutate(date = as_date(paste(year, month, "01", sep = "-"))) %>%
    #   ggplot(aes(x = date, y = suitability_detrended)) +
    #   geom_point()
    
    # add 0s for months not sampled
    ts = left_join(key, ts, by = c("year", "month")) %>%
      mutate(date = as_date(paste(year, month, "01", sep = "-"))) %>%
      mutate(suitability = ifelse(is.na(suitability), 0, suitability)) %>%
      arrange(date)
    
    ## calculate spectral exponent
    se_zeros = append(se_zeros, spec_exp(ts$suitability))
  } 
  else {
    se = append(se, NA)
    se_detrended = append(se_detrended, NA)
    se_zeros = append(se_zeros, NA)
    suitability_trends = rbind(suitability_trends, data.frame(i = i, 
                                                              term = c(NA,NA),
                                                              estimate = c(NA,NA),
                                                              std.error = c(NA,NA),
                                                              statistic = c(NA,NA),
                                                              p.value = c(NA,NA))) 
  }
  
  i = i + 1
  print(i)
}
hist(se)
hist(se_detrended)
hist(se_zeros)

## look at one with high and low autocorrelation - do they look as expected?
## what to do with standrd deviation?


## 1. plot se of suitability on a map - are there spatial patterns? is it just the places with highest suitability?
## 2. look at one with se = 0 - does it look like white noise? 
## 3. what about for another species? how different is this?


## plot for now
temp = df
temp$se = se
temp$se_detrended = se_detrended
temp$se_zeros = se_zeros

temp = dplyr::select(temp, x, y, se, se_detrended, se_zeros)

## add suitability trend 
trends = suitability_trends %>%
  filter(term == "date")
temp$trends = trends$estimate

r_temp = rast(temp)
plot(-r_temp$se)
plot(-r_temp$se_detrended)
plot(-r_temp$se_zeros)
plot(r_temp$trends)

plot(all[[1]])

max(r_temp$trends)

#### NEXT:
### 1. measure latitudinal shift in suitability northern edge for a few species - do we see a trend?
## (to do this, threshold the sdm, turn to presence absence?)
### 2. compare it to latitudinal BBS range shift - does it make sense? 
### 3. write code to calculate trend and spectral exponent of suitability in cells at the range edge 


## plot breeding bird data on top 
bbs <- read.csv("outputs/data-processed/bbs_clean-subset.csv") %>%
  filter(genus_sp == "Passerella iliaca")

bbs <- filter(bbs, TotalAbd != 0)

all_sub = all[[first(which(substr(names(all), 1, 4) == "1980")):nlyr(all)]]


## plot each year in June with occurrences on top
for(x in 1:length(unique(substr(names(all_sub), 1, 4)))) {
  year = paste0(unique(substr(names(all_sub), 1, 4))[x], "07")
  
  bbs_yr = filter(bbs, Year == unique(substr(names(all_sub), 1, 4))[x])
  bbs_yr <- vect(bbs_yr)
  crs(bbs_yr) = crs(all)
  
  i = which(names(all_sub) == year)
  
  plot = ggplot() +
    geom_spatraster(data = all_sub[[i]]) +
    geom_spatvector(data = bbs_yr, size = 0.25, colour = "red") +
    scale_fill_continuous(na.value = "transparent") +
    scale_x_continuous(limits = c(-125,-65)) +
    scale_y_continuous(limits = c(25,50)) +
    labs(title = unique(substr(names(all_sub), 1, 4))[x]) +
    scale_fill_continuous(limits = c(0, 1))
  
  ggsave(plot, path = "outputs/figures/sdm-vs-bbs/", filename = paste0("Passerella_iliaca_", year, ".png"),
         width = 4, height = 3)
}





##############################################
#####              FUNCTIONS:           ######
##############################################
spec_exp = function(noise) {
  
  ## estimate noise colour from a linear regression of power spectrum:
  l <- length(noise)
  dft <- fft(noise)/l
  amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
  amp <- amp[1:(l/2)]	## remove second half of amplitudes (negative half)
  freq <- 1:(l/2)/l ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 
  
  ## create periodogram data by squaring amplitude of FFT output
  spectral <- data.frame(freq = freq, power = amp^2)
  
  ggplot2::ggplot(data = spectral, ggplot2::aes(x = freq, y = power)) + ggplot2::geom_line() +
    ggplot2::scale_y_log10() + ggplot2::scale_x_log10() + ggplot2::geom_smooth(method = "lm") +
    ggplot2::theme_minimal()
  
  true_colour <- lm(data = spectral, log(power) ~ log(freq))
  
  ## save 

  return(as.numeric(true_colour$coefficients[2]))
  
}
