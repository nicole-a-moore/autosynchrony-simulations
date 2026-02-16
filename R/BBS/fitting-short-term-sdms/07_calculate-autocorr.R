## calculate autocorrelation of suitability for each species from sdm predictions
library(tidyverse)
library(terra)

###########################################
##          measure autocorrelation      ##
###########################################
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
  
  ## read in shape representing the buffer around species range edge 
  buffer = vect(paste0("outputs/data-processed/range-edge-buffers/", str_replace_all(cur_sp, " ", "_"), "_buffer.shp"))
  
  ## get suitability predictions
  ## read in sdm predictions 
  pred = rast(files[sp])
  
  ## get rid of data before 1980
  pred = pred[[which(names(pred) == "date_1980_4"):nlyr(pred)]]
  
  ## change extent
  pred = crop(pred, c(-170, -50, 24.5, 88))

  # ggplot() +
  #   geom_spatraster(data = pred[[1]]) +
  #   scale_fill_continuous(limits = c(0,1), na.value = "transparent")
  
  ## crop by buffer
  pred <- mask(pred, buffer)
  #plot(pred[[1]])
  
  ## for each buffer polygon
  poly = 1
  xy = c()
  while(poly <= length(buffer)) {
    ## crop by buffer
    cur_rast = crop(pred, buffer[poly,])
    #plot(cur_rast)
    
    ## convert to xy data frame
    cur_rast = as.data.frame(cur_rast, xy = TRUE)
    
    ## create column for leading/trailing
    cur_rast$buffer_edge = buffer$range_edge[poly]
    
    ## join to all data
    xy = rbind(xy, cur_rast)
    
    poly = poly + 1
  }
  
  ## assign each cell to a longitudinal bin
  bin = c()
  for(i in 1:nrow(xy)) {
    bin = append(bin, ifelse(length(which(xy$x[i] >= lon_bins$min & xy$x[i] <= lon_bins$max)) == 0, 
                             NA, 
                             which(xy$x[i] >= lon_bins$min & xy$x[i] <= lon_bins$max)))
  }
  xy$bin = bin
  xy = left_join(xy, lon_bins, by = "bin")
  
  ## get rid of data outside of longitudinal bands
  xy <- xy[!is.na(xy$Longitude_bin),]
  
  ## calculate spectral exponent of each time series 
  #####################################################
  ## calculate monthly climatology in each cell
  april = which(substr(colnames(xy), 11,11) == "4")
  april_clims = sapply(1:nrow(xy), FUN = function(x) {
    mean(as.numeric(xy[x,april]), na.rm = T)
  })
  may = which(substr(colnames(xy), 11,11) == "5")
  may_clims = sapply(1:nrow(xy), FUN = function(x) {
    mean(as.numeric(xy[x,may]), na.rm = T)
  })
  june = which(substr(colnames(xy), 11,11) == "6")
  june_clims = sapply(1:nrow(xy), FUN = function(x) {
    mean(as.numeric(xy[x,june]), na.rm = T)
  })
  july = which(substr(colnames(xy), 11,11) == "7")
  july_clims = sapply(1:nrow(xy), FUN = function(x) {
    mean(as.numeric(xy[x,july]), na.rm = T)
  })
  
  i = 1
  se_detrended = suitability_trends = sd_suitability = c()
  while(i <= nrow(xy)) {
    
    ## get suitability
    ts = as.numeric(xy[i,3:(ncol(xy)-5)])
    
    ts = data.frame(date = str_split_fixed(colnames(xy)[3:(ncol(xy)-5)], "_", 2)[,2], suitability = ts,
                    year = str_split_fixed(colnames(xy)[3:(ncol(xy)-5)], "_", 3)[,2], 
                    month = str_split_fixed(colnames(xy)[3:(ncol(xy)-5)], "_", 3)[,3]) %>%
      mutate(date = as_date(paste0(year, "-0", month, "-01")))
    
    ## plot
    # ts %>%
    #   ggplot(aes(x = date, y = suitability)) +
    #   geom_point()
    
    ## if not all na 
    if(!all(is.na(ts$suitability))) {
      
      ## subtract climatology, linearly detrend, and calculate spectral exponent
      ## get climatology 
      clim = rep(c(april_clims[i],may_clims[i],june_clims[i], july_clims[i]), length.out = length(ts$date))
      
      ## plot
      # data.frame(date = ts$date, clim = clim) %>%
      #   ggplot(aes(x = date, y = clim)) +
      #   geom_line() +
      #   scale_y_continuous(limits = c(0, 1))
      
      ## subtract from suitability
      ts$suitability_sdetrended = ts$suitability - clim
      
      ## plot
      # ts%>%
      #   ggplot(aes(x = date, y = suitability_detrended)) +
      #   geom_point() 

      ## linearly detrend by fitting linear model and extracting residuals
      ## save slope and intercept
      lm = lm(suitability_sdetrended ~ date,
              data = ts)
      suitability_trends = rbind(suitability_trends, tidy(lm) %>% mutate(i = i)) 
      
      ## plot
      ts%>%
        ggplot(aes(x = date, y = suitability_sdetrended)) +
        geom_point() +
        geom_abline(slope = lm$coefficients[[2]], intercept = lm$coefficients[[1]], colour = "red")
      
      ## set residuals as suitability
      ts$suitability_ldetrended = lm$residuals
      
      ## plot
      ts%>%
        ggplot(aes(x = date, y = suitability_ldetrended)) +
        geom_point() 
      
      ## calculate spectral exponent
      se_detrended = append(se_detrended, spec_exp(ts$suitability_ldetrended))
      
      # ## plot
      # ts%>%
      #   mutate(date = as_date(paste(year, month, "01", sep = "-"))) %>%
      #   ggplot(aes(x = date, y = suitability_detrended)) +
      #   geom_point()
      
      ## calculate sd of suitability around mean 
      sd_suitability = append(sd_suitability, sd(ts$suitability_ldetrended, na.rm = F))
      
    } 
    else {
      se_detrended = append(se_detrended, NA)
      sd_suitability = append(sd_suitability, NA)
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
  # hist(se_detrended)
  # hist(sd_suitability)
  # hist(suitability_trends$estimate)
  
  ## plot 
  results = xy
  results$se_detrended = se_detrended
  results$sd_suitability = sd_suitability
  results$suitability_trends = suitability_trends$estimate[which(suitability_trends$term == "(Intercept)")]
  
  results = dplyr::select(results, x, y, se_detrended, sd_suitability, suitability_trends)
  
  r_temp = rast(results)
  plot(-r_temp$se_detrended)
  plot(r_temp$sd_suitability)
  plot(r_temp$suitability_trends)

  ## save 
  results <- left_join(results, select(xy, x, y, buffer_edge, Longitude_bin, min, max))
  write.csv(results, paste0("outputs/data-processed/suitability-autocorrelation", str_replace_all(cur_sp, " ", "_"), "_autocorrelation.csv"), 
            row.names = F)
  
  ## go to next species 
  sp = sp + 1
  print(sp)
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
  