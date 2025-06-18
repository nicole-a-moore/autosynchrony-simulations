######################################################
##               try making hist plots              ##
######################################################
## try for northern cardinal

## read in clean bbs data
bbs =  read.csv("outputs/data-processed/bbs_clean-subset.csv")

## filter to northern cardinal  
nc <- filter(bbs, AOU == 5930) 

## now, count extinctions, establishments and whether cell is new establishment or not
## split by routes
split = nc %>%
  group_by(route) %>%
  group_split(.)

list = lapply(split, FUN = function(df) {
  ts = df$TotalAbd
  
  ## get rid of NAs but add back later
  ts = ts[which(!is.na(ts))]
  
  zero_count = 0
  ext = 0
  est = 0
  max_zero_count <- 0
  zero_count_vec <- c()
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
    zero_count_vec <- append(zero_count_vec, zero_count)
  }
  
  if(any(ts[which(!is.na(ts))[1:2]] >= 1)) {
    new_est = FALSE
  }
  else {
    new_est = TRUE
  }
  
  ext = zero_count_vec 
  ext[which(ext != 1)] = 0
  ext[1] = 0
  
  est = zero_count_vec 
  shifted = zero_count_vec[2:length(ts)]
  est = ifelse(shifted - est[1:((length(ts))-1)] < 0, 1, 0) 
  est = c(0, est)
  
  ## save results in df 
  ext_est = df %>%
    filter(!is.na(TotalAbd)) %>%
    mutate(est = est, ext = ext)
  
  ## add back NA counts 
  year = df$Year[which(!is.na(df$TotalAbd))]
  ext_est = full_join(ext_est, df) %>% 
    arrange(Year) %>%
    mutate(n_ts = length(ts))
  
  return(ext_est)
})

## bind them all 
ext_est <- bind_rows(list)

## group by year, calculate latitude of centroid and add column for latitude relative to centroid
centroid = ext_est %>%
  filter(!is.na(TotalAbd)) %>%
  group_by(Year) %>%
  mutate(centroid = round(weighted.mean(Latitude, w = TotalAbd, na.rm = T), digits = 2)) %>% 
  ungroup() %>%
  select(Year, centroid) %>%
  distinct() 

## sum number of 'state changes' per latitudinal bin relative to range centroid 
ext_est = ext_est %>%
  left_join(., centroid) %>%
  mutate(rel_lat = Latitude - centroid) %>%
  mutate(rel_lat = round(rel_lat, digits = 0)) %>% ## bind rel lat
  filter(!is.na(TotalAbd)) %>% ## get rid of NA
  group_by(rel_lat) %>%
  mutate(n_samples_rel_lat = n()) %>% ## calculate number of times a bbs route in that latitudinal bin was sampled 
  mutate(state_changes = sum(ext, est, na.rm = T)) %>% 
  mutate(ext_sum = sum(ext, na.rm = T)) %>% 
  mutate(est_sum = sum(est, na.rm = T)) %>% 
  mutate(occ = length(which(TotalAbd != 0))) %>% 
  ungroup() %>%
  mutate(freq_state_changes = state_changes/n_samples_rel_lat,
         ext_rate = ext_sum/n_samples_rel_lat,
         est_rate = est_sum/n_samples_rel_lat,
         occ_freq = occ/n_samples_rel_lat)

hist(ext_est$rel_lat)

## frequency of state changes 
ext_est %>%
  select(rel_lat, freq_state_changes) %>% 
  distinct() %>% 
  ggplot(aes(x = rel_lat, y = freq_state_changes)) +
  geom_col(position = "identity", alpha = 0.7) +
  labs(x = "Latitude relative to centroid of abundance", y = "Frequency of state changes (changes/unit time)")
## increases towards leading edge (since cardinal expanded northern range edge)

## frequency of occupancy 
ext_est %>%
  select(rel_lat, occ_freq) %>% 
  distinct() %>% 
  ggplot(aes(x = rel_lat, y = occ_freq)) +
  geom_col(position = "identity", alpha = 0.7) +
  labs(x = "Latitude relative to centroid of abundance", y = "Frequency of occupancy")
## tapers towards leading edge but not trailing edge
## perhaps since trailing edge of cardinal in BBS is not really a real trailing edge 


## split state changes into ext/est
ext_est %>%
  select(rel_lat, ext_rate, est_rate) %>% 
  distinct() %>% 
  gather(key = "type", value = "rate", c(ext_rate, est_rate)) %>%
  ggplot(aes(x = rel_lat, y = rate, fill = type)) +
  geom_col(position = "identity", alpha = 0.5) +
  labs(x = "Latitude relative to centroid of abundance", y = "Frequency of state changes (changes/unit time)")

## split by time period (before 1980, after 1980)


## bind them all 
ext_est <- bind_rows(list)
ext_est$period = ifelse(ext_est$Year <= 1980, "stable", "shifting")

df_shift = ext_est %>% filter(period == "shifting")
df_stable = ext_est %>% filter(period == "stable")

## group by year, calculate latitude of centroid and add column for latitude relative to centroid
centroid = ext_est %>%
  filter(!is.na(TotalAbd)) %>%
  group_by(Year) %>%
  mutate(centroid = round(weighted.mean(Latitude, w = TotalAbd, na.rm = T), digits = 2)) %>% 
  ungroup() %>%
  select(Year, centroid) %>%
  distinct() 

## sum number of 'state changes' per latitudinal bin relative to range centroid 
ext_est = ext_est %>%
  left_join(., centroid) %>%
  mutate(rel_lat = Latitude - centroid) %>%
  mutate(rel_lat = round(rel_lat, digits = 0)) %>% ## bind rel lat
  filter(!is.na(TotalAbd)) %>% ## get rid of NA
  group_by(rel_lat, period) %>%
  mutate(n_samples_rel_lat = n()) %>% ## calculate number of times a bbs route in that latitudinal bin was sampled 
  mutate(state_changes = sum(ext, est, na.rm = T)) %>% 
  mutate(ext_sum = sum(ext, na.rm = T)) %>% 
  mutate(est_sum = sum(est, na.rm = T)) %>% 
  mutate(occ = length(which(TotalAbd != 0))) %>% 
  ungroup() %>%
  mutate(freq_state_changes = state_changes/n_samples_rel_lat,
         ext_rate = ext_sum/n_samples_rel_lat,
         est_rate = est_sum/n_samples_rel_lat,
         occ_freq = occ/n_samples_rel_lat)

ext_est %>%
  select(rel_lat, freq_state_changes, period) %>% 
  distinct() %>% 
  ggplot(aes(x = rel_lat, y = freq_state_changes, fill = period)) +
  geom_col(position = "identity", alpha = 0.5) +
  labs(x = "Latitude relative to centroid of abundance", y = "Frequency of state changes (changes/unit time)")

## frequency of occupancy 
ext_est %>%
  select(rel_lat, occ_freq, period) %>% 
  distinct() %>% 
  ggplot(aes(x = rel_lat, y = occ_freq, fill = period)) +
  geom_col(position = "identity", alpha = 0.7) +
  labs(x = "Latitude relative to centroid of abundance", y = "Frequency of occupancy")

ext_est %>%
  select(rel_lat, ext_rate, est_rate, period) %>% 
  distinct() %>% 
  gather(key = "type", value = "rate", c(ext_rate, est_rate)) %>%
  ggplot(aes(x = rel_lat, y = rate, fill = period)) +
  geom_col(position = "identity", alpha = 0.5) +
  facet_wrap(~type) +
  labs(x = "Latitude relative to centroid of abundance", y = "Frequency of state changes (changes/unit time)")


ext_est %>%
  select(rel_lat, ext_rate, est_rate, period) %>% 
  distinct() %>% 
  gather(key = "type", value = "rate", c(ext_rate, est_rate)) %>%
  ggplot(aes(x = rel_lat, y = rate, fill = type)) +
  geom_col(position = "identity", alpha = 0.5) +
  facet_wrap(~period) +
  labs(x = "Latitude relative to centroid of abundance", y = "Frequency of state changes (changes/unit time)")

