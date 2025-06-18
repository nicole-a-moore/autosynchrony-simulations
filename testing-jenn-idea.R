## =test Jenn's idea


## set up filenames
files = list.files("outputs/data-processed/range-shifts", full.names = T)
files <- files[which(str_detect(files, "range-shifts.csv"))]

files_stable <- list.files("outputs/data-processed/stable-ranges", full.names = T)
files_stable <- files_stable[which(str_detect(files_stable, "/rep"))]

rep = 10

## get reps
files_sub = files[which(str_detect(files, paste0("rep", rep, "_")))]
files_stable_sub = files_stable[which(str_detect(files_stable, paste0("rep", rep, "_")))]

files_sub[4] = "outputs/data-processed/range-shifts/rep4_p1_b1_icp0.1_range-shifts.csv"

## make one big file of all simulation results for this rep
all_ranges <- c()
for(i in 1:length(files_sub)) {
  all_ranges = rbind(all_ranges, read.csv(files_sub[i]))
}

## add stable ranges 
stable_ranges <- c()
for(i in 1:length(files_stable_sub)) {
  stable_ranges = rbind(stable_ranges, read.csv(files_stable_sub[i]))
}
stable_ranges$period = "stable"

all_ranges$period = "shifting"

## fix time variable 
all_ranges$t = all_ranges$t + 499

## combine
all_ranges <- rbind(stable_ranges, all_ranges)


threshold = 1

split = all_ranges %>%
  group_by(x, y, p, beta) %>%
  group_split(.)

list = lapply(split, FUN = function(df) {
  ts = df$Nt
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
        if(max_zero_count >= threshold) {
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
        if(max_zero_count >= threshold) {
          ext = ext + 1
        }
      }
    }
    else if(ts[i] >= 1 & max_zero_count >= threshold) {
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
    else if(ts[i] >= 1 & max_zero_count <= threshold) {
      zero_count = 0
      max_zero_count = 0
    }
    zero_count_vec <- append(zero_count_vec, zero_count)
  }
  
  ext = zero_count_vec 
  ext[which(ext != 1)] = 0
  ext[1] = 0
  
  est = zero_count_vec 
  shifted = zero_count_vec[2:length(unique(all_ranges$t))]
  est = ifelse(shifted - est[1:(length(unique(all_ranges$t))-1)] < 0, 1, 0) 
  est = c(0, est)
  
  return(data.frame(x = unique(df$x), y = unique(df$y), 
                    p = unique(df$p), beta = unique(df$beta), 
                    t = 1:nrow(df),
                    est = est, ext = ext, 
                    Nt = ts, period = df$period))
}
)

## bind them all 
df <- bind_rows(list)

df <- filter(df, t > 400) ## get rid of times when stable range wasn't yet stable

df_shift = df %>% filter(period == "shifting")
df_stable = df %>% filter(period == "stable")

tmax_shift =  df_shift %>%
  filter(Nt != 0) %>%
  group_by(beta, p) %>%
  mutate(t_max = max(t) - min(t) + 1) %>%
  select(t_max, beta, p, period, t) %>%
  distinct()

tmax_stable =  df_stable %>%
  filter(Nt != 0) %>%
  group_by(beta, p) %>%
  mutate(t_max = max(t) - min(t) + 1) %>%
  select(t_max, beta, p, period, t) %>%
  distinct()

tmax = rbind(tmax_shift, tmax_stable)

## group by time, calculate emax, and add column for latitude relative to centroid
centroid = df %>%
  group_by(t, y, beta, p) %>%
  summarize(Nt = sum(Nt)) %>%
  group_by(t, beta, p) %>%
  mutate(centroid = round(weighted.mean(y, w = Nt), digits = 0)) %>% 
  ungroup() %>%
  select(t, centroid, beta, p) %>%
  distinct() %>%
  left_join(., tmax)

## sum number of 'state changes' per latitudinal bin relative to range centroid 
df = df %>%
  left_join(., centroid) %>%
  mutate(rel_lat = centroid - y) %>%
  filter(t < 750) %>%
  group_by(rel_lat, p, beta, period) %>%
  mutate(state_changes = sum(ext, est, na.rm = T)) %>% 
  mutate(ext_sum = sum(ext, na.rm = T)) %>% 
  mutate(est_sum = sum(est, na.rm = T)) %>% 
  mutate(occ = length(which(Nt != 0))) %>% 
  ungroup() %>%
  mutate(t_max = ifelse(t_max < 250, t_max, 250)) %>%
  mutate(freq_state_changes = state_changes/t_max,
         ext_rate = ext_sum/t_max,
         est_rate = est_sum/t_max,
         occ_freq = occ/t_max)

df %>%
  select(rel_lat, freq_state_changes, beta, p, period) %>% 
  distinct() %>% 
  mutate(p = ifelse(p == 1, "Synchrony", "Asynchrony"),
         beta = ifelse(beta == 1, 'Autocorrelated noise', "White noise")) %>% 
  ggplot(aes(x = rel_lat, y = freq_state_changes, fill = period)) +
  geom_col(position = "identity", alpha = 0.7) +
  facet_grid(beta~p) +
  labs(x = "Latitude relative to centroid of abundance", y = "Frequency of state changes (changes/unit time)")

### MUST MAKE SURE EMAX NOT OUT OF GRID!!

## frequency of occupancy 
df %>%
  select(rel_lat, occ_freq, beta, p, period) %>% 
  distinct() %>% 
  mutate(p = ifelse(p == 1, "Synchrony", "Asynchrony"),
         beta = ifelse(beta == 1, 'Autocorrelated noise', "White noise")) %>% 
  ggplot(aes(x = rel_lat, y = occ_freq, fill = period)) +
  geom_col(position = "identity", alpha = 0.7) +
  facet_grid(beta~p) +
  labs(x = "Latitude relative to centroid of abundance", y = "Frequency of occupancy")






