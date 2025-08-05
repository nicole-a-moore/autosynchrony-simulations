library(terra)


## plot mean growth rate at latitude L over time 
for(t in 1:2000) {
  if(t == 1) {
    r_df = data.frame(t = t, lat = 1:100, 
                      r = sapply(1:100, FUN = function(x) { mean(array[x,1:10,t]) }))
  }
  r_df = rbind(r_df, 
               data.frame(t = t, lat = 1:100, 
                          r = sapply(1:100, FUN = function(x) { mean(array[x,1:10,t]) })))
}

r_df %>%
  filter(t == 796) %>%
  ggplot(aes(x = lat, y = r)) +
  geom_line() +
  scale_y_continuous(limits = c(-2, 4))

# plot mean growth rate at latitude L over moving windows
for(t in 50:1950) {
  if(t == 50) {
    r_df = data.frame(t = t, lat = 1:100,
                      r = sapply(1:100, FUN = function(x) { mean(lattice_r_array[x,1:10,(t-49):(t+49)]) }))
  }
  r_df = rbind(r_df,
               data.frame(t = t, lat = 1:100,
                          r = sapply(1:100, FUN = function(x) { mean(lattice_r_array[x,1:10,(t-49):(t+49)]) })))
}


## apply threshold and calculate 5th and 95th percentiles
r_df <- r_df %>%
  group_by(t) %>%
  mutate(p75 = quantile(.$lat[which(r >= 0)], 0.05),
         p95 = quantile(.$lat[which(r >= 0)], 0.95),
         peak_r = mean(.$lat[which(r == max(r))])) %>%
  ungroup()

r_df %>%
  filter(t == 990) %>%
  ggplot(aes(x = lat, y = r)) +
  geom_line() +
  scale_y_continuous(limits = c(-2, 4)) +
  #geom_smooth() +
  geom_vline(aes(xintercept = p5), colour = "red", linetype = 2) +
  geom_vline(aes(xintercept = p95), colour = "red", linetype = 2) +
  geom_vline(aes(xintercept = peak_r), colour = "red")

#### need to find a way to estimate latitude of niche edge (i.e., lat of 0 growth rate) from the growth rate data 
## moving window mean, first and last where r > 0?


## try with an env grid
env_grid = rast("outputs/data-processed/env-grids/range-shift-grid1_p1_beta1_r1.2_K100_d0.1_icp0.1_L2000.tif")

array = as.array(env_grid)

# plot mean growth rate at latitude L over moving windows
for(t in 50:1950) {
  if(t == 50) {
    r_df = data.frame(t = t, lat = 1:100,
                      r = sapply(1:100, FUN = function(x) { mean(array[x,1:10,(t-49):(t+49)]) }))
  }
  else {
    r_df = rbind(r_df,
               data.frame(t = t, lat = 1:100,
                          r = sapply(1:100, FUN = function(x) { mean(array[x,1:10,(t-49):(t+49)]) })))
  }
  
}


## apply threshold and calculate 5th and 95th percentiles
r_df <- r_df %>%
  group_by(t) %>%
  mutate(p95 = quantile(.$lat[which(r >= 0)], 0.95),
         p5 = quantile(.$lat[which(r >= 0)], 0.05),
         peak_r = mean(.$lat[which(r == max(r))])) %>%
  ungroup()

r_df %>%
  filter(t == 1200) %>%
  ggplot(aes(x = lat, y = r)) +
  geom_line() +
  scale_y_continuous(limits = c(-2, 4)) +
  #geom_smooth() +
  geom_vline(aes(xintercept = p5), colour = "red", linetype = 2) +
  geom_vline(aes(xintercept = p95), colour = "red", linetype = 2) +
  geom_vline(aes(xintercept = peak_r), colour = "red") +
  geom_hline(aes(yintercept = 0), colour = "red")

## plot beside range edges from simulation 

r_df <- r_df %>%
  select(t, p5, p95, peak_r) %>%
  distinct()



## add niche parameters 
params %>%
  mutate(period = ifelse(t < 500, "stable", "shifting")) %>%
  gather(key = "parameter", value = "measurement", c("max_y", "min_y", "mean_y", "q95_y", "q5_y",
                                                     "mean_max_y", "mean_min_y", "abd_centroid", 
                                                     "ext_edge", "est_edge")) %>%
  filter(parameter %in% c("q5_y", "q95_y", "abd_centroid")) %>% 
  mutate(period = paste(period, parameter)) %>%
  ggplot(aes(x = t, y = measurement, colour = parameter)) +
  #geom_line(data = line, inherit.aes = F, aes(x = t, y = Nt)) +
  geom_line() + 
  theme_bw() +
  #geom_smooth(method = "lm", aes(group = period)) +
  scale_x_continuous(limits = c(0, 1500)) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Time", y = "Latitude", colour = "Range edge") +
  geom_line(data = r_df, inherit.aes = F, aes(x = t, y = 100 - peak_r)) +
  geom_line(data = r_df, inherit.aes = F, aes(x = t, y = 100 - p95)) +
  geom_line(data = r_df, inherit.aes = F, aes(x = t, y = 100 - p5)) 

## range edge over time across latitude: estimate yearly from bbs surveys
## change in r over time: project suitability using temperature (+ precipitation) at a daily scale and measure change in latitude of the 95th and 5th percentile of cells above x threshold in suitability 
## even if no directional climate change occurs (like in simulations), could still test hypotheses by looking at dynamics of edges vs. r?


## plot lag over time
temp = params %>%
  left_join(., r_df) 


temp2 = params %>%
  left_join(., r_df) 

temp3 = params %>%
  left_join(., r_df) 

temp$p = 1
temp2$p = 0
temp3$p = "0_0"

rbind(temp, temp2) %>%
  rbind(., temp3) %>%
  ggplot(aes(x = t, y = (100-peak_r) - q95_y, group = p,  colour = p)) + 
  geom_line() +
  geom_line(aes(y = (100-peak_r) - q5_y))


rbind(temp, temp2) %>%
  rbind(., temp3) %>%
  mutate(period = ifelse(t < 500, "stable", "shifting")) %>%
  gather(key = "parameter", value = "measurement", c("max_y", "min_y", "mean_y", "q95_y", "q5_y",
                                                     "mean_max_y", "mean_min_y", "abd_centroid", 
                                                     "ext_edge", "est_edge")) %>%
  filter(parameter %in% c("q5_y", "q95_y", "abd_centroid")) %>% 
  mutate(period = paste(period, parameter), 
         group = paste(parameter, p)) %>%
  ggplot(aes(x = t, y = measurement, colour = p, group = group)) +
  #geom_line(data = line, inherit.aes = F, aes(x = t, y = Nt)) +
  geom_line() + 
  theme_bw() +
  #geom_smooth(method = "lm", aes(group = period)) +
  scale_x_continuous(limits = c(0, 2000)) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Time", y = "Latitude", colour = "Range edge") 


