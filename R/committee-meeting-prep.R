## make plot of r during first 100 shifting time period
for(i in 500:600) {
  ras = rast(lattice_r_array[,,i])
  #plot(ras)
  
  gg = ras %>%
    ggplot() +
    geom_spatraster(data = ras) +
    theme_bw() +
    scale_fill_viridis_c(limits = c(-1.01, 2.06)) +
    labs(fill = "r") +
    theme_void() 
  
  ggsave(gg, path = "outputs/figures/didactic/r-over-time_variable", 
         filename = paste0("r-across-grid_", i, ".png"), width = 2, height = 9)
  
}


## make range shift gifs
library(magick)
library(magrittr)
library(tidyverse)

## order files by time
files = list.files("outputs/figures/didactic/r-over-time_variable", full.names = TRUE)

t = str_split_fixed(files, "\\_", n = 3)[,3]
t = str_split_fixed(t, ".png", n = 2)[,1]
files = files[order(as.numeric(t))]

## make a gif
frames = files %>% 
  image_read() %>% # reads each path file
  image_join() 

## make sure none failed
frames %>% # joins images
  image_apply(function(img) image_background(img, "white")) %>%
  image_animate(fps = 25) %>% # animates, can opt for number of loops
  image_write(paste0("outputs/figures/didactic/r-over-time_variable.gif")) # write to current dir




plot_shift_params_modified <- function(all_sims, 
                              p = c(),
                              beta = c(), 
                              d = c(), 
                              d_dist = c(), 
                              shift_rate = c(), 
                              icp = c(), 
                              K = c(), 
                              sigma = c(),
                              rep = c(), param = c("q5_y", "q95_y")) {
  # p = c()
  # beta = c()
  # d = c()
  # d_dist = c()
  # shift_rate = c()
  # icp = c()
  # K = c()
  # sigma = c()
  # rep = c()
  # param = c("q5_y", "q95_y", )
  
  if(length(p) == 0) {
    p_arg = unique(all_sims$p)
  } else {
    p_arg = p
  }
  if(length(beta) == 0) {
    beta_arg = unique(all_sims$beta)
  }else {
    beta_arg = beta
  }
  if(length(d) == 0) {
    d_arg = unique(all_sims$d)
  }else {
    d_arg = d
  }
  if(length(d_dist) == 0) {
    d_dist_arg = unique(all_sims$d_dist)
  }else {
    d_dist_arg = d_dist
  }
  if(length(shift_rate) == 0) {
    shift_rate_arg = unique(all_sims$shift_rate)
  }else {
    shift_rate_arg = shift_rate
  }
  if(length(icp) == 0) {
    icp_arg = unique(all_sims$icp)
  }else {
    icp_arg = icp
  }
  if(length(K) == 0) {
    K_arg = unique(all_sims$K)
  }else {
    K_arg = K
  }
  if(length(sigma) == 0) {
    sigma_arg = unique(all_sims$sigma)
  }else {
    sigma_arg = sigma
  }
  if(length(rep) == 0) {
    rep_arg = unique(all_sims$rep)
  }else {
    rep_arg = rep
  }
  
  ## plot range shift parameters for specified simulations
  params <- all_sims %>%
    filter(p %in% p_arg, beta %in% beta_arg, d %in% d_arg, shift_rate %in% shift_rate_arg,
           icp %in% icp_arg, d_dist %in% d_dist_arg, K %in% K_arg, rep %in% rep_arg, 
           sigma %in% sigma_arg) 
  
  ## get max time step
  n_ts = max(params$t)
  
  ## plot
  plot = params %>%
    arrange(beta, p, icp, K, d, d_dist, rep, t) %>%
    mutate(period = ifelse(t < 500, "stable", "shifting")) %>%
    mutate(colour = ifelse(parameter %in% c("max_y","q95_y","mean_max_y"),
                           "cold", "hot"),
           linetype = ifelse(parameter %in% c("max_y","min_y"),"Max/min occupied row", 
                                              ifelse(parameter %in% c("q95_y","q5_y"),
                                                     "95th/5th percentile occupied row", 
                                                     "Mean max/min occupied row"))) %>%
    filter(parameter %in% param) %>% 
    ggplot(aes(x = t, y = measurement, colour = colour, group = parameter)) +
    scale_colour_manual(values = c("blue", "red")) +
    scale_linetype_manual(values = c("dashed", "dotted", "solid")) +
    geom_line(aes(linetype = linetype)) + 
    theme_bw() +
    #geom_smooth(method = "lm", aes(group = period)) +
    scale_x_continuous(limits = c(0, 2000)) +
    scale_y_continuous(limits = c(0, 300)) +
    labs(x = "Time", y = "Latitude", linetype = "Measurement/ntype")  +
    theme(legend.position = "none")
  
  return(plot)
}


## make plot of range edge position over time 
for(i in 500:2000) {
  
  all_sims_sub = all_sims %>%
    filter(t <= i)
  
  p = plot_shift_params_modified(all_sims_sub, 
                    p = 0,
                    beta = 1,
                    d = 0.1,
                    d_dist = 4,
                    shift_rate = 0.2,
                    icp = 0.7,
                    K = 200,
                    sigma = 0.2,
                    rep = 1, param = c("max_y","min_y","mean_min_y",
                                        "q95_y","q5_y","mean_max_y"))
  
  ## save
  ggsave(p, path = "outputs/figures/didactic/edge-over-time", filename = paste0("t", i, ".png"),
         width = 4, height = 3)
  print(i)
  
}

## make gif 

## order files by time
files = list.files("outputs/figures/didactic/edge-over-time", full.names = TRUE)

t = str_split_fixed(files, "time/t", n = 2)[,2]
t = str_split_fixed(t, ".png", n = 2)[,1]
files = files[order(as.numeric(t))]

files = files[1:501]

## make a gif
frames = files %>% 
  image_read() %>% # reads each path file
  image_join() 

## make sure none failed
frames %>% # joins images
  image_apply(function(img) image_background(img, "white")) %>%
  image_animate(fps = 25) %>% # animates, can opt for number of loops
  image_write(paste0("outputs/figures/didactic/edge-over-time.gif")) # write to current dir

## make gif of range shift simulation
## order files by time
files = list.files("outputs/figures/p1_b0_icp0.7_K200_d0.1_r_d-dist7_sigma0.2_shift-rate0.2_rep1", full.names = TRUE)

t = str_split_fixed(files, "rep1/t", n = 2)[,2]
t = str_split_fixed(t, ".png", n = 2)[,1]
files = files[order(as.numeric(t))]

## make a gif
frames = files %>% 
  image_read() %>% # reads each path file
  image_join() 

## make sure none failed
frames %>% # joins images
  image_apply(function(img) image_background(img, "white")) %>%
  image_animate(fps = 25) %>% # animates, can opt for number of loops
  image_write(paste0("outputs/figures/didactic/white-sync.gif")) # write to current dir


## make gif of cardinal range edge measurements
files = list.files("/Users/nikkimoore/Documents/autosynchrony-simulations/outputs/figures/range-edge-shifts_after-filtering/species-with-data/Cardinalis_cardinalis", full.names = TRUE)

t = str_split_fixed(files, "cardinalis_", n = 2)[,2]
t = str_split_fixed(t, ".png", n = 2)[,1]
files = files[order(as.numeric(t))]

## make a gif
frames = files %>% 
  image_read() %>% # reads each path file
  image_join() 

## make sure none failed
frames %>% # joins images
  image_apply(function(img) image_background(img, "white")) %>%
  image_animate(fps = 2) %>% # animates, can opt for number of loops
  image_write(paste0("outputs/figures/didactic/cardinal.gif")) # write to current dir



## plot cardinal sdm projections 
library(terra)
library(tidyterra)

r = rast("/Volumes/NIKKI/sdms/Cardinalis_cardinalis/predictions/Cardinalis_cardinalis_predictions_1966-2024.tif")

labs = names(r)
yrs = str_split_fixed(labs, "\\_", 3)[,2]
months = str_split_fixed(labs, "\\_", 3)[,3]
months = ifelse(months == "4", "April", ifelse(months == "5", "May",
                                               ifelse(months == "6", "June",
                                                      ifelse(months == "7", "July", "NA"))))

for(i in 1:nlyr(r)) {
  p = ggplot() +
    geom_spatraster(data = r[[i]]) +
    scale_fill_continuous(na.value = "white", limits = c(0, 1)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    scale_y_continuous(limits = c(85, 25)) +
    scale_x_continuous(limits = c(-180, -40)) +
    labs(fill = "Climate\nsuitability", title = paste0(months[i], ", ", yrs[i]))
  
  ## save
  ggsave(p, path = "outputs/figures/didactic/suitability-over-time", filename = paste0("t", i, ".png"),
         width = 7, height = 4)
  
  print(i)
}



## make gif of cardinal range edge measurements
files = list.files("outputs/figures/didactic/suitability-over-time", full.names = TRUE)

t = str_split_fixed(files, "/t", n = 2)[,2]
t = str_split_fixed(t, ".png", n = 2)[,1]
files = files[order(as.numeric(t))]

## make a gif
frames = files %>% 
  image_read() %>% # reads each path file
  image_join() 

## make sure none failed
frames %>% # joins images
  image_apply(function(img) image_background(img, "white")) %>%
  image_animate(fps = 10) %>% # animates, can opt for number of loops
  image_write(paste0("outputs/figures/didactic/cardinal-suitability.gif")) # write to current dir


### add suitability edge lines 
temp = all_results %>% 
  filter(genus_sp == "Cardinalis_cardinalis") %>%
  mutate(Month = str_split_fixed(year, "\\.", 2)[,2]) %>%
  mutate(Month = ifelse(Month == "4", "April", ifelse(Month == "5", "May",
                                                       ifelse(Month == "6", "June",
                                                              ifelse(Month == "7", "July", "NA")))))

for(i in 1:nlyr(r)) {
  p = ggplot() +
    geom_spatraster(data = r[[i]]) +
    scale_fill_continuous(na.value = "white", limits = c(0, 1)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    scale_y_continuous(limits = c(85, 25)) +
    scale_x_continuous(limits = c(-180, -40)) +
    labs(fill = "Climate\nsuitability", title = paste0(months[i], ", ", yrs[i])) +
    geom_segment(data = temp, inherit.aes = F,
                 aes(x = max, y = 60, yend = 25), colour = "black") +
    geom_segment(data = temp, inherit.aes = F,
                 aes(x = min, y = 60, yend = 25), colour = "black") +
    theme(panel.grid = element_blank())+
    geom_segment(data = filter(temp, Year == yrs[i], Month == months[i]), inherit.aes = F,
                 aes(x = max, xend = min, y = p95), colour = "darkgrey") 
  
  
  ## save
  ggsave(p, path = "outputs/figures/didactic/suitability-over-time_edges", filename = paste0("t", i, ".png"),
         width = 7, height = 4)
  
  print(i)
}

## make gif of cardinal range edge measurements
files = list.files("outputs/figures/didactic/suitability-over-time_edges", full.names = TRUE)

t = str_split_fixed(files, "/t", n = 2)[,2]
t = str_split_fixed(t, ".png", n = 2)[,1]
files = files[order(as.numeric(t))]

## make a gif
frames = files %>% 
  image_read() %>% # reads each path file
  image_join() 

## make sure none failed
frames %>% # joins images
  image_apply(function(img) image_background(img, "white")) %>%
  image_animate(fps = 10) %>% # animates, can opt for number of loops
  image_write(paste0("outputs/figures/didactic/cardinal-suitability_edges.gif")) # write to current dir


### plot both together
bb_card = bbs_edges %>%
  filter(genus_sp == "Cardinalis cardinalis") %>%
  filter(range_edge == "p95")

for(i in 57:nlyr(r)) {
  p = ggplot() +
    geom_spatraster(data = r[[i]]) +
    scale_fill_continuous(na.value = "white", limits = c(0, 1)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    scale_y_continuous(limits = c(85, 25)) +
    scale_x_continuous(limits = c(-180, -40)) +
    labs(fill = "Climate\nsuitability", title = paste0(months[i], ", ", yrs[i])) +
    geom_segment(data = temp, inherit.aes = F,
                 aes(x = max, y = 60, yend = 25), colour = "black") +
    geom_segment(data = temp, inherit.aes = F,
                 aes(x = min, y = 60, yend = 25), colour = "black") +
    theme(panel.grid = element_blank())+
    geom_segment(data = filter(temp, Year == yrs[i], Month == months[i]), inherit.aes = F,
                 aes(x = max, xend = min, y = p95), colour = "darkgrey") + 
    geom_segment(data = filter(bb_card, Year == yrs[i]), inherit.aes = F,
                 aes(x = max, xend = min, y = latitude), colour = "red") +  
    geom_sf(size = 0.1, data = filter(bbs_filtered, TotalAbd != 0 & Year == yrs[i]), inherit.aes = FALSE, colour = "red", 
            alpha = 0.1) 
  
  ## save
  ggsave(p, path = "outputs/figures/didactic/everything-over-time", filename = paste0("t", i, ".png"),
         width = 7, height = 4)
  
  print(i)
}


## make gif of cardinal range edge measurements
files = list.files("outputs/figures/didactic/everything-over-time", full.names = TRUE)

t = str_split_fixed(files, "/t", n = 2)[,2]
t = str_split_fixed(t, ".png", n = 2)[,1]
files = files[order(as.numeric(t))]

## make a gif
frames = files %>% 
  image_read() %>% # reads each path file
  image_join() 

## make sure none failed
frames %>% # joins images
  image_apply(function(img) image_background(img, "white")) %>%
  image_animate(fps = 10) %>% # animates, can opt for number of loops
  image_write(paste0("outputs/figures/didactic/everything-over-time.gif")) # write to current dir






## plot the range edge shift and suitability shift in each longitude band for a given species:
bbs_edges %>%
  filter(genus_sp == "Cardinalis cardinalis") %>%
  filter(range_edge == "p95") %>%
  ggplot(aes(x = Year, y = latitude, colour = range_edge)) +
  geom_line(colour = "red") +
  geom_smooth(method = "lm", se = F, colour = "red") +
  labs(x = "Year", y = "Latitude", title =  "Cardinalis cardinalis") +
  geom_line(data = all_results %>% filter(genus_sp == "Cardinalis_cardinalis", year >= 1980),
            aes(x = year, y = p95), colour = "grey", inherit.aes = F) +
  geom_smooth(data = all_results %>% filter(genus_sp ==  "Cardinalis_cardinalis", year >= 1980),
              method = "lm", aes(x = year, y = p95), inherit.aes = F, colour = "grey", se = F) +
  facet_wrap(~min, nrow = 1)  + ## this orders from west to east 
  scale_y_continuous(limits = c(27, 83)) +
  theme(strip.background = element_blank(), # Removes the background box
        strip.text = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect()) 




