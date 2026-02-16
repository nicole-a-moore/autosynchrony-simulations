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
files = list.files("outputs/figures/range-shift-snapshots/p0_b0_icp0.7_K200_d0.1_r_d-dist7_sigma0.5_shift-rate0.2_rep1", full.names = TRUE)

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
  image_write(paste0("outputs/figures/didactic/white-async.gif")) # write to current dir


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
  image_animate(fps = 25) %>% # animates, can opt for number of loops
  image_write(paste0("outputs/figures/didactic/cardinal.gif")) # write to current dir




