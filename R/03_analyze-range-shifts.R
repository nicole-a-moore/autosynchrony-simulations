## analyze range shift 
## plot 
library(tidyverse)
select = dplyr::select

all_ranges = read.csv(paste0("outputs/data-processed/range-shifts/", 
                             "rep", 1, "_p", 0, "_b", 1, "_icp", 0.1, 
                             "_range-shifts.csv"))

all_ranges = mutate(all_ranges, cell = paste0(x, "_", y)) 


calculate_shift_params <- function(all_ranges) {
  
  ## for each time step, calculate: 
  ## LEADING EDGE
  ## - y position of furthest occupied cell 
  ## - mean y position of furthest 10 occupied cells 
  ## - 95% of occupancy
  ## TRAILING EDGE
  ## - y position of furthest occupied cell  
  ## - mean y position of furthest 10 occupied cells 
  ## - 5% of occupancy
  ## CENTROID
  ## - mean y value occupancy
  ## - abuncance weighted mean y value of occupancy
  
  params = all_ranges %>%
    filter(Nt != 0) %>%
    group_by(t) %>%
    summarize(max_y = max(y), 
              q95_y = quantile(y, c(0.95)),
              min_y = min(y),
              q5_y = quantile(y, c(0.05)),
              abd_centroid = weighted.mean(.$y, weights = .$Nt),
              mean_y = mean(y)) 
  
  params2 = all_ranges %>%
    filter(Nt != 0) %>% 
    group_by(t) %>%
    arrange(-y) %>%
    slice(1:10) %>%
    summarize(mean_max_y = mean(y))
  
  params3 = all_ranges %>%
    filter(Nt != 0) %>% 
    group_by(t) %>%
    arrange(y) %>%
    slice(1:10) %>%
    summarize(mean_min_y = mean(y))
  
  params <- left_join(params, params2) %>% left_join(., params3)
  
  params %>%
    gather(key = "parameter", value = "measurement", c("max_y", "min_y", "mean_y", "q95_y", "q5_y",
                                                       "mean_max_y", "mean_min_y", "abd_centroid")) %>%
    filter(!parameter %in% c("max_y", "min_y", "q5_y", "q95_y")) %>% 
    ggplot(aes(x = t, y = measurement, colour = parameter)) +
    geom_line() + 
    theme_bw() +
    geom_smooth(method = "lm")
  
  ggsave(plot, path = "outputs/figures/range-shift-simulations/", 
         filename = paste0("/rep", rep, "_p", p, "_b", beta, "_icp", icp, "_snapshot", t, ".png"), 
         width = 2, height = 4)
  
  
}



all_ranges %>%
  select(t, N_global) %>%
  distinct() %>%
  ggplot(aes(x = t, y = N_global)) + geom_line()

all_ranges %>%
  dplyr::select(t, N_ext_local) %>%
  distinct() %>%
  ggplot(aes(x = t, y = N_ext_local)) + geom_line()

all_ranges %>%
  select(-x, -y) %>%
  filter(cell == unique(.$cell)[]) %>%
  ggplot(aes(x = t, y = Nt, colour = cell)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none")




