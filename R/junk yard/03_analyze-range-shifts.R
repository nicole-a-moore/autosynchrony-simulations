## analyze range shift 
## plot 
library(tidyverse)
select = dplyr::select


###########################################
##               functions               ##
###########################################
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
  
  ## then calculate lags for each by subtracting position of 95%, 5% and max suitability
  
  ## get ps and betas
  ps = c(0,1)
  betas = c(0,1)
  rep = unique(all_ranges$rep)
  
  for(p in ps) {
    for(beta in betas) {
      
      all_ranges_sub <- all_ranges[which(all_ranges$p == p & all_ranges$beta == beta),]
        
      params = all_ranges_sub %>%
        filter(Nt != 0) %>%
        group_by(t) %>%
        summarize(max_y = max(y), 
                  q95_y = quantile(y, c(0.95)),
                  min_y = min(y),
                  q5_y = quantile(y, c(0.05)),
                  abd_centroid = weighted.mean(y, Nt),
                  mean_y = mean(y)) 
      
      params2 = all_ranges_sub %>%
        filter(Nt != 0) %>% 
        group_by(t) %>%
        arrange(-y) %>%
        slice(1:10) %>%
        summarize(mean_max_y = mean(y))
      
      params3 = all_ranges_sub %>%
        filter(Nt != 0) %>% 
        group_by(t) %>%
        arrange(y) %>%
        slice(1:10) %>%
        summarize(mean_min_y = mean(y))
      
      mid_occ_lat = all_ranges_sub[which(all_ranges_sub$Nt != 0),]
      mid_occ_lat = mean(mid_occ_lat$y, na.rm = TRUE)
      
      params4 = all_ranges_sub %>%
        filter(est == 1) %>%
        group_by(t) %>%
        filter(y <= mid_occ_lat) %>%
        summarize(est_edge = quantile(y, c(0.05), na.rm = T))
      params5 = all_ranges_sub %>%
        filter(ext == 1) %>%
        group_by(t) %>%
        filter(y > mid_occ_lat) %>%
        summarize(ext_edge = quantile(y, c(0.95), na.rm = T))
      
      params <- left_join(params, params2) %>% left_join(., params3) %>% left_join(., params4) %>%
        left_join(., params5)
      
      line = data.frame(t = c(1:1999), Nt = c(rep(75, 499), c(-max(all_ranges_sub$shift_rate)*1:1500 + 75)))
      
      plot = params %>%
        mutate(period = ifelse(t < 500, "stable", "shifting")) %>%
        gather(key = "parameter", value = "measurement", c("max_y", "min_y", "mean_y", "q95_y", "q5_y",
                                                           "mean_max_y", "mean_min_y", "abd_centroid", 
                                                           "ext_edge", "est_edge")) %>%
        filter(parameter %in% c("q5_y", "q95_y", "abd_centroid")) %>% 
        mutate(period = paste(period, parameter)) %>%
        ggplot(aes(x = t, y = measurement, colour = parameter)) +
        geom_line() + 
        theme_bw() +
        #geom_smooth(method = "lm", aes(group = period)) +
        geom_line(data = line, inherit.aes = F, aes(x = t, y = Nt)) +
        scale_x_continuous(limits = c(0, 1500)) +
        scale_y_continuous(limits = c(0, 100)) +
        labs(x = "Time", y = "Latitude", colour = "Range edge")
      
      plot2 = params %>%
        mutate(period = ifelse(t < 500, "stable", "shifting")) %>%
        gather(key = "parameter", value = "measurement", c("max_y", "min_y", "mean_y", "q95_y", "q5_y",
                                                           "mean_max_y", "mean_min_y", "abd_centroid", 
                                                           "ext_edge", "est_edge")) %>%
        filter(parameter %in% c("ext_edge", "est_edge")) %>% 
        mutate(period = paste(period, parameter)) %>%
        ggplot(aes(x = t, y = measurement, colour = parameter)) +
        geom_line() + 
        theme_bw() +
        #geom_smooth(method = "lm", aes(group = period)) +
        geom_line(data = line, inherit.aes = F, aes(x = t, y = Nt)) +
        scale_x_continuous(limits = c(0, 1500)) +
        scale_y_continuous(limits = c(0, 100)) +
        labs(x = "Time", y = "Latitude", colour = "Range edge")
      
      ggsave(plot, path = "outputs/figures/range-shift-simulations/", 
             filename = paste0("rep", rep, "_p", 
                               p, "_b", 
                               beta, "_icp0.1", 
                               "_shift", ".png"), 
             width = 5, height = 3)
      # ggsave(plot2, path = "outputs/figures/range-shift-simulations/", 
      #        filename = paste0("rep", rep, "_p", 
      #                          p, "_b", 
      #                          beta, "_icp0.1", 
      #                          "_shift_est-ext-edge", ".png"), 
      #        width = 5, height = 3)
      
    }
  }

}

calculate_ext_est <- function(all_ranges, threshold) {
  ## for each cell, calculate no. ext and est. over time 
  ## split by cell
  split = all_ranges %>%
    group_by(x, y, p, beta) %>%
    group_split(.)
  
  list = lapply(split, FUN = function(df) {
    ts = df$Nt
    zero_count = 0
    ext = 0
    est = 0
    max_zero_count <- 0
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
    }
    
    if(any(ts[which(!is.na(ts))[1:2]] >= 1)) {
      new_est = FALSE
    }
    else {
      new_est = TRUE
    }
    
    ext_est = df %>%
      mutate(n_est = est, n_ext = ext, new_est = new_est, n_ts = length(ts), p = df$p, beta = df$beta)
    
    return(ext_est)
  })
  
  ## bind them all 
  df <- bind_rows(list)
  
  ## problem: no. of local extinctions/establishments is lower when range goes completely extinct before 1500 time steps
  ## solution: calculate rate of extinction/establishment by dividing it by maximum
  
  df <- df %>%
    group_by(beta, p, t) %>%
    filter(!all(Nt == 0)) %>%
    group_by(beta, p) %>%
    mutate(t_max = max(t),
           mean_rate_ext = mean(n_ext/t_max),
           mean_rate_est = mean(n_est/t_max)) %>%
    ungroup() %>%
    select(x, y, n_ext, n_est, p, beta, t_max, mean_rate_ext, mean_rate_est) %>%
    distinct()
  
  ## plot histogram
  # df %>%
  #   mutate(p = paste0("p", p), beta = paste0("beta", beta)) %>%
  #   ggplot(aes(x = n_ext/t_max)) +
  #   geom_histogram() +
  #   facet_wrap(p~beta) +
  #   labs(x = "Local extinction rate", y = "Frequency") +
  #   geom_vline(aes(xintercept = mean_rate_ext), colour = "red")
  # 
  # df %>%
  #   mutate(p = paste0("p", p), beta = paste0("beta", beta)) %>%
  #   ggplot(aes(x = n_est/t_max)) +
  #   geom_histogram() +
  #   facet_wrap(p~beta) +
  #   labs(x = "Local establishment rate", y = "Frequency") +
  #   geom_vline(aes(xintercept = mean_rate_est), colour = "red")
  
  return(df)
}

est_ext_grid <- function(all_ranges, threshold) {
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
    shifted = zero_count_vec[2:nrow(df)]
    est = ifelse(shifted - est[1:(nrow(df) - 1)] < 0, 1, 0) 
    est = c(0, est)
    
    return(data.frame(x = unique(df$x), y = unique(df$y), 
                      p = unique(df$p), beta = unique(df$beta), 
                      t = 1:nrow(df),
                      est = est, ext = ext))
    
  })
  ## bind them all 
  df <- bind_rows(list)
  
  ## left join to other stats
  all_ranges = left_join(all_ranges, df)
  
  return(all_ranges)
}



## try measuring mean latitude of cells in which the species was gained 
## for each cell, calculate no. ext and est. over time 
## split by cell







## problem: no. of local extinctions/establishments is lower when range goes completely extinct before 1500 time steps
## solution: calculate rate of extinction/establishment by dividing it by maximum

df <- df %>%
  group_by(beta, p, t) %>%
  filter(!all(Nt == 0)) %>%
  group_by(beta, p) %>%
  mutate(t_max = max(t),
         mean_rate_ext = mean(n_ext/t_max),
         mean_rate_est = mean(n_est/t_max)) %>%
  ungroup() %>%
  select(x, y, n_ext, n_est, p, beta, t_max, mean_rate_ext, mean_rate_est) %>%
  distinct()

## plot histogram
# df %>%
#   mutate(p = paste0("p", p), beta = paste0("beta", beta)) %>%
#   ggplot(aes(x = n_ext/t_max)) +
#   geom_histogram() +
#   facet_wrap(p~beta) +
#   labs(x = "Local extinction rate", y = "Frequency") +
#   geom_vline(aes(xintercept = mean_rate_ext), colour = "red")
# 
# df %>%
#   mutate(p = paste0("p", p), beta = paste0("beta", beta)) %>%
#   ggplot(aes(x = n_est/t_max)) +
#   geom_histogram() +
#   facet_wrap(p~beta) +
#   labs(x = "Local establishment rate", y = "Frequency") +
#   geom_vline(aes(xintercept = mean_rate_est), colour = "red")







source("R/functions/plot_range.R")

## set up filenames
files = list.files("outputs/data-processed/range-shifts/all-ranges", full.names = T)

files_stable <- list.files("outputs/data-processed/stable-ranges/all-ranges", full.names = T)

## get reps
reps = str_split_fixed(files, "outputs/data-processed/range-shifts/all-ranges/range-shifts_rep", 2)[,2]
reps = unique(as.numeric(str_split_fixed(reps, "_p", 2)[,1]))

## for each rep
r = 1
while(r <= length(unique(reps))) {
  
  files_sub = files[which(str_detect(files, paste0("rep", reps[r], "_")))]
  files_stable_sub = files_stable[which(str_detect(files_stable, paste0("rep", reps[r], "_")))]
  #files_stable_sub = files_stable_sub[which(str_detect(files_stable_sub,"d0.08"))]
    
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

  ## add est and extinction grid
  all_ranges = est_ext_grid(all_ranges, threshold = 1)
  
  ## calculate range edge parameters and make plots of edges over time
  calculate_shift_params(all_ranges = all_ranges)
  
  ## plot ranges in each time step 
  #plot_range(range_shift = all_ranges, path = "outputs/figures/range-shifts")
  
  ## calculate extinction and establishment rates
  # df <- calculate_ext_est(all_ranges, threshold = 5)
  # 
  # df %>%
  #   ggplot(aes(x = p, y = n_ext/t_max, colour = beta)) +
  #   geom_point(position = position_jitterdodge()) +
  #   facet_wrap(~beta) +
  #   labs(y = "Local extinction rate", x = "Synchrony") +
  #   geom_point(aes(y = mean_rate_ext), colour = "red")
  # 
  # df %>%
  #   ggplot(aes(x = p, y = n_est/t_max, colour = beta)) +
  #   geom_point(position = position_jitterdodge()) +
  #   facet_wrap(~beta) +
  #   labs(y = "Local establishment rate", x = "Synchrony") +
  #   geom_point(aes(y = mean_rate_est), colour = "red")
  
  r = r + 1
}



## MEASURE AND PLOT LAG OVER TIME
## - how? how to measure niche edge from simulations even?

## FIX SCRIPT SO SHIFTING AND NON SHIFTING PERIODS HAPPEN TOGETHER 

## TRY WITH INTERMEDIATE LEVELS OF SYNCHRONY (p=0.5)

## MAKE RANGE WIDER SO THAT EXTINCTION DOESN'T HAPPEN AS QUICKLY WHEN p=1 beta=1










df %>%
  ggplot(aes(x = new_est)) +
  geom_bar() +
  theme(legend.position = "none") +
  labs(x = "No. of routes", y = "New estabishment?")

## plot range edge (5th and 95th percentile of occupied cells) vs. centre points 

## calculate row of 5th and 95th percentile of occupied cells during each time step
df <- df %>%
  filter(Nt != 0) %>%
  group_by(t) %>%
  mutate(p95 = quantile(y, c(0.95)),
         p5 = quantile(y, c(0.05))) %>%
  ungroup() %>%
  mutate(range_position = ifelse(y >= p95, "leading-edge",
                                 ifelse(y <= p5, "trailing-edge", "centre")))

plot(t[[1]])

## calculate extinction and establishment within those cells 
df %>%
  gather(key = "measure", value = "count", c(n_ext, n_est)) %>%
  mutate(measure = ifelse(measure == "n_ext", "Extinctions", "Establishments")) %>%
  select(x, y, measure, count, range_position) %>%
  distinct() %>%
  ggplot(aes(x = count, fill = range_position)) +
  geom_histogram(position = "dodge") +
  facet_wrap(~measure) + 
  labs(x = "No. of local events", y = "Frequency")





#### garbage 
## calculate local ext and est for each grid cell
split = all_ranges %>%
  group_by(x, y) %>%
  group_split(.)

list = lapply(split, FUN = function(df) {
  ext_est <- all_ranges %>%
    filter(x == unique(df$x), y == unique(df$y)) %>%
    arrange(t) %>%  
    group_by(x, y, temp = cumsum(c(0,diff(Nt)) != 0)) %>%
    mutate(cum.zeros = ifelse(Nt, NA_integer_, row_number())) %>%
    ungroup() %>% 
    mutate(estab = ifelse(is.na(cum.zeros) & lag(cum.zeros, 1) >= 2, 1, NA)) %>%
    select(-temp) %>%
    mutate(n_ext = ifelse(nrow(.) == length(which(Nt == 0)), 0, length(which(cum.zeros == 2))),
           n_est = length(which(estab == 1))) %>%
    mutate(new_est = ifelse(first(Nt) == 0 & any(Nt) != 0, 1, NA))
  
  return(ext_est)
})

## bind them all 
list <- bind_rows(list)

## plot histogram 
list %>%
  select(x, y, n_ext, n_est) %>%
  distinct() %>% 
  ggplot(aes(x = n_est)) +
  geom_histogram()

ggsave(path = "outputs/figures/range-shift-simulations", 
       filename = paste0("rep", unique(all_ranges$rep), "_p", unique(all_ranges$p), "_b", unique(all_ranges$beta), "_icp", 
                         unique(all_ranges$icp), "_n-ext", ".png"), width = 4, height = 3)

all_ranges %>%
  filter(p==1, beta==1, t == 1) %>% View



