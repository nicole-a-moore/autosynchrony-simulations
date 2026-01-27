## analyze range shift simulation data
library(tidyverse)
select = dplyr::select


###########################################
##               functions               ##
###########################################
calculate_shift_params <- function(simulation) {
  
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
  
  params = simulation %>%
    filter(Nt != 0) %>%
    group_by(t) %>%
    summarize(max_y = max(y), 
              q95_y = quantile(y, c(0.95)),
              min_y = min(y),
              q5_y = quantile(y, c(0.05)),
              abd_centroid = weighted.mean(y, Nt),
              mean_y = mean(y)) 
  
  params2 = simulation %>%
    filter(Nt != 0) %>% 
    group_by(t) %>%
    arrange(-y) %>%
    slice(1:10) %>%
    summarize(mean_max_y = mean(y))
  
  params3 = simulation %>%
    filter(Nt != 0) %>% 
    group_by(t) %>%
    arrange(y) %>%
    slice(1:10) %>%
    summarize(mean_min_y = mean(y))
  
  mid_occ_lat = simulation[which(simulation$Nt != 0),]
  mid_occ_lat = mean(mid_occ_lat$y, na.rm = TRUE)
  
  params <- left_join(params, params2, by = "t") %>% left_join(., params3, by = "t") 
  
  ## add back info about simulation:
  info = simulation %>%
    select(p, beta, shift_rate) %>%
    unique() 
  
  params$p = info$p
  params$beta = info$beta
  params$shift_rate = info$shift_rate
  
  ## gather different parameters 
  params <- params %>%
    gather(key = "parameter", value = "measurement", c("max_y", "min_y", "mean_y", "q95_y", "q5_y",
                                                     "mean_max_y", "mean_min_y", "abd_centroid")) 

  ## if there is data past the shifting period
  if(max(params$t) > 500) {
    
    rates <- params %>%
      mutate(period = ifelse(t < 500, "stable", "shifting")) %>%
      filter(period == "shifting") %>% 
      group_by(parameter) %>%
      do(broom::tidy(lm(measurement ~ t, data = .), conf.int = TRUE)) %>%
      filter(term == "t") %>%
      rename(shift_rate_estimate = estimate) %>%
      ungroup() %>%
      select(-term)
    
    params = left_join(params, rates, by = c("parameter"))
    
    
    return(params)
  }
  else {
    params$shift_rate_estimate = params$std.error = params$statistic = params$p.value = NA
    params$conf.low = params$conf.high = NA  
  }
  
  return(params)
}

plot_shift_params <- function(all_sims, 
                              p = c(),
                              beta = c(), 
                              d = c(), 
                              d_dist = c(), 
                              shift_rate = c(), 
                              icp = c(), 
                              K = c(), 
                              r_max = c(), 
                              rep = c(), param = "q95_y") {
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
  if(length(r_max) == 0) {
    r_max_arg = unique(all_sims$r_max)
  }else {
    r_max_arg = r_max
  }
  if(length(rep) == 0) {
    rep_arg = unique(all_sims$rep)
  }else {
    rep_arg = rep
  }
  
  ## plot range shift parameters for specified simulations
  params <- all_sims %>%
    filter(p %in% p_arg, beta %in% beta_arg, d %in% d_arg, shift_rate %in% shift_rate_arg,
           icp %in% icp_arg, d_dist %in% d_dist_arg, K %in% K_arg, r_max %in% r_max_arg, rep %in% rep_arg) 
  
  ## get max time step
  n_ts = max(params$t)
  
  ## plot
  if(n_ts > 500) {
    line = data.frame(t = 1:n_ts, Nt = c(rep(41, length.out = 500), shift_rate*(1:(n_ts-500)) + 41))
    
    plot = params %>%
      arrange(beta, p, icp, K, r_max, d, d_dist, rep, t) %>%
      mutate(period = ifelse(t < 500, "stable", "shifting")) %>%
      mutate(beta = as.character(beta)) %>%
      filter(parameter == param) %>% 
      mutate(period = paste(period)) %>%
      mutate(group = paste(beta, p, icp, K, r_max, d, d_dist, rep),
             p_beta = paste0("p = ", p, ", beta = ", beta)) %>% 
      ggplot(aes(x = t, y = measurement, colour = beta, group = group)) +
      geom_line(aes(group = group)) + 
      geom_line(data = line, inherit.aes = F, aes(x = t, y = Nt, group = shift_rate)) +
      theme_bw() +
      #geom_smooth(method = "lm", aes(group = period)) +
      scale_x_continuous(limits = c(0, 2000)) +
      scale_y_continuous(limits = c(0, 300)) +
      labs(x = "Time", y = "Latitude", colour = "p_beta") +
      facet_wrap(p~icp) 
  } else {
    plot = params %>%
      arrange(beta, p, icp, K, r_max, d, d_dist, rep, t) %>%
      mutate(period = ifelse(t < 500, "stable", "shifting")) %>%
      mutate(beta = as.character(beta)) %>%
      filter(parameter == param) %>% 
      mutate(period = paste(period)) %>%
      mutate(group = paste(beta, p, icp, K, r_max, d, d_dist, rep),
             p_beta = paste0("p = ", p, ", beta = ", beta)) %>% 
      ggplot(aes(x = t, y = measurement, colour = beta, group = group)) +
      geom_line(aes(group = group)) + 
      theme_bw() +
      #geom_smooth(method = "lm", aes(group = period)) +
      scale_x_continuous(limits = c(0, 2000)) +
      scale_y_continuous(limits = c(0, 300)) +
      labs(x = "Time", y = "Latitude", colour = "p_beta") +
      facet_wrap(p~icp) 
  }
  
  return(plot)
}

source("R/functions/plot_range.R")

############################################################
##               process simulation results               ##
############################################################
## process results from each simulation separately, measuring range edge parameters over time
## set up filenames
dir = "outputs/data-processed/range-shift-simulations/sim-results"
files = list.files(dir, full.names = T)
folders = list.files(dir)

## for each file
f = 1
while(f <= length(files)) {
  
  print(paste0("On simulation no. ", f))
  
  ## read in file
  cur_file = read.csv(files[f])
  
  ## calculate position of range edge in each time step
  ## estimate rate at which edge shifts with linear regression
  shift_data = calculate_shift_params(simulation = cur_file)
  
  ## get rep, dipsersal distance, dispersal proportion, icp, K, rmax
  ## and add info to shift data
  split = str_split_fixed(folders[f], "_", 10)
  
  shift_data$icp = substr(split[,4], 4,6)
  shift_data$K = substr(split[,5], 2,5)
  shift_data$d = substr(split[,6], 2,4)
  shift_data$r_max = substr(split[,7], 2,3)
  shift_data$d_dist = substr(split[,8], 7,7)
  shift_data$rep = unique(as.numeric(str_split_fixed(substr(split[,10], 4, nchar(split[,10])), ".csv", 2)[,1]))
  
  ## save shift data 
  write.csv(shift_data, paste0("outputs/data-processed/range-shift-simulations/cluster/shift-data/shift-data_",
                               folders[f]), row.names = F)

  f = f + 1
}


#########################################################
##               plot simulation results               ##
#########################################################
#### PLOT SIMULATION RESULTS TOGETHER 

## read in and combine all simulation shift data
files = list.files("outputs/data-processed/range-shift-simulations/cluster/shift-data", full.names = T)

all_sims <- c()
for(i in 1:length(files)) {
  print(paste0("On file no. ", i))
  all_sims = rbind(all_sims, read.csv(files[i]))
}


#### PLOT RANGES OVER TIME FOR SOME REPS x PARAM COMBINATIONS
plot_shift_params(all_sims, beta = c(0,1), p = c(0,1))








## NEXT: 
## add back trailing edge?
## measure and plot r edge over time 











#### garbage 
#############################################
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



split = all_ranges %>%
  filter(t %in% c(500, 1750), Nt != 0) %>%
  group_by(p, beta, t) %>%
  summarise(lat_pos = max(y)) %>% 
  group_by(t) %>%
  group_split()

df1 = split[[1]]
df2 = split[[2]]

df1 <- rename(df1, "lat_pos_500" = lat_pos) %>%
  select(-t)
df2 <- rename(df2, "lat_pos_1750" = lat_pos) %>%
  select(-t)

df = left_join(df1, df2)

df %>%
  mutate(surv_resurv_rate = (lat_pos_1750 - lat_pos_500)/(1750-500),
         p = as.character(p)) %>%
  ggplot(aes(y = surv_resurv_rate, x = beta, colour = p)) +
  geom_point()





# calculate_ext_est <- function(all_ranges, threshold) {
#   ## for each cell, calculate no. ext and est. over time 
#   ## split by cell
#   split = all_ranges %>%
#     group_by(x, y, p, beta) %>%
#     group_split(.)
#   
#   list = lapply(split, FUN = function(df) {
#     ts = df$Nt
#     zero_count = 0
#     ext = 0
#     est = 0
#     max_zero_count <- 0
#     for(i in 1:length(ts)) {
#       if(i == 1) {
#         if(ts[i] == 0) {
#           zero_count = 1
#         }
#       }
#       else if(is.na(ts[i]))  {
#         if(zero_count > max_zero_count) {
#           max_zero_count = zero_count
#         }
#         zero_count = 0
#         if(i == length(ts)) {
#           if(max_zero_count >= threshold) {
#             ext = ext + 1
#           }
#         }
#       }
#       else if(ts[i] == 0) {
#         zero_count = zero_count + 1
#         if(zero_count > max_zero_count) {
#           max_zero_count = zero_count
#         }
#         if(i == length(ts)) {
#           if(max_zero_count >= threshold) {
#             ext = ext + 1
#           }
#         }
#       }
#       else if(ts[i] >= 1 & max_zero_count >= threshold) {
#         if(!any(ts[which(!is.na(ts[1:(i-1)]))] != 0) & first(ts) == 0) {
#           ext = 0
#         }
#         else {
#           ext = ext + 1
#         }
#         est = est + 1
#         zero_count = 0
#         max_zero_count = 0
#       }
#       else if(ts[i] >= 1 & max_zero_count <= threshold) {
#         zero_count = 0
#         max_zero_count = 0
#       }
#     }
#     
#     if(any(ts[which(!is.na(ts))[1:2]] >= 1)) {
#       new_est = FALSE
#     }
#     else {
#       new_est = TRUE
#     }
#     
#     ext_est = df %>%
#       mutate(n_est = est, n_ext = ext, new_est = new_est, n_ts = length(ts), p = df$p, beta = df$beta)
#     
#     return(ext_est)
#   })
#   
#   ## bind them all 
#   df <- bind_rows(list)
#   
#   ## problem: no. of local extinctions/establishments is lower when range goes completely extinct before 1500 time steps
#   ## solution: calculate rate of extinction/establishment by dividing it by maximum
#   
#   df <- df %>%
#     group_by(beta, p, t) %>%
#     filter(!all(Nt == 0)) %>%
#     group_by(beta, p) %>%
#     mutate(t_max = max(t),
#            mean_rate_ext = mean(n_ext/t_max),
#            mean_rate_est = mean(n_est/t_max)) %>%
#     ungroup() %>%
#     select(x, y, n_ext, n_est, p, beta, t_max, mean_rate_ext, mean_rate_est) %>%
#     distinct()
#   
#   ## plot histogram
#   # df %>%
#   #   mutate(p = paste0("p", p), beta = paste0("beta", beta)) %>%
#   #   ggplot(aes(x = n_ext/t_max)) +
#   #   geom_histogram() +
#   #   facet_wrap(p~beta) +
#   #   labs(x = "Local extinction rate", y = "Frequency") +
#   #   geom_vline(aes(xintercept = mean_rate_ext), colour = "red")
#   # 
#   # df %>%
#   #   mutate(p = paste0("p", p), beta = paste0("beta", beta)) %>%
#   #   ggplot(aes(x = n_est/t_max)) +
#   #   geom_histogram() +
#   #   facet_wrap(p~beta) +
#   #   labs(x = "Local establishment rate", y = "Frequency") +
#   #   geom_vline(aes(xintercept = mean_rate_est), colour = "red")
#   
#   return(df)
# }

# est_ext_grid <- function(all_ranges, threshold) {
#   split = all_ranges %>%
#     group_by(x, y, p, beta) %>%
#     group_split(.)
#   
#   list = lapply(split, FUN = function(df) {
#     ts = df$Nt
#     zero_count = 0
#     ext = 0
#     est = 0
#     max_zero_count <- 0
#     zero_count_vec <- c()
#     for(i in 1:length(ts)) {
#       if(i == 1) {
#         if(ts[i] == 0) {
#           zero_count = 1
#         }
#       }
#       else if(is.na(ts[i]))  {
#         if(zero_count > max_zero_count) {
#           max_zero_count = zero_count
#         }
#         zero_count = 0
#         if(i == length(ts)) {
#           if(max_zero_count >= threshold) {
#             ext = ext + 1
#           }
#         }
#       }
#       else if(ts[i] == 0) {
#         zero_count = zero_count + 1
#         if(zero_count > max_zero_count) {
#           max_zero_count = zero_count
#         }
#         if(i == length(ts)) {
#           if(max_zero_count >= threshold) {
#             ext = ext + 1
#           }
#         }
#       }
#       else if(ts[i] >= 1 & max_zero_count >= threshold) {
#         if(!any(ts[which(!is.na(ts[1:(i-1)]))] != 0) & first(ts) == 0) {
#           ext = 0
#         }
#         else {
#           ext = ext + 1
#         }
#         est = est + 1
#         zero_count = 0
#         max_zero_count = 0
#       }
#       else if(ts[i] >= 1 & max_zero_count <= threshold) {
#         zero_count = 0
#         max_zero_count = 0
#       }
#       zero_count_vec <- append(zero_count_vec, zero_count)
#     }
#     
#     ext = zero_count_vec 
#     ext[which(ext != 1)] = 0
#     ext[1] = 0
#     
#     est = zero_count_vec 
#     shifted = zero_count_vec[2:nrow(df)]
#     est = ifelse(shifted - est[1:(nrow(df) - 1)] < 0, 1, 0) 
#     est = c(0, est)
#     
#     return(data.frame(x = unique(df$x), y = unique(df$y), 
#                       p = unique(df$p), beta = unique(df$beta), 
#                       t = 1:nrow(df),
#                       est = est, ext = ext))
#     
#   })
#   ## bind them all 
#   df <- bind_rows(list)
#   
#   ## left join to other stats
#   all_ranges = left_join(all_ranges, df)
#   
#   return(all_ranges)
# }



## calculate mean shift rate and compare to mean climate velocity
shift_params <- calculate_mean_shift_params(all_ranges = all_ranges)

shift_params %>%
  mutate(clim_velo = unique(all_ranges$shift_rate)) %>%
  filter(parameter != "abd_centroid") %>%
  mutate(p_beta = paste(p, beta, sep = "_")) %>%
  ggplot(aes(x = parameter, y = shift_rate, colour = p_beta)) +
  geom_boxplot() +
  geom_point(aes(x = parameter, y = unique(all_ranges$shift_rate)))


## plot ranges in each time step 
plot_range(range_shift = all_ranges, path = "outputs/figures/range-shifts_dispersal_new")

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

