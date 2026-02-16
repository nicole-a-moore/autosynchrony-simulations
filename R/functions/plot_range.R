## creates plots of a set of stable ranges and saves them to a folder 
plot_range = function(path,
                      p,
                      beta,
                      d,
                      d_dist,
                      shift_rate,
                      icp,
                      K,
                      sigma,
                      rep) {


  ## read in file
  range_shift = read.csv(paste0("outputs/data-processed/range-shift-simulations/cluster_both-edges/sim-results/",
                             "range-shifts_p", p, "_b", beta, "_icp", icp, "_K", K, "_d", d, "_r_d-dist", d_dist,
                             "_sigma", sigma, "_shift-rate", shift_rate, "_rep", rep, ".csv"))
  
  
  nrow = length(unique(range_shift$y))
  ncol = length(unique(range_shift$x))

  n_ts = length(unique(range_shift$t))
  
  ## calculate location of max suitability over time
  ## distance at which r < 0 
  ## for these parameters, 41
  max_suit = data.frame(time = 1:n_ts,
                        max_suit = c(rep(41, length.out = 500), shift_rate*(1:(n_ts-500)) + 41))
  
  time = 500
  while(time <= 2000) {
    
    ## filter to time t
    cur = range_shift %>%
      filter(t == time) 
    
    ## calculate metrics
    params = cur %>%
      filter(Nt != 0) %>%
      summarize(q95_y = quantile(y, c(0.95)),
                q5_y = quantile(y, c(0.05)),
                max_y = max(y, na.rm = T),
                min_y = min(y, na.rm = T))
    
    params2 = cur %>%
      filter(Nt != 0) %>%
      arrange(-y) %>%
      slice_head(n = 5) %>%
      mutate(meanmax_y = mean(y, na.rm = T)) 
    
    params3 = cur %>%
      filter(Nt != 0) %>%
      arrange(y) %>%
      slice_head(n = 5) %>%
      mutate(meanmim_y = mean(y, na.rm = T))
    
    params$meanmim_y = unique(params3$meanmim_y)
    params$meanmax_y = unique(params2$meanmax_y)
    
    mid_occ_lat = cur[which(cur$Nt != 0),]
    mid_occ_lat = mean(mid_occ_lat$y, na.rm = TRUE)
    
    if(!length(which(cur$Nt == 0)) == length(cur$Nt)) {
      ## plot raster
      plot = cur %>% 
        ## make 0 into NA
        mutate(layer = ifelse(Nt == 0, NA, Nt)) %>%
        ggplot(aes(x = x, y = y, fill = layer)) +
        geom_raster() +
        coord_fixed() +
        theme_void() +
        labs(fill = "") +
        scale_fill_continuous(limits = c(0, 300), type = "viridis", na.value = "white",
                              alpha = 0.8) +
        theme(panel.border = element_rect(colour = "black", fill = NA)) +
        geom_hline(yintercept = params$q95_y, colour = "blue",linetype = "dashed", alpha = 1, size = 0.5) +
        geom_hline(yintercept = params$q5_y, colour = "red", linetype = "dashed", alpha = 1, size = 0.5) +
        geom_hline(yintercept = params$max_y, colour = "blue", alpha = 1, size = 0.5) +
        geom_hline(yintercept = params$min_y, colour = "red", alpha = 1, size = 0.5) +
        geom_hline(yintercept = params$meanmax_y, colour = "blue", linetype = "dotted", alpha = 1, size = 0.5) +
        geom_hline(yintercept = params$meanmin_y, colour = "red", linetype = "dotted", alpha = 1, size = 0.5) 
      
      
      dir =  paste0(path, "/p", p, "_b", beta, "_icp", icp, "_K", K, "_d", d, "_r_d-dist", d_dist,
                    "_sigma", sigma, "_shift-rate", shift_rate, "_rep", rep)
      if(dir.exists(dir) == FALSE) {
        dir.create(path = dir)
      }
      
      filename = paste0(dir, "/t", time, ".png")
      ## save plot
      ggsave(plot, filename = filename, height = 7, width = 2)
    }
    else {
      time = length(unique(range_shift$t))
    }
    
    print(time)
    time = time + 1
  }
  
}

