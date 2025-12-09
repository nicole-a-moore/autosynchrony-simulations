## creates plots of a set of stable ranges and saves them to a folder 
plot_range = function(range_shift, path) {
  
  nrow = length(unique(range_shift$y))
  ncol = length(unique(range_shift$x))
  
  ## get ps and betas
  ps = unique(all_ranges$p)
  betas = unique(all_ranges$beta)
  
  for(p in ps) {
    for(beta in betas) {
      
      range_shift_sub = range_shift[which(range_shift$p == p & range_shift$beta == beta),]
      
      n_ts = length(unique(range_shift_sub$t))
      shift_rate = unique(range_shift_sub$shift_rate)
      
      ## calculate location of max suitability over time
      ## distance at which r < 0 
      ## for these parameters, 41
      max_suit = data.frame(time = 1:n_ts,
                            max_suit = c(rep(41, length.out = 500), shift_rate*(1:(n_ts-500)) + 41))
      
      time = 1
      while(time <= n_ts) {
        
        ## filter to time t
        cur = range_shift_sub %>%
          filter(t == time) 
        
        ## calculate metrics
        params = cur %>%
          filter(Nt != 0) %>%
          summarize(q95_y = quantile(y, c(0.95)),
                    q5_y = quantile(y, c(0.05)),
                    abd_centroid = weighted.mean(y, Nt))
        
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
            scale_fill_continuous(limits = c(0, 200), type = "viridis", na.value = "white") +
            theme(panel.border = element_rect(colour = "black", fill = NA)) +
            geom_hline(yintercept = params$q95_y, colour = "forestgreen") +
            geom_hline(yintercept = max_suit$max_suit[which(max_suit$time == time)], colour = "red") 
            
            
          dir =  paste0(path, "/rep", unique(range_shift_sub$rep), "_p", 
                        unique(range_shift_sub$p), "_b", 
                        unique(range_shift_sub$beta), "_icp", 
                        unique(range_shift_sub$icp))
          if(dir.exists(dir) == FALSE) {
            dir.create(path = dir)
          }
          
          filename = paste0(dir, "/t", time, "_rep", unique(range_shift_sub$rep), "_p", 
                            unique(range_shift_sub$p), "_b", 
                            unique(range_shift_sub$beta), "_icp", 
                            unique(range_shift_sub$icp), "_shift", ".png")
          ## save plot
          ggsave(plot, filename = filename, width = 2, height = 4)
        }
        else {
          time = length(unique(range_shift_sub$t))
        }

        print(time)
        time = time + 1
      }
    }
  }
  
}

