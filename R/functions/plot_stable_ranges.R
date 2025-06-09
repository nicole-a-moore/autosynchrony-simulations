## creates plots of a set of stable ranges and saves them to a folder 
plot_stable_ranges = function(stable_ranges, 
                             path, filename) {
  
  ## get dimensions
  first = stable_ranges[[1]]
  nrow = nrow(first)
  ncol = ncol(first)
  
  for(range in 1:length(stable_ranges)) {
    ## create raster 
    ras <- raster(nrow = nrow, ncols = ncol, xmn = 0, xmx = ncol, ymn = 0, ymx = nrow)
    values(ras) <- c(t(stable_ranges[[range]]))
    
    if(length(which(values(ras) == 0)) == length(values(ras))) {
      ## plot raster
      plot = ras %>% 
        rasterToPoints() %>%
        as.data.frame() %>%
        ## make 0 into NA
        mutate(layer = ifelse(layer == 0, NA, layer)) %>%
        ggplot(aes(x = x, y = y, fill = layer)) +
        geom_raster() +
        coord_fixed() +
        theme_void() +
        labs(fill = "") +
        #scale_fill_continuous(limits = c(0, 200), type = "viridis", na.value = "white") +
        theme(panel.border = element_rect(colour = "black", fill = NA))
    }
    else {
      ## plot raster
      plot = ras %>% 
        rasterToPoints() %>%
        as.data.frame() %>%
        ## make 0 into NA
        mutate(layer = ifelse(layer == 0, NA, layer)) %>%
        ggplot(aes(x = x, y = y, fill = layer)) +
        geom_raster() +
        coord_fixed() +
        theme_void() +
        labs(fill = "") +
        scale_fill_continuous(limits = c(0, 200), type = "viridis", na.value = "white") +
        theme(panel.border = element_rect(colour = "black", fill = NA))
    }
    
    newname = strsplit(filename, ".rds")[[1]][1]
    newname = strsplit(newname, "stable-ranges/")[[1]][2]
    
    ## save plot
    ggsave(plot, filename = paste0(path, "/range_", range, "_", newname, ".png"), width = 2, height = 4)
  }
}

