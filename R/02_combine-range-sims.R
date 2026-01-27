## script that combines range shift simulation results into one csv 
library(tidyverse)

## make a folder to store results
dir = "outputs/data-processed/range-shift-simulations_dispersal/sim-results"
if(!dir.exists(dir)) {
  dir.create(dir, recursive = T)
}

## go through each folder 
paths = list.files("outputs/data-processed/range-shift-simulations_dispersal/raw-sims", full.names = T)
folders = list.files("outputs/data-processed/range-shift-simulations_dispersal/raw-sims")

f = 1
while(f <= length(folders)) {
  
  ## get rep folders within folder 
  rep_folders = list.files(paths[f], full.names = T)

  ## for each rep
  for(rep in 1:length(rep_folders)) {
    cur_files = list.files(rep_folders[rep], full.names = T)
    r = str_split_fixed(rep_folders[rep], "/rep", 2)[,2]
    
    all_ranges = c()
    i=1
    while(i <= length(cur_files)) {
      all_ranges = rbind(all_ranges, read.csv(cur_files[i]))
      i = i + 1
    }
    
    if(!is.null(nrow(all_ranges))) {
      write.csv(all_ranges, 
                paste0(dir, "/range-shifts_", folders[f], "_rep", r, ".csv"),
                row.names = FALSE)
    }

  }
  
  print(f)
  f = f + 1
}
