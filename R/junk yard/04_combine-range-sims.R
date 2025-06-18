## script that combines simulation results into one csv 
library(tidyverse)

d = 0.1 # proportion of offspring dispersing
icp = 0.1 ## intraspecific competition parameter

if(!dir.exists("outputs/data-processed/stable-ranges/all-ranges")) {
  dir.create("outputs/data-processed/stable-ranges/all-ranges", recursive = T)
}

for(p in c(0, 1)) {
  for(beta in c(0, 1)) {
    ## organize stable ranges 
    stable_files = list.files(paste0("outputs/data-processed/stable-ranges"), full.names = T)
    stable_files = stable_files[str_detect(stable_files, paste0("p", p, "_b", beta, "_icp", icp, "_d", d))]
    
    stable_files = list.files(stable_files, full.names = T)
    
    ## for each rep
    for(rep in 1:length(stable_files)) {
      cur_files = list.files(stable_files[rep], full.names = T)
      r = str_split_fixed(stable_files[rep], "/rep", 2)[,2]
      
      all_ranges = c()
      i=1
      while(i <= length(cur_files)) {
        all_ranges = rbind(all_ranges, read.csv(cur_files[i]))
        i = i + 1
      }
      write.csv(all_ranges, paste0("outputs/data-processed/stable-ranges/all-ranges/all-ranges_rep", r, "_p", p, "_b", beta, "_icp", icp, "_d", d),
                row.names = FALSE)
    }
  }
}


if(!dir.exists("outputs/data-processed/range-shifts/all-ranges")) {
  dir.create("outputs/data-processed/range-shifts/all-ranges", recursive = T)
}

for(p in c(0, 1)) {
  for(beta in c(0, 1)) {
    ## organize stable ranges 
    stable_files = list.files(paste0("outputs/data-processed/range-shifts"), full.names = T)
    stable_files = stable_files[str_detect(stable_files, paste0("p", p, "_b", beta, "_icp", icp, "_d", d))]
    
    stable_files = list.files(stable_files, full.names = T)
    
    ## for each rep
    for(rep in 1:length(stable_files)) {
      cur_files = list.files(stable_files[rep], full.names = T)
      r = str_split_fixed(stable_files[rep], "/rep", 2)[,2]
      
      all_ranges = c()
      i=1
      while(i <= length(cur_files)) {
        all_ranges = rbind(all_ranges, read.csv(cur_files[i]))
        i = i + 1
      }
      write.csv(all_ranges, paste0("outputs/data-processed/range-shifts/all-ranges/range-shifts_rep", r, "_p", p, "_b", beta, "_icp", icp, "_d", d),
                row.names = FALSE)
    }
  }
}
