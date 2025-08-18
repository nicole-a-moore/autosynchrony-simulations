## script that combines range shift simulation results into one csv 
library(tidyverse)

d = 0.1 # proportion of offspring dispersing
icp = 0.1 ## intraspecific competition parameter

if(!dir.exists("outputs/data-processed/range-shift-simulations/all-ranges")) {
  dir.create("outputs/data-processed/range-shift-simulations/all-ranges", recursive = T)
}

for(p in c(0, 1)) {
  for(beta in c(0, 1)) {
    ## organize range shift files 
    files = list.files(paste0("outputs/data-processed/range-shift-simulations"), full.names = T)
    files = files[str_detect(files, paste0("p", p, "_b", beta, "_icp", icp, "_d", d))]
    # 
    # files_dispersal = files[str_detect(files, "dispersal")]
    # files_dispersal = list.files(files_dispersal, full.names = T)
    
    files = files[!str_detect(files, "dispersal")]
    files = list.files(files, full.names = T)
    
    ## for each rep
    for(rep in 1:length(files)) {
      cur_files = list.files(files[rep], full.names = T)
      r = str_split_fixed(files[rep], "/wt_rep", 2)[,2]
      
      all_ranges = c()
      i=1
      while(i <= length(cur_files)) {
        all_ranges = rbind(all_ranges, read.csv(cur_files[i]))
        print(i)
        i = i + 1
      }
      write.csv(all_ranges, paste0("outputs/data-processed/range-shift-simulations/all-ranges/range-shifts_rep", r, "_p", p, "_b", beta, "_icp", icp, "_d", d, ".csv"),
                row.names = FALSE)
    }
    
    # ## for each rep
    # for(rep in 1:length(files_dispersal)) {
    #   cur_files = list.files(files_dispersal[rep], full.names = T)
    #   r = str_split_fixed(files_dispersal[rep], "/rep", 2)[,2]
    #   
    #   all_ranges = c()
    #   i=1
    #   while(i <= length(cur_files)) {
    #     all_ranges = rbind(all_ranges, read.csv(cur_files[i]))
    #     i = i + 1
    #   }
    #   write.csv(all_ranges, paste0("outputs/data-processed/range-shift-simulations/all-ranges/range-shifts_rep", r, "_p", p, "_b", beta, "_icp", icp, "_d", d, "_dispersal"),
    #             row.names = FALSE)
    # }
  }
}
