#######################################
###      setting up environment      ## 
#######################################
.libPaths( c( "~/projects/def-jsunday/nikkim/pop-sims/packages" , .libPaths() ) )
library(tidyverse)
library(parallel)
library(foreach)
library(doParallel)
library(raster)
library(terra)
theme_set(theme_bw())

## set up parallel processing:
starts <- rep(100, 40)
fx <- function(nstart) kmeans(Boston, 4, nstart=nstart)
numCores <- detectCores() ## detect cores
numCores

registerDoParallel(numCores)  ## use multicore, set to the number of our cores

source("R/02_simulate-range-shifts.R")

betas <- c(0, 1)
ps <- c(0, 1)

for(p in ps) {
  for(beta in betas) {
    filename = paste0("outputs/data-processed/stable-ranges/stable-ranges_p", p, 
                      "_beta", beta, "_r1.2_K100_d0.1_icp0.1_L500_reps10.rds")
    ## read in stable ranges
    stable_ranges = readRDS(filename)
    
    ## call function to simulate range shifts 
    simulate_range_shifts(ps = ps, betas = beta, stable_ranges = stable_ranges)
  }
}

