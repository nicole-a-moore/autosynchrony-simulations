#######################################
###      setting up environment      ## 
#######################################
.libPaths( c( "~/projects/def-jsunday/nikkim/packages" , .libPaths() ) )
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

source("R/01_simulate-range-shifts.R")

betas <- c(0, 1)
ps <- c(0, 1)

for(p in ps) {
  for(beta in betas) {
    ## call function to simulate range shifts 
    simulate_range_shifts(p = p, beta = beta)
  }
}

