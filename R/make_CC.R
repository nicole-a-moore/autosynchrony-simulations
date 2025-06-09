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

## read the command line arguments
# command_args <- commandArgs(trailingOnly = TRUE)
# beta = as.numeric(command_args[1])
# p = as.numeric(command_args[2])

## read in functions
source("R/01_generate-stable-ranges_d.R")

# generate stable ranges
for(i in 0:1) {
  for(t in 0:1) {
    filename = generate_stable_ranges(beta = i, p = t, ncol = 10, nrow = 100, reps = 10)
  }
}
