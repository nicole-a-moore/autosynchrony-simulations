#######################################
###      setting up environment      ## 
#######################################
.libPaths( c( "~/projects/def-jsunday/nikkim/sdms/packages" , .libPaths() ) )
library(tidyverse)
library(parallel)
library(foreach)
library(doParallel)
library(raster)
library(terra)

## set up parallel processing:
starts <- rep(100, 40)
fx <- function(nstart) kmeans(Boston, 4, nstart=nstart)
numCores <- detectCores() ## detect cores
numCores

registerDoParallel(numCores)  ## use multicore, set to the number of our cores

source("01_simulate-range-shifts_cluster.R")

## list all parameters to vary
beta <- c(0, 1)
p <- c(0, 1)
icp <- c(0.7)
d <- c(0.2)
d_dist <- c(1,4,7)
sigma <- c(0.2, 0.5)
K <- c(200)
shift_rate <- c(0.2)

## create list of all combinations of parameters to vary 
sims = expand.grid(beta = beta,
            p = p,
            icp = icp,
            d = d, 
            d_dist = d_dist,
            sigma = sigma,
            K = K,
            shift_rate = shift_rate)

## for each combination of parameters 
foreach(row = 1:nrow(sims)) %dopar% {
  simulate_range_shifts(beta = sims$beta[row],
                        p = sims$p[row],
                        icp = sims$icp[row],
                        d = sims$d[row], 
                        d_dist = sims$d_dist[row],
                        sigma = sims$sigma[row],
                        K = sims$K[row],
                        shift_rate = sims$shift_rate[row],
                        r = 1,
                        L = 2000, 
                        nrow = 300,
                        ncol = 10,
                        path = "output_new")
}


