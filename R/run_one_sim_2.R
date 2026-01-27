.libPaths( c( "~/projects/def-jsunday/nikkim/sdms/packages" , .libPaths() ) )
library(terra)
library(tidyverse)
terraOptions(threads = 1)

## load simulation function
source("01_simulate-range-shifts_cluster.R")

## list all parameters to vary
beta <- c(0, 1)
p <- c(0, 1)
icp <- c(0.1)
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

## pick this job's row
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## get this job's params
params <- sims[task_id, ]

## output directory
outdir <- file.path(paste0("output/simulations/p", params$p, "_b", params$beta, "_icp", params$icp, "_K", params$K, "_d", 
                           params$d, "_r", params$r, "_d-dist", params$d_dist, "_sigma", params$sigma, "_shift-rate", params$shift_rate))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## run one simulation
simulate_range_shifts(
  beta = params$beta,
  p = params$p,
  icp = params$icp,
  d = params$d, 
  d_dist = params$d_dist,
  sigma = params$sigma,
  K = params$K,
  shift_rate = params$shift_rate,
  r = 1,
  L = 2000, 
  nrow = 300,
  ncol = 10,
  path = outdir
)