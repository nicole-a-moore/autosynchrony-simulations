#######################################
###      setting up environment      ## 
#######################################
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

## read in functions
source("R/01_generate-stable-ranges.R")
source("R/functions/plot_stable_ranges.R")
source("R/02_simulate-range-shifts.R")

## generate stable ranges with low autocorr (beta = 0)
filename = generate_stable_ranges(beta = 0, ncol = 10, nrow = 100, reps = 10)

## read in the ranges
stable_ranges <- readRDS(filename)

## plot the stable ranges
plot_stable_ranges(stable_ranges, filename = filename, path = "outputs/figures/stable-range-simulations")

## generate more with beta = 1 (high autocorr)
filename = generate_stable_ranges(beta = 1, ncol = 10, nrow = 100, reps = 10)

## read in the ranges
stable_ranges <- readRDS(filename)

## plot the stable ranges
plot_stable_ranges(stable_ranges, filename = filename, path = "outputs/figures/stable-range-simulations")

#######
## simulate range shifts 
simulate_range_shifts(stable_ranges, filename = filename, path = "outputs/data-processed/range-shifts/")

