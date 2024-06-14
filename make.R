#######################################
###      setting up environment      ## 
#######################################
library(tidyverse)
library(parallel)
library(foreach)
library(doParallel)
library(raster)
theme_set(theme_bw())

## set up parallel processing:
starts <- rep(100, 40)
fx <- function(nstart) kmeans(Boston, 4, nstart=nstart)
numCores <- detectCores() ## detect cores
numCores

registerDoParallel(numCores)  ## use multicore, set to the number of our cores


source("R/01_generate-stable-ranges.R")
source("R/functions/plot_stable_ranges.R")
source("R/02_simulate-range-shifts.R")

## generate stable ranges 
stable_ranges = generate_stable_ranges(ncol = 10, nrow = 100, reps = 10)

## save the stable ranges:
saveRDS(stable_ranges, "outputs/data-processed/stable_ranges.rds")
#stable_ranges <- readRDS("outputs/data-processed/stable_ranges.rds")

## plot the stable ranges
plot_stable_ranges(stable_ranges, path = "outputs/figures/stable-range-simulations")

## simulate range shifts 
simulate_range_shifts(stable_ranges, path = "outputs/data-processed/range-shifts/")

