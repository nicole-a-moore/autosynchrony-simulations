## try measuring climate velocity for each cell of simulation and comparing to shift rate 
library(terra)
library(tidyverse)

env_grid = rast("outputs/data-processed/env-grids/range-shift-grid1_p0_beta0_r1.2_K100_d0.1_icp0.1_L2000.tif")

array = as.array(env_grid)

nrow = dim(array)[1]
ncol = dim(array)[2]
ntime = dim(array)[3]

#ntime = 700

ar2d = array[,,1]
# for each cell, calculate temporal trend 
for(y in 1:nrow) {
  for(x in 1:ncol) {
    df = data.frame(x = 501:ntime, y = array[y,x,501:ntime])
    lm = lm(y ~ x, data = df)
    
   # df %>% ggplot(aes(x = x, y = y)) + geom_point()
    
    ar2d[y, x] = lm$coefficients[2]
  }
}

hist(ar2d)

plot(rast(ar2d))

library(climetrics)

r = stack("outputs/data-processed/env-grids/range-shift-grid1_p0_beta0_r1.2_K100_d0.1_icp0.1_L2000.tif")

# dvel = dVelocity(r, t1 = 601, t2 = 701, ny = 100)
# 
# plot(dvel)
# 
# mean(values(dvel), na.rm = T)


gvel = gVelocity(r[[501:1251]])

plot(gvel)

hist(gvel)

mean(values(gvel))
