## playing with climate velocity 

### try measuring shift in position of maximum suitability 
## first try with simulations and see if the shift in latitude of max suitability measured is what was specified by me 
library(tidyverse)
select = dplyr::select

## read in growth rate grid over time 
env = brick("outputs/data-processed/env-grids/range-shift-grid2_p1_beta0_r1.2_K100_d0.1_icp0.1_L2000.tif")

plot(env[[1]])

## convert to df
env = as.data.frame(rasterToPoints(env))

## get lat position of max suitability in each layer
max_suit = sapply(3:ncol(env), FUN = function(ind) {env$y[which(env[,ind] == max(env[,ind]))]})
plot(x = 1:2000, y = max_suit)

## yes, it is 

## now try measuring climate velocity 

## on a grid with no noise, just shift
library(VoCC)
library(raster)
rast = raster::brick(lattice_E_it_array, env[[1]])

ttrend = tempTrend(r = rast,
                   th = 0.25*nlayers(rast) ## set minimum # obs. to 1/4 time series length
)
plot(ttrend)

## use function to calculate spatial gradient in mean daily air temperature for each location:
spgrad = spatGrad(r = rast, 
                  projected = TRUE) ## our raster is projected to a coordinate system

## plot it: 
plot(spgrad)

## calculate gradient based climate velocity:
gvocc = gVoCC(tempTrend = ttrend, spatGrad = spgrad)

## plot it:
plot(gvocc)

mean(values(gvocc$voccMag))
shift_rate


## revelation: duh, r =/= temperature

## first, make an array with temperature that shifts at a given rate
## then, use curve to translate to suitability 

## cellular lattice of microcosms
lattice_N_it = array(dim = c(nrow,ncol,L))
lattice_r = matrix(ncol = ncol, nrow = nrow)
lattice_E_it = matrix(ncol = ncol, nrow = nrow)

## higher proportion dispersing = less pronounced effect of suitability gradient
lattice_r[1:nrow,1:ncol] <- r ## start with growth rate = max growth rate 
lattice_N_it[1:(nrow/2),1:ncol,1] <- K/2 ## start with population size = carrying capacity / 2 in 1/2 of grid
lattice_N_it[(nrow/2):nrow,1:ncol,1] <- 0 ## start with population size = 0 in other half 

## position optimum temperature (20deg) conditions as row 75 on the lattice (Emax)
## create temperature gradient 
lattice_E_it
y = 0.25*(1:100) + 1.25

lattice_E_it[,1:10] = y

## replicate latitudinal gradient L times 
lattice_E_it_array <- replicate(L, lattice_E_it)

## after stable period of 500 time steps, shift temperatures by shift_rate 
plot(x = 1:nrow, y = lattice_E_it_array[1:nrow,1,1])

for(t in 501:2000) {
  lattice_E_it_array[,1:10,t] = (0.25*(1:100) + 1.25) + 0.025*(t-500)
  ## 0.25 = spatial gradient
  ## 0.025 = temporal warming rate
  ## clim velocity  = 0.025/0.25 = 0.1
}

plot(x = 1:nrow, y = lattice_E_it_array[1:nrow,1,1000])

## read in growth rate grid over time 
env = brick("outputs/data-processed/env-grids/range-shift-grid1_p1_beta1_r1.2_K100_d0.1_icp0.1_L2000.tif")

## measure climate velocity 
rast = raster::brick(lattice_E_it_array[,,501:2000], env[[1]])

ttrend = tempTrend(r = rast,
                   th = 0.25*nlayers(rast) ## set minimum # obs. to 1/4 time series length
)
plot(ttrend)
mean(values(ttrend$slpTrends))

## use function to calculate spatial gradient in mean daily air temperature for each location:
spgrad = spatGrad(r = rast, 
                  projected = TRUE) ## our raster is projected to a coordinate system

## plot it: 
plot(spgrad)
mean(values(spgrad$Grad))

## calculate gradient based climate velocity:
gvocc = gVoCC(tempTrend = ttrend, spatGrad = spgrad)

## plot it:
plot(gvocc)

theoretical_climvel = 0.025/0.25
mean(values(gvocc$voccMag))

## make tpc: (from Huey suboptimal is optimal)
## plot a curve!
r = 1.2+0.25 # max suitability
alpha = 1 ## rise rate steepness
beta = 1 ## decline rate steepness
Trmax = 20
## range of body temps 
Tb = seq(from = 0, to = 40, by = 0.1)

## calculate fitness
exponent = -exp(beta*(Tb-Trmax)-8)-alpha*(Tb-Trmax)^2
wb = r*exp(1)^exponent - 0.25

## plot
df <- data.frame(relative_fitness = wb, body_temperature = Tb)

df %>%
  ggplot(aes(x = body_temperature, y = relative_fitness)) +
  geom_line()


## use tpc to translate temps into suitability
lattice_r_array = lattice_E_it_array
lattice_r_array = -exp(beta*(lattice_E_it_array-Trmax)-8)-alpha*(lattice_r_array-Trmax)^2
lattice_r_array = r*exp(1)^lattice_r_array - 0.25

plot(x = 1:nrow, y = lattice_r_array[1:nrow,1,600])


## now add noise 
## create noise time series for each cell, L time steps long
lattice_ac_it = array(dim = c(nrow,ncol,L))

# generate noise of given synchrony and autocorrelation
noise <- generate_noise(beta = beta,
                        p = p,
                        n_ts = nrow*ncol, ## number of cells
                        L1 = 500,
                        L2 = L - 500)

## verify the degree of autocorrelation and synchrony: 
## calculate mean measured spectral exponent for each time series
beta_star =  mean(sapply(noise[[2]], cbind))

## measure cross correlation:
## first-difference the time series values
ts_all <- data.frame(sapply(noise[[1]], cbind))
diffs <- ts_all[2:(nrow(ts_all)-1),] - ts_all[1:(nrow(ts_all)-2),]

cors <- c()
i=1
while(i <= ncol(diffs)) {
  for(n in ncol(diffs)) {
    if(n != i) {
      cors <- append(cors, cor(diffs[,i], diffs[,n], method = 'pearson'))
    }
  }
  i = i + 1
}
p_star = mean(cors)

## assign noise to each cell in the lattice
df <- data.frame(sapply(noise[[1]], cbind))
n = 1
for(a in 1:nrow) {
  for(b in 1:ncol) {
    lattice_ac_it[a,b,1:L] <- df[,n]
    n = n + 1
  }
}

## let noise affect r for each cell in the lattice
lattice_noise_array = (lattice_E_it_array + lattice_ac_it)
#plot(x = 1:nrow, y = lattice_noise_array[1:nrow,1,805])

## use tpc to translate temps into suitability
lattice_r_array = lattice_noise_array
lattice_r_array = -exp(beta*(lattice_noise_array-Trmax)-8)-alpha*(lattice_noise_array-Trmax)^2
lattice_r_array = r*exp(1)^lattice_r_array - 0.25

plot(x = 1:nrow, y = lattice_r_array[1:nrow,1, 700])

plot(x = 1:2000, y = lattice_r_array[50,1,1:2000])


## now try with Berkeley Earth data and hypothetical species response curve 


## get lat position of max suitability in each layer
## convert to df
env_df = as.data.frame(rasterToPoints(brick(lattice_r_array)))
max_suit = sapply(3:ncol(env_df), FUN = function(ind) {mean(env_df$y[which(env_df[,ind] == max(env_df[,ind]))])})
plot(x = 1:2000, y = max_suit)

0.025/0.25

## measure clim velocity after variation was added
## measure climate velocity 
rast = raster::brick(lattice_noise_array[,,501:2000], env[[1]])

ttrend = tempTrend(r = rast,
                   th = 0.25*nlayers(rast) ## set minimum # obs. to 1/4 time series length
)
plot(ttrend)
mean(values(ttrend$slpTrends))
#0.2499

## use function to calculate spatial gradient in mean daily air temperature for each location:
spgrad = spatGrad(r = rast, 
                  projected = TRUE) ## our raster is projected to a coordinate system

## plot it: 
plot(spgrad)
mean(values(spgrad$Grad))
## 0.2852

## calculate gradient based climate velocity:
gvocc = gVoCC(tempTrend = ttrend, spatGrad = spgrad)

## plot it:
plot(gvocc)

theoretical_climvel = 0.025/0.25
mean(values(gvocc$voccMag))
## 0.106


