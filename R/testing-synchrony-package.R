install.packages("synchrony")
library(synchrony)

## generate a community of synchronoues time series
comm <- correlated.matrix(nspecies = 10, ntimes = 100, rho = 0.7)[[4]]

comm 

df = data.frame(time = 1:100, sp1 = comm[,1], sp2 = comm[,2])

df %>%
  filter(time %in% 1:50) %>%
  ggplot(aes(x = time, y = sp1)) + 
  geom_line() +
  geom_line(aes(y = sp2), colour = "red")

meancorr(data = comm, nrands = 999, alternative = "two.tailed", type = 1, quiet = T)

## read in generated noise and test whether synchrony is as specified
p1 = rast("outputs/data-processed/env-grids/grid1_p1_beta0_r1.2_K100_d0.1_icp0.1_L500_reps10.tif")
p1_df = as.data.frame(p1, xy = TRUE)
p0 = rast("outputs/data-processed/env-grids/shift_grid10_p0_beta0_r1.2_K100_d0.1_icp0.1_L1500_reps100.tif")
p0_df = as.data.frame(p0, xy = TRUE)

p1_est = meancorr(data = p1_df[,3:11], nrands = 999, alternative = "two.tailed", type = 1, quiet = T)
p0_est = meancorr(data = p0_df[,3:11], nrands = 999, alternative = "two.tailed", type = 1, quiet = T)
## it is :-)

p1_temp = data.frame(time = 1:(ncol(p1_df)-2), sp1 = as.numeric(as.character(p1_df[1,3:502])), 
                     sp2 = as.numeric(as.character(p1_df[2,3:502])))

p1_temp %>%
  #filter(time %in% 1:50) %>%
  ggplot(aes(x = time, y = sp1)) + 
  geom_line() +
  geom_line(aes(y = sp2), colour = "red")

## measure phase synchrony
phase1 <- phase.sync(p1_temp$sp1, p1_temp$sp2, nrands = 999, quiet = TRUE) 
## 0 phase difference at any time step 

## compute concurrency (whether time series peaks and troughs at the same time)
peaks1 <- peaks(p1_temp$sp1, p1_temp$sp2, nrands = 999, type = 1, quiet = TRUE)
## proportion of common peaks: 0.99

p0_temp = data.frame(time = 1:(ncol(p0_df)-2), sp1 = as.numeric(as.character(p0_df[1,3:502])), 
                     sp2 = as.numeric(as.character(p0_df[2,3:502])))

p0_temp %>%
  #filter(time %in% 1:50) %>%
  ggplot(aes(x = time, y = sp1)) + 
  geom_line() +
  geom_line(aes(y = sp2), colour = "red")

phase0 <- phase.sync(p0_temp$sp1, p0_temp$sp2, nrands = 999, quiet = TRUE) 
## 0 phase difference at any time step 

## compute concurrency (whether time series peaks and troughs at the same time)
peaks0 <- peaks(p0_temp$sp1, p0_temp$sp2, nrands = 999, type = 1, quiet = TRUE)
## proportion of common peaks: 0.44  

## generate some fake time series and test
source("R/functions/generate_noise.R")
noise = generate_noise(beta = 1,
                      p = 0, 
                      n_ts = 100,
                      L = 1000)

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
p_star

noise = noise[[1]]

df = data.frame(time = 1:1000, sp1 = noise[[1]], sp2 = noise[[2]])
matrix = cbind(noise[[1]], noise[[2]])
matrix2 = sapply(noise, FUN = cbind)

df %>%
  filter(time %in% 1:100) %>%
  ggplot(aes(x = time, y = sp1)) + 
  geom_line() +
  geom_line(aes(y = sp2), colour = "red")

## calculate statistics 
meancorr(data = matrix2, nrands = 999, alternative = "two.tailed", type = 1, quiet = T) 
## 

phase <- phase.sync(df$sp1, df$sp2, nrands = 999, quiet = TRUE) 
## 0 phase difference at any time step 

## compute concurrency (whether time series peaks and troughs at the same time)
peaks <- peaks(df$sp1, df$sp2, nrands = 999, type = 1, quiet = TRUE)



df = data.frame(time = 1:1000, sp1 = rep(c(1,-1), length.out = 1000), sp2 = rep(c(-1, 1), length.out = 1000))
peaks(df$sp1, df$sp2, nrands = 999, type = 1, quiet = TRUE)

df = data.frame(time = 1:1000, sp1 = c(rep(c(-0.5,0.5), length.out = 1000)[1:900], rep(0, length.out = 100)), 
                sp2 = rep(c(-1, 1), length.out = 1000))
peaks(df$sp1, df$sp2, nrands = 999, type = 1, quiet = TRUE)
## measure of concurrency does not take into account the amount of difference between values of the time series

df = data.frame(time = 1:1000, sp1 = rep(c(0.1, 0.2), length.out = 1000), 
                sp2 = rep(c(-1, 1), length.out = 1000))
peaks(df$sp1, df$sp2, nrands = 999, type = 1, quiet = TRUE)
