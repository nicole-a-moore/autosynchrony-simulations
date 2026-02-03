library(terra)

r = rast("outputs/data-processed/env-grids/range-shift-grid_both-edges_rep1_p0_b0_icp0.7_K200_d0.2_r1_d-dist1_sigma0.5_shift-rate0.2.tif")

plot(r[[300]])

lattice_r_array = as.array(r)

test = rast(arr)


t = unlist(c(r[150,1,]))

df =data.frame(t = t, x = 1:2000)

df %>%
  ggplot(aes(x = x, y = t)) +
  geom_line()

mean(t)

hist(t[1:1000]-0.25)

## calculate position of r edge over time 
r = rast("outputs/data-processed/env-grids/range-shift-grid1_p1_beta0_r2_K100_d0.1_icp0.1_L2000.tif")

## convert to xy data frame
pred = as.data.frame(r, xy = TRUE)
p95 = c()
for(i in 3:(ncol(pred))) {
  cur = as.data.frame(pred[,c(1:2, i)])
  cur = cur[which(!is.na(cur[,3])),]
  p95 = append(p95, modi::weighted.quantile(cur$y, w = cur[,3], 0.5))
}

#plot(x = 1:2000, y = p95)

data.frame(x = 1:2000, y = p95) %>%
  filter(x > 500) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line()


r = rast("outputs/data-processed/env-grids/range-shift-grid1_p1_beta0_r2_K100_d0.1_icp0.1_L2000.tif")

plot(r[[1]])

t = unlist(c(r[300,1,]))

df =data.frame(t = t, x = 1:2000)

df %>%
  ggplot(aes(x = x, y = t)) +
  geom_line()

mean(t)


##### calculate mean shift in edge of suitability over 1500 time steps where there is shift 
## see if it matches shift rate of 0.2 cells per time step

r = rast("outputs/data-processed/env-grids/range-shift-grid1_p1_beta2.5_r2_K100_d0.1_icp0.1_L2000.tif")


plot(r[[600]])


ts = r[200,1,1:2000]

df = data.frame(x = 1:2000, y = as.numeric(ts))

df %>%
  ggplot(aes(x = x, y = y)) +
  geom_line()


p75 = c()
mean = c()
for(i in 500:nlyr(r)) {
  df = as.data.frame(r[[i]], xy = T) 
  
  p75 = append(p75, modi::weighted.quantile(df$y, w = df[,3], 0.05))
  mean = append(mean, weighted.mean(df$y, w = df[,3]))
  
  print(i)
}

df = data.frame(t = 500:2000, p75 = 300 - p75, mean = mean)

df %>%
  ggplot(aes(x = t, y = p75)) +
  geom_point()


df


df %>%
  lm(data=., x~t)




r = rast("outputs/data-processed/env-grids/range-shift-grid1_p0_beta2.5_r2_K100_d0.1_icp0.1_L2000.tif")

plot(r[[900]])

t = unlist(c(r[150,1,]))

df =data.frame(t = t, x = 1:2000)

df %>%
  ggplot(aes(x = x, y = t)) +
  geom_line()

mean(t)

df$t = df$t*0.25

t = 1250
N = 9
K=100
new_sizes = c()
while(t <= 2000) {
  
  if(t == 1250) {
    N = N
  }
  else {
    N = new_sizes[t - 1250]
  }
 
  ## when population size is above carrying capacity, density feedback becomes negative
  ## if growth rate is also negative, population will grow
  ## solution: change the model to apply monotonic negative feedback
  ## if growth rate r is negative, flip sign of density dependent term so that when r < 0, population declines even when Nt > K
  effective_r = ifelse(df$t[t] >= 0, (1 - as.complex(N/K)^icp),
                       abs(1 - as.complex(N/K)^icp))
  ## add allee effect
  A = 10
  effective_r = effective_r*(N/A-1)
  
  new_size = round(as.numeric(N*exp(df$t[t]*effective_r)))
  
  if(is.na(new_size) || new_size < 0 || is.infinite(new_size)) {
    new_size = 0
  }
  
  ## add demographic stochasticity by sampling # offspring from a poisson distribution where mean depends on N, r, icp, and K
  new_size = sample(rpois(new_size, n = 1000), size = 1)
  
  new_sizes = append(new_sizes,new_size)
  
  t = t + 1
}

plot(new_sizes)




## try something

ts_1 = df

ts_0 = df 

ts_0 %>%
  ggplot(aes(x = x, y = t)) +
  geom_line(data = ts_1,colour = 'red') +
  geom_line() 

ts_0 %>%
  filter(x>500, x<550) %>%
  ggplot(aes(x = x, y = t)) +
  geom_line(data = ts_1 %>% filter(x>500, x<550),colour = 'red') +
  geom_line() 


ts_0_sub = ts_0$t[1:2000]
ts_1_sub = ts_1$t[1:2000]

t = 1
N0 = N1 = 1
new_sizes_0 = new_sizes_1 = c()
K=200
while(t < length(ts_1_sub)) {
  
  effective_r = abs(1 - as.complex(N/K)^icp)
  
  new_size_0 = floor(as.numeric(N0*exp(ts_0_sub[t]*effective_r)))
  new_size_1 = floor(as.numeric(N1*exp(ts_1_sub[t]*effective_r)))
  
  ## add demographic stochasticity by sampling # offspring from a poisson distribution where mean depends on N, r, icp, and K
  # new_size_0 = sample(rpois(new_size_0, n = 1000), size = 1)
  # new_size_1 = sample(rpois(new_size_1, n = 1000), size = 1)

  ## add dispersal
  # new_size_0 = new_size_0 + 1
  # new_size_1 = new_size_1 + 1
  
  new_sizes_0 = append(new_sizes_0, new_size_0)
  new_sizes_1 = append(new_sizes_1, new_size_1)
  
  t = t + 1
}

new_sizes_0
new_sizes_1

min(ts_0$t)
max(ts_0$t)





