## simulate range shift: 
## 1. generate a stable range for 500 time steps
## 2. begin shifting the suitability gradient
## this time: only simulate leading edge (assume trailing edge is out of simulation)

simulate_range_shifts <- function(p,
                                  beta,
                                  r = 1.2, # maximum intrinsic rate of increase 
                                  K = 100, # mean carrying capacity 
                                  d = 0.1, # proportion of offspring dispersing
                                  icp = 0.1, ## intraspecific competition parameter
                                  L = 2000, # number of time steps to simulate
                                  reps = 10, # number of replicates per combination of parameters
                                  nrow, # number of rows in species range matrix 
                                  ncol, # number of columns in species range matrix
                                  path = "outputs/data-processed/range-shift-simulations", # set path
                                  shift_rate = 0.1
) {
  
  r = 0.5 # maximum intrinsic rate of increase 
  K = 30 # mean carrying capacity 
  d = 0.1 # proportion of offspring dispersing
  icp = 0.1 ## intraspecific competition parameter
  L = 2000 # number of time steps to simulate
  reps = 10 # number of replicates per combination of parameters
  nrow=300 # number of rows in species range matrix 
  ncol=10 # number of columns in species range matrix
  path = "outputs/data-processed/range-shift-simulations" # set path
  shift_rate = 0.2
  
  print("Starting!")
  
  ## read in function to generate time series of noise:
  source("R/functions/generate_noise.R")
  
  ## create folder for output 
  path = paste0(path, "/p", p, "_b", beta, "_icp", icp, "_d", d)
  
  if(!dir.exists(path)) {
    dir.create(path, recursive = T)
  }
  
  ######################################################
  ###      set additional simulation parameters       ## 
  ######################################################
  
  ## count number of populations:
  n_pops = nrow*ncol # number of grid cells
  
  #################################
  ###      run simulations       ## 
  #################################
  
  ## for each replicate
  foreach(rep = 1:reps) %dopar% {
    
    print(paste0("On replicate ", rep, " with parameters beta = ", beta, " p = ", p))
    
    new_path = paste0(path, "/rep", rep, "_leading-edge_FAST")
    
    ## create subfolder 
    if(!dir.exists(new_path)) {
      dir.create(new_path, recursive = T)
    }
    
    t = 1 
    while(t <= L) {
      
      if(t == 1) {
        
        ## first, make an array with temperature that shifts at a given rate
        ## then, use thermal performance curve to translate to suitability 
        
        ## cellular lattice of microcosms
        lattice_N_it = array(dim = c(nrow,ncol,L))
        lattice_r = matrix(ncol = ncol, nrow = nrow)
        lattice_E_it = matrix(ncol = ncol, nrow = nrow)
        
        ## higher proportion dispersing = less pronounced effect of suitability gradient
        lattice_r[1:nrow,1:ncol] <- r ## start with growth rate = max growth rate 
        lattice_N_it[(nrow - 50 + 1):nrow,1:ncol,1] <- K/2 ## start with population size = carrying capacity / 2 in first 50 rows of grid
        lattice_N_it[1:(nrow - 50),1:ncol,1] <- 0 ## start with population size = 0 in other half 
        
        ## create temperature gradient 
        ## position optimum temperature (20deg) conditions as row 25 on the lattice (Emax)
        spat_grad = 0.2 #2.5*shift_rate ## spat grad
        temp_trend = spat_grad*shift_rate ## temporal trend 
        clim_velo = temp_trend/spat_grad
      
        lattice_E_it[,1:ncol] = 20 + (1:nrow - (nrow - 25)) * spat_grad
        
        ## replicate latitudinal gradient L times 
        lattice_E_it_array <- replicate(L, lattice_E_it)
        # plot(x = 1:nrow, y = lattice_E_it_array[1:nrow,1,1])
        
        ## after stable period of 500 time steps, shift temperatures by shift_rate 
        for(time in 501:L) {
          lattice_E_it_array[,1:ncol,time] = 20 + (1:nrow - (nrow - 25)) * spat_grad + temp_trend*(time-500)
        }
        #plot(x = 1:nrow, y = lattice_E_it_array[1:nrow,1,1501])
        
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
        
        ## make tpc: (from Huey suboptimal is optimal)
        ## plot a curve!
        # r_max = r + 0.25 ## max suitability
        # alpha_tpc = 0.3 ## rise rate steepness
        # beta_tpc = 0.3 ## decline rate steepness
        # Trmax = 20
        
        r_max = r + 0.25 ## max suitability
        alpha_tpc = 0.1 ## rise rate steepness
        beta_tpc = 0.1 ## decline rate steepness
        Trmax = 20
        
        # plot the curve
        # range of body temps
        Tb = seq(from = 0, to = 40, by = 0.1)

        ## calculate fitness
        exponent = -exp(beta_tpc*(Tb-Trmax)-8)-alpha_tpc*(Tb-Trmax)^2
        wb = r_max*exp(1)^exponent - 0.25

        # data.frame(r = wb, temperature = Tb) %>%
        #  ggplot(aes(x = temperature, y = r)) +
        #  geom_line() +
        #   geom_hline(yintercept = 0, colour = "red") +
        #   theme_bw()
        
        ## make it a leading edge 
        wb[Tb > Trmax] = r_max - 0.25
        
        data.frame(r = wb, temperature = Tb) %>%
          ggplot(aes(x = temperature, y = r)) +
          geom_line() +
          geom_hline(yintercept = 0, colour = "red") +
          theme_bw()

        # ## let noise affect r for each cell in the lattice
        # lattice_noise_array = (lattice_E_it_array + lattice_ac_it*3) ## scale sd to 3
        # #plot(x = 1:nrow, y = lattice_noise_array[1:nrow,1,500])
        # #plot(x = 1:L, y = lattice_noise_array[75,1,1:L])
        # 
        # # max(lattice_noise_array[75,1,1:500])
        # # min(lattice_noise_array[75,1,1:500])
        # 
        # ## use tpc to translate temps into suitability
        # lattice_r_array = lattice_noise_array
        # lattice_r_array = -exp(beta_tpc*(lattice_noise_array-Trmax)-8)-alpha_tpc*(lattice_noise_array-Trmax)^2
        # lattice_r_array = r_max*exp(1)^lattice_r_array - 0.5
        # #plot(x = 1:nrow, y = lattice_r_array[1:nrow,1, 700])
        # #plot(x = 1:L, y = lattice_r_array[50,1,1:L])
        
        # ## make it a leading edge 
        # lattice_r_array[lattice_noise_array > Trmax] = r_max - 0.25
        # #plot(x = 1:nrow, y = lattice_r_array[1:nrow,1, 700])
        
        # plot(x = 1000:1600, y = lattice_r_array[75,1,1000:1600])
        # plot(x = 200:300, y = lattice_r_array[200:300,1,150])
        
        ## use tpc to translate temperature to growth rate
        lattice_r_array = lattice_E_it_array
        lattice_r_array = -exp(beta_tpc*(lattice_r_array-Trmax)-8)-alpha_tpc*(lattice_r_array-Trmax)^2
        lattice_r_array = r_max*exp(1)^lattice_r_array - 0.25
        #plot(x = 1:nrow, y = lattice_r_array[1:nrow,1,500])
        
        ## make it a leading edge 
        lattice_r_array[lattice_E_it_array > Trmax] = r_max - 0.25
        #plot(x = 1:nrow, y = lattice_r_array[1:nrow,1, 700])

        ## add noise to growth rate
        lattice_r_array = lattice_r_array + lattice_ac_it*0.25
        #plot(x = 1:nrow, y = lattice_r_array[1:nrow,1,500])

        ## save environmental array as raster
        filename =  paste0("outputs/data-processed/env-grids/range-shift-grid", rep, "_p", p, "_beta", beta, "_r", r, "_K", K, "_d", 
                           d, "_icp", icp, "_L", L, ".tif")
        writeRaster(rast(lattice_r_array), filename, overwrite = TRUE)
        
        ## measure clim velocity after variation was added
        rast = raster::brick(lattice_r_array[,,501:L], xmn = 0, xmx = 10, ymn = 0, ymx = 300)
        ttrend = tempTrend(r = rast, th = 0.25*nlayers(rast)) ## set minimum # obs. to 1/4 time series length
        #plot(ttrend)
        temp_trend_est = mean(values(ttrend$slpTrends))
        
        ## see how much variation there is in temp trend 
        
        ## use function to calculate spatial gradient in mean daily air temperature for each location:
        spgrad = spatGrad(r = rast, projected = TRUE) ## our raster is projected to a coordinate system
        #plot(spgrad)
        spat_grad_est = mean(values(spgrad$Grad))
        
        ## calculate gradient based climate velocity:
        gvocc = gVoCC(tempTrend = ttrend, spatGrad = spgrad)
        #plot(gvocc)
        clim_velo_est = mean(values(gvocc$voccMag))
        
        ## get latitudinal velocity component 
        # cos_angle = cos(gvocc$voccAng)
        # lat_velo = gvocc$voccMag*cos_angle
        #plot(lat_velo)
        #plot(lat_velo < 0)
        
        ## okay - there are negative local velocities, just like in Berkeley Earth 
        #sd(values(lat_velo))
        
        ## get rid of unnecessary objects 
        rm("noise", "diffs","lattice_ac_it", "ts_all", "lattice_E_it_array")
      }
      
      ######################################
      ###    run population simulations   ## 
      ######################################
      ## run the population simulations
      ## save species range at time t as a data frame:
      matrix = lattice_N_it[,,t]
      df <- expand.grid(y = 1:nrow, x = 1:ncol) 
      pts <- as.data.frame(transform(df, z = matrix[as.matrix(df)]))
      pts$x = pts$x - 0.5
      pts$y = nrow - (pts$y - 0.5)
      pts = rename(pts, "Nt" = "z")
      pts$t = t
      pts$rep = rep
      pts$p = p
      pts$beta = beta
      pts$shift_rate = shift_rate
      pts$temp_trend_est = temp_trend_est
      pts$spat_grad_est = spat_grad_est
      pts$clim_velo_est = clim_velo_est
      pts$N_global = sum(lattice_N_it[,,t])
      pts$N_ext_local = length(which(lattice_N_it[,,t] == 0))
      
      write.csv(pts, paste0(new_path, "/p", p, "_b", beta, "_icp", icp, "_d", d, "_t", t,
                            "_range-shift_rep", rep, ".csv"),
                row.names = FALSE)
      
      ## get grid of growth rates at time t
      lattice_r_curr <- lattice_r_array[,,t]
      
      ## get grid of current population size at time t
      lattice_N_curr <- lattice_N_it[,,t]
      
      ##### DISPERSE #####
      ## fixed proportion of individuals have same probability of dispersing into any patch
      curr_allpops <- lattice_N_curr
      curr_allpops_new <- curr_allpops
      ## loop through grid cells
      x = 1
      while (x <= ncol) {
        y = 1
        while(y <= nrow) {
          ## if all pops are dead, stop
          if(length(which(curr_allpops == 0)) == n_pops) {
            y = n_pops + 1
            i = L
          }
          else {
            # ## get number of dispersers 
            # dispersers = curr_allpops[y,x]*d
            # ## figure out where each one disperses 
            # ## individuals can disperse into neighbouring 8 cells
            # othercells <- expand.grid(x = c(x-1,x,x+1), y = c(y-1,y,y+1)) 
            # othercells = othercells[!(x == othercells$x & y == othercells$y),]
            # othercells$num = 1:8
            # 
            # ## randomly sample from 1:8
            # sample <- data.frame(num = sample(1:8, dispersers, replace = TRUE))
            
            ## get number of dispersers 
            dispersers = round(curr_allpops[y,x]*d, 0)
            
            ## figure out where each one disperses 
            # ## individuals can disperse into neighbouring 24 cells
            # othercells <- expand.grid(x = c(x-2,x-1,x,x+1,x+2), y = c(y-2,y-1,y,y+1,y+2)) 
            # othercells = othercells[!(x == othercells$x & y == othercells$y),]
            # othercells$num = 1:24
            # 
            # ## randomly sample from 1:24
            # sample <- data.frame(num = sample(1:24, dispersers, replace = TRUE))
            
            ### LONG-DISTANCE DISPERSAL
            ## individuals can disperse into any of the neighbouring 230 cells
            othercells <- expand.grid(x = c(seq(from = x - 5, to = x-1), x, seq(from = x+1, to = x+5)), 
                                      y = c(seq(from = y - 10, to = y-1), y, seq(from = y+1, to = y+10)))
            othercells = othercells[!(x == othercells$x & y == othercells$y),]
            othercells$num = 1:230

            ## randomly sample 
            sample <- data.frame(num = sample(1:230, dispersers, replace = TRUE))
            
            othercells <- left_join(sample, othercells, by = "num") %>%
              filter(!(y <= 0 | x <= 0 | y > nrow | x > ncol))
            
            ## subtract dispersers
            curr_allpops_new[y,x] <- curr_allpops_new[y,x] - dispersers
            
            if(nrow(othercells) != 0) {
              ## subtract dispersers from population size, add them to population size at their new home   
              v = 1
              while(v <= nrow(othercells)) {
                curr_allpops_new[othercells$y[v], othercells$x[v]] = 
                  curr_allpops_new[othercells$y[v], othercells$x[v]] + 1
                v = v + 1
              }
            }
            y = y + 1
          }
        }
        x = x + 1
      }
      
      
      ##### REPRODUCE #####
      # for each subpopulation, calculate the new population size after reproduction:
      new_sizes <- matrix(ncol = ncol, nrow = nrow)
      r_curr = lattice_r_array[,,t]
      x = 1
      while (x <= ncol) {
        y = 1
        while(y <= nrow) {
          N = curr_allpops_new[y,x]

          ## when population size is above carrying capacity, density feedback becomes negative
          ## if growth rate is also negative, population will grow
          ## solution: change the model to apply monotonic negative feedback
          ## if growth rate r is negative, flip sign of density dependent term so that when r < 0, population declines even when Nt > K
          effective_r = ifelse(r_curr[y,x] >= 0, (1 - as.complex(N/K)^icp),
                               abs(1 - as.complex(N/K)^icp))

          new_size = round(as.numeric(N*exp(r_curr[y,x]*effective_r)))

          if(is.na(new_size) || new_size < 0 || is.infinite(new_size)) {
            new_size = 0
          }

          new_sizes[y,x] = new_size

          ## add demographic stochasticity by sampling # offspring from a poisson distribution where mean depends on N, r, icp, and K
          new_sizes[y,x] = sample(rpois(new_sizes[y,x], n = 1000), size = 1)
          
          # N = curr_allpops_new[y,x]
          # 
          # new_sizes[y,x] = round(as.numeric(N*exp(r_curr[y,x]*(1 - as.complex(N/K)^icp))))
          # 
          # ## add demographic stochasticity
          # new_sizes[y,x] = sample(rpois(new_sizes[y,x], n = 1000), size = 1)
           
          y = y + 1
        }
        x = x + 1
      }
      ## if pop size < 0, set to 0
      new_sizes[which(new_sizes < 0)] <- 0
      new_sizes[which(is.na(new_sizes))] <- 0
      
      ## update population matrix:
      if(t != L) {
        lattice_N_it[,,t+1] = new_sizes
      }
      
      ## if global population size when t == 500 is less than 500 individuals, restart from the beginning
      if(t == 500 && sum(lattice_N_it[,,t]) < 500) {
        t = 1
        N_global <- c()
        N_ext_local <- 0
      }
      else {
        t = t + 1
      }
    }
    
  t
  }
 }



