## simulate species range shifts under climate change 
simulate_range_shifts <- function(stable_ranges,
                                  ps = c(0, 1), # synchrony parameters
                                  betas = c(0, 1), # autocorrelation parameters (spectral exponents)
                                  r = 1.2, # maximum intrinsic rate of increase 
                                  K = 100, # mean carrying capacity 
                                  d = 0.1, # proportion of offspring dispersing
                                  icp = 0.1, ## intraspecific competition parameter
                                  L = 1500, # number of time steps to simulate
                                  path = "outputs/data-processed/range-shifts" # set path
) {
  
  print("Starting!")
  
  ## get number of rows and columns in species range matrix 
  nrow = nrow(stable_ranges)
  ncol = ncol(stable_ranges)
  
  ## get number of ranges
  reps = 1
  
  ## read in function to generate time series of noise:
  source("R/functions/generate_noise.R")
  
  ######################################################
  ###      set additional simulation parameters       ## 
  ######################################################
  
  ## parameters that specify shape of latitudinal variation in intrinsic rate of increase (r)
  ## shape is sigmoidal, as was done in Mustin et al. 2013
  h = 1 ## half-saturation constant; defines the distance at which Eit = 0.5
  s = -3 ## shape parameter; defines direction (negative = negative slope) and shape (s > 1 gives signmoid)
  Emax = 1 ## set Emax to 1
  
  ## count number of populations:
  n_pops = nrow*ncol # number of grid cells
  
  #################################
  ###      run simulations       ## 
  #################################
  
  rep = 10
  
  ## loop through amounts of autocorrelation and synchrony
  p_curr = 1
  for(p in ps) {
    
    beta_curr = 1
    for(beta in betas) {
      print(paste0("On replicate ", rep, " with parameters beta = ", beta, " p = ", p))
      
      ## cellular lattice of microcosms
      lattice_N_it = array(dim = c(nrow,ncol,L))
      
      ## higher proportion dispersing = less pronounced effect of suitability gradient
      lattice_N_it[1:nrow,1:ncol,1] <- K/2 ## start with population size = carrying capacity / 2
      
      lattice_E_it = matrix(ncol = ncol, nrow = nrow)
      
      ## position optimum climatic conditions as row 25 on the lattice (Emax)
      opt = 25
      lattice_E_it[opt,] = Emax
      
      ## assume that conditions decline sigmoidally away from this optimum in both directions
      lattice_E_it[1:(opt-1),] = Emax*rev((1:(opt-1))^s/ ((1:(opt-1))^s + h^s))
      lattice_E_it[(opt+1):nrow,] = Emax*((1:nrow)^s / ((1:nrow)^s + h^s))[1:(nrow-opt)]
      #plot(x = 1:nrow, y = lattice_E_it[1:nrow,1])
      
      ## read in env raster
      filename =  paste0("outputs/data-processed/env-grids/shift_grid", rep, "_p", p, "_beta", beta, "_r", r, "_K", 100, "_d", 
                         d, "_icp", icp, "_L", L, "_reps", 100, ".tif")
      lattice_r_array = as.array(rast((filename)))
      
      ## get initial species range and set as range at t = 1
      lattice_N_it[,,1] = ceiling(stable_ranges/2)
      
      ## cut pop size in half
      
      ########################################################################
      ###    run population simulations under shifting climate gradient     ## 
      ########################################################################
      
      ## shift climate optimum by 0.33 cells per time step
      lattice_E_it_array <- replicate(1500, lattice_E_it) 
      #plot(x = 1:nrow, y = lattice_E_it_array[1:nrow,1,1])
      q = 1
      new_opt = 0
      while(q <= L & new_opt < nrow-1) {
        step = 0
        new_opt = opt + step
        
        ## shift optimum by "step"
        ## shift the optimum climatic conditions on the lattice
        lattice_E_it_array[new_opt,,q] = Emax
        ## assume that conditions decline sigmoidally away from this optimum in both directions
        lattice_E_it_array[1:((new_opt)-1),,q] = Emax*rev((1:((new_opt)-1))^s/ 
                                                            ((1:((new_opt)-1))^s + h^s))
        lattice_E_it_array[((new_opt)+1):nrow,,q] = Emax*((1:nrow)^s / 
                                                            ((1:nrow)^s + h^s))[1:(nrow-(new_opt))]
        # plot(x = 1:nrow, y = lattice_E_it_array[1:nrow,1,q])
        
        q = q + 1
        
        ## if the optimum has reached the edge of the sampling bounds
        if(new_opt == nrow-1) {
          ## end simulation and set the rest of the env. values to NA 
          lattice_E_it_array[,,q:L] = NA
        }
        
      }
      #plot(x = 1:nrow, y = lattice_E_it_array[1:nrow,1,100])
      
      ## run the population simulations
      t = 1 
      all_ranges <- c()
      while(t <= L) {
        
        ## save species range at time t as a data frame:
        matrix = lattice_N_it[,,t]
        df <- expand.grid(y = 1:nrow, x = 1:ncol) 
        pts <- as.data.frame(transform(df, z = matrix[as.matrix(df)]))
        pts$x = pts$x - 0.5
        pts$y = 100 - (pts$y - 0.5)
        pts = rename(pts, "Nt" = "z")
        pts$t = t
        pts$rep = rep
        pts$p = p
        pts$beta = beta
        pts$N_global = sum(lattice_N_it[,,t])
        pts$N_ext_local = length(which(lattice_N_it[,,t] == 0))
        
        all_ranges = rbind(all_ranges, pts)
        
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
              ## get number of dispersers 
              dispersers = curr_allpops[y,x]*d
              ## figure out where each one disperses 
              ## individuals can disperse into neighbouring 8 cells
              othercells <- expand.grid(x = c(x-1,x,x+1), y = c(y-1,y,y+1)) 
              othercells = othercells[!(x == othercells$x & y == othercells$y),]
              othercells$num = 1:8
              
              ## randomly sample from 1:8
              sample <- data.frame(num = sample(1:8, dispersers, replace = TRUE))
              
              othercells <- left_join(sample, othercells, by = "num") %>%
                filter(!(y == 0 | x == 0 | y > nrow | x > ncol))
              
              ## subtract dispersers
              curr_allpops_new[y,x] <- curr_allpops_new[y,x] - dispersers
              
              if(length(othercells) != 0) {
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
            
            new_sizes[y,x] = round(as.numeric(N*exp(r_curr[y,x]*(1 - as.complex(N/K)^icp))))
            
            ## add demographic stochasticity
            new_sizes[y,x] = sample(rpois(new_sizes[y,x], n = 1000), size = 1)
            
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
        
        t = t + 1
      }
      
      ## save data frame containing species range data at each time step:
      write.csv(all_ranges, paste0(path, 
                                   "/rep", rep, "_p", p, "_b", beta, "_icp", icp, 
                                   "_range-shifts_noshift.csv"), 
                row.names = FALSE)
      
      beta_curr = beta_curr + 1
    }
    p_curr = p_curr + 1
  }
  beta_curr = 1

}