## simulate species range shifts under climate change 
simulate_range_shifts <- function(stable_ranges,
                                  ps = c(0, 1), # synchrony parameters
                                  betas = c(0, 1), # autocorrelation parameters (spectral exponents)
                                  r = 1.2, # maximum intrinsic rate of increase 
                                  K = 100, # mean carrying capacity 
                                  d = 0.1, # proportion of offspring dispersing
                                  icp = 0.1, ## intraspecific competition parameter
                                  L = 1500, # number of time steps to simulate
                                  path = "outputs/data-processed/range-shifts", # set path
                                  shift_rate = 100/1000
) {
  
  print("Starting!")
  
  ## get number of rows and columns in species range matrix 
  nrow = nrow(stable_ranges[[1]])
  ncol = ncol(stable_ranges[[1]])
  
  ## get number of ranges
  reps = length(stable_ranges)
  
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
  
  ## for each replicate
  all <- foreach(rep = 1:reps) %dopar% {
    
    ## loop through amounts of autocorrelation and synchrony
    p_curr = 1
    for(p in ps) {
      
      beta_curr = 1
      for(beta in betas) {
        print(paste0("On replicate ", rep, " with parameters beta = ", beta, " p = ", p))
        
        ## cellular lattice of microcosms
        lattice_N_it = array(dim = c(nrow,ncol,L))
        lattice_r = matrix(ncol = ncol, nrow = nrow)
        lattice_E_it = matrix(ncol = ncol, nrow = nrow)
        
        ## higher proportion dispersing = less pronounced effect of suitability gradient
        lattice_r[1:nrow,1:ncol] <- r ## start with growth rate = max growth rate 
        lattice_N_it[1:nrow,1:ncol,1] <- K/2 ## start with population size = carrying capacity / 2
        
        ## position optimum climatic conditions as row 25 on the lattice (Emax)
        opt = 25
        lattice_E_it[opt,] = Emax
        
        # ## assume that conditions decline sigmoidally away from this optimum in both directions
        # lattice_E_it[1:(opt-1),] = Emax*rev((1:(opt-1))^s/ ((1:(opt-1))^s + h^s))
        # lattice_E_it[(opt+1):nrow,] = Emax*((1:nrow)^s / ((1:nrow)^s + h^s))[1:(nrow-opt)]
        # #plot(x = 1:nrow, y = lattice_E_it[1:nrow,1])
        
        #### TEST
        ## try letting growth rate be negative instead of asymptoting at 0 
        lattice_E_it[1:(opt-1),] = Emax*2*rev((1:(opt-1))^s/ ((1:(opt-1))^s + h^s)) - 0.25
        lattice_E_it[(opt+1):nrow,] = Emax*2*((1:nrow)^s / ((1:nrow)^s + h^s))[1:(nrow-opt)] - 0.25
        # plot(x = 1:nrow, y = lattice_E_it[1:nrow,1])
        
        ## replicate latitudinal gradient L times 
        lattice_E_it_array <- replicate(L, lattice_E_it)
        
        ## create noise time series for each cell during stable conditions, L time steps long
        lattice_ac_it = array(dim = c(nrow,ncol,L))
        
        ## shift climate optimum by 0.33 cells per time step
        lattice_E_it_array <- replicate(1500, lattice_E_it) 
        #plot(x = 1:nrow, y = lattice_E_it_array[1:nrow,1,1])
        q = 1
        new_opt = 0
        while(q <= L) {
          ## if optimum is still on grid
          if(new_opt < nrow-1) {
            step = floor(shift_rate*q)
            new_opt = opt + step
            
            ## shift optimum by "step"
            ## shift the optimum climatic conditions on the lattice
            lattice_E_it_array[new_opt,,q] = Emax
            ## assume that conditions decline sigmoidally away from this optimum in both directions
            lattice_E_it_array[1:((new_opt)-1),,q] = Emax*2*rev((1:((new_opt)-1))^s/ 
                                                                  ((1:((new_opt)-1))^s + h^s)) - 0.25
            lattice_E_it_array[((new_opt)+1):nrow,,q] = Emax*2*((1:nrow)^s / 
                                                                  ((1:nrow)^s + h^s))[1:(nrow-(new_opt))] - 0.25
            # plot(x = 1:nrow, y = lattice_E_it_array[1:nrow,1,q])
          }
          else if(new_opt == (nrow - 1)) {
            step = floor(shift_rate*q)
            new_opt = opt + step

            ## shift the optimum climatic conditions on the lattice
            lattice_E_it_array[new_opt,,q] = Emax
            ## assume that conditions decline sigmoidally away from this optimum in both directions
            lattice_E_it_array[1:(new_opt - 1),,q] = (Emax*2*rev((1:((new_opt)-1))^s/ 
                                                          ((1:((new_opt)-1))^s + h^s)) - 0.25)
            # plot(x = 1:nrow, y = lattice_E_it_array[1:nrow,1,q])
          }
          ## otherwise
          else {
            step = floor(shift_rate*q)
            new_opt = opt + step
            
            ## shift optimum by "step"
            ## assume that conditions decline sigmoidally away from this optimum in both directions
            lattice_E_it_array[1:nrow,,q] = (Emax*2*rev((1:((new_opt)-1))^s/ 
                                                                  ((1:((new_opt)-1))^s + h^s)) - 0.25)[1:nrow]
            # plot(x = 1:nrow, y = lattice_E_it_array[1:nrow,1,q])
          }
          
          q = q + 1
        }
        #plot(x = 1:nrow, y = lattice_E_it_array[1:nrow,1,1499])
        
        ## generate noise of given synchrony and autocorrelation
        noise <- generate_noise(beta = beta,
                                p = p,
                                n_ts = nrow*ncol, ## number of cells
                                L = L)
        
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
        lattice_r_array = (lattice_E_it_array + lattice_ac_it[,,1:L])*r
        # plot(x = 1:nrow, y = lattice_r_array[1:nrow,1,1])
        
        ## save environmental array as raster
        filename =  paste0("outputs/data-processed/env-grids/shift_grid", rep, "_p", p, "_beta", beta, "_r", r, "_K", K, "_d", 
                           d, "_icp", icp, "_L", L, "_reps", reps, ".tif")
        writeRaster(rast(lattice_r_array), filename, overwrite = TRUE)
        
        ## remove unneeded objects from env
        rm("lattice_r", "lattice_E_it_array", "lattice_ac_it")
        
        ## get initial species range and set as range at t = 1
        lattice_N_it[,,1] = stable_ranges[[rep]]
        
        ########################################################################
        ###    run population simulations under shifting climate gradient     ## 
        ########################################################################
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
          pts$shift_rate = shift_rate
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
                                     "_range-shifts.csv"), 
                  row.names = FALSE)
        
        beta_curr = beta_curr + 1
      }
      p_curr = p_curr + 1
    }
    beta_curr = 1
}


}
