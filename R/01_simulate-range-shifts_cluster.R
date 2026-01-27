## simulate range shift: 
## 1. generate a stable range for 500 time steps
## 2. begin shifting the suitability gradient
simulate_range_shifts <- function(p,
                                  beta,
                                  r, # maximum intrinsic rate of increase 
                                  K, # mean carrying capacity 
                                  d, # proportion of offspring dispersing
                                  d_dist,
                                  icp, ## intraspecific competition parameter
                                  L, # number of time steps to simulate
                                  nrow, # number of rows in species range matrix 
                                  ncol, # number of columns in species range matrix
                                  path, # set path
                                  sigma,
                                  shift_rate ## set shift rate
) {
  
  rep = 1
  
  print("Starting!")
  
  ## read in function to generate time series of noise:
  source("functions/generate_noise.R")
  
  ## create folder for output 
  if(!dir.exists(path)) {
    dir.create(path, recursive = T)
  }
  
  ######################################################
  ###      set additional simulation parameters       ## 
  ######################################################
  
  ## parameters that specify shape of latitudinal variation in intrinsic rate of increase (r)
  ## shape is sigmoidal, as was done in Mustin et al. 2013
  h = 20 ## half-saturation constant; defines the distance at which Eit = 0.5
  s = -3 ## shape parameter; defines direction (negative = negative slope) and shape (s > 1 gives signmoid)
  Emax = r ## set Emax to max r
  
  ## count number of populations:
  n_pops = nrow*ncol # number of grid cells
  
  #################################
  ###      run simulations       ## 
  #################################
  
  print(paste0("On replicate ", rep, " with parameters p", p, "_b", beta, "_icp", icp, "_K", K, "_d", d, "_r", r, "_d-dist", 
               d_dist, "_sigma", sigma, "_shift-rate", shift_rate))
  
  new_path = paste0(path, "/rep", rep)
  
  ## create subfolder 
  if(!dir.exists(new_path)) {
    dir.create(new_path, recursive = T)
  }
  
  t = 1 
  while(t <= L) {
    
    if(t == 1) {
      
      ## cellular lattice of microcosms
      lattice_N_it = array(dim = c(nrow,ncol,L))
      lattice_r = matrix(ncol = ncol, nrow = nrow)
      lattice_E_it = matrix(ncol = ncol, nrow = nrow)
      
      ## higher proportion dispersing = less pronounced effect of suitability gradient
      lattice_r[1:nrow,1:ncol] <- r ## start with growth rate = max growth rate 
      lattice_N_it[1:(nrow/3),1:ncol,1] <- K/2 ## start with population size = carrying capacity / 2 in 1/3 of grid
      lattice_N_it[(nrow/3):nrow,1:ncol,1] <- 0 ## start with population size = 0 in other half 
      
      ## position optimum climatic conditions as row 1 on the lattice (Emax)
      opt = 1
      lattice_E_it[opt,1:10] = Emax + 0.05
      
      # ## assume that conditions decline sigmoidally away from this optimum in both directions
      # lattice_E_it[1:(opt-1),] = Emax*2*rev((1:(opt-1))^s/ ((1:(opt-1))^s + h^s)) 
      # lattice_E_it[(opt+1):nrow,] = Emax*2*((1:nrow)^s / ((1:nrow)^s + h^s))[1:(nrow-opt)] 
      # #plot(x = 1:nrow, y = lattice_E_it[1:nrow,1])
      
      ## assume that conditions decline sigmoidally away from this optimum in one direction
      lattice_E_it[1:opt,] = Emax + 0.05
      lattice_E_it[(opt+1):nrow,] = (Emax + 0.05)*((1:nrow)^s / ((1:nrow)^s + h^s))[1:(nrow-opt)] 
      #plot(x = 1:nrow, y = lattice_E_it[1:nrow,1])
      
      ## make growth rate plateau at -0.05
      lattice_E_it = lattice_E_it - 0.05
      #plot(x = 1:nrow, y = lattice_E_it[1:nrow,1])
      
      ## get distance where r = 0 
      first(which(lattice_E_it[,1] <= 0))
      ## row 41
      
      ## replicate latitudinal gradient L times 
      lattice_E_it_array <- replicate(L, lattice_E_it)
      
      ## after stable period of 500 time steps, shift climate optimum by shift_rate 
      #plot(x = 1:nrow, y = lattice_E_it_array[1:nrow,1,1])
      q = 1
      new_opt = 0
      while(q <= (L-500)) {
        ## if optimum is still on grid
        if(new_opt < nrow-1) {
          step = floor(shift_rate*q)
          new_opt = opt + step
          
          ## shift optimum by "step"
          ## shift the optimum climatic conditions on the lattice
          ## assume that conditions decline sigmoidally away from this optimum in one directions
          lattice_E_it_array[1:new_opt,,(q+500)] = Emax + 0.05
          lattice_E_it_array[(new_opt+1):nrow,,(q+500)] = (Emax + 0.05)*((1:nrow)^s /
                                                                           ((1:nrow)^s + h^s))[1:(nrow-(new_opt))]
          # plot(x = 1:nrow, y = lattice_E_it_array[1:nrow,1,q])
          ## make growth rate plateau at -0.05
          lattice_E_it_array[1:nrow,1:ncol,q+500] = lattice_E_it_array[1:nrow,1:ncol,q+500] - 0.05
        }
        ## otherwise
        else {
          step = floor(shift_rate*q)
          new_opt = opt + step
          
          ## shift optimum by "step"
          # ## assume that conditions decline sigmoidally away from this optimum in both directions
          # lattice_E_it_array[1:nrow,,(q+500)] = (Emax*2*rev((1:((new_opt)))^s/ 
          #                                                     ((1:((new_opt)))^s + h^s)))[1:nrow]
          ## assume that conditions decline sigmoidally away from this optimum in one directions
          lattice_E_it_array[1:nrow,,(q+500)] = Emax
          # plot(x = 1:nrow, y = lattice_E_it_array[1:nrow,1,(q+500)])
        }
        
        q = q + 1
      }
      #plot(x = 1:nrow, y = lattice_E_it_array[1:nrow,1,1000])
      
      ## create noise time series for each cell, L time steps long
      lattice_ac_it = array(dim = c(nrow,ncol,L))
      
      # generate noise of given synchrony and autocorrelation
      noise <- generate_noise(beta = beta,
                              p = p,
                              n_ts = nrow*ncol, ## number of cells
                              L1 = 500,
                              L2 = L - 500,
                              sigma = sigma)
      
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
      lattice_r_array = (lattice_E_it_array + lattice_ac_it)
      #plot(x = 1:nrow, y = lattice_r_array[1:nrow,1,500])
      #plot(x = 1:2000, y = lattice_r_array[150,1,1:2000])
      
      ## save environmental array as raster
      filename =  paste0("output/env-grids/range-shift-grid_rep", rep, "_p", p, "_b", beta, "_icp", icp, "_K", K, "_d", 
                         d, "_r", r, "_d-dist", d_dist, "_sigma", sigma, "_shift-rate", shift_rate,  ".tif")
      writeRaster(rast(lattice_r_array), filename, overwrite = TRUE)
      
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
    pts$y = pts$y - 0.5
    pts = rename(pts, "Nt" = "z")
    pts$t = t
    pts$rep = rep
    pts$p = p
    pts$beta = beta
    pts$shift_rate = shift_rate
    pts$N_global = sum(lattice_N_it[,,t])
    pts$N_ext_local = length(which(lattice_N_it[,,t] == 0))
    
    write.csv(pts, paste0(new_path, "/range-shift_t", t, ".csv"), 
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
          ## get number of dispersers 
          dispersers = curr_allpops[y,x]*d
          
          ## calculate number of cells individuals can dispersal to
          ncell_disp = (d_dist*2 + 1)^2 - 1
          
          ## find their cell numbers 
          othercells <- expand.grid(x = c(seq(from = x - d_dist, to = x - 1), x, seq(from = x + 1, to = x + d_dist)),
                                    y = c(seq(from = y - d_dist, to = y - 1), y, seq(from = y + 1, to = y + d_dist)))
          othercells = othercells[!(x == othercells$x & y == othercells$y),]
          othercells$num = 1:ncell_disp
          
          ## randomly sample
          sample <- data.frame(num = sample(1:ncell_disp, dispersers, replace = TRUE))
          
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
        effective_r = r_curr[y,x]*ifelse(r_curr[y,x] >= 0, (1 - as.complex(N/K)^icp),
                                         abs(1 - as.complex(N/K)^icp))
        
        new_size = floor(as.numeric(N*exp(effective_r)))
        
        if(is.na(new_size) || new_size < 0 || is.infinite(new_size)) {
          new_size = 0
        }
        
        new_sizes[y,x] = new_size
        
        ## add demographic stochasticity by sampling # offspring from a poisson distribution where mean depends on N, r, icp, and K
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
    
    ## if global population size when t == 500 is less than 500 individuals, restart from the beginning
    if(t == 500 && sum(lattice_N_it[,,t]) < 500) {
      t = 1
      N_globl <- c()
      N_ext_local <- 0
    }
    else {
      t = t + 1
    }
  }
}
