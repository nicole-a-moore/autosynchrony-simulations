## function to generate initial species ranges under stable climate conditions 
generate_stable_ranges <- function(p, # synchrony parameter 
                                  beta, # autocorrelation parameter (spectral exponent)
                                  r = 1.2, # maximum intrinsic rate of increase 
                                  K = 100, # mean carrying capacity 
                                  d = 0.1, # proportion of offspring dispersing
                                  icp = 0.1, ## intraspecific competition parameter
                                  L = 500, # number of time steps to simulate
                                  nrow, # number of rows in species range matrix 
                                  ncol, # number of columns in species range matrix
                                  reps = 100, # number of ranges to generate 
                                  path = "outputs/data-processed/stable-ranges" # set path
) {
  
  ## read in function to generate time series of noise:
  source("R/functions/generate_noise.R")
  
  ## create folder for output 
  path = paste0(path, "/", "rep", rep, "_p", p, "_b", beta, "_icp", icp, "_d", d, "/")
  
  if(!dir.exists(path)) {
    dir.create(path, recursive = T)
  }
  
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
  
  ## loop through replicates 
  foreach (rep = 1:reps)  %dopar% {
 
    ########################################################################
    ###      run population simulations under stable climate conditions   ## 
    ########################################################################
    t = 1 
    N_global <- c()
    N_ext_local <- 0
    while(t < L) {
      
      if(t == 1) {
        
        ## cellular lattice of microcosms
        lattice_N_it = matrix(ncol = ncol, nrow = nrow)
        lattice_r = matrix(ncol = ncol, nrow = nrow)
        lattice_E_it = matrix(ncol = ncol, nrow = nrow)
        
        ## higher proportion dispersing = less pronounced effect of suitability gradient
        lattice_r[1:nrow,1:ncol] <- r ## start with growth rate = max growth rate 
        lattice_N_it[1:nrow,1:ncol] <- K/2 ## start with population size = carrying capacity / 2
        
        ## position optimum climatic conditions as row 25 on the lattice (Emax)
        opt = 25
        lattice_E_it[opt,] = Emax
        
        ## assume that conditions decline sigmoidally away from this optimum in both directions
        # lattice_E_it[1:(opt-1),] = Emax*rev((1:(opt-1))^s/ ((1:(opt-1))^s + h^s))
        # lattice_E_it[(opt+1):nrow,] = Emax*((1:nrow)^s / ((1:nrow)^s + h^s))[1:(nrow-opt)]
        #plot(x = 1:nrow, y = lattice_E_it[1:nrow,1])
        
        #### TEST
        ## try letting growth rate be negative instead of asympoting at 0 
        lattice_E_it[1:(opt-1),] = Emax*2*rev((1:(opt-1))^s/ ((1:(opt-1))^s + h^s)) - 0.25
        lattice_E_it[(opt+1):nrow,] = Emax*2*((1:nrow)^s / ((1:nrow)^s + h^s))[1:(nrow-opt)] - 0.25
        # plot(x = 1:nrow, y = lattice_E_it[1:nrow,1])
        
        ## create noise time series for each cell during stable conditions, L time steps long
        lattice_ac_it = array(dim = c(nrow,ncol,L))
        
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
        
        lattice_r_array = (replicate(L, lattice_E_it) + lattice_ac_it[,,1:L])*r
        
        ## check that mean r decreases
        # m = function (x){ mean(lattice_r_array[x,1:ncol,1:L]) }
        # means = sapply(1:60, FUN = m)
        # plot(x = 1:60, y = means)
        
        ## save environmental array as raster
        filename =  paste0("outputs/data-processed/env-grids/grid", rep, "_p", p, "_beta", beta, "_r", r, "_K", K, "_d", 
                           d, "_icp", icp, "_L", L, "_reps", reps, ".tif")
        writeRaster(rast(lattice_r_array), filename, overwrite = TRUE)
      }
      
      ## save data frame containing species range data at time step t:
      matrix = lattice_N_it[,]
      df <- expand.grid(y = 1:nrow, x = 1:ncol) 
      pts <- as.data.frame(transform(df, z = matrix[as.matrix(df)]))
      pts$x = pts$x - 0.5
      pts$y = 100 - (pts$y - 0.5)
      pts = rename(pts, "Nt" = "z")
      pts$t = t
      pts$rep = rep
      pts$p = p
      pts$beta = beta
      pts$shift_rate = 0
      pts$N_global = sum(lattice_N_it[,])
      pts$N_ext_local = length(which(lattice_N_it[,] == 0))
      
      write.csv(pts, paste0(path, "/rep", rep, "_p", p, "_b", beta, "_icp", icp, "_d", d, "_t", t, 
                                   "_stable-ranges.csv"), 
                row.names = FALSE)
      
      ## get grid of growth rates at time t
      lattice_r_curr <- lattice_E_it + lattice_ac_it[,,t]*r
      
      ## get grid of current population size at time t
      lattice_N_curr <- lattice_N_it[,]
      
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
      lattice_N_it = new_sizes
      
      ## save stats:
      ## calculate global population size (sum of all pop sizes)
      N_global <- append(N_global, sum(lattice_N_it[,]))
      ## calculate local population extinction frequency (how many pops go extinct per time step)
      N_ext_local = N_ext_local + length(which(new_sizes == 0))
      
      ## if global population size is less than 500 individuals, restart from the beginning
      if(t == (L-1) && sum(lattice_N_it[,]) < 500) {
        t = 1
        N_global <- c()
        N_ext_local <- 0
      }
      else {
        t = t + 1
      }
    }
  }
  
  return()
}