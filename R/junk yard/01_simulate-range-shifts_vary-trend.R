## simulate range shift: 
## 1. generate a stable range for 500 time steps
## 2. begin shifting the suitability gradient

simulate_range_shifts <- function(p,
                                  beta,
                                  icc,
                                  between_var,
                                  r = 1.2, # maximum intrinsic rate of increase 
                                  K = 100, # mean carrying capacity 
                                  d = 0.1, # proportion of offspring dispersing
                                  icp = 0.1, ## intraspecific competition parameter
                                  L = 2000, # number of time steps to simulate
                                  reps = 10, # number of replicates per combination of parameters
                                  nrow, # number of rows in species range matrix 
                                  ncol, # number of columns in species range matrix
                                  path = "outputs/data-processed/range-shift-simulations", # set path
                                  shift_rate = 0.05
) {
  
  library(tidyverse)
  library(raster)
  library(terra)
  select = dplyr::select
  
  r = 1.2 # maximum intrinsic rate of increase 
  K = 100 # mean carrying capacity 
  d = 0.1 # proportion of offspring dispersing
  icp = 0.1 ## intraspecific competition parameter
  L = 2000 # number of time steps to simulate
  reps = 10 # number of replicates per combination of parameters
  nrow=105 # number of rows in species range matrix 
  ncol=12 # number of columns in species range matrix
  path = "outputs/data-processed/range-shift-simulations" # set path
  shift_rate = 0.1
  between_var = 0
  icc = 1
  
  print("Starting!")
  
  ## read in function to generate time series of noise:
  source("R/functions/generate_noise.R")
  
  ## create folder for output 
  path = paste0(path, "/p", p, "_b", beta, "_icp", icp, "_d", d, "fastcc")
  
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
    
    new_path = paste0(path, "/wt_rep", rep)
    
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
        lattice_N_it[ceiling(nrow/2+1):nrow,1:ncol,1] <- K/2 ## start with population size = carrying capacity / 2 in 1/2 of grid
        lattice_N_it[1:ceiling(nrow/2),1:ncol,1] <- 0 ## start with population size = 0 in other half 
        
        ## create temperature gradient 
        ## position optimum temperature (20deg) conditions as row 85 on the lattice (Emax)
        spat_grad = 0.2 #2.5*shift_rate ## spat grad
        temp_trend = spat_grad*shift_rate ## temporal trend 
        clim_velo = temp_trend/spat_grad
      
        lattice_E_it[,1:ncol] = spat_grad*(1:nrow) + (20-spat_grad*85)
        
        ## replicate latitudinal gradient L times 
        lattice_E_it_array <- replicate(L, lattice_E_it)
        # plot(x = 1:nrow, y = lattice_E_it_array[1:nrow,1,1])
        
        ## generate grid of temporal warming trends with specific ICC
        ## ICC = 1 = full synchrony within neighbourhoods = all are warming at shift rate
        ## ICC = 0 = asynchrony = sample trends from randomly such that mean trend = shift rate
        ## measure using the intraclass correlation coefficient (ICC)
        ## to calculate ICC, group adjacent grid cells into groups of 9 neighbours
        ## calculate ICC to measure similarity within groups 
        
        ## function to generate data with specific ICC 
        generate_icc_data <- function(n_groups, n_per_group, icc, between_var, mu) {
          # Total number of observations
          N <- n_groups * n_per_group
          
          # Between-group variance
          sigma2_u <- between_var
          
          # Compute within-group variance to match target ICC
          sigma2_e <- sigma2_u * (1 - icc) / icc
          
          # Group-level random effects
          group_effects <- rnorm(n_groups, mean = 0, sd = sqrt(sigma2_u))
          
          # Assign group IDs
          group_ids <- rep(1:n_groups, each = n_per_group)
          
          # Individual-level noise
          epsilon <- rnorm(N, mean = 0, sd = sqrt(sigma2_e))
          
          # Simulated response
          y <- mu + group_effects[group_ids] + epsilon
          
          # Return data.frame
          data.frame(
            y = y,
            group = factor(group_ids)
          )
        }
        
        ## calculate how many neighbourhoods to generate
        ## convert to grid 
        rast = rast(lattice_r)
        values(rast) = 1:n_pops
        
        ## function to get 3x3 neighbourhood (center + 8 neighbors)
        groups <- lapply(1:ncell(rast), function(cell) {
          adj <- adjacent(rast, cells = cell, directions = 8, include = TRUE)
          return(c(adj))
        })
        
        
        ## create matrix of right size
        mat <- matrix(1:(nrow * ncol), nrow = nrow, byrow = TRUE)
        
        ## define stride (3x3 blocks)
        block_size <- 3
        
        ## initialize list to store center indices
        center_indices <- c()
        
        ## loop over non-overlapping 3x3 blocks
        for (row in seq(2, nrow - 1, by = block_size)) {
          for (col in seq(2, ncol - 1, by = block_size)) {
            ## get index of center cell
            center_idx <- mat[row, col]
            center_indices <- c(center_indices, center_idx)
          }
        }
        
        ## generate data with specific ICC 
        trends = generate_icc_data(n_groups = length(center_indices), n_per_group = 9, 
                                   icc = icc, ## within-neighbourhood similarity
                                   mu = 0, ## mean value
                                   between_var = between_var) ## between-neighbourhood variance
        
        ## set mean to equal temp_trend
        trends$y = trends$y - mean(trends$y, na.rm = TRUE) + temp_trend
        mean(trends$y)
        var(trends$y)
        
        ## make sure icc is right 
        ## pivot wider
        trends <- trends %>%
          group_by(group) %>%
          mutate(obs = row_number()) %>%
          ungroup() %>% 
          pivot_wider(
            id_cols = group,
            names_from = obs,
            values_from = y,
            names_prefix = "obs_"
          ) %>%
          select(-group)
        
        ## calculate icc 
        icc_star = irr::icc(trends, model = "oneway", type = "agreement", unit = "single")
        icc_star$value
        
        ## filter to keep only non-overlapping 3x3 neighborhoods (i.e., 9 cells) 
        groups <- groups[which(sapply(groups, function(x) { return(x[1])}) %in% center_indices)]
        groups = do.call(rbind, groups)
        
        ## make long dataframe 
        groups = as.data.frame(groups) %>%
          mutate(group = 1:nrow(.)) %>%
          pivot_longer(
            cols = starts_with("V"),       # or cols = -group
            names_to = "col_num",
            values_to = "cell_num"
          ) %>%
          select(-col_num)
        
        trends = as.data.frame(trends) %>%
          mutate(group = 1:nrow(.)) %>%
          pivot_longer(
            cols = starts_with("obs"),       # or cols = -group
            names_to = "col_num",
            values_to = "trend"
          ) %>%
          select(-col_num)
        
        trends$cell_num = groups$cell_num
        
        ## make raster containing correct trends for each neighbourhood
        for(i in 1:nrow(trends)) {
          values(rast)[trends$cell_num[i]] = trends$trend[i]
        }
        #plot(rast)
        
        ## convert raster to matrix of trends 
        lattice_trends <- as.array(rast)[,,1]
        
        ## after stable period of 500 time steps, shift temperatures by trend 
        for(row in 1:nrow) {
          for(col in 1:ncol) {
            for(time in 501:2000) {
              lattice_E_it_array[row,col,time] = (spat_grad*row + (20-spat_grad*85)) + lattice_trends[row,col]*(time-500)
            }
          }
        }
        #plot(x = 1:nrow, y = lattice_E_it_array[1:nrow,1,500]) ## at end of stable period, all are same
        #plot(x = 1:nrow, y = lattice_E_it_array[1:nrow,1,2000]) ## diff trends lead to diff values at end of simulation
        #plot(x = 1:2000, y = lattice_E_it_array[1,7,]) ## across space, trends differ :-)
         
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
        r_max = r + 0.25 ## max suitability
        alpha_tpc = 0.3 ## rise rate steepness
        beta_tpc = 0.3 ## decline rate steepness
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

        ## now, we are focusing on what happens at the leading edge 
        ## we can't simulate a really latitudinally-broad range because it would take a long time, but we want to avoid suitability becoming negative (i.e., trailing edge catching up to leading edge) at the leading edge during the simulation
        ## let's modify the tpc so that growth rate plateaus at temperatures hotter than Topt
        
        wb[Tb >= Trmax] = r_max - 0.25
        
        data.frame(r = wb, temperature = Tb) %>%
          ggplot(aes(x = temperature, y = r)) +
          geom_line() +
          geom_hline(yintercept = 0, colour = "red") +
          theme_bw()
        
        ## use tpc to translate temperature to growth rate 
        lattice_r_array = lattice_E_it_array
        lattice_r_array = -exp(beta_tpc*(lattice_r_array-Trmax)-8)-alpha_tpc*(lattice_r_array-Trmax)^2
        lattice_r_array = r_max*exp(1)^lattice_r_array - 0.25
        lattice_r_array[which(lattice_E_it_array >= Trmax)] = r_max - 0.25
        #plot(x = 1:nrow, y = lattice_r_array[1:nrow,1,600])
        
        ## add noise to growth rate
        lattice_r_array = lattice_r_array + lattice_ac_it 
        #plot(x = 1:2000, y = lattice_r_array[1,1,1:2000])
        
        ## save environmental array as raster
        filename =  paste0("outputs/data-processed/env-grids/range-shift-grid", 
                           rep, "_p", p, "_beta", beta, "_r", r, "_K", K, "_d", 
                           d, "_icp", icp, "_L", L, ".tif")
        writeRaster(rast(lattice_r_array), filename, overwrite = TRUE)
        
        # ## measure clim velocity after variation was added
        # rast = raster::brick(lattice_r_array[,,501:L],xmn = 0, xmx = 10, ymn = 0, ymx = 100)
        # ttrend = tempTrend(r = rast, th = 0.25*nlayers(rast)) ## set minimum # obs. to 1/4 time series length
        # #plot(ttrend)
        # temp_trend_est =mean(values(ttrend$slpTrends))
        # #0.025
        # 
        # ## use function to calculate spatial gradient in mean daily air temperature for each location:
        # spgrad = spatGrad(r = rast, projected = TRUE) ## our raster is projected to a coordinate system
        # #plot(spgrad)
        # spat_grad_est = mean(values(spgrad$Grad))
        # ## 0.257
        # 
        # ## calculate gradient based climate velocity:
        # gvocc = gVoCC(tempTrend = ttrend, spatGrad = spgrad)
        # #plot(gvocc)
        # clim_velo_est = mean(values(gvocc$voccMag))
        # ## 0.102
        
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
      pts$mean_temp_trend = temp_trend
      pts$icc = icc_star$value
      pts$between_var = between_var
      # pts$temp_trend_est = temp_trend_est
      # pts$spat_grad_est = spat_grad_est
      # pts$clim_velo_est = clim_velo_est
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
            dispersers = curr_allpops[y,x]*d
            ## figure out where each one disperses 
            ## individuals can disperse into neighbouring 24 cells
            othercells <- expand.grid(x = c(x-2,x-1,x,x+1,x+2), y = c(y-2,y-1,y,y+1,y+2)) 
            othercells = othercells[!(x == othercells$x & y == othercells$y),]
            othercells$num = 1:24
            
            ## randomly sample from 1:24
            sample <- data.frame(num = sample(1:24, dispersers, replace = TRUE))
            
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
    
  }
 }



