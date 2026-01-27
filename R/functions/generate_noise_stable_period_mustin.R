## function to create multiple time series of 1/f noise 
## the time series that are returned are of length L
## have spectral exponent equal to the parameter beta 
## and vary in their synchrony according to the synchrony parameter, p 
## when p = 1, all are synchronous 
## when p = 0, all are asynchronous
## method adapted from Lögdberg & Wennergren and Cuddington & Yodzis
generate_noise <- function(beta, ## spectral exponent, between 0 (low autocorelation) and 1 (high autocorrelation)
                           p, ## synchrony parameter, between 0 (low synchrony) and 1 (high synchrony)
                           n_ts, ## the number of unique time series to create
                           L1, ## the length of each time series 
                           L2 
){
  
  ## set a shared random element of phase:
  randf <- runif(524288, 0, 2*pi) 
  
  ## first, create time series with beta = 0 p = 1 to generate stable range
  ## create a series of 1/fB noise as described in Cuddington and Yodzis
  n_stable = rnorm(524288, mean = 0, sd = 1) ## random numbers with 0 mean and 1 sd 
  f = 1:524288 ## frequencies 
  a_stable <- 1/(f^(0/2))*exp(1)^n_stable ### amplitudes = 1/f^beta/2 * tiny random component
  
  ## now create noise of specified colour & synchrony 
  ## create a series of 1/fB noise as described in Cuddington and Yodzis
  n = rnorm(524288, mean = 0, sd = 1) ## random numbers with 0 mean and 1 sd 
  f = 1:524288 ## frequencies 
  a <- 1/(f^(beta/2))*exp(1)^n ### amplitudes = 1/f^beta/2 * tiny random component
  
  ## if patchy:
  if(p == 2) {

    p = 1
    
    ## for each requested time series:
    list_ts <- list()
    list_colours <- list()
    x = 1
    while(x <= n_ts) {
      
      ## apply phase shift for every 5 time series 
      if(x %in% seq(6, n_ts, by = 5)) {
        p = ifelse(p == 1, 0, 1)
        if(x %in% seq(51, n_ts, by = 50)) {
          p = ifelse(p == 1, 0, 1)
        }
      }

      phases_stable <- randf + (1-0)*runif(1, 0, 2*pi)
      phases <- randf + (1-p)*pi
      
      ## calculate wave coeffs and inverse dft
      complex_stable <- a_stable*cos(phases_stable) + a_stable*sin(phases_stable) ## complex coefficients
      
      dft_stable <- fft(complex_stable, inverse = T) ## inverse fast fourier transform the coefficients to get the temporal noise series
      noise_stable = as.numeric(dft_stable[1:L1]) ## crop the noise series to first L - 500 points
      
      complex <- a*cos(phases) + a*sin(phases) ## complex coefficients
      
      dft <- fft(complex, inverse = T) ## inverse fast fourier transform the coefficients to get the temporal noise series
      noise_shifting = as.numeric(dft[1:L2]) ## crop the noise series to first L - 500 points
      
      ## remove mean and change variance to 1:
      noise_stable <- noise_stable*1/sqrt(var(noise_stable))*sqrt(1^2)
      noise_stable <- noise_stable - mean(noise_stable)
      #plot(x = 1:L1, y = noise_stable)
      noise_shifting <- noise_shifting*1/sqrt(var(noise_shifting))*sqrt(1^2)
      noise_shifting <- noise_shifting - mean(noise_shifting)
      #noise_shifting(x = 1:L2, y = noise_shifting)
      noise = c(noise_stable, noise_shifting)
      #plot(x = 1:(L1+L2), y = noise)
      
      ## remove mean from stable and shifting period separately
      noise[1:L1] = noise[1:L1] - mean(noise[1:L1])
      noise[(L1+1):(L2+L1)] = noise[(L1+1):(L2+L1)] - mean(noise[(L1+1):(L2+L1)])
      #plot(x = 1:(L1+L2), y = noise)
      
      ## estimate noise colour from a linear regression of power spectrum:
      l <- length(noise_shifting)
      dft <- fft(noise_shifting)/l
      amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
      amp <- amp[1:(l/2)]	## remove second half of amplitudes (negative half)
      freq <- 1:(l/2)/l ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 
      
      ## create periodogram data by squaring amplitude of FFT output
      spectral <- data.frame(freq = freq, power = amp^2)
      
      # ggplot2::ggplot(data = spectral, ggplot2::aes(x = freq, y = power)) + ggplot2::geom_line() +
      #   ggplot2::scale_y_log10() + ggplot2::scale_x_log10() + ggplot2::geom_smooth(method = "lm") +
      #   ggplot2::theme_minimal()
      
      true_colour <- lm(data = spectral, log(power) ~ log(freq))
      
      ## save 
      list_ts[[x]] <- noise
      list_colours[[x]] <- as.numeric(true_colour$coefficients[2])
      
      x = x + 1
    }
    
    p = 2
    
  }
  
  else {
    ## for each requested time series:
    list_ts <- list()
    list_colours <- list()
    x = 1
    while(x <= n_ts) {
      
      ## apply phase shift
      phases_stable <- randf + (1-0)*runif(1, 0, 2*pi)
      phases <- randf + (1-p)*runif(1, 0, 2*pi)
      
      ## calculate wave coeffs and inverse dft
      complex_stable <- a_stable*cos(phases_stable) + a_stable*sin(phases_stable) ## complex coefficients
      
      dft_stable <- fft(complex_stable, inverse = T) ## inverse fast fourier transform the coefficients to get the temporal noise series
      noise_stable = as.numeric(dft_stable[1:L1]) ## crop the noise series to first L - 500 points
      
      complex <- a*cos(phases) + a*sin(phases) ## complex coefficients
      
      dft <- fft(complex, inverse = T) ## inverse fast fourier transform the coefficients to get the temporal noise series
      noise_shifting = as.numeric(dft[1:L2]) ## crop the noise series to first L - 500 points
      
      ## remove mean and change variance to 1:
      noise_stable <- noise_stable*1/sqrt(var(noise_stable))*sqrt(1^2)
      noise_stable <- noise_stable - mean(noise_stable)
      #plot(x = 1:L1, y = noise_stable)
      noise_shifting <- noise_shifting*1/sqrt(var(noise_shifting))*sqrt(1^2)
      noise_shifting <- noise_shifting - mean(noise_shifting)
      #noise_shifting(x = 1:L2, y = noise_shifting)
      noise = c(noise_stable, noise_shifting)
      #plot(x = 1:(L1+L2), y = noise)
      
      ## remove mean from stable and shifting period separately
      noise[1:L1] = noise[1:L1] - mean(noise[1:L1])
      noise[(L1+1):(L2+L1)] = noise[(L1+1):(L2+L1)] - mean(noise[(L1+1):(L2+L1)])
      #plot(x = 1:(L1+L2), y = noise)
      
      # ## make noise less intense by lowering standard deviation to 0.2
      noise = noise*1/sqrt(var(noise))*sqrt(0.2^2)
      #plot(x = 1:(L1+L2), y = noise)
      
      ## estimate noise colour from a linear regression of power spectrum:
      l <- length(noise_shifting)
      dft <- fft(noise_shifting)/l
      amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
      amp <- amp[1:(l/2)]	## remove second half of amplitudes (negative half)
      freq <- 1:(l/2)/l ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 
      
      ## create periodogram data by squaring amplitude of FFT output
      spectral <- data.frame(freq = freq, power = amp^2)
      
      # ggplot2::ggplot(data = spectral, ggplot2::aes(x = freq, y = power)) + ggplot2::geom_line() +
      #   ggplot2::scale_y_log10() + ggplot2::scale_x_log10() + ggplot2::geom_smooth(method = "lm") +
      #   ggplot2::theme_minimal()
      
      true_colour <- lm(data = spectral, log(power) ~ log(freq))
      
      ## save 
      list_ts[[x]] <- noise
      list_colours[[x]] <- as.numeric(true_colour$coefficients[2])
      
      x = x + 1
    }
  }
  
  return(list(list_ts, list_colours))
}

