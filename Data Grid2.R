library(tidyverse)
library(deSolve)
library(phaseR)
library(ggplot2)
library(foreach)
library(doSNOW)
library(doParallel)
library(viridis)
library(ggpmisc)
library(patchwork)

# the threshold for extinction -- stored as a global variable
THRESH <<- 1e-6


#### model ####
LV_multispecies_competition <- function(t, x, parameters) {
  
  with(as.list(parameters), {
    
    # zero out species below the threshold    
    x[x<THRESH] <- 0
    
    # Lotka-Volterra in matrix notation
    dxdt <- x*(r - A%*%x)
    
    # zero out species below the threshold 
    x[x < THRESH] <- 0
    
    # equivalent to the mathematical equation:
    # dxi/dt <- ri*xi - aii*xi^2 - sum(a_ij * xi * xj)
    
    # return rate of change
    return(list(dxdt))
    
  })
  
}

##HOI model
HOI <- function(t, x, parameters){
  with(as.list(c(x, parameters)), {
    
    #prevent infinite growing
    x[x < THRESH] <- 0
    x[x > 1/THRESH] <- 1/THRESH
    
    Bm <- A
    for (i in 1:nrow(A)){
      Bm[i,] <- B[i,,] %*% x
    }
    
    dxdt <- x * (r - (A + Bm) %*% x)
    
    #prevent infinite growing
    x[x < THRESH] <- 0
    x[x > 1/THRESH] <- 1/THRESH
    
    
    return(list(dxdt))
  })
}


# generate a tensor B matrix giving the HOIs
get_B <- function(n, mu = 0, sigma = 1, B_diag = 0, use_norm = TRUE){
  
  # choose suitable distribution
  if(use_norm){
    B <- round(array(rnorm(n * n * n, mean = mu, sd = sigma), c(n, n, n)), 2)
  }else{
    B <- round(array(rlnorm(n * n * n, meanlog = mu, sdlog = sigma), c(n, n, n)), 2)
  }
  
  # this zeros out effects of a species on itself, thus
  # only allowing HOIs among three unique species
  for (i in 1:n){
    for (j in 1:n){
      if (i == j){
        B[i,j,] <- 0
        B[i,,j] <- 0
        B[,i,j] <- 0
      }
    }
    
    # this zeros out effects when the two competitors are the same,
    diag(B[i,,]) <- 0
    # adds in a cubic self-regulation term
    B[i,i,i] <- B_diag
    # this ensures the matrix is symmetric with respect to the competitors,
    # ie.,g the xjxk is equal to xkxj
    B[i,,] <- (B[i, ,] + t(B[i,,])) / 2
    
  }
  
  return(B)
}

# test the stability after the time has reached to the end
Comm_statbility <- function(final_state, sta_func, sta_params){
  living_species <- which(final_state > THRESH)
  if(length(living_species) == 0){
    return(NA)
  }
  disturbance_state <- final_state
  for(i in living_species){
    disturbance <- runif(1, -0.05, 0.05) * final_state[i]
    disturbance_state[i] <- max(THRESH, final_state[i] + disturbance)
  }
  teststa_times <- seq(0, 2e4, by = 200)
  teststa_result <- ode(y = disturbance_state, times = teststa_times, func = sta_func, parms = sta_params, method = 'vode') %>%
    data.frame()
  teststa_state <- as.numeric(teststa_result[nrow(teststa_result), -1])
  
  living_after_disturbance <- which(teststa_state > THRESH)
  comparison_dist_result <-  setequal(living_species, living_after_disturbance)
  
  if(comparison_dist_result){
    relative_diff <- abs(teststa_state[living_species] - final_state[living_species]) / final_state[living_species]
    max_diff <- max(relative_diff)
    stability_result <- max_diff < 1e-3
  } else {
    stability_result <- FALSE
  }
  return(stability_result)
}


set.seed(555)
#set the parameters
if(exists("doParallel") && getDoParRegistered()) {
  stopImplicitCluster()
}

par_grid_orig <- expand.grid(nsim = 1:100, 
                             n_start = 30, 
                             A_mean = seq(0.05, 1 , length = 10),
                             A_sd = seq(0.05, 0.1, length = 5),
                             B_mean = seq(0.05, 1 , length = 10),
                             B_sd = seq(0.05, 0.1 , length = 5),
                             B_diag = c(0,0.5, 1, 1.5))

par_grid_new <- expand.grid(nsim = 1:100, 
                            n_start = 30, 
                            A_mean = seq(0.05, 0.5 , length = 10),
                            A_sd = seq(0.1, 0.3, length = 5),
                            B_mean = seq(0.05, 0.5 , length = 10),
                            B_sd = seq(0.1, 0.3 , length = 5),
                            B_diag = c(0,0.5, 1, 1.5))


par_grid <- par_grid_new %>% anti_join(par_grid_orig, by = names(par_grid_new))

cl <- makeCluster(64)
registerDoParallel(cl)
clusterExport(cl, c("LV_multispecies_competition", "THRESH", "HOI", "get_B", "Comm_statbility"))

system.time({result <- foreach(i = 1:nrow(par_grid), .combine = rbind, .packages = c("deSolve", "dplyr", "tidyr"), .inorder = FALSE) %dopar% {
  new_par_grid <- par_grid[i, ]
  
  # # specify the seed -- allows us to all have the same results
  set.seed(555 + i)
  # sample the number of species
  nspp <- new_par_grid$n_start 
  
  # matrix of interaction coefficients and generate the B matrix
  
  # Calculate log-scale parameters for A so that the longormal distributions have the target mean and variance
  sdlogA <- sqrt(log(1 + (new_par_grid$A_sd^2 / new_par_grid$A_mean^2)))
  meanlogA <- log(new_par_grid$A_mean) - (sdlogA^2) / 2
  
  # Calculate log-scale parameters for B so that the longormal distributions have the target mean and variance
  sdlogB <- sqrt(log(1 + (new_par_grid$B_sd^2 / new_par_grid$B_mean^2)))
  meanlogB <- log(new_par_grid$B_mean) - (sdlogB^2) / 2
  
  # generate all of the random matrices
  A <- matrix(rlnorm(nspp^2, meanlog = meanlogA, sdlog =  sdlogA), nrow = nspp, ncol = nspp)
  Blnorm <- get_B(nspp, mu = meanlogB, sigma = sdlogB, B_diag = new_par_grid$B_diag, use_norm = FALSE)
  Bnorm <- get_B(nspp, mu = new_par_grid$B_mean, sigma = new_par_grid$B_sd, B_diag = new_par_grid$B_diag, use_norm = TRUE)
  
  
  # assign a diagonal to A
  diag(A) <- 1
  
  # specify growth rates
  r <- 1
  
  # parameters and initial conditions
  parameters_2D = list(A = A, r = r)
  parameters_3Dnorm = list(A = A, B = Bnorm, r = r) 
  parameters_3Dlnorm = list(A = A, B = Blnorm, r = r) 
  
  # specify initial abundances
  state <- rep(0.1, nspp)
  
  # Time step and time increment for calculating the rate of change)
  times <- seq(0, 1e6, by = 1000)
  
  results2D <- result3Dnorm <- result3Dlnorm <- NULL
  
  # run the numerical solution, store as a tibble. Added a 10 second timeout to prevent things hanging
  result2D <- R.utils::withTimeout(ode(y = state, times = times, func = LV_multispecies_competition, parms = parameters_2D, method = 'vode') %>% data.frame() %>% as_tibble(), timeout = 10, onTimeout = "silent")
  result3Dnorm <- R.utils::withTimeout(ode(y = state, times = times, func = HOI, parms = parameters_3Dnorm, method = 'vode') %>% data.frame() %>% as_tibble(), timeout = 10, onTimeout = "silent")
  result3Dlnorm <- R.utils::withTimeout(ode(y = state, times = times, func = HOI, parms = parameters_3Dlnorm, method = 'vode') %>% data.frame() %>% as_tibble(), timeout = 10, onTimeout = "silent")
  
  # if we successfully ran the pairwise and at least one 3D, store the results
  if(!is.null(result2D) & !(is.null(result3Dnorm) & is.null(result3Dlnorm))){
    # clean up the result by assigning extinct species to zero
    result2D[result2D<THRESH] <- 0
    
    # Getting the final state of the data
    final_state_2D <- as.numeric(result2D[nrow(result2D), -1])
    
    # running the stability function
    stability_result_2D <- Comm_statbility(final_state = final_state_2D,
                                           sta_func = LV_multispecies_competition,
                                           sta_params = parameters_2D)
    
    # importing the data
    n_end_2D <- sum(final_state_2D > THRESH)
  }else{
    # otherwise, store NAs
    n_end_2D <- NA
    stability_result_2D <- NA
  }
  
  # if we successfully ran both the simulations, store the results
  if(!is.null(result3Dnorm)){
    # clean up the result by assigning extinct species to zero
    result3Dnorm[result3Dnorm<THRESH] <- 0
    
    # Getting the final state of the data
    final_state_3Dnorm <- as.numeric(result3Dnorm[nrow(result3Dnorm), -1])
    
    # running the stability function
    stability_result_3Dnorm <- Comm_statbility(final_state = final_state_3Dnorm,
                                               sta_func = HOI,
                                               sta_params = parameters_3Dnorm)
    
    # importing the data
    n_end_3Dnorm <- sum(final_state_3Dnorm > THRESH)
  }else{
    # otherwise, store NAs
    n_end_3Dnorm <- NA
    stability_result_3Dnorm <- NA
  }	
  
  # if we successfully ran both the simulations, store the results
  if(!is.null(result3Dlnorm)){
    # clean up the result by assigning extinct species to zero
    result3Dlnorm[result3Dlnorm<THRESH] <- 0
    
    # Getting the final state of the data
    final_state_3Dlnorm <- as.numeric(result3Dlnorm[nrow(result3Dlnorm), -1])
    
    # running the stability function
    stability_result_3Dlnorm <- Comm_statbility(final_state = final_state_3Dlnorm,
                                                sta_func = HOI,
                                                sta_params = parameters_3Dlnorm)
    
    # importing the data
    n_end_3Dlnorm <- sum(final_state_3Dlnorm > THRESH)
  }else{
    # otherwise, store NAs
    n_end_3Dlnorm <- NA
    stability_result_3Dlnorm <- NA
  }
  
  
  rbind(
    data.frame(
      nsim = new_par_grid$nsim,
      n_start = nspp,
      A_mean = new_par_grid$A_mean,
      A_sd = new_par_grid$A_sd,
      B_mean = new_par_grid$B_mean,
      B_sd = new_par_grid$B_sd,
      B_diag = new_par_grid$B_diag,
      n_end = n_end_2D,
      stability_result = stability_result_2D,
      model = "Pairwise"
    ),
    data.frame(
      nsim = new_par_grid$nsim,
      n_start = nspp,
      A_mean = new_par_grid$A_mean,
      A_sd = new_par_grid$A_sd,
      B_mean = new_par_grid$B_mean,
      B_sd = new_par_grid$B_sd,
      B_diag = new_par_grid$B_diag,
      n_end = n_end_3Dnorm,
      stability_result = stability_result_3Dnorm,
      model = "HOI - norm"
    ),
    data.frame(
      nsim = new_par_grid$nsim,
      n_start = nspp,
      A_mean = new_par_grid$A_mean,
      A_sd = new_par_grid$A_sd,
      B_mean = new_par_grid$B_mean,
      B_sd = new_par_grid$B_sd,
      B_diag = new_par_grid$B_diag,
      n_end = n_end_3Dlnorm,
      stability_result = stability_result_3Dlnorm,
      model = "HOI - lognorm"
    )
  )
  
  
}})


write_csv(result, "~/Dropbox/Temp/new_combinations_norm_and_lognorm.csv")