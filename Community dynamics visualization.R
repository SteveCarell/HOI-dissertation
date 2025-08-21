library(tidyverse)
library(deSolve)
library(phaseR)

# the threshold for extinction -- stored as a global variable
THRESH <<- 1e-8


#### model ####
LV_multispecies_competition <- function(t, x, parameters) {
  
  with(as.list(parameters), {
    
    # zero out species below the threshold    
    x[x<THRESH] <- 0
    
    # Lotka-Volterra in matrix notation
    dxdt <- r*x*(1 - A%*%x)
    
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
      }
    }
    
    # this zeros out effects when the two competitors are the same,
    diag(B[i,,]) <- 0
    # adds in a cubic self-regulation term
    B[i,i,i] <- 0
    # this ensures the matrix is symmetric with respect to the competitors,
    # ie.,g the xjxk is equal to xkxj
    B[i,,] <- (B[i, ,] + t(B[i,,])) / 2
    
  }
  
  return(B)
}

# visualization example for community dynamics
community_visualization <- function(model_func, model_params, initial_state, 
                                        total_time, perturbation_time, time_step,
                                        apply_perturbation = TRUE) {
  
  if (apply_perturbation) {
    times_before <- seq(0, perturbation_time, by = time_step)
    times_after <- seq(0, total_time - perturbation_time, by = time_step)
    result <- ode(y = initial_state, times = times_before, func = model_func, parms = model_params) %>%
      data.frame() %>% as_tibble()
    result[result < THRESH] <- 0
    final_state <- as.numeric(result[nrow(result), -1])
    living_species_indices <- which(final_state > THRESH)
    disturbance_state <- final_state
    if (length(living_species_indices) > 0) {
      for (i in living_species_indices) {
        disturbance <- runif(1, -0.15, 0.15) * final_state[i]
        disturbance_state[i] <- max(THRESH, final_state[i] + disturbance)}}
    
    result_recovery <- ode(y = disturbance_state, times = times_after, func = model_func, parms = model_params) %>%
      data.frame() %>% as_tibble()
    time_shift <- result$time[nrow(result)]
    result_recovery$time <- result_recovery$time + time_shift
    disturbance_row <- result[nrow(result), ]
    disturbance_row[-1] <- as.list(disturbance_state)
    combined_results <- rbind(result, disturbance_row, result_recovery)
    
  } else {
    full_times <- seq(0, total_time, by = time_step)
    combined_results <- ode(y = initial_state, times = full_times, func = model_func, parms = model_params) %>%
      data.frame() %>% as_tibble()
    combined_results[combined_results < THRESH] <- 0}
  final_abundances <- as.numeric(combined_results[nrow(combined_results), -1])
  surviving_species_indices <- which(final_abundances > THRESH)
  surviving_species_names <- paste0('X', surviving_species_indices)
  plot_data <- combined_results %>%
    gather(species, abundance, -time) %>%
    filter(species %in% surviving_species_names)
  return(plot_data)
}

#define the parameters
set.seed(10)
nspp <- 20

par_grid <- expand.grid(n_start = 20, 
                        A_mean = 0.2,
                        A_sd = 0.2,
                        B_mean = 0.2,
                        B_sd = 0.2
)

sdlogA <- sqrt(log(1 + (par_grid$A_sd^2 / par_grid$A_mean^2)))
meanlogA <- log(par_grid$A_mean) - (sdlogA^2) / 2

sdlogB <- sqrt(log(1 + (par_grid$B_sd^2 / par_grid$B_mean^2)))
meanlogB <- log(par_grid$B_mean) - (sdlogB^2) / 2
A <- matrix(rlnorm(nspp^2, meanlog = meanlogA, sdlog =  sdlogA), nrow = nspp, ncol = nspp)
Blnorm <- get_B(nspp, mu = meanlogB, sigma = sdlogB, B_diag = par_grid$B_diag, use_norm = FALSE)
diag(A) <- 1
r <- runif(nspp, 0, 1)
state <- rep(0.01, nspp)
total_time <- 8000
perturb_time <- 5000
time_step <- 100
params_2D <- list(r = r, A = A)
params_3D <- list(r = r, A = A, B = Blnorm)

pairwise_pert <- community_visualization(model_func = LV_multispecies_competition, model_params = params_2D, initial_state = state,
                                             total_time = total_time, perturbation_time = perturb_time, time_step = time_step,
                                             apply_perturbation = TRUE) %>% mutate(model_type = 'Pairwise', perturbation = 'With Perturbation')

HOI_pert <- community_visualization(model_func = HOI, model_params = params_3D, initial_state = state,
                                        total_time = total_time, perturbation_time = perturb_time, time_step = time_step,
                                        apply_perturbation = TRUE) %>% mutate(model_type = 'HOI', perturbation = 'With Perturbation')

#bind the plot data
all_HOI_data <- bind_rows(pairwise_pert, HOI_pert)

all_HOI_data$model_type <- factor(all_HOI_data$model_type,
                                  levels = c('Pairwise', 'HOI'))

# plot data
vline_data <- data.frame(perturbation = 'With Perturbation', xint = 5000)
final_plot <- ggplot(all_HOI_data, aes(x = time, y = abundance, color = species)) +
  geom_line(linewidth = 1.0, alpha = 0.8) + 
  facet_grid(~ model_type) +
  geom_vline(data = vline_data, aes(xintercept = xint), inherit.aes = FALSE, linetype = 'dashed', color = "black") +
  geom_text(data = vline_data, aes(x = xint, y = Inf, label = 'Perturbation'), inherit.aes = FALSE, hjust = 1.1, vjust = 1.5, color = "black", size = 10) +
  labs(x = 'Time', y = 'Species abundance') +
  theme_classic(base_size = 30) +
  theme(legend.position = 'none')

print(final_plot)