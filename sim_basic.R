# This simulation is for examining variance estimates under a basic exogenous LIM model.

library(igraph)
library(foreach)
library(doParallel)
library(dplyr)
library(tidyr)
library(broom)

source('functions/data_generators.R')
source('functions/covariate_functions.R')
source('functions/response_functions.R')
source('functions/existing_estimators.R')
source('functions/proposed_estimators.R')
source('functions/variance_estimators.R')
source('functions/precompute_matrices.R')

run_sim = function(param, g, variance_factor, n_reps, n_cores, pid) {
  # param: row of parameters
  
  # creates a list containing the data
  covariate_fns_for_response=list(
    frac_nbh = fraction_trt_nbrs,
    num_nbh = number_trt_nbrs
  )
  data = generate_covariate_data(g, covariate_fns_for_response)
  
  # generate response
  data$y = linear_response(data$w, data$x_obs, param, noise_sd=1.5)
  
  #v = data %>% linear_variance_estimate(eta, Delta, vars=NULL)
  return(c(
    dm=data %>% difference_in_means,
    adjusted=data %>% linear_adjustment,
    var_est=data %>% linear_variance_estimate(variance_factor)
  ))
}


g = sample_smallworld(dim=1, size=2000, nei=5, p=0.1)

param = list(alpha_trt = 1, beta_trt = c(2, 0.5), alpha_ctrl = 0, beta_ctrl = c(0.5, 0.2))

n_reps = 2000
n_cores = 48

registerDoParallel(cores=n_cores)



covariate_fns_for_estimator=list(
  frac_nbh = fraction_trt_nbrs,
  num_nbh = number_trt_nbrs
)

vf = precompute_variance(g, covariate_fns_for_estimator, n_boot_reps=500, n_cores=n_cores)
vf

print('Running simulation...')
# Run simulation
start = proc.time()
estimates = foreach(rep = 1:n_reps, .combine=rbind) %dopar% {
  run_sim(param, g, vf, n_reps, n_cores, pid=NULL)
} %>% data.frame
print(proc.time())

head(estimates)
var(estimates$adjusted)
mean(estimates$var_est)
var(estimates$var_est)
