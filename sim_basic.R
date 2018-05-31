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

run_sim = function(param, g, Gamma, Delta, n_reps, n_cores, pid) {
  # param: row of parameters
  
  # creates a list containing the data
  covariate_fns_for_response=list(
    frac_nbh = fraction_trt_nbrs,
    num_nbh = number_trt_nbrs
  )
  data = generate_covariate_data(g, covariate_fns_for_response)
  
  # generate response
  data$y = linear_response(data$w, data$x_obs, param, noise_sd=3)
  
  return(c(
    dm=data %>% difference_in_means,
    adjusted=data %>% linear_adjustment,
    #naive_var_est=data %>% linear_variance_estimate,
    var_est=data %>% linear_variance_estimate(Gamma, Delta, vars=NULL)
  ))
}


g = sample_smallworld(dim=1, size=2000, nei=5, p=0.1)

param = list(alpha_trt = 1, beta_trt = c(1, 2), alpha_ctrl = 0, beta_ctrl = c(0.5, 1))

n_reps = 1000
n_cores = 48

registerDoParallel(cores=n_cores)



covariate_fns_for_estimator=list(
  frac_nbh = fraction_trt_nbrs,
  num_nbh = number_trt_nbrs
)

matrices = precompute_matrices(g, covariate_fns_for_estimator, n_boot_reps=n_reps, n_cores=n_cores)
Gamma = matrices$Gamma
Delta = matrices$Delta

print('Running simulation...')
# Run simulation
start = proc.time()
estimates = foreach(rep = 1:n_reps, .combine=rbind) %dopar% {
  run_sim(param, g, Gamma, Delta, n_reps, n_cores, pid=NULL)
} %>% data.frame
print(proc.time())

var(estimates$adjusted)
mean(estimates$var_est)
var(estimates$var_est)
