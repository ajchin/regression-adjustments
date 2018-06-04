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

run_sim = function(param, g, variance_factor, n_reps, n_cores, pid=NULL) {
  # param: row of parameters
  
  # creates a list containing the data
  covariate_fns_for_response=list(
    frac_nbh = fraction_trt_nbrs,
    num_nbh = number_trt_nbrs
  )
  data = generate_covariate_data(g, covariate_fns_for_response)
  
  # generate response
  data$y = linear_response(data$w, data$x_obs, param)
  
  return(c(
    pid=pid,
    dm=data %>% difference_in_means,
    adjusted=data %>% linear_adjustment,
    var_est=data %>% linear_variance_estimate(variance_factor)
  ))
}


load('data/caltech.Rdata')


n_reps = 200
n_cores = 56

registerDoParallel(cores=n_cores)

covariate_fns_for_estimator=list(
  frac_nbh = fraction_trt_nbrs,
  num_nbh = number_trt_nbrs
)
start = proc.time()
vf = precompute_variance(g, covariate_fns_for_estimator, n_boot_reps=200, n_cores=n_cores)
vf
print(proc.time() - start)


params = purrr::cross(list(
  beta0 = list(c(0, 0, 0), c(0, 0.1, 0.05)),
  beta1 = list(c(1, 0, 0), c(1, 0, 0.1), c(1, 0.5, 0), c(1, 0.5, 0.1)),
  noise_sd = c(1)
))

param = params[[2]]

print('Running simulation...')
# Run simulation
start = proc.time()
estimates = foreach(i = 1:length(params), .combine=rbind) %do% {
  param = params[[i]]
  print(unlist(param))
  foreach(rep = 1:n_reps, .combine=rbind) %dopar% {
    run_sim(param, g, vf, n_reps, n_cores, pid = i)
  }
} %>% data.frame
print(proc.time())

d_bar = mean(degree(g))

truth = sapply(params, function(p) {with(p, beta1[1] - beta0[1] + beta1[2] + d_bar * beta1[3])})

estimates %>% mutate(
  truth = truth[pid],
  accepts = abs(adjusted - truth) / sqrt(var_est) < qnorm(0.95)
) %>%
  group_by(pid) %>% 
  summarise(mean(accepts))

estimates %>% mutate(truth=truth[pid]) %>% group_by(pid) %>% summarise(truth=first(truth), mean(adjusted))


head(estimates)
var(estimates$adjusted)
mean(estimates$var_est)
var(estimates$var_est)
