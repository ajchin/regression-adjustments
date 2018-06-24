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

run_sim = function(param, g, n_reps, n_cores, pid=NULL) {
  # param: row of parameters
  
  # creates a list containing the data
  covariate_fns_for_response=list(
    num_nbh = number_trt_nbrs,
    frac_nbh = fraction_trt_nbrs
    #frac_nbh2 = fraction_trt_nbrs2
  )
  data = generate_covariate_data(g, covariate_fns_for_response)
  
  # generate response
  data$y = nonlinear_response(data$w, data$x_obs, param)
  
  
  
  n_folds = 2
  n = length(data$y)
  fold_ids = sample(rep(1:n_folds, ceiling(n / n_folds))[1:n])
  
  return(c(
    pid=pid,
    dm=data %>% difference_in_means,
    dm_var_est = data %>% dm_variance_estimate,
    hajek=data %>% hajek(g, 'frac_nbh', threshold=0.75),
    ols=data %>% linear_adjustment,
    ols_var_est = data %>% linear_variance_estimate(),
    gam=data %>% gam_crossfit(n_folds, fold_ids),
    gam_var_est = data %>% gam_boot(n_folds, fold_ids, g, covariate_fns_for_response, B=50)
  ))
}

param = list(noise_sd = 1)

load('data/stanford.Rdata')



n_reps = 200
n_cores = 50

registerDoParallel(cores=n_cores)

# Compute truth
data = generate_covariate_data(g, covariate_fns_for_response)
truth1 = foreach(i = 1:2000, .combine=c) %dopar% {nonlinear_response(rep(1, vcount(g)), data$x_trt, param)}
truth0 = foreach(i = 1:2000, .combine=c) %dopar% {nonlinear_response(rep(0, vcount(g)), data$x_ctrl, param)}
true_ATE = mean(truth1) - mean(truth0)
true_ATE


start = proc.time()
results = foreach(rep = 1:n_reps, .combine=rbind) %dopar% {
  run_sim(param, g, n_reps, n_cores, pid = i)
} %>% data.frame
print(proc.time() - start)
write.csv(results, file='results/sim_nonlinear.csv')







#OLS
run_sim_ols = function(param, g, vf, n_reps, n_cores, pid=NULL) {
  # param: row of parameters
  
  # creates a list containing the data
  covariate_fns_for_response=list(
    num_nbh = number_trt_nbrs,
    frac_nbh = fraction_trt_nbrs
    #frac_nbh2 = fraction_trt_nbrs2
  )
  data = generate_covariate_data(g, covariate_fns_for_response)
  
  # generate response
  data$y = nonlinear_response(data$w, data$x_obs, param)
  
  
  return(c(
    ols=data %>% linear_adjustment,
    ols_var_est = data %>% linear_variance_estimate(vf)
  ))
}

vf = precompute_variance(
  g, 
  list(frac_nbh=fraction_trt_nbrs, num_nbh=number_trt_nbrs), 
  n_boot_reps=500, 
  n_cores=n_cores
)

start = proc.time()
results_ols = foreach(rep = 1:n_reps, .combine=rbind) %dopar% {
  run_sim_ols(param, g, vf, n_reps, n_cores, pid = i)
} %>% data.frame
print(proc.time() - start)



results_all = cbind(results, results_ols)
results_all %>% 
  gather(estimator, estimate, dm, hajek, gam, ols) %>% 
  mutate(true_ATE=true_ATE, bias=estimate-true_ATE) %>% 
  group_by(estimator) %>%
  summarise(mean_estimate=mean(estimate), bias=mean(bias), bias_pct = mean(abs(bias/true_ATE)), se=sd(estimate)) %>%
  xtable(digits=3) %>% print.xtable(include.rownames=FALSE)

with(results_all, sqrt(mean(gam_var_est)) / sd(gam))
with(results_all, sqrt(mean(ols_var_est)) / sd(ols))

mean(results$dm)
var(results$dm)
mean(hajek)
mean(results$gam)
var(results$gam)
mean(results$gam_var_est)

