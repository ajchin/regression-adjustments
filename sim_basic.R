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
    dm_var_est = data %>% dm_variance_estimate,
    adjusted=data %>% linear_adjustment,
    var_est=data %>% linear_variance_estimate(variance_factor)
  ))
}


load('data/caltech.Rdata')


n_reps = 1000
n_cores = 56

registerDoParallel(cores=n_cores)

covariate_fns_for_estimator=list(
  frac_nbh = fraction_trt_nbrs,
  num_nbh = number_trt_nbrs
)
start = proc.time()
vf = precompute_variance(g_fb,  covariate_fns_for_estimator, n_boot_reps=200, n_cores=n_cores)
vf
print(proc.time() - start)


params = purrr::cross(list(
  beta0 = list(c(0, 0, 0), c(0, 0, 0.05), c(0, 0.1, 0), c(0, 0.1, 0.05)),
  beta1 = list(c(1, 0, 0), c(1, 0, 0.1), c(1, 0.5, 0), c(1, 0.5, 0.1)),
  noise_sd = c(1, 3)
))

param = params[[1]]

print('Running simulation...')
# Run simulation
start = proc.time()
estimates = foreach(i = 1:length(params), .combine=rbind) %do% {
  param = params[[i]]
  print(unlist(param))
  foreach(rep = 1:n_reps, .combine=rbind) %dopar% {
    run_sim(param, g_fb, vf, n_reps, n_cores, pid = i)
  }
} %>% data.frame
print(proc.time())
write.csv(estimates, file='results/sim_basic.csv')

d_bar = mean(degree(g_fb))

truth = sapply(params, function(p) {with(p, beta1[1] - beta0[1] + beta1[2] + d_bar * beta1[3])})

df_dm = estimates %>% 
  select(pid, estimate=dm, var_est=dm_var_est) %>%
  mutate(
    truth=truth[pid],
    accepts=abs(estimate - truth) / sqrt(var_est) < qnorm(0.95)
  ) %>% 
  group_by(pid) %>% 
  summarise(
    truth = first(truth),
    mean=mean(estimate),
    bias=mean - truth,
    sd=sd(estimate),
    avg_sd_est = mean(sqrt(var_est)),
    sd_ratio = avg_sd_est / sd,
    coverage=mean(accepts)
  ) %>%
  ungroup

df_adjusted = estimates %>% 
  select(pid, estimate=adjusted, var_est=var_est) %>%
  mutate(
    truth=truth[pid],
    accepts=abs(estimate - truth) / sqrt(var_est) < qnorm(0.95)
  ) %>% 
  group_by(pid) %>% 
  summarise(
    truth = first(truth),
    mean=mean(estimate),
    bias=mean - truth,
    sd=sd(estimate),
    avg_sd_est = mean(sqrt(var_est)),
    sd_ratio = avg_sd_est / sd,
    coverage=mean(accepts)
  ) %>%
  ungroup

cbind(
  df_dm %>% select(pid, truth, bias_dm=bias, sd_dm = sd, sd_ratio_dm=sd_ratio, coverage_dm=coverage),
  df_adjusted %>% select(bias_adj=bias, sd_adj = sd, sd_ratio_adj=sd_ratio, coverage_adj=coverage)
) %>% filter(pid < 17) %>%
  select(pid, truth, bias_dm, bias_adj, sd_dm, sd_adj, sd_ratio_dm, sd_ratio_adj, coverage_dm, coverage_adj) %>% 
  xtable(digits=3) %>% print.xtable(include.rownames=FALSE)

foreach(p = params[1:16], .combine=rbind) %do% {c(p$beta0[2:3], p$beta1[2:3])} %>% xtable %>% print.xtable(include.rownames=FALSE)



summarised = rbind(
  estimates %>% select(pid, estimate=dm, var_est=dm_var_est) %>% mutate(estimator='dm'),
  estimates %>% select(pid, estimate=adjusted, var_est=var_est) %>% mutate(estimator='adjusted')
) %>% mutate(
  truth=truth[pid],
  accepts=abs(estimate - truth) / sqrt(var_est) < qnorm(0.95)
) %>% 
  group_by(pid, estimator) %>% 
  summarise(
    truth = first(truth),
    mean=mean(estimate),
    bias=mean - truth,
    sd=sd(estimate),
    avg_sd_est = mean(sqrt(var_est)),
    sd_ratio = avg_sd_est / sd,
    coverage=mean(accepts)
  ) %>%
  ungroup

summarised %>% 
  select(estimator, bias, sd, sd_ratio, coverage)

tmp = cbind(summarised %>% filter(estimator == 'adjusted'), summarised %>% filter(estimator == 'dm'))
