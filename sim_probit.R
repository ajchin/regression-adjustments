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

run_sim = function(param, g, vf1, vf2, n_reps, n_cores, pid=NULL) {
  # param: row of parameters
  
  # creates a list containing the data
  covariate_fns_for_response=list(
    frac_nbh = fraction_trt_nbrs,
    frac_nbh2 = fraction_trt_nbrs2
  )
  data = generate_covariate_data(g, covariate_fns_for_response)
  
  # generate response
  data$y = dynamic_time_response(data$w, g, param)
  
  vars1 = c('frac_nbh')
  vars2 = c('frac_nbh', 'frac_nbh2')
  
  return(c(
    pid=pid,
    dm=data %>% difference_in_means,
    hajek=data %>% hajek(g, 'frac_nbh', threshold=0.75, p_design=0.5),
    ols1=data %>% linear_adjustment(vars=vars1),
    #var_est1=data %>% linear_variance_estimate(vf1, vars=vars1),
    ols2=data %>% linear_adjustment(vars=vars2),
    #var_est2=data %>% linear_variance_estimate(vf2, vars=vars2),
    lr1 = data %>% lr_crossfit(n_folds=3, vars=vars1),
    lr2 = data %>% lr_crossfit(n_folds=3, vars=vars2)
  ))
}


#load('data/caltech.Rdata')
set.seed(2018)
#g_sm_large = sample_smallworld(dim=1, size=5000, nei=10, p=0.1)
#save(g_sm_large, file='data/smallworld_large.Rdata')
load('data/smallworld_large.Rdata')


n_reps = 500
n_cores = 48

registerDoParallel(cores=n_cores)
# 
# start = proc.time()
# vf1 = precompute_variance(
#   g_sm, 
#   list(frac_nbh=fraction_trt_nbrs), 
#   n_boot_reps=500, 
#   n_cores=n_cores
# )
# vf1
# 
# vf2 = precompute_variance(
#   g_sm, 
#   list(frac_nbh=fraction_trt_nbrs, frac_nbh2=fraction_trt_nbrs2), 
#   n_boot_reps=500, 
#   n_cores=n_cores
# )
# vf2
# print(proc.time() - start)

params = purrr::cross(list(
  b_intercept = -0.5,
  b_direct = 1,
  b_spill = c(0, 0.25, 0.5, 0.75, 1),#0.1, 0.2, 0.3, 0.4, 0.5),
  max_t = c(2, 4),
  is_probit = TRUE,
  noise_sd = 1 # c(1, 3)
))


print('Running simulation...')
# Run simulation
start = proc.time()
foreach(i = 1:length(params), .combine=rbind) %do% {
  param = params[[i]]
  print(unlist(param))
  estimates = foreach(rep = 1:n_reps, .combine=rbind) %dopar% {
    run_sim(param, g_sm_large, vf1, vf2, n_reps, n_cores, pid = i)
  } %>% data.frame
  write.table(estimates, file='results/sim_lim/results_lr.csv', append=TRUE, col.names=FALSE, row.names=FALSE, sep=',')
  print(proc.time())
}



# results = read.csv('results/sim_lim/results_lr.csv', header=FALSE)
# names(results) = c('pid', 'dm', 'hajek', 'ols1', 'ols2', 'lr1', 'lr2')
# true_ATE
# results %>% 
#   group_by(pid) %>% 
#   summarise_each(funs(mean(., na.rm = TRUE)), dm, hajek, ols1, ols2, lr1, lr2) %>%
#   left_join(true_ATE) %>% gather(k, v, dm, hajek, ols1, ols2, lr1, lr2) %>%
#   mutate(bias_sq = (ATE - v)^2) %>% ggplot(aes(pid, bias_sq, group=k, colour=k)) + geom_point() + geom_line() + theme_bw()
# results %>% 
#   group_by(pid) %>% 
#   summarise_each(funs(var(., na.rm = TRUE)), dm, hajek, ols1, ols2, lr1, lr2)

# 