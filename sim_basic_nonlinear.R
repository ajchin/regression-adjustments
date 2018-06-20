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
    num_nbh = number_trt_nbrs
    #frac_nbh = fraction_trt_nbrs,
    #frac_nbh2 = fraction_trt_nbrs2
  )
  data = generate_covariate_data(g, covariate_fns_for_response)
  
  # generate response
  data$y = basic_nonlinear_response(data$w, data$x_obs, param)

  return(c(
    pid=pid,
    dm=data %>% difference_in_means,
    dm_var_est = data %>% dm_variance_estimate,
    adjusted=data %>% lr_crossfit(n_folds=2)
    #var_est=data %>% linear_variance_estimate(variance_factor)
  ))
}


load('data/stanford.Rdata')


n_reps = 200
n_cores = 56

registerDoParallel(cores=n_cores)

covariate_fns_for_estimator=list(
  frac_nbh = fraction_trt_nbrs,
  frac_nbh2 = fraction_trt_nbrs2
)

# start = proc.time()
# vf = precompute_variance(g_sm_large,  covariate_fns_for_estimator, n_boot_reps=200, n_cores=n_cores)
# vf
# print(proc.time() - start)


# params = purrr::cross(list(
#   beta0 = list(c(0, 0, 0), c(0, 0, 0.05), c(0, 0.1, 0), c(0, 0.1, 0.05)),
#   beta1 = list(c(1, 0, 0), c(1, 0, 0.1), c(1, 0.5, 0), c(1, 0.5, 0.1)),
#   noise_sd = c(1, 3)
# ))

# 
# param = list(
#   beta0 = c(-2, 0, 0.05), beta1 = c(1, 0.05, 0.1), noise_sd = 1
# )

param = list(noise_sd = 1)

# COMPUTE TRUTH
data = generate_covariate_data(g, covariate_fns_for_response)
truth1 = foreach(i = 1:500, .combine=c) %dopar% {binary_response(rep(1, 5000), data$x_trt, param)}
truth0 = foreach(i = 1:500, .combine=c) %dopar% {binary_response(rep(0, 5000), data$x_ctrl, param)}
true_ATE = mean(truth1) - mean(truth0)
true_ATE

results = foreach(rep = 1:n_reps, .combine=rbind) %dopar% {
  run_sim(param, g_sm_large, variance_factor = NULL, n_reps, n_cores, pid = i)
} %>% data.frame
hist(results$adjusted, breaks=30)
mean(results$adjusted)
var(results$adjusted)



df = data.frame(w = data$w, y = data$y, x_obs = data$x_obs, x_trt = data$x_trt, x_ctrl = data$x_ctrl)

lr_for_boot = function(df, ids) {
  df = df[ids,]
  x_obs = df[,3:4]
  x_trt = df[,5:6]
  x_ctrl = df[,7:8]
  names(x_obs) = c('frac_nbh', 'frac_nbh2')
  names(x_ctrl) = names(x_obs)
  names(x_trt) = names(x_obs)
  data_boot = list(y = df$y, w = df$w, x_obs=x_obs, x_trt = x_trt, x_ctrl = x_ctrl)
  lr_crossfit(data_boot, n_folds=2)
}

boot_results = boot(df, lr_for_boot, R=200)
boot_results$t %>% var


w = data$w

# create fold ids
n = length(data$y)
fold_ids = sample(rep(1:n_folds, ceiling(n / n_folds))[1:n])

.get_models = function(data, fold_ids) {
  
  # train models
  models = foreach (fold = 1:n_folds) %do% {
    foreach (group = c(0, 1), .final = function(x) setNames(x, c('fit0', 'fit1'))) %do% {
      id_train = fold_ids != fold & w == group
      X_train = data$x_obs[id_train,, drop=FALSE] %>% select(one_of(vars))
      y_train = data$y[id_train]
      glm(y_train ~ ., data = X_train, family='binomial')
    }
  }
  
  models
}

.predict = function(models, data, fold_ids) {
  
  predictions = foreach (fold = 1:n_folds, .combine = rbind) %do% {
    id_test = fold_ids == fold 
    glo_trt = 1 * (predict(models[[fold]]$fit1, newdata=data$x_trt[id_test,,drop=FALSE], type='response') > 0.5)
    glo_ctrl = 1 * (predict(models[[fold]]$fit0, newdata=data$x_ctrl[id_test,,drop=FALSE], type='response') < 0.5)
    obs_trt = 1 * (predict(models[[fold]]$fit1, newdata=data$x_obs[id_test,,drop=FALSE], type='response') > 0.5)
    obs_ctrl = 1 * (predict(models[[fold]]$fit0, newdata=data$x_obs[id_test,,drop=FALSE], type='response') < 0.5)
    data.frame(glo_trt=glo_trt, glo_ctrl=glo_ctrl, obs_trt=obs_trt, obs_ctrl=obs_ctrl)
  }
  
  #
  adjust_glo = with(predictions, mean(glo_trt) - mean(glo_ctrl))
  adjust_obs = with(predictions, mean(obs_trt[w==1] - mean(obs_ctrl[w==0])))
  
  adjustment = adjust_glo - adjust_obs
  adjustment
  
}

models = .get_models(data, fold_ids)
dm = mean(data$y[w==1]) - mean(data$y[w==0])
dm
.predict(models, data, fold_ids)
adjustments = foreach(i = 1:100, .combine=c) %dopar% {
  newdata = generate_covariate_data(g, covariate_fns_for_response)
  .predict(models, newdata, fold_ids)
}
var(adjustments)

# RSS
(1 /  length(data$y)) * (sum((data$y[w==1] - predictions$obs_trt[w==1])^2) + sum((data$y[w==0] - predictions$obs_ctrl[w==0])^2))



print('Running simulation...')
# Run simulation
start = proc.time()
estimates = foreach(i = 1:length(params), .combine=rbind) %do% {
  param = params[[i]]
  print(unlist(param))
  results = foreach(rep = 1:n_reps, .combine=rbind) %dopar% {
    run_sim(param, g_fb, vf, n_reps, n_cores, pid = i)
  } %>% data.frame
} %>% data.frame
print(proc.time())
write.csv(estimates, file='results/sim_basic_nonlinear.csv')

d_bar = mean(degree(g_fb))

truth = sapply(params, function(p) {with(p, beta1[1] - beta0[1] + beta1[2] + beta1[3])})

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
