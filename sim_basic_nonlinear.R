run_sim = function(param, g, n_reps, n_cores, pid=NULL) {
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
  
  
  n_folds = 2
  n = length(data$y)
  fold_ids = sample(rep(1:n_folds, ceiling(n / n_folds))[1:n])
  
  return(c(
    pid=pid,
    dm=data %>% difference_in_means,
    dm_var_est = data %>% dm_variance_estimate,
    adjusted=data %>% loess_crossfit(n_folds, fold_ids)
    #var_est=data %>% loess_boot(n_folds, fold_ids, g, covariate_fns_for_response, B=50)
    #var_est=data %>% linear_variance_estimate(variance_factor)
  ))
}

param = list(noise_sd = 1)

load('data/smallworld_large.Rdata')


n_reps = 200
n_cores = 56

registerDoParallel(cores=n_cores)

# Compute truth
data = generate_covariate_data(g, covariate_fns_for_response)
truth1 = foreach(i = 1:500, .combine=c) %dopar% {basic_nonlinear_response(rep(1, 5000), data$x_trt, param)}
truth0 = foreach(i = 1:500, .combine=c) %dopar% {basic_nonlinear_response(rep(0, 5000), data$x_ctrl, param)}
true_ATE = mean(truth1) - mean(truth0)
true_ATE



results = foreach(rep = 1:n_reps, .combine=rbind) %dopar% {
  run_sim(param, g, n_reps, n_cores, pid = i)
} %>% data.frame

mean(results$dm)
var(results$dm)
mean(results$adjusted)
var(results$adjusted)


# TEST VARIANCE ESTIMATE

# creates a list containing the data
covariate_fns_for_response=list(
  num_nbh = number_trt_nbrs
  #frac_nbh = fraction_trt_nbrs,
  #frac_nbh2 = fraction_trt_nbrs2
)
data = generate_covariate_data(g, covariate_fns_for_response)
data$y = basic_nonlinear_response(data$w, data$x_obs, param)
n_folds = 2
n = length(data$y)
fold_ids = sample(rep(1:n_folds, ceiling(n / n_folds))[1:n])
var_est=data %>% loess_boot(n_folds, fold_ids, g, covariate_fns_for_response, B=50)
var_est
















df = data.frame(x = data$x_obs$num_nbh, y = data$y, w = as.factor(data$w))
df %>% ggplot(aes(x, y, group = w, colour=w)) + geom_point()
fit = loess(y ~ x, data=df, span=.5)
plot(fit)
lines(predict(fit), x = df$x)

ggplot(df) + geom_point(aes(x,y)) + geom_line(data=data.frame(x=df$x, y=predict(fit)), aes(x, y))

loess_crossfit(data)
