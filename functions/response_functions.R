
# Plain linear model with intercept, frac_nbh, num_nbh
linear_response = function(w, x, p) {
  n = length(w)
  x = cbind(1, as.matrix(x))
  
  y =  x %*% p$beta0 * (w == 0) + x %*% p$beta1 * (w == 1) + rnorm(n, sd=p$noise_sd)
  y
}




# just passes a linear model through logistic function
binary_response = function(w, x, p) {
  .logistic = function(x) {1 / (1 + exp(-x))}
  linear_response(w, x, p) %>% (function(x) {1 * (.logistic(x) > 0.5)})
}

dynamic_time_response = function(w, g, param) {
  
  b_intercept = param$b_intercept
  b_direct = param$b_direct
  b_spill = param$b_spill
  max_t = param$max_t
  is_probit = param$is_probit
  noise_sd = param$noise_sd
  
  adj = as_adj(g)
  n_peers = degree(g)
  n = length(w)
  y = rep(0, n)
  
  for (t in 1:max_t) {
    avg_nbr_y = as.vector(adj %*% y / n_peers)
    y = b_intercept + b_direct * w + b_spill * avg_nbr_y + rnorm(n, sd=noise_sd)
    if (is_probit) y = 1*(y > 0)
  }
  
  return(y)
}


# nonlinear in frac_trt_nbrs
basic_nonlinear_response = function(w, x, p) {
  n = length(w)
  
  # direct effect
  y = rnorm(n, mean=3, sd=0.5) * w
  
  # number of treated neighbors
  num_nbh = x$num_nbh
  y = y + num_nbh * 0.1 + 10 / (1 + 0.001 * exp(-1 * (num_nbh - 20)))
  
  return(y + rnorm(n, sd=p$noise_sd))
}


# nonlinear in both num_trt_nbrs and frac_trt_nbrs
nonlinear_response = function(w, x, param, noise_sd = 1) {
  n = length(w)
  
  # intercept and (heterogeneous) direct effect #TODO FIX MEAN = 1
  y = -5 + 2 * rnorm(n, mean=2, sd=2) * data$w
  
  # number of treated neighbors
  num_nbh = x$num_nbh
  y = y + num_nbh * 0.03 + 1 / (1 + 0.001 * exp(-0.03 * (num_nbh - 300)))
  
  # fraction of treated neighbors
  frac_nbh = x$frac_nbh
  y = y + 10 / (3 + exp(-8 * (frac_nbh - 0.4)))
  
  return(y + rnorm(n, sd=noise_sd))
}
