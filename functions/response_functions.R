
# Plain linear model
linear_response = function(w, x, param, noise_sd = 1) {
  n = length(w)
  x = as.matrix(x)
  
  x_trt = data$x_obs %>% filter(data$w == 1)
  x_ctrl = data$x_ctrl %>% filter(data$w == 0)
  
  y = (param$alpha_ctrl + x %*% param$beta_ctrl) * (w == 0) + 
    (param$alpha_trt +x %*% param$beta_trt) * (w == 1) + rnorm(n, sd=noise_sd)
  return(y)
}

# a la Manski; rom JCI paper
dynamic_time_response = function(w, g, param, noise_sd = 1) {
  
  b_intercept = param$b_intercept
  b_direct = param$b_direct
  b_spill = param$b_spill
  max_t = param$max_t
  is_probit = param$is_probit
  
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


# nonlinear in both num_trt_nbrs and frac_trt_nbrs
nonlinear_response = function(w, x, param, noise_sd = 1) {
  n = length(w)
  
  # intercept and (heterogeneous) direct effect #TODO FIX MEAN = 1
  y = -5 + 2 * rnorm(n, mean=0, sd=2) * data$w
  
  # number of treated neighbors
  num_nbh = x$num_nbh
  y = y + num_nbh * 0.03 + 2 / (1 + 0.001 * exp(-0.03 * (num_nbh - 300)))
  
  # fraction of treated neighbors
  frac_nbh = x$frac_nbh
  y = y + 25 / (3 + exp(-8 * (frac_nbh - 0.4)))
  
  return(y + rnorm(n, sd=noise_sd))
}
