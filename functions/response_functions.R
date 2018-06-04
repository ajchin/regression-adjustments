
# Plain linear model with intercept, frac_nbh, num_nbh
linear_response = function(w, x, p) {
  n = length(w)
  x = cbind(1, as.matrix(x))
  
  y =  x %*% p$beta0 * (w == 0) + x %*% p$beta1 * (w == 1) + rnorm(n, sd=p$noise_sd)
  y
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
