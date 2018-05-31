# Computes Gamma and Delta for variance estimation


# Gamma is comprised of the variance/covariance of sample means
.compute_gamma = function(moments) {
  B = length(moments)
  xbar0 = foreach(m = moments) %do% {m$xbar0}
  xbar1 = foreach(m = moments) %do% {m$xbar1}
  Exbar0 = Reduce(`+`, xbar0) / B
  Exbar1 = Reduce(`+`, xbar1) / B
  
  # manual variance/covariance functions
  XXt0 = foreach(m = moments) %do% {x = m$xbar0; (x - Exbar0) %*% t(x - Exbar0)}
  XXt1 = foreach(m = moments) %do% {y = m$xbar1; (y - Exbar1) %*% t(y - Exbar1)}
  XXt01 = foreach(m = moments) %do% {x = m$xbar0; y = m$xbar1; (x - Exbar0) %*% t(y - Exbar1)}
  
  Gamma00 = Reduce(`+`, XXt0) / B
  Gamma11 = Reduce(`+`, XXt1) / B
  Gamma01 = -Reduce(`+`, XXt01) / B
  
  Gamma = rbind(
    cbind(Gamma00, Gamma01),
    cbind(Gamma01, Gamma11)
  )  
  
  colnames(Gamma) = rep(colnames(Gamma00), 2)
  rownames(Gamma) = rep(colnames(Gamma00), 2)
  Gamma
}

# Delta is comprised of the average of inverse sample covariance matrices
.compute_delta = function(moments) {
  B = length(moments)
  
  S0 = foreach(m = moments) %do% {solve(m$S0)}
  S1 = foreach(m = moments) %do% {solve(m$S1)}
  Delta00 = Reduce(`+`, S0) / B
  Delta11 = Reduce(`+`, S1) / B
  
  Delta = rbind(
    cbind(Delta00, matrix(rep(0, 4), nrow=2)),
    cbind(matrix(rep(0, 4), nrow=2), Delta11)
  )
  colnames(Delta) = rep(colnames(Delta00), 2)
  rownames(Delta) = rep(colnames(Delta00), 2)
  Delta
}

# Compute Gamma and Delta by Monte Carlo
precompute_matrices = function(g, covariate_fns, n_boot_reps, n_cores) {
  registerDoParallel(cores=n_cores) 
  moments = foreach(i = 1:n_boot_reps) %dopar% {
    data = generate_covariate_data(g, covariate_fns)
    x1 = data$x_obs %>% filter(data$w == 1) %>% as.matrix
    x0 = data$x_obs %>% filter(data$w == 0) %>% as.matrix
    
    # return sample means and covariances
    list(
      xbar1 = apply(x1, 2, mean),
      xbar0 = apply(x0, 2, mean),
      S1 = t(x1) %*% x1,
      S0 = t(x0) %*% x0
    )
  }
  
  Gamma = .compute_gamma(moments)
  Delta = .compute_delta(moments)
  
  list(Gamma=Gamma, Delta=Delta)
}
