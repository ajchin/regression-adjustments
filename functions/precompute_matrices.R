# Computes Gamma and Delta for variance estimation


# # # Gamma is comprised of the variance/covariance of sample means
# # .compute_gamma = function(moments) {
# #   B = length(moments)
# #   xbar0 = foreach(m = moments) %do% {m$xbar0}
# #   xbar1 = foreach(m = moments) %do% {m$xbar1}
# #   Exbar0 = Reduce(`+`, xbar0) / B
# #   Exbar1 = Reduce(`+`, xbar1) / B
# #   
# #   # manual variance/covariance functions
# #   XXt0 = foreach(m = moments) %do% {x = m$xbar0; (x - Exbar0) %*% t(x - Exbar0)}
# #   XXt1 = foreach(m = moments) %do% {y = m$xbar1; (y - Exbar1) %*% t(y - Exbar1)}
# #   XXt01 = foreach(m = moments) %do% {x = m$xbar0; y = m$xbar1; (x - Exbar0) %*% t(y - Exbar1)}
# #   
# #   Gamma00 = Reduce(`+`, XXt0) / B
# #   Gamma11 = Reduce(`+`, XXt1) / B
# #   Gamma01 = -Reduce(`+`, XXt01) / B
# #   
# #   Gamma = rbind(
# #     cbind(Gamma00, Gamma01),
# #     cbind(Gamma01, Gamma11)
# #   )  
# #   
# #   colnames(Gamma) = rep(colnames(Gamma00), 2)
# #   rownames(Gamma) = rep(colnames(Gamma00), 2)
# #   Gamma
# # }
# 
# .compute_eta = function(moments) {
#   B = length(moments)
#   
#   eta0 = foreach(m = moments) %do% {m$eta0}
#   eta1 = foreach(m = moments) %do% {m$eta1}
#   eta0 = Reduce(`+`, eta0) / B
#   eta1 = Reduce(`+`, eta1) / B
#   as.vector(eta0 + eta1)
# }
# 
# # Delta is comprised of the average of inverse sample covariance matrices
# .compute_delta = function(moments) {
#   B = length(moments)
#   
#   Delta0 = foreach(m = moments) %do% {m$Delta0}
#   Delta1 = foreach(m = moments) %do% {m$Delta1}
#   Delta00 = Reduce(`+`, Delta0) / B
#   Delta11 = Reduce(`+`, Delta1) / B
#   
#   Delta = rbind(
#     cbind(Delta00, 0 * Delta00),
#     cbind(0 * Delta00, Delta11)
#   )
#   colnames(Delta) = rep(colnames(Delta00), 2)
#   rownames(Delta) = rep(colnames(Delta00), 2)
#   Delta
# }
# 
# # Compute Gamma and Delta by Monte Carlo
# precompute_matrices = function(g, covariate_fns, n_boot_reps, n_cores) {
#   registerDoParallel(cores=n_cores) 
#   moments = foreach(i = 1:n_boot_reps) %dopar% {
#     data = generate_covariate_data(g, covariate_fns)
#     x1 = data$x_obs %>% filter(data$w == 1) %>% as.matrix
#     x0 = data$x_obs %>% filter(data$w == 0) %>% as.matrix
#     
#     #xbar1 = apply(x1, 2, mean)
#     #xbar1_mat = matrix(rep(xbar1, each=nrow(x1)), nrow=nrow(x1))
#     #xbar0 = apply(x0, 2, mean)
#     #xbar0_mat = matrix(rep(xbar0, each=nrow(x0)), nrow=nrow(x0))
#     x1 = cbind(1, x1)
#     x0 = cbind(1, x0)
#     #S1inv = solve(t(x1 - xbar1_mat) %*% (x1 - xbar1_mat))
#     #S0inv = solve(t(x0 - xbar0_mat) %*% (x0 - xbar0_mat))
#     S1inv = solve(t(x1) %*% x1) 
#     S0inv = solve(t(x0) %*% x0)
#     # return sample means and covariances
#     list(
#       #eta1 = t(xbar1) %*% S1inv %*% xbar1,
#       #eta0 = t(xbar0) %*% S0inv %*% xbar0,
#       Delta1 = S1inv,
#       Delta0 = S0inv
#     )
#   }
#   
#  # Gamma = .compute_gamma(moments)
#  # eta = .compute_eta(moments)
#   Delta = .compute_delta(moments)
#   
#   list(eta=NULL, Delta=Delta)
# }

precompute_variance = function(g, covariate_fns, n_boot_reps, n_cores) {
  registerDoParallel(cores=n_cores) 
  moments = foreach(i = 1:n_boot_reps) %dopar% {
    data = generate_covariate_data(g, covariate_fns)
    w = data$w
    
    x0 = cbind(1, data$x_obs[w==0,]) %>% as.matrix
    x1 = cbind(1, data$x_obs[w==1,]) %>% as.matrix
    
    # return sample covariances
    list(
      S1inv = solve(t(x1) %*% x1),
      S0inv = solve(t(x0) %*% x0)
    )
  }
  
  B = length(moments)

  S0inv = foreach(m = moments) %do% {m$S0inv} %>% (function(M) {Reduce(`+`, M) / B})
  S1inv = foreach(m = moments) %do% {m$S1inv} %>% (function(M) {Reduce(`+`, M) / B})

  
  data = generate_covariate_data(g, covariate_fns)
  omega0 = c(1, colMeans(data$x_ctrl))
  omega1 = c(1, colMeans(data$x_trt))
  
  (t(omega0) %*% S0inv %*% omega0 + t(omega1) %*% S1inv %*% omega1) %>% as.vector
}
