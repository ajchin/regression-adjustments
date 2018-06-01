
tmp =  foreach(i = 1:500, .combine=rbind) %dopar% {
  data = generate_covariate_data(g, covariate_fns=list(frac_nbh=fraction_trt_nbrs, num_nbh=number_trt_nbrs))
  data$y = linear_response(data$w, data$x_obs, param, noise_sd=3)
  y1 = data$y[data$w==1]
  x1 = data$x_obs %>% filter(data$w == 1) %>% as.matrix
  beta1=lm(y1 ~ x1) %>% coef
  omega1 = c(1, colMeans(data$x_trt))
  
  y0 = data$y[data$w==0]
  x0 = data$x_obs %>% filter(data$w == 0) %>% as.matrix

  beta0 = lm(y0 ~ x0) %>% coef
  omega0 = c(1, colMeans(data$x_ctrl))
  
  pred = sum(beta1 * omega1) - sum(beta0 * omega0)
  
  
  s1 = solve(t(cbind(1, x1)) %*% cbind(1, x1))
  s0 = solve(t(cbind(1, x0)) %*% cbind(1, x0))
  v1 = t(omega1) %*% s1 %*% omega1
  v0 = t(omega0) %*% s0 %*% omega0
  c(beta0, beta1, v0, v1, pred)
}
