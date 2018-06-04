# Contains variance estimators for regression estimators


linear_variance_estimate = function(data, variance_factor, vars=NULL) {
  df = with(data, data.frame(y, x_obs %>% select(one_of(vars))))
  df = scale(df, center=TRUE, scale=FALSE) %>% as.data.frame
  n = nrow(df)
  #formula = if (is.null(vars)) 'y ~ .' else paste('y ~ ', paste(vars, collapse=' + '), sep='')
  lm_trt = lm(y ~ ., data=df %>% filter(data$w == 1))
  lm_ctrl = lm(y ~ ., data=df %>% filter(data$w == 0))
  resids_trt = lm_trt %>% augment %>% .$.resid
  resids_ctrl = lm_ctrl %>% augment %>% .$.resid
  sigma_hat_sq = (1 / n) * (sum(resids_trt^2) + sum(resids_ctrl^2))
  

  #pi = mean(data$w == 1)
  #resid_term = 1 / (n * pi * (1 - pi))
  sigma_hat_sq * (variance_factor)
  # 
  # omega = c(1, colMeans(data$x_ctrl), 1, colMeans(data$x_trt))
  # delta_term = sigma_hat_sq * t(omega) %*% Delta %*% omega
  # 
  # c(s2=sigma_hat_sq, r=resid_term, eta=eta, delta=delta_term)
  #sigma_hat_sq * (resid_term + eta + delta_term)
}

