# Contains variance estimators for regression estimators


linear_variance_estimate = function(data, Gamma, Delta, vars=NULL) {
  df = with(data, data.frame(y, x_obs))
  df = scale(df, center=TRUE, scale=FALSE) %>% as.data.frame
  n = nrow(df)
  formula = if (is.null(vars)) 'y ~ .' else paste('y ~ ', paste(vars, collapse=' + '), sep='')
  lm_trt = lm(formula, data=df %>% filter(data$w == 1))
  lm_ctrl = lm(formula, data=df %>% filter(data$w == 0))
  resids_trt = lm_trt %>% augment %>% .$.resid
  resids_ctrl = lm_ctrl %>% augment %>% .$.resid
  sigma_hat_sq = (1 / n) * (sum(resids_trt^2) + sum(resids_ctrl^2))
  
  pi_trt = sum(data$w == 1)
  pi_ctrl = sum(data$w == 0)
  resid_term = (1 / n) * (pi_trt / pi_ctrl + pi_ctrl / pi_trt)
  
  beta_trt = lm_trt %>% tidy %>% .$estimate %>% .[-1]
  beta_ctrl = lm_ctrl %>% tidy %>% .$estimate %>% .[-1]
  beta = c(beta_ctrl, beta_trt)
  gamma_term = t(beta) %*% Gamma %*% beta

  omega = c(apply(data$x_ctrl, 2, mean), apply(data$x_trt, 2, mean))
  delta_term = sigma_hat_sq * t(omega) %*% Delta %*% omega
  
  resid_term + gamma_term + delta_term
}
