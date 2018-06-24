# Contains variance estimators for regression estimators
source('functions/proposed_estimators.R')

dm_variance_estimate = function(data) {
  w = data$w
  y = data$y
  var(y[w==1]) / sum(w==1) + var(y[w==0]) / sum(w==0)
}


linear_variance_estimate = function(data, variance_factor, vars=NULL) {
  if (is.null(vars)) vars = names(data$x_obs)
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

gam_boot = function(data, n_folds, fold_ids, g, covariate_fns, B, vars=NULL) {
  
  if (is.null(vars)) vars = names(data$x_obs)
  
  w = data$w
  n = length(data$y)
  
  # train models
  models = foreach (fold = 1:n_folds) %do% {
    foreach (group = c(0, 1), .final = function(x) setNames(x, c('fit0', 'fit1'))) %do% {
      id_train = fold_ids != fold & w == group
      X_train = data$x_obs[id_train,, drop=FALSE] %>% select(one_of(vars))
      y_train = data$y[id_train]
      gam(y_train ~ ., data=X_train)
      #loess(y_train ~ ., data = X_train, control=loess.control(surface="direct"))
    }
  }
  
  predictions = foreach (fold = 1:n_folds, .combine = rbind) %do% {
    id_test = fold_ids == fold 
    obs_trt = predict(models[[fold]]$fit1, newdata=data$x_obs[id_test,,drop=FALSE])
    obs_ctrl = predict(models[[fold]]$fit0, newdata=data$x_obs[id_test,,drop=FALSE])
    data.frame(obs_trt=obs_trt, obs_ctrl=obs_ctrl, y=data$y[id_test], w=w[id_test])
  }
  
  # bootstrap variance estimate
  residuals = with(predictions, w * (y - obs_trt) + (1 - w) * (y - obs_ctrl))
  
  tau_boot = foreach(b = 1:B, .combine=c) %dopar% {
    data = generate_covariate_data(g, covariate_fns)
    resids_boot = sample(residuals, size=n, replace=TRUE)
    
    # PREDICT NEW X THEN GENERATE Y THEN COMPUTE ESTIMSATOR
    x = data$x_obs
    mu_boot = foreach (fold = 1:n_folds, .combine=rbind) %do% {
      id_test = fold_ids == fold
      trt = predict(models[[fold]]$fit1, newdata=data$x_obs[id_test,,drop=FALSE])
      ctrl = predict(models[[fold]]$fit0, newdata=data$x_obs[id_test,,drop=FALSE])
      data.frame(trt=trt, ctrl=ctrl, w = data$w[id_test])
    }
    data$y = with(mu_boot, w * trt + (1 - w) * ctrl + resids_boot)
    data$w = mu_boot$w
    gam_crossfit(data, n_folds, fold_ids)
  }
  
  var(tau_boot)
}

loess_boot = function(data, n_folds, fold_ids, g, covariate_fns, B) {
  w = data$w
  n = length(data$y)
  
  # train models
  models = foreach (fold = 1:n_folds) %do% {
    foreach (group = c(0, 1), .final = function(x) setNames(x, c('fit0', 'fit1'))) %do% {
      id_train = fold_ids != fold & w == group
      df = data.frame(y = data$y[id_train], x = data$x_obs$num_nbh[id_train])
      loess(y ~ x, data = df, control=loess.control(surface="direct"))
    }
  }
  
  predictions = foreach (fold = 1:n_folds, .combine = rbind) %do% {
    id_test = fold_ids == fold 
    obs_trt = predict(models[[fold]]$fit1, newdata=data$x_obs$num_nbh[id_test])
    obs_ctrl = predict(models[[fold]]$fit0, newdata=data$x_obs$num_nbh[id_test])
    data.frame(obs_trt=obs_trt, obs_ctrl=obs_ctrl, y=data$y[id_test], w=w[id_test])
  }
  
  # bootstrap variance estimate
  residuals = with(predictions, w * (y - obs_trt) + (1 - w) * (y - obs_ctrl))
  
  tau_boot = foreach(b = 1:B, .combine=c) %dopar% {
    data = generate_covariate_data(g, covariate_fns)
    resids_boot = sample(residuals, size=n, replace=TRUE)
    
    # PREDICT NEW X THEN GENERATE Y THEN COMPUTE ESTIMSATOR
    x = data$x_obs$num_nbh
    mu_boot = foreach (fold = 1:n_folds, .combine=rbind) %do% {
      id_test = fold_ids == fold
      trt = predict(models[[fold]]$fit1, newdata=x[id_test])
      ctrl = predict(models[[fold]]$fit0, newdata=x[id_test])
      data.frame(trt=trt, ctrl=ctrl, w = data$w[id_test])
    }
    data$y = with(mu_boot, w * trt + (1 - w) * ctrl + resids_boot)
    data$w = mu_boot$w
    loess_crossfit(data, n_folds, fold_ids)
  }
  
  var(tau_boot)
}

lr_boot = function(data, n_folds, fold_ids, g, covariate_fns, B) {
  w = data$w
  n = length(data$y)
  vars = names(data$x_obs)
  
  # train models
  models = foreach (fold = 1:n_folds) %do% {
    foreach (group = c(0, 1), .final = function(x) setNames(x, c('fit0', 'fit1'))) %do% {
      id_train = fold_ids != fold & w == group
      X_train = data$x_obs[id_train,, drop=FALSE] %>% select(one_of(vars))
      y_train = data$y[id_train]
      glm(y_train ~ ., data = X_train, family='binomial')
    }
  }
  
  # predictions = foreach (fold = 1:n_folds, .combine = rbind) %do% {
  #   id_test = fold_ids == fold 
  #   obs_trt = predict(models[[fold]]$fit1, newdata=data$x_obs[id_test,,drop=FALSE], type='response')
  #   obs_ctrl = predict(models[[fold]]$fit0, newdata=data$x_obs[id_test,,drop=FALSE], type='response')
  #   data.frame(obs_trt=obs_trt, obs_ctrl=obs_ctrl, y=data$y[id_test], w=w[id_test])
  # }
  
  # bootstrap variance estimate
  #residuals = with(predictions, w * (y - obs_trt) + (1 - w) * (y - obs_ctrl))
  
  tau_boot = foreach(b = 1:B, .combine=c) %dopar% {
    data = generate_covariate_data(g, covariate_fns)
    #resids_boot = sample(residuals, size=n, replace=TRUE)
    
    # PREDICT NEW X THEN GENERATE Y THEN COMPUTE ESTIMSATOR
    x = data$x_obs
    y_boot = foreach (fold = 1:n_folds, .combine=rbind) %do% {
      id_test = fold_ids == fold
      mu_trt = predict(models[[fold]]$fit1, newdata=x[id_test,,drop=FALSE], type='response')
      mu_ctrl = predict(models[[fold]]$fit0, newdata=x[id_test,,drop=FALSE], type='response')
      y_trt = rbinom(length(mu_trt), size=1, prob=mu_trt)
      y_ctrl = rbinom(length(mu_ctrl), size=1, prob=mu_ctrl)
      data.frame(trt=y_trt, ctrl=y_ctrl, w = data$w[id_test])
    }
    data$y = with(y_boot, w * trt + (1 - w) * ctrl)# + resids_boot)
    data$w = mu_boot$w
    lr_crossfit(data, fold_ids, n_folds)
  }
  
  var(tau_boot)
  
}
