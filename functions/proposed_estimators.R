# Contains code for all proposed adjusted estimators


linear_adjustment = function(data, vars=NULL) {
  if (is.null(vars)) vars = names(data$x_obs)
  w = data$w
  y0 = data$y[w==0]
  y1 = data$y[w==1]
  x0 = data$x_obs %>% select(one_of(vars)) %>% filter(w==0)
  x1 = data$x_obs %>% select(one_of(vars)) %>% filter(w==1)
  
  beta0 = lm(y0 ~ ., data=x0) %>% coef
  beta1 = lm(y1 ~ ., data=x1) %>% coef
  omega0 = c(1, colMeans(data$x_ctrl %>% select(one_of(vars))))
  omega1 = c(1, colMeans(data$x_trt %>% select(one_of(vars))))
  
  sum(omega1 * beta1) - sum(omega0 * beta0)
  
  #formula = if (is.null(vars)) 'y ~ .' else paste('y ~ ', paste(vars, collapse=' + '), sep='')
}


loess_crossfit = function(data) {
  w = data$w
  n = length(data$y)
  n_folds = 2
  fold_ids = sample(rep(1:n_folds, ceiling(n / n_folds))[1:n])
  
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
    glo_trt = predict(models[[fold]]$fit1, newdata=data$x_trt$num_nbh[id_test])
    glo_ctrl = predict(models[[fold]]$fit0,  newdata=data$x_ctrl$num_nbh[id_test])
    obs_trt = predict(models[[fold]]$fit1, newdata=data$x_obs$num_nbh[id_test])
    obs_ctrl = predict(models[[fold]]$fit0, newdata=data$x_obs$num_nbh[id_test])
    data.frame(glo_trt=glo_trt, glo_ctrl=glo_ctrl, obs_trt=obs_trt, obs_ctrl=obs_ctrl)
  }
  
  dm = mean(data$y[w==1]) - mean(data$y[w==0])
  adjust_glo = with(predictions, mean(glo_trt) - mean(glo_ctrl))
  adjust_obs = with(predictions, mean(obs_trt[w==1] - mean(obs_ctrl[w==0])))
  
  dm + adjust_glo - adjust_obs
}

lr_crossfit = function(data, n_folds = 3, vars = NULL) {
  
  if (is.null(vars)) vars = names(data$x_obs)
  
  w = data$w
  
  # create fold ids
  n = length(data$y)
  fold_ids = sample(rep(1:n_folds, ceiling(n / n_folds))[1:n])
  
  # train models
  models = foreach (fold = 1:n_folds) %do% {
    foreach (group = c(0, 1), .final = function(x) setNames(x, c('fit0', 'fit1'))) %do% {
      id_train = fold_ids != fold & w == group
      X_train = data$x_obs[id_train,, drop=FALSE] %>% select(one_of(vars))
      y_train = data$y[id_train]
      glm(y_train ~ ., data = X_train, family='binomial')
    }
  }
  
  predictions = foreach (fold = 1:n_folds, .combine = rbind) %do% {
    id_test = fold_ids == fold 
    glo_trt = 1 * (predict(models[[fold]]$fit1, newdata=data$x_trt[id_test,,drop=FALSE], type='response') > 0.5)
    glo_ctrl = 1 * (predict(models[[fold]]$fit0, newdata=data$x_ctrl[id_test,,drop=FALSE], type='response') < 0.5)
    obs_trt = 1 * (predict(models[[fold]]$fit1, newdata=data$x_obs[id_test,,drop=FALSE], type='response') > 0.5)
    obs_ctrl = 1 * (predict(models[[fold]]$fit0, newdata=data$x_obs[id_test,,drop=FALSE], type='response') < 0.5)
    data.frame(glo_trt=glo_trt, glo_ctrl=glo_ctrl, obs_trt=obs_trt, obs_ctrl=obs_ctrl)
  }
  
  dm = mean(data$y[w==1]) - mean(data$y[w==0])
  adjust_glo = with(predictions, mean(glo_trt) - mean(glo_ctrl))
  adjust_obs = with(predictions, mean(obs_trt[w==1] - mean(obs_ctrl[w==0])))
  
  dm + adjust_glo - adjust_obs
}






# fit bart within specific fold and group
.bart_fit_within_group = function(data, fold_id, fold, group, n_cores) {
  # fold is an integer from 1 to n_folds
  # fold_id is the global list of fold assignments
  # group is 0 (ctrl) or 1 (trt)
  if (!(group %in% c(0, 1))) stop('group must be 0 or 1')
  
  x_train = data$x_obs %>% filter(fold_id != fold, data$w == group)
  y_train = data$y[fold_id != fold & data$w == group]
  if (n_cores == 1) {
    fit = wbart(x_train, y_train)
  } else {
    fit = mc.wbart(x_train, y_train, mc.cores=n_cores)
  }
  x_test = if (group == 1) data$x_trt else data$x_ctrl
  x_test = x_test %>% filter(fold_id == fold) %>% as.matrix
  predict(fit, newdata=x_test, mc.cores=n_cores) %>% apply(2, mean)
}

bart_crossfit = function(data, n_folds=5, n_cores=32) {
  # create fold ids
  n = length(data$y)
  fold_id = sample(rep(1:n_folds, ceiling(n / n_folds))[1:n])
  #folds = createFolds(data$y, k=n_folds, list=TRUE, returnTrain=TRUE)
  
  pred_all_trt = foreach(fold = 1:n_folds, .combine=c) %do% {
    .bart_fit_within_group(data, fold_id, fold, group=1, n_cores=n_cores)
  }
  
  pred_all_ctrl = foreach(fold = 1:n_folds, .combine=c) %do% {
    .bart_fit_within_group(data, fold_id, fold, group=0, n_cores=n_cores)
  }
  
  mean(pred_all_trt) - mean(pred_all_ctrl)
}

