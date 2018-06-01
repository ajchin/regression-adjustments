# Contains code for all proposed adjusted estimators


linear_adjustment = function(data, vars=NULL) {
  w = data$w
  y0 = data$y[w==0]
  y1 = data$y[w==1]
  x0 = data$x_obs[w==0,]
  x1 = data$x_obs[w==1,]
  
  beta0 = lm(y0 ~ ., data=x0) %>% coef
  beta1 = lm(y1 ~ ., data=x1) %>% coef
  omega0 = c(1, colMeans(data$x_ctrl))
  omega1 = c(1, colMeans(data$x_trt))
  
  sum(omega1 * beta1) - sum(omega0 * beta0)
  
  #formula = if (is.null(vars)) 'y ~ .' else paste('y ~ ', paste(vars, collapse=' + '), sep='')
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
