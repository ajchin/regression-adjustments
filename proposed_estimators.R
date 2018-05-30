# Contains code for all proposed adjusted estimators


.adjust_within_group = function(z, y, formula, pred_z) {
  z_scaled = scale(z, center=TRUE, scale=FALSE)
  shift = attr(z_scaled, 'scaled:center')
  df = data.frame(y, z_scaled %>% data.frame)
  fit = lm(formula, df)
  predict(fit, newdata=pred_z - shift)
}

linear_adjustment = function(data, vars=NULL) {
  y = data$y
  w = data$w
  
  formula = if (is.null(vars)) 'y ~ .' else paste('y ~ ', paste(vars, collapse=' + '), sep='')
  pred_all_trt = .adjust_within_group(data$z_obs %>% filter(w == 1), y[w==1], formula, data$z_trt)
  pred_all_ctrl = .adjust_within_group(data$z_obs %>% filter(w == 0), y[w==0], formula, data$z_ctrl)
  mean(pred_all_trt) - mean(pred_all_ctrl)
}


# fit bart within specific fold and group
.bart_fit_within_group = function(data, fold_id, fold, group, n_cores) {
  # fold is an integer from 1 to n_folds
  # fold_id is the global list of fold assignments
  # group is 0 (ctrl) or 1 (trt)
  if (!(group %in% c(0, 1))) stop('group must be 0 or 1')
  
  x_train = data$z_obs %>% filter(fold_id != fold, data$w == group)
  y_train = data$y[fold_id != fold & data$w == group]
  if (n_cores == 1) {
    fit = wbart(x_train, y_train)
  } else {
    fit = mc.wbart(x_train, y_train, mc.cores=n_cores)
  }
  x_test = if (group == 1) data$z_trt else data$z_ctrl
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
