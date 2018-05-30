# Functions for creating necessary data frames

.build_obs_covariates = function(covariate_fns, g, w) {
  x_obs = lapply(covariate_fns, function(f) {f(g, w)})
  do.call(data.frame, x_obs)
}

.build_counterfactual_covariates = function(covariate_fns, g, w){
  x = lapply(covariate_fns, function(f) {f(g, w)})
  #names(z) = sapply(names(z), function(name) {paste(name, suffix, sep='')})
  do.call(data.frame, x)
}

.build_treated_covariates = function(covariate_fns, g) {
  .build_counterfactual_covariates(covariate_fns, g, rep(1, vcount(g)))#, '_T')
}

.build_control_covariates = function(covariate_fns, g) {
  .build_counterfactual_covariates(covariate_fns, g, rep(0, vcount(g)))#, '_C')
}

# Creates a data frame with treatment vector w, and three columns for each 
# specified covariate function: observed, full treatment, and full control.
# Does not include response.

# Creates a list with four elements:
#   w: a vector of treatments
#   z_obs: a data frame of observed covariates
#   z_trt: a data frame of covariate values when all are treated
#   z_ctrl: a data frame of covariate values when all are in control
generate_covariate_data = function(g, covariate_fns) {
  # g is the graph
  # noise is a vector of length vcount(g)
  # covariate functions is a list of functions of the form f(g, w)
  w = rbinom(vcount(g), size=1, prob=0.5)
  x_observed = .build_obs_covariates(covariate_fns, g, w) 
  x_global_trt = .build_treated_covariates(covariate_fns, g)
  x_global_ctrl = .build_control_covariates(covariate_fns, g)
  list(w=w, x_obs=x_observed, x_trt=x_global_trt, x_ctrl=x_global_ctrl)
}