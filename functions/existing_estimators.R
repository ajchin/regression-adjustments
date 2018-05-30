# Contains all functions needed for producing Hajek and DM estimators

difference_in_means = function(data) {
  # Basic DM estimator.  Requires list with entries y and w
  with(data, mean(y[w==1]) - mean(y[w==0]))
}

# For now, handles only x in [0, 1] (i.e. a fraction)
# Exposure probabilities are only for Bernoulli randomization
.indiv_hajek_weight = function(w, x, d, thresh, p_design=0.5) {
  if (w == 0) {
    if (x > 1 - thresh) return(0)
    
    # compute exposure probability 
    p = (1 - p_design) * pbinom(floor(d * (1 - thresh)), size=d, prob=p_design)
    return(1 / p)
  }
  
  # w = 1
  if (x < thresh) return(0)
  p = p_design * (1 - pbinom(floor(d * thresh), size=d, prob=p_design))
  return(1 / p)
}

# Return vector of Hajek weights.  We expose this function because sometimes we want to 
# examine/plot the weights.
hajek_weights = function(data, threshold_var_name, threshold, p_design=0.5) {
  # threshold_var: string name of variable to threshold on.
  #    Assumes the all-treated column has the same name but with '_T' suffix.
  
  w = data$w
  threshold_var = data$x_obs[, threshold_var_name]
  threshold_var_trt = data$x_trt[, threshold_var_name]
  
  # these are unnormalized weights
  wts = sapply(1:nrow(data$x_obs), function(i) {
    .indiv_hajek_weight(w[i], threshold_var[i], threshold_var_trt[i], threshold, p_design)
  })
  
  # now normalize within group
  sum_trt_wt = sum(wts[w == 1])
  sum_ctrl_wt = sum(wts[w == 0])
  return(wts / sum_trt_wt * (w == 1) - wts / sum_ctrl_wt * (w == 0))
}

hajek = function(data, threshold_var_name, threshold, p_design=0.5) {
  wts = hajek_weights(data, threshold_var_name, threshold, p_design)
  return(sum(data$y * wts))
}


# TODO: FIX HAJEK SOMETIMES NaN
# ALLOW FOR HAJEK FOR NUMBER OF NEIGHBORS