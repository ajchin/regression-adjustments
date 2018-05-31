
data = generate_covariate_data(g, covariate_fns_for_response)
data$y = linear_response(data$w, data$x_obs, param, noise_sd=1)

df = data.frame(y=data$y, w = data$w, data$x_obs)

df %>% ggplot(aes(frac_nbh, y, colour=as.factor(w))) + geom_point(alpha=0.1)


test = function(param, g, n_reps, n_cores) {
  # param: row of parameters
  
  # creates a list containing the data
  covariate_fns_for_response=list(
    frac_nbh = fraction_trt_nbrs,
    num_nbh = number_trt_nbrs
  )
  data = generate_covariate_data(g, covariate_fns_for_response)
  
  # generate response
  data$y = linear_response(data$w, data$x_obs, param, noise_sd=1)
  
  y = data$y
  w = data$w
  
  formula = if (is.null(vars)) 'y ~ .' else paste('y ~ ', paste(vars, collapse=' + '), sep='')
  pred_all_trt = .adjust_within_group(data$x_obs %>% filter(w == 1), y[w==1], formula, data$x_trt)
  pred_all_ctrl = .adjust_within_group(data$x_obs %>% filter(w == 0), y[w==0], formula, data$x_ctrl)
  c(mean(pred_all_trt), mean(pred_all_ctrl))

  # df = with(data, data.frame(y, x_obs))
  # lm_trt = lm(y ~ ., data=df %>% filter(data$w == 1))
  # lm_ctrl = lm(y ~ ., data=df %>% filter(data$w == 0))
  
  #c(coef(lm_ctrl), coef(lm_trt))
}


g = sample_smallworld(dim=1, size=2000, nei=5, p=0.1)

param = list(alpha_trt = 1, beta_trt = c(1, 2), alpha_ctrl = 0, beta_ctrl = c(0.5, 1))

n_reps = 500
n_cores = 48

registerDoParallel(cores=n_cores)


print('Running simulation...')
# Run simulation
start = proc.time()
estimates = foreach(rep = 1:n_reps, .combine=rbind) %dopar% {
  test(param, g, n_reps, n_cores)
} %>% data.frame
print(proc.time())

apply(estimates, 2, mean)

var(estimates$adjusted)
mean(estimates$var_est)
var(estimates$var_est)





covariate_fns = list(frac_nbh=fraction_trt_nbrs)
data = generate_covariate_data(g_fb, covariate_fns)

param = list(alpha_trt = 1, beta_trt = 3, alpha_ctrl = 0, beta_ctrl = 1)
data$y = linear_response(data$w, data$x_obs, param, noise_sd = 1)

data %>% linear_adjustment()
data %>% hajek(g_fb, threshold_var_name='frac_nbh', threshold=0.75)

# Regression weights
w = data$w
y0 = data$y[w==0]
y1 = data$y[w==1]
N0 = sum(w == 0)
N1 = sum(w == 1)
X = data$x_obs$frac_nbh
X0 = X[w==0] - mean(X[w==0])
X1 = X[w==1] - mean(X)

solve(t(X0) %*% X0) %*% t(X0) %*% (y0 - mean(y0))
solve(t(X1) %*% X1) %*% t(X1) %*% y1

lm(y ~ ., data=data.frame(y=y0, x=X[w==0]))

t(X0) %*% X0
t(X0 - mean(X0)) %*% (X0 - mean(X0))

omega0 = apply(data$x_ctrl, 2, mean)
omega1 = apply(data$x_trt, 2, mean)
reg_wt0 = - 1 / N0 - t(omega0) %*% solve(t(X0) %*% X0) %*% t(X0) %>% as.vector
reg_wt1 = 1 / N1 + t(omega1) %*% solve(t(X1) %*% X1) %*% t(X1) %>% as.vector

sum(reg_wt1 * data$y[w==1]) - sum(reg_wt0 * data$y[w==0])

