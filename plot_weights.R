# Plots regression and Hajek weights.

library(igraph)
library(foreach)
library(doParallel)
library(dplyr)
library(tidyr)
library(broom)

source('functions/data_generators.R')
source('functions/covariate_functions.R')
source('functions/response_functions.R')
source('functions/existing_estimators.R')
source('functions/proposed_estimators.R')
source('functions/variance_estimators.R')
source('functions/precompute_matrices.R')

set.seed(2010)

load('data/caltech.Rdata')

palette = c("#E69F00", "#56B4E9", "#009E73")


covariate_fns = list(frac_nbh=fraction_trt_nbrs)
data = generate_covariate_data(g_fb, covariate_fns)

param = list(alpha_trt = 3, beta_trt = 0.5, alpha_ctrl = 0, beta_ctrl = 0.2)
data$y = linear_response(data$w, data$x_obs, param, noise_sd = 1)

data %>% linear_adjustment()
data %>% hajek(g_fb, threshold_var_name='frac_nbh', threshold=0.75)

# Regression weights
w = data$w
N0 = sum(w == 0)
N1 = sum(w == 1)
X = data$x_obs$frac_nbh
y0 = data$y[w==0]
y1 = data$y[w==1]
X0 = X[w==0]
X1 = X[w==1]
X0cent = X0 - mean(X0)
X1cent = X1 - mean(X1)
omega0 = apply(data$x_ctrl, 2, mean)
omega1 = apply(data$x_trt, 2, mean)
J0 = diag(N0) - matrix(1/N0, N0, N0)
J1 = diag(N1) - matrix(1/N1, N1, N1)
reg_wt0 = - 1 / N0 - t(omega0 - mean(X0)) %*% solve(t(X0cent) %*% X0cent) %*% t(X0cent) %*% J0 %>% as.vector
reg_wt1 =  1 / N1 + t(omega1 - mean(X1)) %*% solve(t(X1cent) %*% X1cent) %*% t(X1cent) %*% J1 %>% as.vector

sum(reg_wt1 * data$y[w==1]) + sum(reg_wt0 * data$y[w==0])

# Hajek weights
hajek_wts = data %>% hajek_weights(threshold_var_name='frac_nbh', g, threshold=0.75)

df = rbind(
  data.frame(frac_nbh=X, wt=hajek_wts, w=ifelse(data$w, 'treatment', 'control'), estimator='hajek'),
  data.frame(frac_nbh=X[w==0], wt=reg_wt0, w='control', estimator='regression'),
  data.frame(frac_nbh=X[w==1], wt=reg_wt1, w='treatment', estimator='regression')
)

p = df %>% ggplot(aes(frac_nbh, wt, colour=as.factor(w))) + geom_point(shape=0) + 
  geom_hline(yintercept=0, linetype='dashed') + facet_grid(. ~ estimator) + 
  geom_vline(aes(xintercept=xint), data=data.frame(xint=0.25, estimator='hajek'), linetype="dotted", colour="#E69F00") + 
  geom_vline(aes(xintercept=xint), data=data.frame(xint=0.75, estimator='hajek'), linetype="dotted", colour="#56B4E9") + 
  scale_colour_manual(values=palette) + theme_bw() + 
  theme(legend.position='bottom', plot.title = element_text(hjust = 0.5)) + labs(colour = "treatment group") + 
  ggtitle('Estimator weights for the Caltech network') + xlab('proportion of treated neighbors') + ylab('weight')
p

ggsave('figures/weights.png', p, width=6, height=4)

