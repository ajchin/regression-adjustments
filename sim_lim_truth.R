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


calculate_true_ATE = function(param, g, pid) {
  y1 = dynamic_time_response(w=rep(1, vcount(g)), g, param)
  y0 = dynamic_time_response(w=rep(0, vcount(g)), g, param)
  c(pid=pid, ATE=mean(y1) - mean(y0))
}


#load('data/caltech.Rdata')

load('data/smallworld.Rdata')
set.seed(2018)

n_reps = 5000
n_cores = 50

registerDoParallel(cores=n_cores)

params = purrr::cross(list(
  b_intercept = 0,
  b_direct = 1,
  b_spill = c(0, 0.25, 0.5, 0.75, 1),#0.1, 0.2, 0.3, 0.4, 0.5),
  max_t = c(2, 4),
  is_probit = FALSE,
  noise_sd = c(3)
))

print('Calculating true ATE...')
# Calculate true ATE by simulation
true_ATE = foreach(i = 1:length(params), .combine=rbind) %do% {
  param = params[[i]]
  print(unlist(param))
  foreach(rep = 1:n_reps, .combine=rbind, .inorder=FALSE) %dopar% {
    calculate_true_ATE(param, g_sm, pid=i)
  }
} %>% 
data.frame %>%
group_by(pid) %>%
summarise(ATE=mean(ATE))
print(proc.time())

write.csv(true_ATE, file='results/sim_lim/true_ATE_3.csv', row.names=FALSE)
