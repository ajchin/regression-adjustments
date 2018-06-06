params = purrr::cross(list(
  b_intercept = 1,
  b_direct = 1,
  b_spill = c(0, 0.25, 0.5, 0.75, 1),#0.1, 0.2, 0.3, 0.4, 0.5),
  max_t = c(2, 4),
  is_probit = TRUE,
  noise_sd = 1 # c(1, 3)
))

param = params[[5]]
calculate_true_ATE(param, g_sm_large, NULL)

estimates = foreach(rep = 1:n_reps, .combine=rbind) %dopar% {
  run_sim(param, g_sm_large, vf1, vf2, n_reps, n_cores, pid = i)
} %>% data.frame

estimates %>%
  summarise_each(funs(mean(., na.rm = TRUE)), dm, hajek, ols1, ols2, lr1, lr2)
