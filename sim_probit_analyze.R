
results = read.csv('results/sim_lim/results_lr.csv', header=FALSE)
names(results) = c('pid', 'dm', 'hajek', 'ols1', 'ols2', 'lr1', 'lr2')
truth = read.csv('results/sim_lim/true_ATE_probit.csv')

true_ATE
results %>%
  group_by(pid) %>%
  summarise_each(funs(mean(., na.rm = TRUE)), dm, hajek, ols1, ols2, lr1, lr2) %>%
  left_join(true_ATE) %>% gather(k, v, dm, hajek, ols1, ols2, lr1, lr2) %>%
  mutate(bias_sq = (ATE - v)^2) %>% ggplot(aes(pid, bias_sq, group=k, colour=k)) + geom_point() + geom_line() + theme_bw()
results %>%
  group_by(pid) %>%
  summarise_each(funs(var(., na.rm = TRUE)), dm, hajek, ols1, ols2, lr1, lr2)

# 