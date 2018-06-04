results = read.csv('results/sim_lim/results_3.csv', header=FALSE)
colnames(results) = c('pid', 'dm', 'hajek', 'adj1', 'var1', 'adj2', 'var2')

params_df  = params %>% unlist %>% matrix(nrow=length(params), byrow=TRUE) %>% as.data.frame
names(params_df) = c('b_intercept', 'b_direct', 'b_spill', 'max_t', 'is_probit', 'noise_sd')
params_df$pid = 1:nrow(params_df)

truth = read.csv('results/sim_lim/true_ATE_3.csv')
  
#results = params_df %>% left_join(results)

head(results)
summarised = results %>% 
  select(pid, dm, hajek, adj1, adj2) %>% 
  gather(estimator, v, -pid) %>%
  group_by(pid, estimator) %>% 
  summarise(mean=mean(v, na.rm=TRUE), var=var(v, na.rm=TRUE)) %>%
  left_join(truth) %>%
  mutate(abs_bias=abs(mean-ATE), RMSE=sqrt(abs_bias^2 + var), sd=sqrt(var)) %>%
  left_join(params_df)

cbPalette <- c("#E69F00", "#56B4E9", "#009E73",  "#0072B2", "#D55E00", "#CC79A7")
#"#F0E442",

summarised %>% filter(noise_sd == 1) %>% gather(metric, value, abs_bias, sd, RMSE) %>%
  ggplot(aes(b_spill, value, group=estimator, colour=estimator)) + geom_point(shape=0) + geom_line() +
  facet_grid(max_t ~ metric, scales='free_y') + scale_color_manual(values=cbPalette) + theme_bw()

