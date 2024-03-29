.plot_lim_results = function(network, sigma) {
  if (network == 'caltech' & sigma == 1) {
    r_file = 'results/sim_lim/caltech_results_1.csv'
    t_file = 'results/sim_lim/caltech_true_ATE_1.csv'
    title = expression(paste('Caltech network, noise ', sigma, ' = 1'))
    title2 = expression(paste('Caltech, ', sigma, ' = 1'))
  } else if (network == 'caltech' & sigma == 3) {
    r_file = 'results/sim_lim/caltech_results_3.csv'
    t_file = 'results/sim_lim/caltech_true_ATE_3.csv'
    title = expression(paste('Caltech network, noise ', sigma, ' = 3'))
    title2 = expression(paste('Caltech, ', sigma, ' = 3'))
  } else if (network == 'smallworld' & sigma == 1) {
    r_file = 'results/sim_lim/smallworld_results_1.csv'
    t_file = 'results/sim_lim/smallworld_true_ATE_1.csv'
    title = expression(paste('Small-world network, noise ', sigma, ' = 1'))
    title2 = expression(paste('Small-world, ', sigma, ' = 1'))
  } else if (network == 'smallworld' & sigma == 3) {
    r_file = 'results/sim_lim/smallworld_results_3.csv'
    t_file = 'results/sim_lim/smallworld_true_ATE_3.csv'
    title = expression(paste('Small-world network, noise ', sigma, ' = 3'))
    title2 = expression(paste('Small-world, ', sigma, ' = 3'))
  } else stop()
  
  if (network == 'smallworld') {
    breaks = 1:3
    limits = c(0, 3.5)
  } else {
    breaks = 1:5
    limits = c(0, 5.5)
  }
  
  results = read.csv(r_file, header=FALSE)
  colnames(results) = c('pid', 'dm', 'hajek', 'adj1', 'var1', 'adj2', 'var2')
  
  params = purrr::cross(list(
    b_intercept = 0,
    b_direct = 1,
    b_spill = c(0, 0.25, 0.5, 0.75, 1),#0.1, 0.2, 0.3, 0.4, 0.5),
    max_t = c(2, 4),
    is_probit = FALSE,
    noise_sd = 3
  ))
  
  params_df  = params %>% unlist %>% matrix(nrow=length(params), byrow=TRUE) %>% as.data.frame
  names(params_df) = c('b_intercept', 'b_direct', 'b_spill', 'max_t', 'is_probit', 'noise_sd')
  params_df$pid = 1:nrow(params_df)
  
  truth = read.csv(t_file)
  
  #results = params_df %>% left_join(results)
  
  head(results)
  summarised = results %>% 
    select(pid, dm, hajek, adj1, adj2) %>% 
    gather(estimator, v, -pid) %>%
    group_by(pid, estimator) %>% 
    summarise(mean=mean(v, na.rm=TRUE), var=var(v, na.rm=TRUE)) %>%
    left_join(truth) %>%
    mutate(abs_bias=abs(mean-ATE), sd=sqrt(var), RMSE=sqrt(abs_bias^2 + var)) %>%
    left_join(params_df)
  
  #palette = c('adj1'="#E69F00", 'adj2'="#E69F00", 'dm'="#009E73", 'hajek' = "#56B4E9")
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73",  "#0072B2", "#D55E00", "#CC79A7")
  #"#F0E442",
  
  p_rmse = summarised %>% gather(metric, value, abs_bias, sd, RMSE) %>%
    mutate(metric = factor(metric, levels=c('abs_bias', 'sd', 'RMSE'), labels=c('absolute bias', 'standard error', 'RMSE'))) %>% # reorder
    ggplot(aes(b_spill, value, group=estimator, colour=estimator)) + geom_point(shape=0) + geom_line() +
    facet_grid(max_t ~ metric, labeller=
                 label_bquote(rows='steps T' == ~ .(max_t))
    ) + scale_color_manual(values=cbPalette) + theme_bw() +
    scale_y_continuous(limits=limits, breaks=breaks) + 
    theme(
      legend.position="bottom",
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    xlab(expression('spillover effect '*gamma)) +
    ggtitle(title)
  
  
  
  level = 0.1
  cov1 = results %>% 
    select(pid, estimate=adj1, sd=var1) %>% 
    left_join(truth) %>% 
    mutate(bias = abs(estimate - ATE), true_sd = sd(estimate), covers = abs(estimate - ATE) / sd < qnorm(1 - level/2)) %>% 
    group_by(pid) %>% 
    summarise(coverage=mean(covers), mean(sd / true_sd)) %>% 
    mutate(estimator='adj1')
  cov2 = results %>% 
    select(pid, estimate=adj2, sd=var2) %>% 
    left_join(truth) %>% 
    mutate(bias = abs(estimate - ATE), true_sd = sd(estimate), covers = abs(estimate - ATE) / sd < qnorm(1 - level/2)) %>% 
    group_by(pid) %>% 
    summarise(coverage=mean(covers), mean(sd / true_sd)) %>% 
    mutate(estimator='adj2')
  
  p_coverage = rbind(cov1, cov2) %>%
    left_join(params_df) %>%
    ggplot(aes(b_spill, coverage, group=estimator, colour=estimator)) + 
    geom_point(shape=0) + geom_line() +
    scale_y_continuous(limits=c(0, 1)) + geom_hline(yintercept = 1 - level, linetype='dashed') + 
    facet_grid(max_t ~ ., labeller=label_bquote(rows='steps T' == ~ .(max_t))) + 
    scale_color_manual(values=cbPalette) + theme_bw() +
    theme(
      legend.position="bottom",
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    xlab(expression('spillover effect '*gamma)) +
    ggtitle(title2)
  
  list(p_rmse=p_rmse, p_coverage=p_coverage)
}


p1 = .plot_lim_results('smallworld', 1)
p2 = .plot_lim_results('smallworld', 3)
p3 = .plot_lim_results('caltech', 1)
p4 = .plot_lim_results('caltech', 3)

p_all = grid.arrange(p1$p_rmse, p2$p_rmse, p3$p_rmse, p4$p_rmse, nrow=2)
p_all

ggsave(filename='figures/lim_plot.png', p_all, width=8, height=8)




p_cov = grid.arrange(p1$p_coverage, p2$p_coverage, p3$p_coverage, p4$p_coverage, nrow=1)
p_cov

ggsave(filename='figures/lim_coverage.png', p_cov, width=10, height=4)
