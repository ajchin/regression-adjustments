library(ggplot2)
library(gridExtra)

load('data/caltech.Rdata')

covariate_fns=list(
  frac_nbh = fraction_trt_nbrs,
  num_nbh = number_trt_nbrs
)

data = generate_covariate_data(g_fb, covariate_fns)

df = rbind(
  data$x_ctrl %>% mutate(cf = 'control'), # cf = counterfactual
  data$x_trt %>% mutate(cf = 'treatment'),
  data$x_obs %>% mutate(cf = 'observed')
) %>% 
  gather(k, v, frac_nbh, num_nbh)

palette = c(
  "#E69F00", "#56B4E9", "#009E73"
)

.plot = function(df, feature, legend = TRUE) {
  p = df %>% 
    filter(k == feature) %>% 
    ggplot(aes(v, group=cf, fill=cf)) + geom_histogram(bins=30, alpha=0.7, position='identity') + 
    ylab('number of units') + labs(fill = "counterfactual") + 
    scale_fill_manual(values=palette) + theme_bw()
  
  if (feature == 'frac_nbh') p = p + xlab('proportion of treated neighbors')
  if (feature == 'num_nbh') p = p + xlab('number of treated neighbors')
  
  if (legend) {
    p = p + theme(legend.position = c(0.85, 0.85))
  } else {
    p = p + theme(legend.position = 'none')
  }
  
  p
}


p = grid.arrange(
  .plot(df, 'frac_nbh', legend=FALSE), 
  .plot(df, 'num_nbh'), 
  nrow=1,
  top = 'Observed, global control, and global treatment feature distributions'
)


ggsave('figures/cf_distns.png', p, width=10, height=5)
