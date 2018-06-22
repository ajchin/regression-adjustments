library(ggplot2)
library(gridExtra)
load('data/stanford.Rdata')
fns=list(
  num_nbh = number_trt_nbrs,
  frac_nbh = fraction_trt_nbrs
  #frac_nbh2 = fraction_trt_nbrs2
)
data = generate_covariate_data(g, fns)

# generate response
data$y = nonlinear_response(data$w, data$x_obs, param)

frac = data$x_obs$frac_nbh
num = data$x_obs$num_nbh
y = data$y

p_num = data.frame(y, num) %>% 
  ggplot(aes(num, y)) + geom_point(alpha = 0.03, shape=4) + scale_x_log10() + 
  theme_bw() + xlab('number of treated neighbors')
p_frac = data.frame(y, frac) %>% 
  ggplot(aes(frac, y)) + geom_point(alpha = 0.03, shape=4) +
  theme_bw() + xlab('fraction of treated neighbors')
p = grid.arrange(p_frac, p_num, nrow=1, top='Single instance of nonlinear response')
ggsave('figures/nonlinear_response.png', p, width=10, height=5)
