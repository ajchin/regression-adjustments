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

df = data.frame(frac, num, y)

p_x = df %>% ggplot(aes(num, frac)) + 
  geom_point(alpha=0.05, shape=4) + scale_x_log10() + theme_bw() +
  xlab('number of treated neighbors') + ylab('fraction of treated neighbors')
p_num = df %>%
  ggplot(aes(num, y)) + geom_point(alpha = 0.05, shape=4) + scale_x_log10() + 
  theme_bw() + xlab('number of treated neighbors')
p_frac = df %>%
  ggplot(aes(frac, y)) + geom_point(alpha = 0.05, shape=4) +
  theme_bw() + xlab('fraction of treated neighbors')
p = grid.arrange(p_x, p_num,  p_frac, nrow=1, 
                 top='Marginal plots for single instance of nonlinear response on Stanford network')
p
ggsave('figures/nonlinear_response.png', p, width=10, height=4)
