library(ggplot2)
library(GGally)

set.seed(2018)

cai_all = read.dta('cai-data/data/0422allinforawnet.dta')
cai_survey = read.dta('cai-data/data/0422survey.dta')

cai_edgelist = cai_all %>% 
  filter(!is.na(network_id), network_id != 99, id %in% cai_survey$id, network_id %in% cai_survey$id) %>%
  select(id, network_id)

ids = with(cai_edgelist, unique(c(id, network_id)))

g = cai_edgelist %>% 
  as.matrix %>%
  apply(2, as.character) %>%
  graph_from_edgelist(directed=FALSE) # TODO directed

attributes = cai_survey %>% filter(id %in% ids) %>% select(id, takeup_survey, delay, intensive)

id_matcher = match(V(g)$name, attributes$id)
V(g)$y = attributes$takeup_survey[id_matcher]
V(g)$delay = attributes$delay[id_matcher]
V(g)$intensive = attributes$intensive[id_matcher]
V(g)$w = 1*(V(g)$intensive == 1)

y = V(g)$y
w = V(g)$w

covariate_fns=list(
  frac1 = fraction_trt_nbrs,
  frac2 = fraction_trt_nbrs2,
  num1 = number_trt_nbrs,
  num2 = number_trt_nbrs2
)

x_obs=.build_obs_covariates(covariate_fns, g, w)
x_trt=.build_treated_covariates(covariate_fns, g)
x_ctrl=.build_control_covariates(covariate_fns, g)



# point estimate: exploratory
fit = glm(y ~ ., data=data.frame(w, x_obs), family='binomial')

pred_obs = predict(fit, newdata=data.frame(w, x_obs), type='response')
mean(pred_obs[w==1]) - mean(pred_obs[w==0])


data = list(y=y, w=w, x_obs=x_obs, x_trt=x_trt, x_ctrl=x_ctrl)


vf = precompute_variance(g, covariate_fns, n_boot_reps=100, n_cores=n_cores)

df = data.frame(estimator=c('dm', 'hajek1', 'hajek2', 'linear', 'logistic'), 
estimate=c(
  data %>% difference_in_means,
  data %>% hajek(g, 'frac1', threshold=0.75),
  data %>% hajek(g, 'frac2', threshold=0.75),
  data %>% linear_adjustment,
  data %>% lr_crossfit(n_folds=5)
),
se = c(
  0, 
  0, 
  0,
  data %>%linear_variance_estimate(vf) %>% sqrt, 
  0)
) %>% xtable(digits=4) %>% print.xtable(include.rownames=FALSE)





p_scatter = data.frame(y=as.factor(y), x_obs) %>% ggpairs(lower=list(continuous=wrap('points', alpha=0.03, size=0.75))) + theme_bw() + 
  ggtitle('Scatterplot matrix for Cai et al. (2015) variables') + theme(plot.title = element_text(hjust = 0.5))
ggsave('figures/cai_scatter.png', p_scatter, width=8, height=6)
p_scatter




