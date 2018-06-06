create_village = function(village_name, covariate_fns) {
  nodes = cai_survey %>% 
    filter(address == village_name) %>% 
    mutate(trt= delay == 0 & intensive == 1) %>%
    select(id, delay, intensive, trt, takeup_survey)
  edges = cai_all %>% 
    filter(
      address == village_name, network_address == village_name,
      id %in% nodes$id, network_id %in% nodes$id
    )
  nodes = nodes %>% filter(id %in% edges$id | id %in% edges$network_id)
  graph = edges %>% 
    select(id, network_id) %>% 
    as.matrix %>%
    apply(2, as.character) %>% 
    graph_from_edgelist(directed=FALSE) %>%
    simplify # removes bidirectional edges that would be counted twice
  
  w = 1*(nodes$trt)
  structure(
    list(
      name=village_name,
      nodes=nodes,
      edges=edges,
      graph=graph,
      n_nodes=vcount(graph),
      n_edges=ecount(graph),
      y=nodes$takeup_survey,
      w=w,
      x_obs=.build_obs_covariates(covariate_fns, graph, w),
      x_trt=.build_treated_covariates(covariate_fns, graph),
      x_ctrl=.build_control_covariates(covariate_fns, graph)
    ),
    class = 'village'
  )
}





cai_all = read.dta('../seeding/alex/cai-data/data/0422allinforawnet.dta')
cai_survey = read.dta('../seeding/alex/cai-data/data/0422survey.dta')
village_names = unique(cai_all$address)
village_names = village_names[village_names != '']
n_edges = foreach(v = village_names) %dopar% {
  edges = cai_all %>% filter(address == v, network_address == v)
  nrow(edges)
}
village_names = village_names[n_edges > 25]
n = length(village_names)

fns=list(
  frac_nbh = fraction_trt_nbrs,
  frac_nbh2 = fraction_trt_nbrs2,
  num_nbh = number_trt_nbrs
)
villages = foreach(i = 1:n) %do% create_village(village_names[i], fns)

tmp = foreach(v = villages, .combine=c) %dopar% {
  lm(v$y ~ ., data=v$x_obs) %>% glance %>% .$r.squared
}
tmp
v = villages[[1]]

fit = lm(v$y ~ ., data=v$x_obs)

# TODO
# Try and copy Dean's setup from here https://github.com/deaneckles/randomization_inference/blob/master/Randomization%20inference%20in%20networks.ipynb