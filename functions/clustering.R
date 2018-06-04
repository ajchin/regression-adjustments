library(R.matlab)
library(dplyr)
library(doParallel)
library(igraph)

.compute_focal_nodes = function(g, order){
  focal_nodes = c()
  graph = g
  V(graph)$name = V(graph) # name the vertices so the ids are not lost when subgraphing
  while (vcount(graph) > 0) {
    v = sample(V(graph), 1) # v is the index
    focal_nodes = c(focal_nodes, V(graph)[v]$name) # append the actual name
    nbhd = igraph::neighborhood(graph, order=order, nodes=v)[[1]]
    graph = delete_vertices(graph, nbhd)
  }
  focal_nodes
}


.get_distance = function(g, v, focal_node) {
  if (v == focal_node) return(0)
  dist = shortest_paths(g, from=v, to=focal_node)$vpath[[1]] %>% length
  if (dist == 0) return(NA)
  return(dist)
}

.get_closest_focal_node = function(g, v, focal_nodes) {
  dists = sapply(focal_nodes, function(x) {.get_distance(g, v, x)})
  focal_nodes[which.min(dists)]
}

epsilon_net = function(g, order) {
  focal_nodes = .compute_focal_nodes(g, order=order)
  # assign initial cluster centers
  V(g)[focal_nodes]$cluster = focal_nodes
  order = 1
  while(sum(is.na(V(g)$cluster)) > 0) {
    nbhds = igraph::neighborhood(g, order=order, nodes=focal_nodes)
    for (i in sample(1:length(focal_nodes))) {
      nbhd = nbhds[[i]]
      V(g)[nbhd[is.na(nbhd$cluster)]]$cluster = focal_nodes[i]
    }
    n_assigned = sum(!is.na(V(g)$cluster))
    print(sprintf('order %d: assigned %d, not assigned %d', order, n_assigned, vcount(g) - n_assigned))
    order = order + 1
  }
  
  V(g)$cluster
}
