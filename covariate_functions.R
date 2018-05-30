# Defines commonly used covariate generators

# Fraction of treated neighbors
fraction_trt_nbrs = function(g, w) {as.vector(as_adj(g) %*% w) / degree(g)}

# Number of treated neighbors
number_trt_nbrs = function(g, w) {as.vector(as_adj(g) %*% w)}

# Fraction of treated individuals in the 2-hop neighborhood.
fraction_trt_nbrs2 = function(g, w) {
  adj = as_adj(g)
  adj2 = adj %*% adj
  as.vector(adj2 %*% w) / apply(adj2, 2, sum) 
}