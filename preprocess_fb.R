library(R.matlab)
fb_data = readMat('data/Caltech36.mat')
g_fb = graph_from_adjacency_matrix(fb_data$A, mode='undirected')
g_fb = induced_subgraph(g_fb, which(components(g_fb)$membership == 1))
save(g_fb, file='data/caltech.Rdata')

