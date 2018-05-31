library(igraph)
library(foreign)
library(foreach)
library(R.matlab)
library(ggplot2)
library(gridExtra)
library(intergraph)
library(GGally)
library(network)

cai_all = read.dta('~/src/GitHub/seeding/alex/cai-data/data/0422allinforawnet.dta')
source('~/src/GitHub/seeding/alex/simulations/village.R')
village_names = unique(cai_all$address)[1:20]

g_cai = cai_all %>% filter(
  address %in% village_names,
  !is.na(network_id),
  id != 99,
  network_id != 99
) %>%
  dplyr::select(id, network_id) %>%
  as.matrix %>% 
  apply(2, as.character) %>%
  graph_from_edgelist(directed=FALSE) %>% 
  igraph::simplify()


orange = '#E69F00' #56B4E9
light_grey = '#C0C0C0'

p_cai = ggnet2(g_cai, node.size=2, node.color=orange, edge.color=light_grey) + #edge.alpha 
  ggtitle('Subsample of villages from Cai et al. (2015)') + 
  theme(plot.title = element_text(hjust = 0.5))
p_cai




load('data/caltech.Rdata')

p_fb = ggnet2(g_fb, node.size=2, node.color=orange, edge.color=light_grey) + 
  ggtitle('Caltech Facebook network, 2005') + 
  theme(plot.title = element_text(hjust = 0.5))
p_fb

p = grid.arrange(p_cai, p_fb, nrow=1)
ggsave('figures/schematic.png', p)
