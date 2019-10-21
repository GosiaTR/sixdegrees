import numpy as np
import sixdegrees
import networkx as nx
import time

N = 1000
k = 3
p = k/(N-1.)

G = nx.fast_gnp_random_graph(N,p)


tic = time.time()
comp_sizes_sixdegrees = sorted([len(g) for g in sixdegrees.get_components(N,list(G.edges()))])
toc = time.time()

ticnx = time.time()
comp_sizes_nx = sorted([len(g) for g in list(nx.connected_component_subgraphs(G))])
tocnx = time.time()

print(toc-tic,comp_sizes_sixdegrees)
print(tocnx-ticnx,comp_sizes_nx)

