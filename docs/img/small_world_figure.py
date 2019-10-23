from smallworld.draw import draw_network
from smallworld import get_smallworld_graph

import matplotlib.pyplot as pl

# define network parameters
N = 11
k_over_2 = 2
betas = [0.025]

focal_node = 0

fig, ax = pl.subplots(1,1,figsize=(3,3))
ax = [ax]


# scan beta values
for ib, beta in enumerate(betas):

    # generate small-world graphs and draw
    G = get_smallworld_graph(N, k_over_2, beta)
    draw_network(G,k_over_2,focal_node=focal_node,ax=ax[ib])

# show
pl.subplots_adjust(wspace=0.3)
pl.savefig('SW.pdf')
pl.show()
