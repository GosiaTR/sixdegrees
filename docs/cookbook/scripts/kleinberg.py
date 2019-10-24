import sixdegrees as sixd
import netwulf as wulf
from pathlib import Path
import matplotlib.pyplot as plt

B = 10
L = 3
N = B**L
k = 12
mus = [-1,-.5,-.01,0.01]

for imu, mu in enumerate(mus):
    path = Path('data') / ('kleinberg_2d_imu_{}.json'.format(imu))
    if not path.exists():
        G = sixd.to_networkx_graph(*sixd.twoD_random_geometric_kleinberg_network(N,k,mu))
        nw, cfg = wulf.visualize(G,config={'node_fill_color':'#efefef',
                                           'node_collision': False,
                                           'wiggle_nodes': False,
                                           'scale_node_size_by_strength': True,
                                           'node_size': 8,
                                           'node_stroke_width': 2,
                                           'node_stroke_color': '#000000',
                                           })
        wulf.save(str(path), nw, cfg, G)
    else:
        nw, cfg, G = wulf.load(str(path))

    plt.figure()
    fig, ax = wulf.draw_netwulf(nw)

plt.show()
