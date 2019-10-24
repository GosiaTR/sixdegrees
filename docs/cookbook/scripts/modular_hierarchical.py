import sixdegrees as sixd
import netwulf as wulf
from pathlib import Path
import matplotlib.pyplot as plt

B = 10
L = 3
N = B**L
k = 12
mus = [-1,-.5,-.15,1]

for imu, mu in enumerate(mus):
    path = Path('data') / ('modular_hierarchical_imu_{}.json'.format(imu))
    if not path.exists():
        G = sixd.to_networkx_graph(*sixd.modular_hierarchical_network(B,L,k,mu))
        nw, cfg = wulf.visualize(G,config={'node_color':'#dddddd'})
        wulf.save(str(path), nw, cfg, G)
    else:
        nw, cfg, G = wulf.load(str(path))

    plt.figure()
    fig, ax = wulf.draw_netwulf(nw)

plt.show()
