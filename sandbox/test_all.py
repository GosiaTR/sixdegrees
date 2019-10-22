import numpy as np
from matplotlib import pyplot as pl

from sixdegrees import (    
        small_world_network_coord_lists,
        original_small_world_network_coord_lists,
        modular_hierarchical_network_coord_lists,
        random_geometric_kleinberg_network_coord_lists,
        random_geometric_small_world_network_coord_lists,
        twoD_random_geometric_kleinberg_network_coord_lists,
        kleinberg_network_coord_lists,
        )

from sixdegrees.utils import to_sparse_matrix


def get_k_and_k2(A):
    k = np.array(A.sum(axis=0),dtype=float).flatten()
    #print(k)
    return np.mean(k), np.mean(k**2)

if __name__=="__main__":

    NMEAS = 20
    B = 10
    L = 2
    N = B**L
    k = 16
    mus = np.linspace(-1,1,50)

    generators = (
            modular_hierarchical_network_coord_lists,
            small_world_network_coord_lists,
            original_small_world_network_coord_lists,
            random_geometric_kleinberg_network_coord_lists,
            random_geometric_small_world_network_coord_lists,
            twoD_random_geometric_kleinberg_network_coord_lists,
            kleinberg_network_coord_lists,
        )

    all_k = np.zeros((len(generators), len(mus), NMEAS, ))
    all_k2 = np.zeros((len(generators), len(mus), NMEAS, ))

    fig, ax = pl.subplots(1,2,figsize=(8,4),sharex=True)

    for ig, g in enumerate(generators):

        print("now testing:", g.__name__)

        for imu, mu in enumerate(mus):

            for meas in range(NMEAS):

                if ig == 0:
                    A = to_sparse_matrix(*g(B,L,k,mu))
                else:
                    A = to_sparse_matrix(*g(N,k,mu))

                all_k[ig,imu,meas], all_k2[ig,imu,meas] = get_k_and_k2(A)

        _k = all_k[ig].mean(axis=1)
        _k2 = all_k2[ig].mean(axis=1)
        ax[0].plot(mus, _k,label=g.__name__)
        ax[1].plot(mus, _k2)

    ax[0].set_xlabel(r'structural control $\mu$')
    ax[0].set_ylabel(r'mean degree')
    ax[1].set_ylabel(r'second moment degree')

    ax[0].legend()

    fig.tight_layout()

    pl.show()

