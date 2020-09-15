# -*- coding: utf-8 -*-
"""
This module provides convenient wrapper functions for
the C++-functions sampling from all random geometric
small-world network models.
"""

import numpy as np
from _sixdegrees import (
        _random_geometric_small_world_network_coord_lists,
        _random_geometric_small_world_network,
        _twoD_random_geometric_kleinberg_network_coord_lists,
        _twoD_random_geometric_kleinberg_network,
    )

from sixdegrees.kleinberg_helper_functions import (
        get_distance_connection_probability_parameters,
        get_continuous_distance_connection_probability_parameters,
        )

def twoD_random_geometric_kleinberg_network(
                                       N,
                                       k,
                                       mu,
                                       periodic_boundary_conditions=True,
                                       use_largest_component=False,
                                       delete_non_largest_component_nodes=True,
                                       seed=0,
                                       epsilon=1e-9,
                                       return_positions=False,
                                       use_continuous_connection_probability=True,
                                       X=None,
                                       Y=None,
                                       ):

    if X is None:
        X = np.random.rand(N)
    if Y is None:
        Y = np.random.rand(N)

    kappa = 2*(mu-1)
    #kappa = mu - 1

    if not use_continuous_connection_probability:
        C, rmin = get_distance_connection_probability_parameters(
                                                       N, 
                                                       k,
                                                       kappa, 
                                                       epsilon=epsilon, 
                                                       use_periodic_boundary_conditions=periodic_boundary_conditions,
                                                       )
    else:
        C, rmin = get_continuous_distance_connection_probability_parameters(
                                                       N, 
                                                       k,
                                                       kappa, 
                                                       epsilon=epsilon, 
                                                       use_periodic_boundary_conditions=periodic_boundary_conditions,
                                                       )

    _N, _edges = _twoD_random_geometric_kleinberg_network(
                                    N,
                                    C,
                                    rmin,
                                    kappa,
                                    X.tolist(),
                                    Y.tolist(),
                                    periodic_boundary_conditions=periodic_boundary_conditions,
                                    use_largest_component=use_largest_component,
                                    delete_non_largest_component_nodes=delete_non_largest_component_nodes,
                                    seed=seed
                                    )

    if return_positions:
        return _N, _edges, X, Y
    else:
        return _N, _edges


def random_geometric_small_world_network(N,
                                         k,
                                         mu,
                                         use_largest_component=False,
                                         delete_non_largest_component_nodes=True,
                                         seed=0,
                                         return_positions=False,
                                         ):

    beta = k /(N-1-k) * ((N-1-k)/k)**mu

    X = np.random.rand(N) * N
    # positions do not need to be sorted necessarily
    # but it helps knowing that 0, 1, 2,... are 
    # the first, second, third ... nodes in 
    # the interval [0,N]
    
    X.sort()
    _N, _edges = _random_geometric_small_world_network(
                                    N,
                                    k,
                                    beta,
                                    X.tolist(),
                                    use_largest_component=use_largest_component,
                                    delete_non_largest_component_nodes=delete_non_largest_component_nodes,
                                    seed=seed
                                    )

    if return_positions:
        return _N, _edges, X
    else:
        return _N, _edges


def random_geometric_small_world_network_coord_lists(
                                         N,
                                         k,
                                         mu,
                                         use_largest_component=False,
                                         delete_non_largest_component_nodes=True,
                                         seed=0,
                                         return_positions=False,
                                         ):

    beta = k /(N-1-k) * ((N-1-k)/k)**mu
    X = np.random.rand(N) * N

    # positions do not need to be sorted necessarily
    # but it helps knowing that 0, 1, 2,... are 
    # the first, second, third ... nodes in 
    # the interval [0,N]
    
    X.sort() 

    _N, row, col = _random_geometric_small_world_network_coord_lists(
                                    N,
                                    k,
                                    beta,
                                    X.tolist(),
                                    use_largest_component=use_largest_component,
                                    delete_non_largest_component_nodes=delete_non_largest_component_nodes,
                                    seed=seed
                                    )

    if return_positions:
        return _N, row, col, X
    else:
        return _N, row, col

def twoD_random_geometric_kleinberg_network_coord_lists(
                                       N,
                                       k,
                                       mu,
                                       periodic_boundary_conditions=True,
                                       use_largest_component=False,
                                       delete_non_largest_component_nodes=True,
                                       seed=0,
                                       epsilon=1e-9,
                                       return_positions=False,
                                       use_continuous_connection_probability=True,
                                       X = None,
                                       Y = None,
                                       ):

    if X is None:
        X = np.random.rand(N)
    if Y is None:
        Y = np.random.rand(N)

    kappa = 2*(mu-1)
    #kappa = mu - 1

    if not use_continuous_connection_probability:
        C, rmin = get_distance_connection_probability_parameters(
                                                       N, 
                                                       k,
                                                       kappa, 
                                                       epsilon=epsilon, 
                                                       use_periodic_boundary_conditions=periodic_boundary_conditions,
                                                       )
    else:
        C, rmin = get_continuous_distance_connection_probability_parameters(
                                                       N, 
                                                       k,
                                                       kappa, 
                                                       epsilon=epsilon, 
                                                       use_periodic_boundary_conditions=periodic_boundary_conditions,
                                                       )

    _N, row, col = _twoD_random_geometric_kleinberg_network_coord_lists(
                                    N,
                                    C,
                                    rmin,
                                    kappa,
                                    X.tolist(),
                                    Y.tolist(),
                                    periodic_boundary_conditions=periodic_boundary_conditions,
                                    use_largest_component=use_largest_component,
                                    delete_non_largest_component_nodes=delete_non_largest_component_nodes,
                                    seed=seed
                                    )

    if return_positions:
        return _N, row, col, X, Y
    else:
        return _N, row, col


if __name__=="__main__":
    import matplotlib.pyplot as pl
    N = 401
    k = 30
    mu = -1
    kappa = 2*(mu-1)
    _, edges = twoD_random_geometric_kleinberg_network(N,k,mu,periodic_boundary_conditions=False,use_continuous_connection_probability=True)

    print(2*len(edges)/N, k)

    C, rmin = get_continuous_distance_connection_probability_parameters(
                                                   N, 
                                                   k,
                                                   2*(mu-1), 
                                                   epsilon=1e-9, 
                                                   use_periodic_boundary_conditions=False
                                                   )

    def P(r, C, rmin, kappa):
        if r < rmin:
            return 1
        else:
            return C * r**kappa
    r = np.linspace(0,1,10000)
    Ps = [P(_r,C,rmin, kappa) for _r in r]
    pl.plot(r, Ps)
    #pl.xscale('log')
    #pl.yscale('log')

    pl.show()
    #print(edges)
