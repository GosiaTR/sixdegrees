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
        _random_geometric_kleinberg_network_coord_lists,
        _random_geometric_kleinberg_network,
        _twoD_random_geometric_kleinberg_network_coord_lists,
        _twoD_random_geometric_kleinberg_network,
    )

from sixdegrees.kleinberg_helper_functions import (
        get_distance_connection_probability_parameters,
        get_continuous_distance_connection_probability_parameters,
        )

def random_geometric_kleinberg_network(N,
                                       k,
                                       mu,
                                       use_largest_component=False,
                                       delete_non_largest_component_nodes=True,
                                       seed=0,
                                       return_positions=False,
                                       ):
    r"""
    Generate a Kleinberg small-world network from 
    nodes that are uniformly i.i.d. on a ring.

    Parameters
    ==========
    N : int
        Number of nodes.
    k : float
        Mean degree.
    mu : float
        Structural control parameter.
    use_largest_component : bool, default = False
        Whether or not to only return the largest connected component
        of size :math:`N_g\leq N`.
    delete_non_largest_component_nodes : bool, default = True
        If only the largest connected component is returned, relabel
        the nodes to be indexed from zero to :math:`N_g`.
    seed : int, default = 0
        The seed with which to generate this network. Choose ``0``
        for random initialization.
    return_positions : bool, default = False
        Whether or not to return the node positions on the ring.
        
    Returns
    =======
    N : int
        Number of nodes in the generated network (can be smaller than
        the demanded number of nodes if only the largest connected 
        component is returned)
    edges : list of tuple of int
        A list of node index pairs that constitute the edges of the
        generated network.
    X : numpy.ndarray
        Node positions on the ring. Is only returned if 
        ``return_positions = True`` was set.

    Example
    =======

    >>> N, edges = random_geometric_kleinberg_network(N=512,k=7,mu=-0.5)
    >>> print(N)
    512
    """

    X = np.random.rand(N) * N
    # positions do not need to be sorted necessarily
    # but it helps knowing that 0, 1, 2,... are 
    # the first, second, third ... nodes in 
    # the interval [0,N]
    
    X.sort()
    _N, _edges = _random_geometric_kleinberg_network(
                                    N,
                                    k,
                                    mu,
                                    X.tolist(),
                                    use_largest_component=use_largest_component,
                                    delete_non_largest_component_nodes=delete_non_largest_component_nodes,
                                    seed=seed
                                    )

    if return_positions:
        return _N, _edges, X
    else:
        return _N, _edges

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
                                       use_continuous_connection_probability=False,
                                       X=None,
                                       Y=None,
                                       ):
    r"""
    Generate a Kleinberg small-world network from 
    nodes that are uniformly i.i.d. on a square.

    Parameters
    ==========
    N : int
        Number of nodes.
    k : float
        Mean degree.
    mu : float
        Structural control parameter.
    periodic_boundary_conditions : bool, default = True
        Whether or not to apply periodic boundary conditions (2-torus).
    use_largest_component : bool, default = False
        Whether or not to only return the largest connected component
        of size :math:`N_g\leq N`.
    delete_non_largest_component_nodes : bool, default = True
        If only the largest connected component is returned, relabel
        the nodes to be indexed from zero to :math:`N_g`.
    seed : int, default = 0
        The seed with which to generate this network. Choose ``0``
        for random initialization.
    epsilon : float, default = 1e-9
        Minimal distance below which node pairs are always connected.
    return_positions : bool, default = False
        Whether or not to return the node positions on the ring.
    use_continuous_connection_probability : bool, default = False
        Excess connection probability can be redistributed to the whole
        ensemble or to nearest neighbors. If this flag is true, the
        excess probability will be redistributed to the whole ensemble.
    X : numpy.ndarray, default = None
        Node positions in x-direction (:math:`x \in [0,1)`).
    Y : numpy.ndarray, default = None
        Node positions in y-direction (:math:`x \in [0,1)`).
        
    Returns
    =======
    N : int
        Number of nodes in the generated network (can be smaller than
        the demanded number of nodes if only the largest connected 
        component is returned)
    edges : list of tuple of int
        A list of node index pairs that constitute the edges of the
        generated network.
    X : numpy.ndarray
        Node positions on the ring. Is only returned if 
        ``return_positions = True`` was set.
    Y : numpy.ndarray
        Node positions on the ring. Is only returned if 
        ``return_positions = True`` was set.
    """
    if X is None:
        X = np.random.rand(N)
    if Y is None:
        Y = np.random.rand(N)

    kappa = mu-2

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
    r"""
    Generate a small-world network from 
    nodes that are uniformly i.i.d. on a ring.

    Parameters
    ==========
    N : int
        Number of nodes.
    k : float
        Mean degree.
    mu : float
        Structural control parameter.
    use_largest_component : bool, default = False
        Whether or not to only return the largest connected component
        of size :math:`N_g\leq N`.
    delete_non_largest_component_nodes : bool, default = True
        If only the largest connected component is returned, relabel
        the nodes to be indexed from zero to :math:`N_g`.
    seed : int, default = 0
        The seed with which to generate this network. Choose ``0``
        for random initialization.
    return_positions : bool, default = False
        Whether or not to return the node positions on the ring.
        
    Returns
    =======
    N : int
        Number of nodes in the generated network (can be smaller than
        the demanded number of nodes if only the largest connected 
        component is returned)
    edges : list of tuple of int
        A list of node index pairs that constitute the edges of the
        generated network.
    X : numpy.ndarray
        Node positions on the ring. Is only returned if 
        ``return_positions = True`` was set.
    """

    beta = k /(N-1-k) * ((N-1-k)/k)**mu

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


def random_geometric_kleinberg_network_coord_lists(
                                         N,
                                         k,
                                         mu,
                                         use_largest_component=False,
                                         delete_non_largest_component_nodes=True,
                                         seed=0,
                                         return_positions=False,
                                         ):
    r"""
    Generate a Kleinberg small-world network from 
    nodes that are uniformly i.i.d. on a ring.

    Parameters
    ==========
    N : int
        Number of nodes.
    k : float
        Mean degree.
    mu : float
        Structural control parameter.
    use_largest_component : bool, default = False
        Whether or not to only return the largest connected component
        of size :math:`N_g\leq N`.
    delete_non_largest_component_nodes : bool, default = True
        If only the largest connected component is returned, relabel
        the nodes to be indexed from zero to :math:`N_g`.
    seed : int, default = 0
        The seed with which to generate this network. Choose ``0``
        for random initialization.
    return_positions : bool, default = False
        Whether or not to return the node positions on the ring.
        
    Returns
    =======
    N : int
        Number of nodes in the generated network (can be smaller than
        the demanded number of nodes if only the largest connected 
        component is returned)
    row : list of int
        Row coordinates of non-zero entries in the generated adjacency matrix
    col : list of int
        Column coordinates of non-zero entries in the generated adjacency matrix
    X : numpy.ndarray
        Node positions on the ring. Is only returned if 
        ``return_positions = True`` was set.
    """

    X = np.random.rand(N) * N
    # positions do not need to be sorted necessarily
    # but it helps knowing that 0, 1, 2,... are 
    # the first, second, third ... nodes in 
    # the interval [0,N]
    
    X.sort() 

    _N, row, col = _random_geometric_kleinberg_network_coord_lists(
                                    N,
                                    k,
                                    mu,
                                    X.tolist(),
                                    use_largest_component=use_largest_component,
                                    delete_non_largest_component_nodes=delete_non_largest_component_nodes,
                                    seed=seed
                                    )

    if return_positions:
        return _N, row, col, X
    else:
        return _N, row, col

def random_geometric_small_world_network_coord_lists(
                                         N,
                                         k,
                                         mu,
                                         use_largest_component=False,
                                         delete_non_largest_component_nodes=True,
                                         seed=0,
                                         return_positions=False,
                                         ):
    r"""
    Generate a small-world network from 
    nodes that are uniformly i.i.d. on a ring.

    Parameters
    ==========
    N : int
        Number of nodes.
    k : float
        Mean degree.
    mu : float
        Structural control parameter.
    use_largest_component : bool, default = False
        Whether or not to only return the largest connected component
        of size :math:`N_g\leq N`.
    delete_non_largest_component_nodes : bool, default = True
        If only the largest connected component is returned, relabel
        the nodes to be indexed from zero to :math:`N_g`.
    seed : int, default = 0
        The seed with which to generate this network. Choose ``0``
        for random initialization.
    return_positions : bool, default = False
        Whether or not to return the node positions on the ring.
        
    Returns
    =======
    N : int
        Number of nodes in the generated network (can be smaller than
        the demanded number of nodes if only the largest connected 
        component is returned)
    row : list of int
        Row coordinates of non-zero entries in the generated adjacency matrix
    col : list of int
        Column coordinates of non-zero entries in the generated adjacency matrix
    X : numpy.ndarray
        Node positions on the ring. Is only returned if 
        ``return_positions = True`` was set.
    """

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
                                       use_continuous_connection_probability=False,
                                       X = None,
                                       Y = None,
                                       ):
    r"""
    Generate a Kleinberg (power-law) small-world network from 
    nodes that are uniformly i.i.d. on a square.

    Parameters
    ==========
    N : int
        Number of nodes.
    k : float
        Mean degree.
    mu : float
        Structural control parameter.
    periodic_boundary_conditions : bool, default = True
        Whether or not to apply periodic boundary conditions (2-torus).
    use_largest_component : bool, default = False
        Whether or not to only return the largest connected component
        of size :math:`N_g\leq N`.
    delete_non_largest_component_nodes : bool, default = True
        If only the largest connected component is returned, relabel
        the nodes to be indexed from zero to :math:`N_g`.
    seed : int, default = 0
        The seed with which to generate this network. Choose ``0``
        for random initialization.
    epsilon : float, default = 1e-9
        Minimal distance below which node pairs are always connected.
    return_positions : bool, default = False
        Whether or not to return the node positions on the ring.
    use_continuous_connection_probability : bool, default = False
        Excess connection probability can be redistributed to the whole
        ensemble or to nearest neighbors. If this flag is true, the
        excess probability will be redistributed to the whole ensemble.
    X : numpy.ndarray, default = None
        Node positions in x-direction (:math:`x \in [0,1)`).
    Y : numpy.ndarray, default = None
        Node positions in y-direction (:math:`x \in [0,1)`).
        
        
    Returns
    =======
    N : int
        Number of nodes in the generated network (can be smaller than
        the demanded number of nodes if only the largest connected 
        component is returned)
    row : list of int
        Row coordinates of non-zero entries in the generated adjacency matrix
    col : list of int
        Column coordinates of non-zero entries in the generated adjacency matrix
    X : numpy.ndarray
        Node positions on the ring. Is only returned if 
        ``return_positions = True`` was set.
    Y : numpy.ndarray
        Node positions on the ring. Is only returned if 
        ``return_positions = True`` was set.
    """

    if X is None:
        X = np.random.rand(N)
    if Y is None:
        Y = np.random.rand(N)

    kappa = mu-2

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
