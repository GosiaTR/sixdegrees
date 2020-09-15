# -*- coding: utf-8 -*-
"""
This module provides convenient wrapper functions for
the C++-functions sampling from the modified and original
small-world network models.
"""

import numpy as np
from _sixdegrees import (
        _modified_small_world_network_coord_lists,
        _modified_small_world_network,
        _original_small_world_network_coord_lists,
        _original_small_world_network,
    )

def original_small_world_network(
                                 N,
                                 k,
                                 mu,
                                 use_largest_component=False,
                                 delete_non_largest_component_nodes=True,
                                 seed=0,
                                 ):
    r"""
    Generate a small-world network according to the original algorithm
    by Watts and Strogatz.

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
        
    Returns
    =======
    N : int
        Number of nodes in the generated network (can be smaller than
        the demanded number of nodes if only the largest connected 
        component is returned)
    edges : list of tuple of int
        A list of node index pairs that constitute the edges of the
        generated network.

    Example
    =======

    >>> N, edges = original_small_world_network(N=512,k=7,mu=-0.5)
    >>> print(N)
    512
    """

    p_S = 1 / (1+((N-1-k)/k)**mu) # see Eq. (4.29) in B.F. Maier's dissertation, combined with Eq. (2a) in https://arxiv.org/abs/1901.02381
    p_ER = k / (N-1) 

    p_rewire = (1-p_S) / (1-p_ER) # see end of Sec. II.A in https://arxiv.org/abs/1901.02381

    _N, _edges = _original_small_world_network(
                                    N,
                                    k,
                                    p_rewire,
                                    use_largest_component=use_largest_component,
                                    delete_non_largest_component_nodes=delete_non_largest_component_nodes,
                                    seed=seed
                                    )

    return _N, _edges

def original_small_world_network_coord_lists(
                                 N,
                                 k,
                                 mu,
                                 use_largest_component=False,
                                 delete_non_largest_component_nodes=True,
                                 seed=0,
                                 ):
    r"""
    Generate a small-world network according to the original algorithm
    by Watts and Strogatz.

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

    Example
    =======

    >>> N, row, col = small_world_network_coord_lists(B=8,L=3,k=7,mu=-0.5)
    >>> print(N)
    512
    """


    p_S = 1 / (1+((N-1-k)/k)**mu) # see Eq. (4.29) in B.F. Maier's dissertation, combined with Eq. (2a) in https://arxiv.org/abs/1901.02381
    p_ER = k / (N-1) 

    p_rewire = (1-p_S) / (1-p_ER) # see end of Sec. II.A in https://arxiv.org/abs/1901.02381

    _N, row, col = _original_small_world_network_coord_lists(
                                    N,
                                    k,
                                    p_rewire,
                                    use_largest_component=use_largest_component,
                                    delete_non_largest_component_nodes=delete_non_largest_component_nodes,
                                    seed=seed
                                    )
    return _N, row, col

def small_world_network(
                                 N,
                                 k,
                                 mu,
                                 use_largest_component=False,
                                 delete_non_largest_component_nodes=True,
                                 seed=0,
                                 ):
    r"""
    Generate a small-world network according to the modified model
    definition where edges are drawn based on a distance-dependent
    probability.

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
        
    Returns
    =======
    N : int
        Number of nodes in the generated network (can be smaller than
        the demanded number of nodes if only the largest connected 
        component is returned)
    edges : list of tuple of int
        A list of node index pairs that constitute the edges of the
        generated network.

    Example
    =======

    >>> N, edges = small_world_network(N=512,k=7,mu=-0.5)
    >>> print(N)
    512
    """


    beta = k /(N-1-k) * ((N-1-k)/k)**mu

    _N, _edges = _modified_small_world_network(
                                    N,
                                    k,
                                    beta,
                                    use_largest_component=use_largest_component,
                                    delete_non_largest_component_nodes=delete_non_largest_component_nodes,
                                    seed=seed
                                    )

    return _N, _edges

def small_world_network_coord_lists(
                                 N,
                                 k,
                                 mu,
                                 use_largest_component=False,
                                 delete_non_largest_component_nodes=True,
                                 seed=0,
                                 ):
    r"""
    Generate a small-world network according to the modified model
    definition where edges are drawn based on a distance-dependent
    probability.

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

    Example
    =======

    >>> N, row, col = small_world_network_coord_lists(B=8,L=3,k=7,mu=-0.5)
    >>> print(N)
    512
    """

    beta = k /(N-1-k) * ((N-1-k)/k)**mu

    _N, row, col = _modified_small_world_network_coord_lists(
                                    N,
                                    k,
                                    beta,
                                    use_largest_component=use_largest_component,
                                    delete_non_largest_component_nodes=delete_non_largest_component_nodes,
                                    seed=seed
                                    )
    return _N, row, col
