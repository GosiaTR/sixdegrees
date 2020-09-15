# -*- coding: utf-8 -*-
"""
This module provides convenient wrapper functions for
the C++-functions sampling from the modular-hierarchical
small-world network model.
"""

import numpy as np
from _sixdegrees import (
        _modular_hierarchical_network_coord_lists,
        _modular_hierarchical_network,
    )

def modular_hierarchical_network(B,
                                 L,
                                 k,
                                 mu,
                                 use_largest_component=False,
                                 delete_non_largest_component_nodes=True,
                                 seed=0,
                                 ):
    r"""
    Generate a self-similar modular hierarchical random network
    with :math:`N=B^L` nodes.

    Parameters
    ==========
    B : int
        Number of submodules per module.
    L : int
        Number of hierarchical layers (height of hierarchical tree).
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

    >>> N, edges = modular_hierarchical_network(B=8,L=3,k=7,mu=-0.5)
    >>> print(N)
    512
    """

    xi = B**mu

    _N, _edges = _modular_hierarchical_network(
                                    B,
                                    L,
                                    k,
                                    xi,
                                    use_largest_component=use_largest_component,
                                    delete_non_largest_component_nodes=delete_non_largest_component_nodes,
                                    seed=seed
                                    )

    return _N, _edges


def modular_hierarchical_network_coord_lists(
                                 B,
                                 L,
                                 k,
                                 mu,
                                 use_largest_component=False,
                                 delete_non_largest_component_nodes=True,
                                 seed=0,
                                 ):
    r"""
    Generate a self-similar modular hierarchical random network
    with :math:`N=B^L` nodes.

    Parameters
    ==========
    B : int
        Number of submodules per module.
    L : int
        Number of hierarchical layers (height of hierarchical tree).
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

    >>> N, row, col = modular_hierarchical_network_coord_lists(B=8,L=3,k=7,mu=-0.5)
    >>> print(N)
    512
    """

    xi = B**mu
    _N, row, col = _modular_hierarchical_network_coord_lists(
                                    B,
                                    L,
                                    k,
                                    xi,
                                    use_largest_component=use_largest_component,
                                    delete_non_largest_component_nodes=delete_non_largest_component_nodes,
                                    seed=seed
                                    )

    return _N, row, col

