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

