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
