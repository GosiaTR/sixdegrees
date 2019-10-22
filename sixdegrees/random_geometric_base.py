import numpy as np
from _sixdegrees import (
        _random_geometric_small_world_network_coord_lists,
        _random_geometric_small_world_network,
        _random_geometric_kleinberg_network_coord_lists,
        _random_geometric_kleinberg_network,
    )

def random_geometric_kleinberg_network(N,
                                       k,
                                       mu,
                                       use_giant_component=False,
                                       delete_non_giant_component_nodes=True,
                                       seed=0,
                                       return_positions=False,
                                       ):

    X = np.random.rand() * N
    X.sort()
    _N, _edges = _random_geometric_kleinberg_network(
                                    N,
                                    k,
                                    mu,
                                    X.tolist(),
                                    use_giant_component=use_giant_component,
                                    delete_non_giant_component_nodes=delete_non_giant_component_nodes,
                                    seed=seed
                                    )

    if return_posititons:
        return _N, _edges, X
    else:
        return _N, _edges


def random_geometric_small_world_network(N,
                                         k,
                                         beta,
                                         use_giant_component=False,
                                         delete_non_giant_component_nodes=True,
                                         seed=0,
                                         return_positions=False,
                                         ):

    X = np.random.rand() * N
    X.sort()
    _N, _edges = _random_geometric_small_world_network(
                                    N,
                                    k,
                                    beta,
                                    X.tolist(),
                                    use_giant_component=use_giant_component,
                                    delete_non_giant_component_nodes=delete_non_giant_component_nodes,
                                    seed=seed
                                    )

    if return_posititons:
        return _N, _edges, X
    else:
        return _N, _edges


def random_geometric_kleinberg_network_coord_lists(
                                         N,
                                         k,
                                         mu,
                                         use_giant_component=False,
                                         delete_non_giant_component_nodes=True,
                                         seed=0,
                                         return_positions=False,
                                         ):

    X = np.random.rand() * N
    X.sort()
    _N, row, col = _random_geometric_small_world_network_coord_lists(
                                    N,
                                    k,
                                    mu,
                                    X,
                                    use_giant_component=use_giant_component,
                                    delete_non_giant_component_nodes=delete_non_giant_component_nodes,
                                    seed=seed
                                    )

    if return_posititons:
        return _N, row, col, X
    else:
        return _N, row, col

def random_geometric_small_world_network_coord_lists(
                                         N,
                                         k,
                                         beta,
                                         use_giant_component=False,
                                         delete_non_giant_component_nodes=True,
                                         seed=0,
                                         return_positions=False,
                                         ):

    X = np.random.rand() * N
    X.sort()
    _N, row, col = _random_geometric_small_world_network_coord_lists(
                                    N,
                                    k,
                                    beta,
                                    X,
                                    use_giant_component=use_giant_component,
                                    delete_non_giant_component_nodes=delete_non_giant_component_nodes,
                                    seed=seed
                                    )

    if return_posititons:
        return _N, row, col, X
    else:
        return _N, row, col

