# -*- coding: utf-8 -*-
"""
This module provides convenient functions to
convert network data to more convenient
network objects.
"""

import numpy as np
import networkx as nx
import scipy.sparse as sprs

def to_networkx_graph(N, edges):
    """
    Convert a generated network to a networkx Graph object.

    Parameters
    ==========
    N : int
        Number of nodes.
    edges : list of tuple of int
        A list of node index pairs that constitute the edges of the
        generated network.

    Returns
    =======
    G : networkx.Graph
        Same network as a networkx Graph object.
    """

    G = nx.Graph()
    G.add_nodes_from(range(N))
    G.add_edges_from(edges)

    return G

def to_sparse_matrix(N, row, col, generator='csr'):
    """
    Convert a generated network to a scipy sparse matrix.

    Parameters
    ==========
    N : int
        Number of nodes.
    row : list of int
        Row coordinates of non-zero entries in the adjacency matrix
    col : list of int
        Column coordinates of non-zero entries in the adjacency matrix
    generator : str or scipy.sparse, default = 'csr'
        What to use to generate the adjacency matrix. Could be 'csr',
        'csc', or any sparse matrix format from ``scipy.sparse``.

    Returns
    =======
    A : sparse matrix
        Adjacency matrix in ``scipy.sparse`` format.
    """

    if generator == 'csr':
        generator = sprs.csr_matrix
    elif generatro == 'csc':
        generator = sprs.csc_matrix

    A = generator((np.ones_like(row),(row,col)), shape=(N,N))

    return A


