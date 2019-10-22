import numpy as np
import networkx as nx
import scipy.sparse as sprs

def to_networkx_graph(N, edges):

    G = nx.Graph()
    G.add_nodes_from(range(N))
    G.add_edges_from(edges)

    return G

def to_sparse_matrix(N, row, col, X=None, Y=None, generator='csr'):

    if generator == 'csr':
        generator = sprs.csr_matrix

    A = generator((np.ones_like(row),(row,col)), shape=(N,N))

    return A


