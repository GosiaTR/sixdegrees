# -*- coding: utf-8 -*-
"""
This module provides functions to sample 
from the generalized small-world models
where only node-pair distances have to 
be provided.
"""

import numpy as np
from _sixdegrees import find_first_random_edge

def _is_sorted(arr):
    """ Returns `True` if the array is sorted ascendingly and `False if it isn't."""
    i = 0     
    while (i+1 < arr.shape[0]):
        if arr[i] > arr[i+1]:
            return False
        i += 1
    return True

def generalized_kleinberg_network(
                            N,
                            k,
                            all_node_pairs,
                            node_pair_distances = None,
                            node_pair_probabilities = None,
                            kappa = None,
                            seed = None,
                            ):

    mmax = (N*(N-1))//2
    assert(all_node_pairs.shape == (mmax, 2))

    if node_pair_probabilities:
        assert(len(node_pair_probabilities) == len(all_node_pairs))

    if (kappa is None and node_pair_probabilities is None) or (kappa is not None and node_pair_probabilities is not None):
        raise ValueError("either provide node pair probabilities or kappa")

    if kappa is not None:
        assert(kappa <= 0)
        assert(node_pair_distances is not None)
        assert(len(node_pair_distances) == len(all_node_pairs))

        # if kappa and distances are set, compute probabilities as a power-law
        node_pair_probabilities = node_pair_distances**kappa
    
    # seed the engine
    np.random.seed(seed)

    # compute the prefactor of the connection probability such that
    # eventually the total sum of produced
    # edges equals the expected number of edges
    number_of_edges = N*k/2.0
    prefactor = number_of_edges / node_pair_probabilities.sum()
    node_pair_probabilities *= prefactor

    # find the indices that sort the connection probabilities descendingly
    ndx = (np.argsort(node_pair_probabilities))[::-1]

    # sort node pairs and their corresponding existence probabilities
    all_node_pairs = all_node_pairs[ndx,:]
    node_pair_probabilities = node_pair_probabilities[ndx,]

    # compute the excess edges per node pair as the cumulative sum of p-1.0 and find the first
    # node pair where this cumulated number of excess edges is lower than 0
    first_random_edge = find_first_random_edge(node_pair_probabilities)

    # set all node pairs below this "radius" to be connected
    short_range = all_node_pairs[:first_random_edge,:]

    # find indices of the node pairs where edges are set
    # according to a probability p < 1.0
    long_range_ndx = np.where(np.random.rand(mmax-first_random_edge) < \
                              node_pair_probabilities[first_random_edge:]
                             )[0]
    long_range_ndx += first_random_edge

    # get those node pairs
    long_range = all_node_pairs[long_range_ndx,:]

    # concatenate short-ranged and long-ranged edges
    edges = np.vstack((short_range, long_range))

    return edges

def generalized_small_world_network(
                            N,
                            k,
                            beta,
                            all_node_pairs,
                            node_pair_distances,
                            seed = None,
                            ):

    assert(N>1)
    mmax = (N*(N-1))//2
    assert(all_node_pairs.shape == (mmax, 2))
    assert(node_pair_distances is not None)
    assert(len(node_pair_distances) == len(all_node_pairs))
    assert(beta >= 0.0)
    assert(beta <= 1.0)

    # seed the engine
    np.random.seed(seed)

    pS = 1.0 / (1.0 + beta * (N-1-k)/k)
    pL = beta * pS

    # edges equals the expected number of edges
    number_of_edges = short_range_radius = int(N*k/2.0)

    # find the indices that sort the connection probabilities descendingly
    ndx = np.argsort(node_pair_distances)

    # sort node pairs and their corresponding existence probabilities
    all_node_pairs = all_node_pairs[ndx,:]
    node_pair_distances = node_pair_distances[ndx,]

    # set all node pairs below this "radius" to be connected with probability pS
    short_range_ndx = np.where(np.random.rand(short_range_radius) < pS)[0]
    short_range = all_node_pairs[short_range_ndx,:]

    # find indices of the node pairs where edges are set
    # according to the long-range probability pL
    long_range_ndx = np.where(np.random.rand(mmax-short_range_radius) < pL)[0]
    long_range_ndx += short_range_radius

    # get those node pairs
    long_range = all_node_pairs[long_range_ndx,:]

    # concatenate short-ranged and long-ranged edges
    edges = np.vstack((short_range, long_range))

    return edges

def generalized_categorical_kleinberg_network(
                            N,
                            k,
                            all_node_pairs,
                            node_pair_categories,
                            category_probabilities,
                            seed = None,
                            ):


    assert(N>1)
    mmax = (N*(N-1))//2
    print(all_node_pairs.shape, (mmax,2))
    assert(all_node_pairs.shape == (mmax, 2))
    assert(len(node_pair_categories) == len(all_node_pairs))

    # check if probabilities are sorted
    assert(_is_sorted(category_probabilities[::-1]))

    # get the node pair counts of the categories
    unique, m = np.unique(node_pair_categories, return_counts=True)
    m = m[np.argsort(unique)]

    # check that all categories are actually present in the node_pair_categories,
    # otherwise there will be problems in computing the excess probability
    assert(np.all(np.sort(unique) == np.arange(len(category_probabilities),dtype=int)))

    # properly norm the probabilities such that one obtains a mean degree of k
    prefactor = N*k/2.0 / category_probabilities.dot(m)
    category_probabilities *= prefactor

    # redistribute probability if the smaller categories exceed probability of p = 1
    if category_probabilities[0] > 1.0:
        l = 0
        excess_edges = (category_probabilities[0]-1.0) * m[0]
        while excess_edges > 0:
            l += 1
            if (l == len(m)):
                raise ValueError("There's not enough node pairs to redistribute all excess probability.")
            excess_edges += (category_probabilities[l]-1.0) * m[l]

        excess_edges -= (category_probabilities[l]-1.0)*m[l];

        # copy probabilities
        p = np.array(category_probabilities,dtype=float)

        # set short-range probabilities to 1
        p[:l] = 1.0

        # add the remaining excess probability to the critical category
        # all remaining categories will be connected as before
        p[l] += excess_edges / m[l]
    else:
        p = category_probabilities

    np.random.seed(seed)

    # find node pair indices according to their category connection probability
    create_ndx = np.where(np.random.rand(mmax) < p[node_pair_categories])[0]

    return all_node_pairs[create_ndx,:]



if __name__ == "__main__":
    N = 100
    k = 10
    kappas = [-2,-1.5,-1,-0.5,0]
    X = np.random.rand(N)

    all_node_pairs = np.array([ (i,j) for i in range(N-1) for j in range(i+1,N) ],dtype=int)
    all_distances = np.array([ np.min([abs(X[i]-X[j]),abs(float(1.0)-abs((X[j]-X[i])))]) for i in range(N-1) for j in range(i+1,N)])

    for kappa in kappas:
        edges = generalized_kleinberg_network(N, k, all_node_pairs, node_pair_distances = all_distances, kappa=kappa)
        print(edges.shape)
        print("kappa =", kappa, "; expected k =", k, "; measured k =", 2.0*len(edges)/N)
        beta = (k/(N-1-k)) * ((N-1-k)/k)**(kappa+1)
        edges = generalized_small_world_network(N, k, beta, all_node_pairs, all_distances)
        print(edges.shape)
        print("kappa =", kappa, "; expected k =", k, "; measured k =", 2.0*len(edges)/N)

        # some two-block stochastick block model
        _N = 30
        _k = 20
        B = _N // 2
        _all_node_pairs = np.array([ (i,j) for i in range(_N-1) for j in range(i+1,_N) ],dtype=int)
        categorical_distances =  np.array([ 0 if i // B == j // B else 1 for i in range(_N-1) for j in range(i+1,_N) ],dtype=int)
        xi = B**(kappa+1.0)
        probabilities = (xi/B)**np.array([0.0,1.0])

        edges = generalized_categorical_kleinberg_network(_N, _k, _all_node_pairs, categorical_distances, probabilities)
        print("kappa =", kappa, "; expected k =", _k, "; measured k =", 2.0*len(edges)/_N)

        #for i in range(_N-1):
        #    for j in range(i+1,_N):
        #        print(i,j, i % B, j % B, i % B == j % B)
