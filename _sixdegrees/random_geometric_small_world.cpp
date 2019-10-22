/* 
 * The MIT License (MIT)
 * Copyright (c) 2016, Benjamin Maier
 *
 * Permission is hereby granted, free of charge, to any person 
 * obtaining a copy of this software and associated documentation 
 * files (the "Software"), to deal in the Software without 
 * restriction, including without limitation the rights to use, 
 * copy, modify, merge, publish, distribute, sublicense, and/or 
 * sell copies of the Software, and to permit persons to whom the 
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall 
 * be included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NON-
 * INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS 
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN 
 * AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF 
 * OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
 * IN THE SOFTWARE.
 */

#include "Utilities.h"
#include "random_geometric_small_world.h"

#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <set>
#include <utility>
#include <random>
#include <cmath>
#include <numeric>
#include <random>
#include <ctime>
#include <tuple>
#include <assert.h>

using namespace std;

tuple < size_t, vector <size_t>, vector<size_t> > random_geometric_small_world_coord_lists(
        size_t N,
        double k,
        double p,
        vector < double > r,
        bool use_largest_component,
        bool delete_non_largest_component_nodes,
        size_t seed
        )
{
    vector < set < size_t > * > G;

    G = random_geometric_small_world_neighbor_set(N,k,p,r,use_largest_component,seed);

    vector < size_t > rows;
    vector < size_t > cols;
    size_t new_N = neighbor_set_to_coord_lists(G,rows,cols,use_largest_component,delete_non_largest_component_nodes);
    
    return make_tuple(new_N,rows,cols);
}
    
pair < size_t, vector < pair < size_t, size_t > > > random_geometric_small_world_edge_list(
        size_t N,
        double k,
        double p,
        vector < double > r,
        bool use_largest_component,
        bool delete_non_largest_component_nodes,
        size_t seed
        )
{
    vector < set < size_t > * > G;

    G = random_geometric_small_world_neighbor_set(N,k,p,r,use_largest_component,seed);

    vector < pair < size_t, size_t > > edge_list;
    size_t new_N = neighbor_set_to_edge_list(G,edge_list,use_largest_component,delete_non_largest_component_nodes);

    return make_pair(new_N,edge_list);
}

vector < set < size_t > * > random_geometric_small_world_neighbor_set(
        size_t N,
        double k,
        double p,
        vector < double > &r,
        bool use_largest_component,
        size_t seed
        )
{

    assert(k>0);
    assert(N>1);
    assert(p <= 1.);
    assert(p >= 0.);
    assert(k < N);
    assert((r.size()) == N or (r.size() == 0));

    double p0 = k / (k - p*k + p*(N-1.0));
    double p1 = p0 * p;

    double radius = 0.5*k*double(N)/(N-1.0);

    bool list_was_empty = (r.size() == 0);

    //initialize random generators
    mt19937_64 generator;
    if (seed == 0)
        randomly_seed_engine(generator);
    else
        generator.seed(seed);

    uniform_real_distribution<double> random_number(0., 1.);
    uniform_int_distribution<size_t> random_node(0, N-1);
    
    vector < set < size_t > * > G;
    for (size_t node = 0; node < N; ++node)
    {
        G.push_back( new set < size_t >() );
        if (list_was_empty)
            r.push_back( random_number(generator)*double(N) );
    }

    if (list_was_empty)
        sort(r.begin(),r.end());

    // loop over all pairs and draw according to the right probability
    // (this is a lazy slow algorithm running in O(N^2) time
    size_t SR_edges = 0;
    for (size_t i = 0; i < N-1; ++i)
    {
        for (size_t j = i+1; j < N; ++j)
        {
            double distance = fabs(r[j] - r[i]);
            double probability = p1;

            if ((distance <= radius) or ((double(N) - distance) <= radius))
            {
                probability = p0;
                SR_edges++;
            }

            if (random_number(generator) < probability)
            {
                G[i]->insert(j);
                G[j]->insert(i);
            }
        }
    }

    if (use_largest_component)
    {
        get_largest_component(G);
        return G;
    }
    else
    {
        return G;
    }

}
