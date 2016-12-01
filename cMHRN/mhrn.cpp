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
#include "mhrn.h"

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

using namespace std;

pair < vector <size_t>, vector<size_t> > fast_mhrn_coord_lists(
        size_t B,
        size_t L,
        double k,
        double xi,
        bool use_giant_component,
        size_t seed
        )
{
    vector < set < size_t > * > G = fast_mhrn_neighbor_set(B,L,k,xi,use_giant_component,seed);
    size_t N = pow(B,L);
    vector < size_t > rows;
    vector < size_t > cols;

    for(size_t u = 0; u < N; u++)
        for( auto const& v: *G[u] )
        {
            rows.push_back(u);
            cols.push_back(v);
        }
    
    return make_pair(rows,cols);
}

vector < pair < size_t, size_t > > fast_mhrn_edge_list(
        size_t B,
        size_t L,
        double k,
        double xi,
        bool use_giant_component,
        size_t seed
        )
{
    size_t N = pow(B,L);
    vector < set < size_t > * > G = fast_mhrn_neighbor_set(B,L,k,xi,use_giant_component,seed);
    vector < pair < size_t, size_t > > edge_list;

    for(size_t u = 0; u < N; u++)
        for( auto const& v: *G[u] )
        {
            if (u<v)
            {
                edge_list.push_back( make_pair(u,v) );
            }
        }
    
    return edge_list;
}

vector < set < size_t > * > fast_mhrn_neighbor_set(
        size_t B,
        size_t L,
        double k,
        double xi,
        bool use_giant_component,
        size_t seed
        )
{
    //initialize random generators
    default_random_engine generator(seed);
    uniform_real_distribution<double> uni_distribution(0.,1.);
    //binomial_distribution<int> distribution(9,0.5)
    
    size_t N = pow(B,L);

    vector < set < size_t > * > G;
    for(size_t node=0; node<N; node++)
    {
        G.push_back( new set <size_t> );
    }

    double p1;
    if ( xi == 1.0 ) 
        p1 = k / double(B-1) / double(L);
    else
        p1 = k / double(B-1) * (1.0-xi) / (1.0-pow(xi,L));

    if (p1>1.0)
        throw domain_error("The lowest layer connection probability is >1.0, meaning that either xi is too small or k is too large.");

    vector < double > p(L);
    for (size_t l=0; l<L; l++)
        p[l] = p1 * pow(xi/double(B),l);

    //add ER graphs in lowest layer
    for (size_t start_node = 0; start_node<N; start_node += B)
    {
        add_random_subgraph(B,p1,G,generator,uni_distribution,start_node);
    }

    for(size_t l=2;l<=L; l++)
    {
        binomial_distribution<size_t> binom( size_t(0.5*pow(B,(L+l-1))*(B-1)),p[l-1]);
        size_t current_m_l = binom(generator);

        for (size_t m = 0; m< current_m_l; m++)
        {
            size_t B_l = pow(B,l);
            size_t B_lm1 = pow(B,l-1);
            while (true)
            {
                size_t w = size_t(uni_distribution(generator)*N);
                size_t b = w / B_l;
                size_t b_lower = w / B_lm1;
                size_t v;
                uniform_int_distribution<> randint(b*B_l,(b+1)*B_l-1);

                while (true)
                {
                    v = randint(generator);
                    if (b_lower != v / B_lm1)
                        break;
                }

                pair <size_t,size_t> edge = get_sorted_pair(w,v);
                if ( G[w]->find(v) == G[w]->end() )
                {
                    G[w]->insert(v);
                    G[v]->insert(w);
                    break;
                }
            }
        }
    }

    if (use_giant_component)
    {
        G = get_giant_component(G);
        return G;
    }
    else
    {
        return G;
    }

}

vector < pair < size_t, size_t > > fast_gnp(
        size_t N_,
        size_t p,
        size_t start_node,
        size_t seed
        )
{
    size_t N = N_ + start_node;

    //initialize random generators
    default_random_engine generator(seed);
    uniform_real_distribution<double> uni_distribution(0.,1.);

    vector < set < size_t > * > G;
    for(size_t node=0; node<N; node++)
        G.push_back( new set <size_t> );

    add_random_subgraph(N,p,G,generator,uni_distribution,start_node);
    vector < pair < size_t, size_t > > edge_list;

    for(size_t u = 0; u < N; u++)
        for( auto const& v: *G[u] )
        {
            if (u<v)
            {
                edge_list.push_back( make_pair(u,v) );
            }
        }
    
    return edge_list;
}