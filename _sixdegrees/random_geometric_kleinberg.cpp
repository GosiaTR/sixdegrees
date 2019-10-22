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
#include "random_geometric_kleinberg.h"

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

tuple < size_t, vector <size_t>, vector<size_t> > random_geometric_kleinberg_coord_lists(
        size_t N,
        double k,
        double mu,
        vector < double > r,
        bool use_largest_component,
        bool delete_non_largest_component_nodes,
        bool use_theory_algorithm,
        size_t seed,
        double epsilon
        )
{
    vector < set < size_t > * > G;

    if (use_theory_algorithm)
        G = theoretical_random_geometric_kleinberg_neighbor_set(N,k,mu,r,use_largest_component,seed,epsilon);
    else
        G = random_geometric_kleinberg_neighbor_set(N,k,mu,r,use_largest_component,seed);

    vector < size_t > rows;
    vector < size_t > cols;
    size_t new_N = neighbor_set_to_coord_lists(G,rows,cols,use_largest_component,delete_non_largest_component_nodes);

    return make_tuple(new_N,rows,cols);
}
    
pair < size_t, vector < pair < size_t, size_t > > > random_geometric_kleinberg_edge_list(
        size_t N,
        double k,
        double mu,
        vector < double > r,
        bool use_largest_component,
        bool delete_non_largest_component_nodes,
        bool use_theory_algorithm,
        size_t seed,
        double epsilon
        )
{
    vector < set < size_t > * > G;

    if (use_theory_algorithm)
        G = theoretical_random_geometric_kleinberg_neighbor_set(N,k,mu,r,use_largest_component,seed,epsilon);
    else
        G = random_geometric_kleinberg_neighbor_set(N,k,mu,r,use_largest_component,seed);

    vector < pair < size_t, size_t > > edge_list;
    size_t new_N = neighbor_set_to_edge_list(G,edge_list,use_largest_component,delete_non_largest_component_nodes);
    
    return make_pair(new_N,edge_list);
}

vector < set < size_t > * > random_geometric_kleinberg_neighbor_set(
        size_t N,
        double k,
        double mu,
        vector < double > &r,
        bool use_largest_component,
        size_t seed
        )
{

    assert(k>0);
    assert(N>1);
    assert(mu <= 1.);
    assert(k < N);
    assert((r.size() == N) or (r.size()==0));

    bool list_was_empty = (r.size() == 0);

    double kappa = -(1.0-mu);

    assert(kappa <= 0.0);

    //initialize random generators
    mt19937_64 generator;
    if (seed == 0)
        randomly_seed_engine(generator);
    else
        generator.seed(seed);
   
    uniform_real_distribution<double> random_number(0., 1.);
    
    vector < set < size_t > * > G;
    for (size_t node = 0; node < N; ++node)
    {
        G.push_back( new set < size_t >() );
        if (list_was_empty)
            r.push_back( random_number(generator)*double(N) );
    }

    vector < edge_distance > distances;

    // loop over all pairs and compute
    // (this is a lazy slow algorithm running in O(N^2) time
    for (size_t i = 0; i < N-1; ++i)
    {
        for (size_t j = i+1; j < N; ++j)
        {
            double distance = fabs(r[j] - r[i]);
            distance = min(distance, N-distance);
            distances.push_back(make_tuple( make_pair(i,j), distance, 0.0));
        }
    }

    // sort so one can properly assign the short-range edges
    sort(distances.begin(), distances.end(), compare_distance);

    // sum up in order to properly norm the distribution.
    // start the summation from the last entries as these will produce
    // the smallest weights and hence the whole sum will have less accumulation error.
    double norm = 0.0;
    for(vector<edge_distance>::reverse_iterator edge = distances.rbegin(); edge != distances.rend(); ++edge)
    {
        double weight = pow(get<1>(*edge), kappa);
        get<2>(*edge) = weight; // save r^kappa for later use
        norm += weight;
    }

    // compute prefactor
    double C = 0.5 * double(N) * k  / norm;

    // iterate through edges
    auto edge = distances.begin();
    double excess_edges = C * get<2>(*edge) - 1.0;
    size_t l = 0;


    // all distances where the total number of excess edges still exceeds zero 
    // are naturally produced.
    while ((excess_edges > 0.0) and (edge != distances.end()))
    {
        size_t const &i = (get<0>(*edge)).first;
        size_t const &j = (get<0>(*edge)).second;
        G[i]->insert(j);
        G[j]->insert(i);

        ++edge;
        ++l;
        if (edge != distances.end())
            excess_edges += C * get<2>(*edge) - 1.0;
    }

    // all other edges are only produced with probability C*r^kappa
    while (edge != distances.end())
    {
        const double &probability = C * get<2>(*edge);

        if (random_number(generator) < probability)
        {
            size_t const &i = (get<0>(*edge)).first;
            size_t const &j = (get<0>(*edge)).second;
            G[i]->insert(j);
            G[j]->insert(i);
        }

        ++edge;
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

// Compares a tuple of edges according to the distance of their comprising nodes
bool compare_distance(
                    const edge_distance &p1,
                    const edge_distance &p2
                 )
{
    return (get<1>(p1)) < (get<1>(p2));
}

// ============================= THESE ARE THE FLAWED METHODS THAT NEED SOME ADDITIONAL THINKING =============================

vector < set < size_t > * > theoretical_random_geometric_kleinberg_neighbor_set(
        size_t N,
        double k,
        double mu,
        vector < double > &r,
        bool use_largest_component,
        size_t seed,
        double epsilon
        )
{

    assert(k>0);
    assert(N>1);
    assert(mu <= 1.);
    assert(k < N);
    assert(epsilon >= 0.);
    assert((r.size() == N) or (r.size()==0));

    bool list_was_empty = (r.size() == 0);

    // convert to proper power-law convention ~ r^(-kappa)
    double kappa = -(mu - 1.0);

    // the distribution of distances of value r in the interval [0,N/2]
    // is given as
    //                         f(r) = 2/N.
    // Hence the distribution of possible edges of distance r in the
    // interval [0,N/2] is
    //                 M(r) = N/2(N-1) f(r) = N-1.
    // Therefore, on average, there is a single produced edge in the 
    // interval [0,1/N-1].
    //
    // We just arbitrarily say that we want to reduce the number of erroneously 
    // assigned edges to one every 1000 samples and therefore we set
    //                 epsilon = 1/(N-1) * 1/1000
    //
    // (Remember that all distances smaller than epsilon are automati-
    // cally assigned an edge).
    //
    if (epsilon == 0.0)
        epsilon = 1e-3/(N-1);

    // this function might possibly change epsilon if it's too small
    // for the evaluation to converge
    double Rmin = find_R_min(N, k, kappa, epsilon);

    double C;
    if (kappa != 1.0)
        C = (1.0-kappa)*(N/(N-1.0)*k/2.0 - epsilon) / (pow(0.5*N,1-kappa) - pow(epsilon,1-kappa));
    else
        C = (N/(N-1.0)*k/2.0 - epsilon) / (log(0.5*N) - log(epsilon));

    double p_Erdos_Renyi = k/(N-1.0);

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

    // sort position so one knows that the node position 0 is the
    // first in [0,N] and so forth
    if (list_was_empty)
        sort(r.begin(), r.end());

    // loop over all pairs and draw according to the right probability
    // (this is a lazy slow algorithm running in O(N^2) time
    for (size_t i = 0; i < N-1; ++i)
    {
        for (size_t j = i+1; j < N; ++j)
        {
            double distance = fabs(r[j] - r[i]);
            double probability;

            if (mu < 1.0)
            {
                if ((distance <= Rmin) or ((double(N) - distance) <= Rmin))
                    probability = 1.0;
                else
                    probability = C*pow(distance,-kappa);
            }
            else
            {
                probability = p_Erdos_Renyi;
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

double f1(const size_t &N, const double &k, const double &kappa, const double &epsilon, double R)
{
    return ( N/(N-1.0)*k/2.0 - epsilon) * (pow(R,1-kappa)-pow(epsilon,1-kappa)) / (pow(0.5*N,1-kappa) - pow(epsilon,1-kappa)) + (R-epsilon);
}

double f1P(const size_t &N, const double &k, const double &kappa, const double &epsilon, double R)
{
    return ( N/(N-1.0)*k/2.0 - epsilon) * (pow(R,-kappa)) / (1.0-kappa) / (pow(0.5*N,1-kappa) - pow(epsilon,1-kappa)) + 1.0;
}

double f2(const size_t &N, const double &k, const double &epsilon, double R)
{
    return ( N/(N-1.0)*k/2.0 - epsilon) * (log(R)-log(epsilon)) / (log(0.5*N)-log(epsilon)) + (R-epsilon);
}

double f2P(const size_t &N, const double &k, const double &epsilon, double R)
{
    return ( N/(N-1.0)*k/2.0 - epsilon) * 1.0/R / (log(0.5*N)-log(epsilon)) + 1.0;
}

// BEWARE: THIS FUNCTION MAY POSSIBLY CHANGE EPSILON IF IT'S TOO SMALL.
double find_R_min(const size_t &N, const double &k, const double &kappa, double &epsilon)
{
    // for the rationale of choosing this initial value see the comment in the
    // network sampling function. Basically, we assume that R is large
    // enough such that at least a single edge in the network is definitely assigned.
    double x1 = 1/double(N-1);
    double x, fx, fx1;

    // machine eps for newton-raphson
    double eps = 1e-14;
    do
    {
        x = x1;
        if (kappa != 1.0)
        {
            fx = f1(N,k,kappa,epsilon,x);
            fx1 = f1P(N,k,kappa,epsilon,x);
        }
        else
        {
            fx = f2(N,k,epsilon,x);
            fx1 = f2P(N,k,epsilon,x);
        }
           
        x1 = x - (fx/fx1);

    } while(fabs(x1-x) >= eps);

    return x1;
}

