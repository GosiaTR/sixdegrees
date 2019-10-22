/* 
 * The MIT License (MIT)
 * Copyright (c) 2018, Benjamin Maier
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

#include "twoD_random_geometric_kleinberg.h"
#include "Utilities.h"
#include <math.h>


using namespace std;

double dist(double const &x, double const &y)
{
    return sqrt(x*x+y*y);
}

tuple < size_t, vector <size_t>, vector<size_t> > twoD_random_geometric_kleinberg_coord_lists(
        size_t N, 
        double prefactor,
        double rmin,
        double kappa,
        vector < double > & X,
        vector < double > & Y,
        bool PBC,
        bool use_largest_component,
        bool delete_non_largest_component_nodes,
        size_t seed
        )
{
    vector < set < size_t > * > G = 
                twoD_random_geometric_kleinberg_neighbor_set(
                                N, 
                                prefactor,
                                rmin,
                                kappa,
                                X,
                                Y,
                                PBC,
                                use_largest_component,
                                seed
                                );

    vector < size_t > rows;
    vector < size_t > cols;
    size_t new_N = neighbor_set_to_coord_lists(G,rows,cols,use_largest_component,delete_non_largest_component_nodes);

    return make_tuple(new_N,rows,cols);
}

pair < size_t, vector < pair < size_t, size_t > > > twoD_random_geometric_kleinberg_edge_list(
        size_t N, 
        double prefactor,
        double rmin,
        double kappa,
        vector < double > & X,
        vector < double > & Y,
        bool PBC,
        bool use_largest_component,
        bool delete_non_largest_component_nodes,
        size_t seed
        )
{
    vector < set < size_t > * > G = 
                twoD_random_geometric_kleinberg_neighbor_set(
                                N, 
                                prefactor,
                                rmin,
                                kappa,
                                X,
                                Y,
                                PBC,
                                use_largest_component,
                                seed
                                );

    vector < pair < size_t, size_t > > edge_list;
    size_t new_N = neighbor_set_to_edge_list(G,edge_list,use_largest_component,delete_non_largest_component_nodes);

    return make_pair(new_N,edge_list);
}

vector < set < size_t > * > twoD_random_geometric_kleinberg_neighbor_set(
                            size_t N, 
                            double prefactor,
                            double rmin,
                            double kappa,
                            vector < double > & X,
                            vector < double > & Y,
                            bool PBC,
                            bool use_largest_component,
                            size_t seed
                            )
{

    assert(kappa <= 0);
    assert(X.size() == N);
    assert(Y.size() == N);
    assert(prefactor > 0);
    assert(rmin >= 0.0);
    assert(rmin <= 1.0);

    vector < set < size_t > * > G;
    for(size_t node=0; node<N; node++)
    {
        G.push_back( new set <size_t> );
    }

    mt19937_64 generator;
    if (seed == 0)
        randomly_seed_engine(generator);
    else
        generator.seed(seed);
    uniform_real_distribution<double> uni_random(0.,1.);

    for(size_t source = 0; source < N-1; ++source)
    {
        double xs = X[source];
        double ys = Y[source];

        for (size_t target = source+1; target < N; ++target)
        {
            double xt = X[target];
            double yt = Y[target];

            if (PBC)
            {
                if (xs-xt > 0.5)
                    xt += 1;
                else if (xs-xt < -0.5)
                    xt -= 1;

                if (ys-yt > 0.5)
                    yt += 1;
                else if (ys-yt < -0.5)
                    yt -= 1;
            }

            double r = dist(xs-xt,ys-yt);

            if (r <= rmin) 
            {
                const size_t &s = source;
                const size_t &t = target;
                G[t]->insert(s);
                G[s]->insert(t);
            }
            else
            {
                double p = prefactor * pow(r, kappa);

                if (uni_random(generator) < p)
                {
                    const size_t &s = source;
                    const size_t &t = target;
                    G[t]->insert(s);
                    G[s]->insert(t);
                }
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
