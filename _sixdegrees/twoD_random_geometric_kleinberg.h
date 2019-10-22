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
#ifndef __TWOD_RGG_KLEIN__
#define __TWOD_RGG_KLEIN__

#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <set>
#include <utility>
#include <cmath>
#include <numeric>
#include <random>
#include <ctime>
#include <cstdlib>
#include <tuple>

using namespace std;

double dist(double const &x, double const &y);

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
        );

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
        );

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
                            );
#endif
