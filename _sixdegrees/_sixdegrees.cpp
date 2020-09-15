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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "Utilities.h"
#include "ssmh.h"
#include "kleinberg.h"
#include "modified_small_world.h"
#include "original_small_world.h"
#include "random_geometric_small_world.h"
#include "twoD_random_geometric_kleinberg.h"
#include "find_first_random_edge.h"

using namespace std;
namespace py = pybind11;

PYBIND11_MODULE(_sixdegrees, m) {
    m.doc() = R"pydoc(
        Generate generalized small-world networks, including self-similar modular hierarchical and modified Kleinberg networks.

        .. currentmodule:: _sixdegrees

        Lattice small-world generators
        ------------------------------

        .. autosummary::
            :toctree: _generate

            _modified_small_world_network
            _modified_small_world_network_coord_lists
            _original_small_world_network
            _original_small_world_network_coord_lists

        Modular hierarchical generators
        -------------------------------

        .. autosummary::
            :toctree: _generate

            _modular_hierarchical_network
            _modular_hierarchical_network_coord_lists

        Lattice Kleinberg generators
        ----------------------------

        .. autosummary::
            :toctree: _generate

            kleinberg_network
            kleinberg_network_coord_lists

            
        Random geometric generators
        ---------------------------

        .. autosummary::
            :toctree: _generate

            _random_geometric_small_world_network
            _random_geometric_small_world_network_coord_lists
            _twoD_random_geometric_kleinberg_network
            _twoD_random_geometric_kleinberg_network_coord_lists

        Helper functions
        ----------------

        .. autosummary::
            :toctree: _generate

            fast_gnp
            get_components
            get_kleinberg_connection_probability
    )pydoc";

    m.def("_modular_hierarchical_network", &fast_ssmh_edge_list, "Returns a self-similar modular hierarchical network as an edge list. If you want to compare it to a 1d Kleinberg network, be reminded that mu = ",
            py::arg("B"),
            py::arg("L"),
            py::arg("k"),
            py::arg("xi"),
            py::arg("use_largest_component") = false,
            py::arg("delete_non_largest_component_nodes") = true,
            py::arg("allow_probability_redistribution") = true,
            py::arg("seed") = 0
            );

    m.def("_modular_hierarchical_network_coord_lists", &fast_ssmh_coord_lists, "Returns a self-similar modular hierarchical network as lists of adjacency matrix coordinates.",
            py::arg("B"),
            py::arg("L"),
            py::arg("k"),
            py::arg("xi"),
            py::arg("use_largest_component") = false,
            py::arg("delete_non_largest_component_nodes") = true,
            py::arg("allow_probability_redistribution") = true,
            py::arg("seed") = 0
            );

    m.def("kleinberg_network", &kleinberg_edge_list, "Returns a 1d Kleinberg network (with periodic boundary conditions) as edge list. Connection probability of two nodes u and v is ~ d(u,v)^(mu-1) where d(u,v) is the pair's lattice distance. If you want to map from an ssmh, bear in mind that N=B^L and mu=log(xi)/log(B).",
            py::arg("N"),
            py::arg("k"),
            py::arg("mu"),
            py::arg("use_largest_component") = false,
            py::arg("delete_non_largest_component_nodes") = true,
            py::arg("seed") = 0
            );

    m.def("kleinberg_network_coord_lists", &kleinberg_coord_lists, "Returns a 1d Kleinberg network (with periodic boundary conditions) as lists of adjacency matrix coordinates. Connection probability of two nodes u and v is ~ d(u,v)^(mu-1) where d(u,v) is the pair's lattice distance. If you want to map from an ssmh, bear in mind that N=B^L and mu=log(xi)/log(B).",
            py::arg("N"),
            py::arg("k"),
            py::arg("mu"),
            py::arg("use_largest_component") = false,
            py::arg("delete_non_largest_component_nodes") = true,
            py::arg("seed") = 0
            );

    m.def("_twoD_random_geometric_kleinberg_network", &twoD_random_geometric_kleinberg_edge_list, 
            R"pydoc(Returns a 2d Kleinberg network (with or without periodic boundary conditions) as an edge list, where node positions are uniformly distributed in the real-valued unit square [0,1]^2. You need to provide two lists X, Y of N iid random numbers drawn uniformly from [0,1]. Connection probability of two nodes u and v is ~ d(u,v)^kappa where d(u,v) is the pair's Euclidean distance.)pydoc",
            py::arg("N"),
            py::arg("prefactor"),
            py::arg("rmin"),
            py::arg("kappa"),
            py::arg("X"),
            py::arg("Y"),
            py::arg("periodic_boundary_conditions") = true,
            py::arg("use_largest_component") = false,
            py::arg("delete_non_largest_component_nodes") = true,
            py::arg("seed") = 0
            );

    m.def("_twoD_random_geometric_kleinberg_network_coord_lists", &twoD_random_geometric_kleinberg_coord_lists, 
            R"pydoc(Returns a 2d Kleinberg network (with or without periodic boundary conditions) as lists of adjacency matrix coordinates, where node positions are uniformly distributed in the real-valued unit square [0,1]^2. You need to provide two lists X, Y of N iid random numbers drawn uniformly from [0,1]. Connection probability of two nodes u and v is ~ d(u,v)^kappa where d(u,v) is the pair's Euclidean distance.)pydoc",
            py::arg("N"),
            py::arg("prefactor"),
            py::arg("rmin"),
            py::arg("kappa"),
            py::arg("X"),
            py::arg("Y"),
            py::arg("periodic_boundary_conditions") = true,
            py::arg("use_largest_component") = false,
            py::arg("delete_non_largest_component_nodes") = true,
            py::arg("seed") = 0
            );


    m.def("_original_small_world_network", &original_small_world_edge_list, "Returns a Watts-Strogatz small world network as a pair of (number_of_nodes, edge_list). Note that at p = 1, the generated networks are not equal to Erdos-Renyi graphs. The degree k has to be an even integer.",
            py::arg("N"),
            py::arg("k"),
            py::arg("p"),
            py::arg("use_largest_component") = false,
            py::arg("delete_non_largest_component_nodes") = true,
            py::arg("seed") = 0
            );

    m.def("_original_small_world_network_coord_lists", &original_small_world_coord_lists, "Returns a Watts-Strogatz small world network as lists of adjacency matrix coordinates. Note that at p = 1, the generated networks are not equal to Erdos-Renyi graphs. The degree k has to be an even integer",
            py::arg("N"),
            py::arg("k"),
            py::arg("X"),
            py::arg("use_largest_component") = false,
            py::arg("delete_non_largest_component_nodes") = true,
            py::arg("seed") = 0
            );

    m.def("_modified_small_world_network", &modified_small_world_edge_list, "Returns a variant of the Watts-Strogatz small world network model as a pair of (number_of_nodes, edge_list). In this variant, at p = 1, the generated networks are equal to Erdos-Renyi graphs. The degree k has to be an even integer.",
            py::arg("N"),
            py::arg("k"),
            py::arg("beta"),
            py::arg("use_largest_component") = false,
            py::arg("delete_non_largest_component_nodes") = true,
            py::arg("use_fast_algorithm") = true,
            py::arg("seed") = 0
            );

    m.def("_modified_small_world_network_coord_lists", &modified_small_world_coord_lists, "Returns a variant of the Watts-Strogatz small world network model as lists of adjacency matrix coordinates. In this variant, at p = 1, the generated networks are equal to Erdos-Renyi graphs. The degree k has to be an even integer",
            py::arg("N"),
            py::arg("k"),
            py::arg("beta"),
            py::arg("use_largest_component") = false,
            py::arg("delete_non_largest_component_nodes") = true,
            py::arg("use_fast_algorithm") = true,
            py::arg("seed") = 0
            );

    m.def("_random_geometric_small_world_network", &random_geometric_small_world_edge_list, 
            R"pydoc(Returns a variant of the Watts-Strogatz small world network model as a pair of (number_of_nodes, edge_list). In this variant, at p = 1, the generated networks are equal to Erdos-Renyi graphs and at p = 0, they are equal to random geometric graphs of dimension d = 1. The mean degree k can be a float.)pydoc",
            py::arg("N"),
            py::arg("k"),
            py::arg("beta"),
            py::arg("X"),
            py::arg("use_largest_component") = false,
            py::arg("delete_non_largest_component_nodes") = true,
            py::arg("seed") = 0
            );

    m.def("_random_geometric_small_world_network_coord_lists", &random_geometric_small_world_coord_lists, 
            R"pydoc(Returns a variant of the Watts-Strogatz small world network model as lists of adjacency matrix coordinates. In this variant, at p = 1, the generated networks are equal to Erdos-Renyi graphs and at p = 0, they are equal to random geometric graphs of dimension d = 1. The mean degree k can be a float.)pydoc",
            py::arg("N"),
            py::arg("k"),
            py::arg("beta"),
            py::arg("X"),
            py::arg("use_largest_component") = false,
            py::arg("delete_non_largest_component_nodes") = true,
            py::arg("seed") = 0
            );

    m.def("fast_gnp", &fast_gnp, "Returns a G(N,p) random graph in O(N+m) time as described by Batagelj and Brandes (in edge list format).",
            py::arg("N"),
            py::arg("p"),
            py::arg("start_node") = 0,
            py::arg("seed") = 0
           );

    m.def("get_components", &get_components_from_edgelist, "Get a list of sets. Each list entry is a set of nodes which makes up one component of the graph.",
            py::arg("N"),
            py::arg("edge_list")
          );

    m.def("get_kleinberg_connection_probability", &get_kleinberg_pmf, "get a list of connection probabilities by neighbor distance",
            py::arg("N"),
            py::arg("k"),
            py::arg("mu")
        );

    m.def("find_first_random_edge", &find_first_random_edge, "Given a ``numpy.ndarray`` of probabilities in descending order, find the first entry `j` for which :math:`\\sum_i^j(p_i-1)<0`.",
            py::arg("probabilities")
        );

}
