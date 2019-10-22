
About this project
==================

`sixdegrees` is a C++/Python-package to sample small-world networks
from models of different flavor - including a modular hierarchical
model with nested community structure. The name is a reference to the
well-known phenomenon of `six degrees of separation`_ which states
that in social networks it might typically take no more than six steps 
to traverse
from one node to another. For problems regarding this statement,
read `this article by D. Watts`_.

Quick example
-------------

Generate a small-world network, a 1D-Kleinberg-network, and a 
self-similar modular hierarchical network with similar properties:

.. code:: python

   import sixdegrees

   B = number_of_submodules_per_modules = 10
   L = number_of_hierarchy_layers = 3

   N = number_of_nodes = B**L
   k = mean_degree = 10
   mu = structural_control_parameter = -0.5

   _, edge_list_SW = sixdegrees.small_world_network(N, k, mu)
   _, edge_list_KL = sixdegrees.kleinberg_network(N, k, mu)
   _, edge_list_MH = sixdegrees.modular_hierarchical_network(B, L, k, mu)


Install
-------

If you get compiling errors, make sure that `pybind11`_ is installed.

.. code:: bash

   git clone https://github.com/benmaier/sixdegrees
   pip install ./sixdegrees

Note that a C++11-compiler has to be installed on the system before
installing ``sixdegrees``.

.. _six degrees of separation: https://en.wikipedia.org/wiki/Six_degrees_of_separation 
.. _this article by D. Watts: https://medium.com/@duncanjwatts/how-small-is-the-world-really-736fa21808ba
.. _pybind11: https://github.com/pybind/pybind11
