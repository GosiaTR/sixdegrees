# -*- coding: utf-8 -*-
"""
``sixdegrees`` is a package that allows to sample from various
small-world network models in a standardized manner.

Each model generator is provided a number of nodes `N`,
a mean degree `k` and a structural control parameter 
:math:`-\infty \leq mu \leq 1`.

At :math:`\mu = 1`, each model equals a G(N,p) Erdos-Renyi
random graph (with :math:`p = k/(N-1)`).

At :math:`\mu = 0`, a node's short-range degree equals its
long-range degree. This coincided with the resolution 
threshold for community detection.

At :math:`\mu = -\infty`, the network sample will resemble 
either lattices or random geometric graphs. (Relaistically,
this will already happen for :math:`\mu=-2`.
"""

from .metadata import __version__
from _sixdegrees import *

from .random_geometric_api import *
from .modular_hierarchical_api import *

from .small_world_api import *
