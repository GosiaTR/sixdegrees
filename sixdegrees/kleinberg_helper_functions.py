import numpy as np
from scipy.spatial.distance import pdist
from scipy.integrate import quad
from scipy.optimize import brentq

def f(r):
    s = r**2

    if s <= 1:
        result = -4*np.sqrt(s) + np.pi + s
    else:
        result = -2 + 4 * np.arcsin(1/np.sqrt(s)) + 4 * np.sqrt(s-1)\
                 - np.pi - s
        
    return 2 * r * result

def f_PBC(r):
    if r <= 0.5:
        return 2*np.pi*r
    else:
        return 2*np.pi*r - 8*r*np.arccos(0.5/r)

def get_distance_connection_probability_parameters(N, 
                                                   k,
                                                   kappa, 
                                                   epsilon=1e-9, 
                                                   use_periodic_boundary_conditions=True,
                                                   ):

    assert(kappa <= 0)
    np.seterr(all='warn')
    import warnings
    warnings.filterwarnings('error')

    if use_periodic_boundary_conditions:
        distance_density = f_PBC
        rmax = 1/np.sqrt(2)
    else:
        distance_density = f
        rmax = np.sqrt(2)

    f1 = lambda r : distance_density(r)
    f2 = lambda r : distance_density(r) * r**(kappa)

    while True:
        try:
            Ca = (k/(N-1.0) - quad(f1, 0, epsilon)[0] ) / quad(f2, epsilon, rmax)[0]

            if (Ca * epsilon**(kappa)) <= 1.0:
                rmin = epsilon
            else:
                f3 = lambda rm : quad(lambda r: distance_density(r) * (Ca * r**(kappa) - 1.0), epsilon, rm)[0]
                rmin = brentq(f3, epsilon+epsilon/10., rmax)

            warnings.filterwarnings('default')
            return Ca, rmin
        except Warning:
            epsilon *= 2
