import numpy as np
from scipy.spatial.distance import pdist
from scipy.integrate import quad
from scipy.optimize import brentq

def f(r):
    """Distance probability distribution of pairs of uniformly distributed random points in [0,1]^2."""
    s = r**2
    result = np.zeros_like(s)

    ndx1 = np.where(s<=1)
    s1 = s[ndx1]
    result[ndx1] = -4*np.sqrt(s1) + np.pi + s1
    ndx2 = np.where(r>1)
    s2 = s[ndx2]
    result[ndx2] = -2 + 4 * np.arcsin(1/np.sqrt(s2)) + 4 * np.sqrt(s2-1)\
                   - np.pi - s2
        
    return 2 * r * result

def f_(r):
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

def redistribute_edges(bins, geo_counts):
    """Given a list of probabilities and a list of edge counts, redistribute probability such that
    the total amount of generated edges is constant."""
    missing_edges = 0.0
    for i, b in enumerate(bins):
        m_count = geo_counts[i]
        if b > 1.0:
            missing_edges += (b - 1.0) * m_count
            bins[i] = 1.0
        else:
            bins[i] += missing_edges / m_count
            missing_edges = 0.0
            if bins[i] > 1.0:
                missing_edges = (bins[i] - 1.0) * m_count
                bins[i] = 1.0
        
        if missing_edges == 0.0:
            break
        
    return bins
            
def get_distance_connection_probability_parameters(N, 
                                                   k,
                                                   kappa, 
                                                   epsilon=1e-9, 
                                                   use_periodic_boundary_conditions=False,
                                                   ):

    np.seterr(all='warn')
    import warnings
    warnings.filterwarnings('error')

    if use_periodic_boundary_conditions:
        distance_density = f_PBC
        rmax = 1/np.sqrt(2)
    else:
        distance_density = f_
        rmax = np.sqrt(2)

    f1 = lambda r : distance_density(r)
    f2 = lambda r : distance_density(r) * r**(-kappa)

    while True:
        try:
            Ca = (k/(N-1.0) - quad(f1, 0, epsilon)[0] ) / quad(f2, epsilon, rmax)[0]

            if (Ca * epsilon**(-kappa)) <= 1.0:
                rmin = epsilon
            else:
                f3 = lambda rm : quad(lambda r: distance_density(r) * (Ca * r**(-kappa) - 1.0), epsilon, rm)[0]
                rmin = brentq(f3, epsilon+epsilon/10., rmax)

            warnings.filterwarnings('default')
            return Ca, rmin
        except Warning:
            epsilon *= 2
