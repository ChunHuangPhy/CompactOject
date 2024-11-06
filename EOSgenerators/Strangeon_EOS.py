import numpy as np
import math
from scipy import optimize
from TOVsolver.unit import g_cm_3, dyn_cm_2

def compute_EOS(n, theta): 
    """
    Compute the energy density and pressure based on the given parameters.

    Args:
        n (array): An array of n values. Input values of baryon number densities.
        theta (array): An array representing the parameters [epsilon, Nq, ns].
        epsilon: the depth of the potential well; MeV;
        Nq: the number of quarks in a strangeon; 
        ns: the number density of baryons at the surface of the star; fm^-3
        
    Returns:
        tuple: Arrays of energy densities in units of gcm3 and pressures in units of dyncm2.
    """
    
    Nq, epsilon, ns = theta
    
    A12 = 6.2
    A6 = 8.4 
    mq = 300 
    """
    mq: the mass of the quark in this EOS.
    A12 and A6 are fixed throughout the calculation.
    """
   
    sigma = np.sqrt(A6 / (2 * A12)) * (Nq / (3 * ns)) 
   
    energy_density = 2 * epsilon * (A12 * sigma**4 * n**5 - A6 * sigma**2 * n**3) + n * Nq * mq
    pressure = 4 * epsilon * (2 * A12 * sigma**4 * n**5 - A6 * sigma**2 * n**3)
    
    return energy_density , pressure