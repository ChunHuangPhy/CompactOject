# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 14:23:37 2024

@author: wlyua
"""
import numpy as np



c = 3e10  
G = 6.67428e-8  
Msun = 1.989e33 
dyncm2_to_MeVfm3 = 1./(1.6022e33)  
gcm3_to_MeVfm3 = 1./(1.7827e12)  
oneoverfm_MeV = 197.33 



def Strangeon_compute_EOS(n, theta): 
    """
    Compute the energy density and pressure for Strangeon matter based on the given parameters.

    Args:
        n (array): An array of n values. Input values of strangeon number densities.
        theta (array): An array representing the parameters [epsilon, ns].
        epsilon: the depth of the potential well;
        ns: the number density of baryons at the surface of the star.
        Nq: the number of quarks in a strangeon; An integer input.
        
    Returns:
        tuple: Arrays of energy densities and pressures.
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
    
    return energy_density*G/c**2/gcm3_to_MeVfm3, pressure*G/c**4/dyncm2_to_MeVfm3


