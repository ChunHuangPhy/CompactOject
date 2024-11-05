import numpy as np
from TOVsolver.unit import g_cm_3, dyn_cm_2, MeV, fm

def MITbag_compute_EOS(B): 
    """
    Compute the energy density and pressure based on the given parameters.

    Args:
        B: Input value of bag constant; MeVfm^-3
        
    Returns:
        tuple: Arrays of energy densities in units of gcm^3 and pressures in units of dyncm^2.
    """
    
    B_cgs = B * (MeV / (fm)**3) # converting input to cgs
    energy_density  = np.linspace(4 * B_cgs, 10 * B_cgs, 1000) # cgs
    # epsilon has a minimum value of 4B so that pressure >= 0
    
    pressure = ((energy_density / 3) - (4 * B_cgs / 3))
    
    return energy_density, pressure