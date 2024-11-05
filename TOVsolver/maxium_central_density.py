import numpy as np
from scipy.interpolate import interp1d

from TOVsolver import unit
from TOVsolver.main import OutputMR

def maxium_central_density(energy_density, pressure, central_densitys=np.logspace(14.3, 15.6, 5) * unit.g_cm_3, num2=30):
    """Outputs the maxium central density of a stable EoS
    Args:
        energy_density (numpy 1Darray): Density of EoS
        pressure (numpy 1Darray): Pressure of EoS
        central_densitys (numpy 1Darray): The range of central density
        num2 (int): The number of segments in the density interval of second search. At least 5 points
    Returns:
        density (float): The maxium central density, in unit.g_cm_3
    """
    ############## Below is the main part of two searches for peaks
    ######## First search
    Ms = [-1,-2] # Meaningless initialization, only to ensure that the following 'if (i>0) and (Ms [-1]<=Ms [-2])' statements can run properly
    store_d_range = [central_densitys[-2], central_densitys[-1]] # When the following loop does not output a result, initialization here will have its meaning
    # Find the extremum point within the predetermined range of central density and return the central density near the extremum point
    for i, rho in enumerate(central_densitys):
        M, R = OutputMR('', energy_density, pressure, [rho])[0]
        Ms.append(M)
        if (i>0) and (Ms[-1] <= Ms[-2]):
            store_d_range = [central_densitys[i-2], central_densitys[i]] 
            break
    Ms_larg = Ms[-1] # Used to store the mass corresponding to the maximum density during the first peak search
    
    ######## Second search
    Ms = [-1,-2] # Reinitialize
    # Update and refine the central density range, and ultimately return the central density of the extremum points
    store_d_range = np.geomspace(store_d_range[0], store_d_range[1], num2) # At least 5 points
    # Note that the first and last points have already been calculated, so there is no need to traverse them
    for i, rho in enumerate(store_d_range[1:-1]):
        M, R = OutputMR('', energy_density, pressure, [rho])[0]
        Ms.append(M)
        if Ms[-1] <= Ms[-2]:    # Note that due to the second peak search refining the interval, the result is generally not obtained at the first point.
                                # Therefore, initializing Ms=[-1, -2] is acceptable
            return store_d_range[1:][i-1]
    # When the above peak search fails, continue comparing the last point
    if Ms_larg < Ms[-1]:
        return store_d_range[-2]
    else:
        return store_d_range[-1]