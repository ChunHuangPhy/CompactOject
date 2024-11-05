import numpy as np
import math
from scipy import optimize
from TOVsolver.unit import g_cm_3, dyn_cm_2

def compute_EOS(rho, theta):
    """
    Calculate the pressure of a neutron star based on density using a piecewise polytropic equation of state (EOS).

    This function computes the pressure (`pres`) as a function of density (`rho`) by applying different polytropic indices
    (`gamma1`, `gamma2`, `gamma3`) within specified density thresholds (`rho_t1`, `rho_t2`). The EOS is defined in three
    distinct regions:
    
    - **Low-density region:** `rho <= rho_t1`
    - **Intermediate-density region:** `rho_t1 < rho <= rho_t2`
    - **High-density region:** `rho > rho_t2`

    Parameters
    ----------
    rho : array-like
        An array of density values (in cgs units) at which to calculate the pressure.
    
    theta : array-like, length 5
        A list or tuple containing the EOS parameters in the following order:
        
        - `gamma1` (float): Polytropic index for the low-density region.
        - `gamma2` (float): Polytropic index for the intermediate-density region.
        - `gamma3` (float): Polytropic index for the high-density region.
        - `rho_t1` (float): Density threshold between the low and intermediate-density regions (in cgs units).
        - `rho_t2` (float): Density threshold between the intermediate and high-density regions (in cgs units).
    
    Returns
    -------
    pres : ndarray
        An array of pressure values (in cgs units) corresponding to the input density values.
    """
    gamma1, gamma2, gamma3, rho_t1, rho_t2 = theta
    c = 2.99792458E10 # cgs
    G = 6.6730831e-8 # cgs
    rho_ns = 267994004080000.03 #cgs
    rho_t = 4.3721E11*G/c**2
    P_t = 7.7582E29* G / c**4
    
    P_ts, k = np.zeros(3), np.zeros(3)
    P_ts[0] = P_t
    k[0] = P_t / ((rho_t / rho_ns)**gamma1)
    P_ts[1] = k[0] * rho_t1**gamma1
    k[1] = P_ts[1] / (rho_t1**gamma2)
    P_ts[2] = k[1] * rho_t2**gamma2
    k[2] = P_ts[2] / (rho_t2**gamma3)

    # Calculate the pressure for the entire input array `rho`
    pres = np.where(rho <= rho_t1, k[0] * rho**gamma1,
                    np.where(rho <= rho_t2, k[1] * rho**gamma2, k[2] * rho**gamma3))

    return pres
