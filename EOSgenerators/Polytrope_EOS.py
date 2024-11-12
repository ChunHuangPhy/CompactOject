import numpy as np
from scipy.optimize import root

from TOVsolver.unit import g_cm_3, dyn_cm_2

dyncm2_to_MeVfm3 = 1. / (1.6022e33)
gcm3_to_MeVfm3 = 1. / (1.7827e12)
oneoverfm_MeV = 197.33
rho_ns = 267994004080000.03

def compute_EOS(rhos, theta,
              rho_t_start=4.3721E11*g_cm_3, P_t_start=7.7582E29*dyn_cm_2):
    """
    Calculate the pressure of a neutron star based on density using a piecewise polytropic equation of state (EOS).

    This function computes the pressure (`pres`) as a function of density (`rho`) by applying different polytropic indices
    (`gamma_1`, `gamma_2`, ..., `gamma_n`) within specified density thresholds (`rho_t_start`, `rho_t_1`, `rho_t_2`, ..., `rho_t_n-1`).
    The EOS is defined in n distinct regions:
    
    - **Low-density region:** `rho_t_start < rho <= rho_t_1`
    - **Intermediate-density region:** `rho_t_1 < rho <= rho_t_2`, ..., `rho_t_k < rho <= rho_t_k+1`, ..., `rho_t_n-2 < rho <= rho_t_n-1`
    - **High-density region:** `rho > rho_t_n-1`

    Parameters
    ----------
    rho : array-like
        An array of density values (in cgs units) at which to calculate the pressure.
    
    theta : array-like, length `2 * n - 1`
        A list or tuple containing the EOS parameters in the following order:
        
        - `gamma_1` (float): Polytropic index for the low-density region.
        - `gamma_2` to `gamma_2_n-1` (float): Polytropic index for the intermediate-density region.
        - `gamma_n` (float): Polytropic index for the high-density region.
        - `rho_t_1` (float): Density threshold between the low and intermediate-density regions (in cgs units).
        - `rho_t_2` to `rho_t_n-2` (float): Density for the intermediate regions (in cgs units).
        - `rho_t_n-1` (float): Density threshold between the intermediate and high-density regions (in cgs units).
    
    rho_t_start : float
        Start point of the density for polytropic EOS

    P_t_start : float
        Start point of the pressure for polytropic EOS
    Returns
    -------
    pres : ndarray
        An array of pressure values (in cgs units) corresponding to the input density values.
    """
    
    n = len(theta) // 2 + 1
    gammas = theta[:n]
    rho_ts = theta[n:]
    
    P_ts, ks = np.zeros((2, n))
    rho_ts = np.insert(rho_ts, 0, rho_t_start)

    # Calculate the values of P_t and k, based on the parameter gamma and rho_t
    for i in range(len(ks)):
        if i == 0:
            P_ts[i] = P_t_start
            ks[i] = P_ts[i] / (rho_ts[i]**gammas[i])
        else:
            P_ts[i] = ks[i-1] * (rho_ts[i]**gammas[i-1])
            ks[i] = P_ts[i] / (rho_ts[i]**gammas[i])

    # Construct judgement criteria for each section
    conds = [np.array(rhos > rho_ts[i]) & np.array(rhos <= rho_ts[i+1]) for i in range(n-1)]
    conds.append(rhos > rho_ts[-1])
    
    # Build polytrope functions to compute the pressure for each section
    functions = [lambda rho, k=ks[i], g=gammas[i]: k * (rho ** g) for i in range(n)]
    
    # Calculate pressure by using np.piecewise, based on the conds and functions we defined above.
    pres = np.piecewise(rhos, conds, functions)

    return pres

def fun_gamma_max(rho2, rho1, p1):
    """Outputs the maximum gamma for given densities and pressure at both ends of a section of polytropic function.
    Args:
        rho2 (float): density at the end point.
        rho1 (float): density at the start point.
        p1 (float): pressure at the start point.
    Returns:
        gamma_max (float): the maximum gamma.
    """
    def fun(x):
        a1 = np.log(rho2/p1/x)
        a2 = np.log(rho2/rho1)
        return a1 / a2 - x

    gamma_max = root(fun, [3/4]).x[0]
    return gamma_max

def eos_core_pp(gammas,rho_ts,rho_t,rho, P_t):
    P_ts = np.zeros(len(gammas))
    k = np.zeros(len(gammas))
    P_ts[0] = P_t
    k[0] = P_t / ((rho_t / rho_ns)**gammas[0])
    P_ts[1] = k[0] * rho_ts[0]**gammas[0]
    k[1] = P_ts[1] / (rho_ts[0]**gammas[1])
    P_ts[2] = k[1] * rho_ts[1]**gammas[1]
    k[2] = P_ts[2] / (rho_ts[1]**gammas[2])
    pres_ts = P_ts[1::] / dyncm2_to_MeVfm3

    if rho <= rho_ts[0]:
        pres = k[0] * rho**gammas[0]
    
    if rho_ts[0]< rho <= rho_ts[1]:
        #pres = k[1] * rho**gammas[1]
        pres = k[1] * (rho)**gammas[1] - (k[1] * (rho_ts[0])**gammas[1] - k[0] * rho_ts[0]**gammas[0])

    if rho_ts[1] < rho:
        pres = k[2] * (rho)**gammas[2] - (k[2] * (rho_ts[1])**gammas[2] - (k[1] * (rho_ts[1])**gammas[1] - (k[1] * (rho_ts[0])**gammas[1] - k[0] * rho_ts[0]**gammas[0])))

    return pres
