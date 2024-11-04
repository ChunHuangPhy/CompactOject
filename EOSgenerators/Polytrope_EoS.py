import numpy as np
from TOVsolver.unit import g_cm_3, dyn_cm_2

def polytrope(rhos, theta,
              rho_t_start=4.3721E11*g_cm_3, P_t_start=7.7582E29*dyn_cm_2):
    """Generate pressure from given densities and given parameters,

    Args:
        rhos (array): An array that consists of the initial values of sigma, omega, rho, and chemical
        potential obtained from the initial_values function.
        theta (array): An array representing the parameters used to determine a polytrope EoS.
        In this case, theta is composed of gammas and densities of transition in order.
        rho_t_start (fload): Start point of the density for polytrope EoS
        P_t_start (fload): Start point of the pressure for polytrope EoS

    Returns:
        pres (array): EOS ingredient, pressure in dyn/cm2

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
