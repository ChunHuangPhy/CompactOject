import TOVsolver.constant as constant
from TOVsolver.EoS_import import EOS_import
import numpy as np

def PolyInterpolate(eps_crust, pressure_crust):
    """Polytrope connecting crust part of equation of state to core part
    
    Args:
        eps_crust (array): the energy density of crust EoS in MeV/fm3, times a G/c**2 factor
        pres_crust (array): the pressure from crust EoS model in MeV/fm3, times a G/c**4 factor
        
    Returns:
        eps_combine (float): EOS ingredient, combined crust and inter-crust part  energy density in
        MeV/fm3, times a G/c**2 factor
        pres_combine (float): EOS ingredient, combined crust and inter-crust part pressure in MeV/fm3,
        times a G/c**2 factor
        
    """
    c = constant.c
    G = constant.G
    
    xs_polytro = np.logspace(11.7, 14.28, num=1000, base=10)
    ys_polytropic = (4.75764e29 + 9.04238e13 * xs_polytro**(4/3)) * G / c**4
    xs_polytropic = xs_polytro * G / c**2

    eps_combine = np.append(eps_crust, xs_polytropic)
    pres_combine = np.append(pressure_crust, ys_polytropic)
    return eps_combine, pres_combine