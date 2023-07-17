import TOVsolver.constant as constant
from TOVsolver.EoS_import import EOS_import
import numpy as np

def PolyInterpolate(eps_crust, pressure_crust):
    
    c = constant.c
    G = constant.G
    
    xs_polytro = np.logspace(11.7, 14.28, num=1000, base=10)
    ys_polytropic = (4.75764e29 + 9.04238e13 * xs_polytro**(4/3)) * G / c**4
    xs_polytropic = xs_polytro * G / c**2

    eps_combine = np.append(eps_crust, xs_polytropic)
    pres_combine = np.append(pressure_crust, ys_polytropic)
    return eps_combine, pres_combine