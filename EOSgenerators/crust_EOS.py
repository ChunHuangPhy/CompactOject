import numpy as np
from TOVsolver.unit import g_cm_3, dyn_cm_2


def PolyInterpolate(eps_crust, pressure_crust):
    """Polytrope connecting crust part of equation of state to core part

    Args:
        eps_crust (array): the energy density of crust EoS
        pres_crust (array): the pressure from crust EoS model

    Returns:
        eps_combine (float): EOS ingredient, combined crust and inter-crust part energy density
        pres_combine (float): EOS ingredient, combined crust and inter-crust part pressure
    """

    xs_polytro = np.logspace(11.7, 14.28, num=1000, base=10)
    ys_polytro = 4.75764e29 + 9.04238e13 * xs_polytro ** (4 / 3)

    xs_polytropic = xs_polytro * g_cm_3
    ys_polytropic = ys_polytro * dyn_cm_2

    eps_combine = np.append(eps_crust, xs_polytropic)
    pres_combine = np.append(pressure_crust, ys_polytropic)
    return eps_combine, pres_combine
