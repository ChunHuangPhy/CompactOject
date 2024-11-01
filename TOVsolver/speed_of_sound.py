from TOVsolver.constant import c, G
from TOVsolver.unit import g_cm_3, dyn_cm_2


def speed_of_sound_calc(density, pressure):
    """Function that calculates the speed of sound by taking the gradient of the euqation of state.
    Args:
        density (array): numpy 1Darray.
        pressure (array): numpy 1Darray.

    Returns:
        speed_of_sound (array): numpy 1Darray.
    """

    # tzzhou: migrating to new units
    density = density / c**2 * G / g_cm_3
    pressure = pressure / c**4 * G / dyn_cm_2

    drho = density[1:] - density[:-1]
    dp = pressure[1:] - pressure[:-1]
    speed_of_sound = dp / drho

    mask = density[:-1] > 1.5e-14
    d2 = density[:-1][mask]
    C_s = speed_of_sound[mask]

    return C_s, d2
