import numpy as np

def speed_of_sound_calc(density, pressure):


    """Function that calculates the speed of sound by taking the gradient of the euqation of state.
    Args:
        density (array): numpy 1Darray. 
        pressure (array): numpy 1Darray.

    Returns:
        speed_of_sound (array): numpy 1Darray.
    """

    speed_of_sound = np.gradient(density,pressure)
    return speed_of_sound
