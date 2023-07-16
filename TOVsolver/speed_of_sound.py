import numpy as np
from TOVsolver.constant import c,G

def speed_of_sound_calc(density, pressure):


    """Function that calculates the speed of sound by taking the gradient of the euqation of state.
    Args:
        density (array): numpy 1Darray. 
        pressure (array): numpy 1Darray.

    Returns:
        speed_of_sound (array): numpy 1Darray.
    """

    speed_of_sound = []
    #density = density*c**2/G
    #pressure = pressure*c**4/G
    
    for i in range(0,len(density)-1):
        speed_of_sound.append((pressure[i+1]-pressure[i])/(density[i+1]-density[i]))
    d2 = []
    C_s= []
    #eps2 = []
    for i in range(0,len(speed_of_sound)):
        if density[i]> 1.5e-14:
            d2.append(density[i])
            C_s.append(speed_of_sound[i])
    return C_s,d2
