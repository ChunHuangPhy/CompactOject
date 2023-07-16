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
    density = density*G/c**2
    pressure = pressure*G/c**4
    
    for i in range(0,len(density)-1):
        speed_of_sound.append((density[i+1]-density[i])/(pressure[i+1]-pressure[i]))
    p2 = []
    C_s= []
    #eps2 = []
    for i in range(0,len(speed_of_sound)):
        if pressure[i]> 2e-14:
            p2.append(pressure[i])
            C_s.append(speed_of_sound[i])
    return C_s
