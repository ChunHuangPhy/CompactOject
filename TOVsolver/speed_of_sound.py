import numpy as np

def speed_of_sound_calc(density, pressure):

    speed_of_sound = np.gradient(density,pressure)
    return speed_of_sound
