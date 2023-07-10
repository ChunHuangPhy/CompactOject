# Python packages
import numpy as np
import math
from scipy.interpolate import UnivariateSpline
from scipy.constants import pi
from scipy.integrate import odeint, ode
from matplotlib import pyplot
from scipy import optimize
from itertools import repeat

# Import files
import TOV_solver

# Global Variables
dyncm2_to_MeVfm3 = 1./(1.6022e33)
gcm3_to_MeVfm3 = 1./(1.7827e12)
oneoverfm_MeV = 197.33

energy_density, pressure = EOS_import()
RFSU2R = []
MFSU2R = []
density = numpy.logspace(14.3, 15.6, 50)
#This following step is to make a dicision whether the EOS ingredients is always increase. We can do that outsie of this main to the 
#EOS import.
#if   all(x<y for x, y in zip(eps_total_poly[:], eps_total_poly[[1:])) and all(x<y for x, y in zip(pres_total_poly[j][:], pres_total_poly[j][1:])):
for i in range(len(density)):
    try:
        RFSU2R = TOV_solver.solveTOV(density[i], energy_density, pressure)[1]
        MFSU2R = TOV_solver.solveTOV(density[i], energy_density, pressure)[0]
    #This is sentense is for avoiding the outflow of the result, like when solveTOV blow up because of ill EOS, we need to stop
    except OverflowError as e:
        break