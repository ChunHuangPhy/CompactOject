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
pi = np.pi

RFSU2R = []
MFSU2R = []
density = numpy.logspace(14.3, 15.6, 50)
if   all(x<y for x, y in zip(eps_total_poly[j][:], eps_total_poly[j][1:])) and all(x<y for x, y in zip(pres_total_poly[j][:], pres_total_poly[j][1:])):
    for i in range(len(density)):
        try:
            RFSU2R[j].append(solveTOV(density[i], eps_total_poly[j], pres_total_poly[j])[1])
            MFSU2R[j].append(solveTOV(density[i], eps_total_poly[j], pres_total_poly[j])[0])
        except OverflowError as e:
            break
        if i > 20 and solveTOV(density[i], eps_total_poly[j], pres_total_poly[j])[0] - solveTOV(density[i-1], eps_total_poly[j], pres_total_poly[j])[0]< 0:
            break
    if len(RFSU2R[j]) == False:
        break
    else:
        del RFSU2R[j][-1]
        del MFSU2R[j][-1]