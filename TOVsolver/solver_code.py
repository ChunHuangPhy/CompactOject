# Python packages
import numpy as np
import math
from scipy.interpolate import UnivariateSpline
from scipy.constants import pi
from scipy.integrate import odeint, ode
from matplotlib import pyplot
from scipy import optimize
from itertools import repeat

def m1_from_mc_m2(mc, m2):
    m2 = np.array(m2)
    num1 = (2. / 3.)**(1. / 3.) * mc**5.
    denom1 = ((9 * m2**7. * mc**5. + np.sqrt(3.) * 
              np.sqrt(abs(27 * m2**14. * mc**10. - 
                         4. * m2**9. * mc**15.)))**(1. / 3.))
    denom2 = 2.**(1. / 3.) * 3.**(2. / 3.) * m2**3.
    return num1 / denom1 + denom1 / denom2

def TOV(r, y,inveos):
    """a function that packing the whole TOV equations set

    Args:
        r (float): raius as integrate varible
        y (psudo-varible): containing pressure, mass, h and b as intergarte varibles
        to solve out the TOV equation
        inveos: the invert of the eos, pressure and energy density relation to integrate
        and interpolate.

    Returns:
        Mass (array): The array that contains all the Stars' masses, in M_sun as a 
        Units.
        Radius (array): The array that contains all the Stars's radius, in km.
        Tidal Deformability (array): The array that contains correpsonding Tidal property, 
        These are dimension-less.
    """
    pres, m,h,b = y
    
    #eps = 10**inveos(np.log10(pres))
    eps = inveos(pres)
    dpdr = -(eps + pres) * (m + 4.*pi*r**3. * pres)
    dpdr = dpdr/(r*(r - 2.*m))
    dmdr = 4.*pi*r**2.0 * eps
    dhdr = b
    dfdr = 2. * np.power(1. - 2. * m / r, -1) * h * \
        (-2. * np.pi * (5. * eps + 9. * pres + (eps + pres)**2. /
                        (pres )) + 3. /np.power(r,2) + 2. *
            np.power(1. - 2. * m / r,-1) * np.power(m / np.power(r,2) +
         4. * np.pi * r * pres,2)) \
        + 2. * b / r * np.power(1. - 2. * m / r, -1) * \
        (-1. + m / r + 2. * np.pi * np.power(r,2) * (eps - pres))

    return np.array([dpdr, dmdr, dhdr, dfdr])

def tidal_deformability(y2, Mns, Rns):
    """Compute Tidal deformability from y2, neutron star mass and raius

    Args:
        y2 (array): midiate varrible that computing tidal
        Mns (array): neutron star mass in g/cm3
        Rns (array): neutron star radius in cm.

    Returns:
        tidal_def (array): neutron star tidal deformability with unit-less.
    """
    C = Mns / Rns
    Eps = 4. * C**3. * (13. - 11. * y2 + C * (3. * y2 - 2.) +
                        2. * C**2. * (1. + y2)) + \
        3. * (1. - 2. * C)**2. * (2. - y2 + 2. * C * (y2 - 1.)) * \
        np.log(1. - 2. * C) + 2. * C * (6. - 3. * y2 + 3. * C * (5. * y2 - 8.))
    tidal_def = 16. / (15. * Eps) * (1. - 2. * C)**2. *\
        (2. + 2. * C * (y2 - 1.) - y2)

    return tidal_def

# Function solves the TOV equation, returning mass and radius
def solveTOV(center_rho, energy_density, pressure):
    """Solve TOV equation from given Equation of state in the neutron star 
    core density range

    Args:
        center_rho(array): This is the energy density here is fixed in main
        that is np.logspace(14.3, 15.6, 50)
        energy_density (array): Desity array of the neutron star EoS, in MeV/fm^{-3}
        Please check the Test_EOS.csv, The conversion from the g/cm3 to here is 
        (g/cm3)*G/c**2, that should will convert it to here in example.
        
        pressure (array): Pressure array of neutron star EoS, also in nautral unit
        with MeV/fm^{-3}, still please check the Test_EOS.csv, the conversion is 
        (dyn/cm3)*G/c**4.

    Returns:
        Mass (array): The array that contains all the Stars' masses, in M_sun as a 
        Units.
        Radius (array): The array that contains all the Stars's radius, in km.
        Tidal Deformability (array): The array that contains correpsonding Tidal property, 
        These are dimension-less.
    """
    #eos = UnivariateSpline(np.log10(eps), np.log10(pres), k=1, s=0)
    #inveos = UnivariateSpline(np.log10(pres), np.log10(eps), k=1, s=0)
    #We could change this to Double Log Interpolation。


    c = 3e10
    G = 6.67428e-8
    Msun = 1.989e33

    eos = UnivariateSpline(energy_density, pressure, k=3, s=0)
    inveos = UnivariateSpline(pressure, energy_density, k=3, s=0)
    Pmin = pressure[20]
    r = 4.441e-16
    dr = 10.
    center_rho = center_rho * G/c**2.
    #pcent = 10**eos(np.log10(rhocent))
    pcent = eos(center_rho)
    
    P0 = pcent - (2.*pi/3.)*(pcent + center_rho) *(3.*pcent + center_rho)*r**2.
    m0 = 4./3. *pi *center_rho*r**3.
    h0 = r**2.
    b0 = 2. * r
    stateTOV = np.array([P0, m0, h0,b0])
    sy = ode(TOV, None).set_integrator("dopri5")
    
    #have been modified from Irida to this integrator
    sy.set_initial_value(stateTOV , r).set_f_params(inveos)
    
    while sy.successful() and stateTOV[0]>Pmin:
        stateTOV = sy.integrate(sy.t+dr)
        dpdr, dmdr, dhdr, dfdr = TOV(sy.t+dr, stateTOV, inveos)
        dr = 0.46 * (1./stateTOV[1] * dmdr - 1./stateTOV[0]*dpdr)**(-1.)
    Mb = stateTOV[1]
    Rns = sy.t
    y = Rns * stateTOV[3] /stateTOV[2] 
    tidal = tidal_deformability(y, Mb, Rns)
    
    return Mb*c**2./G/Msun, Rns/1e5,tidal

