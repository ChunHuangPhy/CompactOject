def TOV(r, y, inveos):
    pressure, mass = y
    #eps = 10**inveos(numpy.log10(pres))
    energy_density = inveos(pressure)
    dpdr = -(energy_density + pressure) * (mass + 4.*pi*r**3. * pressure)
    dpdr = dpdr/(r*(r - 2.*mass))
    dmdr = 4.*pi*r**2.0 * energy_density
    return numpy.array([dpdr, dmdr])

# Function solves the TOV equation, returning mass and radius
def solveTOV(center_rho, energy_density, pressure):
    #eos = UnivariateSpline(numpy.log10(eps), numpy.log10(pres), k=1, s=0)
    #inveos = UnivariateSpline(numpy.log10(pres), numpy.log10(eps), k=1, s=0)
    #We could change this to Double Log Interpolation。

    eos = UnivariateSpline(energy_density, pressure, k=3, s=0)
    inveos = UnivariateSpline(pressure, energy_density, k=3, s=0)
    Pmin = pres[20]
    r = 4.441e-16
    dr = 10.
    center_rho = center_rho * G/c**2.
    #pcent = 10**eos(numpy.log10(rhocent))
    pcent = eos(center_rho)

    P0 = pcent - (2.*pi/3.)*(pcent + center_rho) *(3.*pcent + center_rho)*r**2.
    m0 = 4./3. *pi *center_rho*r**3.
    stateTOV = numpy.array([P0, m0])
    sy = ode(TOV, None).set_integrator("dopri5")
    
    #have been modified from Irida to this integrator
    sy.set_initial_value(stateTOV , r).set_f_params(inveos)

    while sy.successful() and stateTOV[0]>Pmin:
        stateTOV = sy.integrate(sy.t+dr)
        dpdr, dmdr = TOV(sy.t+dr, stateTOV, inveos)
        dr = 0.46 * (1./stateTOV[1] * dmdr - 1./stateTOV[0]*dpdr)**(-1.)
    radius = stateTOV[1]*c**2./G/Msun #units of km
    mass = sy.t/1e5 #units of solar mass
    return radius, mass
