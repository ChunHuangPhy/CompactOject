def m1_from_mc_m2(mc, m2):
    m2 = numpy.array(m2)
    num1 = (2. / 3.)**(1. / 3.) * mc**5.
    denom1 = ((9 * m2**7. * mc**5. + numpy.sqrt(3.) * 
              numpy.sqrt(abs(27 * m2**14. * mc**10. - 
                         4. * m2**9. * mc**15.)))**(1. / 3.))
    denom2 = 2.**(1. / 3.) * 3.**(2. / 3.) * m2**3.
    return num1 / denom1 + denom1 / denom2

def TOV(r, y,inveos):
    pres, m,h,b = y
    
    #eps = 10**inveos(numpy.log10(pres))
    eps = inveos(pres)
    dpdr = -(eps + pres) * (m + 4.*pi*r**3. * pres)
    dpdr = dpdr/(r*(r - 2.*m))
    dmdr = 4.*pi*r**2.0 * eps
    dhdr = b
    dfdr = 2. * numpy.power(1. - 2. * m / r, -1) * h * \
        (-2. * numpy.pi * (5. * eps + 9. * pres + (eps + pres)**2. /
                        (pres )) + 3. /numpy.power(r,2) + 2. *
            numpy.power(1. - 2. * m / r,-1) * numpy.power(m / numpy.power(r,2) +
         4. * numpy.pi * r * pres,2)) \
        + 2. * b / r * numpy.power(1. - 2. * m / r, -1) * \
        (-1. + m / r + 2. * numpy.pi * numpy.power(r,2) * (eps - pres))

    return numpy.array([dpdr, dmdr, dhdr, dfdr])

def tidal_deformability(y2, Mns, Rns):

    C = Mns / Rns
    Eps = 4. * C**3. * (13. - 11. * y2 + C * (3. * y2 - 2.) +
                        2. * C**2. * (1. + y2)) + \
        3. * (1. - 2. * C)**2. * (2. - y2 + 2. * C * (y2 - 1.)) * \
        numpy.log(1. - 2. * C) + 2. * C * (6. - 3. * y2 + 3. * C * (5. * y2 - 8.))
    tidal_def = 16. / (15. * Eps) * (1. - 2. * C)**2. *\
        (2. + 2. * C * (y2 - 1.) - y2)

    return tidal_def

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
    h0 = r**2.
    b0 = 2. * r
    stateTOV = numpy.array([P0, m0, h0,b0])
    sy = ode(TOV, None).set_integrator("dopri5")
    
    #have been modified from Irida to this integrator
    sy.set_initial_value(stateTOV , r).set_f_params(inveos)
    
    while sy.successful() and stateTOV[0]>Pmin:
        stateTOV = sy.integrate(sy.t+dr)
        dpdr, dmdr = TOV(sy.t+dr, stateTOV, inveos)
        dr = 0.46 * (1./stateTOV[1] * dmdr - 1./stateTOV[0]*dpdr)**(-1.)
    Mb = stateTOV[1]
    Rns = sy.t
    y = Rns * stateTOV[3] /stateTOV[2] 
    tidal = tidal_deformability(y, Mb, Rns)
    
    return Mb*c**2./G/Msun, Rns/1e5,tidal
