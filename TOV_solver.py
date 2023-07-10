def TOV(r, y, inveos):
    pres, m = y
    
    #eps = 10**inveos(numpy.log10(pres))
    eps = inveos(pres)
    dpdr = -(eps + pres) * (m + 4.*pi*r**3. * pres)
    dpdr = dpdr/(r*(r - 2.*m))
    dmdr = 4.*pi*r**2.0 * eps
    
    return numpy.array([dpdr, dmdr])

def solveTOV(rhocent, eps, pres):
    

    #eos = UnivariateSpline(numpy.log10(eps), numpy.log10(pres), k=1, s=0)
    #inveos = UnivariateSpline(numpy.log10(pres), numpy.log10(eps), k=1, s=0)
    #We could change this to Double Log Interpolation。
    eos = UnivariateSpline(eps, pres, k=3, s=0)
    inveos = UnivariateSpline(pres, eps, k=3, s=0)
    Pmin = pres[20]
    
    r = 4.441e-16
    dr = 10.
    
    rhocent = rhocent * G/c**2.
    #pcent = 10**eos(numpy.log10(rhocent))
    pcent = eos(rhocent)
    P0 = pcent - (2.*pi/3.)*(pcent + rhocent) *(3.*pcent + rhocent)*r**2.
    m0 = 4./3. *pi *rhocent*r**3.
    stateTOV = numpy.array([P0, m0])
    
    sy = ode(TOV, None).set_integrator('dopri5')
    
  
    #have been modified from Irida to this integrator
    sy.set_initial_value(stateTOV , r).set_f_params(inveos)
    
    while sy.successful() and stateTOV[0]>Pmin:
        stateTOV = sy.integrate(sy.t+dr)
        dpdr, dmdr = TOV(sy.t+dr, stateTOV, inveos)
        dr = 0.46 * (1./stateTOV[1] * dmdr - 1./stateTOV[0]*dpdr)**(-1.)

    return stateTOV[1]*c**2./G/Msun, sy.t/1e5  
