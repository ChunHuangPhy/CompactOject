import numpy as np

def polytrope(rho, theta):
    gamma1, gamma2, gamma3, rho_t1, rho_t2 = theta
    c = 2.99792458E10 # cgs
    G = 6.6730831e-8 # cgs
    rho_ns = 267994004080000.03 #cgs
    rho_t = 4.3721E11*G/c**2
    P_t = 7.7582E29* G / c**4
    
    P_ts, k = np.zeros(3), np.zeros(3)
    P_ts[0] = P_t
    k[0] = P_t / ((rho_t / rho_ns)**gamma1)
    P_ts[1] = k[0] * rho_t1**gamma1
    k[1] = P_ts[1] / (rho_t1**gamma2)
    P_ts[2] = k[1] * rho_t2**gamma2
    k[2] = P_ts[2] / (rho_t2**gamma3)

    # Calculate the pressure for the entire input array `rho`
    pres = np.where(rho <= rho_t1, k[0] * rho**gamma1,
                    np.where(rho <= rho_t2, k[1] * rho**gamma2, k[2] * rho**gamma3))

    return pres
