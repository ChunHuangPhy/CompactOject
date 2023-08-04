import numpy as np
import math
from matplotlib import pyplot
from NumbaMinpack import lmdif, hybrd, minpack_sig
from numba import njit, jit, cfunc

c = 3e10
G = 6.67428e-8
Msun = 1.989e33

dyncm2_to_MeVfm3 = 1./(1.6022e33)
gcm3_to_MeVfm3 = 1./(1.7827e12)
oneoverfm_MeV = 197.33

m_e = 2.5896 * 10**-3
m_mu = 0.53544
m_n = 4.7583690772
m_p = 4.7583690772

J_B = 1/2.
b_B = 1

m_l = [m_e, m_mu]
m_b = [m_p, m_n]

Matrix_b = np.array([[1., 1., 1/2., 1., 1., 1.],[1., 0., -1/2., 1., 1., 1.]])

Matrix_l = np.array([[0., -1., 1/2.],[0., -1., 1/2.]])

# defining the function to optimize with numba-minpack

@cfunc(minpack_sig)
def myfunc(x, fvec, args):
    """iterate the the sigma, omega, rho term and chemical potential of electron and neutron at 
    any given density(fast version)
    
    Args:
        x (array): initial sigma omega rho and chemical potential from initial_values function
        fvec (array): iterated sigma sigma omega rho and chemical potential for output.
        args (array): paramters of determine a RMF model in lagrangian, here we have 10 parameters.

    Returns:
        sigma (float): sigma term in lagrangian
        omega (float): omega term in lagrangian
        rho_03 (float): rho term in lagrangian
        mu_n (float): chemical potential of neutron matter
        mu_e (float): chemical potential of electron portion
        
    """
    m_sig = args[0]
    m_w = args[1]
    m_rho = args[2]
    g_sigma = args[3]
    g_omega = args[4]
    g_rho = args[5]
    kappa = args[6]
    lambda_0 = args[7]
    zeta = args[8]
    Lambda_w = args[9]
    rho = args[10]

    m_e = 2.5896 * 10**-3
    m_mu = 0.53544
    m_n = 4.7583690772
    m_p = 4.7583690772

    J_B = 1/2.
    b_B = 1

    m_l = np.array([m_e, m_mu])
    m_b = np.array([m_p, m_n])

    Matrix_b = np.array([[1., 1., 1/2., 1., 1., 1.],[1., 0., -1/2., 1., 1., 1.]])

    Matrix_l = np.array([[0., -1., 1/2.],[0., -1., 1/2.]])
    
    sigma = x[0]
    omega=x[1]
    rho_03 = x[2]
    mu_n = x[3]
    mu_e = x[4]

    rho_B_list = []
    rho_SB_list = []
    q_list = []
    
    m_eff = m_n - (g_sigma*sigma)

    for i in range(len(Matrix_b)):    
        mu_b = Matrix_b[i,0] * mu_n - Matrix_b[i, 1] * mu_e
        
        E_fb = mu_b - g_omega * omega - g_rho * rho_03 * Matrix_b[i,2]
        
        k_fb_sq = E_fb**2 - m_eff**2
        if k_fb_sq < 0:
            k_fb_sq = 0.0
            E_fb = m_eff
        
        k_fb = math.sqrt(k_fb_sq)
        
        rho_B = ((2*J_B) + 1)*b_B*k_fb**3 / (6.*math.pi**2)
        rho_SB = (m_eff/(2.*math.pi**2))*(E_fb*k_fb - \
        (m_eff**(2))*np.log((E_fb + k_fb)/m_eff))
        
        rho_B_list.append(rho_B)
        rho_SB_list.append(rho_SB)
        
        Q_B = ((2.*J_B) + 1.)*Matrix_b[i,1]*k_fb**3 / (6.*math.pi**2)
        q_list.append(Q_B)
        
    for j in range(len(Matrix_l)):        
        mu_l = Matrix_l[j,0]*mu_n - Matrix_l[i,1]*mu_e
        E_fl = mu_l
        
        k_fl_sq = E_fl**2 - m_l[j]**2
        
        if k_fl_sq < 0.0:
            k_fl_sq = 0.0
    
        k_fl = math.sqrt(k_fl_sq)
        
        Q_L = ((2.*J_B) + 1.)*Matrix_l[i,1]*(k_fl**3) / (6.*(math.pi**2))
        q_list.append(Q_L)
        
    sum1 = 0.0
    sum2 = 0.0
    sum3 = 0.0
    sum4 = 0.0
    sum5 = 0.0
    for i in range(len(rho_SB_list)):
        sum1 = sum1 + rho_SB_list[i] * Matrix_b[:,3][i]
    for i in range(len(rho_B_list)):
        sum2 = sum2 + rho_B_list[i] * Matrix_b[:,4][i]
    for i in range(len(rho_B_list)):
        sum3 = sum3 + rho_B_list[i] * Matrix_b[:,5][i] * Matrix_b[:,2][i]
    for i in range(len(rho_B_list)):
        sum4 = sum4 + rho_B_list[i]
    for i in range(len(q_list)):
        sum5 = sum5 + q_list[i]

    fvec[0] = (sigma*(m_sig**2)/g_sigma - sum1 +(kappa*(g_sigma*sigma)**2)/2.\
              + (lambda_0*(g_sigma*sigma)**3)/6.)**2
    fvec[1] = (omega*(m_w**2)/g_omega - sum2 + (zeta*(g_omega*omega)**3)/6.\
              + 2.*Lambda_w*g_omega*omega*(rho_03*g_rho)**2)**2
    fvec[2] = (rho_03*(m_rho**2)/g_rho - sum3 + 2.*Lambda_w*g_rho*rho_03*(omega*g_omega)**2)**2
    
    fvec[3] = (rho - sum4)**2
    fvec[4] = sum5**2

    
funcptr = myfunc.address # address in memory to myfunc

# defining the initial values and E-P functions with numba
@jit
def initial_values(rho, theta):
    """Outputs the the sigma, omega, rho term and chemical potential of electron and neutron at 
    given initial density (fast version)
    
    Args:
        rho (float): given nuclear density
        theta (array): paramters of determine a RMF model in lagrangian, here we have 10 parameters.

    Returns:
        sigma (float): sigma term in lagrangian
        omega (float): omega term in lagrangian
        rho_03 (float): rho term in lagrangian
        mu_n (float): chemical potential of neutron matter
        mu_e (float): chemical potential of electron portion
    """
    m_sig, m_w, m_rho, g_sigma, g_omega, g_rho, kappa, lambda_0, zeta, Lambda_w = theta

    m_e = 2.5896 * 10**-3
    m_mu = 0.53544
    m_n = 4.7583690772
    m_p = 4.7583690772
    
    rho_0 = 0.1505

    sigma = g_sigma*rho/(m_sig**2)
    rho_03 = -g_rho*rho/(2.*(m_rho**2))
    omega = rho*((((m_w**2)/g_omega)+\
    (2.*Lambda_w*((g_rho*rho_03)**2)*g_omega))**(-1.))
    m_eff = m_n-(g_sigma*sigma)
    mu_n = m_eff + g_omega*omega + g_rho*rho_03*Matrix_b[1, 2]
    mu_e = 0.12*m_e*(rho/rho_0)**(2/3.)
    
    return sigma, omega, rho_03, mu_n, mu_e

@njit
def Energy_density_Pressure(x, rho, theta):
    """Generate pressure and energy density two EOS ingredient from given RMF term and given parameters,
    (fast version)
    
    Args:
        x (array): An array that consists of the initial values of sigma, omega, rho, and chemical 
        potential obtained from the initial_values function.
        rho (float): The central density from which the computation of the equation of state begins.
        theta (array): An array representing the parameters used to determine a RMF model in the 
        Lagrangian. In this case, the RMF model is defined by 10 parameters.


    Returns:
        energy_density (float): EOS ingredient, energy density in g/cm3
        pressure (float): EOS ingredient, pressure in dyn/cm3
    """
    sigma, omega, rho_03, mu_n, mu_e = x
    
    m_sig, m_w, m_rho, g_sigma, g_omega, g_rho, kappa, lambda_0, zeta, Lambda_w = theta

    m_e = 2.5896 * 10**-3
    m_mu = 0.53544
    m_n = 4.7583690772
    m_p = 4.7583690772

    J_B = 1/2.
    b_B = 1

    m_l = np.array([m_e, m_mu])
    m_b = np.array([m_p, m_n])
    
    energy_b = 0
    energy_l = 0
    multi = 0
    
    m_eff = m_n - (g_sigma*sigma)

    for i in range(len(Matrix_b)):        
        mu_b = Matrix_b[i,0]*mu_n - Matrix_b[i, 1]*mu_e
        
        E_fb = mu_b - g_omega*omega - g_rho*rho_03*Matrix_b[i,2]
        
        k_fb_sq = E_fb**2 - m_eff**2
        if k_fb_sq < 0:
            k_fb_sq = 0.0
            E_fb = m_eff
          
        k_fb = math.sqrt(k_fb_sq)

        rho_B = ((2.*J_B) + 1.)*b_B*(k_fb**3) / (6.*math.pi**2)
        
        multi = multi + mu_b*rho_B
        energy_baryon = (1/(8.*(math.pi**2)))*(k_fb*(E_fb**3)\
        + (k_fb**3)*E_fb - np.log((k_fb + E_fb)/m_eff)*m_eff**4)

        energy_b = energy_b + energy_baryon
        
    for j in range(len(Matrix_l)):
        
        mu_l = Matrix_l[i, 0]*mu_n - Matrix_l[j, 1]*mu_e
        
        k_fl_sq = mu_l**2 - m_l[j]**2
        if k_fl_sq < 0.0:
            k_fl_sq = 0.0
        k_fl = math.sqrt(k_fl_sq)
        
        rho_l = k_fl**3 / (3.*math.pi**2)
        
        multi = multi + mu_l*rho_l
        energy_lepton = (1/(8.*(math.pi**2)))*(k_fl*(mu_l**3)\
        +mu_l*(k_fl**3)-(m_l[j]**4)*np.log((k_fl+mu_l)/m_l[j]))
        
        energy_l = energy_l + energy_lepton
        
    sigma_terms =  0.5*((sigma*m_sig)**2) + (kappa*((g_sigma*sigma)**3))/6.\
                    + (lambda_0*((g_sigma*sigma)**4))/24.
        
    omega_terms = 0.5*((omega*m_w)**2) +(zeta*((g_omega*omega)**4))/8.
        
    rho_terms = 0.5*((rho_03*m_rho)**2)+ + 3.*Lambda_w*((g_rho*rho_03*g_omega*omega)**2)
    
    energy_density = energy_b + energy_l + sigma_terms \
    + omega_terms + rho_terms
        
    Pressure = multi - energy_density

    return energy_density, Pressure

# define a function that computes the EoS 

@njit
def compute_EOS(eps_crust, pres_crust, theta):
    """Generate core part equation of state, main function, from RMF model (fast version)
    
    Args:
        eps_crust (array): the energy density of crust EoS in MeV/fm3, times a G/c**2 factor
        pres_crust (array): the pressure from crust EoS model in MeV/fm3, times a G/c**4 factor
        theta (array): An array representing the parameters used to determine a RMF model in the 
        Lagrangian. In this case, the RMF model is defined by 10 parameters.

    Returns:
        energy_density (float): EOS ingredient, energy density in g/cm3
        pressure (float): EOS ingredient, pressure in dyn/cm3
        
    """
    dt = 0.05
    rho_0 = 0.1505
    
    x_init = np.array(initial_values(0.1*rho_0, theta))
    Energy = []
    Pressure = []
    for i in range(1, 125):

        rho = i*dt*rho_0
        args = np.append(theta, rho)
        xsol, fvec, success, info = lmdif(funcptr, x_init, 5, args)
        Re = Energy_density_Pressure(x_init, rho, theta)

        Energy.append(Re[0]*oneoverfm_MeV/gcm3_to_MeVfm3)
        Pressure.append(Re[1]*oneoverfm_MeV/dyncm2_to_MeVfm3)

        x_init = xsol
    Energy = np.array(Energy)
    Pressure = np.array(Pressure)
    
    end = 0
    for i in range(0,len(Energy)-1):
        if Energy[i]*G/c**2 > max(eps_crust) and i >18:
            end = i+2
            break
        end += 1
    ep = Energy[end::] * G/c**2
    pr = Pressure[end::] * G/c**4
    
    return ep, pr

