import TOVsolver.constant as constant
import TOVsolver.solver_code as TOV_solver
import TOVsolver.EoS_import as EoS_import
import EOSgenerators.crust_EOS as crust
import EOSgenerators.RMF_EOS as RMF
import TOVsolver.main as main

import numpy as np
import math

oneoverfm_MeV = constant.oneoverfm_MeV
c = constant.c
G = constant.G

def MRlikihood_kernel(eps_crust,pres_crust,x,theta):
    kernel = x
    g_sigma, g_omega,g_rho, kappa, lambda_0, zeta, Lambda_w, d1 = theta
    m_sig = 495 / oneoverfm_MeV
    m_w = 3.96544
    m_rho = 3.86662
    theta = np.array([m_sig, m_w, m_rho, g_sigma, g_omega, g_rho, kappa, lambda_0, zeta, Lambda_w])
    ep, pr = RMF.compute_EOS(eps_crust, pres_crust, theta)

    eps_total = np.hstack((eps_crust,ep))
    pres_total = np.hstack((pres_crust,pr))
    
    if d1 ==0 :
        likelihood = -1e101
    else:
        if   all(x<y for x,y in zip(eps_total[:], eps_total[1:])) and all(x<y for x, y in zip(pres_total[:], pres_total[1:])):
            MR = main.OutputMRpoint(d1,eps_total,eps_total).T
        if len(MR[0]) == False:
            likelihood = -1e101
        else:
            likelihood = np.log(kernel.evaluate((MR[0], MR[1])))
    if likelihood <= -1e101:
        return -1e101
    else:
        return likelihood

def MRlikihood_Gaussian(eps_crust,pres_crust,x,theta):
    Mvalue, Rvalue, Mwidth,Rwidth = x
    g_sigma, g_omega,g_rho, kappa, lambda_0, zeta, Lambda_w, d1 = theta
    m_sig = 495 / oneoverfm_MeV
    m_w = 3.96544
    m_rho = 3.86662
    the = np.array([m_sig, m_w, m_rho, g_sigma, g_omega, g_rho, kappa, lambda_0, zeta, Lambda_w])
    ep, pr = RMF.compute_EOS(eps_crust, pres_crust, the)

    eps_total = np.hstack((eps_crust,ep))
    pres_total = np.hstack((pres_crust,pr))
    
    sigma_x = Rwidth
    sigma_y = Mwidth
    
    if d1 ==0 :
        likelihood = -1e101
    else:
        if   all(x<y for x,y in zip(eps_total[:], eps_total[1:])) and all(x<y for x, y in zip(pres_total[:], pres_total[1:])):
            MR = main.OutputMRpoint(d1,eps_total,eps_total).T
        if len(MR[0]) == False:
            likelihood = -1e101
        else:
            fx = 1/(sigma_x*sigma_y*(np.sqrt(2*np.pi))**2)*np.exp(-np.power(MR[0][0]-Rvalue, 2.)/(2*np.power(sigma_x,2.))-np.power(MR[1][0]-Mvalue, 2.)/(2*np.power(sigma_y,2.)))
            likelihood = np.log(fx)
    if likelihood <= -1e101:
        return -1e101
    else:
        return likelihood
    
def Masslikihood_Gaussian(eps_crust,pres_crust,x,theta):
    Mvalue, Mwidth = x
    g_sigma, g_omega,g_rho, kappa, lambda_0, zeta, Lambda_w, d1 = theta
    m_sig = 495 / oneoverfm_MeV
    m_w = 3.96544
    m_rho = 3.86662
    theta = np.array([m_sig, m_w, m_rho, g_sigma, g_omega, g_rho, kappa, lambda_0, zeta, Lambda_w])
    ep, pr = RMF.compute_EOS(eps_crust, pres_crust, theta)

    eps_total = np.hstack((eps_crust,ep))
    pres_total = np.hstack((pres_crust,pr))
    
    sigma_y = Mwidth
    
    if d1 ==0 :
        likelihood = -1e101
    else:
        if   all(x<y for x,y in zip(eps_total[:], eps_total[1:])) and all(x<y for x, y in zip(pres_total[:], pres_total[1:])):
            MR = main.OutputMRpoint(d1,eps_total,eps_total).T
        if len(MR[0]) == False:
            likelihood = -1e101
        else:
            if MR[0][1]>= Mvalue:
                likelihood = 0
            else:    
                fx = 1/(sigma_y*(np.sqrt(2*np.pi))**2)*np.exp(-np.power(MR[1][0]-Mvalue, 2.)/(2*np.power(sigma_y,2.)))
                likelihood = np.log(fx)
    if likelihood <= -1e101:
        return -1e101
    else:
        return likelihood


def Kliklihood(theta,K_low,K_up):
    g_sigma, g_omega,g_rho, kappa, lambda_0, zeta, lambda_w, d1 = theta
    
    m_sig = 495 / oneoverfm_MeV
    m_w = 3.96544
    m_rho = 3.86662
    
    m_b = 939
    E_per =-16.28
    m_sig =495/oneoverfm_MeV
    m_eff = (0.55+np.random.random_sample()*0.1)*m_b
    kappa_one = kappa*oneoverfm_MeV
    m_sig_one = m_sig*oneoverfm_MeV
    m_w_one = m_w*oneoverfm_MeV
    rho = 0.1505* (197.33**3 )
    k_f = (0.5*3.*rho*math.pi**2)**(1/3.)
    E_f = (k_f**2 + m_eff**2)**(1/2)

    g_w_omega = E_per - E_f + m_b
    sigma = m_b - m_eff
    m_sig_star_sqr_g_sig_sqr = (m_sig_one**2/g_sigma**2)+ (kappa_one*sigma) + 0.5*lambda_0*sigma**2
    m_w_star_sqr = m_w_one**2 + 0.5*zeta*(g_omega**2)*(g_w_omega**2)
    # Components o f rho_s_prime
    term1 = k_f/E_f
    term2 = E_f**2 + (2*m_eff**2)
    term3 = 3*m_eff**2
    term4 = np.log((k_f+E_f)/m_eff)
    rho_s_prime = (1/math.pi**2)*(term1*term2 - term3*term4)
    dM_drho = -(m_eff/E_f)*(1./(m_sig_star_sqr_g_sig_sqr+rho_s_prime))
    M_Ef = m_eff/E_f
    dEf_drho = math.pi**2/(2*k_f*E_f)
    dW0_drho = g_omega**2/m_w_star_sqr
    
    K = 9.*rho*(dEf_drho+dW0_drho+(M_Ef*dM_drho))
    center = (K_low + K_up)/2
    width = (K_up - K_low)/2
    p_K = -0.5*abs( center - K)**10./width**10.
    return p_K

def Jliklihood(theta,J_low,J_up):
    g_sigma, g_omega,g_rho, kappa, lambda_0, zeta, Lambda_w, d1 = theta
    
    m_sig = 495 / oneoverfm_MeV
    m_w = 3.96544
    m_rho = 3.86662
    
    m_b = 939
    E_per =-16.28
    m_sig =495/oneoverfm_MeV
    m_eff = (0.55+np.random.random_sample()*0.1)*m_b
    kappa_one = kappa*oneoverfm_MeV
    m_sig_one = m_sig*oneoverfm_MeV
    m_w_one = m_w*oneoverfm_MeV
    rho = 0.1505* (197.33**3 )
    k_f = (0.5*3.*rho*math.pi**2)**(1/3.)
    E_f = (k_f**2 + m_eff**2)**(1/2)
    g_w_omega = E_per - E_f + m_b
    J_0 = k_f**2/(6*E_f)
    J_1 = (1/8.)*(rho*g_rho**2)/(m_rho**2+(2*Lambda_w*(g_w_omega**2)*g_rho**2))
    J = J_0 + J_1
    
    center = (J_low + J_up)/2
    width = (J_up - J_low)/2
    p_J = -0.5*abs( center - J)**10./width**10.
    return p_J

def Lliklihood(theta,L_low,L_up):
    g_sigma, g_omega,g_rho, kappa, lambda_0, zeta, Lambda_w, d1 = theta
    
    m_sig = 495 / oneoverfm_MeV
    m_w = 3.96544
    m_rho = 3.86662
    
    m_b = 939
    E_per =-16.28
    m_sig =495/oneoverfm_MeV
    m_eff = (0.55+np.random.random_sample()*0.1)*m_b
    kappa_one = kappa*oneoverfm_MeV
    m_sig_one = m_sig*oneoverfm_MeV
    m_w_one = m_w*oneoverfm_MeV
    rho = 0.1505* (197.33**3 )
    k_f = (0.5*3.*rho*math.pi**2)**(1/3.)
    E_f = (k_f**2 + m_eff**2)**(1/2)
    J_0 = k_f**2/(6*E_f)
    J_1 = (1/8.)*(rho*g_rho**2)/(m_rho**2+(2*Lambda_w*(g_w_omega**2)*g_rho**2))
    J= J_0 + J_1 
    g_w_omega = E_per - E_f + m_b
    m_star_sqr_E_f_sqr = m_eff**2/E_f**2 
    rho_m_star = 3*rho/m_eff 
    sigma = m_b - m_eff
    m_sig_star_sqr_g_sig_sqr = (m_sig**2/g_sigma**2)+\
    (kappa*sigma) + 0.5*lambda_0*sigma**2
    
    # Components o f rho_s_prime
    term1 = k_f/E_f
    term2 = E_f**2 + (2*m_eff**2)
    term3 = 3*m_eff**2
    term4 = np.log((k_f+E_f)/m_eff)
    rho_s_prime = (1/math.pi**2)*(term1*term2 - term3*term4)
    dM_drho = -(m_eff/E_f)*(1./(m_sig_star_sqr_g_sig_sqr+rho_s_prime))
    m_w_star_sqr = m_w**2 + 0.5*zeta*(g_omega**2)*(g_w_omega**2)
    L_0 = J_0*(1+(m_star_sqr_E_f_sqr*(1. -(rho_m_star*dM_drho))))
    L_1 = 3.*J_1*(1.- 32.*(g_omega**2/m_w_star_sqr)*g_w_omega*Lambda_w*J_1)
    L= L_0 + L_1
    
    center = (L_low + L_up)/2
    width = (L_up - L_low)/2
    p_L = -0.5*abs( center - J)**10./width**10.
    return p_L
