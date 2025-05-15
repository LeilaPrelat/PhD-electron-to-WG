
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila
EELS
see paper#228 Eqs. 3
see paper#149 Eqs. 25
"""
import numpy as np
from scipy.integrate import quad
from global_constants import constants 

#%%

hb,c,alpha,me_c2_eV = constants()
aux = hb*c
epsi1, epsi3 = 1, 1

#%%
# print('Define the EELS')

def Fresnel_coefficient(omegac,u,d,mode,epsi2):
    """    
    Parameters
    ----------
    omegac : omega/c in 1/microns  
    u: k_parallel/(omega/c) 
    d: thickness of the plane in microns
    mode: s o p
    epsi2: permittivity of medium 2 
    (add small imaginary part to epsilon2)
    Returns
    -------
    r_{123}^mode from paper 370 Eq. C2
    """
    k = omegac
 
    def epsilon(i):
        if i == 1:
            return epsi1
        elif i == 2: 
            return epsi2 
                                
        elif i == 3:
            return epsi3 
     
    ## for python add +1j*0 inside kz as \sqrt{ .. + 1j*0}  
    
    # def k(i):
    #     return np.sqrt(epsilon(i))*k
            
    def kz(i): ## divided by k, but this will cancel on the fresnel coeff
        value = np.sqrt(epsilon(i) - u**2  + 1j*0) ## u: k_parallel/(omega/c)
        if np.imag(value) <= 0: ## has to be positive for convergency 
            value = - value
        return value
    
    def rij(i,j,mode):
        rij_s =  (kz(i) - kz(j))/(kz(i) + kz(j))
        rij_p =  (epsilon(j)*kz(i) - epsilon(i)*kz(j))/(epsilon(j)*kz(i) + epsilon(i)*kz(j))
        
        if mode == 's':
            return rij_s
        else:
            return rij_p
        

    def tij(i,j,mode):
        tij_s =  2*kz(i)/(kz(i) + kz(j))
        tij_p =  2*np.sqrt(epsilon(i)*epsilon(j))*kz(i)/(epsilon(j)*kz(i) + epsilon(i)*kz(j))
        
        if mode == 's':
            return tij_s
        else:
            return tij_p
    
    def r123(mode):
        
        num = tij(1,2,mode)*tij(2,1,mode)*rij(2,3,mode)*np.exp(2*1j*kz(2)*k*d) ## kz is defined divided by k
        den = 1 - rij(2,1,mode)*rij(2,3,mode)*np.exp(2*1j*kz(2)*k*d)           ## kz is defined divided by k
        
        r123 =  rij(1,2,mode) + num/den
        return r123

    return  r123(mode)


## EELS under QE approximation
def EELS_QE(energy,kx_norm_k,ze,d,beta,epsi2):
    """    
    Parameters
    ----------
    energy : hbar*omega in eV
    kx_norm_k: k_x/(omega/c)
    ze: position of electron in microns
    d: thickness of the plane in microns
    beta: v/c
    epsi2: permittivity of medium 2 
    Returns
    -------
    Re(EELS) from paper 228 Eq. 3
    divided by L/c and in Gaussian units
    (dimensionless) without integration
    using the QE approx (kz = i*k_par)
    """
    # epsi2 = epsilon(hbw,material)
    L = 1
    
    omegac = energy/(aux)
    k = omegac
    
    u = np.sqrt(kx_norm_k**2 + (1/beta)**2)
    
    ## integration variable u = k_par_over_omegac (dimensionless)
    r123_s =   np.imag(Fresnel_coefficient(omegac,u,d,'s',epsi2))
    r123_p =   np.imag(Fresnel_coefficient(omegac,u,d,'p',epsi2))

    denominator = np.sqrt(u**2 - (1/beta)**2 + 1j*0)
    function_s =   np.exp(-2*u*k*ze)*r123_s/denominator
    function_p =   np.exp(-2*u*k*ze)*r123_p/denominator
    # print(denominator, np.exp(-2*u*k*ze))
    
    # print(r123_s,r123_p,function_s,function_p)
    
    # print(Integral_s,Integral_p)

    factor_Gamma_norm_Lc  = alpha*2*L/(np.pi*beta**2) ## without L/c multipliying

    return  (function_p+function_s)*factor_Gamma_norm_Lc

# EELS per unit lenght
def EELS_no_QE(energy,kx_norm_k,ze,d,beta,epsi2):
    """    
    Parameters
    ----------
    energy : hbar*omega in eV
    kx_norm_k: k_x/(omega/c)
    ze: position of electron in microns
    d: thickness of the plane in microns
    beta: v/c
    epsi2: permittivity of medium 2 
    Returns
    -------
    Re(EELS) from paper 149 Eq. 25
    divided by Leff/c and in Gaussian units
    (dimensionless) without integration
    no QE approximation made
    """
    # epsi2 = epsilon(hbw,material)
    L = 1
    
    omegac = energy/(aux)
    k = omegac
    
    def epsilon(i):
        if i == 1:
            return epsi1
        elif i == 2: 
            return epsi2 
                                
        elif i == 3:
            return epsi3 
    
    u = np.sqrt(kx_norm_k**2 + (1/beta)**2)
    ## for python add +1j*0 inside kz as \sqrt{ .. + 1j*0}  
    # kx_norm_k = np.sqrt(u**2 - (1/beta)**2 + 1j*0)
    # def k(i):
    #     return np.sqrt(epsilon(i))*k
            
    def kz(i): ## divided by k, but this will cancel on the fresnel coeff
        value = np.sqrt(epsilon(i) - u**2  + 1j*0) ## u: k_parallel/(omega/c)
        if np.imag(value) <= 0: ## has to be positive for convergency 
            value = - value
        return value
    

    ## integration variable u = k_par_over_omegac (dimensionless)
    r123_s =  Fresnel_coefficient(omegac,u,d,'s',epsi2)
    r123_p =  Fresnel_coefficient(omegac,u,d,'p',epsi2)

    factor_s = (kx_norm_k*beta/kz(1))**2
    final_function =  np.real(kz(1)*np.exp(2*1j*kz(1)*k*ze)*(r123_s*factor_s - r123_p/epsi1))/(u**2)
    # final_function =  np.real(kz(1)*np.exp(2*1j*kz(1)*k*ze)*(r123_s*factor_s - r123_p/epsi1))/(u*kx_norm_k)
    
    factor_Gamma_norm_Lc  = alpha*2*L/(np.pi*beta**2) ## without L/c multipliying
    
    return final_function*factor_Gamma_norm_Lc/omegac

# EELS per unit lenght
def EELS_integrated_over_k_par_QE(energy,ze,d,beta,epsi2):
    """    
    Parameters
    ----------
    energy : hbar*omega in eV
    ze: position of electron in microns
    d: thickness of the plane in microns
    beta: v/c
    epsi2: permittivity of medium 2 
    Returns
    -------
    Re(EELS) from paper 228 Eq. 3
    divided by Leff/c and in Gaussian units
    (dimensionless)
    using the QE approx (kz = i*k_par)
    """
    # epsi2 = epsilon(hbw,material)
    L = 1
    omegac = energy/(aux)
    k = omegac
    
    ## integration variable u = k_par_over_omegac (dimensionless)
    r123_s = lambda u : np.imag(Fresnel_coefficient(omegac,u,d,'s',epsi2))
    r123_p = lambda u : np.imag(Fresnel_coefficient(omegac,u,d,'p',epsi2))
    
    function_s = lambda u : np.real(1/np.sqrt(u**2 - (1/beta)**2 + 1j*0 ))*np.exp(-2*u*k*ze)*r123_s(u)
    function_p = lambda u : np.real(1/np.sqrt(u**2 - (1/beta)**2 + 1j*0 ))*np.exp(-2*u*k*ze)*r123_p(u)
    
    limit1 = 1.001*(1/beta) ## integral from omega/v
    
    limit2 = 1.3*(1/beta)   ## already zero for this upper limit
    limit2 = np.real(np.sqrt(epsi2)) ## inside light cone
 
    Integral_s = quad(function_s, limit1, limit2)[0]
    Integral_p = quad(function_p, limit1, limit2)[0]
    
    # print(Integral_s,Integral_p)

    factor_Gamma_norm_Lc  = alpha*2*L/(np.pi*beta**2) ## without L/c multipliying

    return  (Integral_s + Integral_p)*factor_Gamma_norm_Lc 


# EELS per unit lenght
def EELS_integrated_over_k_par_no_QE(energy,ze,d,beta,epsi2):
    """    
    Parameters
    ----------
    energy : hbar*omega in eV
    ze: position of electron in microns
    d: thickness of the plane in microns
    beta: v/c
    epsi2: permittivity of medium 2 
    Returns
    -------
    Re(EELS) from paper 149 Eq. 25
    divided by L/c and in Gaussian units
    (dimensionless) integrated over kx
    """
    # epsi2 = epsilon(hbw,material)
    L = 1
    omegac = energy/(aux)
    k = omegac
    
    u = lambda kx_norm_k : np.sqrt(kx_norm_k**2 + (1/beta)**2) 

    ## for python add +1j*0 inside kz as \sqrt{ .. + 1j*0}  
    
    argument = lambda kx_norm_k : np.sqrt(epsi1 - u(kx_norm_k)**2  + 1j*0)
    kz1 = lambda kx_norm_k : argument(kx_norm_k) if np.imag(argument(kx_norm_k))>0  else  - argument(kx_norm_k) ## u: k_parallel/(omega/c)
    
    # if np.imag(kz1) <= 0:
    #     kz1 = - kz1
    
    

    ## integration variable u = k_par_over_omegac (dimensionless)
    r123_s = lambda qx:  Fresnel_coefficient(omegac,u(qx),d,'s',epsi2)
    r123_p = lambda qx: Fresnel_coefficient(omegac,u(qx),d,'p',epsi2)

 

    
    final_function = lambda qx : kz1(qx)*np.exp(2*1j*kz1(qx)*k*ze)*(r123_s(qx)*(qx*beta/kz1(qx))**2 - r123_p(qx)/epsi1)/(u(qx)**2)
    final_function_re = lambda qx : np.real(final_function(qx))


    factor_Gamma_norm_Lc  = alpha*2*L/(np.pi*beta**2) ## without L/c multipliying
    
    
    limit1 = 0.001*omegac ## variable is qx integral from 0
    
    # limit2 = 1.3*(1/beta)   ## already zero for this upper limit
    # limit2 = np.real(np.sqrt(epsi2)) ## inside light cone
 
    limit2 = 50*omegac
    
    limit1 = 0.001 ## variable is qx integral from 0
    limit2 = 50
    
    
    Integral = quad(final_function_re, limit1, limit2)[0]
 
    
    return  Integral*factor_Gamma_norm_Lc 


# EELS integrated over y-coordinate
# the other EELS are per unit lenght
def EELS_integrated_over_electron_trayectory(energy,b,d,beta,epsi2):
    """    
    Parameters
    ----------
    energy : hbar*omega in eV
    b: minimum position of electron along z coordinate in microns
    (most close to the plane)
    d: thickness of the plane in microns
    beta: v/c
    epsi2: permittivity of medium 2 
    Returns
    -------
    Re(EELS) from paper 149 Eq. 25
    divided by L0 and in Gaussian units
    (seconds/microns) integrated over
    electron trajectory and over k_parallel
    (see notes)
    """
    # epsi2 = epsilon(hbw,material)

    omegac = energy/(aux)
    k = omegac
    
    u = lambda kx_norm_k : np.sqrt(kx_norm_k**2 + (1/beta)**2) ## u: k_parallel/(omega/c)

    ## for python add +1j*0 inside kz as \sqrt{ .. + 1j*0}  
    
    argument = lambda kx_norm_k : np.sqrt(epsi1 - u(kx_norm_k)**2  + 1j*0) ## kz = sqrt(k^2 - k_parallel^2)
    kz1 = lambda kx_norm_k : argument(kx_norm_k) if np.imag(argument(kx_norm_k))>0  else  - argument(kx_norm_k) 
    
    # if np.imag(kz1) <= 0:
    #     kz1 = - kz1
    
    

    ## integration variable u = k_par_over_omegac (dimensionless)
    r123_s = lambda qx:  Fresnel_coefficient(omegac,u(qx),d,'s',epsi2)
    r123_p = lambda qx: Fresnel_coefficient(omegac,u(qx),d,'p',epsi2)

 

    
    final_function = lambda qx : kz1(qx)*np.exp(2*1j*kz1(qx)*k*b)*(r123_s(qx)*(qx*beta/kz1(qx))**2 - r123_p(qx)/epsi1)/(u(qx)**(5/2))
    final_function_re = lambda qx : np.real(final_function(qx))

    gamma_e = 1/(np.sqrt(epsi1-beta**2))
    me_over_hb = 8.648539271254356*1e-9 ## seconds/microns^2
    
    omega_sqrt = np.sqrt(omegac*c) 
    factor_Gamma_norm_L0  = alpha*2*np.sqrt(gamma_e)*np.sqrt(me_over_hb)/(np.pi*beta*omega_sqrt) ## Gamma/L0 in unis of seconds/microns
    ## the sqrt(omega) comes from the fact that in the integral I am using dimensionless variables instead of k_par, I use k_par/k
    
    limit1 = 0.001*omegac ## variable is qx integral from0 omega/v
    
    # limit2 = 1.3*(1/beta)   ## already zero for this upper limit
    # limit2 = np.real(np.sqrt(epsi2)) ## inside light cone
 
    limit2 = 50*omegac
    
    
    Integral = quad(final_function_re, limit1, limit2)[0]
    
    return  Integral*factor_Gamma_norm_L0 
 



# function before is without integratation over kx/k
# EELS integrand over y-coordinate: Leff(k_par)
# the other EELS are per unit lenght
def EELS_integrand_over_electron_trayectory(kx_norm_k,energy,b,d,beta,epsi2):
    """    
    Parameters
    ----------
    kx_norm_k : kx/k
    energy : hbar*omega in eV
    b: minimum position of electron along z coordinate in microns
    (most close to the plane)
    d: thickness of the plane in microns
    beta: v/c
    epsi2: permittivity of medium 2 
    Returns
    -------
    Re(EELS) from paper 149 Eq. 25
    divided by L0 and in Gaussian units
    (seconds/microns) using Leff 
    as a function of kx/k
    (see notes)
    """
    # epsi2 = epsilon(hbw,material)

    omegac = energy/(aux)
    k = omegac
    
    u =   np.sqrt(kx_norm_k**2 + (1/beta)**2) ## u: k_parallel/(omega/c)

    ## for python add +1j*0 inside kz as \sqrt{ .. + 1j*0}  
    
    argument =  np.sqrt(epsi1 - u**2  + 1j*0) ## kz = sqrt(k^2 - k_parallel^2)
    kz1 =   argument if np.imag(argument)>0  else  - argument 
    
    # if np.imag(kz1) <= 0:
    #     kz1 = - kz1
    
    

    ## integration variable u = k_par_over_omegac (dimensionless)
    r123_s =  Fresnel_coefficient(omegac,u,d,'s',epsi2)
    r123_p =  Fresnel_coefficient(omegac,u,d,'p',epsi2)

 

    
    final_function =   kz1*np.exp(2*1j*kz1*k*b)*(r123_s*(kx_norm_k*beta/kz1)**2 - r123_p/epsi1)/(u**(5/2))
    final_function_re =   np.real(final_function)

    gamma_e = 1/(np.sqrt(epsi1-beta**2))
    me_over_hb = 8.648539271254356*1e-9 ## seconds/microns^2
    
    omega_sqrt = np.sqrt(omegac*c) 
    factor_Gamma_norm_L0  = alpha*2*np.sqrt(gamma_e)*np.sqrt(me_over_hb)/(np.pi*beta*omega_sqrt) ## Gamma/L0 in unis of seconds/microns
    ## the sqrt(omega) comes from the fact that in the integral I am using dimensionless variables instead of k_par, I use k_par/k
    
    limit1 = 0.001*omegac ## variable is qx integral from0 omega/v
    
    # limit2 = 1.3*(1/beta)   ## already zero for this upper limit
    # limit2 = np.real(np.sqrt(epsi2)) ## inside light cone
 
    limit2 = 50*omegac
 
    
    return  final_function_re*factor_Gamma_norm_L0 
 