
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila
EELS
see paper#228 Eqs. 3
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

def EELS_function(energy,u,ze,d,beta,epsi2):
    """    
    Parameters
    ----------
    energy : hbar*omega in eV
    u: k_parallel/(omega/c)
    ze: position of electron in microns
    d: thickness of the plane in microns
    beta: v/c
    epsi2: permittivity of medium 2 
    Returns
    -------
    Re(EELS) from paper 228 Eq. 3
    divided by L/c and in Gaussian units
    (dimensionless) without integration
    """
    # epsi2 = epsilon(hbw,material)
    L = 1
    
    omegac = energy/(aux)
    k = omegac
    
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

    return  (function_s + function_p)*factor_Gamma_norm_Lc



def EELS_integrated_over_k_par(energy,ze,d,beta,epsi2):
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
    divided by L/c and in Gaussian units
    (dimensionless)
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

#%%

 