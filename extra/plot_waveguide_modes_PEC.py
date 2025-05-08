
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila
 dispersion relation
 of waveguide modes
 see paper#370 Eqs. D1
"""

import numpy as np
import matplotlib.pyplot as plt
 
# from permitivity_epsilon import epsilon
 
#%%
hb = 6.58211899*10**(-16)     ### Planck constant hbar in eV*s
c =  2.99792458*10**(14)      ### light velocity in micron/seg
aux = hb*c

material = 'Si'
# material = 'Ge'
epsi1,epsi3 = 1,1

if material == 'Si':
    epsi2 = 2.13
elif material == 'Ge':
    epsi2 = 17.38

epsi2 = 12
d_microns = 1*1e-3
d = d_microns
    
omegac_WG = (np.pi/d)/np.sqrt(epsi2-1) ## omega_WG/c

#%%
print('Define the eqs of the dispersion relation')

def TM_disp_relation(k_par_over_omegac_WG,omega_omega_WG):
    """    
    Parameters
    ----------
    k_parallel*c/omega_WG : k_parallel divided by omega_{WG}/c  
    omega/omega_WG : omega/omega_WG  
    Returns
    -------
    Eq. TM paper 370 for PEC
    """
    # epsi2 = epsilon(hbw,material)
    
    omegac = omega_omega_WG*omegac_WG
    k = omegac
    
    eta = epsi2
    
    k_par = k_par_over_omegac_WG*omegac_WG
    chi = d*np.sqrt(k_par**2 - k**2)
    chi_prime = d*np.sqrt(k**2*epsi2 - k_par**2)
    
    eq = eta*(chi/chi_prime) - np.tan(chi_prime/2)
    
    # k_par_inf_limit = k
    # k_par_sup_limit = k*np.sqrt(epsi2)
    
    # if k_par_inf_limit <= k_par <= k_par_sup_limit:
    if k < k_par <  k*np.sqrt(epsi2):
    # if k_par/np.sqrt(epsi2) <= k <= k_par:     
        return  eq
    else:
        return np.nan

def TE_disp_relation(k_par_over_omegac_WG,omega_omega_WG):
    """    
    Parameters
    ----------
    k_parallel*c/omega_WG : k_parallel divided by omega_{WG}/c  
    omega/omega_WG : omega/omega_WG   
    Returns
    -------
    Eq. TE paper 370 for PEC
    """
    # epsi2 = epsilon(hbw,material)
    
    omegac = omega_omega_WG*omegac_WG
    k = omegac
    
    eta = 1
    
    k_par = k_par_over_omegac_WG*omegac_WG
    chi = d*np.sqrt(k_par**2 - k**2)
    chi_prime = d*np.sqrt(k**2*epsi2 - k_par**2)
    
    eq = eta*(chi/chi_prime) + 1/np.tan(chi_prime/2)
    
    # k_par_inf_limit = k
    # k_par_sup_limit = k*np.sqrt(epsi2)
    
    # if k_par_inf_limit <= k_par <= k_par_sup_limit:
    if k < k_par <  k*np.sqrt(epsi2):
    # if k_par/np.sqrt(epsi2) <= k <= k_par:     
        return  eq
    else:
        return np.nan
 
#%%
tamfig = [4, 3]
tamletra = 13
tamtitle  = 10
tamnum = tamletra
tamlegend = tamletra
labelpady = 3
labelpadx = 2
pad = 2.5
mk = 1
ms = 1
lw = 1.5
hp = 0.5
length_marker = 1.5
dpi = 500

#%%
N = 100
print('Plot the dispersion relation')

labely =  r'Frequency $\omega/\omega_{\text{WG}}$'
labelx = r'Parallel wave vector $k_\parallel c/\omega_{\text{WG}}$'

listx = np.linspace(0.1,8,N)
listy = np.linspace(0.1,4,N)

X, Y = np.meshgrid(listx,listy)
f_TM = np.vectorize(TM_disp_relation)
Z_TM = f_TM(X, Y)

f_TE = np.vectorize(TE_disp_relation)
Z_TE = f_TE(X, Y)

#%%

import matplotlib as mpl

limits = [np.min(listx) , np.max(listx),np.min(listy) , np.max(listy)]
cmap = plt.cm.RdBu  # define the colormap
n_color = 21
vmin1, vmax1 = np.min(Z_TM), np.max(Z_TM)
bounds =   np.logspace(np.log10(vmin1), np.log10(vmax1), n_color) 
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

plt.figure(figsize=tamfig)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
im_TM = plt.imshow(Z_TM, extent = limits, cmap=cmap, aspect='auto', interpolation = 'bicubic',origin = 'lower' ) 
cbar = plt.colorbar(im_TM, fraction=0.046, pad=0.04, orientation = 'vertical' )
plt.title('TM for PEC',fontsize=tamtitle)
plt.show()


plt.figure(figsize=tamfig)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
im_TE = plt.imshow(Z_TE, extent = limits, cmap=cmap, aspect='auto', interpolation = 'bicubic',origin = 'lower' ) 
cbar = plt.colorbar(im_TE, fraction=0.046, pad=0.04, orientation = 'vertical' )
plt.title('TE for PEC',fontsize=tamtitle)
plt.show()

#%%

print('Finding the dispersion relation')

from scipy.optimize import minimize 

list_freq_norm = np.linspace(0.1,4,N) # omega/omega_WG
 

mode = 's'
solutions_Eq1_mode_s = []
solution_eq = []

x00 = 1e-6
omegac = list_freq_norm[0]
cond_inicial = [omegac]

for omega_omega_WG  in list_freq_norm: 
    
    def solve_Eq1(k_par_over_omegac_WG):
        return TE_disp_relation(k_par_over_omegac_WG,omega_omega_WG )
    
    res = minimize(solve_Eq1, cond_inicial, method='Nelder-Mead', tol=1e-11, 
                     options={'maxiter':1050})
  #        print(res.message)
    if res.message == 'Optimization terminated successfully.':
        solutions_Eq1_mode_s.append(res.x[0])
 
        solution_eq.append(solve_Eq1(res.x[0] ))
        
        cond_inicial =  res.x[0] 
        
        print(res)
        
 #%%   
 