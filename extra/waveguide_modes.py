
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila
 dispersion relation
 of waveguide modes
 see paper#370 Eqs. D1
 reproduce figure S6
"""
real_material = 0
from global_constants import constants 
import numpy as np
import matplotlib.pyplot as plt
 
if real_material == 1:
    from permitivity_epsilon import epsilon as epsilon2
    label_png = '_real'
else:
    material = 'Si'
    # material = 'Ge'
    material = 'generic'
    
    if material == 'Si':
        epsilon2 = 2.13
    elif material == 'Ge':
        epsilon2 = 17.38
    else:
        epsilon2 = 12
    label_png = ''
    Im_epsi2 = 1e-1 
    
#%%
hb,c,alpha,me_c2_eV = constants()
aux = hb*c
epsi1,epsi3 = 1,1

d_microns = 1*1e-3
d = d_microns
    
omegac_WG = (np.pi/d)/np.sqrt(epsilon2-1) ## omega_WG/c

#%%
print('Define the eqs of the dispersion relation')

def Eq1_disp_relation(k_par_over_omegac_WG,omega_omega_WG,d,mode,Im_epsi2):
    """    
    Parameters
    ----------
    k_parallel*c/omega_WG : k_parallel divided by omega_{WG}/c  
    omega/omega_WG : omega/omega_WG  
    d: thickness of the plane in microns
    mode: s (TE) o p (TM)
    Im_epsi2: Im(epsilon_2)
    Returns
    -------
    Eqs. (D1)a paper 370
    """
    omegac = omega_omega_WG*omegac_WG
    k = omegac
    hbw = omegac*aux
 
    
    if real_material == 1:
        epsi2 = epsilon2(hbw,material) 
    else:
        epsi2 = epsilon2 + 1j*Im_epsi2
     
    if mode == 's':
        eta = 1
    else:
        eta = epsi2
    
    k_par = k_par_over_omegac_WG*omegac_WG
    chi = d*np.sqrt(k_par**2 - k**2)
    chi_prime = d*np.sqrt(k**2*epsi2 - k_par**2  )
    
    eq = eta*(chi/chi_prime) - np.tan(chi_prime/2)
    
    # k_par_inf_limit = k
    # k_par_sup_limit = k*np.sqrt(epsi2)
    
    # if k_par_inf_limit <= k_par <= k_par_sup_limit:
    if k < k_par <  k*np.sqrt(epsi2):
    # if k_par/np.sqrt(epsi2) <= k <= k_par:     
        return  np.abs(eq)
    else:
        return np.nan


def Eq2_disp_relation(k_par_over_omegac_WG,omega_omega_WG,d,mode,Im_epsi2):
    """    
    Parameters
    ----------
    k_parallel*c/omega_WG : k_parallel divided by omega_{WG}/c  
    omega/omega_WG : omega/omega_WG 
    d: thickness of the plane in microns
    mode: s (TE) o p (TM)
    Im_epsi2: Im(epsilon_2)
    Returns
    -------
    Eqs. (D1)b paper 370
    """
    omegac = omega_omega_WG*omegac_WG
    k = omegac
    hbw = omegac*aux
    
    if real_material == 1:
        epsi2 = epsilon2(hbw,material) 
    else:
        epsi2 = epsilon2 + 1j*Im_epsi2
    
    if mode == 's':
        eta = 1
    else:
        eta = epsi2
        
    k_par = k_par_over_omegac_WG*omegac_WG    
    chi = d*np.sqrt(k_par**2 - k**2)
    chi_prime = d*np.sqrt(k**2*epsi2 - k_par**2)
    
    eq = eta*(chi/chi_prime) + 1/np.tan(chi_prime/2)
    
    # k_par_inf_limit = k
    # k_par_sup_limit = k*np.sqrt(epsi2)
    
    if k < k_par <  k*np.sqrt(epsi2):
    # if k_par/np.sqrt(epsi2) <= k <= k_par:     
        return np.abs(eq)
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

print('Plot color map of the dispersion relation as 1/|Eq| to see the solutions more clear')

labely =  r'Frequency $\omega/\omega_{\text{WG}}$'
labelx = r'Parallel wave vector $k_\parallel c/\omega_{\text{WG}}$'

N = 100
listx = np.linspace(0.1,8,N)
listy = np.linspace(0.1,4,N)

def TE_modes_1(k_par_over_omegac_WG,omega_omega_WG):
    value = Eq1_disp_relation(k_par_over_omegac_WG,omega_omega_WG,d,'s',Im_epsi2)
    return 1/value

def TE_modes_2(k_par_over_omegac_WG,omega_omega_WG):
    value = Eq2_disp_relation(k_par_over_omegac_WG,omega_omega_WG,d,'s',Im_epsi2)
    return 1/value

def TM_modes_1(k_par_over_omegac_WG,omega_omega_WG):
    value = Eq1_disp_relation(k_par_over_omegac_WG,omega_omega_WG,d,'p',Im_epsi2)
    return 1/value

def TM_modes_2(k_par_over_omegac_WG,omega_omega_WG):
    value = Eq2_disp_relation(k_par_over_omegac_WG,omega_omega_WG,d,'p',Im_epsi2)
    return 1/value

title = r'1/|Eq(D1a)| for TE  $d = %i$ nm' %(d*1e3) 
title = r'1/|Eq(D1b)| for TE  $d = %i$ nm' %(d*1e3) 

title = r'1/|Eq(D1a)| for TM  $d = %i$ nm' %(d*1e3) 
title = r'1/|Eq(D1b)| for TM  $d = %i$ nm' %(d*1e3) 

X, Y = np.meshgrid(listx,listy)
f_TM1 = np.vectorize(TM_modes_1)
Z_TM1 = f_TM1(X, Y)

f_TM2 = np.vectorize(TM_modes_2)
Z_TM2 = f_TM2(X, Y)

f_TE1 = np.vectorize(TE_modes_1)
Z_TE1 = f_TE1(X, Y)

f_TE2 = np.vectorize(TE_modes_2)
Z_TE2 = f_TE2(X, Y)

#%%

Z_TM = (np.array(Z_TE1)/np.nanmax(Z_TE1))

import matplotlib as mpl

limits = [np.min(listx) , np.max(listx),np.min(listy) , np.max(listy)]
cmap = plt.cm.RdBu  # define the colormap
n_color = 21
vmin1, vmax1 = np.nanmin(Z_TM), np.nanmax(Z_TM)
bounds =   np.linspace(vmin1, 0.01, n_color) 
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

plt.figure(figsize=tamfig)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
im_TM = plt.imshow(Z_TM, extent = limits, cmap=cmap, aspect='auto', interpolation = 'bicubic',origin = 'lower' ,norm = norm) 
cbar = plt.colorbar(im_TM, fraction=0.046, pad=0.04, orientation = 'vertical' )
plt.title(title,fontsize=tamtitle)
plt.show()

#%%

print('Finding modes of the dispersion relation')

from scipy.optimize import minimize 

## I think red curves are TM and blue curves are TE curves

# omega/omega_WG
list_freq_norm_TM1_1 = np.linspace(1,2.7,N)  # second red in S6
list_freq_norm_TM1_2 = np.linspace(3,4,N)    # fourth red in S6

list_freq_norm_TM2_1 = np.linspace(0.01,2.7,N)  # first red in S6
list_freq_norm_TM2_2 = np.linspace(2,3.5,N)     # third red in S6


list_freq_norm_TE1_1 = np.linspace(1,2.7,N)  # second blue in S6
list_freq_norm_TE1_2 = np.linspace(3,4,N)    # fourth blue in S6

list_freq_norm_TE2_1 = np.linspace(0.01,2.7,N)  # first blue in S6
list_freq_norm_TE2_2 = np.linspace(2,3.5,N)     # third blue in S6

solutions_TM1_1 = []
solution_eq_TM1_1 = []

mode = 'p'
x00 = 1e-6
omegac = list_freq_norm_TM1_1[0]
cond_inicial = [omegac/omegac_WG]

for omega_omega_WG  in list_freq_norm_TM1_1: 
    
    def solve_Eq1(k_par_over_omegac_WG):
        return Eq1_disp_relation(k_par_over_omegac_WG,omega_omega_WG,d,mode)
    
    res = minimize(solve_Eq1, cond_inicial, method='Nelder-Mead', tol=1e-11, 
                     options={'maxiter':1050})
  #        print(res.message)
    if res.message == 'Optimization terminated successfully.':
        solution_eq_TM1_1.append(res.x[0])
 
        solutions_TM1_1.append(solve_Eq1(res.x[0] ))
        
        cond_inicial =  res.x[0] 
        
        print(res)
        
 #%%   
 