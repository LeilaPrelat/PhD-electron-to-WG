
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila
EELS
see paper#228 Eqs. 3
see paper#149 Eqs. 25
for real materials (\epsilon(\omega))
"""
import numpy as np
import matplotlib.pyplot as plt
from global_constants import constants 
import os
from EELS import EELS_QE, EELS_no_QE
from permittivity_epsilon import epsilon as epsilon2

    
create_data = 1      ## run data for the color maps 
 
real_units = 1 ## Gamma in real units or normalized by c (then Gamma dimensionless)
    
zoom = 0 # zoom to the dispersion relation
   
label_png = '_real'
material = 'Si'     ## default
# material = 'Ge'
 
pwd = os.path.dirname(__file__) 
path_save =  os.path.join(pwd,'plots_EELS_comparison')

#%%
hb,c,alpha,me_c2_eV = constants()
aux = hb*c
epsi1, epsi3 = 1, 1

#omegac_WG = (np.pi/d)/np.sqrt(np.real(epsilon2)-1) ## omega_WG/c

d_microns = 0.1 
d = d_microns
L = 1 # lenght of propagation in microns 
# propagation<infty if im(epsi2)!=0

    
## list of electron energies from jga notes 2025-04-30 ##
ind = 0
list_Ee_electron = [30 , 100 , 200]   ## keV
Ee_electron_keV = list_Ee_electron[ind]
Ee_electron = Ee_electron_keV*1e3
label_Ee = '_Ee%i' %(ind+1)

beta = np.sqrt( 1- (1 + Ee_electron/me_c2_eV)**(-2) )  ## beta = v/c
gamma_e = 1/np.sqrt(1-epsi1*beta**2)

N = 100
if zoom == 0:
    list_energy_eV = np.linspace(0.01,50,N)
else:
    list_energy_eV = np.linspace(0.01,4.5,N)
list_ze_nm =  np.linspace(0.1,200,N)
list_ze_nm =  np.logspace(-1,2,N)

energy_0 = 1

energy_0 = 10
omegac_0 = energy_0/aux
ze_0 = 10*1e-3 ## microns 
ze_0 = 25*1e-3 ## microns 
ze_0 = 50*1e-3 ## microns 

epsi2 = epsilon2(energy_0,material) 

limit1 = 1.001*(1/beta) ## integral from omega/v
limit2 = np.real(np.sqrt(epsi2)) ## inside light cone
    
list_u =  np.linspace(limit1,1.4*(1/beta),N) ## integration of EELS is from k_parallel = \omega/v
list_u =  np.linspace(0.001*omegac_0,45,N) ## integration of EELS is from k_parallel = 0 to infty

beta = 0.9

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

print('1-Plot the EELS without integration as a function of k_parallel, energy = %i eV, ze = %i nm' %(energy_0,ze_0*1e3))

list_EELS1_re = []
list_EELS1_im = []

list_EELS2_re = []
list_EELS2_im = []

for u in list_u: 
    epsi2 = epsilon2(energy_0,material) 
    value1 = EELS_no_QE(energy_0,u,ze_0,d,beta,epsi2)    ## paper 149
    value2 = EELS_QE(energy_0,u,ze_0,d,beta,epsi2)       ## paper 228
    # value = Fresnel_coefficient(omegac,u,d,mode,Im_epsi2)
    list_EELS1_re.append(np.real(value1))
    list_EELS1_im.append(np.imag(value1))
    
    list_EELS2_re.append(np.real(value2))
    list_EELS2_im.append(np.imag(value2))
    
#%% 

labelx = r'Parallel wave vector normalized $k_\parallel/(\omega/c)$'
labely = r'$\Gamma_{\parallel}(k_\parallel) c/L$'
title1 = r'EELS for $d = %.1f$ $\mu$m, $\epsilon_2 = \epsilon_{%s}(\omega)$' %(d,material)
title2 = r'$\hbar\omega = %.1f$ eV, $z_{\text{e}} = %i$ nm, $v = %.2fc$,' %(energy_0,ze_0*1e3,beta)

os.chdir(path_save)
plt.figure(figsize=tamfig)
plt.title(title1 + '\n' + title2,fontsize=tamtitle)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
plt.plot(list_u, np.array(list_EELS1_re) ,'.-',lw = 1.5, label = 'paper #149')
plt.plot(list_u, np.array(list_EELS2_re) ,'.-',lw = 1.5, label = 'paper #228' )
aux_eje_y = np.linspace(np.min(list_EELS2_re),np.max(list_EELS2_re),N)
plt.plot(np.ones(N)*(1/beta)*1.001 , aux_eje_y,'--',color='black',lw = 1.5 ) ## lower limit of integration of EELS
plt.yscale('log')
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1) 
plt.savefig('comparisonEELS.png', format='png',bbox_inches='tight',pad_inches = 0.3, dpi=dpi)  
plt.show() 

#%%

 
