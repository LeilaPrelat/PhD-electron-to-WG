
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila
EELS
see paper#228 Eqs. 3
see paper#149 Eqs. 25
integrated over electron's trajectory
for real materials (\epsilon(\omega))
using Leff(k_par) = L0*f(k_par)
plot EELS/L0
"""
import numpy as np
import matplotlib.pyplot as plt
from global_constants import constants 
import os
from EELS import EELS_integrated_over_electron_trayectory, EELS_integrand_over_electron_trayectory
from permittivity_epsilon import epsilon as epsilon2

real_units = 1     ## Gamma in real units or normalized by c (then Gamma dimensionless)
    
label_png = '_real' 
material = 'Si'     ## default
# material = 'Ge'  
zoom = 1

pwd = os.path.dirname(__file__) 
path_save =  os.path.join(pwd,'plots_EELS')

#%%
hb,c,alpha,me_c2_eV = constants()
aux = hb*c
epsi1, epsi3 = 1, 1

d_microns = 0.2 # microns
d = d_microns
    
## list of electron energies from jga notes 2025-04-30 ##
ind = 2
list_Ee_electron = [30 , 100 , 200]   ## keV
Ee_electron_keV = list_Ee_electron[ind]
Ee_electron = Ee_electron_keV*1e3
label_Ee = label_png + '_d%inm_Ee%i' %(d*1e3,ind+1)

beta = np.sqrt( 1- (1 + Ee_electron/me_c2_eV)**(-2) )  ## beta = v/c
gamma_e = 1/np.sqrt(1-epsi1*beta**2)

N = 100
if d == 0.2:
    if zoom == 0:
        list_energy_eV = np.linspace(0.1,10,N) ## cutoff energy
    else:
        list_energy_eV = np.linspace(0.1,2,N) ## cutoff energy
else:
    if zoom == 0:
        list_energy_eV = np.linspace(0.1,20,N) ## cutoff energy
    else:
        list_energy_eV = np.linspace(0.1,4,N) ## cutoff energy

list_b_nm = [50,100,200]
list_b_nm = [10,50,80]

total_label = material + label_png + label_Ee  + 'zoom%i' %(zoom)

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

energy0 = 50
omegac0 = energy0/aux
list_qx = np.linspace(1e-4*omegac0, 30*omegac0, N)
# list_qx = np.linspace(0.001, 50, N)
# list_qx = np.logspace(-3, 2, N)


print('1-Plot the EELS with Leff as a function of kx/k, for different b, for fix energy = %.2f eV'  %(energy0 ) )

list_EELS_2_re_tot = []
list_EELS_2_im_tot = []
for b_nm in list_b_nm:
    b = b_nm*1e-3
    list_EELS_2_re = []
    list_EELS_2_im = []
    for qx in list_qx: 
        epsi2 = epsilon2(energy0,material) 
        value = EELS_integrand_over_electron_trayectory(qx,energy0,b,d,beta,epsi2)
        # value = Fresnel_coefficient(omegac,u,d,mode,Im_epsi2)
        list_EELS_2_re.append(np.real(value))
        list_EELS_2_im.append(np.imag(value))
        
    list_EELS_2_re_tot.append(list_EELS_2_re)
    list_EELS_2_im_tot.append(list_EELS_2_im)


labelx = r'$k_x/k$'
labely = r'$\Gamma_{\parallel}(k_x)/L_0$ (s/$\mu$m)'

title = r'EELS for $\hbar\omega = %.2f$ eV, $h = %.1f$ $\mu$m, $\epsilon_2 = \epsilon_{%s}(\omega)$, $v = %.2fc$' %(energy0,d,material,beta)
 
plt.figure(figsize=tamfig)
plt.title(title,fontsize=tamtitle)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
for j in range(len(list_b_nm)):
    plt.plot(list_qx, np.array(list_EELS_2_re_tot[j]) ,'.-',lw = 1.5,label = r'$b = %i$ nm' %(list_b_nm[j]) )
 
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1) 
# plt.xscale('log')
plt.yscale('log')
plt.show() 

#%%

print('2-Plot the EELS integrated over k_par and the trajectory as a function of energy, for different b')

list_EELS_re_tot = []
list_EELS_im_tot = []
for b_nm in list_b_nm:
    b = b_nm*1e-3
    list_EELS_re = []
    list_EELS_im = []
    for eV in list_energy_eV: 
        epsi2 = epsilon2(eV,material) 
        value = EELS_integrated_over_electron_trayectory(eV,b,d,beta,epsi2)
        # value = Fresnel_coefficient(omegac,u,d,mode,Im_epsi2)
        list_EELS_re.append(np.real(value))
        list_EELS_im.append(np.imag(value))
        
    list_EELS_re_tot.append(list_EELS_re)
    list_EELS_im_tot.append(list_EELS_im)

#%% 

labelx = r'Electron energy $\hbar\omega$ (eV)'
labely = r'$\Gamma_{\parallel}/L_0$ (s/$\mu$m)'

title = r'EELS for $h = %.1f$ $\mu$m, $\epsilon_2 = \epsilon_{%s}(\omega)$, $v = %.2fc$' %(d,material,beta)
 
plt.figure(figsize=tamfig)
plt.title(title,fontsize=tamtitle)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
for j in range(len(list_b_nm)):
    plt.plot(list_energy_eV, np.array(list_EELS_re_tot[j]) ,'.-',lw = 1.5,label = r'$b = %i$ nm' %(list_b_nm[j]) )
 
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1) 
label_figure = 'EELS_tot_' + material + label_Ee
# plt.xscale('log')
os.chdir(path_save)
plt.savefig(label_figure + '.png', format='png',bbox_inches='tight',pad_inches = 0.04, dpi=dpi)  
plt.show() 
