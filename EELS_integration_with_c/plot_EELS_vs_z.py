
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila
calculate the EELS of a rectangular waveguide
integrated over the energy of a mode
and over z (total EELS)
as a function of theta/V0 (silutions for a fixed bmin)
"""
import numpy as np
import os
import matplotlib.pyplot as plt
# import concurrent.futures

## fit with a lorentzian EELS vs energy to find the width of a mode and integrate over energy later 

path_basic = os.getcwd()
path_data_EELS_integrated = os.path.join(path_basic, 'EELS_integrated/P_intengrand_vs_z')

compare_with_approx = 1 ## compare the results with the approximation wtih V(z)/V0 = a*z/W + b

#%%

tamfig = [4,3]
tamfig2 = [4,3.3]
tamletra = 16
tamtitle  = 8
tamnum = tamletra
tamlegend = tamletra
labelpady = 2
labelpadx = 3
pad = 3
mk = 1
ms = 2.5
hp = 0.3
length_marker = 0
dpi = 500
lw = 1.5

deltax,deltay = 3,8

values = tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad,deltax,deltay

#%%
print('1-Define the parameters')

############## BEM ############################################################## 
w = 500
Ee_electron_keV = 200
Ee_electron = Ee_electron_keV*1e3
position_of_the_electron = [0,0.4,0.6]
ind_xe=0
xe_norm_w = position_of_the_electron[ind_xe]
N=500
##### potential from c++ code  ##################################################
wp_norm_w = 0.5 ## w'/w
d_norm_w = 0.2  ## d/w distance between the side wires normalized to w
h_norm_w = 0.4   ## aspect ratio height/w
Nc = 400   ## discretization points.
x0 = xe_norm_w    ## V(xe,z) with z from z0 to z1
 
# x=0 is in the middle of the center waveguide
# z=0 is at the interface
                             
bmin_norm_w = 50/w
z_min_val_norm_w = bmin_norm_w + h_norm_w
###################################################################################
h=h_norm_w*w
 
   
label_txt = '_wp%.2f_h%.2f_d%.2f_xe%.1f' %(wp_norm_w,h_norm_w,d_norm_w,xe_norm_w) ## all normalize to width w
label_Ee = '_Ee%ikeV' %(Ee_electron_keV) 

energy = 0.8
V0 = 0.01331 ## eV
theta = 0.25858 ## mrad
 
all_parameters = '_energy%.4feV_w%inm_Ee%ikeV_wp%.2f_h%.2f_d%.2f_xe%.2f_V0_%.5feV_theta%.5fmrad_zmin%.2f' %(energy,w,Ee_electron_keV,wp_norm_w,h_norm_w,d_norm_w,xe_norm_w,V0,theta,z_min_val_norm_w)
total_name1 = 'P_integrand_vs_z%s' %(all_parameters)
total_name2 = 'P_integrand_vs_z_approx%s' %(all_parameters)

#%%
os.chdir(path_data_EELS_integrated)
print('Plot EELS vs z')

labelx=r'$z/W$'
labely=r'$\Gamma(\omega,x,z)$ (1/eV $\cdot$ nm)' 
#labely=r'$\Gamma$ (s)'

title1 = r'$W_{\text{p}}/W$ = %.2f, $h/W$ = %.2f, $d/W$ = %.2f, x/W = %.1f, $E_{\rm e}$ = %i keV' % (wp_norm_w,h_norm_w,d_norm_w,x0,Ee_electron_keV)
title2 = r'$z_{\text{min}} = %i$ nm, $W$ = %i nm, $\hbar\omega = %.2f$ eV, $V_0$ = %.4f eV, $\theta$ = %.4f mrad' %(z_min_val_norm_w*w,w, energy, V0, theta)

plt.figure(figsize=tamfig)
plt.title(title1 + '\n' + title2,fontsize=tamtitle)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady) 
table_P_integrand_over_energy = np.loadtxt(total_name1 + '.dat', delimiter='\t', skiprows=1, encoding=None)
table_P_integrand_over_energy2 = np.transpose(table_P_integrand_over_energy)
listz = table_P_integrand_over_energy2[0]
listP = table_P_integrand_over_energy2[1]
 
plt.plot(listz, listP ,'.-' ,label = r'$\phi$ from BEM' )
 
table_P_integrand_over_energy_approx = np.loadtxt(total_name2 + '.dat', delimiter='\t', skiprows=1, encoding=None)
table_P_integrand_over_energy2_approx = np.transpose(table_P_integrand_over_energy_approx)
listz_approx = table_P_integrand_over_energy2_approx[0]
listP_approx = table_P_integrand_over_energy2_approx[1]

plt.plot(listz_approx, listP_approx ,'--',label = r'$\phi$ linear' )
plt.xticks(np.arange(0.5,1.2+0.1,0.1),['', '0.6', '', '0.8', '', '1.0' , '', '1.2'])
plt.xlim([0.49,0.6])
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1)
plt.yscale('log')

plt.savefig('EELS_vs_z' + all_parameters + '.png',bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)
plt.savefig('EELS_vs_z' + all_parameters + '.pdf',bbox_inches='tight',pad_inches = 0.01, format='pdf', dpi=dpi)
