
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila
EELS
see paper#228 Eqs. 3
see paper#149 Eqs. 25
integrated over electron's trajectory
for real materials (\epsilon(\omega))
and over the frequency from 0 to omega_f : plot total EELS vs omega_f

using Leff(k_par) = L0*f(k_par)

use the data from Fig.2b (integration over kx) and integrate it over omega
"""
import numpy as np
import matplotlib.pyplot as plt
from global_constants import constants 
import os
from scipy.interpolate import CubicSpline,interp1d
 

create_data = 0     ## run data for the color maps 

label_png = '_real'
material = 'Si'   ## default
# material = 'Ge'  
zoom = 0

delta = 0.2 ## extra loss for the imaginary part of the permittivity
pwd = os.path.dirname(__file__) 
path_save =  os.path.join(pwd,'plots_EELS_permittivity_extra_loss_%.2f'%(delta))

#%%
hb,c,alpha,me_c2_eV = constants()
aux = hb*c
epsi1, epsi3 = 1, 1

d_microns = 0.2 # microns
d = d_microns
    
## list of electron energies from jga notes 2025-04-30 ##
ind = 2
list_Ee_electron = [30, 100 , 200]   ## keV . dont use the first value 
Ee_electron_keV = list_Ee_electron[ind]
Ee_electron = Ee_electron_keV*1e3
label_Ee = '_d%inm_Ee%i' %(d_microns*1e3,ind+1)

beta = np.sqrt( 1- (1 + Ee_electron/me_c2_eV)**(-2) )  ## beta = v/c
gamma_e = 1/np.sqrt(1-epsi1*beta**2)

energy_initial = 0.1 ## eV
step_energy = 0.0001
# step_energy = 0.00001
N = 1000
# if zoom == 0:
#     list_upper_eV_limit = np.arange(energy_initial,10,step_energy) ## cutoff energy
# else:
#     list_upper_eV_limit = np.linspace(energy_initial,2,N) ## cutoff energy

list_b_nm = [0,10,50]

 
total_label = material + label_png + label_Ee  + 'zoom%i' %(zoom)

list_upper_eV_limit_tot =  np.arange(energy_initial,50,step_energy)

#%%
 
tamfig = [3.45, 3]
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

print('1-Load data from Fig. 2b): EELS integrated over k_par and the trajectory as a function of energy, for different b')

os.chdir(path_save)    
list_energy_eV = np.loadtxt('list_energy_' + total_label + '.txt' , delimiter='\t',  skiprows = 1, encoding=None)

gamma_e = 1/(np.sqrt(epsi1-beta**2))
me_over_hb = 8.648539271254356*1e-9 ## seconds/microns^2

aux2 = np.sqrt(c)*hb ## add hb because the integral is over energy and not over frequency ------------------------------------- IMPORTANT 
factor_Gamma_norm_L0  = alpha*2*np.sqrt(gamma_e)*np.sqrt(me_over_hb)/(np.pi*beta*aux2) ## Gamma/L0 in unis of seconds/microns


title = r'EELS for $h = %.1f$ $\mu$m, $\epsilon_2 = \epsilon_{%s}(\omega)$, $v = %.2fc$' %(d,material,beta)
header = title + r', $\beta$ = %.2f. Re(EELS) from paper 149 Eq. 25 divided by L0 in Gaussian units (1/microns), integrated over k_par, trajectory, and eV' %(beta)


if create_data == 1:
    list_EELS_re_int_tot  = []
    os.chdir(path_save)
    for j, b_nm in enumerate(list_b_nm):
    
        list_EELS_re = np.loadtxt('list_EELS_over_L0_vs_energy_' + total_label + '_b%inm.txt' %(b_nm) , delimiter='\t', skiprows = 1, encoding=None)
     
    
        list_EELS_re_int = []
        
        #Build spline once
        spline = CubicSpline(list_energy_eV, list_EELS_re)
        f_lin = interp1d(list_energy_eV, list_EELS_re, kind="linear", fill_value="extrapolate")

       #  # Integrate directly using spline.integrate
       #  list_EELS_re_int = [
       #     spline.integrate(energy_initial, eV) / hb
       #     for eV in list_upper_eV_limit_tot
       # ]
        # Compute integrals for each upper limit
        list_EELS_re_int = []
        for eV in list_upper_eV_limit_tot:
            # Generate energy grid between initial and eV
            energy_grid = np.arange(energy_initial, eV, step_energy)
            values = f_lin(energy_grid)
            integral = np.trapz(values, energy_grid) / hb  # trapezoidal integration
            list_EELS_re_int.append(integral)


        ## DIVIDE BY HBAR BECAUSE THE INTEGRAL IS OVER ENERGY AND IT SHOULD BE OVER FREQUENCY ## 
        # print(b_nm, list_EELS_re_int)
   
    
        list_EELS_re_int_tot.append(list_EELS_re_int)
     
    for j in range(len(list_b_nm)):
        b_nm = list_b_nm[j]
        np.savetxt('list_EELS_over_L0_' + total_label + '_b%inm.txt' %(b_nm) ,list_EELS_re_int_tot[j], fmt='%.10f', delimiter='\t', header = header, encoding=None)
    np.savetxt('list_eV_upper_limit_' + total_label + '.txt' ,list_upper_eV_limit_tot, fmt='%.10f', delimiter='\t', header = header, encoding=None)

else:
    
    list_EELS_re_int_tot = []
    os.chdir(path_save)
    for j in range(len(list_b_nm)):
        b_nm = list_b_nm[j]    
        listy = np.loadtxt('list_EELS_over_L0_' + total_label + '_b%inm.txt' %(b_nm) , delimiter='\t', skiprows = 1)
        list_EELS_re_int_tot.append(listy)
    list_energy_eV = np.loadtxt('list_eV_upper_limit_' + total_label + '.txt' , delimiter='\t', skiprows = 1)

#%% 

print('2-Fig. 2c) Plot the EELS integrated over k_par and over omega_j, for different b')

from mycolorpy import colorlist as mcp
color1 = mcp.gen_color(cmap="hot",n=5)

labelx = r'Cutoff energy $\hbar\omega_f$ (eV)'
labely = r'$\Gamma/L_0$ (1/$\mu$m)'

index_cut = np.argmin(np.abs(np.array(list_upper_eV_limit_tot)-10))


plt.figure(figsize=tamfig)
# plt.title(title,fontsize=tamtitle)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)

for j in range(len(list_b_nm)):
    plt.plot(list_upper_eV_limit_tot[0:index_cut], np.array(list_EELS_re_int_tot[j][0:index_cut]) ,'-',color = color1[j],lw = 1.5,label = r'$b = %i$ nm' %(list_b_nm[j]) )
 
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.xticks(np.arange(0,12,2))
plt.yticks(np.arange(0,80,10), [ "0", "", "20", "", "40", "", "60", ""])
#plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1) 
label_figure = 'EELS_int_energy_' + total_label
# plt.xscale('log')
os.chdir(path_save)
plt.xlim(-0.5,10.5)
plt.ylim(-4,75)
plt.savefig(label_figure + '.png', format='png',bbox_inches='tight',pad_inches = 0.04, dpi=dpi)  
plt.savefig(label_figure + '.pdf', format='pdf',bbox_inches='tight',pad_inches = 0.04, dpi=dpi)  
plt.show() 

#%% 

print('3-Fig. 2d) Calculate the total integral from 0 to infinity to use it as normalization')
   
 
# list_EELS_norm_tot  = []
# os.chdir(path_save)
# for j in range(len(list_b_nm)):
#     b_nm = list_b_nm[j]    
 
 
 
#     list_EELS_re_interp_new = np.interp(list_upper_eV_limit_tot,list_upper_eV_limit, list_EELS_re_int_tot[j])  
    
#     list_EELS_norm_tot.append(np.sum(list_EELS_re_interp_new)*step_energy/hb)

plt.figure(figsize=tamfig)
# plt.title(title,fontsize=tamtitle)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
for j in range(len(list_b_nm)):
    maxi = np.max(list_EELS_re_int_tot[j])
    norm = list_EELS_re_int_tot[j][-1]
    norm = maxi
    plt.plot(list_upper_eV_limit_tot[0:index_cut], np.array(list_EELS_re_int_tot[j][0:index_cut])/maxi ,'-',color = color1[j],lw = 1.5,label = r'$b = %i$ nm' %(list_b_nm[j]) )
 
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.xticks(np.arange(0,12,2))
plt.xlim(-0.5,10.5)
plt.ylim(-0.06,1.1)
# plt.yticks(np.arange(0,180,20))
#plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1) 
label_figure = 'EELS_int_energy_norm' + total_label
# plt.xscale('log')

os.chdir(path_save)
plt.savefig(label_figure + '.png', format='png',bbox_inches='tight',pad_inches = 0.04, dpi=dpi)  
plt.savefig(label_figure + '.pdf', format='pdf',bbox_inches='tight',pad_inches = 0.04, dpi=dpi)  
plt.show() 

