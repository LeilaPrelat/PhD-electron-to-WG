
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila
calculate the EELS of a rectangular waveguide
integrated over the energy of a mode (\int\omega dEELS/dz)
as a function of h/W
"""
import numpy as np
import os
import matplotlib.pyplot as plt
from EELS_integrated import find_width_of_peak
# import concurrent.futures
from scipy.interpolate import interp1d

## fit with a lorentzian EELS vs energy to find the width of a mode and integrate over energy later 

path_basic = os.getcwd()
path_data = os.path.join(path_basic, 'bem2D_files_fixed_z')
path_data_created = os.path.join(path_basic, 'dEELSdz_integrated_over_omega')


from global_constants import constants 
hb,c,alpha,me_c2_eV = constants()
aux = hb*c

#%%

tamfig = [4.5,3.5]
tamletra = 14
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
# //     N   = the total number of parametrization points is N+50
# //     a   = total width of wg     || x
# //     h   = total thickness of wg  || y
# //     s   = rounding radius of wg (s<t/2)

print('Define the parameters')

############ parameters for bem2d ############
a=500 ## width 300 400 500
s=a/10 ## 10% of the width
list_h = np.arange(210,500,10)
# list_h = [300]
w_microns=a*1e-3

N = 400
Ee_electron_keV = 100
Ee_electron = Ee_electron_keV*1e3
label_Ee = '_Ee%ikeV' %(Ee_electron_keV)

step2=0.001 #integration over freq
mode=1
ze=50
list_modes=[1,2]
plot_figure = 0

position_of_the_electron = [0,40,80]
ind_xe=0
xe = position_of_the_electron[ind_xe]

if xe == 40:  
    labelxe ='_xe40W'
    xe_real = 0.4*a
elif xe == 80:
    labelxe ='_xe80W'
    xe_real = 0.8*a
else:
    labelxe ='_xe0'
    xe_real = 0
    

#%%

print('1-Plot EELS vs energy to fit it and find the width for all (V0,theta), 6-later use the width to integrate over energy')
#print('IMPORTANT: we are assuming the position of the peak DOES vary with V0,theta')

header0 = 'h/W    omega_j*W/c    omega_left*W/c    omega_right*W/c    \int\d\\Delta\omega_jdP/dz for j = %i' %(mode)
title1 = r'$W$ = %i nm, $s$ = %i nm, $z_{\text{e}} = %i$ nm, $x_{\text{e}} = %i$ nm, $E_{\rm e}$ = %i keV' %(a,s,ze,xe_real,Ee_electron_keV)

list_dP_integrated_over_mode_tot = []
list_omega_j_tot = []
list_h_fit_tot = []
for mode in list_modes: 
    
    os.chdir(path_data)
    
    list_dP_integrated_over_mode = []
    list_omega_j = []
    list_omega_left = []
    list_omega_right = []
    list_h_fit = []
    
    index_mode = -mode  ## if mode = 1, is the highest one because is the last element of the list (the array is sorted from min to max)
    for h in list_h: 
        
        name1 = 'EELS_epsi12_N%i_a%inm_h%inm_ze%inm_Ee%ikeV' %(N,a,h,ze,Ee_electron_keV) + labelxe + '.dat'
        tabla1 = np.loadtxt(name1, delimiter=' ',dtype=None)
        tabla1_2 = np.transpose(tabla1)
        
        listeV = tabla1_2[0]
        listEELS = tabla1_2[4]
        
        listeV_correct = []
        listEELS_correct = []
        for k in range(len(listEELS)):
            if listEELS[k]>0: 
                listeV_correct.append(listeV[k])
                listEELS_correct.append(listEELS[k])
        
        listeV_correct_2 = np.array(listeV_correct)*w_microns/aux
        listx = listeV_correct_2
        
        
 
        listy = listEELS_correct
        EELS_vs_energy = interp1d(listx, listy) ############### interpolation of P(omega) ###############
    
        print('2-Find the wmin and wmax of the mode = %i by fitting it with a Lorenztian' %(mode))
        
        try:
            ######### find peaks and sort them by the highest to minimum (then we identify each mode according to its amplitude) by fitting the curve with a lorentzian
            x_left, x_right, x0_fit,function_lorentzian, fit_success = find_width_of_peak(listx,listy,index_mode,title1,plot_figure)
            # Integration by summation (small steps)
            if fit_success:
            
            # x_integrate[0]>listx[0]: ## interpolation can be made
                x_integrate = np.arange(x_left, x_right + step2, step2)
                integration_over_energy = np.sum(EELS_vs_energy(x_integrate))*(step2)
                result_sum = integration_over_energy*a
                list_dP_integrated_over_mode.append(result_sum)
                list_omega_j.append(x0_fit)
                list_omega_left.append(x_left)
                list_omega_right.append(x_right)
                list_h_fit.append(h)
                print(h,integration_over_energy,function_lorentzian)
            
        # else: ## we approximate to the first point of listx if the difference is small (1e-6)
        #     if np.abs(x_integrate[0]-listx[0])<1e-6:
            # x_integrate = np.arange(listx[0], x_right + step2, step2)
            # integration_over_energy = np.sum(EELS_vs_energy(x_integrate))*(step2)

        
            # list_dP_integrated_over_mode.append(integration_over_energy)
            # list_omega_j.append(x0_fit)
            # list_omega_left.append(x_left)
            # list_omega_right.append(x_right)
            # list_h_fit.append(h)
            
        except (ValueError, TypeError) as e:
            print("Skipping due to error:", e)
            continue
        

        # except TypeError: 
        #     continue
            
    print('3-Save data of P integrated over energy as a function of h/W for mode =%i' %(mode))
    
    os.chdir(path_data_created)
    
    table = [list_h_fit, list_omega_j, list_omega_left, list_omega_right, list_dP_integrated_over_mode]
    table_P_integrated_over_energy = np.transpose(table)
    np.savetxt('dPdz_integrated_over_mode%i' %(mode) + label_Ee +'.txt', table_P_integrated_over_energy, fmt='%.10f', delimiter='\t', header = header0, encoding=None)
    

    list_dP_integrated_over_mode_tot.append(list_dP_integrated_over_mode)
    list_omega_j_tot.append(list_omega_j)
    list_h_fit_tot.append(list_h_fit)
    
#%%

labelx=r'$h/W$'

print('4-Plot P integrate over energy as a function of h/W')

label_png1 = 'freq_res_over_mode_w%inm' %(a) + label_Ee + labelxe  +'.png'
label_png2 = 'dPdz_integrated_over_mode_w%inm' %(a) + label_Ee + labelxe +'.png'

labely1=r'$\omega_j W/c$'
labely2=r'$W \, dP_j/dz$'

plt.figure(figsize=tamfig)
plt.title(title1,fontsize=tamtitle)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely1,fontsize=tamletra,labelpad =labelpady)
for j in range(len(list_modes)):
    mode = list_modes[j]
    list_h_norm_tot = np.array(list_h_fit_tot[j])/a
    plt.plot(list_h_norm_tot, list_omega_j_tot[j],'.-',label = 'mode = %i' %(mode))
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1)
# plt.yscale('log')
# plt.xlim([0,3.1])
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.tight_layout()
os.chdir(path_data_created)
plt.savefig(label_png1,bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)
plt.show()


labely1=r'$\hbar\omega_j$ (eV)'
plt.figure(figsize=tamfig)
plt.title(title1,fontsize=tamtitle)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely1,fontsize=tamletra,labelpad =labelpady)
for j in range(len(list_modes)):
    mode = list_modes[j]
    list_h_norm_tot = np.array(list_h_fit_tot[j])/a
    plt.plot(list_h_norm_tot, np.array(list_omega_j_tot[j])*aux/w_microns,'.-',label = 'mode = %i' %(mode))
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1)
# plt.yscale('log')
# plt.xlim([0,3.1])
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.tight_layout()
os.chdir(path_data_created)
plt.savefig(label_png1,bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)
plt.show()

#%%

plt.figure(figsize=tamfig)
plt.title(title1,fontsize=tamtitle)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely2,fontsize=tamletra,labelpad =labelpady)
for j in range(len(list_modes)):
    mode = list_modes[j]
    list_h_norm_tot = np.array(list_h_fit_tot[j])/a
    plt.plot(list_h_norm_tot, list_dP_integrated_over_mode_tot[j],'.-',label = 'mode = %i' %(mode))
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1)
# plt.yscale('log')
# plt.xlim([0,3.1])
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.tight_layout()
os.chdir(path_data_created)
plt.savefig(label_png2,bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)
plt.show()






