
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila
plot the EELS of a rectangular waveguide
integrated over the energy of a mode
and over z (total EELS)
as a function of theta/V0
"""
import numpy as np
import os
import matplotlib.pyplot as plt

## fit with a lorentzian EELS vs energy to find the width of a mode and integrate over energy later 

path_basic = os.getcwd()
path_data = os.path.join(path_basic, 'bem_files_EELS')
path_data_created = os.path.join(path_basic, 'EELS_integrated_over_z')

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

print('1-Define the parameters')

############## BEM ############################################################## 
w = 500
Ee_electron_keV = 200
Ee_electron = Ee_electron_keV*1e3
position_of_the_electron = [0,0.4,0.6]
ind_xe=0
xe = position_of_the_electron[ind_xe]
N=500
##### potential from c++ code  ##################################################
wp_norm_w = 0.5 ## w'/w
d_norm_w = 0.2  ## d/w distance between the side wires normalized to w
h_norm_w = 0.4   ## aspect ratio height/w
Nc = 400   ## discretization points.
epsilon=9 ## permittivity = 9
x0 = xe    ## V(xe,z) with z from z0 to z1
z0 = h_norm_w*1.001    ## z0/w ## start outside the waveguide
z1 = 2          ## z1/w 
nz = 200
# x=0 is in the middle of the center waveguide
# z=0 is at the interface
                            
label_Ee = '_Ee%ikeV' %(Ee_electron_keV)
label_txt = '_wp%.2f_h%.2f_d%.2f_xe%.1f' %(wp_norm_w,h_norm_w,d_norm_w,xe) ## all normalize to width w
if xe == 0.4:  
    labelxe ='_xe40W'
    xe_real = xe*w
elif xe == 0.8:
    labelxe ='_xe80W'
    xe_real = xe*w
else:
    labelxe ='_xe0'
    xe_real = 0
 
# mode = 1
# index_mode = -mode ## if mode = 1, is the highest one because is the last element of the list (the array is sorted from min to max)
 
# if dd == 5: 
#     list_mode = [1,2,3]
# else:
#     list_mode = [1,2]

# if a == 300: 
#     bmin_vals0 = 0.2
#     bmin_vals0 = 0.1
# elif a == 400: 
#     bmin_vals0 = 0.15    ## same bmin = 60 nm
 
# if bmin_vals0 == 0.2:
#     list_mode = [1, 2]
# else:
#     list_mode = [1, 2, 3]    



mode0 = 1
bmin_vals0 = 50/w
list_bmin = [0.1, 0.15, 0.2]
list_bmin = [bmin_vals0]
list_mode = [1,2] 

#%%

print('2-Plot total EELS from the code plot_EELS_integrated_1mode_vs_theta.py for different modes')

labelx=r'$\theta$ (mrad)'
labely='Total EELS' 
title1 = r'$W_{\text{p}}/W$ = %.2f, $h/W$ = %.2f, $d/W$ = %.2f, x/W = %.1f, $E_{\rm e}$ = %i keV' % (wp_norm_w,h_norm_w,d_norm_w,x0,Ee_electron_keV)
title2 = r'$b_{\text{min}} = %i$ nm' %(bmin_vals0*w)
# fig, ax1 = plt.subplots(figsize=tamfig)
# #ax1.title(title1,fontsize=tamtitle)
# ax1.set_xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
# ax1.set_ylabel(labely,fontsize=tamletra,labelpad =labelpady)
# ax1.plot(list_theta, list_Pvalue ,'.-' ,label = r'$b_{\text{min}}/W = %.2f$' %(bmin_vals0))


plt.figure(figsize=tamfig)
plt.title(title1 + ', ' + title2,fontsize=tamtitle)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
os.chdir(path_data_created)
for mode in list_mode: 
    info_label = '_mode%i_w%inm' %(mode,w)  + label_Ee + label_txt + '_bmin%.2f' %(bmin_vals0)
    
    table_P_integrated_over_energy = np.loadtxt('P_integrated_over_z_over' + info_label + '.txt', delimiter='\t', skiprows=1, encoding=None)
    table_P_integrated_over_energy2 = np.transpose(table_P_integrated_over_energy)
    listV0 = table_P_integrated_over_energy2[0]
    list_theta = table_P_integrated_over_energy2[1]
    list_Pvalue = table_P_integrated_over_energy2[2]


    plt.plot(list_theta, list_Pvalue ,'.-' ,label = r'$n = %i$' %(mode))
# ax1.plot(x_filtered, y_filtered ,'.')

# Create second x-axis sharing the same y-axis
# ax2 = ax1.twiny()
# Set ticks for the top axis
# top_ticks =  np.arange(0.1,1.1,0.1)
# ax2.set_xticks(top_ticks)
# ax2.set_xlabel(r"$V_0$ (eV)",fontsize=tamletra,labelpad =labelpadx)
# ax2.set_xlim(listV0[0], listV0[-1])

# ax1.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
# plt.plot(peak_x, peak_y, "x")
# plt.plot(x_left_peak_value, np.ones(len(x_left_peak_value))*0, "x",color = 'blue')
# plt.plot(x_right_peak_value, np.ones(len(x_right_peak_value))*0, "x",color = 'red')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1)

plt.savefig( 'totalEELS_integrated_vs_modes' + '_w%inm' %(w) + label_Ee + label_txt + '_bmin%.2f' %(bmin_vals0) + '.png',bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)


#%%
    
print('3-Plot total EELS from the code plot_EELS_integrated_1mode_vs_theta.py for different bmin')

title3 = r'$n = %i$' %(mode0)
# fig, ax1 = plt.subplots(figsize=tamfig)
#ax1.title(title1,fontsize=tamtitle)
# ax1.set_xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
# ax1.set_ylabel(labely,fontsize=tamletra,labelpad =labelpady)
# ax1.plot(list_theta, list_Pvalue ,'.-' ,label = r'$b_{\text{min}}/W = %.2f$' %(bmin_vals0))

plt.figure(figsize=tamfig)
plt.title(title1 + ', ' + title3,fontsize=tamtitle)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
os.chdir(path_data_created)
for bmin in list_bmin: 
    info_label = '_mode%i_w%inm' %(mode0,w)  + label_Ee + label_txt + '_bmin%.2f' %(bmin_vals0)
    table_P_integrated_over_energy = np.loadtxt('P_integrated_over_z_over'  + info_label + '.txt', delimiter='\t', skiprows=1, encoding=None)
    table_P_integrated_over_energy2 = np.transpose(table_P_integrated_over_energy)
    listV0 = table_P_integrated_over_energy2[0]
    list_theta = table_P_integrated_over_energy2[1]
    list_Pvalue = table_P_integrated_over_energy2[2]


    plt.plot(list_theta, list_Pvalue ,'.-' ,label =  r'$b_{\text{min}} = %i$ nm' %(bmin*w))
# ax1.plot(x_filtered, y_filtered ,'.')

# Create second x-axis sharing the same y-axis
# ax2 = ax1.twiny()
# Set ticks for the top axis
# top_ticks =  np.arange(0.1,1.1,0.1)
# ax2.set_xticks(top_ticks)
# ax2.set_xlabel(r"$V_0$ (eV)",fontsize=tamletra,labelpad =labelpadx)
# ax2.set_xlim(listV0[0], listV0[-1])

# ax1.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
# plt.plot(peak_x, peak_y, "x")
# plt.plot(x_left_peak_value, np.ones(len(x_left_peak_value))*0, "x",color = 'blue')
# plt.plot(x_right_peak_value, np.ones(len(x_right_peak_value))*0, "x",color = 'red')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1)

plt.savefig('totalEELS_integrated_vs_bmin' + '_w%inm' %(w) + label_Ee + label_txt + '_bmin%.2f' %(bmin_vals0) + '.png',bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)

#%%

# from scipy.ndimage import gaussian_filter1d

# def exp_decay(x, a, b, c):
#     return a * np.exp(-b * x) + c

# listx_to_fit = np.array(listy_sorted[0:ind_max])

# list_Pvalue = np.array(list_Pvalue)
# #### removing the noise of the data by detecting sharp changes in the derivative ######
# dy = np.gradient(list_Pvalue, listx_to_fit)
# slope_std = np.std(dy)
# mask = np.abs(dy) < 0.1 * slope_std  # Keep where slope is reasonable
# # Filter: keep points where residual is within 2 standard deviations
# x_filtered = listx_to_fit[mask.astype(bool)]
# y_filtered = list_Pvalue[mask.astype(bool)]

# # Apply Savitzky-Golay filter
# y_smooth = savgol_filter(list_Pvalue, window_length=31, polyorder=3)
# y_gauss = gaussian_filter1d(list_Pvalue, sigma=3)

# # Initial parameter guesses: [a, b, offset]
# x1, x2 = listx_to_fit[0], listx_to_fit[10]
# y1, y2 = y_gauss[0], y_gauss[10]
# b0 = np.log(y1 / y2) / (x2 - x1)
# initial_guess = [max(y_gauss), b0,  min(y_gauss)]

# # Perform the fit
# popt, pcov = curve_fit(exp_decay, listy_sorted[0:ind_max], y_gauss, p0=initial_guess)
# a_fit, b_fit, offset_fit = popt
# y_fit = exp_decay(y_gauss, *popt)




