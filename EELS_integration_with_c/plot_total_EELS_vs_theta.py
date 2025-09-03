
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
from scipy.interpolate import make_interp_spline
# import concurrent.futures

## fit with a lorentzian EELS vs energy to find the width of a mode and integrate over energy later 

path_basic = os.getcwd()
path_data_EELS_integrated = os.path.join(path_basic, 'EELS_integrated')

compare_with_approx = 0 ## compare the results with the approximation wtih V(z)/V0 = a*z/W + b

# Function to smooth data
def smooth_curve(x, y, points=200):
    x_new = np.linspace(x.min(), x.max(), points)
    spline = make_interp_spline(x, y, k=3)  # Cubic spline
    y_smooth = spline(x_new)
    return x_new, y_smooth

def moving_average(y, window=5):
    return np.convolve(y, np.ones(window)/window, mode="same")

from scipy.signal import savgol_filter

#%%

tamfig = [4,3]
tamfig2 = [4.3,3.5]
tamletra = 15.25
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
# mode = 1
# index_mode = -mode ## if mode = 1, is the highest one because is the last element of the list (the array is sorted from min to max)
plot_figure = 1    ## plot the lorentzian fitting
step2 = 0.001      ## thinner range to integrate over energy

list_mode = [1,3]
# list_mode = [3]
   
label_txt = '_wp%.2f_h%.2f_d%.2f_xe%.1f' %(wp_norm_w,h_norm_w,d_norm_w,xe_norm_w) ## all normalize to width w
label_Ee = '_Ee%ikeV' %(Ee_electron_keV) 

# Ntot = 3mode0 = 1
bmin_vals0 = 50/w
mode0=1
#list_bmin = [0.1, 0.15, 0.2]
list_bmin = [50/w, 80/w, 100/w]
#list_bmin = [50/w, 80/w, 100/w]
list_mode = [1,2]
list_mode = [1,3]

#%%
os.chdir(path_data_EELS_integrated)
print('2-Plot total EELS from the code plot_EELS_integrated_1mode_vs_theta.py for different modes')

labelx=r'$\theta$ (mrad)'
labely='Number of photons' 
#labely=r'$\Gamma$ (s)'

title1 = r'$W_{\text{p}}/W$ = %.2f, $h/W$ = %.2f, $d/W$ = %.2f, x/W = %.1f, $E_{\rm e}$ = %i keV' % (wp_norm_w,h_norm_w,d_norm_w,x0,Ee_electron_keV)
title2 = r'$b_{\text{min}} = %i$ nm' %(bmin_vals0*w)
# fig, ax1 = plt.subplots(figsize=tamfig)
# #ax1.title(title1,fontsize=tamtitle)
# ax1.set_xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
# ax1.set_ylabel(labely,fontsize=tamletra,labelpad =labelpady)
# ax1.plot(list_theta, list_Pvalue ,'.-' ,label = r'$b_{\text{min}}/W = %.2f$' %(bmin_vals0))


fig, ax1 = plt.subplots(figsize=tamfig2)
#plt.title(title1 + ', ' + title2,fontsize=tamtitle)
plt.title('',fontsize=tamtitle)
ax1.set_ylabel(labely,fontsize=tamletra,labelpad =labelpady)
angle_ticks = np.arange(0.1,2.2,0.4)
V0_ticks = np.arange(0.2,1.2,0.2)

# Create a second x-axis sharing the same y
ax1.set_xticks(angle_ticks)
ax1.set_xlabel(r'$\theta$ (mrad)',fontsize=tamletra,labelpad =labelpadx)
plt.yticks(np.arange(5,35,5))

for mode in list_mode: 
    info_label = '_mode%i_w%inm' %(mode,w)  + label_Ee + label_txt + '_bmin%.2f' %(bmin_vals0)
    
    table_P_integrated_over_energy = np.loadtxt('P_integrated_over_z_over' + info_label + '.txt', delimiter='\t', skiprows=1, encoding=None)
    table_P_integrated_over_energy2 = np.transpose(table_P_integrated_over_energy)
    listV0 = table_P_integrated_over_energy2[0]
    list_theta = table_P_integrated_over_energy2[1]
    list_Pvalue = table_P_integrated_over_energy2[2]

    # Smooth each dataset
    #list_theta_smooth, list_Pvalue_smooth = smooth_curve(list_theta, list_Pvalue)
    
    # Apply Savitzky–Golay filter
    photons_smooth = savgol_filter(list_Pvalue, window_length=9, polyorder=3)
    photons_smooth2 = savgol_filter(photons_smooth, window_length=4, polyorder=3)
 
    if compare_with_approx == 0: 
        plt.plot(list_theta, photons_smooth ,'.-' , label =  r'$n = %i$' %(mode))
    else:
        plt.plot(list_theta, photons_smooth ,'.-' )
    
    if compare_with_approx == 1: 
        table_P_integrated_over_energy_approx = np.loadtxt('P_integrated_over_z_over_approx' + info_label + '.txt', delimiter='\t', skiprows=1, encoding=None)
        table_P_integrated_over_energy2_approx = np.transpose(table_P_integrated_over_energy_approx)
        listV0_approx = table_P_integrated_over_energy2_approx[0]
        list_theta_approx = table_P_integrated_over_energy2_approx[1]
        list_Pvalue_approx = table_P_integrated_over_energy2_approx[2]
    
        list_Pvalue_approx_smooth = smooth_curve(list_theta_approx, list_Pvalue_approx)
        photons_approx_smooth = savgol_filter(list_Pvalue_approx, window_length=9, polyorder=3)
        
        plt.plot(list_theta_approx, photons_approx_smooth ,'-' )
        index = np.argmin(np.abs(list_theta_approx - 0.25858))
        ## 0.25858 used to compare both EELS
        
        
        plt.plot(list_theta_approx[index],photons_approx_smooth[index],'o',color = 'black')
        
# ax1.plot(x_filtered, y_filtered ,'.')

# Create second x-axis sharing the same y-axis
# ax2 = ax1.twiny()
# Set ticks for the top axis
# top_ticks =  np.arange(0.1,1.1,0.1)
# ax2.set_xticks(top_ticks)
# ax2.set_xlabel(r"$V_0$ (eV)",fontsize=tamletra,labelpad =labelpadx)
# ax2.set_xlim(listV0[0], listV0[-1])

# Move ax2 to the bottom
ax2 = ax1.twiny()
# ax2.spines["bottom"].set_position(("outward", 40))
# ax2.xaxis.set_ticks_position("bottom")
# ax2.xaxis.set_label_position("bottom")
ax2.spines["bottom"].set_visible(False)

# Secondary axis data and ticks
ax2.set_xticks(V0_ticks)
ax2.set_xlabel(r'$V_0$ (eV)',fontsize=tamletra,labelpad =labelpadx)


ax1.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
ax2.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
#plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
# plt.plot(peak_x, peak_y, "x")
# plt.plot(x_left_peak_value, np.ones(len(x_left_peak_value))*0, "x",color = 'blue')
# plt.plot(x_right_peak_value, np.ones(len(x_right_peak_value))*0, "x",color = 'red')
ax1.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1)

plt.savefig( 'totalEELS_integrated_vs_modes' + '_w%inm' %(w) + label_Ee + label_txt + '_bmin%.2f' %(bmin_vals0) + '.png',bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)
plt.savefig( 'totalEELS_integrated_vs_modes' + '_w%inm' %(w) + label_Ee + label_txt + '_bmin%.2f' %(bmin_vals0) + '.pdf',bbox_inches='tight',pad_inches = 0.01, format='pdf', dpi=dpi)

#%%
    
print('3-Plot total EELS from the code plot_EELS_integrated_1mode_vs_theta.py for different bmin')

title3 = r'$n = %i$' %(mode0)
# fig, ax1 = plt.subplots(figsize=tamfig)
#ax1.title(title1,fontsize=tamtitle)
# ax1.set_xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
# ax1.set_ylabel(labely,fontsize=tamletra,labelpad =labelpady)
# ax1.plot(list_theta, list_Pvalue ,'.-' ,label = r'$b_{\text{min}}/W = %.2f$' %(bmin_vals0))


fig, ax1 = plt.subplots(figsize=tamfig2)
#plt.title(title1 + ', ' + title2,fontsize=tamtitle)
plt.title('',fontsize=tamtitle)
ax1.set_ylabel(labely,fontsize=tamletra,labelpad =labelpady)
angle_ticks = np.arange(0.1,2.2,0.4)
V0_ticks = np.arange(0.2,1.2,0.2)

# Create a second x-axis sharing the same y
ax1.set_xticks(angle_ticks)
ax1.set_xlabel(r'$\theta$ (mrad)',fontsize=tamletra,labelpad =labelpadx)
plt.yticks(np.arange(5,35,5))
for bmin in list_bmin: 
    info_label = '_mode%i_w%inm' %(mode0,w)  + label_Ee + label_txt + '_bmin%.2f' %(bmin)
    table_P_integrated_over_energy = np.loadtxt('P_integrated_over_z_over'  + info_label + '.txt', delimiter='\t', skiprows=1, encoding=None)
    table_P_integrated_over_energy2 = np.transpose(table_P_integrated_over_energy)
    listV0 = table_P_integrated_over_energy2[0]
    list_theta = table_P_integrated_over_energy2[1]
    list_Pvalue = table_P_integrated_over_energy2[2]

        
    # Apply Savitzky–Golay filter
    photons_smooth = savgol_filter(list_Pvalue, window_length=9, polyorder=3)
    photons_smooth2 = savgol_filter(photons_smooth, window_length=4, polyorder=3)
    
    plt.plot(list_theta, photons_smooth ,'.-' ,label =  r'$b_{\text{min}} = %i$ nm' %(bmin*w))
    
    
    if compare_with_approx == 1: 
        table_P_integrated_over_energy_approx = np.loadtxt('P_integrated_over_z_over_approx' + info_label + '.txt', delimiter='\t', skiprows=1, encoding=None)
        table_P_integrated_over_energy2_approx = np.transpose(table_P_integrated_over_energy_approx)
        listV0_approx = table_P_integrated_over_energy2_approx[0]
        list_theta_approx = table_P_integrated_over_energy2_approx[1]
        list_Pvalue_approx = table_P_integrated_over_energy2_approx[2]
    
        plt.plot(list_theta_approx, list_Pvalue_approx ,'-' )
    
    
# Move ax2 to the bottom
ax2 = ax1.twiny()
# ax2.spines["bottom"].set_position(("outward", 40))
# ax2.xaxis.set_ticks_position("bottom")
# ax2.xaxis.set_label_position("bottom")
ax2.spines["bottom"].set_visible(False)

# Secondary axis data and ticks
ax2.set_xticks(V0_ticks)
ax2.set_xlabel(r'$V_0$ (eV)',fontsize=tamletra,labelpad =labelpadx)


ax1.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
ax2.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
#plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
# plt.plot(peak_x, peak_y, "x")
# plt.plot(x_left_peak_value, np.ones(len(x_left_peak_value))*0, "x",color = 'blue')
# plt.plot(x_right_peak_value, np.ones(len(x_right_peak_value))*0, "x",color = 'red')
ax1.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1)

plt.savefig('totalEELS_integrated_vs_bmin' + '_w%inm' %(w) + label_Ee + label_txt + '_bmin%.2f' %(bmin_vals0) + '.png',bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)
plt.savefig('totalEELS_integrated_vs_bmin' + '_w%inm' %(w) + label_Ee + label_txt + '_bmin%.2f' %(bmin_vals0) + '.pdf',bbox_inches='tight',pad_inches = 0.01, format='pdf', dpi=dpi)

