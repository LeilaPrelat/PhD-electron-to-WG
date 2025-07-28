
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
from fit_peaks import find_width_of_peak
# import concurrent.futures
import matplotlib as mpl
from scipy.interpolate import interp1d

## fit with a lorentzian EELS vs energy to find the width of a mode and integrate over energy later 

path_basic = os.getcwd()
path_data_EELS = os.path.join(path_basic, 'datfiles_EELS_along_z/new_files')
path_data_EELS_integrated = os.path.join(path_basic, 'EELS_integrated')
path_data_zmin = os.path.join(path_basic, 'potential_data')

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
me_c2_eV = 510998.95069  ## me*c**2 in eV

#%%
print('1-Define the parameters')

############## BEM ############################################################## 
w = 400
Ee_electron_keV = 200
Ee_electron = Ee_electron_keV*1e3
position_of_the_electron = [0,0.4,0.6]
ind_xe=0
xe_norm_w = position_of_the_electron[ind_xe]
N=500
##### potential from c++ code  ##################################################
wp_norm_w = 0.5 ## w'/w
d_norm_w = 0.2  ## d/w distance between the side wires normalized to w
h_norm_w = 0.5   ## aspect ratio height/w
Nc = 400   ## discretization points.
x0 = xe_norm_w    ## V(xe,z) with z from z0 to z1
 
# x=0 is in the middle of the center waveguide
# z=0 is at the interface
                             
bmin_norm_w = 50/w ## 80/w ## 100/w
z_min_val_norm_w = bmin_norm_w + h_norm_w
###################################################################################
h=h_norm_w*w
# mode = 1
# index_mode = -mode ## if mode = 1, is the highest one because is the last element of the list (the array is sorted from min to max)
plot_figure = 0    ## plot the lorentzian fitting
step2 = 0.001      ## thinner range to integrate over energy

list_mode = [1]
# list_mode = [3]

if xe_norm_w == 0.4:  
    labelxe ='_xe40W'
    xe_real = xe_norm_w*w
elif xe_norm_w == 0.8:
    labelxe ='_xe80W'
    xe_real = xe_norm_w*w
else:
    labelxe ='_xe0'
    xe_real = 0
    
label_txt = '_wp%.2f_h%.2f_d%.2f_xe%.1f' %(wp_norm_w,h_norm_w,d_norm_w,xe_norm_w) ## all normalize to width w
label_Ee = '_Ee%ikeV' %(Ee_electron_keV) 
 
# list_energy0 = tabla_energy_2[0:-1]

#%%

print('2-Load data from the code plot_bmin_vs_V0_theta.py of theta,V0 for a fix bmin')

os.chdir(path_data_zmin)
## loading data for plot from find_bmin_vs_V0_theta.py
bmin_vals = np.loadtxt('bmin_for_plot' + label_Ee + label_txt + '.txt', delimiter='\t', skiprows = 1)
V0_vals = np.loadtxt('V0_for_plot_for_bmin%.2f' %(bmin_norm_w) + label_Ee + label_txt + '.txt', delimiter='\t', skiprows = 1)
theta_mrad_vals = np.loadtxt('theta_mrad_for_plot_for_bmin%.2f' %(bmin_norm_w) + label_Ee + label_txt + '.txt', delimiter='\t', skiprows = 1)

## loading solutions (V0,theta) such that bmin = bmin_vals0 from find_bmin_vs_V0_theta.py
name_sol = 'V0_theta_mrad_sols_for_bmin%.2f' %(bmin_norm_w) + label_Ee + label_txt + '.txt'
table_solutions = np.loadtxt(name_sol, delimiter='\t', skiprows=1, encoding=None)
table_solutions2 = np.transpose(table_solutions)
listx_sorted, listy_sorted = table_solutions2

tamfig2 = [5, 3]

print('3-Plot b_min as a function of (V0,theta) and see if they match with the contorn plot')

title = r'$W_{\text{p}}/W$ = %.2f, $h/W$ = %.2f, $d/W$ = %.2f, x/W = %.1f, $E_{\rm e}$ = %i keV' % (wp_norm_w,h_norm_w,d_norm_w,x0,Ee_electron_keV)
# another cbar with constant a*1e-3 (real units, in microns)
cte_cbar2 = w*1e-3
vmin1 , vmax1 = np.nanmin(bmin_vals), np.nanmax(bmin_vals)
cmap = plt.cm.RdBu   # define the colormap
bounds1 =   np.logspace(np.log10(0.05), np.log10(vmax1) , 10) 
# bounds1 =   np.logspace(np.log10(vmin1), np.log10(100) , 10) 
norm1 = mpl.colors.BoundaryNorm(bounds1, cmap.N)
norm2 = mpl.colors.BoundaryNorm(bounds1*cte_cbar2, cmap.N) ## second colorbar with real units 


limits1 = [np.min(V0_vals) , np.max(V0_vals),np.min(theta_mrad_vals) , np.max(theta_mrad_vals)]
plt.figure(figsize=tamfig2)
#plt.title(title1 + r', $E_{\text{e}}$ = %i keV' %(Ee_electron_keV),fontsize=tamtitle)
im_show = plt.imshow(bmin_vals, extent = limits1, cmap=cmap, aspect='auto', interpolation = 'bicubic',origin = 'lower' , norm = norm1  ) 
    #plt.clabel(contours, fmt='%.2f', colors='green', fontsize=tamletra, manual=[(0.5, 5) ])  # Label contours
cbar = plt.colorbar(im_show, fraction=0.046, pad=0.18 , format = '%.2f') 
im_show2 = plt.imshow(np.array(bmin_vals)*cte_cbar2, extent = limits1, cmap=cmap, aspect='auto', interpolation = 'bicubic',origin = 'lower' ,norm=norm2  )  ## second colorbar with real units 
cbar2 = plt.colorbar(im_show2, fraction=0.046, pad=0.04, orientation = 'vertical')
cbar.ax.set_title(r'$b_{\text{min}}/W$',fontsize=tamletra-1)
cbar.ax.tick_params(labelsize = tamnum-2, width=0.1, direction="in",which = 'both', length = 2,pad = pad)
cbar2.ax.tick_params(labelsize = tamnum-2, width=0.1, direction="in",which = 'both', length = 2,pad = pad)
cbar2.ax.set_title(r'$b_{\text{min}}$ ($\mu$m)',fontsize=tamletra-1)
plt.xlabel(r'$V_0$ (eV)',fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(r'$\theta$ (mrad)',fontsize=tamletra,labelpad =labelpady)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.plot(listx_sorted,listy_sorted,'--',color = 'green')
plt.savefig('zmin_aux' + label_Ee + label_txt + '.png', format='png',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
plt.show()

header = title
Ntot = len(listx_sorted)
 
#%%

print('4-Plot EELS vs energy to fit it and find the width for all (V0,theta)')
#print('IMPORTANT: we are assuming the position of the peak DOES vary with V0,theta')

for mode in list_mode: 
    index_mode = -mode  ## if mode = 1, is the highest one because is the last element of the list (the array is sorted from min to max)
    
    list_Pvalue = []
    header0 = 'energy (eV)    P(omega), '
    
    for j in range(Ntot): ## Ntot
        V0 = listx_sorted[j]
        theta_mrad = listy_sorted[j]
        theta = theta_mrad*1e-3
        info_bmin = '_V0_%.5feV_theta%.5fmrad_bmin%.3f' %(V0,theta_mrad,bmin_norm_w) ## info from contour plot from plot_bmin_vs_V0_theta.py
        
        all_info_label = '_mode%i_w%inm' %(mode,w) + label_Ee + label_txt + info_bmin
        
        os.chdir(path_data_EELS_integrated)
        name_P_integrated_over_z = 'P_integrated_over_z_w%inm_Ee%ikeV_wp%.2f_h%.2f_d%.2f_xe%.2f' %(w,Ee_electron_keV,wp_norm_w,h_norm_w,d_norm_w,xe_norm_w) + info_bmin + '.dat'
 
        tabla_P = np.loadtxt(name_P_integrated_over_z, delimiter=' ',dtype=None,skiprows=1)
        tabla_P_2 = np.transpose(tabla_P)
        list_energy0_EELS_positive,list_P_integrated_over_z = tabla_P_2
        
        
        if j == 0:
            print('5-Find the dw of the mode = %i by fitting it with a Lorenztian' %(mode))
        

        ######### find peaks and sort them by the highest to minimum (then we identify each mode according to its amplitude) by fitting the curve with a lorentzian
        x_left, x_right, x0_fit, function_lorentzian, fit_success = find_width_of_peak(list_energy0_EELS_positive,list_P_integrated_over_z,index_mode,title,plot_figure)
        # Integration by summation (small steps)
        x_integrate = np.arange(x_left, x_right + step2, step2)
        
        if fit_success:
            ############### interpolation of P(omega) after integration over z ###############
            EELS_vs_energy = interp1d(list_energy0_EELS_positive, list_P_integrated_over_z)
            integration_over_energy = np.sum(EELS_vs_energy(x_integrate))*(step2)
            
            np.savetxt('list_energy_for_P_integration'  + all_info_label + '.txt', x_integrate, fmt='%.10f', delimiter='\t', header = header, encoding=None)
            list_Pvalue.append(function_lorentzian) ## integration using the lorentzian fit 
            
            print(j,V0,theta_mrad,integration_over_energy,function_lorentzian)
    
    
    print('6-Integrate over energy from x_left_peak, x_right_peak as a function of (theta,V0) for a fixed b_min = %.2f' %(bmin_norm_w))
    
    info_label = '_mode%i_w%inm' %(mode,w) + label_Ee + label_txt + '_bmin%.3f' %(bmin_norm_w)
    
    ind_max = len(list_Pvalue)
    table_P_integrated_over_energy = np.transpose([listx_sorted[0:ind_max],listy_sorted[0:ind_max], list_Pvalue])
    header1 = 'V0 (eV)    theta (mrad)    P_mode1, '
    np.savetxt('P_integrated_over_z_over' + info_label + '.txt', table_P_integrated_over_energy, fmt='%.10f', delimiter='\t', header = header1 + header, encoding=None)


#%%