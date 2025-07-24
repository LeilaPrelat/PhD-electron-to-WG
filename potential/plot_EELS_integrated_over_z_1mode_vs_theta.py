
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila
calculate the EELS of a rectangular waveguide
integrated over the energy of a mode
and over z (total EELS)
as a function of theta/V0
"""
import numpy as np
import os
import matplotlib.pyplot as plt
from EELS_integrated import P_integrated_over_z, run_e_out_EELS, find_width_of_peak
# import concurrent.futures
import matplotlib as mpl
from scipy.interpolate import interp1d

## fit with a lorentzian EELS vs energy to find the width of a mode and integrate over energy later 

path_basic = os.getcwd()
path_data = os.path.join(path_basic, 'bem2D_files_along_z')
path_data_created = os.path.join(path_basic, 'EELS_integrated_over_z')
path_data_zmin = os.path.join(path_basic, 'zmin_solutions')

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
h_norm_w = 0.6   ## aspect ratio height/w
Nc = 400   ## discretization points.
epsilon=9 ## permittivity = 9
x0 = xe    ## V(xe,z) with z from z0 to z1
z0 = h_norm_w   ## z0/w ## start outside the waveguide
z1 = 2          ## z1/w 
nz = 200
# x=0 is in the middle of the center waveguide
# z=0 is at the interface
                             
list_z_norm_w, listV_normU = run_e_out_EELS(wp_norm_w,d_norm_w,h_norm_w,Nc,epsilon,x0,z0,z1) ## import z/W and V(z)/V0 from c++ code
zz, V0_over_U_list = run_e_out_EELS(wp_norm_w,d_norm_w,h_norm_w,Nc,epsilon,0,h_norm_w,h_norm_w)
V0_over_U = np.abs(V0_over_U_list[0])

bmin_vals0 = 50/w
z_min_val_norm_w = bmin_vals0 + h_norm_w
###################################################################################
h=h_norm_w*w
# mode = 1
# index_mode = -mode ## if mode = 1, is the highest one because is the last element of the list (the array is sorted from min to max)
plot_figure = 1    ## plot the lorentzian fitting
step2 = 0.001      ## thinner range to integrate over energy

list_mode = [1]
# list_mode = [3]

if xe == 0.4:  
    labelxe ='_xe40W'
    xe_real = xe*w
elif xe == 0.8:
    labelxe ='_xe80W'
    xe_real = xe*w
else:
    labelxe ='_xe0'
    xe_real = 0
    
label_txt = '_wp%.2f_h%.2f_d%.2f_xe%.1f' %(wp_norm_w,h_norm_w,d_norm_w,xe) ## all normalize to width w
label_Ee = '_Ee%ikeV' %(Ee_electron_keV) 
 
ze=50
name_list_energy = 'energy_list_N%i_a%inm_h%inm_ze%inm_Ee%ikeV%s'%(N,w,h,ze,Ee_electron_keV,labelxe) + '.txt'
os.chdir(path_data)
tabla_energy = np.loadtxt(name_list_energy, delimiter=' ',dtype=None)
tabla_energy_2 = np.transpose(tabla_energy)
list_energy0 = tabla_energy_2
# list_energy0 = tabla_energy_2[0:-1]

#%%

print('2-Load data from the code plot_bmin_vs_V0_theta.py of theta,V0 for a fix bmin')

os.chdir(path_data_zmin)
bmin_vals = np.loadtxt('bmin' + label_Ee + label_txt + '.txt', delimiter='\t', skiprows = 1)
V0_vals = np.loadtxt('V0_for_bmin%.2f' %(bmin_vals0) + label_Ee + label_txt + '.txt', delimiter='\t', skiprows = 1)
theta_mrad_vals = np.loadtxt('theta_mrad_for_bmin%.2f' %(bmin_vals0) + label_Ee + label_txt + '.txt', delimiter='\t', skiprows = 1)
## solutions (V0,theta) such that bmin = bmin_vals0
listx_sorted = np.loadtxt('V0_sol_for_bmin%.2f'%(bmin_vals0)  + label_Ee + label_txt + '.txt', delimiter='\t', skiprows = 1)
listy_sorted = np.loadtxt('theta_sol_mrad_for_bmin%.2f'%(bmin_vals0) + label_Ee + label_txt + '.txt', delimiter='\t', skiprows = 1)
 
tamfig2 = [5, 3]

print('3-Plot b_min as a function of (V0,theta) and see if they match with the contorn plot')

title = r'$W_{\text{p}}/W$ = %.2f, $h/W$ = %.2f, $d/W$ = %.2f, x/W = %.1f, $E_{\rm e}$ = %i keV' % (wp_norm_w,h_norm_w,d_norm_w,x0,Ee_electron_keV)
# another cbar with constant a*1e-3 (real units, in microns)
cte_cbar2 = w*1e-3
cmap = plt.cm.RdBu   # define the colormap
bounds1 =   np.logspace(np.log10(0.05), np.log10(100) , 10) 
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
# Ntot = 3
#%%

print('4-Plot EELS vs energy to fit it and find the width for all (V0,theta)')
#print('IMPORTANT: we are assuming the position of the peak DOES vary with V0,theta')

for mode in list_mode: 
    index_mode = -mode  ## if mode = 1, is the highest one because is the last element of the list (the array is sorted from min to max)
    
    list_Pvalue = []
    header0 = 'energy (eV)    P(omega), '
    
    for j in range(Ntot): 

        V0 = listx_sorted[j]
        theta_mrad = listy_sorted[j]
        theta = theta_mrad*1e-3
        info_bmin = '_V0_%.2feV_theta%.2fmrad_bmin%.2f' %(V0,theta_mrad,bmin_vals0) ## info from contour plot from plot_bmin_vs_V0_theta.py
        
        all_info_label = '_mode%i_w%inm' %(mode,w) + label_Ee + label_txt + info_bmin
        
        os.chdir(path_data)
        list_P_integrated_over_z = []   ## we integrate EELS over z for every theta, V0 and we fit each of those plots with Lorenztian
        
        k=0
        list_energy0_EELS_positive = []
        for energy in list_energy0:
            theta = theta_mrad*1e-3
            
            try:                                                                               
                P_int_value,list_z_norm_a_EELS,P_array = P_integrated_over_z(z_min_val_norm_w, energy, w, h, xe, Ee_electron_keV, N, theta, V0, list_z_norm_w, listV_normU, V0_over_U)
                # print(list_z_norm_a_EELS, P_array)
                list_P_integrated_over_z.append(P_int_value)
                list_energy0_EELS_positive.append(energy)
               # print(k,P_array)
                k=k+1
                # plt.figure()
                # plt.plot(list_z_norm_a_EELS,EELS_array)
            except ValueError: ## some EELS are negative, we omit them 
                pass
            
        table_P_integrated_over_z = np.transpose([list_energy0_EELS_positive, list_P_integrated_over_z])
        os.chdir(path_data_created)
        np.savetxt('P_integrated_over_z' + all_info_label + '.txt', table_P_integrated_over_z, fmt='%.10f', delimiter='\t', header = header0 + header, encoding=None)
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
    
            ######################## OLD [NOT NEEDED] #########################################################################################################
                
        # Find the index of the closest value
        # idx_x_left = (np.abs(list_energy0 - x_left)).argmin()
        # closest_x_left = list_energy0[idx_x_left]
        # idx_x_right = (np.abs(list_energy0 - x_right)).argmin()
        # closest_x_right = list_energy0[idx_x_right]
        
        # def Pintegrated_over_energy(V0,theta_mrad,list_energy):
        #     theta = theta_mrad*1e-3
             
        #     P_integrated_over_energy = 0
        #     for energy in list_energy:
        #         P_int_value = P_integrated_over_z(z_min_val, theta, V0, Ee_electron, bb, ss, dd, energy, a, N, list_z_norm_a, listV_normV0)
         
        #         P_integrated_over_energy = P_int_value + P_integrated_over_energy
            
        #     delta_energy = list_energy[1] - list_energy[0]
         
        #     return P_integrated_over_energy*delta_energy
        
 
        # P_example = Pintegrated_over_energy(V0,theta_mrad,list_energy_over_mode)
    
    #    print("Difference between using the Lorenztian to integrate over energy and using the data:",np.abs(y_integrate2-P_example))
        
        # os.chdir(path_data)
        # table_P_integrated_over_z = np.loadtxt('P_integrated_over_z' + label_Ee + all_info_label, delimiter='\t', skiprows=1)
        # table_P_integrated_over_z2 = np.transpose(table_P_integrated_over_z)
        # list_energy0 = table_P_integrated_over_z2[0]
        # list_P_integrated_over_z = table_P_integrated_over_z2[1]
    
        ######################## OLD #########################################################################################################
    
    print('6-Integrate over energy from x_left_peak, x_right_peak as a function of (theta,V0) for a fixed b_min = %.2f' %(bmin_vals0))
    
    info_label = '_mode%i_w%inm' %(mode,w) + label_Ee + label_txt + '_bmin%.2f' %(bmin_vals0)
    
    ind_max = len(list_Pvalue)
    table_P_integrated_over_energy = np.transpose([listx_sorted[0:ind_max],listy_sorted[0:ind_max], list_Pvalue])
    header1 = 'V0 (eV)    theta (mrad)    P_mode1, '
    np.savetxt('P_integrated_over_z_over' + info_label + '.txt', table_P_integrated_over_energy, fmt='%.10f', delimiter='\t', header = header1 + header, encoding=None)


#%%
