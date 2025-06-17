
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
from EELS_integrated import P_integrated_over_z, dy_cached, find_width_of_peak
# import concurrent.futures
import matplotlib as mpl
from scipy.interpolate import interp1d

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

############ parameters for bem2d ############
a = 400 ## nm
h = 300 ## nm
s = 20  ## nm
############ parameters for c++ ############
bb = h/a ## h/a ## zmin is bb/2 (upper surface of the wg)
ss = s/a
dd = 5   ## distance to the plane
N = 400
Ee_electron_keV = 200
Ee_electron = Ee_electron_keV*1e3
label_Ee = '_Ee%ikeV' %(Ee_electron_keV)

step = 0.005
list_energy0 = np.arange(0.2,1.5 + step,step) ## we will reduce this energy list to the energies around the first peak to integrate over \omega
list_energy0 = np.arange(0.2,3 + step,step)
list_z_norm_a, listV_normV0 = dy_cached(bb,ss,dd,N) ## import z/W and V(z)/V0 from c++ code

# mode = 1
# index_mode = -mode ## if mode = 1, is the highest one because is the last element of the list (the array is sorted from min to max)
plot_figure = 0    ## plot the lorentzian fitting
step2 = 0.001      ## thinner range to integrate over energy

list_mode = [1,2,3]

#%%

print('2-Load data of z_min from the code plot_potential_V_V0_theta_only_above_wg.py')

tamfig2 = [5, 3]
bmin_vals = np.loadtxt('bmin' + label_Ee + '_dd%i_hh%.2f.txt' %(dd,bb),  delimiter='\t', skiprows = 1)
V0_vals = np.loadtxt('V0' + label_Ee + '_dd%i_hh%.2f.txt' %(dd,bb), delimiter='\t', skiprows = 1)
theta_mrad_vals = np.loadtxt('theta_mrad' + label_Ee + '_dd%i_hh%.2f.txt' %(dd,bb), delimiter='\t', skiprows = 1)

# find the index closest to the values we want for bmin 
if a == 300: 
    bmin_vals0 = 0.2
    bmin_vals0 = 0.1
elif a == 400: 
    bmin_vals0 = 0.15    ## same bmin = 60 nm
    bmin_vals0 = 0.1
    bmin_vals0 = 0.2
    
z_min_val = bmin_vals0 + bb/2 
print('3-Find the V0,theta values for a fixed b_min = %.2f' %(bmin_vals0))
print('4-Plot b_min as a function of (V0,theta) and see if they match with the contorn plot')

title1 = r'$W$ = %i nm, $h$ = %i nm, $s$ = %i nm, $d/W = %i$, $E_{\rm e}$ = %i keV' %(a,h,s,dd,Ee_electron_keV)

limits1 = [np.min(V0_vals) , np.max(V0_vals),np.min(theta_mrad_vals) , np.max(theta_mrad_vals)]
cmap = plt.cm.RdBu   # define the colormap

# another cbar with constant a*1e-3 (real units, in microns)
cte_cbar2 = a*1e-3
# Create contour plot
contour_levels = [bmin_vals0]

vmin1 , vmax1 = np.nanmin(bmin_vals), np.nanmax(bmin_vals)
bounds1 =  [vmin1,0.1,bmin_vals0,0.5,1]
bounds1 =   np.logspace(np.log10(0.1), np.log10(100) , 10) 
norm1 = mpl.colors.BoundaryNorm(bounds1, cmap.N)
norm2 = mpl.colors.BoundaryNorm(bounds1*cte_cbar2, cmap.N) ## second colorbar with real units 

plt.figure(figsize=tamfig2)
#plt.title(title1 + r', $E_{\text{e}}$ = %i keV' %(Ee_electron_keV),fontsize=tamtitle)
im_show = plt.imshow(bmin_vals, extent = limits1, cmap=cmap, aspect='auto', interpolation = 'bicubic',origin = 'lower' , norm = norm1  ) 
contours = plt.contour(V0_vals, theta_mrad_vals, bmin_vals, levels=contour_levels, colors='green', linestyles='dashed'  )

# Access the contour data
allsegs = contours.collections[0].get_paths()  # Assuming only one contour collection
list_V0_fix_bmin = []
list_theta_fix_bmin = [] 
for xy  in allsegs[0].vertices:
    x,y = xy
    list_theta_fix_bmin.append(y)
    list_V0_fix_bmin.append(x)
# Subtract the V0 and \theta values of corresponding points on different levels
# Sorte the lists argsort()
index_sorted = np.argsort(list_V0_fix_bmin)
listx_sorted = np.sort(list_V0_fix_bmin)
listy_sorted = [] 
for index in index_sorted:
    listy_sorted.append(list_theta_fix_bmin[index])

#plt.clabel(contours, fmt='%.2f', colors='green', fontsize=tamletra, manual=[(0.5, 5) ])  # Label contours
cbar = plt.colorbar(im_show, fraction=0.046, pad=0.2   ,format = '%.2f') 
im_show2 = plt.imshow(np.array(bmin_vals)*cte_cbar2, extent = limits1, cmap=cmap, aspect='auto', interpolation = 'bicubic',origin = 'lower' ,norm=norm2  )  ## second colorbar with real units 
cbar2 = plt.colorbar(im_show2, fraction=0.046, pad=0.04, orientation = 'vertical')
cbar.ax.set_title(r'$b_{\text{min}}/W$',fontsize=tamletra-1)
cbar.ax.tick_params(labelsize = tamnum-2, width=0.1, direction="in",which = 'both', length = 2,pad = pad)
cbar2.ax.tick_params(labelsize = tamnum-2, width=0.1, direction="in",which = 'both', length = 2,pad = pad)
cbar2.ax.set_title(r'$b_{\text{min}}$ ($\mu$m)',fontsize=tamletra-1)
# plt.xticks(np.arange(0,np.max(V0_vals)+0.5,0.5))
#plt.plot(listx_sorted,listy_sorted,'--',color = 'blue')
plt.xlabel(r'$V_0$ (eV)',fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(r'$\theta$ (mrad)',fontsize=tamletra,labelpad =labelpady)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.savefig('zmin_aux' + label_Ee + '_bmin%.2f_dd%i_hh%.2f.png' %(bmin_vals0,dd,bb), format='png',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
#plt.title('Root x as function of V0 and theta')
plt.show()

header = title1
os.chdir(path_data_created)
np.savetxt('V0_vals_for_bmin%.2f' %(bmin_vals0) + label_Ee + '_dd%i_hh%.2f.txt' %(dd,bb), listx_sorted, fmt='%.10f', delimiter='\t', header = header, encoding=None)
np.savetxt('theta_mrad_vals_for_bmin%.2f' %(bmin_vals0) + label_Ee + '_dd%i_hh%.2f.txt' %(dd,bb), listy_sorted, fmt='%.10f', delimiter='\t', header = header, encoding=None)

Ntot = len(listx_sorted)

#%%

print('5-Plot EELS vs energy to fit it and find the width for all (V0,theta), 6-later use the width to integrate over energy')
#print('IMPORTANT: we are assuming the position of the peak DOES vary with V0,theta')

for mode in list_mode: 
    index_mode = -mode  ## if mode = 1, is the highest one because is the last element of the list (the array is sorted from min to max)
    
    list_Pvalue = []
    header0 = 'energy (eV)    P(omega), '
    
    for j in range(Ntot): 
        V0 = listx_sorted[j]
        theta_mrad = listy_sorted[j]
        theta = theta_mrad*1e-3
        all_info_label = '_mode%i_dd%i_hh%.2f_V0_%.2feV_theta%.2fmrad_bmin%.2f' %(mode,dd,bb,V0,theta_mrad,bmin_vals0)
        os.chdir(path_data)
        list_P_integrated_over_z = []   ## we integrate EELS over z for every theta, V0 and we fit each of those plots with Lorenztian
        for energy in list_energy0:
     
            P_int_value = P_integrated_over_z(z_min_val, theta, V0, Ee_electron, bb, ss, dd, energy, a, N, list_z_norm_a, listV_normV0)
            # print(energy,P_int_value)
            list_P_integrated_over_z.append(P_int_value)
        
        table_P_integrated_over_z = np.transpose([list_energy0, list_P_integrated_over_z])
        os.chdir(path_data_created)
        np.savetxt('P_integrated_over_z' + label_Ee + all_info_label + '.txt', table_P_integrated_over_z, fmt='%.10f', delimiter='\t', header = header0 + header, encoding=None)
        if j == 0:
            print('6-Find the wmin and wmax of the mode = %i by fitting it with a Lorenztian' %(mode))
    
        ######### find peaks and sort them by the highest to minimum (then we identify each mode according to its amplitude) by fitting the curve with a lorentzian
        x_left, x_right = find_width_of_peak(list_energy0,list_P_integrated_over_z,V0,theta_mrad,index_mode,title1,plot_figure)
        # Integration by summation (small steps)
        x_integrate = np.arange(x_left, x_right + step2, step2)
        
        ############### interpolation of P(omega) after integration over z ###############
        EELS_vs_energy = interp1d(list_energy0, list_P_integrated_over_z)
        integration_over_energy = np.sum(EELS_vs_energy(x_integrate))*(step2)
        
        np.savetxt('list_energy_for_P_integration' + label_Ee + all_info_label + '.txt', x_integrate, fmt='%.10f', delimiter='\t', header = header, encoding=None)
        list_Pvalue.append(integration_over_energy)
        
        print(j,V0,theta_mrad,integration_over_energy)
        
        ######################## OLD #########################################################################################################
        
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
    
    print('7-Integrate over energy from x_left_peak, x_right_peak as a function of (theta,V0) for a fixed b_min = %.2f' %(bmin_vals0))
    
    ind_max = len(list_Pvalue)
    table_P_integrated_over_energy = np.transpose([listx_sorted[0:ind_max],listy_sorted[0:ind_max], list_Pvalue])
    header1 = 'V0 (eV)    theta (mrad)    P_mode1, '
    np.savetxt('P_integrated_over_z_over_mode%i' %(mode) + label_Ee + '_dd%i_hh%.2f_bmin%.2f.txt' %(dd,bb,bmin_vals0), table_P_integrated_over_energy, fmt='%.10f', delimiter='\t', header = header1 + header, encoding=None)


#%%
    
# print('8-Plot the result of the integration')

# labelx=r'$\theta$ (mrad)'
# labely='EELS of n = %i' %(mode)

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


 