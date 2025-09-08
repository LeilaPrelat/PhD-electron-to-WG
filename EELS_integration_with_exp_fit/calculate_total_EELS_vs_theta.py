
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila
calculate the EELS of a rectangular waveguide
integrated over z/W as a function of energy
using the exponential fits of the EELS
"""
import numpy as np
import os
import matplotlib.pyplot as plt
# import concurrent.futures
import matplotlib as mpl
from scipy.interpolate import interp1d
from EELS_integrated_over_trajectory import P_integrated_over_z, find_width_of_peak
import pandas as pd
import ast
## fit with a lorentzian EELS vs energy to find the width of a mode and integrate over energy later 

path_basic = os.getcwd()
path_data_zmin = os.path.join(path_basic, 'potential_data')
path_data_EELS_integrated = os.path.join(path_basic, 'EELS_integrated')

#%%

tamfig = [4,3]
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
xe_norm_w = position_of_the_electron[ind_xe]
##### potential from c++ code  ##################################################
wp_norm_w = 0.5 ## w'/w
d_norm_w = 0.2  ## d/w distance between the side wires normalized to w
h_norm_w = 0.4   ## aspect ratio height/w
x0 = xe_norm_w    ## V(xe,z) with z from z0 to z1
 
# x=0 is in the middle of the center waveguide
# z=0 is at the interface
Nc = 400 ## points of the potential code
                        
bmin_norm_w = 50/w ## 80/w ## 100/w
z_min_val_norm_w = bmin_norm_w + h_norm_w
###################################################################################
h=h_norm_w*w
# mode = 1
# index_mode = -mode ## if mode = 1, is the highest one because is the last element of the list (the array is sorted from min to max)
plot_figure = 0    ## plot the lorentzian fitting
step2 = 0.001      ## thinner range to integrate over energy

mode = 1
index_mode = -1 ## maximum peak always. 

#%%

print('Extract info of the modes from data_modes.txt')

# Step 1: Read file, skipping comment lines
try: 
    with open("data_modes.txt", "r") as f:
        lines = [line.strip() for line in f if not line.startswith("##") and line.strip()]
except FileNotFoundError:
    raise('data_modes.txt not found. Copy paste it from /bem2D' )
# Step 2: Split the remaining data line into parts
data_parts = lines[0].split('\t')  # assuming tab-separated

# Step 3: Convert each string like '[500, 200]' into actual Python lists
parsed_parts = [ast.literal_eval(part) for part in data_parts]

# Step 4: Assign names for clarity (customize as needed)
width_height = parsed_parts[0]
w1 = parsed_parts[1]
w2 = parsed_parts[2]
w3 = parsed_parts[3]

# Step 5: Put into DataFrame (optional but handy)
df = pd.DataFrame({
    "width": [width_height[0]],
    "height": [width_height[1]],
    "w1_left": [w1[0]],
    "w1_right": [w1[1]],
    "w1": [w1[2]],
    "w2_left": [w2[0]],
    "w2_right": [w2[1]],
    "w2": [w2[2]],
    "w3_left": [w3[0]],
    "w3_right": [w3[1]],
    "w3": [w3[2]],
})

if width_height == [w,h]:
    x_left_limits = np.array([w1[0], w2[0],w3[0]])
    x_right_limits = np.array([w1[1], w2[1],w3[1]])

list_energy = np.arange(x_left_limits[mode-1], x_right_limits[mode-1], step2)

#%%
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
try: 
    table_solutions = np.loadtxt(name_sol, delimiter='\t', skiprows=1, encoding=None)
except FileNotFoundError:
    raise(name_sol, ' not found. Run run_potential_norm_V0.sh and save the files in /potential_data' )
table_solutions2 = np.transpose(table_solutions)
listx_sorted, listy_sorted = table_solutions2

    
name_potential = 'potential_wp%.1f_d%.1f_h%.1f_N%i_xe%.1f.txt' %(wp_norm_w,d_norm_w,h_norm_w,Nc,xe_norm_w)
try: 
    tabla_sh = np.loadtxt(name_potential,skiprows=1,delimiter=' ') 
except FileNotFoundError:
    raise(name_potential,'not found. Run run_potential_norm_V0.sh and save the files in /potential_data' )

tabla_sh2 = np.transpose(tabla_sh)
list_z_norm_w_sh, listV_normV0_sh = tabla_sh2
    
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
plt.savefig('zmin_aux' + label_Ee + label_txt + '.pdf', format='pdf',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
plt.show()

header = title
Ntot = len(listx_sorted)
# Ntot = 1

 

 
#%%

print('4-Calculate EELS integrated over z vs energy for all (V0,theta)')
#print('IMPORTANT: we are assuming the position of the peak DOES vary with V0,theta')



header0 = 'energy (eV)    P(omega) (1/eV), '
list_Pvalue_tot = []
for j in range(Ntot): ## Ntot
    V0 = listx_sorted[j]
    theta_mrad = listy_sorted[j]
    theta = theta_mrad*1e-3
    
    
    info_bmin = '_V0_%.5feV_theta%.5fmrad_bmin%.2f' %(V0,theta_mrad,bmin_norm_w) ## info from contour plot from plot_bmin_vs_V0_theta.py    
    all_info_label = '_mode%i_w%inm' %(mode,w) + label_Ee + label_txt + info_bmin
    name_P_integrated_over_z1 = 'P_integrated_exp_fit_over_z_w%inm_mode%i_Ee%ikeV_wp%.2f_h%.2f_d%.2f_xe%.2f' %(w,mode,Ee_electron_keV,wp_norm_w,h_norm_w,d_norm_w,xe_norm_w) + info_bmin + '.dat'
 
    list_Pvalue = []
    os.chdir(path_data_zmin)
    for energy in list_energy:
        energy = np.round(energy,5)
        Pvalue = P_integrated_over_z(z_min_val_norm_w, energy, w, h, Ee_electron_keV, theta, V0, list_z_norm_w_sh, listV_normV0_sh, mode)
        list_Pvalue.append(Pvalue)

    os.chdir(path_data_EELS_integrated)
    
    np.savetxt(name_P_integrated_over_z1, list_Pvalue, fmt='%.10f', delimiter='\t', header = header, encoding=None)
 

 
    results = find_width_of_peak(list_energy,list_Pvalue,index_mode,title,plot_figure)
    

    x_left, x_right, x0_fit, function_lorentzian, fit_success = results 
    x_integrate = np.arange(x_left, x_right + step2, step2) ## width of the lorenztian
    ## function_lorentzian is the integration of the fit (the lorentzian)
    bandwidth  = x_right - x_left
    
    if fit_success:
        ############### interpolation of P(omega) after integration over z --> integration over mode ###############
        EELS_vs_energy = interp1d(list_energy, list_Pvalue)
        try: 
            integration_over_energy = np.sum(EELS_vs_energy(x_integrate))*(step2) # Integration by summation (small steps)
        except ValueError:
            integration_over_energy = function_lorentzian
        
 
        list_Pvalue_tot.append(function_lorentzian) ## integration using the lorentzian fit 
        
        print(j,V0,theta_mrad,integration_over_energy,function_lorentzian)

#%%
print('6-EELS integrated over energy from x_left_peak, x_right_peak as a function of (theta,V0) for a fixed b_min = %.2f' %(bmin_norm_w))


os.chdir(path_data_EELS_integrated)

info_label = '_mode%i_w%inm' %(mode,w) + label_Ee + label_txt + '_bmin%.2f' %(bmin_norm_w)

ind_max = len(list_Pvalue_tot)
table_P_integrated_over_energy = np.transpose([listx_sorted[0:ind_max],listy_sorted[0:ind_max], list_Pvalue_tot])
header1 = 'V0 (eV)    theta (mrad)    P_mode1 (photons), '
final_name = 'P_integrated_exp_fit_over_z_over' + info_label + '.txt'

np.savetxt(final_name, table_P_integrated_over_energy, fmt='%.10f', delimiter='\t', header = header1 + header, encoding=None)


#%%