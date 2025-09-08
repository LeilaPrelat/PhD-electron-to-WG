
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

#%%
# //     N   = the total number of parametrization points is N+50
# //     a   = total width of wg     || x
# //     h   = total thickness of wg  || y
# //     s   = rounding radius of wg (s<t/2)

print('1-Define the parameters')

############ parameters for bem2d ############
a = 500 ## nm ## 300 o 400
h = 300 ## nm
s = a*0.1  ## nm
############ parameters for c++ normalized by a ############
bb = h/a ## h/a ## zmin is bb/2 (upper surface of the wg)
ss = s/a
dd = 10/(a*1e-3) ## = d/a  ### a = 400 nm : 12.5 for d = 5 microns / 25 for d = 10 microns / 50 for d = 20 microns 
           ### a = 300 nm : 3.33 for d = 1 microns
           ### a = 400 nm : 2.5 for d = 1 microns
           

dd = 20

N = 400
Ee_electron_keV = 200
label_Ee = '_Ee%ikeV' %(Ee_electron_keV)

# step = 0.005
# step = 0.000713 ## thinner steps
# step = 0.000712899999999999
# list_energy = np.arange(0.2,1.5 + step,step) ## we will reduce this energy list to the energies around the first peak to integrate over \omega
# list_energy = np.arange(0.01,5 + step,step)
 

list_z_norm_a, listV_normV0 = dy_cached(bb,ss,dd,N) ## import z/W and V(z)/V0 from c++ code

# mode = 1
# index_mode = -mode ## if mode = 1, is the highest one because is the last element of the list (the array is sorted from min to max)
plot_figure = 0    ## plot the lorentzian fitting
step2 = 0.001      ## thinner range to integrate over energy

list_mode = [1,2]
# list_mode = [3]

# find the index closest to the values we want for bmin 
if a == 300: 
    bmin_vals0 = 0.15
    # bmin_vals0 = 0.1
    
elif a == 400: ## LETS USES ONLY THIS VALUE
    bmin_vals0 = 0.1    ## same bmin = 60 nm
    # bmin_vals0 = 0.2
    bmin_vals0 = 0.2

bmin_vals0 = 80/a
z_min_val = bmin_vals0 + bb/2 
label_txt = '_dd%.2f_hh%.2f.txt' %(dd,bb)


position_of_the_electron = [0,40,80]
ind_xe=1
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
    

# list_energy0 = tabla_energy_2[0:-1]


#%%

print('2-Load data from the code plot_bmin_vs_V0_theta.py of theta,V0 for a fix bmin')

os.chdir(path_data_zmin)
bmin_vals = np.loadtxt('bmin' + label_Ee + label_txt, delimiter='\t', skiprows = 1)
V0_vals = np.loadtxt('V0_for_bmin%.2f' %(bmin_vals0) + label_Ee + label_txt, delimiter='\t', skiprows = 1)
theta_mrad_vals = np.loadtxt('theta_mrad_for_bmin%.2f' %(bmin_vals0) + label_Ee + label_txt, delimiter='\t', skiprows = 1)
## solutions (V0,theta) such that bmin = bmin_vals0
listx_sorted = np.loadtxt('V0_sol_for_bmin%.2f'%(bmin_vals0)  + label_Ee + label_txt, delimiter='\t', skiprows = 1)
listy_sorted = np.loadtxt('theta_sol_mrad_for_bmin%.2f'%(bmin_vals0) + label_Ee + label_txt, delimiter='\t', skiprows = 1)
 
tamfig2 = [5, 3]

print('3-Plot b_min as a function of (V0,theta) and see if they match with the contorn plot')

title1 = r'$W$ = %i nm, $h$ = %i nm, $s$ = %i nm, $d/W = %i$, $E_{\rm e}$ = %i keV' %(a,h,s,dd,Ee_electron_keV)
# another cbar with constant a*1e-3 (real units, in microns)
cte_cbar2 = a*1e-3
cmap = plt.cm.RdBu   # define the colormap
bounds1 =   np.logspace(np.log10(0.1), np.log10(100) , 10) 
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
# plt.plot(listx_sorted,listy_sorted,'--',color = 'green')
plt.savefig('zmin_aux' + label_Ee + 'bmin%.2f_dd%i_hh%.2f.png' %(bmin_vals0,dd,bb), format='png',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
plt.show()



ze=50
name_list_energy = 'energy_list_N%i_a%inm_h%inm_ze%inm_Ee%ikeV%s'%(N,a,h,ze,Ee_electron_keV,labelxe) + '.txt'
os.chdir(path_data)
tabla_energy = np.loadtxt(name_list_energy, delimiter=' ',dtype=None)
tabla_energy_2 = np.transpose(tabla_energy)
list_energy0 = tabla_energy_2


header = title1
Ntot = len(listx_sorted)
# Ntot = 3
#%%

print('4-Plot EELS vs energy to fit it and find the width for all (V0,theta), 6-later use the width to integrate over energy')
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
        
        k=0
        list_energy0_EELS_positive = []
        for energy in list_energy0:
            
            try:
                P_int_value,list_z_norm_a_EELS,EELS_array = P_integrated_over_z(z_min_val, theta, V0, xe, Ee_electron_keV, bb, dd, energy, a, N, list_z_norm_a, listV_normV0)
                # print(energy,P_int_value)
                list_P_integrated_over_z.append(P_int_value)
                list_energy0_EELS_positive.append(energy)
                # print(k,energy,P_int_value)
                k=k+1
                # plt.figure()
                # plt.plot(list_z_norm_a_EELS,EELS_array)
            except ValueError: ## some EELS are negative, we omit them 
                pass
            
        table_P_integrated_over_z = np.transpose([list_energy0_EELS_positive, list_P_integrated_over_z])
        os.chdir(path_data_created)
        np.savetxt('P_integrated_over_z' + label_Ee + all_info_label + '.txt', table_P_integrated_over_z, fmt='%.10f', delimiter='\t', header = header0 + header, encoding=None)
        if j == 0:
            print('6-Find the wmin and wmax of the mode = %i by fitting it with a Lorenztian' %(mode))
        

        ######### find peaks and sort them by the highest to minimum (then we identify each mode according to its amplitude) by fitting the curve with a lorentzian
        x_left, x_right, x0_fit, function_lorentzian, fit_success = find_width_of_peak(list_energy0_EELS_positive,list_P_integrated_over_z,index_mode,title1,plot_figure)
        # Integration by summation (small steps)
        x_integrate = np.arange(x_left, x_right + step2, step2)
        
        if fit_success:
            ############### interpolation of P(omega) after integration over z ###############
            EELS_vs_energy = interp1d(list_energy0_EELS_positive, list_P_integrated_over_z)
            integration_over_energy = np.sum(EELS_vs_energy(x_integrate))*(step2)
            
            np.savetxt('list_energy_for_P_integration' + label_Ee + all_info_label + '.txt', x_integrate, fmt='%.10f', delimiter='\t', header = header, encoding=None)
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
    
    print('5-Integrate over energy from x_left_peak, x_right_peak as a function of (theta,V0) for a fixed b_min = %.2f' %(bmin_vals0))
    
    ind_max = len(list_Pvalue)
    table_P_integrated_over_energy = np.transpose([listx_sorted[0:ind_max],listy_sorted[0:ind_max], list_Pvalue])
    header1 = 'V0 (eV)    theta (mrad)    P_mode1, '
    np.savetxt('P_integrated_over_z_over_mode%i' %(mode) + label_Ee + '_dd%i_hh%.2f_bmin%.2f.txt' %(dd,bb,bmin_vals0), table_P_integrated_over_energy, fmt='%.10f', delimiter='\t', header = header1 + header, encoding=None)


#%%
