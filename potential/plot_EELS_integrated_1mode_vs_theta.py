
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila
plot the EELS of a rectangular waveguide
integrated over a peak
"""
import numpy as np
import os
import matplotlib.pyplot as plt
from EELS_integrated import P_integrated_over_z, dy_cached
# import concurrent.futures
from scipy.optimize import curve_fit

create_data = 0
path_basic = os.getcwd()
path_data = os.path.join(path_basic, 'bem_files_EELS')

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
a = 300 ## nm
h = 300 ## nm
s = 20  ## nm
############ parameters for c++ ############
bb = h/a ## h/a ## zmin is bb/2 (upper surface of the wg)
ss = 0.1
dd = 5   ## distance to the plane
N = 400
Ee_electron_keV = 200
Ee_electron = Ee_electron_keV*1e3
label_Ee = '_Ee%ikeV' %(Ee_electron_keV)

step = 0.005
list_energy0 = np.arange(0.2,1.505,step) ## we will reduce this energy list to the energies around the first peak to integrate over \omega
list_z_norm_a, listV_normV0 = dy_cached(bb,ss,dd,N) ## import z/W and V(z)/V0 from c++ code

#%%

print('2-Load data of z_min from the code plot_potential_V_V0_theta_only_above_wg.py')

bmin_vals = np.loadtxt('bmin' + label_Ee + '_dd%i_hh%.2f.txt' %(dd,bb),  delimiter='\t', skiprows = 1)
V0_vals = np.loadtxt('V0' + label_Ee + '_dd%i_hh%.2f.txt' %(dd,bb), delimiter='\t', skiprows = 1)
theta_mrad_vals = np.loadtxt('theta_mrad' + label_Ee + '_dd%i_hh%.2f.txt' %(dd,bb), delimiter='\t', skiprows = 1)

# find the index closest to the values we want for bmin 
bmin_vals0 = 0.15
bmin_vals1 = 0.25
bmin_vals3 = 0.5
z_min_val = bmin_vals0 + bb/2 
print('3-Find the V0,theta values for a fixed b_min = %.2f' %(bmin_vals0))
print('4-Plot b_min as a function of (V0,theta) and see if they match with the contorn plot')

title1 = r'$W$ = %i nm, $h$ = %i nm, $s$ = %i nm, $d/W = %i$, $E_{\rm e}$ = %i keV' %(a,h,s,dd,Ee_electron_keV)

limits1 = [np.min(V0_vals) , np.max(V0_vals),np.min(theta_mrad_vals) , np.max(theta_mrad_vals)]
cmap = plt.cm.hot  # define the colormap

# Create contour plot
contour_levels = [bmin_vals0]


plt.figure(figsize=tamfig)
#plt.title(title1 + r', $E_{\text{e}}$ = %i keV' %(Ee_electron_keV),fontsize=tamtitle)
im_show = plt.imshow(bmin_vals, extent = limits1, cmap=cmap, aspect='auto', interpolation = 'bicubic',origin = 'lower'   ) 
contours = plt.contour(V0_vals, theta_mrad_vals,bmin_vals, levels=contour_levels, colors='green', linestyles='dashed')

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

plt.clabel(contours, fmt='%.2f', colors='green',fontsize=tamletra)  # Label contours
cbar = plt.colorbar(im_show, fraction=0.046, pad=0.04   ,format = '%i') 
cbar.ax.set_title(r'$b_{\text{min}}/W$',fontsize=tamletra)
cbar.ax.tick_params(labelsize = tamnum, width=0.1, direction="in",which = 'both', length = 2,pad = pad)
# plt.xticks(np.arange(0,np.max(V0_vals)+0.5,0.5))
plt.plot(listx_sorted,listy_sorted,'-',color = 'blue')
plt.xlabel(r'$V_0$ (eV)',fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(r'$\theta$ (mrad)',fontsize=tamletra,labelpad =labelpady)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.savefig('zmin_aux' + label_Ee + '_dd%i_hh%.2f.png' %(dd,bb), format='png',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
#plt.title('Root x as function of V0 and theta')
plt.show()

header = title1 + r', $E_{\text{e}}$ = %i keV' %(Ee_electron_keV)
os.chdir(path_data)
np.savetxt('V0_vals_for_bmin%.2f' %(bmin_vals0) + label_Ee + '_dd%i_hh%.2f.txt' %(dd,bb), listx_sorted, fmt='%.10f', delimiter='\t', header = header, encoding=None)
np.savetxt('theta_mrad_vals_for_bmin%.2f' %(bmin_vals0) + label_Ee + '_dd%i_hh%.2f.txt' %(dd,bb), listy_sorted, fmt='%.10f', delimiter='\t', header = header, encoding=None)

#%%

print('5-Plot EELS vs energy for a fixed (V0,theta) to fit it and find the width')
print('IMPORTANT: we are assuming the position of the peak does not vary with V0,theta')

V0 = listx_sorted[0]
theta_mrad = listy_sorted[0]
theta = theta_mrad*1e-3
all_info_label = '_dd%i_hh%.2f_V0%.2feV_theta%.2fmrad_bmin%.2f.txt' %(dd,bb,V0,theta_mrad,bmin_vals0)

if create_data == 1:
    list_P_integrated_over_z = []
    for energy in list_energy0:
        print(energy)
        P_int_value = P_integrated_over_z(z_min_val, theta, V0, Ee_electron, bb, ss, dd, energy, a, N, list_z_norm_a, listV_normV0)
        list_P_integrated_over_z.append(P_int_value)
    
    os.chdir(path_data)
    np.savetxt('P_integrated_over_z' + label_Ee + all_info_label, list_P_integrated_over_z, fmt='%.10f', delimiter='\t', header = header, encoding=None)
    np.savetxt('list_energy_of_P' + label_Ee + all_info_label, list_energy0, fmt='%.10f', delimiter='\t', header = header, encoding=None)

else:
    os.chdir(path_data)
    list_P_integrated_over_z = np.loadtxt('P_integrated_over_z' + label_Ee + all_info_label, delimiter='\t', skiprows=1)
    list_energy0 = np.loadtxt('list_energy_of_P' + label_Ee + all_info_label,   delimiter='\t', skiprows=1)
    
#%%

print('6-Find the wmin and wmax of the highest mode by fitting it with a Lorenztian')

# Lorentzian function
def lorentzian(x, A, x0, gamma, offset):
    return A * gamma**2 / ((x - x0)**2 + gamma**2) + offset

ind_max = int(np.argmax(list_P_integrated_over_z))
ind0 = int(ind_max*0.8)
ind1 = int(ind_max*1.1)
x_data = list_energy0[ind0:ind1] ## we fit around the maximum
y_data = list_P_integrated_over_z[ind0:ind1]

# Initial parameter guesses: [A, x0, gamma, offset]
initial_guess = [max(y_data), list_energy0[ind_max], 0.1, min(y_data)]

# Perform the fit
popt, pcov = curve_fit(lorentzian, x_data, y_data, p0=initial_guess)
A_fit, x0_fit, gamma_fit, offset_fit = popt

# Threshold at 1% of amplitude
threshold = offset_fit + 0.01 * A_fit

# Solve for x where Lorentzian equals threshold
delta = np.sqrt((A_fit * gamma_fit**2) / (threshold - offset_fit) - gamma_fit**2)
x_left = x0_fit - delta
x_right = x0_fit + delta
bottom_width = x_right - x_left

title2 = r'$V_0$ = %.2f eV, $\theta$ = %.2f mrad' %(V0,theta_mrad)
 
plot_figure = 1
if plot_figure == 1:
    
    labelx='Photon energy $\hbar\omega$ (eV)'
    labely='EELS per electron (1/eV)'
     
    tamfig = [4.5,3.5]
    tamletra = 13
    tamtitle  = tamletra - 3
    tamnum = tamletra
    labelpady = 3
    labelpadx = 2
    pad = 2.5
    dpi = 500
    
    plt.title(title1 + ', ' + title2,fontsize=tamtitle)
    plt.figure(figsize=tamfig)
    plt.title('',fontsize=tamtitle)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.plot(list_energy0, list_P_integrated_over_z ,'.-' )
    plt.plot(x_data, lorentzian(x_data, *popt), 'r-', label='Lorentzian fit')
    # plt.plot(peak_x, peak_y, "x")
    # plt.plot(x_left_peak_value, np.ones(len(x_left_peak_value))*0, "x",color = 'blue')
    # plt.plot(x_right_peak_value, np.ones(len(x_right_peak_value))*0, "x",color = 'red')
    plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)

#%%

print('7-Integrate over omega from x_left_peak, x_right_peak as a function of (theta,V0) for a fixed b_min = %.2f' %(bmin_vals0))

# Find the index of the closest value
idx_x_left = (np.abs(list_energy0 - x_left)).argmin()
closest_x_left = list_energy0[idx_x_left]
idx_x_right = (np.abs(list_energy0 - x_right)).argmin()
closest_x_right = list_energy0[idx_x_right]


# we are assuming the position of the peak does not vary with V0,theta
# list_energy_over_mode = np.arange(x_left_peak_value[0], x_right_peak_value[0] + 0.005,0.005)
list_energy_over_mode = np.arange(closest_x_left, closest_x_right + step, step)
np.savetxt('list_energy_of_P' + label_Ee + '_dd%i_hh%.2f.txt' %(dd,bb), list_energy_over_mode, fmt='%.10f', delimiter='\t', header = header, encoding=None)

def Pintegrated_over_energy(V0,theta_mrad,list_energy):
    theta = theta_mrad*1e-3
     
    P_integrated_over_energy = 0
    for energy in list_energy:
        P_int_value = P_integrated_over_z(z_min_val, theta, V0, Ee_electron, bb, ss, dd, energy, a, N, list_z_norm_a, listV_normV0)
 
        P_integrated_over_energy = P_int_value + P_integrated_over_energy
    
    delta_energy = list_energy[1] - list_energy[0]
 
    return P_integrated_over_energy*delta_energy


# Integration by summation (small steps)
x_integrate = np.arange(closest_x_left, closest_x_right + step, step)
y_integrate = lorentzian(x_integrate, *popt)
integral_sum = np.sum(y_integrate * step)
P_example = Pintegrated_over_energy(V0,theta_mrad,list_energy_over_mode)
print("Difference between using the Lorenztian to integrate over energy and using the data:",np.abs(integral_sum,P_example))

#%%


list_Pvalue = []
for j in range(len(listx_sorted)): 
    V0 = listx_sorted[j]
    theta_mrad = listy_sorted[j]
    Pvalue = Pintegrated_over_energy(V0,theta_mrad,list_energy_over_mode)
    list_Pvalue.append(Pvalue)
    print(j,V0,theta_mrad,Pvalue)

header = title1 + ', ' + title2  
np.savetxt('P_integrated_over_z_over_1mode' + label_Ee + '_dd%i_hh%.2f_bmin%.2f.txt' %(dd,bb,bmin_vals0), list_Pvalue, fmt='%.10f', delimiter='\t', header = header, encoding=None)

    
    # list_P_integrated_over_z_tot = []
    # with concurrent.futures.ProcessPoolExecutor(max_workers=6) as executor:
    #     futures = [executor.submit(Pintegrated_over_z_vs_theta, V0,theta_mrad) for V0,theta_mrad in zip(listx_sorted,listy_sorted)]
        
    #     for future in concurrent.futures.as_completed(futures):
    #         result = future.result()
    #         list_P_integrated_over_z_tot.append(result)
            # you can also process result immediately here
        

#%%



