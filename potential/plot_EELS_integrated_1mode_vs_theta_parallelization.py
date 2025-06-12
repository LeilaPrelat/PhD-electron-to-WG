
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
import concurrent.futures

create_data = 1
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

plot_figure = 0
list_energy = np.arange(0.2,1.505,0.005)
list_z_norm_a, listV_normV0 = dy_cached(bb,ss,dd,N) ## z/W and V(z)/V0 from c++ code

#%%

print('2-Load data of z_min from the code plot_potential_V_V0_theta_only_above_wg.py')

bmin_vals = np.loadtxt('bmin' + label_Ee + '_dd%i_hh%.2f.txt' %(dd,bb),  delimiter='\t', skiprows = 1)
V0_vals = np.loadtxt('V0' + label_Ee + '_dd%i_hh%.2f.txt' %(dd,bb), delimiter='\t', skiprows = 1)
theta_mrad_vals = np.loadtxt('theta_mrad' + label_Ee + '_dd%i_hh%.2f.txt' %(dd,bb), delimiter='\t', skiprows = 1)

# find the index closest to the values we want for bmin 
bmin_vals0 = 0.2
bmin_vals1 = 0.25
bmin_vals3 = 0.5
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
cbar.ax.set_title(r'$b_{min}/W$',fontsize=tamletra)
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
np.savetxt('V0_vals_for_bmin%.2f' %(bmin_vals0) + label_Ee + '_dd%i_hh%.2f.txt' %(dd,bb), listx_sorted, fmt='%.10f', delimiter='\t', header = header, encoding=None)
np.savetxt('theta_mrad_vals_for_bmin%.2f' %(bmin_vals0) + label_Ee + '_dd%i_hh%.2f.txt' %(dd,bb), listy_sorted, fmt='%.10f', delimiter='\t', header = header, encoding=None)

#%%

np.savetxt('list_energy_of_P' + label_Ee + '_dd%i_hh%.2f.txt' %(dd,bb), list_energy, fmt='%.10f', delimiter='\t', header = header, encoding=None)

def Pintegrated_over_z_vs_theta(V0,theta_mrad,bmin_vals,Ee_electron, bb, ss, dd, a, N, list_z_norm_a, listV_normV0):
    theta = theta_mrad*1e-3
    z_min_val = bmin_vals0 + bb/2  
    list_P_integrated_over_z = []
    for energy in list_energy:
        P_int_value = P_integrated_over_z(z_min_val, theta, V0, Ee_electron, bb, ss, dd, energy, a, N, list_z_norm_a, listV_normV0)
        # print(energy,P_int_value,V0,theta)
        list_P_integrated_over_z.append(P_int_value)
    
    # Build unique filename
    base_path = "E:\Desktop\Leila\EELs_omega_vs_theta\PhD-electron-to-WG-main\potential\EELS_integrated_over_z"
    filename0 = r"width%inm_dd%i_hh%.2f_V0%.2feV_theta%.2fmrad_bmin%.2f.txt" %(a,dd,bb,V0,theta_mrad,bmin_vals0)
    filename = "P_integrated_over_z" + label_Ee + filename0
    full_path = os.path.join(base_path, filename)
    
    # Save result to file during parallelization
    print('saving file for V0=%.2f, theta=%.2f to %s' %(V0, theta_mrad,base_path))
    with open(full_path, "w") as f:
        f.write(str(list_P_integrated_over_z))
        
    # title2 = r'$d/W = %i$, $V_0$ = %.2f eV, $\theta$ = %.2f mrad' %( dd,V0,theta_mrad)
    # header = title1 + ', ' + title2  
    # np.savetxt('P_integrated_over_z' + label_Ee + '_dd%i_hh%.2f_V0%.2feV_theta%.2fmrad_bmin%.2f.txt' %(dd,bb,V0,theta_mrad,bmin_vals0), list_P_integrated_over_z, fmt='%.10f', delimiter='\t', header = header, encoding=None)
    
    # return list_P_integrated_over_z


    return list_P_integrated_over_z

print('5-Integration of P(omega) over z/a as a function of theta,V0 for a fixed b_min = %.2f' %(bmin_vals0))
max_workers = 25
import time

if __name__ == "__main__":
    start = time.time()
 
    ## process dont see global variables, I need to write them again as before ## 
    bmin_vals0 = 0.2 
    Ee_electron_keV = 200
    Ee_electron = Ee_electron_keV*1e3
    
    bb = h/a ## h/a ## zmin is bb/2 (upper surface of the wg)
    ss = 0.1
    dd = 5   ## distance to the plane
    a = 300 ## nm
    N = 400
    list_z_norm_a, listV_normV0 = dy_cached(bb,ss,dd,N) ## z/W and V(z)/V0 from c++ code
    listx_sorted = np.loadtxt('V0_vals_for_bmin%.2f' %(bmin_vals0) + label_Ee + '_dd%i_hh%.2f.txt' %(dd,bb), delimiter='\t', skiprows = 1)
    listy_sorted = np.loadtxt('theta_mrad_vals_for_bmin%.2f' %(bmin_vals0) + label_Ee + '_dd%i_hh%.2f.txt' %(dd,bb), delimiter='\t', skiprows = 1)
    
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(Pintegrated_over_z_vs_theta, V0, theta_mrad,
                                   bmin_vals,Ee_electron, bb, ss, dd, a, N, list_z_norm_a, listV_normV0)
                  for V0, theta_mrad in zip(listx_sorted, listy_sorted)]
        #print(V0,theta_mrad)
        results = [future.result() for future in concurrent.futures.as_completed(futures)]
        print(futures)
 
 
    end = time.time()
    print(f"Finished in {end - start:.2f} seconds")



# if __name__== "__main__": # necessary in windows to parallelize 

#     # list_P_integrated_over_z_tot = []
#     with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
#         results =  list(executor.map(Pintegrated_over_z_vs_theta, listx_sorted, listy_sorted))
 
#     # Save all results after parallelization
#     for V0, theta_mrad, result in results:
#         filename = f"result_V0_{V0}_theta_{theta_mrad}.txt"
#         with open(filename, "w") as f:
#             f.write(str(result))
              
 
            # list_P_integrated_over_z_tot.append(results)
            # you can also process result immediately here
        
    # for j in range(len(listx_sorted)): 
    #     V0 = listx_sorted[j]
    #     theta_mrad = listy_sorted[j]
    #     theta = theta_mrad*1e-3
    #     list_P_integrated_over_z = Pintegrated_over_z_vs_theta(V0,theta,bmin_vals0,bb,dd)
    #     list_P_integrated_over_z_tot.append(list_P_integrated_over_z)

 #%%
 




#%%



