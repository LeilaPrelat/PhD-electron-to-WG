#code d.h : 
#- Models a grounded plane + biased rectangular wire with rounded corners.
#- Uses the boundary element method to compute induced surface charge.
#- Computes electrostatic potential/V0 at any point in space (x,y)
# where: 
#bb = 2, aspect ratio  (side length along y (labelled as b) divided by side length along x).
#ss = 0.1, rounding radius.
#dd = 1, distance to the plane.
#N = 200, discretization points.
### All lengths are given in units of the x side-length, label as "a"

## plot the EELS integrated over z/W
import numpy as np
import matplotlib.pyplot as plt
import os
from EELS_integrated import P_integrated_over_z, P_integrand_over_z, EELS_from_BEM_interpolated, dy_cached

create_data = 1

tamfig = [4.5,3.5]
tamletra = 13
tamtitle  = tamletra - 3
tamnum = tamletra
tamlegend = tamletra
labelpady = 3
labelpadx = 2
pad = 2.5
mk = 1
ms = 1
lw = 1.5
hp = 0.5
length_marker = 1.5
dpi = 500

path_basic = os.getcwd()
path_data = os.path.join(path_basic, 'bem_files_EELS')

#%%
# //     N   = the total number of parametrization points is N+50
# //     a   = total width of wg     || x
# //     h   = total thickness of wg  || y
# //     s   = rounding radius of wg (s<t/2)
print('1-Define the parameters')
############ parameters for bem2d ############
a=300 ## nm
h=300 ## nm
s=20  ## nm
############ parameters for c++ ############
bb = h/a ## h/a ## zmin is bb/2 (upper surface of the wg)
ss = 0.1
dd = 5   ## distance to the plane
N = 400

Nint = 10
Ee_electron_keV = 200 ## keV
Ee_electron = Ee_electron_keV*1e3
me_c2_eV = 510998.95069  ## me*c**2 in eV
beta = np.sqrt( 1- (1 + Ee_electron/me_c2_eV)**(-2) )  ## beta = v/c
energy_eV0 = 0.75
list_energy = np.arange(0.2,1.505,0.005)
list_energy = np.arange(0.2,3.005,0.005)
list_z_norm_a, listV_normV0 = dy_cached(bb,ss,dd,N) ## z/W and V(z)/V0 from c++ code

label_Ee = '_Ee%ikeV' %(Ee_electron_keV)

#%%

print('2-Load data of z_min from the code plot_potential_V_V0_theta_only_above_wg.py')

bmin_vals = np.loadtxt('bmin' + label_Ee + '_dd%i_hh%.2f.txt' %(dd,bb),  delimiter='\t', skiprows = 1)
V0_vals = np.loadtxt('V0' + label_Ee + '_dd%i_hh%.2f.txt' %(dd,bb), delimiter='\t', skiprows = 1)
theta_mrad_vals = np.loadtxt('theta_mrad' + label_Ee + '_dd%i_hh%.2f.txt' %(dd,bb), delimiter='\t', skiprows = 1)

# find the index closest to the values we want for bmin 
bmin_vals0 = 0.15
bmin_vals1 = 0.25
bmin_vals3 = 0.5
print('3-Define the value of b_min = %.2f with the respective (V0,theta)' %(bmin_vals0))
distance_bmin = np.abs(bmin_vals - bmin_vals0)
index = np.unravel_index(np.nanargmin(distance_bmin, axis=None), distance_bmin.shape)
z_min_val = bmin_vals0 + bb/2   
V0 = V0_vals[index[0]] ## eV
theta_mrad = theta_mrad_vals[index[1]]
theta = theta_mrad*1e-3

print('4-Plot b_min as a function of V0,theta and mark the selected point in blue')

title1 = r'$W$ = %i nm, $h$ = %i nm, $s$ = %i nm, $E_{\rm e}$ = %i keV' %(a,h,s,Ee_electron_keV)
title2 = r'$d/W = %i$, $V_0$ = %.2f eV, $\theta$ = %.2f mrad' %(dd,V0,theta_mrad)

limits1 = [np.min(V0_vals) , np.max(V0_vals),np.min(theta_mrad_vals) , np.max(theta_mrad_vals)]
cmap = plt.cm.hot  # define the colormap

# Create contour plot
contour_levels = [0.15,   0.5]

plt.figure(figsize=tamfig)
#plt.title(title1 + r', $E_{\text{e}}$ = %i keV' %(Ee_electron_keV),fontsize=tamtitle)
im_show = plt.imshow(bmin_vals, extent = limits1, cmap=cmap, aspect='auto', interpolation = 'bicubic',origin = 'lower'   ) 
contours = plt.contour(V0_vals, theta_mrad_vals,bmin_vals, levels=contour_levels, colors='green', linestyles='dashed')
plt.clabel(contours, fmt='%.2f', colors='green',fontsize=tamletra)  # Label contours
plt.plot(V0,theta_mrad,'o',color = 'blue')
cbar = plt.colorbar(im_show, fraction=0.046, pad=0.04   ,format = '%i') 
cbar.ax.set_title(r'$b_{\text{min}}/W$',fontsize=tamletra)
cbar.ax.tick_params(labelsize = tamnum, width=0.1, direction="in",which = 'both', length = 2,pad = pad)
# plt.xticks(np.arange(0,np.max(V0_vals)+0.5,0.5))
plt.xlabel(r'$V_0$ (eV)',fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(r'$\theta$ (mrad)',fontsize=tamletra,labelpad =labelpady)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.savefig('zmin_aux' + label_Ee + '_dd%i_hh%.2f.png' %(dd,bb), format='png',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
#plt.title('Root x as function of V0 and theta')
plt.show()

#%%

print('6-Plot P(omega,z) as a function of z/W for a fixed energy')

listz_norm_a_BEM, EELS_interp = EELS_from_BEM_interpolated(energy_eV0,a,h,N)
listz_norm_a_BEM_interp = np.linspace(z_min_val, np.max(listz_norm_a_BEM),int(N*Nint)) 

labelx='z/W'
labely=r'$P(\hbar\omega,z)$'
plt.figure(figsize=tamfig)    
listP = P_integrand_over_z(listz_norm_a_BEM, theta, V0, Ee_electron, bb, ss, dd, energy_eV0, a, N , list_z_norm_a, listV_normV0)
aux_eje_y = np.linspace(np.min(listP), np.max(listP),10)
plt.plot(listz_norm_a_BEM,listP,'.-',label = r'$\hbar\omega = %.2f$ eV' %(energy_eV0))
plt.plot(np.ones(10)*z_min_val,aux_eje_y,'-',color = 'black')
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
os.chdir(path_basic)
plt.savefig('integrandP' + label_Ee + '_dd%i_hh%.2f.png' %(dd,bb), format='png',bbox_inches='tight',pad_inches = 0.09, dpi=dpi)  
plt.show()

#%%

print('7-Integration of P(omega) over z/a')

if create_data == 1:
    
    list_P_integrated_over_z = []
    for energy in list_energy:
        P_int_value = P_integrated_over_z(z_min_val, theta, V0, Ee_electron, bb, ss, dd, energy, a, N, list_z_norm_a, listV_normV0)
        print(energy,P_int_value)
        list_P_integrated_over_z.append(P_int_value)
        
    os.chdir(path_data)
    header = title1 + ', ' + title2 
    np.savetxt('P_integrated_over_z' + label_Ee + '_dd%i_hh%.2f.txt' %(dd,bb), list_P_integrated_over_z, fmt='%.10f', delimiter='\t', header = header, encoding=None)
    np.savetxt('list_energy_of_P' + label_Ee + '_dd%i_hh%.2f.txt' %(dd,bb), list_energy, fmt='%.10f', delimiter='\t', header = header, encoding=None)


else:
    os.chdir(path_data)
    list_P_integrated_over_z = np.loadtxt('P_integrated_over_z' + label_Ee + '_dd%i_hh%.2f.txt' %(dd,bb), delimiter='\t', skiprows = 1)
#    list_energy = np.loadtxt('list_energy_of_P' + label_Ee + '_dd%i_hh%.2f.txt' %(dd,bb), delimiter='\t', skiprows = 1) 

#%%

labelx='Electron energy loss $\hbar\omega$ (eV)'
labely='EELS per electron (1/eV)'

plt.figure(figsize=tamfig)
plt.title(title1 + '\n' + title2,fontsize=tamtitle)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
plt.plot(list_energy,list_P_integrated_over_z,'.-')
plt.xticks(np.arange(0.5,3.5,0.5))
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.tight_layout()
os.chdir(path_basic)
plt.savefig( 'EELS_integrated_over_z' + '_dd%i_hh%.2f.png' %(dd,bb),bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)
plt.show()

