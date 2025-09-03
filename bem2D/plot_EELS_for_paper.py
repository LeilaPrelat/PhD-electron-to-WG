
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila
plot the EELS of a rectangular waveguide
on top of substrate
with epsilon= 12.25 + i*0.2
we increase the losses to decrease the L
of the virtual geometry needed to converge 
    (read notes)

"""
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.signal import find_peaks
from fit_peaks import find_width_of_peak
import pandas as pd
import ast

path_basic = os.getcwd()
os.chdir(path_basic)

from global_constants import constants 
hb,c,alpha,me_c2_eV = constants()
aux = hb*c
epsi1, epsi3 = 1, 1

path_data = os.path.join(path_basic, 'datfiles_EELS')

list_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
               '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

#%%

tamfig = [4.3,3.5]
tamletra = 14
tamtitle  = 6
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
lw = 0.8

deltax,deltay = 3,8

values = tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad,deltax,deltay

# //     N   = the total number of parametrization points is N+50
# //     w   = total width of wg     || x
# //     h   = total thickness of wg  || y
# //     ze  = impact parameter is h/2+ze. ze distance between electron and wg surface  

W=400 # nm
h=400 ## nm

listL = [1000]
N0=500 ### plot vs L for N = N0
L0=700 ## plot vs N for L = L0 for mode = mode0 because mode0 has negative EELS

Ee_electron_keV = 200
Ee_electron = Ee_electron_keV*1e3
me_c2_eV = 510998.95069  ## me*c**2 in eV
beta = np.sqrt( 1- (1 + Ee_electron/me_c2_eV)**(-2) )  ## beta = v/c
gamma_e = 1/np.sqrt(1-epsi1*beta**2) ## gamma lorentz
ze = 50

xticks1=[0.70,0.75,0.80,0.85]
xticks2=[0.30,0.35,0.40]
xticks3=[1.45,1.50,1.55,1.60]

labelx=r'Photon energy $\hbar\omega$ (eV)'
labely=r'd$\Gamma$/d$y$ (1/eV$\cdot$nm)'
title = r'$W$ = %i nm, $h$ = %i nm, $z_{\rm e}$ = %i nm, $E_{\rm e}$ = %i keV, Im($\epsilon$) = 0.2' %(W,h,ze,Ee_electron_keV)
os.chdir(path_data)

#%%

print('Extract info of the modes from data_modes.txt')

# Step 1: Read file, skipping comment lines
with open("data_modes.txt", "r") as f:
    lines = [line.strip() for line in f if not line.startswith("##") and line.strip()]

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

if width_height == [W,h]:
    x_left_limits = np.array([w1[0], w2[0],w3[0]])
    x_right_limits = np.array([w1[1], w2[1],w3[1]])

#%%

print('Plot EELS with N = %i and L = %i for all modes' %(N0,L0))
  
plt.figure(figsize=tamfig)
plt.title(title,fontsize=tamtitle)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)

name3 = 'EELS_comsol2_spectrum_N%i_W%inm_h%inm_L%inm_Ee%ikeV.dat' %(N0,W,h,L0,Ee_electron_keV) ## virtual geometry
#name3 = 'EELS_comsol2_spectrum_N%i_W%inm_h%inm_L%inm_Ee%ikeV_mode%i.dat' %(N0,W,h,L,Ee_electron_keV,mode) ## virtual geometry
## the label "EELS_comsol2.." means the permittivity of Si is 12.25 + i*0.2

#name3 = 'EELS_N%i_W%inm_h%inm_L%inm_R%inm_Ee%ikeV.dat' %(N,W,h,L,R,Ee_electron_keV) ## virtual geometry

try:
    tabla1 = np.loadtxt(name3, delimiter=' ',dtype=None)
    tabla1_2 = np.transpose(tabla1)
except (ValueError, FileNotFoundError):
    tabla1 = np.loadtxt(name3, delimiter='\t',dtype=None)
    tabla1_2 = np.transpose(tabla1)
  
listeV = tabla1_2[0]
listEELS = tabla1_2[1]

listeV_correct = []
listEELS_correct = []

for k in range(len(listEELS)):
    if listEELS[k]>0: 
        listeV_correct.append(listeV[k])
        listEELS_correct.append(listEELS[k])

peaks, _ = find_peaks(listEELS , height=np.max(listEELS)*0.1)
listeV  = np.array(listeV )
listEELS  = np.array(listEELS )

sorted_indices = np.argsort(listeV)
listeV = listeV[sorted_indices]
listEELS = listEELS[sorted_indices]

plt.plot(listeV, listEELS,'.-',lw=lw ,color = list_colors[1] )
# plt.plot(listeV[peaks], listEELS[peaks], "x",color = list_colors[LL])

## IMPORTANT: sort the y-values from minimum to maximum --> mode = -1 is the highest (last one), mode = -2 is the previous one, and so on, .. 
sorted_index = np.argsort( listEELS[peaks])
# listy_peaks_sorted = np.sort(listy_peaks) 
# listx_peaks_sorted = list_energy0[peaks[sorted_index]]
list_index_peaks = peaks[sorted_index]
listx_peaks_sorted = listeV[peaks[sorted_index]]
     
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.ylim(-1e-4*1.3,np.max(listEELS)*1.1)
# plt.yscale('log')

plt.xticks([0.25,0.5,0.75,1,1.25,1.5], ['','0.5','','1.0','','1.5'])
plt.yticks(np.arange(0,0.0008,0.0001),['0.0000','','0.0002','','0.0004','','0.0006',''])
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1)
plt.tight_layout()
plt.savefig('EELS_paper_all_modes.png' ,bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)
plt.savefig('EELS_paper_all_modes.pdf' ,bbox_inches='tight',pad_inches = 0.01, format='pdf', dpi=dpi)
plt.show()

#%%

print('Plot EELS with N = %i and L = %i for highest mode' %(N0,L0))

mode = 1
kk = mode - 1 
   

plt.figure(figsize=tamfig)
plt.title(title,fontsize=tamtitle)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
name3 = 'EELS_comsol2_spectrum_N%i_W%inm_h%inm_L%inm_Ee%ikeV_zoom.dat' %(N0,W,h,L0,Ee_electron_keV) ## virtual geometry
#name3 = 'EELS_comsol2_spectrum_N%i_W%inm_h%inm_L%inm_Ee%ikeV_mode%i.dat' %(N0,W,h,L,Ee_electron_keV,mode) ## virtual geometry
## the label "EELS_comsol2.." means the permittivity of Si is 12.25 + i*0.2

#name3 = 'EELS_N%i_W%inm_h%inm_L%inm_R%inm_Ee%ikeV.dat' %(N,W,h,L,R,Ee_electron_keV) ## virtual geometry

try:
    tabla1 = np.loadtxt(name3, delimiter=' ',dtype=None)
    tabla1_2 = np.transpose(tabla1)
except (ValueError, FileNotFoundError):
    tabla1 = np.loadtxt(name3, delimiter='\t',dtype=None)
    tabla1_2 = np.transpose(tabla1)
  
listeV = tabla1_2[0]
listEELS = tabla1_2[1]

listeV_correct = []
listEELS_correct = []

for k in range(len(listEELS)):
    if listEELS[k]>0: 
        listeV_correct.append(listeV[k])
        listEELS_correct.append(listEELS[k])


listeV  = np.array(listeV )
listEELS  = np.array(listEELS )

sorted_indices = np.argsort(listeV)
listeV = listeV[sorted_indices]
listEELS = listEELS[sorted_indices]

peaks, _ = find_peaks(listEELS , height=np.max(listEELS)*0.1)
              
plt.plot(listeV, listEELS,'.-',lw=lw,color = list_colors[1])
# plt.plot(listeV[peaks], listEELS[peaks], "x",color = list_colors[LL])

## IMPORTANT: sort the y-values from minimum to maximum --> mode = -1 is the highest (last one), mode = -2 is the previous one, and so on, .. 
sorted_index = np.argsort( listEELS[peaks])
 
# listx_peaks_sorted = list_energy0[peaks[sorted_index]]
list_index_peaks = peaks[sorted_index]
listx_peaks_sorted = listeV[peaks[sorted_index]]
listy_peaks_sorted = listEELS[peaks[sorted_index]]
maximum = peaks[sorted_index][-1]

if width_height == [W,h]:
    plt.xlim([x_left_limits[kk],x_right_limits[kk]])
if mode == 1: 
    plt.xticks(xticks1)
elif mode == 2: 
    plt.xticks(xticks2)
else:
    plt.xticks(xticks3)

if mode == 1:
    plt.plot(listeV[maximum], listEELS[maximum], 'o', color = 'black')

plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.ylim(-1e-4*1.3,np.max(listEELS)*1.1)
plt.yticks(np.arange(0,0.0008,0.0001),['0.0000','','0.0002','','0.0004','','0.0006',''])
# plt.yscale('log')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1)
plt.tight_layout()
plt.savefig('EELS_paper_mode%i.png'%(mode),bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)
plt.savefig('EELS_paper_mode%i.pdf'%(mode),bbox_inches='tight',pad_inches = 0.01, format='pdf', dpi=dpi)
plt.show()



