
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila
plot the EELS of a rectangular waveguide
on top of substrate
"""
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline,interp1d
import matplotlib as mpl
from scipy.signal import find_peaks

path_basic = os.getcwd()
os.chdir(path_basic)

from global_constants import constants 
hb,c,alpha,me_c2_eV = constants()
aux = hb*c
epsi1, epsi3 = 1, 1

path_data = os.path.join(path_basic, 'datfiles_EELS')

list_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
               '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

plot_EELS_map = 0
 
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
lw = 0.8

deltax,deltay = 3,8

values = tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad,deltax,deltay

# //     N   = the total number of parametrization points is N+50
# //     w   = total width of wg     || x
# //     h   = total thickness of wg  || y
# //     s   = rounding radius of wg (s<t/2)
# //     ze  = impact parameter is h/2+ze. ze distance between electron and wg surface  

 

W=500 # nm
h=200 ## nm


listN = [500,700,800,1000]
listN = [500,700]
listL = [500,1000]
N0=500
L0=500


energy=1

R=500 ## virtual geometry

Ee_electron_keV = 200
Ee_electron = Ee_electron_keV*1e3
me_c2_eV = 510998.95069  ## me*c**2 in eV
beta = np.sqrt( 1- (1 + Ee_electron/me_c2_eV)**(-2) )  ## beta = v/c
gamma_e = 1/np.sqrt(1-epsi1*beta**2) ## gamma lorentz


xticks=np.arange(0.5,3.5,0.5)

labelx=r'$\omega W/c$'
labely='EELS from BEM2D'
title = r'$W$ = %i nm, $h$ = %i nm, $E_{\rm e}$ = %i keV' %(W,h,Ee_electron_keV)
Nint = 50

#%%

print('Plot EELS map')
if plot_EELS_map == 1:
    N=N0
    L=L0
    tamfig = [5.5, 4]
    
    name1 = 'EELS_along_xz_no_virtual_N%i_W%inm_h%inm_L%inm_Ee%ikeV_energy%ieV.dat' %(N,W,h,L,Ee_electron_keV,energy)
    name2 = 'EELS_along_xz_N%i_W%inm_h%inm_L%inm_R%inm_Ee%ikeV_energy%ieV.dat' %(N,W,h,L,R,Ee_electron_keV,energy) ## virtual geometry
    
    name = name1
    os.chdir(path_data)
    tabla1 = np.loadtxt(name, delimiter=' ',dtype=None)
    tabla1_2 = np.transpose(tabla1)
      
    listeV = tabla1_2[0]
    listq = tabla1_2[1]
    listx = tabla1_2[2]
    listy = tabla1_2[3] ## this is actually z in our coordinates
    # listReEx = tabla1_2[4]
    # listImEx = tabla1_2[5]
    # listReEy = tabla1_2[6]
    # listImEy = tabla1_2[7]
    # listReEz = tabla1_2[8]
    # listImEz = tabla1_2[9]
    listEELS = tabla1_2[4]
    
     
    n_color = 21
    vmin1, vmax1 = np.nanmin(listEELS), np.nanmax(listEELS)
    cmap = plt.cm.hot  # define the colormap
    bounds =   np.logspace( -14, np.log10(vmax1) , n_color)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    
    #
    delta = 0
    labelx='x (nm)'
    labely='z (nm)'
    labelz1 = 'region'
    limits = [np.min(listx) -delta, np.max(listx)+ delta, np.min(listy) -delta, np.max(listy)+ delta]
    
    plt.figure(figsize=tamfig)
    # plt.title(title,fontsize=tamtitle)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    
    tpc = plt.tripcolor(listx,listy, listEELS  )
    plt.colorbar(tpc)
    plt.title('EELS',fontsize=tamletra)
    # plt.scatter(listx,listy,c=listEELS)
    plt.xlim(limits[0],limits[1])
    plt.ylim(limits[2],limits[3])
    # plt.plot(listx,np.ones(len(listx))*ze_2,'--',color='black')
    plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
    # plt.set_aspect(abs((np.min(listR)- np.max(listR))/( np.min(listz)-np.max(listz)))*ratio)
    plt.tight_layout()
    
    # plt.axis('scaled')
    os.chdir(path_data)
    plt.savefig( 'EELS_map.png',bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)
    plt.show()


#%%

print('Plot EELS vs N with L=%i'%(L0))

L=L0
plt.figure(figsize=tamfig)
plt.title(title,fontsize=tamtitle)
plt.xlabel(r'$\hbar\omega$ (eV)',fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(r'EELS with L = %i nm'%(L),fontsize=tamletra,labelpad =labelpady)
for N in listN:
    name3 = 'EELS_N%i_W%inm_h%inm_L%inm_Ee%ikeV.dat' %(N,W,h,L,Ee_electron_keV) ## virtual geometry
    name3 = 'EELS_comsol_spectrum_N%i_W%inm_h%inm_L%inm_Ee%ikeV.dat' %(N,W,h,L,Ee_electron_keV) ## virtual geometry
    #name3 = 'EELS_N%i_W%inm_h%inm_L%inm_R%inm_Ee%ikeV.dat' %(N,W,h,L,R,Ee_electron_keV) ## virtual geometry
    os.chdir(path_data)
    try:
        tabla1 = np.loadtxt(name3, delimiter=' ',dtype=None)
        tabla1_2 = np.transpose(tabla1)
    except ValueError:
        tabla1 = np.loadtxt(name3, delimiter='\t',dtype=None)
        tabla1_2 = np.transpose(tabla1)
      
    if len(tabla1_2)>2:
        listeV = tabla1_2[0]
        listq = tabla1_2[1]
        listx = tabla1_2[2]
        listy = tabla1_2[3] ## this is actually z in our coordinates
        # listReEx = tabla1_2[4]
        # listImEx = tabla1_2[5]
        # listReEy = tabla1_2[6]
        # listImEy = tabla1_2[7]
        # listReEz = tabla1_2[8]
        # listImEz = tabla1_2[9]
        listEELS = tabla1_2[4]
    else:
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
                          
    
    plt.plot(listeV, listEELS,'.-',lw=lw,label = r'$N = %i$' %(N))
    plt.plot(listeV[peaks], listEELS[peaks], "x")


    ## IMPORTANT: sort the y-values from minimum to maximum --> mode = -1 is the highest (last one), mode = -2 is the previous one, and so on, .. 
    sorted_index = np.argsort( listEELS[peaks])
    # listy_peaks_sorted = np.sort(listy_peaks) 
    # listx_peaks_sorted = list_energy0[peaks[sorted_index]]
    list_index_peaks = peaks[sorted_index]

    listx_peaks_sorted = listeV[peaks[sorted_index]]
    
    print(N, W,h,listx_peaks_sorted)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.xticks(xticks)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1)
plt.tight_layout()
plt.savefig( 'EELS_convergency_L%inm.png'%(L),bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)
plt.show()

#%%

print('Plot EELS vs L (length of the virtual geometry) with N=%i'%(N0))

N=N0
plt.figure(figsize=tamfig)
plt.title(title,fontsize=tamtitle)
plt.xlabel(r'$\hbar\omega$ (eV)',fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(r'EELS with N = %i nm'%(N),fontsize=tamletra,labelpad =labelpady)
for L in listL:
    name3 = 'EELS_N%i_W%inm_h%inm_L%inm_Ee%ikeV.dat' %(N,W,h,L,Ee_electron_keV) ## virtual geometry
    name3 = 'EELS_comsol_spectrum_N%i_W%inm_h%inm_L%inm_Ee%ikeV.dat' %(N,W,h,L,Ee_electron_keV) ## virtual geometry
    os.chdir(path_data)
    try:
        tabla1 = np.loadtxt(name3, delimiter=' ',dtype=None)
        tabla1_2 = np.transpose(tabla1)
    except ValueError:
        tabla1 = np.loadtxt(name3, delimiter='\t',dtype=None)
        tabla1_2 = np.transpose(tabla1)
    
    if len(tabla1_2)>2:
        listeV = tabla1_2[0]
        listq = tabla1_2[1]
        listx = tabla1_2[2]
        listy = tabla1_2[3] ## this is actually z in our coordinates
        # listReEx = tabla1_2[4]
        # listImEx = tabla1_2[5]
        # listReEy = tabla1_2[6]
        # listImEy = tabla1_2[7]
        # listReEz = tabla1_2[8]
        # listImEz = tabla1_2[9]
        listEELS = tabla1_2[4]
    else:
        listeV = tabla1_2[0]
        listEELS = tabla1_2[1]

    listeV_correct = []
    listEELS_correct = []
    for k in range(len(listEELS)):
        if listEELS[k]>0: 
            listeV_correct.append(listeV[k])
            listEELS_correct.append(listEELS[k])
            
    plt.plot(listeV, listEELS,'.-',lw=lw,label = r'$L = %i$ nm' %(L))
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.xticks(xticks)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1)
plt.tight_layout()
plt.savefig('EELS_convergency_virtual_geometry_N%i.png'%(N),bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)
plt.show()


#%% double waveguides 
"""
print('Plot EELS vs H (height of the second rectangular waveguide-substrate)')
N=300
listH = [10,500,1000]
plt.figure(figsize=tamfig)
# plt.title(title,fontsize=tamtitle)
plt.xlabel(r'$\hbar\omega$ (eV)',fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(r'EELS with virtual geometry',fontsize=tamletra,labelpad =labelpady)
for H in listH:
    name3 = 'EELS_N%i_W%inm_h%inm_L%inm_H%inm_Ee%ikeV.dat' %(N,W,h,L,H,Ee_electron_keV) ## virtual geometry
    os.chdir(path_data)
    tabla1 = np.loadtxt(name3, delimiter=' ',dtype=None)
    tabla1_2 = np.transpose(tabla1)
      
    
    listeV = tabla1_2[0]
    listq = tabla1_2[1]
    listx = tabla1_2[2]
    listy = tabla1_2[3] ## this is actually z in our coordinates
    # listReEx = tabla1_2[4]
    # listImEx = tabla1_2[5]
    # listReEy = tabla1_2[6]
    # listImEy = tabla1_2[7]
    # listReEz = tabla1_2[8]
    # listImEz = tabla1_2[9]
    listEELS = tabla1_2[4]

    listeV_correct = []
    listEELS_correct = []
    for k in range(len(listEELS)):
        if listEELS[k]>0: 
            listeV_correct.append(listeV[k])
            listEELS_correct.append(listEELS[k])
            
    plt.plot(listeV, listEELS,'.-',label = r'$H = %i$ nm' %(H))
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1)
plt.tight_layout()
plt.savefig( 'EELS_convergency_double_waveguide.png',bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)
plt.show()

"""



