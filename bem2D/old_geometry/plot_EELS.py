
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila
plot the EELS of a rectangular waveguide
"""
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline,interp1d
import matplotlib as mpl

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

listN = [300,400]
listN = [400]

w=400 # nm
h=120 ## nm
s=w*0.1  ## nm. fixed to 10% of w
ze=50 ## nm
w_microns = w*1e-3

Ee_electron_keV = 200
Ee_electron = Ee_electron_keV*1e3
me_c2_eV = 510998.95069  ## me*c**2 in eV
beta = np.sqrt( 1- (1 + Ee_electron/me_c2_eV)**(-2) )  ## beta = v/c
gamma_e = 1/np.sqrt(1-epsi1*beta**2) ## gamma lorentz

energy_mode = 0.7
omegac_mode = energy_mode/(hb*c)
Lp = beta*gamma_e/(2*omegac_mode)

print('Plot EELS at the position of the electron vs energy')

labelx=r'$\omega W/c$'
labely='EELS from BEM2D'
title = r'$W$ = %i nm, $h$ = %i nm, $b$ = %i nm, $E_{\rm e}$ = %i keV' %(w,h,ze,Ee_electron_keV)
Nint = 50
plt.figure(figsize=tamfig)
plt.title(title,fontsize=tamtitle)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
for j in range(len(listN)):
    
    N = listN[j]
    
    os.chdir(path_data)
    
    name1 = 'EELS_N%i_a%inm_h%inm_s%inm_ze%inm.dat' %(N,w,h,s,ze)
    name1 = 'EELS_epsi12_N%i_a%inm_h%inm_ze%inm_Ee%ikeV.dat' %(N,w,h,ze,Ee_electron_keV)
    tabla1 = np.loadtxt(name1, delimiter=' ',dtype=None)
    tabla1_2 = np.transpose(tabla1)
    

    # The output file contains five columns: energy/wavelength, 
    # q, x, y, and EELS probability per electron per nm of path length and  
    # per eV of energy-loss range
    
    # q=w/v for the electron
    
    listeV = tabla1_2[0]
    listq = tabla1_2[1]
    listx = tabla1_2[2]
    listy = tabla1_2[3] ## this is actually z in our coordinates
    listEELS = tabla1_2[4]
    
    listeV_correct = []
    listEELS_correct = []
    for k in range(len(listEELS)):
        if listEELS[k]>0: 
            listeV_correct.append(listeV[k])
            listEELS_correct.append(listEELS[k])
    
    listeV_correct_2 = np.array(listeV_correct)*w_microns/aux
    plt.plot(listeV_correct_2  ,listEELS_correct   ,'.',lw =lw, color = list_colors[j], label = r'$N$ = %i' %(N))
    
    x_int = np.linspace(np.min(listeV),np.max(listeV),int(len(listeV_correct)*Nint))
    try:
        cs = interp1d(listeV_correct, listEELS_correct)
     
        # plt.plot(x_int,cs(x_int) )
    except ValueError:
        continue
    

    name2 = 'EELS_zoom_epsi12_N%i_a%inm_h%inm_ze%inm_Ee%ikeV.dat' %(N,w,h,ze,Ee_electron_keV)
    
    try: 
        print('adding zoom for N = %i' %(N))
        tabla2 = np.loadtxt(name2, delimiter=' ',dtype=None)
        tabla2_2 = np.transpose(tabla2)
        
        listeV_zoom = tabla2_2[0]
        listEELS_zoom = tabla2_2[4]
        
        listeV_correct_zoom  = []
        listEELS_correct_zoom  = []
        for k in range(len(listEELS_zoom)):
            if listEELS_zoom[k]>0: 
                listeV_correct_zoom.append(listeV_zoom[k])
                listEELS_correct_zoom.append(listEELS_zoom[k])
        
        listeV_correct_zoom_2 = np.array(listeV_correct_zoom)*w_microns/aux
        plt.plot(listeV_correct_zoom_2,listEELS_correct_zoom,'.-' ,lw =lw, color = list_colors[j+1] )
 
    except IOError: 
        continue
    

plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1)
# plt.yscale('log')
# plt.xticks(np.arange(0,4,0.5))
# plt.xlim(0,4.2)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.tight_layout()
os.chdir(path_data)
plt.savefig( 'EELS_convergency.png',bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)
plt.show()


#%%

tamfig = [3, 4]

tabla1 = np.loadtxt('EELS.dat', delimiter=' ',dtype=None)
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

print('Plot EELS along z for a fixed energy (choose a peak)')
list_energy_eV=[0.68,1.11,1.49]
#list_energy_eV=np.arange(0.1,1.6,0.1)
N=300
title = r'$w$ = %i nm, $h$ = %i nm, $s$ = %i nm, $N$ = %i, $E_{\rm e}$ = %i keV' %(w,h,s,N,Ee_electron_keV)
labelx='z coordinate (nm)'
labely='EELS from BEM2D'
Nint = 50
plt.figure(figsize=tamfig)
plt.title(title,fontsize=tamtitle)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
for energy_eV in list_energy_eV:
    os.chdir(path_data)
    name1 = 'EELS_along_z_N%i_a%inm_h%inm_energy%.2feV.dat' %(N,w,h,energy_eV)
 
    
    tabla1 = np.loadtxt(name1, delimiter=' ',dtype=None)
    tabla1_2 = np.transpose(tabla1)
    
    # The output file contains five columns: energy/wavelength, 
    # q, x, y, and EELS probability per electron per nm of path length and  
    # per eV of energy-loss range
    
    # q=w/v for the electron
    
    listeV = tabla1_2[0]
    listq = tabla1_2[1]
    listx = tabla1_2[2]
    listy = tabla1_2[3] ## this is actually z in our coordinates
    listEELS = tabla1_2[4]

     
    # plt.plot(x_int,cs(x_int), label = r'$N$ = %i' %(N))
    plt.plot(listy,listEELS,'.-', label = r'$\hbar\omega$ = %.2f eV' %(energy_eV))

plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1)
# plt.yscale('log')
# plt.xlim([0,3.1])
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.tight_layout()
os.chdir(path_data)
plt.savefig( 'EELS_along_z.png',bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)
plt.show()
