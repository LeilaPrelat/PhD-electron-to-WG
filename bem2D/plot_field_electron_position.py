
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila
plot the geometry of a rectangular waveguide
"""
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl

path_basic = os.getcwd()
path_data = os.path.join(path_basic, 'datfiles_EELS')
os.chdir(path_data)
 
#%%

tamfig = [3, 4]
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

a=200 
h=400 
s=20
N=300
ze=50
ze_2 = ze + h/2

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

#%%

field = r'Re(Ez)'

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
plt.title(field,fontsize=tamletra)
# plt.scatter(listx,listy,c=listEELS)
plt.xlim(limits[0],limits[1])
plt.ylim(limits[2],limits[3])
# plt.plot(listx,np.ones(len(listx))*ze_2,'--',color='black')
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
# plt.set_aspect(abs((np.min(listR)- np.max(listR))/( np.min(listz)-np.max(listz)))*ratio)
plt.tight_layout()

# plt.axis('scaled')
os.chdir(path_data)
plt.savefig( 'field.png',bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)
plt.show()

 