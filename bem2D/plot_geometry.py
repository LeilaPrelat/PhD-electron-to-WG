
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila
plot the geometry of a rectangular waveguide
"""
import numpy as np
import os
import matplotlib.pyplot as plt

path_basic = os.getcwd()
path_data = os.path.join(path_basic, 'geometry')
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

tabla1 = np.loadtxt('geometry.dat', delimiter=' ',dtype=None)
tabla1_2 = np.transpose(tabla1)

listR = tabla1_2[0]
listz = tabla1_2[1]
list_region = tabla1_2[2]

#%%

zcolor=[]
for zz in list_region: 
    if zz==22: 
        zcolor.append('tab:orange')
    else:
        zcolor.append('tab:green')
#        print('orange')

#%%
#
delta = 10
labelx='x (nm)'
labely='z (nm)'
labelz1 = 'region'
limits = [np.min(listR) -delta, np.max(listR)+ delta, np.min(listz) -delta, np.max(listz)+ delta]

plt.figure(figsize=tamfig)
# plt.title(title,fontsize=tamtitle)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
plt.scatter(listR,listz,c=zcolor)
plt.xlim(limits[0],limits[1])
plt.ylim(limits[2],limits[3])
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
# plt.set_aspect(abs((np.min(listR)- np.max(listR))/( np.min(listz)-np.max(listz)))*ratio)
plt.tight_layout()
# plt.axis('scaled')
os.chdir(path_data)
plt.savefig( 'geometry.png',bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)
plt.show()

 