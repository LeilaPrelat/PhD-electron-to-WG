
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila
EELS
see paper#228 Eqs. 3
see paper#149 Eqs. 25

dEELS/dy
"""
import numpy as np
import matplotlib.pyplot as plt
from global_constants import constants 
import matplotlib as mpl
import matplotlib.ticker as ticker
import os
from EELS import EELS_no_QE, EELS_QE, EELS_integrated_over_kx_no_QE
from permittivity_epsilon import epsilon as epsilon2

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'$10^{{{}}}$'.format(b)

formatt = ticker.FuncFormatter(fmt)
    
create_data = 1      ## run data for the color maps 
if_real_material = 1 ## if epsilon_2 is constant or not. if is real material, epsi2 = epsilon(omega)

if if_real_material == 1: 
    real_units = 1 ## Gamma in real units or normalized by c (then Gamma dimensionless)
else:
    real_units = 0
    
zoom = 0 # zoom in energy (if = 1 better definition of the modes in the dispersion relation) 

if if_real_material == 1:
    
    label_png = '_real'
    material = 'Si'     ## default
    # material = 'Ge'
else:
    material = 'Si'
    
    if material == 'Si':
        Re_epsi2 = 2.13
        Re_epsi2 = 15
    elif material == 'Ge':
        Re_epsi2 = 17.38
    else:
        Re_epsi2 = 12    ## default if if_real_material == 0
    label_png = ''
    Im_epsi2 = 1e-1  ## value of EELS depends on this Im_epsi2

delta = 1e-1  ## extra loss for the imaginary part of the permittivity
pwd = os.path.dirname(__file__) 
path_save =  os.path.join(pwd,'plots_EELS_permittivity_extra_loss_%.2f'%(delta))

#%%
hb,c,alpha,me_c2_eV = constants()
aux = hb*c
epsi1, epsi3 = 1, 1

#omegac_WG = (np.pi/d)/np.sqrt(np.real(epsilon2)-1) ## omega_WG/c

d_microns = 0.2
d = d_microns
L = 1 # lenght of propagation in microns 
# propagation<infty if im(epsi2)!=0

## list of electron energies from jga notes 2025-04-30 ##
ind = 2
list_Ee_electron = [30 , 100 , 200]   ## keV
Ee_electron_keV = list_Ee_electron[ind]
Ee_electron = Ee_electron_keV*1e3
label_Ee = '_Ee%i' %(ind+1)

beta = np.sqrt( 1- (1 + Ee_electron/me_c2_eV)**(-2) )  ## beta = v/c
gamma_e = 1/np.sqrt(1-epsi1*beta**2)

N = 150
N = 100
list_ze_nm =  np.linspace(0.1,200,N)
if ind == 1:
    list_ze_nm =  np.logspace(-1,np.log10(200),N)
    list_ze_nm =  np.linspace(0.1,200,N)
elif ind == 2:
    list_ze_nm =  np.logspace(-1,np.log10(250),N)
    list_ze_nm =  np.linspace(0.1,250,N)

    
energy_0 = 1
energy_0 = 10
omegac_0 = energy_0/aux
ze_0 = 10*1e-3 ## microns 
ze_0 = 25*1e-3 ## microns 
ze_0 = 50*1e-3 ## microns 
ze0_nm = ze_0*1e3

if if_real_material == 1:
    epsi2 = epsilon2(energy_0,delta,material) 
    Re_epsi2 = np.real(epsi2)
else:
    epsi2 = Re_epsi2 +  1j*Im_epsi2

# delta_num = (epsi1 - epsi2)*(epsi3 - epsi2)
# delta_den = (epsi1 + epsi2)*(epsi3 + epsi2)
# lambda_p_microns = 4*np.pi*d*(np.log(delta_num/delta_den))**(-1)
# kp_microns = 2*np.pi/lambda_p_microns
# kp_microns_norm = kp_microns/omegac_0
# kx_microns_norm = np.sqrt(kp_microns_norm**2 - 1/beta**2 +1j*0  )
omegac_WG = np.real((np.pi/d)/np.sqrt(epsi2-1+1j*0)) ## omega_WG/c

limit1 = 1.001*(1/beta) ## integral from omega/v
limit2 = np.real(np.sqrt(epsi2)) ## inside light cone
    
list_u =  np.linspace(limit1,1.4*(1/beta),N) ## integration of EELS is from k_parallel = \omega/v
list_u =  np.linspace(limit1,50*omegac_0,N) ## integration of EELS is from k_parallel = \omega/v
# list_u =  np.linspace(kx_microns_norm*1e-5,kx_microns_norm*100,N) ## integration of EELS is from k_parallel = \omega/v
list_u =  np.linspace(omegac_0*1e-5,omegac_0*50,N) ## integration of EELS is from k_parallel = \omega/v
# list_u =  np.linspace(1e-6,30,N) ## integration of EELS is from k_parallel = \omega/v

if if_real_material == 1:
    if d_microns == 0.2:
        if zoom == 0:
            if ze0_nm == 10:
                list_energy_eV_2 = np.linspace(1e-1,10,N)   ##  absorption part
                list_k_parallel = np.linspace(1e-1,80,N)    ## 
                
                list_energy_eV_2 = np.linspace(1e-1,10,N)   ##  absorption part
                list_k_parallel = np.linspace(1e-1,80,N)    ## 
                
            else:
                list_energy_eV_2 = np.linspace(1e-1,10,N)   ##  absorption part
                list_k_parallel = np.linspace(1e-1,80,N)   
            
        else:
            
            list_energy_eV_2 = np.linspace(1e-2,2,N)   ##  absorption part
            list_k_parallel = np.linspace(1e-2,20,N)    ## 
    else:
        if zoom == 0:
            if ze0_nm == 10:
                list_energy_eV_2 = np.linspace(1e-1,10,N)   ##  absorption part
                list_k_parallel = np.linspace(1e-1,80,N)    ## 
        
                
            else:
                list_energy_eV_2 = np.linspace(1e-1,6,N)   ##  absorption part
                list_k_parallel = np.linspace(1e-1,60,N)    ## 
            
        else:
            
            list_energy_eV_2 = np.linspace(1e-2,4,N)   ##  absorption part
            list_k_parallel = np.linspace(1e-2,40,N)    ## 
else:
    omegac_WG = (np.pi/d)/np.sqrt(np.real(Re_epsi2)-1) ## omega_WG/c
    energy_WG = omegac_WG*aux
    
    list_energy_eV_2 = np.linspace(0.01,4*energy_WG,N) ## paper 370. Fig S6
    list_k_parallel = np.linspace(0.01,8*omegac_WG,N)  ## paper 370. Fig S6
    
#%%

#tamfig = [3.5, 3]
tamfig = [4.5, 3]
tamletra = 13
tamtitle  = 10
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

#%%

print('1-Verification: Plot the EELS without integration as a function of kx, energy = %i eV, ze = %i nm' %(energy_0,ze_0*1e3))

list_EELS_re = []
list_EELS_im = []
for kx_norm in list_u: 
    if if_real_material == 1:
        epsi2 = epsilon2(energy_0,delta,material) 
        
    u = np.sqrt(kx_norm**2 + 1/beta**2 )
    if u >= 1/beta:
        
        value = EELS_no_QE(energy_0,kx_norm,ze_0,d,beta,epsi2)
    else:
        value = 1e-9 ## small number  (cannot use zero because I use log scale)
    # value = Fresnel_coefficient(omegac,u,d,mode,Im_epsi2)
    list_EELS_re.append(np.real(value))
    list_EELS_im.append(np.imag(value))

#%% 

labelx = r'Parallel wave vector normalized $k_x/(\omega/c)$'
if real_units == 0: 
    labely = r'$\Gamma_{\parallel}(k_\parallel) c/L$'
else:
    labely = r'$\Gamma_{\parallel}(k_\parallel)/L$ (s/$\mu$m)'

if if_real_material == 0:
    title1 = r'EELS for $h = %.1f$ $\mu$m, Re($\epsilon_2$) = %i, Im($\epsilon_2$) = %.1e' %(d,Re_epsi2,Im_epsi2)
else:
    title1 = r'EELS for $h = %.1f$ $\mu$m, $\epsilon_2 = \epsilon_{%s}(\omega)$' %(d,material)
title2 = r'$\hbar\omega = %.1f$ eV, $z_{\text{e}} = %i$ nm, $v = %.2fc$,' %(energy_0,ze_0*1e3,beta)

if if_real_material == 0:
    list_EELS_re = np.array(list_EELS_re)/np.max(list_EELS_re)
    ## depends on the value of delta, so I normalized
    ## to the maximum because the value is arbitrary 
    labely = r'$\Gamma_{\parallel}(k_\parallel)/\Gamma_{\text{max}}$'

 
plt.figure(figsize=tamfig)
plt.title(title1 + '\n' + title2,fontsize=tamtitle)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
plt.plot(list_u, np.array(list_EELS_re) ,'.-',lw = 1.5 )
aux_eje_y = np.linspace(np.min(list_EELS_re),np.max(list_EELS_re),N)
plt.plot(np.ones(N)*omegac_WG ,aux_eje_y,'--',color='black',lw = 1.5 ) ## lower limit of integration of EELS
# plt.plot(list_u, np.array(list_EELS_im) ,'.-',lw = 1.5 )
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
#plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1) 
plt.show() 

#%%

print('2-Plot the EELS without integration as a function of (k_paralllel, energy), ze = %i nm, beta = %.2f'%(ze_0*1e3,beta))

def EELS_color_map(k_par,energy):
    omegac = energy/aux
    u = k_par/omegac
    if if_real_material == 1:
        epsi2 = epsilon2(energy,delta,material)    
    
    if u >= 1/beta: ## kx has to be positive. The formula from paper #149 ("EELS_no_QE") is defined for positive kx
                    ## 2025/04/15 comments. see #issue02
                    
        kx_norm_k = np.sqrt(u**2-(1/beta)**2)
        value = EELS_no_QE(energy,kx_norm_k,ze_0,d,beta,epsi2)
    else:
        value = 1e-9 ## small number  (cannot use zero because I use log scale)
    return np.real(value)

listx, listy = list_k_parallel, list_energy_eV_2 
X, Y = np.meshgrid(listx,listy)
f_EELS = np.vectorize(EELS_color_map)
Z_EELS = f_EELS(X, Y)

#%%

labelx = r'Parallel wave vector $k_\parallel$ (1/$\mu$m)'
# labelx = r'Wave vector $k_x$ (1/$\mu$m)'
labely = r'Electron energy loss $\hbar\omega$ (eV)'
labelz = r'$\text{d}\Gamma_{\parallel}(k_\parallel)c/\text{d}y$'

if if_real_material == 0:
    Z_EELS = np.array(Z_EELS)/np.max(Z_EELS)
    labelz = r'$\Gamma_{\parallel}(k_\parallel)$ (a.u.)'
    # units of EELS arbitrary because it is not integrated, is divided by L and multiplied by c
    # so i normalized by the maximum

saturation_number_inf = -5
saturation_number_up = -2
delta = 1e-8 ## vmin can be zero
limits = [np.nanmin(listx) , np.nanmax(listx),np.nanmin(listy) , np.nanmax(listy)]
cmap = plt.cm.hot  # define the colormap
n_color = 21
vmin1, vmax1 = np.nanmin(Z_EELS), np.nanmax(Z_EELS)
vmin = vmin1 + delta
vmax = vmax1 + delta
bounds =   np.logspace(  np.log10(vmin1), np.log10(vmax1) , n_color)
bounds =   np.logspace( saturation_number_inf, saturation_number_up , n_color)

norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

maxlog = int(np.ceil( np.log10( np.abs(vmax) )))
minlog = int(np.ceil( np.log10( np.abs(vmin) ))) 
ticks_z = [(10.0**x) for x in np.linspace(minlog,maxlog-1, int(np.sign(vmax1)*maxlog) -  int(np.sign(vmin1)*minlog) ) ]
ticks_z = [(10.0**x) for x in np.linspace(saturation_number_inf,saturation_number_up, int(saturation_number_up) -  int(saturation_number_inf) +1 ) ]
 

plt.figure(figsize=tamfig)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
if if_real_material == 0:
    im_EELS = plt.imshow(Z_EELS, extent = limits, cmap=cmap, aspect='auto', interpolation = 'bicubic',origin = 'lower'  ) 
    if material == 'generic':
        plt.xticks(np.arange(1,8,1))
        plt.yticks(np.arange(0.1,0.8,0.1))
else:
    im_EELS = plt.imshow(Z_EELS, extent = limits, cmap=cmap, aspect='auto', interpolation = 'bicubic',origin = 'lower'  ,norm = norm   ) 
    
    # plt.xticks(np.arange(1,7,1))
    # plt.yticks(np.arange(0.1,0.9,0.1))
    # plt.yticks(np.arange(0.1,3,0.5))
    
#cbar = plt.colorbar(im_EELS, fraction=0.2, pad=0.04, location = 'top'  ,format = formatt  )
cbar = plt.colorbar(im_EELS, fraction=0.2, pad=0.04 , format = formatt  )

plt.plot(np.array(listy)/(aux*beta),np.array(listy),'--',color = 'green')   ## electron velocity 
plt.plot(np.array(listx),np.array(listx)*aux,'-',color = 'green')            ## light cone 1. omega = k_par/c
if zoom == 1:
    if if_real_material==0:     
        
        plt.plot(np.array(listx),np.array(listx)*aux/np.sqrt(Re_epsi2),'-',color = 'green') ## light cone 2. omega = k_par/\sqrt{\epsilon_2}c
    else:
        list_re_epsi2 = []
        for y in listy: 
            list_re_epsi2.append(np.real(np.sqrt(epsilon2(y,delta,material))))
    
        plt.plot(np.array(listx),np.array(listx)*aux/np.array(list_re_epsi2),'-',color = 'green') ## light cone 2. omega = k_par/\sqrt{\epsilon_2}

cbar.set_ticks(ticks_z)
plt.xlim(np.nanmin(listx) , np.nanmax(listx))
plt.ylim(np.nanmin(listy) , np.nanmax(listy))
if zoom == 0:
    plt.xticks(np.arange(10,90,10))
# if zoom == 0:
#      plt.text(11, 4,r"$k_\parallel = \omega/v$",color = 'green', rotation=52,fontsize = tamletra)
#      plt.text(11, 8,r"$k_\parallel = \omega/c$",color = 'green', rotation=52,fontsize = tamletra)
# #     plt.xticks(np.arange(10,90,10),['','20','','40','','60','','80'])
# #     plt.yticks(np.arange(1,11,1),['','2','','4','','6','','8','','10'])

# else:
#       plt.text(0.5, 0.5,r"$k_\parallel = \omega/v$",color = 'green', rotation=52,fontsize = tamletra)
#       plt.text(0.5, 1,r"$k_\parallel = \omega/c$",color = 'green', rotation=52,fontsize = tamletra)
#     plt.xticks(np.arange(5,35,5),['5','10','15','20','25','30' ])
#     plt.yticks(np.arange(0.5,3.5,0.5),['0.5','1.0','1.5','2.0','2.5','3.0' ])
# plt.xticks(np.arange(10,60,10))
cbar.ax.set_title(labelz,fontsize=tamletra,pad = pad+8)
cbar.ax.tick_params(labelsize = tamnum, width=0.1, direction="in",which = 'both', length = 2,pad = pad)
data_figure = title1 + r', $z_{\text{e}} = %i$ nm, $\beta$ = %.2f' %(ze_0*1e3,beta)
#plt.title(data_figure,fontsize=tamtitle)
os.chdir(path_save)
label_figure ='disp_relation_ze%inm_beta%.2f_h%inm' %(ze_0*1e3,beta,d*1e3) + material + label_png  
np.savetxt("info_of_" + label_figure + ".txt", [data_figure], fmt='%s')
plt.savefig(label_figure + '.png', format='png',bbox_inches='tight',pad_inches = 0.04, dpi=dpi)  
plt.show()

#%%

print('3-Verification: Plot the EELS integrated over kx as function of energy for a fixed ze')

def EELS_integrated_over_k_par_color_map(energy,ze_nm):
    ze = ze_nm*1e-3
    if if_real_material == 1:
        epsi2 = epsilon2(energy,delta,material) 
    return np.real(EELS_integrated_over_kx_no_QE(energy,ze,d,beta,epsi2))

list_EELS_int_re = []
for energy in list_energy_eV_2: 
    value = EELS_integrated_over_k_par_color_map(energy,ze0_nm)
    list_EELS_int_re.append(np.real(value))
 
#%% 
tamfig = [4, 3]
labelx = r'Electron energy loss $\hbar\omega$ (eV)'
labely = r'$\text{d}\Gamma_{\parallel}c/\text{d}y$'

if if_real_material == 0:
    title1 = r'EELS for $h = %.1f$ $\mu$m, Re($\epsilon_2$) = %i, Im($\epsilon_2$) = %.1e' %(d,Re_epsi2,Im_epsi2)
else:
    title1 = r'EELS for $h = %.1f$ $\mu$m, $\epsilon_2 = \epsilon_{%s}(\omega)$' %(d,material)
title2 = r'$z_{\text{e}} = %i$ nm, $v = %.2fc$,' %(ze_0*1e3,beta)

if if_real_material == 0:
    list_EELS_re = np.array(list_EELS_re)/np.max(list_EELS_re)
    ## depends on the value of delta, so I normalized
    ## to the maximum because the value is arbitrary 
    labely = r'$\Gamma_{\parallel}(k_\parallel)/\Gamma_{\text{max}}$'

 
plt.figure(figsize=tamfig)
plt.title(title1 + '\n' + title2,fontsize=tamtitle)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
plt.plot(list_energy_eV_2, np.array(list_EELS_int_re) ,'.-',lw = 1.5 )
# plt.plot(list_u, np.array(list_EELS_im) ,'.-',lw = 1.5 )
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.yticks(np.arange(0,0.015,0.0025))
#plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1) 
os.chdir(path_save)
data_figure = title2
total_label = label_png + label_Ee + '_zoom%i_ze%inm'  %(zoom,ze0_nm) 
label_figure = 'EELS_2D' + material + total_label
np.savetxt("info_of_" + label_figure + ".txt", [data_figure], fmt='%s')
plt.savefig(label_figure + '.png', format='png',bbox_inches='tight',pad_inches = 0.04, dpi=dpi)  
plt.show() 

#%%

os.chdir(pwd)
import time

start = time.time()

print('4-Plot the EELS integrated over kx as function of (energy,ze)')

listx2, listy2 = list_energy_eV_2 , list_ze_nm

total_label = label_png + label_Ee + '_zoom%i'  %(zoom) 
if create_data == 1:

    X, Y = np.meshgrid(listx2,listy2)
    f_EELS_2 = np.vectorize(EELS_integrated_over_k_par_color_map)
    Z_EELS_2 = f_EELS_2(X, Y)
    
    os.chdir(path_save)
#    header = title1 + r', $\beta$ = %.2f. Re(EELS) from paper 228 Eq. 3 divided by L/c and in Gaussian units (dimensionless) ' %(beta)
    header = title1 + r', $\beta$ = %.2f. Re(EELS) from paper 149 Eq. 25 divided by L/c and in Gaussian units (dimensionless) ' %(beta)
    np.savetxt('Z_EELS_norm_%s'%(material)  + total_label + '.txt' ,Z_EELS_2, fmt='%.10f', delimiter='\t', header = header, encoding=None)
    np.savetxt('list_x_energy_eV_%s'%(material)  + total_label + '.txt',listx2, fmt='%.10f', delimiter='\t', header = header, encoding=None)
    np.savetxt('list_y_ze_nm_%s'%(material)  + total_label+ '.txt' ,listy2, fmt='%.10f', delimiter='\t', header = header, encoding=None)

else:
    os.chdir(path_save)
    Z_EELS_2 = np.loadtxt('Z_EELS_norm_%s'%(material)  + total_label + '.txt', delimiter='\t', skiprows = 1)
    listx2 = np.loadtxt('list_x_energy_eV_%s'%(material)  + total_label + '.txt' , delimiter='\t', skiprows = 1)
    listy2 = np.loadtxt('list_y_ze_nm_%s'%(material)  + total_label + '.txt', delimiter='\t', skiprows = 1)
    
    # Z_EELS_2 = np.transpose(Z_EELS_2)
    # listx2 = np.transpose(listx2)
    # listy2 = np.transpose(listy2)
    
#%%
labelx = r'Electron energy loss $\hbar\omega$ (eV)'
labely = r'Electron-plane distance $z_{\rm e}$ (nm)'
labelz = r'$\text{d}\Gamma_{\parallel}c/\text{d}y$'

if if_real_material == 0:
    Z_EELS_2 = np.array(Z_EELS_2)/np.max(Z_EELS_2)
    labelz = r'$\Gamma_{\parallel}$ (a.u.)'


limits2 = [np.min(listx2) , np.max(listx2),np.min(listy2) , np.max(listy2)]
delta2 = 1e-2 ## vmin can be zero
delta2 = 0 ## vmin can be zero
cmap = plt.cm.hot  # define the colormap
n_color = 45
vmin12, vmax12 = np.nanmin(Z_EELS_2), np.nanmax(Z_EELS_2)
vmin2 = vmin12 + delta2
vmax2 = vmax12 + delta2
if if_real_material == 0:   
    bounds2 =   np.logspace(np.log10(vmin12+1e-4), np.log10(vmax12), n_color) 
else:
    bounds2 =   np.logspace(np.log10(vmin12 + 1e-3), np.log10(vmax12)  , n_color) 
    bounds2 =   np.logspace(np.log10(vmin2), np.log10(vmax2) , n_color) 
    bounds2 =   np.linspace(vmin12, vmax12 - delta2 , n_color) 
    # bounds =   np.logspace(np.log10(vmin1 +1e-2), np.log10(vmax1)  , n_color) 
norm2 = mpl.colors.BoundaryNorm(bounds2, cmap.N)

plt.figure(figsize=tamfig)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
im_EELS_2 = plt.imshow(Z_EELS_2, extent = limits2, cmap=cmap, aspect='auto', interpolation = 'bicubic',origin = 'lower' ,norm = norm2  ) 
cbar = plt.colorbar(im_EELS_2, fraction=0.046, pad=0.04   ,format = '%.3f')
cbar.ax.set_title(labelz,fontsize=tamletra)
cbar.ax.tick_params(labelsize = tamnum , width=0.1, length = 0,pad = 2)
# plt.xticks(np.arange(0,3.5,0.5))
# plt.yticks(np.arange(0,175,25))
# plt.xlim(np.min(listx2) , np.max(listx2))
# plt.ylim(np.min(listy2) , np.max(listy2))
plt.plot(listx2,np.ones(len(listx2))*ze0_nm)
plt.yscale('log')
data_figure = title1 + r', $\beta$ = %.2f' %(beta)
#plt.title(data_figure,fontsize=tamtitle)
os.chdir(path_save)
label_figure = 'EELS_' + material + total_label
np.savetxt("info_of_" + label_figure + ".txt", [data_figure], fmt='%s')
plt.savefig(label_figure + '.png', format='png',bbox_inches='tight',pad_inches = 0.04, dpi=dpi)  
plt.show()
        
 #%%   
 
time.sleep(1)
end = time.time()
if create_data == 1:
    print("Total runtime of the program is ", end - start ," seconds")

