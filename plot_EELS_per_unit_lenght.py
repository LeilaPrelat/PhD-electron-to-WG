
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
from EELS import EELS_no_QE,EELS_QE, EELS_integrated_over_k_par_no_QE
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
    
zoom = 1 # zoom to the dispersion relation

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

pwd = os.path.dirname(__file__) 
path_save =  os.path.join(pwd,'plots_EELS')

#%%
hb,c,alpha,me_c2_eV = constants()
aux = hb*c
epsi1, epsi3 = 1, 1

#omegac_WG = (np.pi/d)/np.sqrt(np.real(epsilon2)-1) ## omega_WG/c

d_microns = 0.1
d = d_microns
L = 1 # lenght of propagation in microns 
# propagation<infty if im(epsi2)!=0

    
## list of electron energies from jga notes 2025-04-30 ##
ind = 1
list_Ee_electron = [30 , 100 , 200]   ## keV
Ee_electron_keV = list_Ee_electron[ind]
Ee_electron = Ee_electron_keV*1e3
label_Ee = '_Ee%i' %(ind+1)

beta = np.sqrt( 1- (1 + Ee_electron/me_c2_eV)**(-2) )  ## beta = v/c
gamma_e = 1/np.sqrt(1-epsi1*beta**2)

N = 100
if zoom == 0:
    list_energy_eV = np.linspace(0.01,100,N)
else:
    list_energy_eV = np.linspace(0.01,4.5,N)
list_ze_nm =  np.linspace(0.1,200,N)

if ind == 1:
    list_ze_nm =  np.logspace(-1,np.log10(200),N)
elif ind == 2:
    list_ze_nm =  np.logspace(-1,np.log10(250),N)
if zoom == 1:
    list_energy_eV = np.logspace(-2,0,N)
else:
    list_energy_eV = np.logspace(-2,1,N)

energy_0 = 1

energy_0 = 5
omegac_0 = energy_0/aux
ze_0 = 10*1e-3 ## microns 
ze_0 = 25*1e-3 ## microns 
ze_0 = 10*1e-3 ## microns 

if if_real_material == 1:
    epsi2 = epsilon2(energy_0,material) 
    Re_epsi2 = np.real(epsi2)
else:
    epsi2 = Re_epsi2 +  1j*Im_epsi2

# delta_num = (epsi1 - epsi2)*(epsi3 - epsi2)
# delta_den = (epsi1 + epsi2)*(epsi3 + epsi2)
# lambda_p_microns = 4*np.pi*d*(np.log(delta_num/delta_den))**(-1)
# kp_microns = 2*np.pi/lambda_p_microns
# kp_microns_norm = kp_microns/omegac_0

limit1 = 1.001*(1/beta) ## integral from omega/v
limit2 = np.real(np.sqrt(epsi2)) ## inside light cone
    
list_u =  np.linspace(limit1,1.4*(1/beta),N) ## integration of EELS is from k_parallel = \omega/v
list_u =  np.linspace(limit1,50*omegac_0,N) ## integration of EELS is from k_parallel = \omega/v
list_u =  np.linspace(0.001*omegac_0,50*omegac_0,N) ## integration of EELS is from k_parallel = \omega/v
# list_u =  np.linspace(0.00001,40,N) ## integration of EELS is from k_parallel = \omega/v

#%%
 
tamfig = [4, 3]
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

print('1-Plot the EELS without integration as a function of k_parallel, energy = %i eV, ze = %i nm' %(energy_0,ze_0*1e3))

list_EELS_re = []
list_EELS_im = []
for u in list_u: 
    if if_real_material == 1:
        epsi2 = epsilon2(energy_0,material) 
    value = EELS_no_QE(energy_0,u,ze_0,d,beta,epsi2)
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
plt.plot(np.ones(N)*(1/beta)*1.001 ,aux_eje_y,'--',color='black',lw = 1.5 ) ## lower limit of integration of EELS
# plt.plot(list_u, np.array(list_EELS_im) ,'.-',lw = 1.5 )
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
#plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1) 
plt.show() 

#%%

## we want to excite all the modes with an electron velocity of almost 1, so we see the dispersion relation
## with the EELS 
beta0 = beta
print('2-Plot the EELS without integration as a function of (k_paralllel, energy), ze = %i nm, beta = %.2f'%(ze_0*1e3,beta0))

if if_real_material == 1:
    list_energy_eV_2 = np.linspace(0.01,0.8,int(N*2)) ## non absorption part
    list_energy_eV_2 = np.linspace(0.01,3,int(N*2))   ##  absorption part
    list_k_parallel = np.linspace(0.01,6,int(N*2))    ## 
    
    if zoom == 0:
        list_energy_eV_2 = np.linspace(0.001,10,int(N*2))   ##  absorption part
        list_k_parallel = np.linspace(0.001,50,int(N*2))    ## 
        
        list_energy_eV_2 = np.logspace(-1,np.log10(4),int(N*2))   ##  absorption part
        list_k_parallel = np.logspace(-1,np.log10(40),int(N*2))    ## 
        
        list_energy_eV_2 = np.linspace(1e-1,10,int(N*2))   ##  absorption part
        list_k_parallel = np.linspace(1e-1,40,int(N*2))    ## 
        
    else:
        list_energy_eV_2 = np.linspace(0.001,4,int(N*2))  ##  absorption part
        list_k_parallel = np.linspace(0.001,55,int(N*2))    ## 
        
        list_energy_eV_2 = np.logspace(-1,np.log10(10),int(N*2))   ##  absorption part
        list_k_parallel = np.logspace(-1,np.log10(60),int(N*2))    ## 
        
        list_energy_eV_2 = np.linspace(1e-2,10,int(N*2))   ##  absorption part
        list_k_parallel = np.linspace(1e-3,60,int(N*2))    ## 
        
else:
    omegac_WG = (np.pi/d)/np.sqrt(np.real(Re_epsi2)-1) ## omega_WG/c
    energy_WG = omegac_WG*aux
    
    list_energy_eV_2 = np.linspace(0.01,4*energy_WG,int(N*2)) ## paper 370. Fig S6
    list_k_parallel = np.linspace(0.01,8*omegac_WG,int(N*2))  ## paper 370. Fig S6

def EELS_color_map(k_parallel,energy):
    omegac = energy/aux
    u = k_parallel/omegac
    # print(beta)
    if if_real_material == 1:
        epsi2 = epsilon2(energy,material)         
    return np.real(EELS_no_QE(energy,u,ze_0,d,beta0,epsi2))

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

list_k_parallel_for_disp_relation = []
for j in range(len(listx)): 
    energy = listy[j]
    x = listx[j]
    omegac = energy/aux
    list_k_parallel_for_disp_relation.append(np.sqrt(x**2 +  (omegac/beta)**2))


saturation_number = -5
delta = 1e-8 ## vmin can be zero
limits = [np.nanmin(list_k_parallel_for_disp_relation) , np.nanmax(list_k_parallel_for_disp_relation),np.nanmin(listy) , np.nanmax(listy)]
cmap = plt.cm.hot  # define the colormap
n_color = 21
vmin1, vmax1 = np.nanmin(Z_EELS), np.nanmax(Z_EELS)
vmin = vmin1 + delta
vmax = vmax1 + delta
if if_real_material == 0:   
    bounds =   np.logspace(np.log10(vmin1+1e-4), np.log10(vmax1), n_color) 
else:
    bounds =   np.logspace(np.log10(vmin1 + 1e-3), np.log10(vmax1)  , n_color) 
    bounds =   np.logspace( saturation_number, np.log10(vmax1) , n_color) 

    # bounds =   np.logspace(-5, np.log10(vmax1) , n_color) 
    # bounds =   np.linspace(1e-16, vmax1 , n_color) 
    # bounds =   np.logspace(np.log10(vmin1 +1e-2), np.log10(vmax1)  , n_color) 
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

maxlog = int(np.ceil( np.log10( np.abs(vmax) )))
minlog = int(np.ceil( np.log10( np.abs(vmin) ))) 
# ticks_z = [(10.0**x) for x in np.linspace(saturation_number,maxlog-1, int(np.sign(vmax1)*maxlog) - int(saturation_number) ) ]
# ticks_z = [(10.0**x) for x in np.linspace(-8, 0, 9) ]

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
    
cbar = plt.colorbar(im_EELS, fraction=0.046, pad=0.04, orientation = 'vertical'    )


plt.plot(np.array(listy)/(aux*beta0),np.array(listy),'-',color = 'green')   ## electron velocity 

plt.plot(np.array(list_k_parallel_for_disp_relation),np.array(list_k_parallel_for_disp_relation)*aux,'-',color = 'green')            ## light cone 1. omega = k_par/c
# if if_real_material==0:     
    
#     plt.plot(np.array(listx),np.array(listx)*aux/np.sqrt(Re_epsi2),'-',color = 'green') ## light cone 2. omega = k_par/\sqrt{\epsilon_2}c
# else:
#     list_re_epsi2 = []
#     for y in listy: 
#         list_re_epsi2.append(np.real(np.sqrt(epsilon2(y,material))))

#     plt.plot(np.array(listx),np.array(listx)*aux/np.array(list_re_epsi2),'-',color = 'green') ## light cone 2. omega = k_par/\sqrt{\epsilon_2}

# cbar.set_ticks(ticks_z)
plt.xlim(np.nanmin(list_k_parallel_for_disp_relation) , np.nanmax(list_k_parallel_for_disp_relation))
plt.ylim(np.nanmin(listy) , np.nanmax(listy))
# plt.xticks(np.arange(10,60,10))
# if zoom == 0:
#     plt.text(11, 4,r"$k_\parallel = \omega/v$",color = 'green', rotation=52,fontsize = tamletra)
#     plt.xticks(np.arange(10,90,10),['','20','','40','','60','','80'])
#     plt.yticks(np.arange(1,11,1),['','2','','4','','6','','8','','10'])

# else:
#     # plt.text(0.5, 0.5,r"$k_\parallel = \omega/v$",color = 'green', rotation=52,fontsize = tamletra)
#     plt.xticks(np.arange(5,35,5),['5','10','15','20','25','30' ])
#     plt.yticks(np.arange(0.5,3.5,0.5),['0.5','1.0','1.5','2.0','2.5','3.0' ])
# plt.xticks(np.arange(10,60,10))
cbar.ax.set_title(labelz,fontsize=tamletra)
cbar.ax.tick_params(labelsize = tamnum, width=0.1, length = 0,pad = 2)
data_figure = title1 + r', $z_{\text{e}} = %i$ nm, $\beta$ = %.2f' %(ze_0*1e3,beta0)
#plt.title(data_figure,fontsize=tamtitle)
os.chdir(path_save)
label_figure ='disp_relation_ze%inm_beta%.2f_' %(ze_0*1e3,beta0) + material + label_png  
np.savetxt("info_of_" + label_figure + ".txt", [data_figure], fmt='%s')
plt.savefig(label_figure + '.png', format='png',bbox_inches='tight',pad_inches = 0.04, dpi=dpi)  
plt.show()

#%%

os.chdir(pwd)
import time

start = time.time()

print('3-Plot the EELS integrated over k_parallel as function of (energy,ze)')

def EELS_integrated_over_k_par_color_map(energy,ze_nm):
    ze = ze_nm*1e-3
    if if_real_material == 1:
        epsi2 = epsilon2(energy,material) 
    return np.real(EELS_integrated_over_k_par_no_QE(energy,ze,d,beta,epsi2))

listx2, listy2 = list_energy_eV , list_ze_nm

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
delta = 1e-8 ## vmin can be zero
cmap = plt.cm.hot  # define the colormap
n_color = 21
vmin12, vmax12 = np.nanmin(Z_EELS_2), np.nanmax(Z_EELS_2)
vmin2 = vmin12 + delta
vmax2 = vmax12 + delta
if if_real_material == 0:   
    bounds2 =   np.logspace(np.log10(vmin12+1e-4), np.log10(vmax12), n_color) 
else:
    bounds2 =   np.logspace(np.log10(vmin12 + 1e-3), np.log10(vmax12)  , n_color) 
    bounds2 =   np.logspace(np.log10(vmin2), np.log10(vmax2) , n_color) 
    # bounds =   np.logspace(np.log10(vmin1 +1e-2), np.log10(vmax1)  , n_color) 
norm2 = mpl.colors.BoundaryNorm(bounds2, cmap.N)

plt.figure(figsize=tamfig)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
im_EELS_2 = plt.imshow(Z_EELS_2, extent = limits2, cmap=cmap, aspect='auto', interpolation = 'bicubic',origin = 'lower'   ) 
cbar = plt.colorbar(im_EELS_2, fraction=0.046, pad=0.04, orientation = 'vertical'  )
cbar.ax.set_title(labelz,fontsize=tamletra)
cbar.ax.tick_params(labelsize = tamnum , width=0.1, length = 0,pad = 2)
# plt.xticks(np.arange(0,3.5,0.5))
# plt.yticks(np.arange(0,175,25))
# plt.xlim(np.min(listx2) , np.max(listx2))
# plt.ylim(np.min(listy2) , np.max(listy2))
# plt.yscale('log')
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

print("Total runtime of the program is ", end - start ," seconds")

