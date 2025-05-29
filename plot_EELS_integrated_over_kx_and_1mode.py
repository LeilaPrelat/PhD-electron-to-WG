
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila
EELS
see paper#228 Eqs. 3
see paper#149 Eqs. 25
integrated over electron's trajectory
for real materials (\epsilon(\omega))
and over the frequency from \Delta\omega_TM : plot total EELS vs b (electron-plane distance)

using Leff(k_par) = L0*f(k_par)
"""
import numpy as np
import matplotlib.pyplot as plt
from global_constants import constants 
import os
from EELS import EELS_integrated_over_electron_trayectory
from EELS import Fresnel_coefficient
from permittivity_epsilon import epsilon as epsilon2
from mycolorpy import colorlist as mcp
from scipy.signal import find_peaks
    
label_png = '_real' 
material = 'Si'     ## default
# material = 'Ge'  EELS_integrated_over_electron_trayectory
zoom = 0
create_data = 1

delta = 1e-1 ## extra loss for the imaginary part of the permittivity
pwd = os.path.dirname(__file__) 
path_save =  os.path.join(pwd,'plots_EELS_permittivity_extra_loss_%.2f'%(delta))

#%%
hb,c,alpha,me_c2_eV = constants()
aux = hb*c
epsi1, epsi3 = 1, 1

d_microns = 0.2 # microns
d = d_microns
    
## list of electron energies from jga notes 2025-04-30 ##
ind = 1
list_Ee_electron = [30 , 100 , 200]   ## keV
Ee_electron_keV = list_Ee_electron[ind]
Ee_electron = Ee_electron_keV*1e3
label_Ee = '_d%inm_Ee%i' %(d*1e3,ind+1)

beta = np.sqrt( 1- (1 + Ee_electron/me_c2_eV)**(-2) )  ## beta = v/c
gamma_e = 1/np.sqrt(1-epsi1*beta**2)

N = 100
if d == 0.2:
    if zoom == 0:
        list_energy_eV = np.linspace(0.1,10,N) ## cutoff energy
    else:
        list_energy_eV = np.linspace(0.1,2,N) ## cutoff energy
else:
    if zoom == 0:
        list_energy_eV = np.linspace(0.1,20,N) ## cutoff energy
    else:
        list_energy_eV = np.linspace(0.1,4,N) ## cutoff energy

list_b_nm  = np.linspace(0,50,51)
list_b_nm  = np.linspace(0,30,31)
total_label = material + label_png + label_Ee  + 'zoom%i' %(zoom)
# list_b_nm = [0,10,50]
# list_b_nm = [0,10,30]
Nint = 50 ## N interpolation factor
title = r'EELS for $h = %.1f$ $\mu$m, $\epsilon_2 = \epsilon_{%s}(\omega)$, $v = %.2fc$' %(d,material,beta) 
plot_figure = 0
list_modes = [-1,-2,-3]
# list_modes = [-1 ]

#%%
 
tamfig = [3.45, 3]
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
color1 = mcp.gen_color(cmap="hot",n=len(list_b_nm)+1)

#%% 

def find_FWHW(mode):
    
    print('1-Plot the EELS integrated over k_par and the trajectory as a function of energy, for different b')
    
    list_EELS_re_tot = []
    list_EELS_im_tot = []
    for b_nm in list_b_nm:
        b = b_nm*1e-3
        list_EELS_re = []
        list_EELS_im = []
        for eV in list_energy_eV: 
            epsi2 = epsilon2(eV,delta,material) 
            value = EELS_integrated_over_electron_trayectory(eV,b,d,beta,epsi2)
            # value = Fresnel_coefficient(omegac,u,d,mode,Im_epsi2)
            list_EELS_re.append(np.real(value))
            list_EELS_im.append(np.imag(value))
            
        list_EELS_re_tot.append(list_EELS_re)
        list_EELS_im_tot.append(list_EELS_im)
        
    labelx = r'Electron energy loss $\hbar\omega$ (eV)'
    labely = r'$\Gamma_{\parallel}/L_0$ (fs/$\mu$m)'
    

    
    print('1b- Find the FWHM of the highest mode for each b to integrate over it later')
    
    if plot_figure == 1:
        plt.figure(figsize=tamfig)
        #plt.title(title,fontsize=tamtitle)
        plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
        plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
        
    x_FHM1_tot = []  ## lower limit of integration 
    x_FHM2_tot = []  ## upper limit of integration 
    for j in range(len(list_b_nm)):
        listy =  np.array(list_EELS_re_tot[j])*1e15 ## fentosecond
        
        
        ## interpolation
        x_int = np.linspace(np.min(list_energy_eV),np.max(list_energy_eV),int(len(list_energy_eV)*Nint))
        y_int = np.interp(x_int, list_energy_eV, listy)
        if plot_figure == 1:
            plt.plot(x_int, y_int ,'-',color = color1[j],lw = 1.5,label = r'$b = %i$ nm' %(list_b_nm[j]) )
        
        ## find peaks
        peaks, _ = find_peaks(y_int, height=0)
    
        sorted_listy = np.sort(y_int[peaks])        ## sort the y-values from minimum to maximum  --> mode = -1 is the highest (last one), mode = -2 is the previous one, and so on, .. 
        sorted_index = np.argsort(y_int[peaks]) 
        sorted_listx = x_int[peaks[sorted_index]]
        
            
        if plot_figure == 1:
            plt.plot(sorted_listx, sorted_listy, "x")
        
       
        y_max = sorted_listy[mode]
        x_max = sorted_listx[mode]
        
        ## moving the curve down, then the x value that crosses the FWHM are when the curve changes signs 
        ynew = np.array(y_int) - y_max/2
        aux_eje_y = np.linspace(np.min(listy),np.max(listy),10)
    
        for j in range(len(ynew)-1): 
            if np.sign(ynew[j]*ynew[j+1]) == -1: ## crossing the zero
                mean_x =  (x_int[j+1]+ x_int[j])/2
                print(j,mean_x)
                if x_int[j] < x_max:        
                    x_FHM1_tot.append(mean_x)                    
                    if plot_figure == 1:
                        plt.plot(np.ones(10)*mean_x, aux_eje_y, "--")
                else: 
                    x_FHM2_tot.append(mean_x)
                        
                    if plot_figure == 1:
                        plt.plot(np.ones(10)*mean_x, aux_eje_y, "-")
        
    if plot_figure == 1:
        plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
        plt.xticks(np.arange(0,12,2))
        plt.yticks(np.arange(0,30,5))
        #plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1) 
        # label_figure = 'EELS_tot_' + total_label
        # plt.xscale('log')
        os.chdir(path_save)
        # plt.savefig(label_figure + '.png', format='png',bbox_inches='tight',pad_inches = 0.04, dpi=dpi)  
        plt.show() 
    
    return list_b_nm, x_FHM1_tot, x_FHM2_tot

#%%

print('2-Plot the EELS integrated over k_par, over the trajectory, and over energy close to the highest peak for different b')

#EELS integrated over y-coordinate and over frequency
def EELS_double_integral_as_sum(lower_eV_limit,upper_eV_limit,delta_hbw,b,d,beta):
    """    
    Parameters
    ----------
    lower_eV_limit : lower limit of the integral over energy in eV (omega_f)
    upper_eV_limit : upper limit of the integral over energy in eV (omega_f)
    delta_hbw : discretization of the omega_f in eV
    b: minimum position of electron along z coordinate in microns
    (most close to the plane)
    d: thickness of the plane in microns
    beta: v/c
    Returns
    -------
    Re(EELS) from paper 149 Eq. 25
    divided by L0 and in Gaussian units
    (1/microns) integrated over
    electron trajectory Leff(k_par)
    integrate over k_parallel
    and over energy
    (see notes)
    """
    # epsi2 = epsilon(hbw,material)

    Integral = 0
    N = 150

    # list_eV = np.logspace(-1,upper_eV_limit,N)
    
    list_eV = np.arange(lower_eV_limit, upper_eV_limit + delta_hbw ,delta_hbw)
    N = int(len(list_eV))
    Nk = 150
    for k in range(N-1):
        eV = list_eV[k]
        omegac = eV/(aux)
        
        list_kx_norm_k = np.linspace(1e-4*omegac,30*omegac,Nk)
        
        delta_energy =  list_eV[k+1] - list_eV[k]
        
        for j in range(Nk-1):
            kx_norm_k = list_kx_norm_k[j]
            
            delta_qx = list_kx_norm_k[j+1] - list_kx_norm_k[j]
            
            qx = kx_norm_k
            u =  np.sqrt(kx_norm_k**2 + (1/beta)**2) ## u = k_parallel/k with kx variable and ky = \omega/v
        
            ## for python add +1j*0 inside kz as \sqrt{ .. + 1j*0}  
            
            argument =  np.sqrt(epsi1 - u**2  + 1j*0) ## kz = sqrt(k^2 - k_parallel^2)
            kz1 = argument if np.imag(argument)>0  else  - argument  
            
            # if np.imag(kz1) <= 0:
            #     kz1 = - kz1
            epsi2 =  epsilon2(eV,delta,material) ## epsilon2 depends on frequency so i cannot put it as argument
            
            ## integration variable u = k_par_over_omegac (dimensionless)
            r123_s =  Fresnel_coefficient(omegac,u,d,'s',epsi2)
            r123_p = Fresnel_coefficient(omegac,u,d,'p',epsi2)
        
            final_function =  kz1*np.exp(2*1j*kz1*omegac*b)*(r123_s*(qx*beta/kz1)**2 - r123_p/epsi1)/(u**(5/2)*np.sqrt(omegac))
            final_function_re =  np.real(final_function)
        
            gamma_e = 1/(np.sqrt(epsi1-beta**2))
            me_over_hb = 8.648539271254356*1e-9 ## seconds/microns^2
            
            aux2 = np.sqrt(c)*hb ## add hb because the integral is over energy and not over frequency ------------------------------------- IMPORTANT 
            factor_Gamma_norm_L0  = alpha*2*np.sqrt(gamma_e)*np.sqrt(me_over_hb)/(np.pi*beta*aux2) ## Gamma/L0 in unis of seconds/microns
            ## the sqrt(omega) comes from the fact that in the integral I am using dimensionless variables instead of k_par, I use k_par/k --> part with omega is inside function
         
            Integral0 = final_function_re*factor_Gamma_norm_L0
            
            Integral = Integral + Integral0*delta_qx*delta_energy
        
    return  Integral ,N

#%%

labely = r'$\Gamma_{\text{TM}_\nu}/L_0$ (1/$\mu$m)'
labelx = r'Electron-plane distance $b$ (nm)'

plt.figure(figsize=tamfig)
# plt.title(title,fontsize=tamtitle)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)

for mode in list_modes: ## two highest modes

    if create_data == 1:
        
        list_b_nm, x_FHM1_tot, x_FHM2_tot = find_FWHW(mode)

        
        delta_hbw = 0.005
        list_EELS_re_tot = []
        list_EELS_im_tot = []
        
        k = 0
        for b_nm in list_b_nm :
            b = b_nm*1e-3
         
            
            ev_initial = x_FHM1_tot[k]
            eV_final = x_FHM2_tot[k]
                
            value,Nvalue = EELS_double_integral_as_sum(ev_initial,eV_final,delta_hbw,b,d,beta)
            print(b_nm, k, value, Nvalue)
            # value = Fresnel_coefficient(omegac,u,d,mode,Im_epsi2)
         
            k = k + 1
            
            list_EELS_re_tot.append(np.real(value))
            list_EELS_im_tot.append(np.imag(value))
        
        os.chdir(path_save)
        header = title + r', $mode$ = %i. Re(EELS) from paper 149 Eq. 25 divided by L0, in Gaussian units, integrated over kx and trajectory Leff' %(mode)
        np.savetxt('Pj_list_b_nm' + total_label + '_mode%i.txt' %(mode), list_b_nm, fmt='%.10f', delimiter='\t', header = header, encoding=None)
        np.savetxt('Pj_eV_FHM1' + total_label + '_mode%i.txt' %(mode), x_FHM1_tot, fmt='%.10f', delimiter='\t', header = header, encoding=None)
        np.savetxt('Pj_eV_FHM2' + total_label+ '_mode%i.txt'%(mode), x_FHM2_tot, fmt='%.10f', delimiter='\t', header = header, encoding=None)
        np.savetxt('Pj_value' + total_label+ '_mode%i.txt'%(mode), list_EELS_re_tot, fmt='%.10f', delimiter='\t', header = header, encoding=None)
        
    else:
        os.chdir(path_save)
        list_b_nm = np.loadtxt('Pj_list_b_nm' + total_label + '_mode%i.txt' %(mode),  delimiter='\t', skiprows = 1, encoding=None)
        list_EELS_re_tot = np.loadtxt('Pj_value' + total_label+ '_mode%i.txt'%(mode),delimiter='\t', skiprows = 1, encoding=None)
    
    plt.plot(list_b_nm , np.array(list_EELS_re_tot) ,'-',lw = 1.5  ,label = r'$\nu = %i$' %(mode))
    plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
    label_figure = 'Pj_int_energy_' + total_label
    plt.xticks(np.arange(0,60,10))
    # plt.xscale('log')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1) 
    plt.savefig(label_figure + '.png', format='png',bbox_inches='tight',pad_inches = 0.045, dpi=dpi)  
    plt.show() 




