
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila
EELS
see paper#228 Eqs. 3
see paper#149 Eqs. 25
integrated over electron's trajectory
for real materials (\epsilon(\omega))
and over the frequency : total EELS
"""
import numpy as np
import matplotlib.pyplot as plt
from global_constants import constants 
import os
from EELS import Fresnel_coefficient
from scipy.integrate import dblquad
from permittivity_epsilon import epsilon as epsilon2

create_data = 1      ## run data for the color maps 
real_units = 1       ## Gamma in real units or normalized by c (then Gamma dimensionless)
    
label_png = '_real'
material = 'Si'   ## default
# material = 'Ge'  
zoom = 1

pwd = os.path.dirname(__file__) 
path_save =  os.path.join(pwd,'plots_EELS')

#%%
hb,c,alpha,me_c2_eV = constants()
aux = hb*c
epsi1, epsi3 = 1, 1

d_microns = 0.2 # microns
d = d_microns
    
## list of electron energies from jga notes 2025-04-30 ##
ind = 1
list_Ee_electron = [30, 100 , 200]   ## keV . dont use the first value 
Ee_electron_keV = list_Ee_electron[ind]
Ee_electron = Ee_electron_keV*1e3
label_Ee = '_Ee%i' %(ind+1)

beta = np.sqrt( 1- (1 + Ee_electron/me_c2_eV)**(-2) )  ## beta = v/c
gamma_e = 1/np.sqrt(1-epsi1*beta**2)

N = 100
if zoom == 0:
    list_upper_eV_limit = np.linspace(0.1,10,N) ## cutoff energy
else:
    list_upper_eV_limit = np.linspace(0.1,2,N) ## cutoff energy

list_b_nm = [10,50,80]

 
total_label = material + label_png + label_Ee  + 'zoom%i' %(zoom)

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

#EELS integrated over y-coordinate and over frequency
def EELS_integrated_over_electron_trayectory_and_energy(upper_eV_limit,b,d,beta):
    """    
    Parameters
    ----------
    upper_eV_limit : upper limit of the integral over energy in eV
    b: minimum position of electron along z coordinate in microns
    (most close to the plane)
    d: thickness of the plane in microns
    beta: v/c
    Returns
    -------
    Re(EELS) from paper 149 Eq. 25
    divided by L0 and in Gaussian units
    (1/microns) integrated over
    electron trajectory, over k_parallel
    and over energy
    (see notes)
    """
    # epsi2 = epsilon(hbw,material)


    
    u = lambda kx_norm_k : np.sqrt(kx_norm_k**2 + (1/beta)**2) ## u = k_parallel/k with kx variable and ky = \omega/v

    ## for python add +1j*0 inside kz as \sqrt{ .. + 1j*0}  
    
    argument = lambda kx_norm_k : np.sqrt(epsi1 - u(kx_norm_k)**2  + 1j*0) ## kz = sqrt(k^2 - k_parallel^2)
    kz1 = lambda kx_norm_k : argument(kx_norm_k) if np.imag(argument(kx_norm_k))>0  else  - argument(kx_norm_k)  
    
    # if np.imag(kz1) <= 0:
    #     kz1 = - kz1
    epsi2 = lambda eV: epsilon2(eV,material)
    omegac = lambda eV: eV/(aux)
    ## integration variable u = k_par_over_omegac (dimensionless)
    r123_s = lambda qx,eV:  Fresnel_coefficient(omegac(eV),u(qx),d,'s',epsi2(eV))
    r123_p = lambda qx,eV: Fresnel_coefficient(omegac(eV),u(qx),d,'p',epsi2(eV))

    final_function = lambda qx,eV : kz1(qx)*np.exp(2*1j*kz1(qx)*omegac(eV)*b)*(r123_s(qx,eV)*(qx*beta/kz1(qx))**2 - r123_p(qx,eV)/epsi1)/(u(qx)**(5/2)*np.sqrt(omegac(eV)))
    final_function_re = lambda qx,eV : np.real(final_function(qx,eV))

    gamma_e = 1/(np.sqrt(epsi1-beta**2))
    me_over_hb = 8.648539271254356*1e-9 ## seconds/microns^2
    
    aux2 = np.sqrt(c)*hb ## add hb because the integral is over energy and not over frequency ------------------------------------- IMPORTANT 
    factor_Gamma_norm_L0  = alpha*2*np.sqrt(gamma_e)*np.sqrt(me_over_hb)/(np.pi*beta*aux2) ## Gamma/L0 in unis of seconds/microns
    ## the sqrt(omega) comes from the fact that in the integral I am using dimensionless variables instead of k_par, I use k_par/k --> part with omega is inside function
    ## because I am integrating over energy
    
    ############ integration over k_parallel ##################
    k_par1 = lambda eV: 1e-4*omegac(eV) ## integration from 0
    
    # limit2 = 1.3*(1/beta)   ## already zero for this upper limit
    # limit2 = np.real(np.sqrt(epsi2)) ## inside light cone
 
    k_par2 = lambda eV: 30*omegac(eV)
    ###########################################################
    ############ integration over energy ######################
    eV1, eV2 = 0.001, upper_eV_limit
    
    Integral = dblquad(final_function_re, eV1, eV2, k_par1, k_par2)[0]
    
    return  Integral*factor_Gamma_norm_L0 



#EELS integrated over y-coordinate and over frequency
def EELS_double_integral_as_sum(upper_eV_limit,b,d,beta):
    """    
    Parameters
    ----------
    upper_eV_limit : upper limit of the integral over energy in eV
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
    delta_hbw = 0.005
    list_eV = np.arange(0.01, upper_eV_limit + delta_hbw ,delta_hbw)
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
            epsi2 =  epsilon2(eV,material)
            
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

print('1-Plot the EELS integrated over k_par, over the trajectory, and over energy as a function of the upper limit of integration, for different b')

labelx = r'Upper integration limit $\hbar\omega_{\text{f}}$ (eV)'
labelx = r'Cutoff energy $\hbar\omega_f$ (eV)'
labely = r'$\Gamma_{\parallel}/L_0$ (1/$\mu$m)'

title = r'EELS for $h = %.1f$ $\mu$m, $\epsilon_2 = \epsilon_{%s}(\omega)$, $v = %.2fc$' %(d,material,beta)

if create_data == 1:
    
    list_EELS_re_tot = []
    list_EELS_im_tot = []
    for b_nm in list_b_nm:
        b = b_nm*1e-3
        list_EELS_re = []
        list_EELS_im = []
        k = 0
        for eV in list_upper_eV_limit:  
            
            value,Nvalue = EELS_double_integral_as_sum(eV,b,d,beta)
            print(k, value,Nvalue)
            # value = Fresnel_coefficient(omegac,u,d,mode,Im_epsi2)
            list_EELS_re.append(np.real(value))
            list_EELS_im.append(np.imag(value))
            k = k + 1
            
        list_EELS_re_tot.append(list_EELS_re)
        list_EELS_im_tot.append(list_EELS_im)
        
    os.chdir(path_save)
#    header = title1 + r', $\beta$ = %.2f. Re(EELS) from paper 228 Eq. 3 divided by L/c and in Gaussian units (dimensionless) ' %(beta)
    header = title + r', $\beta$ = %.2f. Re(EELS) from paper 149 Eq. 25 divided by L0 in Gaussian units (1/microns), integrated over k_par, trajectory, and eV' %(beta)
    for j in range(len(list_b_nm)):
        b_nm = list_b_nm[j]
        
        np.savetxt('list_EELS_over_L0_%s'%(material)  + total_label + '_b%inm.txt' %(b_nm) ,list_EELS_re_tot[j], fmt='%.10f', delimiter='\t', header = header, encoding=None)
    np.savetxt('list_eV_upper_limit_%s'%(material)  + total_label+ '.txt' ,list_upper_eV_limit, fmt='%.10f', delimiter='\t', header = header, encoding=None)

else:
    list_EELS_re_tot = []
    os.chdir(path_save)
    for j in range(len(list_b_nm)):
        b_nm = list_b_nm[j]    
        listy = np.loadtxt('list_EELS_over_L0_%s'%(material)  + total_label + '_b%inm.txt' %(b_nm) , delimiter='\t', skiprows = 1)
        
    list_EELS_re_tot.append(listy)
    list_energy_eV = np.loadtxt('list_eV_upper_limit_%s'%(material)  + total_label+ '.txt' , delimiter='\t', skiprows = 1)

#%% 

plt.figure(figsize=tamfig)
plt.title(title,fontsize=tamtitle)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
for j in range(len(list_b_nm)):
    plt.plot(list_upper_eV_limit, np.array(list_EELS_re_tot[j]) ,'.-',lw = 1.5,label = r'$b = %i$ nm' %(list_b_nm[j]) )
 
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1) 
label_figure = 'EELS_int_energy_' + total_label
# plt.xscale('log')
os.chdir(path_save)
plt.savefig(label_figure + '.png', format='png',bbox_inches='tight',pad_inches = 0.04, dpi=dpi)  
plt.show() 



