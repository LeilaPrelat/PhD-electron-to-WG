
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
    
fit the EELS with an exponential 

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

path_data2 = os.path.join(path_basic, 'datfiles_EELS/along_z')

#%%

tamfig = [4.5,3.5]
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

W=500 # nm
h=200 ## nm

listL = [700,1000,1500]
N0=500 ### plot vs L for N = N0
L0=1000 ## plot vs N for L = L0 for mode = mode0 because mode0 has negative EELS
mode0 = 2 

list_modes = [1,2,3]
energy=1

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

print('Plot EELS vs L (length of the virtual geometry) with N=%i for each mode'%(N0))

kk = 0
for mode in list_modes:    

    plt.figure(figsize=tamfig)
    plt.title(title,fontsize=tamtitle)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    LL = 0
    for L in listL:
        name3 = 'EELS_N%i_W%inm_h%inm_L%inm_Ee%ikeV.dat' %(N0,W,h,L,Ee_electron_keV) ## virtual geometry
        name3 = 'EELS_comsol2_spectrum_N%i_W%inm_h%inm_L%inm_Ee%ikeV_zoom.dat' %(N0,W,h,L,Ee_electron_keV) ## virtual geometry
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
            # else:
            #     listeV_correct.append(listeV[k])
            #     listEELS_correct.append(-listEELS[k])
        
        peaks, _ = find_peaks(listEELS , height=np.max(listEELS)*0.1)
        listeV  = np.array(listeV )
        listEELS  = np.array(listEELS )
        
        sorted_indices = np.argsort(listeV)
        listeV = listeV[sorted_indices]
        listEELS = listEELS[sorted_indices]
                      
        plt.plot(listeV_correct, listEELS_correct,'.-',lw=lw,color = list_colors[LL], label = r'$L$ = %.1f $\mu$m' %(L*1e-3))
        # plt.plot(listeV[peaks], listEELS[peaks], "x",color = list_colors[LL])
    
        ## IMPORTANT: sort the y-values from minimum to maximum --> mode = -1 is the highest (last one), mode = -2 is the previous one, and so on, .. 
        sorted_index = np.argsort( listEELS[peaks])
        # listy_peaks_sorted = np.sort(listy_peaks) 
        # listx_peaks_sorted = list_energy0[peaks[sorted_index]]
        list_index_peaks = peaks[sorted_index]
        listx_peaks_sorted = listeV[peaks[sorted_index]]
        
        if L == 700: 
            try: 
               x_left, x_right, x0_fit, function_lorentzian, fit_success = find_width_of_peak(listeV_correct,listEELS_correct,-mode,title,0)
               print("mode = %i, L = %i nm, xleft = %.3f eV, xright = %.3f eV, x_res = %.3f, N = %i, W = %i nm, h = %i nm" %(mode,L,x_left,x_right,x0_fit,N0,W,h))
            except TypeError: 
                continue
    
        if width_height == [W,h]:
            plt.xlim([x_left_limits[kk],x_right_limits[kk]])
        if mode == 1: 
            plt.xticks(xticks1)
        elif mode == 2: 
            plt.xticks(xticks2)
        else:
            plt.xticks(xticks3)
        
        if mode == 1:
            plt.plot(listeV[peaks[0]], listEELS[peaks[0]], 'o', color = 'black')
            
    
        LL = LL + 1

    kk = kk + 1
    plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
    plt.ylim(-1e-4*1.3,np.max(listEELS)*1.1)
# plt.yscale('log')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1)
    plt.tight_layout()
    plt.savefig('EELS_convergency_virtual_geometry_mode%i_N%i.png'%(mode,N0),bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)
    plt.savefig('EELS_convergency_virtual_geometry_mode%i_N%i.pdf'%(mode,N0),bbox_inches='tight',pad_inches = 0.01, format='pdf', dpi=dpi)
    plt.show()

#%%

print('Plot EELS vs L (length of the virtual geometry) with N=%i for all modes'%(N0))
os.chdir(path_data)

plt.figure(figsize=tamfig)
plt.title(title,fontsize=tamtitle)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
LL = 0
for L in listL:
    name3 = 'EELS_N%i_W%inm_h%inm_L%inm_Ee%ikeV.dat' %(N0,W,h,L,Ee_electron_keV) ## virtual geometry
    name3 = 'EELS_comsol2_spectrum_N%i_W%inm_h%inm_L%inm_Ee%ikeV_zoom.dat' %(N0,W,h,L,Ee_electron_keV) ## virtual geometry
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
                  
    plt.plot(listeV, listEELS,'.-',lw=lw,color = list_colors[LL], label = r'$L$ = %.1f $\mu$m' %(L*1e-3))
    # plt.plot(listeV[peaks], listEELS[peaks], "x",color = list_colors[LL])

    ## IMPORTANT: sort the y-values from minimum to maximum --> mode = -1 is the highest (last one), mode = -2 is the previous one, and so on, .. 
    sorted_index = np.argsort( listEELS[peaks])
    # listy_peaks_sorted = np.sort(listy_peaks) 
    # listx_peaks_sorted = list_energy0[peaks[sorted_index]]
    list_index_peaks = peaks[sorted_index]
    listx_peaks_sorted = listeV[peaks[sorted_index]]
    
    ### study of a negative EELS as a function of L, N in plot_EELS_correct_negative.py 
    EELS_neg_energy = 0.35
    index_EELS_neg = int(np.where(listeV == 0.35)[0])
    EELS_neg_value = listEELS[index_EELS_neg]
   
    
    if L == 700: 
        plt.plot([0.35], [EELS_neg_value],'x',color = 'black')
        try: 
           x_left, x_right, x0_fit, function_lorentzian, fit_success = find_width_of_peak(listeV_correct,listEELS_correct,-mode,title,0)
           print("mode = %i, L = %i nm, xleft = %.3f eV, xright = %.3f eV, x_res = %.3f, N = %i, W = %i nm, h = %i nm" %(mode,L,x_left,x_right,x0_fit,N0,W,h))
        except TypeError: 
            continue

 
    LL = LL + 1

kk = kk + 1
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.ylim(-1e-4*1.3,np.max(listEELS)*1.1)
# plt.yscale('log')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1)
plt.tight_layout()
plt.savefig('EELS_convergency_virtual_geometry_N%i.png'%(N0),bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)
plt.savefig('EELS_convergency_virtual_geometry_N%i.pdf'%(N0),bbox_inches='tight',pad_inches = 0.01, format='pdf', dpi=dpi)
plt.show()


#%%

## mode = 2 has negative parts 
print('Plot EELS vs N (number of points) with L=%i for mode=%i'%(L0,mode0))

listN = [500,700]
mode=mode0
L=L0

plt.figure(figsize=tamfig)
plt.title(title,fontsize=tamtitle)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
LL = 0
for N in listN:
    name3 = 'EELS_N%i_W%inm_h%inm_L%inm_Ee%ikeV.dat' %(N0,W,h,L,Ee_electron_keV) ## virtual geometry
    name3 = 'EELS_comsol2_spectrum_N%i_W%inm_h%inm_L%inm_Ee%ikeV_mode%i.dat' %(N,W,h,L0,Ee_electron_keV,mode0) ## virtual geometry
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
                  
    plt.plot(listeV, listEELS,'.-',lw=lw,color = list_colors[LL], label = r'$N$ = %i' %(N))
    # plt.plot(listeV[peaks], listEELS[peaks], "x",color = list_colors[LL])

    ## IMPORTANT: sort the y-values from minimum to maximum --> mode = -1 is the highest (last one), mode = -2 is the previous one, and so on, .. 
    sorted_index = np.argsort( listEELS[peaks])
    # listy_peaks_sorted = np.sort(listy_peaks) 
    # listx_peaks_sorted = list_energy0[peaks[sorted_index]]
    list_index_peaks = peaks[sorted_index]
    listx_peaks_sorted = listeV[peaks[sorted_index]]
    
    plt.xticks(xticks2)
        
    LL = LL + 1

plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
# plt.yscale('log')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1)
plt.tight_layout()
plt.savefig('EELS_convergency_virtual_geometry_mode%i_L%i.png'%(mode0,L0),bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)
plt.savefig('EELS_convergency_virtual_geometry_mode%i_L%i.pdf'%(mode0,L0),bbox_inches='tight',pad_inches = 0.01, format='pdf', dpi=dpi)
plt.show()


#%% 

os.chdir(path_data2)

print('Plot EELS vs L (length of the virtual geometry) with N=%i as a function of z, fitted with a exponential' %(N0))


from scipy.optimize import curve_fit

def fit_exponential(z, amp, k_par_times_w):
    z_norm_w = z/W
    fit = amp*np.exp(-2*k_par_times_w*(z_norm_w-h/W)) 
    return fit 


ind_L = 1
L2 = listL[ind_L] ## convergency
energy_mode = [0.78,1] # see modes_data.txt
energy, mode = energy_mode
labelx2 = r'$z$ (nm)' 


name3 = 'EELS_along_z_N%i_W%inm_h%inm_L%inm_Ee%ikeV_xe0.00_energy%.4f.dat' %(N0,W,h,L2,Ee_electron_keV,energy) ## virtual geometry
#name3 = 'EELS_comsol2_spectrum_N%i_W%inm_h%inm_L%inm_Ee%ikeV_mode%i.dat' %(N0,W,h,L,Ee_electron_keV,mode) ## virtual geometry

#name3 = 'EELS_N%i_W%inm_h%inm_L%inm_R%inm_Ee%ikeV.dat' %(N,W,h,L,R,Ee_electron_keV) ## virtual geometry

try:
    tabla1 = np.loadtxt(name3, delimiter=' ',dtype=None)
    tabla1_2 = np.transpose(tabla1)
except (ValueError, FileNotFoundError):
    tabla1 = np.loadtxt(name3, delimiter='\t',dtype=None)
    tabla1_2 = np.transpose(tabla1)
  
listeV = tabla1_2[0]
list_z = tabla1_2[3]
listEELS = tabla1_2[4]

x_data = list_z[1:-1]
y_data = listEELS[1:-1]

init_conditions = (listEELS[1], 0.001)
# init_conditions = (max(y_data)-min(y_data), 2, min(y_data))
# ------------------------
# Fit curve
# ------------------------
popt, pcov = curve_fit(fit_exponential, x_data, y_data, p0=init_conditions)
A_fit, k_par_times_W_fit = popt
x_fit = np.linspace(min(x_data), max(x_data), 1000)
y_fit = fit_exponential(x_fit, *popt)

info = 'w%inm_h%inm_N%i_Ee%ikeV_mode%i_energy%.4feV' %(W,h,N,Ee_electron_keV,mode,energy)


############### interpolation of P(omega) after integration over z ###############
EELS_vs_z = interp1d(list_z, listEELS)
            

title = r'$W$ = %i nm, $h$ = %i nm, L = %i nm, $E_{\rm e}$ = %i keV, Im($\epsilon$) = 0.2' %(W,h,L2,Ee_electron_keV)
plt.figure(figsize=tamfig)
plt.title(title,fontsize=tamtitle)
plt.xlabel(labelx2,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
plt.plot(list_z[1:-1], listEELS[1:-1],'.-',lw=lw,color = list_colors[ind_L], label = r'$\hbar\omega$ = %.2f eV' %(energy))
plt.plot(x_fit, y_fit, 'r-', label="Exponential fit")
plt.plot([h+ze],[EELS_vs_z(h+ze)],'o',color = 'black')
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1)
plt.tight_layout()
plt.savefig('EELS_vs_z_mode%i_N%i.png'%(mode,N0),bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)
plt.savefig('EELS_vs_z_mode%i_N%i.pdf'%(mode,N0),bbox_inches='tight',pad_inches = 0.01, format='pdf', dpi=dpi)
plt.show()

#%% 


print('Plot the parameters of the fit as as function of energy')

bmin = 50 ## nm
zmin = bmin + h

step2=0.001 #integration over freq
z_integrate = np.arange(x_left, x_right + step2, step2)





list_energy = np.arange(0.7,0.84 + 0.001,0.001)
list_integration_over_z = []
list_A_fit = []
list_k_par_times_W_fit = []
for energy in list_energy:
    
    name3 = 'EELS_along_z_N%i_W%inm_h%inm_L%inm_Ee%ikeV_xe0.00_energy%.4f.dat' %(N0,W,h,L2,Ee_electron_keV,energy) ## virtual geometry
    #name3 = 'EELS_comsol2_spectrum_N%i_W%inm_h%inm_L%inm_Ee%ikeV_mode%i.dat' %(N0,W,h,L,Ee_electron_keV,mode) ## virtual geometry

    #name3 = 'EELS_N%i_W%inm_h%inm_L%inm_R%inm_Ee%ikeV.dat' %(N,W,h,L,R,Ee_electron_keV) ## virtual geometry

    try:
        tabla1 = np.loadtxt(name3, delimiter=' ',dtype=None)
        tabla1_2 = np.transpose(tabla1)
    except (ValueError, FileNotFoundError):
        tabla1 = np.loadtxt(name3, delimiter='\t',dtype=None)
        tabla1_2 = np.transpose(tabla1)
      
    listeV = tabla1_2[0]
    list_z = tabla1_2[3]
    listEELS = tabla1_2[4]

    x_data = list_z[1:-1]
    y_data = listEELS[1:-1]

    init_conditions = (listEELS[1], 0.001)
    # init_conditions = (max(y_data)-min(y_data), 2, min(y_data))
    # ------------------------
    # Fit curve
    # ------------------------
    popt, pcov = curve_fit(fit_exponential, x_data, y_data, p0=init_conditions)
    A_fit, k_par_times_W_fit = popt
 
    
    
    
    
    function_lorentzian = np.sum(fit_exponential(z_integrate, A_fit, k_par_times_W_fit))*step2
    
    
    list_k_par_times_W_fit.append(k_par_times_W_fit)
    list_A_fit.append(A_fit)
    print(pcov)
    
    
    
print("Save fit parameters as a function of energy")
name_file = "fitted_parameters_%s_bmin%inm.txt" %(info,bmin)
header = 'energy (eV)     A_fit     kappa*w     int_over_z'
table = np.column_stack((list_energy, list_A_fit, list_k_par_times_W_fit))
# Save parameters
np.savetxt(name_file, table, fmt='%.10f', delimiter='\t', header = header, encoding=None)

#%%
title = r'$W$ = %i nm, $h$ = %i nm, L = %i nm, $E_{\rm e}$ = %i keV, Im($\epsilon$) = 0.2' %(W,h,L2,Ee_electron_keV)
plt.figure(figsize=[4,3.5])
plt.title(title,fontsize=tamtitle)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.ylabel('Exponential fit',fontsize=tamletra,labelpad =labelpady)
plt.plot(list_energy, list_A_fit, label=r"$A$ from fit")
plt.plot(list_energy, list_k_par_times_W_fit, label=r"$\kappa\cdot W$ from fit")
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in",which = 'both', pad = pad)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=1)
plt.tight_layout()
plt.xticks(xticks1)
plt.yscale("log")
plt.savefig('fit_parameters_vs_energy_mode%i_N%i.png'%(mode,N0),bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)
plt.savefig('fit_parameters_vs_z_mode%i_N%i.pdf'%(mode,N0),bbox_inches='tight',pad_inches = 0.01, format='pdf', dpi=dpi)
plt.show()






